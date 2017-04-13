Table of Contents
=================
   * [Geno imputation](#geno-imputation)
   * [Pipeline](#pipeline)
      * [Setup and software dependencies](#setup-and-software-dependencies)
      * [Common folder tree](#common-folder-tree)
      * [Split collection files.](#split-collection-files)
      * [Automatically summarize raw data in folder tree.](#automatically-summarize-raw-data-in-folder-tree)
      * [Prepare marker map files.](#prepare-marker-map-files)
         * [Convert annotation file to map file accepted by ioSNP.py](#convert-annotation-file-to-map-file-accepted-by-iosnppy)
      * [Convert raw-data to Plink text input files](#convert-raw-data-to-plink-text-input-files)
         * [Affymetrix 55K](#affymetrix-55k)
         * [Illumina](#illumina)
      * [Convert Plink text input files to binary files](#convert-plink-text-input-files-to-binary-files)
      * [Merge converted plink files](#merge-converted-plink-files)
      * [QC of converted raw data <strong>before</strong> imputation.](#qc-of-converted-raw-data-before-imputation)
         * [Suggestions for filtering:](#suggestions-for-filtering)
      * [Checklist: Validate that Tim's <a href="https://en.wikipedia.org/wiki/Tacit_knowledge">tacit knowledge</a> is accounted for!](#checklist-validate-that-tims-tacit-knowledge-is-accounted-for)
      * [Convert to alphaimpute format](#convert-to-alphaimpute-format)
      * [Documentation and scripts for imputation.](#documentation-and-scripts-for-imputation)

Created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc)

# Geno imputation
Welcome to the Geno Imputation github repository. We want this to b a common code base for developing the genotype reference for Geno imputed with AlphaImpute. Our goal is to write scripts that both document and does the convertion from Chip raw-data and imputation. As long as the software dependencies are met and the setup code block is edited we should be able to run the scripts on the HPC clusters in Ås / Edinburg / Oslo and on Geno's production server.

# Pipeline
Link to other .md docs,, or have everything here with a big TOC at the top. 

## Setup and software dependencies
First set up some user-specific variables pointing to a clone of this repository, the [snptranslate](https://github.com/timknut/snptranslate) repository, raw genotype data from the ftp server ftpgeno.geno.no. The code below dependends on the Bash version 4, Python 2 with numpy (for snptranslate) and R with a few packages (data.table,knitr).

```bash
#user-specific variable
prefix=/mnt/users/gjuvslan/geno/geno_imputation          #git clone of https://github.com/argju/geno_imputation
ftpgeno=/mnt/users/gjuvslan/geno/geno_imputation/ftpgeno #raw genotype data from ftpgeno.geno.no:/avlgeno/Raw_Data_Files
snptranslatepath=/mnt/users/gjuvslan/geno/snptranslate/  #git clone of https://github.com/timknut/snptranslate.git
export PATH=$PATH:$snptranslatepath

#check required software
set -e
echo "#### Checking dependency Python 2 with Numpy." 
python -c "import sys; assert sys.version_info[0]==2; import numpy as np"
echo "#### Checking dependency Plink 1.9."
plink --version | grep v1.9
echo "#### Checking dependency R with packages data.table, knitr, ggplot2, DT."
R -e 'library(data.table); library(knitr); library(ggplot2); library(DT)'
set +e

#check raw genotype data
Nrawfiles=$(ls -1 $ftpgeno/Raw_Data_Files/FinalReport* $ftpgeno/Raw_Data_Files/*.calls.txt | wc -l)
if [ "$Nrawfiles" -ge "90" ]
then
    echo "Found all unzipped raw genotype files in $ftpgeno"
else
    echo "Did not find all raw genotype files in $ftpgeno."
    ls -1 $ftpgeno/Raw_Data_Files/FinalReport*.txt $ftpgeno/Raw_Data_Files/Batch*.txt
    echo "Run scripts/download_raw_data to get all the raw data."
    exit 1
fi
```

## Common folder tree

```bash
# Code to from the raw data from ftpgeno.geno.no:/avlgeno/Raw_Data_Files to the common folder tree
cd genotype_rawdata
mkdir -p affymetrix25k illumina54k_v1 illumina54k_v2 illumina54k_v2/collections illumina777k affymetrix54k
ln -s -t affymetrix25k $ftpgeno/Raw_Data_Files/FinalReport_25k.txt
ln -s -t illumina54k_v1 $ftpgeno/Raw_Data_Files/FinalReport_54kV1*
ln -s -t illumina54k_v2 $ftpgeno/Raw_Data_Files/FinalReport_54kV2*
ln -s -t illumina54k_v2/ $ftpgeno/Raw_Data_Files/Nordic_54k*
ln -s -t illumina54k_v1/ $ftpgeno/Raw_Data_Files/Swedish_54k_ed1.txt 
mv illumina54k_v2/FinalReport_54kV2_collection* illumina54k_v2/collections
ln -s -t illumina777k $ftpgeno/Raw_Data_Files/FinalReport_777k*
ln -s -t illumina777k/ $ftpgeno/Raw_Data_Files/Nordic_HDexchange_201110.txt
ln -s -t affymetrix54k/ $ftpgeno/Raw_Data_Files/*.calls.txt 
cd ..
```

## Automatically summarize raw data in folder tree. 

The bash code below extracts key information from the raw data files, and this is collected into a full report in [reports/genotypereport.Rmd](reports/genotypereport.Rmd) .

```bash
cd genotype_rawdata

#grep Illumina FinalReports for headers, normalize whitespace and create table
grep -A 8 -m 1 "\[Header" illumina*/*/Final* illumina*/Final* | grep -v -e "\[" -e "--" > tmp
sed -i -e s/Processing\\t/Processing" "/g -e s/Num\\t/Num" "/g -e s/Total\\t/Total" "/g tmp
sed -i -e s/GSGT\\t/GSGT" "/g -e s/2010\\t/2010" "/g -e s/54\\t/54" "/g tmp
sed -i -e s/bovinehd-manifest-b.bpm/bovinehd_manifest_b.bpm/g -e 's/\r$//' tmp
cat tmp | sed -e s/-/\\t/g -e s/\\t\\t/\\t/g > illumina_headers

## grep Illumina FinalReports for ids, normalize whitespace and create table of filenames and ids
# 1. files in GenomeMatrix format
matrixfiles=illumina54k_v2/collections/Final*" "illumina777k/FinalReport_777k_apr2015.txt" "illumina777k/FinalReport_777k_jun2015.txt" "illumina54k_v2/FinalReport_54kV2_nov2011_ed1.txt" "illumina777k/FinalReport_777k.txt
grep -A 1 -m 1 "\[Data" $matrixfiles | grep -v -e "\[Data" -e "--" > tmp 
sed -i -e s/FinalReport_777k.txt-2402/FinalReport_777k.txt-\\t2402/g tmp #FinalReport_777k.txt lacks tab before first ID
cat tmp | sed -e s/-//g -e s/[[:space:]]/\\t/g | awk '{for(i=2;i<=NF;i++) print $1"\t"$i}' > illumina_ids
rm illumina_formats
for file in $matrixfiles; do echo -e $file"\t"Genomematrix >> illumina_formats ; done

# 2. files in GenomeList format (~10 min)
listfiles=illumina54k_v1/Final*" "illumina54k_v2/FinalReport_54kV2_feb2011_ed1.txt" "illumina54k_v2/FinalReport_54kV2_genoskan.txt" "illumina54k_v2/FinalReport_54kV2_ed1.txt" "illumina777k/FinalReport_777k_jan2015.txt" "illumina54k_v1/Swedish_54k_ed1.txt" "illumina54k_v2/Nordic*.txt" "illumina777k/Nordic_HDexchange_201110.txt" "affymetrix25k/FinalReport_25k.txt
time for file in $listfiles; do grep -q $file illumina_ids_list && echo 'Skipping (already in illumina_ids_list)' $file || ./id.R $file >> illumina_ids_list;done
cat illumina_ids_list >> illumina_ids
for file in $listfiles; do echo -e $file"\t"Genomelist >> illumina_formats ; done

## Affymetrix reports
grep -E "time-str|chip-type|cel-count" affymetrix54k/*.calls.txt > tmp
sed -i -e s/:#%/\\t/g -e s/affymetrix-algorithm-param-apt-//g -e s/=/\\t/ tmp
sed -e s/time-str/Date/g -e s/opt-chip-type/Chip/g -e s/opt-cel-count/Nsamples/g tmp > affymetrix_headers
for file in $(ls affymetrix54k/*.calls.txt)
do 
  echo "Counting SNPs in $file"
  echo -n -e "$file\tNum SNPs\t$(grep AX-[[:digit:]]* $file | cut -f 1 | wc -l)\n" >> affymetrix_headers
done

grep probeset affymetrix54k/*.calls.txt | awk '{for (i=2; i<=NF; i++) print $1"\t"$i}' > tmp
sed -e s/:probeset_id//g tmp > affymetrix_ids

cd ..
```

## Prepare marker map files.  
ioSNP.py will create the plink map file, which is a better solution than doing it manually in R.

### Convert annotation file to map file accepted by ioSNP.py
Annotation files for 50Kv1&v2 and 777K Illumina chips was [downloaded from SNPchimp](http://bioinformatics.tecnoparco.org/SNPchimp/index.php/download/download-cow-data) and the tab-separated gzipped raw files are in [genotype_rawdata/marker_mapfiles/snpchimp/](genotype_rawdata/marker_mapfiles/snpchimp/). The SNP positions refer to a particular genom assembly and choosing "Native platform" will not give consistent SNP positions across chips, we should therefore use the "UMD3.1" positions when creating plink files, see [snpchimp.Rmd](genotype_rawdata/marker_mapfiles/snpchimp/snpchimp.Rmd) / [snpchimp.pdf](genotype_rawdata/marker_mapfiles/snpchimp/snpchimp.pdf) for details on the SNPchimp data.

```bash
cd genotype_rawdata/marker_mapfiles/
#marker map files for Illumina chips with Native platform positions
zgrep Bov_Illu50Kv1 snpchimp/illumina_50Kv1_50Kv2_777K_native.tsv.gz | gawk '{print $5"\t"$7"\t"0"\t"$6}' > illumina50Kv1_annotationfile_native.map
zgrep Bov_Illu50Kv2 snpchimp/illumina_50Kv1_50Kv2_777K_native.tsv.gz | gawk '{print $5"\t"$7"\t"0"\t"$6}' > illumina50Kv2_annotationfile_native.map
zgrep Bov_IlluHD snpchimp/illumina_50Kv1_50Kv2_777K_native.tsv.gz | gawk '{print $5"\t"$7"\t"0"\t"$6}' > illumina777K_annotationfile_native.map

#marker map files for Illumina chips with UMD3.1 dbSNP positions
zgrep Bov_Illu50Kv1 snpchimp/illumina_50Kv1_50Kv2_777K_UMD3.1.tsv.gz | gawk '{print $5"\t"$7"\t"0"\t"$6}' | sed s/^99/0/ > illumina50Kv1_annotationfile_umd3_1.map
zgrep Bov_Illu50Kv2 snpchimp/illumina_50Kv1_50Kv2_777K_UMD3.1.tsv.gz | gawk '{print $5"\t"$7"\t"0"\t"$6}' | sed s/^99/0/ > illumina50Kv2_annotationfile_umd3_1.map
zgrep Bov_IlluHD snpchimp/illumina_50Kv1_50Kv2_777K_UMD3.1.tsv.gz | gawk '{print $5"\t"$7"\t"0"\t"$6}' | sed s/^99/0/ > illumina777K_annotationfile_umd3_1.map

#marker map file for Affymetrix 50K chip
cut -f 1-4 affy50k_annotation_final_list_20160715.txt | grep -v NA > affymetrix50K.map

cd ../..
```


## Convert raw-data to Plink text input files

Convert the Affymetrix and Illumina genotype report files to Plink standard [text input files (.ped)](https://www.cog-genomics.org/plink2/formats#ped) files paired up with matching [.map](https://www.cog-genomics.org/plink2/formats#map) files.

### Affymetrix 55K
Use snptranslate-script from https://github.com/timknut/snptranslate/blob/master/seqreport_edit.py to convert from Affymetrix format to Geno format. Then use ioSNP.py to convert from Geno to Plink text input format.

```bash

#Convert Affymetrix 50K files
markermap_affy=genotype_rawdata/marker_mapfiles/affy50k_annotation_merged_20170119.txt
mkdir -p genotype_data/plink_txt genotype_data/plink_bin genotype_data/reports
for affycall in `gawk '{print $1}' genotype_rawdata/affymetrix_headers | sort | uniq`
do
    pedfile=genotype_data/plink_txt/$(basename $affycall).ped
    if [ -e $pedfile ]
    then
       echo "Skipping conversion of: "$affycall", plink input file exists: "$pedfile
    else
        echo -e "Converting Affymetrix report file: "$affycall"\tMarkermap: "$markermap_affy
        time seqreport_edit.py -m $markermap_affy -o genotype_data/plink_txt/$(basename $affycall) -r genotype_data/reports/$(basename $affycall).seqreport genotype_rawdata/$affycall
        ioSNP.py -i genotype_data/plink_txt/$(basename $affycall) -n Geno -o genotype_data/plink_txt/$(basename $affycall).ped -u Plink --output2 genotype_data/plink_txt/$(basename $affycall).map -m genotype_rawdata/marker_mapfiles/affymetrix50K.map
    fi
done
```
### Illumina

Use ioSNP.py to convert all Finalreport files.

```bash
snptranslatepath=/mnt/users/gjuvslan/geno/snptranslate/
export PATH=$PATH:$snptranslatepath

## Assign marker position maps for each chip
declare -A markermap=([illumina54k_v1]=illumina50Kv1_annotationfile_umd3_1.map [illumina54k_v2]=illumina50Kv2_annotationfile_umd3_1.map [illumina777k]=illumina777K_annotationfile_umd3_1.map)

## Loop over chips and convert all genotype report files to .ped format 
for chip in "${!markermap[@]}"
do
    for report in `grep $chip genotype_rawdata/illumina_formats | gawk '{print $1}'`
    do
        pedfile=genotype_data/plink_txt/$(basename $report).ped
        if [ -e $pedfile ]
        then
            echo "Skipping conversion of: "$report", plink input file exists: "$pedfile
        else
            echo -e "Converting Illumina report file: "$report"\tMarkermap: "${markermap[$chip]}  
            informat=$(grep $report[[:space:]] genotype_rawdata/illumina_formats | gawk '{print $2}')
            ioSNP.py -i genotype_rawdata/$report -n $informat -o $pedfile".tmp" -u Plink --output2 genotype_data/plink_txt/$(basename $report).map -m genotype_rawdata/marker_mapfiles/${markermap[$chip]} && mv $pedfile".tmp" $pedfile
        fi
    done
done
```

## Convert Plink text input files to binary files

```bash
# convert all .ped and .map files to Plink binary files
for chip in illumina54k_v1 illumina54k_v2 illumina777k affymetrix54k
do
    for file in `grep -h $chip genotype_rawdata/illumina_formats genotype_rawdata/affymetrix_headers | gawk '{print $1}' | sort | uniq`
    do
        plinkbin=genotype_data/plink_bin/$(basename $file)
        if [[ -e $plinkbin".bed" && -e $plinkbin".bim" && -e $plinkbin".fam" ]]
        then
	    echo "Skipping conversion of: "$(basename $file)", plink binaries exist: "$binfile".bim/.bed/.fam"
        else
            time plink --cow --file genotype_data/plink_txt/$(basename $file) --out genotype_data/plink_bin/$(basename $file)
        fi
    done 
done

#parse plink logs to extract table with number of SNPs and individuals
for file in genotype_data/plink_bin/*.log
do
    echo -n $(echo $file | sed s/.log//)
    grep single-pass $file | gawk -F "[( ,)]" '{print "\t"$6"\t"$9}'
done
```

## Update pedigree and generate reports in plink

* Update IDs and pedigree to GenoId, using mapping file ped_genoid_all.txt (separate repo: [geno_imputation_idmapping](git@bitbucket.org:genogit/geno_imputation_idmapping.git))
* Write reports on missingness, HW, MAF and heterozygosity

```bash
cd genotype_data
## per File id mapping (remove, update ids, update pedigree, update sex)
idmap=$ftpgeno/Id_Raw_Data_Files/ped_genoid_all.txt
mkdir -p plink_bin_updateid plink_bin_updateped updates plink_bin_reports
for chip in illumina54k_v1 illumina54k_v2 illumina777k affymetrix54k
do
    for file in $(grep $chip $idmap | cut -f 2 | uniq)
    do
        grep "$file.*Remove" $idmap | gawk '{print "F0\t"$3}' | sort |uniq > updates/$file.remove 
        grep "$file.*Impute" $idmap | gawk -v F="$file" '{print "F0\t"$3"\t"F"\t"$4}' | sort | uniq > updates/$file.ids
        plink --cow --bfile plink_bin/$file --remove updates/$file.remove --update-ids updates/$file.ids --make-bed --out plink_bin_updateid/$file
        grep "$file.*Impute" $idmap | gawk -v F="$file" '{print F"\t"$4"\t"$5"\t"$6}' | sort | uniq > updates/$file.parents
        grep "$file.*Impute" $idmap | gawk -v F="$file" '{print F"\t"$4"\t"$8}' | sort | uniq > updates/$file.sex
        plink --cow --bfile plink_bin_updateid/$file --update-parents updates/$file.parents --update-sex updates/$file.sex --make-bed --out plink_bin_updateped/$file
	plink --cow --bfile plink_bin_updateped/$file --hardy --missing --het --freqx --nonfounders --out plink_bin_reports/$file
    done
done

##check for errors and warnings on id updates
grep -i Error plink_bin_update*/*.log
grep -i Warning plink_bin_update*/*.log
grep -i Note plink_bin_update*/*.log

cd ..
```

## Merge converted plink files into one file per chip
See the [plink documentation](https://www.cog-genomics.org/plink2/data#merge) for more information.
Given a file that list file prefixes that are to be merged with a plink binary file-set, the following command would work.
```bash
cd genotype_data
mkdir -p plink_merged_chip
for chip in affymetrix25k illumina54k_v1 illumina54k_v2 illumina777k affymetrix54k
do
    grep $chip ../genotype_rawdata/illumina_formats | cut -f 1 | sed s/$chip/plink_bin/g | sed s/collections//g > $chip.files
    grep $chip ../genotype_rawdata/affymetrix_headers | cut -f 1 | sort | uniq | sed s/$chip/plink_bin/g >> $chip.files
    tail -n +1 $chip.files > $chip.merge
    plink --cow --bfile $(head -1 $chip.files) --merge-list $chip.merge --out plink_merged_chip/$chip
done
```

## QC of converted raw data **before** imputation. 
### Suggestions for filtering:
* Missingess per animal. 90 % 
* Missingness per SNP 90 %
* HWE p < 1e-7
* MAF < 0.01
* Mendelian error filtering per SNP and animal. (Although AlphaImpute do a good job at this.)
* Heterozygosity per animal
Example plink command: `plink --bfile infile --cow --hwe 1e-7 --maf 0.01 --geno 0.90 --mind 0.90`

## Checklist: Validate that Tim's [tacit knowledge](https://en.wikipedia.org/wiki/Tacit_knowledge) is accounted for!
- [ ] **Flipping**: Illumina 50K and 777k should be flipped to the forward strand. She lists in chip-folders on the Dropbox Geno_sharing. Flipping is done by Paolo with `plink --flip`
- [ ] **Remap 25k markers on affy 50k** Convertion list `final_corrected_Affy25K_on_50K.txt` in Dropbox.
- [ ] **Correct subsetting of Affy 50k markers.** See list in Affy50k Dropbox folder. `final_list_to_Paolo_may_2016_affy50k_markers2_keep`
- [ ] **Markers wrongly positioned because of assembly problems** See `assembly_fail_markers_to_exclude_*` in affy50k and illumina 777k folders in Dropbox.
- [ ] **Have Paolo recieved and included the latest Affymetrix data?**
- [ ] **Keep most recent or higher density sample** when same sample is genotyped on same or several platforms. 


## Convert to alphaimpute format

See Paolos [script](https://github.com/timknut/geno_imputation/blob/master/scripts/plink2_alphaimpute.Rmd). In addition to this, do a chromsomoe loop in plink and do something like `for i in seq 1 29; do plink --recode A --chr $i --cow --bfile prefix --out prefix.$i.raw; done`

## Documentation and scripts for imputation.
Paolo: TODO
