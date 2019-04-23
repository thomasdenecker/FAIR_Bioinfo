##------------------------------------------------------------------------------
## FAIR script
## Author: T. Denecker & C. Toffano-Nioche
## Affiliation: I2BC
## Aim: A workflow to process RNA-Seq.
## Organism : O. tauri
## Date: Jan 2019
## Step :
## 1- Create tree structure
## 2- Download data from SRA
## 3- Run workflow (snakemake)
##------------------------------------------------------------------------------


echo "=============================================================="
echo "Creation of tree structure"
echo "=============================================================="

mkdir Project
mkdir Project/samples
mkdir Project/annotations
mkdir Project/bowtie2
mkdir Project/fastqc
mkdir Project/genome
mkdir Project/graphics
mkdir Project/htseq
mkdir Project/reference
mkdir Project/samtools

echo "=============================================================="
echo "Download data from SRA"
echo "=============================================================="

# Get accession number (comment / uncomment to change methods)
 Accession="$(cut -f7 conditions.txt | tail -n +2)" # ascp method
# Accession="$(cut -f1 conditions.txt | tail -n +2)" # fastq dump method
# Accession="$(cut -f6 conditions.txt | tail -n +2)" # wget method

cd Project/samples

IFS=$'\n'       # make newlines the only separator
for j in $(tail -n +2 ../../conditions.txt)
do

    access=$( echo "$j" |cut -f7 )
    id=$( echo "$j" |cut -f1 )

    echo "--------------------------------------------------------------"
    echo ${id}
    echo "--------------------------------------------------------------"

    ascp -QT --file-checksum=md5 --file-manifest=text --file-manifest-path=. -l 300m -P33001 -i /home/miniconda3/etc/asperaweb_id_dsa.openssh ${access} .

    md5_local="$(md5sum $id.fastq.gz | cut -d' ' -f1)"
    echo $md5_local

    if grep -q $md5_local *.txt
    then
        echo "Done"
    else
        rm $id.fastq.gz
        access=$( echo "$j" |cut -f6 )
        wget ${access} # wget method
    fi
    rm *.txt
done

cd ../..

echo "=============================================================="
echo "Download annotations"
echo "=============================================================="

wget https://raw.githubusercontent.com/thomasdenecker/FAIR_Bioinfo/master/Data/O.tauri_annotation.gff -P Project/annotations

echo "=============================================================="
echo "Download genome"
echo "=============================================================="

wget https://raw.githubusercontent.com/thomasdenecker/FAIR_Bioinfo/master/Data/O.tauri_genome.fna -P Project/genome

echo "=============================================================="
echo "Snakemake"
echo "=============================================================="

snakemake

echo "=============================================================="
echo "Create unique count table"
echo "=============================================================="

Rscript R-code/countTable.R

echo "=============================================================="
echo "Shiny app"
echo "=============================================================="

R -e "shiny::runApp('R-code', port=4444)"
