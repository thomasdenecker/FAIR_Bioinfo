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

# Get accession number
SRA_Accession="$(cut -f1 conditions.txt | tail -n +2)"

cd Project/samples

for SRA in ${SRA_Accession}
do
    echo "=============================================================="
    echo ${SRA}
    echo "=============================================================="
    fastq-dump ${SRA}
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
