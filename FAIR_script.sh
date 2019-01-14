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

# nom du fichier contenant le génome de référence
genome=GCF_000214015.3_version_140606_genomic.fna
# nom du fichier contenant les annotations
annotations=GCF_000214015.3_version_140606_genomic_DUO2.gff

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
