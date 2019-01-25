##------------------------------------------------------------------------------
## FAIR script
## Author: T. Denecker & C. Toffano-Nioche
## Affiliation: I2BC
## Aim: A workflow to process RNA-Seq.
## Organism: O. tauri
## Date: Jan 2019
## Step :
## 1- Create tree structure
## 2- Download data from SRA
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

cd Project/samples

IFS=$'\n'       # make newlines the only separator
for j in $(tail -n +2 ../../conditions.txt)
do
    
    # Get important information from the line
    access=$( echo "${j}" | cut -f6 )
    id=$( echo "${j}" | cut -f1 )
    md5=$( echo "${j}" | cut -f7 )

    echo "--------------------------------------------------------------"
    echo ${id}
    echo "--------------------------------------------------------------"

    # Download file
    wget ${access} # wget method

    # Get md5 of downloaded file
    md5_local="$(md5sum ${id}.fastq.gz | cut -d' ' -f1)"
    echo ${md5_local}
    
    # Test md5 
    if [ "${md5_local}" == "${md5}" ]
    then
        echo "Done"
    else
        echo "Nope"
        exit 1
    fi
done

cd ../..
exit 0
