#!/bin/bash
# les parametres pour le planificateur slurm commencent par "#SBATCH" :
#SBATCH -o slurm.%N.%j.out  # fichier de redirection de l'output standard renommé avec l'identifiant du job
#SBATCH -e slurm.%N.%j.err  # fichier de redirection de l'erreu standard renommé avec l'identifiant du job
#SBATCH --partition fast    # choix des queues/noeuds de calcul ("fast" => le calcul sera coupé après n heures)
#SBATCH --mem 10GB          # réservation de RAM

# lancement du docker et du script FAIR_script.sh
singularity exec -B $PWD:/home/rstudio fair_bioinfo.simg bash ./FAIR_script.sh
