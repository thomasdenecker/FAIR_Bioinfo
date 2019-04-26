# Welcome !

**Bienvenue à FAIR_bioinfo**

Vous trouverez ici des communications réalisées lors des sessions FAIR_bioinfo. Les communications sont en français. Tout le contenu présenté existe déjà en anglais sur internet. Nous proposons donc ici des ressources pour les francophones.

*You will find here some communications made during the I2BC Bioinformatics Club. Communications will be mainly in French. All the content presented also exists in English on the Internet. Therefore, we propose here resources for Francophones.*

**Informations pratiques**
- Quand ? : le dernier vendredi après midi de chaque mois (sauf juillet à définir), rdv 12h30
- Durée ? : 1h30 (questions incluses)
- Lieu ? : Salle de conférence A.Kalogeropoulos, b. 400, campus Orsay

**Objectifs**

L'objectif est de proposer et d'utiliser un panel d'outils permettant la réalisation d'un projet complet de bio-info en partant de rien et aboutissant à la création d'un conteneur (technologie Docker). Le partage, la valorisation et l'analyse dynamique des données seront inclus dans le panel.
FAIR correspond à l'acronyme anglais "Findable, Accessible, Interoperable, & Reusable", initialement défini pour les données mais que nous détournons ici pour leurs protocoles d'analyse.
Le projet support est une étude "d'expression différentielle de gènes" à partir de données RNAseq d'O.tauri.

**Pré-requis**

Quasi-rien ... Savoir taper sur un clavier ?

**Contact**

- Thomas DENECKER (<thomas.denecker@gmail.com>)
- Claire Toffano-Nioche (<claire.toffano-nioche@u-psud.fr>)

# Communications orales

## La mémoire du code

**Nouveaux outils pour la session**
- Git
- Github

**Programme :**
- Présentation du dépot local Git et du dépôt distant Github
- Mise en pratique par un exemple simple

**Date :** 26/11/2018

**Orateur :** Thomas Denecker

**Slides :** [Git/Github](https://thomasdenecker.github.io/Club-Bioinfo/docs/git-github.html)

**Wiki  :** [Git](https://github.com/thomasdenecker/FAIR_Bioinfo/wiki/Git), [Github](https://github.com/thomasdenecker/FAIR_Bioinfo/wiki/Github) et [Markdown](https://github.com/thomasdenecker/FAIR_Bioinfo/wiki/Markdown)

## Ce n'est pas de la magie

**Nouveaux outils pour la session**
- Terminal

**Programme :**
- Présentation du projet et de l'application finale
- Présentation des données
- Présentation du workflow de l'analyse
- Ouvrir un terminal (Mac Os, Ubuntu ou Windows)
- Initialisation du projet sur Git/Github
- Premier script shell : création de l'arborescence du projet
- Amélioration du script shell : récupération des données

**Date :** 25/01/2019

**Orateur :** Thomas Denecker

**Slides :** [Session 2](https://thomasdenecker.github.io/FAIR_Bioinfo/docs/session2.html)

## Installer et jouer avec les outils d'analyse

**Nouveaux outils pour la session**
- conda
- quelques outils d'analyse pour des données RNAseq

**Programme :**
- Présentation des outils du workflow (FastQC, Bowtie2, Samtools, HTseq-Count)
- Installation de FastQC à la main
- Présentation de conda
- Script pour l'installation de FastQC avec conda
- Amélioration du script d'installation : ajout des autres outils (travail à la maison)
- Amélioration du script d'analyse : workflow d'analyse à faire en boucle pour chaque échantillon (à faire tourner à la maison)

**Date :** 22/02/2019

**Orateur :** Thomas Denecker & Claire Toffano

**Slides :** [Session 3](https://thomasdenecker.github.io/FAIR_Bioinfo/docs/session3.html)

## Une virée en mer

**Nouveaux outils pour la session**
- Docker
- Cloud IFB

**Programme :**
- Présentation de docker
- Ecriture du dockerfile
- Utilisation du docker avec le script d'analyse
- Partager le docker : Docker Hub
- Pourquoi ? utilisation dans un cloud (ex. IFB)

**Date :** 29/03/2019

**Orateur :** Thomas Denecker

**Slides :** [Session 4](https://thomasdenecker.github.io/FAIR_Bioinfo/docs/session4.html)

## I've got the power !

**Programme :**

_La parallélisation_
   - Snakemake
   - Comparaison avec le script.sh

_Le cluster de calcul_

   - C'est quoi?
   - Cluster de l'IFB / I2BC
   - Singularity ↔ Docker
   - Exemple sur l'IFB et l'I2BC

**Date :** 26/04/2019

**Orateur :** Thomas Denecker

**Slides :** [Session 5](https://thomasdenecker.github.io/FAIR_Bioinfo/docs/session5.html)
