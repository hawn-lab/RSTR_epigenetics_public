#!/bin/bash

##### Installation and setup to run ATAC_pipeline.sh ######
##### Linux or Mac only #####

########################################
## Update OS
########################################
sudo apt update --yes
sudo apt upgrade --yes

########################################
## Python v3 or greater 
########################################
sudo apt-get install python3

########################################
## Conda
########################################
# Install anaconda
sudo mkdir ~/tmp
cd ~/tmp
#Change to correct URL if not using Linux 64-bit
sudo curl -O https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh
sudo bash Anaconda3-2019.03-Linux-x86_64.sh
### Save
source ~/.bashrc

########################################
## git
########################################
sudo apt install git

########################################
## Conda packages
########################################
# Install packages in conda
### Configure channel priority
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

### Install packages
conda install -c conda-forge fastqc bowtie2 samtools
conda install adapterremoval genrich bedtools
conda update --all

########################################
## Git packages
########################################
# Install packages in git
sudo apt install git-all
mkdir ~/applications
cd ~/applications
git clone https://github.com/Ensembl/ensembl-vep.git

################# END ##################