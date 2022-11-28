# OCRA #
#### Author: Kaustubh Hakim ####

## 1. Overview ##

OCRA (Ocean Chemistry with Reaktoro And beyond) is a code to simulate water ocean chemistry. This code calculates ocean pH and the amounts of aqueous, gaseous and mineral phases with a special focus on pure Ca, Mg and Fe carbonate systems. This code makes use of Reaktoro v2 (https://reaktoro.org/intro.html) for modelling ocean chemistry.

## 2. Installation ##

OCRA is written in Python 3 and requires a Reaktoro conda environment: https://reaktoro.org/installation/installation-using-conda.html

If conda is not installed, please follow conda installation instructions: https://conda.io/projects/conda/en/stable/user-guide/install/index.html

Before using OCRA, please follow the steps below to install required libraries. 

Add conda-forge channel to find Reaktoro

> conda config --add channels conda-forge

Create new conda environment and install Reaktoro and its dependencies

> conda create --name rkt reaktoro

Activate conda environment

> conda activate rkt

Install OCRA dependencies that are not part of Reaktoro installation

> conda install pandas scipy astropy matplotlib

Further information on Reaktoro can be found here https://reaktoro.org/intro.html

Download OCRA directly or using git

> git clone https://github.com/kaustubhhakim/ocra.git

Go to ocra directory

> cd ocra/

## 3. Running OCRA ##

### Example ###

After the installation, you can test if OCRA is working by running the following command. This command generates and saves three figures and three csv tables in the ocra directory.

> python plots_example.py

### All figures ###

The following command generates and saves tables in CSV format and figures in PDF format. These figures are included in Hakim et al. (2023).

> python plots_paper.py

## 4. References ##

Hakim et al. (2023)
