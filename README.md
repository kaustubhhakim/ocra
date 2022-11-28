# OCRA #
#### Author: Kaustubh Hakim ####

# Overview #

OCRA (Ocean Chemistry with Reaktoro And beyond) is a code to simulate water ocean chemistry. This code calculates ocean pH and the amounts of aqueous, gaseous and mineral phases with a special focus on pure Ca, Mg and Fe carbonate systems. This code makes use of Reaktoro v2 for modelling ocean chemistry.

# Installation #

The code is written in Python. Before runnning OCRA on a personal computer, Python 3 needs to be installed. Additionally, the following libraries are required. Please follow the steps below. 

New conda environment

> conda create --name ocra reaktoro

Activate conda environment

> conda activate ocra

Add conda-forge channel to find Reaktoro

> conda config --add channels conda-forge

Reaktoro installation

> conda install reaktoro

Install dependencies that are not part of Reaktoro

> conda install pandas scipy astropy matplotlib

Further information on Reaktoro can be found here https://reaktoro.org/intro.html

Download ocra

> git clone https://github.com/kaustubhhakim/ocra.git

# Running OCRA #

## Example ##

After the installation, you can test if OCRA is working by running the following command. This command generates and saves three figures and three csv tables.

> python plots_example.py

## All figures ##

The following command generates and saves tables in CSV format and figures in PDF format. These figures are included in Hakim et al. (2023).

> python plots_paper.py

# References #

Hakim et al. (2023)
