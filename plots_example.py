#!/usr/bin/env python
# coding: utf-8

# # OCRA: Ocean Chemistry with Reacktoro And beyond
# ### This Python code implements Reaktoro software to calculate ocean chemistry
# 
# ## Reference: Hakim et al. (2023) ApJL
# 
# ### plots_example.py # contains functions to make Example plots to test OCRA

# Import ocra
from ocra import *

# Plot ocean pH as a function of PCO2
## DIV decides the carbonate system: Ca, Mg or Fe
pH_PCO2(DIV='Ca') 

# Plot Ca-CCD as a function of PCO2 and T
## nSiO2 = 1 includes silica and 0 excludes silica
## numQ1 is the number of steps in x-axis
## numQ2 is the number of steps in y-axis
## totnum is the number of steps in z-axis
CaCCD_PCO2_T(nSiO2 = 1, totnum = 10, numQ1 = 20, numQ2 = 20) 

# Plot stable phases as a function of PCO2
## DIV decides the carbonate system: Ca, Mg or Fe
## nSiO2 = 1 includes silica and 0 excludes silica
## Temp is temperature in kelvin
phases_PCO2(DIV='Ca', nSiO2 = 1, Temp = 310)
