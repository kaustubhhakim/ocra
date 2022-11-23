#!/usr/bin/env python
# coding: utf-8

# # Example plots using OCRA
# ## Hakim et al. (2023) ApJL
# 
# This Python script implements Reaktoro software to calculate ocean chemistry

from ocra import *

# Plot ocean pH as a function of PCO2
plot_pH(DIV='Ca')

# Plot Ca-CCD as a function of PCO2 and T
CaCCD_PCO2_T(nSiO2 = 1, totnum = 10, numQ1 = 20, numQ2 = 20)

# Plot stable phases as a function of PCO2
phases_PCO2(DIV='Ca', nSiO2 = 1, Temp = 310)
