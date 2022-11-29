#!/usr/bin/env python
# coding: utf-8

# # OCRA: Ocean Chemistry with Reacktoro And beyond
# ### This Python code implements Reaktoro software to calculate ocean chemistry
# 
# ## Reference: Hakim et al. (2023) ApJL
# 
# ### store.py # contains analytical functions and dictionary objects to store data


# Import libraries

import numpy as np
from astropy.constants import R


# Set global constants

numden = 1000 # 1 mol/kg = 1000 mol/m3

totH2O = 55.5 # moles
totN2 = 0.1 # moles --- produces pressure equivalent to 0.8 bar

low_cutoff = 1e-10 # mol/m3


# Analytical functions and scaling relations

def ocean_depth(totP = 1): # Leroy and Parthiot (1998) 
    '''
    Returns ocean depth in km as a function of ocean pressure in bar
    '''
    num = 97.266 * totP - 2.512e-3 * totP**2 + 2.28e-7 * totP**3 - 1.8e-11 * totP**4
    den = 9.7803 + 1.1e-5 * totP
    return num/den / 1000 # m to km

def ocean_pH_upp(PCO2, nDIV, logK9, logK16):
    '''
    Returns analytical upper limit of ocean pH in the presence of carbonates
    '''
    pH = -0.5 * (np.log10(PCO2) + logK16 + logK9 + np.log10(nDIV/numden))
    return pH

def ocean_pH_low(PCO2, logK3):
    '''
    Returns analytical lower limit of ocean pH in the absence of carbonates
    '''
    pH = -0.5 * (np.log10(PCO2) + logK3)
    return pH

def weath_scaling(PCO2, T, beta = 0.3, Ea = 31e3, delT = 13.7):
    '''
    Returns scaling in the divalent cation number density due to weathering
    '''
    PCO20 = 0.3e-3 # bar
    T0 = 288 # k
    
    if beta == 0:
        scaling = 1
    elif beta == 0.3:
        scaling = (PCO2 / PCO20)**beta * np.exp((T-T0)/delT) 
    else:
        scaling = (PCO2 / PCO20)**beta * np.exp(-Ea/R.value * (1/T - 1/T0)) 
    
    return scaling


# Custom dictionary objects to save chemical species as number density [dm^-3]

def chem_dict1(totnum):
    '''
    Returns 1D Chemical Dictionary Object to be used for Plotting
    '''
    chem = {'Ca+2': np.zeros(totnum), 'Mg+2': np.zeros(totnum), 
            'Fe+2': np.zeros(totnum), 'H+' : np.zeros(totnum),  
            'Na+' : np.zeros(totnum), 'K+' : np.zeros(totnum), 
            'OH-': np.zeros(totnum), 'HCO3-': np.zeros(totnum), 
            'CO3-2': np.zeros(totnum), 'Cl-': np.zeros(totnum), 
            'CO2(aq)': np.zeros(totnum), 'CO2(g)': np.zeros(totnum), 
            'PCO2': np.zeros(totnum), 'SiO2(aq)': np.zeros(totnum),
            'ALK': np.zeros(totnum), 'DIC': np.zeros(totnum), 
            'Calcite': np.zeros(totnum), 'Magnesite': np.zeros(totnum), 
            'Siderite': np.zeros(totnum), 'Dolomite': np.zeros(totnum), 
            'Wollastonite': np.zeros(totnum), 'Lime': np.zeros(totnum), 
            'Clino-Enstatite': np.zeros(totnum), 'Ferrosilite': np.zeros(totnum), 
            'Antigorite': np.zeros(totnum), 'Fayalite': np.zeros(totnum), 
            'SiO2(a)': np.zeros(totnum), 'Quartz': np.zeros(totnum), 
            'Chalcedony': np.zeros(totnum), 'Cristobalite': np.zeros(totnum), 
            'Coesite': np.zeros(totnum), 'Diopside': np.zeros(totnum), 
            'pH': np.zeros(totnum)}

    return chem

def chem_dict2(numQ,totnum):
    '''
    Returns 2D Chemical Dictionary Object to be used for Plotting
    '''
    chem = {'Ca+2': np.zeros((numQ,totnum)), 'Mg+2': np.zeros((numQ,totnum)), 
            'Fe+2': np.zeros((numQ,totnum)), 'H+' : np.zeros((numQ,totnum)),  
            'Na+' : np.zeros((numQ,totnum)), 'K+' : np.zeros((numQ,totnum)), 
            'OH-': np.zeros((numQ,totnum)), 'HCO3-': np.zeros((numQ,totnum)), 
            'CO3-2': np.zeros((numQ,totnum)), 'Cl-': np.zeros((numQ,totnum)), 
            'CO2(aq)': np.zeros((numQ,totnum)), 'CO2(g)': np.zeros((numQ,totnum)), 
            'PCO2': np.zeros((numQ,totnum)), 'SiO2(aq)': np.zeros((numQ,totnum)),
            'ALK': np.zeros((numQ,totnum)), 'DIC': np.zeros((numQ,totnum)), 
            'Calcite': np.zeros((numQ,totnum)), 'Magnesite': np.zeros((numQ,totnum)), 
            'Siderite': np.zeros((numQ,totnum)), 'Dolomite': np.zeros((numQ,totnum)), 
            'Wollastonite': np.zeros((numQ,totnum)), 'Lime': np.zeros((numQ,totnum)), 
            'Clino-Enstatite': np.zeros((numQ,totnum)), 'Ferrosilite': np.zeros((numQ,totnum)), 
            'Antigorite': np.zeros((numQ,totnum)), 'Fayalite': np.zeros((numQ,totnum)), 
            'SiO2(a)': np.zeros((numQ,totnum)), 'Quartz': np.zeros((numQ,totnum)), 
            'Chalcedony': np.zeros((numQ,totnum)), 'Cristobalite': np.zeros((numQ,totnum)), 
            'Coesite': np.zeros((numQ,totnum)), 'Diopside': np.zeros((numQ,totnum)), 
            'pH': np.zeros((numQ,totnum))}

    return chem

def chem_dict3(numQ1,numQ2,totnum):
    '''
    Returns 3D Chemical Dictionary Object to be used for Plotting
    '''
    chem = {'Ca+2': np.zeros((numQ1,numQ2,totnum)), 'Mg+2': np.zeros((numQ1,numQ2,totnum)), 
            'Fe+2': np.zeros((numQ1,numQ2,totnum)), 'H+' : np.zeros((numQ1,numQ2,totnum)),  
            'Na+' : np.zeros((numQ1,numQ2,totnum)), 'K+' : np.zeros((numQ1,numQ2,totnum)), 
            'OH-': np.zeros((numQ1,numQ2,totnum)), 'HCO3-': np.zeros((numQ1,numQ2,totnum)), 
            'CO3-2': np.zeros((numQ1,numQ2,totnum)), 'Cl-': np.zeros((numQ1,numQ2,totnum)), 
            'CO2(aq)': np.zeros((numQ1,numQ2,totnum)), 'CO2(g)': np.zeros((numQ1,numQ2,totnum)), 
            'PCO2': np.zeros((numQ1,numQ2,totnum)), 'SiO2(aq)': np.zeros((numQ1,numQ2,totnum)),
            'ALK': np.zeros((numQ1,numQ2,totnum)), 'DIC': np.zeros((numQ1,numQ2,totnum)), 
            'Calcite': np.zeros((numQ1,numQ2,totnum)), 'Magnesite': np.zeros((numQ1,numQ2,totnum)), 
            'Siderite': np.zeros((numQ1,numQ2,totnum)), 'Dolomite': np.zeros((numQ1,numQ2,totnum)), 
            'Wollastonite': np.zeros((numQ1,numQ2,totnum)), 'Lime': np.zeros((numQ1,numQ2,totnum)), 
            'Clino-Enstatite': np.zeros((numQ1,numQ2,totnum)), 'Ferrosilite': np.zeros((numQ1,numQ2,totnum)), 
            'Antigorite': np.zeros((numQ1,numQ2,totnum)), 'Fayalite': np.zeros((numQ1,numQ2,totnum)), 
            'SiO2(a)': np.zeros((numQ1,numQ2,totnum)), 'Quartz': np.zeros((numQ1,numQ2,totnum)), 
            'Chalcedony': np.zeros((numQ1,numQ2,totnum)), 'Cristobalite': np.zeros((numQ1,numQ2,totnum)), 
            'Coesite': np.zeros((numQ1,numQ2,totnum)), 'Diopside': np.zeros((numQ1,numQ2,totnum)), 
            'pH': np.zeros((numQ1,numQ2,totnum))}

    return chem


# Save chemical species in 1D dictionary objects in units of number density [dm^-3] 

def save_chems1_Ca(state, PCO2, chems1, j):
    '''
    Returns chems1 dictionary object for Ca by updating chems1[j]
    '''    
    chems1['Ca+2'][j] = numden*state.speciesAmount('Ca+2')[0]
    chems1['H+'][j] = numden*state.speciesAmount('H+')[0]
    chems1['OH-'][j] = numden*state.speciesAmount('OH-')[0]
    chems1['CO3-2'][j] = numden*state.speciesAmount('CO3-2')[0]
    chems1['HCO3-'][j] = numden*state.speciesAmount('HCO3-')[0]
    chems1['SiO2(aq)'][j] = numden*state.speciesAmount('SiO2(aq)')[0]

    chems1['CO2(aq)'][j] = numden*state.speciesAmount('CO2(aq)')[0]
    chems1['CO2(g)'][j] = numden*state.speciesAmount('CO2(g)')[0]
    chems1['PCO2'][j] = PCO2

    chems1['Calcite'][j] = numden*state.speciesAmount('Calcite')[0]
    chems1['Wollastonite'][j] = numden*state.speciesAmount('Wollastonite')[0]

    chems1['pH'][j] = -np.log10(state.speciesAmount('H+')[0])
    
    return chems1

def save_chems1_Mg(state, PCO2, chems1, j):
    '''
    Returns chems1 dictionary object for Mg by updating chems1[j]
    '''  
    chems1['Mg+2'][j] = numden*state.speciesAmount('Mg+2')[0]
    chems1['H+'][j] = numden*state.speciesAmount('H+')[0]
    chems1['OH-'][j] = numden*state.speciesAmount('OH-')[0]
    chems1['CO3-2'][j] = numden*state.speciesAmount('CO3-2')[0]
    chems1['HCO3-'][j] = numden*state.speciesAmount('HCO3-')[0]
    chems1['SiO2(aq)'][j] = numden*state.speciesAmount('SiO2(aq)')[0]

    chems1['CO2(aq)'][j] = numden*state.speciesAmount('CO2(aq)')[0]
    chems1['CO2(g)'][j] = numden*state.speciesAmount('CO2(g)')[0]
    chems1['PCO2'][j] = PCO2

    chems1['Magnesite'][j] = numden*state.speciesAmount('Magnesite')[0]
    chems1['Clino-Enstatite'][j] = numden*state.speciesAmount('Clino-Enstatite')[0]

    chems1['pH'][j] = -np.log10(state.speciesAmount('H+')[0])
            
    return chems1

def save_chems1_Fe(state, PCO2, chems1, j):
    '''
    Returns chems1 dictionary object for Fe by updating chems1[j]
    '''  
    chems1['Fe+2'][j] = numden*state.speciesAmount('Fe+2')[0]
    chems1['H+'][j] = numden*state.speciesAmount('H+')[0]
    chems1['OH-'][j] = numden*state.speciesAmount('OH-')[0]
    chems1['CO3-2'][j] = numden*state.speciesAmount('CO3-2')[0]
    chems1['HCO3-'][j] = numden*state.speciesAmount('HCO3-')[0]
    chems1['SiO2(aq)'][j] = numden*state.speciesAmount('SiO2(aq)')[0]

    chems1['CO2(aq)'][j] = numden*state.speciesAmount('CO2(aq)')[0]
    chems1['CO2(g)'][j] = numden*state.speciesAmount('CO2(g)')[0]
    chems1['PCO2'][j] = PCO2

    chems1['Siderite'][j] = numden*state.speciesAmount('Siderite')[0]
    chems1['Fayalite'][j] = numden*state.speciesAmount('Fayalite')[0]

    chems1['pH'][j] = -np.log10(state.speciesAmount('H+')[0])
            
    return chems1


# Save chemical species in 2D dictionary objects in units of number density [dm^-3] 

def save_chems2_an_Ca(PCO2, logK3, logK9, logK16, nDIV_fixed, nDIV_numerical, chems2_an, chems2_san, i, j):
    '''
    Returns chems2 dictionary objects for Ca by updating chems2_an[i][j] and chems2_san[i][j] 
    '''
    if i == 0: 
        chems2_an['pH'][i][j] = ocean_pH_low(PCO2, logK3)
    else:
        chems2_san['pH'][i][j] = ocean_pH_upp(PCO2, nDIV_numerical, logK9, logK16)
        chems2_an['pH'][i][j] = ocean_pH_upp(PCO2, nDIV_fixed, logK9, logK16)
        
    return chems2_an, chems2_san

def save_chems2_Ca(state, PCO2, chems2, i, j):
    '''
    Returns chems2 dictionary object for Ca by updating chems2[i][j] 
    '''
    chems2['Ca+2'][i][j] = numden*state.speciesAmount('Ca+2')[0]
    chems2['H+'][i][j] = numden*state.speciesAmount('H+')[0]
    chems2['OH-'][i][j] = numden*state.speciesAmount('OH-')[0]
    chems2['CO3-2'][i][j] = numden*state.speciesAmount('CO3-2')[0]
    chems2['HCO3-'][i][j] = numden*state.speciesAmount('HCO3-')[0]
    chems2['SiO2(aq)'][i][j] = numden*state.speciesAmount('SiO2(aq)')[0]

    chems2['CO2(aq)'][i][j] = numden*state.speciesAmount('CO2(aq)')[0]
    chems2['CO2(g)'][i][j] = numden*state.speciesAmount('CO2(g)')[0]
    chems2['PCO2'][i][j] = PCO2

    chems2['Calcite'][i][j] = numden*state.speciesAmount('Calcite')[0]
    chems2['Wollastonite'][i][j] = numden*state.speciesAmount('Wollastonite')[0]
    chems2['Quartz'][i][j] = numden*state.speciesAmount('Quartz')[0]

    chems2['pH'][i][j] = -np.log10(state.speciesAmount('H+')[0])
    
    return chems2

def save_chems2_Mg(state, PCO2, chems2, i, j):
    '''
    Returns chems2 dictionary object for Mg by updating chems2[i][j] 
    '''
    chems2['Mg+2'][i][j] = numden*state.speciesAmount('Mg+2')[0]
    chems2['H+'][i][j] = numden*state.speciesAmount('H+')[0]
    chems2['OH-'][i][j] = numden*state.speciesAmount('OH-')[0]
    chems2['CO3-2'][i][j] = numden*state.speciesAmount('CO3-2')[0]
    chems2['HCO3-'][i][j] = numden*state.speciesAmount('HCO3-')[0]
    chems2['SiO2(aq)'][i][j] = numden*state.speciesAmount('SiO2(aq)')[0]

    chems2['CO2(aq)'][i][j] = numden*state.speciesAmount('CO2(aq)')[0]
    chems2['CO2(g)'][i][j] = numden*state.speciesAmount('CO2(g)')[0]
    chems2['PCO2'][i][j] = PCO2

    chems2['Magnesite'][i][j] = numden*state.speciesAmount('Magnesite')[0]
    chems2['Clino-Enstatite'][i][j] = numden*state.speciesAmount('Clino-Enstatite')[0]
    chems2['Quartz'][i][j] = numden*state.speciesAmount('Quartz')[0]

    chems2['pH'][i][j] = -np.log10(state.speciesAmount('H+')[0])
    
    return chems2

def save_chems2_Fe(state, PCO2, chems2, i, j):
    '''
    Returns chems2 dictionary object for Fe by updating chems2[i][j] 
    '''
    chems2['Fe+2'][i][j] = numden*state.speciesAmount('Fe+2')[0]
    chems2['H+'][i][j] = numden*state.speciesAmount('H+')[0]
    chems2['OH-'][i][j] = numden*state.speciesAmount('OH-')[0]
    chems2['CO3-2'][i][j] = numden*state.speciesAmount('CO3-2')[0]
    chems2['HCO3-'][i][j] = numden*state.speciesAmount('HCO3-')[0]
    chems2['SiO2(aq)'][i][j] = numden*state.speciesAmount('SiO2(aq)')[0]

    chems2['CO2(aq)'][i][j] = numden*state.speciesAmount('CO2(aq)')[0]
    chems2['CO2(g)'][i][j] = numden*state.speciesAmount('CO2(g)')[0]
    chems2['PCO2'][i][j] = PCO2

    chems2['Siderite'][i][j] = numden*state.speciesAmount('Siderite')[0]
    chems2['Fayalite'][i][j] = numden*state.speciesAmount('Fayalite')[0]
    chems2['Quartz'][i][j] = numden*state.speciesAmount('Quartz')[0]

    chems2['pH'][i][j] = -np.log10(state.speciesAmount('H+')[0])
    
    return chems2


# Save chemical species in 3D dictionary objects in units of number density [dm^-3] 

def save_chems3_Ca(state, PCO2, chems3, i, j, k):
    '''
    Returns chems3 dictionary object for Ca by updating chems3[k][i][j]
    '''    
    chems3['Ca+2'][k][i][j] = numden*state.speciesAmount('Ca+2')[0]
    chems3['H+'][k][i][j] = numden*state.speciesAmount('H+')[0]
    chems3['OH-'][k][i][j] = numden*state.speciesAmount('OH-')[0]
    chems3['CO3-2'][k][i][j] = numden*state.speciesAmount('CO3-2')[0]
    chems3['HCO3-'][k][i][j] = numden*state.speciesAmount('HCO3-')[0]
    chems3['SiO2(aq)'][k][i][j] = numden*state.speciesAmount('SiO2(aq)')[0]

    chems3['CO2(aq)'][k][i][j] = numden*state.speciesAmount('CO2(aq)')[0]
    chems3['CO2(g)'][k][i][j] = numden*state.speciesAmount('CO2(g)')[0]
    chems3['PCO2'][k][i][j] = PCO2

    chems3['Calcite'][k][i][j] = numden*state.speciesAmount('Calcite')[0]
    chems3['Wollastonite'][k][i][j] = numden*state.speciesAmount('Wollastonite')[0]
    chems3['Quartz'][k][i][j] = numden*state.speciesAmount('Quartz')[0]

    chems3['pH'][k][i][j] = -np.log10(state.speciesAmount('H+')[0])
    
    return chems3

def save_chems3_Mg(state, PCO2, chems3, i, j, k):
    '''
    Returns chems3 dictionary object for Mg by updating chems3[k][i][j]
    '''    
    chems3['Mg+2'][k][i][j] = numden*state.speciesAmount('Mg+2')[0]
    chems3['H+'][k][i][j] = numden*state.speciesAmount('H+')[0]
    chems3['OH-'][k][i][j] = numden*state.speciesAmount('OH-')[0]
    chems3['CO3-2'][k][i][j] = numden*state.speciesAmount('CO3-2')[0]
    chems3['HCO3-'][k][i][j] = numden*state.speciesAmount('HCO3-')[0]
    chems3['SiO2(aq)'][k][i][j] = numden*state.speciesAmount('SiO2(aq)')[0]

    chems3['CO2(aq)'][k][i][j] = numden*state.speciesAmount('CO2(aq)')[0]
    chems3['CO2(g)'][k][i][j] = numden*state.speciesAmount('CO2(g)')[0]
    chems3['PCO2'][k][i][j] = PCO2

    chems3['Magnesite'][k][i][j] = numden*state.speciesAmount('Magnesite')[0]
    chems3['Clino-Enstatite'][k][i][j] = numden*state.speciesAmount('Clino-Enstatite')[0]
    chems3['Quartz'][k][i][j] = numden*state.speciesAmount('Quartz')[0]

    chems3['pH'][k][i][j] = -np.log10(state.speciesAmount('H+')[0])
    return chems3

def save_chems3_Fe(state, PCO2, chems3, i, j, k):
    '''
    Returns chems3 dictionary object for Fe by updating chems3[k][i][j]
    '''    
    chems3['Fe+2'][k][i][j] = numden*state.speciesAmount('Fe+2')[0]
    chems3['H+'][k][i][j] = numden*state.speciesAmount('H+')[0]
    chems3['OH-'][k][i][j] = numden*state.speciesAmount('OH-')[0]
    chems3['CO3-2'][k][i][j] = numden*state.speciesAmount('CO3-2')[0]
    chems3['HCO3-'][k][i][j] = numden*state.speciesAmount('HCO3-')[0]
    chems3['SiO2(aq)'][k][i][j] = numden*state.speciesAmount('SiO2(aq)')[0]

    chems3['CO2(aq)'][k][i][j] = numden*state.speciesAmount('CO2(aq)')[0]
    chems3['CO2(g)'][k][i][j] = numden*state.speciesAmount('CO2(g)')[0]
    chems3['PCO2'][k][i][j] = PCO2

    chems3['Siderite'][k][i][j] = numden*state.speciesAmount('Siderite')[0]
    chems3['Fayalite'][k][i][j] = numden*state.speciesAmount('Fayalite')[0]
    chems3['Quartz'][k][i][j] = numden*state.speciesAmount('Quartz')[0]

    chems3['pH'][k][i][j] = -np.log10(state.speciesAmount('H+')[0])
    return chems3

