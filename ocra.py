#!/usr/bin/env python
# coding: utf-8

# # OCRA: Ocean Chemistry with Reacktoro And beyond
# ### This Python code implements Reaktoro software to calculate ocean chemistry
# 
# ## Reference: Hakim et al. (2023) ApJL
# 
# ### okra.py # contains functions to calculate ocean pH, CCD and stable phases

# In[ ]:


# Import libraries

import numpy as np
import pandas as pd

from store import *
from solve import *
from output import *


# In[ ]:


# Calculate pH as a function of PCO2 for Ca, Mg or Fe systems

def pH_PCO2(DIV = 'Ca', totP = 1, Temp = 288, nDIV = 1, nSiO2 = 1, totnum = 100,
            plot_flag = True, table_flag = True):
    '''
    Returns ocean pH as a function of PCO2 [bar] for Ca, Mg or Fe carbonate systems
    '''
        
    if DIV != 'Ca' and DIV != 'Mg' and DIV != 'Fe':
        print('Error: Enter DIV = "Ca" or "Mg" or "Fe"')

    PCO2s = np.logspace(-8, -0.5, num=totnum) # bar

    ####################################
    
    if DIV == 'Ca':

        numQ = 3
        betas = np.array([-1, 0, 0.3]) 

        chems2 = chem_dict2(numQ, totnum)

        system, specs, solver = setup_Ca()

        for i in range(numQ):
            j = 0

            beta = betas[i]
            while j < totnum:

                PCO2 = PCO2s[j]

                if beta == -1:
                    addDIVtot = 0
                    addSiO2   = 0
                elif beta == 0:
                    addDIVtot = 1e-2
                    addSiO2   = 0
                else:
                    addDIVtot = nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden
                    addSiO2   = nSiO2 * addDIVtot

                state = solve_Ca(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
                
                chems2 = save_chems2_Ca(state, PCO2, chems2, i, j)

                j = j + 1

        output_pH_PCO2(PCO2s, chems2, DIV = DIV, plot_flag = plot_flag, table_flag = table_flag)
    
    ####################################
    
    elif DIV == 'Mg':

        numQ = 3
        betas = np.array([-1, 0, 0.3]) 

        chems2 = chem_dict2(numQ, totnum)

        system, specs, solver = setup_Mg()

        for i in range(numQ):
            j = 0

            beta = betas[i]
            while j < totnum:

                PCO2 = PCO2s[j]

                if beta == -1:
                    addDIVtot = 0
                    addSiO2   = 0
                elif beta == 0:
                    addDIVtot = 1e-2
                    addSiO2   = 0
                else:
                    addDIVtot = nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden
                    addSiO2   = nSiO2 * addDIVtot

                state = solve_Mg(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
                
                chems2 = save_chems2_Mg(state, PCO2, chems2, i, j)

                j = j + 1

        output_pH_PCO2(PCO2s, chems2, DIV = DIV, plot_flag = plot_flag, table_flag = table_flag)
        
    ####################################
    
    elif DIV == 'Fe':

        numQ = 3
        betas = np.array([-1, 0, 0.3]) 

        chems2 = chem_dict2(numQ, totnum)

        system, specs, solver = setup_Fe()

        for i in range(numQ):
            j = 0

            beta = betas[i]
            while j < totnum:

                PCO2 = PCO2s[j]

                if beta == -1:
                    addDIVtot = 0
                    addSiO2   = 0
                elif beta == 0:
                    addDIVtot = 1e-2
                    addSiO2   = 0
                else:
                    addDIVtot = nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden
                    addSiO2   = nSiO2 * addDIVtot

                state = solve_Fe(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
                
                chems2 = save_chems2_Fe(state, PCO2, chems2, i, j)

                j = j + 1

        output_pH_PCO2(PCO2s, chems2, DIV = DIV, plot_flag = plot_flag, table_flag = table_flag)

    return


# In[ ]:


# Calculate analytical and numerical solutions of ocean pH

def pH_PCO2_an(DIV = 'Ca', totP = 1, Temp = 288, nDIV = 1, nSiO2 = 1, totnum = 100, nDIV_fixed = 1,
           plot_flag = True, table_flag = True):
    '''
    Returns analytical/numerical solutions of ocean pH as a function of PCO2 [bar]
    '''    
    if DIV != 'Ca':
        print('Error: Enter DIV = "Ca"')

    PCO2s = np.logspace(-8, -0.5, num=totnum) # bar
    
    ####################################
    
    if DIV == 'Ca':
        
        numQ = 2
        addDIVtots = numden*np.linspace(0, 1e-2, num=numQ) # in units of 1000 * mol/dm3 = 1 mol/m3
        chems2 = chem_dict2(numQ,totnum)
        chems2_san = chem_dict2(numQ,totnum)
        chems2_an = chem_dict2(numQ,totnum)

        system, specs, solver = setup_Ca()
        
        logK3, logK9, logK16 = setup_an_Ca(Temp, totP)
        
        for i in range(numQ):
            j = 0
            addDIVtot = addDIVtots[i] / numden
            addSiO2   = 0

            while j < totnum:

                PCO2 = PCO2s[j]

                state = solve_Ca(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
                
                chems2 = save_chems2_Ca(state, PCO2, chems2, i, j)
                
                chems2_an, chems2_san = save_chems2_an_Ca(PCO2, logK3, logK9, logK16, nDIV_fixed,
                                                          chems2['Ca+2'][i][j], chems2_an, chems2_san, i, j)

                j = j + 1

        output_pH_PCO2_an(PCO2s, chems2, chems2_an, chems2_san, DIV = DIV, nDIV_fixed = nDIV_fixed,
                          plot_flag = plot_flag, table_flag = table_flag)
    
    return


# In[ ]:


# Calculate ocean pH as a function of local pressure

def pH_P(DIV = 'Ca', PCO2 = 0.3e-3, Temp = 288, nDIV = 1, nSiO2 = 1, totnum = 100,
           plot_flag = True, table_flag = True):
    '''
    Returns ocean pH as a function of P [bar]
    '''  
    if DIV != 'Ca':
        print('Error: Enter DIV = "Ca"')

    totPs = np.logspace(0, 3, num=totnum) # bar

    ####################################
    
    if DIV == 'Ca':
        
        numQ = 3
        betas = np.array([-1, 0, 0.3]) 

        chems2 = chem_dict2(numQ, totnum)

        system, specs, solver = setup_Ca()

        for i in range(numQ):
            j = 0

            beta = betas[i]
            while j < totnum:

                totP = totPs[j]

                if beta == -1:
                    addDIVtot = 0
                    addSiO2   = 0
                elif beta == 0:
                    addDIVtot = 1e-2
                    addSiO2 = 0
                else:
                    addDIVtot = nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden
                    addSiO2   = nSiO2 * addDIVtot

                state = solve_Ca(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
                
                chems2 = save_chems2_Ca(state, PCO2, chems2, i, j)

                j = j + 1
                
        output_pH_P(totPs, chems2, DIV = DIV, plot_flag = plot_flag, table_flag = table_flag)

    return


# In[ ]:


# Calculate ocean pH as a function of temperature

def pH_T(DIV = 'Ca', PCO2 = 0.3e-3, totP = 1, nDIV = 1, nSiO2 = 1, totnum = 100,
           plot_flag = True, table_flag = True):
    '''
    Returns ocean pH as a function of Temp [K]
    '''  
    if DIV != 'Ca':
        print('Error: Enter DIV = "Ca"')

    Temps = np.linspace(273.16, 372.16, num=totnum) # bar

    ####################################
    
    if DIV == 'Ca':

        numQ = 3
        betas = np.array([-1, 0, 0.3]) 

        chems2 = chem_dict2(numQ, totnum)

        system, specs, solver = setup_Ca()

        for i in range(numQ):
            j = 0

            beta = betas[i]
            while j < totnum:

                Temp = Temps[j]

                if beta == -1:
                    addDIVtot = 0
                    addSiO2   = 0
                elif beta == 0:
                    addDIVtot = 1e-2
                    addSiO2   = 0
                else:
                    addDIVtot = nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden
                    addSiO2   = nSiO2 * addDIVtot

                state = solve_Ca(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
                
                chems2 = save_chems2_Ca(state, PCO2, chems2, i, j)

                j = j + 1

        output_pH_T(Temps, chems2, DIV = DIV, plot_flag = plot_flag, table_flag = table_flag)
    
    return


# In[ ]:


# Calculate Ca-CCD as a function of PCO2 and T

def CaCCD_PCO2_T(beta = 0.3, nSiO2 = 1, nDIV = 1, totnum = 10, numQ1 = 10, numQ2 = 10,
                plot_flag = True, table_flag = True):
    '''
    Returns Ca-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    if totnum > 10:
        print('Please be patient. A high-resolution figure is being generated.')
        
    Temps = np.linspace(273.16, 372.16, num=numQ1) # Temperature in K
    PCO2s = np.logspace(-8, -0.5, num=numQ2) # surface CO2 pressure in bar
    totPs = np.logspace(0, np.log10(5000), num=totnum)

    chems3 = chem_dict3(numQ1, numQ2, totnum)

    CCDs = np.zeros((numQ1, numQ2))
    
    system, specs, solver = setup_Ca()
                
    for k in range(len(Temps)):
        Temp = Temps[k] 

        for i in range(len(PCO2s)):
            PCO2 = PCO2s[i]

            addDIVtot = nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden

            addSiO2 = nSiO2 * addDIVtot

            j = 0
            while j < totnum:

                totP = totPs[j]

                state = solve_Ca(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
                
                chems3 = save_chems3_Ca(state, PCO2, chems3, i, j, k)

                j = j + 1

    ####################################

    for k in range(len(Temps)):
        for i in range(len(PCO2s)):
            nCarb_surf = chems3['Calcite'][k][i][0]
            if nCarb_surf < low_cutoff:
                CCDs[k][i] = 1e-3 # 0 # km
            else:
                CCDs[k][i] = 100 # km
            j = 0
            while j < totnum:
                if chems3['Calcite'][k][i][j] < 0.001 * nCarb_surf:
                    CCDs[k][i] = ocean_depth(totPs[j])
                    j = totnum
                else:
                    j = j + 1
    
    output_CaCCD_PCO2_T(PCO2s, Temps, CCDs, beta = beta, nSiO2 = nSiO2, nDIV = nDIV,
                        plot_flag = plot_flag, table_flag = table_flag)
        
    return


# In[ ]:


# Calculate Mg-CCD as a function of PCO2 and T

def MgCCD_PCO2_T(beta = 0.3, nSiO2 = 1, nDIV = 1, totnum = 10, numQ1 = 10, numQ2 = 10,
                plot_flag = True, table_flag = True):
    '''
    Returns Mg-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    if totnum > 10:
        print('Please be patient. A high-resolution figure is being generated.')
        
    Temps = np.linspace(273.16, 372.16, num=numQ1) # Temperature in K
    PCO2s = np.logspace(-8, -0.5, num=numQ2) # surface CO2 pressure in bar
    totPs = np.logspace(0, np.log10(5000), num=totnum)

    chems3 = chem_dict3(numQ1, numQ2, totnum)

    CCDs = np.zeros((numQ1, numQ2))
    
    system, specs, solver = setup_Mg()
                
    for k in range(len(Temps)):
        Temp = Temps[k] 

        for i in range(len(PCO2s)):
            PCO2 = PCO2s[i]

            addDIVtot = nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden

            addSiO2 = nSiO2 * addDIVtot

            j = 0
            while j < totnum:

                totP = totPs[j]

                state = solve_Mg(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
                
                chems3 = save_chems3_Mg(state, PCO2, chems3, i, j, k)

                j = j + 1

    ####################################

    for k in range(len(Temps)):
        for i in range(len(PCO2s)):
            nCarb_surf = chems3['Magnesite'][k][i][0]
            if nCarb_surf < low_cutoff:
                CCDs[k][i] = 1e-3 # 0 # km
            else:
                CCDs[k][i] = 100 # km
            j = 0
            while j < totnum:
                if chems3['Magnesite'][k][i][j] < 0.001 * nCarb_surf:
                    CCDs[k][i] = ocean_depth(totPs[j])
                    j = totnum
                else:
                    j = j + 1
    
    output_MgCCD_PCO2_T(PCO2s, Temps, CCDs, beta = beta, nSiO2 = nSiO2, nDIV = nDIV,
                        plot_flag = plot_flag, table_flag = table_flag)
        
    return


# In[ ]:


# Calculate Fe-CCD as a function of PCO2 and T

def FeCCD_PCO2_T(beta = 0.3, nSiO2 = 1, nDIV = 1, totnum = 10, numQ1 = 10, numQ2 = 10,
                plot_flag = True, table_flag = True):
    '''
    Returns Fe-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    if totnum > 10:
        print('Please be patient. A high-resolution figure is being generated.')
        
    Temps = np.linspace(273.16, 372.16, num=numQ1) # Temperature in K
    PCO2s = np.logspace(-8, -0.5, num=numQ2) # surface CO2 pressure in bar
    totPs = np.logspace(0, np.log10(5000), num=totnum)

    chems3 = chem_dict3(numQ1, numQ2, totnum)

    CCDs = np.zeros((numQ1, numQ2))
    
    system, specs, solver = setup_Fe()
                
    for k in range(len(Temps)):
        Temp = Temps[k] 

        for i in range(len(PCO2s)):
            PCO2 = PCO2s[i]

            addDIVtot = nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden

            addSiO2 = nSiO2 * addDIVtot

            j = 0
            while j < totnum:

                totP = totPs[j]

                state = solve_Fe(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
                
                chems3 = save_chems3_Fe(state, PCO2, chems3, i, j, k)

                j = j + 1

    ####################################

    for k in range(len(Temps)):
        for i in range(len(PCO2s)):
            nCarb_surf = chems3['Siderite'][k][i][0]
            if nCarb_surf < low_cutoff:
                CCDs[k][i] = 1e-3 # 0 # km
            else:
                CCDs[k][i] = 100 # km
            j = 0
            while j < totnum:
                if chems3['Siderite'][k][i][j] < 0.001 * nCarb_surf:
                    CCDs[k][i] = ocean_depth(totPs[j])
                    j = totnum
                else:
                    j = j + 1
    
    output_FeCCD_PCO2_T(PCO2s, Temps, CCDs, beta = beta, nSiO2 = nSiO2, nDIV = nDIV,
                        plot_flag = plot_flag, table_flag = table_flag)
        
    return


# In[ ]:


# Calculate stable phases as a function of PCO2

def phases_PCO2 (DIV = 'Ca', Temp = 298, totP = 1, beta = 0.3, nDIV = 1, nSiO2 = 1, totnum = 100,
                table_flag = True, plot_flag = True):
    '''
    Returns stable phases as a function of PCO2 [bar]
    '''  
    if DIV != 'Ca' and DIV != 'Mg' and DIV != 'Fe':
        print('Error: Enter DIV = "Ca" or "Mg" or "Fe"')
    
    PCO2s = np.logspace(-8, -0.5, num=totnum) # bar
    chems1 = chem_dict1(totnum)

    ####################################
    
    if DIV == 'Ca':
        
        system, specs, solver = setup_Ca()

        j = 0
        while j < totnum:

            PCO2 = PCO2s[j]
            
            addDIVtot = nDIV * weath_scaling(PCO2,Temp,beta=beta) / numden
            addSiO2 = nSiO2 * addDIVtot

            state = solve_Ca(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
                
            chems1 = save_chems1_Ca(state, PCO2, chems1, j)

            j = j + 1
    
        df = pd.DataFrame({
            'Ca++': chems1['Ca+2'],
            'Calcite': chems1['Calcite'],
            'Silicates': chems1['Wollastonite'],
        }, index=PCO2s)
        
        output_phases_PCO2(df, DIV = DIV, beta = beta, nDIV = nDIV, nSiO2 = nSiO2,
                           table_flag = table_flag, plot_flag = plot_flag)
        
    ####################################
    
    elif DIV == 'Mg':
        
        system, specs, solver = setup_Mg()

        j = 0
        while j < totnum:

            PCO2 = PCO2s[j]
            
            addDIVtot = nDIV * weath_scaling(PCO2,Temp,beta=beta) / numden
            addSiO2 = nSiO2 * addDIVtot

            state = solve_Mg(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
                
            chems1 = save_chems1_Mg(state, PCO2, chems1, j)

            j = j + 1

        df = pd.DataFrame({
            'Mg++': chems1['Mg+2'],
            'Magnesite': chems1['Magnesite'],
            'Silicates': 2*chems1['Clino-Enstatite'],
        }, index=PCO2s)
        
        output_phases_PCO2(df, DIV = DIV, beta = beta, nDIV = nDIV, nSiO2 = nSiO2,
                           table_flag = table_flag, plot_flag = plot_flag)
    
    ####################################
    
    elif DIV == 'Fe':
        
        system, specs, solver = setup_Fe()

        j = 0
        while j < totnum:

            PCO2 = PCO2s[j]
            
            addDIVtot = nDIV * weath_scaling(PCO2,Temp,beta=beta) / numden
            addSiO2 = nSiO2 * addDIVtot

            state = solve_Fe(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP)
                
            chems1 = save_chems1_Fe(state, PCO2, chems1, j)

            j = j + 1

        df = pd.DataFrame({
            'Fe++': chems1['Fe+2'],
            'Siderite': chems1['Siderite'],
            'Silicates': 2*chems1['Fayalite'],
        }, index=PCO2s)
        
        output_phases_PCO2(df, DIV = DIV, beta = beta, nDIV = nDIV, nSiO2 = nSiO2,
                           table_flag = table_flag, plot_flag = plot_flag)
        
    return

