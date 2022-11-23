#!/usr/bin/env python
# coding: utf-8

# # OCRA: Ocean Chemistry with Reacktoro And beyond
# ## Hakim et al. (2023) ApJL
# 
# This Python script implements Reaktoro software to calculate ocean chemistry

# In[ ]:


# Import libraries

from reaktoro import *

import matplotlib.pyplot as plt
plt.rcParams['pcolor.shading'] ='nearest'
import matplotlib.gridspec as gridspec
from matplotlib import ticker

import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

from astropy.constants import R, sigma_sb

from matplotlib.cm import get_cmap

col1 = get_cmap('Dark2').colors  # type: matplotlib.colors.ListedColormap
col2 = get_cmap('Set1').colors
col3 = get_cmap('Set3').colors
colors = col1 + col2 + col3


# In[ ]:


# Set global constants

numden = 1000 # 1 mol/kg = 1000 mol/m3

totH2O = 55.5 # moles
totN2 = 0.1 # moles --- produces pressure equivalent to 0.8 bar

low_cutoff = 1e-10 # mol/m3


# In[ ]:


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


# In[ ]:


# Custom dictionary objects to store chemical species as number density in m^-3

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


# In[ ]:


# Calculate Ca-CCD as a function of PCO2 and T
def CaCCD_PCO2_T(beta = 0.3, nSiO2 = 1, nDIV = 1, totnum = 10, numQ1 = 10, numQ2 = 10,
                plot_flag = True, table_flag = True):
    '''
    Returns tabulated values of Ca-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    if totnum > 10:
        print('Please be patient. A high-resolution figure is being generated.')
        
    Temps = np.linspace(273.16, 372.16, num=numQ1) # Temperature in K
    PCO2s = np.logspace(-8, -0.5, num=numQ2) # surface CO2 pressure in bar
    totPs = np.logspace(0, np.log10(5000), num=totnum)

    chems = chem_dict3(numQ1, numQ2, totnum)

    CCDs = np.zeros((numQ1, numQ2))
    
    ####################################
    
    db = SupcrtDatabase('supcrtbl')

    solution = AqueousPhase(['H2O(aq)','CO2(aq)', 'HCO3-', 'CO3-2', 'H+', 'OH-', 'Ca+2', 'SiO2(aq)'])
    solution.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond('CO2(aq)'),
    ))

    gases = GaseousPhase(['CO2(g)', 'N2(g)'])
    gases.setActivityModel(ActivityModelPengRobinson())

    minerals = MineralPhases(['Calcite', 'Wollastonite', 'Quartz'])

    system = ChemicalSystem(db, solution, gases, minerals)

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.fugacity('CO2')

    solver = EquilibriumSolver(specs)
        
    ####################################
                
    for k in range(len(Temps)):
        Temp = Temps[k] 

        for i in range(len(PCO2s)):
            PCO2 = PCO2s[i]

            addDIVtot = nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden

            addSiO2 = nSiO2 * addDIVtot

            j = 0
            while j < totnum:

                totP = totPs[j]

                state = ChemicalState(system)
                state.setTemperature(Temp, 'K')
                state.setPressure(totP, 'bar')
                state.set('H2O(aq)', totH2O - addDIVtot, 'mol')     # add ~ one kg of water
                state.set('N2(g)', totN2, 'mol')
                state.set('HCO3-', 2*addDIVtot, 'mol')
                state.set('Ca+2', addDIVtot, 'mol')
                state.set('SiO2(aq)', addSiO2, 'mol')

                conditions = EquilibriumConditions(specs)
                conditions.temperature(state.temperature())
                conditions.pressure(state.pressure())
                conditions.fugacity('CO2', PCO2, 'bar')

                result = solver.solve(state, conditions)

                assert result.optima.succeeded

                chems['Ca+2'][k][i][j] = numden*state.speciesAmount('Ca+2')[0]
                chems['H+'][k][i][j] = numden*state.speciesAmount('H+')[0]
                chems['OH-'][k][i][j] = numden*state.speciesAmount('OH-')[0]
                chems['CO3-2'][k][i][j] = numden*state.speciesAmount('CO3-2')[0]
                chems['HCO3-'][k][i][j] = numden*state.speciesAmount('HCO3-')[0]
                chems['SiO2(aq)'][k][i][j] = numden*state.speciesAmount('SiO2(aq)')[0]

                chems['CO2(aq)'][k][i][j] = numden*state.speciesAmount('CO2(aq)')[0]
                chems['CO2(g)'][k][i][j] = numden*state.speciesAmount('CO2(g)')[0]
                chems['PCO2'][k][i][j] = PCO2

                chems['Calcite'][k][i][j] = numden*state.speciesAmount('Calcite')[0]
                chems['Wollastonite'][k][i][j] = numden*state.speciesAmount('Wollastonite')[0]
                chems['Quartz'][k][i][j] = numden*state.speciesAmount('Quartz')[0]

                chems['pH'][k][i][j] = -np.log10(state.speciesAmount('H+')[0])

                j = j + 1

    ####################################

    for k in range(len(Temps)):
        for i in range(len(PCO2s)):
            nCarb_surf = chems['Calcite'][k][i][0]
            if nCarb_surf < low_cutoff:
                CCDs[k][i] = 1e-3 # 0 # km
            else:
                CCDs[k][i] = 100 # km
            j = 0
            while j < totnum:
                if chems['Calcite'][k][i][j] < 0.001 * nCarb_surf:
                    CCDs[k][i] = ocean_depth(totPs[j])
                    j = totnum
                else:
                    j = j + 1
    
    if table_flag == True:
        df = pd.DataFrame(CCDs, index=PCO2s, columns=Temps)
        if nSiO2 == 0:
            df.to_csv('figA3a.csv')
        elif nSiO2 == 1:
            df.to_csv('fig3a.csv')
        #df.to_csv('CaCCD_PCO2_T_'+str(int(beta*100))+'_nSiO2'+str(int(nSiO2))+'.csv')
        
    if plot_flag == True:
        plot_CaCCD_PCO2_T(PCO2s, Temps, CCDs, beta = beta, nSiO2 = nSiO2, nDIV = nDIV, 
                          totnum = totnum, numQ1 = numQ1, numQ2 = numQ2)
        
    return


# In[ ]:


# Plot Ca-CCD as a function of PCO2 and T
def plot_CaCCD_PCO2_T(PCO2s, Temps, CCDs, 
                       beta = 0.3, nSiO2 = 1, nDIV = 1, totnum = 10, numQ1 = 10, numQ2 = 10):
    '''
    Returns plots of Ca-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    
    if beta == 0:
        
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.15)

        ax.set_xscale('log')
        ax.set_yscale('linear')

        ax.set_xlabel(r'$P_{\rm CO_2}$ [bar]', fontsize=16)
        ax.set_ylabel(r'$T$ [K]', fontsize=16)
        
        ax.text(-0.1, 1.01, '(a)', transform=ax.transAxes, fontsize=16, va='top', ha='right')

        # contourf

        levels = np.array([1, 2, 4, 10, 20, 40]) 

        cf = ax.contourf(PCO2s, Temps, CCDs, levels = levels, extend = 'both', locator=ticker.LogLocator(), 
                          cmap = plt.get_cmap('viridis'))

        fig.colorbar(cf, ax=ax, label='CCD [km]')

        ax.scatter(0.3e-3, 288, color='gray', marker='o')

        ax.text(1e-5, 338, 'Carbon Cycle', c='black',
         ha='left', va='center'
        )
        
        ax.text(1e-8, 360, 'No Carbon Cycle \n(too little CO$_2$)', c='white',
         ha='left', va='center'
        )

        ax.text(1e-3, 280, 'No Carbon Cycle \n(too acidic)', c='white',
         ha='left', va='center'
        )

        ax.text(1e-5, 365, r'$n_{\rm Ca,tot}$ = %.1f m$^{-3}$'%nDIV, c='gray',
         ha='left', va='center'
        )

        ax.set_title(r'Ca-CCD', fontsize=18)
        # plt.savefig('CCD_Ca_beta0_nSiO2%s.pdf'%(int(nSiO2)), bbox_inches='tight')
        
    elif beta == 0.3:
        
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.15)

        ax.set_xscale('log')
        ax.set_yscale('linear')

        ax.set_xlabel(r'$P_{\rm CO_2}$ [bar]', fontsize=16)
        ax.set_ylabel(r'$T$ [K]', fontsize=16)
        
        ax.text(-0.1, 1.01, '(a)', transform=ax.transAxes, fontsize=16, va='top', ha='right')

        # contourf

        levels = np.array([1, 2, 4, 10, 20, 40]) 

        cf = ax.contourf(PCO2s, Temps, CCDs, levels = levels, extend = 'both', locator=ticker.LogLocator(), 
                          cmap = plt.get_cmap('viridis')) 

        PCO2sA, TempsA = np.meshgrid(PCO2s, Temps)

        nCarb = nDIV * weath_scaling(PCO2sA, TempsA, beta=beta)

        levels = np.logspace(-2, 2, num=3)  

        ax.contour(PCO2sA, TempsA, nCarb, levels = levels, locator=ticker.LogLocator(), 
                          colors='gray')

        fig.colorbar(cf, ax=ax, label='CCD [km]')

        ax.scatter(0.3e-3, 288, color='gray', marker='o')

        ax.text(1e-3, 340, r'$n_{Ca, tot} = 100$ m$^{-3}$', c='gray',
         ha='left', va='center', rotation=-28
        )

        ax.text(1e-5, 298, r'$n_{Ca, tot} = 1$ m$^{-3}$', c='gray',
         ha='left', va='center', rotation=-28
        )

        ax.text(1e-3, 320, 'Carbon Cycle', c='black',
         ha='left', va='center'
        )
        
        if nSiO2 != 0:
            ax.text(2e-8, 360, 'No Carbon Cycle \n(cations consumed \nby silicates)', c='white',
             ha='left', va='center'
            )
        
        ax.text(2e-8, 280, 'No Carbon Cycle \n(too little CO$_2$)', c='white',
         ha='left', va='center'
        )

        ax.text(1e-3, 280, 'No Carbon Cycle \n(too acidic)', c='white',
         ha='left', va='center'
        )

        ax.text(5e-4, 310, r'$n_{\rm Ca, tot} = f_{\rm W}(P_{\rm CO_2}, T)$', c='gray',
         ha='left', va='center'
        )

        ax.set_title(r'Ca-CCD', fontsize=18)
        if nSiO2 == 0:
            plt.savefig('figA3a.pdf', bbox_inches='tight')
        elif nSiO2 == 1:
            plt.savefig('fig3a.pdf', bbox_inches='tight')
        # plt.savefig('CCD_Ca_beta30_nSiO2%s.pdf'%(int(nSiO2)), bbox_inches='tight')
        
    return


# In[ ]:


# Calculate Mg-CCD as a function of PCO2 and T
def MgCCD_PCO2_T(beta = 0.3, nSiO2 = 1, nDIV = 1, totnum = 10, numQ1 = 10, numQ2 = 10,
                plot_flag = True, table_flag = True):
    '''
    Returns tabulated values of Mg-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    if totnum > 10:
        print('Please be patient. A high-resolution figure is being generated.')
        
    Temps = np.linspace(273.16, 372.16, num=numQ1) # Temperature in K
    PCO2s = np.logspace(-8, -0.5, num=numQ2) # surface CO2 pressure in bar
    totPs = np.logspace(0, np.log10(5000), num=totnum)

    chems = chem_dict3(numQ1, numQ2, totnum)

    CCDs = np.zeros((numQ1, numQ2))
    
    ####################################
    
    db = SupcrtDatabase('supcrtbl')

    solution = AqueousPhase(['H2O(aq)','CO2(aq)', 'HCO3-', 'CO3-2', 'H+', 'OH-', 'Mg+2', 'SiO2(aq)'])
    solution.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond('CO2(aq)'),
    ))

    gases = GaseousPhase(['CO2(g)', 'N2(g)'])
    gases.setActivityModel(ActivityModelPengRobinson())

    minerals = MineralPhases(['Magnesite', 'Clino-Enstatite', 'Quartz'])

    system = ChemicalSystem(db, solution, gases, minerals)

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.fugacity('CO2')
    
    solver = EquilibriumSolver(specs)
        
    ####################################
                
    for k in range(len(Temps)):
        Temp = Temps[k] 

        for i in range(len(PCO2s)):
            PCO2 = PCO2s[i]

            addDIVtot = nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden

            addSiO2 = nSiO2 * addDIVtot

            j = 0
            while j < totnum:

                totP = totPs[j]

                state = ChemicalState(system)
                state.setTemperature(Temp, 'K')
                state.setPressure(totP, 'bar')
                state.set('H2O(aq)', totH2O - addDIVtot, 'mol')     # add ~ one kg of water
                state.set('N2(g)', totN2, 'mol')
                state.set('HCO3-', 2*addDIVtot, 'mol')
                state.set('Mg+2', addDIVtot, 'mol')
                state.set('SiO2(aq)', addSiO2, 'mol')

                conditions = EquilibriumConditions(specs)
                conditions.temperature(state.temperature())
                conditions.pressure(state.pressure())
                conditions.fugacity('CO2', PCO2, 'bar')

                result = solver.solve(state, conditions)

                assert result.optima.succeeded

                chems['Mg+2'][k][i][j] = numden*state.speciesAmount('Mg+2')[0]
                chems['H+'][k][i][j] = numden*state.speciesAmount('H+')[0]
                chems['OH-'][k][i][j] = numden*state.speciesAmount('OH-')[0]
                chems['CO3-2'][k][i][j] = numden*state.speciesAmount('CO3-2')[0]
                chems['HCO3-'][k][i][j] = numden*state.speciesAmount('HCO3-')[0]
                chems['SiO2(aq)'][k][i][j] = numden*state.speciesAmount('SiO2(aq)')[0]

                chems['CO2(aq)'][k][i][j] = numden*state.speciesAmount('CO2(aq)')[0]
                chems['CO2(g)'][k][i][j] = numden*state.speciesAmount('CO2(g)')[0]
                chems['PCO2'][k][i][j] = PCO2

                chems['Magnesite'][k][i][j] = numden*state.speciesAmount('Magnesite')[0]
                chems['Clino-Enstatite'][k][i][j] = numden*state.speciesAmount('Clino-Enstatite')[0]
                chems['Quartz'][k][i][j] = numden*state.speciesAmount('Quartz')[0]

                chems['pH'][k][i][j] = -np.log10(state.speciesAmount('H+')[0])

                j = j + 1

    ####################################

    for k in range(len(Temps)):
        for i in range(len(PCO2s)):
            nCarb_surf = chems['Magnesite'][k][i][0]
            if nCarb_surf < low_cutoff:
                CCDs[k][i] = 1e-3 # 0 # km
            else:
                CCDs[k][i] = 100 # km
            j = 0
            while j < totnum:
                if chems['Magnesite'][k][i][j] < 0.001 * nCarb_surf:
                    CCDs[k][i] = ocean_depth(totPs[j])
                    j = totnum
                else:
                    j = j + 1
    
    if table_flag == True:
        df = pd.DataFrame(CCDs, index=PCO2s, columns=Temps)
        if nSiO2 == 0:
            df.to_csv('figA3b.csv')
        elif nSiO2 == 1:
            df.to_csv('fig3b.csv')
        #df.to_csv('MgCCD_PCO2_T_'+str(int(beta*100))+'_nSiO2'+str(int(nSiO2))+'.csv')
    
    if plot_flag == True:
        plot_MgCCD_PCO2_T(PCO2s, Temps, CCDs, beta = beta, nSiO2 = nSiO2, nDIV = nDIV, 
                          totnum = totnum, numQ1 = numQ1, numQ2 = numQ2)
        
    return


# In[ ]:


# Plot Mg-CCD as a function of PCO2 and T
def plot_MgCCD_PCO2_T(PCO2s, Temps, CCDs, 
                       beta = 0.3, nSiO2 = 1, nDIV = 1, totnum = 10, numQ1 = 10, numQ2 = 10):
    '''
    Returns plots of Mg-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    
    if beta == 0:
        
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.15)

        ax.set_xscale('log')
        ax.set_yscale('linear')

        ax.set_xlabel(r'$P_{\rm CO_2}$ [bar]', fontsize=16)
        ax.set_ylabel(r'$T$ [K]', fontsize=16)
        
        ax.text(-0.1, 1.01, '(b)', transform=ax.transAxes, fontsize=16, va='top', ha='right')

        # contourf

        levels = np.array([1, 2, 4, 10, 20, 40]) 

        cf = ax.contourf(PCO2s, Temps, CCDs, levels = levels, extend = 'both', locator=ticker.LogLocator(), 
                          cmap = plt.get_cmap('viridis'))

        fig.colorbar(cf, ax=ax, label='CCD [km]')

        ax.scatter(0.3e-3, 288, color='gray', marker='o')

        ax.text(1e-5, 338, 'Carbon Cycle', c='black',
         ha='left', va='center'
        )

        ax.text(1e-8, 360, 'No Carbon Cycle \n(too little CO$_2$)', c='white',
         ha='left', va='center'
        )

        ax.text(1e-3, 280, 'No Carbon Cycle \n(too acidic)', c='white',
         ha='left', va='center'
        )

        ax.text(1e-5, 365, r'$n_{\rm Mg,tot}$ = %.1f m$^{-3}$'%nDIV, c='gray',
         ha='left', va='center'
        )

        ax.set_title(r'Mg-CCD', fontsize=18)
        # plt.savefig('CCD_Mg_beta0_nSiO2%s.pdf'%(int(nSiO2)), bbox_inches='tight')
        
    elif beta == 0.3:
        
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.15)

        ax.set_xscale('log')
        ax.set_yscale('linear')

        ax.set_xlabel(r'$P_{\rm CO_2}$ [bar]', fontsize=16)
        ax.set_ylabel(r'$T$ [K]', fontsize=16)
        
        ax.text(-0.1, 1.01, '(b)', transform=ax.transAxes, fontsize=16, va='top', ha='right')

        # contourf

        levels = np.array([1, 2, 4, 10, 20, 40]) 

        cf = ax.contourf(PCO2s, Temps, CCDs, levels = levels, extend = 'both', locator=ticker.LogLocator(), 
                          cmap = plt.get_cmap('viridis')) 

        PCO2sA, TempsA = np.meshgrid(PCO2s, Temps)

        nCarb = nDIV * weath_scaling(PCO2sA, TempsA, beta=beta)

        levels = np.logspace(-2, 2, num=3)  

        ax.contour(PCO2sA, TempsA, nCarb, levels = levels, locator=ticker.LogLocator(), 
                          colors='gray')

        fig.colorbar(cf, ax=ax, label='CCD [km]')

        ax.scatter(0.3e-3, 288, color='gray', marker='o')

        ax.text(1e-3, 340, r'$n_{Mg, tot} = 100$ m$^{-3}$', c='gray',
         ha='left', va='center', rotation=-28
        )

        ax.text(1e-5, 298, r'$n_{Mg, tot} = 1$ m$^{-3}$', c='gray',
         ha='left', va='center', rotation=-28
        )

        ax.text(1e-3, 320, 'Carbon Cycle', c='black',
         ha='left', va='center'
        )

        if nSiO2 != 0:
            ax.text(2e-8, 360, 'No Carbon Cycle \n(cations consumed \nby slicates)', c='white',
             ha='left', va='center'
            )
        
        ax.text(2e-8, 280, 'No Carbon Cycle \n(too little CO$_2$)', c='white',
         ha='left', va='center'
        )

        ax.text(1e-3, 280, 'No Carbon Cycle \n(too acidic)', c='white',
         ha='left', va='center'
        )

        ax.text(5e-4, 310, r'$n_{\rm Mg, tot} = f_{\rm W}(P_{\rm CO_2}, T)$', c='gray',
         ha='left', va='center'
        )

        ax.set_title(r'Mg-CCD', fontsize=18)
        if nSiO2 == 0:
            plt.savefig('figA3b.pdf', bbox_inches='tight')
        elif nSiO2 == 1:
            plt.savefig('fig3b.pdf', bbox_inches='tight')
        # plt.savefig('CCD_Mg_beta30_nSiO2%s.pdf'%(int(nSiO2)), bbox_inches='tight')
        
    return


# In[ ]:


# Calculate Fe-CCD as a function of PCO2 and T
def FeCCD_PCO2_T(beta = 0.3, nSiO2 = 1, nDIV = 1, totnum = 10, numQ1 = 10, numQ2 = 10,
                plot_flag = True, table_flag = True):
    '''
    Returns tabulated values of Fe-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    if totnum > 10:
        print('Please be patient. A high-resolution figure is being generated.')
        
    Temps = np.linspace(273.16, 372.16, num=numQ1) # Temperature in K
    PCO2s = np.logspace(-8, -0.5, num=numQ2) # surface CO2 pressure in bar
    totPs = np.logspace(0, np.log10(5000), num=totnum)

    chems = chem_dict3(numQ1, numQ2, totnum)

    CCDs = np.zeros((numQ1, numQ2))
    
    ####################################
    
    db = SupcrtDatabase('supcrtbl')

    solution = AqueousPhase(['H2O(aq)','CO2(aq)', 'HCO3-', 'CO3-2', 'H+', 'OH-', 'Fe+2', 'SiO2(aq)'])
    solution.setActivityModel(chain(
        ActivityModelHKF(),
        ActivityModelDrummond('CO2(aq)'),
    ))

    gases = GaseousPhase(['CO2(g)', 'N2(g)'])
    gases.setActivityModel(ActivityModelPengRobinson())

    minerals = MineralPhases(['Siderite', 'Fayalite', 'Quartz'])

    system = ChemicalSystem(db, solution, gases, minerals)

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.fugacity('CO2')

    solver = EquilibriumSolver(specs)
        
    ####################################
                
    for k in range(len(Temps)):
        Temp = Temps[k] 

        for i in range(len(PCO2s)):
            PCO2 = PCO2s[i]

            addDIVtot = nDIV * weath_scaling(PCO2, Temp, beta=beta) / numden

            addSiO2 = nSiO2 * addDIVtot

            j = 0
            while j < totnum:

                totP = totPs[j]

                state = ChemicalState(system)
                state.setTemperature(Temp, 'K')
                state.setPressure(totP, 'bar')
                state.set('H2O(aq)', totH2O - addDIVtot, 'mol')     # add ~ one kg of water
                state.set('N2(g)', totN2, 'mol')
                state.set('HCO3-', 2*addDIVtot, 'mol')
                state.set('Fe+2', addDIVtot, 'mol')
                state.set('SiO2(aq)', addSiO2, 'mol')

                conditions = EquilibriumConditions(specs)
                conditions.temperature(state.temperature())
                conditions.pressure(state.pressure())
                conditions.fugacity('CO2', PCO2, 'bar')

                result = solver.solve(state, conditions)

                assert result.optima.succeeded

                chems['Fe+2'][k][i][j] = numden*state.speciesAmount('Fe+2')[0]
                chems['H+'][k][i][j] = numden*state.speciesAmount('H+')[0]
                chems['OH-'][k][i][j] = numden*state.speciesAmount('OH-')[0]
                chems['CO3-2'][k][i][j] = numden*state.speciesAmount('CO3-2')[0]
                chems['HCO3-'][k][i][j] = numden*state.speciesAmount('HCO3-')[0]
                chems['SiO2(aq)'][k][i][j] = numden*state.speciesAmount('SiO2(aq)')[0]

                chems['CO2(aq)'][k][i][j] = numden*state.speciesAmount('CO2(aq)')[0]
                chems['CO2(g)'][k][i][j] = numden*state.speciesAmount('CO2(g)')[0]
                chems['PCO2'][k][i][j] = PCO2

                chems['Siderite'][k][i][j] = numden*state.speciesAmount('Siderite')[0]
                chems['Fayalite'][k][i][j] = numden*state.speciesAmount('Fayalite')[0]
                chems['Quartz'][k][i][j] = numden*state.speciesAmount('Quartz')[0]

                chems['pH'][k][i][j] = -np.log10(state.speciesAmount('H+')[0])

                j = j + 1

    ####################################

    for k in range(len(Temps)):
        for i in range(len(PCO2s)):
            nCarb_surf = chems['Siderite'][k][i][0]
            if nCarb_surf < low_cutoff:
                CCDs[k][i] = 1e-3 # 0 # km
            else:
                CCDs[k][i] = 100 # km
            j = 0
            while j < totnum:
                if chems['Siderite'][k][i][j] < 0.001 * nCarb_surf:
                    CCDs[k][i] = ocean_depth(totPs[j])
                    j = totnum
                else:
                    j = j + 1
    
    if table_flag == True:
        df = pd.DataFrame(CCDs, index=PCO2s, columns=Temps)
        if nSiO2 == 0:
            df.to_csv('figA3c.csv')
        elif nSiO2 == 1:
            df.to_csv('fig3c.csv')
        #df.to_csv('FeCCD_PCO2_T_'+str(int(beta*100))+'_nSiO2'+str(int(nSiO2))+'.csv')
    
    if plot_flag == True:
        plot_FeCCD_PCO2_T(PCO2s, Temps, CCDs, beta = beta, nSiO2 = nSiO2, nDIV = nDIV, 
                          totnum = totnum, numQ1 = numQ1, numQ2 = numQ2)
        
    return


# In[ ]:


# Plot Fe-CCD as a function of PCO2 and T
def plot_FeCCD_PCO2_T(PCO2s, Temps, CCDs, 
                       beta = 0.3, nSiO2 = 1, nDIV = 1, totnum = 10, numQ1 = 10, numQ2 = 10):
    '''
    Returns plots of Fe-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    
    if beta == 0:
        
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.15)

        ax.set_xscale('log')
        ax.set_yscale('linear')

        ax.set_xlabel(r'$P_{\rm CO_2}$ [bar]', fontsize=16)
        ax.set_ylabel(r'$T$ [K]', fontsize=16)
        
        ax.text(-0.1, 1.01, '(c)', transform=ax.transAxes, fontsize=16, va='top', ha='right')

        # contourf

        levels = np.array([1, 2, 4, 10, 20, 40]) 

        cf = ax.contourf(PCO2s, Temps, CCDs, levels = levels, extend = 'both', locator=ticker.LogLocator(), 
                          cmap = plt.get_cmap('viridis'))

        fig.colorbar(cf, ax=ax, label='CCD [km]')

        ax.scatter(0.3e-3, 288, color='gray', marker='o')

        ax.text(1e-5, 338, 'Carbon Cycle', c='black',
         ha='left', va='center'
        )

        ax.text(1e-5, 365, r'$n_{\rm Fe,tot}$ = %.1f m$^{-3}$'%nDIV, c='gray',
         ha='left', va='center'
        )

        ax.set_title(r'Fe-CCD', fontsize=18)
        # plt.savefig('CCD_Fe_beta0_n%s_nSiO2%s.pdf'%(int(nDIV), int(nSiO2)), bbox_inches='tight')
        
    elif beta == 0.3:
        
        fig, ax = plt.subplots()
        plt.subplots_adjust(bottom=0.15)

        ax.set_xscale('log')
        ax.set_yscale('linear')

        ax.set_xlabel(r'$P_{\rm CO_2}$ [bar]', fontsize=16)
        ax.set_ylabel(r'$T$ [K]', fontsize=16)
        
        ax.text(-0.1, 1.01, '(c)', transform=ax.transAxes, fontsize=16, va='top', ha='right')

        # contourf

        levels = np.array([1, 2, 4, 10, 20, 40]) 

        cf = ax.contourf(PCO2s, Temps, CCDs, levels = levels, extend = 'both', locator=ticker.LogLocator(), 
                          cmap = plt.get_cmap('viridis')) 

        PCO2sA, TempsA = np.meshgrid(PCO2s, Temps)

        nCarb = nDIV * weath_scaling(PCO2sA, TempsA, beta=beta)

        levels = np.logspace(-2, 2, num=3)  

        ax.contour(PCO2sA, TempsA, nCarb, levels = levels, locator=ticker.LogLocator(), 
                          colors='gray')

        fig.colorbar(cf, ax=ax, label='CCD [km]')

        ax.scatter(0.3e-3, 288, color='gray', marker='o')

        ax.text(1e-3, 340, r'$n_{Fe, tot} = 100$ m$^{-3}$', c='gray',
         ha='left', va='center', rotation=-28
        )

        ax.text(1e-5, 298, r'$n_{Fe, tot} = 1$ m$^{-3}$', c='gray',
         ha='left', va='center', rotation=-28
        )

        ax.text(1e-3, 320, 'Carbon Cycle', c='black',
         ha='left', va='center'
        )

        if nSiO2 != 0:
            ax.text(2e-8, 360, 'No Carbon Cycle \n(cations consumed \nby slicates)', c='white',
             ha='left', va='center'
            )

        ax.text(5e-4, 310, r'$n_{\rm Fe, tot} = f_{\rm W}(P_{\rm CO_2}, T)$', c='gray',
         ha='left', va='center'
        )

        ax.set_title(r'Fe-CCD', fontsize=18)
        if nSiO2 == 0:
            plt.savefig('figA3c.pdf', bbox_inches='tight')
        elif nSiO2 == 1:
            plt.savefig('fig3c.pdf', bbox_inches='tight')
        # plt.savefig('CCD_Fe_beta30_n%s_nSiO2%s.pdf'%(int(nDIV), int(nSiO2)), bbox_inches='tight')
        
    return


# In[ ]:


# Calculate and plot pH as a function of PCO2 for Ca, Mg or Fe systems
def plot_pH(DIV = 'Ca', totP = 1, Temp = 288, nDIV = 1, nSiO2 = 1, totnum = 100,
           plot_flag = True, table_flag = True):
    '''
    Returns tabulated values and plots of ocean pH as a function of PCO2 [bar] for Ca, Mg or Fe carbonate systems
    '''
        
    if DIV != 'Ca' and DIV != 'Mg' and DIV != 'Fe':
        print('Error: Enter DIV = "Ca" or "Mg" or "Fe"')

    PCO2s = np.logspace(-8, -0.5, num=totnum) # bar
    
    ####################################
    ####################################
    
    if DIV == 'Ca':

        numQ = 3
        betas = np.array([-1, 0, 0.3]) 

        chems2 = chem_dict2(numQ, totnum)

        db = SupcrtDatabase('supcrtbl')

        solution = AqueousPhase(['H2O(aq)','CO2(aq)', 'HCO3-', 'CO3-2', 'H+', 'OH-', 'Ca+2', 'SiO2(aq)'])
        solution.setActivityModel(chain(
            ActivityModelHKF(),
            ActivityModelDrummond('CO2(aq)'),
        ))

        gases = GaseousPhase(['CO2(g)', 'N2(g)'])
        gases.setActivityModel(ActivityModelPengRobinson())

        minerals = MineralPhases(['Calcite', 'Wollastonite', 'Quartz'])

        system = ChemicalSystem(db, solution, gases, minerals)

        specs = EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()
        specs.fugacity('CO2')

        solver = EquilibriumSolver(specs)

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

                state = ChemicalState(system)
                state.setTemperature(Temp, 'K')
                state.setPressure(totP, 'bar')
                state.set('H2O(aq)', totH2O - addDIVtot, 'mol')     # add ~ one kg of water
                state.set('N2(g)', totN2, 'mol')
                state.set('HCO3-', 2*addDIVtot, 'mol')
                state.set('Ca+2', addDIVtot, 'mol')
                state.set('SiO2(aq)', addSiO2, 'mol')

                conditions = EquilibriumConditions(specs)
                conditions.temperature(state.temperature())
                conditions.pressure(state.pressure())
                conditions.fugacity('CO2', PCO2, 'bar')

                result = solver.solve(state, conditions)

                assert result.optima.succeeded

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

                j = j + 1

        ####################################
        
        if table_flag == True:
            
            df = pd.DataFrame({
                'pH lower limit': chems2['pH'][0],
                'pH upper limit': chems2['pH'][1],
                'pH weathering': chems2['pH'][2],
            }, index=PCO2s)
            df.to_csv('fig2a.csv')
            #df.to_csv('pH_PCO2_Ca.csv')
            
        ####################################
            
        if plot_flag == True:

            fig1 = plt.figure(constrained_layout=False,figsize=(5,5))
            spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)

            p11  = fig1.add_subplot(spec1[0,0])
            
            p11.text(-0.1, 1.01, '(a)', transform=p11.transAxes, fontsize=16, va='top', ha='right')

            p11.fill_between(PCO2s, chems2['pH'][1], chems2['pH'][0], facecolor='blue', alpha=0.3)
            p11.fill_between(PCO2s, chems2['pH'][0], 4, facecolor='red', alpha=0.1)
            p11.fill_between(PCO2s, 11, chems2['pH'][1], facecolor='red', alpha=0.1)
            p11.plot(PCO2s, chems2['pH'][1], lw=5,c='black',ls='-')
            p11.plot(PCO2s, chems2['pH'][2], lw=3,c='cyan',ls=':', 
                     label=r'$n_{\rm Ca,tot}$ = $f_{\rm W} (P_{\rm CO_2})$')

            p11.scatter(0.3e-3, 8.1, color='green', marker='x')

            p11.set_ylim([4,11])

            p11.set_xscale('log')

            p11.set_xlabel(r'$P_{\rm CO_2}$ [bar]',fontsize=16)
            p11.set_ylabel(r'Ocean pH', fontsize=18)
            p11.tick_params(axis='y')

            p11.text(1e-6, 9.6, 'Carbon Cycle', c='black',
                     ha='left', va='center', rotation=-35
                    )

            p11.text(5e-7, 9, 'No Carbon Cycle', c='blue',
                     ha='left', va='center', rotation=-35
                    )

            p11.text(3e-5, 7.7, 'Modern \nEarth pH', c='green',
                     ha='left', va='center', rotation=0
                    )

            p11.text(1e-4, 10, 'Forbidden', c='red',
                     ha='left', va='center', rotation=0
                    )

            p11.text(1e-7, 5, 'Forbidden', c='red',
                     ha='left', va='center', rotation=0
                    )

            p11.legend(fontsize=10, bbox_to_anchor=(0.05,0.3), borderaxespad=0, loc='lower left')

            p11.set_title(r'Ca', fontsize=16)
            plt.savefig('fig2a.pdf', bbox_inches='tight')
            # plt.savefig('pH_PCO2_Ca.pdf', bbox_inches='tight')
    
    ####################################
    ####################################
    
    elif DIV == 'Mg':

        numQ = 3
        betas = np.array([-1, 0, 0.3]) 

        chems2 = chem_dict2(numQ, totnum)

        db = SupcrtDatabase('supcrtbl')

        solution = AqueousPhase(['H2O(aq)','CO2(aq)', 'HCO3-', 'CO3-2', 'H+', 'OH-', 'Mg+2', 'SiO2(aq)'])
        solution.setActivityModel(chain(
            ActivityModelHKF(),
            ActivityModelDrummond('CO2(aq)'),
        ))

        gases = GaseousPhase(['CO2(g)', 'N2(g)'])
        gases.setActivityModel(ActivityModelPengRobinson())

        minerals = MineralPhases(['Magnesite', 'Clino-Enstatite', 'Quartz'])

        system = ChemicalSystem(db, solution, gases, minerals)

        specs = EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()
        specs.fugacity('CO2')

        solver = EquilibriumSolver(specs)

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

                state = ChemicalState(system)
                state.setTemperature(Temp, 'K')
                state.setPressure(totP, 'bar')
                state.set('H2O(aq)', totH2O - addDIVtot, 'mol')     # add ~ one kg of water
                state.set('N2(g)', totN2, 'mol')
                state.set('HCO3-', 2*addDIVtot, 'mol')
                state.set('Mg+2', addDIVtot, 'mol')
                state.set('SiO2(aq)', addSiO2, 'mol')

                conditions = EquilibriumConditions(specs)
                conditions.temperature(state.temperature())
                conditions.pressure(state.pressure())
                conditions.fugacity('CO2', PCO2, 'bar')

                result = solver.solve(state, conditions)

                assert result.optima.succeeded

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

                j = j + 1

        ####################################
        
        if table_flag == True:
            
            df = pd.DataFrame({
                'pH lower limit': chems2['pH'][0],
                'pH upper limit': chems2['pH'][1],
                'pH weathering': chems2['pH'][2],
            }, index=PCO2s)
            df.to_csv('fig2b.csv')
            #df.to_csv('pH_PCO2_Mg.csv')
            
        ####################################
        
        if plot_flag == True:

            fig1 = plt.figure(constrained_layout=False,figsize=(5,5))
            spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)

            p11  = fig1.add_subplot(spec1[0,0])
            
            p11.text(-0.1, 1.01, '(b)', transform=p11.transAxes, fontsize=16, va='top', ha='right')

            p11.fill_between(PCO2s, chems2['pH'][1], chems2['pH'][0], facecolor='blue', alpha=0.3)
            p11.fill_between(PCO2s, chems2['pH'][0], 4, facecolor='red', alpha=0.1)
            p11.fill_between(PCO2s, 11, chems2['pH'][1], facecolor='red', alpha=0.1)
            p11.plot(PCO2s, chems2['pH'][1], lw=5,c='black',ls='-')
            p11.plot(PCO2s, chems2['pH'][2], lw=3,c='cyan',ls=':', 
                     label=r'$n_{\rm Mg,tot}$ = $f_{\rm W} (P_{\rm CO_2})$')

            p11.scatter(0.3e-3, 8.1, color='green', marker='x')

            p11.set_ylim([4,11])

            p11.set_xscale('log')

            p11.set_xlabel(r'$P_{\rm CO_2}$ [bar]',fontsize=16)
            p11.set_ylabel(r'Ocean pH', fontsize=18)
            p11.tick_params(axis='y')

            p11.text(1e-6, 9.8, 'Carbon Cycle', c='black',
                     ha='left', va='center', rotation=-35
                    )

            p11.text(5e-7, 9, 'No Carbon Cycle', c='blue',
                     ha='left', va='center', rotation=-35
                    )

            p11.text(3e-5, 7.7, 'Modern \nEarth pH', c='green',
                     ha='left', va='center', rotation=0
                    )

            p11.text(1e-4, 10, 'Forbidden', c='red',
                     ha='left', va='center', rotation=0
                    )

            p11.text(1e-7, 5, 'Forbidden', c='red',
                     ha='left', va='center', rotation=0
                    )

            p11.legend(fontsize=10, bbox_to_anchor=(0.05,0.3), borderaxespad=0, loc='lower left')

            p11.set_title(r'Mg', fontsize=16)
            plt.savefig('fig2b.pdf', bbox_inches='tight')
            # plt.savefig('pH_PCO2_Mg.pdf', bbox_inches='tight')
        
    ####################################
    ####################################
    
    elif DIV == 'Fe':

        numQ = 3
        betas = np.array([-1, 0, 0.3]) 

        chems2 = chem_dict2(numQ, totnum)

        db = SupcrtDatabase('supcrtbl')

        solution = AqueousPhase(['H2O(aq)','CO2(aq)', 'HCO3-', 'CO3-2', 'H+', 'OH-', 'Fe+2', 'SiO2(aq)'])
        solution.setActivityModel(chain(
            ActivityModelHKF(),
            ActivityModelDrummond('CO2(aq)'),
        ))

        gases = GaseousPhase(['CO2(g)', 'N2(g)'])
        gases.setActivityModel(ActivityModelPengRobinson())

        minerals = MineralPhases(['Siderite', 'Fayalite', 'Quartz'])
        
        system = ChemicalSystem(db, solution, gases, minerals)

        specs = EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()
        specs.fugacity('CO2')

        solver = EquilibriumSolver(specs)

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

                state = ChemicalState(system)
                state.setTemperature(Temp, 'K')
                state.setPressure(totP, 'bar')
                state.set('H2O(aq)', totH2O - addDIVtot, 'mol')     # add ~ one kg of water
                state.set('N2(g)', totN2, 'mol')
                state.set('HCO3-', 2*addDIVtot, 'mol')
                state.set('Fe+2', addDIVtot, 'mol')
                state.set('SiO2(aq)', addSiO2, 'mol')

                conditions = EquilibriumConditions(specs)
                conditions.temperature(state.temperature())
                conditions.pressure(state.pressure())
                conditions.fugacity('CO2', PCO2, 'bar')

                result = solver.solve(state, conditions)

                assert result.optima.succeeded

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

                j = j + 1

        ####################################
        
        if table_flag == True:
            
            df = pd.DataFrame({
                'pH lower limit': chems2['pH'][0],
                'pH upper limit': chems2['pH'][1],
                'pH weathering': chems2['pH'][2],
            }, index=PCO2s)
            df.to_csv('fig2c.csv')
            #df.to_csv('pH_PCO2_Fe.csv')
            
        ####################################
        
        if plot_flag == True:
        
            fig1 = plt.figure(constrained_layout=False,figsize=(5,5))
            spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)

            p11  = fig1.add_subplot(spec1[0,0])
            
            p11.text(-0.1, 1.01, '(c)', transform=p11.transAxes, fontsize=16, va='top', ha='right')

            p11.fill_between(PCO2s, chems2['pH'][1], chems2['pH'][0], facecolor='blue', alpha=0.3)
            p11.fill_between(PCO2s, chems2['pH'][0], 4, facecolor='red', alpha=0.1)
            p11.fill_between(PCO2s, 11, chems2['pH'][1], facecolor='red', alpha=0.1)
            p11.plot(PCO2s, chems2['pH'][1], lw=5,c='black',ls='-')
            p11.plot(PCO2s, chems2['pH'][2], lw=3,c='cyan',ls=':', 
                     label=r'$n_{\rm Fe,tot}$ = $f_{\rm W} (P_{\rm CO_2})$')

            p11.scatter(0.3e-3, 8.1, color='green', marker='x')

            p11.set_ylim([4,11])

            p11.set_xscale('log')

            p11.set_xlabel(r'$P_{\rm CO_2}$ [bar]',fontsize=16)
            p11.set_ylabel(r'Ocean pH', fontsize=18)
            p11.tick_params(axis='y')

            p11.text(1e-6, 8.5, 'Carbon Cycle', c='black',
                     ha='left', va='center', rotation=-35
                    )

            p11.text(1e-7, 8.2, 'No Carbon Cycle', c='blue',
                     ha='left', va='center', rotation=-35
                    )

            p11.text(1e-4, 7.7, 'Modern \nEarth pH', c='green',
                     ha='left', va='center', rotation=0
                    )

            p11.text(1e-4, 10, 'Forbidden', c='red',
                     ha='left', va='center', rotation=0
                    )

            p11.text(1e-7, 5, 'Forbidden', c='red',
                     ha='left', va='center', rotation=0
                    )

            p11.legend(fontsize=10, bbox_to_anchor=(0.05,0.3), borderaxespad=0, loc='lower left')

            p11.set_title(r'Fe', fontsize=16)
            plt.savefig('fig2c.pdf', bbox_inches='tight')
            # plt.savefig('pH_PCO2_Fe.pdf', bbox_inches='tight')

    return


# In[ ]:


# Calculate and plot analytical and numerical solutions of ocean pH
def plot_pH_an(DIV = 'Ca', totP = 1, Temp = 288, nDIV = 1, nSiO2 = 1, totnum = 100, nDIV_fixed = 1,
           plot_flag = True, table_flag = True):
    '''
    Returns tabulated values and plots of analytical/numerical solutions of ocean pH as a function of PCO2 [bar]
    '''    
    if DIV != 'Ca':
        print('Error: Enter DIV = "Ca"')

    PCO2s = np.logspace(-8, -0.5, num=totnum) # bar
    
    ####################################
    ####################################
    
    if DIV == 'Ca':
        
        numQ = 2
        addDIVtots = numden*np.linspace(0, 1e-2, num=numQ) # in units of 1000 * mol/dm3 = 1 mol/m3
        pH_limits = chem_dict2(numQ,totnum)
        pH_limits_san = chem_dict2(numQ,totnum)
        pH_limits_an = chem_dict2(numQ,totnum)

        db = SupcrtDatabase('supcrtbl')

        solution = AqueousPhase(['H2O(aq)','CO2(aq)', 'HCO3-', 'CO3-2', 'H+', 'OH-', 'Ca+2', 'SiO2(aq)'])
        solution.setActivityModel(chain(
            ActivityModelHKF(),
            ActivityModelDrummond('CO2(aq)'),
        ))

        gases = GaseousPhase(['CO2(g)', 'N2(g)'])
        gases.setActivityModel(ActivityModelPengRobinson())

        minerals = MineralPhases(['Calcite', 'Wollastonite', 'Quartz'])

        system = ChemicalSystem(db, solution, gases, minerals)

        specs = EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()
        specs.fugacity('CO2')

        solver = EquilibriumSolver(specs)
        
        rxn9 = db.reaction('Ca+2 + CO3-2 = Calcite')
        logK9 = rxn9.props(Temp, 'K', totP, 'bar').lgK[0]
        
        rxn15 = db.reaction('CO2(g) + H2O(aq) = 2*H+ + CO3-2')
        logK16 = rxn15.props(Temp, 'K', totP, 'bar').lgK[0]
        
        rxn3= db.reaction('CO2(g) + H2O(aq) = H+ + HCO3-')
        logK3 = rxn3.props(Temp, 'K', totP, 'bar').lgK[0]
        
        for i in range(numQ):
            j = 0
            addDIVtot = addDIVtots[i] / numden
            addSiO2   = 0

            while j < totnum:

                PCO2 = PCO2s[j]

                state = ChemicalState(system)
                state.setTemperature(Temp, 'K')
                state.setPressure(totP, 'bar')
                state.set('H2O(aq)', totH2O - addDIVtot, 'mol')     # add ~ one kg of water
                state.set('N2(g)', totN2, 'mol')
                state.set('HCO3-', 2*addDIVtot, 'mol')
                state.set('Ca+2', addDIVtot, 'mol')
                state.set('SiO2(aq)', addSiO2, 'mol')

                conditions = EquilibriumConditions(specs)
                conditions.temperature(state.temperature())
                conditions.pressure(state.pressure())
                conditions.fugacity('CO2', PCO2, 'bar')

                result = solver.solve(state, conditions)

                assert result.optima.succeeded

                pH_limits['Ca+2'][i][j] = numden*state.speciesAmount('Ca+2')[0]
                pH_limits['H+'][i][j] = numden*state.speciesAmount('H+')[0]
                pH_limits['OH-'][i][j] = numden*state.speciesAmount('OH-')[0]
                pH_limits['CO3-2'][i][j] = numden*state.speciesAmount('CO3-2')[0]
                pH_limits['HCO3-'][i][j] = numden*state.speciesAmount('HCO3-')[0]
                pH_limits['SiO2(aq)'][i][j] = numden*state.speciesAmount('SiO2(aq)')[0]

                pH_limits['CO2(aq)'][i][j] = numden*state.speciesAmount('CO2(aq)')[0]
                pH_limits['CO2(g)'][i][j] = numden*state.speciesAmount('CO2(g)')[0]
                pH_limits['PCO2'][i][j] = PCO2
                pH_limits['Calcite'][i][j] = numden*state.speciesAmount('Calcite')[0]
                pH_limits['Wollastonite'][i][j] = numden*state.speciesAmount('Wollastonite')[0]
                pH_limits['Quartz'][i][j] = numden*state.speciesAmount('Quartz')[0]

                pH_limits['pH'][i][j] = -np.log10(state.speciesAmount('H+')[0])
                
                if i == 0: 
                    pH_limits_an['pH'][i][j] = ocean_pH_low(PCO2, logK3)
                else:
                    pH_limits_san['pH'][i][j] = ocean_pH_upp(PCO2, pH_limits['Ca+2'][i][j], logK9, logK16)
                    pH_limits_an['pH'][i][j] = ocean_pH_upp(PCO2, nDIV_fixed, logK9, logK16)

                j = j + 1

        ####################################
        
        if table_flag == True:
            
            df = pd.DataFrame({
                'pH lower numerical': pH_limits['pH'][0],
                'pH lower analytical': pH_limits_an['pH'][0],
                'pH upper numerical': pH_limits['pH'][1],
                'pH upper analytical': pH_limits_an['pH'][1],
                'pH upper semi-analy.': pH_limits_san['pH'][1],
            }, index=PCO2s)
            df.to_csv('figA1.csv')
            #df.to_csv('pHan_PCO2_Ca.csv')
            
        ####################################
        
        if plot_flag == True:

            fig1 = plt.figure(constrained_layout=False,figsize=(5,5))
            spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)

            p11  = fig1.add_subplot(spec1[0,0])

            p11.fill_between(PCO2s, pH_limits['pH'][1], pH_limits['pH'][0], facecolor='blue', alpha=0.3)
            p11.fill_between(PCO2s, pH_limits['pH'][0], 4, facecolor='red', alpha=0.1)
            p11.fill_between(PCO2s, 11, pH_limits['pH'][1], facecolor='red', alpha=0.1)
            p11.plot(PCO2s, pH_limits['pH'][1], lw=5,c='black',ls='-', label = 'Up (numerical)')
            p11.plot(PCO2s, pH_limits_san['pH'][1], lw=3,c='orange',ls=':', label = 'Up (semi-analytical)')
            p11.plot(PCO2s, pH_limits_an['pH'][1], lw=3,c='orange',ls='--', 
                     label = r'Up (ana., $n_{\rm Ca^{2+}}$ = %d m$^{-3}$)'%nDIV_fixed)

            p11.plot(PCO2s, pH_limits['pH'][0], lw=5,c='blue',ls='-', label = 'Low (numerical)')
            p11.plot(PCO2s, pH_limits_an['pH'][0], lw=3,c='cyan',ls='--', label = 'Low (analytical)')

            p11.scatter(0.3e-3, 8.1, color='green', marker='x')

            p11.set_ylim([4,11])

            p11.set_xscale('log')

            p11.set_xlabel(r'$P_{\rm CO_2}$ [bar]',fontsize=16)
            p11.set_ylabel(r'Ocean pH', fontsize=18)
            p11.tick_params(axis='y')

            p11.text(1e-6, 9.6, 'Carbon Cycle', c='black',
                     ha='left', va='center', rotation=-35
                    )

            p11.text(5e-7, 9, 'No Carbon Cycle', c='blue',
                     ha='left', va='center', rotation=-35
                    )

            p11.text(3e-5, 7.7, 'Modern \nEarth pH', c='green',
                     ha='left', va='center', rotation=0
                    )

            p11.text(1e-4, 10, 'Forbidden', c='red',
                     ha='left', va='center', rotation=0
                    )

            p11.legend(fontsize=10, bbox_to_anchor=(0.05,0.05), borderaxespad=0, loc='lower left')

            p11.set_title(r'Ca', fontsize=16)
            plt.savefig('figA1.pdf', bbox_inches='tight')
            # plt.savefig('pHan_PCO2_Ca.pdf', bbox_inches='tight')
    
    return


# In[ ]:


# Calculate and plot ocean pH as a function of local pressure
def plot_pH_P(DIV = 'Ca', PCO2 = 0.3e-3, Temp = 288, nDIV = 1, nSiO2 = 1, totnum = 100,
           plot_flag = True, table_flag = True):
    '''
    Returns tabulated values and plots of ocean pH as a function of P [bar]
    '''  
    if DIV != 'Ca':
        print('Error: Enter DIV = "Ca"')

    totPs = np.logspace(0, 3, num=totnum) # bar
    
    ####################################
    ####################################
    
    if DIV == 'Ca':
        
        numQ = 3
        betas = np.array([-1, 0, 0.3]) 

        chems2 = chem_dict2(numQ, totnum)

        db = SupcrtDatabase('supcrtbl')

        solution = AqueousPhase(['H2O(aq)','CO2(aq)', 'HCO3-', 'CO3-2', 'H+', 'OH-', 'Ca+2', 'SiO2(aq)'])
        solution.setActivityModel(chain(
            ActivityModelHKF(),
            ActivityModelDrummond('CO2(aq)'),
        ))

        gases = GaseousPhase(['CO2(g)', 'N2(g)'])
        gases.setActivityModel(ActivityModelPengRobinson())

        minerals = MineralPhases(['Calcite', 'Wollastonite', 'Quartz'])

        system = ChemicalSystem(db, solution, gases, minerals)

        specs = EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()
        specs.fugacity('CO2')

        solver = EquilibriumSolver(specs)

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

                state = ChemicalState(system)
                state.setTemperature(Temp, 'K')
                state.setPressure(totP, 'bar')
                state.set('H2O(aq)', totH2O - addDIVtot, 'mol')     # add ~ one kg of water
                state.set('N2(g)', totN2, 'mol')
                state.set('HCO3-', 2*addDIVtot, 'mol')
                state.set('Ca+2', addDIVtot, 'mol')
                state.set('SiO2(aq)', addSiO2, 'mol')

                conditions = EquilibriumConditions(specs)
                conditions.temperature(state.temperature())
                conditions.pressure(state.pressure())
                conditions.fugacity('CO2', PCO2, 'bar')

                result = solver.solve(state, conditions)

                assert result.optima.succeeded

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

                j = j + 1
                
        ####################################
        
        if table_flag == True:
            
            df = pd.DataFrame({
                'pH lower limit': chems2['pH'][0],
                'pH upper limit': chems2['pH'][1],
                'pH weathering': chems2['pH'][2],
            }, index=totPs)
            df.to_csv('figA2a.csv')
            #df.to_csv('pH_P_Ca.csv')
            
        ####################################
        
        if plot_flag == True:

            fig1 = plt.figure(constrained_layout=False,figsize=(5,5))
            spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)

            p11  = fig1.add_subplot(spec1[0,0])
            
            p11.text(-0.1, 1.01, '(a)', transform=p11.transAxes, fontsize=16, va='top', ha='right')

            p11.fill_between(totPs, chems2['pH'][1], chems2['pH'][0], facecolor='blue', alpha=0.3)
            p11.fill_between(totPs, chems2['pH'][0], 4, facecolor='red', alpha=0.1)
            p11.fill_between(totPs, 11, chems2['pH'][1], facecolor='red', alpha=0.1)
            p11.plot(totPs, chems2['pH'][1], lw=5,c='black',ls='-')
            p11.plot(totPs, chems2['pH'][2], lw=3,c='cyan',ls=':', 
                     label=r'$n_{\rm Ca,tot}$ = $f_{\rm W} (P_{\rm CO_2})$')

            p11.scatter(1, 8.1, color='green', marker='x')

            p11.set_ylim([4,11])

            p11.set_xscale('log')

            p11.set_xlabel(r'$P$ [bar]',fontsize=16)
            p11.set_ylabel(r'Ocean pH', fontsize=18)
            p11.tick_params(axis='y')

            p11.text(10, 8.5, 'Carbon Cycle', c='black',
                     ha='left', va='center', rotation=0
                    )

            p11.text(10, 7.8, 'No Carbon Cycle', c='blue',
                     ha='left', va='center', rotation=0
                    )

            p11.text(1, 7.8, 'Modern \nEarth pH', c='green',
                     ha='left', va='center', rotation=0
                    )

            p11.text(10, 10, 'Forbidden', c='red',
                     ha='left', va='center', rotation=0
                    )

            p11.text(10, 5, 'Forbidden', c='red',
                     ha='left', va='center', rotation=0
                    )

            p11.legend(fontsize=10, bbox_to_anchor=(0.05,0.3), borderaxespad=0, loc='lower left')

            p11.set_title(r'Ca', fontsize=16)
            plt.savefig('figA2a.pdf', bbox_inches='tight')
            # plt.savefig('pH_P_Ca.pdf', bbox_inches='tight')

    return


# In[ ]:


# Calculate and plot ocean pH as a function of temperature
def plot_pH_T(DIV = 'Ca', PCO2 = 0.3e-3, totP = 1, nDIV = 1, nSiO2 = 1, totnum = 100,
           plot_flag = True, table_flag = True):
    '''
    Returns tabulated values and plots of ocean pH as a function of Temp [K]
    '''  
    if DIV != 'Ca':
        print('Error: Enter DIV = "Ca"')

    Temps = np.linspace(273.16, 372.16, num=totnum) # bar
    
    ####################################
    ####################################
    
    if DIV == 'Ca':

        numQ = 3
        betas = np.array([-1, 0, 0.3]) 

        chems2 = chem_dict2(numQ, totnum)

        db = SupcrtDatabase('supcrtbl')

        solution = AqueousPhase(['H2O(aq)','CO2(aq)', 'HCO3-', 'CO3-2', 'H+', 'OH-', 'Ca+2', 'SiO2(aq)'])
        solution.setActivityModel(chain(
            ActivityModelHKF(),
            ActivityModelDrummond('CO2(aq)'),
        ))

        gases = GaseousPhase(['CO2(g)', 'N2(g)'])
        gases.setActivityModel(ActivityModelPengRobinson())

        minerals = MineralPhases(['Calcite', 'Wollastonite', 'Quartz'])

        system = ChemicalSystem(db, solution, gases, minerals)

        specs = EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()
        specs.fugacity('CO2')

        solver = EquilibriumSolver(specs)

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

                state = ChemicalState(system)
                state.setTemperature(Temp, 'K')
                state.setPressure(totP, 'bar')
                state.set('H2O(aq)', totH2O - addDIVtot, 'mol')     # add ~ one kg of water
                state.set('N2(g)', totN2, 'mol')
                state.set('HCO3-', 2*addDIVtot, 'mol')
                state.set('Ca+2', addDIVtot, 'mol')
                state.set('SiO2(aq)', addSiO2, 'mol')

                conditions = EquilibriumConditions(specs)
                conditions.temperature(state.temperature())
                conditions.pressure(state.pressure())
                conditions.fugacity('CO2', PCO2, 'bar')

                result = solver.solve(state, conditions)

                assert result.optima.succeeded

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

                j = j + 1

        ####################################
        
        if table_flag == True:
            
            df = pd.DataFrame({
                'pH lower limit': chems2['pH'][0],
                'pH upper limit': chems2['pH'][1],
                'pH weathering': chems2['pH'][2],
            }, index=Temps)
            df.to_csv('figA2b.csv')
            #df.to_csv('pH_T_Ca.csv')
            
        ####################################
        
        if plot_flag == True:

            fig1 = plt.figure(constrained_layout=False,figsize=(5,5))
            spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)

            p11  = fig1.add_subplot(spec1[0,0])
            
            p11.text(-0.1, 1.01, '(b)', transform=p11.transAxes, fontsize=16, va='top', ha='right')

            p11.fill_between(Temps, chems2['pH'][1], chems2['pH'][0], facecolor='blue', alpha=0.3)
            p11.fill_between(Temps, chems2['pH'][0], 4, facecolor='red', alpha=0.1)
            p11.fill_between(Temps, 11, chems2['pH'][1], facecolor='red', alpha=0.1)
            p11.plot(Temps, chems2['pH'][1], lw=5,c='black',ls='-')
            p11.plot(Temps, chems2['pH'][2], lw=3,c='cyan',ls=':', 
                     label=r'$n_{\rm Ca,tot}$ = $f_{\rm W} (P_{\rm CO_2})$')

            p11.scatter(288, 8.1, color='green', marker='x')

            p11.set_ylim([4,11])

            p11.set_xlabel(r'$T$ [K]',fontsize=16)
            p11.set_ylabel(r'Ocean pH', fontsize=18)
            p11.tick_params(axis='y')

            p11.text(320, 8.5, 'Carbon Cycle', c='black',
                     ha='left', va='center', rotation=0
                    )

            p11.text(320, 7.8, 'No Carbon Cycle', c='blue',
                     ha='left', va='center', rotation=0
                    )

            p11.text(285, 7.7, 'Modern \nEarth pH', c='green',
                     ha='left', va='center', rotation=0
                    )

            p11.text(280, 10, 'Forbidden', c='red',
                     ha='left', va='center', rotation=0
                    )

            p11.text(280, 5, 'Forbidden', c='red',
                     ha='left', va='center', rotation=0
                    )

            p11.legend(fontsize=10, bbox_to_anchor=(0.05,0.3), borderaxespad=0, loc='lower left')

            p11.set_title(r'Ca', fontsize=16)
            plt.savefig('figA2b.pdf', bbox_inches='tight')
            # plt.savefig('pH_T_Ca.pdf', bbox_inches='tight')
    
    return


# In[ ]:


# Calculate stable phases as a function of PCO2
def phases_PCO2 (DIV = 'Ca', Temp = 298, totP = 1, beta = 0.3, nDIV = 1, nSiO2 = 1, totnum = 100,
                table_flag = True, plot_flag = True):
    '''
    Returns tabulated values of stable phases as a function of PCO2 [bar]
    '''  
    if DIV != 'Ca' and DIV != 'Mg' and DIV != 'Fe':
        print('Error: Enter DIV = "Ca" or "Mg" or "Fe"')
    
    PCO2s = np.logspace(-8, -0.5, num=totnum) # bar
    chems = chem_dict1(totnum)

    ####################################
    ####################################
    
    if DIV == 'Ca':
        
        db = SupcrtDatabase('supcrtbl')

        solution = AqueousPhase(['H2O(aq)','CO2(aq)', 'HCO3-', 'CO3-2', 'H+', 'OH-', 'Ca+2', 'SiO2(aq)'])
        solution.setActivityModel(chain(
            ActivityModelHKF(),
            ActivityModelDrummond('CO2(aq)'),
        ))

        gases = GaseousPhase(['CO2(g)', 'N2(g)'])
        gases.setActivityModel(ActivityModelPengRobinson())

        minerals = MineralPhases(['Calcite', 'Wollastonite', 'Quartz'])

        system = ChemicalSystem(db, solution, gases, minerals)

        specs = EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()
        specs.fugacity('CO2')

        solver = EquilibriumSolver(specs)

        j = 0
        while j < totnum:

            PCO2 = PCO2s[j]
            
            addDIVtot = nDIV * weath_scaling(PCO2,Temp,beta=beta) / numden
            addSiO2 = nSiO2 * addDIVtot

            state = ChemicalState(system)
            state.setTemperature(Temp, 'K')
            state.setPressure(totP, 'bar')
            state.set('H2O(aq)', totH2O - addDIVtot, 'mol')     # add ~ one kg of water
            state.set('N2(g)', totN2, 'mol')
            state.set('HCO3-', 2*addDIVtot, 'mol')
            state.set('Ca+2', addDIVtot, 'mol')
            state.set('SiO2(aq)', addSiO2, 'mol')

            conditions = EquilibriumConditions(specs)
            conditions.temperature(state.temperature())
            conditions.pressure(state.pressure())
            conditions.fugacity('CO2', PCO2, 'bar')

            result = solver.solve(state, conditions)

#            assert result.optima.succeeded

            chems['Ca+2'][j] = numden*state.speciesAmount('Ca+2')[0]
            chems['H+'][j] = numden*state.speciesAmount('H+')[0]
            chems['OH-'][j] = numden*state.speciesAmount('OH-')[0]
            chems['CO3-2'][j] = numden*state.speciesAmount('CO3-2')[0]
            chems['HCO3-'][j] = numden*state.speciesAmount('HCO3-')[0]
            chems['SiO2(aq)'][j] = numden*state.speciesAmount('SiO2(aq)')[0]

            chems['CO2(aq)'][j] = numden*state.speciesAmount('CO2(aq)')[0]
            chems['CO2(g)'][j] = numden*state.speciesAmount('CO2(g)')[0]
            chems['PCO2'][j] = PCO2

            chems['Calcite'][j] = numden*state.speciesAmount('Calcite')[0]
            chems['Wollastonite'][j] = numden*state.speciesAmount('Wollastonite')[0]

            chems['pH'][j] = -np.log10(state.speciesAmount('H+')[0])


            j = j + 1
    
        df = pd.DataFrame({
            'Ca++': chems['Ca+2'],
            'Calcite': chems['Calcite'],
            'Silicates': chems['Wollastonite'],
        }, index=PCO2s)
        
        if table_flag == True:
            if nSiO2 == 0:
                df.to_csv('figA4a.csv')
            elif nSiO2 == 1:
                df.to_csv('fig4a.csv')
            #df.to_csv('phases_PCO2_'+DIV+'_beta'+str(int(beta*100))+'_nSiO2'+str(int(nSiO2))+'.csv')
            
        if plot_flag == True:
            plot_phases_PCO2(df, DIV = DIV, Temp = Temp, totP = totP,  beta = beta, 
                      nDIV = nDIV, nSiO2 = nSiO2, totnum = totnum)
        
    ####################################
    ####################################
    
    elif DIV == 'Mg':
        
        db = SupcrtDatabase('supcrtbl')

        solution = AqueousPhase(['H2O(aq)','CO2(aq)', 'HCO3-', 'CO3-2', 'H+', 'OH-', 'Mg+2', 'SiO2(aq)'])
        solution.setActivityModel(chain(
            ActivityModelHKF(),
            ActivityModelDrummond('CO2(aq)'),
        ))

        gases = GaseousPhase(['CO2(g)', 'N2(g)'])
        gases.setActivityModel(ActivityModelPengRobinson())

        minerals = MineralPhases(['Magnesite', 'Clino-Enstatite', 'Quartz'])

        system = ChemicalSystem(db, solution, gases, minerals)

        specs = EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()
        specs.fugacity('CO2')

        solver = EquilibriumSolver(specs)

        j = 0
        while j < totnum:

            PCO2 = PCO2s[j]
            
            addDIVtot = nDIV * weath_scaling(PCO2,Temp,beta=beta) / numden
            addSiO2 = nSiO2 * addDIVtot

            state = ChemicalState(system)
            state.setTemperature(Temp, 'K')
            state.setPressure(totP, 'bar')
            state.set('H2O(aq)', totH2O - addDIVtot, 'mol')     # add ~ one kg of water
            state.set('N2(g)', totN2, 'mol')
            state.set('HCO3-', 2*addDIVtot, 'mol')
            state.set('Mg+2', addDIVtot, 'mol')
            state.set('SiO2(aq)', addSiO2, 'mol')

            conditions = EquilibriumConditions(specs)
            conditions.temperature(state.temperature())
            conditions.pressure(state.pressure())
            conditions.fugacity('CO2', PCO2, 'bar')

            result = solver.solve(state, conditions)

#            assert result.optima.succeeded

            chems['Mg+2'][j] = numden*state.speciesAmount('Mg+2')[0]
            chems['H+'][j] = numden*state.speciesAmount('H+')[0]
            chems['OH-'][j] = numden*state.speciesAmount('OH-')[0]
            chems['CO3-2'][j] = numden*state.speciesAmount('CO3-2')[0]
            chems['HCO3-'][j] = numden*state.speciesAmount('HCO3-')[0]
            chems['SiO2(aq)'][j] = numden*state.speciesAmount('SiO2(aq)')[0]

            chems['CO2(aq)'][j] = numden*state.speciesAmount('CO2(aq)')[0]
            chems['CO2(g)'][j] = numden*state.speciesAmount('CO2(g)')[0]
            chems['PCO2'][j] = PCO2

            chems['Magnesite'][j] = numden*state.speciesAmount('Magnesite')[0]
            chems['Clino-Enstatite'][j] = numden*state.speciesAmount('Clino-Enstatite')[0]

            chems['pH'][j] = -np.log10(state.speciesAmount('H+')[0])

            j = j + 1

        df = pd.DataFrame({
            'Mg++': chems['Mg+2'],
            'Magnesite': chems['Magnesite'],
            'Silicates': 2*chems['Clino-Enstatite'],
        }, index=PCO2s)
        
        if table_flag == True:
            if nSiO2 == 0:
                df.to_csv('figA4b.csv')
            elif nSiO2 == 1:
                df.to_csv('fig4b.csv')
            #df.to_csv('phases_PCO2_'+DIV+'_beta'+str(int(beta*100))+'_nSiO2'+str(int(nSiO2))+'.csv')
            
        if plot_flag == True:
            plot_phases_PCO2(df, DIV = DIV, Temp = Temp, totP = totP,  beta = beta, 
                      nDIV = nDIV, nSiO2 = nSiO2, totnum = totnum)
    
    ####################################
    ####################################
    
    elif DIV == 'Fe':
        
        db = SupcrtDatabase('supcrtbl')

        solution = AqueousPhase(['H2O(aq)','CO2(aq)', 'HCO3-', 'CO3-2', 'H+', 'OH-', 'Fe+2', 'SiO2(aq)'])
        solution.setActivityModel(chain(
            ActivityModelHKF(),
            ActivityModelDrummond('CO2(aq)'),
        ))

        gases = GaseousPhase(['CO2(g)', 'N2(g)'])
        gases.setActivityModel(ActivityModelPengRobinson())

        minerals = MineralPhases(['Siderite', 'Fayalite', 'Quartz'])
        
        system = ChemicalSystem(db, solution, gases, minerals)

        specs = EquilibriumSpecs(system)
        specs.temperature()
        specs.pressure()
        specs.fugacity('CO2')

        solver = EquilibriumSolver(specs)

        j = 0
        while j < totnum:

            PCO2 = PCO2s[j]
            
            addDIVtot = nDIV * weath_scaling(PCO2,Temp,beta=beta) / numden
            addSiO2 = nSiO2 * addDIVtot

            state = ChemicalState(system)
            state.setTemperature(Temp, 'K')
            state.setPressure(totP, 'bar')
            state.set('H2O(aq)', totH2O - addDIVtot, 'mol')     # add ~ one kg of water
            state.set('N2(g)', totN2, 'mol')
            state.set('HCO3-', 2*addDIVtot, 'mol')
            state.set('Fe+2', addDIVtot, 'mol')
            state.set('SiO2(aq)', addSiO2, 'mol')

            conditions = EquilibriumConditions(specs)
            conditions.temperature(state.temperature())
            conditions.pressure(state.pressure())
            conditions.fugacity('CO2', PCO2, 'bar')

            result = solver.solve(state, conditions)

#            assert result.optima.succeeded

            chems['Fe+2'][j] = numden*state.speciesAmount('Fe+2')[0]
            chems['H+'][j] = numden*state.speciesAmount('H+')[0]
            chems['OH-'][j] = numden*state.speciesAmount('OH-')[0]
            chems['CO3-2'][j] = numden*state.speciesAmount('CO3-2')[0]
            chems['HCO3-'][j] = numden*state.speciesAmount('HCO3-')[0]
            chems['SiO2(aq)'][j] = numden*state.speciesAmount('SiO2(aq)')[0]

            chems['CO2(aq)'][j] = numden*state.speciesAmount('CO2(aq)')[0]
            chems['CO2(g)'][j] = numden*state.speciesAmount('CO2(g)')[0]
            chems['PCO2'][j] = PCO2

            chems['Siderite'][j] = numden*state.speciesAmount('Siderite')[0]
            chems['Fayalite'][j] = numden*state.speciesAmount('Fayalite')[0]

            chems['pH'][j] = -np.log10(state.speciesAmount('H+')[0])

            j = j + 1

        df = pd.DataFrame({
            'Fe++': chems['Fe+2'],
            'Siderite': chems['Siderite'],
            'Silicates': 2*chems['Fayalite'],
        }, index=PCO2s)
        
        if table_flag == True:
            if nSiO2 == 0:
                df.to_csv('figA4c.csv')
            elif nSiO2 == 1:
                df.to_csv('fig4c.csv')
            #df.to_csv('phases_PCO2_'+DIV+'_beta'+str(int(beta*100))+'_nSiO2'+str(int(nSiO2))+'.csv')
            
        if plot_flag == True:
            plot_phases_PCO2(df, DIV = DIV, Temp = Temp, totP = totP,  beta = beta, 
                      nDIV = nDIV, nSiO2 = nSiO2, totnum = totnum)
        
    return


# In[ ]:


# Plot stable phases as a function of PCO2
def plot_phases_PCO2 (df, DIV = 'Ca', Temp = 298, totP = 1, beta = 0.3, 
                      nDIV = 1, nSiO2 = 1, totnum = 100):
    '''
    Returns plots of stable phases as a function of PCO2 [bar]
    '''  
    if DIV == 'Ca':
        
        p11 = df.plot.area(stacked=True, color=colors[::3])
        
        p11.text(-0.1, 1.01, '(a)', transform=p11.transAxes, fontsize=16, va='top', ha='right')
        
        p11.set_xscale('log')

        p11.set_xlabel(r'$P_{\rm CO_2}$ [bar]',fontsize=16)
        p11.set_ylabel(r'$n$ [m$^{-3}$]', fontsize=18)
        p11.tick_params(axis='y')           

        if beta == 0:
            p11.text(1e-8, 1.4, r'$n_{\rm Ca,tot}$ = %.1f m$^{-3}$'%nDIV, c='black',
             ha='left', va='center'
            )
            p11.set_ylim([0,2])
        
        elif beta == 0.3:
            p11.text(1e-8, 5, r'$n_{\rm Ca,tot}$ = $f_{\rm W}(P_{\rm CO_2})$', c='black',
             ha='left', va='center'
            )
            p11.set_yscale('log')
            p11.set_ylim([0.1,50])

        p11.legend(fontsize=10, bbox_to_anchor=(0,1), borderaxespad=0, loc='upper left')

        p11.set_title(r'Ca Partitioning', fontsize=16)
        
        if nSiO2 == 0:
            plt.savefig('figA4a.pdf',bbox_inches='tight')
        elif nSiO2 == 1:
            plt.savefig('fig4a.pdf',bbox_inches='tight')
        #plt.savefig('phases_PCO2_Ca_beta%s_nSiO2%s.pdf'%(int(beta*100), int(nSiO2)),bbox_inches='tight')
        
    elif DIV == 'Mg':
        
        p11 = df.plot.area(stacked=True, color=colors[::3])
        
        p11.text(-0.1, 1.01, '(b)', transform=p11.transAxes, fontsize=16, va='top', ha='right')
        
        p11.set_xscale('log')

        p11.set_xlabel(r'$P_{\rm CO_2}$ [bar]',fontsize=16)
        p11.set_ylabel(r'$n$ [m$^{-3}$]', fontsize=18)
        p11.tick_params(axis='y')           

        if beta == 0:
            p11.text(1e-8, 1.4, r'$n_{\rm Mg,tot}$ = %.1f m$^{-3}$'%nDIV, c='black',
             ha='left', va='center'
            )
            p11.set_ylim([0,2])
        
        elif beta == 0.3:
            p11.text(1e-8, 5, r'$n_{\rm Mg,tot}$ = $f_{\rm W}(P_{\rm CO_2})$', c='black',
             ha='left', va='center'
            )
            p11.set_yscale('log')
            p11.set_ylim([0.1,50])

        p11.legend(fontsize=10, bbox_to_anchor=(0,1), borderaxespad=0, loc='upper left')

        p11.set_title(r'Mg Partitioning', fontsize=16)
        
        if nSiO2 == 0:
            plt.savefig('figA4b.pdf',bbox_inches='tight')
        elif nSiO2 == 1:
            plt.savefig('fig4b.pdf',bbox_inches='tight')
        # plt.savefig('phases_PCO2_Mg_beta%s_nSiO2%s.pdf'%(int(beta*100), int(nSiO2)),bbox_inches='tight')    
        
    elif DIV == 'Fe':

        p11 = df.plot.area(stacked=True, color=colors[::3])
        
        p11.text(-0.1, 1.01, '(c)', transform=p11.transAxes, fontsize=16, va='top', ha='right')
        
        p11.set_xscale('log')

        p11.set_xlabel(r'$P_{\rm CO_2}$ [bar]',fontsize=16)
        p11.set_ylabel(r'$n$ [m$^{-3}$]', fontsize=18)
        p11.tick_params(axis='y')           

        if beta == 0:
            p11.text(1e-8, 1.4, r'$n_{\rm Fe,tot}$ = %.1f m$^{-3}$'%nDIV, c='black',
             ha='left', va='center'
            )
            p11.set_ylim([0,2])
        
        elif beta == 0.3:
            p11.text(1e-8, 5, r'$n_{\rm Fe,tot}$ = $f_{\rm W}(P_{\rm CO_2})$', c='black',
             ha='left', va='center'
            )
            p11.set_yscale('log')
            p11.set_ylim([0.1,50])
            
        p11.legend(fontsize=10, bbox_to_anchor=(0,1), borderaxespad=0, loc='upper left')

        p11.set_title(r'Fe Partitioning', fontsize=16)
        
        if nSiO2 == 0:
            plt.savefig('figA4c.pdf',bbox_inches='tight')
        elif nSiO2 == 1:
            plt.savefig('fig4c.pdf',bbox_inches='tight')
        #plt.savefig('phases_PCO2_Fe_beta%s_nSiO2%s.pdf'%(int(beta*100), int(nSiO2)),bbox_inches='tight')
        
    return

