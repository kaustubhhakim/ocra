#!/usr/bin/env python
# coding: utf-8

# # OCRA: Ocean Chemistry with Reacktoro And beyond
# ### This Python code implements Reaktoro software to calculate ocean chemistry
# 
# ## Reference: Hakim et al. (2023) ApJL
# 
# ### output.py # contains functions to output plots and csv tables

# In[ ]:


import numpy as np
import pandas as pd

from store import *

import matplotlib.pyplot as plt
# plt.rcParams['pcolor.shading'] ='nearest'
import matplotlib.gridspec as gridspec
from matplotlib import ticker

from matplotlib.cm import get_cmap

col1 = get_cmap('Dark2').colors  # type: matplotlib.colors.ListedColormap
col2 = get_cmap('Set1').colors
col3 = get_cmap('Set3').colors
colors = col1 + col2 + col3


# In[ ]:


# Plot analytical and numerical limits of pH as a function of PCO2

def output_pH_PCO2(PCO2s, chems2, DIV = 'Ca', plot_flag = True, table_flag = True):
    '''
    Returns tables and plots of ocean pH as a function of PCO2 [bar] for Ca, Mg or Fe carbonate systems
    '''
    if DIV == 'Ca':
        
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


# Plot numerical and analytical limits of pH as a function of PCO2

def output_pH_PCO2_an(PCO2s, chems2, chems2_an, chems2_san, DIV = 'Ca', nDIV_fixed = 1,
                      plot_flag = True, table_flag = True):
    '''
    Returns tables and plots of numerical and analytical limits of ocean pH as a function of PCO2 [bar] 
    '''
    if DIV == 'Ca':
        
        ####################################
        
        if table_flag == True:
            
            df = pd.DataFrame({
                'pH lower numerical': chems2['pH'][0],
                'pH lower analytical': chems2_an['pH'][0],
                'pH upper numerical': chems2['pH'][1],
                'pH upper analytical': chems2_an['pH'][1],
                'pH upper semi-analy.': chems2_san['pH'][1],
            }, index=PCO2s)
            df.to_csv('figA1.csv')
            #df.to_csv('pHan_PCO2_Ca.csv')
            
        ####################################
        
        if plot_flag == True:

            fig1 = plt.figure(constrained_layout=False,figsize=(5,5))
            spec1 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig1)

            p11  = fig1.add_subplot(spec1[0,0])

            p11.fill_between(PCO2s, chems2['pH'][1], chems2['pH'][0], facecolor='blue', alpha=0.3)
            p11.fill_between(PCO2s, chems2['pH'][0], 4, facecolor='red', alpha=0.1)
            p11.fill_between(PCO2s, 11, chems2['pH'][1], facecolor='red', alpha=0.1)
            p11.plot(PCO2s, chems2['pH'][1], lw=5,c='black',ls='-', label = 'Up (numerical)')
            p11.plot(PCO2s, chems2_san['pH'][1], lw=3,c='orange',ls=':', label = 'Up (semi-analytical)')
            p11.plot(PCO2s, chems2_an['pH'][1], lw=3,c='orange',ls='--', 
                     label = r'Up (ana., $n_{\rm Ca^{2+}}$ = %d m$^{-3}$)'%nDIV_fixed)

            p11.plot(PCO2s, chems2['pH'][0], lw=5,c='blue',ls='-', label = 'Low (numerical)')
            p11.plot(PCO2s, chems2_an['pH'][0], lw=3,c='cyan',ls='--', label = 'Low (analytical)')

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


# Plot pH as a function of local pressure

def output_pH_P(totPs, chems2, DIV = 'Ca', plot_flag = True, table_flag = True):
    '''
    Returns tables and plots of ocean pH as a function of P [bar]
    '''
    if DIV == 'Ca':

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


# Plot pH as a function of temperature

def output_pH_T(Temps, chems2, DIV = 'Ca', plot_flag = True, table_flag = True):
    '''
    Returns tables and plots of ocean pH as a function of Temp [K]
    '''
    if DIV == 'Ca':

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


# Output Ca-CCD as a function of PCO2 and T

def output_CaCCD_PCO2_T(PCO2s, Temps, CCDs, beta = 0.3, nSiO2 = 1, nDIV = 1, 
                        table_flag = True, plot_flag = True):
    '''
    Returns tables and plots of Ca-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    if table_flag == True:
        df = pd.DataFrame(CCDs, index=PCO2s, columns=Temps)
        if nSiO2 == 0:
            df.to_csv('figA3a.csv')
        elif nSiO2 == 1:
            df.to_csv('fig3a.csv')
        #df.to_csv('CaCCD_PCO2_T_'+str(int(beta*100))+'_nSiO2'+str(int(nSiO2))+'.csv')
        
    if plot_flag == True:
        plot_CaCCD_PCO2_T(PCO2s, Temps, CCDs, beta = beta, nSiO2 = nSiO2, nDIV = nDIV)

    return

def plot_CaCCD_PCO2_T(PCO2s, Temps, CCDs, beta = 0.3, nSiO2 = 1, nDIV = 1):
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


# Plot Mg-CCD as a function of PCO2 and T

def output_MgCCD_PCO2_T(PCO2s, Temps, CCDs, beta = 0.3, nSiO2 = 1, nDIV = 1, 
                        table_flag = True, plot_flag = True):
    '''
    Returns tables and plots of Mg-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    if table_flag == True:
        df = pd.DataFrame(CCDs, index=PCO2s, columns=Temps)
        if nSiO2 == 0:
            df.to_csv('figA3b.csv')
        elif nSiO2 == 1:
            df.to_csv('fig3b.csv')
        #df.to_csv('MgCCD_PCO2_T_'+str(int(beta*100))+'_nSiO2'+str(int(nSiO2))+'.csv')
        
    if plot_flag == True:
        plot_MgCCD_PCO2_T(PCO2s, Temps, CCDs, beta = beta, nSiO2 = nSiO2, nDIV = nDIV)

    return

def plot_MgCCD_PCO2_T(PCO2s, Temps, CCDs, beta = 0.3, nSiO2 = 1, nDIV = 1):
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





# In[ ]:


# Plot Fe-CCD as a function of PCO2 and T

def output_FeCCD_PCO2_T(PCO2s, Temps, CCDs, beta = 0.3, nSiO2 = 1, nDIV = 1, 
                        table_flag = True, plot_flag = True):
    '''
    Returns tables and plots of Fe-CCD [km] as a function of PCO2 [bar] and Temp [K]
    '''
    if table_flag == True:
        df = pd.DataFrame(CCDs, index=PCO2s, columns=Temps)
        if nSiO2 == 0:
            df.to_csv('figA3c.csv')
        elif nSiO2 == 1:
            df.to_csv('fig3c.csv')
        #df.to_csv('FeCCD_PCO2_T_'+str(int(beta*100))+'_nSiO2'+str(int(nSiO2))+'.csv')
        
    if plot_flag == True:
        plot_FeCCD_PCO2_T(PCO2s, Temps, CCDs, beta = beta, nSiO2 = nSiO2, nDIV = nDIV)

    return

def plot_FeCCD_PCO2_T(PCO2s, Temps, CCDs, beta = 0.3, nSiO2 = 1, nDIV = 1):
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


def output_phases_PCO2(df, DIV = 'Ca', beta = 0.3, nDIV = 1, nSiO2 = 1, table_flag = True, plot_flag = True):
    '''
    Returns tables and plots of stable phases as a function of PCO2 [bar]
    '''  
    if DIV == 'Ca':
        
        if table_flag == True:
            if nSiO2 == 0:
                df.to_csv('figA4a.csv')
            elif nSiO2 == 1:
                df.to_csv('fig4a.csv')
            #df.to_csv('phases_PCO2_'+DIV+'_beta'+str(int(beta*100))+'_nSiO2'+str(int(nSiO2))+'.csv')
            
        if plot_flag == True:
            plot_phases_PCO2(df, DIV = DIV, beta = beta, nDIV = nDIV, nSiO2 = nSiO2)
            
    elif DIV == 'Mg':
        
        if table_flag == True:
            if nSiO2 == 0:
                df.to_csv('figA4b.csv')
            elif nSiO2 == 1:
                df.to_csv('fig4b.csv')
            #df.to_csv('phases_PCO2_'+DIV+'_beta'+str(int(beta*100))+'_nSiO2'+str(int(nSiO2))+'.csv')
            
        if plot_flag == True:
            plot_phases_PCO2(df, DIV = DIV, beta = beta, nDIV = nDIV, nSiO2 = nSiO2)
            
    elif DIV == 'Fe':

        if table_flag == True:
            if nSiO2 == 0:
                df.to_csv('figA4c.csv')
            elif nSiO2 == 1:
                df.to_csv('fig4c.csv')
            #df.to_csv('phases_PCO2_'+DIV+'_beta'+str(int(beta*100))+'_nSiO2'+str(int(nSiO2))+'.csv')
            
        if plot_flag == True:
            plot_phases_PCO2(df, DIV = DIV, beta = beta, nDIV = nDIV, nSiO2 = nSiO2)
            
    return


# In[ ]:


# Plot stable phases as a function of PCO2

def plot_phases_PCO2 (df, DIV = 'Ca', beta = 0.3, nDIV = 1, nSiO2 = 1):
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

