#!/usr/bin/env python
# coding: utf-8

# # OCRA: Ocean Chemistry with Reacktoro And beyond
# ### This Python code implements Reaktoro software to calculate ocean chemistry
# 
# ## Reference: Hakim et al. (2023) ApJL
# 
# ### plots_paper.py # contains functions to make plots available in the published paper

# In[1]:


from ocra import *


# In[2]:


# Plot ocean pH as a function of PCO2 for the Ca-system
pH_PCO2(DIV='Ca')


# In[3]:


# Plot ocean pH as a function of PCO2 for the Mg-system
pH_PCO2(DIV='Mg')


# In[4]:


# Plot ocean pH as a function of PCO2 for the Fe-system
pH_PCO2(DIV='Fe')


# In[2]:


# Plot numerical and analytical solutions of ocean pH
pH_PCO2_an(DIV='Ca', nDIV_fixed = 1)


# In[6]:


# Plot ocean pH as a function of P
pH_P(DIV='Ca')


# In[7]:


# Plot ocean pH as a function of T
pH_T(DIV='Ca')


# In[11]:


# Plot Ca-CCD as a function of PCO2 and T with no silicates
CaCCD_PCO2_T(nSiO2 = 0, totnum = 20, numQ1 = 100, numQ2 = 100)


# In[12]:


# Plot Ca-CCD as a function of PCO2 and T with silicates
CaCCD_PCO2_T(totnum = 20, numQ1 = 100, numQ2 = 100)


# In[13]:


# Plot Mg-CCD as a function of PCO2 and T with no silicates
MgCCD_PCO2_T(beta = 0.3, totnum = 20, numQ1 = 100, numQ2 = 100)


# In[14]:


# Plot Mg-CCD as a function of PCO2 and T with silicates
MgCCD_PCO2_T(totnum = 20, numQ1 = 100, numQ2 = 100)


# In[9]:


# Plot Fe-CCD as a function of PCO2 and T with no silicates
FeCCD_PCO2_T(nSiO2 = 0, totnum = 20, numQ1 = 100, numQ2 = 100)


# In[10]:


# Plot Fe-CCD as a function of PCO2 and T with silicates
FeCCD_PCO2_T(totnum = 20, numQ1 = 100, numQ2 = 100)


# In[3]:


# Plot stable phases as a function of PCO2 for the Ca-system with no silicates
phases_PCO2(DIV='Ca', nSiO2 = 0, Temp = 310)


# In[4]:


# Plot stable phases as a function of PCO2 for the Ca-system with silicates
phases_PCO2(DIV='Ca', nSiO2 = 1, Temp = 310)


# In[5]:


# Plot stable phases as a function of PCO2 for the Mg-system with no silicates
phases_PCO2(DIV='Mg', nSiO2 = 0, Temp = 310)


# In[6]:


# Plot stable phases as a function of PCO2 for the Mg-system with silicates
phases_PCO2(DIV='Mg', nSiO2 = 1, Temp = 310)


# In[7]:


# Plot stable phases as a function of PCO2 for the Fe-system with no silicates
phases_PCO2(DIV='Fe', nSiO2 = 0, Temp = 310)


# In[8]:


# Plot stable phases as a function of PCO2 for the Fe-system with silicates
phases_PCO2(DIV='Fe', nSiO2 = 1, Temp = 310)


# In[ ]:




