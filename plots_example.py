#!/usr/bin/env python
# coding: utf-8

# # OCRA: Ocean Chemistry with Reacktoro And beyond
# ### This Python code implements Reaktoro software to calculate ocean chemistry
# 
# ## Reference: Hakim et al. (2023) ApJL
# 
# ### plots_example.py # contains functions to make Example plots to test OCRA

# In[1]:


from ocra import *


# In[2]:


# Plot ocean pH as a function of PCO2
pH_PCO2(DIV='Ca')


# In[4]:


# Plot Ca-CCD as a function of PCO2 and T
CaCCD_PCO2_T(nSiO2 = 1, totnum = 10, numQ1 = 20, numQ2 = 20)


# In[3]:


# Plot stable phases as a function of PCO2
phases_PCO2(DIV='Ca', nSiO2 = 1, Temp = 310)


# In[ ]:




