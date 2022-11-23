#!/usr/bin/env python
# coding: utf-8

# # Plots available in the paper using OCRA
# ## Hakim et al. (2023) ApJL
# 
# This Python script implements Reaktoro software to calculate ocean chemistry

# In[ ]:


from ocra import *


# In[ ]:


# Plot ocean pH as a function of PCO2 for the Ca-system
plot_pH(DIV='Ca')


# In[ ]:


# Plot ocean pH as a function of PCO2 for the Mg-system
plot_pH(DIV='Mg')


# In[ ]:


# Plot ocean pH as a function of PCO2 for the Fe-system
plot_pH(DIV='Fe')


# In[ ]:


# Plot numerical and analytical solutions of ocean pH
plot_pH_an(DIV='Ca', nDIV_fixed = 1)


# In[ ]:


# Plot ocean pH as a function of P
plot_pH_P(DIV='Ca')


# In[ ]:


# Plot ocean pH as a function of T
plot_pH_T(DIV='Ca')


# In[ ]:


# Plot Ca-CCD as a function of PCO2 and T with no silicates
CaCCD_PCO2_T(nSiO2 = 0, totnum = 20, numQ1 = 100, numQ2 = 100)


# In[ ]:


# Plot Ca-CCD as a function of PCO2 and T with silicates
CaCCD_PCO2_T(totnum = 20, numQ1 = 100, numQ2 = 100)


# In[ ]:


# Plot Mg-CCD as a function of PCO2 and T with no silicates
MgCCD_PCO2_T(nSiO2 = 0, totnum = 20, numQ1 = 100, numQ2 = 100)


# In[ ]:


# Plot Mg-CCD as a function of PCO2 and T with silicates
MgCCD_PCO2_T(totnum = 20, numQ1 = 100, numQ2 = 100)


# In[ ]:


# Plot Fe-CCD as a function of PCO2 and T with no silicates
FeCCD_PCO2_T(nSiO2 = 0, totnum = 20, numQ1 = 100, numQ2 = 100)


# In[ ]:


# Plot Fe-CCD as a function of PCO2 and T with silicates
FeCCD_PCO2_T(totnum = 20, numQ1 = 100, numQ2 = 100)


# In[ ]:


# Plot stable phases as a function of PCO2 for the Ca-system with no silicates
phases_PCO2(DIV='Ca', nSiO2 = 0, Temp = 310)


# In[ ]:


# Plot stable phases as a function of PCO2 for the Ca-system with silicates
phases_PCO2(DIV='Ca', nSiO2 = 1, Temp = 310)


# In[ ]:


# Plot stable phases as a function of PCO2 for the Mg-system with no silicates
phases_PCO2(DIV='Mg', nSiO2 = 0, Temp = 310)


# In[ ]:


# Plot stable phases as a function of PCO2 for the Mg-system with silicates
phases_PCO2(DIV='Mg', nSiO2 = 1, Temp = 310)


# In[ ]:


# Plot stable phases as a function of PCO2 for the Fe-system with no silicates
phases_PCO2(DIV='Fe', nSiO2 = 0, Temp = 310)


# In[ ]:


# Plot stable phases as a function of PCO2 for the Fe-system with silicates
phases_PCO2(DIV='Fe', nSiO2 = 1, Temp = 310)

