#!/usr/bin/env python
# coding: utf-8

# # OCRA: Ocean Chemistry with Reacktoro And beyond
# ### This Python code implements Reaktoro software to calculate ocean chemistry
# 
# ## Reference: Hakim et al. (2023) ApJL
# 
# ### solve.py # contains functions to setup and solve chemical systems using Reaktoro

# Import libraries

from reaktoro import *
from store import *

# Setup Reaktoro to solve ocean chemistry

def setup_an_Ca(Temp, totP):
    '''
    Returns chemical analytical setup with system, specs, solver for Ca
    '''
    db = SupcrtDatabase('supcrtbl')

    rxn9 = db.reaction('Ca+2 + CO3-2 = Calcite')
    logK9 = rxn9.props(Temp, 'K', totP, 'bar').lgK[0]

    rxn15 = db.reaction('CO2(g) + H2O(aq) = 2*H+ + CO3-2')
    logK16 = rxn15.props(Temp, 'K', totP, 'bar').lgK[0]

    rxn3= db.reaction('CO2(g) + H2O(aq) = H+ + HCO3-')
    logK3 = rxn3.props(Temp, 'K', totP, 'bar').lgK[0]
    
    return logK3, logK9, logK16

def setup_Ca():
    '''
    Returns chemical setup with system, specs, solver for Ca
    '''
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
    
    return system, specs, solver

def setup_Mg():
    '''
    Returns chemical setup with system, specs, solver for Mg
    '''
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
    
    return system, specs, solver

def setup_Fe():
    '''
    Returns chemical setup with system, specs, solver for Fe
    '''
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
    
    return system, specs, solver


# Solve Reaktoro ocean chemistry system

def solve_Ca(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP):
    '''
    Returns state for Ca
    '''
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

    return state

def solve_Mg(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP):
    '''
    Returns state for Mg
    '''
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

    return state

def solve_Fe(system, specs, solver, addDIVtot, addSiO2, PCO2, Temp, totP):
    '''
    Returns state for Fe
    '''
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

    return state
