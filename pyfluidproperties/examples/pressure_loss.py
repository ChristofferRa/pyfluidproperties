# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 20:55:36 2023

@author: Christoffer
"""

from pyfluidproperties import properties as prop
from pyfluidproperties import unit_conv
from dimensionless_numbers import dimensionless_numbers

import numpy as np

def main():
    
    ### Inputs ###
    pressure =      10.3            # bar(a)
    temperature =   250             # C
    
    pipe_diameter = 80*1e-3         # m
    pipe_length =   100             # m
    roughness =     .05*1e-3        # m
    
    mass_flow =     .5              # kg/s 
    
    ### Init objects ###
    fluid = prop(fluid = "H2O", unit_system = 'SI_bar_kJ')
    dim_num = dimensionless_numbers(fluid)
    
    ### Calculations ###
    fluid.update(p = pressure, T = temperature)
    
    pipe_area = pipe_diameter**2*np.pi/4
    fluid_velocity = mass_flow/fluid.rho/pipe_area
    
    re_d = dim_num.reynolds_number(fluid_velocity, pipe_diameter)
    fric = friction_factor(pipe_diameter, roughness, re_d)
    
    pressure_loss = fric*pipe_length/pipe_diameter*.5*fluid.rho*fluid_velocity**2
    
    ### Output ###
    print(fluid)
    print(f'Fluid velocity = {fluid_velocity:.1f} m/s, pressure loss = {(pressure_loss*1e-5):.2f} bar')

    
def friction_factor(pipe_diameter,roughness, reynolds_number):
    """
    Returns the darcy friction factor using 
    Swamee-Jain approximation (turbulent flow)

    Parameters
    ----------
    pipe_diameter:          float   pipe diameter (m)
    roughness :             float   roughness eps (m)
    reynolds_number :       float   Reynolds number (-)

    Returns
    -------
    float   darcy friction factor (-)

    """
    
    reynolds_laminar = 2000  # (-)
    reynold_turbulent = 4000 # (-)
    
    if reynolds_number < reynolds_laminar:
        friction_factor = 64/ reynolds_number
    else:
        # Swamee-Jain
        relative_roughness = roughness/pipe_diameter
        friction_factor = 0.25/((np.log10((relative_roughness/3.7)+(5.74/(reynolds_number**.9))))**2)
        if reynolds_number <= reynold_turbulent:
            print('Warning, Critical regime. Flow is unsteady, Laminar friction factor assumed -> higher pressure loss')
    return friction_factor

if __name__ == "__main__":
    main()
