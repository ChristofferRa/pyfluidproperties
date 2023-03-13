# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 21:30:37 2023

@author: Christoffer
"""

from pyfluidproperties import properties as prop
from pyfluidproperties import unit_conv
import numpy as np

def main():
    
    ### Inputs ###
    pressure =      10.3*1e5        # Pa
    temperature =   273.15 + 150    # K
    
    pipe_diameter = 80*1e-3         # m
    pipe_length =   100             # m
    roughness =     .05*1e-3        # m
    
    mass_flow =     20              # kg/s 
    
    ### Init objects ###
    fluid = prop(fluid = "H2O", unit_system = 'SI')
    dim_num = dimensionless_numbers(fluid)
    
    ### Calculations ###
    fluid.update(p = pressure, T = temperature)
    
    pipe_area = pipe_diameter**2*np.pi/4
    fluid_velocity = mass_flow/fluid.rho/pipe_area
    
    re_d = dim_num.reynolds_number(fluid_velocity, pipe_diameter)
    
    ### Output ###
    print(fluid)
    print(f'Fluid velocity = {fluid_velocity:.1f} m/s, Reynolds number = {(re_d):.0f} -')

class dimensionless_numbers:
    """
    Class for calculating dimensionless numbers
    """
    
    
    def __init__(self, fluid_object):
        """
        Initilize dimensionless object 

        Parameters
        ----------
        fluid_object : pyfluidproperties     fluid object containing fluid properties.

        Returns
        -------
        None.

        """
        # Set upp unit-conversions 
        #self.unit_system = unit_conv
        
        # Save fluid properties
        self.fluid_properties = fluid_object
        
    def reynolds_number(self, fluid_velocity ,characteristic_length):
        """
        Reynolds number
        https://en.wikipedia.org/wiki/Reynolds_number

        Parameters
        ----------
        fluid_velocity :        float   fluid velocity
        characteristic_length : float   characteristic length.

        Returns
        -------
        float   Reynolds number.

        """
        
        return self.fluid_properties.rho * fluid_velocity * characteristic_length / self.fluid_properties.my

if __name__ == "__main__":
    main()
