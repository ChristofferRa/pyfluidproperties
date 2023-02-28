# -*- coding: utf-8 -*-
"""
Subclass to fluid_prop, implements IAPWSIF-97.

"Properties of Water and Steam
IAPWS Industrial Formulation 1997"

This class implements functions from the different modules sorted by IAPWS-releases.
The class functions checks that the input is within the valid region (find_region* functions)
and then uses the appropiate functions for that region. The modules only contains functions
for the property-equations presented in these releases. These are either used directly or 
in calculations to determine the desired property.

The class may be used either as an object, that contain all the basic properties given
pressure-temperature, pressure-quality, pressure-enthalpy or pressure-entropy using the
correspondig update functions. Or it can be used directly for each property by using 
its built in static methods. 

Basic properties for iapwsif97 class:
    Description                         Unit        Letter
    Pressure                            (Pa)        p
    Temperature                         (K)         T,t
    Specific Volume                     (m^3/kg)    v
    Specific Enthalpy                   (J/kg)      h
    Specific Internal Energy            (J/kg)      u
    Specific Entropy                    (J/kg/K)    s
    Specific Isobaric heat capacity     (J/kg/K)    cp
    Specific Isochoric heat capacity    (J/kg/K)    cv
    Speed of sound                      (m/s)       w 
    Vapor mass fraction                 (-)         x
    Dynamic Viscosity                   (Pa*s)      my
    Thermal Conductivity                (W/m/K)     tc
    Surface Tension                     (N/m)       st

List of functions and code structure for iapwsif97 class:
    Global settings set-functions
       set_rho_r3_itt(itt)
       
   Updating functions
       update_pt(p, t)
       update_px(p, x)
       update_tx(t, x)
       update_ph(p, h)
       update_ps(p, s)
       update_hs(h, s)
   
   _pt - properties
       v_pt(p, t)
       h_pt(p, t)
       u_pt(p, t)
       s_pt(p, t)
       cp_ptv
       cv_pt(p, t)
       w_pt(p, t)
       rho_pt(p, t)
       my_pt(p, t)
       tc_pt(p, t)
   
   _ph - properties
       v_ph(p, h)
       u_ph(p, h)
       s_ph(p, h)
       cp_ph(p, h)
       cv_ph(p, h)
       w_ph(p, h)
       rho_ph(p, h)
       my_ph(p, h)
       tc_ph(p, h)
   
   _ps - properties
       v_ps(p, s)
       h_ps(p, s)
       u_ps(p, s)
       cp_ps(p, s)
       cv_ps(p, s)
       w_ps(p, s)
       rho_ps(p, s)
       my_ps(p, s)
       tc_ps(p, s)
   
   _hs - properties
       v_hs(h, s)
       u_hs(h, s)
       cp_hs(h, s)
       cv_hs(h, s)
       w_hs(h, s)
       rho_hs(h, s)
       my_hs(h, s)
       tc_hs(h, s)
       
    _px - properties
       v_px(p, x)
       h_px(p, x)
       u_px(p, x)
       s_px(p, x)
       rho_px(p, x)
    
    _tx- properties
       v_tx(t, x)
       h_tx(t, x)
       u_tx(u, x)
       s_tx(t, x)
       rho_tx(t, x)
       
   Temperature, pressure, enthalpy backward functions
       T_ph(p, h)
       T_ps(p, s)
       T_hs(h, s)
       
       Tsat_p(p)
       Tsat_hs(h, s)
       
       psat_t(t)
       psat_s(s)
       psat_hs(h, s)
       
       p_hs(h, s)
       
       h_prho(p, rho)
       
   Vapor Mass fractions
       x_ph(p, h)
       x_ps(p, s)
       x_hs(h, s)
   
   Vapor Volume fractions
       vx_ph(p, h)
       vx_ps(p, s)
       vx_hs(h, s)
   
   Saturation properties, Liquid and Vapor
       x pressure
           hL_p(p),     hV_p(p)
           sL_p(p),     sV_p(p)
           rhoV_p(p),   rhoL_p(p)
           uV_p(p),     uL_p(p)
           CpV_p(p),    CpL_p(p)
           CvV_p(p),    CvL_p(p)
           wV_p(p),     wL_p(p)
           tcL_p(p),    tcV_p(p)
           myL_p(p),    myV_p(p)
       
       x temperature
           hL_T(t),     hV_T(t)
           sL_T(t),     sV_T(t)
           vL_T(t),     vV_T(t)
           uL_T(t),     uV_T(t)
           CpL_T(t),    CpV_T(t)
           CvL_T(t),    CvV_T(t)
           wL_T(t),     wV_T(t)
           rhoV_T(t),   rhoL_Ts(t)
           tcL_T(t),    tcV_T(t)
           myL_T(t),    myV_T(t)
           
    Misc
        st_t(t)
        st_p(p)
       
Known bugs:
    - 
Not yet implemented/todo:
    -
Future uppgrades:
    - Change itteration method to ITP or Illinois to speed upp?
    - Tsat_hs, lower accuracy for pressures close to p_crit. mainly in region 3 p>16MPa, double itterations, possible to do it a better way?
    - high accuracy-mode, itterate instead of using backwards equations. (gives higher accuracy) use backward equation for initial guess. implemented for _pt functions in region 3...
    - make function-call in class for meta-stable region 2...
    - make function for gamma/kappa and calculate two-phase values correct
    - add two-phase correlation for speed of sound
    - update_hs function
Bugs:
    - find_region_hs, uses fitted polynominals for region boundaries. som inacuracies may show h-s point as outside region when not.

@author: Christoffer Rappmann, christoffer.rappmann@gmail.com
"""
# Parent class
from fluid_prop import fluid_prop

# sub functions
# Global Constants
import iapwsif97_globals as global_property
# Main iapwsif97-implementation
import iapwsif97_main as if97
# Iapwsif97 suplementary release 5, IAPWS SR5-05(2016)
import iapwsif97_vpt3 as if97_vpt3
# Iapwsif97 suplementary release 3, IAPWS SR3-03(2014)
import iapwsif97_tps3_tph3_vph3_vps3 as if97_tvphps
# Iapwsif97 suplementary release 2 and 4, IAPWS SR2-01(2014) and SR4-04(2014)
import iapwsif97_phs as if97_phs
# iapwsf08, viscosity of ordinary water, R12-08
import iapwsf08_my as f08
# iapwsf11, thermal conductivity of ordinary water, R15-11
import iapwsf11_tc as f11
# iapwsR1-76, surface tension
import iapwsR176_st as r176
# Helper functions
import iapwsif97_helper_functions as aux

# Dependencies
import numpy as np

#### Global Settings ###
rho_r3_itt = False # Chose to use v_r3 backward equations or use itteration

### Class ###
class iapwsif97(fluid_prop):
    """
    Subclass to fluid_prop, implements IAPWSIF-97
    
    Properties of Water and Steam
    IAPWS Industrial Formulation 1997
    
    @author: Christoffer Rappmann, christoffer.rappmann@gmail.com
    """
    
    def __init__(self, p = 1*1e5, T = 20+273.15):
        """
        Initialize fluid properties
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        T:      double  temperature (K).

        Returns
        -------
        None.

        """
        super().__init__(p,T)
                
        iapwsif97.update_pt(self,p, T)    
      
    def __str__(self):
        """
        Lists basic properties of object as a string
        """
        return(f'\n\nFluid: Water (H2O)\n\nPressure,\t\tp = {self.p*10**-6:.6f}\tMPa\n' +
               f'Temperature,\tT = {self.T:.6f}\tK\n\n' +
               f'region =\t{self.region}\n\n'
               f'Specific Volume,\t\t\t\t\t\tv  =\t\t{self.v:.11f}\t\t(m^3/kg)\n' +
               f'Specific Enthalpy,\t\t\t\t\t\th  =\t\t{self.h*10**-3:.11f}\t(kJ/kg)\n' +
               f'Specific Internal Energy,\t\t\t\tu  =\t\t{self.u*10**-3:.11f}\t(kJ/kg)\n' +
               f'Specific Entropy,\t\t\t\t\t\ts  =\t\t{self.s*10**-3:.11f}\t\t(kJ/kg/K)\n' +
               f'Specific Isobaric heat capacity,\t\tcp =\t\t{self.cp*10**-3:.11f}\t(kJ/kg/K)\n' +
               f'Specific Isochoric heat capacity,\t\tcv =\t\t{self.cv*10**-3:.11f}\t\t(kJ/kg/K)\n' +
               f'Speed of sound,\t\t\t\t\t\t\tw  =\t\t{self.w:.11f}\t\t(m/s)\n' + 
               f'Vapor mass fraction,\t\t\t\t\tx  =\t\t{self.x:.11f}\t\t(-)\n' +
               f'Dynamic Viscosity,\t\t\t\t\t\tmy =\t\t{self.my:.11f}\t\t(Pa*s)\n' +
               f'Thermal Conductivity,\t\t\t\t\ttc =\t\t{self.tc:.11f}\t\t(W/m/K)\n' +
               f'Surface Tension,\t\t\t\t\t\tst =\t\t{self.st:.11f}\t\t(N/m)\n')
    
    
    def set_rho_r3_itt(itt):
        """
        Set method to calculate rho/v in region 3
        
        True = use itteration (higher accuracy)
        False = use v_pt function from SR5-05 suplementary release (default)
        
        This is a global setting, and will affect all objects?
        
        Parameters
        ----------
        itt:      boolean     -

        Returns
        -------
        None.
        """
        #Region 3 uses p-rho to calculate properties, no backward eqation for
        #rho available i original iapwsif97, however a suplementary release is 
        #available, SR5-05, that defines a backward equation for v_pt in region 3.
        #The accuracy of this equation is lower, therefore a option to use itteration
        #is provided using this setting. v_pt-function is used to provide a initial
        #guess for the itterations.
            
        global rho_r3_itt
        
        rho_r3_itt = itt
    
    ### Object updating functions ####
    
    def update_pt(self, p, T):
        """
        Update fluid properties
        IAPWS-IF97
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        T:      double  temperature (K).

        Returns
        -------
        None.
        """
        
        self.p = p                      # Pressure (Pa)
        self.T = T                      # Temperature (K)
        
        region = if97.find_region(p, T)
        self.region = region
        
        if region == -1:
            # If outside applicable region, set value to -1
            print(f'update_pt, out of bounds. p = {p*10**-6:.3f} MPa, T = {T:.3f} K')
            self.v = np.nan
            self.h = np.nan
            self.u = np.nan
            self.s = np.nan
            self.cp = np.nan
            self.cv = np.nan
            self.w = np.nan
            self.x = np.nan
        elif region == 1:
            self.v = if97.v_r1(p, T)    # Specific volume                       (m^3/kg)
            self.h = if97.h_r1(p, T)    # Specific enthalpy                     (J/kg)
            self.u = if97.u_r1(p, T)    # Specific internal energy              (J/kg)
            self.s = if97.s_r1(p, T)    # Specific entropy                      (J/kg/K)
            self.cp = if97.cp_r1(p, T)  # Specific isobaric heat capacity       (J/kg/K)
            self.cv = if97.cv_r1(p, T)  # Specific isochoric heat capacity      (J/kg/K)
            self.w = if97.w_r1(p, T)    # Speed of sound                        (m/s)
            self.x = 0.0                # Mass fraction of steam                (-)
        elif region == 2:
            self.v = if97.v_r2(p, T)
            self.h = if97.h_r2(p, T)
            self.u = if97.u_r2(p, T)
            self.s = if97.s_r2(p, T)
            self.cp = if97.cp_r2(p, T)
            self.cv = if97.cv_r2(p, T)
            self.w = if97.w_r2(p, T)
            self.x = 1.0
        elif region == 3:
            # Region defined by rho,T. rho needs to be determined
            # by backwards eq. for v_pt defined in Supp-VPT3-2016 
            # (IAPWS SR5-05(2016)) or by itteration
            
            rho = 1/if97_vpt3.v_pt(p, T)

            if rho_r3_itt: # if this option is set, itterate a more precise rho
                rho = if97.rho_r3(p, T, rho)
                
            self.v = 1/rho
            self.h = if97.h_r3(rho, T)
            self.u = if97.u_r3(rho, T)
            self.s = if97.s_r3(rho, T)
            self.cp = if97.cp_r3(rho, T)
            self.cv = if97.cv_r3(rho, T)
            self.w = if97.w_r3(rho, T)
            
            # Check if super critical, steam or liquid since region 3
            # cointains all states.
            x = np.nan
            h_crit = global_property.h_crit # J/kg, critical enthalpy
            
            if self.h > h_crit:
                x = 1.0
            elif p < if97_tvphps.p_3sat_h(self.h):
                x = 1.0
            else:
                x = 0.0
                
            self.x = x
        elif region == 4:
            # Saturation properties. 
            # Will not be entered this way
            print('Error, region 4')
        elif region == 5:
            self.v = if97.v_r5(p, T)
            self.h = if97.h_r5(p, T)
            self.u = if97.u_r5(p, T)
            self.s = if97.s_r5(p, T)
            self.cp = if97.cp_r5(p, T)
            self.cv = if97.cv_r5(p, T)
            self.w = if97.w_r5(p, T)
            self.x = 1.0
            
        if region == -1:
            self.rho = np.nan
            self.my = np.nan
            self.tc = np.nan
            self.st = np.nan
        else:
            if self.v != 0.0 and self.v != np.nan:
                self.rho = 1/self.v
            else:
                self.rho = np.nan
            self.my = iapwsif97.my_pt(p,T)
            self.tc = iapwsif97.tc_pt(p,T)
            self.st = np.nan # surface tension only between phases, needs to be two-phase
        
    def update_px(self, p, x):
        """
        Update fluid properties
        IAPWS-IF97
        
        Saturated properties, region 4
        2-phase
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        x:      double  mass-fraction steam (-).

        Returns
        -------
        None.
        """
        
        ### Constants ###
        p_crit =                global_property.p_crit                  # Pa, Critical pressure
        p_r3_lower_boundary =   global_property.p_r3_lower_boundary     # Pa, Region 1/2-3 boundary
        
        # Initialize properties
        v_l = np.nan    # Liquid Specific volume
        h_l = np.nan    # Liquid Specific enthalpy                     
        u_l = np.nan    # Liquid Specific internal energy              
        s_l = np.nan    # Liquid Specific entropy                      
        v_v = np.nan    # Vapor Specific volume                       
        h_v = np.nan    # Vapor Specific enthalpy                     
        u_v = np.nan    # Vapor Specific internal energy              
        s_v = np.nan    # Vapor Specific entropy        
        
        # mass fraction must be within 0.0<=x<=1.0
        if x >= 0.0 and x <= 1.0:
            T = if97.t_sat_p_r4(p)
            self.region = 4
        
            self.p = p
            self.T = T
            self.x = x
                      
            # Check if pressure is below region 3 boundary
            if p < p_r3_lower_boundary:
                # Use region 1 and 2 to calculate saturation properties
                
                # Saturated liquid properties
                v_l = if97.v_r1(p, T)    # Liquid Specific volume                       (m^3/kg)
                h_l = if97.h_r1(p, T)    # Liquid Specific enthalpy                     (J/kg)
                u_l = if97.u_r1(p, T)    # Liquid Specific internal energy              (J/kg)
                s_l = if97.s_r1(p, T)    # Liquid Specific entropy                      (J/kg/K)
                
                # Saturated steam properties
                v_v = if97.v_r2(p, T)    # Vapor Specific volume                       (m^3/kg)
                h_v = if97.h_r2(p, T)    # Vapor Specific enthalpy                     (J/kg)
                u_v = if97.u_r2(p, T)    # Vapor Specific internal energy              (J/kg)
                s_v = if97.s_r2(p, T)    # Vapor Specific entropy                      (J/kg/K)
        
            # pressure is above region 3 lower boundary, check if its below critical pressure
            elif p < p_crit:
                # Use region 3 for saturation properties    
                #h_l and h_v is determined from their functions and then used
                # to determine densities for the liquid and vapor respectevily, then both
                # rho and T is known and thus the rest of the properties can be calculated.
                
                h_l = iapwsif97.hL_p(p)
                h_v = iapwsif97.hV_p(p)
            
                # Saturated liquid properties
                v_l = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_l),h_l) # use p from p_3sat_h to ensure region 3 otherwise error
                u_l = if97.u_r3(1.0/v_l, T)
                s_l = if97.s_r3(1.0/v_l, T)
                # Saturated steam properties
                v_v = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_v),h_v) # # use p from p_3sat_h to ensure region 3 otherwise error
                u_v = if97.u_r3(1.0/v_v, T)
                s_v = if97.s_r3(1.0/v_v, T)
            
            else:
                print("Pressure above critical pressure")
                self.x = np.nan

            # Mixed 2-phase properties
            self.v = x*v_v + (1-x)*v_l
            self.h = x*h_v + (1-x)*h_l
            self.u = x*u_v + (1-x)*u_l
            self.s = x*s_v + (1-x)*s_l
            
            self.st = iapwsif97.st_t(T)
            
        else:
            print(f'x ={x:.2f} must be within 0.0 <= x <= 1.0')
            
            self.region = -1
        
            self.p = p
            self.T = np.nan
            self.x = np.nan
            self.v = np.nan
            self.h = np.nan
            self.u = np.nan
            self.s = np.nan
            self.st = np.nan
                       
        # N/A for 2-phase
        self.cp = np.nan
        self.cv = np.nan
        self.w  = np.nan # Not implemented, not part of iapwsif97
        self.my = np.nan # Not defined for two-phase
        self.tc = np.nan # Not defined for two-phase
        
        if self.v != 0.0 and self.v != np.nan:
            self.rho = x*(1/v_v) + (1-x)*(1/v_l)
        else:
            self.rho = np.nan
            
    def update_tx(self, T, x):
        """
        Update fluid properties
        IAPWS-IF97
        
        Saturated properties, region 4
        2-phase
    
        Parameters
        ----------
        T:      double  temperature (K).
        x:      double  mass-fraction steam (-).

        Returns
        -------
        None.
        """
        # Initialize properties
        v_l = np.nan    # Liquid Specific volume
        h_l = np.nan    # Liquid Specific enthalpy                     
        u_l = np.nan    # Liquid Specific internal energy              
        s_l = np.nan    # Liquid Specific entropy                      
        v_v = np.nan    # Vapor Specific volume                       
        h_v = np.nan    # Vapor Specific enthalpy                     
        u_v = np.nan    # Vapor Specific internal energy              
        s_v = np.nan    # Vapor Specific entropy        
        
        # mass fraction must be within 0.0<=x<=1.0
        if x >= 0.0 and x <= 1.0:
            # Saturation pressure
            p = iapwsif97.psat_t(T)

            self.region = 4
        
            self.p = p
            self.T = T
            self.x = x
        
            if np.isnan(p): # Valid region checked inside p_sat_t...
                u_l = np.nan  
            
            # Check if temperature is below region 3 boundary
            elif T < global_property.T_boundary_13:
                # Use region 1 and 2 to calculate saturation properties
                
                # Saturated liquid properties
                v_l = if97.v_r1(p, T)    # Liquid Specific volume                       (m^3/kg)
                h_l = if97.h_r1(p, T)    # Liquid Specific enthalpy                     (J/kg)
                u_l = if97.u_r1(p, T)    # Liquid Specific internal energy              (J/kg)
                s_l = if97.s_r1(p, T)    # Liquid Specific entropy                      (J/kg/K)
                
                # Saturated steam properties
                v_v = if97.v_r2(p, T)    # Vapor Specific volume                       (m^3/kg)
                h_v = if97.h_r2(p, T)    # Vapor Specific enthalpy                     (J/kg)
                u_v = if97.u_r2(p, T)    # Vapor Specific internal energy              (J/kg)
                s_v = if97.s_r2(p, T)    # Vapor Specific entropy                      (J/kg/K)
      
            # temperature is above region 3 lower boundary, check if its below critical temperature
            elif T < global_property.T_crit:
                # Use region 3 for saturation properties    
                #h_l and h_v is determined from their functions and then used
                # to determine densities for the liquid and vapor respectevily, then both
                # rho and T is known and thus the rest of the properties can be calculated.
                
                h_l = iapwsif97.hL_T(T)
                h_v = iapwsif97.hV_T(T)
            
                # Saturated liquid properties
                v_l = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_l),h_l) # use p from p_3sat_h to ensure region 3 otherwise error
                u_l = if97.u_r3(1.0/v_l, T)
                s_l = if97.s_r3(1.0/v_l, T)
                # Saturated steam properties
                v_v = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_v),h_v) # # use p from p_3sat_h to ensure region 3 otherwise error
                u_v = if97.u_r3(1.0/v_v, T)
                s_v = if97.s_r3(1.0/v_v, T)
            else:
                print(f'Temperature, T = {T:.3f} K, above critical temperature')
                u_l = np.nan
                
            # Mixed 2-phase properties
            self.v = x*v_v + (1-x)*v_l
            self.h = x*h_v + (1-x)*h_l
            self.u = x*u_v + (1-x)*u_l
            self.s = x*s_v + (1-x)*s_l
            
            self.st = iapwsif97.st_t(T)
            
        else:
            print(f'x ={x:.2f} must be within 0.0 <= x <= 1.0')
            
            self.region = -1
        
            self.p = p
            self.T = np.nan
            self.x = np.nan
            self.v = np.nan
            self.h = np.nan
            self.u = np.nan
            self.s = np.nan
            self.st = np.nan
                       
        # N/A for 2-phase
        self.cp = np.nan
        self.cv = np.nan
        self.w  = np.nan # Not implemented, not part of iapwsif97
        self.my = np.nan # Not defined for two-phase
        self.tc = np.nan # Not defined for two-phase
        
        if self.v != 0.0 and self.v != np.nan:
            self.rho = x*(1/v_v) + (1-x)*(1/v_l)
        else:
            self.rho = np.nan
            

    def update_ph(self, p, h):
        """
        Update fluid properties using p,h
        IAPWS-IF97
        
        Determines applicable region and calculates temperature or 
        vapor-mass-fraction. Then calls class update-function update_pt or 
        update_px
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        h:      double  enthalpy (J/kg).

        Returns
        -------
        None.
        """  
        region = if97.find_region_ph(p,h)
        
        if region != 4:
            # Region 1,2,3,5, Not saturated
            iapwsif97.update_pt(self,p,iapwsif97.T_ph(p,h))
        else:
            # Region 4, saturation properties, 2-phase            
    
            # Uppdate properties using p,x
            iapwsif97.update_px(self, p, iapwsif97.x_ph(p,h))
            
    def update_ps(self, p, s):
        """
        Update fluid properties using p,s
        IAPWS-IF97
        
        Determines applicable region and calculates temperature using
        backward function T_ps. Then calls class update-function
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        s:      double  entropy (J/kg/K).
    
        Returns
        -------
        None.
        """  
        region = if97.find_region_ps(p, s)
        
        if region != 4:
            # Region 1,2,3,5, Not saturated
            iapwsif97.update_pt(self,p,iapwsif97.T_ps(p,s))
        else:
            # Region 4, saturation properties, 2-phase
            
            # Uppdate properties using p,x
            iapwsif97.update_px(self, p, iapwsif97.x_ps(p,s))
        
    def update_hs(self, h, s):
        """
        Update fluid properties using h,s
        IAPWS-IF97
        
        Determines applicable region and calculates temperature using
        backward function P_hs and T_hs. Then calls class update-function
    
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).
    
        Returns
        -------
        None.
        """  
        
        p = iapwsif97.p_hs(h, s)
        T = iapwsif97.T_hs(h, s)
        
        region = if97.find_region_pt(p, T)
        
        if region != 4:
            # Region 1,2,3,5, Not saturated
            iapwsif97.update_pt(self,p,T)
        else:
            # Region 4, saturation properties, 2-phase
            
            # Uppdate properties using p,x
            iapwsif97.update_px(self, p, iapwsif97.x_ph(p,h))
    
        
    ### Properties as a function of p,T ###
        
    @staticmethod
    def v_pt(p, T):
        """
        Specific volume IAPWS-IF97
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        T:      double  temperature (K).

        Returns
        -------
        double  specific volume (m^3/kg)
        """
        
        region = if97.find_region(p, T)
        
        if region == -1:
            print(f'v_pt, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T-273.15):.3f} C')
            v = np.nan
        elif region == 1:
            v = if97.v_r1(p, T)
        elif region == 2:
            v = if97.v_r2(p, T)
        elif region == 3:
            v = if97_vpt3.v_pt(p, T)
            if rho_r3_itt: # if this option is set, itterate a more precise rho
                v = 1/if97.rho_r3(p, T, 1/v)
        elif region == 4:
            # Saturation properties. 
            # Will not be entered this way
            print('Error, region 4')
        elif region == 5:
            v = if97.v_r5(p, T)
        
        return v
    
    @staticmethod
    def h_pt(p, T):
        """
        Specific enthalpy IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        T:      double  temperature (K).
        
        Returns
        -------
        double  specific enthalpy (J/kg)    
        """
        region = if97.find_region(p, T)
        
        if region == -1:
            print(f'h_pt, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T-273.15):.3f} C')
            h = np.nan
        elif region == 1:
            h = if97.h_r1(p, T)
        elif region == 2:
            h = if97.h_r2(p, T)
        elif region == 3:
            rho = 1 / if97_vpt3.v_pt(p, T)
            if rho_r3_itt: # if this option is set, itterate a more precise rho
                rho = if97.rho_r3(p, T, rho)
            h = if97.h_r3(rho, T)
        elif region == 4:
            # Saturation properties. 
            # Will not be entered this way
            print('Error, region 4')
        elif region == 5:
            h = if97.h_r5(p, T)
        
        return h
    
    @staticmethod
    def u_pt(p, T):
        """
        Specific internal energy IAPWS-IF97
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        T:      double  temperature (K).

        Returns
        -------
        double  specific internal energy (J/kg)    
        """
        
        region = if97.find_region(p, T)
        
        if region == -1:
            print(f'u_pt, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T-273.15):.3f} C')
            u = np.nan
        elif region == 1:
            u = if97.u_r1(p, T)
        elif region == 2:
            u = if97.u_r2(p, T)
        elif region == 3:
            rho = 1 / if97_vpt3.v_pt(p, T)
            if rho_r3_itt: # if this option is set, itterate a more precise rho
                rho = if97.rho_r3(p, T, rho)
            u = if97.u_r3(rho, T)
        elif region == 4:
            # Saturation properties. 
            # Will not be entered this way
            print('Error, region 4')
        elif region == 5:
            u = if97.u_r5(p, T)
        
        return u
    
    @staticmethod
    def s_pt(p, T):
        """
        Specific entropy IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        T:      double  temperature (K).
        
        Returns
        -------
        double  specific entropy (J/kg/K)
        """
    
        region = if97.find_region(p, T)
        
        if region == -1:
            print(f's_pt, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T-273.15):.3f} C')
            s = np.nan
        elif region == 1:
            s = if97.s_r1(p, T)
        elif region == 2:
            s = if97.s_r2(p, T)
        elif region == 3:
            rho = 1 / if97_vpt3.v_pt(p, T)
            if rho_r3_itt: # if this option is set, itterate a more precise rho
                rho = if97.rho_r3(p, T, rho)
            s = if97.s_r3(rho, T)
        elif region == 4:
            # Saturation properties. 
            # Will not be entered this way
            print('Error, region 4')
        elif region == 5:
            s = if97.s_r5(p, T)
        
        return s
    
    @staticmethod
    def cp_pt(p, T):
        """
        Specific isobaric heat capacity IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        T:      double  temperature (K).
        
        Returns
        -------
        double  specific isobaric heat capacity (J/kg/K)
        """
    
        region = if97.find_region(p, T)

        if region == -1:
            print(f'cp_pt, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T-273.15):.3f} C')
            cp = np.nan
        elif region == 1:
            cp = if97.cp_r1(p, T)
        elif region == 2:
            cp = if97.cp_r2(p, T)
        elif region == 3:
            rho = 1 / if97_vpt3.v_pt(p, T)
            if rho_r3_itt: # if this option is set, itterate a more precise rho
                rho = if97.rho_r3(p, T, rho)
            cp = if97.cp_r3(rho, T)
        elif region == 4:
            # Saturation properties. 
            # Will not be entered this way
            print('Error, region 4')
        elif region == 5:
            cp = if97.cp_r5(p, T)
        
        return cp
    
    @staticmethod
    def cv_pt(p, T):
        """
        Specific isochoric heat capacity IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        T:      double  temperature (K).
        
        Returns
        -------
        double  specific isochoric heat capacity (J/kg/K)
        """
        
        region = if97.find_region(p, T)
        
        if region == -1:
            print(f'cv_pt, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T-273.15):.3f} C')
            cv = np.nan
        elif region == 1:
            cv = if97.cv_r1(p, T)
        elif region == 2:
            cv = if97.cv_r2(p, T)
        elif region == 3:
            rho = 1 / if97_vpt3.v_pt(p, T)
            if rho_r3_itt: # if this option is set, itterate a more precise rho
                rho = if97.rho_r3(p, T, rho)
            cv = if97.cv_r3(rho, T)
        elif region == 4:
            # Saturation properties. 
            # Will not be entered this way
            print('Error, region 4')
        elif region == 5:
            cv = if97.cv_r5(p, T)
        
        return cv
    
    @staticmethod
    def w_pt(p, T):
        """
        Speed of sound IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        T:      double  temperature (K).
        
        Returns
        -------
        double  Speed of sound (m/s)
        """
        
        region = if97.find_region(p, T)
        
        if region == -1:
            print(f'w_pt, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T-273.15):.3f} C')
            w = np.nan
        elif region == 1:
            w = if97.w_r1(p, T)
        elif region == 2:
            w = if97.w_r2(p, T)
        elif region == 3:
            rho = 1 / if97_vpt3.v_pt(p, T)
            if rho_r3_itt: # if this option is set, itterate a more precise rho
                rho = if97.rho_r3(p, T, rho)
            w = if97.w_r3(rho, T)
        elif region == 4:
            # Saturation properties. 
            # Will not be entered this way
            print('Error, region 4')
        elif region == 5:
            w = if97.w_r5(p, T)
        
        return w
    
    @staticmethod
    def rho_pt(p, T):
        """
        Density IAPWS-IF97
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        T:      double  temperature (K).

        Returns
        -------
        double  density (kg/m^3)
        """
        
        region = if97.find_region(p, T)
        
        if region == -1:
            print(f'rho_pt, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T-273.15):.3f} C')
            v = np.nan
        elif region == 1:
            v = if97.v_r1(p, T)
        elif region == 2:
            v = if97.v_r2(p, T)
        elif region == 3:
            v = if97_vpt3.v_pt(p, T)
            if rho_r3_itt: # if this option is set, itterate a more precise v
                v = 1/if97.rho_r3(p, T, 1/v)
        elif region == 4:
            # Saturation properties. 
            # Will not be entered this way
            print('Error, region 4')
        elif region == 5:
            v = if97.v_r5(p, T)
        
        if v != 0.0 and v != np.nan:
            rho = 1/v
        else:
            rho = np.nan
            
        return rho
    
    @staticmethod
    def my_pt(p,T):
        """
        Viscosity IAPWSF08
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        T:      double  temperature (K).

        Returns
        -------
        double  viscosity (Pa*s)
        """
        my = 0.0
                
        # Check if within valid region
        if f08.region_my_pt(p, T) :
            
            rho = iapwsif97.rho_pt(p,T)
            
            # Check if critical enhancement is needed or not...
            if f08.crit_enhance_trho(T,rho):
                # Derivatives used in critical enhancement is based on iapws95, iapwsif97 has a lower accuracy, espacially on the derivative
                # These will therefore not match with validation tables in the iapwsf08 release. The release however suggest to skip the
                # critical enhancement all together for greater speed. It is still used in this implementation but will give a warning.
                
                # critical enhancement is only active within region 3.
                #p_rho = lambda rho: if97.p_r3(rho, T) 
                
                # determine derivative drho/dp at constant temperature
                dp = 100 # Pa, discrete length of derivative
                drhodpT = (iapwsif97.rho_pt(p+dp/2,T)-iapwsif97.rho_pt(p-dp/2,T))/dp
                
                # drho/dp at reference temperature
                #Tr = global_property.Tref_f08 #K
                #drhodpTr = (iapwsif97.rho_pt(p+dp/2,Tr)-iapwsif97.rho_pt(p-dp/2,Tr))/dp
                
                drhodpTr = f11.sigma_tr(rho) # use function provided in iapwsf11 instead
                
                my = f08.my_trhosigma(T, rho, drhodpT, drhodpTr)
                
            else:
                my = f08.my_trho(T, rho)
                    
        else:
            print(f'my_pt, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T):.3f} K')
            my = np.nan
            
        return my
    
    @staticmethod
    def tc_pt(p,T):
        """
        Thermal conductivity IAPWSF11
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        T:      double  temperature (K).

        Returns
        -------
        double  thermal conductivity (W/m/K)
        """
        tc = 0.0
        
        # Check if within valid region
        if f11.region_tc_pt(p, T) :

            rho = iapwsif97.rho_pt(p,T)
            cp = iapwsif97.cp_pt(p, T)
            cv = iapwsif97.cv_pt(p, T)
            my = f08.my_trho(T, rho) # calculate my without critical enhancement as is done in iapwsif11 verifications

            # determine derivative drho/dp at constant temperature
            dp = 100 # Pa, discrete length of derivative
            
            if p == global_property.p_r123_upper_boundary or p == global_property.p_r5_upper_boundary and T > global_property.T_boundary_25:
                # On region boundary, do backward derivative
                drhodpT = (iapwsif97.rho_pt(p,T)-iapwsif97.rho_pt(p-dp,T))/dp
            elif p == global_property.p_min:
                # On region boundary, do forward derivative
                drhodpT = (iapwsif97.rho_pt(p+dp,T)-iapwsif97.rho_pt(p,T))/dp
            else:
                # Centeral derivative
                drhodpT = (iapwsif97.rho_pt(p+dp/2,T)-iapwsif97.rho_pt(p-dp/2,T))/dp
                
            # drho/dp at reference temperature
            drhodpTr = f11.sigma_tr(rho)
                
            tc = f11.lambda_tc(T, rho, cp, cv, my, drhodpT, drhodpTr)
                    
        else:
            print(f'tc_pt, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T):.3f} K')
            tc = np.nan
            
        return tc
        
    
    ### Properties as a function of p,h ###
    
    @staticmethod
    def v_ph(p,h):
        """
        Specific volume IAPWS-IF97
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        h:      double  enthalpy (J/kg).

        Returns
        -------
        double  specific volume (m^3/kg)
        """
        region = if97.find_region_ph(p, h)
        
        T = 0.0
        
        if region == -1:
            print(f'v_ph, out of bounds. p = {p*10**-6:.3f} MPa, h = {h*10**-3:.3f} kJ/kg')
            v = np.nan
        elif region == 1:
            T = if97.t_ph_r1(p, h)
            v = if97.v_r1(p, T)
        elif region == 2:
            T = if97.t_ph_r2(p, h)
            v = if97.v_r2(p, T)
        elif region == 3:
            v = if97_tvphps.v_3_ph(p, h)
        elif region == 4:
            # Saturation properties.
            # Determine mass fration vapor
            x = iapwsif97.x_ph(p, h)
            
            T = if97.t_sat_p_r4(p) # # Saturation temperature (K)
            
            # Get liquid and vapor properties
            v_l = if97.v_r1(p, T)    # Liquid Specific volume          (m^3/kg)
            v_v = if97.v_r2(p, T)    # Vapor Specific volume           (m^3/kg)
            
            # Calculate mix
            v = x*v_v + (1-x)*v_l
        elif region == 5:
            T = if97.t_ph_r5(p, h)
            v = if97.v_r5(p, T)
        
        return v
    
    @staticmethod
    def u_ph(p,h):
        """
        Specific internal energy IAPWS-IF97
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        h:      double  enthalpy (J/kg).

        Returns
        -------
        double  specific internal energy (J/kg)    
        """
        region = if97.find_region_ph(p, h)
        
        T = 0.0
        
        if region == -1:
            print(f'u_ph, out of bounds. p = {p*10**-6:.3f} MPa, h = {h*10**-3:.3f} kJ/kg')
            u = np.nan
        elif region == 1:
            T = if97.t_ph_r1(p, h)
            u = if97.u_r1(p, T)
        elif region == 2:
            T = if97.t_ph_r2(p, h)
            u = if97.u_r2(p, T)
        elif region == 3:
            v = if97_tvphps.v_3_ph(p, h)
            T = if97_tvphps.t_3_ph(p, h)
            
            u = if97.u_r3(1/v, T)
        elif region == 4:
            # Saturation properties.
            # Determine mass fration vapor
            x = iapwsif97.x_ph(p, h)
            
            T = if97.t_sat_p_r4(p) # # Saturation temperature (K)
            
            # Get liquid and vapor properties
            u_l = if97.u_r1(p, T)    # Liquid Specific internal energy   (J/kg)
            u_v = if97.u_r2(p, T)    # Liquid Specific internal energy   (J/kg)
            
            # Calculate mix
            u = x*u_v + (1-x)*u_l
        elif region == 5:
            T = if97.t_ph_r5(p, h)
            u = if97.u_r5(p, T)
        
        return u
    
    @staticmethod
    def s_ph(p,h):
        """
        Specific entropy IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        h:      double  enthalpy (J/kg).
        
        Returns
        -------
        double  specific entropy (J/kg/K)
        """
        region = if97.find_region_ph(p, h)
        
        T = 0.0
        
        if region == -1:
            print(f's_ph, out of bounds. p = {p*10**-6:.3f} MPa, h = {h*10**-3:.3f} kJ/kg')
            s = np.nan
        elif region == 1:
            T = if97.t_ph_r1(p, h)
            s = if97.s_r1(p, T)
        elif region == 2:
            T = if97.t_ph_r2(p, h)
            s = if97.s_r2(p, T)
        elif region == 3:
            v = if97_tvphps.v_3_ph(p, h)
            T = if97_tvphps.t_3_ph(p, h)
            
            s = if97.s_r3(1/v, T)
        elif region == 4:
            # Saturation properties.
            # Determine mass fration vapor
            x = iapwsif97.x_ph(p, h)
            
            T = if97.t_sat_p_r4(p) # # Saturation temperature (K)
            
            # Get liquid and vapor properties
            s_l = if97.s_r1(p, T)    # Liquid Specific entropy        (J/kg/K)
            s_v = if97.s_r2(p, T)    # Vapor Specific entropy         (J/kg/K)
            
            # Calculate mix
            s = x*s_v + (1-x)*s_l
        elif region == 5:
            T = if97.t_ph_r5(p, h)
            s = if97.s_r5(p, T)
        
            
        return s
    
    @staticmethod
    def cp_ph(p,h):
        """
        Specific isobaric heat capacity IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        h:      double  enthalpy (J/kg).
        
        Returns
        -------
        double  specific isobaric heat capacity (J/kg/K)
        """
        region = if97.find_region_ph(p, h)
        
        T = 0.0
        
        if region == -1:
            print(f'cp_ph, out of bounds. p = {p*10**-6:.3f} MPa, h = {h*10**-3:.3f} kJ/kg')
            cp = np.nan
        elif region == 1:
            T = if97.t_ph_r1(p, h)
            cp = if97.cp_r1(p, T)
        elif region == 2:
            T = if97.t_ph_r2(p, h)
            cp = if97.cp_r2(p, T)
        elif region == 3:
            v = if97_tvphps.v_3_ph(p, h)
            T = if97_tvphps.t_3_ph(p, h)
            
            cp = if97.cp_r3(1/v, T)
        elif region == 4:
            # Saturation properties.
            cp = 0.0 # N/A for 
            print('Warning cp not applicable for saturated mixtures')
        elif region == 5:
            T = if97.t_ph_r5(p, h)
            cp = if97.cp_r5(p, T)
    
        return cp
    
    @staticmethod
    def cv_ph(p,h):
        """
        Specific isochoric heat capacity IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        h:      double  enthalpy (J/kg).
        
        Returns
        -------
        double  specific isochoric heat capacity (J/kg/K)
        """
        region = if97.find_region_ph(p, h)
        
        T = 0.0
        
        if region == -1:
            print(f'cv_ph, out of bounds. p = {p*10**-6:.3f} MPa, h = {h*10**-3:.3f} kJ/kg')
            cv = np.nan
        elif region == 1:
            T = if97.t_ph_r1(p, h)
            cv = if97.cv_r1(p, T)
        elif region == 2:
            T = if97.t_ph_r2(p, h)
            cv = if97.cv_r2(p, T)
        elif region == 3:
            v = if97_tvphps.v_3_ph(p, h)
            T = if97_tvphps.t_3_ph(p, h)
            
            cv = if97.cv_r3(1/v, T)
        elif region == 4:
            # Saturation properties.
            cv = np.nan
            print('Warning cv not applicable for saturated mixtures')
        elif region == 5:
            T = if97.t_ph_r5(p, h)
            cv = if97.cv_r5(p, T)
        
        return cv
    
    @staticmethod
    def w_ph(p,h):
        """
        Speed of sound IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        h:      double  enthalpy (J/kg).
        
        Returns
        -------
        double  Speed of sound (m/s)
        """
        region = if97.find_region_ph(p, h)
        
        T = 0.0
        
        if region == -1:
            print(f'w_ph, out of bounds. p = {p*10**-6:.3f} MPa, h = {h*10**-3:.3f} kJ/kg')
            w = np.nan
        elif region == 1:
            T = if97.t_ph_r1(p, h)
            w = if97.w_r1(p, T)
        elif region == 2:
            T = if97.t_ph_r2(p, h)
            w = if97.w_r2(p, T)
        elif region == 3:
            v = if97_tvphps.v_3_ph(p, h)
            T = if97_tvphps.t_3_ph(p, h)
            
            w = if97.w_r3(1/v, T)
        elif region == 4:
            # Saturation properties.
            w = np.nan
            print('Warning, speed of sound for 2-phase not implemented')
        elif region == 5:
            T = if97.t_ph_r5(p, h)
            w = if97.w_r5(p, T)
            
        return w
    
    @staticmethod
    def rho_ph(p,h):
        """
        Density IAPWS-IF97
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        h:      double  enthalpy (J/kg).

        Returns
        -------
        double  density (kg/m^3)
        """
        
        v = iapwsif97.v_ph(p, h)
        if v != 0.0 and v != np.nan:
            rho = 1/v
        else:
            rho = 0.0
            
        return rho
    
    @staticmethod
    def my_ph(p,h):
        """
        Viscosity IAPWSF08
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        h:      double  enthalpy  (kJ/kg).

        Returns
        -------
        double  viscosity (Pa*s)
        """
        my = 0.0
        
        region = if97.find_region_ph(p, h)
        
        if region !=4:
            T = iapwsif97.T_ph(p,h)
            my = iapwsif97.my_pt(p, T)
        else:
            print('Warning, Viscosity is not defined for two-phase fluid')
            my = np.nan
            
        return my
    
    @staticmethod
    def tc_ph(p,h):
        """
        Thermal conductivity IAPWSF11
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        h:      double  enthalpy  (kJ/kg).

        Returns
        -------
        double  Thermal conductivity (W/m/K)
        """
        tc = 0.0
        
        region = if97.find_region_ph(p, h)
        
        if region !=4:
            T = iapwsif97.T_ph(p,h)
            tc = iapwsif97.tc_pt(p, T)
        else:
            print('Warning, Thermal conductivity is not defined for two-phase fluid')
            tc = np.nan
        
        return tc
    
    ### Properties as a function of p,s ###
    
    @staticmethod
    def v_ps(p,s):
        """
        Specific volume IAPWS-IF97
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        s:      double  entropy (J/kg/K).

        Returns
        -------
        double  specific volume (m^3/kg)
        """
        region = if97.find_region_ps(p, s)
        
        T = 0.0
        
        if region == -1:
            print(f'v_ps, out of bounds. p = {p*10**-6:.3f} MPa, s = {s*10**-3:.3f} kJ/kg/K')
            v = np.nan
        elif region == 1:
            T = if97.t_ps_r1(p, s)
            v = if97.v_r1(p, T)
        elif region == 2:
            T = if97.t_ps_r2(p, s)
            v = if97.v_r2(p, T)
        elif region == 3:
            v = if97_tvphps.v_3_ps(p, s)
        elif region == 4:
            # Saturation properties.
            # Determine mass fration vapor
            x = iapwsif97.x_ps(p, s)
            
            T = if97.t_sat_p_r4(p) # # Saturation temperature (K)
            
            # Get liquid and vapor properties
            v_l = if97.v_r1(p, T)    # Liquid Specific volume          (m^3/kg)
            v_v = if97.v_r2(p, T)    # Vapor Specific volume           (m^3/kg)
            
            # Calculate mix
            v = x*v_v + (1-x)*v_l
        elif region == 5:
            T = if97.t_ps_r5(p, s)
            v = if97.v_r5(p, T)
        
        return v
    
    
    @staticmethod
    def u_ps(p,s):
        """
        Specific internal energy IAPWS-IF97
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        s:      double  entropy (J/kg/K).

        Returns
        -------
        double  specific internal energy (J/kg)    
        """
        region = if97.find_region_ps(p, s)
        
        T = 0.0
        
        if region == -1:
            print(f'u_ps, out of bounds. p = {p*10**-6:.3f} MPa, s = {s*10**-3:.3f} kJ/kg/K')
            u = np.nan
        elif region == 1:
            T = if97.t_ps_r1(p, s)
            u = if97.u_r1(p, T)
        elif region == 2:
            T = if97.t_ps_r2(p, s)
            u = if97.u_r2(p, T)
        elif region == 3:
            v = if97_tvphps.v_3_ps(p, s)
            T = if97_tvphps.t_3_ps(p, s)
            
            u = if97.u_r3(1/v, T)
        elif region == 4:
            # Saturation properties.
            # Determine mass fration vapor
            x = iapwsif97.x_ps(p, s)
            
            T = if97.t_sat_p_r4(p) # # Saturation temperature (K)
            
            # Get liquid and vapor properties
            u_l = if97.u_r1(p, T)    # Liquid Specific internal energy              (J/kg)
            u_v = if97.u_r2(p, T)    # Vapor Specific internal energy              (J/kg)
            
            # Calculate mix
            u = x*u_v + (1-x)*u_l
        elif region == 5:
            T = if97.t_ps_r5(p, s)
            u = if97.u_r5(p, T)
        
        return u
    
    @staticmethod
    def h_ps(p,s):
        """
        Specific enthalpy IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  specific enthalpy (J/kg)
        """
        region = if97.find_region_ps(p, s)
        
        T = 0.0
        
        if region == -1:
            print(f'h_ps, out of bounds. p = {p*10**-6:.3f} MPa, s = {s*10**-3:.3f} kJ/kg/K')
            h = np.nan
        elif region == 1:
            T = if97.t_ps_r1(p, s)
            h = if97.h_r1(p, T)
        elif region == 2:
            T = if97.t_ps_r2(p, s)
            h = if97.h_r2(p, T)
        elif region == 3:
            v = if97_tvphps.v_3_ps(p, s)
            T = if97_tvphps.t_3_ps(p, s)
            
            h = if97.h_r3(1/v, T)
        elif region == 4:
            # Saturation properties.
            # Determine mass fration vapor
            x = iapwsif97.x_ps(p, s)
            
            T = if97.t_sat_p_r4(p) # # Saturation temperature (K)
            
            # Get liquid and vapor properties
            h_l = if97.h_r1(p, T)    # Liquid Specific enthalpy                     (J/kg)
            h_v = if97.h_r2(p, T)    # Vapor Specific enthalpy                     (J/kg)
            
            # Calculate mix
            h = x*h_v + (1-x)*h_l
        elif region == 5:
            T = if97.t_ps_r5(p, s)
            h = if97.h_r5(p, T)
            
        return h
    
    @staticmethod
    def cp_ps(p,s):
        """
        Specific isobaric heat capacity IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  specific isobaric heat capacity (J/kg/K)
        """
        region = if97.find_region_ps(p, s)
        
        T = 0.0
        
        if region == -1:
            print(f'cp_ps, out of bounds. p = {p*10**-6:.3f} MPa, s = {s*10**-3:.3f} kJ/kg/K')
            cp = np.nan
        elif region == 1:
            T = if97.t_ps_r1(p, s)
            cp = if97.cp_r1(p, T)
        elif region == 2:
            T = if97.t_ps_r2(p, s)
            cp = if97.cp_r2(p, T)
        elif region == 3:
            v = if97_tvphps.v_3_ps(p, s)
            T = if97_tvphps.t_3_ps(p, s)
            
            cp = if97.cp_r3(1/v, T)
        elif region == 4:
            # Saturation properties.
            cp = np.nan # N/A for 
            print('Warning cp not applicable for saturated mixtures')
        elif region == 5:
            T = if97.t_ps_r5(p, s)
            cp = if97.cp_r5(p, T)
    
        return cp
    
    @staticmethod
    def cv_ps(p,s):
        """
        Specific isochoric heat capacity IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  specific isochoric heat capacity (J/kg/K)
        """
        region = if97.find_region_ps(p, s)
        
        T = 0.0
        
        if region == -1:
            print(f'cv_ps, out of bounds. p = {p*10**-6:.3f} MPa, s = {s*10**-3:.3f} kJ/kg/K')
            cv = np.nan
        elif region == 1:
            T = if97.t_ps_r1(p, s)
            cv = if97.cv_r1(p, T)
        elif region == 2:
            T = if97.t_ps_r2(p, s)
            cv = if97.cv_r2(p, T)
        elif region == 3:
            v = if97_tvphps.v_3_ps(p, s)
            T = if97_tvphps.t_3_ps(p, s)
            
            cv = if97.cv_r3(1/v, T)
        elif region == 4:
            # Saturation properties.
            cv = np.nan # N/A for 
            print('Warning cv not applicable for saturated mixtures')
        elif region == 5:
            T = if97.t_ps_r5(p, s)
            cv = if97.cv_r5(p, T)
        
        return cv
    
    @staticmethod
    def w_ps(p,s):
        """
        Speed of sound IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  Speed of sound (m/s)
        """
        region = if97.find_region_ps(p, s)
        
        T = 0.0
        
        if region == -1:
            print(f'w_ps, out of bounds. p = {p*10**-6:.3f} MPa, s = {s*10**-3:.3f} kJ/kg/K')
            w = np.nan
        elif region == 1:
            T = if97.t_ps_r1(p, s)
            w = if97.w_r1(p, T)
        elif region == 2:
            T = if97.t_ps_r2(p, s)
            w = if97.w_r2(p, T)
        elif region == 3:
            v = if97_tvphps.v_3_ps(p, s)
            T = if97_tvphps.t_3_ps(p, s)
            
            w = if97.w_r3(1/v, T)
        elif region == 4:
            # Saturation properties.
            w = np.nan
            print('Warning, speed of sound for 2-phase not implemented')
        elif region == 5:
            T = if97.t_ps_r5(p, s)
            w = if97.w_r5(p, T)
            
        return w
    
    @staticmethod
    def rho_ps(p,s):
        """
        Density IAPWS-IF97
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        s:      double  entropy (J/kg/K).

        Returns
        -------
        double  density (kg/m^3)
        """
        
        v = iapwsif97.v_ps(p, s)
        
        if v != 0.0 and v != np.nan:
            rho = 1/v
        else:
            rho = 0.0
            
        return rho
    
    @staticmethod
    def my_ps(p,s):
        """
        Viscosity IAPWSF08
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        s:      double  entropy  (kJ/kg).

        Returns
        -------
        double  viscosity (Pa*s)
        """
        my = 0.0
        
        region = if97.find_region_ps(p, s)
        
        if region !=4:
            T = iapwsif97.T_ps(p,s)
            my = iapwsif97.my_pt(p, T)
        else:
            print('Warning, Viscosity is not defined for two-phase fluid')
            my = np.nan
            
        return my

    @staticmethod
    def tc_ps(p,s):
        """
        Thermal conductivity IAPWSF11
    
        Parameters
        ----------
        p:      double  pressure (Pa).
        s:      double  entropy  (kJ/kg).

        Returns
        -------
        double  Thermal conductivity (W/m/K)
        """
        tc = 0.0
        
        region = if97.find_region_ps(p, s)
        
        if region !=4:
            T = iapwsif97.T_ps(p,s)
            tc = iapwsif97.tc_pt(p, T)
        else:
            print('Warning, Thermal conductivity is not defined for two-phase fluid')
            tc = np.nan
            
        return tc
    
    ### Properties as a function of h,s ###
        
    @staticmethod
    def v_hs(h, s):
        """
        Specific volume IAPWS-IF97
    
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).

        Returns
        -------
        double  specific volume (m^3/kg)
        """
        
        p = iapwsif97.p_hs(h, s)
        T = iapwsif97.T_hs(h, s)
        
        return iapwsif97.v_pt(p, T)

    @staticmethod
    def u_hs(h, s):
        """
        Specific internal energy IAPWS-IF97
    
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).

        Returns
        -------
        double  specific internal energy (J/kg)    
        """
        
        p = iapwsif97.p_hs(h, s)
        T = iapwsif97.T_hs(h, s)
        
        return iapwsif97.u_pt(p, T)
    
    @staticmethod
    def cp_hs(h, s):
        """
        Specific isobaric heat capacity IAPWS-IF97
        
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  specific isobaric heat capacity (J/kg/K)
        """
        
        p = iapwsif97.p_hs(h, s)
        T = iapwsif97.T_hs(h, s)
        
        return iapwsif97.cp_pt(p, T)
    
    @staticmethod
    def cv_hs(h, s):
        """
        Specific isochoric heat capacity IAPWS-IF97
        
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  specific isochoric heat capacity (J/kg/K)
        """

        p = iapwsif97.p_hs(h, s)
        T = iapwsif97.T_hs(h, s)
        
        return iapwsif97.cv_pt(p, T)
    
    @staticmethod
    def w_hs(h, s):
        """
        Speed of sound IAPWS-IF97
        
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  Speed of sound (m/s)
        """

        p = iapwsif97.p_hs(h, s)
        T = iapwsif97.T_hs(h, s)
        
        return iapwsif97.w_pt(p, T)
    
    @staticmethod
    def rho_hs(h, s):
        """
        Density IAPWS-IF97
    
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).

        Returns
        -------
        double  density (kg/m^3)
        """
        
        p = iapwsif97.p_hs(h, s)
        T = iapwsif97.T_hs(h, s)
        
        return iapwsif97.rho_pt(p, T)
    
    @staticmethod
    def my_hs(h, s):
        """
        VIscosity IAPWSF08
    
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).

        Returns
        -------
        double  viscosity (Pa*s)
        """
        
        p = iapwsif97.p_hs(h, s)
        T = iapwsif97.T_hs(h, s)
        
        return iapwsif97.my_pt(p, T)
    
    @staticmethod
    def tc_hs(h, s):
        """
        Thermal conductivity IAPWSF11
    
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).

        Returns
        -------
        double  Thermal conductivity (W/m/K)
        """
        
        p = iapwsif97.p_hs(h, s)
        T = iapwsif97.T_hs(h, s)
        
        return iapwsif97.tc_pt(p, T)
    
    ### Properties as a function of p,x ###
    
    @staticmethod
    def v_px(p,x):
        """
        Specific volume for pressure and vapor mass fraction p,x

        Parameters
        ----------
        p:      double  pressure (Pa).
        x:      double  vapor mass fraction (-)

        Returns
        -------
        double  Specific volume (m^3/kg)

        """
        v = np.nan
        
        # mass fraction must be within 0.0<=x<=1.0
        if x >= 0.0 and x <= 1.0:
            vV = iapwsif97.vV_p(p)
            vL = iapwsif97.vL_p(p)
        
            if vV != np.nan and vL != np.nan:
                v = x*vV + (1-x)*vL
        else:
            print(f'x ={x:.2f} must be within 0.0 <= x <= 1.0')
            
        return v
        
    @staticmethod
    def h_px(p,x):
        """
        Enthalpy for pressure and vapor mass fraction p,x

        Parameters
        ----------
        p:      double  pressure (Pa).
        x:      double  vapor mass fraction (-)

        Returns
        -------
        double  Enthalpy (J/kg)

        """
        h = np.nan
        
        # mass fraction must be within 0.0<=x<=1.0
        if x >= 0.0 and x <= 1.0:
            hV = iapwsif97.hV_p(p)
            hL = iapwsif97.hL_p(p)
        
            if hV != np.nan and hL != np.nan:
                h = x*hV + (1-x)*hL
        else:
            print(f'x ={x:.2f} must be within 0.0 <= x <= 1.0')
            
        return h
    
    @staticmethod
    def u_px(p,x):
        """
        Internal energy for pressure and vapor mass fraction p,x

        Parameters
        ----------
        p:      double  pressure (Pa).
        x:      double  vapor mass fraction (-)

        Returns
        -------
        double  Internal energy (J/kg)

        """
        u = np.nan
        
        # mass fraction must be within 0.0<=x<=1.0
        if x >= 0.0 and x <= 1.0:
            uV = iapwsif97.uV_p(p)
            uL = iapwsif97.uL_p(p)
        
            if uV != np.nan and uL != np.nan:
                u = x*uV + (1-x)*uL
        else:
            print(f'x ={x:.2f} must be within 0.0 <= x <= 1.0')
            
        return u
    
    @staticmethod
    def s_px(p,x):
        """
        Entropy for temperature and vapor mass fraction T,x

        Parameters
        ----------
        p:      double  pressure (Pa).
        x:      double  vapor mass fraction (-)

        Returns
        -------
        double  Entropy (J/kg/K)

        """
        s = np.nan
        
        # mass fraction must be within 0.0<=x<=1.0
        if x >= 0.0 and x <= 1.0:
            sV = iapwsif97.sV_p(p)
            sL = iapwsif97.sL_p(p)
        
            if sV != np.nan and sL != np.nan:
                s = x*sV + (1-x)*sL
        else:
            print(f'x ={x:.2f} must be within 0.0 <= x <= 1.0')
            
        return s
    
    @staticmethod
    def rho_px(p,x):
        """
        Density for temperature and vapor mass fraction T,x

        Parameters
        ----------
        p:      double  pressure (Pa).
        x:      double  vapor mass fraction (-)

        Returns
        -------
        double  Density (kg/m^3)

        """
        rho = np.nan
        
        # mass fraction must be within 0.0<=x<=1.0
        if x >= 0.0 and x <= 1.0:
            rhoV = iapwsif97.rhoV_p(p)
            rhoL = iapwsif97.rhoL_p(p)
        
            if rhoV != np.nan and rhoL != np.nan:
                rho = x*rhoV + (1-x)*rhoL
        else:
            print(f'x ={x:.2f} must be within 0.0 <= x <= 1.0')
            
        return rho
    
    
    ### Properties as a function of T,x ###
    
    @staticmethod
    def v_tx(T,x):
        """
        Specific volume for temperature and vapor mass fraction T,x

        Parameters
        ----------
        T:      double  temperature (K).
        x:      double  vapor mass fraction (-)

        Returns
        -------
        double  Specific volume (m^3/kg)

        """
        v = np.nan
        
        # mass fraction must be within 0.0<=x<=1.0
        if x >= 0.0 and x <= 1.0:
            vV = iapwsif97.vV_T(T)
            vL = iapwsif97.vL_T(T)
        
            if vV != np.nan and vL != np.nan:
                v = x*vV + (1-x)*vL
        else:
            print(f'x ={x:.2f} must be within 0.0 <= x <= 1.0')
            
        return v
    
    @staticmethod
    def h_tx(T,x):
        """
        Enthalpy for temperature and vapor mass fraction T,x

        Parameters
        ----------
        T:      double  temperature (K).
        x:      double  vapor mass fraction (-)

        Returns
        -------
        double  Enthalpy (J/kg)

        """
        h = np.nan
        
        # mass fraction must be within 0.0<=x<=1.0
        if x >= 0.0 and x <= 1.0:
            hV = iapwsif97.hV_T(T)
            hL = iapwsif97.hL_T(T)
        
            if hV != np.nan and hL != np.nan:
                h = x*hV + (1-x)*hL
        else:
            print(f'x ={x:.2f} must be within 0.0 <= x <= 1.0')
            
        return h
    
    @staticmethod
    def u_tx(T,x):
        """
        Internal energy for temperature and vapor mass fraction T,x

        Parameters
        ----------
        T:      double  temperature (K).
        x:      double  vapor mass fraction (-)

        Returns
        -------
        double  Internal energy (J/kg)

        """
        u = np.nan
        
        # mass fraction must be within 0.0<=x<=1.0
        if x >= 0.0 and x <= 1.0:
            uV = iapwsif97.uV_T(T)
            uL = iapwsif97.uL_T(T)
        
            if uV != np.nan and uL != np.nan:
                u = x*uV + (1-x)*uL
        else:
            print(f'x ={x:.2f} must be within 0.0 <= x <= 1.0')
            
        return u
    
    @staticmethod
    def s_tx(T,x):
        """
        Entropy for temperature and vapor mass fraction T,x

        Parameters
        ----------
        T:      double  temperature (K).
        x:      double  vapor mass fraction (-)

        Returns
        -------
        double  Entropy (J/kg/K)

        """
        s = np.nan
        
        # mass fraction must be within 0.0<=x<=1.0
        if x >= 0.0 and x <= 1.0:
            sV = iapwsif97.sV_T(T)
            sL = iapwsif97.sL_T(T)
        
            if sV != np.nan and sL != np.nan:
                s = x*sV + (1-x)*sL
        else:
            print(f'x ={x:.2f} must be within 0.0 <= x <= 1.0')
            
        return s
    
    @staticmethod
    def rho_tx(T,x):
        """
        Density for temperature and vapor mass fraction T,x

        Parameters
        ----------
        T:      double  temperature (K).
        x:      double  vapor mass fraction (-)

        Returns
        -------
        double  Density (kg/m^3)

        """
        rho = np.nan
        
        # mass fraction must be within 0.0<=x<=1.0
        if x >= 0.0 and x <= 1.0:
            rhoV = iapwsif97.rhoV_T(T)
            rhoL = iapwsif97.rhoL_T(T)
        
            if rhoV != np.nan and rhoL != np.nan:
                rho = x*rhoV + (1-x)*rhoL
        else:
            print(f'x ={x:.2f} must be within 0.0 <= x <= 1.0')
            
        return rho 
    
    ### Temperature, pressure backward functions ###
    
    @staticmethod
    def T_ph(p, h):
        """
        Temperature given p,h IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        h:      double  enthalpy (J/kg).
        
        Returns
        -------
        double  Temperature (K)
        """

        region = if97.find_region_ph(p, h)

        if region == -1:
            print(f'T_ph, out of bounds. p = {p*10**-6:.3f} MPa, h = {h*10**-3:.3f} kJ/kg')
            T = np.nan
        elif region == 1:
            T = if97.t_ph_r1(p, h)
        elif region == 2:
            T = if97.t_ph_r2(p, h)
        elif region == 3:
            # IAPWS SR3-03(2014)
            T = if97_tvphps.t_3_ph(p, h)
        elif region == 4:
            # Saturation properties. 
            T = if97.t_sat_p_r4(p)
        elif region == 5:
            # No backward equation for region 5, itterate
            T = if97.t_ph_r5(p, h)
        
        return T
    
    @staticmethod
    def T_ps(p, s):
        """
        Temperature given p,s IAPWS-IF97
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  Temperature (K)
        """
        
        T = np.nan
        region = if97.find_region_ps(p, s)
        
        if region == -1:
            print(f'T_ps, out of bounds. p = {p*10**-6:.3f} MPa, s = {s*10**-3:.3f} kJ/kg/K')
        elif region == 1:
            T = if97.t_ps_r1(p, s)
        elif region == 2:
            T = if97.t_ps_r2(p, s)
        elif region == 3:
            # IAPWS SR3-03(2014)
            T = if97_tvphps.t_3_ps(p, s)
        elif region == 4:
            # Saturation properties. 
            T = if97.t_sat_p_r4(p)
        elif region == 5:
            # No backward equation for region 5, itterate
            T = if97.t_ps_r5(p, s)
        
        return T
    
    @staticmethod
    def T_hs(h, s):
        """
        Temperature given p,s IAPWS-IF97
        
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  Temperature (K)
        """
        
        T = np.nan
        
        if if97_phs.find_region_hs(h, s) != -1:
            p = iapwsif97.p_hs(h, s)
            T = iapwsif97.T_ph(p, h)
        else:
            print(f'T_hs, out of bounds. h = {h*10**-3:.3f} MPa, s = {s*10**-3:.3f} kJ/kg/K')
        
        return T
    
    @staticmethod
    def Tsat_p(p):
        """
        Saturation Temperature given p 
        IAPWS-IF97 Region 4
        
        Parameters
        ----------
        p:      double  pressure (Pa).
        
        Returns
        -------
        double  Temperature (K)
        """
        
        T = np.nan
        
        # Check validity range
        if p <= global_property.p_crit and p >= global_property.p_min:
            T = if97.t_sat_p_r4(p)
        else:
            print(f'Pressure, p = {(p*10**-6):.3f} MPa outside of valid range {(global_property.p_min*10**-6):.3f} MPa <= p <= {(global_property.p_crit*10**-6):.3f} MPa for saturation properties')
                    
        return T
    
    @staticmethod
    def Tsat_hs(h,s):
        """
        Saturation Temperature given h,s 
        IAPWS-IF97 Region 4
        
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  Temperature (K)
        """
        
        region = if97_phs.find_region_hs(h, s)
        
        if region == 4:
            if s < global_property.s_h_bis_623:
                x_p = lambda p: iapwsif97.x_ps(p, s) - iapwsif97.x_ph(p, h)
                # Sensitive to upper pressure guess in itteration, assume psat_s as upper limit works fine
                p = aux.bisection(x_p, global_property.p_min, iapwsif97.psat_s(s), 0, global_property.itt_tol, global_property.itt_max)
                
                T = iapwsif97.Tsat_p(p)
            else:
                T = if97_phs.T_sat_hs(h,s)
                #if97_tvphps.p_3sat_s
            #T = -1
            #print(f'Specified s < {global_property.s_h_bis_623*10**-3:.3f} kJ/kg/K')
            
        else:
            T = np.nan
            print('Specified h,s not within the saturated region')
        
        return T
    
    @staticmethod
    def psat_t(T):
        """
        Saturation pressure given T 
        IAPWS-IF97 Region 4
        
        Parameters
        ----------
        T:      double  Temperature.
        
        Returns
        -------
        double  pressure (Pa)
        """
        
        p = 0.0
        
        # Check validity range
        if T <= global_property.T_crit and T >= global_property.T_boundary_1:
            p = if97.p_sat_t_r4(T)
        else:
            print(f'Temperature, T = {(T):.2f} K outside of valid range {(global_property.T_boundary_1):.2f} K <= T <= {(global_property.T_crit):.2f} K for saturation properties')
            p = np.nan
            
        return p
    
    """
    Several possible pressures for saturated steam. se p-h diagram.
    Leave out or give warning for saturated vapor
    @staticmethod
    def psat_h(h):
        
        Saturation pressure given h
        IAPWS-IF97 Region 4
        
        Parameters
        ----------
        s:      double  Enthalpy (J/kg).
        
        Returns
        -------
        double  pressure (Pa)
        
        p = 0.0
    
        if h < global_property.h_r3_lower:
            # Region 1
            hsat_r1_p = lambda p: if97.h_r1(p,if97.t_sat_p_r4(p))
            
            # Itterate, first using  a hybrid of the bisection method the 
            # secant method. in this way it always converges and is way faster 
            # then using just the bisection method.
            
            # Intermeadiate tolerance, determines how many itterations with bisection method that is to be used
            tol_inter = 100
            
            p = aux.hybrid_bisec_secant(hsat_r1_p, 1, global_property.p_r3_lower_boundary, h, tol_inter)
            
        elif h <= global_property.h_r3_upper:
            # Region 3
            p = if97_tvphps.p_3sat_h(h)
        else:
            # Region 2
            hsat_r2_p = lambda p: if97.h_r2(p,if97.t_sat_p_r4(p))
            
            # Itterate, first using  a hybrid of the bisection method the 
            # secant method. in this way it always converges and is way faster 
            # then using just the bisection method.
            
            # Intermeadiate tolerance, determines how many itterations with bisection method that is to be used
            tol_inter = 100
            p = aux.hybrid_bisec_secant(hsat_r2_p, global_property.p_r3_lower_boundary, 1, h, tol_inter)
            
        return p"""
    
    @staticmethod
    def psat_s(s):
        """
        Saturation pressure given s
        IAPWS-IF97 Region 4
        
        Parameters
        ----------
        s:      double  Entropy (J/kg/K).
        
        Returns
        -------
        double  pressure (Pa)
        """
        
        p = 0.0
    
        if s < global_property.s_r3_lower:
            # Region 1
            ssat_r1_p = lambda p: if97.s_r1(p,if97.t_sat_p_r4(p))
            
            # Itterate, first using  a hybrid of the bisection method the 
            # secant method. in this way it always converges and is way faster 
            # then using just the bisection method.
            
            # Intermeadiate tolerance, determines how many itterations with bisection method that is to be used
            tol_inter = 100
            
            p = aux.hybrid_bisec_secant(ssat_r1_p, 1, global_property.p_r3_lower_boundary, s, tol_inter)
            
        elif s <= global_property.s_r3_upper:
            # Region 3
            p = if97_tvphps.p_3sat_s(s)
        elif s <= global_property.s_h_bis_273:
            # Region 2
            ssat_r2_p = lambda p: if97.s_r2(p,if97.t_sat_p_r4(p))
            
            # Itterate, first using  a hybrid of the bisection method the 
            # secant method. in this way it always converges and is way faster 
            # then using just the bisection method.
            
            # Intermeadiate tolerance, determines how many itterations with bisection method that is to be used
            tol_inter = 100
            p = aux.hybrid_bisec_secant(ssat_r2_p, global_property.p_r3_lower_boundary, 1, s, tol_inter)
        else:
            print(f'Entropy, s = {s*10**-3:.2f}, above saturated vapor entropy in triple point')
            p = np.nan
            
        return p
    
    @staticmethod
    def psat_hs(h,s):
        """
        Saturation Pressure given h,s 
        IAPWS-IF97 Region 4
        
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  Pressure (Pa)
        """
        
        region = if97_phs.find_region_hs(h, s)
        
        if region == 4:
            if s < global_property.s_h_bis_623:
                x_p = lambda p: iapwsif97.x_ps(p, s) - iapwsif97.x_ph(p, h)
                # Sensitive to upper pressure guess in itteration, assume psat_s as upper limit works fine
                p = aux.bisection(x_p, global_property.p_min, iapwsif97.psat_s(s), 0, global_property.itt_tol, global_property.itt_max)
                
            else:
                T = if97_phs.T_sat_hs(h,s)
                p = iapwsif97.psat_t(T)
        else:
            p = np.nan
            print('Specified h,s not within the saturated region')
        
        return p
    
    
    @staticmethod
    def p_hs(h,s):
        """
        Pressure given h,s IAPWS-IF97 SR2-01 and SR4-04
        
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  pressure (Pa)
        """
        p = np.nan
        region = if97_phs.find_region_hs(h, s)
        
        if region == -1:
            print(f'p_hs, Out of bounds. h = {h*1e-3:.3f} kJ/kg, s = {s*1e-3:.3f} kJ/kg')
            p = np.nan
        elif region == 1:
            p = if97_phs.p_hs_r1(h,s)
        elif region == 2:
            p = if97_phs.p_hs_r2(h,s)
        elif region == 3:
            p = if97_phs.p_hs_r3(h,s)
        elif region == 4:
            p = iapwsif97.psat_hs(h,s)
        elif region == 5:
            s_p = lambda p: if97.s_r5(p, if97.t_ph_r5(p, h))
            p = aux.bisection(s_p, global_property.p_min, global_property.p_r5_upper_boundary, s, global_property.itt_tol, global_property.itt_max)
            
        return p
    
    @staticmethod
    def h_prho(p,rho):
        """
        Enthalpy for pressure and density p,rho

        Parameters
        ----------
        p:      double  pressure (Pa).
        rho:    double  density (kg/m^3)

        Returns
        -------
        double  Enthalpy (J/kg)

        """
        # Make sure given p is within applicable region
        if p <= global_property.p_r5_upper_boundary and p>= global_property.p_min and rho >= iapwsif97.rho_ph(p, if97.h_r5(p,global_property.T_boundary_5)) and rho <= iapwsif97.rho_ph(p, if97.h_r1(p,global_property.T_boundary_1)):
            # Region 1-5 possible
            rho_h = lambda h: iapwsif97.rho_ph(p, h)
        
            h = aux.bisection(rho_h, if97.h_r1(p,global_property.T_boundary_1), if97.h_r5(p,global_property.T_boundary_5), rho, global_property.itt_tol, global_property.itt_max)
        elif p <= global_property.p_r123_upper_boundary and p> global_property.p_r5_upper_boundary and rho >= iapwsif97.rho_ph(p, if97.h_r2(p,global_property.T_boundary_25)) and rho <= iapwsif97.rho_ph(p, if97.h_r1(p,global_property.T_boundary_1)):
            # Region 1-4 posible
            rho_h = lambda h: iapwsif97.rho_ph(p, h)
        
            h = aux.bisection(rho_h, if97.h_r1(p,global_property.T_boundary_1), if97.h_r2(p,global_property.T_boundary_25), rho, global_property.itt_tol, global_property.itt_max)
        else:
            h = np.nan
            print(f'h_prho, out of bounds. p = {p*10**-6:.3f} MPa, rho = {rho:.3f} kg/m^3')
        return h
    
       
    ### Mass fractions ####
    
    @staticmethod
    def x_ph(p,h):
        """
        Vapor mass fraction given p,h
        
        Parameters
        ----------
        Parameters
        ----------
        p:      double  pressure (Pa).
        h:      double  enthalpy (J/kg).
        
        Returns
        -------
        double  vapor-mass-fraction (-)
        """

        # Region 4, saturation properties, 2-phase           
        
        p_crit = global_property.p_crit              # Pa, Critical pressure
        
        if p <= p_crit:
            # Get saturation enthalpies and calculate vapor-mass-fraction, x
            h_l = iapwsif97.hL_p(p) # Saturated liquid enthalpy
            h_v = iapwsif97.hV_p(p) # Saturated steam enthalpy
            
            if h_l != 0.0 and h_v != 0.0 and h_l != np.nan and h_v != np.nan:
                x = (h - h_l)/(h_v - h_l)
            elif h_v >= 1.0:
                x = 1.0
            else:
                x = 0.0
        else:
            # Above critical pressure
            print("Warning, Pressure above critical pressure")
            x = 1.0
        
        # mass fraction should always between 0-1
        if x < 0.0:
            x = 0.0
        elif x > 1.0:
            x = 1.0

        return x
    
    @staticmethod
    def x_ps(p,s):
        """
        Vapor mass fraction given p,s
        
        Parameters
        ----------
        Parameters
        ----------
        p:      double  pressure (Pa).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  vapor-mass-fraction (-)
        """

        # Region 4, saturation properties, 2-phase           
        
        p_crit =  global_property.p_crit              # Pa, Critical pressure
        
        if p <= p_crit:
            # Get saturation enthalpies and calculate vapor-mass-fraction, x
            s_l = iapwsif97.sL_p(p) # Saturated liquid enthalpy
            s_v = iapwsif97.sV_p(p) # Saturated steam enthalpy
            
            if s_l != 0.0 and s_v != 0.0 and s_l != np.nan and s_v != np.nan:
                x = (s - s_l)/(s_v - s_l)
            elif s_v >= 1.0:
                x = 1.0
            else:
                x = 0.0
        else:
            # Above critical pressure
            print("Warning, Pressure above critical pressure")
            x = 1.0
        
        # mass fraction should always between 0-1
        if x < 0.0:
            x = 0.0
        elif x > 1.0:
            x = 1.0

        return x
    
    @staticmethod
    def x_hs(h,s):
        """
        Vapor mass fraction given h,s
        
        Parameters
        ----------
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  vapor-mass-fraction (-)
        """

        # Region 4, saturation properties, 2-phase           
        
        p_crit =  global_property.p_crit              # Pa, Critical pressure
        
        region = if97_phs.find_region_hs(h, s)
        
        if region == 4:
            if s < global_property.s_h_bis_623:
                x_p = lambda p: iapwsif97.x_ps(p, s) - iapwsif97.x_ph(p, h)
                # Sensitive to upper pressure guess in itteration, assume psat_s as upper limit works fine
                p = aux.bisection(x_p, global_property.p_min, iapwsif97.psat_s(s), 0, global_property.itt_tol, global_property.itt_max)
                
            else:
                T = if97_phs.T_sat_hs(h,s)
                p = iapwsif97.psat_t(T)
                                
            if p <= p_crit:
                x = iapwsif97.x_ph(p,h)
            
            else:
                # Above critical pressure
                print("Warning, Pressure above critical pressure")
                x = 1.0
        else:
            if region == 1:
                x = 0.0
            elif region == 2:
                x = 1.0
            else:
                # Region 3, check if above or below saturation line (p-T)
                T = iapwsif97.T_hs(h,s)
                p_sat = iapwsif97.psat_t(T)
                p = iapwsif97.p_hs(h,s)
                
                if p < p_sat:
                    x = 1.0
                else:
                    x = 0.0
        
        # mass fraction should always between 0-1
        if x < 0.0:
            x = 0.0
        elif x > 1.0:
            x = 1.0

        return x
    
    ### Volume vapor fractions ####
    
    @staticmethod
    def vx_ph(p,h):
        """
        Vapor volume fraction given p,h
        
        Parameters
        ----------
        Parameters
        ----------
        p:      double  pressure (Pa).
        h:      double  enthalpy (J/kg).
        
        Returns
        -------
        double  vapor-volume-fraction (-)
        """
        
        # Mass fraction
        x = iapwsif97.x_ph(p,h)
        
        # Densities
        rhoL = iapwsif97.rhoL_p(p)
        rhoV = iapwsif97.rhoV_p(p)
        
        if x != np.nan and rhoL != np.nan and rhoV != np.nan:
            vx = x/(x+(1-x)*(rhoV/rhoL)) 
        else:
            vx = np.nan
        
        return vx
    
    @staticmethod
    def vx_ps(p,s):
        """
        Vapor volume fraction given p,h
        
        Parameters
        ----------
        Parameters
        ----------
        p:      double  pressure (Pa).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  vapor-volume-fraction (-)
        """
        
        # Mass fraction
        x = iapwsif97.x_ps(p,s)
        
        # Densities
        rhoL = iapwsif97.rhoL_p(p)
        rhoV = iapwsif97.rhoV_p(p)
        
        if x != np.nan and rhoL != np.nan and rhoV != np.nan:
            vx = x/(x+(1-x)*(rhoV/rhoL)) 
        else:
            vx = np.nan

        return vx
    
    @staticmethod
    def vx_hs(h,s):
        """
        Vapor volume fraction given h,s
        
        Parameters
        ----------
        Parameters
        ----------
        h:      double  enthalpy (J/kg).
        s:      double  entropy (J/kg/K).
        
        Returns
        -------
        double  vapor-volume-fraction (-)
        """
        
        # Pressure
        p = iapwsif97.p_hs(h, s)
        
        # Mass fraction
        x = iapwsif97.x_hs(h,s)
        
        # Densities
        rhoL = iapwsif97.rhoL_p(p)
        rhoV = iapwsif97.rhoV_p(p)
        
        if x != np.nan and rhoL != np.nan and rhoV != np.nan:
            vx = x/(x+(1-x)*(rhoV/rhoL)) 
        else:
            vx = np.nan

        return vx
    
    
    ### Saturation properties, Liquid and Vapor ###
    ### _pressure ###
    
    @staticmethod
    def hL_p(p):
        """
        Saturated liquid enthalpy for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  Enthalpy (J/kg)

        """
        h_l = 0.0 # (J/kg)
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            h_l = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 1  to calculate liquid saturation properties
            h_l = if97.h_r1(p, T)    # Liquid Specific enthalpy                     (J/kg)

        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties

            # Both saturated liquid and saturated vapor properties is defined by 
            # region 3, itteration in p-h space is done with p_3sat_h(h) which is
            # the saturation line p-h between region 3 and 4. 

            int_tol = 1000
            h_l = aux.hybrid_bisec_secant(if97_tvphps.p_3sat_h, global_property.h_r3_lower , global_property.h_crit , p, int_tol)
            
        else:
            print(f'Warning, Pressure, p = {p*10**-6:.3f} MPa, above critical pressure')
            h_l = np.nan
                        
        return h_l
    
    @staticmethod
    def hV_p(p):
        """
        Saturated vapor enthalpy for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  Enthalpy (J/kg)

        """
        h_v = 0.0 # (J/kg)
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            h_v = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary: # Check if pressure is below region 3 boundary
            # Use region 2  to calculate vapor saturation properties
            h_v = if97.h_r2(p, T)    # Vapor Specific enthalpy                     (J/kg)
            
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties
            
            # Both saturated liquid and saturated vapor properties is defined by 
            # region 3, itteration in p-h space is done with p_3sat_h(h) which is
            # the saturation line p-h between region 3 and 4. 

            int_tol = 1000
            h_v = aux.hybrid_bisec_secant(if97_tvphps.p_3sat_h, global_property.h_crit, global_property.h_r3_upper , p, int_tol)
            
        else:
            print(f'Warning, Pressure, p = {p*10**-6:.3f} MPa, above critical pressure')
            h_v = np.nan
                        
        return h_v
    
    @staticmethod
    def sL_p(p):
        """
        Saturated liquid entropy for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  Entropy (J/kg/K)

        """
        s_l = 0.0 # (J/kg)
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            s_l = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 1  to calculate liquid saturation properties
            s_l = if97.s_r1(p, T)    # Liquid Specific enthalpy                     (J/kg)
      
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties
            
            # Both saturated liquid and saturated vapor properties is defined by 
            # region 3, itteration in p-s space is done with p_3sat_s(s) which is
            # the saturation line p-s between region 3 and 4. 
             
            int_tol = 75
            s_l = aux.hybrid_bisec_secant(if97_tvphps.p_3sat_s, global_property.s_r3_lower, global_property.s_crit, p, int_tol)
            
        else:
            print(f'Pressure, p = {p*10**-6:.3f} MPa, above critical pressure')
            s_l = np.nan
                        
        return s_l

    @staticmethod
    def sV_p(p):
        """
        Saturated vapor entropy for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  Entropy (J/kg/K)

        """
        s_v = 0.0 # (J/kg)
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            s_v = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 1  to calculate liquid saturation properties
            s_v = if97.s_r2(p, T)    # Vapor Specific enthalpy                     (J/kg)
      
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties
            
            # Both saturated liquid and saturated vapor properties is defined by 
            # region 3, itteration in p-s space is done with p_3sat_s(s) which is
            # the saturation line p-s between region 3 and 4. 
            
            int_tol = 75
            s_v = aux.hybrid_bisec_secant(if97_tvphps.p_3sat_s, global_property.s_crit, global_property.s_r3_upper, p, int_tol)
                
        else:
            print(f'Pressure, p = {p*10**-6:.3f} MPa, above critical pressure')
            s_v = np.nan
                        
        return s_v
    
    @staticmethod
    def vL_p(p):
        """
        Saturated liquid specific volume for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  specific volume (m^3/kg)

        """
        # Init
        v_l = 0.0
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            v_l = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 1  to calculate liquid saturation properties
            v_l = if97.v_r1(p, T)    # Liquid Specific volume (m^3/kg)
      
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_l = iapwsif97.hL_p(p)
            
            # Determine v_l using v_ph function for region 3
            v_l = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_l),h_l) # use p from p_3sat_h to ensure region 3 otherwise error  
        else:
            print(f'Pressure, p = {p*10**-6:.3f} MPa, above critical pressure')
            v_l = np.nan
            
        return v_l
    
    @staticmethod
    def vV_p(p):
        """
        Saturated vapor specific volume for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  specific volume (m^3/kg)

        """
        v_v = 0.0 # (J/kg)
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            v_v = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            
            # Use region 2  to calculate vapor saturation properties
            v_v = if97.v_r2(p, T)    # Vapor Specific volume (m^3/kg)
            
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_v = iapwsif97.hV_p(p)
            
            # Determine v_l using v_ph function for region 3
            v_v = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_v),h_v) # use p from p_3sat_h to ensure region 3 otherwise error
            
        else:
            print(f'Pressure, p = {p*10**-6:.3f} MPa, above critical pressure')
            v_v = np.nan
            
        return v_v
    
    @staticmethod
    def uL_p(p):
        """
        Saturated liquid internal energy for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  Internal energy (J/kg)

        """
        # Init
        u_l = 0.0
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            u_l = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 1  to calculate liquid saturation properties
            u_l = if97.u_r1(p, T)    # Liquid internal energy                    (J/kg)
      
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_l = iapwsif97.hL_p(p)
            
            # Determine v_l using v_ph function for region 3, v_l is used 
            # for region 3 inout (p, rho)
            v_l = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_l),h_l) # use p from p_3sat_h to ensure region 3 otherwise error
            u_l = if97.u_r3(1.0/v_l, T)
            
        else:
            print(f'Pressure, p = {p*10**-6:.3f} MPa, above critical pressure')
            u_l = np.nan
            
        return u_l
    
    @staticmethod
    def uV_p(p):
        """
        Saturated vapor internal energy for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  Internal energy (J/kg)

        """
        u_v = 0.0 # (J/kg)
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            u_v = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 2  to calculate vapor saturation properties
            u_v = if97.u_r2(p, T)    # Vapor internal energy   (J/kg)
            
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_v = iapwsif97.hV_p(p)
            
            # Determine v_l using v_ph function for region 3
            v_v = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_v),h_v) # use p from p_3sat_h to ensure region 3 otherwise error
            u_v = if97.u_r3(1.0/v_v, T)
            
        else:
            print(f'Pressure, p = {p*10**-6:.3f} MPa, above critical pressure')
            u_v = np.nan
                        
        return u_v
    
    @staticmethod
    def cpL_p(p):
        """
        Saturated liquid isobaric heat capacity for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  isobaric heat capacity (J/kg/K)

        """
        # Init
        cp_l = 0.0
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            cp_l = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 1  to calculate liquid saturation properties
            cp_l = if97.cp_r1(p, T)    # Liquid isobaric heat capacity (J/kg/K)
      
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_l = iapwsif97.hL_p(p)
            
            # Determine v_l using v_ph function for region 3, v_l is used 
            # for region 3 inout (p, rho)
            v_l = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_l),h_l) # use p from p_3sat_h to ensure region 3 otherwise error
            cp_l = if97.cp_r3(1.0/v_l, T)
            
        else:
            print(f'Pressure, p = {p*10**-6:.3f} MPa, above critical pressure')
            cp_l = np.nan
            
        return cp_l
    
    @staticmethod
    def cpV_p(p):
        """
        Saturated vapor isobaric heat capacity for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  isobaric heat capacity (J/kg/K)

        """
        cp_v = 0.0 # (J/kg)
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            cp_v = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            
            # Use region 2  to calculate vapor saturation properties
            cp_v = if97.cp_r2(p, T)    # Vapor isobaric heat capacity (J/kg/K)
            
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_v = iapwsif97.hV_p(p)
            
            # Determine v_l using v_ph function for region 3
            v_v = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_v),h_v) # use p from p_3sat_h to ensure region 3 otherwise error
            cp_v = if97.cp_r3(1.0/v_v, T)
            
        else:
            print(f'Pressure, p = {p*10**-6:.3f} MPa, above critical pressure')
            cp_v = np.nan
                        
        return cp_v
    
    @staticmethod
    def cvL_p(p):
        """
        Saturated liquid isochoric heat capacity for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  isochoric heat capacity (J/kg/K)

        """
        # Init
        cv_l = 0.0
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            cv_l = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 1  to calculate liquid saturation properties
            cv_l = if97.cv_r1(p, T)    # Liquid isochoric heat capacity (J/kg/K)
      
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_l = iapwsif97.hL_p(p)
            
            # Determine v_l using v_ph function for region 3, v_l is used 
            # for region 3 inout (p, rho)
            v_l = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_l),h_l) # use p from p_3sat_h to ensure region 3 otherwise error
            cv_l = if97.cv_r3(1.0/v_l, T)
            
        else:
            print(f'Pressure, p = {p*10**-6:.3f} MPa, above critical pressure')
            cv_l = np.nan
            
        return cv_l
    
    @staticmethod
    def cvV_p(p):
        """
        Saturated vapor isochoric heat capacity for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  isochoric heat capacity (J/kg/K)

        """
        cv_v = 0.0 # (J/kg)
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            cv_v = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            
            # Use region 2  to calculate vapor saturation properties
            cv_v = if97.cv_r2(p, T)    # Vapor isochoric heat capacity (J/kg/K)
            
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_v = iapwsif97.hV_p(p)
            
            # Determine v_l using v_ph function for region 3
            v_v = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_v),h_v) # use p from p_3sat_h to ensure region 3 otherwise error
            cv_v = if97.cv_r3(1.0/v_v, T)
            
        else:
            print(f'Pressure, p = {p*10**-6:.3f} MPa, above critical pressure')
            cv_v = np.nan
                        
        return cv_v
    
    @staticmethod
    def wL_p(p):
        """
        Saturated liquid speed of sound for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  speed of sound (m/s)

        """
        # Init
        w_l = 0.0
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            w_l = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 1  to calculate liquid saturation properties
            w_l = if97.w_r1(p, T)    # speed of sound (m/s)
      
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_l = iapwsif97.hL_p(p)
            
            # Determine v_l using v_ph function for region 3, v_l is used 
            # for region 3 inout (p, rho)
            v_l = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_l),h_l) # use p from p_3sat_h to ensure region 3 otherwise error
            w_l = if97.w_r3(1.0/v_l, T)
            
        else:
            print(f'Pressure, p = {p*10**-6:.3f} MPa, above critical pressure')
            w_l = np.nan
            
        return w_l
    
    @staticmethod
    def wV_p(p):
        """
        Saturated vapor speed of sound for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  speed of sound (m/s)

        """
        w_v = 0.0 # (J/kg)
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            w_v = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:

            # Use region 2  to calculate vapor saturation properties
            w_v = if97.w_r2(p, T)    # speed of sound (m/s)
            
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_v = iapwsif97.hV_p(p)
            
            # Determine v_l using v_ph function for region 3
            v_v = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_v),h_v) # use p from p_3sat_h to ensure region 3 otherwise error
            w_v = if97.w_r3(1.0/v_v, T)
            
        else:
            print(f'Pressure, p = {p*10**-6:.3f} MPa, above critical pressure')
            w_v = np.nan
                        
        return w_v
    
    @staticmethod
    def rhoL_p(p):
        """
        Saturated liquid density for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  Density (kg/m^3)

        """
        return 1/iapwsif97.vL_p(p)
    
    @staticmethod
    def rhoV_p(p):
        """
        Saturated vapor density for pressure p

        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  Density (kg/m^3)

        """
        return 1/iapwsif97.vV_p(p)
    
    @staticmethod
    def tcL_p(p):
        """
        Saturated liquid thermal conductivity IAPWSF11
    
        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  thermal conductivity (W/m/K)
        
        """
        tc_l = 0.0
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            tc_l = np.nan
        # Check if within valid region for tc_pt
        elif f11.region_tc_pt(p, T):

            rho = iapwsif97.rhoL_p(p)
            cp = iapwsif97.cpL_p(p)
            cv = iapwsif97.cvL_p(p)
            my = f08.my_trho(T, rho) # calculate my without critical enhancement as is done in iapwsif11 verifications

            # determine derivative drho/dp at constant temperature
            dp = 100 # Pa, discrete length of derivative
            # one sided derivative to make sure density in liquid-state
            drhodpT = (iapwsif97.rho_pt(p+dp+1,T)-iapwsif97.rho_pt(p+1,T))/dp
                
            # drho/dp at reference temperature
            drhodpTr = f11.sigma_tr(rho)
                
            tc_l = f11.lambda_tc(T, rho, cp, cv, my, drhodpT, drhodpTr)
                    
        else:
            print('tcL_p, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T):.3f} K')
            tc_l = np.nan
            
        return tc_l
    
    @staticmethod
    def tcV_p(p):
        """
        Saturated vapor thermal conductivity IAPWSF11
    
        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  thermal conductivity (W/m/K)
        
        """
        
        tc_v = 0.0
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            tc_v = np.nan
        # Check if within valid region for tc_pt
        elif f11.region_tc_pt(p, T):

            rho = iapwsif97.rhoV_p(p)
            cp = iapwsif97.cpV_p(p)
            cv = iapwsif97.cvV_p(p)
            my = f08.my_trho(T, rho) # calculate my without critical enhancement as is done in iapwsif11 verifications

            # determine derivative drho/dp at constant temperature
            dp = 100 # Pa, discrete length of derivative
            # one sided derivative to make sure density in liquid-state
            drhodpT = (iapwsif97.rho_pt(p-1,T)-iapwsif97.rho_pt(p-dp-1,T))/dp
                
            # drho/dp at reference temperature
            drhodpTr = f11.sigma_tr(rho)
                
            tc_v = f11.lambda_tc(T, rho, cp, cv, my, drhodpT, drhodpTr)
            
        else:
            print(f'tcV_p, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T):.3f} K')
            tc_v = np.nan
            
        return tc_v
    
    @staticmethod
    def myL_p(p):
        """
        Saturated liquid Viscosity IAPWSF08
    
        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  viscosity (Pa*s)
        """
        my_l = 0.0
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            my_l = np.nan
        # Check if within valid region for my_pt
        elif f08.region_my_pt(p, T):

            rho = iapwsif97.rhoL_p(p)
            my_l = f08.my_trho(T, rho) # calculate my without critical enhancement as suggested using iapwsif97
        else:
            print(f'myL_p, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T):.3f} K')
            my_l = np.nan
            
        return my_l
    
    @staticmethod
    def myV_p(p):
        """
        Saturated vapor Viscosity IAPWSF08
    
        Parameters
        ----------
        p:      double  pressure (Pa).

        Returns
        -------
        double  viscosity (Pa*s)
        """
        
        my_v = 0.0
        
        # Saturation temperature
        T = iapwsif97.Tsat_p(p)
        
        if np.isnan(T): # Valid region checked inside tsat_p...
            my_v = np.nan
        # Check if within valid region for my_pt
        elif f08.region_my_pt(p, T):

            rho = iapwsif97.rhoV_p(p)
            my_v = f08.my_trho(T, rho) # calculate my without critical enhancement as suggested using iapwsif97
        else:
            print(f'myV_p, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T):.3f} K')
            my_v = np.nan
            
        return my_v
    
    ### Saturation properties, Liquid and Vapor ###
    ### _temperature ###
    
    @staticmethod
    def hL_T(T):
        """
        Saturated liquid enthalpy for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  Enthalpy (J/kg)

        """
        h_l = 0.0 # (J/kg)
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)
        
        if np.isnan(p): # Valid region checked inside p_sat_t...
            h_l = np.nan
        # Check if temperature is below region 3 boundary
        elif T < global_property.T_boundary_13:
            # Use region 1  to calculate liquid saturation properties
            h_l = if97.h_r1(p, T)    # Liquid Specific enthalpy                     (J/kg)
      
        # temperature is above region 3 lower boundary, check if its below critical temperature
        elif T < global_property.T_crit:
            # Use region 3 for saturation properties
            
            # Both saturated liquid and saturated vapor properties is defined by 
            # region 3, itteration in p-h space is done with p_3sat_h(h) which is
            # the saturation line p-h between region 3 and 4. 

            int_tol = 1000
            h_l = aux.hybrid_bisec_secant(if97_tvphps.p_3sat_h, global_property.h_r3_lower , global_property.h_crit , p, int_tol)
        else:
            print(f'Temperature, T = {T:.3f} K, above critical temperature')
            h_l = np.nan
        
        return h_l
    
    @staticmethod
    def hV_T(T):
        """
        Saturated vapor enthalpy for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  Enthalpy (J/kg)

        """
        h_v = 0.0 # (J/kg)
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)
        
        
        if np.isnan(p): # Valid region checked inside p_sat_t...
            h_v = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 2  to calculate vapor saturation properties
            h_v = if97.h_r2(p, T)    # Vapor Specific enthalpy                     (J/kg)
            
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties
            
            # Both saturated liquid and saturated vapor properties is defined by 
            # region 3, itteration in p-h space is done with p_3sat_h(h) which is
            # the saturation line p-h between region 3 and 4. 

            int_tol = 1000
            h_v = aux.hybrid_bisec_secant(if97_tvphps.p_3sat_h, global_property.h_crit, global_property.h_r3_upper , p, int_tol)
        else:
            print(f'Temperature, T = {T:.3f} K, above critical temperature')
            h_v = np.nan
        
        return h_v
    
    @staticmethod
    def sL_T(T):
        """
        Saturated liquid entropy for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  Entropy (J/kg/K)

        """
        s_l = 0.0 # (J/kg)
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)
        
        if np.isnan(p): # Valid region checked inside p_sat_t...
            s_l = np.nan
        # Check if temperature is below region 3 boundary
        elif T < global_property.T_boundary_13:
            # Use region 1  to calculate liquid saturation properties
            s_l = if97.s_r1(p, T)    # Liquid Specific entropy         (J/kg/K)
      
        # temperature is above region 3 lower boundary, check if its below critical temperature
        elif T < global_property.T_crit:
            # Use region 3 for saturation properties
            
            # Both saturated liquid and saturated vapor properties is defined by 
            # region 3, itteration in p-s space is done with p_3sat_s(s) which is
            # the saturation line p-s between region 3 and 4. 
             
            int_tol = 75
            s_l = aux.hybrid_bisec_secant(if97_tvphps.p_3sat_s, global_property.s_r3_lower, global_property.s_crit, p, int_tol)
        else:
            print(f'Temperature, T = {T:.3f} K, above critical temperature')
            s_l = np.nan
        
              
        return s_l
    
    @staticmethod
    def sV_T(T):
        """
        Saturated vapor entropy for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  Entropy (J/kg/K)

        """
        s_v = 0.0 # (J/kg)
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)
        
        if np.isnan(p): # Valid region checked inside p_sat_t...
            s_v = np.nan
        # Check if pressure is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 2  to calculate vapor saturation properties
            s_v = if97.s_r2(p, T)    # Vapor Specific entropy        (J/kg/K)
            
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties
            
            # Both saturated liquid and saturated vapor properties is defined by 
            # region 3, itteration in p-s space is done with p_3sat_s(s) which is
            # the saturation line p-s between region 3 and 4. 
            
            int_tol = 75
            s_v = aux.hybrid_bisec_secant(if97_tvphps.p_3sat_s, global_property.s_crit, global_property.s_r3_upper, p, int_tol)
            
        else:
            print(f'Temperature, T = {T:.3f} K, above critical temperature')
            s_v = np.nan
                        
        return s_v
    
    @staticmethod
    def vL_T(T):
        """
        Saturated liquid specific volume for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  specific volume (m^3/kg)

        """
        v_l = 0.0 # (m^3/kg)
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)

        if np.isnan(p): # Valid region checked inside p_sat_t...
            v_l = np.nan    
        # Check if temperature is below region 3 boundary
        elif T < global_property.T_boundary_13:
            # Use region 1  to calculate liquid saturation properties
            v_l = if97.v_r1(p, T)    # Liquid Specific volume          (m^3/kg)
      
        # temperature is above region 3 lower boundary, check if its below critical temperature
        elif T < global_property.T_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_l = iapwsif97.hL_T(T)            

            # Determine v_l using v_ph function for region 3
            v_l = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_l),h_l) # use p from p_3sat_h to ensure region 3 otherwise error
        
        else:
            print(f'Temperature, T = {T:.3f} K, above critical temperature')
            v_l = np.nan
                 
        return v_l
    
    @staticmethod
    def vV_T(T):
        """
        Saturated vapor specific volume for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  specific volume (m^3/kg)

        """
        v_v = 0.0 # (J/kg)
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)
        
        if np.isnan(p): # Valid region checked inside p_sat_t...
            v_v = np.nan    
        # Check if temperature is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 2  to calculate vapor saturation properties
            v_v = if97.v_r2(p, T)    # Vapor Specific volume          (m^3/kg)
            
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_v = iapwsif97.hV_T(T)            

            # Determine v_l using v_ph function for region 3
            v_v = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_v),h_v) # use p from p_3sat_h to ensure region 3 otherwise error
            
        else:
            print(f'Temperature, T = {T:.3f} K, above critical temperature')
            v_v = np.nan
            
        return v_v
    
    @staticmethod
    def uL_T(T):
        """
        Saturated liquid internal energy for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  internal energy (kJ/kg)

        """
        u_l = 0.0 # (kJ/kg)
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)

        if np.isnan(p): # Valid region checked inside p_sat_t...
            u_l = np.nan    
        # Check if temperature is below region 3 boundary
        elif T < global_property.T_boundary_13:
            # Use region 1  to calculate liquid saturation properties
            u_l = if97.u_r1(p, T)    # Liquid internal energy          (kJ/kg)
      
        # temperature is above region 3 lower boundary, check if its below critical temperature
        elif T < global_property.T_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_l = iapwsif97.hL_T(T)            

            # Determine v_l using v_ph function for region 3
            v_l = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_l),h_l) # use p from p_3sat_h to ensure region 3 otherwise error
            u_l = if97.u_r3(1.0/v_l, T)
        else:
            print(f'Temperature, T = {T:.3f} K, above critical temperature')
            u_l = np.nan
                 
        return u_l
    
    @staticmethod
    def uV_T(T):
        """
        Saturated vapor internal energy for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  internal energy (kJ/kg)

        """
        u_v = 0.0 # (J/kg)
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)
        
        if np.isnan(p): # Valid region checked inside p_sat_t...
            u_v = np.nan    
        # Check if temperature is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 2  to calculate vapor saturation properties
            u_v = if97.u_r2(p, T)    # Vapor internal energy          (kJ/kg)
            
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_v = iapwsif97.hV_T(T)            

            # Determine v_l using v_ph function for region 3
            v_v = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_v),h_v) # use p from p_3sat_h to ensure region 3 otherwise error
            u_v = if97.u_r3(1.0/v_v, T)
            
        else:
            print(f'Temperature, T = {T:.3f} K, above critical temperature')
            u_v = np.nan
            
        return u_v
    
    @staticmethod
    def cpL_T(T):
        """
        Saturated liquid isobaric heat capacity for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  isobaric heat capacity (J/kg/K)

        """
        cp_l = 0.0 # (kJ/kg)
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)

        if np.isnan(p): # Valid region checked inside p_sat_t...
            cp_l = np.nan    
        # Check if temperature is below region 3 boundary
        elif T < global_property.T_boundary_13:
            # Use region 1  to calculate liquid saturation properties
            cp_l = if97.cp_r1(p, T)    # Liquid isobaric heat capacity  (J/kg/K)
      
        # temperature is above region 3 lower boundary, check if its below critical temperature
        elif T < global_property.T_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_l = iapwsif97.hL_T(T)            

            # Determine v_l using v_ph function for region 3
            v_l = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_l),h_l) # use p from p_3sat_h to ensure region 3 otherwise error
            cp_l = if97.cp_r3(1.0/v_l, T)
        else:
            print(f'Temperature, T = {T:.3f} K, above critical temperature')
            cp_l = np.nan
                 
        return cp_l
    
    @staticmethod
    def cpV_T(T):
        """
        Saturated vapor isobaric heat capacity for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  isobaric heat capacity (J/kg/K)

        """
        cp_v = 0.0 # (J/kg)
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)
        
        if np.isnan(p): # Valid region checked inside p_sat_t...
            cp_v = np.nan    
        # Check if temperature is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 2  to calculate vapor saturation properties
            cp_v = if97.cp_r2(p, T)    # Vapor isobaric heat capacity   (J/kg/K)
            
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_v = iapwsif97.hV_T(T)            

            # Determine v_l using v_ph function for region 3
            v_v = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_v),h_v) # use p from p_3sat_h to ensure region 3 otherwise error
            cp_v = if97.cp_r3(1.0/v_v, T)
            
        else:
            print(f'Temperature, T = {T:.3f} K, above critical temperature')
            cp_v = np.nan
            
        return cp_v
    
    @staticmethod
    def cvL_T(T):
        """
        Saturated liquid isochoric heat capacity for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  isochoric heat capacity (J/kg/K)

        """
        cv_l = 0.0 # (kJ/kg)
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)

        if np.isnan(p): # Valid region checked inside p_sat_t...
            cv_l = np.nan    
        # Check if temperature is below region 3 boundary
        elif T < global_property.T_boundary_13:
            # Use region 1  to calculate liquid saturation properties
            cv_l = if97.cv_r1(p, T)    # Liquid isochoric heat capacity  (J/kg/K)
      
        # temperature is above region 3 lower boundary, check if its below critical temperature
        elif T < global_property.T_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_l = iapwsif97.hL_T(T)            

            # Determine v_l using v_ph function for region 3
            v_l = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_l),h_l) # use p from p_3sat_h to ensure region 3 otherwise error
            cv_l = if97.cv_r3(1.0/v_l, T)
        else:
            print(f'Temperature, T = {T:.3f} K, above critical temperature')
            cv_l = np.nan
                 
        return cv_l
    
    @staticmethod
    def cvV_T(T):
        """
        Saturated vapor isochoric heat capacity for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  isochoric heat capacity (J/kg/K)

        """
        cv_v = 0.0 # (J/kg)
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)
        
        if np.isnan(p): # Valid region checked inside p_sat_t...
            cv_v = np.nan    
        # Check if temperature is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 2  to calculate vapor saturation properties
            cv_v = if97.cv_r2(p, T)    # Vapor isobaric heat capacity   (J/kg/K)
            
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_v = iapwsif97.hV_T(T)            

            # Determine v_l using v_ph function for region 3
            v_v = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_v),h_v) # use p from p_3sat_h to ensure region 3 otherwise error
            cv_v = if97.cv_r3(1.0/v_v, T)
            
        else:
            print(f'Temperature, T = {T:.3f} K, above critical temperature')
            cv_v = np.nan
            
        return cv_v
    
    @staticmethod
    def wL_T(T):
        """
        Saturated liquid speed of sound for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  speed of sound (m/s)

        """
        w_l = 0.0 # (kJ/kg)
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)

        if np.isnan(p): # Valid region checked inside p_sat_t...
            w_l = np.nan    
        # Check if temperature is below region 3 boundary
        elif T < global_property.T_boundary_13:
            # Use region 1  to calculate liquid saturation properties
            w_l = if97.w_r1(p, T)    # Liquid speed of sound  (m/s)
      
        # temperature is above region 3 lower boundary, check if its below critical temperature
        elif T < global_property.T_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_l = iapwsif97.hL_T(T)            

            # Determine v_l using v_ph function for region 3
            v_l = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_l),h_l) # use p from p_3sat_h to ensure region 3 otherwise error
            w_l = if97.w_r3(1.0/v_l, T)
        else:
            print(f'Temperature, T = {T:.3f} K, above critical temperature')
            w_l = np.nan
                 
        return w_l
    
    @staticmethod
    def wV_T(T):
        """
        Saturated vapor speed of sound for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  speed of sound (m/s)

        """
        w_v = 0.0 # (J/kg)
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)
        
        if np.isnan(p): # Valid region checked inside p_sat_t...
            w_v = np.nan    
        # Check if temperature is below region 3 boundary
        elif p < global_property.p_r3_lower_boundary:
            # Use region 2  to calculate vapor saturation properties
            w_v = if97.w_r2(p, T)    # Vapor isobaric heat capacity   (J/kg/K)
            
        # pressure is above region 3 lower boundary, check if its below critical pressure
        elif p < global_property.p_crit:
            # Use region 3 for saturation properties, done by
            # reusing existing hL-function
            
            h_v = iapwsif97.hV_T(T)            

            # Determine v_l using v_ph function for region 3
            v_v = if97_tvphps.v_3_ph(if97_tvphps.p_3sat_h(h_v),h_v) # use p from p_3sat_h to ensure region 3 otherwise error
            w_v = if97.w_r3(1.0/v_v, T)
            
        else:
            print(f'Temperature, T = {T:.3f} K, above critical temperature')
            w_v = np.nan
            
        return w_v

    @staticmethod
    def rhoL_T(T):
        """
        Saturated liquid density for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  Density (kg/m^3)

        """
        return 1/iapwsif97.vL_T(T)
    
    @staticmethod
    def rhoV_T(T):
        """
        Saturated vapor density for temperature T

        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  Density (kg/m^3)

        """
        return 1/iapwsif97.vV_T(T)
    
    @staticmethod
    def tcL_T(T):
        """
        Saturated liquid thermal conductivity IAPWSF11
    
        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  thermal conductivity (W/m/K)
        
        """
        
        tc_l = 0.0
          
        # Saturation pressure
        p = iapwsif97.psat_t(T)
        
        if np.isnan(p): # Valid region checked inside p_sat_t...
            tc_l = np.nan    
        #Check if valid region for tc
        elif f11.region_tc_pt(p, T):

            rho = iapwsif97.rhoL_p(p)
            cp = iapwsif97.cpL_p(p)
            cv = iapwsif97.cvL_p(p)
            my = f08.my_trho(T, rho) # calculate my without critical enhancement as is done in iapwsif11 verifications

            # determine derivative drho/dp at constant temperature
            dp = 100 # Pa, discrete length of derivative
            # one sided derivative to make sure density in liquid-state
            drhodpT = (iapwsif97.rho_pt(p+dp+1,T)-iapwsif97.rho_pt(p+1,T))/dp
                
            # drho/dp at reference temperature
            drhodpTr = f11.sigma_tr(rho)
                
            tc_l = f11.lambda_tc(T, rho, cp, cv, my, drhodpT, drhodpTr)
        else:
            print(f'myV_p, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T):.3f} K')
            tc_l = np.nan
            
        return tc_l
    
    @staticmethod
    def tcV_T(T):
        """
        Saturated vapor thermal conductivity IAPWSF11
    
        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  thermal conductivity (W/m/K)
        
        """
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)
        
        tc_v = 0.0
        
        if np.isnan(p): # Valid region for T checked inside p_sat_t...
            tc_v = np.nan    
        # Check if valid region for tc 
        elif f11.region_tc_pt(p, T):

            rho = iapwsif97.rhoV_p(p)
            cp = iapwsif97.cpV_p(p)
            cv = iapwsif97.cvV_p(p)
            my = f08.my_trho(T, rho) # calculate my without critical enhancement as is done in iapwsif11 verifications

            # determine derivative drho/dp at constant temperature
            dp = 100 # Pa, discrete length of derivative
            # one sided derivative to make sure density in liquid-state
            drhodpT = (iapwsif97.rho_pt(p-1,T)-iapwsif97.rho_pt(p-dp-1,T))/dp
                
            # drho/dp at reference temperature
            drhodpTr = f11.sigma_tr(rho)
                
            tc_v = f11.lambda_tc(T, rho, cp, cv, my, drhodpT, drhodpTr)   
        else:
            print(f'myV_p, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T):.3f} K')
            tc_v = np.nan
            
        return tc_v
    
    @staticmethod
    def myL_T(T):
        """
        Saturated liquid Viscosity IAPWSF08
    
        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  viscosity (Pa*s)
        """
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)
        
        my_l = 0.0
        
        if np.isnan(p): # Valid region for T checked inside p_sat_t...
            my_l = np.nan    
        # Check if valid region for tc 
        elif f08.region_my_pt(p, T):

            rho = iapwsif97.rhoL_p(p)
            my_l = f08.my_trho(T, rho) # calculate my without critical enhancement as suggested using iapwsif97
        else:
            print(f'myV_p, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T):.3f} K')
            my_l = np.nan    
            
        return my_l
    
    @staticmethod
    def myV_T(T):
        """
        Saturated vapor Viscosity IAPWSF08
    
        Parameters
        ----------
        T:      double  temperature (K).

        Returns
        -------
        double  viscosity (Pa*s)
        """
        
        # Saturation pressure
        p = iapwsif97.psat_t(T)
        
        my_v = 0.0
          
        if np.isnan(p): # Valid region for T checked inside p_sat_t...
            my_v = np.nan    
        # Check if valid region for tc 
        elif f08.region_my_pt(p, T):

            rho = iapwsif97.rhoV_p(p)
            my_v = f08.my_trho(T, rho) # calculate my without critical enhancement as suggested using iapwsif97
        else:
            print(f'myV_t, out of bounds. p = {p*10**-6:.3f} MPa, T = {(T):.3f} K')
            my_v = np.nan
            
        return my_v

### MISC ###
    @staticmethod
    def st_t(T):
        """
        Surface tension as a function of temperature IAPWS R1-76
        Between liquid and vapor phases of water

        Parameters
        ----------
        t : double  temperature (K).
        
        Returns
        -------
        double  surface tension (N/m)
        
        """   
        st_t = 0.0
        
        if r176.region_st_t(T):
            st_t = r176.st_t(T)
        else:
            st_t = np.nan
            
        return st_t
    
    @staticmethod
    def st_p(p):
        """
        Surface tension saturated liquid IAPWS R1-76
        Between liquid and vapor phases of water

        Parameters
        ----------
        p:      double  pressure (Pa).
        
        Returns
        -------
        double  surface tension (N/m)
        
        """   
        st_t = 0.0
        
        
        T = iapwsif97.Tsat_p(p)
        
        if r176.region_st_t(T):
            st_t = r176.st_t(T)
        else:
            st_t = np.nan
            
        return st_t

if __name__ == "__main__":
    main()
