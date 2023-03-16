# -*- coding: utf-8 -*-
"""
IAPWSIF97 Itterative methods
Backwards equations, for rho_pt (Region 3), t_ph and t_ps using itteration for higher accuracy

List of functions and code structure:
    rho_r3(p,T, rho_init = 0.0)
    t_itt(func, p,hs, T_init, initial_guess_accuracy = .025)
    t_ph_r1_itt(p,h)
    t_ph_r2_itt(p,h)
    t_ph_r3_itt(p,h)
    t_ps_r1_itt(p,h)
    t_ps_r2_itt(p,h)
    t_ps_r3_itt(p,h)
    
t_p*_r3_itt not implemented yet, just a dummy function.

@author: Christoffer Rappmann, christoffer.rappmann@gmail.com
"""

# sub functions
#
# Helper functions
from .... import utils as aux
# Global Constants
from .iapwsif97_globals import iapwsif_globals as global_property
# t_ph_r* and t_ps_r*
from .iapwsif97_main import t_ph_r1, t_ph_r2, t_ps_r1, t_ps_r2, h_r1, h_r2, s_r1, s_r2 , p_r3
from .iapwsif97_tps3_tph3_vph3_vps3 import t_3_ph, t_3_ps



def rho_r3(p,T, rho_init = 0.0, initial_guess_accuracy = .025):
    """
    Density for region 3 IAPWS-IF97
    Itteration...
    
    Parameters
    ----------
    p:            double              pressure (Pa).
    T:            double              temperature (K).
    rho_init:     double              initial guess density (kg/m^3) (optional).
    
    Returns
    -------
    double  density (kg/m^3)

    """
    
    # check if initial guess is provided...
    if rho_init == 0.0:
        # set a range bigger than the range of possible densities in region 3...
        rho_0 = 1.0 # kg/m^3
        rho_1 = 1000.0 # kg/m^3
    else:
        # set a range +-accuracy around initial guess
        rho_0 = (1-initial_guess_accuracy)*rho_init
        rho_1 = (1+initial_guess_accuracy)*rho_init
        
    p_rho = lambda rho: p_r3(rho, T)
    rho = aux.bisection(p_rho, rho_0, rho_1, p, global_property.itt_tol_r3, global_property.itt_max_r3)
       
    return rho

def t_itt(func, p,hs, T_init, initial_guess_accuracy = .025):
    """

    Itteration for t_ph and t_ps...
    
    Parameters
    ----------
    func:       function        function to itterate over, h_pt_r* * = 1-5
    p:          double          pressure (Pa).
    hs:          double         Enthalpy (J/kg) or entropy
    T_init:     double          initial guess temperature (K).
    initial_guess_accuracy double accuracy of initial guess
    
    Returns
    -------
    double  density (kg/m^3)

    """
    
    # set a range +-accuracy around initial guess
    T_0 = (1-initial_guess_accuracy)*T_init
    T_1 = (1+initial_guess_accuracy)*T_init
        
    hs_T = lambda T: func(p, T)
    T = aux.bisection(hs_T, T_0, T_1, hs, global_property.itt_tol, global_property.itt_max)
       
    return T

def t_ph_r1_itt(p,h):
    """
    Temperature from p,h region 1 IAPWS-IF97
    Itterations...
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    h:      double  enthalpy (J/kg).

    Returns
    -------
    double  Temperature (K)

    """
    
    # Initial guess from backward equation
    T_init = t_ph_r1(p,h)
    
    return t_itt(h_r1, p, h, T_init)

def t_ph_r2_itt(p,h):
    """
    Temperature from p,h region 2 IAPWS-IF97
    Itterations...
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    h:      double  enthalpy (J/kg).

    Returns
    -------
    double  Temperature (K)

    """
    
    # Initial guess from backward equation
    T_init = t_ph_r2(p,h)
    
    return t_itt(h_r2, p, h, T_init)

def t_ph_r3_itt(p,h):
    """
    Temperature from p,h region 3 IAPWS-IF97
    Itterations...
    
    Parameters
    ----------
    rho:      double  density (kg/m3).
    h:      double  enthalpy (J/kg).

    Returns
    -------
    double  Temperature (K)

    """
    
    # To be implemented...
    
    # Initial guess from backward equation
    T_init = t_3_ph(p,h)
    
    return T_init

def t_ps_r1_itt(p,s):
    """
    Temperature from p,s region 1 IAPWS-IF97
    Itterations...
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    s:      double  entropy (J/kg/K).

    Returns
    -------
    double  Temperature (K)

    """
    
    # Initial guess from backward equation
    T_init = t_ps_r1(p,s)
    
    return t_itt(s_r1, p, s, T_init)

def t_ps_r2_itt(p,s):
    """
    Temperature from p,s region 2 IAPWS-IF97
    Itterations...
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    s:      double  entropy (J/kg/K).

    Returns
    -------
    double  Temperature (K)

    """
    
    # Initial guess from backward equation
    T_init = t_ps_r2(p,s)
    
    return t_itt(s_r2, p, s, T_init)

def t_ps_r3_itt(p,s):
    """
    Temperature from p,s region 3 IAPWS-IF97
    Itterations...
    
    Parameters
    ----------
    rho:      double  density (kg/m3).
    h:      double  entropy (J/kg/K).

    Returns
    -------
    double  Temperature (K)

    """
    
    # To be implemented...
    
    # Initial guess from backward equation
    T_init = t_3_ps(p,s)
    
    return T_init


    