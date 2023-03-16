# -*- coding: utf-8 -*-
"""
Module containing functions implemented from equations presentend in the 
SR2-01(2014) and SR4-04(2014) suplementary release of IAPWSIF97.

SR2-01 (Region 1 and 2)
SR4-04 (Region 3 and 4 and region boundaries)

List of functions and code structure:
    Determine region
        find_region_hs(h,s)
    Region boundaries
        h_max_s(s)
        h_max_s_r5(s)
        h_prim_s_r1(s)
        h_prim_s_r3(s)
        h_s_b13(s)
        h_bis_s_r2(s)
        h_bis_s_r23(s)
        T_hs_b23(h,s)
    Region 1
        p_hs_r1(h,s)
        h_b2ab(s)
    Region 2
        p_hs_r2a(h,s)
        p_hs_r2b(h,s)
        p_hs_r2c(h,s)
        p_hs_r2(h,s)
    Region 3
        p_hs_r3a(h,s)
        p_hs_r3b(h,s)
        p_hs_r3(h,s)
    Region 4
        T_sat_hs(h,s)

@author: Christoffer Rappmann, christoffer.rappmann@gmail.com
"""

from .iapwsif97_globals import iapwsif_globals as global_property
from .iapwsif97_main import p_t_b23
from .... import utils as aux
import numpy as np


###########################
###  Determine Region   ###
###########################

def find_region_hs(h,s):
    """
    The IAPWS Industrial Formulation 1997 consists of a set of equations for different regions
    which cover the following range of validity:
    273.15 K  <= T <= 1073.15 K, p <= 100 MPa if97.s_r1(p*10**6,T_boundary_1)*10**-3
    1073.15 K < T  <= 2273.15 K, p <= 50 MPa, region 5, not implemented
    
    Which corresponds 
    -0.0002 s'(273.15) < s <= 11.921 (s_pt(0.000611*10**6,1073.15)) kJ/kg/K
    -0.042 h_prim_s_r1(s'(273.15)) < h <= 4160.82 (h_pt(100*10**6,1073.15) kJ/kg
        
    This method determine the applicable region given conditions h,s
    the boundaries needs to be converted to enthalpy-lines for a given entropy

    Parameters
    ----------
    h : double      Enthalpy (kJ/kg).
    s : double      Entropy (kJ/kg/K).

    Returns
    -------
    region : int    Applicable region (-)
    """
    
    ### Fixed boundaries ###
    s_r1_min = global_property.s_r1_min # J/kg/K
    s_r2_max = global_property.s_r2_max # J/kg/K
    s_r4_max = global_property.s_h_bis_273 # J/kg/K
    s_r5_min = global_property.s_r5_min # J/kg/K
    s_r5_max = global_property.s_r5_max # J/kg/K
    
    # Entropy boundaries, between min-max ranges check enthalpy with b13 respectively b23 eq.
    s_boundary_13_1 = global_property.s_h_b13_min # kJ/kg/K
    s_boundary_13_2 = global_property.s_h_b13_max # kJ/kg/K
    
    s_boundary_32_1 = global_property.s_b23_min # kJ/kg/K
    s_boundary_32_2 = global_property.s_b23_max # kJ/kg/K
    
    s_crit = global_property.s_crit # kJ/kg/K
    s_585 = global_property.s_h_bis_585 # kJ/kg/K

    # init
    region = 0
        
    # Check which region should be used
    if h < h_max_s(s) and s < s_r2_max and h > h_611pa(s) and s > s_r1_min:
        # Within region 1-4
        if s < s_boundary_13_1:
            # Region 1 and 4 entropy-range
            # Use equation h_prim_s_r1
            if h>h_prim_s_r1(s):
                region = 1
            else:
                region = 4
        elif s <= s_boundary_13_2:
            # Region 1,3 and 4 entropy range
            # Use equation h_s_b13 and h_prim_s_r1
            if h>h_s_b13(s):
                region = 3
            elif h>h_prim_s_r1(s):
                region = 1
            else:
                region = 4
        elif s < s_boundary_32_1:
            # Region 3 and 4 entropy-range
            # Use eq. h_prim_s_r3
            if s <= s_crit:
                # h_prim_s_r3
                if h > h_prim_s_r3(s):
                    region = 3
                else:
                    region = 4
            else:
                # h_bis_s_r23
                if h > h_bis_s_r23(s):
                    region = 3
                else:
                    region = 4
        elif s <= s_boundary_32_2:
            # Region 2,3 and 4 entropy-range
            # Use eq. T_hs_b23, p_t_b23 and h_bis_s_r23
            # Se ch. 4.6 in SR4-04
            if h > h_bis_s_r23(s):
                T = T_hs_b23(h, s)
                if p_hs_r2c(h,s) > p_t_b23(T):
                    region = 3
                else:
                    region = 2
            else:
                region = 4
                    
        elif s <= s_r4_max:
            # Region 2 and 4 entropy range
            # Use eq. h_bis_s_r2
            if s < s_585:
                # h_bis_s_r23
                if h > h_bis_s_r23(s):
                    region = 2
                else:
                    region = 4
            else:
                # h_bis_s_r2
                if h > h_bis_s_r2(s):
                    region = 2
                else:
                    region = 4
        else:
            # Region 2
            region = 2
    elif h < h_max_s_r5(s) and s < s_r5_max and h > h_611pa(s) and s > s_r5_min:
        # Region 5
        region = 5
    else:
        # Out of bounds
        region = -1
            
    return region


###########################
### Region boundaries   ###
### SR4-04              ###
###########################

# Polynomials of region outside boundaries
# Simplifies outer boundary check, provides enough accuracy since its only a
# check that the input is inside the valid range.
# see figure 3 in SR4-04

def h_611pa(s):
    """
    Isobar for 0.000611 MPa, h(s)
    
    Divided into two polynomials to be accurate enough for region 1 and region 5 
    boundaries at the same time.
    
    Validity-range
    0.0 < s <= 11.92 kJ/kg/K
    
    Parameters
    ----------
    s:      double  Entropy (J/kg/K)
    
    Returns
    -------
    double Enthalpy (J/kg)
    """
    h_611pa_poly_r14 = [-7.21856844e+02,  
                     2.71392499e+02,  
                     7.36098774e-03, 
                    -5.38041352e-06,  
                     1.49098281e-09, 
                    -1.75257150e-13,  
                     7.35042671e-18]

    #h_611pa = lambda s: aux.gen_poly(s, h_611pa_poly)

    h_611pa_poly_r25 = [ 6.49829046e+03,
                         2.28558669e+02,
                         6.12522141e-02,
                        -3.24671405e-05,
                         8.14208042e-09,
                        -1.02356442e-12,
                         6.12915243e-17,
                        -1.35024719e-21]
    
    if s < global_property.s_h_bis_273:
        h = aux.gen_poly(s, h_611pa_poly_r14)
    else:
        h = aux.gen_poly(s, h_611pa_poly_r25)   
    
    return h

"""
 Isobar for 50 MP, h(s)
 
 Validity-range
 6.5 < s <= 8.66886 kJ/kg/K
 
 Parameters
 ----------
 s:      double  Entropy (J/kg/K)

 Returns
 -------
 double Enthalpy (J/kg)
"""   
h_50mpa_poly = [-3.68056906e+06,
                2.54104728e+03,
                -4.07265943e-01,
                3.01230582e-05]

h_50mpa = lambda s: aux.gen_poly(s, h_50mpa_poly)

"""
 Isobar for 100 MP, h(s)
 
 Validity-range
 0.0 < s <= 6.04 kJ/kg/K
 
 Parameters
 ----------
 s:      double  Entropy (J/kg/K)

 Returns
 -------
 double Enthalpy (J/kg)
"""  
""" 
h_100mpa_poly = [ 9.91420620e+04,
                     2.69024560e+02,
                     3.71511690e-02,
                     2.85428071e-06]

h_100mpa_poly = [ 9.84758019e+04,
                   2.66405274e+02,
                   5.07118530e-02,
                  -1.08623071e-05,
                   5.42184579e-09,
                  -9.26633815e-13,
                   5.73095879e-17]
"""
h_100mpa_poly = [ 9.76286067e+04,
                  2.75537849e+02,
                  2.93034662e-02,
                  9.31250120e-06,
                 -3.90438968e-09,
                  1.32085216e-12,
                 -2.13646569e-16,
                  1.29026742e-20]

h_100mpa = lambda s: aux.gen_poly(s, h_100mpa_poly)

"""
 Isotherm for 1073.15 K, h(s)
 
 Validity-range
 6.04 < s <= 8.502 kJ/kg/K
 
 Parameters
 ----------
 s:      double  Entropy (J/kg/K)

 Returns
 -------
 double Enthalpy (J/kg)
"""  
h_1073k_poly = [    -1.67604322e+07,  
                     7.55036188e+03, 
                     -9.10564053e-01,  
                     3.66822309e-05]

h_1073k = lambda s: aux.gen_poly(s, h_1073k_poly)

"""
 Isotherm for 2273.15 K, h(s)
 
 Validity-range
 8.66886 < s <= 13.905 kJ/kg/K
 
 Parameters
 ----------
 s:      double  Entropy (J/kg/K)

 Returns
 -------
 double Enthalpy (J/kg)
"""  
h_2273k_poly = [ 4.82746449e+06,
                 8.65502721e+02,
                -1.09772902e-01,
                 6.16435458e-06,
                -1.29313133e-10]

h_2273k = lambda s: aux.gen_poly(s, h_2273k_poly)


def h_max_s(s):
    """
    Combines h_100mpa, h_1073k and h_r2_max (double)
    to cover entire range (region 1-3) with one function.
 
    Validity-range
    0.0 < s <= 11.92 kJ/kg/K
 
    Parameters
    ----------
    s:      double  Entropy (J/kg/K)

    Returns
    ------
    double Enthalpy (J/kg)
    """
    
    if s <= 6.04*10**3:
        h = h_100mpa(s)
    elif s <= 8.502*10**3:
        h = h_1073k(s)
    else:
        h = global_property.h_r2_max
    return h

def h_max_s_r5(s):
    """
    Combines h_50mpa, h_2273k and h_r5_max (double)
    to cover entire range for region 5 with one function.
 
    Validity-range
    0.0 < s <= 11.92 kJ/kg/K
 
    Parameters
    ----------
    s:      double  Entropy (J/kg/K)

    Returns
    ------
    double Enthalpy (J/kg)
    """
    
    if s <= 8.66886*10**3:
        h = h_50mpa(s)
    elif s <= global_property.s_r5_max:
        h = h_2273k(s)
    else:
        h = global_property.h_r5_max
    return h


# Liquid
def h_prim_s_r1(s):
    """
    Region 1-4 boundary eq. 3 (SR4-04)
    Saturated liquid line from the triple-point temperature 273.16 K to 623.15 K 
    
    Validity range:
        s'(273.15) <= s <= s'(623.15)
        where:
        s'(273.15) = -1.545 495 919 * 10**-4 kJ/kg/K
        s'(623.15) = 3.778 281 340 kJ/kg/K (s_r3_lower)
    
    Parameters
    ----------
    s:      double  Entropy (J/kg/K)

    Returns
    -------
    double Enthalpy (J/kg)

    """
    #Reducing quantaties
    h_star = 1700 * 10**3 # J/kg
    s_star = 3.8 * 10**3 # J/kg/K
    
    #eta = h/h_star
    sigma = s/s_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 7, 8, 12, 12, 14, 14, 16, 20, 20, 22, 24, 28, 32, 32]
    J = [14, 36, 3, 16, 0, 5, 4, 36, 4, 16, 24, 18, 24, 1, 4, 2, 4, 1, 22, 10, 12, 28, 8, 3, 0, 6, 8]
    n_array = [ 0.332171191705237,
                0.611217706323496*10**-3,
               -0.882092478906822*10**1,
               -0.455628192543250,
               -0.263483840850452*10**-4,
               -0.223949661148062*10**2,
               -0.428398660164013*10**1,
               -0.616679338856916,
               -0.146823031104040*10**2,
                0.284523138727299*10**3,
               -0.113398503195444*10**3,
                0.115671380760859*10**4,
                0.395551267359325*10**3,
               -0.154891257229285*10**1,
                0.194486637751291*10**2,
               -0.357915139457043*10**1,
               -0.335369414148819*10**1,
               -0.664426796332460,
                0.323321885383934*10**5,
                0.331766744667084*10**4,
               -0.223501257931087*10**5,
                0.573953875852936*10**7,
                0.173226193407919*10**3,
               -0.363968822121321*10**-1,
                0.834596332878346*10**-6,
                0.503611916682674*10**1,
                0.655444787064505*10**2]
    
    eta_prim = 0
    for i,j,n in zip(I, J, n_array):
        eta_prim += n * (sigma - 1.09)**i * (sigma + 0.366*10**-4)**j
        
    return eta_prim*h_star

def h_prim_s_r3(s):
    """
    Region 3-4 boundary (3a) eq. 4 (SR4-04)
    Saturated liquid line from the 623.15 K to critical point
    
    Validity range:
        s'(623.15) <= s < s_crit
        where:
        s'(623.15) = 3.778 281 340 kJ/kg/K (s_r3_lower)
        s_crit =     4.412 021 482 kJ/kg/K
    
    Parameters
    ----------
    s:      double  Entropy (J/kg/K)

    Returns
    -------
    double Enthalpy (J/kg)

    """
    #Reducing quantaties
    h_star = 1700 * 10**3 # J/kg
    s_star = 3.8 * 10**3 # J/kg/K
    
    #eta = h/h_star
    sigma = s/s_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [0, 0, 0,  0,  2, 3,  4, 4,  5,  5,  6, 7, 7,  7,  10, 10, 10, 32, 32]
    J = [1, 4, 10, 16, 1, 36, 3, 16, 20, 36, 4, 2, 28, 32, 14, 32, 36, 0,  6]
    n_array = [ 0.822673364673336,
                0.181977213534479,
               -0.112000260313624*10**-1,
               -0.746778287048033*10**-3,
               -0.179046263257381,
                0.424220110836657*10**-1,
               -0.341355823438768,
               -0.209881740853565*10**1,
               -0.822477343323596*10**1,
               -0.499684082076008*10**1,
                0.191413958471069,
                0.581062241093136*10**-1,
               -0.165505498701029*10**4,
                0.158870443421201*10**4,
               -0.850623535172818*10**2,
               -0.317714386511207*10**5,
               -0.945890406632871*10**5,
               -0.139273847088690*10**-5,
                0.631052532240980]
    
    eta_prim = 0
    for i,j,n in zip(I, J, n_array):
        eta_prim += n * (sigma - 1.09)**i * (sigma + 0.366*10**-4)**j
        
    return eta_prim*h_star

def h_s_b13(s):
    """
    Region 1-3 boundary eq. 7 (SR4-04)
    enthalpy as a function of entropy for the isotherm
    T = 623.15 K from the saturated liquid line up to 100 MPa.
    
    Validity range:
        s(100 MPa,623.15 K) <= s <= s'(623.15)
        where:
        s'(623.15) = 3.778 281 340 kJ/kg/K (s_r3_lower)
        s(100 MPa,623.15 K) = 3.397 782 955 kJ/kg/K
    
    Parameters
    ----------
    s:      double  Entropy (J/kg/K)

    Returns
    -------
    double Enthalpy (J/kg)

    """
    #Reducing quantaties
    h_star = 1700 * 10**3 # J/kg
    s_star = 3.8 * 10**3 # J/kg/K
    
    #eta = h/h_star
    sigma = s/s_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [0, 1, 1, 3, 5, 6]
    J = [0, -2, 2, -12, -4, -3, ]
    n_array = [ 0.913965547600543,
               -0.430944856041991*10**-4,
                0.603235694765419*10**2,
                0.117518273082168*10**-17,
                0.220000904781292,
               -0.690815545851641*10**2]
    
    eta = 0
    for i,j,n in zip(I, J, n_array):
        eta += n * (sigma - 0.884)**i * (sigma - 0.864)**j
        
    return eta*h_star

# Vapor

def h_bis_s_r2(s):
    """
    Region 2-4 boundary (2a, 2b) eq. 5 (SR4-04)
    Saturated vapor line from the triple point to 5.85 kJ/kg/K
    
    Validity range:
        5.85 kJ/kg/K <= s <= s''(273.15)
        where:
        s''(273.15) = 9.155 759 395 kJ/kg/K
    
    Parameters
    ----------
    s:      double  Entropy (J/kg/K)

    Returns
    -------
    double Enthalpy (J/kg)

    """
    #Reducing quantaties
    h_star = 2800 * 10**3 # J/kg
    s_star_1 = 5.21 * 10**3 # J/kg/K
    s_star_2 = 9.2 * 10**3 # J/kg/K
    
    sigma_1 = s/s_star_1
    sigma_2 = s/s_star_2
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [1, 1, 2, 2, 4, 4, 7, 8, 8, 10, 12, 12, 18, 20, 24, 28, 28, 28, 28, 28, 32, 32, 32, 32, 32, 36, 36, 36, 36, 36]
    J = [8, 24, 4, 32, 1, 2, 7, 5, 12, 1, 0, 7, 10, 12, 32, 8, 12, 20, 22, 24, 2, 7, 12, 14, 24, 10, 12, 20, 22, 28]
    n_array = [-0.524581170928788*10**3,
               -0.926947218142218*10**7,
               -0.237385107491666*10**3,
                0.210770155812776*10**11,
               -0.239494562010986*10**2,
                0.221802480294197*10**3,
               -0.510472533393438*10**7,
                0.124981396109147*10**7,
                0.200008436996201*10**10,
               -0.815158509791035*10**3,
               -0.157612685637523*10**3,
               -0.114200422332791*10**11,
                0.662364680776872*10**16,
               -0.227622818296144*10**19,
               -0.171048081348406*10**32,
                0.660788766938091*10**16,
                0.166320055886021*10**23,
               -0.218003784381501*10**30,
               -0.787276140295618*10**30,
                0.151062329700346*10**32,
                0.795732170300541*10**7,
                0.131957647355347*10**16,
               -0.325097068299140*10**24,
               -0.418600611419248*10**26,
                0.297478906557467*10**35,
               -0.953588761745473*10**20,
                0.166957699620939*10**25,
               -0.175407764869978*10**33,
                0.347581490626396*10**35,
               -0.710971318427851*10**39]
    
    eta_bis = 0
    for i,j,n in zip(I, J, n_array):
        eta_bis += n * (sigma_1**-1 - 0.513)**i * (sigma_2 - 0.524)**j
    eta_bis = np.exp(eta_bis)
        
    return eta_bis*h_star

def h_bis_s_r23(s):
    """
    Region 2,3-4 boundary (2c,3b) eq. 6 (SR4-04)
    Saturated vapor line in the entropy range
    
    Validity range:
        s_crit <= s <= 5.85 kJ/kg/K
        where:
        s_crit =     4.412 021 482 kJ/kg/K
        
    Parameters
    ----------
    s:      double  Entropy (J/kg/K)

    Returns
    -------
    double Enthalpy (J/kg)

    """
    h_star = 2800 * 10**3 # J/kg
    s_star = 5.9 * 10**3 # J/kg/K
    
    sigma = s/s_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [0, 0, 0, 1, 1,  5,  6 , 7,  8, 8,  12, 16, 22, 22, 24, 36]
    J = [0, 3, 4, 0, 12, 36, 12, 16, 2, 20, 32, 36, 2,  32, 7,  20]
    n_array = [ 0.104351280732769*10**1,
               -0.227807912708513*10**1,
                0.180535256723202*10**1,
                0.420440834792042,
               -0.105721244834660*10**6,
                0.436911607493884*10**25,
               -0.328032702839753*10**12,
               -0.678686760804270*10**16,
                0.743957464645363*10**4,
               -0.356896445355761*10**20,
                0.167590585186801*10**32,
               -0.355028625419105*10**38,
                0.396611982166538*10**12,
               -0.414716268484468*10**41,
                0.359080103867382*10**19,
               -0.116994334851995*10**41]
    
    eta_bis = 0
    for i,j,n in zip(I, J, n_array):
        eta_bis += n * (sigma - 1.02)**i * (sigma - 0.726)**j
    eta_bis = eta_bis**4
        
    return eta_bis*h_star

def T_hs_b23(h,s):
    """
    Region 2-3 boundary eq. 8 (SR4-04)
    Backward equation for temperature
    
    The range of validity of the equation T_hs_b23 is from the saturated vapor line x=1 up
    to 100 MPa in the entropy range:
    
        s_min_b23 <= s <= s_max_b23
        where:
        s_min_b23 = 5.048 096 828 kJ/kg/K
        s_max_b23 = 5.260 578 707 kJ/kg/K
        
    and in enthalpy range:
        h_min_b23 <= h <= h_max_b23
        where
        h_min_b23 = 2.563 592 004*10**3 kJ/kg
        h_max_b23 = 2.812 942 061*10**3 kJ//kg
    
    Parameters
    ----------
    h :     double  Enthalpy (J/kg)
    s:      double  Entropy (J/kg/K)

    Returns
    -------
    double Temperature (K)

    """
    #Reducing quantaties
    h_star = 3000 * 10**3 # J/kg
    s_star = 5.3 * 10**3 # J/kg/K
    T_star = 900 # K
    
    #pi = p/p_star
    eta = h/h_star
    sigma = s/s_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [-12, -10, -8, -4, -3, -2, -2, -2, -2,   0,  1,  1, 1, 3, 3, 5, 6, 6, 8, 8, 8, 12, 12, 14, 14]
    J = [10, 8, 3, 4, 3, -6, 2, 3, 4, 0, -3, -2, 10, -2, -1, -5, -6, -3, -8, -2, -1, -12, -1, -12, 1]
    n_array = [ 0.629096260829810*10**-3,
               -0.823453502583165*10**-3,
                0.515446951519474*10**-7,
               -0.117565945784945*10**1,
                0.348519684726192*10**1,
               -0.507837382408313*10**-11,
               -0.284637670005479*10**1,
               -0.236092263939673*10**1,
                0.601492324973779*10**1,
                0.148039650824546*10**1,
                0.360075182221907*10**-3,
               -0.126700045009952*10**-1,
               -0.122184332521413*10**7,
                0.149276502463272,
                0.698733471798484,
               -0.252207040114321*10**-1,
                0.147151930985213*10**-1,
               -0.108618917681849*10**1,
               -0.936875039816322*10**-3,
                0.819877897570217*10**2,
               -0.182041861521835*10**3,
                0.261907376402688*10**-5,
               -0.291626417025961*10**5,
                0.140660774926165*10**-4,
                0.783237062349385*10**7]
    
    omega = 0
    for i,j,n in zip(I, J, n_array):
        omega += n * (eta - 0.727)**i * (sigma - 0.864)**j
        
    return omega*T_star

################
### Region 1 ###
### SR2-02   ###
################

def p_hs_r1(h,s):
    """
    Backwards equation p_hs_r1(h,s), region 1. eq. 1
    Gives pressure from h,s input.
    
    Parameters
    ----------
    h :     double  Enthalpy (J/kg)
    s:      double  Entropy (J/kg/K)

    Returns
    -------
    double Pressure (Pa)

    """
    #Reducing quantaties
    p_star = 100 * 10**6 # Pa
    h_star = 3400 * 10**3 # J/kg
    s_star = 7.6 * 10**3 # J/kg/K
    
    #pi = p/p_star
    eta = h/h_star
    sigma = s/s_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [0, 0, 0, 0, 0, 0, 0, 0,  1, 1, 1, 1, 2, 2, 2,  3, 4, 4, 5]
    J = [0, 1, 2, 4, 5, 6, 8, 14, 0, 1, 4, 6, 0, 1, 10, 4, 1, 4, 0]
    n_array = [ -0.691997014660582,
                -0.183612548787560*10**2,
                -0.928332409297335*10**1,
                 0.659639569909906*10**2,
                -0.162060388912024*10**2,
                 0.450620017338667*10**3,
                 0.854680678224170*10**3,
                 0.607523214001162*10**4,
                 0.326487682621856*10**2,
                -0.269408844582931*10**2,
                -0.319947848334300*10**3,
                -0.928354307043320*10**3,
                 0.303634537455249*10**2,
                -0.650540422444146*10**2,
                -0.430991316516130*10**4,
                -0.747512324096068*10**3,
                 0.730000345529245*10**3,
                 0.114284032569021*10**4,
                -0.436407041874559*10**3]
    
    pi = 0
    for i,j,n in zip(I, J, n_array):
        pi += n * (eta + 0.05)**i * (sigma + 0.05)**j
        
    return pi*p_star

################
### Region 2 ###
### SR2-02   ###
################

def h_b2ab(s):
    """
    Backwards equation, boundary calculation betwwen region 2a/2b
    Approximatios isobaric line, p = 4.0 MPa
    Region 2 IAPWS-IF97
    Chapter 6.1, eq. 2
    
    Parameters
    ----------
    s:      double  Entropy (J/kg/K).

    Returns
    -------
    double  Enthalpy (J/kg)

    """
    #Reducing quantaties
    h_star = 1*10**3 # 1 kJ/kg
    s_star = 1*10**3 # 1 kJ/kg/K
        
    # Numerical coefficients from table 19
    n = [-0.349898083432139*10**4, 0.257560716905876*10**4, -0.421073558227969*10**3, 0.276349063799944*10**2]
    
    sigma = s/s_star
    eta = n[0] + n[1]*sigma + n[2]*sigma**2 + n[3]*sigma**3
    
    return eta * h_star

def p_hs_r2a(h,s):
    """
    Backwards equation p_hs_r2a(h,s), region 2a. eq. 3
    Gives pressure from h,s input.
    
    Parameters
    ----------
    h :     double  Enthalpy (J/kg)
    s:      double  Entropy (J/kg/K)

    Returns
    -------
    double Pressure (Pa)

    """
    #Reducing quantaties
    p_star = 4 * 10**6 # Pa
    h_star = 4200 * 10**3 # J/kg
    s_star = 12 * 10**3 # J/kg/K
     
    eta = h/h_star
    sigma = s/s_star
    
    #    1, 2, 3, 4,  5,  6,  7, 8, 9, 10,11,12,13, 14, 15, 16, 17,18, 19, 20,21,22,23,24, 25, 26,27, 28,29, 30, 31
    I = [0, 0, 0, 0,  0,  0,  1, 1, 1, 1, 1, 1, 1,  1,  1,  1,  2, 2,  2,  3, 3, 3, 3, 3,  4,  5, 5,  6, 7]
    J = [1, 3, 6, 16, 20, 22, 0, 1, 2, 3, 5, 6, 10, 16, 20, 22, 3, 16, 20, 0, 2, 3, 6, 16, 16, 3, 16, 3, 1]
    n_array = [-0.182575361923032*10**-1,
               -0.125229548799536,
                0.592290437320145,
                0.604769706185122*10**1,
                0.238624965444474*10**3,
               -0.298639090222922*10**3,
                0.512250813040750*10**-1,
               -0.437266515606486,
                0.413336902999504,
               -0.516468254574773*10**1,
               -0.557014838445711*10**1,
                0.128555037824478*10**2,
                0.114144108953290*10**2,
               -0.119504225652714*10**3,
               -0.284777985961560*10**4,
                0.431757846408006*10**4,
                0.112894040802650*10**1,
                0.197409186206319*10**4,
                0.151612444706087*10**4,
                0.141324451421235*10**-1,
                0.585501282219601,
               -0.297258075863012*10**1,
                0.594567314847319*10**1,
               -0.623656565798905*10**4,
                0.965986235133332*10**4,
                0.681500934948134*10**1,
               -0.633207286824489*10**4,
               -0.558919224465760*10**1,
                0.400645798472063*10**-1
                   ]
    
    pi = 0
    for i,j,n in zip(I, J, n_array):
        pi += n * (eta - 0.5)**i * (sigma - 1.2)**j
        
    return pi**4*p_star


def p_hs_r2b(h,s):
    """
    Backwards equation p_hs_r2b(h,s), region 2b. eq. 4
    Gives pressure from h,s input.
    
    Parameters
    ----------
    h :     double  Enthalpy (J/kg)
    s:      double  Entropy (J/kg/K)

    Returns
    -------
    double Pressure (Pa)

    """
    #Reducing quantaties
    p_star = 100 * 10**6 # Pa
    h_star = 4100 * 10**3 # J/kg
    s_star = 7.9 * 10**3 # J/kg/K
    
    eta = h/h_star
    sigma = s/s_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,  2, 2, 2,  3, 3, 3, 3,  4, 4,  5, 5,  6, 6, 6,  7, 7,  8, 8, 8,  8,  12, 14]
    J = [0, 1, 2, 4, 8, 0, 1, 2, 3, 5, 12, 1, 6, 18, 0, 1, 7, 12, 1, 16, 1, 12, 1, 8, 18, 1, 16, 1, 3, 14, 18, 10, 16]
    n_array = [ 0.801496989929495*10**-1,
               -0.543862807146111,
                0.337455597421283,
                0.890555451157450*10**1,
                0.313840736431485*10**3,
                0.797367065977789,
               -0.121616973556240*10**1,
                0.872803386937477*10**1,
               -0.169769781757602*10**2,
               -0.186552827328416*10**3,
                0.951159274344237*10**5,
               -0.189168510120494*10**2,
               -0.433407037194840*10**4,
                0.543212633012715*10**9,
                0.144793408386013,
                0.128024559637516*10**3,
               -0.672309534071268*10**5,
                0.336972380095287*10**8,
               -0.586634196762720*10**3,
               -0.221403224769889*10**11,
                0.171606668708389*10**4,
               -0.570817595806302*10**9,
               -0.312109693178482*10**4,
               -0.207841384633010*10**7,
                0.305605946157786*10**13,
                0.322157004314333*10**4,
                0.326810259797295*10**12,
               -0.144104158934487*10**4,
                0.410694867802691*10**3,
                0.109077066873024*10**12,
               -0.247964654258893*10**14,
                0.188801906865134*10**10,
               -0.123651009018773*10**15]
    
    pi = 0
    for i,j,n in zip(I, J, n_array):
        pi += n * (eta - 0.6)**i * (sigma - 1.01)**j
        
    return pi**4*p_star


def p_hs_r2c(h,s):
    """
    Backwards equation p_hs_r2c(h,s), region 2c. eq. 5
    Gives pressure from h,s input.
    
    Parameters
    ----------
    h :     double  Enthalpy (J/kg)
    s:      double  Entropy (J/kg/K)

    Returns
    -------
    double Pressure (Pa)

    """
    #Reducing quantaties
    p_star = 100 * 10**6 # Pa
    h_star = 3500 * 10**3 # J/kg
    s_star = 5.9 * 10**3 # J/kg/K
    
    eta = h/h_star
    sigma = s/s_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,  2, 2, 2, 2,  2,  3, 3, 3, 3,  3,  4,  5, 5, 5, 5,  6, 6,  10, 12, 16]
    J = [0, 1, 2, 3, 4, 8, 0, 2, 5, 8, 14, 2, 3, 7, 10, 18, 0, 5, 8, 16, 18, 18, 1, 4, 6, 14, 8, 18, 7,  7,  10]
    n_array = [ 0.112225607199012,
               -0.339005953606712*10**1,
               -0.320503911730094*10**2,
               -0.197597305104900*10**3,
               -0.407693861553446*10**3,
                0.132943775222331*10**5,
                0.170846839774007*10**1,
                0.373694198142245*10**2,
                0.358144365815434*10**4,
                0.423014446424664*10**6,
               -0.751071025760063*10**9,
                0.523446127607898*10**2,
               -0.228351290812417*10**3,
               -0.960652417056937*10**6,
               -0.807059292526074*10**8,
                0.162698017225669*10**13,
                0.772465073604171,
                0.463929973837746*10**5,
               -0.137317885134128*10**8,
                0.170470392630512*10**13,
               -0.251104628187308*10**14,
                0.317748830835520*10**14,
                0.538685623675312*10**2,
               -0.553089094625169*10**5,
               -0.102861522421405*10**7,
                0.204249418756234*10**13,
                0.273918446626977*10**9,
               -0.263963146312685*10**16,
               -0.107890854108088*10**10,
               -0.296492620980124*10**11,
               -0.111754907323424*10**16]

    pi = 0
    for i,j,n in zip(I, J, n_array):
        pi += n * (eta - 0.7)**i * (sigma - 1.1)**j
        
    return pi**4*p_star


def p_hs_r2(h,s):
    """
    Finds correct subregion for backwards equations 2a-c
    and returns result of applicable eqation
    
    For pressures below 4 MPa (h_b2ab(s)) -> region a
    if entropy s >= 5.85 kJ/kg/K -> region b
    s < 5.85 -> region c
    
    Parameters
    ----------
    h:      double  Enthalpy (J/kg).
    s :     double  Entropy (J/kg/k)

    Returns
    -------
    double Pressure (Pa)

    """
    
    if h <= h_b2ab(s):
        p = p_hs_r2a(h,s)
    elif s >= global_property.s_h_bis_585:
        p = p_hs_r2b(h,s)
    else:
        p = p_hs_r2c(h,s)
    
    return p

################
### Region 3 ###
### SR4-04   ###
################

def p_hs_r3a(h,s):
    """
    Backwards equation p_hs_r3a(h,s), region 3a. eq. 1
    Gives pressure from h,s input.
    
    Parameters
    ----------
    h :     double  Enthalpy (J/kg)
    s:      double  Entropy (J/kg/K)

    Returns
    -------
    double Pressure (Pa)

    """
    #Reducing quantaties
    p_star = 99 * 10**6 # Pa
    h_star = 2300 * 10**3 # J/kg
    s_star = 4.4 * 10**3 # J/kg/K
    
    eta = h/h_star
    sigma = s/s_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 6, 7, 8, 10, 10, 14, 18, 20, 22, 22, 24, 28, 28, 32, 32]
    J = [0, 1, 5, 0, 3, 4, 8, 14, 6, 16, 0, 2, 3, 0, 1, 4, 5, 28, 28, 24, 1, 32, 36, 22, 28, 36, 16, 28, 36, 16, 36, 10, 28]
    n_array = [ 0.770889828326934*10**1,
               -0.260835009128688*10**2,
                0.267416218930389*10**3,
                0.172221089496844*10**2,
               -0.293542332145970*10**3,
                0.614135601882478*10**3,
               -0.610562757725674*10**5,
               -0.651272251118219*10**8,
                0.735919313521937*10**5,
               -0.116646505914191*10**11,
                0.355267086434461*10**2,
               -0.596144543825955*10**3,
               -0.475842430145708*10**3,
                0.696781965359503*10**2,
                0.335674250377312*10**3,
                0.250526809130882*10**5,
                0.146997380630766*10**6,
                0.538069315091534*10**20,
                0.143619827291346*10**22,
                0.364985866165994*10**20,
               -0.254741561156775*10**4,
                0.240120197096563*10**28,
               -0.393847464679496*10**30,
                0.147073407024852*10**25,
               -0.426391250432059*10**32,
                0.194509340621077*10**39,
                0.666212132114896*10**24,
                0.706777016552858*10**34,
                0.175563621975576*10**42,
                0.108408607429124*10**29,
                0.730872705175151*10**44,
                0.159145847398870*10**25,
                0.377121605943324*10**41]

    pi = 0
    for i,j,n in zip(I, J, n_array):
        pi += n * (eta - 1.01)**i * (sigma - 0.750)**j
        
    return pi*p_star

def p_hs_r3b(h,s):
    """
    Backwards equation p_hs_r3b(h,s), region 3b. eq. 2
    Gives pressure from h,s input.
    
    Parameters
    ----------
    h :     double  Enthalpy (J/kg)
    s:      double  Entropy (J/kg/K)

    Returns
    -------
    double Pressure (Pa)

    """
    #Reducing quantaties
    p_star = 16.6 * 10**6 # Pa
    h_star = 2800 * 10**3 # J/kg
    s_star = 5.3 * 10**3 # J/kg/K
    
    eta = h/h_star
    sigma = s/s_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [-12, -12, -12, -12, -12, -10, -10, -10, -10, -8, -8, -6, -6, -6, -6, -5, -4, -4, -4, -3, -3, -3, -3, -2, -2, -1, 0, 2, 2, 5, 6, 8, 10, 14, 14]
    J = [ 2,   10,  12,  14,  20,  2,   10,  14,  18,  2,  8,  2,  6,  7,  8,  10, 4,  5,  8,  1,  3,  5,  6,  0,  1,  0, 3, 0, 1, 0, 1, 1, 1,  3,  7]
    n_array = [ 0.125244360717979*10**-12,
               -0.126599322553713*10**-1,
                0.506878030140626*10**1,
                0.317847171154202*10**2,
               -0.391041161399932*10**6,
               -0.975733406392044*10**-10,
               -0.186312419488279*10**2,
                0.510973543414101*10**3,
                0.373847005822362*10**6,
                0.299804024666572*10**-7,
                0.200544393820342*10**2,
               -0.498030487662829*10**-5,
               -0.102301806360030*10**2,
                0.552819126990325*10**2,
               -0.206211367510878*10**3,
               -0.794012232324823*10**4,
                0.782248472028153*10**1,
               -0.586544326902468*10**2,
                0.355073647696481*10**4,
               -0.115303107290162*10**-3,
               -0.175092403171802*10**1,
                0.257981687748160*10**3,
               -0.727048374179467*10**3,
                0.121644822609198*10**-3,
                0.393137871762692*10**-1,
                0.704181005909296*10**-2,
               -0.829108200698110*10**2,
               -0.265178818131250,
                0.137531682453991*10**2,
               -0.522394090753046*10**2,
                0.240556298941048*10**4,
               -0.227361631268929*10**5,
                0.890746343932567*10**5,
               -0.239234565822486*10**8,
                0.568795808129714*10**10]
    
    pi = 0
    for i,j,n in zip(I, J, n_array):
        pi += n * (eta - 0.681)**i * (sigma - 0.792)**j
    pi = 1/pi
        
    return pi*p_star

def p_hs_r3(h,s):
    """
    Finds correct subregion for backwards equations 3a-b
    and returns result of applicable eqation
    
    region 3a and 3b seperated by critical entropy line, s_crit. 
    s<=s_crit = 3a
    s>s_crit = 3b
    s
    Parameters
    ----------
    h:      double  Enthalpy (J/kg).
    s :     double  Entropy (J/kg/k)

    Returns
    -------
    double Pressure (Pa)

    """
    
    if s<= global_property.s_crit:
        p = p_hs_r3a(h,s)
    else:
        p = p_hs_r3b(h,s)
    
    return p

################
### Region 4 ###
### SR4-04   ###
################

def T_sat_hs(h,s):
    """
    Backward equation for saturation temperature  eq. 9 (SR4-04)
    given h,s
    
    The range of validity:
        s >= s''(623.15K)
        where:
        s''(623.15K) = 5.210 887 825 kJ/kg/K
        
    Parameters
    ----------
    h :     double  Enthalpy (J/kg)
    s:      double  Entropy (J/kg/K)

    Returns
    -------
    double Temperature (K)

    """
    #Reducing quantaties
    h_star = 2800 * 10**3 # J/kg
    s_star = 9.2 * 10**3 # J/kg/K
    T_star = 550 # K
    
    eta = h/h_star
    sigma = s/s_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 8, 10, 10, 12, 14, 14, 16, 16, 18, 18, 18, 20, 28]
    J = [0, 3, 12, 0, 1, 2, 5, 0, 5, 8, 0, 2, 3, 4, 0, 1, 1, 2, 4, 16, 6, 8, 22, 1, 20, 36, 24, 1, 28, 12, 32, 14, 22, 36, 24, 36]
    n_array = [ 0.179882673606601,
               -0.267507455199603,
                0.116276722612600*10**1,
                0.147545428713616,
               -0.512871635973248,
                0.421333567697984,
                0.563749522189870,
                0.429274443819153,
               -0.335704552142140*10**1,
                0.108890916499278*10**2,
               -0.248483390456012,
                0.304153221906390,
               -0.494819763939905,
                0.107551674933261*10**1,
                0.733888415457688*10**-1,
                0.140170545411085*10**-1,
               -0.106110975998808,
                0.168324361811875*10**-1,
                0.125028363714877*10**1,
                0.101316840309509*10**4,
               -0.151791558000712*10**1,
                0.524277865990866*10**2,
                0.230495545563912*10**5,
                0.249459806365456*10**-1,
                0.210796467412137*10**7,
                0.366836848613065*10**9,
               -0.144814105365163*10**9,
               -0.179276373003590*10**-2,
                0.489955602100459*10**10,
                0.471262212070518*10**3,
               -0.829294390198652*10**11,
               -0.171545662263191*10**4,
                0.355777682973575*10**7,
                0.586062760258436*10**12,
               -0.129887635078195*10**8,
                0.317247449371057*10**11]
    
    omega = 0
    for i,j,n in zip(I, J, n_array):
        omega += n * (eta - 0.119)**i * (sigma - 1.07)**j
        
    return omega*T_star

if __name__ == "__main__":
    main()
