# -*- coding: utf-8 -*-
"""
Module containing functions implemented from equations presentend in the R15-11
Release of IAPWS Formulation 2011 for the Thermal conductivity of Ordinary Water Substance.

List of functions and code structure:
    Determine within valid region
        region_tc_pt(p,t)
    sub-functions/ equation factors
        lambda_0(t_bar)
        lambda_1(t_bar,rho_bar)
        lambda_2(t_bar, rho_bar, cp_bar, my_bar, kappa, sigma_t_bar, sigma_tr_bar)
    main functions/equations
        lambda_trho(t, rho)
        lambda_tc(t, rho, cp, cv, my, sigma_t, sigma_tr)  
    suplementary functions
        sigma_tr(rho)

Comments regarding verification tables 7-9. None of the verification tables match perfectly
however the error is well within the uncertainty for the formulation. This is because of 
small uncertainties in the itterations and discrete derivative. The last verification
data, 647.35K rho = 322 kg/m^3, this is near the critical point, the differens here is much 
larger because the viscosity used in the verification, in the presented table the viscosity
used neglect the critical enhancement for the viscoity as suggesten using IAPWSIF97. 
However in this implementation the critical enhancement is still used near the critical region.

@author: Christoffer Rappmann, christoffer.rappmann@gmail.com
"""
import numpy as np

from .iapwsif97 import iapwsif97_globals as global_property

# Reference constants
t_star      = global_property.T_crit    # K         reference temperature
rho_star    = global_property.rho_crit  # kg/m^3    reference density
p_star      = global_property.p_crit    # Pa        reference pressure
my_star     = 1.00*10**-6               # Pa*s      reference viscosity
lambda_star = 1.00*10**-3               # W/K/m     reference thermal conductivity
R           = global_property.R_f11     # J/kg/K    specific gas constant eq-6, not the same value for R used in if97...

def region_tc_pt(p,t):
    """
    Determines if given pressure-temperature is within the applicable region.
    There is only one region in iapwsf11.

    Parameters
    ----------
    p :     double     pressure (Pa)
    t :     double     temperature (K)

    Returns
    -------
    boolean,        within applicable region (-)
    """
    # Tm(p) is assumed to be the same as the triple point temperature since
    # IAPWSif97 is only valid over 273.15, then 273.16 is low enough.
    region = False
    
    if p <= 0:
        region = False
    elif p < global_property.pt:
        if t >= global_property.T_r1_f11_lower and t <= global_property.T_r1_f11_upper:
            region = True
        else:
            region = False
    elif p <= global_property.p_r2_f11_upper:
        if t >= global_property.Tm_f11 and t <= global_property.T_r2_f11_upper:
            region = True
        else:
            region = False
    elif p <= global_property.p_r3_f11_upper:
        if t >= global_property.Tm_f11 and t <= global_property.T_r3_f11_upper:
            region = True
        else:
            region = False
    elif p <= global_property.p_r4_f11_upper:
        if t >= global_property.Tm_f11 and t <= global_property.T_r4_f11_upper:
            region = True
        else:
            region = False
    elif p <= global_property.p_r5_f11_upper:
        if t >= global_property.Tm_f11 and t <= global_property.T_r5_f11_upper:
            region = True
        else:
            region = False
    elif p <= global_property.p_r6_f11_upper:
        if t >= global_property.Tm_f11 and t <= global_property.T_r6_f11_upper:
            region = True
        else:
            region = False
    else:
        region = False
    
    return region
    
def lambda_0(t_bar):
    """
    Thermalconductivity in the dilute-gas limit, eq 16
    First factor in eq. 15.

    Parameters
    ----------
    t_bar :     double    dimensionless temperature (-)

    Returns
    -------
    double dimensionless thermal conductivity (-)
    """
    
    # Coefficients from table 1
    Lk = [2.443221*10**-3, 1.323095*10**-2, 6.770357*10**-3, -3.454586*10**-3, 4.096266*10**-4]
    
    lambda_sum = 0
    for k,L in enumerate(Lk):
        lambda_sum = lambda_sum + L/(t_bar**k)
        
    lambda_0 = (t_bar**.5)/lambda_sum
    
    return lambda_0

def lambda_1(t_bar,rho_bar):
    """
    Thermal conductivity due to finite density, eq 17
    Second factor in eq. 15.

    Parameters
    ----------
    t_bar :     double    dimensionless temperature (-)
    rho_bar:    double    dimensionless density (-)

    Returns
    -------
    double dimensionless thermal conductivity (-)
    """
    
    # Coefficients from Table 2
    Lij_array = np.array([[ 1.60397357, -0.646013523, 0.111443906, 0.102997357, -0.0504123634, 0.00609859258],
                          [ 2.33771842, -2.78843778,  1.53616167, -0.463045512,  0.0832827019,-0.00719201245],
                          [ 2.19650529, -4.54580785,  3.55777244, -1.40944978,   0.275418278, -0.0205938816],
                          [-1.21051378,  1.60812989, -0.621178141, 0.0716373224, 0.0,          0.0],
                          [-2.7203370,   4.57586331, -3.18369245,  1.1168348,   -0.19268305,   0.012913842]])
    
    lambda_sum = 0
    for i in range(0,len(Lij_array[:,0])):
        
        lambda_sum2 = 0
        for j,L in enumerate(Lij_array[i,:]):
            lambda_sum2 = lambda_sum2 + L*(rho_bar-1)**j
            
        lambda_sum = lambda_sum + ((1/t_bar)-1)**i*lambda_sum2
    
    lambda_1 = np.exp(rho_bar * lambda_sum)
    
    return lambda_1

def lambda_2(t_bar, rho_bar, cp_bar, my_bar, kappa, sigma_t_bar, sigma_tr_bar):
    """
    Critical enhancement of thermal conductivity, eq 18-24
    Second factor in eq. 15.

    Parameters
    ----------
    t_bar :         double    dimensionless temperature (-)
    rho_bar:        double    dimensionless density (-)
    cp_bar:         double    dimensionless isobaric heat capacity (-)
    my_bar:         double    dimensionless viscosity (-)
    kappa:          double    heat capacity ratio (-)
    sigma_t_bar:    double    dimensionless density-pressure derivative at constant t (-)
    sigma_tr_bar:   double    dimensionless density-pressure derivative at t-ref (-)

    Returns
    -------
    double dimensionless thermal conductivity (-)
    """
    # Critical region constants
    cap_lambda  = 177.8514  # (-)
    qd_bar      = 1/.4      # (nm)
    nu          = 0.630     # (-)
    gamma       = 1.239     # (-)
    xi0       = 0.13      # (nm)
    cap_gamma0  = 0.06      # (-)
    t_barr     = 1.5        # (-)
    
    # fot note 2 on page 12
    if sigma_t_bar > 10**13 or sigma_t_bar < 0:
        sigma_t_bar = 10**13
        
    if cp_bar > 10**13 or cp_bar < 0:
        cp_bar = 10**13
        
    # eq-23
    delta_chi_bar = rho_bar * (sigma_t_bar-sigma_tr_bar*t_barr/t_bar)
    
    if delta_chi_bar < 0:
        delta_chi_bar = 0
        
    # eq-22
    xi = xi0 * (delta_chi_bar/cap_gamma0)**(nu/gamma)
    
    # eq-20
    y = qd_bar * xi
    
    # eq-21, set to zero when y is small to avoid truncation error
    if y <1.2*10**-7:
        z = 0.0
    else:
        # eq-19
        z = 2/(np.pi*y)*(((1-1/kappa)*np.arctan(y)+1/kappa*y)-(1-np.exp(-1/(1/y+y**2/(3*rho_bar**2)))))
     
    # eq-18
    lambda_2 = cap_lambda*(rho_bar*cp_bar*t_bar)/my_bar*z
    
    return lambda_2

def lambda_trho(t, rho):
    """
    Thermal conductivity, eq 15.
    critical enhancement term lambda_2 asumed equal to 1

    Parameters
    ----------
    t :             double    temperature (K)
    rho:            double    density (kg/m^3)
   
    Returns
    -------
    double      thermal conductivity (-)
    """   
    
    # dimensionless quantities eq. 7-13
    t_bar           = t/t_star
    rho_bar         = rho/rho_star
    
    # eq 15 where lambda_1 = 1
    return lambda_0(t_bar)*lambda_1(t_bar, rho_bar)*lambda_star
    
def lambda_tc(t, rho, cp, cv, my, sigma_t, sigma_tr):
    """
    Thermal conductivity, eq 15.
    with critical enhancement term

    Parameters
    ----------
    t :             double    temperature (K)
    rho:            double    density (kg/m^3)
    cp:             double    isobaric heat capacity (J/kg/K)
    cv:             double    isochoric heat capacity (J/kg/K)
    my:             double    dimensionless viscosity (Pa*s)
    sigma_t:        double    density-pressure derivative at constant t (kg/m^3/Pa)
    sigma_tr:       double    density-pressure derivative at t-ref (kg/m^3/Pa)

    Returns
    -------
    double      thermal conductivity (-)
    """ 
    
    # dimensionless quantities eq. 7-13
    t_bar           = t/t_star
    rho_bar         = rho/rho_star
    cp_bar          = cp/R
    my_bar          = my/my_star
    sigma_t_bar     = sigma_t*(p_star/rho_star)
    sigma_tr_bar    = sigma_tr*(p_star/rho_star)
    
    kappa           = cp/cv
    
    lambda2 = lambda_2(t_bar, rho_bar, cp_bar, my_bar, kappa, sigma_t_bar, sigma_tr_bar)
    
    
    # eq 15 
    return (lambda_0(t_bar)*lambda_1(t_bar, rho_bar)+lambda2)*lambda_star

def sigma_tr(rho):
    """
    drho/dp derivative approximation at constant Tr
    eq-25
    is to be used in combination IAPWSIF97 for faster computingspeeds
    instead of using IAPWS95.

    Parameters
    ----------
    rho:            double    density (kg/m^3)

    Returns
    -------
    double      drho/dp derivative
    """
    
    # dimensionless quantities eq. 8
    rho_bar         = rho/rho_star
    
    # Coefficients from table 6
    #               j =    0                  1                  2                  3                   4  
    aij_array = np.array([[6.53786807199516,  6.52717759281799,  5.35500529896124,  1.55225959906681,   1.11999926419994],
                         [-5.61149954923348, -6.30816983387575, -3.96415689925446,  0.464621290821181,  0.595748562571649],
                         [ 3.39624167361325,  8.08379285492595,  8.91990208918795,  8.93237374861479,   9.88952565078920],
                         [-2.27492629730878, -9.82240510197603, -12.0338729505790, -11.0321960061126,  -10.3255051147040],
                         [ 10.2631854662709,  12.1358413791395,  9.19494865194302,  6.16780999933360,   4.66861294457414],
                         [ 1.97815050331519, -5.54349664571295, -2.16866274479712, -0.965458722086812, -0.503243546373828]])
    # Coefficients from eq-26, rho_bar<range_rho_bar for each j
    range_rho = np.array([ 0.0,               0.310559006,       0.776397516,       1.242236025,       1.863354037])
    
    # Determine range to be used from eq-26
    j = 0
    
    for k,rho_range in enumerate(range_rho):
        if rho_bar > rho_range:
            j = k

    # eq-25
    sigma_sum = 0
    for i,a in enumerate(aij_array[:,j]):
        sigma_sum = sigma_sum + a*rho_bar**i
    sigma_tr_bar = 1/sigma_sum
 
    return sigma_tr_bar/(p_star/rho_star)
        
if __name__ == "__main__":
    main()

    
    
    
    
    
    