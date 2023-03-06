# -*- coding: utf-8 -*-
"""
Module containing functions implemented from equations presentend in the R12-08 
release of IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance.

List of functions and code structure:
    Determine within valid region
        region_my_pt(p,t)
        crit_enhance_trho(t,rho)
    sub-functions/ equation factors
        my_0(t_bar)
        my_1(t_bar,rho_bar)
        my_2(t_bar, rho_bar, sigma_t_bar, sigma_tr_bar)
    main functions/equations
        my_trho(t,rho)
        my_trhosigma(t,rho, sigma_t, sigma_tr)        

@author: Christoffer Rappmann, christoffer.rappmann@gmail.com
"""

import numpy as np

from .iapwsif97 import iapwsif97_globals as global_property

# Reference constants
t_star      = global_property.T_crit    # K         reference temperature
rho_star    = global_property.rho_crit  # kg/m^3    reference density
p_star      = global_property.p_crit    # Pa        reference pressure
my_star     = 1.00*10**-6               # Pa*s      reference viscosity


        
def region_my_pt(p,t):
    """
    Determines if given pressure-temperature is within the applicable region.
    There is only one region in iapwsf08.

    Parameters
    ----------
    p :     double     pressure (Pa)
    t :     double     temperature (K)

    Returns
    -------
    boolean within applicable region (-)

    """
    # Tm(p) is assumed to be the same as the triple point temperature since
    # IAPWSif97 is only valid over 273.15, then 273.16 is low enough.
    region = False
    
    if p <= 0:
        region = False
    elif p < global_property.pt:
        if t >= global_property.T_r1_f08_lower and t <= global_property.T_r1_f08_upper:
            region = True
        else:
            region = False
    elif p <= global_property.p_r2_f08_upper:
        if t >= global_property.Tm_f08 and t <= global_property.T_r2_f08_upper:
            region = True
        else:
            region = False
    elif p <= global_property.p_r3_f08_upper:
        if t >= global_property.Tm_f08 and t <= global_property.T_r3_f08_upper:
            region = True
        else:
            region = False
    elif p <= global_property.p_r4_f08_upper:
        if t >= global_property.Tm_f08 and t <= global_property.T_r4_f08_upper:
            region = True
        else:
            region = False
    elif p <= global_property.p_r5_f08_upper:
        if t >= global_property.Tm_f08 and t <= global_property.T_r5_f08_upper:
            region = True
        else:
            region = False
    else:
        region = False
        
    return region

def crit_enhance_trho(t,rho):
    """
    Determines if given pressure-temperature is within the region where
    the critical enhancement term in iapwsf08 is significant

    Parameters
    ----------
    p :     double     pressure (Pa)
    t :     double     temperature (-)

    Returns
    -------
    boolean within significant region (-)

    """
    region = False
    
    if t < global_property.T_f08_upper and t > global_property.T_f08_lower and rho < global_property.rho_f08_upper and rho > global_property.rho_f08_lower:
        region = True
    else:
        region = False
    
    return region

def my_0(t_bar):
    """
    Viscosity in the dilute-gas limit, eq 11
    First factor in eq. 10.

    Parameters
    ----------
    t_bar :     double    dimensionless temperature (-)

    Returns
    -------
    double dimensionless viscosity (-)

    """
    
    # table 1 coeficients
    H = [1.67752, 2.20462, 0.6366564, -.241605]
    
    # eq 11
    htsum = 0
    for i,h in enumerate(H):
        htsum = htsum + h/t_bar**i
        
    return 100*(t_bar**.5)/htsum

def my_1(t_bar,rho_bar):
    """
    Viscosity due to finite density, eq 12
    Second factor in eq. 10.

    Parameters
    ----------
    t_bar :     double    dimensionless temperature (-)
    rho_bar:    double    dimensionless density (-)

    Returns
    -------
    double dimensionless viscosity (-)

    """
    
    # table 2
    ijH_array = np.array([[0, 1, 2, 3, 0, 1, 2, 3, 5, 0, 1, 2, 3, 4, 0, 1, 0, 3, 4, 3, 5],
                          [0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4, 4, 5, 6, 6],
                          [ 5.20094*10**-1,
                            8.50895*10**-2,
                           -1.08374,
                           -2.89555*10**-1,
                            2.22531*10**-1,
                            9.99115*10**-1,
                            1.88797,
                            1.26613,
                            1.20573*10**-1,
                           -2.81378*10**-1,
                           -9.06851*10**-1,
                           -7.72479*10**-1,
                           -4.89837*10**-1,
                           -2.57040*10**-1,
                            1.61913*10**-1, 
                            2.57399*10**-1,
                           -3.25372*10**-2,
                            6.98452*10**-2,
                            8.72102*10**-3,
                           -4.35673*10**-3,
                           -5.93264*10**-4]])
    sum1 = 0
    for i,j,h in zip(ijH_array[0,:],ijH_array[1,:],ijH_array[2,:]):
        sum1 = sum1 + (1/t_bar-1)**i*h*(rho_bar-1)**j
    
    return np.exp(rho_bar*sum1)

def my_2(t_bar, rho_bar, sigma_t_bar, sigma_tr_bar):
    """
    Critical enhancement of viscosity, eq 14-21
    Second factor in eq. 10.

    Parameters
    ----------
    t_bar :     double    dimensionless temperature (-)
    rho_bar:    double    dimensionless density (-)
    sigma_t_bar:    double    dimensionless density-pressure derivative at constant t (-)
    sigma_tr_bar:    double   dimensionless density-pressure derivative at t-ref (-)

    Returns
    -------
    double dimensionless viscosity (-)

    """
    # Critical region constants
    x_my        = .068  # -
    qc          = 1/1.9 # nm
    qd          = 1/1.1 # nm
    nu          = .630  # -
    gamma       = 1.239 # -
    xi0       = 0.13  # nm
    cap_gamma0   = 0.06  # -
    t_barr      = 1.5   # -
    
    # eq-21
    delta_chi_bar = rho_bar * (sigma_t_bar-sigma_tr_bar*t_barr/t_bar)
    
    if delta_chi_bar < 0:
        delta_chi_bar = 0
        
    # eq-20
    xi = xi0 * (delta_chi_bar/cap_gamma0)**(nu/gamma)
    # eq-17
    psi_d = np.arccos((1+qd**2*xi**2)**-.5)
    # eq-19
    w = np.abs((qc*xi-1)/(qc*xi+1))**.5*np.tan(psi_d/2)
    # eq-18
    if qc*xi > 1:
        l = np.log((1+w)/(1-w))
    else:
        l = 2*np.arctan(np.abs(w))
        # eq-15 eq-16
    if xi < 0:
        y = None
        print('error in iapwsf08 critical enhancement...')
    elif xi <= .3817016416:
        y = (1/5)*qc*xi*(qd*xi)**5*(1-qc*xi+(qc*xi)**2-(765/504)*(qd*xi)**2)
    else:
        y = (1/12)*np.sin(3*psi_d)-(1/(4*qc*xi))*np.sin(2*psi_d)+(1/(qc*xi)**2)*(1-(5/4)*(qc*xi)**2)*np.sin(psi_d)-(1/(qc*xi)**3)*((1-(3/2)*(qc*xi)**2)*psi_d-np.abs((qc*xi)**2-1)**(3/2)*l)
    
    # eq-14
    return np.exp(x_my*y)   

def my_trho(t,rho):
    """
    Viscosity according to eq. 10
    critical enhancement term my2 asumed equal to 1

    Parameters
    ----------
    t :             double    temperature (K)
    rho:            double    density (kg/m^3)
    Returns
    -------
    double viscosity (Pa*s)

    """
    
    # dimensionless quantities
    t_bar           = t/t_star
    rho_bar         = rho/rho_star
    
    #equation 10 where my2 = 1 assumed
    return my_0(t_bar)*my_1(t_bar,rho_bar)*my_star
     

def my_trhosigma(t,rho, sigma_t, sigma_tr):
    """
    Viscosity according to eq. 10
    with critical enhancement term

    Parameters
    ----------
    t :             double    temperature (K)
    rho:            double    density (kg/m^3)
    sigma_t:        double    density-pressure derivative at constant t (kg/m^3/Pa)
    sigma_tr:       double    density-pressure derivative at t-ref (kg/m^3/Pa)
    
    Returns
    -------
    double viscosity (Pa*s)

    """
    
    # dimensionless quantities
    t_bar           = t/t_star
    rho_bar         = rho/rho_star
    sigma_t_bar     = sigma_t*(p_star/rho_star)
    sigma_tr_bar    = sigma_tr*(p_star/rho_star)
    
    #equation 10
    my2 = my_2(t_bar,rho_bar, sigma_t_bar, sigma_tr_bar)
    
    return my_0(t_bar)*my_1(t_bar,rho_bar)*my2*my_star

if __name__ == "__main__":
    main()
        