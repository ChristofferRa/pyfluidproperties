# -*- coding: utf-8 -*-
"""
Module containing functions implemented from equations presentend in the R7-97 release of IAPWSIF97.

List of functions and code structure:
    
    Determine region and region boundaries
        p_t_b23(T, n = n_array_b23)
        t_p_b23(p, n = n_array_b23)
        find_region(p,T)
        find_region_ph(p,h)
        find_region_ps(p,s)
        
    Region 1 
        Equations of state
            gamma_r1(pi,tau, ijn=ijn_array_r1)
            gamma_r1_pi(pi,tau,  ijn=ijn_array_r1)
            gamma_r1_pipi(pi,tau, ijn=ijn_array_r1)
            gamma_r1_tau(pi,tau, ijn=ijn_array_r1)
            gamma_r1_tautau(pi,tau, ijn=ijn_array_r1)
            gamma_r1_pitau(pi,tau, ijn=ijn_array_r1)
        Properties
            v_r1(p,T)
            u_r1(p,T)
            s_r1(p,T)
            h_r1(p,T)
            cp_r1(p,T)
            cv_r1(p,T)
            w_r1(p,T)
        Backward equations
            t_ph_r1(p,h)
            t_ps_r1(p,s)
    
    Region 2
        Equations of state
            gamma_ideal_r2(pi,tau,jn = jn_array_ideal_r2)
            gamma_ideal_r2_pi(pi,tau, jn = jn_array_ideal_r2)
            gamma_ideal_r2_pipi(pi,tau, jn = jn_array_ideal_r2)
            gamma_ideal_r2_tau(pi, tau, jn = jn_array_ideal_r2)
            gamma_ideal_r2_tautau(pi, tau, jn = jn_array_ideal_r2)
            gamma_ideal_r2_pitau(pi, tau, jn = jn_array_ideal_r2)
            gamma_residual_r2(pi, tau, ijn = ijn_array_residual_r2)
            gamma_residual_r2_pi(pi, tau, ijn = ijn_array_residual_r2)
            gamma_residual_r2_pipi(pi, tau, ijn = ijn_array_residual_r2)
            gamma_residual_r2_tau(pi, tau, ijn = ijn_array_residual_r2)
            gamma_residual_r2_tautau(pi, tau, ijn = ijn_array_residual_r2)
            gamma_residual_r2_pitau(pi, tau, ijn = ijn_array_residual_r2)
            
        Properties
            v_r2(p,T)
            u_r2(p,T)
            s_r2(p,T)
            h_r2(p,T)
            cp_r2(p,T)
            cv_r2(p,T)
            w_r2(p,T)
        
            v_r2_metastable(p,T)
            u_r2_metastable(p,T)
            s_r2_metastable(p,T)
            h_r2_metastable(p,T)
            cp_r2_metastable(p,T)
            cv_r2_metastable(p,T)
            w_r2_metastable(p,T)
            
        Backward equations
            p_b2bc(h)
            h_b2bc(p)
            t_ph_r2a(p,h)
            t_ph_r2b(p,h)
            t_ph_r2c(p,h)
            t_ph_r2(p,h)
            t_ps_r2a(p,s)
            t_ps_r2b(p,s)
            t_ps_r2c(p,s)
            t_ps_r2(p,s)
            
    
    Region 3
        Equations of state
            phi_r3(delta, tau, ijn = ijn_array_r3,n_r3_1 = n_r3_1)
            phi_r3_delta(delta, tau, ijn = ijn_array_r3,n_r3_1 = n_r3_1)
            phi_r3_deltadelta(delta, tau, ijn = ijn_array_r3,n_r3_1 = n_r3_1)
            phi_r3_tau(delta, tau, ijn = ijn_array_r3,n_r3_1 = n_r3_1)
            phi_r3_tautau(delta, tau, ijn = ijn_array_r3,n_r3_1 = n_r3_1)
            phi_r3_deltatau(delta, tau, ijn = ijn_array_r3,n_r3_1 = n_r3_1)
            
        Properties
            p_r3(rho,T)
            u_r3(rho,T)
            s_r3(rho,T)
            h_r3(rho,T)
            cp_r3(rho,T)
            cv_r3(rho,T)
            w_r3(rho,T)
            
            rho_r3(p,T, rho_init = 0.0)
            
        Backward equations
            -
    
    Region 4
        Equations of state
        -
        
        Properties
            p_sat_t_r4(t)
            t_sat_p_r4(p)
        
        Backward equations
        -
        
    Region 5
        Equations of state
            gamma_ideal_r5(pi,tau,jn = jn_array_ideal_r5)
            gamma_ideal_r5_pi(pi,tau,jn = jn_array_ideal_r5)
            gamma_ideal_r5_pipi(pi,tau,jn = jn_array_ideal_r5)
            gamma_ideal_r5_tau(pi,tau,jn = jn_array_ideal_r5)
            gamma_ideal_r5_tautau(pi,tau,jn = jn_array_ideal_r5)
            gamma_ideal_r5_pitau(pi,tau,jn = jn_array_ideal_r5)
            gamma_residual_r5(pi, tau, ijn = ijn_array_residual_r5)
            gamma_residual_r5_pi(pi, tau, ijn = ijn_array_residual_r5)
            gamma_residual_r5_pipi(pi, tau, ijn = ijn_array_residual_r5)
            gamma_residual_r5_tau(pi, tau, ijn = ijn_array_residual_r5)
            gamma_residual_r5_tautau(pi, tau, ijn = ijn_array_residual_r5)
            gamma_residual_r5_pitau(pi, tau, ijn = ijn_array_residual_r5)
        
        Properties
            v_r5(p,T)
            u_r5(p,T)
            s_r5(p,T)
            h_r5(p,T)
            cp_r5(p,T)
            cv_r5(p,T)
            w_r5(p,T)
        
        Backward equations
            t_ph_r5(p,h)
            t_ps_r5(p,s)

@author: Christoffer Rappmann, christoffer.rappmann@gmail.com
"""
import numpy as np
from .iapwsif97_tps3_tph3_vph3_vps3 import p_3sat_h, p_3sat_s

# Helper functions
from . import iapwsif97_helper_functions as aux
# Global Constants
from . import iapwsif97_globals as global_property

# Specific gas constant
R = global_property.R #J/kg/K      Eq. 1

# Critical parameters
T_crit =    global_property.T_crit      # K          Eq. 2
p_c =       global_property.p_crit      # MPa        Eq. 3
rho_c =     global_property.rho_crit    # kg/m^3     Eq. 4

"""
    Auxiliary Equation for the 
    Boundary between Regions 2 and 3, 
    Chapter 4
"""
#Reducing quantaties
p_star_b23 = 1 * 10**6 #Pa
T_star_b23 = 1 # K

# Table 1, Numerical values of the coefficients of the B23-equation
n_array_b23 = np.array([ 0.34805185628969*10**3,
                        -0.11671859879975*10,
                         0.10192970039326*10**-2,
                         0.57254459862746*10**3,
                         0.13918839778870*10**2])

def p_t_b23(T, n = n_array_b23):
    """
    Boundary equation for region 2-3,
    Chapter 4, Eq. 5
    
    Parameters
    ----------
    T:     double       temperature (K)
    n      np array     coefficients from table 1

    Returns
    -------
    p :     double      corresponding pressure
    """
    
    teta = T/T_star_b23
    return (n[0]+n[1]*teta+n[2]*teta**2)*p_star_b23

def t_p_b23(p, n = n_array_b23):
    """
    Boundary equation for region 2-3,
    Chapter 4, Eq. 5
    
    Parameters
    ----------
    p:     double       pressure (Pa)
    n      np array     coefficients from table 1

    Returns
    -------
    T :     double      corresponding temperature
    """
    
    pi = p / p_star_b23
    return (n[3]+((pi-n[4])/n[2])**.5)*T_star_b23

"""
Determine applicable region
"""

def find_region(p,T):
    """
    The IAPWS Industrial Formulation 1997 consists of a set of equations for different regions
    which cover the following range of validity:
    273.15 K  <= T <= 1073.15 K, p <= 100 MPa
    1073.15 K < T  <= 2273.15 K, p <= 50 MPa
        
    This method determine the applicable region given conditions p,T

    Parameters
    ----------
    p : double      Pressure (Pa).
    T : double      Temperature (K).

    Returns
    -------
    region : int    Applicable region (-)
    """
        
    region = 0
        
    # Check which region should be used
    
    if np.isnan(p) or np.isnan(T):
        # Out of applicable region
        region = -1
    elif p <= global_property.p_r123_upper_boundary:
        if T >= global_property.T_boundary_1 and T <= global_property.T_boundary_25:
            # Regions 1-3
            if T >= global_property.T_boundary_1 and T <= global_property.T_boundary_13:
                if p >= p_sat_t_r4(T):
                    # Region 1
                    region = 1
                else:
                    # Region 2
                    region = 2
            else:
                # Regions 2-3
                if p < p_t_b23(T):
                    # Region 2
                    region = 2
                else:
                    # Region 3
                    region = 3
        elif T <= global_property.T_boundary_5 and p <= global_property.p_r5_upper_boundary:
            # Corresponds to region 5
            region = 5
        else:
            # Out of range
            region = -1
    else:
        # Out of range
        region = -1
            
    return region

def find_region_ph(p,h):
    """
    The IAPWS Industrial Formulation 1997 consists of a set of equations for different regions
    which cover the following range of validity:
    273.15 K  <= T <= 1073.15 K, p <= 100 MPa if97.h_r1(p*10**6,T_boundary_1)*10**-3
    1073.15 K < T  <= 2273.15 K, p <= 50 MPa, region 5
        
    This method determine the applicable region given conditions p,h
    the boundaries needs to be converted to enthalpy-lines for a given pressure

    Parameters
    ----------
    p : double      Pressure (Pa).
    h : double      Enthalpy (kJ/kg).

    Returns
    -------
    region : int    Applicable region (-)
    """
   
    ### Fixed boundaries ###
    # Temperature boundaries
    T_boundary_13 = global_property.T_boundary_13 # K
    T_boundary_25 = global_property.T_boundary_25 # K
    T_boundary_1 = global_property.T_boundary_1 # K
    T_boundary_5 = global_property.T_boundary_5 # K
    
    # Pressure
    p_boundary_123_upper = global_property.p_r123_upper_boundary
    p_boundary_134 = global_property.p_r3_lower_boundary # Pa 
    p_boundary_5 = global_property.p_r5_upper_boundary # Pa

    region = 0
        
    # Check which region should be used
    
    if np.isnan(p) or np.isnan(h):
        # Out of applicable region
        region = -1
    elif p <= p_boundary_123_upper:
        # In applicable pressure region
        if p <= p_boundary_134:
            # Region 1, 4, 2, 5
            if h < h_r1(p,T_boundary_1) or h > h_r5(p,T_boundary_5):
                # Out of applicable region
                region = -1
            elif h < h_r1(p,(t_sat_p_r4(p))):
                # Region 1
                region = 1
            elif h > h_r2(p,(t_sat_p_r4(p))) and h <= h_r2(p,T_boundary_25):
                # Region 2
                region = 2
            elif h > h_r2(p,T_boundary_25):
                # Region 5
                region = 5
            else:
                # Region 4, saturation properties
                region = 4
        elif p <= p_boundary_5:
            # Region 1, 3, 4, 2, 5
            if h < h_r1(p,T_boundary_1) or h > h_r5(p,T_boundary_5):
                # Out of applicable region
                region = -1
            elif h < h_r1(p,T_boundary_13):
                # Region 1
                region = 1
            elif h > h_r2(p, (t_p_b23(p))) and h <= h_r2(p,T_boundary_25):
                # Region 2
                region = 2
            elif h > h_r2(p,T_boundary_25):
                # Region 5
                region = 5
            elif p > p_3sat_h(h):
                # Region 3
                region = 3
            else:
                # Region 4
                region = 4
        else:
            # Region 1, 3, 4, 2
            if h < h_r1(p,T_boundary_1) or h > h_r2(p,T_boundary_25):
                # Out of applicable region
                region = -1
            elif h < h_r1(p,T_boundary_13):
                # Region 1
                region = 1
            elif h > h_r2(p, (t_p_b23(p))) and h <= h_r2(p,T_boundary_25):
                # Region 2
                region = 2
            elif p > p_3sat_h(h):
                # Region 3
                region = 3
            else:
                # Region 4
                region = 4
    else:
        # Out of range
        region = -1
            
    return region

def find_region_ps(p,s):
    """
    The IAPWS Industrial Formulation 1997 consists of a set of equations for different regions
    which cover the following range of validity:
    273.15 K  <= T <= 1073.15 K, p <= 100 MPa if97.s_r1(p*10**6,T_boundary_1)*10**-3
    1073.15 K < T  <= 2273.15 K, p <= 50 MPa, region 5, not implemented
        
    This method determine the applicable region given conditions p,h
    the boundaries needs to be converted to entropy-lines for a given pressure

    Parameters
    ----------
    p : double      Pressure (Pa).
    s : double      Entropy (kJ/kg/K).

    Returns
    -------
    region : int    Applicable region (-)
    """
    
    ### Fixed boundaries ###
    # Temperature boundaries
    T_boundary_13 = global_property.T_boundary_13 # K
    T_boundary_25 = global_property.T_boundary_25 # K
    T_boundary_1 = global_property.T_boundary_1 #K
    T_boundary_5 = global_property.T_boundary_5 # K
    
    # Pressure
    p_boundary_134 = global_property.p_r3_lower_boundary # Pa 
    p_boundary_5 = global_property.p_r5_upper_boundary # Pa
      
    region = 0
        
    # Check which region should be used
    
    if np.isnan(p) or np.isnan(s):
        # Out of applicable region
        region = -1
    elif p <= global_property.p_r123_upper_boundary:
        # In applicable pressure region
        if p <= p_boundary_134:
            # Region 1, 4, 2
            if s < s_r1(p,T_boundary_1) or s > s_r5(p,T_boundary_5):
                # Out of applicable region              
                region = -1
            elif s < s_r1(p,(t_sat_p_r4(p))):
                # Region 1
                region = 1
            elif s > s_r2(p,(t_sat_p_r4(p))) and s <= s_r2(p,T_boundary_25):
                # Region 2
                region = 2
            elif s > s_r2(p,T_boundary_25):
                # Region 5
                region = 5
            else:
                # Region 4, saturation properties
                region = 4
                
        elif p <= p_boundary_5:
            # Region 1, 3, 4, 2
            if s < s_r1(p,T_boundary_1) or s > s_r5(p,T_boundary_5):
                # Out of applicable region
                region = -1
            elif s < s_r1(p,T_boundary_13):
                # Region 1
                region = 1
            elif s > s_r2(p, (t_p_b23(p))) and s <= s_r2(p,T_boundary_25):
                # Region 2
                region = 2
            elif s > s_r2(p,T_boundary_25):
                # Region 5
                region = 5
            elif p > p_3sat_s(s):
                # Region 3
                region = 3
            else:
                # Region 4
                region = 4
        else:
            # Region 1, 3, 4, 2
            if s < s_r1(p,T_boundary_1) or s > s_r2(p,T_boundary_25):
                # Out of applicable region
                region = -1
            elif s < s_r1(p,T_boundary_13):
                # Region 1
                region = 1
            elif s > s_r2(p, (t_p_b23(p))) and s <= s_r2(p,T_boundary_25):
                # Region 2
                region = 2
            elif p > p_3sat_s(s):
                # Region 3
                region = 3
            else:
                # Region 4
                region = 4
    else:
        # Out of range
        region = -1
            
    return region

"""
    Region 1 
    IAPWS-IF97, Chapter 5

    Validity range:
    273.15 K <= T <= 623.15 K . 
    p_s(T)<= p <= 100 MPa
"""

#Reducing quantaties
p_star_r1 = 16.53*10**6          #MPa
T_star_r1 = 1386                 #K

#Gibbs free energy
# Table 2, Numerical values of the coefficients and exponents
ijn_array_r1 = np.array([[  0,  0,  0, 0, 0, 0, 0, 0, 1,  1,  1 , 1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  4,  4,  4,  5,  8,   8,  21,  23,  29,  30,  31,  32],
                        [  -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0,  1,  3,  -3, 0,  1,  3,  17, -4, 0,  6, -5, -2,  10, -8, -11, -6, -29, -31, -38, -39, -40, -41],
                        [0.14632971213167,
                    -0.84548187169114,
                    -0.37563603672040*10,
                     0.33855169168385*10,
                    -0.95791963387872,
                     0.15772038513228,
                    -0.16616417199501*10**-1,
                     0.81214629983568*10**-3,
                     0.28319080123804*10**-3,
                    -0.60706301565874*10**-3,
                    -0.18990068218419*10**-1,
                    -0.32529748770505*10**-1,
                    -0.21841717175414*10**-1,
                    -0.52838357969930*10**-4,
                    -0.47184321073267*10**-3,
                    -0.30001780793026*10**-3,
                     0.47661393906987*10**-4,
                    -0.44141845330846*10**-5,
                    -0.72694996297594*10**-15,
                    -0.31679644845054*10**-4,
                    -0.28270797985312*10**-5,
                    -0.85205128120103*10**-9,
                    -0.22425281908000*10**-5,
                    -0.65171222895601*10**-6,
                    -0.14341729937924*10**-12,
                    -0.40516996860117*10**-6,
                    -0.12734301741641*10**-8,
                    -0.17424871230634*10**-9,
                    -0.68762131295531*10**-18,
                     0.14478307828521*10**-19,
                     0.26335781662795*10**-22,
                    -0.11947622640071*10**-22,
                     0.18228094581404*10**-23,
                    -0.93537087292458*10**-25]])

def gamma_r1(pi,tau, ijn=ijn_array_r1):
    """
    Gibs free energy eq. 7, g(p,T)
    Region 1 of IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    ijn_table_r1 np array       I,J,n coefficients

    Returns
    -------
    gamma :   double    g(p,T), Gibbs free energy
    """
    
    gamma = 0
    for x,y,z in zip(ijn[0,:],ijn[1,:],ijn[2,:]):
        gamma += (z*(7.1-pi)**x)*((tau-1.222)**y)
    return gamma

def gamma_r1_pi(pi,tau,  ijn=ijn_array_r1):
    """
    Gibs free energy eq. 7 derivative dg(p,T)/dpi
    Region 1 of IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    ijn_table_r1 np array       I,J,n coefficients

    Returns
    -------
    gamma_pi :   double    dg(p,T)/dpi, Gibbs free energy pi derivative
    """

    # g(p,T), Gibbs free energy eq. 7
    gamma = 0
    for x,y,z in zip(ijn[0,:],ijn[1,:],ijn[2,:]):
        gamma += (-z*x*(7.1-pi)**(x-1))*((tau-1.222)**y)
    return gamma

def gamma_r1_pipi(pi,tau, ijn=ijn_array_r1):
    """
    Gibs free energy eq. 7 derivative d^2g(p,T)/dpi^2
    Region 1 of IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    ijn_table_r1 np array       I,J,n coefficients

    Returns
    -------
    gamma_pipi :   double    d^2g(p,T)/dpi^2, Gibbs free energy pipi derivative
    """
    
    gamma = 0
    for x,y,z in zip(ijn[0,:],ijn[1,:],ijn[2,:]):
        gamma += (z*x*(x-1)*(7.1-pi)**(x-2))*((tau-1.222)**y)
        
    return gamma

def gamma_r1_tau(pi,tau, ijn=ijn_array_r1):
    """
    Gibs free energy eq. 7 derivative dg(p,T)/dtau
    Region 1 of IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    ijn_table_r1 np array       I,J,n coefficients

    Returns
    -------
    gamma_tau :   double    dg(p,T)/dtau, Gibbs free energy tau derivative
    """
    
    gamma = 0
    for x,y,z in zip(ijn[0,:],ijn[1,:],ijn[2,:]):
        gamma += (z*(7.1-pi)**x)*(y*(tau-1.222)**(y-1))
        
    return gamma

def gamma_r1_tautau(pi,tau, ijn=ijn_array_r1):
    """
    Gibs free energy eq. 7 derivative d^2g(p,T)/dtau^2
    Region 1 of IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    ijn_table_r1 np array       I,J,n coefficients

    Returns
    -------
    gamma_tautau :   double    d^2g(p,T)/dtau^2, Gibbs free energy tautau derivative
    """
    
    gamma = 0
    for x,y,z in zip(ijn[0,:],ijn[1,:],ijn[2,:]):
        gamma += (z*(7.1-pi)**x)*(y*(y-1)*(tau-1.222)**(y-2))
        
    return gamma

def gamma_r1_pitau(pi,tau, ijn=ijn_array_r1):
    """
    Gibs free energy eq. 7 derivative d^2g(p,T)/dpi/dtau
    Region 1 of IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    ijn_table_r1 np array       I,J,n coefficients

    Returns
    -------
    gamma_pitau :   double    d^2g(p,T)/dpi/dtau, Gibbs free energy pitau derivative
    """
    
    gamma = 0
    for x,y,z in zip(ijn[0,:],ijn[1,:],ijn[2,:]):
        gamma += (-z*x*(7.1-pi)**(x-1))*(y*(tau-1.222)**(y-1))
        
    return gamma

def v_r1(p,T):
    """
    Specific volume for region 1 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    T:      double  temperature (K).

    Returns
    -------
    double  specific volume (m^3/kg)

    """
    pi = p/p_star_r1
    tau = T_star_r1/T
   
    return pi*gamma_r1_pi(pi, tau)*R*T/p

def u_r1(p,T):
    """
    Specific internal energy for region 1 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    T:      double  temperature (K).

    Returns
    -------
    double  specific internal energy (J/kg)

    """
    pi = p/p_star_r1
    tau = T_star_r1/T
    
    return (tau*gamma_r1_tau(pi,tau)-pi*gamma_r1_pi(pi, tau))*R*T

def s_r1(p,T):
    """
    Specific entropy for region 1 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    T:      double  temperature (K).

    Returns
    -------
    double  specific entropy (J/kg/K)

    """
    pi = p/p_star_r1
    tau = T_star_r1/T

    return (tau*gamma_r1_tau(pi, tau)-gamma_r1(pi, tau))*R   

def h_r1(p,T):
    """
    Specific enthalpy for region 1 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    T:      double  temperature (K).

    Returns
    -------
    double  specific enthalpy (J/kg)

    """
    pi = p/p_star_r1
    tau = T_star_r1/T
    
    return tau*gamma_r1_tau(pi, tau)*R*T

def cp_r1(p,T):
    """
    Specific isobaric heat capacity for region 1 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    T:      double  temperature (K).

    Returns
    -------
    double  specific isobaric heat capacity (J/kg/K)

    """
    pi = p/p_star_r1
    tau = T_star_r1/T
    
    return -1*(tau**2)*gamma_r1_tautau(pi, tau)*R

def cv_r1(p,T):
    """
    Specific isochoric heat capacity for region 1 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    T:      double  temperature (K).

    Returns
    -------
    double  specific isochoric heat capacity (J/kg/K)

    """
    pi = p/p_star_r1
    tau = T_star_r1/T
    
    return (-1*(tau**2)*gamma_r1_tautau(pi, tau)+(gamma_r1_pi(pi, tau)-tau*gamma_r1_pitau(pi, tau))**2/gamma_r1_pipi(pi, tau))*R

def w_r1(p,T):
    """
    Speed of sound for region 1 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    T:      double  temperature (K).

    Returns
    -------
    double  Speed of sound (m/s)

    """
    pi = p/p_star_r1
    tau = T_star_r1/T
    
    return (gamma_r1_pi(pi, tau)**2/((gamma_r1_pi(pi, tau)-tau*gamma_r1_pitau(pi, tau))**2/(tau**2*gamma_r1_tautau(pi, tau))-gamma_r1_pipi(pi, tau))*R*T)**.5

def t_ph_r1(p,h):
    """
    Backwards equation for 
    Temperature from p,h region 1 IAPWS-IF97
    Chapter 5, eq. 11
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    h:      double  enthalpy (J/kg).

    Returns
    -------
    double  Temperature (K)

    """
    
    #Reducing quantaties
    T_star = 1              #K
    p_star = 1 * 10**6      #Pa
    h_star = 2500 * 10**3   #J/kg
    
    pi = p/p_star
    eta = h/h_star
    
    # Table 6, Numerical values of the coefficients and exponents
    #             1, 2, 3 ,4 ,5 , 6 , 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
    I = np.array([0, 0, 0, 0, 0,  0,  1, 1, 1, 1,  1,  1,  1,  2,  2,  3,  3,  4,  5,  6])
    J = np.array([0, 1, 2, 6, 22, 32, 0, 1, 2, 3,  4,  10, 32, 10, 32, 10, 32, 32, 32, 32])
    
    n_array = np.array([-0.23872489924521*10**3,
                        0.40421188637945*10**3,
                        0.11349746881718*10**3,
                        -0.58457616048039*10,
                        -0.15285482413140*10**-3,
                        -0.10866707695377*10**-5,
                        -0.13391744872602*10**2,
                        0.43211039183559*10**2,
                        -0.54010067170506*10**2,
                        0.30535892203916*10**2,
                        -0.65964749423638*10,
                        0.93965400878363*10**-2,
                        0.11573647505340*10**-6,
                        -0.25858641282073*10**-4,
                        -0.40644363084799*10**-8,
                        0.66456186191635*10**-7,
                        0.80670734103027*10**-10,
                        -0.93477771213947*10**-12,
                        0.58265442020601*10**-14,
                        -0.15020185953503*10**-16])
    
    # backward equation, eq 11
    theta = 0
    for x,y,z in zip(I,J,n_array):
        theta = theta + z*pi**x*(eta+1)**y

    T = theta*T_star
    
    return T

def t_ps_r1(p,s):
    """
    Backwards equation for 
    Temperature from p,s region 1 IAPWS-IF97
    Chapter 5, eq. 13
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    s:      double  entropy (J/kg/K).

    Returns
    -------
    double  Temperature (K)

    """
    
    #Reducing quantaties
    T_star = 1              #K
    p_star = 1 * 10**6      #Pa
    s_star = 1*10**3        #J/kg/K
    
    pi = p/p_star
    sigma = s/s_star
    
    # Table 6, Numerical values of the coefficients and exponents
    #             1, 2, 3 ,4 ,5 , 6 , 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
    I = np.array([0, 0, 0, 0, 0,  0,  1, 1, 1, 1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  4])
    J = np.array([0, 1, 2, 3, 11, 31, 0, 1, 2, 3,  12, 31, 0,  1,  2,  9,  31, 10, 32, 32])
    
    n_array = np.array([0.17478268058307*10**3,
                        0.34806930892873*10**2,
                        0.65292584978455*10,
                        0.33039981775489,
                        -0.19281382923196*10**-6,
                        -0.24909197244573*10**-22,
                        -0.26107636489332,
                        0.22592965981586,
                        -0.64256463395226*10**-1,
                        0.78876289270526*10**-2,
                        0.35672110607366*10**-9,
                        0.17332496994895*10**-23,
                        0.56608900654837*10**-3,
                        -0.32635483139717*10**-3,
                        0.44778286690632*10**-4,
                         -0.51322156908507*10**-9,
                         -0.42522657042207*10**-25,
                         0.26400441360689*10**-12,
                         0.78124600459723*10**-28,
                         -0.30732199903668*10**-30])

    
    # backward equation, eq 11
    theta = 0
    for x,y,z in zip(I,J,n_array):
        theta = theta + z*pi**x*(sigma+2)**y

    T = theta*T_star
    
    return T


"""
    Region 2
    IAPWS-IF97, Chapter 6

    Validity range:
    6.1 Basic eqs,
    273.15 K <= T <= 623.15 K  and 0 < p <ps(T) (eq.30) 
    623.15 K <= T <= 863.15 K  and 0 < p <p(T) (eq.5) 
    863.15 K <= T <= 1073.15 K  and 0 < p < 100 MPa
    
    6.2 Metastable region, 
    Eq. 18 is valid in the metastable-vapor region from the saturated-vapor 
    line to the 5% equilibrium moisture line (determined from the equilibrium 
    h_prim and h_bis values calculated at the given pressure) for pressures 
    from the triple-point pressure, see Eq. 9, up to 10 MPa
"""
#Reducing quantaties
p_star_r2 = 1*10**6          #Pa
T_star_r2 = 540              #K

# Chapter 6.1 basic equation
# Gibbs free energy, ideal part region 2
# Table 10, Numerical values of the coefficients and exponents
jn_array_ideal_r2 = np.array([[0, 1, -5, -4, -3, -2, -1, 2, 3],
                              [-0.96927686500217*10,
                     0.10086655968018*10**2,
                    -0.56087911283020*10**-2,
                     0.71452738081455*10**-1,
                    -0.40710498223928,
                     0.14240819171444*10,
                    -0.43839511319450*10,
                    -0.28408632460772,
                     0.21268463753307*10**-1]])

def gamma_ideal_r2(pi,tau,jn = jn_array_ideal_r2):
    """
    Gibs free energy eq. 16, g(p,T), ideal part
    Region 2 of IAPWS-IF97
    Chapter 6, Table 13
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    jn_array:        np-array of int     coefficient J, n of table 10

    Returns
    -------
    gamma :   double    g(p,T), Gibbs free energy
    """
    gamma_ideal = np.log(pi)
    for j,n in zip(jn[0,:], jn[1,:]):
        gamma_ideal = gamma_ideal + n*tau**j
    
    return gamma_ideal

def gamma_ideal_r2_pi(pi,tau, jn = jn_array_ideal_r2):
    """
    Gibs free energy eq. 16, dg(p,T)/dpi, ideal part
    Region 2 of IAPWS-IF97
    Chapter 6, Table 13
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    jn_array:        np-array of int     coefficient J, n of table 10
    
    Returns
    -------
    gamma :   double    dg(p,T)/dpi, Gibbs free energy
    """
    
    return 1/pi
        
def gamma_ideal_r2_pipi(pi,tau, jn = jn_array_ideal_r2):
    """
    Gibs free energy eq. 16, d^2g(p,T)/dpi^2, ideal part
    Region 2 of IAPWS-IF97
    Chapter 6, Table 13
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    jn_array:        np-array of int     coefficient J, n of table 10

    Returns
    -------
    gamma :   double    d^2g(p,T)/dpi^2, Gibbs free energy
    """
    
    return -1/(pi**2)

def gamma_ideal_r2_tau(pi, tau, jn = jn_array_ideal_r2):
    """
    Gibs free energy eq. 16, dg(p,T)/dtau, ideal part
    Region 2 of IAPWS-IF97
    Chapter 6, Table 13
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    jn_array:        np-array of int     coefficient J, n of table 10

    Returns
    -------
    gamma :   double    dg(p,T)/dtau, Gibbs free energy
    """
   
    gamma_ideal = 0
    for j,n in zip(jn[0,:], jn[1,:]):
        gamma_ideal = gamma_ideal + n*j*tau**(j-1)
    
    return gamma_ideal

def gamma_ideal_r2_tautau(pi, tau, jn = jn_array_ideal_r2):
    """
    Gibs free energy eq. 16, d^2g(p,T)/dtau^2, ideal part
    Region 2 of IAPWS-IF97
    Chapter 6, Table 13
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    jn_array:        np-array of int     coefficient J, n of table 10

    Returns
    -------
    gamma :   double    d^2g(p,T)/dtau^2, Gibbs free energy
    """
   
    gamma_ideal = 0
    for j,n in zip(jn[0,:], jn[1,:]):
        gamma_ideal = gamma_ideal + n*j*(j-1)*tau**(j-2)
    
    return gamma_ideal

def gamma_ideal_r2_pitau(pi, tau, jn = jn_array_ideal_r2):
    """
    Gibs free energy eq. 16, d^2g(p,T)/dpi/dtau, ideal part
    Region 2 of IAPWS-IF97
    Chapter 6, Table 13
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)

    Returns
    -------
    gamma :   double    d^2g(p,T)/dpi/dtau, Gibbs free energy
    """
    
    return 0

#Gibbs free energy, residual part region 2
# Table 11, Numerical values of the coefficients and exponents
ijn_array_residual_r2 = np.array([[1, 1, 1, 1, 1, 2, 2, 2,  2,  2,  3,  3,  3,  3,  3,  4,  4,  4,  5,  6,  6,  6,  7,  7,  7,  8,  8,  9,  10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24],
                              [0, 1, 2, 3, 6, 1, 2, 4,  7,  36, 0,  1,  3,  6,  35, 1,  2,  3,  7,  3,  16, 35, 0,  11, 25, 8,  36, 13, 4,  10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58],
                              [-0.17731742473213*10**-2,
                               -0.17834862292358*10**-1,
                               -0.45996013696365*10**-1,   
                               -0.57581259083432*10**-1,
                               -0.50325278727930*10**-1,
                               -0.33032641670203*10**-4,
                               -0.18948987516315*10**-3,
                               -0.39392777243355*10**-2,
                               -0.43797295650573*10**-1,
                               -0.26674547914087*10**-4,
                                0.20481737692309*10**-7,
                                0.43870667284435*10**-6,
                               -0.32277677238570*10**-4,
                               -0.15033924542148*10**-2,
                               -0.40668253562649*10**-1,
                               -0.78847309559367*10**-9,
                                0.12790717852285*10**-7,
                                0.48225372718507*10**-6,
                                0.22922076337661*10**-5,
                               -0.16714766451061*10**-10,
                               -0.21171472321355*10**-2,
                               -0.23895741934104*10**2,
                               -0.59059564324270*10**-17,
                               -0.12621808899101*10**-5,
                               -0.38946842435739*10**-1,
                                0.11256211360459*10**-10,
                               -0.82311340897998*10,
                                0.19809712802088*10**-7,
                                0.10406965210174*10**-18,
                               -0.10234747095929*10**-12,
                               -0.10018179379511*10**-8,
                               -0.80882908646985*10**-10,
                                0.10693031879409,
                               -0.33662250574171,
                                0.89185845355421*10**-24,
                                0.30629316876232*10**-12,
                               -0.42002467698208*10**-5,
                               -0.59056029685639*10**-25,
                                0.37826947613457*10**-5,
                               -0.12768608934681*10**-14,
                                0.73087610595061*10**-28,
                                0.55414715350778*10**-16,
                               -0.94369707241210*10**-6]])

def gamma_residual_r2(pi, tau, ijn = ijn_array_residual_r2):
    """
    Gibs free energy eq. 17, g(p,T), residual part
    Region 2 of IAPWS-IF97
    Chapter 6, Table 14
    
    Parameters
    ----------
    pi:             double              reduced pressure (-)
    tau:            double              reduced temperature (-)    
    ijn_table_r2    np array            I,J,n coefficients of table 11

    Returns
    -------
    gamma :   double    g(p,T), Gibbs free energy
    """
    gamma_res = 0
    for i,j,n in zip(ijn[0,:], ijn[1,:], ijn[2,:]):
        gamma_res = gamma_res + n*pi**i*(tau-.5)**j
    
    return gamma_res

def gamma_residual_r2_pi(pi, tau, ijn = ijn_array_residual_r2):
    """
    Gibs free energy eq. 17, dg(p,T)/dpi, residual part
    Region 2 of IAPWS-IF97
    Chapter 6, Table 14
    
    Parameters
    ----------
    pi:             double              reduced pressure (-)
    tau:            double              reduced temperature (-)
    ijn_table_r2    np array            I,J,n coefficients of table 11

    Returns
    -------
    gamma :   double    dg(p,T)/dpi, Gibbs free energy
    """
    
    gamma_res = 0
    for i,j,n in zip(ijn[0,:], ijn[1,:], ijn[2,:]):
        gamma_res = gamma_res + n*i*pi**(i-1)*(tau-.5)**j
    
    return gamma_res
   
def gamma_residual_r2_pipi(pi, tau, ijn = ijn_array_residual_r2):
    """
    Gibs free energy eq. 17, d^2g(p,T)/dpi^2, residual part
    Region 2 of IAPWS-IF97
    Chapter 6, Table 14
    
    Parameters
    ----------
    pi:             double              reduced pressure (-)
    tau:            double              reduced temperature (-)
    ijn_table_r2    np array            I,J,n coefficients of table 11

    Returns
    -------
    gamma :   double    d^2g(p,T)/dpi^2, Gibbs free energy
    """
    
    gamma_res = 0
    for i,j,n in zip(ijn[0,:], ijn[1,:], ijn[2,:]):
        gamma_res = gamma_res + n*i*(i-1)*pi**(i-2)*(tau-.5)**j
    
    return gamma_res

def gamma_residual_r2_tau(pi, tau, ijn = ijn_array_residual_r2):
    """
    Gibs free energy eq. 17, dg(p,T)/dtau, residual part
    Region 2 of IAPWS-IF97
    Chapter 6, Table 14
    
    Parameters
    ----------
    pi:             double              reduced pressure (-)
    tau:            double              reduced temperature (-)
    ijn_table_r2    np array            I,J,n coefficients of table 11

    Returns
    -------
    gamma :   double    dg(p,T)/dtau, Gibbs free energy
    """
    
    gamma_res = 0
    for i,j,n in zip(ijn[0,:], ijn[1,:], ijn[2,:]):
        gamma_res = gamma_res + n*pi**i*j*(tau-.5)**(j-1)
    
    return gamma_res
   
def gamma_residual_r2_tautau(pi, tau, ijn = ijn_array_residual_r2):
    """
    Gibs free energy eq. 17, d^2g(p,T)/dtau^2, residual part
    Region 2 of IAPWS-IF97
    Chapter 6, Table 14
    
    Parameters
    ----------
    pi:             double              reduced pressure (-)
    tau:            double              reduced temperature (-)
    ijn_table_r2    np array            I,J,n coefficients of table 11

    Returns
    -------
    gamma :   double    d^2g(p,T)/dtau^2, Gibbs free energy
    """
        
    gamma_res = 0
    for i,j,n in zip(ijn[0,:], ijn[1,:], ijn[2,:]):
        gamma_res = gamma_res + n*pi**i*j*(j-1)*(tau-.5)**(j-2)
    
    return gamma_res
   
def gamma_residual_r2_pitau(pi, tau, ijn = ijn_array_residual_r2):
    """
    Gibs free energy eq. 17, d^2g(p,T)/dpi/dtau, residual part
    Region 2 of IAPWS-IF97
    Chapter 6, Table 14
    
    Parameters
    ----------
    pi:             double              reduced pressure (-)
    tau:            double              reduced temperature (-)
    ijn_table_r2    np array            I,J,n coefficients of table 11

    Returns
    -------
    gamma :   double    d^2g(p,T)/dpi/dtau, Gibbs free energy
    """
        
    gamma_res = 0
    for i,j,n in zip(ijn[0,:], ijn[1,:], ijn[2,:]):
        gamma_res = gamma_res + n*i*pi**(i-1)*j*(tau-.5)**(j-1)
    
    return gamma_res
   
# Relations of thermodynamic properties
# Table 12 Chapter 6.1
   
def v_r2(p,T):
    """
    Specific volume for region 2 IAPWS-IF97
    Chapter 6, Table 12
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).

    Returns
    -------
    double  specific volume (m^3/kg)

    """
    pi = p/p_star_r2
    tau = T_star_r2/T
   
    return pi*(gamma_ideal_r2_pi(pi, tau) + gamma_residual_r2_pi(pi, tau))*R*T/p

def u_r2(p,T):
    """
    Specific internal energy for region 2 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific internal energy (J/kg)

    """
    pi = p/p_star_r2
    tau = T_star_r2/T
    
    return (tau*(gamma_ideal_r2_tau(pi, tau)+gamma_residual_r2_tau(pi, tau))-pi*(gamma_ideal_r2_pi(pi, tau)+gamma_residual_r2_pi(pi, tau)))*R*T

def s_r2(p,T):
    """
    Specific entropy for region 2 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific entropy (J/kg/K)

    """
    pi = p/p_star_r2
    tau = T_star_r2/T

    return (tau*(gamma_ideal_r2_tau(pi, tau)+gamma_residual_r2_tau(pi, tau))-(gamma_ideal_r2(pi, tau)+gamma_residual_r2(pi, tau)))*R   

def h_r2(p,T):
    """
    Specific enthalpy for region 2 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific enthalpy (J/kg)

    """
    pi = p/p_star_r2
    tau = T_star_r2/T
    
    return tau*(gamma_ideal_r2_tau(pi, tau)+gamma_residual_r2_tau(pi, tau))*R*T

def cp_r2(p,T):
    """
    Specific isobaric heat capacity for region 2 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific isobaric heat capacity (J/kg/K)

    """
    pi = p/p_star_r2
    tau = T_star_r2/T
    
    return -1*(tau**2)*(gamma_ideal_r2_tautau(pi, tau)+gamma_residual_r2_tautau(pi, tau))*R

def cv_r2(p,T):
    """
    Specific isochoric heat capacity for region 2 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific isochoric heat capacity (J/kg/K)

    """
    pi = p/p_star_r2
    tau = T_star_r2/T
    
    return (-1*(tau**2)*(gamma_ideal_r2_tautau(pi, tau)+gamma_residual_r2_tautau(pi, tau))-((1+pi*gamma_residual_r2_pi(pi, tau)-tau*pi*gamma_residual_r2_pitau(pi, tau))**2/(1-pi**2*gamma_residual_r2_pipi(pi, tau))))*R

def w_r2(p,T):
    """
    Speed of sound for region 2 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  Speed of sound (m/s)

    """
    pi = p/p_star_r2
    tau = T_star_r2/T
    
    return (((1+2*pi*gamma_residual_r2_pi(pi, tau)+pi**2*gamma_residual_r2_pi(pi, tau)**2)/((1-pi**2*gamma_residual_r2_pipi(pi, tau))+((1+pi*gamma_residual_r2_pi(pi, tau)-tau*pi*gamma_residual_r2_pitau(pi, tau))**2/(tau**2*(gamma_ideal_r2_tautau(pi, tau)+gamma_residual_r2_tautau(pi, tau))))))*R*T)**.5


# Chapter 6.2 Supplementary Equation for the Metastable-Vapor Region
# Gibbs free energy, ideal part region 2
# Table 10, Numerical values of the coefficients and exponents
jn_array_metastable_ideal_r2 = jn_array_ideal_r2.copy()
jn_array_metastable_ideal_r2[1,0] = -0.96937268393049*10
jn_array_metastable_ideal_r2[1,1] = 0.10087275970006*10**2

# Gibbs free energy, ideal part region 2
# Table 16, Numerical values of the coefficients and exponents
ijn_array_residual_metastable_r2 = np.array([[   1, 1, 1, 1,  2, 2, 2,  3, 3,  4, 4,  5,  5],
                                    [   0, 2, 5, 11, 1, 7, 16, 4, 16, 7, 10, 9,  10],
                                    [-0.73362260186506*10**-2,
                         -0.88223831943146*10**-1,
                         -0.72334555213245*10**-1,
                         -0.40813178534455*10**-2,
                          0.20097803380207*10**-2,
                         -0.53045921898642*10**-1,
                         -0.76190409086970*10**-2,
                         -0.63498037657313*10**-2,
                         -0.86043093028588*10**-1,
                         0.75321581522770*10**-2,
                         -0.79238375446139*10**-2,
                         -0.22888160778447*10**-3,
                         -0.26456501482810*10**-2]])

# Metastable equations same as ordinary region 2 but with I,J and n as above.
def v_r2_metastable(p,T,jn_ideal_metastable = jn_array_metastable_ideal_r2,ijn_residual_metastable = ijn_array_residual_metastable_r2):
    """
    Specific volume for region 2 IAPWS-IF97
    Chapter 6, Table 12
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).

    Returns
    -------
    double  specific volume (m^3/kg)

    """
    pi = p/p_star_r2
    tau = T_star_r2/T
   
    return pi*(gamma_ideal_r2_pi(pi, tau,jn_ideal_metastable) + gamma_residual_r2_pi(pi, tau,ijn_residual_metastable))*R*T/p

def u_r2_metastable(p,T,jn_ideal_metastable = jn_array_metastable_ideal_r2,ijn_residual_metastable = ijn_array_residual_metastable_r2):
    """
    Specific internal energy for region 2 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific internal energy (J/kg)

    """
    pi = p/p_star_r2
    tau = T_star_r2/T
    
    return (tau*(gamma_ideal_r2_tau(pi, tau,jn_ideal_metastable)+gamma_residual_r2_tau(pi, tau,ijn_residual_metastable))-pi*(gamma_ideal_r2_pi(pi, tau,jn_ideal_metastable)+gamma_residual_r2_pi(pi, tau,ijn_residual_metastable)))*R*T

def s_r2_metastable(p,T,jn_ideal_metastable = jn_array_metastable_ideal_r2,ijn_residual_metastable = ijn_array_residual_metastable_r2):
    """
    Specific entropy for region 2 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific entropy (J/kg/K)

    """
    pi = p/p_star_r2
    tau = T_star_r2/T

    return (tau*(gamma_ideal_r2_tau(pi, tau,jn_ideal_metastable)+gamma_residual_r2_tau(pi, tau,ijn_residual_metastable))-(gamma_ideal_r2(pi, tau,jn_ideal_metastable)+gamma_residual_r2(pi, tau,ijn_residual_metastable)))*R   

def h_r2_metastable(p,T,jn_ideal_metastable = jn_array_metastable_ideal_r2,ijn_residual_metastable = ijn_array_residual_metastable_r2):
    """
    Specific enthalpy for region 2 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific enthalpy (J/kg)

    """
    pi = p/p_star_r2
    tau = T_star_r2/T
    
    return tau*(gamma_ideal_r2_tau(pi, tau,jn_ideal_metastable)+gamma_residual_r2_tau(pi, tau,ijn_residual_metastable))*R*T

def cp_r2_metastable(p,T,jn_ideal_metastable = jn_array_metastable_ideal_r2,ijn_residual_metastable = ijn_array_residual_metastable_r2):
    """
    Specific isobaric heat capacity for region 2 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific isobaric heat capacity (J/kg/K)

    """
    pi = p/p_star_r2
    tau = T_star_r2/T
    
    return -1*(tau**2)*(gamma_ideal_r2_tautau(pi, tau,jn_ideal_metastable)+gamma_residual_r2_tautau(pi, tau,ijn_residual_metastable))*R

def cv_r2_metastable(p,T,jn_ideal_metastable = jn_array_metastable_ideal_r2,ijn_residual_metastable = ijn_array_residual_metastable_r2):
    """
    Specific isochoric heat capacity for region 2 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific isochoric heat capacity (J/kg/K)

    """
    pi = p/p_star_r2
    tau = T_star_r2/T
    
    return (-1*(tau**2)*(gamma_ideal_r2_tautau(pi, tau,jn_ideal_metastable)+gamma_residual_r2_tautau(pi, tau,ijn_residual_metastable))-((1+pi*gamma_residual_r2_pi(pi, tau,ijn_residual_metastable)-tau*pi*gamma_residual_r2_pitau(pi, tau,ijn_residual_metastable))**2/(1-pi**2*gamma_residual_r2_pipi(pi, tau,ijn_residual_metastable))))*R

def w_r2_metastable(p,T,jn_ideal_metastable = jn_array_metastable_ideal_r2,ijn_residual_metastable = ijn_array_residual_metastable_r2):
    """
    Speed of sound for region 2 IAPWS-IF97
    Chapter 5, Table 4
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  Speed of sound (m/s)

    """
    pi = p/p_star_r2
    tau = T_star_r2/T
    
    return (((1+2*pi*gamma_residual_r2_pi(pi, tau,ijn_residual_metastable)+pi**2*gamma_residual_r2_pi(pi, tau,ijn_residual_metastable)**2)/((1-pi**2*gamma_residual_r2_pipi(pi, tau,ijn_residual_metastable))+((1+pi*gamma_residual_r2_pi(pi, tau,ijn_residual_metastable)-tau*pi*gamma_residual_r2_pitau(pi, tau,ijn_residual_metastable))**2/(tau**2*(gamma_ideal_r2_tautau(pi, tau,jn_ideal_metastable)+gamma_residual_r2_tautau(pi, tau,ijn_residual_metastable))))))*R*T)**.5

### Backwards equations for region 2, Chapt. 6.3 ###
# Equations (20), pi_b2bc(h) and (21), eta_b2bc(p) give the boundary line 
# between subregions 2b and 2c from the saturation state at T = 554.485 K and 
# ps = 6.546 70 MPa to T = 1019.32 K and p = 100 MPa. 
# Subregion 2a is defined for all pressures below 4 MPa

def p_b2bc(h):
    """
    Backwards equation, boundary calculation betwwen region 2b/2c
    Approximatios isentropic line, s = 5.5 kJ/kg/K
    Region 2 IAPWS-IF97
    Chapter 6.3, eq. 20
    
    Parameters
    ----------
    h:      double  enthalpy (J/kg).

    Returns
    -------
    double  Pressure (Pa)

    """
    p_star = 1*10**6 # 1 MPa
    h_star = 1*10**3 # 1 kJ/kg  
    
    # Numerical coefficients from table 19
    n = np.array([0.90584278514723*10**3, -0.67955786399241, 0.12809002730136*10**-3])
    
    eta = h/h_star
    pi = n[0] + n[1] * eta + n[2] * eta**2
    
    return pi * p_star

def h_b2bc(p):
    """
    Backwards equation, boundary calculation betwwen region 2b/2c
    Approximatios isentropic line, s = 5.85 kJ/kg/K
    Region 2 IAPWS-IF97
    Chapter 6.3, eq. 21
    
    Parameters
    ----------
    p:      double  Pressure (Pa).

    Returns
    -------
    double  Enthalpy (J/kg)

    """
    
    p_star = 1*10**6 # 1 MPa
    h_star = 1*10**3 # 1 kJ/kg 
    
    # Numerical coefficients from table 19
    n = np.array([0.12809002730136*10**-3, 0.26526571908428*10**4, 0.45257578905948*10])
    
    pi = p/p_star
    eta = n[1] + ((pi-n[2])/(n[0]))**.5
    
    return eta * h_star

def t_ph_r2a(p,h):
    
    # Table 20, Numerical values of the coefficients and exponents for backwards 
    # eq. T_ph region 2a

    #                     1, 2, 3, 4, 5, 6,  7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34 
    I_r2a_backwards_ph = [0, 0, 0, 0, 0, 0,  1, 1, 1, 1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  4,  4,  4,  5,  5,  5,  6,  6,  7]
    J_r2a_backwards_ph = [0, 1, 2, 3, 7, 20, 0, 1, 2, 3,  7,  9,  11, 18, 44, 0,  2,  7,  36, 38, 40, 42, 44, 24, 44, 12, 32, 44, 32, 36, 42, 34, 44, 28]
    n_r2a_backwards_ph = [0.10898952318288*10**4,
                          0.84951654495535*10**3,
                         -0.10781748091826*10**3,
                          0.33153654801263*10**2,
                         -0.74232016790248*10,
                          0.11765048724356*10**2,
                          0.18445749355790*10,
                         -0.41792700549624*10,
                          0.62478196935812*10,
                         -0.17344563108114*10**2,
                         -0.20058176862096*10**3,
                          0.27196065473796*10**3,
                         -0.45511318285818*10**3,
                          0.30919688604755*10**4,
                          0.25226640357872*10**6,
                         -0.61707422868339*10**-2,
                         -0.31078046629583,
                          0.11670873077107*10**2,
                          0.12812798404046*10**9,
                         -0.98554909623276*10**9,
                          0.28224546973002*10**10,
                         -0.35948971410703*10**10,
                          0.17227349913197*10**10,
                         -0.13551334240775*10**5,
                          0.12848734664650*10**8,
                          0.13865724283226*10,
                          0.23598832556514*10**6,
                         -0.13105236545054*10**8,
                          0.73999835474766*10**4,
                         -0.55196697030060*10**6,
                          0.37154085996233*10**7,
                          0.19127729239660*10**5,
                         -0.41535164835634*10**6,
                         -0.62459855192507*10**2]
    
    p_star = 1 * 10**6 #Pa
    h_star = 2000 * 10**3 #J/kg
    T_star = 1 #K
    
    pi = p/p_star
    eta = h/h_star
    
    # Backwards equation, Eq 22
    teta2a = 0
    for i,j,n in zip(I_r2a_backwards_ph, J_r2a_backwards_ph, n_r2a_backwards_ph):
        teta2a = teta2a + n*pi**i*(eta-2.1)**j
    
    return teta2a*T_star


def t_ph_r2b(p,h):
    
    # Table 21, Numerical values of the coefficients and exponents for backwards 
    # eq. T_ph region 2b

    #                     1, 2, 3, 4,  5,  6,  7,  8,  9, 10,11,12, 13, 14, 15, 16, 17,18,19,20, 21, 22,23,24, 25, 26, 27, 28,29, 30, 31, 32, 33, 34, 35, 36, 37, 38 
    I_r2b_backwards_ph = [0, 0, 0, 0,  0,  0,  0,  0,  1, 1, 1, 1,  1,  1,  1,  1,  2, 2, 2,  2,  3, 3, 3,  3,  4, 4,  4,  4, 4,  4,  5,  5,  5,  6,  7,  7,  9,  9]
    J_r2b_backwards_ph = [0, 1, 2, 12, 18, 24, 28, 40, 0, 2, 6, 12, 18, 24, 28, 40, 2, 8, 18, 40, 1, 2, 12, 24, 2, 12, 18, 24,28, 40, 18, 24, 40, 28, 2,  28, 1,  40]
    n_r2b_backwards_ph = [ 0.14895041079516*10**4,
                           0.74307798314034*10**3,
                          -0.97708318797837*10**2,
                           0.24742464705674*10,
                          -0.63281320016026,
                           0.11385952129658*10,
                          -0.47811863648625,
                           0.85208123431544*10**-2,
                           0.93747147377932,
                           0.33593118604916*10,
                           0.33809355601454*10,
                           0.16844539671904,
                           0.73875745236695,
                          -0.47128737436186,
                           0.15020273139707,
                          -0.21764114219750*10**-2,
                          -0.21810755324761*10**-1,
                          -0.10829784403677,
                          -0.46333324635812*10**-1,
                           0.71280351959551*10**-4,
                           0.11032831789999*10**-3,
                           0.18955248387902*10**-3,
                           0.30891541160537*10**-2,
                           0.13555504554949*10**-2,
                           0.28640237477456*10**-6,
                          -0.10779857357512*10**-4,
                          -0.76462712454814*10**-4,
                           0.14052392818316*10**-4,
                          -0.31083814331434*10**-4,
                          -0.10302738212103*10**-5,
                           0.28217281635040*10**-6,
                           0.12704902271945*10**-5,
                           0.73803353468292*10**-7,
                          -0.11030139238909*10**-7,
                          -0.81456365207833*10**-13,
                          -0.25180545682962*10**-10,
                          -0.17565233969407*10**-17,
                           0.86934156344163*10**-14]
    
    p_star = 1 * 10**6 #Pa
    h_star = 2000 * 10**3 #J/kg
    T_star = 1 #K
    
    pi = p/p_star
    eta = h/h_star
    
    # Backwards equation, Eq 23
    teta2b = 0
    for i,j,n in zip(I_r2b_backwards_ph, J_r2b_backwards_ph, n_r2b_backwards_ph):
        teta2b = teta2b + n*(pi-2)**i*(eta-2.6)**j
    
    return teta2b*T_star


def t_ph_r2c(p,h):
    
    # Table 22, Numerical values of the coefficients and exponents for backwards 
    # eq. T_ph region 2c

    #                     1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,12,13,14,15,16,17,18,19, 20, 21, 22, 23 
    I_r2c_backwards_ph = [-7, -7, -6, -6, -5, -5, -2, -2, -1, -1, 0, 0, 1, 1, 2, 6, 6, 6, 6,  6,  6,  6,  6]
    J_r2c_backwards_ph = [0,  4,  0,  2,  0,  2,  0,  1,  0,  2,  0, 1, 4, 8, 4, 0, 1, 4, 10, 12, 16, 20, 22]
    n_r2c_backwards_ph =[-0.32368398555242*10**13,
                       0.73263350902181*10**13,
                       0.35825089945447*10**12,
                      -0.58340131851590*10**12,
                      -0.10783068217470*10**11,
                       0.20825544563171*10**11,
                       0.61074783564516 *10**6,
                       0.85977722535580*10**6,
                      -0.25745723604170*10**5,
                       0.31081088422714*10**5,
                       0.12082315865936*10**4,
                       0.48219755109255*10**3,
                       0.37966001272486*10,
                      -0.10842984880077*10**2,
                      -0.45364172676660*10**-1,
                       0.14559115658698*10**-12,
                       0.11261597407230*10**-11,
                      -0.17804982240686*10**-10,
                       0.12324579690832*10**-6,
                      -0.11606921130984*10**-5,
                       0.27846367088554*10**-4,
                      -0.59270038474176*10**-3,
                       0.12918582991878*10**-2]
    
    p_star = 1 * 10**6 #Pa
    h_star = 2000 * 10**3 #J/kg
    T_star = 1 #K
    
    pi = p/p_star
    eta = h/h_star
    
    # Backwards equation, Eq 24
    teta2c = 0
    for i,j,n in zip(I_r2c_backwards_ph, J_r2c_backwards_ph, n_r2c_backwards_ph):
        teta2c = teta2c + n*(pi+25)**i*(eta-1.8)**j
    
    return teta2c*T_star

def t_ph_r2(p,h):
    """
    Finds correct subregion for backwards equations 2a-c
    and returns result of applicable eqation.
    
    For pressures below 4 MPa -> region a
    if entropy s >= 5.85 kJ/kg/K -> region b
    s < 5.85 -> region c
    the isoentropic line is calculated based on enthalpy using function p_b2bc 
    returning the corresponding pressure
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    h :     double  Enthalpy (J/kg)

    Returns
    -------
    double Temperature (K)

    """
    
    if p <= 4*10**6:
        T = t_ph_r2a(p,h)
    elif p < p_b2bc(h):
        T = t_ph_r2b(p,h)
    else:
        T = t_ph_r2c(p,h)
    
    return T

def t_ps_r2a(p,s):
    
    # Table 25, Numerical values of the coefficients and exponents for backwards 
    # eq. T_ps region 2a

    #                     1,    2,    3,    4,    5,    6,    7,     8,     9,     10,   11,   12,   13,   14,   15,   16,   17,   18,  19,  20,  21,  22,   23,   24,   25,   26,  27,  28,  29,   30, 31, 32, 33, 34, 35, 36, 37,  38,  39,  40,  41,  42,  43,   44,   45,  46 
    I_r2a_backwards_ps = [-1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.25, -1.25, -1.25, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -.75, -.75, -.5, -.5, -.5, -.5, -.25, -.25, -.25, -.25, .25, .25, .25, .25,  .5, .5, .5, .5, .5, .5, .5, .75, .75, .75, .75, 1.0, 1.0, 1.25, 1.25, 1.5, 1.5]
    J_r2a_backwards_ps = [-24,  -23,  -19,  -13,  -11,  -10,  -19,   -15,   -6,    -26,  -21,  -17,  -16,  -9,   -8,   -15,  -14,  -26, -13, -9,  -7,  -27,  -25,  -11,  -6,   1,   4,   8,   11,   0,  1,  5,  6,  10, 14, 16, 0,   4,   9,   17,  7,   18,  3,    15,   5,   18]
    n_r2a_backwards_ps =[-0.39235983861984*10**6,
                          0.51526573827270*10**6,
                          0.40482443161048*10**5,
                         -0.32193790923902*10**3,
                          0.96961424218694*10**2,
                         -0.22867846371773*10**2,
                         -0.44942914124357*10**6,
                         -0.50118336020166*10**4,
                          0.35684463560015,
                          0.44235335848190*10**5,
                         -0.13673388811708*10**5,
                          0.42163260207864*10**6,
                          0.22516925837475*10**5,
                          0.47442144865646*10**3,
                         -0.14931130797647*10**3,
                         -0.19781126320452*10**6,
                         -0.23554399470760*10**5,
                         -0.19070616302076*10**5,
                          0.55375669883164*10**5,
                          0.38293691437363*10**4,
                         -0.60391860580567*10**3,
                          0.19363102620331*10**4,
                          0.42660643698610*10**4,
                         -0.59780638872718*10**4,
                         -0.70401463926862*10**3,
                          0.33836784107553*10**3,
                          0.20862786635187*10**2,
                          0.33834172656196*10**-1,
                         -0.43124428414893*10**-4,
                          0.16653791356412*10**3,
                         -0.13986292055898*10**3,
                         -0.78849547999872,
                          0.72132411753872*10**-1,
                         -0.59754839398283*10**-2,
                         -0.12141358953904*10**-4,
                          0.23227096733871*10**-6,
                         -0.10538463566194*10**2,
                          0.20718925496502*10,
                         -0.72193155260427*10**-1,
                          0.20749887081120*10**-6,
                         -0.18340657911379*10**-1,
                          0.29036272348696*10**-6,
                          0.21037527893619,
                          0.25681239729999*10**-3,
                         -0.12799002933781*10**-1,
                         -0.82198102652018*10**-5]
    
    p_star = 1 * 10**6 #Pa
    T_star = 1 #K
    s_star = 2*10**3 #J/kg/K
    
    pi = p/p_star
    sigma = s/s_star
    
    # Backwards equation, Eq 25
    teta2a = 0
    for i,j,n in zip(I_r2a_backwards_ps, J_r2a_backwards_ps, n_r2a_backwards_ps):
        teta2a += n*(pi)**i*(sigma-2)**j
    
    return teta2a*T_star

def t_ps_r2b(p,s):
    
    # Table 26, Numerical values of the coefficients and exponents for backwards 
    # eq. T_ps region 2b

    #                     1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44 
    I_r2b_backwards_ps = [-6, -6, -5, -5, -4, -4, -4, -3, -3, -3, -3, -2, -2, -2, -2, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5]
    J_r2b_backwards_ps = [0,  11, 0,  11, 0,  1,  11, 0,  1,  11, 12, 0,  1,  6,  10, 0,  1,  5,  8,  9,  0, 1, 2, 4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 0, 1, 5, 0, 1, 3, 0, 1, 0, 1, 2]
    n_r2b_backwards_ps =[ 0.31687665083497*10**6,
                          0.20864175881858*10**2,
                         -0.39859399803599*10**6,
                         -0.21816058518877*10**2,
                          0.22369785194242*10**6,
                         -0.27841703445817*10**4,
                          0.99207436071480*10,
                         -0.75197512299157*10**5,
                          0.29708605951158*10**4,
                         -0.34406878548526*10,
                          0.38815564249115,
                          0.17511295085750*10**5,
                         -0.14237112854449*10**4,
                          0.10943803364167*10,
                          0.89971619308495,
                         -0.33759740098958*10**4,
                          0.47162885818355*10**3,
                         -0.19188241993679*10,
                          0.41078580492196,
                         -0.33465378172097,
                          0.13870034777505*10**4,
                         -0.40663326195838*10**3,
                          0.41727347159610*10**2,
                          0.21932549434532*10,
                         -0.10320050009077*10,
                          0.35882943516703,
                          0.52511453726066*10**-2,
                          0.12838916450705*10**2,
                         -0.28642437219381*10,
                          0.56912683664855,
                         -0.99962954584931*10**-1,
                         -0.32632037778459*10**-2,
                          0.23320922576723*10**-3,
                         -0.15334809857450,
                          0.29072288239902*10**-1,
                          0.37534702741167*10**-3,
                          0.17296691702411*10**-2,
                         -0.38556050844504*10**-3,
                         -0.35017712292608*10**-4,
                         -0.14566393631492*10**-4,
                          0.56420857267269*10**-5,
                          0.41286150074605*10**-7,
                         -0.20684671118824*10**-7,
                          0.16409393674725*10**-8]
    
    p_star = 1 * 10**6 #Pa
    T_star = 1 #K
    s_star = 0.7853*10**3 #J/kg/K
    
    pi = p/p_star
    sigma = s/s_star
    
    # Backwards equation, Eq 26
    teta2b = 0
    for i,j,n in zip(I_r2b_backwards_ps, J_r2b_backwards_ps, n_r2b_backwards_ps):
        teta2b += n*(pi)**i*(10-sigma)**j
    
    return teta2b*T_star

def t_ps_r2c(p,s):
    
    # Table 27, Numerical values of the coefficients and exponents for backwards 
    # eq. T_ps region 2c

    #                     1,  2,  3,  4, 5, 6, 7, 8, 9, 10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30 
    I_r2c_backwards_ps = [-2, -2, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7]
    J_r2c_backwards_ps = [ 0,  1,  0, 0, 1, 2, 3, 0, 1, 3, 4, 0, 1, 2, 0, 1, 5, 0, 1, 4, 0, 1, 2, 0, 1, 0, 1, 3, 4, 5]
    n_r2c_backwards_ps =[ 0.90968501005365*10**3,
                          0.24045667088420*10**4,
                         -0.59162326387130*10**3,
                          0.54145404128074*10**3,
                         -0.27098308411192*10**3,
                          0.97976525097926*10**3,
                         -0.46966772959435*10**3,
                          0.14399274604723*10**2,
                         -0.19104204230429*10**2,
                          0.53299167111971*10,
                         -0.21252975375934*10**2,
                         -0.31147334413760,
                          0.60334840894623,
                         -0.42764839702509*10**-1,
                          0.58185597255259*10**-2,
                         -0.14597008284753*10**-1,
                          0.56631175631027*10**-2,
                         -0.76155864584577*10**-4,
                          0.22440342919332*10**-3,
                         -0.12561095013413*10**-4,
                          0.63323132660934*10**-6,
                         -0.20541989675375*10**-5,
                          0.36405370390082*10**-7,
                         -0.29759897789215*10**-8,
                          0.10136618529763*10**-7,
                          0.59925719692351*10**-11,
                         -0.20677870105164*10**-10,
                         -0.20874278181886*10**-10,
                          0.10162166825089*10**-9,
                         -0.16429828281347*10**-9]
    
    p_star = 1 * 10**6 #Pa
    T_star = 1 #K
    s_star = 2.9251*10**3 #J/kg/K
    
    pi = p/p_star
    sigma = s/s_star
    
    # Backwards equation, Eq 27
    teta2b = 0
    for i,j,n in zip(I_r2c_backwards_ps, J_r2c_backwards_ps, n_r2c_backwards_ps):
        teta2b += n*(pi)**i*(2-sigma)**j
    
    return teta2b*T_star

def t_ps_r2(p,s):
    """
    Finds correct subregion for backwards equations 2a-c
    and returns result of applicable eqation
    
    For pressures below 4 MPa -> region a
    if entropy s >= 5.85 kJ/kg/K -> region b
    s < 5.85 -> region c
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    s :     double  Entropy (J/kg/k)

    Returns
    -------
    double Temperature (K)

    """
    
    if p <= 4*10**6:
        T = t_ps_r2a(p,s)
    elif s >= 5.85*10**3:
        T = t_ps_r2b(p,s)
    else:
        T = t_ps_r2c(p,s)
    
    return T

"""
    Region 3
    IAPWS-IF97, Chapter 7

    Validity range:
    623.15 K <= T <= (eq. 6) K  and (eq. 5) <= p <= 100 MPa 

"""

#Reducing quantaties
rho_star_r3 = rho_c        #kg/m^3
T_star_r3 = T_crit            #K

#Helmholtz free energy, residual part region 3
# Table 30, Numerical values of the coefficients and exponents
#                1, 2, 3, 4, 5,  6,  7,  8, 9, 10, 11, 12,13,14,15,16, 17, 18,19,20,21, 22, 23,24,25,26, 27,28,29, 30,31,32,33, 34,35,36, 37,38,39,40, 41, 42, 43
n_r3_1 = 0.10658070028513*10
ijn_array_r3 = np.array([[0, 0, 0, 0, 0,  0,  0,  1, 1, 1,  1,  2, 2, 2, 2, 2,  2,  3, 3, 3, 3,  3,  4, 4, 4, 4,  5, 5, 5,  6, 6, 6, 7,  8, 9, 9,  10,10,11],
                         [0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26,2,  26,2, 26, 0, 1, 26],
                         [      -0.15732845290239*10**2,
                        0.20944396974307*10**2,
                       -0.76867707878716*10,
                        0.26185947787954*10,
                       -0.28080781148620*10,
                        0.12053369696517*10,
                       -0.84566812812502*10**-2,
                       -0.12654315477714*10,
                       -0.11524407806681*10,
                        0.88521043984318,
                       -0.64207765181607,
                        0.38493460186671,
                       -0.85214708824206,
                        0.48972281541877*10,
                       -0.30502617256965*10,
                        0.39420536879154*10**-1,
                        0.12558408424308,
                       -0.27999329698710,
                        0.13899799569460*10,
                       -0.20189915023570*10,
                       -0.82147637173963*10**-2,
                       -0.47596035734923,
                        0.43984074473500*10**-1,
                       -0.44476435428739,
                        0.90572070719733,
                        0.70522450087967,
                        0.10770512626332,
                       -0.32913623258954,
                       -0.50871062041158,
                       -0.22175400873096*10**-1,
                        0.94260751665092*10**-1,
                        0.16436278447961,
                       -0.13503372241348*10**-1,
                       -0.14834345352472*10**-1,
                        0.57922953628084*10**-3,
                        0.32308904703711*10**-2,
                        0.80964802996215*10**-4,
                       -0.16557679795037*10**-3,
                       -0.44923899061815*10**-4]])

def phi_r3(delta, tau, ijn = ijn_array_r3,n_r3_1 = n_r3_1):
    """
    Helmholtz free energy eq. 28, f(p,T)/RT
    Region 3 of IAPWS-IF97
    Chapter 7, Table 32
    
    Parameters
    ----------
    rho:            double              reduced density (-)
    tau:            double              reduced temperature (-)
    ijn_table_r3    np array            I,J,n coefficients of table 30
    n_r3_1          double              n(i=1) coefficient of table 30

    Returns
    -------
    phi :   double    f(p,T)/RT, Helmholtz free energy
    """  
    phi = n_r3_1*np.log(delta)
    for i,j,n in zip(ijn[0,:],ijn[1,:],ijn[2,:]):
        phi += n*delta**i*tau**j

    return phi

def phi_r3_delta(delta, tau, ijn = ijn_array_r3,n_r3_1 = n_r3_1):
    """
    Helmholtz free energy eq. 28, d/d delta (f(p,T)/RT)
    Region 3 of IAPWS-IF97
    Chapter 7, Table 32
    
    Parameters
    ----------
    rho:             double              reduced density (-)
    tau:            double              reduced temperature (-)
    ijn_table_r3    np array            I,J,n coefficients of table 30
    n_r3_1          double              n(i=1) coefficient of table 30


    Returns
    -------
    phi :   double    d/d delta (f(p,T)/RT), Helmholtz free energy
    """
    
    phi = n_r3_1/delta
    for i,j,n in zip(ijn[0,:],ijn[1,:],ijn[2,:]):
        phi += n*i*delta**(i-1)*tau**j

    return phi


def phi_r3_deltadelta(delta, tau, ijn = ijn_array_r3,n_r3_1 = n_r3_1):
    """
    Helmholtz free energy eq. 28, d^2/d delta^2 (f(p,T)/RT)
    Region 3 of IAPWS-IF97
    Chapter 7, Table 32
    
    Parameters
    ----------
    rho:             double              reduced density (-)
    tau:            double              reduced temperature (-)    
    ijn_table_r3    np array            I,J,n coefficients of table 30
    n_r3_1          double              n(i=1) coefficient of table 30

    Returns
    -------
    phi :   double    d^2/d delta^2 (f(p,T)/RT), Helmholtz free energy
    """
    
    phi = -n_r3_1/(delta**2)
    for i,j,n in zip(ijn[0,:],ijn[1,:],ijn[2,:]):
        phi += n*i*(i-1)*delta**(i-2)*tau**j

    return phi

def phi_r3_tau(delta, tau, ijn = ijn_array_r3,n_r3_1 = n_r3_1):
    """
    Helmholtz free energy eq. 28, d/d delta (f(p,T)/RT)
    Region 3 of IAPWS-IF97
    Chapter 7, Table 32
    
    Parameters
    ----------
    rho:             double              reduced density (-)
    tau:            double              reduced temperature (-)
    ijn_table_r3    np array            I,J,n coefficients of table 30
    n_r3_1          double              n(i=1) coefficient of table 30


    Returns
    -------
    phi :   double    d/d tau (f(p,T)/RT), Helmholtz free energy
    """
    
    phi = 0
    for i,j,n in zip(ijn[0,:],ijn[1,:],ijn[2,:]):
        phi += n*delta**i*j*tau**(j-1)

    return phi

def phi_r3_tautau(delta, tau, ijn = ijn_array_r3,n_r3_1 = n_r3_1):
    """
    Helmholtz free energy eq. 28, d^2/d delta^2 (f(p,T)/RT)
    Region 3 of IAPWS-IF97
    Chapter 7, Table 32
    
    Parameters
    ----------
    rho:             double              reduced density (-)
    tau:            double              reduced temperature (-)    
    ijn_table_r3    np array            I,J,n coefficients of table 30
    n_r3_1          double              n(i=1) coefficient of table 30

    Returns
    -------
    phi :   double    d^2/d tau^2 (f(p,T)/RT), Helmholtz free energy
    """
    
    phi = 0
    for i,j,n in zip(ijn[0,:],ijn[1,:],ijn[2,:]):
        phi += n*delta**i*j*(j-1)*tau**(j-2)

    return phi

def phi_r3_deltatau(delta, tau, ijn = ijn_array_r3,n_r3_1 = n_r3_1):
    """
    Helmholtz free energy eq. 28, d^2/d delta^2 (f(p,T)/RT)
    Region 3 of IAPWS-IF97
    Chapter 7, Table 32
    
    Parameters
    ----------
    rho:             double              reduced density (-)
    tau:            double              reduced temperature (-)    
    ijn_table_r3    np array            I,J,n coefficients of table 30
    n_r3_1          double              n(i=1) coefficient of table 30

    Returns
    -------
    phi :   double    d^2/(d delta d tau) (f(p,T)/RT), Helmholtz free energy
    """
    
    phi = 0
    for i,j,n in zip(ijn[0,:],ijn[1,:],ijn[2,:]):
        phi += n*i*delta**(i-1)*j*tau**(j-1)

    return phi

# Relations of thermodynamic properties
# Table 31 Chapter 7
   
def p_r3(rho,T):
    """
    Pressure for region 3 IAPWS-IF97
    Chapter 7, Table 31
    
    Parameters
    ----------
    rho:            double              density (kg/m^3).
    T:              double              temperature (K).
    ijn_table_r3    np array            I,J,n coefficients of table 30
    n_r3_1          double              n(i=1) coefficient of table 30

    Returns
    -------
    double  pressure (Pa)

    """
    delta = rho/rho_star_r3
    tau = T_star_r3/T
   
    return delta*phi_r3_delta(delta, tau)*rho*R*T

def u_r3(rho,T):
    """
    Specific internal energy for region 3 IAPWS-IF97
    Chapter 7, Table 31
    
    Parameters
    ----------
    rho:            double              density (kg/m^3).
    T:              double              temperature (K).
    ijn_table_r3    np array            I,J,n coefficients of table 30
    n_r3_1          double              n(i=1) coefficient of table 30

    
    Returns
    -------
    double  specific internal energy (J/kg)

    """
    delta = rho/rho_star_r3
    tau = T_star_r3/T
    
    return tau*phi_r3_tau(delta, tau)*R*T

def s_r3(rho,T):
    """
    Specific entropy for region 3 IAPWS-IF97
    Chapter 7, Table 31
    
    Parameters
    ----------
    rho:            double              density (kg/m^3).
    T:              double              temperature (K).
    ijn_table_r3    np array            I,J,n coefficients of table 30
    n_r3_1          double              n(i=1) coefficient of table 30
    
    Returns
    -------
    double  specific entropy (J/kg/K)

    """
    delta = rho/rho_star_r3
    tau = T_star_r3/T

    return (tau*phi_r3_tau(delta, tau)-phi_r3(delta, tau))*R

def h_r3(rho,T):
    """
    Specific enthalpy for region 3 IAPWS-IF97
    Chapter 7, Table 31
    
    Parameters
    ----------
    rho:            double              density (kg/m^3).
    T:              double              temperature (K).
    ijn_table_r3    np array            I,J,n coefficients of table 30
    n_r3_1          double              n(i=1) coefficient of table 30

    
    Returns
    -------
    double  specific enthalpy (J/kg)

    """
    delta = rho/rho_star_r3
    tau = T_star_r3/T
    
    return (tau*phi_r3_tau(delta, tau) + delta*phi_r3_delta(delta, tau))*R*T

def cp_r3(rho,T):
    """
    Specific isobaric heat capacity for region 3 IAPWS-IF97
    Chapter 7, Table 31
    
    Parameters
    ----------
    rho:            double              density (kg/m^3).
    T:              double              temperature (K).
    ijn_table_r3    np array            I,J,n coefficients of table 30
    n_r3_1          double              n(i=1) coefficient of table 30
    
    Returns
    -------
    double  specific isobaric heat capacity (J/kg/K)

    """
    delta = rho/rho_star_r3
    tau = T_star_r3/T
    
    return (-1*tau**2*phi_r3_tautau(delta, tau)+(((delta*phi_r3_delta(delta, tau)-delta*tau*phi_r3_deltatau(delta, tau))**2)/(2*delta*phi_r3_delta(delta, tau)+delta**2*phi_r3_deltadelta(delta, tau))))*R

def cv_r3(rho,T):
    """
    Specific isochoric heat capacity for region 3 IAPWS-IF97
    Chapter 7, Table 31
    
    Parameters
    ----------
    rho:            double              density (kg/m^3).
    T:              double              temperature (K).
    ijn_table_r3    np array            I,J,n coefficients of table 30
    n_r3_1          double              n(i=1) coefficient of table 30
    
    Returns
    -------
    double  specific isochoric heat capacity (J/kg/K)

    """
    delta = rho/rho_star_r3
    tau = T_star_r3/T
    
    return (-tau**2*phi_r3_tautau(delta, tau))*R

def w_r3(rho,T):
    """
    Speed of sound for region 3 IAPWS-IF97
    Chapter 7, Table 31
    
    Parameters
    ----------
    rho:            double              density (kg/m^3).
    T:              double              temperature (K).
    ijn_table_r3    np array            I,J,n coefficients of table 30
    n_r3_1          double              n(i=1) coefficient of table 30
    
    Returns
    -------
    double  Speed of sound (m/s)

    """
    delta = rho/rho_star_r3
    tau = T_star_r3/T
    
    return ((2*delta*phi_r3_delta(delta, tau)+(delta**2)*phi_r3_deltadelta(delta, tau)-((delta*phi_r3_delta(delta, tau)-delta*tau*phi_r3_deltatau(delta, tau))**2/((tau**2)*phi_r3_tautau(delta, tau))))*R*T)**.5

def rho_r3(p,T, rho_init = 0.0):
    """
    Density for region 3 IAPWS-IF97
    Itteration...
    
    Parameters
    ----------
    p:            double              pressure (Pa).
    T:            double              temperature (K).
    p:            double              initial guess density (kg/m^3) (optional).
    
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
        # set a range of +-2.5% of initial guess
        rho_0 = .975*rho_init
        rho_1 = 1.025*rho_init
        
    p_rho = lambda rho: p_r3(rho, T)
    rho = aux.bisection(p_rho, rho_0, rho_1, p, global_property.itt_tol, global_property.itt_max)
   
    return rho

"""
    Region 4
    Saturation pressure
    
    IAPWS-IF97, Chapter 8
    
    Validity range:
    273.15 K <= T >= 647.096 K . 
    
"""
def p_sat_t_r4(t):
    """
    Function calculating region 4 or IAPWS-IF97, Saturation pressure
    Chapter 8
    
    Eq. 30
    
    Validity range:
    Triple point to critical point
    273.15 K <= T >= 647.096 K . 
    
    Parameters
    ----------
    t_s:    double    temperature (K).

    Returns
    -------
    p_s :   double    saturation pressure (Pa) 

    """
    
    # Reducing quantaties
    p_star = 1 * 10**6 #1 MPa
    t_star = 1 # K
    
    #Table-34
    #Numerical values of the coefficients of the dimensionless saturation equations
    n_array = np.array([0.11670521452767*10**4,     #0
                        -0.72421316703206*10**6,    #1
                        -0.17073846940092*10**2,    #2
                        0.12020824702470*10**5,     #3
                        -0.32325550322333*10**7,    #4
                        0.14915108613530*10**2,     #5
                        -0.48232657361591*10**4,    #6
                        0.40511340542057*10**6,     #7
                        -0.23855557567849,          #8
                        0.65017534844798*10**3])    #9
    
    # Eq. 29b
    theta = t/t_star + n_array[8]/(t/t_star-n_array[9])
    
    # Coeficients for eq. 30
    A = theta**2+n_array[0]*theta+n_array[1]
    B = n_array[2]*theta**2+n_array[3]*theta+n_array[4]
    C = n_array[5]*theta**2+n_array[6]*theta+n_array[7]
    
    # Eq. 30
    p_s =  ((2*C/(-B+(B**2-4*A*C)**.5))**4)*p_star
    
    return p_s

def t_sat_p_r4(p):
    """
    Function calculating region 4 or IAPWS-IF97, Saturation temperature (backward equation) 
    Chapter 8
    
    eq. 31
    
    Validity range:
    from 611 Pa to critical pressure
    611.213 Pa <= p <= 22.064 MPa
    
    Parameters
    ----------
    p_s:    double    pressure (Pa).

    Returns
    -------
    t_s :   double    saturation temperature (K) 

    """
    # Reducing quantaties
    p_star = 1 * 10**6 #1 MPa
    t_star = 1 # K
    
    #Table-34
    #Numerical values of the coefficients of the dimensionless saturation equations
    n_array = np.array([0.11670521452767*10**4,     #0
                        -0.72421316703206*10**6,    #1
                        -0.17073846940092*10**2,    #2
                        0.12020824702470*10**5,     #3
                        -0.32325550322333*10**7,    #4
                        0.14915108613530*10**2,     #5
                        -0.48232657361591*10**4,    #6
                        0.40511340542057*10**6,     #7
                        -0.23855557567849,          #8
                        0.65017534844798*10**3])    #9
    
    # Eq. 29a
    beta = (p/p_star)**(1/4)
    
    # Coeficients for eq. 31
    E = beta**2 + n_array[2]*beta+n_array[5]
    F = n_array[0]*beta**2 + n_array[3]*beta + n_array[6]
    G = n_array[1]*beta**2 + n_array[4]*beta + n_array[7]
    
    D =  2*G/(-F-(F**2-4*E*G)**.5)
    
    # Eq. 31
    t_s = (n_array[9]+D-((n_array[9]+D)**2-4*(n_array[8]+n_array[9]*D))**.5)/2*t_star
    
    return t_s


"""
    Region 5
    
    IAPWS-IF97, Chapter 9
    
    Validity range:
    0 MPa < p <= 50 MPa
    1073.15 K <= T >= 2273.15 K.
    
    
"""

#Reducing quantaties
p_star_r5 = 1*10**6          #MPa
T_star_r5 = 1000             #K

# Gibbs free energy, ideal part region 5
# Table 37, Numerical values of the coefficients and exponents
#                   1, 2 , 3 , 4 , 5 , 6 , 7, 8, 9
J_ideal_r5 = np.array([0, 1, -3, -2, -1, 2])
n_ideal_r5 = np.array([-0.13179983674201*10**2,
                     0.68540841634434*10,
                    -0.24805148933466*10**-1,
                     0.36901534980333,
                    -0.31161318213925*10,
                    -0.32961626538917])
jn_array_ideal_r5 = np.array([[0, 1, -3, -2, -1, 2],
                              [-0.13179983674201*10**2,
                     0.68540841634434*10,
                    -0.24805148933466*10**-1,
                     0.36901534980333,
                    -0.31161318213925*10,
                    -0.32961626538917]])

def gamma_ideal_r5(pi,tau,jn = jn_array_ideal_r5):
    """
    Gibs free energy eq. 33, g(p,T), ideal part
    Region 5 of IAPWS-IF97
    Chapter 9, Table 40
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    jn_array:        np-array of int     coefficient J, n of table 37

    Returns
    -------
    gamma :   double    g(p,T), Gibbs free energy
    """
    
    gamma_ideal = np.log(pi)
    for j,n in zip(jn[0,:], jn[1,:]):
        gamma_ideal = gamma_ideal + n*tau**j
    
    return gamma_ideal

def gamma_ideal_r5_pi(pi,tau,jn = jn_array_ideal_r5):
    """
    Gibs free energy eq. 33, dg(p,T)/dpi, ideal part
    Region 5 of IAPWS-IF97
    Chapter 9, Table 40
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    jn_array:        np-array of int     coefficient J, n of table 37

    Returns
    -------
    gamma :   double    dg(p,T)/dpi, Gibbs free energy
    """
    
    return 1/pi
        
def gamma_ideal_r5_pipi(pi,tau,jn = jn_array_ideal_r5):
    """
    Gibs free energy eq. 33, d^2g(p,T)/dpi^2, ideal part
    Region 5 of IAPWS-IF97
    Chapter 9, Table 40
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    jn_array:        np-array of int     coefficient J, n of table 37

    Returns
    -------
    gamma :   double    d^2g(p,T)/dpi^2, Gibbs free energy
    """
    
    return -1/(pi**2)

def gamma_ideal_r5_tau(pi,tau,jn = jn_array_ideal_r5):
    """
    Gibs free energy eq. 33, dg(p,T)/dtau, ideal part
    Region 5 of IAPWS-IF97
    Chapter 9, Table 40
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    jn_array:        np-array of int     coefficient J, n of table 37

    Returns
    -------
    gamma :   double    dg(p,T)/dtau, Gibbs free energy
    """
   
    gamma_ideal = 0
    for j,n in zip(jn[0,:], jn[1,:]):
        gamma_ideal = gamma_ideal + n*j*tau**(j-1)
    
    return gamma_ideal

def gamma_ideal_r5_tautau(pi,tau,jn = jn_array_ideal_r5):
    """
    Gibs free energy eq. 33, d^2g(p,T)/dtau^2, ideal part
    Region 5 of IAPWS-IF97
    Chapter 9, Table 40
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    jn_array:        np-array of int     coefficient J, n of table 37

    Returns
    -------
    gamma :   double    d^2g(p,T)/dtau^2, Gibbs free energy
    """
   
    gamma_ideal = 0
    for j,n in zip(jn[0,:], jn[1,:]):
        gamma_ideal = gamma_ideal + n*j*(j-1)*tau**(j-2)
    
    return gamma_ideal

def gamma_ideal_r5_pitau(pi,tau,jn = jn_array_ideal_r5):
    """
    Gibs free energy eq. 33, d^2g(p,T)/dpi/dtau, ideal part
    Region 5 of IAPWS-IF97
    Chapter 9, Table 40
    
    Parameters
    ----------
    pi:     double              reduced pressure (-)
    tau:    double              reduced temperature (-)
    jn_array:        np-array of int     coefficient J, n of table 37

    Returns
    -------
    gamma :   double    d^2g(p,T)/dpi/dtau, Gibbs free energy
    """
    
    return 0

#Gibbs free energy, residual part region 5
# Table 38, Numerical values of the coefficients and exponents
#                         1, 2, 3, 4, 5, 6, 7, 8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43
ijn_array_residual_r5 = np.array([[1, 1, 1, 2, 2, 3],
                                  [1, 2, 3, 3, 9, 7],
                                  [ 0.15736404855259*10**-2,
                           0.90153761673944*10**-3,
                          -0.50270077677648*10**-2,
                           0.22440037409485*10**-5,
                          -0.41163275453471*10**-5,
                           0.37919454822955*10**-7]])

def gamma_residual_r5(pi, tau, ijn = ijn_array_residual_r5):
    """
    Gibs free energy eq. 34, g(p,T), residual part
    Region 5 of IAPWS-IF97
    Chapter 9, Table 41
    
    Parameters
    ----------
    pi:             double              reduced pressure (-)
    tau:            double              reduced temperature (-)
    ijn    np array            I,J,n coefficients of table 38

    Returns
    -------
    gamma :   double    g(p,T), Gibbs free energy
    """
    
    gamma_res = 0
    for i,j,n in zip(ijn[0,:], ijn[1,:], ijn[2,:]):
        gamma_res = gamma_res + n*pi**i*tau**j
    
    return gamma_res

def gamma_residual_r5_pi(pi, tau, ijn = ijn_array_residual_r5):
    """
    Gibs free energy eq. 34, dg(p,T)/dpi, residual part
    Region 5 of IAPWS-IF97
    Chapter 9, Table 41
    
    Parameters
    ----------
    pi:             double              reduced pressure (-)
    tau:            double              reduced temperature (-)
    ijn    np array            I,J,n coefficients of table 38

    Returns
    -------
    gamma :   double    dg(p,T)/dpi, Gibbs free energy
    """
    
    gamma_res = 0
    for i,j,n in zip(ijn[0,:], ijn[1,:], ijn[2,:]):
        gamma_res = gamma_res + n*i*pi**(i-1)*tau**j
    
    return gamma_res
   
def gamma_residual_r5_pipi(pi, tau, ijn = ijn_array_residual_r5):
    """
    Gibs free energy eq. 34, d^2g(p,T)/dpi^2, residual part
    Region 5 of IAPWS-IF97
    Chapter 9, Table 41
    
    Parameters
    ----------
    pi:             double              reduced pressure (-)
    tau:            double              reduced temperature (-)
    ijn    np array            I,J,n coefficients of table 38

    Returns
    -------
    gamma :   double    d^2g(p,T)/dpi^2, Gibbs free energy
    """
    
    gamma_res = 0
    for i,j,n in zip(ijn[0,:], ijn[1,:], ijn[2,:]):
        gamma_res = gamma_res + n*i*(i-1)*pi**(i-2)*tau**j
    
    return gamma_res

def gamma_residual_r5_tau(pi, tau, ijn = ijn_array_residual_r5):
    """
    Gibs free energy eq. 34, dg(p,T)/dtau, residual part
    Region 5 of IAPWS-IF97
    Chapter 9, Table 41
    
    Parameters
    ----------
    pi:             double              reduced pressure (-)
    tau:            double              reduced temperature (-)
    ijn    np array            I,J,n coefficients of table 38

    Returns
    -------
    gamma :   double    dg(p,T)/dtau, Gibbs free energy
    """
    
    gamma_res = 0
    for i,j,n in zip(ijn[0,:], ijn[1,:], ijn[2,:]):
        gamma_res = gamma_res + n*pi**i*j*tau**(j-1)
    
    return gamma_res
   
def gamma_residual_r5_tautau(pi, tau, ijn = ijn_array_residual_r5):
    """
    Gibs free energy eq. 34, d^2g(p,T)/dtau^2, residual part
    Region 5 of IAPWS-IF97
    Chapter 9, Table 41
    
    Parameters
    ----------
    pi:             double              reduced pressure (-)
    tau:            double              reduced temperature (-)
    ijn    np array            I,J,n coefficients of table 38

    Returns
    -------
    gamma :   double    d^2g(p,T)/dtau^2, Gibbs free energy
    """
        
    gamma_res = 0
    for i,j,n in zip(ijn[0,:], ijn[1,:], ijn[2,:]):
        gamma_res = gamma_res + n*pi**i*j*(j-1)*tau**(j-2)
    
    return gamma_res
   
def gamma_residual_r5_pitau(pi, tau, ijn = ijn_array_residual_r5):
    """
    Gibs free energy eq. 34, d^2g(p,T)/dpi/dtau, residual part
    Region 5 of IAPWS-IF97
    Chapter 9, Table 41
    
    Parameters
    ----------
    pi:             double              reduced pressure (-)
    tau:            double              reduced temperature (-)
    ijn             np array            I,J,n coefficients of table 38

    Returns
    -------
    gamma :   double    d^2g(p,T)/dpi/dtau, Gibbs free energy
    """
        
    gamma_res = 0
    for i,j,n in zip(ijn[0,:], ijn[1,:], ijn[2,:]):
        gamma_res = gamma_res + n*i*pi**(i-1)*j*tau**(j-1)
    
    return gamma_res 

# Relations of thermodynamic properties
# Table 39 Chapter 9
   
def v_r5(p,T):
    """
    Specific volume for region 5 IAPWS-IF97
    Chapter 9, Table 39
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).

    Returns
    -------
    double  specific volume (m^3/kg)

    """
    pi = p/p_star_r5
    tau = T_star_r5/T
   
    return pi*(gamma_ideal_r5_pi(pi, tau) + gamma_residual_r5_pi(pi, tau))*R*T/p

def u_r5(p,T):
    """
    Specific internal energy for region 5 IAPWS-IF97
    Chapter 9, Table 39
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific internal energy (J/kg)

    """
    pi = p/p_star_r5
    tau = T_star_r5/T
    
    return (tau*(gamma_ideal_r5_tau(pi, tau)+gamma_residual_r5_tau(pi, tau))-pi*(gamma_ideal_r5_pi(pi, tau)+gamma_residual_r5_pi(pi, tau)))*R*T

def s_r5(p,T):
    """
    Specific entropy for region 5 IAPWS-IF97
    Chapter 9, Table 39
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific entropy (J/kg/K)

    """
    pi = p/p_star_r5
    tau = T_star_r5/T

    return (tau*(gamma_ideal_r5_tau(pi, tau)+gamma_residual_r5_tau(pi, tau))-(gamma_ideal_r5(pi, tau)+gamma_residual_r5(pi, tau)))*R   

def h_r5(p,T):
    """
    Specific enthalpy for region 5 IAPWS-IF97
    Chapter 9, Table 39
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific enthalpy (J/kg)

    """
    pi = p/p_star_r5
    tau = T_star_r5/T
    
    return tau*(gamma_ideal_r5_tau(pi, tau)+gamma_residual_r5_tau(pi, tau))*R*T

def cp_r5(p,T):
    """
    Specific isobaric heat capacity for region 5 IAPWS-IF97
    Chapter 9, Table 39
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific isobaric heat capacity (J/kg/K)

    """
    pi = p/p_star_r5
    tau = T_star_r5/T
    
    return -1*(tau**2)*(gamma_ideal_r5_tautau(pi, tau)+gamma_residual_r5_tautau(pi, tau))*R

def cv_r5(p,T):
    """
    Specific isochoric heat capacity for region 5 IAPWS-IF97
    Chapter 9, Table 39
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  specific isochoric heat capacity (J/kg/K)

    """
    pi = p/p_star_r5
    tau = T_star_r5/T
    
    return (-1*(tau**2)*(gamma_ideal_r5_tautau(pi, tau)+gamma_residual_r5_tautau(pi, tau))-((1+pi*gamma_residual_r5_pi(pi, tau)-tau*pi*gamma_residual_r5_pitau(pi, tau))**2/(1-pi**2*gamma_residual_r5_pipi(pi, tau))))*R

def w_r5(p,T):
    """
    Speed of sound for region 5 IAPWS-IF97
    Chapter 9, Table 39
    
    Parameters
    ----------
    p:              double              pressure (Pa).
    T:              double              temperature (K).
    
    Returns
    -------
    double  Speed of sound (m/s)

    """
    pi = p/p_star_r5
    tau = T_star_r5/T
    
    return (((1+2*pi*gamma_residual_r5_pi(pi, tau)+pi**2*gamma_residual_r5_pi(pi, tau)**2)/((1-pi**2*gamma_residual_r5_pipi(pi, tau))+((1+pi*gamma_residual_r5_pi(pi, tau)-tau*pi*gamma_residual_r5_pitau(pi, tau))**2/(tau**2*(gamma_ideal_r5_tautau(pi, tau)+gamma_residual_r5_tautau(pi, tau))))))*R*T)**.5

def t_ph_r5(p,h):
    """
    Backwards function for 
    Temperature from p,h region 5 IAPWS-IF97
    
    No function presented in iapwsif97, use itteration instead
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    h:      double  enthalpy (J/kg).

    Returns
    -------
    double  Temperature (K)

    """
    T = 0.0
    
    # Range for itterations, must be within region 5...
    t0 = global_property.T_boundary_25
    t1 = global_property.T_boundary_5
    
    # Itterate
    h_t = lambda t: h_r5(p,t)
    T = aux.bisection(h_t, t0, t1, h, global_property.itt_tol, global_property.itt_max)
    
    return T

def t_ps_r5(p,s):
    """
    Backwards function for 
    Temperature from p,s region 5 IAPWS-IF97
    
    No function presented in iapwsif97, use itteration instead
    
    Parameters
    ----------
    p:      double  pressure (Pa).
    s:      double  entropy (J/kg/K).

    Returns
    -------
    double  Temperature (K)

    """
    T = 0.0
    
    # Range for itterations, must be within region 5...
    t0 = global_property.T_boundary_25
    t1 = global_property.T_boundary_5
    
    # Itterate
    s_t = lambda t: s_r5(p,t)
    T = aux.bisection(s_t, t0, t1, s, global_property.itt_tol, global_property.itt_max)
    
    return T

if __name__ == "__main__":
    main()