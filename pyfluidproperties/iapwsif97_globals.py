# -*- coding: utf-8 -*-
"""
Global constants/properties used within the different modules in the iapwsif-97 and iapwsf08 implementations.

@author: Christoffer Rappmann, christoffer.rappmann@gmail.com
"""
##########################
###  Global Settings   ###
##########################

# Tolerance and max itterations for where itteration is needed
itt_tol = 1e-6                                      # -
itt_max = 100

##########################
### Physical constants ###
##########################

### Specific gas constant ###
R =                     0.461526*10**3              # J/kg/K    Specific gas constant used in iapwsif97 (Eq. 1)
R_f11 =                 0.46151805*10**3            # J/kg/K    Specific gas constant used in iapwsf11 (eq-6)


### Critical properties ###
p_crit =                22.064 * 10**6              # Pa, Critical pressure
T_crit =                647.096                     # K          Eq. 2
h_crit =                2.087546845*10**6           # J/kg, critical enthalpy, from IAPWS SR4-04(2014)
s_crit =                4.41202148223476*10**3      # J/kg/K critical entropy, from IAPWS SR4-04(2014)
rho_crit =              322                         # kg/m^3     Eq. 4

### Triple point properties ###
pt =                    611.657                     # Pa, Pressure at triplepoint of water, https://en.wikipedia.org/wiki/Triple_point
Tt =                    273.16                      # K, Temperature at triplepoint of water, https://en.wikipedia.org/wiki/Triple_point

##################
### iapwsif-97 ###
##################

### Region Boundaries ###

# Pressures
p_min =                 0.000611*10**6              # Pa, lower pressure boundary (SR4-04 fig. 1)
p_r123_upper_boundary = 100*10**6                   # Pa
p_r5_upper_boundary =   50*10**6                    # Pa
p_r3_lower_boundary =   16.529164252604478*10**6    # Pa, Region 1/2-3 boundary, if97.p_sat_t_r4(T_boundary_13)

# Temperatures
T_boundary_1 =          273.15                      # K
T_boundary_13 =         623.15                      # K
T_boundary_25 =         1073.15                     # K
T_boundary_5 =          2273.15                     # K

# Enthalpies
h_r2_max =              4160.815448611362*10**3     # J/kg, Upper bound for region 2, h_pt(0.000611*10**6,1073.15)
h_r5_max =              7376.980141800623*10**3     # J/kg

h_r3_lower =            1670.8582182746243*10**3    # J/kg, lowest liquid saturation enthalyp at p_r3_lower, used to narrow itteration window, calculated from region boundaries
h_r3_upper =            2563.7465232886752*10**3    # J/kg, highest liquid saturation enthalyp at p_r3_lower, used to narrow itteration window, calculated from region boundaries

h_b23_min =             2.563592004*10**6           # J/kg, Min enthalpy range of validity for b23 equation (SR4-04 ch. 4.6), h''(623.15K)
h_b23_max =             2.812942061*10**6           # J/kg, Max enthalpy range of validity for b23 equation (SR4-04 ch. 4.6), h''(100MPa. 863.15K)

# Entropies
s_r1_min =              -1.545495919*10**-1         # J/kg/K, lower bound for region 1, s'(273.15)
s_r2_max  =             11.921657919039976*10**3    # J/kg/K, Upper bound for region 2, s_pt(0.000611*10**6,1073.15)
s_r5_min =              6.522642312284054*10**3     # J/kg/K, lower bound for region 5, s_pt(50*10**6,1073.15)
s_r5_max =              13.90511670393939*10**3    # J/kg/K, Upper bound for region 5, s_pt(0.000611*10**6, 2273.15)

s_r3_lower =            3.778281339544306*10**3     # J/kg/K, lowest liquid saturation entropy at p_r3_lower  
s_r3_upper =            5.21133005419787*10**3      # J/kg/K, highest vapor saturation entropy at p_r3_lower

s_h_b13_min =           3.397782955 * 10**3         # J/kg/K, s(100 MPa,623.15 K), min range of validity h_b13 eq.
s_h_b13_max =           s_r3_lower                  # J/kg/K, s'(623.15), max range of validity h_b13 eq. 

s_h_bis_585 =           5.85*10**3                  # J/kg/K, break-point between h_bis_s_r23 and h_bis_s_r2 in iapwsif97_phs (SR4-04)
s_h_bis_273 =           9.155759395*10**3           # J/kg/K, Saturated vapor entropy in triple point
s_h_bis_623 =           5.210887825*10**3           # J/kg/K, s''(623.15K) 
s_h_b13 =               3.397782955*10**3           # J/kg/K, Boundary point region 1-3 s(100 MPa,623.15 K)

s_b23_min =             5.048096828*10**3           # J/kg/K, Min entropy range of validity for b23 equation (SR4-04 ch. 4.6)
s_b23_max =             5.260578707*10**3           # J/kg/K, Max entropy range of validity for b23 equation (SR4-04 ch. 4.6)

##################
###  iapwsf08  ###
##################

#### Range of validity ###
# Values for pressure and temperatures from eq-9 in iapwsf08 (R12-08)
# Pressures
p_r2_f08_upper =        300*10**6                   # Pa

p_r3_f08_lower =        p_r2_f08_upper                  # Pa
p_r3_f08_upper =        350*10**6                   # Pa

p_r4_f08_lower =        p_r3_f08_upper                  # Pa
p_r4_f08_upper =        500*10**6                   # Pa

p_r5_f08_lower =        p_r4_f08_upper                  # Pa
p_r5_f08_upper =        1000*10**6                  # Pa

# Temperatures
T_r1_f08_lower =        Tt                          # K
T_r1_f08_upper =        1173.15                     # K
 
T_r2_f08_upper =        1173.15                     # K

T_r3_f08_upper =        873.15                      # K

T_r4_f08_upper =        433.15                      # K

T_r5_f08_upper =        373.15                      # K

Tm_f08 =                Tt                          # K, Pressure dependent temperature melting-point, 
                                                    # assumed to be the same as triple point since this value always will be higher. 
                                                    # IAPWS-if97 is only valid to 273.15 and there is therfore no need to go lower 
                                                    # in temperature
Tref_f08 =              1.5*647.096                 # K, Reference temperature in table 3 of iapwsf08

# Critical enhancement term significant region eq-13 in iapwsf08 (R12-08)
T_f08_lower =           645.91                      # K
T_f08_upper =           650.77                      # K

rho_f08_lower =         245.8                       # kg/m^3
rho_f08_upper =         405.3                       # kg/m^3

##################
###  iapwsf11  ###
##################

#### Range of validity ###
# Values for pressure and temperatures from eq-14 in iapwsf11 (R15-11)
# Pressures
p_r2_f11_upper =        100*10**6                  # Pa

p_r3_f11_lower =        p_r2_f11_upper             # Pa
p_r3_f11_upper =        250*10**6                  # Pa

p_r4_f11_lower =        p_r3_f11_upper             # Pa
p_r4_f11_upper =        687*10**6                  # Pa

p_r5_f11_lower =        p_r4_f11_upper             # Pa
p_r5_f11_upper =        785*10**6                  # Pa

p_r6_f11_lower =        p_r5_f11_upper             # Pa
p_r6_f11_upper =        1000*10**6                 # Pa
# Temperatures
T_r1_f11_lower =        Tt                         # K
T_r1_f11_upper =        1173.15                    # K
 
T_r2_f11_upper =        1173.15                    # K

T_r3_f11_upper =        874.0                      # K

T_r4_f11_upper =        573.0                      # K

T_r5_f11_upper =        403.0                      # K

T_r6_f11_upper =        348.0                      # K

Tm_f11 =                Tt                          # K, Pressure dependent temperature melting-point, 
                                                    # assumed to be the same as triple point since this value always will be higher. 
                                                    # IAPWS-if97 is only valid to 273.15 and there is therfore no need to go lower 
                                                    # in temperature
Tref_f11 =              Tref_f08                    # K, Reference temperature in table 3 of iapwsf08
