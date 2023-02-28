# -*- coding: utf-8 -*-
"""
Validation of iapwsif97_class according to validation tables in IAPWS-IF97 releases and suplementary releases.
Validation is done using the highest available function, for example if a class-functions exists this is the
to be used. If no class-function exist but there is a sub-function implemented this will be used directly.

In this way the class i validated (for used functions) as well as the sub-functions called directly or the
ones used by the class.

@author: Christoffer Rappmann, christoffer.rappman@gmail.com
"""

from iapwsif97_class import iapwsif97 as w_prop
from iapwsif97_class import if97, if97_tvphps, if97_phs, f08 

import numpy as np
from tabulate import tabulate

def main():
    
    print('Validation of iapwsif97_class according to validation tables in IAPWS-IF97 releases')
    
    # Settings
    owerwrite_verification = False
    show_table_console = False
    write_table_html = True
    
    # Initiate object
    w = w_prop() 
    
    
    # Verification matrix
    # 0    1    2  3  4  5  6  7   8  9  10    11  12  13  14
    # rel, tbl, p, t, v, h, u, s, cp, w, rho, ,my, th, st, reg
    verification_tbls = np.zeros((0,15))
    
    ##############
    ### Inputs ###
    ##############
    
    # input data matrix
    #    0           1     2        3            4       5               6           8
    #    rel         tbl,  type,    p,           t,      h,              s,          rho
    verification_input_data = np.array([
        ['R7-97',    5,    'pt',    3*10**6,     300,    0,              0,          0],
        ['R7-97',    5,    'pt',    80*10**6,    300,    0,              0,          0],
        ['R7-97',    5,    'pt',    3*10**6,     500,    0,              0,          0],
        ['R7-97',    7,    'ph',    3*10**6,     0,      500*10**3,      0,          0],
        ['R7-97',    7,    'ph',    80*10**6,    0,      500*10**3,      0,          0],
        ['R7-97',    7,    'ph',    80*10**6,    0,      1500*10**3,     0,          0],
        ['R7-97',    9,    'ps',    3*10**6,     0,      0,              .5*10**3,   0],
        ['R7-97',    9,    'ps',    80*10**6,    0,      0,              .5*10**3,   0],
        ['R7-97',    9,    'ps',    80*10**6,    0,      0,              3*10**3,    0],
        ['R7-97',    15,   'pt',    .0035*10**6, 300,    0,              0,          0],
        ['R7-97',    15,   'pt',    .0035*10**6, 700,    0,              0,          0],
        ['R7-97',    15,   'pt',    30*10**6,    700,    0,              0,          0],
        ['R7-97',    18,   'ptm',   1*10**6,     450,    0,              0,          0],
        ['R7-97',    18,   'ptm',   1*10**6,     440,    0,              0,          0],
        ['R7-97',    18,   'ptm',   1.5*10**6,   450,    0,              0,          0],
        ['R7-97',    24,   'ph',    .001*10**6,  0,      3000*10**3,     0,          0],
        ['R7-97',    24,   'ph',    3*10**6,     0,      3000*10**3,     0,          0],
        ['R7-97',    24,   'ph',    3*10**6,     0,      4000*10**3,     0,          0],
        ['R7-97',    24,   'ph',    5*10**6,     0,      3500*10**3,     0,          0],
        ['R7-97',    24,   'ph',    5*10**6,     0,      4000*10**3,     0,          0],
        ['R7-97',    24,   'ph',    25*10**6,    0,      3500*10**3,     0,          0],
        ['R7-97',    24,   'ph',    40*10**6,    0,      2700*10**3,     0,          0],
        ['R7-97',    24,   'ph',    60*10**6,    0,      2700*10**3,     0,          0],
        ['R7-97',    24,   'ph',    60*10**6,    0,      3200*10**3,     0,          0],
        ['R7-97',    29,   'ps',    .1*10**6,    0,      0,              7.5*10**3,  0],
        ['R7-97',    29,   'ps',    .1*10**6,    0,      0,              8*10**3,    0],
        ['R7-97',    29,   'ps',    2.5*10**6,   0,      0,              8*10**3,    0],
        ['R7-97',    29,   'ps',    8*10**6,     0,      0,              6*10**3,    0],
        ['R7-97',    29,   'ps',    8*10**6,     0,      0,              7.5*10**3,  0],
        ['R7-97',    29,   'ps',    90*10**6,    0,      0,              6*10**3,    0],
        ['R7-97',    29,   'ps',    20*10**6,    0,      0,              5.75*10**3, 0],
        ['R7-97',    29,   'ps',    80*10**6,    0,      0,              5.25*10**3, 0],
        ['R7-97',    29,   'ps',    80*10**6,    0,      0,              5.75*10**3, 0],
        ['R7-97',    33,   'prho',  0,           650,    0,              0,          500],
        ['R7-97',    33,   'prho',  0,           650,    0,              0,          200],
        ['R7-97',    33,   'prho',  0,           750,    0,              0,          500],
        ['R7-97',    35,   'psat',  0,           300,    0,              0,          0],
        ['R7-97',    35,   'psat',  0,           500,    0,              0,          0],
        ['R7-97',    35,   'psat',  0,           600,    0,              0,          0],
        ['R7-97',    36,   'tsat',  .1*10**6,    0,      0,              0,          0],
        ['R7-97',    36,   'tsat',  1*10**6,     0,      0,              0,          0],
        ['R7-97',    36,   'tsat',  10*10**6,    0,      0,              0,          0],
        ['R7-97',    42,   'pt',    .5*10**6,   1500,    0,              0,          0],
        ['R7-97',    42,   'pt',    30*10**6,   1500,    0,              0,          0],
        ['R7-97',    42,   'pt',    30*10**6,   2000,    0,              0,          0],
        ['SR2-01',   3,    'hs',    0,          0,       .001*10**3,     0,          0],
        ['SR2-01',   3,    'hs',    0,          0,       90*10**3,       0,          0],
        ['SR2-01',   3,    'hs',    0,          0,       1500*10**3,     3.4*10**3,  0],
        ['SR2-01',   9,    'hs',    0,          0,       2800*10**3,     6.5*10**3,  0],
        ['SR2-01',   9,    'hs',    0,          0,       2800*10**3,     9.5*10**3,  0],
        ['SR2-01',   9,    'hs',    0,          0,       4100*10**3,     9.5*10**3,  0],
        ['SR2-01',   9,    'hs',    0,          0,       2800*10**3,     6*10**3,    0],
        ['SR2-01',   9,    'hs',    0,          0,       3600*10**3,     6*10**3,    0],
        ['SR2-01',   9,    'hs',    0,          0,       3600*10**3,     7*10**3,    0],
        ['SR2-01',   9,    'hs',    0,          0,       2800*10**3,     5.1*10**3,  0],
        ['SR2-01',   9,    'hs',    0,          0,       2800*10**3,     5.8*10**3,  0],
        ['SR2-01',   9,    'hs',    0,          0,       3400*10**3,     5.8*10**3,  0],
        ['SR3-03',   5,    'ph',    20*10**6,   0,       1700*10**3,     0,          0],
        ['SR3-03',   5,    'ph',    50*10**6,   0,       2000*10**3,     0,          0],
        ['SR3-03',   5,    'ph',    100*10**6,  0,       2100*10**3,     0,          0],
        ['SR3-03',   5,    'ph',    20*10**6,   0,       2500*10**3,     0,          0],
        ['SR3-03',   5,    'ph',    50*10**6,   0,       2400*10**3,     0,          0],
        ['SR3-03',   5,    'ph',    100*10**6,  0,       2700*10**3,     0,          0],
        ['SR3-03',   8,    'phv',   20*10**6,   0,       1700*10**3,     0,          0],
        ['SR3-03',   8,    'phv',   50*10**6,   0,       2000*10**3,     0,          0],
        ['SR3-03',   8,    'phv',   100*10**6,  0,       2100*10**3,     0,          0],
        ['SR3-03',   8,    'phv',   20*10**6,   0,       2500*10**3,     0,          0],
        ['SR3-03',   8,    'phv',   50*10**6,   0,       2400*10**3,     0,          0],
        ['SR3-03',   8,    'phv',   100*10**6,  0,       2700*10**3,     0,          0],
        ['SR3-03',   12,   'ps',    20*10**6,   0,       0,              3.8*10**3,  0],
        ['SR3-03',   12,   'ps',    50*10**6,   0,       0,              3.6*10**3,  0],
        ['SR3-03',   12,   'ps',    100*10**6,  0,       0,              4.0*10**3,  0],
        ['SR3-03',   12,   'ps',    20*10**6,   0,       0,              5.0*10**3,  0],
        ['SR3-03',   12,   'ps',    50*10**6,   0,       0,              4.5*10**3,  0],
        ['SR3-03',   12,   'ps',    100*10**6,  0,       0,              5.0*10**3,  0],
        ['SR3-03',   15,   'psv',    20*10**6,  0,       0,              3.8*10**3,  0],
        ['SR3-03',   15,   'psv',    50*10**6,  0,       0,              3.6*10**3,  0],
        ['SR3-03',   15,   'psv',    100*10**6, 0,       0,              4.0*10**3,  0],
        ['SR3-03',   15,   'psv',    20*10**6,  0,       0,              5.0*10**3,  0],
        ['SR3-03',   15,   'psv',    50*10**6,  0,       0,              4.5*10**3,  0],
        ['SR3-03',   15,   'psv',    100*10**6, 0,       0,              5.0*10**3,  0],
        ['SR3-03',   18,   'psath3', 0,         0,       1700*10**3,     0,          0],
        ['SR3-03',   18,   'psath3', 0,         0,       2000*10**3,     0,          0],
        ['SR3-03',   18,   'psath3', 0,         0,       2400*10**3,     0,          0],
        ['SR3-03',   20,   'psats3', 0,         0,       0,              3.8*10**3,  0],
        ['SR3-03',   20,   'psats3', 0,         0,       0,              4.2*10**3,  0],
        ['SR3-03',   20,   'psats3', 0,         0,       0,              5.2*10**3,  0],
        ['SR4-04',   5,    'hs',     0,         0,       1700*10**3,     3.8*10**3,  0],
        ['SR4-04',   5,    'hs',     0,         0,       2000*10**3,     4.2*10**3,  0],
        ['SR4-04',   5,    'hs',     0,         0,       2100*10**3,     4.3*10**3,  0],
        ['SR4-04',   5,    'hs',     0,         0,       2600*10**3,     5.1*10**3,  0],
        ['SR4-04',   5,    'hs',     0,         0,       2400*10**3,     4.7*10**3,  0],
        ['SR4-04',   5,    'hs',     0,         0,       2700*10**3,     5.0*10**3,  0],
        ['SR4-04',   11,   'hprim1', 0,         0,       0,              1.0*10**3,  0],
        ['SR4-04',   11,   'hprim1', 0,         0,       0,              2.0*10**3,  0],
        ['SR4-04',   11,   'hprim1', 0,         0,       0,              3.0*10**3,  0],
        ['SR4-04',   11,   'hprim3', 0,         0,       0,              3.8*10**3,  0],
        ['SR4-04',   11,   'hprim3', 0,         0,       0,              4.0*10**3,  0],
        ['SR4-04',   11,   'hprim3', 0,         0,       0,              4.2*10**3,  0],
        ['SR4-04',   18,   'hbis2',  0,         0,       0,              7.0*10**3,  0],
        ['SR4-04',   18,   'hbis2',  0,         0,       0,              8.0*10**3,  0],
        ['SR4-04',   18,   'hbis2',  0,         0,       0,              9.0*10**3,  0],
        ['SR4-04',   18,   'hbis23', 0,         0,       0,              5.5*10**3,  0],
        ['SR4-04',   18,   'hbis23', 0,         0,       0,              5.0*10**3,  0],
        ['SR4-04',   18,   'hbis23', 0,         0,       0,              4.5*10**3,  0],
        ['SR4-04',   24,   'hb13',   0,         0,       0,              3.7*10**3,  0],
        ['SR4-04',   24,   'hb13',   0,         0,       0,              3.6*10**3,  0],
        ['SR4-04',   24,   'hb13',   0,         0,       0,              3.5*10**3,  0],
        ['SR4-04',   26,   'Tb23',   0,         0,       2600*10**3,     5.1*10**3,  0],
        ['SR4-04',   26,   'Tb23',   0,         0,       2700*10**3,     5.15*10**3,  0],
        ['SR4-04',   26,   'Tb23',   0,         0,       2800*10**3,     5.2*10**3,  0],
        ['SR4-04',   29,   'tsaths', 0,         0,       1800*10**3,     5.3*10**3,  0],
        ['SR4-04',   29,   'tsaths', 0,         0,       2400*10**3,     6.0*10**3,  0],
        ['SR4-04',   29,   'tsaths', 0,         0,       2500*10**3,     5.5*10**3,  0],
        ['SR5-05',    5,   'pt',     50*10**6,  630,     0,              0,          0],
        ['SR5-05',    5,   'pt',     80*10**6,  670,     0,              0,          0],
        ['SR5-05',    5,   'pt',     50*10**6,  710,     0,              0,          0],
        ['SR5-05',    5,   'pt',     80*10**6,  750,     0,              0,          0],
        ['SR5-05',    5,   'pt',     20*10**6,  630,     0,              0,          0],
        ['SR5-05',    5,   'pt',     30*10**6,  650,     0,              0,          0],
        ['SR5-05',    5,   'pt',     26*10**6,  656,     0,              0,          0],
        ['SR5-05',    5,   'pt',     30*10**6,  670,     0,              0,          0],
        ['SR5-05',    5,   'pt',     26*10**6,  661,     0,              0,          0],
        ['SR5-05',    5,   'pt',     30*10**6,  675,     0,              0,          0],
        ['SR5-05',    5,   'pt',     26*10**6,  671,     0,              0,          0],
        ['SR5-05',    5,   'pt',     30*10**6,  690,     0,              0,          0],
        ['SR5-05',    5,   'pt',     23.6*10**6,649,     0,              0,          0],
        ['SR5-05',    5,   'pt',     24*10**6,  650,     0,              0,          0],
        ['SR5-05',    5,   'pt',     23.6*10**6,652,     0,              0,          0],
        ['SR5-05',    5,   'pt',     24*10**6,  654,     0,              0,          0],
        ['SR5-05',    5,   'pt',     23.6*10**6,653,     0,              0,          0],
        ['SR5-05',    5,   'pt',     24*10**6,  655,     0,              0,          0],
        ['SR5-05',    5,   'pt',     23.5*10**6,655,     0,              0,          0],
        ['SR5-05',    5,   'pt',     24*10**6,  660,     0,              0,          0],
        ['SR5-05',    5,   'pt',     23*10**6,  660,     0,              0,          0],
        ['SR5-05',    5,   'pt',     24*10**6,  670,     0,              0,          0],
        ['SR5-05',    5,   'pt',     22.6*10**6,646,     0,              0,          0],
        ['SR5-05',    5,   'pt',     23*10**6,  646,     0,              0,          0],
        ['SR5-05',    5,   'pt',     22.6*10**6,648.6,   0,              0,          0],
        ['SR5-05',    5,   'pt',     22.8*10**6,649.3,   0,              0,          0],
        ['SR5-05',    5,   'pt',     22.6*10**6,649,     0,              0,          0],
        ['SR5-05',    5,   'pt',     22.8*10**6,649.7,   0,              0,          0],
        ['SR5-05',    5,   'pt',     22.6*10**6,649.1,   0,              0,          0],
        ['SR5-05',    5,   'pt',     22.8*10**6,649.9,   0,              0,          0],
        ['SR5-05',    5,   'pt',     22.6*10**6,649.4,   0,              0,          0],
        ['SR5-05',    5,   'pt',     22.8*10**6,650.2,   0,              0,          0],
        ['SR5-05',    5,   'pt',     21.1*10**6,640,     0,              0,          0],
        ['SR5-05',    5,   'pt',     21.8*10**6,643,     0,              0,          0],
        ['SR5-05',    5,   'pt',     21.1*10**6,644,     0,              0,          0],
        ['SR5-05',    5,   'pt',     21.8*10**6,648,     0,              0,          0],
        ['SR5-05',    5,   'pt',     19.1*10**6,635,     0,              0,          0],
        ['SR5-05',    5,   'pt',     20*10**6,  638,     0,              0,          0],
        ['SR5-05',    5,   'pt',     17*10**6,  626,     0,              0,          0],
        ['SR5-05',    5,   'pt',     20*10**6,  640,     0,              0,          0],
        ['SR5-05',    13,  'pt',     21.5*10**6,644.6,   0,              0,          0],
        ['SR5-05',    13,  'pt',     22*10**6,  646.1,   0,              0,          0],
        ['SR5-05',    13,  'pt',     22.5*10**6,648.6,   0,              0,          0],
        ['SR5-05',    13,  'pt',     22.3*10**6,647.9,   0,              0,          0],
        ['SR5-05',    13,  'pt',     22.15*10**6, 647.5, 0,              0,          0],
        ['SR5-05',    13,  'pt',     22.3*10**6,648.1,   0,              0,          0],
        ['SR5-05',    13,  'pt',     22.11*10**6,648,    0,              0,          0],
        ['SR5-05',    13,  'pt',     22.3*10**6,649,     0,              0,          0],
        ['SR5-05',    13,  'pt',     22*10**6,  646.84,  0,              0,          0],
        ['SR5-05',    13,  'pt',     22.064*10**6, 647.05, 0,            0,          0],
        ['SR5-05',    13,  'pt',     22*10**6,  646.89,  0,              0,          0],
        ['SR5-05',    13,  'pt',     22.064*10**6, 647.15, 0,            0,          0],
        ['R12-08',    4,   'visc',   0.0,       298.15,  0,              0,          998],
        ['R12-08',    4,   'visc',   0.0,       298.15,  0,              0,          1200],
        ['R12-08',    4,   'visc',   0.0,       373.15,  0,              0,          1000],
        ['R12-08',    4,   'visc',   0.0,       433.15,  0,              0,          1],
        ['R12-08',    4,   'visc',   0.0,       433.15,  0,              0,          1000],
        ['R12-08',    4,   'visc',   0.0,       873.15,  0,              0,          1],
        ['R12-08',    4,   'visc',   0.0,       873.15,  0,              0,          100],
        ['R12-08',    4,   'visc',   0.0,       873.15,  0,              0,          600],
        ['R12-08',    4,   'visc',   0.0,       1173.15, 0,              0,          1],
        ['R12-08',    4,   'visc',   0.0,       1173.15, 0,              0,          100],
        ['R12-08',    4,   'visc',   0.0,       1173.15, 0,              0,          400],
        ['R15-11',    7,   'pt',     20*10**6,  620,     0,              0,          0],
        ['R15-11',    7,   'pt',     50*10**6,  620,     0,              0,          0],
        ['R15-11',    8,   'pt',     .3*10**6,  650,     0,              0,          0],
        ['R15-11',    8,   'pt',     50*10**6,  800,     0,              0,          0],
        ['R1-76',     1,   'pt',     1.0*10**5, 0.01001+273.15,  0,         0,          0.0],
        ['R1-76',     1,   'pt',     1.0*10**5, 50+273.15,  0,         0,          0.0],
        ['R1-76',     1,   'pt',     1.0*10**5, 100+273.15,  0,         0,          0.0],
        ['R1-76',     1,   'pt',     1.0*10**5, 150+273.15,  0,         0,          0.0],
        ['R1-76',     1,   'pt',     1.0*10**5, 200+273.15,  0,         0,          0.0],
        ['R1-76',     1,   'pt',     1.0*10**5, 250+273.15,  0,         0,          0.0],
        ['R1-76',     1,   'pt',     1.0*10**5, 300+273.15,  0,         0,          0.0],
        ['R1-76',     1,   'pt',     1.0*10**5, 350+273.15,  0,         0,          0.0],
        ['R1-76',     1,   'pt',     1.0*10**5, 370+273.15,  0,         0,          0.0],
        ['R12-08',    5,   'trho_r3',   0.0,    647.35,  0,              0,          222],
        ['R12-08',    5,   'trho_r3',   0.0,    647.35,  0,              0,          272],
        ['R12-08',    5,   'trho_r3',   0.0,    647.35,  0,              0,          322],
        ['R12-08',    5,   'trho_r3',   0.0,    647.35,  0,              0,          372],
        ['R12-08',    5,   'trho_r3',   0.0,    647.35,  0,              0,          422],
        ['R15-11',    9,   'trho_r3',   0.0,    647.35,  0,              0,          222],
        ['R15-11',    9,   'trho_r3',   0.0,    647.35,  0,              0,          322]
        ])

    ##################################
    #### Verification calculations ###
    ##################################
    
    for i in range(0, len(verification_input_data)):
        
        ### Class-functions
        
        # p-T
        if verification_input_data[i,2] == 'pt':
            p = float(verification_input_data[i,3])
            T = float(verification_input_data[i,4])
            
            w.update_pt(p, T)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], w.p*10**-6, w.T, w.v, w.h*10**-3, w.u*10**-3, w.s*10**-3, w.cp*10**-3, w.w, w.rho, w.my, w.tc, w.st, w.region]],axis = 0)
        
        # p-h
        elif verification_input_data[i,2] == 'ph':
            p = float(verification_input_data[i,3])
            h = float(verification_input_data[i,5])
            
            w.update_ph(p, h)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], w.p*10**-6, w.T, 0, h*10**-3, 0, 0, 0, 0, 0, 0,0,0,  w.region]],axis = 0)
        
        # p-s  
        elif verification_input_data[i,2] == 'ps':
            p = float(verification_input_data[i,3])
            s = float(verification_input_data[i,6])
            
            w.update_ps(p, s)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], w.p*10**-6, w.T, 0, 0, 0, s*10**-3, 0, 0, 0, 0,0,0,  w.region]],axis = 0)
        
        # h-s
        if verification_input_data[i,2] == 'hs':
            h = float(verification_input_data[i,5])
            s = float(verification_input_data[i,6])
            
            region = if97_phs.find_region_hs(h, s)
            p_hs = w.p_hs(h,s)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], p_hs*10**-6, 0, 0, h*10**-3, 0, s*10**-3, 0, 0, 0, 0,0,0,  region]],axis = 0)
        
        # Specific volume, p-h
        elif verification_input_data[i,2] == 'phv':
            p = float(verification_input_data[i,3])
            h = float(verification_input_data[i,5])
            
            region = if97.find_region_ph(p, h)
            v_ph = w.v_ph(p, h)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], p*10**-6, 0, v_ph, h*10**-3, 0, 0, 0, 0, 0, 0,0,0,  region]],axis = 0)
        
        # Specific volume, p-s
        elif verification_input_data[i,2] == 'psv':
            p = float(verification_input_data[i,3])
            s = float(verification_input_data[i,6])
            
            region = if97.find_region_ps(p, s)
            v_ps = w.v_ps(p, s)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], p*10**-6, 0, v_ps, 0, s*10**-3, 0, 0, 0, 0, 0,0,0,  region]],axis = 0)
        
        # saturation pressure, T
        elif verification_input_data[i,2] == 'psat':
            T = float(verification_input_data[i,4])
            
            region = 4
            p_sat = w.psat_t(T)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], p_sat*10**-6, T, 0, 0, 0, 0, 0, 0, 0, 0,0,0,  region]],axis = 0)

        # saturation temperature, p
        elif verification_input_data[i,2] == 'tsat':
            p = float(verification_input_data[i,3])
            
            region = 4
            T_sat = w.Tsat_p(p)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0], verification_input_data[i,1], p*10**-6, T_sat, 0, 0, 0, 0, 0, 0, 0, 0,0,0, region]],axis = 0)
        
        # saturation temperature, hs
        elif verification_input_data[i,2] == 'tsaths':
            h = float(verification_input_data[i,5])
            s = float(verification_input_data[i,6])
            
            region = if97_phs.find_region_hs(h, s)
            t_sat = w.Tsat_hs(h, s)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], 0, t_sat, h*10**-3, 0, s*10**-3, 0, 0, 0, 0, 0,0,0, region]],axis = 0)
        
        # trho_r3
        elif verification_input_data[i,2] == 'trho_r3':
            T = float(verification_input_data[i,4])
            rho = float(verification_input_data[i,7])
            
            p = if97.p_r3(rho,T)
            if if97.find_region(p, T) == 3:
                w_prop.set_rho_r3_itt(True) # itterate rho for region 3 to ger higher accuracy

                my = w.my_pt(p, T)
                tc = w.tc_pt(p, T)
                
                w_prop.set_rho_r3_itt(False)
            else:
                
                my = 0.0
                tc = 0.0
                print('Warning, trho_r3, not region 3...')
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], p*10**-6, T, 0, 0, 0, 0, 0, 0, rho,  my, tc,0,region]],axis = 0)
        
        ### property sub functions
        
        # p-T metastable (region 2)
        elif verification_input_data[i,2] == 'ptm':
            p = float(verification_input_data[i,3])
            T = float(verification_input_data[i,4])

            region = 2
             
            v_m = if97.v_r2_metastable(p,T)
            h_m = if97.h_r2_metastable(p,T)
            u_m = if97.u_r2_metastable(p,T)
            s_m = if97.s_r2_metastable(p,T)
            cp_m = if97.cp_r2_metastable(p,T)
            w_m = if97.w_r2_metastable(p,T)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], p*10**-6, T, v_m, h_m*10**-3, u_m*10**-3, s_m*10**-3, cp_m*10**-3, w_m, 0, 0,0,0,  region]],axis = 0)
        
        # p-rho           
        elif verification_input_data[i,2] == 'prho':
            T = float(verification_input_data[i,4])
            rho = float(verification_input_data[i,7])
            
            region = 3
            
            p_3 = if97.p_r3(rho, T)
            h_3 = if97.h_r3(rho, T)
            u_3 = if97.u_r3(rho, T)
            s_3 = if97.s_r3(rho, T)
            cp_3 = if97.cp_r3(rho, T)
            w_3 = if97.w_r3(rho, T)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], p_3*10**-6, T, 0, h_3*10**-3, u_3*10**-3, s_3*10**-3, cp_3*10**-3, w_3, rho, 0,0,0,  region]],axis = 0)
        
        # saturation pressure, h
        elif verification_input_data[i,2] == 'psath3':
            h = float(verification_input_data[i,5])
            
            region = 3
            p_sat = if97_tvphps.p_3sat_h(h)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], p_sat*10**-6, 0, 0, h*10**-3, 0, 0, 0, 0, 0, 0,0,0,  region]],axis = 0)
            
        # saturation pressure, s
        elif verification_input_data[i,2] == 'psats3':
            s = float(verification_input_data[i,6])
            
            region = 3
            p_sat = if97_tvphps.p_3sat_s(s)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], p_sat*10**-6, 0, 0, 0, s*10**-3, 0, 0, 0, 0, 0,0,0,  region]],axis = 0)
        
        ### Region boundary sub functions
        
        # h_prim_s_r1
        elif verification_input_data[i,2] == 'hprim1':
            s = float(verification_input_data[i,6])
            
            region = -1
            h = if97_phs.h_prim_s_r1(s)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], 0, 0, h*10**-3, 0, s*10**-3, 0, 0, 0, 0,  0,0,0, region]],axis = 0)

        # h_prim_s_r3
        elif verification_input_data[i,2] == 'hprim3':
            s = float(verification_input_data[i,6])
            
            region = -1
            h = if97_phs.h_prim_s_r3(s)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], 0, 0, h*10**-3, 0, s*10**-3, 0, 0, 0, 0,  0, 0,0,region]],axis = 0)

        # h_bis_s_r2
        elif verification_input_data[i,2] == 'hbis2':
            s = float(verification_input_data[i,6])
            
            region = -1
            h = if97_phs.h_bis_s_r2(s)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], 0, 0, h*10**-3, 0, s*10**-3, 0, 0, 0, 0,  0,0,0, region]],axis = 0)

        # h_bis_s_r23
        elif verification_input_data[i,2] == 'hbis23':
            s = float(verification_input_data[i,6])
            
            region = -1
            h = if97_phs.h_bis_s_r23(s)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], 0, 0, h*10**-3, 0, s*10**-3, 0, 0, 0, 0,  0,0,0, region]],axis = 0)

        # h_s_b13
        elif verification_input_data[i,2] == 'hb13':
            s = float(verification_input_data[i,6])
            
            region = -1
            h = if97_phs.h_s_b13(s)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], 0, 0, h*10**-3, 0, s*10**-3, 0, 0, 0, 0,  0,0,0, region]],axis = 0)

        # T_hs_b23
        elif verification_input_data[i,2] == 'Tb23':
            h = float(verification_input_data[i,5])
            s = float(verification_input_data[i,6])
            
            region = -1
            T = if97_phs.T_hs_b23(h,s)
            
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], 0, T, h*10**-3, 0, s*10**-3, 0, 0, 0, 0,  0,0,0, region]],axis = 0)
           
        # Viscosity
        elif verification_input_data[i,2] == 'visc':
            T = float(verification_input_data[i,4])
            rho = float(verification_input_data[i,7])
            
            """
            #Many of the validation data-points is outside iapwsif97 valid region
            #hence the need to access the function directly.
            
            rho_p = lambda pp: w.rho_pt(pp, T)
            p = aux.bisection(rho_p,10,99*10**6,rho,0.000001,1000)
            region = if97.find_region(p, T)
            
            my = w.my_pt(p, T)
            """
            my = f08.my_trho(T, rho)
            verification_tbls = np.append(verification_tbls,[[  verification_input_data[i,0],verification_input_data[i,1], p*10**-6, T, 0, 0, 0, 0, 0, 0, rho,  my, 0,0,region]],axis = 0)
        
    ###########################
    ### Verification checks ###
    ###########################
    
    # Save verification tables
    if owerwrite_verification:
        np.save('verification_tbl_data', verification_tbls)
        
    # Load verification tables and check if correct.
    verification_tbls_original = np.load('verification_tbl_data.npy')
    
    if np.array_equal(verification_tbls, verification_tbls_original):
        print('\nAll verification tables correct!')
    else:
        print('\nverification tables not correct!')
        
        for i in range(0,len(verification_tbls)):
            if (verification_tbls == verification_tbls_original)[i,:].all() != True:
                print(f'\n\nrad = {i}\n{(verification_tbls == verification_tbls_original)[i,:]}\n')
                print(verification_tbls[i,:])
                print(verification_tbls_original[i,:])

    #####################################
    ### Show results for manual check ###
    #####################################
    
    # Print results...
    header = ['Release', 'Table', 'p (MPa)', 'T (K)', 'v (m^3/kg)', 'h (kJ/kg)', 'u (kJ/kg)', 's (kJ/kg/K)', 'cp (kJ/kg/K)', 'w (m/s)', 'rho (kg/m^3)','viscosity (Pa*s)','thermal cond. (W/m)','surface tens. (N/m)','Region']
    fformat =['',        ''    ,'.13f'   ,'.9f'    ,'.9e'        ,'.9e'       ,'.9e'       ,'.9e'         ,'.9e'          ,'.9e'     ,'.9e'         ,'.9e', '.9e','.9e', '.9e']
    # Print table to console    
    if show_table_console:
        print(tabulate(verification_tbls, headers = header, tablefmt = 'plain', numalign='right',floatfmt=fformat))
    
    # Print table to html-file
    if write_table_html:
        file = 'IAPWSIF97_Validation.html'
        with open(file, 'w') as f:
            f.write(tabulate(verification_tbls, headers = header, tablefmt = 'html', numalign='right',floatfmt=fformat))
        print(f'\nVerification data printed to: {file}')

    print('\nThere may be warnings printed out, this need not be a fault. One case will be a warning for the viscosity, which uses another iapws realase with a different validity range. thus the tested point might be outside of this region. \n\n-Verifications from SR5-05 regardning the boundary equations are not included because these will indirectly be tested when the correct region is selected for other verification tables. \n- Verifications for the viscosity with critical enhancement does not match exatly due to the table presented is meant to be used in conjuction with IAPWS95, IAPWSIF97 has a lower accuracy\n- Verifications for the thermal conductivity near the critical point. table 9, does not match exactly. Probably due to the method used to calculate the discrete derivative is different.\n\nValidation should be done with the rho_r3_itt setting set to False (Default) everywhere except for trho_r3.')

if __name__ == "__main__":
    main()

