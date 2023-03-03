# -*- coding: utf-8 -*-
"""
Validation of iapwsif97_class by comparison to NIST database.
https://webbook.nist.gov/cgi/fluid.cgi?ID=C7732185&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm&Type=IsoTherm&RefState=DEF&Action=Page
    
Fluid states from 0.1 to 100 MPa and 300 to 2000 K has been copied from the NIST database into a csv-file
which is then compared to calculated states using the iapwsif97 class.
Saturation states for liquid and vapor 0.1 to 22 Mpa and 373 to 640 K is also provided and compared in a seperate csv-file

The aim is to compare all functions in the main class. 
All functions except: vx_ph(p, h), vx_ps(p, s), vx_hs(h, s) since these are not part of the NIST data.

@author: Christoffer Rappmann, christoffer.rappman@gmail.com
"""

from pyfluidproperties.iapwsif97_class import iapwsif97 as w_prop
import numpy as np
import time

def main():
    print('Validation of IAPWSIF97-class, comparing results to NIST(IAPWSIF95)\n\n')
    
    # Initiate object
    w = w_prop(10*10**6, 100+273.15)
    
    # iapwsif97 results matrix
    
    # Properties
    # 0  1  2    3  4  5  6  7   8   9  10  11  12  13
    # T, p, rho, v, u, h, s, Cv, Cp, w, JT, my, tc, phase
    iapwsif_results_pt = np.zeros((0,14))
    iapwsif_results_pt_static = np.zeros((0,14))
    iapwsif_results_ph = np.zeros((0,14))
    iapwsif_results_ph_static = np.zeros((0,14))
    iapwsif_results_ps = np.zeros((0,14))
    iapwsif_results_ps_static = np.zeros((0,14))
    iapwsif_results_hs_static = np.zeros((0,14))
    
    iapwsif_results_h_prho = np.zeros((0,14))
    
    # Saturation properties
    # 0  1  2    3  4  5  6  7   8   9  10  11  12  13
    # T, p, rho, v, u, h, s, Cv, Cp, w, JT, my, tc, phase
    iapwsif_results_sat_vap_p = np.zeros((0,14))
    iapwsif_results_sat_vap_p_static = np.zeros((0,14))
    iapwsif_results_sat_vap_T = np.zeros((0,14))
    iapwsif_results_sat_vap_T_static = np.zeros((0,14))
    
    iapwsif_results_sat_vap_px = np.zeros((0,14))
    iapwsif_results_sat_vap_tx = np.zeros((0,14))
    # 0  1  2    3  4  5  6  7   8   9  10  11  12  13  14
    # T, p, rho, v, u, h, s, Cv, Cp, w, JT, my, tc, st, phase
    iapwsif_results_sat_liq_p = np.zeros((0,15))
    iapwsif_results_sat_liq_p_static = np.zeros((0,15))
    iapwsif_results_sat_liq_T = np.zeros((0,15))
    iapwsif_results_sat_liq_T_static = np.zeros((0,15))

    iapwsif_results_sat_liq_px = np.zeros((0,15))
    iapwsif_results_sat_liq_tx = np.zeros((0,15))
    
    iapwsif_results_twophase_p = np.zeros((0,14))
    iapwsif_results_twophase_px = np.zeros((0,14))
    iapwsif_results_twophase_tx = np.zeros((0,14))
    iapwsif_results_twophase_sat_hs = np.zeros((0,14))
    
    iapwsif_results_sat_p = np.zeros((0,15))
    
    iapwsif_results_twophase_x = np.zeros((0,3))
    
    ##############
    ### States ###
    ##############
    """
    Fluid states from NIST, 
    https://webbook.nist.gov/cgi/fluid.cgi?ID=C7732185&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm&Type=IsoTherm&RefState=DEF&Action=Page
    
    Range is specified on nist website and the corresponding table is copied to a csv-file with the following columns
    Temperature (K),Pressure (MPa),Density (kg/m3),Volume (m3/kg),Internal Energy (kJ/kg),Enthalpy (kJ/kg),Entropy (J/g*K),Cv (J/g*K),Cp (J/g*K),Sound Spd. (m/s),Joule-Thomson (K/MPa),Viscosity (Pa*s),Therm. Cond. (W/m*K),Phase
    """
    
    print('Fluid states from NIST,\nhttps://webbook.nist.gov/cgi/fluid.cgi?ID=C7732185&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm&Type=IsoTherm&RefState=DEF&Action=Page\nRange is specified on nist website and the corresponding table is copied to a csv-file with the following columns\nTemperature (K),Pressure (MPa),Density (kg/m3),Volume (m3/kg),Internal Energy (kJ/kg),Enthalpy (kJ/kg),Entropy (J/g*K),Cv (J/g*K),Cp (J/g*K),Sound Spd. (m/s),Joule-Thomson (K/MPa),Viscosity (Pa*s),Therm. Cond. (W/m*K),Phase')
    nist_header = np.array(['Temperature (K)','Pressure (MPa)','Density (kg/m3)','Volume (m3/kg)','Internal Energy (kJ/kg)','Enthalpy (kJ/kg)','Entropy (kJ/g*K)','Cv (kJ/g*K)','Cp (kJ/g*K)','Sound Spd. (m/s)','Joule-Thomson (K/MPa)','Viscosity (Pa*s)','Therm. Cond. (W/m*K)','Phase'])
    # Ordinary data points
    data_file = 'fluid_states_IAPWSIF95_nist_0to2000K_01to100MPa.csv'
    nist_data = np.genfromtxt(data_file, delimiter=',',skip_header=True)
    # Saturation data pressure
    data_file_sat_liq = 'fluid_states_IAPWSIF95_nist_saturation_liquid_p.csv'
    data_file_sat_vap = 'fluid_states_IAPWSIF95_nist_saturation_vapor_p.csv'
    nist_data_sat_liq = np.genfromtxt(data_file_sat_liq, delimiter=',',skip_header=True)
    nist_data_sat_vap = np.genfromtxt(data_file_sat_vap, delimiter=',',skip_header=True)
    
    print(f'\nCalculating fluid properties for all fluid states provided in:\n - {data_file}\n - {data_file_sat_liq}\n - {data_file_sat_vap}\n ')
    
    ##################################
    #### Verification calculations ###
    ##################################
    print('\nCalculating properties using _pt, _ph, _ps, _hs and _prho functions....')
    tic = time.perf_counter()
    for i in range(0, len(nist_data)):
        
        # update_pt & *_pt
        p = nist_data[i,1]*10**6 # Pa
        T = nist_data[i,0] # K
        
        w.update_pt(p, T)
                
        iapwsif_results_pt = np.append(iapwsif_results_pt, [[w.T, w.p*10**-6, 1/w.v, w.v, w.u*10**-3, w.h*10**-3, w.s*10**-3, w.cv*10**-3, w.cp*10**-3, w.w, np.nan, w.my, w.tc, np.nan]], axis = 0)
        iapwsif_results_pt_static = np.append(iapwsif_results_pt_static, [[T, p*10**-6, w_prop.rho_pt(p, T), w_prop.v_pt(p,T), w_prop.u_pt(p,T)*10**-3, w_prop.h_pt(p, T)*10**-3, w_prop.s_pt(p, T)*10**-3, w_prop.cv_pt(p, T)*10**-3, w_prop.cp_pt(p, T)*10**-3, w_prop.w_pt(p, T), np.nan, w_prop.my_pt(p, T), w_prop.tc_pt(p, T), np.nan]], axis = 0)
        
        # update_ph & *_ph
        h = nist_data[i,5]*10**3 # J/kg  
            
        w.update_ph(p, h)
        iapwsif_results_ph = np.append(iapwsif_results_ph, [[w.T, w.p*10**-6, 1/w.v, w.v, w.u*10**-3, w.h*10**-3, w.s*10**-3, w.cv*10**-3, w.cp*10**-3, w.w, np.nan, w.my, w.tc, np.nan]], axis = 0)
        iapwsif_results_ph_static = np.append(iapwsif_results_ph_static, [[w_prop.T_ph(p, h), p*10**-6, w_prop.rho_ph(p, h), w_prop.v_ph(p,h), w_prop.u_ph(p,h)*10**-3, h*10**-3, w_prop.s_ph(p, h)*10**-3, w_prop.cv_ph(p, h)*10**-3, w_prop.cp_ph(p, h)*10**-3, w_prop.w_ph(p, h), np.nan, w_prop.my_ph(p, h), w_prop.tc_ph(p, h), np.nan]], axis = 0)
         
        # update_ps  & *_ps
        s = nist_data[i,6]*10**3 # J/kg  
            
        w.update_ps(p, s)
        iapwsif_results_ps = np.append(iapwsif_results_ps, [[w.T, w.p*10**-6, 1/w.v, w.v, w.u*10**-3, w.h*10**-3, w.s*10**-3, w.cv*10**-3, w.cp*10**-3, w.w, np.nan, w.my, w.tc, np.nan]], axis = 0)
        iapwsif_results_ps_static = np.append(iapwsif_results_ps_static, [[w_prop.T_ps(p, s), p*10**-6, w_prop.rho_ps(p, s), w_prop.v_ps(p,s), w_prop.u_ps(p,s)*10**-3, w_prop.h_ps(p, s)*10**-3, s*10**-3, w_prop.cv_ps(p, s)*10**-3, w_prop.cp_ps(p, s)*10**-3, w_prop.w_ps(p, s), np.nan, w_prop.my_ps(p, s), w_prop.tc_ps(p, s), np.nan]], axis = 0)
        
        
        # *_hs
        iapwsif_results_hs_static = np.append(iapwsif_results_hs_static, [[w_prop.T_hs(h, s), w_prop.p_hs(h, s)*10**-6, w_prop.rho_hs(h, s), w_prop.v_hs(h,s), w_prop.u_hs(h,s)*10**-3, h*10**-3, s*10**-3, w_prop.cv_hs(h, s)*10**-3, w_prop.cp_hs(h, s)*10**-3, w_prop.w_hs(h, s), np.nan, w_prop.my_hs(h, s), w_prop.tc_hs(h, s), np.nan]], axis = 0)
        
        # h_phro
        rho = nist_data[i,2]
        iapwsif_results_h_prho = np.append(iapwsif_results_h_prho, [[np.nan, p*10**-6, rho, np.nan, np.nan, w_prop.h_prho(p, rho)*10**-3, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]], axis = 0)
        
    # Saturation properties
    print('Calculating saturation properties using _px, _Tx, L_p, V_p, L_T, V_T and sat_s functions....')
    for j in range(0, len(nist_data_sat_liq)):
        # update_px & *L_p
        p = nist_data_sat_liq[j,1]*10**6
                
        w.update_px(p, 0.0)
        iapwsif_results_sat_liq_p = np.append(iapwsif_results_sat_liq_p, [[w.T, w.p*10**-6, 1/w.v, w.v, w.u*10**-3, w.h*10**-3, w.s*10**-3, w.cv*10**-3, w.cp*10**-3, w.w, np.nan, w.my, w.tc, w.st, np.nan]], axis = 0)
        iapwsif_results_sat_liq_p_static = np.append(iapwsif_results_sat_liq_p_static, [[w_prop.Tsat_p(p), p*10**-6, w_prop.rhoL_p(p), w_prop.vL_p(p), w_prop.uL_p(p)*10**-3, w_prop.hL_p(p)*10**-3, w_prop.sL_p(p)*10**-3, w_prop.cvL_p(p)*10**-3, w_prop.cpL_p(p)*10**-3, w_prop.wL_p(p), np.nan, w_prop.myL_p(p), w_prop.tcL_p(p), w_prop.st_p(p), np.nan]], axis = 0)
        
        #*L_T
        T = nist_data_sat_liq[j,0]
                
        w.update_tx(T, 0.0)
        iapwsif_results_sat_liq_T = np.append(iapwsif_results_sat_liq_T, [[w.T, w.p*10**-6, 1/w.v, w.v, w.u*10**-3, w.h*10**-3, w.s*10**-3, w.cv*10**-3, w.cp*10**-3, w.w, np.nan, w.my, w.tc, w.st, np.nan]], axis = 0)
        iapwsif_results_sat_liq_T_static = np.append(iapwsif_results_sat_liq_T_static, [[T, w_prop.psat_t(T)*10**-6, w_prop.rhoL_T(T), w_prop.vL_T(T), w_prop.uL_T(T)*10**-3, w_prop.hL_T(T)*10**-3, w_prop.sL_T(T)*10**-3, w_prop.cvL_T(T)*10**-3, w_prop.cpL_T(T)*10**-3, w_prop.wL_T(T), np.nan, w_prop.myL_T(T), w_prop.tcL_T(T), w_prop.st_t(T), np.nan]], axis = 0)
        
        #h_tx(t, x), h_px(p, x), s_tx(t, x), s_px(p, x)      
        iapwsif_results_sat_liq_px = np.append(iapwsif_results_sat_liq_px, [[np.nan, np.nan, np.nan, np.nan, np.nan, w_prop.h_px(p, 0.0)*10**-3, w_prop.s_px(p, 0.0)*10**-3, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]], axis = 0)
        iapwsif_results_sat_liq_tx = np.append(iapwsif_results_sat_liq_tx, [[np.nan, np.nan, np.nan, np.nan, np.nan, w_prop.h_tx(T, 0.0)*10**-3, w_prop.s_tx(T, 0.0)*10**-3, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]], axis = 0)
        
        #psat_s(s)
        s = nist_data_sat_liq[j,6]*10**3
        iapwsif_results_sat_p = np.append(iapwsif_results_sat_p, [[np.nan, w_prop.psat_s(s)*10**-6, np.nan, np.nan, np.nan, np.nan, s*10**-3, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]], axis = 0)
        
    for j in range(0, len(nist_data_sat_vap)):    
        # update_px & *V_p
        p = nist_data_sat_vap[j,1]*10**6
        
        w.update_px(p, 1.0)
        iapwsif_results_sat_vap_p = np.append(iapwsif_results_sat_vap_p, [[w.T, w.p*10**-6, 1/w.v, w.v, w.u*10**-3, w.h*10**-3, w.s*10**-3, w.cv*10**-3, w.cp*10**-3, w.w, np.nan, w.my, w.tc, np.nan]], axis = 0)
        iapwsif_results_sat_vap_p_static = np.append(iapwsif_results_sat_vap_p_static, [[w_prop.Tsat_p(p), p*10**-6, w_prop.rhoV_p(p), w_prop.vV_p(p), w_prop.uV_p(p)*10**-3, w_prop.hV_p(p)*10**-3, w_prop.sV_p(p)*10**-3, w_prop.cvV_p(p)*10**-3, w_prop.cpV_p(p)*10**-3, w_prop.wV_p(p), np.nan, w_prop.myV_p(p), w_prop.tcV_p(p), np.nan]], axis = 0)
        
        # update_Tx & *V_T
        T = nist_data_sat_vap[j,0]
        w.update_tx(T, 1.0)
        iapwsif_results_sat_vap_T = np.append(iapwsif_results_sat_vap_T, [[w.T, w.p*10**-6, 1/w.v, w.v, w.u*10**-3, w.h*10**-3, w.s*10**-3, w.cv*10**-3, w.cp*10**-3, w.w, np.nan, w.my, w.tc, np.nan]], axis = 0)
        iapwsif_results_sat_vap_T_static = np.append(iapwsif_results_sat_vap_T_static, [[T, w_prop.psat_t(T)*10**-6, w_prop.rhoV_T(T), w_prop.vV_T(T), w_prop.uV_T(T)*10**-3, w_prop.hV_T(T)*10**-3, w_prop.sV_T(T)*10**-3, w_prop.cvV_T(T)*10**-3, w_prop.cpV_T(T)*10**-3, w_prop.wV_T(T), np.nan, w_prop.myV_T(T), w_prop.tcV_T(T), np.nan]], axis = 0)
      
        #h_tx(t, x), h_px(p, x), s_tx(t, x), s_px(p, x)      
        iapwsif_results_sat_vap_px = np.append(iapwsif_results_sat_vap_px, [[np.nan, np.nan, np.nan, np.nan, np.nan, w_prop.h_px(p, 1.0)*10**-3, w_prop.s_px(p, 1.0)*10**-3, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]], axis = 0)
        iapwsif_results_sat_vap_tx = np.append(iapwsif_results_sat_vap_tx, [[np.nan, np.nan, np.nan, np.nan, np.nan, w_prop.h_tx(T, 1.0)*10**-3, w_prop.s_tx(T, 1.0)*10**-3, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]], axis = 0)
        
    # Two-phase properties
    
    print('Calculating two-phase properties using _px and _Tx functions....')
    x = .5 # massfraction to validate
    
    # Check that liq and sat data contains equal amount of data-points
    if len(nist_data_sat_liq) == len(nist_data_sat_vap):
        for k in range(0, len(nist_data_sat_vap)):
            p = nist_data_sat_liq[k,1]*10**6
            
            w.update_px(p,x)
            iapwsif_results_twophase_p = np.append(iapwsif_results_twophase_p, [[w.T, w.p*10**-6, w.rho, w.v, w.u*10**-3, w.h*10**-3, w.s*10**-3, w.cv*10**-3, w.cp*10**-3, w.w, np.nan, w.my, w.tc, np.nan]], axis = 0)
        
            #h_px(p, x), s_px(p, x)      
            iapwsif_results_twophase_px = np.append(iapwsif_results_twophase_px, [[np.nan, p*10**-6, np.nan, np.nan, np.nan, w_prop.h_px(p, x)*10**-3, w_prop.s_px(p, x)*10**-3, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]], axis = 0)
        
            T = nist_data_sat_liq[k,0]
            #h_tx(t, x), s_tx(t, x)
            iapwsif_results_twophase_tx = np.append(iapwsif_results_twophase_tx, [[T, np.nan, np.nan, np.nan, np.nan, w_prop.h_tx(T, x)*10**-3, w_prop.s_tx(T, x)*10**-3, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]], axis = 0)
            
            h = ((1-x)*nist_data_sat_liq[k,5] + (x)*nist_data_sat_vap[k,5])*10**3
            s = ((1-x)*nist_data_sat_liq[k,6] + (x)*nist_data_sat_vap[k,6])*10**3

            # Tsat_hs(h, s), psat_hs(h, s)
            iapwsif_results_twophase_sat_hs = np.append(iapwsif_results_twophase_sat_hs, [[w_prop.Tsat_hs(h, s), w_prop.psat_hs(h, s)*10**-6, np.nan, np.nan, np.nan, h*10**-3, s*10**-3, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]], axis = 0)
            
            # x_ph, x_ps, x_hs
            iapwsif_results_twophase_x = np.append(iapwsif_results_twophase_x, [[w_prop.x_ph(p, h), w_prop.x_ps(p, s), w_prop.x_hs(h, s)]], axis = 0)
    
    toc = time.perf_counter()
    
    ###########################
    ### Verification checks ###
    ###########################
    
    print('\n\t##################################################################' +
          '\n\t###            Verification of property functions              ###' +
          '\n\t###        - Objects: update_pt, update_ph, update_ps          ###' + 
          '\n\t###        - Static functions: *_pt, *_ph, *_ps, *_hs          ###' +
          '\n\t##################################################################\n')
    # update_pt & *_pt
    results_diff_pt = diff(nist_data,iapwsif_results_pt)
    print_max_diff('update_pt', results_diff_pt, nist_data, nist_header)
    
    results_diff_pt_static = diff(nist_data,iapwsif_results_pt_static)
    print_max_diff('*_pt', results_diff_pt_static, nist_data, nist_header)
    
    # update_ph  & *_ph
    results_diff_ph = diff(nist_data,iapwsif_results_ph)
    print_max_diff('update_ph', results_diff_ph, nist_data, nist_header)
    
    results_diff_ph_static = diff(nist_data,iapwsif_results_ph_static)
    print_max_diff('*_ph', results_diff_ph_static, nist_data, nist_header)
    
    # update_ps  & *_ps
    results_diff_ps = diff(nist_data,iapwsif_results_ps)
    print_max_diff('update_ps', results_diff_ps, nist_data, nist_header)
    
    results_diff_ps_static = diff(nist_data,iapwsif_results_ps_static)
    print_max_diff('*_ps', results_diff_ps_static, nist_data, nist_header)
    
    # *_hs
    results_diff_hs_static = diff(nist_data,iapwsif_results_hs_static)
    print_max_diff('*_hs', results_diff_hs_static, nist_data, nist_header)
    
    # h_prho
    results_diff_h_prho = diff(nist_data, iapwsif_results_h_prho)
    print_max_diff('h_prho', results_diff_h_prho, nist_data, nist_header)
    
    print('\n\t###################################################################' +
          '\n\t###        Verification of Saturation property functions        ###' +
          '\n\t###                       objects             static func.      ###'
          '\n\t###   - Saturated Liquid: update_px(p,x=0.0), *L_p,*L_T,*_*x    ###' + 
          '\n\t###   - Saturated Vapor:  update_px(p,x=1.0), *V_p,*V_T *_*x    ###' +
          '\n\t###   - Saturation pressure:                  psat_s            ###'
          '\n\t###################################################################\n')
    # update_px(x=0)
    results_diff_sat_liq_p = diff(nist_data_sat_liq, iapwsif_results_sat_liq_p)
    print_max_diff('Saturated Liquid, update_px(x=0.0)', results_diff_sat_liq_p, nist_data_sat_liq, nist_header)
    
    # *L_p
    results_diff_sat_liq_p_stat = diff(nist_data_sat_liq,iapwsif_results_sat_liq_p_static)
    print_max_diff('Saturated Liquid, *L_p - functions', results_diff_sat_liq_p_stat, nist_data_sat_liq, nist_header)
    
    # update_px(x=1.0)
    results_diff_sat_vap_p = diff(nist_data_sat_vap,iapwsif_results_sat_vap_p)
    print_max_diff('Saturated Vapor, update_px(x=1.0)', results_diff_sat_vap_p, nist_data_sat_vap, nist_header)
    
    # *V_p
    results_diff_sat_vap_p_stat = diff(nist_data_sat_vap,iapwsif_results_sat_vap_p_static)
    print_max_diff('Saturated Vapor, *V_p - functions', results_diff_sat_vap_p_stat, nist_data_sat_vap, nist_header)
    
    # update_Tx(x=0)
    results_diff_sat_liq_T = diff(nist_data_sat_liq, iapwsif_results_sat_liq_T)
    print_max_diff('Saturated Liquid, update_Tx(x=0.0)', results_diff_sat_liq_T, nist_data_sat_liq, nist_header)
    
    # *L_T
    results_diff_sat_liq_T = diff(nist_data_sat_liq,iapwsif_results_sat_liq_T_static)
    print_max_diff('Saturated Liquid, *L_T - functions', results_diff_sat_liq_T, nist_data_sat_liq, nist_header)
    
    # update_Tx(x=1.0)
    results_diff_sat_vap_T = diff(nist_data_sat_vap,iapwsif_results_sat_vap_T)
    print_max_diff('Saturated Vapor, update_Tx(x=1.0)', results_diff_sat_vap_T, nist_data_sat_vap, nist_header)
    
    # *V_T
    results_diff_sat_vap_T = diff(nist_data_sat_vap,iapwsif_results_sat_vap_T_static)
    print_max_diff('Saturated Vapor, *V_T - functions', results_diff_sat_vap_T, nist_data_sat_vap, nist_header)
    
    # psat_s
    results_diff_sat_s = diff(nist_data_sat_liq, iapwsif_results_sat_p)
    print_max_diff('Saturation pressure, psat_s', results_diff_sat_s, nist_data_sat_liq, nist_header)
    
    #h_px(p, x), s_px(p, x)
    results_diff_sat_px = diff(nist_data_sat_liq,iapwsif_results_sat_liq_px)
    print_max_diff('Saturated Liquid h_px(p, x=0.0), s_px(p, x=0.0)', results_diff_sat_px, nist_data_sat_liq, nist_header)
    
    #h_tx(p, x), s_tx(p, x)
    results_diff_sat_tx = diff(nist_data_sat_liq,iapwsif_results_sat_liq_tx)
    print_max_diff('Saturated Liquid h_tx(p, x=0.0), s_tx(p, x=0.0)', results_diff_sat_tx, nist_data_sat_liq, nist_header)
    
    #h_px(p, x = 1.0), s_px(p, x = 1.0)
    results_diff_sat_vap_px = diff(nist_data_sat_vap,iapwsif_results_sat_vap_px)
    print_max_diff('Saturated Vapor h_px(p, x=1.0), s_px(p, x=1.0)', results_diff_sat_vap_px, nist_data_sat_vap, nist_header)
    
    #h_tx(p, x = 1.0), s_tx(p, x = 1.0)
    results_diff_sat_vap_tx = diff(nist_data_sat_vap,iapwsif_results_sat_vap_tx)
    print_max_diff('Saturated Vapor h_tx(p, x=1.0), s_tx(p, x=1.0)', results_diff_sat_vap_tx, nist_data_sat_vap, nist_header)
    
    
    print('\n\t###################################################################' +
          '\n\t###         Verification of Two phase proerty functions         ###' +
          '\n\t###         - Objects: update_px(p,x=0.5)                       ###' + 
          '\n\t###         - Static functions: *_*x(p,x=0.5), *sat_hs          ###' +
          '\n\t###                             x_ph, x_ps, x_hs                ###'
          '\n\t###################################################################\n')
    
    #two-phase
    #update_px(p, x = 0.5)
    results_diff_twophase_p = diff((nist_data_sat_liq[:,[0,1,2,3,4,5,6,7,8,9,10,11,12,13]]*(1-x)+nist_data_sat_vap*(x)),iapwsif_results_twophase_p)
    print_max_diff('Two phase update_px(p, x = 0.5)', results_diff_twophase_p, nist_data_sat_vap, nist_header)
    
    #h_px(p, x = 0.5), s_px(p, x = 0.5)
    results_diff_twophase_px = diff((nist_data_sat_liq[:,[0,1,2,3,4,5,6,7,8,9,10,11,12,13]]*(1-x)+nist_data_sat_vap*(x)),iapwsif_results_twophase_px)
    print_max_diff('Two phase h_px(p, x=0.5), s_px(p, x=0.5)', results_diff_twophase_px, nist_data_sat_vap, nist_header)
    
    #h_tx(p, x = 0.5), s_tx(p, x = 0.5)
    results_diff_twophase_tx = diff((nist_data_sat_liq[:,[0,1,2,3,4,5,6,7,8,9,10,11,12,13]]*(1-x)+nist_data_sat_vap*(x)),iapwsif_results_twophase_tx)
    print_max_diff('Two phase h_tx(p, x=0.5), s_tx(p, x=0.5)', results_diff_twophase_tx, nist_data_sat_vap, nist_header)
    
    #Tsat_hs(h, s), psat_hs(h, s)
    results_diff_twophase_sat_hs = diff((nist_data_sat_liq[:,[0,1,2,3,4,5,6,7,8,9,10,11,12,13]]*(1-x)+nist_data_sat_vap*(x)),iapwsif_results_twophase_sat_hs)
    print_max_diff('Two phase Tsat_hs, psat_hs', results_diff_twophase_sat_hs, nist_data_sat_vap, nist_header)
    
    # x_ph, x_ps, x_hs
    results_diff_twophase_x = diff(np.ones((len(iapwsif_results_twophase_x),3))*x,iapwsif_results_twophase_x)
    max_diff_x = np.nanmax(results_diff_twophase_x, axis = 0)
    print(f'Tested functions: x_ph, x_ps, x_hs with x = {x:.2f}\n\nx_ph(p, h)\tx_ps(p, s)\tx_hs(h, s)\n{max_diff_x[0]:.3f} %\t\t{max_diff_x[1]:.3f} %\t\t{max_diff_x[2]:.3f} %' + 
          '\n__________________________________________________________________________')
            
    print(f'\n\nValidation run completed!\n{len(nist_data)+len(nist_data_sat_vap)+len(nist_data_sat_liq)} states calculated in {(toc-tic):.1f} s\n\n')
    
def diff(a,b):
    """
    Compares np-arrays and calculates relative difference in each position
    
    Parameters
    ----------
    a : np-array of double,     array a
    b : np-array of double,     array b
    
    Returns
    -------
    np-array of double,     relative difference
    """
    
    return np.divide(np.absolute(np.subtract(a,b)),a)

def print_max_diff(desc,a,data,header):
    """
    Determine max-value in array, print and return value and position
    
    Parameters
    ----------
    desc : string,  description
    a : np-array of double,     array a
    
    Returns
    -------
    double, [int, int],     relative difference, position in array
    """
    
    max_diff = np.nanmax(a)
    max_diff_pos = np.where(a == max_diff)
    max_diff_prop = np.nanmax(a, axis = 0)
    
    print(f'Tested functions: {desc}\n' + 
          f'Maximum error: {max_diff*100:.2f} % for "{header[max_diff_pos[1]][0]}" @ p = {data[max_diff_pos[0],1][0]} MPa, T = {data[max_diff_pos[0],0][0]} K')
    # T, p, rho, v, u, h, s, Cv, Cp, w, JT, my, tc, phase
    prop_table = ['T', 'p', 'rho', 'v', 'u', 'h', 's', 'Cv', 'Cp', 'w', 'JT', 'my', 'tc', 'st']
    
    for i in range(0,14):
        if not np.isnan(max_diff_prop[i]):
            print(f'{prop_table[i]}\t\t',end = '')
    print('\n', end = '')
    for j in range(0,14):
        if not np.isnan(max_diff_prop[j]):
            print(f'{max_diff_prop[j]*100:.2f} %\t',end = '')
    print('\n__________________________________________________________________________')
    return max_diff,max_diff_pos

if __name__ == "__main__":
    main()