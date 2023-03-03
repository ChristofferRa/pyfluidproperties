# -*- coding: utf-8 -*-
"""
Module containing functions implemented from equations presentend in the SR3-03(2014)
suplementary release of IAPWSIF97.

List of functions and code structure:
    Region boundaries
        h_3ab(p)
    Sub regions, region 3
        t_3a_ph(p,h)
        t_3b_ph(p,h)
        
        v_3a_ph(p,h)
        v_3b_ph(p,h)
        
        t_3a_ps(p,s)
        t_3b_ps(p,s)
        
        v_3a_ps(p,s)
        v_3b_ps(p,s)
        
        p_3sat_h(h)
        p_3sat_s(s)
        
    Determine region
        find_region_ph(p,h)
        find_region_ps(p,s)
        
    Properties region 3
        t_3_ph(p,h)
        v_3_ph(p,h)
        t_3_ps(p,s)
        v_3_ps(p,s)
        
@author: Christoffer Rappmann, christoffer.rappmann@gmail.com  
"""
# Global Constants
from . import iapwsif97_globals as global_property

import numpy as np
            
def h_3ab(p):
    """
    Backwards equation, boundary calculation betwwen region 3a/3b
    Approximates critical isentropic line, s_c = 4.412 021 482 234 76 kJ/kg/K
    Region 3 IAPWS-IF97
    Chapter 3.2, eq. 1
    
    Valid from critical point (22.064 MPa) to 100 MPa
    
    Parameters
    ----------
    p:      double  Pressure (Pa).

    Returns
    -------
    double  Enthalpy (J/kg)

    """
    
    p_star = 1*10**6 # 1 MPa
    h_star = 1*10**3 # 1 kJ/kg 
    
    # Numerical coefficients from table 2
    n = np.array([0.201464004206875*10**4, 0.374696550136983*10**1, -0.219921901054187*10**-1, 0.875131686009950*10**-4,])
    
    pi = p/p_star
    eta = n[0] + n[1]*pi + n[2]*pi**2 + n[3]*pi**3
    
    return eta * h_star

def t_3a_ph(p,h):
    """
    Backwards equation T_3a(p,h), subregion 3a. eq. 2
    Gives temperature from p,h input.
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    h :     double  Enthalpy (J/kg)

    Returns
    -------
    double Temperature (K)

    """
    t_star = 760 #K
    p_star = 100 * 10**6 # Pa
    h_star = 2300 * 10**3 # J/kg
    
    pi = p/p_star
    eta = h/h_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [-12, -12, -12, -12, -12, -12, -12, -12, -10, -10, -10, -8, -8, -8, -8, -5, -3, -2, -2, -2, -1, -1, 0,  0,  1,  3,  3,  4,  4,  10, 12]
    J = [0,   1,   2,   6,   14,  16,  20,  22,  1,   5,   12,  0,  2,  4,  10,  2,  0,  1,  3,  4,  0,  2, 0,  1,  1,  0,  1,  0,  3,  4,  5]
    n_array = [ -0.133645667811215*10**-6, 
                 0.455912656802978*10**-5,
                -0.146294640700979*10**-4,
                 0.639341312970080*10**-2,
                 0.372783927268847*10**3,
                -0.718654377460447*10**4,
                 0.573494752103400*10**6,
                -0.267569329111439*10**7,
                -0.334066283302614*10**-4,
                -0.245479214069597*10**-1,
                 0.478087847764996*10**2,
                 0.764664131818904*10**-5,
                 0.128350627676972*10**-2,
                 0.171219081377331*10**-1,
                -0.851007304583213*10**1,
                -0.136513461629781*10**-1,
                -0.384460997596657*10**-5,
                 0.337423807911655*10**-2,
                -0.551624873066791,
                 0.729202277107470,
                -0.992522757376041*10**-2,
                -0.119308831407288,
                 0.793929190615421,
                 0.454270731799386,
                 0.209998591259910,
                -0.642109823904738*10**-2,
                -0.235155868604540*10**-1,
                 0.252233108341612*10**-2,
                -0.764885133368119*10**-2,
                 0.136176427574291*10**-1,
                -0.133027883575669*10**-1]
    
    teta = 0
    for i,j,n in zip(I, J, n_array):
        teta += n * (pi + 0.240)**i * (eta - 0.615)**j
        
    return teta*t_star

def t_3b_ph(p,h):
    """
    Backwards equation T_3b(p,h), subregion 3b. eq. 3
    Gives temperature from p,h input.
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    h :     double  Enthalpy (J/kg)

    Returns
    -------
    double Temperature (K)

    """
    t_star = 860 #K
    p_star = 100 * 10**6 # Pa
    h_star = 2800 * 10**3 # J/kg
    
    pi = p/p_star
    eta = h/h_star
    
    #    1,   2,    3,   4,   5,   6,   7,   8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,27,28,29,30,31,32,33
    I = [-12, -12, -10, -10, -10, -10, -10, -8, -8, -8, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2, -1, -1, -1, -1, -1, -1, 0, 0, 1, 3, 5, 6, 8]
    J = [ 0,   1,   0,   1,   5,   10,  12,  0,  1,  2,  4,  10, 0,  1,  2,  0,  1,  5,  0,  4,  2,  4,  6, 10, 14, 16, 0, 2, 1, 1, 1, 1, 1]
    n_array = [ 0.323254573644920*10**-4,
               -0.127575556587181*10**-3,
               -0.475851877356068*10**-3,
                0.156183014181602*10**-2,
                0.105724860113781,
               -0.858514221132534*10**2,
                0.724140095480911*10**3,
                0.296475810273257*10**-2,
               -0.592721983365988*10**-2,
               -0.126305422818666*10**-1,
               -0.115716196364853,
                0.849000969739595*10**2,
               -0.108602260086615*10**-1,
                0.154304475328851*10**-1,
                0.750455441524466*10**-1,
                0.252520973612982*10**-1,
               -0.602507901232996*10**-1,
               -0.307622221350501*10**1,
               -0.574011959864879*10**-1,
                0.503471360939849*10**1,
               -0.925081888584834,
                0.391733882917546*10**1,
               -0.773146007130190*10**2,
                0.949308762098587*10**4,
               -0.141043719679409*10**7,
                0.849166230819026*10**7,
                0.861095729446704,
                0.323346442811720,
                0.873281936020439,
               -0.436653048526683,
                0.286596714529479,
               -0.131778331276228,
                0.676682064330275*10**-2]
    
    teta = 0
    for i,j,n in zip(I, J, n_array):
        teta += n * (pi + 0.298)**i * (eta - 0.720)**j
        
    return teta*t_star


def v_3a_ph(p,h):
    """
    Backwards equation v_3a(p,h), subregion 3a. eq. 4
    Gives specific volume from p,h input.
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    h :     double  Enthalpy (J/kg)

    Returns
    -------
    double Specific volume (m^3/kg)

    """
    v_star = 0.0028 #m^3/kg
    p_star = 100 * 10**6 # Pa
    h_star = 2100 * 10**3 # J/kg
    
    pi = p/p_star
    eta = h/h_star
    
    #    1,   2,    3,   4,   5,   6,   7,   8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,27,28,29,30,31,32,33
    I = [-12, -12, -12, -12, -10, -10, -10, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2, -1, -1, -1, -1, 0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 8]
    J = [ 6,   8,   12,  18,  4,   7,   10,  5,  12, 3,  4,  22, 2,  3,  7,  3,  16, 0,  1,  2,  3, 0, 1, 0, 1, 2, 0, 2, 0, 2, 2, 2]
    n_array = [ 0.529944062966028*10**-2,
               -0.170099690234461,
                0.111323814312927*10**2,
               -0.217898123145125*10**4,
               -0.506061827980875*10**-3,
                0.556495239685324,
               -0.943672726094016*10**1,
               -0.297856807561527,
                0.939353943717186*10**2,
                0.192944939465981*10**-1,
                0.421740664704763,
               -0.368914126282330*10**7,
               -0.737566847600639*10**-2,
               -0.354753242424366,
               -0.199768169338727*10**1,
                0.115456297059049*10**1,
                0.568366875815960*10**4,
                0.808169540124668*10**-2,
                0.172416341519307,
                0.104270175292927*10**1,
               -0.297691372792847,
                0.560394465163593,
                0.275234661176914,
               -0.148347894866012,
               -0.651142513478515*10**-1,
               -0.292468715386302*10**1,
                0.664876096952665*10**-1,
                0.352335014263844*10**1,
               -0.146340792313332*10**-1,
               -0.224503486668184*10**1,
                0.110533464706142*10**1,
               -0.408757344495612*10**-1]
    
    teta = 0
    for i,j,n in zip(I, J, n_array):
        teta += n * (pi + 0.128)**i * (eta - 0.727)**j
        
    return teta*v_star

def v_3b_ph(p,h):
    """
    Backwards equation v_3b(p,h), subregion 3b. eq. 5
    Gives specific volume from p,h input.
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    h :     double  Enthalpy (J/kg)

    Returns
    -------
    double Specific volume (m^3/kg)

    """
    v_star = 0.0088 #m^3/kg
    p_star = 100 * 10**6 # Pa
    h_star = 2800 * 10**3 # J/kg
    
    pi = p/p_star
    eta = h/h_star
    
    #    1,   2,    3,   4,   5,   6,   7,   8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,27,28,29,30,31,32,33
    I = [-12, -12, -8, -8, -8, -8, -8, -8, -6, -6, -6, -6, -6, -6, -4, -4, -4, -3, -3, -2, -2, -1, -1, -1, -1, 0, 1, 1, 2, 2]
    J = [ 0,   1,   0,  1,  3,  6,  7,  8,  0,  1,  2,  5,  6,  10, 3,  6,  10, 0,  2,  1,  2,  0,  1,  4,  5, 0, 0, 1, 2, 6]
    n_array = [ -0.225196934336318*10**-8,
                0.140674363313486*10**-7,
                0.233784085280560*10**-5,
               -0.331833715229001*10**-4,
                0.107956778514318*10**-2,
               -0.271382067378863,
                0.107202262490333*10,
               -0.853821329075382,
               -0.215214194340526*10**-4,
                0.769656088222730*10**-3,
               -0.431136580433864*10**-2,
                0.453342167309331,
               -0.507749535873652,
               -0.100475154528389*10**3,
               -0.219201924648793,
               -0.321087965668917*10,
                0.607567815637771*10**3,
                0.557686450685932*10**-3,
                0.187499040029550,
                0.905368030448107*10**-2,
                0.285417173048685,
                0.329924030996098*10**-1,
                0.239897419685483,
                0.482754995951394*10,
               -0.118035753702231*10**2,
                0.169490044091791,
               -0.179967222507787*10**-1,
                0.371810116332674*10**-1,
               -0.536288335065096*10**-1,
                0.160697101092520*10]
    
    teta = 0
    for i,j,n in zip(I, J, n_array):
        teta += n * (pi + 0.0661)**i * (eta - 0.720)**j
        
    return teta*v_star

def t_3a_ps(p,s):
    """
    Backwards equation T_3a(p,s), subregion 3a. eq. 6
    Gives temperature from p,s input.
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    s :     double  Entropy (J/kg/K)

    Returns
    -------
    double Temperature (K)

    """
    t_star = 760 #K
    p_star = 100 * 10**6 # Pa
    s_star = 4.4 * 10**3 # J/kg
    
    pi = p/p_star
    sigma = s/s_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [-12, -12, -10, -10, -10, -10, -8, -8, -8, -8, -6, -6, -6, -5, -5, -5, -4, -4, -4, -2, -2, -1, -1, 0, 0, 0, 1, 2, 2, 3, 8, 8, 10]
    J = [ 28,  32,  4,   10,  12,  14,  5,  7,  8,  28, 2,  6,  32, 0,  14, 32, 6,  10, 36, 1,  4,  1,  6, 0, 1, 4, 0, 0, 3, 2, 0, 1, 2]
    n_array = [ 0.150042008263875*10**10,
               -0.159397258480424*10**12,
                0.502181140217975*10**-3,
               -0.672057767855466*10**2,
                0.145058545404456*10**4,
               -0.823889534888890*10**4,
               -0.154852214233853,
                0.112305046746695*10**2,
               -0.297000213482822*10**2,
                0.438565132635495*10**11,
                0.137837838635464*10**-2,
               -0.297478527157462*10**1,
                0.971777947349413*10**13,
               -0.571527767052398*10**-4,
                0.288307949778420*10**5,
               -0.744428289262703*10**14,
                0.128017324848921*10**2,
               -0.368275545889071*10**3,
                0.664768904779177*10**16,
                0.449359251958880*10**-1,
               -0.422897836099655*10**1,
               -0.240614376434179,
               -0.474341365254924*10**1,
                0.724093999126110,
                0.923874349695897,
                0.399043655281015*10**1,
                0.384066651868009*10**-1,
               -0.359344365571848*10**-2,
               -0.735196448821653,
                0.188367048396131,
                0.141064266818704*10**-3,
               -0.257418501496337*10**-2,
                0.123220024851555*10**-2]
        
    teta = 0
    for i,j,n in zip(I, J, n_array):
        teta += n * (pi + 0.240)**i * (sigma - 0.703)**j
        
    return teta*t_star

def t_3b_ps(p,s):
    """
    Backwards equation T_3a(p,s), subregion 3b. eq. 7
    Gives temperature from p,s input.
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    s :     double  Entropy (J/kg/K)

    Returns
    -------
    double Temperature (K)

    """
    t_star = 860 #K
    p_star = 100 * 10**6 # Pa
    s_star = 5.3 * 10**3 # J/kg/K
    
    pi = p/p_star
    sigma = s/s_star
    
    #    1,    2,   3,   4,   5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [-12, -12, -12, -12, -8, -8, -8, -6, -6, -6, -5, -5, -5, -5, -5, -4, -3, -3, -2,  0,  2,  3,  4,  5,  6,  8,  12, 14]
    J = [ 1,   3,   4,   7,   0,  1,  3,  0,  2,  4,  0,  1,  2,  4,  6,  12, 1,  6,  2,  0,  1,  1,  0,  24, 0,  3,  1,  2]
    n_array = [ 0.527111701601660,
               -0.401317830052742*10**2,
                0.153020073134484*10**3,
               -0.224799398218827*10**4,
               -0.193993484669048,
               -0.140467557893768*10**1,
                0.426799878114024*10**2,
                0.752810643416743,
                0.226657238616417*10**2,
               -0.622873556909932*10**3,
               -0.660823667935396,
                0.841267087271658,
               -0.253717501764397*10**2,
                0.485708963532948*10**3,
                0.880531517490555*10**3,
                0.265015592794626*10**7,
               -0.359287150025783,
               -0.656991567673753*10**3,
                0.241768149185367*10**1,
                0.856873461222588,
                0.655143675313458,
               -0.213535213206406,
                0.562974957606348*10**-2,
               -0.316955725450471*10**15,
               -0.699997000152457*10**-3,
                0.119845803210767*10**-1,
                0.193848122022095*10**-4,
               -0.215095749182309*10**-4]
    
    teta = 0
    for i,j,n in zip(I, J, n_array):
        teta += n * (pi + 0.760)**i * (sigma - 0.818)**j
        
    return teta*t_star

def v_3a_ps(p,s):
    """
    Backwards equation v_3a(p,s), subregion 3a. eq. 8
    Gives temperature from p,s input.
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    s :     double  Entropy (J/kg/K)

    Returns
    -------
    double Specific volume (m^3/kg)

    """
    
    v_star = 0.0028 #m^3/kg
    p_star = 100 * 10**6 # Pa
    s_star = 4.4 * 10**3 # J/kg
    
    pi = p/p_star
    sigma = s/s_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [-12, -12, -12, -10, -10, -10, -10, -8, -8, -8, -8, -6, -5, -4, -3, -3, -2, -2, -1, -1, 0, 0, 0, 1, 2, 4, 5, 6]
    J = [ 10,  12,  14,  4,   8,   10,  20,  5,  6,  14, 16, 28, 1,  5,  2,  4,  3,  8,  1,  2, 0, 1, 3, 0, 0, 2, 2, 0]
    n_array = [ 0.795544074093975*10**2,
               -0.238261242984590*10**4,
                0.176813100617787*10**5,
               -0.110524727080379*10**-2,
               -0.153213833655326*10**2,
                0.297544599376982*10**3,
               -0.350315206871242*10**8,
                0.277513761062119,
               -0.523964271036888,
               -0.148011182995403*10**6,
                0.160014899374266*10**7,
                0.170802322663427*10**13,
                0.246866996006494*10**-3,
                0.165326084797980*10**1,
               -0.118008384666987,
                0.253798642355900*10**1,
                0.965127704669424,
               -0.282172420532826*10**2,
                0.203224612353823,
                0.110648186063513*10**1,
                0.526127948451280,
                0.277000018736321,
                0.108153340501132*10**1,
               -0.744127885357893*10**-1,
                0.164094443541384*10**-1,
               -0.680468275301065*10**-1,
                0.257988576101640*10**-1,
               -0.145749861944416*10**-3]
    
    teta = 0
    for i,j,n in zip(I, J, n_array):
        teta += n * (pi + 0.187)**i * (sigma - 0.755)**j
        
    return teta*v_star

def v_3b_ps(p,s):
    """
    Backwards equation v_3b(p,s), subregion 3b. eq. 9
    Gives temperature from p,s input.
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    s :     double  Entropy (J/kg/K)

    Returns
    -------
    double Specific volume (m^3/kg)

    """
    
    v_star = 0.0088 #m^3/kg
    p_star = 100 * 10**6 # Pa
    s_star = 5.3 * 10**3 # J/kg
    
    pi = p/p_star
    sigma = s/s_star
    
    #    1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [-12, -12, -12, -12, -12, -12, -10, -10, -10, -10, -8, -5, -5, -5, -4, -4, -4, -4, -3, -2, -2, -2, -2, -2, -2, 0, 0, 0, 1, 1, 2]
    J = [ 0,   1,   2,   3,   5,   6,   0,   1,   2,   4,   0,  1,  2,  3,  0,  1,  2,  3,  1,  0,  1,  2,  3,  4,  12,0, 1, 2, 0, 2, 2]
    n_array = [ 0.591599780322238*10**-4,
               -0.185465997137856*10**-2,
                0.104190510480013*10**-1,
                0.598647302038590*10**-2,
               -0.771391189901699,
                0.172549765557036*10**1,
               -0.467076079846526*10**-3,
                0.134533823384439*10**-1,
               -0.808094336805495*10**-1,
                0.508139374365767,
                0.128584643361683*10**-2,
               -0.163899353915435*10**1,
                0.586938199318063*10**1,
               -0.292466667918613*10**1,
               -0.614076301499537*10**-2,
                0.576199014049172*10**1,
               -0.121613320606788*10**2,
                0.167637540957944*10**1,
               -0.744135838773463*10**1,
                0.378168091437659*10**-1,
                0.401432203027688*10**1,
                0.160279837479185*10**2,
                0.317848779347728*10**1,
               -0.358362310304853*10**1,
               -0.115995260446827*10**7,
                0.199256573577909,
               -0.122270624794624,
               -0.191449143716586*10**2,
               -0.150448002905284*10**-1,
                0.146407900162154*10**2,
               -0.327477787188230*10**1]
    
    teta = 0
    for i,j,n in zip(I, J, n_array):
        teta += n * (pi + 0.298)**i * (sigma - 0.816)**j
        
    return teta*v_star

def p_3sat_h(h):
    """
    Determines region boundary between region 3 and two-phase region 4.
    For a given enthalpy this function returns the pressure above which is 
    defined as region 3 and below region 4
    
    Eq. 10

    Parameters
    ----------
    h :     double  Enthalpy (J/kg)

    Returns
    -------
    double Pressure (Pa)

    """
    p_star = 22 * 10**6     # Pa
    h_star = 2600 * 10**3   # J/kg
    
    eta = h / h_star
    
    #   1,  2, 3, 4, 5,  6, 7, 8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [0, 1, 1, 1, 1,  5, 7, 8,  14, 20, 22, 24, 28, 36]
    J = [0, 1, 3, 4, 36, 3, 0, 24, 16, 16, 3,  18, 8,  24]
    n_array = [ 0.600073641753024,
               -0.936203654849857*10**1,
                0.246590798594147*10**2,
               -0.107014222858224*10**3,
               -0.915821315805768*10**14,
               -0.862332011700662*10**4,
               -0.235837344740032*10**2,
                0.252304969384128*10**18,
               -0.389718771997719*10**19,
               -0.333775713645296*10**23,
                0.356499469636328*10**11,
               -0.148547544720641*10**27,
                0.330611514838798*10**19,
                0.813641294467829*10**38]
    
    pi = 0
    for i,j,n in zip(I, J, n_array):
        pi += n * (eta - 1.02)**i * (eta - 0.608)**j
        
    return pi*p_star

def p_3sat_s(s):
    """
    Determines region boundary between region 3 and two-phase region 4.
    For a given entropy this function returns the pressure above which is 
    defined as region 3 and below region 4
    
    Eq. 11

    Parameters
    ----------
    s :     double  Entropy (J/kg/K)

    Returns
    -------
    double Pressure (Pa)

    """
    p_star = 22 * 10**6     # Pa
    s_star = 5.2 * 10**3   # J/kg
    
    sigma = s / s_star
    
    #   1,  2, 3, 4, 5,  6, 7, 8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31
    I = [0, 1, 1, 4, 12, 12, 16, 24, 28, 32]
    J = [0, 1, 32, 7, 4, 14, 36, 10, 0, 18]
    n_array = [ 0.639767553612785,
               -0.129727445396014*10**2,
               -0.224595125848403*10**16,
                0.177466741801846*10**7,
                0.717079349571538*10**10,
               -0.378829107169011*10**18,
               -0.955586736431328*10**35,
                0.187269814676188*10**24,
                0.119254746466473*10**12,
                0.110649277244882*10**37]
    
    pi = 0
    for i,j,n in zip(I, J, n_array):
        pi += n * (sigma - 1.03)**i * (sigma - 0.699)**j
        
    return pi*p_star

def find_region_ph(p,h):
    """
    find correct sub-region, 3a or 3b given p,h
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    h :     double  Enthalpy (J/kg)

    Returns
    -------
    str region (-)

    """
    
    p_crit = global_property.p_crit             # Critical pressure
    h_crit = global_property.h_crit             # J/kg, critical enthalpy, from IAPWS SR4-04(2014)
    p_r3_upper_boundary = global_property.p_r123_upper_boundary
    p_r3_lower_boundary = global_property.p_r3_lower_boundary
    
    region = None
    if p < p_crit and p >= p_r3_lower_boundary :
        # For pressures lower than the critical pressure, 
        # boundary defined as saturation line
        #print(f'p_3sat_h = {p_3sat_h(h)*10**-6} (MPa)')
        #
        if h <= h_crit and p >= p_3sat_h(h):
            # h>h_crit -> enthalpy is a saturated liquid, x = 0
            # for pressures above p_sat, region a (liquid)
            region = '3a'
        elif h >= h_crit and p >= p_3sat_h(h):
            # h>h_crit -> enthalpy is a saturated steam, x = 1
            # for pressures above p_sat, region b
            # region b
            region = '3b'
        else:
            # pressure in the saturation region, wrong main region, 
            # should be region 4 not 3. give error
            print('Error, saturation conditions for p,h given in t_3_ph. region 4 should be used')
        
    elif p >= p_crit and p <= p_r3_upper_boundary:
        # For pressures higher than the critical pressure, 
        # boundary defined as critical isentropic line estimated by h_3ab
        #print(f'h_3ab(p) = {h_3ab(p)*10**-3} (kJ/kg)')
        if h <= h_3ab(p):
            # region a
            # lower h -> lower temp -> region a
            region = '3a'
        else:
            # region b
            region = '3b'   
            
    return region
    
def find_region_ps(p,s):
    """
    find correct sub-region, 3a or 3b given p,s
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    s :     double  Entropy (J/kg/K)

    Returns
    -------
    str region (-)

    """
    
    p_crit = global_property.p_crit             # Critical pressure
    s_crit = global_property.s_crit             #J/kg/K
    
    p_r3_upper_boundary = global_property.p_r123_upper_boundary
    p_r3_lower_boundary = global_property.p_r3_lower_boundary
    
    region = None
    if p < p_crit and p >= p_r3_lower_boundary :
        # For pressures lower than the critical pressure, 
        # boundary defined as saturation line
        #print(f'p_3sat_s = {p_3sat_s(s)*10**-6} (MPa)')
        #
        if s < s_crit and p > p_3sat_s(s):
            # s>s_crit -> entropy is a saturated liquid, x = 0
            # for pressures above p_sat, region a (liquid)
            region = '3a'
        elif s > s_crit and p > p_3sat_s(s):
            # s>s_crit -> entropy is a saturated steam, x = 1
            # for pressures above p_sat, region b
            # region b
            region = '3b'
        else:
            # pressure in the saturation region, wrong main region, 
            # should be region 4 not 3. give error
            print('Error, saturation conditions for p,h given in t_3_ph. region 4 should be used')
    elif p >= p_crit and p <= p_r3_upper_boundary:
        # For pressures higher than the critical pressure, 
        # boundary defined as critical isentropic line
        if s <= s_crit:
            # region a
            region = '3a'
        else:
            # region b
            region = '3b'
            
    return region


def t_3_ph(p,h):
    """
    Backwards equation T_3(p,h), identifies correct sub-region, a or b
    and uses the corresponding eq. to return a temperature.
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    h :     double  Enthalpy (J/kg)

    Returns
    -------
    double Temperature (K)

    """
    region = find_region_ph(p, h)
    
    t = None
    if region == '3a':
        # region a
        t = t_3a_ph(p,h)
    elif region == '3b':
        # region b
        t = t_3b_ph(p,h)
            
    return t

def v_3_ph(p,h):
    """
    Backwards equation v_3(p,h), identifies correct sub-region, a or b
    and uses the corresponding eq.
    Gives specific volume from p,h input.
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    h :     double  Enthalpy (J/kg)

    Returns
    -------
    double Specific volume (m^3/kg)
    """
    
    region = find_region_ph(p, h)
    
    v = None
    if region == '3a':
        # region a
        v = v_3a_ph(p,h)
    elif region == '3b':
        # region b
        v = v_3b_ph(p,h)
            
    return v
    
def t_3_ps(p,s):
    """
    Backwards equation T_3(p,s), identifies correct sub-region, a or b
    and uses the corresponding eq. to return a temperature.
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    s :     double  Entropy (J/kg/K)

    Returns
    -------
    double Temperature (K)

    """
    
    region = find_region_ps(p, s)
    
    t = None
    if region == '3a':
        # region a
        t = t_3a_ps(p,s)
    elif region == '3b':
        # region b
        t = t_3b_ps(p,s)
            
    return t

def v_3_ps(p,s):
    """
    Backwards equation v_3(p,s), identifies correct sub-region, a or b
    and uses the corresponding eq.
    Gives specific volume from p,s input.
    
    Parameters
    ----------
    p:      double  Pressure (Pa).
    h :     double  Entropy (J/kg/K)

    Returns
    -------
    double Specific volume (m^3/kg)
    """
    
    region = find_region_ps(p, s)
    
    v = None
    if region == '3a':
        # region a
        v = v_3a_ps(p,s)
    elif region == '3b':
        # region b
        v = v_3b_ps(p,s)
            
    return v

if __name__ == "__main__":
    main()
