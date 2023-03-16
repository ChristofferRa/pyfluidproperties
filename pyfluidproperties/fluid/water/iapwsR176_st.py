# -*- coding: utf-8 -*-
"""
Module containing functions implemented from equations presentend in the R1-76
Revised Release on Surface Tension of Ordinary Water Substance.

Surface tension in the interface between liquid av vapor phase of water

List of functions and code structure:
    region_st_t(t)
    st_t(t)

@author: Christoffer Rappmann, christoffer.rappmann@gmail.com
"""

from .iapwsif97 import iapwsif_globals as global_property

def region_st_t(t):
    """
    Determines if given temperature is within the applicable region.
    There is only one region in iapwsf R1-76.

    Parameters
    ----------
    t :     double     temperature (K)

    Returns
    -------
    boolean,     within applicable region (-)
    """
    
    #       Critical temperature  (K)       Triple point temperature (K)
    if t <= global_property.T_crit and t >= global_property.Tt:
        region = True
    else:
        region = False
    
    return region


def st_t(t):
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
    
    # Reference constants
    tcrit      = global_property.T_crit    # K         reference temperature
    cap_b   = 235.8*10**-3               # N/m       reference surface tension
    b       = -0.625                    # -         -
    my      = 1.256                     # -         -
    
    # Dimensionless quantities
    tau = 1 - t/tcrit
    
    st = cap_b * tau**my * (1 + b * tau)
    
    return st

if __name__ == "__main__":
    main()