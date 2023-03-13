# -*- coding: utf-8 -*-
"""
Available unit systems for pyfluidproperties

@author: Christoffer Rappmann, christoffer.rappmann@gmail.com
"""

class unit_sys:
    """
    Class containing all available unit systems for the pyfluidproperties package.
    Use when specifing unit systems for easier typing.
    """
    
    # Name of fluid        # ID         # Name
    SI =                   'SI'         # SI Units
    US =                   'US'         # Imperial units / US
    SI_bar_kJ =            'SI_bar_kJ'  # SI units with bar(a) for pressure, degC for temprature and kJ for energy
    SI_MPa_kJ =            'SI_MPa_kJ'  # SI units with MPa for pressure, degC for temperature and kJ for energy