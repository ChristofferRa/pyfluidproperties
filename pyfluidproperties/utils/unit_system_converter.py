# -*- coding: utf-8 -*-
"""
Unit conversion for pyfluidproperties

Converts property to and from SI-units used by pyfluidproperties

List of functions and code structure for unit_conv class:
    __init__(self, unit_system)
    get_unit_system()
    get_units(self)
    convert(self,value,**kwargs)
    set_units(self)
    set_conversion_factors(slef)
    get_var_conversion_factors(self,prop,x)
    

@author: Christoffer Rappmann, christoffer.rappmann@gmail.com
"""

class unit_conv():
    """
    Unit conversions
    
    Converts units l,m,e p,T,v,h,u,s,cp/cv,w,my,tc,st from chosen unit-system to si 
    and vice versa.
    
    """
    
    def __init__(self, unit_system):
        """
        Initilize unit_conv using chosen unit-system.

        Parameters
        ----------
        unit_system :   string      unit system

        Returns
        -------
        None.

        """
        
        # Save unit-system to use
        self.unit_system = unit_system
        
        # Get units and conversion factors given unit-system
        unit_conv.set_units(self)
        unit_conv.set_conversion_factors(self)
        
    @staticmethod
    def get_unit_systems():
        """
        Return list with implemented unit-systems

        Returns
        -------
        unit_systems_list : list of strings    list of unit systems

        """
        unit_system_list = ['SI', 'US', 'SI_bar_kJ', 'SI_MPa_kJ']
        return unit_system_list
        
    def get_units(self):
        """
        Returns dict with used units for choosen unit system

        Returns
        -------
        Dict        used units.

        """
            
        return self.units
           
    def convert(self,value,**kwargs):
        """
        Convert submitted value to or from SI-units using choosen unit system

        Parameters
        ----------
        value :     double          value to convert.
        prop =:     double          which property to convert
        to_from=:   string          "to_SI" or "From_SI"
        
        Returns
        -------
        converted_unit : double new value after conversion

        """
        # Set defaults
        to_from = 'to_SI' # Convert to SI-units as deafult
        prop = ''
        # Read  and set kwargs, ignore onces not defined
        for key, arg in kwargs.items():
            if key == 'prop':
                prop = arg
            elif key == 'to_from':
                to_from = arg
                
        if prop in self.conversion_factors or prop == 'T':
            # Convert to SI units from choosen unit system
            if to_from == 'to_SI':
                # Check if property is temperature, handled differently
                if prop == 'T':
                    converted_unit = unit_conv.get_var_conversion_factors(self, 'T_to', value)
                    # do conversion 
                else:
                    converted_unit = self.conversion_factors[prop]*value
                    # Convert to SI units from choosen unit system
            elif to_from == 'from_SI':
                # Check if property is temperature, handled differently
                if prop == 'T':
                    converted_unit = unit_conv.get_var_conversion_factors(self, 'T_from', value)
                # do conversion 
                else:
                    converted_unit = 1/self.conversion_factors[prop]*value
            else:
                raise Exception(f'"{to_from}" is not a valid input, choose "from_SI" or "to_SI"')
        else:
            raise Exception(f'"{prop}" is not a valid input, valid properties: {self.conversion_factors.keys()}')
            
        return converted_unit
    
    def set_units(self):
        """
        Set which units to use given chosen unit system chosen.

        Returns
        -------
        None.

        """
        units = {'SI': {'l': 'm',
                          'm': 'kg',
                          'E': 'J',
                          'p': 'Pa',
                          'T': 'K',
                          'v': 'm^3/kg',
                          'h': 'J/kg',
                          'u': 'J/kg',
                          's': 'J/kg/C',
                          'cp': 'J/kg/C',
                          'cv': 'J/kg/C',
                          'w': 'm/s',
                          'rho': 'kg/m^3',
                          'x': '-',
                          'my': 'Pa*s',
                          'tc': 'W/m/K',
                          'st': 'N/m'}, 
               'US': {'l': 'inch',
                          'm': 'lbs',
                          'E': 'Btu',
                          'p': 'Psi',
                          'T': 'degF',
                          'v': 'ft^3/lbs',
                          'h': 'Btu/lbm',
                          'u': 'Btu/lbm',
                          's': 'Btu/lbm/R',
                          'cp': 'Btu/lbm/R',
                          'cv': 'Btu/lbm/R',
                          'w': 'ft/s',
                          'rho': 'lbs/ft^3',
                          'x': '-',
                          'my': 'lb*s/ft^2',
                          'tc': 'Btu/h/ft/F',
                          'st': 'lb/inch'}, 
               'SI_bar_kJ': {'l': 'm',
                          'm': 'kg',
                          'E': 'kJ',
                          'p': 'bar(a)',
                          'T': 'degC',
                          'v': 'm^3/kg',
                          'h': 'kJ/kg',
                          'u': 'kJ/kg',
                          's': 'kJ/kg/C',
                          'cp': 'kJ/kg/C',
                          'cv': 'kJ/kg/C',
                          'w': 'm/s',
                          'rho': 'kg/m^3',
                          'x': '-',
                          'my': 'Pa*s',
                          'tc': 'W/m/K',
                          'st': 'N/m'}, 
               'SI_MPa_kJ': {'l': 'm',
                          'm': 'kg',
                          'E': 'kJ',
                          'p': 'MPa',
                          'T': 'degC',
                          'v': 'm^3/kg',
                          'h': 'kJ/kg',
                          'u': 'kJ/kg',
                          's': 'kJ/kg/C',
                          'cp': 'kJ/kg/C',
                          'cv': 'kJ/kg/C',
                          'w': 'm/s',
                          'rho': 'kg/m^3',
                          'x': '-',
                          'my': 'Pa*s',
                          'tc': 'W/m/K',
                          'st': 'N/m'}}
        
        self.units = units[self.unit_system]
        
    def set_conversion_factors(self):
        """
        Set which conversion factors to use given unit_system choosen

        Returns
        -------
        None.

        """
        conversion_factors = {'SI': {'l':       1,                          # m
                                     'm':       1,                          # kg
                                     'E':       1,                          # J
                                     'p':       1,                          # Pa
                                     'v':       1,                          # m^3/kg
                                     'h':       1,                          # J/kg
                                     'u':       1,                          # J/kg
                                     's':       1,                          # J/kg/C
                                     'cp':      1,                          # J/kg/C
                                     'cv':      1,                          # J/kg/C
                                     'w':       1,                          # m/s
                                     'rho':     1,                          # kg/m^3
                                     'x':       1,                          # -
                                     'my':      1,                          # Pa s
                                     'tc':      1,                          # W/m/K
                                     'st':      1},                         # N/m
                              'US': {'l':       0.0254,                     # inch
                                     'm':       0.45359237,                 # lbs
                                     'E':       1055.0558526,               # Btu
                                     'p':       0.0689475729*1e5,           # Psi
                                     'v':       0.3048**3/0.45359237,       # ft^3/lbs
                                     'h':       1055.0558526/0.45359237,         # Btu/lbm
                                     'u':       1055.0558526/0.45359237,         # Btu/lbm
                                     's':       1055.0558526/0.45359237/(5/9),   # Btu/lbm/R
                                     'cp':      1055.0558526/0.45359237/(5/9),   # Btu/lbm/R
                                     'cv':      1055.0558526/0.45359237/(5/9),   # Btu/lbm/R
                                     'w':       0.3048,                     # ft/s
                                     'rho':     1/0.0624279606,               # lbs/ft^3
                                     'x':       1,                          # -
                                     'my':      1/0.0208854342,               # lb s / ft^2
                                     'tc':      1/0.5777893165,               # Btu/h/ft/F
                                     'st':      1/0.0057101471},              # lb/inch
                              'SI_bar_kJ': {'l': 1,                         # m
                                     'm':       1,                          # kg
                                     'E':       1e3,                        # kJ 
                                     'p':       1e5,                        # bar
                                     'v':       1,                          # m^3/kg
                                     'h':       1e3,                        # kJ/kg
                                     'u':       1e3,                        # kJ/kg
                                     's':       1e3,                        # kJ/kg/C
                                     'cp':      1e3,                        # kJ/kg/C
                                     'cv':      1e3,                        # kJ/kg/C
                                     'w':       1,                          # m/s
                                     'rho':     1,                          # kg/m^3
                                     'x':       1,                          # -
                                     'my':      1,                          # Pa s
                                     'tc':      1,                          # W/m/K
                                     'st':      1},                          # N/m
                                'SI_MPa_kJ': {'l': 1,                         # m
                                     'm':       1,                          # kg
                                     'E':       1e3,                        # kJ 
                                     'p':       1e6,                        # bar
                                     'v':       1,                          # m^3/kg
                                     'h':       1e3,                        # kJ/kg
                                     'u':       1e3,                        # kJ/kg
                                     's':       1e3,                        # kJ/kg/C
                                     'cp':      1e3,                        # kJ/kg/C
                                     'cv':      1e3,                        # kJ/kg/C
                                     'w':       1,                          # m/s
                                     'rho':     1,                          # kg/m^3
                                     'x':       1,                          # -
                                     'my':      1,                          # Pa s
                                     'tc':      1,                          # W/m/K
                                     'st':      1}}                         # N/m
        
        self.conversion_factors = conversion_factors[self.unit_system]
            
    def get_var_conversion_factors(self, prop, x):
        """
        Return conversion for variable conversion factors

        Parameters
        ----------
        prop :  string   property to get conversion for .
        x :     double   value to convert.
        
        Returns
        -------
        double  converted property value

        """
        variable_conversion_factors = {'T_to': {'SI':           x,                  # K
                                                'US':           (x-32)*5/9+273.15,  # F
                                                'SI_bar_kJ':    x+273.15,           # C
                                                'SI_MPa_kJ':    x+273.15},          # C
                                       'T_from': {'SI':         x,                  # K
                                                  'US':         ((x-273.15)*9/5)+32,# F
                                                  'SI_bar_kJ':  x-273.15,           # C
                                                  'SI_MPa_kJ':  x-273.15},          # C
                                               }
        return variable_conversion_factors[prop][self.unit_system]
        
if __name__ == "__main__":
    main()
