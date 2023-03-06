# -*- coding: utf-8 -*-
"""
Main class for pyfluidproperties library. 

Initialize fluid and unit-system.
Access corresponding fluid class and do conversions

List of functions and code structure for pyfluidproperties class:
    Init functions:
        __init__(self, **kwargs)
        init_fluid(self, p = 1.013*1e5, t = 20+273.15)
    
    Update All properties:
        update(self, **kwargs)
    
    Update Single properties:
        get_h(self, **kwargs)
        get_v(self, **kwargs)
        get_u(self, **kwargs)
        get_s(self, **kwargs)
        get_cp(self, **kwargs)
        get_cv(self, **kwargs)
        get_w(self, **kwargs)
        get_rho(self, **kwargs)
        get_my(self, **kwargs)
        get_tc(self, **kwargs)
        get_st(self, **kwargs)
        get_x(self, **kwargs)
        get_vx(self, **kwargs)
        get_T(self, **kwargs)
        get_p(self, **kwargs)
    
    Misc:
        set_properties(self,properties_dict)
        get_fluids()
        get_unit_systems()

@author: Christoffer Rappmann, christoffer.rappmann@gmail.com
"""

from .fluid.water import iapwsif97
from .unit_system_converter import unit_conv
import numpy as np

class properties:
    """
    Main class for pyfluidproperties library. 

    Initialize fluid and unit-system.
    Access corresponding fluid class and do conversions 
    """
    
    def __init__(self, **kwargs):
        """
        Initialize pyfluidproperties object. 
        
        Choose fluid and unit system

        Parameters
        ----------
        fluid :         string      fluid 
        unit_system :   string      unit system si/us

        Returns
        -------
        None.
        """
         
        # Set defaults
        unit_system = 'SI' # Convert to SI-units as deafult
        fluid = None
        
        # Read  and set kwargs, ignore onces not defined
        for key, arg in kwargs.items():
            if key == 'fluid':
                fluid = arg
            elif key == 'unit_system':
                unit_system = arg
        
        # Check for enough and valid inputs, otherwhise abort...
        if fluid not in properties.get_fluids():
            raise Exception(f'"{fluid}" is not a valid input for "fluid = "')
            
        if unit_system not in properties.get_unit_systems():
            raise Exception(f'"{unit_system}" is not a valid input for "unit_system ="')
        
        # If all good, set fluid and unit_system
        self.fluid = fluid
        self.unit_system = unit_system
        
        # Init converter class and fluid class
        self.unit_converter = unit_conv(unit_system)
        properties.init_fluid(self)
            
    
    def init_fluid(self, p = 1.013*1e5, t = 20+273.15):
        """
        Returning object for chosen fluid when initializing the object.
        Optional initiation properties, p-T

        Parameters
        ----------
        p : double, optional    Pressure (Pa), default 1 Atm
        T : double, optional    Temperature (K), defualt 20 degC

        Returns
        -------
        fluid_class     fluid properties class

        """
        
        if self.fluid == 'H2O':
            self.fluid_class = iapwsif97(p, t)
        elif self.fluid == 'N2':
            self.fluid_class = None
        else:
            print('No valid fluid specified')
            self.fluid_class = None
        
        # Set default properties
        properties.set_properties(self, {'p': self.fluid_class.p,
                                                'T': self.fluid_class.T,
                                                'v': self.fluid_class.v,
                                                'h': self.fluid_class.h,
                                                'u': self.fluid_class.u,
                                                's': self.fluid_class.s,
                                                'cp': self.fluid_class.cp,
                                                'cv': self.fluid_class.cv,
                                                'w': self.fluid_class.w,
                                                'rho': self.fluid_class.rho,
                                                'x': self.fluid_class.x,
                                                'my': self.fluid_class.my,
                                                'tc': self.fluid_class.tc,
                                                'st': self.fluid_class.st})
        
    def update(self, **kwargs):
        """
        Update function, converts and updates fluid class
        Saves converted properties to this class

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        T :     double,     Temperature (valid unit in selected unit system)
        h:      double      Entropy     (valid unit in selected unit system)
        s:      double      Entropy     (valid unit in selected unit system)
        x:      double      Mass-fraction steam (-).

        Returns
        -------
        None.

        """
        # init
        p = np.nan
        T = np.nan
        h = np.nan
        s = np.nan
        x = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'T':
                T = self.unit_converter.convert(arg, prop = 'T', to_from = 'to_SI')
            elif key == 'h':
                h = self.unit_converter.convert(arg, prop = 'h', to_from = 'to_SI')
            elif key == 's':
                s = self.unit_converter.convert(arg, prop = 's', to_from = 'to_SI')
            elif key == 'x':
                x = arg
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            if 'T' in kwargs.keys():
                # p-T
                self.fluid_class.update_pt(p, T)
            elif 'h' in kwargs.keys():
                # p-h
                self.fluid_class.update_ph(p, h)
            elif 's' in kwargs.keys():
                # p-s
                self.fluid_class.update_ps(p, s)
            elif 'x' in kwargs.keys():
                # p-x
                self.fluid_class.update_px(p, x)
        elif 'h' in kwargs.keys() and 's' in kwargs.keys():
            # h-s
            self.fluid_class.update_hs(h, s)
        elif 'T' in kwargs.keys() and 'x' in kwargs.keys():
            # T-x
            self.fluid_class.update_tx(T,x)
        else:
            raise Exception('Not enough input arguments or inputarguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        
        # Convert fluid-class properties from SI to specified units system
        properties.set_properties(self, {'p': self.fluid_class.p,
                                                'T': self.fluid_class.T,
                                                'v': self.fluid_class.v,
                                                'h': self.fluid_class.h,
                                                'u': self.fluid_class.u,
                                                's': self.fluid_class.s,
                                                'cp': self.fluid_class.cp,
                                                'cv': self.fluid_class.cv,
                                                'w': self.fluid_class.w,
                                                'rho': self.fluid_class.rho,
                                                'x': self.fluid_class.x,
                                                'my': self.fluid_class.my,
                                                'tc': self.fluid_class.tc,
                                                'st': self.fluid_class.st})
           
    def get_h(self, **kwargs):
        """
        calls fluid class for enthalpy h and returns enthalpy in choosen units.

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        T :     double,     Temperature (valid unit in selected unit system)
        s:      double      Entropy     (valid unit in selected unit system)
        rho:    double      Density     (valid unit in selected unit system)
        x:      double      Mass-fraction steam (-)

        Returns
        -------
        Enthalpy (valid unit in selected unit system).

        """
        # init
        p = np.nan
        T = np.nan
        h = np.nan
        s = np.nan
        x = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'T':
                T = self.unit_converter.convert(arg, prop = 'T', to_from = 'to_SI')
            elif key == 's':
                s = self.unit_converter.convert(arg, prop = 's', to_from = 'to_SI')
            elif key == 'rho':
                rho = self.unit_converter.convert(arg, prop = 'rho', to_from = 'to_SI')
            elif key == 'x':
                x = arg
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: p, T, s, rho')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            if 'T' in kwargs.keys():
                # p-T
                h = self.fluid_class.h_pt(p, T)
            elif 's' in kwargs.keys():
                # p-s
                h = self.fluid_class.h_ps(p, s)
            elif 'x' in kwargs.keys():
                # p-x
                h = self.fluid_class.h_px(p, x)
            elif 'rho' in kwargs.keys():
                # p-rho
                h = self.fluid_class.h_prho(p, rho)
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-s, p-rho, p-x, T-x')
        elif 'T' in kwargs.keys() and 'x' in kwargs.keys():
            # T-x
            h = self.fluid_class.h_tx(T, x)
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-s, p-rho, p-x, T-x')
        
        # Convert from SI to specified units system
        return self.unit_converter.convert(h,prop = 'h', to_from = 'from_SI')
       
    def get_v(self, **kwargs):
        """
        calls fluid class for specific volume and returns results in choosen units.

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        T :     double,     Temperature (valid unit in selected unit system)
        h:      double      Entropy     (valid unit in selected unit system)
        s:      double      Entropy     (valid unit in selected unit system)
        x:      double      Mass-fraction steam (-)

        Returns
        -------
        Specific volume (valid unit in selected unit system).

        """
        # init
        p = np.nan
        T = np.nan
        h = np.nan
        s = np.nan
        x = np.nan
        
        v = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'T':
                T = self.unit_converter.convert(arg, prop = 'T', to_from = 'to_SI')
            elif key == 'h':
                h = self.unit_converter.convert(arg, prop = 'h', to_from = 'to_SI')
            elif key == 's':
                s = self.unit_converter.convert(arg, prop = 's', to_from = 'to_SI')
            elif key == 'x':
                x = arg
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: p, T, h, s')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            if 'T' in kwargs.keys():
                # p-T
                v = self.fluid_class.v_pt(p, T)
            elif 'h' in kwargs.keys():
                # p-h
                v = self.fluid_class.v_ph(p, h)
            elif 's' in kwargs.keys():
                # p-s
                v = self.fluid_class.v_ps(p, s)
            elif 'x' in kwargs.keys():
                # p-x
                v = self.fluid_class.v_px(p, x)
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        elif 'h' in kwargs.keys() and 's' in kwargs.keys():
            # h-s
            v = self.fluid_class.v_hs(h, s)
        elif 'T' in kwargs.keys() and 'x' in kwargs.keys():
            # T-x
            v = self.fluid_class.v_tx(T, x)
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        
        # Convert from SI to specified units system
        return self.unit_converter.convert(v,prop = 'v', to_from = 'from_SI')
    
    def get_u(self, **kwargs):
        """
        calls fluid class for specific internal energy and returns results in choosen units.

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        T :     double,     Temperature (valid unit in selected unit system)
        h:      double      Entropy     (valid unit in selected unit system)
        s:      double      Entropy     (valid unit in selected unit system)
        x:      double      Mass-fraction steam (-)

        Returns
        -------
        Specific internal energy (valid unit in selected unit system).

        """
        # init
        p = np.nan
        T = np.nan
        h = np.nan
        s = np.nan
        x = np.nan
        
        u = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'T':
                T = self.unit_converter.convert(arg, prop = 'T', to_from = 'to_SI')
            elif key == 'h':
                h = self.unit_converter.convert(arg, prop = 'h', to_from = 'to_SI')
            elif key == 's':
                s = self.unit_converter.convert(arg, prop = 's', to_from = 'to_SI')
            elif key == 'x':
                x = arg
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: p, T, h, s, x')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            if 'T' in kwargs.keys():
                # p-T
                u = self.fluid_class.u_pt(p, T)
            elif 'h' in kwargs.keys():
                # p-h
                u = self.fluid_class.u_ph(p, h)
            elif 's' in kwargs.keys():
                # p-s
                u = self.fluid_class.u_ps(p, s)
            elif 'x' in kwargs.keys():
                # p-x
                u = self.fluid_class.u_px(p, x)
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        elif 'h' in kwargs.keys() and 's' in kwargs.keys():
            # h-s
            u = self.fluid_class.u_hs(h, s)
        elif 'T' in kwargs.keys() and 'x' in kwargs.keys():
            # T-x
            u = self.fluid_class.u_tx(T, x)
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        
        # Convert from SI to specified units system
        return self.unit_converter.convert(u,prop = 'u', to_from = 'from_SI')
    
    def get_s(self, **kwargs):
        """
        calls fluid class for entropy and returns result in choosen units.

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        T :     double,     Temperature (valid unit in selected unit system)
        h:      double      Entropy     (valid unit in selected unit system)
        x:      double      Mass-fraction steam (-)

        Returns
        -------
        Entropy (valid unit in selected unit system).

        """
        # init
        p = np.nan
        T = np.nan
        h = np.nan
        s = np.nan
        x = np.nan
                
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'T':
                T = self.unit_converter.convert(arg, prop = 'T', to_from = 'to_SI')
            elif key == 'h':
                h = self.unit_converter.convert(arg, prop = 'h', to_from = 'to_SI')
            elif key == 'x':
                x = arg
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: p, T, h, rho')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            if 'T' in kwargs.keys():
                # p-T
                s = self.fluid_class.s_pt(p, T)
            elif 'h' in kwargs.keys():
                # p-h
                s = self.fluid_class.s_ph(p, h)
            elif 'x' in kwargs.keys():
                # p-x
                s = self.fluid_class.s_px(p, x)
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, h-s, p-x, T-x')
        elif 'T' in kwargs.keys() and 'x' in kwargs.keys():
            # T-x
            s = self.fluid_class.s_tx(T, x)
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, h-s, p-x, T-x')
        
        # Convert from SI to specified units system
        return self.unit_converter.convert(s,prop = 's', to_from = 'from_SI')
    
    def get_cp(self, **kwargs):
        """
        calls fluid class for isobaric heat capcity and returns results in choosen units.

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        T :     double,     Temperature (valid unit in selected unit system)
        h:      double      Entropy     (valid unit in selected unit system)
        s:      double      Entropy     (valid unit in selected unit system)
        x:      double      Mass-fraction steam (-)

        Returns
        -------
        Isobaric heat capacity (valid unit in selected unit system).

        """
        # init
        p = np.nan
        T = np.nan
        h = np.nan
        s = np.nan
        x = np.nan
        
        cp = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'T':
                T = self.unit_converter.convert(arg, prop = 'T', to_from = 'to_SI')
            elif key == 'h':
                h = self.unit_converter.convert(arg, prop = 'h', to_from = 'to_SI')
            elif key == 's':
                s = self.unit_converter.convert(arg, prop = 's', to_from = 'to_SI')
            elif key == 'x':
                x = arg
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: p, T, h, s, x')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            if 'T' in kwargs.keys():
                # p-T
                cp = self.fluid_class.cp_pt(p, T)
            elif 'h' in kwargs.keys():
                # p-h
                cp = self.fluid_class.cp_ph(p, h)
            elif 's' in kwargs.keys():
                # p-s
                cp = self.fluid_class.cp_ps(p, s)
            elif 'x' in kwargs.keys():
                # p-x
                if x == 0.0:
                    cp = self.fluid_class.cpL_p(p)
                elif x == 1.0:
                    cp = self.fluid_class.cpV_p(p)
                else:
                    cp = np.nan
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        elif 'h' in kwargs.keys() and 's' in kwargs.keys():
            # h-s
            cp = self.fluid_class.cp_hs(h, s)
        elif 'T' in kwargs.keys() and 'x' in kwargs.keys():
            # T-x
            if x == 0.0:
                cp = self.fluid_class.cpL_T(T)
            elif x == 1.0:
                cp = self.fluid_class.cpV_T(T)
            else:
                cp = np.nan
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        
        # Convert from SI to specified units system
        return self.unit_converter.convert(cp,prop = 'cp', to_from = 'from_SI')
    
    def get_cv(self, **kwargs):
        """
        calls fluid class for isochoric heat capcity and returns results in choosen units.

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        T :     double,     Temperature (valid unit in selected unit system)
        h:      double      Entropy     (valid unit in selected unit system)
        s:      double      Entropy     (valid unit in selected unit system)
        x:      double      Mass-fraction steam (-)

        Returns
        -------
        Isochoric heat capacity (valid unit in selected unit system).

        """
        # init
        p = np.nan
        T = np.nan
        h = np.nan
        s = np.nan
        x = np.nan
        
        cv = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'T':
                T = self.unit_converter.convert(arg, prop = 'T', to_from = 'to_SI')
            elif key == 'h':
                h = self.unit_converter.convert(arg, prop = 'h', to_from = 'to_SI')
            elif key == 's':
                s = self.unit_converter.convert(arg, prop = 's', to_from = 'to_SI')
            elif key == 'x':
                x = arg
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: p, T, h, s, x')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            if 'T' in kwargs.keys():
                # p-T
                cv = self.fluid_class.cv_pt(p, T)
            elif 'h' in kwargs.keys():
                # p-h
                cv = self.fluid_class.cv_ph(p, h)
            elif 's' in kwargs.keys():
                # p-s
                cv = self.fluid_class.cv_ps(p, s)
            elif 'x' in kwargs.keys():
                # p-x
                if x == 0.0:
                    cv = self.fluid_class.cvL_p(p)
                elif x == 1.0:
                    cv = self.fluid_class.cvV_p(p)
                else:
                    cv = np.nan
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        elif 'h' in kwargs.keys() and 's' in kwargs.keys():
            # h-s
            cv = self.fluid_class.cv_hs(h, s)
        elif 'T' in kwargs.keys() and 'x' in kwargs.keys():
            # T-x
            if x == 0.0:
                cv = self.fluid_class.cvL_T(T)
            elif x == 1.0:
                cv = self.fluid_class.cvV_T(T)
            else:
                cv = np.nan
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        
        # Convert from SI to specified units system
        return self.unit_converter.convert(cv,prop = 'cv', to_from = 'from_SI')
    
    def get_w(self, **kwargs):
        """
        calls fluid class for speed of sound and returns results in choosen units.

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        T :     double,     Temperature (valid unit in selected unit system)
        h:      double      Entropy     (valid unit in selected unit system)
        s:      double      Entropy     (valid unit in selected unit system)
        x:      double      Mass-fraction steam (-)

        Returns
        -------
        None.

        """
        # init
        p = np.nan
        T = np.nan
        h = np.nan
        s = np.nan
        x = np.nan
        
        w = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'T':
                T = self.unit_converter.convert(arg, prop = 'T', to_from = 'to_SI')
            elif key == 'h':
                h = self.unit_converter.convert(arg, prop = 'h', to_from = 'to_SI')
            elif key == 's':
                s = self.unit_converter.convert(arg, prop = 's', to_from = 'to_SI')
            elif key == 'x':
                x = arg
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: p, T, h, s, x')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            if 'T' in kwargs.keys():
                # p-T
                w = self.fluid_class.w_pt(p, T)
            elif 'h' in kwargs.keys():
                # p-h
                w = self.fluid_class.w_ph(p, h)
            elif 's' in kwargs.keys():
                # p-s
                w = self.fluid_class.w_ps(p, s)
            elif 'x' in kwargs.keys():
                # p-x
                if x == 0.0:
                    w = self.fluid_class.wL_p(p)
                elif x == 1.0:
                    w = self.fluid_class.wV_p(p)
                else:
                    w = np.nan
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        elif 'h' in kwargs.keys() and 's' in kwargs.keys():
            # h-s
            w = self.fluid_class.w_hs(h, s)
        elif 'T' in kwargs.keys() and 'x' in kwargs.keys():
            # T-x
            if x == 0.0:
                w = self.fluid_class.wL_T(T)
            elif x == 1.0:
                w = self.fluid_class.wV_T(T)
            else:
                w = np.nan
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        
        # Convert from SI to specified units system
        return self.unit_converter.convert(w,prop = 'w', to_from = 'from_SI')
    
    def get_rho(self, **kwargs):
        """
        calls fluid class for density and returns results in choosen units.

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        T :     double,     Temperature (valid unit in selected unit system)
        h:      double      Entropy     (valid unit in selected unit system)
        s:      double      Entropy     (valid unit in selected unit system)
        x:      double      Mass-fraction steam (-)

        Returns
        -------
        None.

        """
        # init
        p = np.nan
        T = np.nan
        h = np.nan
        s = np.nan
        x = np.nan
        
        rho = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'T':
                T = self.unit_converter.convert(arg, prop = 'T', to_from = 'to_SI')
            elif key == 'h':
                h = self.unit_converter.convert(arg, prop = 'h', to_from = 'to_SI')
            elif key == 's':
                s = self.unit_converter.convert(arg, prop = 's', to_from = 'to_SI')
            elif key == 'x':
                x = arg
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: p, T, h, s, x')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            if 'T' in kwargs.keys():
                # p-T
                rho = self.fluid_class.rho_pt(p, T)
            elif 'h' in kwargs.keys():
                # p-h
                rho = self.fluid_class.rho_ph(p, h)
            elif 's' in kwargs.keys():
                # p-s
                rho = self.fluid_class.rho_ps(p, s)
            elif 'x' in kwargs.keys():
                # p-x
                rho = self.fluid_class.rho_px(p, x)
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        elif 'h' in kwargs.keys() and 's' in kwargs.keys():
            # h-s
            rho = self.fluid_class.rho_hs(h, s)
        elif 'T' in kwargs.keys() and 'x' in kwargs.keys():
            # T-x
            rho = self.fluid_class.rho_tx(T, x)
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        
        # Convert from SI to specified units system
        return self.unit_converter.convert(rho,prop = 'rho', to_from = 'from_SI')
    
    
    def get_my(self, **kwargs):
        """
        calls fluid class for dynamic viscosity and returns results in choosen units.

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        T :     double,     Temperature (valid unit in selected unit system)
        h:      double      Entropy     (valid unit in selected unit system)
        s:      double      Entropy     (valid unit in selected unit system)
        x:      double      Mass-fraction steam (-)

        Returns
        -------
        None.

        """
        # init
        p = np.nan
        T = np.nan
        h = np.nan
        s = np.nan
        x = np.nan
        
        my = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'T':
                T = self.unit_converter.convert(arg, prop = 'T', to_from = 'to_SI')
            elif key == 'h':
                h = self.unit_converter.convert(arg, prop = 'h', to_from = 'to_SI')
            elif key == 's':
                s = self.unit_converter.convert(arg, prop = 's', to_from = 'to_SI')
            elif key == 'x':
                x = arg
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: p, T, h, s, x')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            if 'T' in kwargs.keys():
                # p-T
                my = self.fluid_class.my_pt(p, T)
            elif 'h' in kwargs.keys():
                # p-h
                my = self.fluid_class.my_ph(p, h)
            elif 's' in kwargs.keys():
                # p-s
                my = self.fluid_class.my_ps(p, s)
            elif 'x' in kwargs.keys():
                # p-x
                if x == 0.0:
                    my = self.fluid_class.myL_p(p)
                elif x == 1.0:
                    my = self.fluid_class.myV_p(p)
                else:
                    my = np.nan
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        elif 'h' in kwargs.keys() and 's' in kwargs.keys():
            # h-s
            my = self.fluid_class.my_hs(h, s)
        elif 'T' in kwargs.keys() and 'x' in kwargs.keys():
            # T-x
            if x == 0.0:
                my = self.fluid_class.myL_T(T)
            elif x == 1.0:
                my = self.fluid_class.myV_T(T)
            else:
                my = np.nan
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        
        # Convert from SI to specified units system
        return self.unit_converter.convert(my,prop = 'my', to_from = 'from_SI')
      
    def get_tc(self, **kwargs):
        """
        calls fluid class for thermal conductivity and returns results in choosen units.

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        T :     double,     Temperature (valid unit in selected unit system)
        h:      double      Entropy     (valid unit in selected unit system)
        s:      double      Entropy     (valid unit in selected unit system)
        x:      double      Mass-fraction steam (-)

        Returns
        -------
        None.

        """
        # init
        p = np.nan
        T = np.nan
        h = np.nan
        s = np.nan
        x = np.nan
        
        tc = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'T':
                T = self.unit_converter.convert(arg, prop = 'T', to_from = 'to_SI')
            elif key == 'h':
                h = self.unit_converter.convert(arg, prop = 'h', to_from = 'to_SI')
            elif key == 's':
                s = self.unit_converter.convert(arg, prop = 's', to_from = 'to_SI')
            elif key == 'x':
                x = arg
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: p, T, h, s, x')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            if 'T' in kwargs.keys():
                # p-T
                tc = self.fluid_class.tc_pt(p, T)
            elif 'h' in kwargs.keys():
                # p-h
                tc = self.fluid_class.tc_ph(p, h)
            elif 's' in kwargs.keys():
                # p-s
                tc = self.fluid_class.tc_ps(p, s)
            elif 'x' in kwargs.keys():
                # p-x
                if x == 0.0:
                    tc = self.fluid_class.tcL_p(p)
                elif x == 1.0:
                    tc = self.fluid_class.tcV_p(p)
                else:
                    tc = np.nan
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        elif 'h' in kwargs.keys() and 's' in kwargs.keys():
            # h-s
            tc = self.fluid_class.tc_hs(h, s)
        elif 'T' in kwargs.keys() and 'x' in kwargs.keys():
            # T-x
            if x == 0.0:
                tc = self.fluid_class.tcL_T(T)
            elif x == 1.0:
                tc = self.fluid_class.tcV_T(T)
            else:
                tc = np.nan
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-T, p-h, p-s, h-s, p-x, T-x')
        
        # Convert from SI to specified units system
        return self.unit_converter.convert(tc,prop = 'tc', to_from = 'from_SI')
       
    def get_st(self, **kwargs):
        """
        calls fluid class for surface-tension and returns results in choosen units.

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        T :     double,     Temperature (valid unit in selected unit system)

        Returns
        -------
        None.

        """
        # init
        p = np.nan
        T = np.nan
        
        st = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'T':
                T = self.unit_converter.convert(arg, prop = 'T', to_from = 'to_SI')
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: p, T')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            # _p
            st = self.fluid_class.st_p(p)
        elif 'T' in kwargs.keys():
            # _T
            st = self.fluid_class.st_t(T)
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p, T')
        
        # Convert from SI to specified units system
        return self.unit_converter.convert(st,prop = 'st', to_from = 'from_SI')
    
    
    def get_x(self, **kwargs):
        """
        calls fluid class for vapor mass fraction and returns results in choosen units.

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        h:      double      Entropy     (valid unit in selected unit system)
        s:      double      Entropy     (valid unit in selected unit system)

        Returns
        -------
        None.

        """
        # init
        p = np.nan
        h = np.nan
        s = np.nan
        
        x = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'h':
                h = self.unit_converter.convert(arg, prop = 'h', to_from = 'to_SI')
            elif key == 's':
                s = self.unit_converter.convert(arg, prop = 's', to_from = 'to_SI')
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: p, h, s')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            if 'h' in kwargs.keys():
                # p-h
                x = self.fluid_class.x_ph(p, h)
            elif 's' in kwargs.keys():
                # p-s
                x = self.fluid_class.x_ps(p, s)
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-h, p-s, h-s')
        elif 'h' in kwargs.keys() and 's' in kwargs.keys():
            # h-s
            x = self.fluid_class.x_hs(h, s)
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-h, p-s, h-s')
        
        # Convert from SI to specified units system
        return x
    
    def get_vx(self, **kwargs):
        """
        calls fluid class for vapor volume fraction and returns results in choosen units.

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        h:      double      Entropy     (valid unit in selected unit system)
        s:      double      Entropy     (valid unit in selected unit system)

        Returns
        -------
        None.

        """
        # init
        p = np.nan
        h = np.nan
        s = np.nan
        
        x = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'h':
                h = self.unit_converter.convert(arg, prop = 'h', to_from = 'to_SI')
            elif key == 's':
                s = self.unit_converter.convert(arg, prop = 's', to_from = 'to_SI')
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: p, h, s')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            if 'h' in kwargs.keys():
                # p-h
                x = self.fluid_class.vx_ph(p, h)
            elif 's' in kwargs.keys():
                # p-s
                x = self.fluid_class.vx_ps(p, s)
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-h, p-s, h-s')
        elif 'h' in kwargs.keys() and 's' in kwargs.keys():
            # h-s
            x = self.fluid_class.vx_hs(h, s)
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-h, p-s, h-s')
        
        # Convert from SI to specified units system
        return x
    
    def get_T(self, **kwargs):
        """
        calls fluid class for temperature and returns results in choosen units.

        Parameters
        ----------
        kwargs:
        p :     double,     Pressure    (valid unit in selected unit system)
        h:      double      Entropy     (valid unit in selected unit system)
        s:      double      Entropy     (valid unit in selected unit system)
        x:      double      Mass-fraction steam (-)

        Returns
        -------
        None.

        """
        # init
        p = np.nan
        h = np.nan
        s = np.nan
        
        T = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'p':
                p = self.unit_converter.convert(arg, prop = 'p', to_from = 'to_SI')
            elif key == 'h':
                h = self.unit_converter.convert(arg, prop = 'h', to_from = 'to_SI')
            elif key == 's':
                s = self.unit_converter.convert(arg, prop = 's', to_from = 'to_SI')
            elif key == 'x':
                x = arg
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: p, h, s, x')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'p' in kwargs.keys():
            if 'h' in kwargs.keys():
                # p-h
                T = self.fluid_class.T_ph(p, h)
            elif 's' in kwargs.keys():
                # p-s
                T = self.fluid_class.T_ps(p, s)
            elif 'x' in kwargs.keys():
                # p-x
                T = self.fluid_class.Tsat_p(p)
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-h, p-s, h-s, p-x')
        elif 'h' in kwargs.keys() and 's' in kwargs.keys():
            # h-s
            T = self.fluid_class.T_hs(h, s)
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: p-h, p-s, h-s, p-x')
        
        # Convert from SI to specified units system
        return self.unit_converter.convert(T,prop = 'T', to_from = 'from_SI')
    
    def get_p(self, **kwargs):
        """
        calls fluid class for pressure and returns results in choosen units.

        Parameters
        ----------
        kwargs:
        T :     double,     Temperature (valid unit in selected unit system)
        h:      double      Entropy     (valid unit in selected unit system)
        s:      double      Entropy     (valid unit in selected unit system)
        x:      double      Mass-fraction steam (-)

        Returns
        -------
        None.

        """
        # init
        T = np.nan
        h = np.nan
        s = np.nan
        x = np.nan
        
        p = np.nan
        
        # Read convert and set kwargs, ignore onces not defined
        # Inputs converted to SI for use by fluid-class
        for key, arg in kwargs.items():
            if key == 'T':
                T = arg
            elif key == 'h':
                h = self.unit_converter.convert(arg, prop = 'h', to_from = 'to_SI')
            elif key == 's':
                s = self.unit_converter.convert(arg, prop = 's', to_from = 'to_SI')
            elif key == 'x':
                x = arg
            else:
                raise Exception(f'{key} is not a valid keyword, Valid keywords: T, h, s, x')
                
        # Identify valid combination of kwargs and call corresponding fluid-class function
        if 'x' in kwargs.keys():
            if 'T' in kwargs.keys():
                # T
                p = self.fluid_class.psat_t(T)
            elif 'h' in kwargs.keys():
                if 's' in kwargs.keys():
                    # h-s-sat
                    p = self.fluid_class.psat_hs(h, s)
                else:
                    raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: T-x, s-x, h-s, h-s-x')
            elif 's' in kwargs.keys():
                # s
                p = self.fluid_class.psat_s(s)
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: T-x, s-x, h-s, h-s-x')
        elif 'h' in kwargs.keys():
            if 's' in kwargs.keys():
                # h-s
                p = self.fluid_class.p_hs(h, s)
            else:
                raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: h-s')
        else:
            raise Exception('Not enough input arguments or input arguments not a valid combination. Valid combinations: h-s')
        
        # Convert from SI to specified units system
        return self.unit_converter.convert(p,prop = 'p', to_from = 'from_SI')
       
    def set_properties(self,properties_dict):
        """
        Convert and set each property after calculation
        All other functions uses SI-units. This function converts incomming
        units to specified unit-systems and stores them in "self"
        
        Parameters
        ----------
        properties_dict : dict of doubles    Dict with all property values in SI-units

        Returns
        -------
        None.

        """
        converted_dict = {'p': np.nan,
                          'T': np.nan,
                          'v': np.nan,
                          'h': np.nan,
                          'u': np.nan,
                          's': np.nan,
                          'cp': np.nan,
                          'cv': np.nan,
                          'w': np.nan,
                          'rho': np.nan,
                          'x': np.nan,
                          'my': np.nan,
                          'tc': np.nan,
                          'st': np.nan}
        
        for key in properties_dict:
            converted_dict[key] = self.unit_converter.convert(properties_dict[key],prop = key, to_from = 'from_SI')
        
        self.p = converted_dict['p']
        self.T = converted_dict['T']
        self.v = converted_dict['v']
        self.h = converted_dict['h']
        self.u = converted_dict['u']
        self.s = converted_dict['s']
        self.cp = converted_dict['cp']
        self.cv = converted_dict['cv']
        self.w = converted_dict['w']
        self.x = converted_dict['x']
        self.rho = converted_dict['rho']
        self.my = converted_dict['my']
        self.tc = converted_dict['tc']
        self.st = converted_dict['st']
    
    
    @staticmethod
    def get_fluids():
        """
        give list with implemented fluids

        Returns
        -------
        fluid_list : list of strings    list of fluids

        """
        fluid_list = ['H2O', 'N2']
        return fluid_list
    
    @staticmethod
    def get_unit_systems():
        """
        give list with implemented unit-systems

        Returns
        -------
        unit_systems_list : list of strings    list of unit systems

        """
        
        return unit_conv.get_unit_systems()
    
    def __str__(self):
        """
        Lists basic properties of object as a string
        """
        
        # Units used:
        units = self.unit_converter.get_units()
        
        return(f'\n\nFluid:\t\t\tWater (H2O)\nUnit system:\t{self.unit_converter.unit_system}\n\n' +
               f'Pressure,\t\tp = {self.p:.6f}\t{units["p"]}\n' +
               f'Temperature,\tT = {self.T:.6f}\t{units["T"]}\n\n' +
               f'region =\t{self.fluid_class.region}\n\n'
               f'Specific Volume,\t\t\t\t\t\tv  =\t\t{self.v:.11f}\t\t{units["v"]}\n' +
               f'Specific Enthalpy,\t\t\t\t\t\th  =\t\t{self.h:.11f}\t{units["h"]}\n' +
               f'Specific Internal Energy,\t\t\t\tu  =\t\t{self.u:.11f}\t{units["u"]}\n' +
               f'Specific Entropy,\t\t\t\t\t\ts  =\t\t{self.s:.11f}\t\t{units["s"]}\n' +
               f'Specific Isobaric heat capacity,\t\tcp =\t\t{self.cp:.11f}\t{units["cp"]}\n' +
               f'Specific Isochoric heat capacity,\t\tcv =\t\t{self.cv:.11f}\t\t{units["cv"]}\n' +
               f'Speed of sound,\t\t\t\t\t\t\tw  =\t\t{self.w:.11f}\t\t{units["w"]}\n' + 
               f'Vapor mass fraction,\t\t\t\t\tx  =\t\t{self.x:.11f}\t\t(-)\n' +
               f'Dynamic Viscosity,\t\t\t\t\t\tmy =\t\t{self.my:.11f}\t\t{units["my"]}\n' +
               f'Thermal Conductivity,\t\t\t\t\ttc =\t\t{self.tc:.11f}\t\t{units["tc"]}\n' +
               f'Surface Tension,\t\t\t\t\t\tst =\t\t{self.st:.11f}\t\t{units["st"]}\n')
    
    
if __name__ == "__main__":
    main()
