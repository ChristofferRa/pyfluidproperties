# pyfluidproperties
Thermodynamic properties of fluids in python

pyfluidproperties is a library of fluid properties for different fluids.
All fluid-class:es is extension of the fluid_prop class and share the same syntax and basic properties. 
Using the same syntax for all fluids it is easy to change fluid in a calculation.

## Available fluids:
- Water (IAPWSIF-97 extended with SR2-01,SR3-03,SR4-04,SR5-05, R12-08, R15-11 etc.)

At the moment the only available fluid is water, but the aim is to extend the library with more fluids.
Ex. Nitrogen, Air, Hydrogen etc.

### Basic properties:
      Description                         Unit        Letter
    - Pressure                            (Pa)        p
    - Temperature                         (K)         T,t
    - Specific Volume                     (m^3/kg)    v
    - Specific Enthalpy                   (J/kg)      h   
    - Specific Internal Energy            (J/kg)      u
    - Specific Entropy                    (J/kg/K)    s
    - Specific Isobaric heat capacity     (J/kg/K)    cp
    - Specific Isochoric heat capacity    (J/kg/K)    cv
    - Speed of sound                      (m/s)       w 
    - Vapor mass fraction                 (-)         x
    - Dynamic Viscosity                   (Pa*s)      my
    - Thermal Conductivity                (W/m/K)     tc
    - Surface Tension                     (N/m)       st

### Syntax:
Use either as an object containing all the properties or access them individually using the provided static functions. 
The syntax for static functions is inspired by xSteam excel-library (https://xsteam.sourceforge.net/).

#### As and object:
       update_pt(p, t)
       update_px(p, x)
       update_ph(p, h)
       update_ps(p, s)
       
#### Static functions:
       *_pt - properties
       v_pt(p, t)
       h_pt(p, t)
       u_pt(p, t)
       s_pt(p, t)
       cp_ptv
       cv_pt(p, t)
       w_pt(p, t)
       rho_pt(p, t)
       my_pt(p, t)
       tc_pt(p, t)
   
       *_ph - properties
       v_ph(p, h)
       u_ph(p, h)
       s_ph(p, h)
       cp_ph(p, h)
       cv_ph(p, h)
       w_ph(p, h)
       rho_ph(p, h)
       my_ph(p, h)
       tc_ph(p, h)
   
       *_ps - properties
       v_ps(p, s)
       h_ps(p, s)
       u_ps(p, s)
       cp_ps(p, s)
       cv_ps(p, s)
       w_ps(p, s)
       rho_ps(p, s)
       my_ps(p, s)
       tc_ps(p, s)
   
       *_hs - properties
       v_hs(h, s)
       u_hs(h, s)
       cp_hs(h, s)
       cv_hs(h, s)
       w_hs(h, s)
       rho_hs(h, s)
       my_hs(h, s)
       tc_hs(h, s)
   
       Temperature, pressure, enthalpy backward functions
       T_ph(p, h)
       T_ps(p, s)
       T_hs(h, s)
       
       Tsat_p(p)
       Tsat_hs(h, s)
       
       psat_t(t)
       psat_s(s)
       psat_hs(h, s)
       
       p_hs(h, s)
       
       h_tx(t, x)
       h_px(p, x)
       h_prho(p, rho)
       
       s_tx(t, x)
       s_px(p, x)
   
       Vapor Mass fractions
       x_ph(p, h)
       x_ps(p, s)
       x_hs(h, s)
   
       Vapor Volume fractions
       vx_ph(p, h)
       vx_ps(p, s)
       vx_hs(h, s)
   
       Saturation properties, Liquid and Vapor
           x pressure
           hL_p(p),     hV_p(p)
           sL_p(p),     sV_p(p)
           rhoV_p(p),   rhoL_p(p)
           uV_p(p),     uL_p(p)
           CpV_p(p),    CpL_p(p)
           CvV_p(p),    CvL_p(p)
           wV_p(p),     wL_p(p)
           tcL_p(p),    tcV_p(p)
           myL_p(p),    myV_p(p)
       
           x temperature
           hL_T(t),     hV_T(t)
           sL_T(t),     sV_T(t)
           vL_T(t),     vV_T(t)
           uL_T(t),     uV_T(t)
           CpL_T(t),    CpV_T(t)
           CvL_T(t),    CvV_T(t)
           wL_T(t),     wV_T(t)
           rhoV_T(t),   rhoL_Ts(t)
           tcL_T(t),    tcV_T(t)
           myL_T(t),    myV_T(t)
           
        Misc
        st_t(t)
        st_p(p)
   
