# Verification and Validation of pyfluidproperties
This folder contains modules for verification and validation of fluid properties

## Water (IAPWSIF97)
### Verification tables
Comparison to tables provided in IAPWS-releases is done with "verification_iapws.py"
Releases included:
- IAPWS R7-97 (2012), IF97 main relase
- IAPWS SR2-01 (2014), Backward equations for p_hs region 1 and 2
- IAPWS SR3-03 (2014), Backward equations for T_ph, v_ph, T_ps, v_ps in region 3
- IAPWS SR4-04 (2014), Backward equations for p_hs for region 3 and T_sat_hs in region 4
- IAPWS SR5-05 (2016), Backward equations for v_pt in region 3
- IAPWS R12-08 (2008), Viscosity
- IAPWS R1-76 (2014), Surface Tension
- IAPWS R15-11 (2011), Thermal conductivity

### Validation against NIST
About 2200 fluid states spread out over the applicability range of IAPWSIF97, from 0.1 to 100 MPa and 300 to 2000 K, is compared to NIST (IAPWSIF95), the near critical region is excluded in this comparison. The maximum error for each property is calculated and presented below. 
| Input properties/combination  | $$T$$    | $$p$$    | $$\rho$$  | $$\nu$$	| $$u$$	   | $$h$$    | $$s$$    | $$C_v$$   | $$C_p$$   | $$w$$    | $$\mu$$   | $$\lambda$$   |
| ------                        |------|------|------|------|------|------|------|------|------|------|------|------|
| Pressure-Temperature (p-T)    |0.00 %|0.00 %|0.07 %|0.07 %|0.03 %|0.03 %|0.02 %|0.66 %|0.80 %|0.52 %|0.11 %|0.14 %|	
| Pressure-Enthalpy (p-h)       |0.02 %|0.00 %|0.07 %|0.07 %|0.08 %|0.08 %|0.08 %|0.67 %|0.76 %|0.55 %|0.12 %|0.13 %|
| Pressure-Entropy (p-s)        |0.02 %|0.00 %|0.07 %|0.07 %|0.04 %|0.03 %|0.04 %|0.67 %|0.76 %|0.57 %|0.14 %|0.13 %|
| Enthalpy-Entropy (h-s)        |0.02 %|1.75 %|0.09 %|0.09 %|0.08 %|0.00 %|0.00 %|0.67 %|0.84 %|0.60 %|0.14 %|0.14 %|
| Pressure - Mass fraction Vapor (p-x) x = 0.0|0.00 %|0.00 %|0.07 %|0.07 %|0.02 %|0.02 %|0.02 %|2.39 %|1.27 %|2.10 %|0.06 %|0.10 %|0.05 %	
| Pressure - Mass fraction Vapor (p-x) x = 1.0|0.00 %|0.00 %|0.12 %|0.12 %|0.04 %|0.04 %|0.03 %|0.59 %|0.80 %|0.11 %|0.01 %|0.31 %|
| Temperature - Mass fraction Vapor (T-x) x = 0.0|0.00 %|0.02 %|0.07 %|0.07 %|0.02 %|0.02 %|0.01 %|2.38 %|1.35 %|2.08 %|0.07 %|0.10 %|0.07 %|
| Temperature - Mass fraction Vapor (T-x) x = 1.0|0.00 %|0.02 %|0.15 %|0.15 %|0.04 %|0.05 %|0.04 %|0.58 %|0.82 %|0.11 %|0.01 %|0.32 %|

Accuracy data for all functions in the IAPWIF97 implementation, run "validation/validation_nist.py".
#### Heat maps, comparison to IAPWSIF97
A even more comprehensive comparison to NIST-data (IAPWSIF95) has been performed. Where 250 000 fluid states evenly spread in the range 300 K < T < 1070, 0.1< p < 100 MPa is calculated. A heat map for specific volume ($$\nu$$)  and isobaric heat capacity ($$C_p$$).  

![Heat map for density](heat_maps/nist_comparison_to_IAPWSIF95_2.png)
![Heat map for density Zoomed](heat_maps/nist_comparison_to_IAPWSIF95_2_near_crit.png)
![Heat map for specific volume](heat_maps/nist_comparison_to_IAPWSIF95_3.png)
![Heat map for specific volume Zoomed](heat_maps/nist_comparison_to_IAPWSIF95_3_near_crit.png)
![Heat map for internal energy](heat_maps/nist_comparison_to_IAPWSIF95_4.png)
![Heat map for internal energy Zoomed](heat_maps/nist_comparison_to_IAPWSIF95_4_near_crit.png)
![Heat map for enthalpy](heat_maps/nist_comparison_to_IAPWSIF95_5.png)
![Heat map for enthalpy Zoomed](heat_maps/nist_comparison_to_IAPWSIF95_5_near_crit.png)
![Heat map for entropy](heat_maps/nist_comparison_to_IAPWSIF95_6.png)
![Heat map for entropy Zoomed](heat_maps/nist_comparison_to_IAPWSIF95_6_near_crit.png)
![Heat map for isochoric heat capacity](heat_maps/nist_comparison_to_IAPWSIF95_7.png)
![Heat map for isochoric heat capacity Zoomed](heat_maps/nist_comparison_to_IAPWSIF95_7_near_crit.png)
![Heat map for isobaric heat capacity](heat_maps/nist_comparison_to_IAPWSIF95_8.png)
![Heat map for isobaric heat capacity Zoomed](heat_maps/nist_comparison_to_IAPWSIF95_8_near_crit.png)
![Heat map for speed of sound](heat_maps/nist_comparison_to_IAPWSIF95_9.png)
![Heat map for speed of sound Zoomed](heat_maps/nist_comparison_to_IAPWSIF95_9_near_crit.png)
![Heat map for viscosity](heat_maps/nist_comparison_to_IAPWSIF95_11.png)
![Heat map for viscosity Zoomed](heat_maps/nist_comparison_to_IAPWSIF95_11_near_crit.png)
![Heat map for isobaric heat capacity](heat_maps/nist_comparison_to_IAPWSIF95_12.png)
![Heat map for isobaric heat capacity Zoomed](heat_maps/nist_comparison_to_IAPWSIF95_12_near_crit.png)
