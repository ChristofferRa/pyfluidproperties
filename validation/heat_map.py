# -*- coding: utf-8 -*-
"""
This module contains functions for calculating accuracy data given a csv-file
with pre calculated fluid states to compare with. Designed to use NIST-data
https://webbook.nist.gov/cgi/fluid.cgi?ID=C7732185&TUnit=K&PUnit=MPa&DUnit=kg%2Fm3&HUnit=kJ%2Fkg&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm&Type=IsoTherm&RefState=DEF&Action=Page

Copy the table from the link above to a csv-file to get the fluid states

A function for turning the accuracy data into "heat maps" is provided.

The main-function reads a specified .csv, generates accuracy data and makes a 
heat map for each property.

@author: Christoffer Rappmann, christoffer.rappman@gmail.com
"""
import numpy as np
import matplotlib.pyplot as plt
import time

from pyfluidproperties import properties, fluid_id, unit_sys
from pyfluidproperties.fluid.water.iapwsif97 import p_t_b23 as pb23
from pyfluidproperties.fluid.water.iapwsif97 import t_p_b23 as tb23

def main():
    ### Input/Output files
    input_data = 'nist_data_without_sat.csv'                        # Fluid states as csv-file, see "nist_header" below for columns
    output_data_file = 'fluid_states_diff.npy'          # File to save accuracy data to for faster loading
    figure_file_name = 'nist_comparison_to_IAPWSIF95.png'  # Name for saved figures
    
    ### Settings
    recalculate = False # Recalculate accuracy data (True) or load .npy-file (False)
    
    ### Read (NIST) data from csv-file
    nist_header = np.array(['Temperature (K)','Pressure (MPa)','Density (kg/m3)','Specific Volume (m3/kg)','Internal Energy (kJ/kg)','Enthalpy (kJ/kg)','Entropy (kJ/g*K)','Cv (kJ/g*K)','Cp (kJ/g*K)','Sound Spd. (m/s)','Joule-Thomson (K/MPa)','Viscosity (Pa*s)','Therm. Cond. (W/m*K)','Phase'])
    nist_data = np.genfromtxt(input_data, delimiter=',',skip_header=True)
    
    ### Load or generate accuracy data
    if recalculate:
        results_diff_pt = generate_diff_file(nist_data, output_data_file)
    else:
        results_diff_pt = np.load(output_data_file)
    
    ### Generate and save heat-map
    plot_ranges_max = [0, 0, .15, .15, .1, .1, .1, 3.0, 3.0, 2.0, 0, .5, .5]
    #plot_ranges_max = [0, 0, 1.0, 1.0, 0.5, 0.5, 0.5, 2.0, 4.0, 2.0, 0, 2.0, 5.0]
    
    for property_pos in range(2,13): # Loop throug all properties
        figure_save_name = f'{figure_file_name[:-4]}_{property_pos}{figure_file_name[-4:]}'
        heat_map(f'Accuracy of {nist_header[property_pos]} using p-T functions compared to NIST (IAPWSIF95)',  results_diff_pt[:,property_pos], nist_data[:,1], nist_data[:,0], nist_header[1], nist_header[0],plot_ranges_max[property_pos],  figure_save_name)

def heat_map(desc, accuracy_data,x_data, y_data, x_label, y_label, plot_range_max, figure_file_name):
    """
    Plots heat-mapp for supplied accuracy_data

    Parameters
    ----------
    desc :          str                 Plot title
    accuracy_data : np-array of float   2-D array with accuracy data
    x_data :        np-array of float   data for x-axis
    y_data :        np-array of float   data for y-axis
    x_label :       str                 x-axis title
    y_label :       str                 y_axis title

    Returns
    -------
    None.

    """
    
    print('Generating heat-mapp...')
    
    ### Sort into 2d array using binning
    
    # Choose number of bins
    number_of_bins_x = 400 # Set number of bins in x-direction
    number_of_bins_y = 500 # Set number of bins in y-direction
    
    # Caluclate bin edges
    x_max = max(x_data)
    x_min = min(x_data)
    y_max = max(y_data)
    y_min = min(y_data)
    
    x_bins_edges = np.linspace(x_min, x_max, number_of_bins_x+1)
    y_bins_edges = np.linspace(y_min, y_max, number_of_bins_y+1)
    
    data_array_2d = np.zeros((number_of_bins_x,number_of_bins_y))-1
    
    # Sort data into 2d-bins
    for x, y, z in zip(x_data, y_data, accuracy_data) :
        #x_bins
        x_pos = 0
        y_pos = 0
        for bin_pos in range(0, number_of_bins_x):
            if x >= x_bins_edges[bin_pos]:
                x_pos = bin_pos
        for bin_pos in range(0, number_of_bins_y):      
            if y >= y_bins_edges[bin_pos]:
                y_pos = bin_pos
                
        # Write value in bin if larger than existing value
        if data_array_2d[x_pos, y_pos] < z*100: 
            data_array_2d[x_pos, y_pos] = z*100
        #print(f'x = {x} x_pos = {x_pos}, y = {y} y_pos = {y_pos}')
        
    ### Region boundaries
    
    fluid = properties(fluid = fluid_id.Water, unit_sys = unit_sys.SI_MPa_kJ)

    # Saturation line
    t_sat = np.linspace(300,647.08, 300)
    p_sat = [fluid.get_p(T=t, x = 0)*1e-6 for t in t_sat]

    # Region boundary 1 - 3 
    r13 = np.linspace(623.15, 623.15,100)
    p13 = np.linspace(fluid.get_p(T=623.15,x = 0)*1e-6,100,100)
    # Region boundary 2 - 3 
    t23 = np.linspace(623.15, tb23(100*1e6), 100)
    p23 = [pb23(T)*1e-6 for T in t23]

    ### Plot Settings
    aspect_ratio = (y_max - y_min)/((x_max-x_min))
    axis_dimensions = [y_min, y_max, x_min, x_max]
    axis_near_critical = [620, 700, 22, 30]
    aspect_ratio_near_crit = (axis_near_critical[1] - axis_near_critical[0])/((axis_near_critical[3]-axis_near_critical[2]))
    max_accuracy = plot_range_max
    min_accuracy = 0.0
    
    ### Plot figure
    fig, ax = plt.subplots(figsize=(12, 9))
    im = ax.imshow(data_array_2d, aspect = aspect_ratio, vmax = max_accuracy, vmin = min_accuracy, origin = 'lower', extent = axis_dimensions)
    
    # Region Boundaries
    ax.plot(r13,p13, 'k-', linewidth=0.4,  label = '623.15')
    ax.plot(t_sat, p_sat, 'k-', linewidth=0.4, label = 'saturation-line')
    ax.plot(t23,p23,'k-', linewidth=0.4, label = 'region 2-3')
    ax.axis([300, 1070, 0, 100])
    
    # Near critical
    fig1, ax1 = plt.subplots(figsize=(12,9))
    x_bin_size = ((x_max-x_min)/number_of_bins_x)
    y_bin_size = ((y_max-y_min)/number_of_bins_y)

    im1 = ax1.imshow(data_array_2d[int((axis_near_critical[2]-x_min)/x_bin_size):int((axis_near_critical[3]-x_min)/x_bin_size), int((axis_near_critical[0]-y_min)/y_bin_size):int((axis_near_critical[1]-y_min)/y_bin_size)], aspect = aspect_ratio_near_crit, origin = 'lower', extent = axis_near_critical)
        
    fig.colorbar(im, ax=ax, label='Accuracy (%)')
    ax.set_title(desc)
    ax.set_ylabel(x_label)
    ax.set_xlabel(y_label)
    
    fig1.colorbar(im1, ax=ax1, label='Accuracy (%)')
    ax1.set_title(desc)
    ax1.set_ylabel(x_label)
    ax1.set_xlabel(y_label)
    
    ## Save plot
    fig.savefig(figure_file_name, bbox_inches='tight', dpi = 300)
    fig1.savefig(f'{figure_file_name[:-4]}_near_crit{figure_file_name[-4:]}', bbox_inches='tight', dpi = 300)

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
    desc :  string,                 description
    a :     np-array of float,      accuracy data
    data :  np-array of float,      fluid states from which array a is based
    header: np-array of str,        property-names
    
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


def generate_diff_file(fluid_states: str, output_data_file: str = ''):
    """
    Compares calculated data using pyfluidproperties with supplied csv-file (NIST-data)
    Returns data-set with difference/accuracy for each fluid state in input_data
    
    Optional: save to .npy for fast loading

    Parameters
    ----------
    fluid_states : np-array  numpy array containing fluid states to be compared.
    output_data_file : str   output npy-file (optional).

    Returns
    -------
    numpy array of floats   array with difference for each property at supplied fluid states.

    """
    from tqdm import tqdm
    
    ### Initiate fluid object
    w = properties(fluid = "H2O", unit_system = 'SI_MPa_kJ', high_accuracy = False)

    print('\nCalculating properties using p-T....')
    
    ### Calculate properties for each fluid state
    
    # Properties
    # 0  1  2    3  4  5  6  7   8   9  10  11  12  13
    # T, p, rho, v, u, h, s, Cv, Cp, w, JT, my, tc, phase
    results_pt = np.zeros((0,14))
    
    tic = time.perf_counter()
    for i in tqdm(range(0, len(fluid_states))):
        
        # update_pt
        p = fluid_states[i,1] # MPa
        T = fluid_states[i,0]-273.15 # C
        
        w.update(p = p, T = T)
                
        results_pt = np.append(results_pt, [[w.T+273.15, w.p, w.rho, w.v, w.u, w.h, w.s, w.cv, w.cp, w.w, np.nan, w.my, w.tc, np.nan]], axis = 0)
    
    toc = time.perf_counter()
    print(f'\n\Calculatiosn completed!\n{len(fluid_states)} states calculated in {(toc-tic):.1f} s\n\n')
    
    ### Calculate difference/accuracy betweed pyfluidporperties and supplied data
    results_diff_pt = diff(fluid_states,results_pt)
    
    ### Save .npy-file
    if output_data_file:
        print(f'Saving calculated accuracy-data to {output_data_file}')
        np.save(output_data_file, results_diff_pt)
    
    print('All done')
    return results_diff_pt
    
if __name__ == "__main__":
    main()