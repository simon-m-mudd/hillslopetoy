"""Main module."""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def ss_nonlinear_elevation(x,S_c = 1.0,C_0 = 0.001 ,D = 0.01 ,rho_ratio =2 ,c =20):
    """This returns the elevations along a transect using the analytical solutions of the hillslope equations. 
    From Roering et al 2007

    Args:
        x (float array): Location of the elevation nodes
        S_c (float): Critical slope (dimensionless)
        C_0 (float): Steady uplift/erosion in m/yr
        D (float): sediment transport coefficient in m^2/yr
        rho_ratio (float): ratio between rock and soil density
        c (float): a scaling constant that sets the elevation of the divide in metres
        
    Returns:
        the elevation along the profile as a numpy array

    Author:
        Simon M Mudd
        
    Date:
        14/05/2020  
    """    

    beta = rho_ratio*C_0
    first_part = - S_c*S_c*0.5/beta
    second_part = np.sqrt( D*D + (2*beta*x/S_c)*(2*beta*x/S_c) )
    third_part = D*np.log(S_c*(second_part+D)*0.5/beta)
    full_z = first_part*(second_part-third_part)+c
    return full_z    
 
def ss_nonlinear_hilltop_curvature(C_0 = 0.001, D = 0.01, rho_ratio =2):
    """This returns the elevations along a transect using the analytical solutions of the hillslope equations. 
    From Roering et al 2007

    Args:
        C_0 (float): Steady uplift/erosion in m/yr
        D (float): sediment transport coefficient in m^2/yr
        rho_ratio (float): ratio between rock and soil density
        
    Returns:
        Curvature of the ridgetop at steady state

    Author:
        Simon M Mudd
        
    Date:
        14/05/2020  
    """    

    curv_analytical = -rho_ratio*C_0/D
    return curv_analytical      
    


def set_profile_locations_constant(minimum_x = -8, maximum_x = 8.01,spacing = 1):
    """This is the most basic function for setting up the profile locations

    Args:
        minimum_x (float): The minimum value of x in metres (divide is at 0)
        maximum_x (float): The maximum value of x in metres  (divide is at 0)
        spacing (float): The spacing of x  in metres 

    Returns:
        A numpy array with the hillslope locations in metres

    Author:
        Simon M Mudd
        
    Date:
        14/05/2020    
    """
    x_locs = np.arange(minimum_x,maximum_x,spacing)
    return x_locs
    

def fit_hilltop_curvature(x,z):
    """This takes the curvature from the profiles and fits a polynomial to minic
    what happens inside LSDTopoTools curvature algorithm

    Args:
        x (array): The array of x locations
        z (array): The array of elevations

    Returns:
        The ridgetop curvature in 1/m

    Author:
        Simon M Mudd
        
    Date:
        14/05/2020    
    """    
    poly_degree = 2
    p = np.polyfit(x_array, z_list, poly_degree)
    curv_meas = 2*p[0]
    return curv_meas
    
def plot_ss_hillslope(x,z,show_figure = False, print_to_file = True, filename = "hillslope_profile.png", ):
    """This prints a hillslope profile
    Args:
        x (array): The array of x locations
        z (array): The array of elevations
        show_figure (bool): If true, show figure
        print_to_file (bool): If true, print to file
        filename (string): Name of file to which the function prints

    Returns:
        Either a shown figure, a printed figure, or both

    Author:
        Simon M Mudd
        
    Date:
        14/05/2020    
    """    
    fig, ax = plt.subplots()
    ax.plot(x, z)

    ax.set(xlabel='distance from hillcrest (m)', ylabel='elevation (m)',
       title='Hillslope elevation')
    ax.grid()

    if print_to_file:
        fig.savefig(filename)

    if show_figure:
        plt.show()
    