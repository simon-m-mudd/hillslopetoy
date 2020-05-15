"""Main module."""

from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import optimize

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

    # This comes from the Roering et al. 2007 solution
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

    # This comes from the Roering et al. 2007 solution
    curv_analytical = -rho_ratio*C_0/D
    return curv_analytical      
    

def set_profile_locations_half_length(half_length = 8,spacing = 1):
    """This sets the profile for a distance either side of the hillcrest

    Args:
        half_length (float): The distance covered by the profile away rom the hillcrest
        spacing (float): The spacing of x  in metres 

    Returns:
        A numpy array with the hillslope locations in metres

    Author:
        Simon M Mudd
        
    Date:
        14/05/2020    
    """
    minimum_x = -half_length
    maximum_x = half_length+0.01
    x_locs = np.arange(minimum_x,maximum_x,spacing)
    return x_locs   
    

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
    

def displace_profile_locations_constant(x, displacement_distance = 0.1):
    """This shifts the profile locations by a fixed amount

    Args:
        x (array): The array of x locations
        displacement_distance (float): How far the profile locations will be displaced in metres. 
        Simulates missing the hillcrest in the DEM

    Returns:
        A numpy array with the hillslope locations in metres

    Author:
        Simon M Mudd
        
    Date:
        14/05/2020    
    """
    x_locs = np.add(x,displacement_distance)
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
    p = np.polyfit(x, z, poly_degree)
    curv_meas = 2*p[0]
    return curv_meas


def calculate_apparent_D(meas_erosion,meas_curv,half_length = 7, spacing = 1, S_c = 1.0, rho_ratio=2):
    """This functions aim is to see whqat the actualk curvature would be if the profile were displaced

    Args:
        meas_erosion (float): the measured erosion rate (in m/yr)
        meas_curv (float): the curvature measured from the DEM
        half_length (float): The distance covered by the profile away rom the hillcrest
        spacing (float): The spacing of x  in metres       
        S_c (float): Critical slope (dimensionless)
        rho_ratio (float): ratio between rock and soil density

    Returns:
        The ridgetop curvature in 1/m

    Author:
        Simon M Mudd
        
    Date:
        15/05/2020    
    """  
    
    # first guess the diffusivity from the measured data
    D_apparent = -rho_ratio*meas_erosion/meas_curv
    print("The apparent D is:"+str(D_apparent)+" or "+str(D_apparent*1000)+" in m^2/kyr")
    
    # now we get the baseline sampling
    x_loc_baseline = set_profile_locations_half_length(half_length = half_length,spacing = spacing)
    
    # now solve the D needed to get the apparent D
    print("Right, going into the optimization loop")
    #print("By the way, the starting C is")
    #first_curv = curvature_difference_function(D_apparent, x_loc_baseline, S_c, meas_erosion, rho_ratio,meas_curv)
    #print("meas: "+str(meas_curv)+" and difference from the test curvature: "+str(first_curv))
    #big = curvature_difference_function(D_apparent*2, x_loc_baseline, S_c, meas_erosion, rho_ratio,meas_curv)
    #small = curvature_difference_function(D_apparent*0.5, x_loc_baseline, S_c, meas_erosion, rho_ratio,meas_curv)
    #print("big:" +str(big)+" and small: "+str(small))
    
    # plot the curvature function
    #D_vals = [D_apparent*0.5,D_apparent,D_apparent*2]
    #C_offset_vals = []
    #for D in D_vals:
    #    C_offset_vals.append( curvature_difference_function(D, x_loc_baseline, S_c, meas_erosion, rho_ratio,meas_curv) )
    #print("C offsets are: ")
    #print(C_offset_vals)
        
    
    # Run the optimisation
    root = optimize.newton(curvature_difference_function,D_apparent,args=(x_loc_baseline,S_c,meas_erosion,rho_ratio,meas_curv))
    #print("The D needed to get the measured curvature is: "+str(root))
    
    # now test
    #z = ss_nonlinear_elevation(x_loc_baseline,S_c = S_c,C_0 = meas_erosion , D=root ,rho_ratio = rho_ratio)
    #app_curv = fit_hilltop_curvature(x_loc_baseline,z)
    
    #print("Measured curvature is: "+str(meas_curv)+ " and the apparent curvature measured from a gridded sample is: "+str(app_curv))
    
    return root
    
    

    
def curvature_difference_function(D, x, S_c, C_0, rho_ratio,meas_curvature):
    """This gets the difference between the measured curvature and the curvature calculated for a given value of D.
    Used with Newton-Raphson optimisation.

    Args:
        D (float): sediment transport coefficient in m^2/yr   
        x (array): The array of x locations
        S_c (float): Critical slope (dimensionless)
        C_0 (float): Steady uplift/erosion in m/yr       
        rho_ratio (float): ratio between rock and soil density
        meas_curvature (float): the measured curvature. When the apparent cur

    Returns:
        The curvature as measured from the profile x 

    Author:
        Simon M Mudd
        
    Date:
        15/05/2020    
    """ 
    
    # first create the x data
    z = ss_nonlinear_elevation(x,S_c = S_c,C_0 = C_0 , D = D ,rho_ratio = rho_ratio)
    fit_curvature = fit_hilltop_curvature(x,z)
    print("fit_curvature is: "+str(fit_curvature))
    curvature_offset = meas_curvature-fit_curvature
    
    return curvature_offset
    

            
    
    
def ratio_of_fit_to_analytical_hilltop_curvature(displacement_distances, half_length = 7, spacing = 1, S_c = 1.0, C_0 = 0.001 ,D=0.01, rho_ratio=2):
    """This takes the curvature from the profiles and fits a polynomial to minic
    what happens inside LSDTopoTools curvature algorithm

    Args:
        displacement_distance (float array or list): Displacement distances: how far away from the divide the gidded divide point is
        half_length (float): The distance covered by the profile away rom the hillcrest
        spacing (float): The spacing of x  in metres       
        S_c (float): Critical slope (dimensionless)
        C_0 (float): Steady uplift/erosion in m/yr
        D (float): sediment transport coefficient in m^2/yr
        rho_ratio (float): ratio between rock and soil density

    Returns:
        A list of the ratio of measured curvature to analytical curvature

    Author:
        Simon M Mudd
        
    Date:
        14/05/2020    
    """ 
    
    # first create the x data
    x_loc_baseline = set_profile_locations_half_length(half_length = half_length,spacing = spacing)
    analytical_curvature = ss_nonlinear_hilltop_curvature(C_0 = C_0, D = D, rho_ratio = rho_ratio)
    
    curv_ratios = []
    
    # now enter a loop where we do displacements
    for displacement in displacement_distances:
        displaced_x_loc = displace_profile_locations_constant(x_loc_baseline, displacement_distance = displacement)
        z = ss_nonlinear_elevation(displaced_x_loc,S_c = S_c,C_0 = C_0 ,D = D ,rho_ratio =rho_ratio)
        
        # get the curvature fit
        meas_curv = fit_hilltop_curvature(displaced_x_loc,z)
        curv_ratios.append(meas_curv/analytical_curvature)
        
    return curv_ratios
        
    
    
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
    