#!/usr/bin/env python
# This is a script that checks uncertainty caused by grid resolution
# on estimates of the diffusion coefficient.
# It was written by Simon Mudd at the University of Edinburgh


import hillslopetoy as ht
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd


def power_fxn(E,a,b):
    C = a*(E**b)
    return C

def get_fit_and_label(xdata,ydata):
    '''
    This is a simple function for extracting a power law fit from data
    You feed in the data and then you get back out both the fit as well
    as a string that can be used to label plots.

    Author: SMM

    Date: 30/11/2020
    '''
    popt, pcov = curve_fit(power_fxn,xdata,ydata, p0=(0.6082,0.5))
    residuals = ydata- power_fxn(xdata, popt[0],popt[1])
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((ydata-np.mean(ydata))**2)
    r_squared_field2 = 1 - (ss_res / ss_tot)
    fit_ydata = power_fxn(xdata, popt[0],popt[1])

    round_p0 = np.around(popt[0],2)
    round_p1 = np.around(popt[1],2)
    round_r2 = np.around(r_squared_field2,2)

    sp0 = str(round_p0)
    sp1 = str(round_p1)
    sr2 = str(round_r2)

    power_string = r'$%s*E^{%s}, r^2 = %s$' % (sp0,sp1,sr2)

    return [fit_ydata,power_string]


def plot_results():
    '''
    This plots the results of the uncertainty analysis.
    You need to run run_tests_2 before you plot the results
    since this function reads a file created by
    run_tests_2

    Author: SMM

    Date: 30/11/2020
    '''
    datadf = pd.read_csv("manny_uncert_6mradius_our_data.csv")

    # first we need a fit
    Erate = datadf['erate (mm/kyr)'].to_numpy()
    naive_D = datadf['naive D (m2/kyr)'].to_numpy()

    median_D = datadf['median D (m2/kyr)'].to_numpy()
    D16 = datadf['16th percentile D (m2/kyr)'].to_numpy()
    D84 = datadf['84th percentile D (m2/kyr)'].to_numpy()
    max_D = datadf['maximum D (m2/kyr)'].to_numpy()
    min_D = datadf['minimum D (m2/kyr)'].to_numpy()

    D16err = np.subtract(median_D,D16)
    D84err = np.subtract(D84,median_D)

    min_err = np.subtract(median_D,min_D)
    max_err = np.subtract(max_D,median_D)

    errors = np.vstack((D16err,D84err))
    errors2= np.vstack((min_err,max_err))

    #print(errors)


    # Get the fit to the truncated field data
    datadf2 = datadf[datadf['erate (mm/kyr)'] < 400]

    # first we need a fit
    Erate2 = datadf2['erate (mm/kyr)'].to_numpy()
    naive_D2 = datadf2['naive D (m2/kyr)'].to_numpy()

    median_D2 = datadf2['median D (m2/kyr)'].to_numpy()
    D162 = datadf2['16th percentile D (m2/kyr)'].to_numpy()
    D842 = datadf2['84th percentile D (m2/kyr)'].to_numpy()

    popt, pcov = curve_fit(power_fxn,Erate2,naive_D2, p0=(0.6082,0.5))

    #print(popt)
    #print(pcov)

    residuals = naive_D2- power_fxn(Erate2, popt[0],popt[1])
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((naive_D2-np.mean(naive_D2))**2)
    r_squared_field2 = 1 - (ss_res / ss_tot)
    fit_D2 = power_fxn(Erate2, popt[0],popt[1])
    #print(r_squared_field2)

    fig, ax = plt.subplots()
    fig.set_size_inches(6, 6)
    ax.scatter(Erate, naive_D,label="$D$ directly from $E$ and $C$", color='g',alpha = 0.5)
    ax.errorbar(Erate, median_D, yerr=errors, fmt='ro',label = "$D$ corrected for grid resolution")


    [fit_ydata,power_string1] = get_fit_and_label(Erate2,naive_D2)
    #power_string1 = power_string1+", exculding 2 outliers"
    ax.plot(Erate2,fit_ydata,linestyle=":",alpha = 0.6,label=power_string1)


    [fit_ydata2,power_string2] = get_fit_and_label(Erate2,median_D2)
    power_string2 = power_string2+", corrected"
    ax.plot(Erate2,fit_ydata2,linestyle=":",alpha = 0.6,label=power_string2)


    [fit_ydata3,power_string3] = get_fit_and_label(Erate,naive_D)
    ax.plot(Erate,fit_ydata3,linestyle="-",alpha = 1,label=power_string3)


    [fit_ydata4,power_string4] = get_fit_and_label(Erate,median_D)
    power_string4 = power_string4+", corrected"
    ax.plot(Erate,fit_ydata4,linestyle="-",alpha = 1,label=power_string4)


    ax.set(xlabel='Erosion rate, $E$ (mm/kyr)', ylabel='Transport coefficient, $D$ (m$^2$/kyr)')
    ax.grid()
    ax.legend()

    ax = ht.axis_styler(ax,axis_style="Normal")

    plt.yscale("log")
    plt.xscale("log")

    filename = "calculated_updated.svg"
    fig.savefig(filename,dpi=300)





def run_tests():
    "I am going to run some hillslope tests for you."

    meas_erosion = 36.5/1000000
    meas_curv =  -0.02536
    apparent_D = ht.calculate_apparent_D(meas_erosion,meas_curv,half_length = 7, spacing = 1, S_c = 1.0, rho_ratio=2)
    print("apparent diffusion coefficient is: ")
    print(apparent_D)
    print(" m^2/kyr")


    data = np.loadtxt("manny_data.csv",delimiter=",",skiprows=1)
    print(data)

    calc_D = []


    for row in data:
        original_data = []
        original_data.append(row[0])
        original_data.append(row[1])
        original_data.append(row[2])
        meas_erosion = row[0]/1000000
        meas_curv = row[1]
        print(str(meas_erosion)+"  "+str(meas_curv))
        apparent_D = ht.calculate_apparent_D(meas_erosion,meas_curv,half_length = 7, spacing = 1, S_c = 1.0, rho_ratio=2)

        # convert D in m^2/kyr
        apparent_D= [element*1000 for element in apparent_D]
        apparent_D.insert(0,original_data)
        calc_D.append(apparent_D)


    np.savetxt("manny_uncert.csv", calc_D, delimiter=",",header="erate (mm/kyr), measured curvature (1/m),naive D (m2/kyr),minimum D (m2/kyr),median D (m2/kyr),ymaximum D (m2/kyr), ,")


def run_tests_2():
    "I am going to run some hillslope tests for you."

    #meas_erosion = 36.5/1000000
    ##meas_curv =  -0.02536
    #apparent_D = ht.calculate_apparent_D_with_noise(meas_erosion,meas_curv,half_length = 7, spacing = 1, S_c = 1.0, rho_ratio=2,topographic_uncert = 0.2, n_iterations = 3)
    #print("apparent diffusion coefficient is: ")
    #print(apparent_D)
    #print(" m^2/kyr")


    data = np.loadtxt("manny_data_v5.csv",delimiter=",",skiprows=1)
    print(data)

    calc_D = []

    for row in data:
        original_data = []
        original_data.append(row[0])
        original_data.append(row[1])
        original_data.append(row[2])
        meas_erosion = row[0]/1000000
        meas_curv = row[1]
        print(str(meas_erosion)+"  "+str(meas_curv))
        apparent_D = ht.calculate_apparent_D_with_noise(meas_erosion,meas_curv,half_length = 6, spacing = 1, S_c = 1.0, rho_ratio=2,topographic_uncert = 0.2, n_iterations = 250)
        apparent_D= [element*1000 for element in apparent_D]
        apparent_D[0:0] = original_data
        print(apparent_D)
        calc_D.append(apparent_D)


    np.savetxt("manny_uncert_v5.csv", calc_D, delimiter=",",header="erate (mm/kyr),measured curvature (1/m),naive D (m2/kyr),minimum D (m2/kyr),16th percentile D (m2/kyr),median D (m2/kyr),84th percentile D (m2/kyr),maximum D (m2/kyr),")

def plot_data():
    plot_results()


if __name__ == "__main__":

    '''
    Hello and welcome to the check curvature function.
    You need a csv file with the erosion rate (mm/kyr), curvature (1/m) and diffusion coefficient (m^2/ky)
    You go into run_tests_2 and update the name of this csv
    The output will be a csv with all the uncertainty (you can rename this, currently manny_uncert_2.csv)
    Once you are finished with that, update the below code to run plot_data()
    On line 38 you need to make sure your filename matches that in line 193
    Have fun.'''

    plot_data()
    #run_tests_2()

