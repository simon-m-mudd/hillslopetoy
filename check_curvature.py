#!/usr/bin/env python



import hillslopetoy as ht
import numpy as np

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
        meas_erosion = row[0]/1000000
        meas_curv = row[1]
        print(str(meas_erosion)+"  "+str(meas_curv))
        apparent_D = ht.calculate_apparent_D(meas_erosion,meas_curv,half_length = 7, spacing = 1, S_c = 1.0, rho_ratio=2)
        apparent_D= [element*1000 for element in apparent_D]
        apparent_D.append(row[2])
        apparent_D.append(row[0])
        apparent_D.append(row[1])
        calc_D.append(apparent_D)
        

    np.savetxt("manny_uncert.csv", calc_D, delimiter=",")

    
def run_tests_2():
    "I am going to run some hillslope tests for you."
    
    meas_erosion = 36.5/1000000
    ##meas_curv =  -0.02536
    #apparent_D = ht.calculate_apparent_D_with_noise(meas_erosion,meas_curv,half_length = 7, spacing = 1, S_c = 1.0, rho_ratio=2,topographic_uncert = 0.2, n_iterations = 3)
    #print("apparent diffusion coefficient is: ")
    #print(apparent_D)
    #print(" m^2/kyr")
    
    
    data = np.loadtxt("manny_data.csv",delimiter=",",skiprows=1)
    #print(data)
    
    calc_D = []
    
    for row in data:
        #print("This is row:" + str(row))
        meas_erosion = row[0]/1000000
        meas_curv = row[1]
        print(str(meas_erosion)+"  "+str(meas_curv))
        apparent_D = ht.calculate_apparent_D_with_noise(meas_erosion,meas_curv,half_length = 7, spacing = 1, S_c = 1.0, rho_ratio=2,topographic_uncert = 0.2, n_iterations = 50)
        apparent_D= [element*1000 for element in apparent_D]
        apparent_D.append(row[2])
        apparent_D.append(row[0])
        apparent_D.append(row[1])
        calc_D.append(apparent_D)
        

    np.savetxt("manny_uncert.csv", calc_D, delimiter=",")    
    
if __name__ == "__main__":
    run_tests_2()
