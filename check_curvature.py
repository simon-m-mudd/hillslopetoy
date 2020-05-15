#!/usr/bin/env python



import hillslopetoy as ht

def run_tests():
    "I am going to run some hillslope tests for you."
    
    meas_erosion = 36.5/1000000
    meas_curv =  -0.02536
    apparent_D = ht.calculate_apparent_D(meas_erosion,meas_curv,half_length = 7, spacing = 1, S_c = 1.0, rho_ratio=2)
    print("apparent diffusion coefficient is: "+str(apparent_D*1000)+" m^2/kyr")


if __name__ == "__main__":
    run_tests()
