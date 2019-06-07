#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 16:10:31 2019

@author: Julio
"""

import numpy as np
from scipy.optimize import least_squares

def lorentzian(ppm, parameters):
    '''
    Inputs:
        ppm  = vector of saturation offsets (ppm)
        
        pars = vector of parameters
        [A1...An, W1...Wn,C1...Cn]
        
    Returns: Single Vector of a Lorentzian
    
    '''
    pars = parameters[0:-1:] # all except the last item
    scaling = parameters[-1]
    
    #pars = np.array(pars)
    pools = int(len(pars)/3)
    freqs = len(ppm)
    
    # 1) Preallocate output
    L = np.zeros( (freqs,pools) )
    for pool in np.arange(0,pools,1):
        P = pars[pool::pools]
        fwhm = P[1]**(2/4)
        L[:,pool] = P[0]*fwhm / (fwhm + (ppm - P[2])**2 )
    
    return np.sum(L,1) + scaling 

def fit_lorentzian(xdata, ydata, x0, LB, UB):
    '''
    Inputs: xdata, ydata, x0, LB, UB
    
    '''
    
    # define function to estimate residuals
    def res(par_guess):
        yhat = lorentzian(xdata, par_guess)
        return yhat - ydata
    
    # fitting
    res_1 = least_squares(res, x0, bounds= (LB,UB), method='trf')
    
    return res_1.x
    