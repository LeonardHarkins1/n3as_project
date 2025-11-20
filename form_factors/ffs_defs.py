#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2025-11-6 8:59:30
# @Version : 0.2
# Form factor defintions


import numpy as np


def z_conformal(Q2, tc):
    """ calculate conformal variable z
    Input:
    Q2 : value of -(pf-pi)^2
    tc : value of branch cut in the Q2 complex plane
    
    Other vairable: 
    t0 : optimal conformal map variable, 
    minimizes |z| within experimental range according to the equation
    t0 = tc * ( 1 - np.sqrt( 1 + Q2max/tc) )
    given in paper to be -0.28 GeV^2 for Q2max = 1 GeV
    """
    t0 = -0.28

    return (np.sqrt(tc + Q2) - np.sqrt(tc - t0))/(np.sqrt(tc + Q2) + np.sqrt(tc - t0))


def axial_conf(Q2, tc, a_list):
    """ calculate axial form factor (conformal map model)
    Input:
    Q2 : value of -(pf-pi)^2
    tc : value of branch cut in the Q2 complex plane
    a_list : list of coefficients of the conformal map expansion
    """
    zz = z_conformal(Q2, tc)

    ff = 0
    
    for nn, aval in enumerate(a_list):
        ff += aval * zz**nn

    return ff



def construct_A_inverse(tc):
    """constructs inverse matrix that solves for a values when multiplied by b, which depends on base parameters of a1-a4
    A: Sum rules, identity matrix for a1-a4, F(0) = gA constraint
    """
    #initialize matrix for system of equations
    A = np.zeros((9,9))
    
    #identity matrix for kown a1-a4
    A[1:5,1:5] = np.eye(4)
    
    #sum rules, k represents column number (0-8)
    k_values = np.arange(9)
    A[5,:] = 1
    A[6, :] = k_values
    A[7,2:]= k_values[2:]*(k_values[2:]-1)
    A[8,3:]= k_values[3:]*(k_values[3:]-1)*(k_values[3:]-2)

    #F_A(0) = gA constraint
    A[0, :] = z_conformal(0, tc)**k_values

    # Calculate the inverse of A
    A_inv = np.linalg.inv(A)
    return A_inv

def construct_a_list(A_inverse, base_params, gA): 
    """ takes A inverse, uses b to find the full a list from the a1-a4 base parameters
    """
    #Initialize the b vector for Aa = b
    b = np.zeros(9)
    b[1:5] = base_params
    b[5:] = 0
    b[0] = gA

    return A_inverse @ b

def axial_conf_err(Q2, tc, a_list, a_std_dv_list, correlation_matrix, gA):
    """ calculate axial form factor errors (conformal map model)
    Input:
    Q2 : value of -(pf-pi)^2
    tc : value of branch cut in the Q2 complex plane
    a_list : list of coefficients of the conformal map expansion
    a_std_dv_list : list of std deviations for each value of coefficient of conformal map expansion
    correlation_matrix : correlation matrix for the coefficients
    """
    N=10**4
    # define covariance matrix from correlation matrix and standard deviations
    D = np.diag(a_std_dv_list)
    covariance_matrix = D @ correlation_matrix @ D

    # Use numpy's multivariate_normal function to draw N samples of a1-a4
    data_points = np.random.multivariate_normal(
        mean=a_list,
        cov=covariance_matrix,
        size=N
    )
    
    A_inverse = construct_A_inverse(tc)
    
    # Apply the function to all 10,000 sampled data points
    ff_list = [axial_conf(Q2, tc, construct_a_list(A_inverse, dat, gA)) for dat in data_points]
    
    return np.std(ff_list)

def axial_dipole(Q2, gA, mA):
    """ calculate axial form factor (dipole model)
    Input:
    Q2 : value of -(pf-pi)^2
    gA : axial charge (F_A(0))
    mA : dipole mass
    """
    ff = gA*((1+Q2/(mA**2))**-2)
    return ff

def axial_dipole_err(Q2, gA, mA,  mA_std_dv):
    """ calculate axial form factor (dipole model)
    Input:
    Q2 : value of -(pf-pi)^2
    gA : axial charge (F_A(0))
    mA : dipole mass
    """
    N=10**4
   
    # Use numpy's normal function to draw N samples of a1-a4
    data_points = np.random.normal(
        loc = mA,
        scale = mA_std_dv,
        size = N
    )
    
    
    # Apply the function to all 10,000 sampled data points
    ff_list = [axial_dipole(Q2, gA, dat) for dat in data_points]
    
    return np.std(ff_list)
