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



def construct_a_list(base_params, tc, gA):
    """takes in the list of a1-a4 and constructs a0, a5 - a8 using the sum rules and F_A(0) = gA constraint
    base_params: known a1-a4
    A: Sum rules matrix
    """
    #initialize matrix for system of equations
    A = np.zeros((9,9))
    #identity matrix for kown a1-a4
    A[0:4,0:4] = np.eye(4)
    
    #sum rules, k represents column number (0-8)
    k_values = np.arange(9)
    A[4,:] = 1
    A[5, :] = k_values
    A[6,2:]= k_values[2:]*(k_values[2:]-1)
    A[7,3:]= k_values[3:]*(k_values[3:]-1)*(k_values[3:]-2)

    #F_A(0) = gA constraint
    A[8, :] = z_conformal(0, tc)**k_values
    
    #Initialize the b vector
    b = np.zeros(9)
    b[:4] = base_params
    b[4:8] = 0
    b[8] = gA

    full_a_list = np.linalg.solve(A,b)
    return full_a_list


def axial_conf_err(Q2, tc, a_list, a_std_dv_list, correlation_matrix, gA):
    """ calculate axial form factor errors (conformal map model)
    Input:
    Q2 : value of -(pf-pi)^2
    tc : value of branch cut in the Q2 complex plane
    a_list : list of coefficients of the conformal map expansion
    a_std_dv_list : list of std deviations for each value of coefficient of conformal map expansion
    correlation_matrix : correlation matrix for the coefficients
    """
    N=10**3
    # define covariance matrix from correlation matrix and standard deviations
    D = np.diag(a_std_dv_list)
    covariance_matrix = D @ correlation_matrix @ D

    # Use numpy's multivariate_normal function to draw N samples
    data_points = np.random.multivariate_normal(
        mean=a_list,
        cov=covariance_matrix,
        size=N
    )

    
    # Apply the function to all 10,000 sampled data points
    ff_list = [axial_conf(Q2, tc, construct_a_list(dat, tc, gA)) for dat in data_points]
    
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
