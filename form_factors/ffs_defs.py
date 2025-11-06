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



def axial_dipole(Q2, gA, mA):
    """ calculate axial form factor (conformal map model)
    Input:
    Q2 : value of -(pf-pi)^2
    gA : axial charge (F_A(0))
    mA : dipole mass
    """
    ff = -gA*((1+Q2/(mA**2))**-2)

    return ff
