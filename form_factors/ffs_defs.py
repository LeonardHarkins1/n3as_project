#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2025-10-30 21:46:30
# @Version : 0.1
# Form factor defintions


import numpy as np


def z_conformal(Q2, tc):
    """ calculate conformal variable z
    Input:
    Q2 : value of -(pf-pi)^2
    tc : value of branch cut in the Q2 complex plane
    """
    return np.sqrt(Q2 - tc)


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

    ff = gA*(1+Q2/mA)

    return ff
