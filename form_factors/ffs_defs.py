
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2026-15-5 2:25
# @Version : 0.3
# Form factor defintions


import numpy as np
import scipy

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


def construct_a_list(base_params, gA, tc): 
    """constructsa list that solves for a values when multiplied by b, which depends on base parameters of a1-a4
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


    #Initialize the b vector for Aa = b
    b = np.zeros(9)
    b[1:5] = base_params
    b[0] = gA

    return np.linalg.solve(A,b)

def generate_samples(a_list, a_std_dv_list, correlation_matrix, N):
    
    # define covariance matrix from correlation matrix and standard deviations
    D = np.diag(a_std_dv_list)
    covariance_matrix = D @ correlation_matrix @ D

    # Use numpy's multivariate_normal function to draw N samples of a1-a4
    data_points = np.random.multivariate_normal(
        mean=a_list,
        cov=covariance_matrix,
        size=N
    )
    return data_points

def axial_conf_err(Q2, tc, data_points):
    """ calculate axial form factor errors (conformal map model)
    Input:
    Q2 : value of -(pf-pi)^2
    tc : value of branch cut in the Q2 complex plane
    a_list : list of coefficients of the conformal map expansion
    a_std_dv_list : list of std deviations for each value of coefficient of conformal map expansion
    correlation_matrix : correlation matrix for the coefficients
    """
    
    #need to only make samples once...
    #this should have the 1000 samples as an input rather than resampling every time for every Q^2
    
    # Apply the function to all 10,000 sampled data points
    ff_list = [axial_conf(Q2, tc, dat) for dat in data_points]
    
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

def pade_G(Q2, M, a_list, b_list):
    """ calculate G(Q^2) using rational function expansion
    Input:
    Q2     : value of momentum transfer squared
    M      : mass scale parameter
    a_list : list of coefficients for the numerator [a0, a1, ... an]
    b_list : list of coefficients for the denominator [b1, b2, ... bn]
             (note: denominator sum starts at k=1, so b_list[0] is b1)
    """
    tau = Q2 / (4 * M**2)
    # Calculate numerator: sum_{k=0}^{n} a_k * tau^k
    num = 0
    for k, ak in enumerate(a_list):
        num += ak * tau**k

    # Calculate denominator: 1 + sum_{k=1}^{n} b_k * tau^k
    den = 1
    for k, bk in enumerate(b_list):
        # b_list starts at index 0 which corresponds to b1, so power is k+1
        den += bk * tau**(k+1)

    return num / den

def dipole_Gd(Q2, Lambda2=0.71):
    """ Standard dipole form factor
    G_d(Q^2) = (1 + Q^2 / Lambda^2)^{-2}

    Input:
    Q2      : momentum transfer squared (GeV^2)
    Lambda2 : dipole mass parameter squared (default 0.71 GeV^2)
    """
    return (1 + Q2 / Lambda2)**(-2)

def calculate_F1(Q2, mN, GE, GM):
    """ calculate Dirac (F1) form factor
    F_1 = (G_E + tau * G_M) / (1 + tau)
    Where tau = Q^2 / (4 * mN^2)

    Input:
    Q2 : value of -(pf-pi)^2
    mN : Nucleon mass
    GE : Electric form factor G_E(Q^2)
    GM : Magnetic form factor G_M(Q^2)
    """
    tau = Q2 / (4 * mN**2)
    return (GE + (tau * GM)) / (1 + tau)

def calculate_F2(Q2, mN, GE, GM):
    """ calculate Pauli (F2) form factor
    F_2 = (G_M - G_E) / (1 + tau)
    Where tau = Q^2 / (4 * mN^2)

    Input:
    Q2 : value of -(pf-pi)^2
    mN : Nucleon mass
    GE : Electric form factor G_E(Q^2)
    GM : Magnetic form factor G_M(Q^2)
    """
    tau = Q2 / (4 * mN**2)
    return (GM - GE) / (1 + tau)

def generate_G_coef_samples(a_list, b_list, a_err, b_err, N):
    
    # define covariance matrix from correlation matrix and standard deviations 
    a_covariance_matrix = np.diag(a_err)**2
    b_covariance_matrix = np.diag(b_err)**2

    # Use numpy's multivariate_normal function to draw N samples of a1-a4
    a_data_points = np.random.multivariate_normal(
        mean=a_list,
        cov=a_covariance_matrix,
        size=N
    )
    b_data_points = np.random.multivariate_normal(
        mean=b_list,
        cov=b_covariance_matrix,
        size=N
    )

    # need to somehow impose constraints here to not get samples with poles, specifically for en
    #could just use one of the constraints to directly solve for Gen, but this seems suspicious - why are they giving values for all the parameters if some directly depend on the others? 
    #alternatively could use a different parameterization...
    #I wonder what Aaron did here
    
    return a_data_points, b_data_points


def pade_G_uncertainty_samples(Q2, M, a_data_points, b_data_points):
    """ calculate G(Q^2) using rational function expansion
    Input:
    Q2     : value of momentum transfer squared
    M      : mass scale parameter
    a_list : list of coefficients for the numerator [a0, a1, ... an]
    b_list : list of coefficients for the denominator [b1, b2, ... bn]
             (note: denominator sum starts at k=1, so b_list[0] is b1)
    """
    tau = Q2 / (4 * M**2)

    # Calculate numerator: vector of tau^k 's
    tau_powers =  tau ** (np.arange(a_data_points.shape[-1]))

    #sum up tau^k * ak 's
    num = np.sum(tau_powers * a_data_points, axis = 1)
    
    # Calculate numerator: vector of tau^k 's
    tau_powers =  tau ** (np.arange(1, b_data_points.shape[-1] + 1))

    #sum up tau^k * bk 's
    den = np.sum(tau_powers * b_data_points, axis = 1)
    
    G_list = num / (1 + den)
    
    return G_list


def F1_uncertainty(Q2, mN, GEp_samples, GMp_samples, GEn_samples, GMn_samples):
    """ calculate uncertainty for Dirac (F1) form factor via error propagation

    Input:
    Q2     : value of -(pf-pi)^2
    mN     : Nucleon mass
    GE     : Electric form factor G_E(Q^2)
    GM     : Magnetic form factor G_M(Q^2)
    GE_unc : uncertainty on G_E
    GM_unc : uncertainty on G_M
    """
    mu_p = 2.793
    mu_n = -1.913
    
    tau = Q2 / (4 * mN**2)
    #initialize sachs vector current form factors
    G_ve = GEp_samples - GEn_samples
    G_vm = GMp_samples - GMn_samples
    
    return np.std((G_ve + tau*G_vm)/(1 + tau))

def F2_uncertainty(Q2, mN, GEp_samples, GMp_samples, GEn_samples, GMn_samples):
    """ calculate uncertainty for Pauli (F2) form factor via error propagation

    Input:
    Q2     : value of -(pf-pi)^2
    mN     : Nucleon mass
    GE     : Electric form factor G_E(Q^2)
    GM     : Magnetic form factor G_M(Q^2)
    GE_unc : uncertainty on G_E
    GM_unc : uncertainty on G_M
    """
    mu_p = 2.793
    mu_n = -1.913
    
    tau = Q2 / (4 * mN**2)
    #initialize sachs vector current form factors
    G_ve = GEp_samples - GEn_samples
    G_vm = GMp_samples - GMn_samples

    return np.std((G_vm - G_ve)/(1 + tau))
        
def g2(Q2, axial_ff, Mn, Mpi):
    
    return 2*(Mn**2)*axial_ff/(Mpi**2+Q2)
    
def A_t(Q2, m_l, m_N, f1, f2, g1, g2):
    term1 =(1+Q2/(4*m_N**2))*(g1**2)
    term2 = -(1-Q2/(4*m_N**2))*(f1**2- (Q2/(4* m_N**2))*f2**2)
    term3 = Q2/(m_N**2)*f1*f2
    term4 = -((m_l**2)/(4*m_N**2))*((f1+f2)**2+(g1+2*g2)**2 - 4*((1+Q2/(4*m_N**2)))*(g2**2))

    A = ((m_l**2 + Q2)/(m_N**2))*(term1+term2+term3+term4)
    return A

def B_t(Q2, m_l, m_N, f1, f2, g1, g2):
    return (Q2/m_N**2)*g1*(f1+f2) 
    
def C_t(Q2, m_l, m_N, f1, f2, g1, g2):
    return 0.25*(g1**2 + f1**2 + (Q2/(4*(m_N**2))*(f2**2)))

def diff_cs(E_v, Q2, m_l, m_N, A, B, C, G_F=1.16637e-5, cosc=0.9746, hbc=0.1973e-13):
    """ 
    Calculate the differential cross section d sigma / dQ^2
    """
    s_u = 4 * E_v * m_N - Q2 - m_l**2
    norm = (m_N * G_F * cosc)**2 / (8 * np.pi * (E_v**2)) * (hbc)**2
    return norm * (A - B * (s_u / m_N**2) + C * (s_u / m_N**2)**2)


def total_cs(E_v, Q2, m_l, m_N, A, B, C, G_F=1.16637e-5, cosc=0.9746, hbc=0.1973e-13):
    """ 
    Calculate the total cross section sigma by integrating 
    d sigma / dQ^2 over Q2 until the kinematic limit. 
    """
    E_thr = (m_l**2 + 2 * m_N * m_l) / (2 * m_N)  # exact threshold
    if E_v < E_thr:
        return 0
    
    #insert Q2 limit here
    Q2_max = E_v/(2*E_v + m_N)*(2*E_v*m_N + m_l**2 + np.sqrt(4*(E_v**2)*(m_N**2) - 4*m_N*(E_v + m_N)*m_l**2 + m_l**4)) - m_l**2
    Q2_min = E_v/(2*E_v + m_N)*(2*E_v*m_N+m_l**2 - np.sqrt(4*(E_v**2)*(m_N**2) - 4*m_N*(E_v+m_N)*m_l**2 + m_l**4)) - m_l**2

    #mask only allows values at indexes that meet the condition
    mask = np.logical_and(Q2 <= Q2_max, Q2 >= Q2_min)

    Q2_bdd = Q2[mask]
    A_bdd = A[mask]
    B_bdd = B[mask]
    C_bdd = C[mask]
    
    s_u = 4 * E_v * m_N - Q2_bdd - m_l**2
    norm = (m_N * G_F * cosc)**2 / (8 * np.pi * (E_v**2)) * (hbc)**2
  
    return scipy.integrate.trapezoid(norm * (A_bdd - B_bdd * (s_u / m_N**2) + C_bdd * (s_u / m_N**2)**2), Q2_bdd)


