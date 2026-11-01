"""Modello a parametri concentrati (3 + 5 GdL, telescopio appeso al mozzo): equazioni di Newton-Eulero non lineari.

Generato da gimbal_mdof.simulatore_concentrato.esporta_codice. Non modificare a mano.
q = [psi_e, psi_m, psi_L, theta_e, theta_h, theta_m, theta_r, theta_p] (rad); x = [q, dq].
M(theta_p) ddq = E tau - h(q, dq);  tau = [tau_psi, tau_theta] (N m).
Telescopio con le componenti nella terna 2 (I_hat, J_hat, K_hat; K_hat = puntamento).
AVVISO: risultato da parametri segnaposto.
"""

import numpy as np
from numpy import sin, cos

PARAMETRI = {
    'J_psi_e': 0.00018494309307983683,
    'J_psi_m': 0.00021224446514033664,
    'J_psi_L': 1.2521778203365463,
    'J_theta_e': 8.941574012585403e-06,
    'J_theta_h': 0.032269662626647655,
    'J_theta_m': 0.0006566945493165637,
    'J_theta_r': 0.0006441475210322294,
    'J_theta_p': 0.0,
    'k_psi_1': 9547.550768128942,
    'k_psi_2': 13090.759416925233,
    'k_theta_1': 1740.5409614753523,
    'k_theta_2': 1533.77362649146,
    'k_theta_3': 9694.189682015414,
    'k_theta_4': 9800.0,
    'c_psi_1': 0.039262275314698954,
    'c_psi_2': 0.05383296855791751,
    'c_theta_1': 0.06364414246590469,
    'c_theta_2': 0.056083544918199904,
    'c_theta_3': 0.35447507577800225,
    'c_theta_4': 0.35834410678688206,
    'b_psi_e': 0.0,
    'b_psi_m': 0.0,
    'b_psi_L': 0.0,
    'b_theta_e': 0.0,
    'b_theta_h': 0.0,
    'b_theta_m': 0.0,
    'b_theta_r': 0.0,
    'b_theta_p': 0.0,
    'I_II': 0.8578208118885499,
    'I_JJ': 0.8662079938719208,
    'I_KK': 0.2562011574375517,
    'I_IJ': 7.126428076259492e-07,
    'I_IK': -0.007507889466143888,
    'I_JK': -1.139346362048e-06,
    'd_I': -0.005197498733951209,
    'd_J': -5.231081386703321e-07,
    'd_K': 0.07312446314306373,
    'm_p': 15.30788838490337,
    'g': 9.80665,
    'J_ax': 0.05331244002435073,
}


def massa(q, p=PARAMETRI):
    psi_e, psi_m, psi_L, theta_e, theta_h, theta_m, theta_r, theta_p = q
    J_psi_e = p['J_psi_e']
    J_psi_m = p['J_psi_m']
    J_psi_L = p['J_psi_L']
    J_theta_e = p['J_theta_e']
    J_theta_h = p['J_theta_h']
    J_theta_m = p['J_theta_m']
    J_theta_r = p['J_theta_r']
    J_theta_p = p['J_theta_p']
    k_psi_1 = p['k_psi_1']
    k_psi_2 = p['k_psi_2']
    k_theta_1 = p['k_theta_1']
    k_theta_2 = p['k_theta_2']
    k_theta_3 = p['k_theta_3']
    k_theta_4 = p['k_theta_4']
    c_psi_1 = p['c_psi_1']
    c_psi_2 = p['c_psi_2']
    c_theta_1 = p['c_theta_1']
    c_theta_2 = p['c_theta_2']
    c_theta_3 = p['c_theta_3']
    c_theta_4 = p['c_theta_4']
    b_psi_e = p['b_psi_e']
    b_psi_m = p['b_psi_m']
    b_psi_L = p['b_psi_L']
    b_theta_e = p['b_theta_e']
    b_theta_h = p['b_theta_h']
    b_theta_m = p['b_theta_m']
    b_theta_r = p['b_theta_r']
    b_theta_p = p['b_theta_p']
    I_II = p['I_II']
    I_JJ = p['I_JJ']
    I_KK = p['I_KK']
    I_IJ = p['I_IJ']
    I_IK = p['I_IK']
    I_JK = p['I_JK']
    d_I = p['d_I']
    d_J = p['d_J']
    d_K = p['d_K']
    m_p = p['m_p']
    g = p['g']
    J_ax = p['J_ax']
    return np.array([
        [J_psi_e, 0, 0, 0, 0, 0, 0, 0],
        [0, J_psi_m, 0, 0, 0, 0, 0, 0],
        [0, 0, (1/2)*I_II*cos(2*theta_p) + (1/2)*I_II + I_IK*sin(2*theta_p) - 1/2*I_KK*cos(2*theta_p) + (1/2)*I_KK + J_ax + J_psi_L, 0, 0, 0, 0, I_IJ*cos(theta_p) + I_JK*sin(theta_p)],
        [0, 0, 0, J_theta_e, 0, 0, 0, 0],
        [0, 0, 0, 0, J_theta_h, 0, 0, 0],
        [0, 0, 0, 0, 0, J_theta_m, 0, 0],
        [0, 0, 0, 0, 0, 0, J_theta_r, 0],
        [0, 0, I_IJ*cos(theta_p) + I_JK*sin(theta_p), 0, 0, 0, 0, I_JJ + J_theta_p],
    ])


def termini(q, dq, p=PARAMETRI):
    """h(q, dq): giroscopici, centrifughi, elastici, dissipativi, gravita'."""
    psi_e, psi_m, psi_L, theta_e, theta_h, theta_m, theta_r, theta_p = q
    dpsi_e, dpsi_m, dpsi_L, dtheta_e, dtheta_h, dtheta_m, dtheta_r, dtheta_p = dq
    J_psi_e = p['J_psi_e']
    J_psi_m = p['J_psi_m']
    J_psi_L = p['J_psi_L']
    J_theta_e = p['J_theta_e']
    J_theta_h = p['J_theta_h']
    J_theta_m = p['J_theta_m']
    J_theta_r = p['J_theta_r']
    J_theta_p = p['J_theta_p']
    k_psi_1 = p['k_psi_1']
    k_psi_2 = p['k_psi_2']
    k_theta_1 = p['k_theta_1']
    k_theta_2 = p['k_theta_2']
    k_theta_3 = p['k_theta_3']
    k_theta_4 = p['k_theta_4']
    c_psi_1 = p['c_psi_1']
    c_psi_2 = p['c_psi_2']
    c_theta_1 = p['c_theta_1']
    c_theta_2 = p['c_theta_2']
    c_theta_3 = p['c_theta_3']
    c_theta_4 = p['c_theta_4']
    b_psi_e = p['b_psi_e']
    b_psi_m = p['b_psi_m']
    b_psi_L = p['b_psi_L']
    b_theta_e = p['b_theta_e']
    b_theta_h = p['b_theta_h']
    b_theta_m = p['b_theta_m']
    b_theta_r = p['b_theta_r']
    b_theta_p = p['b_theta_p']
    I_II = p['I_II']
    I_JJ = p['I_JJ']
    I_KK = p['I_KK']
    I_IJ = p['I_IJ']
    I_IK = p['I_IK']
    I_JK = p['I_JK']
    d_I = p['d_I']
    d_J = p['d_J']
    d_K = p['d_K']
    m_p = p['m_p']
    g = p['g']
    J_ax = p['J_ax']
    return np.array([
        b_psi_e*dpsi_e + c_psi_1*dpsi_e - c_psi_1*dpsi_m + k_psi_1*psi_e - k_psi_1*psi_m,
        b_psi_m*dpsi_m - c_psi_1*dpsi_e + c_psi_1*dpsi_m - c_psi_2*dpsi_L + c_psi_2*dpsi_m - k_psi_1*psi_e + k_psi_1*psi_m - k_psi_2*psi_L + k_psi_2*psi_m,
        -I_II*dpsi_L*dtheta_p*sin(2*theta_p) - I_IJ*dtheta_p**2*sin(theta_p) + 2*I_IK*dpsi_L*dtheta_p*cos(2*theta_p) + I_JK*dtheta_p**2*cos(theta_p) + I_KK*dpsi_L*dtheta_p*sin(2*theta_p) + b_psi_L*dpsi_L + c_psi_2*dpsi_L - c_psi_2*dpsi_m + k_psi_2*psi_L - k_psi_2*psi_m,
        b_theta_e*dtheta_e + c_theta_1*dtheta_e - c_theta_1*dtheta_h + k_theta_1*theta_e - k_theta_1*theta_h,
        b_theta_h*dtheta_h - c_theta_1*dtheta_e + c_theta_1*dtheta_h + c_theta_2*dtheta_h - c_theta_2*dtheta_m + c_theta_4*dtheta_h - c_theta_4*dtheta_p - k_theta_1*theta_e + k_theta_1*theta_h + k_theta_2*theta_h - k_theta_2*theta_m + k_theta_4*theta_h - k_theta_4*theta_p,
        b_theta_m*dtheta_m - c_theta_2*dtheta_h + c_theta_2*dtheta_m + c_theta_3*dtheta_m - c_theta_3*dtheta_r - k_theta_2*theta_h + k_theta_2*theta_m + k_theta_3*theta_m - k_theta_3*theta_r,
        b_theta_r*dtheta_r - c_theta_3*dtheta_m + c_theta_3*dtheta_r - k_theta_3*theta_m + k_theta_3*theta_r,
        (1/2)*I_II*dpsi_L**2*sin(2*theta_p) - I_IK*dpsi_L**2*cos(2*theta_p) - 1/2*I_KK*dpsi_L**2*sin(2*theta_p) + b_theta_p*dtheta_p - c_theta_4*dtheta_h + c_theta_4*dtheta_p - d_I*g*m_p*sin(theta_p) + d_K*g*m_p*cos(theta_p) - k_theta_4*theta_h + k_theta_4*theta_p,
    ])


E = np.array([[-0., -0.],
 [ 1., -0.],
 [-0., -0.],
 [-0., -0.],
 [-0., -0.],
 [-0.,  1.],
 [-0., -0.],
 [-0., -0.]])


def rhs(t, x, tau, p=PARAMETRI):
    """dx/dt per scipy.integrate.solve_ivp; tau(t) restituisce [tau_psi, tau_theta]."""
    q, dq = x[:8], x[8:]
    ddq = np.linalg.solve(massa(q, p), E @ np.asarray(tau(t)) - termini(q, dq, p))
    return np.concatenate([dq, ddq])
