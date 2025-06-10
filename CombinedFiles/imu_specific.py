import numpy as np
from scipy.interpolate import splrep, splev, spalde
from scipy.linalg import inv


def dAdr(a, r, w):
    """
    UNTITLED Summary of this function goes here
    Detailed explanation goes here
    """
    
    N_a = r.shape[1]  # Number of acc triads
    
    R_skew = []
    for k in range(N_a):
        R_skew.append(skew_sym(r[:, k]))
    
    H = np.vstack([-np.vstack(R_skew), np.tile(np.eye(3), (N_a, 1))])
    
    O2 = np.linalg.matrix_power(skew_sym(w), 2)
    h = (O2 @ r).flatten('F')  # Flatten in column-major order like MATLAB
    a_1 = a - h
    A_1 = H.T @ H
    A_1_inv = inv(A_1)
    
    f = np.linalg.solve(A_1, H.T @ a_1)
    
    assert len(f) == 6
    A_2 = np.kron(np.eye(N_a), O2)
    A_3 = np.kron(f.T, -A_1_inv)
    A_4 = np.kron(a_1.T, A_1_inv)
    A_5 = -np.linalg.solve(A_1, H.T @ A_2)
    
    K = commutation_matrix(3*N_a, 6)
    A_6 = np.kron(H.T, np.eye(6)) @ K
    A_7 = np.kron(np.eye(6), H.T)
    A_8 = A_3 @ A_6 + A_3 @ A_7 + A_4 @ K
    
    B_1 = skew_sym(np.array([-1, 0, 0]))
    B_2 = skew_sym(np.array([0, -1, 0]))
    B_3 = skew_sym(np.array([0, 0, -1]))
    
    A_9 = np.vstack([
        np.kron(np.eye(N_a), -B_1),
        np.kron(np.eye(N_a), -B_2),
        np.kron(np.eye(N_a), -B_3),
        np.zeros((9*N_a, 3*N_a))
    ])
    
    J = A_8 @ A_9 + A_5
    
    assert J.shape == (6, 3*N_a)
    
    return f, J


def d_h_d_omega(omega, r):
    """
    D_H_D_OMEGA Summary of this function goes here
    Detailed explanation goes here
    """
    N_a = r.shape[1]
    d_h_d_omega_parts = []
    omega_hat = HatSO3(omega)
    
    for k in range(N_a):
        r_k = r[:, k]
        d_h_d_omega_parts.append(-HatSO3(omega_hat @ r_k) - omega_hat @ HatSO3(r_k))
    
    d_h_d_omega = np.vstack(d_h_d_omega_parts)
    
    return d_h_d_omega


def d_h_d_omega_opt(w, r):
    """
    D_H_D_OMEGA Summary of this function goes here
    Detailed explanation goes here
    """
    N_a = r.shape[1]
    row1 = np.arange(0, 3*N_a, 3)  # 0-based indexing in Python
    row2 = np.arange(1, 3*N_a, 3)
    row3 = np.arange(2, 3*N_a, 3)
    
    d_h_d_omega = np.zeros((3*N_a, 3))
    r1 = r[0, :]
    r2 = r[1, :]
    r3 = r[2, :]
    
    r1w1 = w[0] * r1
    r1w2 = r1 * w[1]
    r1w3 = r1 * w[2]
    
    r2w1 = r2 * w[0]
    r2w2 = w[1] * r2
    r2w3 = r2 * w[2]
    
    r3w1 = r3 * w[0]
    r3w2 = r3 * w[1]
    r3w3 = w[2] * r3
    
    d_h_d_omega[row1, 0] = r2w2 + r3w3
    d_h_d_omega[row2, 0] = r1w2 - 2*r2w1
    d_h_d_omega[row3, 0] = r1w3 - 2*r3w1
    
    d_h_d_omega[row1, 1] = r2w1 - 2*r1w2
    d_h_d_omega[row2, 1] = r1w1 + r3w3
    d_h_d_omega[row3, 1] = r2w3 - 2*r3w2
    
    d_h_d_omega[row1, 2] = r3w1 - 2*r1w3
    d_h_d_omega[row2, 2] = r3w2 - 2*r2w3
    d_h_d_omega[row3, 2] = r1w1 + r2w2
    
    return d_h_d_omega


def get_norm_g_kth():
    """
    Returns normalized gravity constant
    """
    g = 9.8183037  # Lantmateriet, m/s^2
    return g


def get_triad_form(x):
    """
    GET_TRIAD_FORM Summary of this function goes here
    Detailed explanation goes here
    """
    N_sens = x.shape[0]
    assert N_sens % 3 == 0
    N_imu = N_sens // 3
    
    y = x.reshape(3, N_imu, -1, order='F')  # Fortran order like MATLAB
    
    return y


def gravity(lmbda, h):
    """
    function g=gravity(lambda,h)
    
    function for calculation of the local gravity vector, in
    the geographic reference frame (same as tangent plane is 
    stationary).
    
    Based upon the WGS_84 Geodetic and Gravity model. For more 
    info see [pp 222-223,1].
    
    lambda -> Latitude [degrees]
    h -> Altitude [m]
    g 
    
    edit: Isaac Skog, 2006-08-17
    """
    
    # degrees to radians
    lmbda = np.pi/180 * lmbda
    
    gamma = 9.780327 * (1 + 0.0053024*np.sin(lmbda)**2 - 0.0000058*np.sin(2*lmbda)**2)
    
    g = gamma - ((3.0877e-6) - (0.004e-6)*np.sin(lmbda)**2)*h + (0.072e-12)*h**2
    
    g = np.array([0, 0, -g])
    
    return g


def stationaryAcc2rollPitch(u):
    """
    GET_ROLL_PITCH Summary of this function goes here
    Detailed explanation goes here
    """
    f_x = np.mean(u[0, :])
    f_y = np.mean(u[1, :])
    f_z = np.mean(u[2, :])
    
    roll = np.arctan2(-f_x, f_z)
    pitch = np.arctan2(f_y, np.sqrt(f_x**2 + f_z**2))
    
    return roll, pitch


def stationaryAcc2rollPitch_IS(u):
    """
    GET_ROLL_PITCH Summary of this function goes here
    Detailed explanation goes here
    """
    f_u = np.mean(u[0, :])
    f_v = np.mean(u[1, :])
    f_w = np.mean(u[2, :])
    
    roll = np.arctan2(-f_v, -f_w)
    pitch = np.arctan2(f_u, np.sqrt(f_v**2 + f_w**2))
    
    return roll, pitch


def triad_mean(x):
    """
    TRIAD_MEAN mean of triad data from 2D matrix
    y = triad_mean(x)
    """
    N_sens = x.shape[0]
    assert N_sens % 3 == 0
    N_imu = N_sens // 3
    
    y = np.mean(x.reshape(3, N_imu, -1, order='F'), axis=1).reshape(3, -1, order='F')
    
    return y


def triad_norm(x):
    """
    TRIAD_NORM Summary of this function goes here
    Detailed explanation goes here
    """
    N_sens = x.shape[0]
    assert N_sens % 3 == 0
    N_imu = N_sens // 3
    
    y = np.sqrt(np.sum(get_triad_form(x)**2, axis=0)).reshape(N_imu, -1, order='F')
    
    return y


def w2w_dot_splines(t, w):
    """
    W2W_DOT_SPLINES Interpolate w to w_dot using splines
    
    w_dot_interp = w2w_dot_splines(t, w)
    """
    w_dot_interp = np.zeros_like(w)
    
    for i in range(3):
        # Create spline representation
        tck = splrep(t, w[i, :], s=0)  # s=0 for interpolating spline
        
        # Compute derivative by evaluating spline derivative
        w_dot_interp[i, :] = splev(t, tck, der=1)
    
    return w_dot_interp


# Helper functions that need to be implemented or imported
# These functions are referenced but not defined in the original MATLAB file

def skew_sym(v):
    """
    Create skew-symmetric matrix from 3D vector
    Note: This function needs to be implemented based on your specific requirements
    """
    # Placeholder implementation - you'll need to implement this based on your needs
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])


def commutation_matrix(m, n):
    """
    Create commutation matrix
    Note: This function needs to be implemented based on your specific requirements
    """
    # Placeholder implementation - you'll need to implement this based on your needs
    K = np.zeros((m*n, m*n))
    # Implementation would go here
    return K


def HatSO3(v):
    """
    Hat operator for SO(3) - creates skew-symmetric matrix
    Note: This function needs to be implemented based on your specific requirements
    """
    # This is likely the same as skew_sym
    return skew_sym(v)
