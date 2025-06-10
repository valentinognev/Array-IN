import numpy as np
from scipy.linalg import svd


def adj_so3(w):
    """Adjoint representation for SO(3)"""
    adw = hat_so3(w)
    return adw


def ad_so3(R):
    """Adjoint representation for SO(3)"""
    AdR = R
    return AdR


def dual_hat_so3(x):
    """
    Dual hat operator for SO(3)
    w_hat*x = x_dual_hat*w
    """
    x_dual_hat = -hat_so3(x)
    return x_dual_hat


def error_so3(Rhat, Rtrue):
    """
    Calculate errors between estimated and true rotations
    Need to be same length 
    """
    N = Rhat.shape[2] if Rhat.ndim == 3 else 1
    
    # Calculate errors
    err = np.zeros((3, N))
    for n in range(N):
        if Rtrue.ndim == 3:
            err[:, n] = log_so3(inv_so3(Rhat[:, :, n]) @ Rtrue[:, :, n])
        else:
            err[:, n] = log_so3(inv_so3(Rhat[:, :, n]) @ Rtrue)
    
    return err


def exp_so3(w):
    """Exponential map from so(3) to SO(3)"""
    normw = np.linalg.norm(w)
    
    if normw == 0:
        R = np.eye(3)
        return R
    
    w_hat = hat_so3(w)
    
    R = np.eye(3) + np.sin(normw) * w_hat / normw + (1 - np.cos(normw)) * w_hat @ w_hat / (normw**2)
    
    return R


def hat_so3(w):
    """
    Hat operator: converts 3D vector to skew-symmetric matrix
    """
    w_hat = np.array([[0, -w[2], w[1]],
                      [w[2], 0, -w[0]],
                      [-w[1], w[0], 0]])
    
    # Alternative implementation (commented out in original):
    # w_hat = np.zeros((3, 3))
    # w_hat[1, 0] = w[2]
    # w_hat[2, 0] = -w[1]
    # 
    # w_hat[0, 1] = -w[2]
    # w_hat[2, 1] = w[0]
    # 
    # w_hat[0, 2] = w[1]
    # w_hat[1, 2] = -w[0]
    
    return w_hat


def inv_so3(R):
    """Inverse of SO(3) matrix"""
    error_flag = 0
    Rinv = R.T
    return Rinv, error_flag


def log_so3(R):
    """Logarithm map from SO(3) to so(3)"""
    phy = np.arccos((np.trace(R) - 1) / 2)
    if abs(phy) > np.pi:
        raise ValueError('angle supérieur à pi')
    
    if phy == 0:
        w = np.zeros(3)
    elif abs(phy) == np.pi:
        A = (R - np.eye(3)) / 2
        w1 = np.sqrt(-((A[1, 1] + A[2, 2] - A[0, 0]) / 2))
        
        w2 = np.sqrt(-((A[0, 0] + A[2, 2] - A[1, 1]) / 2))
        w3 = np.sqrt(-((A[0, 0] + A[1, 1] - A[2, 2]) / 2))
        
        if w1 != 0:
            if A[0, 1] < 0:
                w2 = -w2
            if A[0, 2] < 0:
                w3 = -w3
        elif w2 != 0:
            if A[1, 2] < 0:
                w3 = -w3
        
        w = np.array([w1, w2, w3]) * phy
    else:
        # on remultiplie par phy pour retrouver le vecteur avec sa norme originale
        w_hat = (R - R.T) / (2 * np.sin(phy)) * phy
        w = vec_so3(w_hat)
    
    error_flag = 0
    return w, error_flag


def normalize_so3(R):
    """Normalize matrix to proper SO(3) using SVD"""
    u, s, v = svd(R)
    Rnorm = u @ v
    return Rnorm


def phi_so3(w):
    """
    Left-Jacobian to SO(3)
    sum_k 1/(k + 1)! ad(w)^k
    """
    normw = np.linalg.norm(w)
    
    # if normw > np.pi/2:
    #     raise ValueError('formula not sure')
    
    if normw > 0:
        adw = adj_so3(w)
        
        Phiw = (np.eye(3) + 
                (1 / (2 * normw**2)) * (4 - normw * np.sin(normw) - 4 * np.cos(normw)) * adw +
                (1 / (2 * normw**3)) * (4 * normw - 5 * np.sin(normw) + normw * np.cos(normw)) * adw @ adw +
                (1 / (2 * normw**4)) * (2 - normw * np.sin(normw) - 2 * np.cos(normw)) * np.linalg.matrix_power(adw, 3) +
                (1 / (2 * normw**5)) * (2 * normw - 3 * np.sin(normw) + normw * np.cos(normw)) * np.linalg.matrix_power(adw, 4))
    else:
        Phiw = np.eye(3)
    
    return Phiw


def vec_so3(w_hat):
    """
    Vec operator: converts skew-symmetric matrix to 3D vector
    Inverse of hat operator
    """
    w = np.array([w_hat[2, 1], w_hat[0, 2], w_hat[1, 0]])
    return w