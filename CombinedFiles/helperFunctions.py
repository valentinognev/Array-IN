import numpy as np
from scipy.linalg import qr, svd, eig
from scipy.spatial.transform import Rotation as R
import warnings


def commutation_matrix(m, n):
    """COMMUTATIONMATRIX Summary of this function goes here
    Detailed explanation goes here
    [m, n] = size(A);
    """
    # [m, n] = size(A);
    I = np.reshape(np.arange(1, m*n + 1), (m, n), order='F')  # initialize a matrix of indices of size(A)
    I = I.T  # Transpose it
    I = I.flatten(order='F')  # vectorize the required indices
    Y = np.eye(m*n)  # Initialize an identity matrix
    Y = Y[I-1, :]  # Re-arrange the rows of the identity matrix (adjust for 0-indexing)
    return Y


def compute_As(r_tot):
    """COMPUTE_AS The matrix for the rotation 
    
    As = compute_As(r_tot)
    Where r_tot is centered.
    """
    assert np.all(np.abs(np.mean(r_tot, axis=1)) < 10*np.finfo(float).eps)
    K = r_tot.shape[1]  # Number of acc triads

    R_skew = np.zeros((3, 3, K))
    for k in range(K):
        R_skew[:, :, k] = skew_sym(r_tot[:, k])
    
    R_square = np.zeros((3, 3))
    for k in range(K):
        R_square = R_square + R_skew[:, :, k].T @ R_skew[:, :, k]
    
    As = np.zeros((3, 3, K))
    for k in range(K):
        As[:, :, k] = np.linalg.solve(R_square, R_skew[:, :, k])
    
    return As, R_square


def compute_projection_matrix(H, Q):
    """COMPUTE_PROJECTION_MATRIX Summary of this function goes here
    Detailed explanation goes here
    """
    H_t_Q_inv = H.T @ np.linalg.inv(Q)
    M = np.linalg.solve(H_t_Q_inv @ H, H_t_Q_inv)
    return M


def factorize_T(T):
    """FACTORIZE_T Factorize scale matrix as T = D*L*Q
    D: Diagonal matrix with scale factors
    L: Upper triangular matrix with unit diagonal. Account for
    non-orthogonalities 
    Q: Rotation matrix 
    """
    R_mat, Q = rq(T)

    S = np.diag(np.diag(R_mat))
    D = np.diag(S)
    L = np.linalg.solve(D, R_mat)

    return D, L, Q


def get_spherical_motion(t, inp):
    """UNTITLED Summary of this function goes here
    Detailed explanation goes here
    """
    
    if inp['phi'] == "sinus":
        phi, phi_dot, phi_dot_2 = get_sinus(t, inp['phi_params'])
    elif inp['phi'] == "linear":
        phi, phi_dot, phi_dot_2 = get_linear(t, inp['phi_params'])
    elif inp['phi'] == "quadratic":
        phi, phi_dot, phi_dot_2 = get_quadratic(t, inp['phi_params'])
    elif inp['phi'] == "poly":
        phi, phi_dot, phi_dot_2 = get_polynomial(t, inp['phi_params'])
    elif inp['phi'] == "constant":
        phi, phi_dot, phi_dot_2 = get_constant(t, inp['phi_params'])
    else:
        raise ValueError("No correct motion for phi")

    if inp['theta'] == "sinus":
        theta, theta_dot, theta_dot_2 = get_sinus(t, inp['theta_params'])
    elif inp['theta'] == "linear":
        theta, theta_dot, theta_dot_2 = get_linear(t, inp['theta_params'])
    elif inp['theta'] == "quadratic":
        theta, theta_dot, theta_dot_2 = get_quadratic(t, inp['theta_params'])
    elif inp['theta'] == "poly":
        theta, theta_dot, theta_dot_2 = get_polynomial(t, inp['theta_params'])
    elif inp['theta'] == "constant":
        theta, theta_dot, theta_dot_2 = get_constant(t, inp['theta_params'])
    else:
        raise ValueError("No correct motion for theta")
    
    m = {
        'phi': phi,
        'phi_dot': phi_dot,
        'phi_dot_2': phi_dot_2,
        'theta': theta,
        'theta_dot': theta_dot,
        'theta_dot_2': theta_dot_2
    }
    return m


def get_constant(t, inp):
    """Get constant motion parameters"""
    A = inp.get('A', 1)
    
    s = A * np.ones_like(t)
    s_dot = np.zeros_like(t)
    s_dot_2 = np.zeros_like(t)
    
    return s, s_dot, s_dot_2


def get_sinus(t, inp):
    """Get sinusoidal motion parameters"""
    A = inp.get('A', 1)
    f = inp.get('f', 1)
    b = inp.get('b', 0)
    
    s = A * np.sin(2*np.pi*f*t) + b
    s_dot = A * np.cos(2*np.pi*f*t) * 2*np.pi*f
    s_dot_2 = -A * np.sin(2*np.pi*f*t) * (2*np.pi*f)**2
    
    return s, s_dot, s_dot_2


def get_linear(t, inp):
    """Get linear motion parameters"""
    A = inp.get('A', 1)
    
    s = A * t
    s_dot = A * np.ones_like(t)
    s_dot_2 = np.zeros_like(t)
    
    return s, s_dot, s_dot_2


def get_quadratic(t, inp):
    """Get quadratic motion parameters"""
    A = inp.get('A', 1)
    
    s = A * t**2
    s_dot = 2 * A * t
    s_dot_2 = 2 * A * np.ones_like(t)
    
    return s, s_dot, s_dot_2


def get_polynomial(t, inp):
    """Get polynomial motion parameters"""
    p = inp['p']
    
    s = np.polyval(p, t)
    p1 = np.polyder(p)
    s_dot = np.polyval(p1, t)
    p2 = np.polyder(p1)
    s_dot_2 = np.polyval(p2, t)
    
    return s, s_dot, s_dot_2


def get_T_components(T):
    """GET_T_COMPONENTS Get the components of the T matrix
    s: scale factors 
    l: angles for non-orthogonalities
    q: rotation vector 
    """
    S, L, Q = factorize_T(T)

    s = np.diag(S)

    l = np.zeros(3)
    l[0] = L[0, 1]
    l[1] = L[0, 2]
    l[2] = L[1, 2]

    q = logSO3(Q)
    return s, l, q


def initialAttitude2Rotm(u):
    """INITIAL_ATTITUDE Summary of this function goes here
    Detailed explanation goes here
    """
    
    f_u = np.mean(u[0, :])
    f_v = np.mean(u[1, :])
    f_w = np.mean(u[2, :])

    roll = np.arctan2(-f_v, -f_w)
    pitch = np.arctan2(f_u, np.sqrt(f_v**2 + f_w**2))
    heading = 0

    # Convert Euler angles to rotation matrix (XYZ convention)
    r = R.from_euler('XYZ', [roll, pitch, heading])
    R_matrix = r.as_matrix()

    return R_matrix


def norm_time(x):
    """NORM_TIME Summary of this function goes here
    Detailed explanation goes here
    """
    norm_v = np.sqrt(np.sum(x**2, axis=0))
    return norm_v


def ppdiff(pp, j=1):
    """PPDIFF Differentiate piecewise polynomial.
    QQ = PPDIFF(PP,J) returns the J:th derivative of a piecewise
    polynomial PP. PP must be on the form evaluated by PPVAL. QQ is a
    piecewise polynomial on the same form. Default value for J is 1.

    Example:
        x = linspace(-pi,pi,9);
        y = sin(x);
        pp = spline(x,y);
        qq = ppdiff(pp);
        xx = linspace(-pi,pi,201);
        plot(xx,cos(xx),'b',xx,ppval(qq,xx),'r')

    See also PPVAL, SPLINE, SPLINEFIT, PPINT

    Author: Jonas Lundgren <splinefit@gmail.com> 2009
    """
    
    # Check diff order
    if not np.isreal(j) or j != int(j) or j < 0:
        raise ValueError('Order of derivative must be a non-negative integer!')

    # Get coefficients
    coefs = pp['coefs'].copy()
    m, n = coefs.shape

    if j == 0:
        # Do nothing
        pass
    elif j < n:
        # Derivative of order J
        D = np.zeros((j, n-j))
        D[0, :] = np.arange(n-j, 0, -1)
        for i in range(1, j):
            D[i, :] = 1
        D = np.cumprod(D, axis=0)
        D = np.prod(D, axis=0)
        
        coefs = coefs[:, :n-j]
        for k in range(n-j):
            coefs[:, k] = D[k] * coefs[:, k]
    else:
        # Derivative kills PP
        coefs = np.zeros((m, 1))

    # Set output
    qq = pp.copy()
    qq['coefs'] = coefs
    qq['order'] = coefs.shape[1]
    
    return qq


def q_minus(q1, q2):
    """UNTITLED17 Summary of this function goes here
    Detailed explanation goes here
    Half the angle 
    """
    # Note: This function appears to work with quaternions
    # The implementation would depend on the specific quaternion library being used
    # This is a placeholder implementation
    warnings.warn("q_minus function needs quaternion library implementation")
    
    # Placeholder implementation - would need proper quaternion operations
    dtheta = np.zeros((len(q1), 3))
    return dtheta


def rq(T):
    """RQ RQ factorization
    Same as QR and R have positive diagonals 
    """
    
    Q, R_temp = qr(np.flipud(T).T)
    Q = np.fliplr(Q)  # Upper triangularize T{1} from the left
    Q = Q @ np.diag(np.diag(np.sign(Q)))  # To not change coordinate system orientation
    Q = Q.T

    R_matrix = T @ Q.T

    return R_matrix, Q


def skew_sym(a):
    """Create skew symmetric matrix from vector"""
    A = np.array([
        [0,     -a[2],  a[1]],
        [a[2],   0,    -a[0]],
        [-a[1],  a[0],   0]
    ])
    return A


def solve_Wahbas_problem(W, V):
    """solve_Wahbas_problem Estimate initial rotation matrix from gravity 
    R = argmin sum_{k,n} || w_{k,n} - R*v_{k,n}||^2
    where R in SO(3).
    """
    assert W.shape == V.shape
    N = W.shape[1]
    K = W.shape[0] // 3
    B = np.zeros((3, 3))
    inds = np.reshape(np.arange(3*K), (3, K), order='F')
    
    for n in range(N):
        for k in range(K):
            kk = inds[:, k]
            B = B + np.outer(W[kk, n], V[kk, n])

    U, _, Vt = svd(B)
    V_svd = Vt.T
    M = np.eye(3)
    M[2, 2] = np.linalg.det(U) * np.linalg.det(V_svd)
    R_matrix = U @ M @ V_svd.T

    return R_matrix


def stochastic_observability(P):
    """STOCHASTIC_OBSERVABILITY Summary of this function goes here
    Detailed explanation goes here
    """
    
    assert P.ndim == 3

    P0 = P[:, :, 0]
    n = P.shape[0]

    assert np.allclose(P0, np.diag(np.diag(P0)))  # Check if diagonal

    F = np.linalg.inv(np.sqrt(P0))

    P_norm = np.zeros_like(P)

    for k in range(P.shape[2]):
        P_k_in = P[:, :, k]
        
        P_k_in_1 = F @ P_k_in @ F

        # Normalize to unit norm for the eigen values 
        P_norm_k = P_k_in_1 / np.trace(P_k_in_1)
        
        P_norm[:, :, k] = P_norm_k

    # eigenshuffle: Consistent sorting for an eigenvalue/vector sequence
    # [Vseq,Dseq] = eigenshuffle(Asequence)
    # Note: eigenshuffle is not a standard function, would need custom implementation
    warnings.warn("eigenshuffle function needs custom implementation")
    
    # Placeholder - compute eigenvalues and eigenvectors for each slice
    V = np.zeros((n, n, P.shape[2]))
    D = np.zeros((n, n, P.shape[2]))
    
    for k in range(P.shape[2]):
        eigenvals, eigenvecs = eig(P_norm[:, :, k])
        D[:, :, k] = np.diag(eigenvals)
        V[:, :, k] = eigenvecs

    return D, V, P_norm


def logSO3(R_matrix):
    """Logarithm map for SO(3) - convert rotation matrix to rotation vector"""
    # Using scipy's Rotation class
    r = R.from_matrix(R_matrix)
    rotvec = r.as_rotvec()
    return rotvec 