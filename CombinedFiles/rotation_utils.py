import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import interp1d
import warnings


def average_rotation(Rs, nb_it_max=20, tol_r=1e-10):
    """Average rotation from multiple rotation matrices
    
    Args:
        Rs: Rotation matrices of shape (3, 3, N), N number of rotations
        nb_it_max: max number of iterations, default 20 
        tol_r: tolerance for deviations, default 1e-10
    
    Returns:
        R_mean: Average rotation 
        list_r: residuals in lie algebra 
    """
    
    number_rotations = Rs.shape[2]
    R_mean = Rs[:, :, 0]  # First approx of R [1]
    
    for nb_it in range(1, nb_it_max + 1):  # [2]
        list_r = np.full((3, number_rotations), np.nan)  # [3]
        for i in range(number_rotations):
            list_r[:, i] = logSO3(R_mean.T @ Rs[:, :, i])
        
        r = np.mean(list_r, axis=1)
        
        print(f"{nb_it}/{nb_it_max}: tol: {np.linalg.norm(r):.3e} / {tol_r:.3e}")
        if np.linalg.norm(r) < tol_r:  # [4]
            break
        R_mean = R_mean @ expSO3(r)  # Update [7]
    # [8]
    
    if nb_it == nb_it_max:
        raise RuntimeError('the maximum number of iteration where reached')
    
    return R_mean, list_r


def interpolate_rotation(t, R, t_inter):
    """Interpolate IMU pos and rotation estimates 
    
    Args:
        t: time points corresponding to R data points 
        R: rotation matrices of shape (3, 3, N)
        t_inter: time points where interpolation should occur
    
    Returns:
        R_inter: interpolated rotation matrices
    """
    
    S = R.shape
    assert len(t) == S[2]
    assert S[0] == 3 and S[1] == 3
    
    # Find the fractional indices using linear interpolation 
    interp_func = interp1d(t, np.arange(len(t)), kind='linear', 
                          bounds_error=False, fill_value=np.nan)
    inds_imu_time = interp_func(t_inter)
    
    R_inter = np.full((3, 3, len(inds_imu_time)), np.nan)
    
    for n in range(len(inds_imu_time)):
        # Extrapolation set to NaN
        if np.isnan(inds_imu_time[n]):
            continue
        
        frac = inds_imu_time[n] - np.floor(inds_imu_time[n])
        if frac > 1:
            warnings.warn("fraction larger than 1")
        elif frac == 0:
            # Same point in time
            R_inter[:, :, n] = R[:, :, int(np.round(inds_imu_time[n]))]
        else:
            # Calculate the rotation vector and scale it
            left_ind = int(np.floor(inds_imu_time[n]))
            right_ind = left_ind + 1
            R_left = R[:, :, left_ind]
            R_right = R[:, :, right_ind]
            
            theta = logSO3(invSO3(R_left) @ R_right)
            R_inter[:, :, n] = R_left @ expSO3(frac * theta)
    
    return R_inter


def my_rotm2eul(R_mat):
    """Rotation matrix to euler angles [roll, pitch, yaw]
    
    Args:
        R_mat: Rotation matrix (3, 3, N) 
    
    Returns:
        E: Euler angles (3, N)
    
    Roll: around x-axis
    Pitch: around y-axis
    Yaw: around z-axis (heading)
    """
    
    # Convert to scipy Rotation object and get euler angles
    # rotm2eul gives [yaw, pitch, roll] intrinsic rotation
    # R = R_z(yaw)*R_y(pitch)*R_z(roll)
    # unwrap: adds 2pi when wrapping 
    # flipud to get in order [roll, pitch, yaw]
    
    if R_mat.ndim == 2:
        # Single rotation matrix
        rot = R.from_matrix(R_mat)
        angles = rot.as_euler('ZYX')  # [yaw, pitch, roll]
        E = np.flip(angles)  # [roll, pitch, yaw]
    else:
        # Multiple rotation matrices
        E = np.zeros((3, R_mat.shape[2]))
        for i in range(R_mat.shape[2]):
            rot = R.from_matrix(R_mat[:, :, i])
            angles = rot.as_euler('ZYX')  # [yaw, pitch, roll]
            E[:, i] = np.flip(angles)  # [roll, pitch, yaw]
        
        # Apply unwrap to handle angle wrapping
        E = np.unwrap(E, axis=1)
    
    return E


def R2w_central_diff(R, t):
    """Rotation matrix 2 angular velocity using central difference 
    
    Args:
        R: Rotation matrices (3, 3, N)
        t: time vector
    
    Returns:
        w: angular velocity in body frame
    
    Based on:
    R_{t+1} = R_{t} exp_SO3(w*t)
    w in body frame
    """
    
    w = np.full((3, len(t)), np.nan)
    
    for n in range(1, len(t) - 1):
        dt = t[n+1] - t[n-1]
        w[:, n] = logSO3(R[:, :, n-1].T @ R[:, :, n+1]) / dt
    
    return w


def rotationMatrixFromTwoUnitVectors(a, b):
    """Find rotation matrix from a to b
    
    Args:
        a: first unit vector
        b: second unit vector
    
    Returns:
        R: rotation matrix that rotates a to b
    """
    
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = np.dot(a, b)
    
    R = np.eye(3) + skew_sym(v) + np.linalg.matrix_power(skew_sym(v), 2) * (1-c)/(s**2)
    
    return R


# Helper functions for SO(3) operations
def logSO3(R):
    """Logarithm map from SO(3) to so(3)"""
    theta = np.arccos((np.trace(R) - 1) / 2)
    if np.abs(theta) < 1e-6:
        return np.array([R[2,1] - R[1,2], R[0,2] - R[2,0], R[1,0] - R[0,1]]) / 2
    else:
        return theta / (2 * np.sin(theta)) * np.array([R[2,1] - R[1,2], R[0,2] - R[2,0], R[1,0] - R[0,1]])


def expSO3(omega):
    """Exponential map from so(3) to SO(3)"""
    theta = np.linalg.norm(omega)
    if theta < 1e-6:
        return np.eye(3) + skew_sym(omega)
    else:
        omega_hat = skew_sym(omega)
        return np.eye(3) + np.sin(theta)/theta * omega_hat + (1-np.cos(theta))/(theta**2) * omega_hat @ omega_hat


def invSO3(R):
    """Inverse of rotation matrix (transpose)"""
    return R.T


def skew_sym(v):
    """Skew symmetric matrix from 3D vector"""
    return np.array([[0, -v[2], v[1]],
                     [v[2], 0, -v[0]],
                     [-v[1], v[0], 0]])
