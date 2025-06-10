import numpy as np
from scipy.linalg import cholesky, inv
from scipy.interpolate import interp1d
import warnings

# Assuming these functions are available from other modules
# from your_so3_module import errorSO3, logSO3, expSO3, invSO3
# from your_matrix_utils import matrix3d2blkdiag

def calculate_trajectory_error(S_hat, S_true, *args):
    """CALCULATE_TRAJECTORY_ERROR Summary of this function goes here
       Detailed explanation goes here"""
    
    if len(args) > 0:
        section = args[0]
    else:
        section = range(S_true['R'].shape[2])
    
    err = {}
    err['R'] = errorSO3(S_hat['R'], S_true['R'][:, :, section])
    err['R_deg'] = np.rad2deg(err['R'])
    
    if 'w' in S_hat:
        err['w'] = S_hat['w'] - S_true['w'][:, section]
        err['w_deg'] = np.rad2deg(err['w'])
    
    err['p'] = S_hat['p'] - S_true['p'][:, section]
    err['v'] = S_hat['v'] - S_true['v'][:, section]
    
    if 'b_a' in S_hat:
        inds_bias = range(S_hat['b_a'].shape[0])
        err['b_a'] = S_hat['b_a'] - S_true['b_a'][inds_bias, section]
    
    if 'b_g' in S_hat:
        err['b_g'] = S_hat['b_g'] - S_true['b_g'][:, section]
        err['b_g_deg'] = np.rad2deg(err['b_g'])
    
    return err


def compensate_covariance(Q_y, T):
    """COMPENSATE_MEASUREMENTS Summary of this function goes here
       Detailed explanation goes here"""
    
    L = cholesky(Q_y, lower=True)
    
    q = np.linalg.solve(T, L)
    
    Q_u = q @ q.T
    
    return Q_u


def compensate_measurements(y, T, b):
    """COMPENSATE_MEASUREMENTS Summary of this function goes here
       Detailed explanation goes here"""
    
    if T.ndim == 3:
        T_diag = matrix3d2blkdiag(T)
    else:
        T_diag = T
    
    b = b.reshape(-1, 1)
    u = np.linalg.solve(T_diag, (y - b))
    
    return u


def compute_error(S, S_ref):
    """COMPUTE_ERROR Summary of this function goes here
       Detailed explanation goes here"""
    
    err = {}
    
    if 'R' in S_ref and 'R' in S:
        try:
            err['R'] = errorSO3(S['R'], S_ref['R'])
        except:
            warnings.warn('Angle Error is too high.')
    
    if 'v' in S_ref and 'v' in S:
        err['v'] = S['v'] - S_ref['v']
    
    if 'p' in S_ref and 'p' in S:
        err['p'] = S['p'] - S_ref['p']
    
    if 'w' in S_ref and 'w' in S:
        err['w'] = S['w'] - S_ref['w']
    
    if 'omega_dot' in S_ref and 'omega_dot' in S:
        err['omega_dot'] = S['omega_dot'] - S_ref['omega_dot']
    
    if 'v_dot' in S_ref and 'v_dot' in S:
        err['v_dot'] = S['v_dot'] - S_ref['v_dot']
    
    if 's' in S_ref and 's' in S:
        err['s'] = S['s'] - S_ref['s']
    
    if 'b_g' in S_ref and 'b_g' in S:
        err['b_g'] = S['b_g'] - S_ref['b_g']
    
    if 'b_s' in S_ref and 'b_s' in S:
        err['b_s'] = S['b_s'] - S_ref['b_s']
    
    if 'T_a' in S_ref and 'T_a' in S:
        err['T_a'] = S['T_a'] - S_ref['T_a']
    
    if 'b_omega_dot' in S_ref and 'b_omega_dot' in S:
        err['b_omega_dot'] = S['b_omega_dot'] - S_ref['b_omega_dot']
    
    return err


def estimate_T_and_b(y, u, Q):
    """ESTIMATE_T_AND_B Summary of this function goes here
       Detailed explanation goes here"""
    
    assert y.shape == u.shape
    
    A = np.zeros((12, 12))
    b = np.zeros((12, 1))
    
    for n in range(y.shape[1]):
        H_n = np.hstack([np.kron(u[:, n].T, np.eye(3)), np.eye(3)])
        Ht_Q_inv = H_n.T @ np.linalg.inv(Q)
        A = A + Ht_Q_inv @ H_n
        b = b + Ht_Q_inv @ y[:, n:n+1]
    
    Tb = np.linalg.solve(A, b)
    
    T = Tb[:9].reshape(3, 3)
    b = Tb[9:12]
    
    return T, b


def get_release_inds(time, release_times, IN_time, T):
    """Get release indices based on time parameters"""
    
    # IN_time = 5
    IN_samples = int(IN_time / T)
    inds_growth = np.zeros((IN_samples, len(release_times)), dtype=int)
    
    for i_t in range(len(release_times)):
        inds_t = np.where(time >= release_times[i_t])[0]
        inds_growth[:, i_t] = inds_t[:IN_samples]
    
    IN_time_array = np.arange(IN_samples) * T
    
    res = {}
    res['inds_growth'] = inds_growth
    res['IN_time_array'] = IN_time_array
    
    return res


def interpolate_pos_and_rotation(S, imu_time, rig_time):
    """INTERPOLATE_POS_AND_ROTATION Interpolate IMU pos and rotation estimates 
    to rig time"""
    
    S['p_rig_time'] = np.zeros((3, len(rig_time)))
    
    for i in range(3):
        interp_func = interp1d(imu_time, S['p'][i, :], kind='cubic', 
                              fill_value='extrapolate')
        S['p_rig_time'][i, :] = interp_func(rig_time)
    
    # Find the fractional indices using linear interpolation
    interp_func = interp1d(imu_time, np.arange(len(imu_time)), kind='linear', 
                          fill_value='extrapolate')
    inds_time = interp_func(rig_time)
    
    R_rig_time = np.zeros((3, 3, len(inds_time)))
    
    for n in range(len(inds_time)):
        frac = inds_time[n] - np.floor(inds_time[n])
        
        if frac > 1:
            warnings.warn("fraction larger than 1")
        elif frac == 0:
            # Same point in time
            R_rig_time[:, :, n] = S['R'][:, :, int(np.round(inds_time[n]))]
        else:
            # Calculate the rotation vector and scale it
            left_ind = int(np.floor(inds_time[n]))
            right_ind = left_ind + 1
            R_left = S['R'][:, :, left_ind]
            R_right = S['R'][:, :, right_ind]
            
            theta = logSO3(invSO3(R_left) @ R_right)
            R_rig_time[:, :, n] = R_left @ expSO3(frac * theta)
    
    S['R_rig_time'] = R_rig_time
    
    return S


def lsq_triad(y, Q):
    """LSQ_TRIAD Weighted Mean of triad
       y = Hu + e , e ~ N(0,Q)
       u = (H'*Q^{-1}*H)^{-1}(H'*Q^{-1}*y)
       Where u is triad"""
    
    assert y.shape[0] % 3 == 0
    assert Q.shape[0] == Q.shape[1]
    
    N = y.shape[0] // 3
    H = np.tile(np.eye(3), (N, 1))
    L = cholesky(Q, lower=True)
    
    t1 = H.T @ np.linalg.inv(L)
    t2 = np.linalg.solve(L, y)
    u = np.linalg.solve(t1 @ t1.T, t1 @ t2)
    Qu = np.linalg.inv(t1 @ t1.T)
    
    return u, Qu


def lsq_triad_naive(y, Q):
    """LSQ_TRIAD Summary of this function goes here
       y = Hu + e , e ~ N(0,Q)
       u = (H'*Q^{-1}*H)^{-1}(H'*Q^{-1}*y)
       Where u is triad"""
    
    assert y.shape[0] % 3 == 0
    assert Q.shape[0] == Q.shape[1]
    
    N = y.shape[0] // 3
    H = np.kron(np.ones((N, 1)), np.eye(3))
    HT_Q_inv = H.T @ np.linalg.inv(Q)
    
    u = np.linalg.solve(HT_Q_inv @ H, HT_Q_inv @ y)
    Qu = np.linalg.inv(HT_Q_inv @ H)
    
    return u, Qu


def rotate_measurements(S_in, R):
    """ROTATE_MEASUREMENTS Summary of this function goes here
       Detailed explanation goes here"""
    
    S_out = {}
    S_out['y'] = R @ S_in['y']
    S_out['Q'] = R @ S_in['Q'] @ R.T
    S_out['Q_inv'] = np.linalg.inv(S_out['Q'])
    
    return S_out


def run_filter(sensorData, initData, my_settings, S_ref, myFilter):
    """RUN_FILTER Run filter and calculate error 
       Detailed explanation goes here"""
    
    res = {}
    res['filt'], res['pred'] = myFilter(sensorData, initData, my_settings)
    
    err = {}
    
    if 'R' in S_ref:
        try:
            err['R'] = errorSO3(res['filt']['R'], S_ref['R'])
        except:
            warnings.warn('Angle Error is too high.')
    
    if 'v' in S_ref:
        err['v'] = res['filt']['v'] - S_ref['v']
    
    if 'p' in S_ref:
        err['p'] = res['filt']['p'] - S_ref['p']
    
    if 'w' in S_ref and 'w' in res['filt']:
        err['w'] = res['filt']['w'] - S_ref['w']
    
    if 'omega_dot' in S_ref:
        err['omega_dot'] = res['pred']['omega_dot'] - S_ref['omega_dot']
    
    if 'v_dot' in S_ref:
        err['v_dot'] = res['pred']['v_dot'] - S_ref['v_dot']
    
    res['err'] = err
    
    return res


def run_filter_w_error(sensorData, initData, my_settings, S_ref, myFilter):
    """RUN_FILTER Run filter and calculate error 
       Detailed explanation goes here"""
    
    res = {}
    res['filt'], res['pred'] = myFilter(sensorData, initData, my_settings)
    
    err = {}
    err['filt'] = compute_error(res['filt'], S_ref)
    err['pred'] = compute_error(res['pred'], S_ref)
    res['err'] = err
    
    return res 