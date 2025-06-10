import numpy as np
from scipy.linalg import cholesky, solve_triangular
from typing import Dict, List, Optional, Any, Tuple
import warnings

def DLGEKFv4(sensor_data: Dict, init: Dict, model: Any, settings: Dict) -> Dict:
    """
    D-LG-EKF-ARRAY - Discrete Lie-Group Extended Kalman Filter for Inertial Navigation
    
    State vector:
    R           - Rotation between body frame and navigation frame
    x[0:3]      - angular velocity in body frame
    x[3:6]      - position in navigation frame
    x[6:9]      - velocity in navigation frame
    
    Covariance matrix:
    P[0:3]      - Rotation 
    P[3:6]      - angular velocity in body frame
    P[6:9]      - position in navigation frame
    P[9:12]     - velocity in navigation frame
    
    Args:
        sensor_data: Dictionary containing sensor measurements
        init: Dictionary containing initial conditions
        model: Model object with required methods
        settings: Dictionary containing algorithm settings
        
    Returns:
        Dictionary containing filtered results
    """
    
    # Get verbose setting
    verbose = settings.get('verbose', False)
    
    if verbose:
        print("Run Discrete Lie-Group EKF Inertial Navigation v4")
    
    # -------------------------------------------------------------------------
    # These must exist
    u_prop = model.get_input(sensor_data)  # acceleration measurements history data in body frame
    Q_prop = model.get_Q()
    Nt = u_prop.shape[1]

    if verbose:
        print(f"Number of time-sample: {Nt}")
        print(f"Size of input vector: {u_prop.shape[0]}")
     
    # -------------------------------------------------------------------------
    # Kalman updates using gyroscopes
    do_gyro_updates = settings.get('do_gyro_updates', True)
    
    if 'gyro_measurements' in sensor_data and do_gyro_updates:
        u_g = sensor_data['gyro_measurements']
        assert u_g.shape == (3, Nt), f"Expected gyro shape (3, {Nt}), got {u_g.shape}"
    else:
        u_g = np.full((3, Nt), np.nan)
        
    if verbose:
        print(f"Number of gyroscope updates: {np.sum(~np.isnan(u_g[0, :]))}")
    
    # -------------------------------------------------------------------------
    # Position updates
    do_position_updates = settings.get('do_position_updates', True)
    
    if 'position_measurements' in sensor_data and do_position_updates:
        p_obs = sensor_data['position_measurements']
        assert p_obs.shape[0] == 3, f"Expected position shape (3, N), got {p_obs.shape}"
        assert p_obs.shape[1] == Nt, f"Expected position shape (3, {Nt}), got {p_obs.shape}"
    else:
        p_obs = np.full((3, Nt), np.nan)
        
    if verbose:
        print(f"Number of position updates: {np.sum(~np.isnan(p_obs[0, :]))}")
    
    # -------------------------------------------------------------------------
    # Rotation updates
    do_rotation_updates = settings.get('do_rotation_updates', True)
    
    if 'rotation_measurements' in sensor_data and do_rotation_updates:
        R_obs = sensor_data['rotation_measurements']
        assert R_obs.shape[0] == 3, f"Expected rotation shape (3, 3, N), got {R_obs.shape}"
        assert R_obs.shape[2] == Nt, f"Expected rotation shape (3, 3, {Nt}), got {R_obs.shape}"
    else:
        R_obs = np.full((3, 3, Nt), np.nan)
        
    if verbose:
        print(f"Number of rotation updates: {np.sum(~np.isnan(R_obs[0, 0, :]))}")
    
    # -------------------------------------------------------------------------
    # Zero velocity updates
    if 'zero_velocity_updates' in sensor_data:
        zupts = sensor_data['zero_velocity_updates']
    else:
        zupts = np.full((3, Nt), np.nan)
        
    if verbose:
        print(f"Number of ZUPTs: {np.sum(~np.isnan(zupts[0, :]))}")
    
    # -------------------------------------------------------------------------
    # Settings for what to save
    save_full_covariances = settings.get('save_full_covariances', False)
    save_pred = settings.get('save_pred', False)
    save_aux_vars = settings.get('save_aux_vars', False)
    save_jacobians = settings.get('save_jacobians', False)
    
    if verbose:
        print(f"Save full covariances: {save_full_covariances}")
        print(f"Save prediction: {save_pred}")
        print(f"Save auxiliary variables: {save_aux_vars}")
        print(f"Save jacobians: {save_jacobians}")
    
    # -------------------------------------------------------------------------
    Nx = model.Nx
    Nw = model.Nw
    
    if verbose:
        print(f"Size Nx: {Nx}")
        print(f"Size Nw: {Nw}")
    
    # Initial aposteriori state
    init_out = model.get_initial_conditions(init)
    R0 = init_out['R0']
    x0 = init_out['x0']
    P0 = init_out['P0']
    
    assert x0.shape[0] == Nx - 3, f"Expected x0 shape ({Nx-3},), got {x0.shape}"
    assert P0.shape == (Nx, Nx), f"Expected P0 shape ({Nx}, {Nx}), got {P0.shape}"
    
    # Initialize arrays
    # n|n-1 (prediction)
    R_pred = np.full((3, 3, Nt), np.nan)
    x_pred = np.full((Nx-3, Nt), np.nan)  # Reduce by 3 for rotation
    if save_full_covariances:
        P_pred = np.full((Nx, Nx, Nt), np.nan)
    std_pred_diag = np.full((Nx, Nt), np.nan)

    # n|n (filtered)
    R_filt = np.zeros((3, 3, Nt))
    x_filt = np.zeros((Nx-3, Nt))  # Reduce by 3 for rotation
    if save_full_covariances:
        P_filt = np.zeros((Nx, Nx, Nt))
    std_filt_diag = np.zeros((Nx, Nt))
    
    # Accelerations    
    omega_dot = np.full((3, Nt), np.nan)
    v_dot = np.full((3, Nt), np.nan)
    s = np.full((3, Nt), np.nan)
    
    # Set initial values (convert to 0-based indexing)
    R_filt[:, :, 0] = R0
    x_filt[:, 0] = x0
    if save_full_covariances:
        P_filt[:, :, 0] = P0
    std_filt_diag[:, 0] = np.sqrt(np.diag(P0))
    
    # Basic checks
    assert not _is_any_nan(R0), "R0 contains NaN values"
    assert not _is_any_nan(x0), "x0 contains NaN values"
    assert not _is_any_nan(P0), "P0 contains NaN values"
    
    # Current state (mise à jour)
    R_maj = R0.copy()
    x_maj = x0.copy()
    P_maj = P0.copy()
    
    # Variables storage 
    logLL = np.full(Nt, np.nan)
    residuals_normalized = [None] * Nt
    
    if save_jacobians:
        K_tot = [None] * Nt
        F_tot = [None] * Nt
        G_tot = [None] * Nt
        F_pre_tot = [None] * Nt
        G_pre_tot = [None] * Nt
    
    if verbose:
        S0 = {
            'R': init_out['R0'],
            'x': init_out['x0'],
            'std': np.sqrt(np.diag(init_out['P0']))
        }
        S_init = model.extract_variables(S0)
        model.print_info(S_init, sensor_data, settings)
    
    # Print initial measurements if verbose
    if verbose and np.any(~np.isnan(u_g)):
        print("Do Kalman updates using gyroscopes")
        print("\tStart of gyroscope triad measurements as updates:")
        for i in range(min(3, u_g.shape[1])):
            print(f"\t{np.rad2deg(u_g[0, i]):10.1f} {np.rad2deg(u_g[1, i]):10.1f} {np.rad2deg(u_g[2, i]):10.1f}")

    if verbose and np.any(~np.isnan(p_obs)):
        print("Do Kalman updates using positions")
        print("\tStart of position updates:")
        for i in range(min(3, p_obs.shape[1])):
            print(f"\t{p_obs[0, i]:10.1f} {p_obs[1, i]:10.1f} {p_obs[2, i]:10.1f}")
            
    if verbose and np.any(~np.isnan(R_obs)):
        print("Do Kalman updates using rotations")
        print(f"\t{R_obs[:, :, 0]}")
    
    if verbose and np.any(~np.isnan(zupts)):
        print("Do Kalman updates using ZUPTs")
    
    w = np.zeros(Nw)  # Set noise to zero
    res_prop = {}
    
    # Main filter loop
    for n in range(1, Nt):  # Convert to 0-based indexing (start from 1, not 2)
        if 'progress' in settings and callable(settings['progress']):
            settings['progress'](n+1, Nt)  # Convert back to 1-based for progress callback
        
        # Propagation
        # Get n-1 values (0-based indexing)
        u_n = u_prop[:, n-1]   # acceleration measurements data in body frame
        
        if n == (Nt - 11):  # Adjust debug condition for 0-based indexing
            pass  # Debug point
                
        # Update angular velocity if available from previous propagation
        if 'omega' in res_prop:
            x_maj[0:3] = res_prop['omega']
            
        res_prop = model.propagate(R_maj, x_maj, u_n, w)
        Omega_n = res_prop['Omega']
        dOmega_de_n = res_prop['dOmega_de']
        dOmega_dw_n = res_prop['dOmega_dw']
        
        if save_jacobians:
            F_pre_tot[n] = dOmega_de_n
            G_pre_tot[n] = dOmega_dw_n
        
        # Store accelerations at n-1
        if 'omega_dot' in res_prop:
            omega_dot[:, n-1] = res_prop['omega_dot']
        if 'v_dot' in res_prop:
            v_dot[:, n-1] = res_prop['v_dot']
        if 's' in res_prop:
            s[:, n-1] = res_prop['s']
        
        # Mean propagation (n|n-1)
        R_prop = R_maj @ _exp_so3(Omega_n[0:3])  # SO(3)  
        x_prop = x_maj + Omega_n[3:]  # Euclidean
        
        # Covariance propagation
        Ad = np.eye(P_maj.shape[0])  # adjoint
        Ad[0:3, 0:3] = _ad_so3(_exp_so3(-Omega_n[0:3]))
        Phi_n = np.eye(P_maj.shape[0])
        Phi_n[0:3, 0:3] = _phi_so3(-Omega_n[0:3])
        
        F_n = Ad + Phi_n @ dOmega_de_n
        G_n = Phi_n @ dOmega_dw_n
        P_prop = F_n @ P_maj @ F_n.T + G_n @ Q_prop @ G_n.T
        
        # Store jacobians
        if save_jacobians:
            F_tot[n] = F_n
            G_tot[n] = G_n
        
        # Basic checks
        assert not _is_any_nan(R_prop), f"R_prop contains NaN at step {n}"
        assert not _is_any_nan(x_prop), f"x_prop contains NaN at step {n}"
        assert not _is_any_nan(P_prop), f"P_prop contains NaN at step {n}"
        
        # Measurement update (n|n)
        Hs = []
        es = []
        Rs = []
        
        # Gyroscope update
        if not np.isnan(u_g[0, n]):
            e, H, R = model.gyroscope_update(u_g[:, n], R_prop, x_prop)
            es.append(e)
            Hs.append(H)
            Rs.append(R)
        
        # Position update
        if not np.isnan(p_obs[0, n]):
            e, H, R = model.position_update(p_obs[:, n], R_prop, x_prop)
            es.append(e)
            Hs.append(H)
            Rs.append(R)
        
        # Rotation update
        if not np.isnan(R_obs[0, 0, n]):
            e, H, R = model.rotation_update(R_obs[:, :, n], R_prop, x_prop)
            es.append(e)
            Hs.append(H)
            Rs.append(R)
        
        # Zero velocity update
        if not np.isnan(zupts[0, n]):
            e, H, R = model.zero_velocity_update(zupts[:, n], R_prop, x_prop)
            es.append(e)
            Hs.append(H)
            Rs.append(R)
        
        if len(es) > 0:
            # Combine all measurements
            H_update = np.vstack(Hs)
            e_update = np.concatenate(es)
            R_update = _block_diag(Rs)
            
            # Innovations and likelihood
            S = H_update @ P_prop @ H_update.T + R_update
            try:
                L = cholesky(S, lower=True)
            except np.linalg.LinAlgError as e:
                print(f"Cholesky decomposition failed at step {n}")
                raise e
            
            e_norm = solve_triangular(L, e_update, lower=True)
            logdetS = 2 * np.sum(np.log(np.diag(L)))
            logLL[n] = np.dot(e_norm, e_norm) + logdetS
            residuals_normalized[n] = e_norm
            
            # Kalman update
            A = solve_triangular(L.T, (P_prop @ H_update.T).T, lower=False).T
            K = solve_triangular(L, A.T, lower=True).T
            m = K @ e_update
            
            # Update mean
            R_maj = R_prop @ _exp_so3(m[0:3])  # SO(3)
            x_maj = x_prop + m[3:]
            
            # Update covariance
            Phi_maj = np.eye(P_prop.shape[0])
            Phi_maj[0:3, 0:3] = _phi_so3(-m[0:3])
            I_KH = np.eye(Nx) - K @ H_update
            P_maj = Phi_maj @ (I_KH @ P_prop @ I_KH.T + K @ R_update @ K.T) @ Phi_maj.T
            
            if save_jacobians:
                K_tot[n] = K
        else:
            # No Kalman updates
            R_maj = R_prop.copy()
            x_maj = x_prop.copy()
            P_maj = P_prop.copy()
        
        # Final checks
        assert not _is_any_nan(R_maj), f"R_maj contains NaN at step {n}"
        assert not _is_any_nan(x_maj), f"x_maj contains NaN at step {n}"
        assert not _is_any_nan(P_maj), f"P_maj contains NaN at step {n}"
        
        # Store prediction values (n|n-1)
        R_pred[:, :, n] = R_prop
        x_pred[:, n] = x_prop
        if save_full_covariances:
            P_pred[:, :, n] = P_prop
        std_pred_diag[:, n] = np.sqrt(np.diag(P_prop))
        
        # Store filtered values (n|n)
        R_filt[:, :, n] = R_maj
        x_filt[:, n] = x_maj
        if save_full_covariances:
            P_filt[:, :, n] = P_maj
        std_filt_diag[:, n] = np.sqrt(np.diag(P_maj))
    
    # Extract variables
    Sfilt = {
        'R': R_filt,
        'x': x_filt,
        'std': std_filt_diag
    }
    if save_full_covariances:
        Sfilt['P'] = P_filt
    if save_jacobians:
        Sfilt['K'] = K_tot
    
    Sprop = {
        'R': R_pred,
        'x': x_pred,
        'std': std_pred_diag
    }
    if save_full_covariances:
        Sprop['P'] = P_pred
    if save_jacobians:
        Sprop['F'] = F_tot
        Sprop['G'] = G_tot
        Sprop['F_pre'] = F_pre_tot
        Sprop['G_pre'] = G_pre_tot
    
    # Build result dictionary
    res = {
        'filt': model.extract_variables(Sfilt)
    }
    
    if save_pred:
        res['pred'] = model.extract_variables(Sprop)
    
    if save_aux_vars:
        aux = {
            'omega_dot': omega_dot,
            'v_dot': v_dot,
            's': s
        }
        res['aux'] = aux
    
    # Log likelihood
    valid_logLL = logLL[~np.isnan(logLL)]
    res['logL'] = {
        'value': -0.5 * np.sum(valid_logLL),
        'parts': logLL,
        'residuals_normalized': residuals_normalized
    }
    
    if save_full_covariances:
        res['tot'] = {
            'filt': Sfilt,
            'pred': Sprop
        }
    
    if 'label' in settings:
        res['label'] = settings['label']
    
    return res


# Helper functions (these would need to be implemented based on your specific Lie group operations)
def _is_any_nan(arr: np.ndarray) -> bool:
    """Check if array contains any NaN values"""
    return np.any(np.isnan(arr))


def _exp_so3(omega: np.ndarray) -> np.ndarray:
    """
    Exponential map for SO(3) - converts angular velocity vector to rotation matrix
    This is a placeholder - implement based on your specific requirements
    """
    # Rodrigues' rotation formula implementation
    angle = np.linalg.norm(omega)
    if angle < 1e-8:
        return np.eye(3) + _skew_symmetric(omega)
    
    axis = omega / angle
    K = _skew_symmetric(axis)
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)


def _skew_symmetric(v: np.ndarray) -> np.ndarray:
    """Create skew-symmetric matrix from 3D vector"""
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])


def _ad_so3(R: np.ndarray) -> np.ndarray:
    """
    Adjoint representation for SO(3)
    This is a placeholder - implement based on your specific requirements
    """
    return R  # Simplified - actual implementation depends on your Lie group library


def _phi_so3(omega: np.ndarray) -> np.ndarray:
    """
    Left Jacobian for SO(3)
    This is a placeholder - implement based on your specific requirements
    """
    angle = np.linalg.norm(omega)
    if angle < 1e-8:
        return np.eye(3)
    
    axis = omega / angle
    K = _skew_symmetric(axis)
    return (np.sin(angle) / angle) * np.eye(3) + ((1 - np.cos(angle)) / angle) * K


def _block_diag(matrices: List[np.ndarray]) -> np.ndarray:
    """Create block diagonal matrix from list of matrices"""
    if not matrices:
        return np.array([])
    
    # Calculate total size
    sizes = [m.shape[0] for m in matrices]
    total_size = sum(sizes)
    
    # Create block diagonal matrix
    result = np.zeros((total_size, total_size))
    start_idx = 0
    
    for matrix in matrices:
        size = matrix.shape[0]
        result[start_idx:start_idx + size, start_idx:start_idx + size] = matrix
        start_idx += size
    
    return result 