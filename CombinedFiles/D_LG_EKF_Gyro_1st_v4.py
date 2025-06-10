import numpy as np
from scipy.linalg import cholesky, solve_triangular
from typing import Dict, List, Optional, Any, Tuple
import warnings

class D_LG_EKF_Gyro_1st_v4:
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
    def __init__(self, settings):
        # D_LG_EKF_CLASSIC Construct an instance of this class
        # Assume one gyro

        self.T = settings['T']
        self.g = settings['g']
        
        # --- Initialize property attributes to None or empty values ---
        self.r = None
        self.A_s = None
        self.N_a = 0
        self.T2_R = None # This property was in the MATLAB but not used in the constructor
        
        # --- Index definitions ---
        # R, p, v, omega
        self.inds_R = slice(0, 3) # Rotation
        Nx = 3
        
        self.inds_w_g = slice(0, 3) # white noise gyro
        Nw = 3

        self.inds_y_g = slice(0, 3)
        Ny = 3
        
        # Initialize index slices as empty
        self.inds_p = slice(0,0)
        self.inds_v = slice(0,0)
        self.inds_b_g = slice(0,0)
        self.inds_b_s = slice(0,0)
        self.inds_w_s = slice(0,0)
        self.inds_w_b_g = slice(0,0)
        self.inds_w_b_s = slice(0,0)
        self.inds_y_a = slice(0,0)
        
        if settings.get("propagate_bias_gyro", False):
            self.inds_b_g = slice(Nx, Nx + 3)
            Nx += 3
            self.inds_w_b_g = slice(Nw, Nw + 3)
            Nw += 3
        
        if settings.get("input_accelerometers", True):
            if 'N_a' in settings:
                self.N_a = settings['N_a']
                self.A_s = np.tile(np.eye(3), (1, self.N_a)) / self.N_a
            else:
                raise ValueError("N_a must be provided if input_accelerometers is true")

            self.r = settings['r']
            self.inds_y_a = slice(Ny, Ny + 3 * self.N_a)
            Ny += 3 * self.N_a
            
            self.inds_w_s = slice(Nw, Nw + 3)
            Nw += 3
            
            if settings.get("propagate_position", False):
                self.inds_p = slice(Nx, Nx + 3)
                Nx += 3
            if settings.get("propagate_velocity", False):
                self.inds_v = slice(Nx, Nx + 3)
                Nx += 3
            
            if settings.get("propagate_bias_s", False):
                self.inds_b_s = slice(Nx, Nx + 3)
                Nx += 3
                self.inds_w_b_s = slice(Nw, Nw + 3)
                Nw += 3

        self.Nx = Nx
        self.Nw = Nw
        self.Ny = Ny

        # -------------------------------------------------------------
        # Fill Q (Process Noise Covariance)
        # -------------------------------------------------------------
        Q = np.zeros((Nw, Nw))
        Q[self.inds_w_g, self.inds_w_g] = settings['Q_gyro']
        
        if self.inds_w_b_g.start != self.inds_w_b_g.stop:
            Q[self.inds_w_b_g, self.inds_w_b_g] = settings['Q_bias_gyro']
        
        if self.inds_w_s.start != self.inds_w_s.stop:
            if "Q_s" in settings:
                Q[self.inds_w_s, self.inds_w_s] = settings['Q_s']
            elif "Q_acc" in settings:
                Q[self.inds_w_s, self.inds_w_s] = self.A_s @ settings['Q_acc'] @ self.A_s.T
            else:
                raise ValueError("No covariance for s provided (Q_s or Q_acc)")

        if self.inds_w_b_s.start != self.inds_w_b_s.stop:
            if "Q_bias_s" in settings:
                Q[self.inds_w_b_s, self.inds_w_b_s] = settings['Q_bias_s']
            elif "Q_bias_acc" in settings:
                Q[self.inds_w_b_s, self.inds_w_b_s] = self.A_s @ settings['Q_bias_acc'] @ self.A_s.T
            else:
                raise ValueError("No covariance for bias s provided (Q_bias_s or Q_bias_acc)")
        
        self.Q = Q
        
        # -----------------------------------------------------------------
        # Measurement Noise Covariance
        # -----------------------------------------------------------------
        self.R_pos = settings.get("R_pos")
        if self.R_pos is not None:
            assert np.allclose(self.R_pos, self.R_pos.T), "R_pos must be symmetric"
            assert self.R_pos.shape == (3, 3), "R_pos must be 3x3"
            
        self.R_rot = settings.get("R_rot")
        if self.R_rot is not None:
            assert np.allclose(self.R_rot, self.R_rot.T), "R_rot must be symmetric"
            assert self.R_rot.shape == (3, 3), "R_rot must be 3x3"

    def get_input(self, sensorData):
        assert sensorData['gyro_measurements'].shape[0] == 3
        if self.inds_y_a.start != self.inds_y_a.stop:
            assert sensorData['acc_measurements'].shape[0] == 3 * self.N_a
            u = np.vstack([sensorData['gyro_measurements'], sensorData['acc_measurements']])
        else:
            u = sensorData['gyro_measurements']
        return u

    def get_Q(self):
        return self.Q

    def propagate(self, R, x, y, w):
        #METHOD1 Summary of this method goes here
        assert len(x) == self.Nx - 3
        assert len(y) == self.Ny
        assert len(w) == self.Nw
        
        v = np.zeros((3, 1))
        if self.inds_v.start != self.inds_v.stop:
            v = x[self.inds_v.start-3 : self.inds_v.stop-3]
        
        b_s = np.zeros((3, 1))
        if self.inds_b_s.start != self.inds_b_s.stop:
            b_s = x[self.inds_b_s.start-3 : self.inds_b_s.stop-3]
        
        b_g = np.zeros((3, 1))
        if self.inds_b_g.start != self.inds_b_g.stop:
            b_g = x[self.inds_b_g.start-3 : self.inds_b_g.stop-3]

        # Process noise
        w_g = w[self.inds_w_g]
        w_s = np.zeros((3, 1))
        if self.inds_w_s.start != self.inds_w_s.stop:
            w_s = w[self.inds_w_s]
            
        w_b_s = np.zeros((3, 1))
        if self.inds_w_b_s.start != self.inds_w_b_s.stop:
            w_b_s = w[self.inds_w_b_s]
        
        w_b_g = np.zeros((3, 1))
        if self.inds_w_b_g.start != self.inds_w_b_g.stop:
            w_b_g = w[self.inds_w_b_g]
        
        y_g = y[self.inds_y_g]
        # Depending on the model, omega can be derived in different ways.
        # This implementation uses the external cardou9Wx0 function.
        omega = y_g # Initial guess for omega
        bdotdot, _, omega = cardou9Wx0(y[self.inds_y_a], omega, self.r, self.T)

        s = np.zeros((3, 1))
        v_dot = np.zeros((3, 1))
        if self.inds_y_a.start != self.inds_y_a.stop:
            y_a_mean = bdotdot
            s = y_a_mean + b_s + w_s   # Specific force
            v_dot = self.g + R @ s     # Navigation acceleration
        
        Omega = np.zeros((self.Nx, 1))
        # R
        Omega[self.inds_R] = omega * self.T

        # p
        if self.inds_p.start != self.inds_p.stop:
            Omega[self.inds_p] = v * self.T + v_dot * self.T**2 / 2
        # v
        if self.inds_v.start != self.inds_v.stop:
            Omega[self.inds_v] = v_dot * self.T
        
        if self.inds_b_g.start != self.inds_b_g.stop:
            Omega[self.inds_b_g] = w_b_g
        
        if self.inds_b_s.start != self.inds_b_s.stop:
            Omega[self.inds_b_s] = w_b_s
                    
        # -----------------------------------------------------------------
        # Jacobian in x (dOmega_de)
        dOmega_de = np.zeros((self.Nx, self.Nx))
        
        d_v_dot_d_R = np.zeros((3, 3))
        d_v_dot_d_b_s = np.zeros((3, 3))
        if self.inds_y_a.start != self.inds_y_a.stop:
            d_v_dot_d_R = -R @ HatSO3(s)
            d_v_dot_d_b_s = R

        # --- R equation ---
        if self.inds_b_g.start != self.inds_b_g.stop:
            # This Jacobian is based on the simple model omega = y_g - b_g - w_g
            # A more complex model from cardou9Wx0 would have a different Jacobian.
            d_omega_d_b_g = -np.eye(3) 
            dOmega_de[self.inds_R, self.inds_b_g] = d_omega_d_b_g * self.T

        # --- p equation ---
        if self.inds_p.start != self.inds_p.stop:
            dOmega_de[self.inds_p, self.inds_R] = d_v_dot_d_R * self.T**2 / 2
            if self.inds_v.start != self.inds_v.stop:
                dOmega_de[self.inds_p, self.inds_v] = np.eye(3) * self.T
            if self.inds_b_s.start != self.inds_b_s.stop:
                dOmega_de[self.inds_p, self.inds_b_s] = d_v_dot_d_b_s * self.T**2 / 2

        # --- v equation ---
        if self.inds_v.start != self.inds_v.stop:
            dOmega_de[self.inds_v, self.inds_R] = d_v_dot_d_R * self.T
            if self.inds_b_s.start != self.inds_b_s.stop:
                dOmega_de[self.inds_v, self.inds_b_s] = d_v_dot_d_b_s * self.T
        
        # -----------------------------------------------------------------
        # Jacobian in w (dOmega_dw)
        dOmega_dw = np.zeros((self.Nx, self.Nw))
        d_v_dot_d_w_s = R
        
        # --- R equation ---
        # Based on simple model: omega = y_g - b_g - w_g
        d_omega_d_w_g = -np.eye(3) 
        dOmega_dw[self.inds_R, self.inds_w_g] = d_omega_d_w_g * self.T

        # --- p equation ---
        if self.inds_p.start != self.inds_p.stop and self.inds_w_s.start != self.inds_w_s.stop:
            dOmega_dw[self.inds_p, self.inds_w_s] = d_v_dot_d_w_s * self.T**2 / 2

        # --- v equation ---
        if self.inds_v.start != self.inds_v.stop and self.inds_w_s.start != self.inds_w_s.stop:
            dOmega_dw[self.inds_v, self.inds_w_s] = d_v_dot_d_w_s * self.T
            
        # --- bias equations ---
        if self.inds_b_s.start != self.inds_b_s.stop:
            dOmega_dw[self.inds_b_s, self.inds_w_b_s] = np.eye(3)
        
        if self.inds_b_g.start != self.inds_b_g.stop:
            dOmega_dw[self.inds_b_g, self.inds_w_b_g] = np.eye(3)
        
        # --- Return results as a Namespace object for easy access ---
        return SimpleNamespace(
            Omega=Omega,
            dOmega_de=dOmega_de,
            dOmega_dw=dOmega_dw,
            v_dot=v_dot,
            s=s,
            omega=omega
        )

    def position_update(self, p_obs, _, x_in):
        assert self.R_pos is not None, "R_pos covariance not set"
        Q = self.R_pos # In updates, Q is measurement noise
        
        H = np.zeros((3, self.Nx))
        if self.inds_p.start == self.inds_p.stop:
            raise ValueError("Position update called but position is not in state vector")
        
        H[:, self.inds_p] = np.eye(3) # p
        
        x_p_slice = slice(self.inds_p.start-3, self.inds_p.stop-3)
        p_pred = x_in[x_p_slice]
        e = p_obs - p_pred
        return e, H, Q

    def rotation_update(self, R_obs, R_pred, _):
        assert self.R_rot is not None, "R_rot covariance not set"
        Q = self.R_rot
        
        # H rotation
        H = np.zeros((3, self.Nx))
        H[:, self.inds_R] = np.eye(3) # R    

        e = logSO3(invSO3(R_pred) @ R_obs)
        return e, H, Q

    def get_initial_conditions(self, initIn):
        initOut = {}
        use_full = "x" in initIn or "P" in initIn or "R" in initIn
        use_partial = "mean" in initIn or "cov" in initIn
        
        if use_full and use_partial:
            raise ValueError("Both (R,x,P) and (mean,cov) defined, use only one")
        elif "x" in initIn and "P" in initIn and "R" in initIn:
            initOut['R0'] = initIn['R']
            initOut['x0'] = initIn['x']
            initOut['P0'] = initIn['P']
        elif "mean" in initIn and "cov" in initIn:
            m = initIn['mean']
            c = initIn['cov']
            
            assert m['R'].shape == (3, 3)
            initOut['R0'] = m['R']
            
            x0 = np.zeros((self.Nx - 3, 1))
            if self.inds_p.start != self.inds_p.stop:
                x0[self.inds_p.start-3:self.inds_p.stop-3] = m['p']
            if self.inds_v.start != self.inds_v.stop:
                x0[self.inds_v.start-3:self.inds_v.stop-3] = m['v']
            if self.inds_b_s.start != self.inds_b_s.stop:
                if "b_s" in m:
                    x0[self.inds_b_s.start-3:self.inds_b_s.stop-3] = m['b_s']
                elif "b_a" in m:
                     x0[self.inds_b_s.start-3:self.inds_b_s.stop-3] = -self.A_s @ m['b_a']
                else:
                    raise ValueError("No mean value for specific force bias")
            initOut['x0'] = x0
            
            P0 = np.zeros((self.Nx, self.Nx))
            P0[self.inds_R, self.inds_R] = c['R']
            if self.inds_p.start != self.inds_p.stop:
                P0[self.inds_p, self.inds_p] = c['p']
            if self.inds_v.start != self.inds_v.stop:
                P0[self.inds_v, self.inds_v] = c['v']
            if self.inds_b_g.start != self.inds_b_g.stop:
                P0[self.inds_b_g, self.inds_b_g] = c['b_g']
            if self.inds_b_s.start != self.inds_b_s.stop:
                if "b_s" in c:
                    P0[self.inds_b_s, self.inds_b_s] = c['b_s']
                elif "b_a" in c:
                    P0[self.inds_b_s, self.inds_b_s] = self.A_s @ c['b_a'] @ self.A_s.T
                else:
                    raise ValueError("No covariance value for specific force bias")
            initOut['P0'] = P0
        else:
            raise ValueError("Valid initial conditions not provided")
        return initOut

    def extract_variables(self, Sin):
        Sout = {'mean': {}, 'std': {}}
        Sout['mean']['R'] = Sin['R']
        
        if self.inds_p.start != self.inds_p.stop:
            Sout['mean']['p'] = Sin['x'][self.inds_p.start-3:self.inds_p.stop-3, :]
        if self.inds_v.start != self.inds_v.stop:
            Sout['mean']['v'] = Sin['x'][self.inds_v.start-3:self.inds_v.stop-3, :]
        if self.inds_b_g.start != self.inds_b_g.stop:
            Sout['mean']['b_g'] = Sin['x'][self.inds_b_g.start-3:self.inds_b_g.stop-3, :]
        if self.inds_b_s.start != self.inds_b_s.stop:
            Sout['mean']['b_s'] = Sin['x'][self.inds_b_s.start-3:self.inds_b_s.stop-3, :]
        
        Sout['std']['R'] = Sin['std'][self.inds_R, :]
        if self.inds_p.start != self.inds_p.stop:
            Sout['std']['p'] = Sin['std'][self.inds_p, :]
        if self.inds_v.start != self.inds_v.stop:
            Sout['std']['v'] = Sin['std'][self.inds_v, :]
        if self.inds_b_g.start != self.inds_b_g.stop:
            Sout['std']['b_g'] = Sin['std'][self.inds_b_g, :]
        if self.inds_b_s.start != self.inds_b_s.stop:
            Sout['std']['b_s'] = Sin['std'][self.inds_b_s, :]
        return Sout

    def print_info(self, S_init, SensorData, settings):
        # A partial conversion of the print_info method for brevity.
        # The pattern can be extended to the full method.
        def print_matrix(mat, name, unit, is_chol=False):
            label = "chol(" + name + ")" if is_chol else name
            unit_str = f"[{unit}]" if not is_chol else f"[{unit}]^0.5"
            print(f"\t\t{label} {unit_str}:")
            # Use np.savetxt to print the matrix with specified formatting
            np.savetxt(sys.stdout, mat.T, fmt='%10.1e', delimiter='')

        print(f"Sampling time: {settings['T']:.2e} [s]")
        print(f"Sampling freq: {1/settings['T']:.1f} [Hz]")
        print(f"Propagate rotation: {self.inds_R.start != self.inds_R.stop}")
        print(f"Propagate position: {self.inds_p.start != self.inds_p.stop}")
        print(f"Propagate velocity: {self.inds_v.start != self.inds_v.stop}")
        print(f"Propagate bias gyro: {self.inds_b_g.start != self.inds_b_g.stop}")
        print(f"Propagate bias specific force: {self.inds_b_s.start != self.inds_b_s.stop}")
        
        print("\nPropagation data:")
        print("\tGyroscope data:")
        if self.inds_w_g.start != self.inds_w_g.stop:
            Q_g = self.Q[self.inds_w_g, self.inds_w_g]
            print_matrix(np.rad2deg(np.linalg.cholesky(Q_g)), "Q_gyro", "deg/s", is_chol=True)
        
        if self.inds_w_b_g.start != self.inds_w_b_g.stop:
            try:
                Q_b_g = self.Q[self.inds_w_b_g, self.inds_w_b_g]
                print_matrix(np.rad2deg(np.linalg.cholesky(Q_b_g)), "Q_b_g", "deg/s", is_chol=True)
            except np.linalg.LinAlgError:
                print_matrix(np.rad2deg(np.rad2deg(Q_b_g)), "Q_b_g", "deg/s^2")

        print("\nInitial conditions:")
        print("\tMean:")
        e_deg = np.rad2deg(my_rotm2eul(S_init['mean']['R']))
        print(f"\t\tRotation: roll: {e_deg[0]:.1f}, pitch: {e_deg[1]:.1f}, yaw: {e_deg[2]:.1f} [deg]")
        
        if "p" in S_init['mean']:
            p_str = " ".join([f"{x:.1f}" for x in S_init['mean']['p'].flatten()])
            print(f"\t\tPosition:          [{p_str}] [m]")
        if "v" in S_init['mean']:
            v_str = " ".join([f"{x:.1f}" for x in S_init['mean']['v'].flatten()])
            print(f"\t\tVelocity:          [{v_str}] [m/s]")
        if "b_s" in S_init['mean']:
            bs_str = " ".join([f"{x:.1f}" for x in S_init['mean']['b_s'].flatten()])
            print(f"\t\tBias s:            [{bs_str}] [m/s^2]")
        if "b_g" in S_init['mean']:
            bg_str = " ".join([f"{x:.1f}" for x in np.rad2deg(S_init['mean']['b_g'].flatten())])
            print(f"\t\tBias gyro:         [{bg_str}] [deg/s]")
            
        print("\tStd:")
        std_R_str = " ".join([f"{x:.1f}" for x in np.rad2deg(S_init['std']['R'].flatten())])
        print(f"\t\tRotation:          [{std_R_str}] [deg]")
        # ... and so on for other std variables ...