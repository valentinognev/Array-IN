import numpy as np
from scipy.linalg import cholesky
from typing import Dict, Optional, Tuple, Union

class D_LG_EKF_Gyro_2nd_v4:
    """
    D_LG_EKF_CLASSIC Summary of this class goes here
    Detailed explanation goes here
    alpha = [omega_dot; s]
    """
    
    def __init__(self, settings: Dict):
        """
        D_LG_EKF_CLASSIC Construct an instance of this class
        Detailed explanation goes here
        Assume one gyro
        """
        self.T = settings['T']
        self.g = settings['g']
        
        # Initialize indices
        self.inds_R = np.arange(3)  # Rotation matrix R_nb, from body to navigation frame
        self.inds_p = None  # Position in navigation frame
        self.inds_v = None  # Velocity in navigation frame
        self.inds_b_g = None  # Bias in gyroscope
        self.inds_b_s = None  # Bias in specific force
        self.inds_w_g = np.arange(3)  # white noise gyro
        self.inds_w_alpha = None  # alpha white noise
        self.inds_w_b_g = None  # Process noise bias gyroscope
        self.inds_w_b_s = None  # Process noise bias specific force
        self.inds_y_g = np.arange(3)
        self.inds_y_a = None
        
        if 'set_T2_R_zero' in settings and settings['set_T2_R_zero']:
            self.T2_R = 0
        else:
            self.T2_R = self.T**2

        # R, p, v, omega
        Nx = 3  # Rotation
        Nw = 3
        Ny = 3
        
        if 'propagate_bias_gyro' in settings and settings['propagate_bias_gyro']:
            self.inds_b_g = np.arange(3) + Nx
            Nx += 3
            self.inds_w_b_g = np.arange(3) + Nw
            Nw += 3
            
        if 'input_accelerometers' not in settings or settings['input_accelerometers']:
            self.N_a = settings['r'].shape[1]
            self.r = settings['r']
            self.A = self.compute_A_non_center(self.r)
            self.A_omega_dot = self.A[0:3, :]
            self.A_s = self.A[3:6, :]
            
            self.inds_y_a = np.arange(3 * self.N_a) + Ny
            Ny += 3 * self.N_a
            
            self.inds_w_alpha = np.arange(6) + Nw
            Nw += 6
            
            if 'propagate_position' in settings and settings['propagate_position']:
                self.inds_p = np.arange(3) + Nx
                Nx += 3
                
            if 'propagate_velocity' in settings and settings['propagate_velocity']:
                self.inds_v = np.arange(3) + Nx
                Nx += 3
                
            if 'propagate_bias_s' in settings and settings['propagate_bias_s']:
                self.inds_b_s = np.arange(3) + Nx
                Nx += 3
                self.inds_w_b_s = np.arange(3) + Nw
                Nw += 3

        self.Nx = Nx
        self.Nw = Nw
        self.Ny = Ny

        # Fill Q
        Q = np.zeros((Nw, Nw))
        Q[np.ix_(self.inds_w_g, self.inds_w_g)] = settings['Q_gyro']
        
        if self.inds_w_b_g is not None:
            Q[np.ix_(self.inds_w_b_g, self.inds_w_b_g)] = settings['Q_bias_gyro']
            
        if self.inds_w_alpha is not None:
            if 'Q_alpha' in settings:
                Q[np.ix_(self.inds_w_alpha, self.inds_w_alpha)] = settings['Q_alpha']
            elif 'Q_acc' in settings:
                Q[np.ix_(self.inds_w_alpha, self.inds_w_alpha)] = self.A @ settings['Q_acc'] @ self.A.T
            else:
                raise ValueError("No covariance for alpha")

        if self.inds_w_b_s is not None:
            if 'Q_bias_s' in settings:
                Q[np.ix_(self.inds_w_b_s, self.inds_w_b_s)] = settings['Q_bias_s']
            elif 'Q_bias_acc' in settings:
                Q[np.ix_(self.inds_w_b_s, self.inds_w_b_s)] = self.A_s @ settings['Q_bias_acc'] @ self.A_s.T
            else:
                raise ValueError("No covariance for bias s")
            
        self.Q = Q
        
        # Update covariance
        if 'R_pos' in settings:
            self.R_pos = settings['R_pos']
            assert np.allclose(self.R_pos, self.R_pos.T)  # Check symmetry
            assert self.R_pos.shape == (3, 3)
            
        if 'R_rot' in settings:
            self.R_rot = settings['R_rot']
            assert np.allclose(self.R_rot, self.R_rot.T)  # Check symmetry
            assert self.R_rot.shape == (3, 3)

    def get_input(self, sensorData: Dict) -> np.ndarray:
        """Get input from sensor data"""
        assert sensorData['gyro_measurements'].shape[0] == 3
        if self.inds_y_a is not None:
            assert sensorData['acc_measurements'].shape[0] == 3 * self.N_a
            return np.vstack((sensorData['gyro_measurements'], sensorData['acc_measurements']))
        else:
            return sensorData['gyro_measurements']

    def get_Q(self) -> np.ndarray:
        """Get process noise covariance matrix"""
        return self.Q

    def propagate(self, R: np.ndarray, x: np.ndarray, y: np.ndarray, w: np.ndarray) -> Dict:
        """
        Propagate the state
        """
        assert len(x) == self.Nx - 3
        assert len(y) == self.Ny
        assert len(w) == self.Nw
        
        v = x[self.inds_v - 3] if self.inds_v is not None else None
        b_s = x[self.inds_b_s - 3] if self.inds_b_s is not None else np.zeros(3)
        b_g = x[self.inds_b_g - 3] if self.inds_b_g is not None else np.zeros(3)
        
        # Process noise
        w_g = w[self.inds_w_g]
        w_alpha = w[self.inds_w_alpha] if self.inds_w_alpha is not None else np.zeros(6)
        w_b_s = w[self.inds_w_b_s] if self.inds_w_b_s is not None else np.zeros(3)
        w_b_g = w[self.inds_w_b_g] if self.inds_w_b_g is not None else np.zeros(3)
        
        if self.inds_y_a is not None:
            y_g = y[self.inds_y_g]
            omega = y_g - b_g - w_g
            y_a = y[self.inds_y_a]
            h = self.HatSO3(omega) @ self.HatSO3(omega) @ self.r
            h = h.flatten()
            alpha = self.A @ (y_a - h) + w_alpha
            
            omega_dot = alpha[0:3]  # Angular acceleration
            s = alpha[3:6] + b_s  # Specific force
            v_dot = self.g + R @ s  # Navigation acceleration
        else:
            omega_dot = np.zeros(3)
            s = np.zeros(3)
            v_dot = np.zeros(3)
        
        Omega = np.zeros(self.Nx)
        # R
        Omega[self.inds_R] = omega * self.T + omega_dot * self.T2_R / 2
        
        # p
        if self.inds_p is not None:
            Omega[self.inds_p] = v * self.T + v_dot * self.T**2 / 2
            
        # v
        if self.inds_v is not None:
            Omega[self.inds_v] = v_dot * self.T
            
        if self.inds_b_g is not None:
            Omega[self.inds_b_g] = w_b_g
            
        if self.inds_b_s is not None:
            Omega[self.inds_b_s] = w_b_s
            
        # Jacobian in x
        dOmega_de = np.zeros((self.Nx, self.Nx))
        
        if self.inds_y_a is not None:
            d_h_d_omega = self.compute_d_h_d_omega(omega)
            d_omega_dot_d_omega = -self.A_omega_dot @ d_h_d_omega
            d_v_dot_d_R = -R @ self.HatSO3(s)
            d_v_dot_d_omega = -R @ self.A_s @ d_h_d_omega
        else:
            d_omega_dot_d_omega = np.zeros((3, 3))
            d_v_dot_d_R = np.zeros((3, 3))
            
        # R equation
        if self.inds_b_g is not None:
            d_omega_d_b_g = -np.eye(3)
            if self.inds_y_a is not None:
                d_omega_dot_d_b_g = d_omega_dot_d_omega @ d_omega_d_b_g
            else:
                d_omega_dot_d_b_g = np.zeros((3, 3))
            dOmega_de[np.ix_(self.inds_R, self.inds_b_g)] = d_omega_d_b_g * self.T + d_omega_dot_d_b_g * self.T2_R / 2
            
        # p equation
        if self.inds_p is not None:
            dOmega_de[np.ix_(self.inds_p, self.inds_R)] = d_v_dot_d_R * self.T**2 / 2
            
            if self.inds_v is not None:
                dOmega_de[np.ix_(self.inds_p, self.inds_v)] = np.eye(3) * self.T
                
            if self.inds_b_g is not None:
                dOmega_de[np.ix_(self.inds_p, self.inds_b_g)] = d_v_dot_d_omega @ d_omega_d_b_g * self.T**2 / 2
                
            if self.inds_b_s is not None:
                d_v_dot_d_b_s = R
                dOmega_de[np.ix_(self.inds_p, self.inds_b_s)] = d_v_dot_d_b_s * self.T**2 / 2
                
        # v equation
        if self.inds_v is not None:
            dOmega_de[np.ix_(self.inds_v, self.inds_R)] = d_v_dot_d_R * self.T
            
            if self.inds_b_g is not None:
                dOmega_de[np.ix_(self.inds_v, self.inds_b_g)] = d_v_dot_d_omega @ d_omega_d_b_g * self.T
                
            if self.inds_b_s is not None:
                dOmega_de[np.ix_(self.inds_v, self.inds_b_s)] = d_v_dot_d_b_s * self.T
                
        # Jacobian in w
        dOmega_dw = np.zeros((self.Nx, self.Nw))
        
        d_omega_dot_d_w_alpha = np.block([np.eye(3), np.zeros((3, 3))])
        d_v_dot_d_w_alpha = np.block([np.zeros((3, 3)), R])
        
        # R equation
        d_omega_d_w_g = -np.eye(3)
        if self.inds_y_a is not None:
            d_omega_dot_d_w_g = d_omega_dot_d_omega @ d_omega_d_w_g
        else:
            d_omega_dot_d_w_g = np.zeros((3, 3))
            
        dOmega_dw[np.ix_(self.inds_R, self.inds_w_g)] = d_omega_d_w_g * self.T + d_omega_dot_d_w_g * self.T2_R / 2
        
        if self.inds_w_alpha is not None:
            dOmega_dw[np.ix_(self.inds_R, self.inds_w_alpha)] = d_omega_dot_d_w_alpha * self.T2_R / 2
            
        # p equation
        if self.inds_p is not None:
            dOmega_dw[np.ix_(self.inds_p, self.inds_w_alpha)] = d_v_dot_d_w_alpha * self.T**2 / 2
            
            if self.inds_w_g is not None:
                dOmega_dw[np.ix_(self.inds_p, self.inds_w_g)] = d_v_dot_d_omega @ d_omega_d_w_g * self.T**2 / 2
                
        # v equation
        if self.inds_v is not None:
            if self.inds_w_g is not None:
                dOmega_dw[np.ix_(self.inds_v, self.inds_w_g)] = d_v_dot_d_omega @ d_omega_d_w_g * self.T
                
            dOmega_dw[np.ix_(self.inds_v, self.inds_w_alpha)] = d_v_dot_d_w_alpha * self.T
            
        # b_s equation
        if self.inds_b_s is not None:
            dOmega_dw[np.ix_(self.inds_b_s, self.inds_w_b_s)] = np.eye(3)
            
        if self.inds_b_g is not None:
            dOmega_dw[np.ix_(self.inds_b_g, self.inds_w_b_g)] = np.eye(3)
            
        return {
            'Omega': Omega,
            'dOmega_de': dOmega_de,
            'dOmega_dw': dOmega_dw,
            'v_dot': v_dot,
            'omega_dot': omega_dot,
            's': s,
            'omega': omega
        }

    def compute_d_h_d_omega(self, w: np.ndarray) -> np.ndarray:
        """Compute derivative of h with respect to omega"""
        row1 = np.arange(0, 3*self.N_a, 3)
        row2 = np.arange(1, 3*self.N_a, 3)
        row3 = np.arange(2, 3*self.N_a, 3)
        d_h_d_omega = np.zeros((3*self.N_a, 3))
        
        r1 = self.r[0, :]
        r2 = self.r[1, :]
        r3 = self.r[2, :]
        
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

    def position_update(self, p_obs: np.ndarray, _, x_in: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Position update"""
        assert self.R_pos is not None
        Q = self.R_pos
        
        H = np.zeros((3, self.Nx))
        H[:, self.inds_p] = np.eye(3)  # p
        
        p_pred = x_in[self.inds_p - 3]
        e = p_obs - p_pred
        
        return e, H, Q

    def rotation_update(self, R_obs: np.ndarray, R_pred: np.ndarray, _) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Rotation update"""
        assert self.R_rot is not None
        Q = self.R_rot
        
        H = np.zeros((3, self.Nx))
        H[:, 0:3] = np.eye(3)  # R
        
        e = self.logSO3(self.invSO3(R_pred) @ R_obs)
        
        return e, H, Q

    def get_initial_conditions(self, initIn: Dict) -> Dict:
        """Get initial conditions"""
        initOut = {}
        use_full = ('x' in initIn or 'P' in initIn or 'R' in initIn)
        use_partial = ('mean' in initIn or 'cov' in initIn)
        
        if use_full and use_partial:
            raise ValueError("Both (R0,x0,P0) and (mean,cov) defined, use only one")
        elif 'x' in initIn and 'P' in initIn:
            initOut['R0'] = initIn['R']
            initOut['x0'] = initIn['x']
            initOut['P0'] = initIn['P']
        elif 'mean' in initIn and 'cov' in initIn:
            m = initIn['mean']
            
            assert m['R'].shape == (3, 3)
            initOut['R0'] = m['R']
            x0 = np.zeros(self.Nx - 3)
            
            if self.inds_p is not None:
                x0[self.inds_p - 3] = m['p']
                
            if self.inds_v is not None:
                x0[self.inds_v - 3] = m['v']
                
            if self.inds_b_s is not None:
                if 'b_s' in m:
                    x0[self.inds_b_s - 3] = m['b_s']
                elif 'b_a' in m:
                    x0[self.inds_b_s - 3] = -self.A_s @ m['b_a']
                else:
                    raise ValueError("No mean value alpha bias")
                    
            initOut['x0'] = x0
            
            c = initIn['cov']
            P0 = np.zeros((self.Nx, self.Nx))
            P0[np.ix_(self.inds_R, self.inds_R)] = c['R']
            
            if self.inds_p is not None:
                P0[np.ix_(self.inds_p, self.inds_p)] = c['p']
            if self.inds_v is not None:
                P0[np.ix_(self.inds_v, self.inds_v)] = c['v']
            if self.inds_b_g is not None:
                P0[np.ix_(self.inds_b_g, self.inds_b_g)] = c['b_g']
            if self.inds_b_s is not None:
                if 'b_s' in m:
                    P0[np.ix_(self.inds_b_s, self.inds_b_s)] = c['b_s']
                elif 'b_a' in m:
                    P0[np.ix_(self.inds_b_s, self.inds_b_s)] = self.A_s @ c['b_a'] @ self.A_s.T
                else:
                    raise ValueError("No cov value alpha bias")
                    
            initOut['P0'] = P0
        else:
            raise ValueError("No initial conditions")
            
        return initOut

    def extract_variables(self, Sin: Dict) -> Dict:
        """Extract variables from state"""
        Sout = {'mean': {}, 'std': {}}
        Sout['mean']['R'] = Sin['R']
        
        if self.inds_p is not None:
            Sout['mean']['p'] = Sin['x'][self.inds_p - 3, :]
        if self.inds_v is not None:
            Sout['mean']['v'] = Sin['x'][self.inds_v - 3, :]
        if self.inds_b_g is not None:
            Sout['mean']['b_g'] = Sin['x'][self.inds_b_g - 3, :]
        if self.inds_b_s is not None:
            Sout['mean']['b_s'] = Sin['x'][self.inds_b_s - 3, :]
            
        Sout['std']['R'] = Sin['std'][self.inds_R, :]
        
        if self.inds_p is not None:
            Sout['std']['p'] = Sin['std'][self.inds_p, :]
        if self.inds_v is not None:
            Sout['std']['v'] = Sin['std'][self.inds_v, :]
        if self.inds_b_g is not None:
            Sout['std']['b_g'] = Sin['std'][self.inds_b_g, :]
        if self.inds_b_s is not None:
            Sout['std']['b_s'] = Sin['std'][self.inds_b_s, :]
            
        return Sout

    def print_info(self, S_init: Dict, SensorData: Dict, settings: Dict) -> None:
        """Print information about the filter"""
        print(f"Sampling time: {settings['T']:.2e} [s]")
        print(f"Sampling freq: {1/settings['T']:.1f} [Hz]")
        print(f"T2 for rotation: {self.T2_R:.1e} [s]")
        print(f"Propagate rotation: {self.inds_R is not None}")
        print(f"Propagate position: {self.inds_p is not None}")
        print(f"Propagate velocity: {self.inds_v is not None}")
        print(f"Propagate bias gyro: {self.inds_b_g is not None}")
        print(f"Propagate bias specific force: {self.inds_b_s is not None}")
        
        print("Propagation data:")
        print("\tAccelerometer data:")
        if self.inds_w_alpha is not None:
            try:
                Q_sqrt = cholesky(self.Q[np.ix_(self.inds_w_alpha, self.inds_w_alpha)])
                print("\t\twhite noise chol(Q_alpha) [mixed units]")
                print(f"\t\t{Q_sqrt.T}")
            except:
                Q_tot = self.Q[np.ix_(self.inds_w_alpha, self.inds_w_alpha)]
                print("\t\twhite noise Q_alpha [mixed units]^2")
                print(f"\t\t{Q_tot.T}")
                
        if self.inds_b_s is not None:
            try:
                Q_sqrt = cholesky(self.Q[np.ix_(self.inds_w_b_s, self.inds_w_b_s)])
                print("\t\tbias specific force: [m/s^2]")
                print(f"\t\t{Q_sqrt.T}")
                print("\t\t[m/s^2]")
            except:
                Q_s = self.Q[np.ix_(self.inds_w_b_s, self.inds_w_b_s)]
                print("\t\tbias specific force: [m/s^2]^2")
                print(f"\t\t{Q_s.T}")
                print("\t\t[m/s^2]^2")
                
        if 'constant_b_s' in settings:
            print("\t\tConstant b_s:          [", end='')
            print(f"{settings['constant_b_s']}", end='')
            print("] []")
            
        print("\t\tStart of accelerometer triad measurements:")
        print(f"\t\t{SensorData['acc_measurements'][0:3, 0:3].T}")
        
        if self.r is not None:
            print("\t\tStart of accelerometer positions:")
            try:
                print(f"\t\t{self.r[0:3, 0:3].T}")
            except:
                pass
            print("\t\tMean accelerometer positions: [", end='')
            print(f"{np.mean(self.r, axis=1)}", end='')
            print("] [m]")
            
        print("\tGyroscope data:")
        if self.inds_b_g is not None:
            try:
                Q_sqrt = np.rad2deg(cholesky(self.Q[np.ix_(self.inds_w_b_g, self.inds_w_b_g)]))
                print("\t\tbias gyro noise chol(Q): [deg/s]")
                print(f"\t\t{Q_sqrt.T}")
                print("\t\t[deg/s]")
            except:
                Q_b_g = np.rad2deg(np.rad2deg(self.Q[np.ix_(self.inds_w_b_g, self.inds_w_b_g)]))
                print("\t\tbias gyro noise Q: [deg/s]^2")
                print(f"\t\t{Q_b_g.T}")
                print("\t\t[deg/s]^2")
                
        print("Initial conditions:")
        print("\tMean:")
        e_deg = np.rad2deg(self.rotm2eul(S_init['mean']['R']))
        print(f"\t\tRotation: roll: {e_deg[0]:.1f}, pitch: {e_deg[1]:.1f}, yaw: {e_deg[2]:.1f} [deg]")
        
        if self.inds_p is not None:
            print("\t\tPosition:          [", end='')
            print(f"{S_init['mean']['p']}", end='')
            print("] [m]")
            
        if self.inds_v is not None:
            print("\t\tVelocity:          [", end='')
            print(f"{S_init['mean']['v']}", end='')
            print("] [m/s]")
            
        if self.inds_b_s is not None:
            print("\t\tBias s:        [", end='')
            print(f"{S_init['mean']['b_s']}", end='')
            print("] [m/s^2]")
            
        if self.inds_b_g is not None:
            print("\t\tBias gyro:        [", end='')
            print(f"{np.rad2deg(S_init['mean']['b_g'])}", end='')
            print("] [deg/s]")
            
        print("\tStd:")
        print("\t\tRotation:          [", end='')
        print(f"{np.rad2deg(S_init['std']['R'])}", end='')
        print("] [deg]")
        
        if self.inds_p is not None:
            print("\t\tPosition:          [", end='')
            print(f"{S_init['std']['p']}", end='')
            print("] [m]")
            
        if self.inds_v is not None:
            print("\t\tVelocity:          [", end='')
            print(f"{S_init['std']['v']}", end='')
            print("] [m/s]")
            
        if self.inds_b_s is not None:
            print("\t\tBias s:        [", end='')
            print(f"{S_init['std']['b_s']}", end='')
            print("] [m/s^2]")
            
        if self.inds_b_g is not None:
            print("\t\tBias gyro:        [", end='')
            print(f"{np.rad2deg(S_init['std']['b_g'])}", end='')
            print("] [deg/s]")

    @staticmethod
    def HatSO3(w: np.ndarray) -> np.ndarray:
        """Convert 3D vector to skew-symmetric matrix"""
        return np.array([
            [0, -w[2], w[1]],
            [w[2], 0, -w[0]],
            [-w[1], w[0], 0]
        ])

    @staticmethod
    def invSO3(R: np.ndarray) -> np.ndarray:
        """Inverse of SO(3) matrix (transpose)"""
        return R.T

    @staticmethod
    def logSO3(R: np.ndarray) -> np.ndarray:
        """Logarithm map of SO(3)"""
        theta = np.arccos((np.trace(R) - 1) / 2)
        if np.abs(theta) < 1e-10:
            return np.zeros(3)
        return theta / (2 * np.sin(theta)) * np.array([
            R[2, 1] - R[1, 2],
            R[0, 2] - R[2, 0],
            R[1, 0] - R[0, 1]
        ])

    @staticmethod
    def rotm2eul(R: np.ndarray) -> np.ndarray:
        """Convert rotation matrix to Euler angles"""
        # Implementation depends on your specific Euler angle convention
        # This is a placeholder - you'll need to implement the specific conversion
        # based on your needs
        pass

    @staticmethod
    def compute_A_non_center(r: np.ndarray) -> np.ndarray:
        """Compute A matrix for non-centered accelerometers"""
        # Implementation depends on your specific requirements
        # This is a placeholder - you'll need to implement the specific computation
        pass