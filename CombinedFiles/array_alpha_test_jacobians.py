import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.transform import Rotation as R
from types import SimpleNamespace

# ==============================================================================
# Helper Functions and Classes (replacing MATLAB's proprietary code)
# ==============================================================================

def expSO3(phi):
    """
    Computes the exponential map for SO(3).
    Converts a 3x1 rotation vector into a 3x3 rotation matrix.
    """
    phi = np.asarray(phi).flatten()
    # Scipy's from_rotvec expects a 1D array
    return R.from_rotvec(phi).as_matrix()

def numeric_jacobian(f, x, h=1e-7):
    """
    Computes the numeric Jacobian of a function f at point x using finite differences.
    """
    x = np.asarray(x).flatten()
    y0 = f(x).flatten()
    J = np.zeros((len(y0), len(x)))
    for i in range(len(x)):
        x_h = x.copy()
        x_h[i] += h
        y_h = f(x_h).flatten()
        J[:, i] = (y_h - y0) / h
    return J

class D_LG_EKF_Array_v4_alpha:
    """
    A placeholder class for the MATLAB model 'D_LG_EKF_Array_v4_alpha'.
    This class mimics the structure and computes dummy outputs to allow
    the main script logic to run.
    """
    def __init__(self, simdata):
        self.simdata = simdata
        self.N_a = 5  # Dummy value for number of accelerometers
        
        # Determine state vector size (Nx) based on flags
        # Base state is rotation (3)
        self.Nx = 3
        if simdata.get('propagate_position', False):
            self.Nx += 3
        if simdata.get('propagate_velocity', False):
            self.Nx += 3
        if simdata.get('propagate_bias_alpha', False):
            self.Nx += 6
        if simdata.get('propagate_bias_gyro', False):
            self.Nx += 3

        # Dummy value for process noise vector size
        self.Nw = 9

    def propagate(self, R_in, e_minus_R, y_a, w):
        """
        A placeholder for the 'propagate' method.
        Returns a namespace object with dummy outputs of the correct shape.
        """
        # In a real implementation, this method would perform state propagation.
        # Here, we generate random data with the expected shapes for demonstration.
        Omega_dim = 3
        
        # Create dummy analytical Jacobians
        # These would be calculated by the actual EKF propagation equations
        F = np.random.randn(Omega_dim, self.Nx) * 0.1
        G = np.random.randn(Omega_dim, self.Nw) * 0.2
        
        # Create a dummy output state
        Omega = np.random.randn(Omega_dim, 1)

        # Use a SimpleNamespace to mimic MATLAB's struct output
        res = SimpleNamespace(
            Omega=Omega,
            dOmega_de=F,
            dOmega_dw=G
        )
        return res

# ==============================================================================
# Main Script
# ==============================================================================

def run_jacobian_test(simdata_dict, test_case_title):
    """
    Runs a single Jacobian test case.
    """
    print(f"--- Running Test Case: {test_case_title} ---")
    simdata = simdata_dict

    # simdata.set_T2_R_zero = True

    Q_alpha_sqrt = np.random.randn(6, 6)
    simdata['Q_alpha'] = Q_alpha_sqrt.T @ Q_alpha_sqrt
    simdata['Q_bias_alpha'] = Q_alpha_sqrt.T @ Q_alpha_sqrt

    Q_bias_gyro_sqrt = np.random.randn(3, 3)
    simdata['Q_bias_gyro'] = Q_bias_gyro_sqrt.T @ Q_bias_gyro_sqrt
    
    # simdata.set_T2_R_zero = True;
    model = D_LG_EKF_Array_v4_alpha(simdata)

    R0 = expSO3(np.random.randn(3, 1))
    y_a = np.random.randn(3 * model.N_a, 1)
    w0 = np.zeros((model.Nw, 1))

    # Important to check jacobian around 0 for SO(3). The jacobian is defined
    # there.
    # e0 = zeros(simdata.Nx,1);
    e0 = np.vstack([np.zeros((3, 1)), np.random.randn(model.Nx - 3, 1)])

    # Define lambda functions to pass to the numeric Jacobian calculator
    f_e = lambda e: model.propagate(R0 @ expSO3(e[:3]), e[3:].reshape(-1, 1), y_a, w0).Omega
    f_w = lambda w: model.propagate(R0 @ expSO3(e0[:3]), e0[3:], y_a, w.reshape(-1, 1)).Omega

    f_e(e0)
    f_w(w0)

    jac_f_e = numeric_jacobian(f_e, e0)
    jac_f_w = numeric_jacobian(f_w, w0)
    res = model.propagate(R0 @ expSO3(e0[:3]), e0[3:], y_a, w0)
    F = res.dOmega_de
    G = res.dOmega_dw
    
    # Calculate and print the norm of the differences
    norm_F_diff = np.linalg.norm(F - jac_f_e)
    norm_G_diff = np.linalg.norm(G - jac_f_w)
    print(f"Norm of difference (F - jac_f_e): {norm_F_diff}")
    print(f"Norm of difference (G - jac_f_w): {norm_G_diff}\n")


    # Plotting
    plt.figure(figsize=(18, 10))
    plt.suptitle(test_case_title, fontsize=16)

    # Difference between analytical and numerical Jacobian F
    plt.subplot(2, 4, 1)
    sns.heatmap(np.log10(np.abs(F - jac_f_e) + 1e-16), vmin=-15, vmax=1)
    plt.title("log10(|F - F_num|)")

    # Analytical Jacobian F
    plt.subplot(2, 4, 2)
    sns.heatmap(F)
    plt.title("Analytical F (J)")

    # Numerical Jacobian F
    plt.subplot(2, 4, 3)
    sns.heatmap(jac_f_e)
    plt.title("Numerical F (J num)")
    
    # Difference between analytical and numerical Jacobian G
    plt.subplot(2, 4, 5)
    sns.heatmap(np.log10(np.abs(G - jac_f_w) + 1e-16), vmin=-15, vmax=1)
    plt.title("log10(|G - G_num|)")

    # Analytical Jacobian G
    plt.subplot(2, 4, 6)
    sns.heatmap(G)
    plt.title("Analytical G")

    # Numerical Jacobian G
    plt.subplot(2, 4, 7)
    sns.heatmap(jac_f_w)
    plt.title("Numerical G (G num)")
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()


# %%
# --- Test Case 1: Propagate all states ---
simdata1 = {
    'T': 0.1,
    'g': np.array([[0], [0], [1]]),
    'r': np.random.randn(3, 3),
    'propagate_position': True,
    'propagate_velocity': True,
    'propagate_bias_alpha': True,
    'propagate_bias_gyro': True,
    'set_T2_R_zero': False
}
run_jacobian_test(simdata1, "Test Case 1: Propagate All States")


# %%
# --- Test Case 2: Propagate only bias states ---
simdata2 = {
    'T': 0.1,
    'g': np.array([[0], [0], [1]]),
    'r': np.random.randn(3, 3),
    'propagate_position': False,
    'propagate_velocity': False,
    'propagate_bias_alpha': True,
    'propagate_bias_gyro': True
}
run_jacobian_test(simdata2, "Test Case 2: Propagate Biases Only")


# %%
# --- Test Case 3: Propagate no states (only rotation) ---
simdata3 = {
    'T': 0.1,
    'g': np.array([[0], [0], [1]]),
    'r': np.random.randn(3, 3),
    'propagate_position': False,
    'propagate_velocity': False,
    'propagate_bias_alpha': False,
    'propagate_bias_gyro': False
}
run_jacobian_test(simdata3, "Test Case 3: Propagate Rotation Only")