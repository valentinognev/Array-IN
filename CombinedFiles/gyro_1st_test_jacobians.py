# Test jacobians
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.linalg import expm
from scipy.optimize import approx_fprime

# Helper function for matrix exponential on SO(3)
def expSO3(omega):
    """Matrix exponential for SO(3) - converts axis-angle to rotation matrix"""
    if omega.ndim == 1:
        omega = omega.reshape(-1, 1)
    
    angle = np.linalg.norm(omega)
    if angle < 1e-8:
        return np.eye(3)
    
    axis = omega / angle
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]], 
                  [-axis[1], axis[0], 0]])
    
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * K @ K

def numeric_jacobian(func, x, eps=1e-8):
    """Compute numerical jacobian of function at point x"""
    return approx_fprime(x, func, eps).T

# Placeholder for the model class - you'll need to implement this based on your actual model
class D_LG_EKF_Gyro_1st_v4:
    def __init__(self, simdata):
        self.simdata = simdata
        self.Ny = simdata.get('Ny', 6)  # Assuming some default values
        self.Nw = simdata.get('Nw', 12)
        self.Nx = simdata.get('Nx', 15)
    
    def propagate(self, R, x, y, w):
        # Placeholder implementation - replace with actual propagation logic
        result = type('obj', (object,), {})()
        result.Omega = np.zeros(self.Nx)  # Placeholder
        result.dOmega_de = np.random.randn(self.Nx, self.Nx)  # Placeholder
        result.dOmega_dw = np.random.randn(self.Nx, self.Nw)  # Placeholder
        return result

# Initialize simulation data structure
simdata = {
    'T': 0.1,
    'g': np.array([0, 0, 1]),
    'N_a': 3,
    'propagate_position': True,
    'propagate_velocity': True,
    'propagate_bias_s': True,
    'propagate_bias_gyro': True,
    'input_accelerometers': True,
    'Ny': 6,  # You may need to adjust these based on your model
    'Nw': 12,
    'Nx': 15
}

Q_gyro_sqrt = np.random.randn(3, 3)
simdata['Q_gyro'] = Q_gyro_sqrt.T @ Q_gyro_sqrt

Q_s_sqrt = np.random.randn(3, 3)
simdata['Q_s'] = Q_s_sqrt.T @ Q_s_sqrt

Q_bias_s_sqrt = np.random.randn(3, 3)
simdata['Q_bias_s'] = Q_bias_s_sqrt.T @ Q_bias_s_sqrt

Q_bias_gyro_sqrt = np.random.randn(3, 3)
simdata['Q_bias_gyro'] = Q_bias_gyro_sqrt.T @ Q_bias_gyro_sqrt

# simdata['set_T2_R_zero'] = True
model = D_LG_EKF_Gyro_1st_v4(simdata)

R0 = expSO3(np.random.randn(3, 1))
y = np.random.randn(model.Ny, 1)
w0 = np.zeros((model.Nw, 1))

# Important to check jacobian around 0 for SO(3). The jacobian is defined
# there.
# e0 = np.zeros(simdata['Nx'])
e0 = np.concatenate([np.zeros(3), np.random.randn(model.Nx - 3)])

def f_e(e):
    return model.propagate(R0 @ expSO3(e[:3]), e[3:], y, w0).Omega

def f_w(w):
    return model.propagate(R0 @ expSO3(e0[:3]), e0[3:], y, w).Omega

f_e(e0)
f_w(w0.flatten())

jac_f_e = numeric_jacobian(f_e, e0)
jac_f_w = numeric_jacobian(f_w, w0.flatten())
res = model.propagate(R0 @ expSO3(e0[:3]), e0[3:], y, w0)
F = res.dOmega_de
G = res.dOmega_dw

print(f"Norm difference F: {np.linalg.norm(F - jac_f_e)}")
print(f"Norm difference G: {np.linalg.norm(G - jac_f_w)}")

# Plotting
plt.figure(figsize=(12, 10))

plt.subplot(2, 3, 1)
sns.heatmap(np.log10(np.abs(F - jac_f_e) + 1e-16), vmin=-13, vmax=0, cmap='viridis')
plt.title("Diff log10")

plt.subplot(2, 3, 2)
sns.heatmap(F, cmap='viridis')
plt.title("F")

plt.subplot(2, 3, 3)
sns.heatmap(jac_f_e, cmap='viridis')
plt.title("F num")

plt.subplot(2, 3, 4)
sns.heatmap(np.log10(np.abs(G - jac_f_w) + 1e-16), vmin=-13, vmax=0, cmap='viridis')
plt.title("Diff for G (log10)")

plt.subplot(2, 3, 5)
sns.heatmap(G, cmap='viridis')
plt.title("G")

plt.subplot(2, 3, 6)
sns.heatmap(jac_f_w, cmap='viridis')
plt.title("G num")

plt.tight_layout()
plt.show() 