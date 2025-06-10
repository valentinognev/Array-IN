#!/usr/bin/env python3
"""
Test jacobians
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.linalg import expm

# Assuming these are custom modules that need to be imported
# from your_module import expSO3, D_LG_EKF_Gyro_2nd_v4, numeric_jacobian

def main():
    # Test jacobians
    simdata = {}
    simdata['T'] = 0.1
    simdata['g'] = np.array([[0], [0], [1]])
    simdata['r'] = np.random.randn(3, 3)
    simdata['propagate_position'] = True
    simdata['propagate_velocity'] = True
    simdata['propagate_bias_s'] = True
    simdata['propagate_bias_gyro'] = True
    simdata['input_accelerometers'] = True
    simdata['set_T2_R_zero'] = True

    Q_gyro_sqrt = np.random.randn(3, 3)
    simdata['Q_gyro'] = Q_gyro_sqrt.T @ Q_gyro_sqrt

    Q_alpha_sqrt = np.random.randn(6, 6)
    simdata['Q_alpha'] = Q_alpha_sqrt.T @ Q_alpha_sqrt

    Q_bias_s_sqrt = np.random.randn(3, 3)
    simdata['Q_bias_s'] = Q_bias_s_sqrt.T @ Q_bias_s_sqrt

    Q_bias_gyro_sqrt = np.random.randn(3, 3)
    simdata['Q_bias_gyro'] = Q_bias_gyro_sqrt.T @ Q_bias_gyro_sqrt
    # simdata['set_T2_R_zero'] = True
    
    model = D_LG_EKF_Gyro_2nd_v4(simdata)

    R0 = expSO3(np.random.randn(3, 1))
    y = np.random.randn(model.Ny, 1)
    w0 = np.zeros((model.Nw, 1))

    # Important to check jacobian around 0 for SO(3). The jacobian is defined
    # there.
    # e0 = np.zeros((simdata['Nx'], 1))
    e0 = np.vstack([np.zeros((3, 1)), np.random.randn(model.Nx - 3, 1)])

    def f_e(e):
        return model.propagate(R0 @ expSO3(e[:3]), e[3:], y, w0).Omega
    
    def f_w(w):
        return model.propagate(R0 @ expSO3(e0[:3]), e0[3:], y, w).Omega

    f_e(e0)
    f_w(w0)

    jac_f_e = numeric_jacobian(f_e, e0)
    jac_f_w = numeric_jacobian(f_w, w0)
    res = model.propagate(R0 @ expSO3(e0[:3]), e0[3:], y, w0)
    F = res.dOmega_de
    G = res.dOmega_dw

    print(f"norm(F-jac_f_e): {np.linalg.norm(F - jac_f_e)}")
    print(f"norm(G-jac_f_w): {np.linalg.norm(G - jac_f_w)}")

    plt.figure()
    sns.heatmap(np.log10(np.abs(F - jac_f_e)), vmin=-13, vmax=0, cmap='viridis')
    plt.title("Diff log10")
    
    plt.figure()
    sns.heatmap(F, cmap='viridis')
    plt.title("J")
    
    plt.figure()
    sns.heatmap(jac_f_e, cmap='viridis')
    plt.title("J num")
    
    plt.figure()
    sns.heatmap(np.log10(np.abs(G - jac_f_w)), vmin=-13, vmax=0, cmap='viridis')
    plt.title("Diff for G (log10)")

    plt.figure(5)
    plt.clf()
    sns.heatmap(G, cmap='viridis')
    plt.title("G")
    
    plt.figure()
    sns.heatmap(jac_f_w, cmap='viridis')
    plt.title("G num")

    plt.show()

if __name__ == "__main__":
    main() 