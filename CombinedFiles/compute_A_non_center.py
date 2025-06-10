import numpy as np

def compute_A_non_center(sensorPos, T):
    # COMPUTE_AS The matrix for the rotation 
    # Based on: Appendix A of "Inertial Navigation using an Inertial sensor array"
    #
    # As = compute_As(sensorPos)
    # Where sensorPos is centered - position relative to centroid.

    K = sensorPos.shape[1] # Number of acc triads

    R_skew = [skew_sym(sensorPos[:,k]) for k in range(K)]
    H = [-np.concatenate(R_skew, axis=0), np.tile(np.eye(3), (K, 1))]

    if T is not None:
        H = matrix3d2blkdiag(T)*H

    A = np.linalg.inv(H.T @ H) @ H.T  # 5b - Inertial Navigation using an Inertial sensor array

    return A



