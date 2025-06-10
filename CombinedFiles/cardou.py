import numpy as np
from scipy.linalg import qr, lstsq
from scipy import linalg
import warnings

def cardou12(acc_b, omega_b_est, pos_b, dt):
    """
    IN:
       acc_b - acceleration measurements
       ombest - omega body estimation
    OUT:
       bdotdot - acceleration of body c.g.
       omega_dot_b - angular acceleration
       omega_b - angular velocity
    """
    
    rb = pos_b   # accelerometer position
    
    npoints = rb.shape[1]
    
    acc = acc_b.flatten()
    evecbase = np.array([[1, 0, 0], 
                         [0, 1, 0], 
                         [0, 0, 1]])
    evec = np.tile(evecbase, (npoints, 1))
    
    mdofs = 3 * npoints
    rvec = np.zeros((3, mdofs))
    for i in range(npoints):
        mat = np.tile(rb[:, i:i+1], (1, 3))
        rvec[:, 3*i:3*i+3] = mat
    
    def CPM(x):
        """Cross Product Matrix"""
        return np.array([[0, -x[2], x[1]], 
                        [x[2], 0, -x[0]], 
                        [-x[1], x[0], 0]])
    
    def Sigmaa(x):
        """Sigma matrix function"""
        return np.array([[0, -x[0], -x[0], 0, x[2], x[1]], 
                        [-x[1], 0, -x[1], x[2], 0, x[0]], 
                        [-x[2], -x[2], 0, x[1], x[0], 0]])
    
    Ap = evec
    R = np.zeros((3, 3*mdofs))
    F = np.zeros((mdofs, 3*mdofs))
    Sigma = np.zeros((6, 3*mdofs))
    
    for i in range(mdofs):
        F[i, 3*i:3*i+3] = evec[i, :]
        R[:, 3*i:3*i+3] = CPM(rvec[:, i])
        Sigma[:, 3*i:3*i+3] = Sigmaa(rvec[:, i]).T
    
    At = F @ R.T
    Ar = F @ Sigma.T
    A = np.hstack([Ap, At, Ar])
    
    res = linalg.lstsq(A, acc)[0]
    
    bdotdot = res[0:3]
    omega_dot_b = res[3:6]
    w0sq = res[6]
    w1sq = res[7]
    w2sq = res[8]
    w1w2 = res[9]
    w0w2 = res[10]
    w0w1 = res[11]
    
    Ws = np.array([[-w1sq-w2sq, w0w1, w0w2], 
                   [w0w1, -w0sq-w2sq, w1w2], 
                   [w0w2, w1w2, -w0sq-w1sq]])
    
    omb = omega_b_est
    if np.linalg.norm(omega_b_est) < 1e-6:
        omb = omega_dot_b * dt
    
    wCANP = calcCANP(Ws, omb).T
    wCAD = calcCAD(Ws, omb).T
    wCAAD = calcCAAD(Ws, omb).T
    wCAAM = calcCAAM(Ws, omb).T
    omega_b = wCAAM
    
    ### test
    omega_b = wCAAM
    ksiCAAM = np.hstack([bdotdot.T, omega_dot_b.T, omega_b[0]**2, omega_b[1]**2, omega_b[2]**2, 
                        omega_b[1]*omega_b[2], omega_b[0]*omega_b[2], omega_b[0]*omega_b[1]])
    testCAAM = (A @ ksiCAAM.T - acc) / acc * 100
    
    omega_b = wCAD
    ksiCAD = np.hstack([bdotdot.T, omega_dot_b.T, omega_b[0]**2, omega_b[1]**2, omega_b[2]**2, 
                       omega_b[1]*omega_b[2], omega_b[0]*omega_b[2], omega_b[0]*omega_b[1]])
    testCAD = (A @ ksiCAD.T - acc) / acc * 100
    
    omega_b = wCANP
    ksiCANP = np.hstack([bdotdot.T, omega_dot_b.T, omega_b[0]**2, omega_b[1]**2, omega_b[2]**2, 
                        omega_b[1]*omega_b[2], omega_b[0]*omega_b[2], omega_b[0]*omega_b[1]])
    testCANP = (A @ ksiCANP.T - acc) / acc * 100
    
    omega_b = wCAAD
    ksiCAAD = np.hstack([bdotdot.T, omega_dot_b.T, omega_b[0]**2, omega_b[1]**2, omega_b[2]**2, 
                        omega_b[1]*omega_b[2], omega_b[0]*omega_b[2], omega_b[0]*omega_b[1]])
    testCAAD = (A @ ksiCAAD.T - acc) / acc * 100
    
    print('')
    
    return bdotdot, omega_dot_b, omega_b


#########################################################################################################
## calcCANP
def calcCANP(Ws, wTA):
    """CANP calculation function"""
    if wTA.shape[1] == 3:
        wTA = wTA.T
    
    mu2 = -(Ws[0,0] + Ws[1,1] + Ws[2,2])
    s12 = Ws[0,1]**2
    s23 = Ws[1,2]**2
    s31 = Ws[2,0]**2
    s13 = Ws[0,1] * Ws[1,2]
    d12 = Ws[0,0] * Ws[1,1]
    mu1 = d12 + Ws[1,1]*Ws[2,2] + Ws[0,0]*Ws[2,2] - s12 - s23 - s31
    mu0 = -d12*Ws[2,2] - 2*s13*Ws[2,0] + s12*Ws[2,2] + s23*Ws[0,0] + s31*Ws[1,1]
    ni2 = mu2/3
    theta2 = ni2**2
    q = mu1/3 - theta2
    
    if q >= 0:
        wCANP = np.array([0, 0, 0])
        return wCANP
    else:
        r = (mu1*ni2 - mu0) / 2 - ni2*theta2
        alpha = np.sqrt(-q)
        beta = alpha**3
    
    if beta <= r:
        wCANP = np.array([0, 0, 0])
        return wCANP
    else:
        lam = 2 * alpha * np.cos(np.arccos(r / beta) / 3) - ni2
        delta = (lam + mu2) / 2
    
    if delta <= 0 or (lam * mu0) > 0:
        wCANP = np.array([0, 0, 0])
        return wCANP
    else:
        wcanpnorm = np.sqrt(delta)
        zeta11 = Ws[0,0] - lam
        zeta22 = Ws[1,1] - lam
        zeta33 = Ws[2,2] - lam
        ksi11 = zeta22 * zeta33 - s23
        ksi22 = zeta33 * zeta11 - s31
        ksi33 = zeta11 * zeta22 - s12
        ksi12 = Ws[1,2] * Ws[2,0] - Ws[0,1] * Ws[2,2]
        ksi23 = Ws[0,1] * Ws[2,0] - Ws[1,2] * Ws[0,0]
        ksi31 = s13 - Ws[2,0] * Ws[1,1]
        adjX = np.array([[ksi11, ksi12, ksi31], 
                        [ksi12, ksi22, ksi23], 
                        [ksi31, ksi23, ksi33]])
        v = adjX @ wTA
        vunit = v / np.linalg.norm(v)
        wCANP = wcanpnorm * vunit.T
        return wCANP


#########################################################################################################
## calcCAD
def calcCAD(Ws, wTA):
    """CAD calculation function"""
    trWs = np.trace(Ws)
    if trWs < 0 and np.sum(np.abs(wTA)) > 0:
        zeta0 = Ws[0,0] - 0.5 * trWs
        zeta1 = Ws[1,1] - 0.5 * trWs
        zeta2 = Ws[2,2] - 0.5 * trWs
        
        def heaviside(x):
            return np.where(x >= 0, 1, 0)
        
        wCAD0 = np.sign(wTA[0]) * heaviside(zeta0) * np.sqrt(np.abs(zeta0))
        wCAD1 = np.sign(wTA[1]) * heaviside(zeta1) * np.sqrt(np.abs(zeta1))
        wCAD2 = np.sign(wTA[2]) * heaviside(zeta2) * np.sqrt(np.abs(zeta2))
        wCAD = np.array([wCAD0, wCAD1, wCAD2])
        return wCAD
    else:
        wCAD = np.array([0, 0, 0])
        return wCAD


#########################################################################################################
## calcCAAD
def calcCAAD(Ws, wTA):
    """CAAD calculation function"""
    if wTA.shape[1] == 3:
        wTA = wTA.T
    
    trWs = np.trace(Ws)
    if trWs < 0 and np.sum(np.abs(wTA)) > 0:
        adjWs = np.linalg.inv(Ws) * np.linalg.det(Ws)
        wCAADnorm = np.sqrt(-0.5 * trWs)
        v = adjWs @ wTA
        vunit = v / np.linalg.norm(v)
        wCAAD = wCAADnorm * vunit.T
        return wCAAD
    else:
        wCAAD = np.array([0, 0, 0])
        return wCAAD


#########################################################################################################
## calcCAAM
def calcCAAM(Ws, wTA):
    """CAAM calculation function"""
    if wTA.shape[0] == 3:
        wTA = wTA.T
    
    trWs = np.trace(Ws)
    if trWs < 0:
        adjWs = np.linalg.inv(Ws) * np.linalg.det(Ws)
        wCAAMnorm = np.sqrt(-0.5 * trWs)
        Xtop = Ws
        Xbot = wTA @ adjWs
        X = np.vstack([Xtop, Xbot])
        Y = np.array([0, 0, 0, (-trWs/2)**3])
        Q, R = qr(X)
        uw = linalg.solve(R, Q.T @ Y.T)
        uwunit = uw / np.linalg.norm(uw)
        wCAAM = wCAAMnorm * uwunit.T
        return wCAAM
    else:
        wCAAM = np.array([0, 0, 0])
        return wCAAM


def cardou9(acc_b, omega_b_est, pos_b, dt):
    """
    Cardou, 2010, Computing the Rigid body acceleration field from nine accelerometer measurements
    
    IN:
       acc_b - acceleration measurements
       ombest - omega body estimation
    OUT:
       bdotdot - acceleration of body c.g.
       omega_dot_b - angular acceleration
       omega_b - angular velocity
    """
    bdotdot = np.array([0.9981, -66.6620, 1.0000])
    omega_dot_b = np.array([2.9999, 2.0002, 0.9999])
    omega = np.array([0.1022, 0.2015, 0.3007])
    
    pos_b = pos_b[:, [0, 1, 3]]
    acc_b = acc_b[[0, 1, 2, 3, 4, 5, 9, 10, 11]]
    # bdotdot=[]; omega_dot_b=[]; omega_b=[];
    rb = pos_b   # accelerometer position
    
    npoints = rb.shape[1]
    
    acc = acc_b.flatten()
    evecbase = np.array([[1, 0, 0], 
                         [0, 1, 0], 
                         [0, 0, 1]])
    evec = np.tile(evecbase, (npoints, 1))
    
    mdofs = 3 * npoints
    rvec = np.zeros((3, mdofs))
    for i in range(npoints):
        mat = np.tile(rb[:, i:i+1], (1, 3))
        rvec[:, 3*i:3*i+3] = mat
    
    def CPM(x):
        """Cross Product Matrix"""
        return np.array([[0, -x[2], x[1]], 
                        [x[2], 0, -x[0]], 
                        [-x[1], x[0], 0]])
    
    def Sigmaa(x):
        """Sigma matrix function"""
        return np.array([[0, -x[0], -x[0], 0, x[2], x[1]], 
                        [-x[1], 0, -x[1], x[2], 0, x[0]], 
                        [-x[2], -x[2], 0, x[1], x[0], 0]])
    
    Ap = evec
    R = np.zeros((3, 3*mdofs))
    F = np.zeros((mdofs, 3*mdofs))
    Sigma = np.zeros((6, 3*mdofs))
    
    for i in range(mdofs):
        F[i, 3*i:3*i+3] = evec[i, :]
        R[:, 3*i:3*i+3] = CPM(rvec[:, i])
        Sigma[:, 3*i:3*i+3] = Sigmaa(rvec[:, i]).T
    
    At = F @ R.T
    Ar = F @ Sigma.T
    
    # Commented out parts that were specific to Python and need adaptation
    # # Ar = Ar.' * rho;
    A = np.hstack([Ap, At, Ar])
    Q, S = qr(A)
    
    Q1 = Q[:, 0:6]
    Q2 = Q[:, 6:9]
    if Q.shape[1] > 9:
        Q3 = Q[:, 9:]
    else:
        Q3 = np.array([])
    
    S11 = S[0:6, 0:6]
    S22 = S[6:9, 6:12]
    S12 = S[0:6, 6:12]
    O3x6 = S[6:9, 0:6]
    if S.shape[1] > 12:
        S32 = S[:, 6:12]
        Onm9x6 = S[:, 0:6]
    else:
        S32 = np.array([])
        Onm9x6 = np.array([])
    
    rankA = np.linalg.matrix_rank(A)
    Q2a = Q2.T @ acc
    
    s11 = S22[0, 0]; s12 = S22[0, 1]; s13 = S22[0, 2]; s14 = S22[0, 3]; s15 = S22[0, 4]; s16 = S22[0, 5]
    s22 = S22[1, 1]; s23 = S22[1, 2]; s24 = S22[1, 3]; s25 = S22[1, 4]; s26 = S22[1, 5]
    s33 = S22[2, 2]; s34 = S22[2, 3]; s35 = S22[2, 4]; s36 = S22[2, 5]
    q7a = Q2a[0]; q8a = Q2a[1]; q9a = Q2a[2]
    
    v0, v1, v2, v3, v4 = getVcoeffs(s11, s12, s13, s14, s15, s16, s22, s23, s24, s25, s26, s33, s34, s35, s36, q7a, q8a, q9a)
    zeta = np.roots([v4, v3, v2, v1, v0])
    
    zeta_0, zeta_1, zeta_2, zeta_3 = quarticEquationRoots(v0, v1, v2, v3, v4)
    zeta2 = np.array([zeta_0, zeta_1, zeta_2, zeta_3]).T
    w3arr = []
    
    for i in range(len(zeta)):
        w3p = np.sqrt(zeta[i])
        w3n = -np.sqrt(zeta[i])
        w3arr.extend([w3p, w3n])
    
    wArr = []
    wdotArr = []
    bddArr = []
    
    for w3 in w3arr:
        u1, u2, u3, u4, u5, u6, u7, u8, u9, u10 = getuvector(s11, s12, s13, s14, s15, s16, s22, s23, s24, s25, s26, s33, s34, s35, s36, q7a, q8a, q9a, w3)
        
        C1 = np.array([[s11, s16, s15*w3, s12, s14*w3],
                      [0, s26, s25*w3, s22, s24*w3],
                      [0, s36, s35*w3, 0, s34*w3],
                      [3*u1, 2*u2, 2*u3, u4, u5],
                      [u2, 2*u4, u5, 3*u7, 2*u8],
                      [u3, u5, 2*u6, u8, 2*u9]])
        
        c2 = -np.array([s13*w3**2 - q7a,
                       s23*w3**2 - q8a,
                       s33*w3**2 - q9a,
                       u6,
                       u9,
                       3*u10])
        
        rankC1 = np.linalg.matrix_rank(C1)
        print(f'rankC1: {rankC1}')
        
        Q_2, S_2 = qr(C1)
        lhs = Q_2.T @ c2
        ww4 = lhs[4] / S_2[4, 4]
        ww3 = (lhs[3] - S_2[3, 4] * ww4) / S_2[3, 3]
        ww2 = (lhs[2] - S_2[2, 3] * ww3 - S_2[2, 4] * ww4) / S_2[2, 2]
        ww1 = (lhs[1] - S_2[1, 2] * ww2 - S_2[1, 3] * ww3 - S_2[1, 4] * ww4) / S_2[1, 1]
        ww0 = (lhs[0] - S_2[0, 1] * ww1 - S_2[0, 2] * ww2 - S_2[0, 3] * ww3 - S_2[0, 4] * ww4) / S_2[0, 0]
        wwi = np.array([ww0, ww1, ww2, ww3, ww4])
        
        if ww0 < 0 or ww3 < 0:
            continue
        
        w1abs = np.real(np.sqrt(ww0))
        w2abs = np.real(np.sqrt(ww3))
        
        if np.isclose(np.abs(np.real(ww2)), w1abs) and np.isclose(np.abs(np.real(ww4)), w2abs):
            wArr.append(np.real([ww2, ww4, w3]))
        elif np.real(ww1) >= 0:
            wArr.extend([[w1abs, w2abs, np.real(w3)], [-w1abs, -w2abs, np.real(w3)]])
        else:
            wArr.extend([[w1abs, -w2abs, np.real(w3)], [-w1abs, w2abs, np.real(w3)]])
    
    wArr = np.unique(np.real(wArr), axis=0)
    
    for i in range(wArr.shape[0]):
        w1 = wArr[i, 0]
        w2 = wArr[i, 1]
        w3 = wArr[i, 2]
        
        ksii = np.array([w1**2, w2**2, w3**2, w2*w3, w1*w3, w1*w2])
        Q1a = Q1.T @ acc
        S12ksii = S12 @ ksii.T
        Q1amS12ksii = Q1a - S12ksii
        
        xpt5 = Q1amS12ksii[5] / S11[5, 5]
        xpt4 = (Q1amS12ksii[4] - S11[4, 5] * xpt5) / S11[4, 4]
        xpt3 = (Q1amS12ksii[3] - S11[3, 4] * xpt4 - S11[3, 5] * xpt5) / S11[3, 3]
        xpt2 = (Q1amS12ksii[2] - S11[2, 3] * xpt3 - S11[2, 4] * xpt4 - S11[2, 5] * xpt5) / S11[2, 2]
        xpt1 = (Q1amS12ksii[1] - S11[1, 2] * xpt2 - S11[1, 3] * xpt3 - S11[1, 4] * xpt4 - S11[1, 5] * xpt5) / S11[1, 1]
        xpt0 = (Q1amS12ksii[0] - S11[0, 1] * xpt1 - S11[0, 2] * xpt2 - S11[0, 3] * xpt3 - S11[0, 4] * xpt4 - S11[0, 5] * xpt5) / S11[0, 0]
        xpti = np.array([xpt0, xpt1, xpt2, xpt3, xpt4, xpt5])
        
        bdotdot = xpti[0:3]
        wdot = xpti[3:6]
        wdotArr.append(wdot)
        bddArr.append(bdotdot)
    
    wdotArr = np.real(wdotArr)
    bddArr = np.real(bddArr)
    
    for i in range(wArr.shape[0]):
        omdotb_ = wdotArr[i, :]
        omb_ = wArr[i, :]
        sb_ = bddArr[i, :]
        
        ksii = np.array([omb_[0]**2, omb_[1]**2, omb_[2]**2, omb_[1] * omb_[2], omb_[0] * omb_[2], omb_[0] * omb_[1]])
        test = Ap @ sb_.T + At @ omdotb_.T + Ar @ ksii.T - acc
        # additional tests or operations can go here \
        test2 = []
        
        for j in range(len(rvec)):
            rb_ = rvec[j, :]
            fb_ = acc[j]
            test2.append(np.dot(evec[j, :], sb_ + np.cross(omdotb_, rb_) + np.cross(omb_, np.cross(omb_, rb_))) - fb_)
    
    print('')
    
    return bdotdot, omega_dot_b, omega_b


def cubicEquationRoots(v0, v1, v2, v3):
    """Cubic equation roots solver"""
    term1 = -v2/(3*v3)
    
    discriminant = 4*(-v2**2 + 3*v1*v3)**3 + (-2*v2**3 + 9*v1*v2*v3 - 27*v0*v3**2)**2
    sqrt_disc = np.sqrt(discriminant + 0j)
    
    cbrt_term = (-2*v2**3 + 9*v1*v2*v3 - 27*v0*v3**2 + sqrt_disc)**(1/3)
    
    term2 = -(2**(1/3)*(-v2**2 + 3*v1*v3))/(3*v3*cbrt_term)
    term3 = cbrt_term/(3*2**(1/3)*v3)
    
    zeta0 = term1 + term2 + term3
    
    omega = np.exp(2j*np.pi/3)  # cube root of unity
    
    zeta1 = term1 + term2*omega**2 + term3*omega
    zeta2 = term1 + term2*omega + term3*omega**2
    
    roots = np.array([zeta0, zeta1, zeta2])
    return roots


def quarticEquationRoots(v0, v1, v2, v3, v4):
    """Quartic equation roots solver"""
    sqrtV = np.sqrt(-4 * (v2**2 - 3 * v1 * v3 + 12 * v0 * v4)**3 + (2 * v2**3 - 9 * v1 * v2 * v3 + 27 * v0 * v3**2 + 27 * v1**2 * v4 - 72 * v0 * v2 * v4)**2 + 0j)
    
    sqrtA2_ = 0j + v3**2 / (2 * v4**2) - (4 * v2) / (3 * v4) - (2**(1/3) * (v2**2 - 3 * v1 * v3 + 12 * v0 * v4)) / (3 * v4 * (2 * v2**3 - 9 * v1 * v2 * v3 + 27 * v0 * v3**2 + 27 * v1**2 * v4 - 72 * v0 * v2 * v4 + sqrtV)**(1/3)) - (2 * v2**3 - 9 * v1 * v2 * v3 + 27 * v0 * v3**2 + 27 * v1**2 * v4 - 72 * v0 * v2 * v4 + sqrtV)**(1/3) / (3 * 2**(1/3) * v4)

    sqrtA3_ = (-(v3**3 / v4**3) + (4 * v2 * v3) / v4**2 - (8 * v1) / v4) / (4 * np.sqrt(v3**2 / (4 * v4**2) - (2 * v2) / (3 * v4) + (2**(1/3) * (v2**2 - 3 * v1 * v3 + 12 * v0 * v4)) / (3 * v4 * (2 * v2**3 - 9 * v1 * v2 * v3 + 27 * v0 * v3**2 + 27 * v1**2 * v4 - 72 * v0 * v2 * v4 + sqrtV)**(1/3)) + (2 * v2**3 - 9 * v1 * v2 * v3 + 27 * v0 * v3**2 + 27 * v1**2 * v4 - 72 * v0 * v2 * v4 + sqrtV)**(1/3) / (3 * 2**(1/3) * v4)))
    
    sqrtA1 = np.sqrt(0j + v3**2 / (4 * v4**2) - (2 * v2) / (3 * v4) + (2**(1/3) * (v2**2 - 3 * v1 * v3 + 12 * v0 * v4)) / (3 * v4 * (2 * v2**3 - 9 * v1 * v2 * v3 + 27 * v0 * v3**2 + 27 * v1**2 * v4 - 72 * v0 * v2 * v4 + sqrtV)**(1/3)) + (2 * v2**3 - 9 * v1 * v2 * v3 + 27 * v0 * v3**2 + 27 * v1**2 * v4 - 72 * v0 * v2 * v4 + sqrtV)**(1/3) / (3 * 2**(1/3) * v4))
    
    zeta0 = -(v3 / (4 * v4)) - (1/2) * sqrtA1 - (1/2) * np.sqrt(sqrtA2_ - sqrtA3_)
    zeta1 = -(v3 / (4 * v4)) - (1/2) * sqrtA1 + (1/2) * np.sqrt(sqrtA2_ - sqrtA3_)
    zeta2 = -(v3 / (4 * v4)) + (1/2) * sqrtA1 - (1/2) * np.sqrt(sqrtA2_ + sqrtA3_)
    zeta3 = -(v3 / (4 * v4)) + (1/2) * sqrtA1 + (1/2) * np.sqrt(sqrtA2_ + sqrtA3_)
    
    return zeta0, zeta1, zeta2, zeta3


def getuvector(s11, s12, s13, s14, s15, s16, s22, s23, s24, s25, s26, s33, s34, s35, s36, q7a, q8a, q9a, w3):
    """Get U vector components"""
    u1 = (2*s11*s26*s35 - 2*s11*s25*s36)*w3
    u2 = (2*s11*s26*s34 + 4*s11*s22*s35 - 2*s11*s24*s36)*w3
    u3 = (4*s11*s26*s33 - 2*s11*s25*s34 + 2*s11*s24*s35 - 4*s11*s23*s36)*w3**2 + 4*q8a*s11*s36 - 4*q9a*s11*s26
    u4 = (4*s11*s22*s34 + 2*s16*s22*s35 - 2*s12*s26*s35 - 2*s15*s22*s36 + 2*s12*s25*s36)*w3
    u5 = (8*s11*s22*s33 - 2*s16*s25*s34 + 2*s15*s26*s34 + 2*s16*s24*s35 - 2*s14*s26*s35 - 2*s15*s24*s36 + 2*s14*s25*s36)*w3**2 - 8*q9a*s11*s22
    u6 = (4*s11*s24*s33 - 2*s16*s25*s33 + 2*s15*s26*s33 - 4*s11*s23*s34 + 2*s16*s23*s35 - 2*s13*s26*s35 - 2*s15*s23*s36 + 2*s13*s25*s36)*w3**3 + \
        (-4*q9a*s11*s24 + 2*q9a*s16*s25 - 2*q9a*s15*s26 + 4*q8a*s11*s34 - 2*q8a*s16*s35 + 2*q7a*s26*s35 + 2*q8a*s15*s36 - 2*q7a*s25*s36)*w3
    u7 = (2*s16*s22*s34 - 2*s12*s26*s34 - 2*s14*s22*s36 + 2*s12*s24*s36)*w3
    u8 = (4*s16*s22*s33 - 4*s12*s26*s33 + 2*s15*s22*s34 - 2*s12*s25*s34 - 2*s14*s22*s35 + 2*s12*s24*s35 - 4*s13*s22*s36 + 4*s12*s23*s36)*w3**2 + \
        (-4*q9a*s16*s22 + 4*q9a*s12*s26 - 4*q8a*s12*s36 + 4*q7a*s22*s36)
    u9 = (4*s15*s22*s33 + 2*s16*s24*s33 - 4*s12*s25*s33 - 2*s14*s26*s33 - 2*s16*s23*s34 + 2*s13*s26*s34 - 4*s13*s22*s35 + 4*s12*s23*s35 + 2*s14*s23*s36 - 2*s13*s24*s36)*w3**3 + \
        (-4*q9a*s15*s22 - 2*q9a*s16*s24 + 4*q9a*s12*s25 + 2*q9a*s14*s26 + 2*q8a*s16*s34 - 2*q7a*s26*s34 - 4*q8a*s12*s35 + 4*q7a*s22*s35 - 2*q8a*s14*s36 + 2*q7a*s24*s36)*w3
    u10 = -2*q9a*s15*s24*w3**2 + 2*q9a*s14*s25*w3**2 + 2*q8a*s15*s34*w3**2 - 2*q7a*s25*s34*w3**2 - 2*q8a*s14*s35*w3**2 + 2*q7a*s24*s35*w3**2 + 2*s15*s24*s33*w3**4 - \
        2*s14*s25*s33*w3**4 - 2*s15*s23*s34*w3**4 + 2*s13*s25*s34*w3**4 + 2*s14*s23*s35*w3**4 - 2*s13*s24*s35*w3**4
    
    return u1, u2, u3, u4, u5, u6, u7, u8, u9, u10 


def getVcoeffs(s11,s12,s13,s14,s15,s16,s22,s23,s24,s25,s26,s33,s34,s35,s36,q7a,q8a,q9a)
    v0=512*s11**2*(q9a**2*(s11*s22**2 + s26*((-s16)*s22 + s12*s26)) +  \
    q9a*(q8a*s16*s22 - 2*q8a*s12*s26 + q7a*s22*s26)*s36 +  \
    q8a*(q8a*s12 - q7a*s22)*s36**2)**2

    v1_0=-2048*q9a**3*s11**4*s22**4*s33 + 4096*q9a**3*s11**3*s16*s22**3*s26*s33 - 4096*q9a**3*s11**3*s12*s22**2*s26**2*s33 - 2048*q9a**3*s11**2*s16**2*s22**2*s26**2*s33 + 4096*q9a**3*s11**2*s12*s16*s22*s26**3*s33 -  \
    2048*q9a**3*s11**2*s12**2*s26**4*s33 + 1024*q9a**3*s11**4*s22**3*s24*s34 - 512*q9a**3*s11**3*s16*s22**3*s25*s34 - 512*q9a**3*s11**3*s15*s22**3*s26*s34 - 1536*q9a**3*s11**3*s16*s22**2*s24*s26*s34 +  \
    1024*q9a**3*s11**3*s12*s22**2*s25*s26*s34 + 512*q9a**3*s11**2*s16**2*s22**2*s25*s26*s34 + 512*q9a**3*s11**3*s14*s22**2*s26**2*s34 + 512*q9a**3*s11**2*s15*s16*s22**2*s26**2*s34 + 1024*q9a**3*s11**3*s12*s22*s24*s26**2*s34 +  \
    512*q9a**3*s11**2*s16**2*s22*s24*s26**2*s34 - 1536*q9a**3*s11**2*s12*s16*s22*s25*s26**2*s34 - 512*q9a**3*s11**2*s12*s15*s22*s26**3*s34 - 512*q9a**3*s11**2*s14*s16*s22*s26**3*s34 - 512*q9a**3*s11**2*s12*s16*s24*s26**3*s34 +  \
    1024*q9a**3*s11**2*s12**2*s25*s26**3*s34 + 512*q9a**3*s11**2*s12*s14*s26**4*s34 - 1024*q8a*q9a**2*s11**4*s22**3*s34**2 + 1536*q8a*q9a**2*s11**3*s16*s22**2*s26*s34**2 - 1024*q8a*q9a**2*s11**3*s12*s22*s26**2*s34**2 -  \
    512*q8a*q9a**2*s11**2*s16**2*s22*s26**2*s34**2 - 512*q7a*q9a**2*s11**3*s22**2*s26**2*s34**2 + 512*q8a*q9a**2*s11**2*s12*s16*s26**3*s34**2 + 512*q7a*q9a**2*s11**2*s16*s22*s26**3*s34**2 - 512*q7a*q9a**2*s11**2*s12*s26**4*s34**2 +  \
    1024*q9a**3*s11**3*s15*s22**4*s35 - 512*q9a**3*s11**3*s16*s22**3*s24*s35 - 1024*q9a**3*s11**3*s12*s22**3*s25*s35 + 512*q9a**3*s11**2*s16**2*s22**3*s25*s35 - 512*q9a**3*s11**3*s14*s22**3*s26*s35 -  \
    1536*q9a**3*s11**2*s15*s16*s22**3*s26*s35 + 1024*q9a**3*s11**3*s12*s22**2*s24*s26*s35 + 512*q9a**3*s11**2*s16**2*s22**2*s24*s26*s35 + 512*q9a**3*s11**2*s12*s16*s22**2*s25*s26*s35 - 512*q9a**3*s11*s16**3*s22**2*s25*s26*s35 +  \
    1536*q9a**3*s11**2*s12*s15*s22**2*s26**2*s35 + 512*q9a**3*s11**2*s14*s16*s22**2*s26**2*s35 + 512*q9a**3*s11*s15*s16**2*s22**2*s26**2*s35 - 1536*q9a**3*s11**2*s12*s16*s22*s24*s26**2*s35 - 1024*q9a**3*s11**2*s12**2*s22*s25*s26**2*s35 +  \
    1024*q9a**3*s11*s12*s16**2*s22*s25*s26**2*s35 - 512*q9a**3*s11**2*s12*s14*s22*s26**3*s35 - 1024*q9a**3*s11*s12*s15*s16*s22*s26**3*s35 + 1024*q9a**3*s11**2*s12**2*s24*s26**3*s35 - 512*q9a**3*s11*s12**2*s16*s25*s26**3*s35 +  \
    512*q9a**3*s11*s12**2*s15*s26**4*s35 + 1024*q8a*q9a**2*s11**3*s16*s22**3*s34*s35 - 2048*q8a*q9a**2*s11**3*s12*s22**2*s26*s34*s35 - 1024*q8a*q9a**2*s11**2*s16**2*s22**2*s26*s34*s35 + 1024*q7a*q9a**2*s11**3*s22**3*s26*s34*s35 +  \
    3072*q8a*q9a**2*s11**2*s12*s16*s22*s26**2*s34*s35 - 1024*q7a*q9a**2*s11**2*s16*s22**2*s26**2*s34*s35 - 2048*q8a*q9a**2*s11**2*s12**2*s26**3*s34*s35 + 1024*q7a*q9a**2*s11**2*s12*s22*s26**3*s34*s35 +  \
    1024*q8a*q9a**2*s11**3*s12*s22**3*s35**2 - 512*q8a*q9a**2*s11**2*s16**2*s22**3*s35**2 - 1024*q7a*q9a**2*s11**3*s22**4*s35**2 - 512*q8a*q9a**2*s11**2*s12*s16*s22**2*s26*s35**2 + 512*q8a*q9a**2*s11*s16**3*s22**2*s26*s35**2 +  \
    1536*q7a*q9a**2*s11**2*s16*s22**3*s26*s35**2 + 1024*q8a*q9a**2*s11**2*s12**2*s22*s26**2*s35**2 - 1024*q8a*q9a**2*s11*s12*s16**2*s22*s26**2*s35**2 - 1536*q7a*q9a**2*s11**2*s12*s22**2*s26**2*s35**2 -  \
    512*q7a*q9a**2*s11*s16**2*s22**2*s26**2*s35**2 + 512*q8a*q9a**2*s11*s12**2*s16*s26**3*s35**2 + 1024*q7a*q9a**2*s11*s12*s16*s22*s26**3*s35**2 - 512*q7a*q9a**2*s11*s12**2*s26**4*s35**2 - 1024*q9a**3*s11**3*s16*s22**3*s23*s36 -  \
    512*q9a**3*s11**3*s15*s22**3*s24*s36 + 512*q9a**3*s11**3*s16*s22**2*s24**2*s36 - 1536*q9a**3*s11**3*s14*s22**3*s25*s36 + 1536*q9a**3*s11**2*s15*s16*s22**3*s25*s36 + 2048*q9a**3*s11**3*s12*s22**2*s24*s25*s36 -  \
    1024*q9a**3*s11**2*s16**2*s22**2*s24*s25*s36 - 1536*q9a**3*s11**2*s12*s16*s22**2*s25**2*s36 + 512*q9a**3*s11*s16**3*s22**2*s25**2*s36 - 1024*q9a**3*s11**3*s13*s22**3*s26*s36 + 512*q9a**3*s11**2*s15**2*s22**3*s26*s36 +  \
    2048*q9a**3*s11**3*s12*s22**2*s23*s26*s36 + 1024*q9a**3*s11**2*s16**2*s22**2*s23*s26*s36 + 1536*q9a**3*s11**3*s14*s22**2*s24*s26*s36 - 512*q9a**3*s11**2*s15*s16*s22**2*s24*s26*s36 - 2048*q9a**3*s11**3*s12*s22*s24**2*s26*s36 -  \
    2560*q9a**3*s11**2*s12*s15*s22**2*s25*s26*s36 + 512*q9a**3*s11**2*s14*s16*s22**2*s25*s26*s36 - 512*q9a**3*s11*s15*s16**2*s22**2*s25*s26*s36 + 2048*q9a**3*s11**2*s12*s16*s22*s24*s25*s26*s36 +  \
    2048*q9a**3*s11**2*s12**2*s22*s25**2*s26*s36 - 1024*q9a**3*s11*s12*s16**2*s22*s25**2*s26*s36 - 1024*q9a**3*s11**2*s14*s15*s22**2*s26**2*s36 + 1024*q9a**3*s11**2*s13*s16*s22**2*s26**2*s36 -  \
    3072*q9a**3*s11**2*s12*s16*s22*s23*s26**2*s36 + 1536*q9a**3*s11**2*s12*s15*s22*s24*s26**2*s36 - 512*q9a**3*s11**2*s14*s16*s22*s24*s26**2*s36 + 512*q9a**3*s11**2*s12*s16*s24**2*s26**2*s36 +  \
    512*q9a**3*s11**2*s12*s14*s22*s25*s26**2*s36 + 1024*q9a**3*s11*s12*s15*s16*s22*s25*s26**2*s36 - 2048*q9a**3*s11**2*s12**2*s24*s25*s26**2*s36 + 512*q9a**3*s11*s12**2*s16*s25**2*s26**2*s36 - 1024*q9a**3*s11**2*s12*s13*s22*s26**3*s36 +  \
    512*q9a**3*s11**2*s14**2*s22*s26**3*s36 + 2048*q9a**3*s11**2*s12**2*s23*s26**3*s36 - 512*q9a**3*s11**2*s12*s14*s24*s26**3*s36 - 512*q9a**3*s11*s12**2*s15*s25*s26**3*s36 - 3072*q8a*q9a**2*s11**3*s16*s22**3*s33*s36 +  \
    6144*q8a*q9a**2*s11**3*s12*s22**2*s26*s33*s36 + 3072*q8a*q9a**2*s11**2*s16**2*s22**2*s26*s33*s36 - 3072*q7a*q9a**2*s11**3*s22**3*s26*s33*s36 - 9216*q8a*q9a**2*s11**2*s12*s16*s22*s26**2*s33*s36 +  \
    3072*q7a*q9a**2*s11**2*s16*s22**2*s26**2*s33*s36 + 6144*q8a*q9a**2*s11**2*s12**2*s26**3*s33*s36 - 3072*q7a*q9a**2*s11**2*s12*s22*s26**3*s33*s36 + 1024*q8a*q9a**2*s11**3*s15*s22**3*s34*s36 +  \
    512*q8a*q9a**2*s11**3*s16*s22**2*s24*s34*s36 - 3072*q8a*q9a**2*s11**3*s12*s22**2*s25*s34*s36 + 512*q8a*q9a**2*s11**2*s16**2*s22**2*s25*s34*s36 + 2048*q7a*q9a**2*s11**3*s22**3*s25*s34*s36 -  \
    2560*q8a*q9a**2*s11**3*s14*s22**2*s26*s34*s36 - 512*q8a*q9a**2*s11**2*s15*s16*s22**2*s26*s34*s36 + 2048*q8a*q9a**2*s11**3*s12*s22*s24*s26*s34*s36 - 1024*q8a*q9a**2*s11**2*s16**2*s22*s24*s26*s34*s36 +  \
    1024*q8a*q9a**2*s11**2*s12*s16*s22*s25*s26*s34*s36 - 1536*q7a*q9a**2*s11**2*s16*s22**2*s25*s26*s34*s36 + 2048*q8a*q9a**2*s11**2*s14*s16*s22*s26**2*s34*s36 + 512*q7a*q9a**2*s11**2*s15*s22**2*s26**2*s34*s36 +  \
    512*q8a*q9a**2*s11**2*s12*s16*s24*s26**2*s34*s36 - 512*q7a*q9a**2*s11**2*s16*s22*s24*s26**2*s34*s36 - 1024*q8a*q9a**2*s11**2*s12**2*s25*s26**2*s34*s36 + 1024*q7a*q9a**2*s11**2*s12*s22*s25*s26**2*s34*s36 -  \
    1536*q8a*q9a**2*s11**2*s12*s14*s26**3*s34*s36 - 512*q7a*q9a**2*s11**2*s14*s22*s26**3*s34*s36 + 1024*q7a*q9a**2*s11**2*s12*s24*s26**3*s34*s36 - 1024*q8a**2*q9a*s11**3*s16*s22**2*s34**2*s36 +  \
    1024*q8a**2*q9a*s11**2*s16**2*s22*s26*s34**2*s36 + 1024*q7a*q8a*q9a*s11**3*s22**2*s26*s34**2*s36 - 1024*q8a**2*q9a*s11**2*s12*s16*s26**2*s34**2*s36 - 1024*q7a*q8a*q9a*s11**2*s16*s22*s26**2*s34**2*s36 +  \
    1024*q7a*q8a*q9a*s11**2*s12*s26**3*s34**2*s36 + 2048*q8a*q9a**2*s11**3*s14*s22**3*s35*s36 - 3072*q8a*q9a**2*s11**3*s12*s22**2*s24*s35*s36 + 512*q8a*q9a**2*s11**2*s16**2*s22**2*s24*s35*s36 +  \
    1024*q7a*q9a**2*s11**3*s22**3*s24*s35*s36 + 2560*q8a*q9a**2*s11**2*s12*s16*s22**2*s25*s35*s36 - 512*q8a*q9a**2*s11*s16**3*s22**2*s25*s35*s36 - 2560*q7a*q9a**2*s11**2*s16*s22**3*s25*s35*s36 -  \
    512*q8a*q9a**2*s11**2*s12*s15*s22**2*s26*s35*s36 - 1536*q8a*q9a**2*s11**2*s14*s16*s22**2*s26*s35*s36 - 512*q8a*q9a**2*s11*s15*s16**2*s22**2*s26*s35*s36 + 512*q7a*q9a**2*s11**2*s15*s22**3*s26*s35*s36 +  \
    1024*q8a*q9a**2*s11**2*s12*s16*s22*s24*s26*s35*s36 - 512*q7a*q9a**2*s11**2*s16*s22**2*s24*s26*s35*s36 - 2048*q8a*q9a**2*s11**2*s12**2*s22*s25*s26*s35*s36 + 2048*q7a*q9a**2*s11**2*s12*s22**2*s25*s26*s35*s36 +  \
    2048*q7a*q9a**2*s11*s16**2*s22**2*s25*s26*s35*s36 + 1024*q8a*q9a**2*s11**2*s12*s14*s22*s26**2*s35*s36 + 2048*q8a*q9a**2*s11*s12*s15*s16*s22*s26**2*s35*s36 + 512*q7a*q9a**2*s11**2*s14*s22**2*s26**2*s35*s36 -  \
    1024*q7a*q9a**2*s11*s15*s16*s22**2*s26**2*s35*s36 - 1024*q8a*q9a**2*s11**2*s12**2*s24*s26**2*s35*s36 + 512*q8a*q9a**2*s11*s12**2*s16*s25*s26**2*s35*s36 - 3072*q7a*q9a**2*s11*s12*s16*s22*s25*s26**2*s35*s36 -  \
    1536*q8a*q9a**2*s11*s12**2*s15*s26**3*s35*s36 + 1024*q7a*q9a**2*s11*s12*s15*s22*s26**3*s35*s36 + 1024*q7a*q9a**2*s11*s12**2*s25*s26**3*s35*s36 + 4096*q8a**2*q9a*s11**3*s12*s22**2*s34*s35*s36 -  \
    4096*q7a*q8a*q9a*s11**3*s22**3*s34*s35*s36 - 4096*q8a**2*q9a*s11**2*s12*s16*s22*s26*s34*s35*s36 + 4096*q7a*q8a*q9a*s11**2*s16*s22**2*s26*s34*s35*s36 + 4096*q8a**2*q9a*s11**2*s12**2*s26**2*s34*s35*s36 -  \
    4096*q7a*q8a*q9a*s11**2*s12*s22*s26**2*s34*s35*s36 - 1024*q8a**2*q9a*s11**2*s12*s16*s22**2*s35**2*s36 + 1024*q7a*q8a*q9a*s11**2*s16*s22**3*s35**2*s36 + 1024*q8a**2*q9a*s11*s12*s16**2*s22*s26*s35**2*s36 +  \
    1024*q7a*q8a*q9a*s11**2*s12*s22**2*s26*s35**2*s36 - 1024*q7a*q8a*q9a*s11*s16**2*s22**2*s26*s35**2*s36 - 1024*q7a**2*q9a*s11**2*s22**3*s26*s35**2*s36 - 1024*q8a**2*q9a*s11*s12**2*s16*s26**2*s35**2*s36 +  \
    1024*q7a**2*q9a*s11*s16*s22**2*s26**2*s35**2*s36 + 1024*q7a*q8a*q9a*s11*s12**2*s26**3*s35**2*s36 - 1024*q7a**2*q9a*s11*s12*s22*s26**3*s35**2*s36 + 1024*q8a*q9a**2*s11**3*s13*s22**3*s36**2 - 512*q8a*q9a**2*s11**2*s15**2*s22**3*s36**2 -  \
    2048*q8a*q9a**2*s11**3*s12*s22**2*s23*s36**2 - 1024*q8a*q9a**2*s11**2*s16**2*s22**2*s23*s36**2 + 1024*q7a*q9a**2*s11**3*s22**3*s23*s36**2 - 1536*q8a*q9a**2*s11**3*s14*s22**2*s24*s36**2 + 512*q8a*q9a**2*s11**2*s15*s16*s22**2*s24*s36**2 +  \
    2048*q8a*q9a**2*s11**3*s12*s22*s24**2*s36**2 - 512*q7a*q9a**2*s11**3*s22**2*s24**2*s36**2 + 2560*q8a*q9a**2*s11**2*s12*s15*s22**2*s25*s36**2 - 512*q8a*q9a**2*s11**2*s14*s16*s22**2*s25*s36**2 +  \
    512*q8a*q9a**2*s11*s15*s16**2*s22**2*s25*s36**2 - 1536*q7a*q9a**2*s11**2*s15*s22**3*s25*s36**2 - 2048*q8a*q9a**2*s11**2*s12*s16*s22*s24*s25*s36**2 + 2048*q7a*q9a**2*s11**2*s16*s22**2*s24*s25*s36**2 -  \
    2048*q8a*q9a**2*s11**2*s12**2*s22*s25**2*s36**2 + 1024*q8a*q9a**2*s11*s12*s16**2*s22*s25**2*s36**2 + 1536*q7a*q9a**2*s11**2*s12*s22**2*s25**2*s36**2 - 1536*q7a*q9a**2*s11*s16**2*s22**2*s25**2*s36**2 +  \
    2048*q8a*q9a**2*s11**2*s14*s15*s22**2*s26*s36**2 - 2048*q8a*q9a**2*s11**2*s13*s16*s22**2*s26*s36**2 + 6144*q8a*q9a**2*s11**2*s12*s16*s22*s23*s26*s36**2 - 2048*q7a*q9a**2*s11**2*s16*s22**2*s23*s26*s36**2 -  \
    3072*q8a*q9a**2*s11**2*s12*s15*s22*s24*s26*s36**2 + 1024*q8a*q9a**2*s11**2*s14*s16*s22*s24*s26*s36**2 + 512*q7a*q9a**2*s11**2*s15*s22**2*s24*s26*s36**2 - 1024*q8a*q9a**2*s11**2*s12*s16*s24**2*s26*s36**2 -  \
    1024*q8a*q9a**2*s11**2*s12*s14*s22*s25*s26*s36**2 - 2048*q8a*q9a**2*s11*s12*s15*s16*s22*s25*s26*s36**2 - 512*q7a*q9a**2*s11**2*s14*s22**2*s25*s26*s36**2 + 1024*q7a*q9a**2*s11*s15*s16*s22**2*s25*s26*s36**2 +  \
    4096*q8a*q9a**2*s11**2*s12**2*s24*s25*s26*s36**2 - 2048*q7a*q9a**2*s11**2*s12*s22*s24*s25*s26*s36**2 - 1024*q8a*q9a**2*s11*s12**2*s16*s25**2*s26*s36**2 + 2048*q7a*q9a**2*s11*s12*s16*s22*s25**2*s26*s36**2 +  \
    3072*q8a*q9a**2*s11**2*s12*s13*s22*s26**2*s36**2 - 1536*q8a*q9a**2*s11**2*s14**2*s22*s26**2*s36**2 - 1024*q7a*q9a**2*s11**2*s13*s22**2*s26**2*s36**2 - 6144*q8a*q9a**2*s11**2*s12**2*s23*s26**2*s36**2 +  \
    3072*q7a*q9a**2*s11**2*s12*s22*s23*s26**2*s36**2 + 1536*q8a*q9a**2*s11**2*s12*s14*s24*s26**2*s36**2 + 512*q7a*q9a**2*s11**2*s14*s22*s24*s26**2*s36**2 - 512*q7a*q9a**2*s11**2*s12*s24**2*s26**2*s36**2 +  \
    1536*q8a*q9a**2*s11*s12**2*s15*s25*s26**2*s36**2 - 1024*q7a*q9a**2*s11*s12*s15*s22*s25*s26**2*s36**2 - 512*q7a*q9a**2*s11*s12**2*s25**2*s26**2*s36**2 - 2048*q8a**2*q9a*s11**3*s12*s22**2*s33*s36**2 -  \
    1024*q8a**2*q9a*s11**2*s16**2*s22**2*s33*s36**2 + 2048*q7a*q8a*q9a*s11**3*s22**3*s33*s36**2 + 6144*q8a**2*q9a*s11**2*s12*s16*s22*s26*s33*s36**2 - 4096*q7a*q8a*q9a*s11**2*s16*s22**2*s26*s33*s36**2 -  \
    6144*q8a**2*q9a*s11**2*s12**2*s26**2*s33*s36**2 + 6144*q7a*q8a*q9a*s11**2*s12*s22*s26**2*s33*s36**2 - 1024*q7a**2*q9a*s11**2*s22**2*s26**2*s33*s36**2 + 2048*q8a**2*q9a*s11**3*s14*s22**2*s34*s36**2 -  \
    3072*q8a**2*q9a*s11**3*s12*s22*s24*s34*s36**2 + 512*q8a**2*q9a*s11**2*s16**2*s22*s24*s34*s36**2 + 1024*q7a*q8a*q9a*s11**3*s22**2*s24*s34*s36**2 + 512*q8a**2*q9a*s11**2*s12*s16*s22*s25*s34*s36**2 -  \
    512*q7a*q8a*q9a*s11**2*s16*s22**2*s25*s34*s36**2 + 1536*q8a**2*q9a*s11**2*s12*s15*s22*s26*s34*s36**2 - 2560*q8a**2*q9a*s11**2*s14*s16*s22*s26*s34*s36**2 - 1536*q7a*q8a*q9a*s11**2*s15*s22**2*s26*s34*s36**2;

    v1_1=  \
    512*q8a**2*q9a*s11**2*s12*s16*s24*s26*s34*s36**2 + 1024*q7a*q8a*q9a*s11**2*s16*s22*s24*s26*s34*s36**2 - 1024*q8a**2*q9a*s11**2*s12**2*s25*s26*s34*s36**2 + 1024*q7a**2*q9a*s11**2*s22**2*s25*s26*s34*s36**2 +  \
    1536*q8a**2*q9a*s11**2*s12*s14*s26**2*s34*s36**2 + 1024*q7a*q8a*q9a*s11**2*s14*s22*s26**2*s34*s36**2 - 2048*q7a*q8a*q9a*s11**2*s12*s24*s26**2*s34*s36**2 + 1024*q8a**3*s11**3*s12*s22*s34**2*s36**2 -  \
    512*q8a**3*s11**2*s16**2*s22*s34**2*s36**2 - 1024*q7a*q8a**2*s11**3*s22**2*s34**2*s36**2 + 512*q8a**3*s11**2*s12*s16*s26*s34**2*s36**2 + 512*q7a*q8a**2*s11**2*s16*s22*s26*s34**2*s36**2 - 512*q7a*q8a**2*s11**2*s12*s26**2*s34**2*s36**2 -  \
    1024*q8a**2*q9a*s11**2*s12*s15*s22**2*s35*s36**2 + 1024*q8a**2*q9a*s11**2*s14*s16*s22**2*s35*s36**2 + 1024*q7a*q8a*q9a*s11**2*s15*s22**3*s35*s36**2 + 512*q8a**2*q9a*s11**2*s12*s16*s22*s24*s35*s36**2 -  \
    1536*q7a*q8a*q9a*s11**2*s16*s22**2*s24*s35*s36**2 + 3072*q8a**2*q9a*s11**2*s12**2*s22*s25*s35*s36**2 - 1024*q8a**2*q9a*s11*s12*s16**2*s22*s25*s35*s36**2 - 5120*q7a*q8a*q9a*s11**2*s12*s22**2*s25*s35*s36**2 +  \
    1024*q7a*q8a*q9a*s11*s16**2*s22**2*s25*s35*s36**2 + 2048*q7a**2*q9a*s11**2*s22**3*s25*s35*s36**2 - 512*q8a**2*q9a*s11**2*s12*s14*s22*s26*s35*s36**2 - 1024*q8a**2*q9a*s11*s12*s15*s16*s22*s26*s35*s36**2 -  \
    512*q7a*q8a*q9a*s11**2*s14*s22**2*s26*s35*s36**2 + 1024*q7a*q8a*q9a*s11*s15*s16*s22**2*s26*s35*s36**2 - 1024*q8a**2*q9a*s11**2*s12**2*s24*s26*s35*s36**2 + 2048*q7a*q8a*q9a*s11**2*s12*s22*s24*s26*s35*s36**2 +  \
    512*q8a**2*q9a*s11*s12**2*s16*s25*s26*s35*s36**2 + 2048*q7a*q8a*q9a*s11*s12*s16*s22*s25*s26*s35*s36**2 - 2560*q7a**2*q9a*s11*s16*s22**2*s25*s26*s35*s36**2 + 1536*q8a**2*q9a*s11*s12**2*s15*s26**2*s35*s36**2 -  \
    2048*q7a*q8a*q9a*s11*s12*s15*s22*s26**2*s35*s36**2 + 512*q7a**2*q9a*s11*s15*s22**2*s26**2*s35*s36**2 - 2048*q7a*q8a*q9a*s11*s12**2*s25*s26**2*s35*s36**2 + 2048*q7a**2*q9a*s11*s12*s22*s25*s26**2*s35*s36**2 +  \
    1024*q8a**3*s11**2*s12*s16*s22*s34*s35*s36**2 - 1024*q7a*q8a**2*s11**2*s16*s22**2*s34*s35*s36**2 - 2048*q8a**3*s11**2*s12**2*s26*s34*s35*s36**2 + 3072*q7a*q8a**2*s11**2*s12*s22*s26*s34*s35*s36**2 -  \
    1024*q7a**2*q8a*s11**2*s22**2*s26*s34*s35*s36**2 - 1024*q8a**3*s11**2*s12**2*s22*s35**2*s36**2 + 2048*q7a*q8a**2*s11**2*s12*s22**2*s35**2*s36**2 - 1024*q7a**2*q8a*s11**2*s22**3*s35**2*s36**2 + 512*q8a**3*s11*s12**2*s16*s26*s35**2*s36**2 -  \
    1024*q7a*q8a**2*s11*s12*s16*s22*s26*s35**2*s36**2 + 512*q7a**2*q8a*s11*s16*s22**2*s26*s35**2*s36**2 - 512*q7a*q8a**2*s11*s12**2*s26**2*s35**2*s36**2 + 1024*q7a**2*q8a*s11*s12*s22*s26**2*s35**2*s36**2 -  \
    512*q7a**3*s11*s22**2*s26**2*s35**2*s36**2 - 1024*q8a**2*q9a*s11**2*s14*s15*s22**2*s36**3 + 1024*q8a**2*q9a*s11**2*s13*s16*s22**2*s36**3 - 3072*q8a**2*q9a*s11**2*s12*s16*s22*s23*s36**3 + 2048*q7a*q8a*q9a*s11**2*s16*s22**2*s23*s36**3 +  \
    1536*q8a**2*q9a*s11**2*s12*s15*s22*s24*s36**3 - 512*q8a**2*q9a*s11**2*s14*s16*s22*s24*s36**3 - 512*q7a*q8a*q9a*s11**2*s15*s22**2*s24*s36**3 + 512*q8a**2*q9a*s11**2*s12*s16*s24**2*s36**3 +  \
    512*q8a**2*q9a*s11**2*s12*s14*s22*s25*s36**3 + 1024*q8a**2*q9a*s11*s12*s15*s16*s22*s25*s36**3 + 512*q7a*q8a*q9a*s11**2*s14*s22**2*s25*s36**3 - 1024*q7a*q8a*q9a*s11*s15*s16*s22**2*s25*s36**3 -  \
    2048*q8a**2*q9a*s11**2*s12**2*s24*s25*s36**3 + 2048*q7a*q8a*q9a*s11**2*s12*s22*s24*s25*s36**3 - 1024*q7a**2*q9a*s11**2*s22**2*s24*s25*s36**3 + 512*q8a**2*q9a*s11*s12**2*s16*s25**2*s36**3 -  \
    2048*q7a*q8a*q9a*s11*s12*s16*s22*s25**2*s36**3 + 1536*q7a**2*q9a*s11*s16*s22**2*s25**2*s36**3 - 3072*q8a**2*q9a*s11**2*s12*s13*s22*s26*s36**3 + 1536*q8a**2*q9a*s11**2*s14**2*s22*s26*s36**3 +  \
    2048*q7a*q8a*q9a*s11**2*s13*s22**2*s26*s36**3 + 6144*q8a**2*q9a*s11**2*s12**2*s23*s26*s36**3 - 6144*q7a*q8a*q9a*s11**2*s12*s22*s23*s26*s36**3 + 1024*q7a**2*q9a*s11**2*s22**2*s23*s26*s36**3 -  \
    1536*q8a**2*q9a*s11**2*s12*s14*s24*s26*s36**3 - 1024*q7a*q8a*q9a*s11**2*s14*s22*s24*s26*s36**3 + 1024*q7a*q8a*q9a*s11**2*s12*s24**2*s26*s36**3 - 1536*q8a**2*q9a*s11*s12**2*s15*s25*s26*s36**3 +  \
    2048*q7a*q8a*q9a*s11*s12*s15*s22*s25*s26*s36**3 - 512*q7a**2*q9a*s11*s15*s22**2*s25*s26*s36**3 + 1024*q7a*q8a*q9a*s11*s12**2*s25**2*s26*s36**3 - 1024*q7a**2*q9a*s11*s12*s22*s25**2*s26*s36**3 -  \
    1024*q8a**3*s11**2*s12*s16*s22*s33*s36**3 + 1024*q7a*q8a**2*s11**2*s16*s22**2*s33*s36**3 + 2048*q8a**3*s11**2*s12**2*s26*s33*s36**3 - 3072*q7a*q8a**2*s11**2*s12*s22*s26*s33*s36**3 + 1024*q7a**2*q8a*s11**2*s22**2*s26*s33*s36**3 -  \
    1024*q8a**3*s11**2*s12*s15*s22*s34*s36**3 + 1024*q8a**3*s11**2*s14*s16*s22*s34*s36**3 + 1024*q7a*q8a**2*s11**2*s15*s22**2*s34*s36**3 - 512*q8a**3*s11**2*s12*s16*s24*s34*s36**3 - 512*q7a*q8a**2*s11**2*s16*s22*s24*s34*s36**3 +  \
    1024*q8a**3*s11**2*s12**2*s25*s34*s36**3 - 1024*q7a*q8a**2*s11**2*s12*s22*s25*s34*s36**3 - 512*q8a**3*s11**2*s12*s14*s26*s34*s36**3 - 512*q7a*q8a**2*s11**2*s14*s22*s26*s34*s36**3 + 1024*q7a*q8a**2*s11**2*s12*s24*s26*s34*s36**3 +  \
    1024*q8a**3*s11**2*s12**2*s24*s35*s36**3 - 2048*q7a*q8a**2*s11**2*s12*s22*s24*s35*s36**3 + 1024*q7a**2*q8a*s11**2*s22**2*s24*s35*s36**3 - 512*q8a**3*s11*s12**2*s16*s25*s35*s36**3 + 1024*q7a*q8a**2*s11*s12*s16*s22*s25*s35*s36**3 -  \
    512*q7a**2*q8a*s11*s16*s22**2*s25*s35*s36**3 - 512*q8a**3*s11*s12**2*s15*s26*s35*s36**3 + 1024*q7a*q8a**2*s11*s12*s15*s22*s26*s35*s36**3 - 512*q7a**2*q8a*s11*s15*s22**2*s26*s35*s36**3 +  \
    1024*q7a*q8a**2*s11*s12**2*s25*s26*s35*s36**3 - 2048*q7a**2*q8a*s11*s12*s22*s25*s26*s35*s36**3 + 1024*q7a**3*s11*s22**2*s25*s26*s35*s36**3 + 1024*q8a**3*s11**2*s12*s13*s22*s36**4 - 512*q8a**3*s11**2*s14**2*s22*s36**4 -  \
    1024*q7a*q8a**2*s11**2*s13*s22**2*s36**4 - 2048*q8a**3*s11**2*s12**2*s23*s36**4 + 3072*q7a*q8a**2*s11**2*s12*s22*s23*s36**4 - 1024*q7a**2*q8a*s11**2*s22**2*s23*s36**4 + 512*q8a**3*s11**2*s12*s14*s24*s36**4 +  \
    512*q7a*q8a**2*s11**2*s14*s22*s24*s36**4 - 512*q7a*q8a**2*s11**2*s12*s24**2*s36**4 + 512*q8a**3*s11*s12**2*s15*s25*s36**4 - 1024*q7a*q8a**2*s11*s12*s15*s22*s25*s36**4 + 512*q7a**2*q8a*s11*s15*s22**2*s25*s36**4 -  \
    512*q7a*q8a**2*s11*s12**2*s25**2*s36**4 + 1024*q7a**2*q8a*s11*s12*s22*s25**2*s36**4 - 512*q7a**3*s11*s22**2*s25**2*s36**4;
    v1=v1_0+v1_1

    v2_0=3072*q9a**2*s11**4*s22**4*s33**2 - 6144*q9a**2*s11**3*s16*s22**3*s26*s33**2 + 6144*q9a**2*s11**3*s12*s22**2*s26**2*s33**2 + 3072*q9a**2*s11**2*s16**2*s22**2*s26**2*s33**2 - 6144*q9a**2*s11**2*s12*s16*s22*s26**3*s33**2 +  \
    3072*q9a**2*s11**2*s12**2*s26**4*s33**2 - 3072*q9a**2*s11**4*s22**3*s24*s33*s34 + 1536*q9a**2*s11**3*s16*s22**3*s25*s33*s34 + 1536*q9a**2*s11**3*s15*s22**3*s26*s33*s34 + 4608*q9a**2*s11**3*s16*s22**2*s24*s26*s33*s34 -  \
    3072*q9a**2*s11**3*s12*s22**2*s25*s26*s33*s34 - 1536*q9a**2*s11**2*s16**2*s22**2*s25*s26*s33*s34 - 1536*q9a**2*s11**3*s14*s22**2*s26**2*s33*s34 - 1536*q9a**2*s11**2*s15*s16*s22**2*s26**2*s33*s34 -  \
    3072*q9a**2*s11**3*s12*s22*s24*s26**2*s33*s34 - 1536*q9a**2*s11**2*s16**2*s22*s24*s26**2*s33*s34 + 4608*q9a**2*s11**2*s12*s16*s22*s25*s26**2*s33*s34 + 1536*q9a**2*s11**2*s12*s15*s22*s26**3*s33*s34 +  \
    1536*q9a**2*s11**2*s14*s16*s22*s26**3*s33*s34 + 1536*q9a**2*s11**2*s12*s16*s24*s26**3*s33*s34 - 3072*q9a**2*s11**2*s12**2*s25*s26**3*s33*s34 - 1536*q9a**2*s11**2*s12*s14*s26**4*s33*s34 + 1024*q9a**2*s11**4*s22**3*s23*s34**2 +  \
    512*q9a**2*s11**4*s22**2*s24**2*s34**2 - 512*q9a**2*s11**3*s15*s22**3*s25*s34**2 - 512*q9a**2*s11**3*s16*s22**2*s24*s25*s34**2 + 512*q9a**2*s11**3*s12*s22**2*s25**2*s34**2 - 1536*q9a**2*s11**3*s16*s22**2*s23*s26*s34**2 -  \
    512*q9a**2*s11**3*s15*s22**2*s24*s26*s34**2 - 512*q9a**2*s11**3*s16*s22*s24**2*s26*s34**2 + 1024*q9a**2*s11**3*s14*s22**2*s25*s26*s34**2 + 512*q9a**2*s11**2*s15*s16*s22**2*s25*s26*s34**2 +  \
    512*q9a**2*s11**2*s16**2*s22*s24*s25*s26*s34**2 - 512*q9a**2*s11**2*s12*s16*s22*s25**2*s26*s34**2 + 512*q9a**2*s11**3*s13*s22**2*s26**2*s34**2 + 1024*q9a**2*s11**3*s12*s22*s23*s26**2*s34**2 +  \
    512*q9a**2*s11**2*s16**2*s22*s23*s26**2*s34**2 + 512*q9a**2*s11**2*s15*s16*s22*s24*s26**2*s34**2 + 512*q9a**2*s11**3*s12*s24**2*s26**2*s34**2 - 512*q9a**2*s11**2*s12*s15*s22*s25*s26**2*s34**2 -  \
    1024*q9a**2*s11**2*s14*s16*s22*s25*s26**2*s34**2 - 512*q9a**2*s11**2*s12*s16*s24*s25*s26**2*s34**2 + 512*q9a**2*s11**2*s12**2*s25**2*s26**2*s34**2 - 512*q9a**2*s11**2*s13*s16*s22*s26**3*s34**2 -  \
    512*q9a**2*s11**2*s12*s16*s23*s26**3*s34**2 - 512*q9a**2*s11**2*s12*s15*s24*s26**3*s34**2 + 1024*q9a**2*s11**2*s12*s14*s25*s26**3*s34**2 + 512*q9a**2*s11**2*s12*s13*s26**4*s34**2 + 2048*q8a*q9a*s11**4*s22**3*s33*s34**2 -  \
    3072*q8a*q9a*s11**3*s16*s22**2*s26*s33*s34**2 + 2048*q8a*q9a*s11**3*s12*s22*s26**2*s33*s34**2 + 1024*q8a*q9a*s11**2*s16**2*s22*s26**2*s33*s34**2 + 1024*q7a*q9a*s11**3*s22**2*s26**2*s33*s34**2 -  \
    1024*q8a*q9a*s11**2*s12*s16*s26**3*s33*s34**2 - 1024*q7a*q9a*s11**2*s16*s22*s26**3*s33*s34**2 + 1024*q7a*q9a*s11**2*s12*s26**4*s33*s34**2 - 1024*q8a*q9a*s11**4*s22**2*s24*s34**3 + 512*q8a*q9a*s11**3*s16*s22**2*s25*s34**3 +  \
    512*q8a*q9a*s11**3*s15*s22**2*s26*s34**3 + 1024*q8a*q9a*s11**3*s16*s22*s24*s26*s34**3 - 512*q8a*q9a*s11**2*s16**2*s22*s25*s26*s34**3 - 1024*q7a*q9a*s11**3*s22**2*s25*s26*s34**3 - 512*q8a*q9a*s11**2*s15*s16*s22*s26**2*s34**3 -  \
    1024*q8a*q9a*s11**3*s12*s24*s26**2*s34**3 + 512*q8a*q9a*s11**2*s12*s16*s25*s26**2*s34**3 + 1024*q7a*q9a*s11**2*s16*s22*s25*s26**2*s34**3 + 512*q8a*q9a*s11**2*s12*s15*s26**3*s34**3 - 1024*q7a*q9a*s11**2*s12*s25*s26**3*s34**3 +  \
    512*q8a**2*s11**4*s22**2*s34**4 - 512*q8a**2*s11**3*s16*s22*s26*s34**4 + 512*q8a**2*s11**3*s12*s26**2*s34**4 - 3072*q9a**2*s11**3*s15*s22**4*s33*s35 + 1536*q9a**2*s11**3*s16*s22**3*s24*s33*s35 +  \
    3072*q9a**2*s11**3*s12*s22**3*s25*s33*s35 - 1536*q9a**2*s11**2*s16**2*s22**3*s25*s33*s35 + 1536*q9a**2*s11**3*s14*s22**3*s26*s33*s35 + 4608*q9a**2*s11**2*s15*s16*s22**3*s26*s33*s35 - 3072*q9a**2*s11**3*s12*s22**2*s24*s26*s33*s35 -  \
    1536*q9a**2*s11**2*s16**2*s22**2*s24*s26*s33*s35 - 1536*q9a**2*s11**2*s12*s16*s22**2*s25*s26*s33*s35 + 1536*q9a**2*s11*s16**3*s22**2*s25*s26*s33*s35 - 4608*q9a**2*s11**2*s12*s15*s22**2*s26**2*s33*s35 -  \
    1536*q9a**2*s11**2*s14*s16*s22**2*s26**2*s33*s35 - 1536*q9a**2*s11*s15*s16**2*s22**2*s26**2*s33*s35 + 4608*q9a**2*s11**2*s12*s16*s22*s24*s26**2*s33*s35 + 3072*q9a**2*s11**2*s12**2*s22*s25*s26**2*s33*s35 -  \
    3072*q9a**2*s11*s12*s16**2*s22*s25*s26**2*s33*s35 + 1536*q9a**2*s11**2*s12*s14*s22*s26**3*s33*s35 + 3072*q9a**2*s11*s12*s15*s16*s22*s26**3*s33*s35 - 3072*q9a**2*s11**2*s12**2*s24*s26**3*s33*s35 +  \
    1536*q9a**2*s11*s12**2*s16*s25*s26**3*s33*s35 - 1536*q9a**2*s11*s12**2*s15*s26**4*s33*s35 - 1024*q9a**2*s11**3*s16*s22**3*s23*s34*s35 + 1536*q9a**2*s11**3*s15*s22**3*s24*s34*s35 - 512*q9a**2*s11**3*s16*s22**2*s24**2*s34*s35 -  \
    1536*q9a**2*s11**3*s14*s22**3*s25*s34*s35 + 512*q9a**2*s11**2*s15*s16*s22**3*s25*s34*s35 + 512*q9a**2*s11**2*s16**2*s22**2*s24*s25*s34*s35 - 512*q9a**2*s11**2*s12*s16*s22**2*s25**2*s34*s35 -  \
    1024*q9a**2*s11**3*s13*s22**3*s26*s34*s35 - 512*q9a**2*s11**2*s15**2*s22**3*s26*s34*s35 + 2048*q9a**2*s11**3*s12*s22**2*s23*s26*s34*s35 + 1024*q9a**2*s11**2*s16**2*s22**2*s23*s26*s34*s35 +  \
    512*q9a**2*s11**3*s14*s22**2*s24*s26*s34*s35 - 2048*q9a**2*s11**2*s15*s16*s22**2*s24*s26*s34*s35 + 512*q9a**2*s11**2*s16**2*s22*s24**2*s26*s34*s35 + 512*q9a**2*s11**2*s12*s15*s22**2*s25*s26*s34*s35 +  \
    1024*q9a**2*s11**2*s14*s16*s22**2*s25*s26*s34*s35 - 512*q9a**2*s11*s15*s16**2*s22**2*s25*s26*s34*s35 - 512*q9a**2*s11*s16**3*s22*s24*s25*s26*s34*s35 + 512*q9a**2*s11*s12*s16**2*s22*s25**2*s26*s34*s35 +  \
    512*q9a**2*s11**2*s14*s15*s22**2*s26**2*s34*s35 + 1024*q9a**2*s11**2*s13*s16*s22**2*s26**2*s34*s35 + 512*q9a**2*s11*s15**2*s16*s22**2*s26**2*s34*s35 - 3072*q9a**2*s11**2*s12*s16*s22*s23*s26**2*s34*s35 +  \
    1536*q9a**2*s11**2*s12*s15*s22*s24*s26**2*s34*s35 - 512*q9a**2*s11**2*s14*s16*s22*s24*s26**2*s34*s35 + 512*q9a**2*s11*s15*s16**2*s22*s24*s26**2*s34*s35 - 512*q9a**2*s11**2*s12*s16*s24**2*s26**2*s34*s35 -  \
    1536*q9a**2*s11**2*s12*s14*s22*s25*s26**2*s34*s35 + 512*q9a**2*s11*s14*s16**2*s22*s25*s26**2*s34*s35 + 512*q9a**2*s11*s12*s16**2*s24*s25*s26**2*s34*s35 - 512*q9a**2*s11*s12**2*s16*s25**2*s26**2*s34*s35 -  \
    1024*q9a**2*s11**2*s12*s13*s22*s26**3*s34*s35 - 512*q9a**2*s11*s12*s15**2*s22*s26**3*s34*s35 - 512*q9a**2*s11*s14*s15*s16*s22*s26**3*s34*s35 + 2048*q9a**2*s11**2*s12**2*s23*s26**3*s34*s35 +  \
    512*q9a**2*s11**2*s12*s14*s24*s26**3*s34*s35 - 512*q9a**2*s11*s12*s15*s16*s24*s26**3*s34*s35 + 512*q9a**2*s11*s12**2*s15*s25*s26**3*s34*s35 - 512*q9a**2*s11*s12*s14*s16*s25*s26**3*s34*s35 +  \
    512*q9a**2*s11*s12*s14*s15*s26**4*s34*s35 - 2048*q8a*q9a*s11**3*s16*s22**3*s33*s34*s35 + 4096*q8a*q9a*s11**3*s12*s22**2*s26*s33*s34*s35 + 2048*q8a*q9a*s11**2*s16**2*s22**2*s26*s33*s34*s35 -  \
    2048*q7a*q9a*s11**3*s22**3*s26*s33*s34*s35 - 6144*q8a*q9a*s11**2*s12*s16*s22*s26**2*s33*s34*s35 + 2048*q7a*q9a*s11**2*s16*s22**2*s26**2*s33*s34*s35 + 4096*q8a*q9a*s11**2*s12**2*s26**3*s33*s34*s35 -  \
    2048*q7a*q9a*s11**2*s12*s22*s26**3*s33*s34*s35 - 1024*q8a*q9a*s11**3*s15*s22**3*s34**2*s35 + 1536*q8a*q9a*s11**3*s16*s22**2*s24*s34**2*s35 - 1024*q8a*q9a*s11**3*s12*s22**2*s25*s34**2*s35 -  \
    512*q8a*q9a*s11**2*s16**2*s22**2*s25*s34**2*s35 + 2048*q7a*q9a*s11**3*s22**3*s25*s34**2*s35 - 1536*q8a*q9a*s11**3*s14*s22**2*s26*s34**2*s35 + 1536*q8a*q9a*s11**2*s15*s16*s22**2*s26*s34**2*s35 -  \
    1536*q8a*q9a*s11**2*s16**2*s22*s24*s26*s34**2*s35 + 1024*q8a*q9a*s11**2*s12*s16*s22*s25*s26*s34**2*s35 + 512*q8a*q9a*s11*s16**3*s22*s25*s26*s34**2*s35 - 1536*q7a*q9a*s11**2*s16*s22**2*s25*s26*s34**2*s35 -  \
    1024*q8a*q9a*s11**2*s12*s15*s22*s26**2*s34**2*s35 + 1536*q8a*q9a*s11**2*s14*s16*s22*s26**2*s34**2*s35 - 512*q8a*q9a*s11*s15*s16**2*s22*s26**2*s34**2*s35 - 512*q7a*q9a*s11**2*s15*s22**2*s26**2*s34**2*s35 +  \
    1536*q8a*q9a*s11**2*s12*s16*s24*s26**2*s34**2*s35 - 1024*q8a*q9a*s11**2*s12**2*s25*s26**2*s34**2*s35 - 512*q8a*q9a*s11*s12*s16**2*s25*s26**2*s34**2*s35 + 2048*q7a*q9a*s11**2*s12*s22*s25*s26**2*s34**2*s35 -  \
    512*q7a*q9a*s11*s16**2*s22*s25*s26**2*s34**2*s35 - 1536*q8a*q9a*s11**2*s12*s14*s26**3*s34**2*s35 + 512*q8a*q9a*s11*s12*s15*s16*s26**3*s34**2*s35 + 512*q7a*q9a*s11*s15*s16*s22*s26**3*s34**2*s35 +  \
    512*q7a*q9a*s11*s12*s16*s25*s26**3*s34**2*s35 - 512*q7a*q9a*s11*s12*s15*s26**4*s34**2*s35 - 1024*q8a**2*s11**3*s16*s22**2*s34**3*s35 + 1024*q8a**2*s11**2*s16**2*s22*s26*s34**3*s35 + 1024*q7a*q8a*s11**3*s22**2*s26*s34**3*s35 -  \
    1024*q8a**2*s11**2*s12*s16*s26**2*s34**3*s35 - 1024*q7a*q8a*s11**2*s16*s22*s26**2*s34**3*s35 + 1024*q7a*q8a*s11**2*s12*s26**3*s34**3*s35 + 1024*q9a**2*s11**3*s13*s22**4*s35**2 + 512*q9a**2*s11**2*s15**2*s22**4*s35**2 -  \
    1024*q9a**2*s11**3*s12*s22**3*s23*s35**2 + 512*q9a**2*s11**2*s16**2*s22**3*s23*s35**2 - 512*q9a**2*s11**3*s14*s22**3*s24*s35**2 - 512*q9a**2*s11**2*s15*s16*s22**3*s24*s35**2 + 512*q9a**2*s11**3*s12*s22**2*s24**2*s35**2 -  \
    1024*q9a**2*s11**2*s12*s15*s22**3*s25*s35**2 + 1024*q9a**2*s11**2*s14*s16*s22**3*s25*s35**2 - 512*q9a**2*s11**2*s12*s16*s22**2*s24*s25*s35**2 + 512*q9a**2*s11**2*s12**2*s22**2*s25**2*s35**2 - 512*q9a**2*s11**2*s14*s15*s22**3*s26*s35**2 -  \
    1536*q9a**2*s11**2*s13*s16*s22**3*s26*s35**2 - 512*q9a**2*s11*s15**2*s16*s22**3*s26*s35**2 + 512*q9a**2*s11**2*s12*s16*s22**2*s23*s26*s35**2 - 512*q9a**2*s11*s16**3*s22**2*s23*s26*s35**2 +  \
    1024*q9a**2*s11**2*s12*s15*s22**2*s24*s26*s35**2 + 512*q9a**2*s11**2*s14*s16*s22**2*s24*s26*s35**2 + 512*q9a**2*s11*s15*s16**2*s22**2*s24*s26*s35**2 - 512*q9a**2*s11**2*s12*s16*s22*s24**2*s26*s35**2 -  \
    512*q9a**2*s11**2*s12*s14*s22**2*s25*s26*s35**2 + 1024*q9a**2*s11*s12*s15*s16*s22**2*s25*s26*s35**2 - 1024*q9a**2*s11*s14*s16**2*s22**2*s25*s26*s35**2 + 512*q9a**2*s11*s12*s16**2*s22*s24*s25*s26*s35**2 -  \
    512*q9a**2*s11*s12**2*s16*s22*s25**2*s26*s35**2 + 1536*q9a**2*s11**2*s12*s13*s22**2*s26**2*s35**2 + 512*q9a**2*s11*s12*s15**2*s22**2*s26**2*s35**2 + 512*q9a**2*s11*s14*s15*s16*s22**2*s26**2*s35**2 +  \
    512*q9a**2*s11*s13*s16**2*s22**2*s26**2*s35**2 - 1024*q9a**2*s11**2*s12**2*s22*s23*s26**2*s35**2 + 1024*q9a**2*s11*s12*s16**2*s22*s23*s26**2*s35**2 - 512*q9a**2*s11**2*s12*s14*s22*s24*s26**2*s35**2 -  \
    1536*q9a**2*s11*s12*s15*s16*s22*s24*s26**2*s35**2 + 512*q9a**2*s11**2*s12**2*s24**2*s26**2*s35**2 - 1024*q9a**2*s11*s12**2*s15*s22*s25*s26**2*s35**2 + 1536*q9a**2*s11*s12*s14*s16*s22*s25*s26**2*s35**2 -  \
    512*q9a**2*s11*s12**2*s16*s24*s25*s26**2*s35**2 + 512*q9a**2*s11*s12**3*s25**2*s26**2*s35**2 - 512*q9a**2*s11*s12*s14*s15*s22*s26**3*s35**2 - 1024*q9a**2*s11*s12*s13*s16*s22*s26**3*s35**2 - 512*q9a**2*s11*s12**2*s16*s23*s26**3*s35**2 +  \
    1024*q9a**2*s11*s12**2*s15*s24*s26**3*s35**2 - 512*q9a**2*s11*s12**2*s14*s25*s26**3*s35**2 + 512*q9a**2*s11*s12**2*s13*s26**4*s35**2 - 2048*q8a*q9a*s11**3*s12*s22**3*s33*s35**2 + 1024*q8a*q9a*s11**2*s16**2*s22**3*s33*s35**2 +  \
    2048*q7a*q9a*s11**3*s22**4*s33*s35**2 + 1024*q8a*q9a*s11**2*s12*s16*s22**2*s26*s33*s35**2 - 1024*q8a*q9a*s11*s16**3*s22**2*s26*s33*s35**2 - 3072*q7a*q9a*s11**2*s16*s22**3*s26*s33*s35**2 -  \
    2048*q8a*q9a*s11**2*s12**2*s22*s26**2*s33*s35**2 + 2048*q8a*q9a*s11*s12*s16**2*s22*s26**2*s33*s35**2 + 3072*q7a*q9a*s11**2*s12*s22**2*s26**2*s33*s35**2 + 1024*q7a*q9a*s11*s16**2*s22**2*s26**2*s33*s35**2 -  \
    1024*q8a*q9a*s11*s12**2*s16*s26**3*s33*s35**2 - 2048*q7a*q9a*s11*s12*s16*s22*s26**3*s33*s35**2 + 1024*q7a*q9a*s11*s12**2*s26**4*s33*s35**2 + 2048*q8a*q9a*s11**3*s14*s22**3*s34*s35**2 -  \
    1024*q8a*q9a*s11**3*s12*s22**2*s24*s34*s35**2 - 512*q8a*q9a*s11**2*s16**2*s22**2*s24*s34*s35**2 - 1024*q7a*q9a*s11**3*s22**3*s24*s34*s35**2 + 1536*q8a*q9a*s11**2*s12*s16*s22**2*s25*s34*s35**2 -  \
    1536*q7a*q9a*s11**2*s16*s22**3*s25*s34*s35**2 - 1536*q8a*q9a*s11**2*s12*s15*s22**2*s26*s34*s35**2 - 1536*q8a*q9a*s11**2*s14*s16*s22**2*s26*s34*s35**2 + 1536*q7a*q9a*s11**2*s15*s22**3*s26*s34*s35**2 +  \
    1024*q8a*q9a*s11**2*s12*s16*s22*s24*s26*s34*s35**2 + 512*q8a*q9a*s11*s16**3*s22*s24*s26*s34*s35**2 + 1536*q7a*q9a*s11**2*s16*s22**2*s24*s26*s34*s35**2 - 1536*q8a*q9a*s11*s12*s16**2*s22*s25*s26*s34*s35**2;
    v2_1= \
    1536*q7a*q9a*s11*s16**2*s22**2*s25*s26*s34*s35**2 + 2048*q8a*q9a*s11**2*s12*s14*s22*s26**2*s34*s35**2 + 1536*q8a*q9a*s11*s12*s15*s16*s22*s26**2*s34*s35**2 - 512*q8a*q9a*s11*s14*s16**2*s22*s26**2*s34*s35**2 -  \
    512*q7a*q9a*s11**2*s14*s22**2*s26**2*s34*s35**2 - 1536*q7a*q9a*s11*s15*s16*s22**2*s26**2*s34*s35**2 - 1024*q8a*q9a*s11**2*s12**2*s24*s26**2*s34*s35**2 - 512*q8a*q9a*s11*s12*s16**2*s24*s26**2*s34*s35**2 -  \
    1024*q7a*q9a*s11**2*s12*s22*s24*s26**2*s34*s35**2 - 512*q7a*q9a*s11*s16**2*s22*s24*s26**2*s34*s35**2 + 1536*q8a*q9a*s11*s12**2*s16*s25*s26**2*s34*s35**2 - 1536*q7a*q9a*s11*s12*s16*s22*s25*s26**2*s34*s35**2 -  \
    1536*q8a*q9a*s11*s12**2*s15*s26**3*s34*s35**2 + 512*q8a*q9a*s11*s12*s14*s16*s26**3*s34*s35**2 + 1536*q7a*q9a*s11*s12*s15*s22*s26**3*s34*s35**2 + 512*q7a*q9a*s11*s14*s16*s22*s26**3*s34*s35**2 +  \
    512*q7a*q9a*s11*s12*s16*s24*s26**3*s34*s35**2 - 512*q7a*q9a*s11*s12*s14*s26**4*s34*s35**2 + 1024*q8a**2*s11**3*s12*s22**2*s34**2*s35**2 + 512*q8a**2*s11**2*s16**2*s22**2*s34**2*s35**2 - 1024*q7a*q8a*s11**3*s22**3*s34**2*s35**2 -  \
    1024*q8a**2*s11**2*s12*s16*s22*s26*s34**2*s35**2 - 512*q8a**2*s11*s16**3*s22*s26*s34**2*s35**2 + 1024*q8a**2*s11**2*s12**2*s26**2*s34**2*s35**2 + 512*q8a**2*s11*s12*s16**2*s26**2*s34**2*s35**2 -  \
    1024*q7a*q8a*s11**2*s12*s22*s26**2*s34**2*s35**2 + 1024*q7a*q8a*s11*s16**2*s22*s26**2*s34**2*s35**2 + 512*q7a**2*s11**2*s22**2*s26**2*s34**2*s35**2 - 1024*q7a*q8a*s11*s12*s16*s26**3*s34**2*s35**2 -  \
    512*q7a**2*s11*s16*s22*s26**3*s34**2*s35**2 + 512*q7a**2*s11*s12*s26**4*s34**2*s35**2 + 1024*q8a*q9a*s11**2*s12*s15*s22**3*s35**3 - 1024*q8a*q9a*s11**2*s14*s16*s22**3*s35**3 - 1024*q7a*q9a*s11**2*s15*s22**4*s35**3 +  \
    512*q8a*q9a*s11**2*s12*s16*s22**2*s24*s35**3 + 512*q7a*q9a*s11**2*s16*s22**3*s24*s35**3 - 1024*q8a*q9a*s11**2*s12**2*s22**2*s25*s35**3 + 1024*q7a*q9a*s11**2*s12*s22**3*s25*s35**3 + 512*q8a*q9a*s11**2*s12*s14*s22**2*s26*s35**3 -  \
    1024*q8a*q9a*s11*s12*s15*s16*s22**2*s26*s35**3 + 1024*q8a*q9a*s11*s14*s16**2*s22**2*s26*s35**3 + 512*q7a*q9a*s11**2*s14*s22**3*s26*s35**3 + 1024*q7a*q9a*s11*s15*s16*s22**3*s26*s35**3 -  \
    512*q8a*q9a*s11*s12*s16**2*s22*s24*s26*s35**3 - 1024*q7a*q9a*s11**2*s12*s22**2*s24*s26*s35**3 - 512*q7a*q9a*s11*s16**2*s22**2*s24*s26*s35**3 + 1024*q8a*q9a*s11*s12**2*s16*s22*s25*s26*s35**3 -  \
    1024*q7a*q9a*s11*s12*s16*s22**2*s25*s26*s35**3 + 1024*q8a*q9a*s11*s12**2*s15*s22*s26**2*s35**3 - 1536*q8a*q9a*s11*s12*s14*s16*s22*s26**2*s35**3 - 1024*q7a*q9a*s11*s12*s15*s22**2*s26**2*s35**3 -  \
    512*q7a*q9a*s11*s14*s16*s22**2*s26**2*s35**3 + 512*q8a*q9a*s11*s12**2*s16*s24*s26**2*s35**3 + 1536*q7a*q9a*s11*s12*s16*s22*s24*s26**2*s35**3 - 1024*q8a*q9a*s11*s12**3*s25*s26**2*s35**3 +  \
    1024*q7a*q9a*s11*s12**2*s22*s25*s26**2*s35**3 + 512*q8a*q9a*s11*s12**2*s14*s26**3*s35**3 + 512*q7a*q9a*s11*s12*s14*s22*s26**3*s35**3 - 1024*q7a*q9a*s11*s12**2*s24*s26**3*s35**3 - 1024*q8a**2*s11**2*s12*s16*s22**2*s34*s35**3 +  \
    1024*q7a*q8a*s11**2*s16*s22**3*s34*s35**3 + 1024*q8a**2*s11*s12*s16**2*s22*s26*s34*s35**3 + 1024*q7a*q8a*s11**2*s12*s22**2*s26*s34*s35**3 - 1024*q7a*q8a*s11*s16**2*s22**2*s26*s34*s35**3 - 1024*q7a**2*s11**2*s22**3*s26*s34*s35**3 -  \
    1024*q8a**2*s11*s12**2*s16*s26**2*s34*s35**3 + 1024*q7a**2*s11*s16*s22**2*s26**2*s34*s35**3 + 1024*q7a*q8a*s11*s12**2*s26**3*s34*s35**3 - 1024*q7a**2*s11*s12*s22*s26**3*s34*s35**3 + 512*q8a**2*s11**2*s12**2*s22**2*s35**4 -  \
    1024*q7a*q8a*s11**2*s12*s22**3*s35**4 + 512*q7a**2*s11**2*s22**4*s35**4 - 512*q8a**2*s11*s12**2*s16*s22*s26*s35**4 + 1024*q7a*q8a*s11*s12*s16*s22**2*s26*s35**4 - 512*q7a**2*s11*s16*s22**3*s26*s35**4 +  \
    512*q8a**2*s11*s12**3*s26**2*s35**4 - 1024*q7a*q8a*s11*s12**2*s22*s26**2*s35**4 + 512*q7a**2*s11*s12*s22**2*s26**2*s35**4 + 3072*q9a**2*s11**3*s16*s22**3*s23*s33*s36 + 1536*q9a**2*s11**3*s15*s22**3*s24*s33*s36 -  \
    1536*q9a**2*s11**3*s16*s22**2*s24**2*s33*s36 + 4608*q9a**2*s11**3*s14*s22**3*s25*s33*s36 - 4608*q9a**2*s11**2*s15*s16*s22**3*s25*s33*s36 - 6144*q9a**2*s11**3*s12*s22**2*s24*s25*s33*s36 +  \
    3072*q9a**2*s11**2*s16**2*s22**2*s24*s25*s33*s36 + 4608*q9a**2*s11**2*s12*s16*s22**2*s25**2*s33*s36 - 1536*q9a**2*s11*s16**3*s22**2*s25**2*s33*s36 + 3072*q9a**2*s11**3*s13*s22**3*s26*s33*s36 -  \
    1536*q9a**2*s11**2*s15**2*s22**3*s26*s33*s36 - 6144*q9a**2*s11**3*s12*s22**2*s23*s26*s33*s36 - 3072*q9a**2*s11**2*s16**2*s22**2*s23*s26*s33*s36 - 4608*q9a**2*s11**3*s14*s22**2*s24*s26*s33*s36 +  \
    1536*q9a**2*s11**2*s15*s16*s22**2*s24*s26*s33*s36 + 6144*q9a**2*s11**3*s12*s22*s24**2*s26*s33*s36 + 7680*q9a**2*s11**2*s12*s15*s22**2*s25*s26*s33*s36 - 1536*q9a**2*s11**2*s14*s16*s22**2*s25*s26*s33*s36 +  \
    1536*q9a**2*s11*s15*s16**2*s22**2*s25*s26*s33*s36 - 6144*q9a**2*s11**2*s12*s16*s22*s24*s25*s26*s33*s36 - 6144*q9a**2*s11**2*s12**2*s22*s25**2*s26*s33*s36 + 3072*q9a**2*s11*s12*s16**2*s22*s25**2*s26*s33*s36 +  \
    3072*q9a**2*s11**2*s14*s15*s22**2*s26**2*s33*s36 - 3072*q9a**2*s11**2*s13*s16*s22**2*s26**2*s33*s36 + 9216*q9a**2*s11**2*s12*s16*s22*s23*s26**2*s33*s36 - 4608*q9a**2*s11**2*s12*s15*s22*s24*s26**2*s33*s36 +  \
    1536*q9a**2*s11**2*s14*s16*s22*s24*s26**2*s33*s36 - 1536*q9a**2*s11**2*s12*s16*s24**2*s26**2*s33*s36 - 1536*q9a**2*s11**2*s12*s14*s22*s25*s26**2*s33*s36 - 3072*q9a**2*s11*s12*s15*s16*s22*s25*s26**2*s33*s36 +  \
    6144*q9a**2*s11**2*s12**2*s24*s25*s26**2*s33*s36 - 1536*q9a**2*s11*s12**2*s16*s25**2*s26**2*s33*s36 + 3072*q9a**2*s11**2*s12*s13*s22*s26**3*s33*s36 - 1536*q9a**2*s11**2*s14**2*s22*s26**3*s33*s36 -  \
    6144*q9a**2*s11**2*s12**2*s23*s26**3*s33*s36 + 1536*q9a**2*s11**2*s12*s14*s24*s26**3*s33*s36 + 1536*q9a**2*s11*s12**2*s15*s25*s26**3*s33*s36 + 3072*q8a*q9a*s11**3*s16*s22**3*s33**2*s36 -  \
    6144*q8a*q9a*s11**3*s12*s22**2*s26*s33**2*s36 - 3072*q8a*q9a*s11**2*s16**2*s22**2*s26*s33**2*s36 + 3072*q7a*q9a*s11**3*s22**3*s26*s33**2*s36 + 9216*q8a*q9a*s11**2*s12*s16*s22*s26**2*s33**2*s36 -  \
    3072*q7a*q9a*s11**2*s16*s22**2*s26**2*s33**2*s36 - 6144*q8a*q9a*s11**2*s12**2*s26**3*s33**2*s36 + 3072*q7a*q9a*s11**2*s12*s22*s26**3*s33**2*s36 - 1024*q9a**2*s11**3*s15*s22**3*s23*s34*s36 -  \
    512*q9a**2*s11**3*s16*s22**2*s23*s24*s34*s36 - 512*q9a**2*s11**3*s15*s22**2*s24**2*s34*s36 + 512*q9a**2*s11**3*s16*s22*s24**3*s34*s36 - 2048*q9a**2*s11**3*s13*s22**3*s25*s34*s36 + 1024*q9a**2*s11**2*s15**2*s22**3*s25*s34*s36 +  \
    3072*q9a**2*s11**3*s12*s22**2*s23*s25*s34*s36 - 512*q9a**2*s11**2*s16**2*s22**2*s23*s25*s34*s36 - 512*q9a**2*s11**3*s14*s22**2*s24*s25*s34*s36 + 1024*q9a**2*s11**2*s15*s16*s22**2*s24*s25*s34*s36 +  \
    1024*q9a**2*s11**3*s12*s22*s24**2*s25*s34*s36 - 1024*q9a**2*s11**2*s16**2*s22*s24**2*s25*s34*s36 - 2048*q9a**2*s11**2*s12*s15*s22**2*s25**2*s34*s36 - 512*q9a**2*s11**2*s14*s16*s22**2*s25**2*s34*s36 +  \
    512*q9a**2*s11*s15*s16**2*s22**2*s25**2*s34*s36 - 512*q9a**2*s11**2*s12*s16*s22*s24*s25**2*s34*s36 + 512*q9a**2*s11*s16**3*s22*s24*s25**2*s34*s36 + 1024*q9a**2*s11**2*s12**2*s22*s25**3*s34*s36 -  \
    512*q9a**2*s11*s12*s16**2*s22*s25**3*s34*s36 + 2560*q9a**2*s11**3*s14*s22**2*s23*s26*s34*s36 + 512*q9a**2*s11**2*s15*s16*s22**2*s23*s26*s34*s36 + 512*q9a**2*s11**2*s15**2*s22**2*s24*s26*s34*s36 -  \
    2048*q9a**2*s11**3*s12*s22*s23*s24*s26*s34*s36 + 1024*q9a**2*s11**2*s16**2*s22*s23*s24*s26*s34*s36 + 512*q9a**2*s11**3*s14*s22*s24**2*s26*s34*s36 - 512*q9a**2*s11**2*s15*s16*s22*s24**2*s26*s34*s36 -  \
    1024*q9a**2*s11**3*s12*s24**3*s26*s34*s36 - 2048*q9a**2*s11**2*s14*s15*s22**2*s25*s26*s34*s36 + 1536*q9a**2*s11**2*s13*s16*s22**2*s25*s26*s34*s36 - 512*q9a**2*s11*s15**2*s16*s22**2*s25*s26*s34*s36 -  \
    1024*q9a**2*s11**2*s12*s16*s22*s23*s25*s26*s34*s36 + 1024*q9a**2*s11**2*s14*s16*s22*s24*s25*s26*s34*s36 - 512*q9a**2*s11*s15*s16**2*s22*s24*s25*s26*s34*s36 + 1536*q9a**2*s11**2*s12*s16*s24**2*s25*s26*s34*s36 +  \
    2560*q9a**2*s11**2*s12*s14*s22*s25**2*s26*s34*s36 - 512*q9a**2*s11*s14*s16**2*s22*s25**2*s26*s34*s36 - 1024*q9a**2*s11**2*s12**2*s24*s25**2*s26*s34*s36 - 512*q9a**2*s11*s12*s16**2*s24*s25**2*s26*s34*s36 +  \
    512*q9a**2*s11*s12**2*s16*s25**3*s26*s34*s36 - 512*q9a**2*s11**2*s13*s15*s22**2*s26**2*s34*s36 - 2048*q9a**2*s11**2*s14*s16*s22*s23*s26**2*s34*s36 - 512*q9a**2*s11**2*s14*s15*s22*s24*s26**2*s34*s36 +  \
    512*q9a**2*s11**2*s13*s16*s22*s24*s26**2*s34*s36 - 512*q9a**2*s11**2*s12*s16*s23*s24*s26**2*s34*s36 + 1024*q9a**2*s11**2*s12*s15*s24**2*s26**2*s34*s36 - 1024*q9a**2*s11**2*s12*s13*s22*s25*s26**2*s34*s36 +  \
    1024*q9a**2*s11**2*s14**2*s22*s25*s26**2*s34*s36 + 512*q9a**2*s11*s12*s15**2*s22*s25*s26**2*s34*s36 + 512*q9a**2*s11*s14*s15*s16*s22*s25*s26**2*s34*s36 + 1024*q9a**2*s11**2*s12**2*s23*s25*s26**2*s34*s36 -  \
    2560*q9a**2*s11**2*s12*s14*s24*s25*s26**2*s34*s36 + 512*q9a**2*s11*s12*s15*s16*s24*s25*s26**2*s34*s36 - 512*q9a**2*s11*s12**2*s15*s25**2*s26**2*s34*s36 + 512*q9a**2*s11*s12*s14*s16*s25**2*s26**2*s34*s36 +  \
    512*q9a**2*s11**2*s13*s14*s22*s26**3*s34*s36 + 1536*q9a**2*s11**2*s12*s14*s23*s26**3*s34*s36 - 1024*q9a**2*s11**2*s12*s13*s24*s26**3*s34*s36 - 512*q9a**2*s11*s12*s14*s15*s25*s26**3*s34*s36 -  \
    2048*q8a*q9a*s11**3*s15*s22**3*s33*s34*s36 - 1024*q8a*q9a*s11**3*s16*s22**2*s24*s33*s34*s36 + 6144*q8a*q9a*s11**3*s12*s22**2*s25*s33*s34*s36 - 1024*q8a*q9a*s11**2*s16**2*s22**2*s25*s33*s34*s36 -  \
    4096*q7a*q9a*s11**3*s22**3*s25*s33*s34*s36 + 5120*q8a*q9a*s11**3*s14*s22**2*s26*s33*s34*s36 + 1024*q8a*q9a*s11**2*s15*s16*s22**2*s26*s33*s34*s36 - 4096*q8a*q9a*s11**3*s12*s22*s24*s26*s33*s34*s36 +  \
    2048*q8a*q9a*s11**2*s16**2*s22*s24*s26*s33*s34*s36 - 2048*q8a*q9a*s11**2*s12*s16*s22*s25*s26*s33*s34*s36 + 3072*q7a*q9a*s11**2*s16*s22**2*s25*s26*s33*s34*s36 - 4096*q8a*q9a*s11**2*s14*s16*s22*s26**2*s33*s34*s36 -  \
    1024*q7a*q9a*s11**2*s15*s22**2*s26**2*s33*s34*s36 - 1024*q8a*q9a*s11**2*s12*s16*s24*s26**2*s33*s34*s36 + 1024*q7a*q9a*s11**2*s16*s22*s24*s26**2*s33*s34*s36 + 2048*q8a*q9a*s11**2*s12**2*s25*s26**2*s33*s34*s36 -  \
    2048*q7a*q9a*s11**2*s12*s22*s25*s26**2*s33*s34*s36 + 3072*q8a*q9a*s11**2*s12*s14*s26**3*s33*s34*s36 + 1024*q7a*q9a*s11**2*s14*s22*s26**3*s33*s34*s36 - 2048*q7a*q9a*s11**2*s12*s24*s26**3*s33*s34*s36 +  \
    2048*q8a*q9a*s11**3*s16*s22**2*s23*s34**2*s36 + 1536*q8a*q9a*s11**3*s15*s22**2*s24*s34**2*s36 - 1024*q8a*q9a*s11**3*s16*s22*s24**2*s34**2*s36 - 512*q8a*q9a*s11**3*s14*s22**2*s25*s34**2*s36 -  \
    1536*q8a*q9a*s11**2*s15*s16*s22**2*s25*s34**2*s36 - 2048*q8a*q9a*s11**3*s12*s22*s24*s25*s34**2*s36 + 1536*q8a*q9a*s11**2*s16**2*s22*s24*s25*s34**2*s36 + 1024*q7a*q9a*s11**3*s22**2*s24*s25*s34**2*s36 +  \
    1024*q8a*q9a*s11**2*s12*s16*s22*s25**2*s34**2*s36 - 512*q8a*q9a*s11*s16**3*s22*s25**2*s34**2*s36 + 512*q7a*q9a*s11**2*s16*s22**2*s25**2*s34**2*s36 - 1024*q8a*q9a*s11**3*s13*s22**2*s26*s34**2*s36 -  \
    512*q8a*q9a*s11**2*s15**2*s22**2*s26*s34**2*s36 - 2048*q8a*q9a*s11**2*s16**2*s22*s23*s26*s34**2*s36 - 1024*q7a*q9a*s11**3*s22**2*s23*s26*s34**2*s36 - 1024*q8a*q9a*s11**3*s14*s22*s24*s26*s34**2*s36 +  \
    2048*q8a*q9a*s11**3*s12*s24**2*s26*s34**2*s36 + 1024*q8a*q9a*s11**2*s12*s15*s22*s25*s26*s34**2*s36 + 1024*q8a*q9a*s11**2*s14*s16*s22*s25*s26*s34**2*s36 + 512*q8a*q9a*s11*s15*s16**2*s22*s25*s26*s34**2*s36 +  \
    1536*q7a*q9a*s11**2*s15*s22**2*s25*s26*s34**2*s36 - 2048*q8a*q9a*s11**2*s12*s16*s24*s25*s26*s34**2*s36 - 2048*q7a*q9a*s11**2*s16*s22*s24*s25*s26*s34**2*s36 + 512*q8a*q9a*s11*s12*s16**2*s25**2*s26*s34**2*s36 -  \
    2048*q7a*q9a*s11**2*s12*s22*s25**2*s26*s34**2*s36 + 512*q7a*q9a*s11*s16**2*s22*s25**2*s26*s34**2*s36 + 512*q8a*q9a*s11**2*s14*s15*s22*s26**2*s34**2*s36 + 1024*q8a*q9a*s11**2*s13*s16*s22*s26**2*s34**2*s36 +  \
    2048*q8a*q9a*s11**2*s12*s16*s23*s26**2*s34**2*s36 + 1024*q7a*q9a*s11**2*s16*s22*s23*s26**2*s34**2*s36 - 512*q8a*q9a*s11**2*s12*s15*s24*s26**2*s34**2*s36 - 512*q8a*q9a*s11**2*s12*s14*s25*s26**2*s34**2*s36 -  \
    512*q8a*q9a*s11*s12*s15*s16*s25*s26**2*s34**2*s36 - 1024*q7a*q9a*s11**2*s14*s22*s25*s26**2*s34**2*s36 - 512*q7a*q9a*s11*s15*s16*s22*s25*s26**2*s34**2*s36 + 3072*q7a*q9a*s11**2*s12*s24*s25*s26**2*s34**2*s36 -  \
    512*q7a*q9a*s11*s12*s16*s25**2*s26**2*s34**2*s36 - 1024*q8a*q9a*s11**2*s12*s13*s26**3*s34**2*s36 - 1024*q7a*q9a*s11**2*s12*s23*s26**3*s34**2*s36 + 512*q7a*q9a*s11*s12*s15*s25*s26**3*s34**2*s36 +  \
    1024*q8a**2*s11**3*s16*s22**2*s33*s34**2*s36 - 1024*q8a**2*s11**2*s16**2*s22*s26*s33*s34**2*s36 - 1024*q7a*q8a*s11**3*s22**2*s26*s33*s34**2*s36 + 1024*q8a**2*s11**2*s12*s16*s26**2*s33*s34**2*s36;
    v2_2= \
    1024*q7a*q8a*s11**2*s16*s22*s26**2*s33*s34**2*s36 - 1024*q7a*q8a*s11**2*s12*s26**3*s33*s34**2*s36 - 1024*q8a**2*s11**3*s15*s22**2*s34**3*s36 + 512*q8a**2*s11**3*s16*s22*s24*s34**3*s36 + 1024*q8a**2*s11**3*s12*s22*s25*s34**3*s36 -  \
    512*q8a**2*s11**2*s16**2*s22*s25*s34**3*s36 + 512*q8a**2*s11**3*s14*s22*s26*s34**3*s36 + 512*q8a**2*s11**2*s15*s16*s22*s26*s34**3*s36 - 1024*q8a**2*s11**3*s12*s24*s26*s34**3*s36 + 512*q8a**2*s11**2*s12*s16*s25*s26*s34**3*s36 -  \
    512*q8a**2*s11**2*s12*s15*s26**2*s34**3*s36 - 2048*q9a**2*s11**3*s14*s22**3*s23*s35*s36 - 1024*q9a**2*s11**3*s13*s22**3*s24*s35*s36 - 512*q9a**2*s11**2*s15**2*s22**3*s24*s35*s36 + 3072*q9a**2*s11**3*s12*s22**2*s23*s24*s35*s36 -  \
    512*q9a**2*s11**2*s16**2*s22**2*s23*s24*s35*s36 + 1024*q9a**2*s11**3*s14*s22**2*s24**2*s35*s36 + 512*q9a**2*s11**2*s15*s16*s22**2*s24**2*s35*s36 - 1024*q9a**2*s11**3*s12*s22*s24**3*s35*s36 -  \
    512*q9a**2*s11**2*s14*s15*s22**3*s25*s35*s36 + 2560*q9a**2*s11**2*s13*s16*s22**3*s25*s35*s36 + 512*q9a**2*s11*s15**2*s16*s22**3*s25*s35*s36 - 2560*q9a**2*s11**2*s12*s16*s22**2*s23*s25*s35*s36 +  \
    512*q9a**2*s11*s16**3*s22**2*s23*s25*s35*s36 + 1536*q9a**2*s11**2*s12*s15*s22**2*s24*s25*s35*s36 - 2048*q9a**2*s11**2*s14*s16*s22**2*s24*s25*s35*s36 - 512*q9a**2*s11*s15*s16**2*s22**2*s24*s25*s35*s36 +  \
    1536*q9a**2*s11**2*s12*s16*s22*s24**2*s25*s35*s36 + 512*q9a**2*s11**2*s12*s14*s22**2*s25**2*s35*s36 - 1024*q9a**2*s11*s12*s15*s16*s22**2*s25**2*s35*s36 + 1024*q9a**2*s11*s14*s16**2*s22**2*s25**2*s35*s36 -  \
    1024*q9a**2*s11**2*s12**2*s22*s24*s25**2*s35*s36 - 512*q9a**2*s11*s12*s16**2*s22*s24*s25**2*s35*s36 + 512*q9a**2*s11*s12**2*s16*s22*s25**3*s35*s36 - 512*q9a**2*s11**2*s13*s15*s22**3*s26*s35*s36 +  \
    512*q9a**2*s11*s15**3*s22**3*s26*s35*s36 + 512*q9a**2*s11**2*s12*s15*s22**2*s23*s26*s35*s36 + 1536*q9a**2*s11**2*s14*s16*s22**2*s23*s26*s35*s36 + 512*q9a**2*s11*s15*s16**2*s22**2*s23*s26*s35*s36 +  \
    1024*q9a**2*s11**2*s14*s15*s22**2*s24*s26*s35*s36 + 512*q9a**2*s11**2*s13*s16*s22**2*s24*s26*s35*s36 - 512*q9a**2*s11*s15**2*s16*s22**2*s24*s26*s35*s36 - 1024*q9a**2*s11**2*s12*s16*s22*s23*s24*s26*s35*s36 -  \
    1536*q9a**2*s11**2*s12*s15*s22*s24**2*s26*s35*s36 - 512*q9a**2*s11**2*s14*s16*s22*s24**2*s26*s35*s36 + 512*q9a**2*s11**2*s12*s16*s24**3*s26*s35*s36 - 2048*q9a**2*s11**2*s12*s13*s22**2*s25*s26*s35*s36 -  \
    512*q9a**2*s11**2*s14**2*s22**2*s25*s26*s35*s36 - 2048*q9a**2*s11*s12*s15**2*s22**2*s25*s26*s35*s36 + 1024*q9a**2*s11*s14*s15*s16*s22**2*s25*s26*s35*s36 - 2048*q9a**2*s11*s13*s16**2*s22**2*s25*s26*s35*s36 +  \
    2048*q9a**2*s11**2*s12**2*s22*s23*s25*s26*s35*s36 + 2048*q9a**2*s11**2*s12*s14*s22*s24*s25*s26*s35*s36 + 1024*q9a**2*s11*s12*s15*s16*s22*s24*s25*s26*s35*s36 + 512*q9a**2*s11*s14*s16**2*s22*s24*s25*s26*s35*s36 -  \
    1024*q9a**2*s11**2*s12**2*s24**2*s25*s26*s35*s36 - 512*q9a**2*s11*s12*s16**2*s24**2*s25*s26*s35*s36 + 2560*q9a**2*s11*s12**2*s15*s22*s25**2*s26*s35*s36 - 3072*q9a**2*s11*s12*s14*s16*s22*s25**2*s26*s35*s36 +  \
    1536*q9a**2*s11*s12**2*s16*s24*s25**2*s26*s35*s36 - 1024*q9a**2*s11*s12**3*s25**3*s26*s35*s36 - 512*q9a**2*s11**2*s13*s14*s22**2*s26**2*s35*s36 - 1024*q9a**2*s11*s14*s15**2*s22**2*s26**2*s35*s36 +  \
    1024*q9a**2*s11*s13*s15*s16*s22**2*s26**2*s35*s36 - 1024*q9a**2*s11**2*s12*s14*s22*s23*s26**2*s35*s36 - 2048*q9a**2*s11*s12*s15*s16*s22*s23*s26**2*s35*s36 + 512*q9a**2*s11**2*s14**2*s22*s24*s26**2*s35*s36 +  \
    1536*q9a**2*s11*s12*s15**2*s22*s24*s26**2*s35*s36 - 512*q9a**2*s11*s14*s15*s16*s22*s24*s26**2*s35*s36 + 1024*q9a**2*s11**2*s12**2*s23*s24*s26**2*s35*s36 - 512*q9a**2*s11**2*s12*s14*s24**2*s26**2*s35*s36 +  \
    512*q9a**2*s11*s12*s15*s16*s24**2*s26**2*s35*s36 + 1024*q9a**2*s11*s12*s14*s15*s22*s25*s26**2*s35*s36 + 3072*q9a**2*s11*s12*s13*s16*s22*s25*s26**2*s35*s36 - 512*q9a**2*s11*s14**2*s16*s22*s25*s26**2*s35*s36 -  \
    512*q9a**2*s11*s12**2*s16*s23*s25*s26**2*s35*s36 - 2560*q9a**2*s11*s12**2*s15*s24*s25*s26**2*s35*s36 + 512*q9a**2*s11*s12*s14*s16*s24*s25*s26**2*s35*s36 + 1024*q9a**2*s11*s12**2*s14*s25**2*s26**2*s35*s36 -  \
    1024*q9a**2*s11*s12*s13*s15*s22*s26**3*s35*s36 + 512*q9a**2*s11*s14**2*s15*s22*s26**3*s35*s36 + 1536*q9a**2*s11*s12**2*s15*s23*s26**3*s35*s36 - 512*q9a**2*s11*s12*s14*s15*s24*s26**3*s35*s36 -  \
    1024*q9a**2*s11*s12**2*s13*s25*s26**3*s35*s36 - 4096*q8a*q9a*s11**3*s14*s22**3*s33*s35*s36 + 6144*q8a*q9a*s11**3*s12*s22**2*s24*s33*s35*s36 - 1024*q8a*q9a*s11**2*s16**2*s22**2*s24*s33*s35*s36 -  \
    2048*q7a*q9a*s11**3*s22**3*s24*s33*s35*s36 - 5120*q8a*q9a*s11**2*s12*s16*s22**2*s25*s33*s35*s36 + 1024*q8a*q9a*s11*s16**3*s22**2*s25*s33*s35*s36 + 5120*q7a*q9a*s11**2*s16*s22**3*s25*s33*s35*s36 +  \
    1024*q8a*q9a*s11**2*s12*s15*s22**2*s26*s33*s35*s36 + 3072*q8a*q9a*s11**2*s14*s16*s22**2*s26*s33*s35*s36 + 1024*q8a*q9a*s11*s15*s16**2*s22**2*s26*s33*s35*s36 - 1024*q7a*q9a*s11**2*s15*s22**3*s26*s33*s35*s36 -  \
    2048*q8a*q9a*s11**2*s12*s16*s22*s24*s26*s33*s35*s36 + 1024*q7a*q9a*s11**2*s16*s22**2*s24*s26*s33*s35*s36 + 4096*q8a*q9a*s11**2*s12**2*s22*s25*s26*s33*s35*s36 - 4096*q7a*q9a*s11**2*s12*s22**2*s25*s26*s33*s35*s36 -  \
    4096*q7a*q9a*s11*s16**2*s22**2*s25*s26*s33*s35*s36 - 2048*q8a*q9a*s11**2*s12*s14*s22*s26**2*s33*s35*s36 - 4096*q8a*q9a*s11*s12*s15*s16*s22*s26**2*s33*s35*s36 - 1024*q7a*q9a*s11**2*s14*s22**2*s26**2*s33*s35*s36 +  \
    2048*q7a*q9a*s11*s15*s16*s22**2*s26**2*s33*s35*s36 + 2048*q8a*q9a*s11**2*s12**2*s24*s26**2*s33*s35*s36 - 1024*q8a*q9a*s11*s12**2*s16*s25*s26**2*s33*s35*s36 + 6144*q7a*q9a*s11*s12*s16*s22*s25*s26**2*s33*s35*s36 +  \
    3072*q8a*q9a*s11*s12**2*s15*s26**3*s33*s35*s36 - 2048*q7a*q9a*s11*s12*s15*s22*s26**3*s33*s35*s36 - 2048*q7a*q9a*s11*s12**2*s25*s26**3*s33*s35*s36 + 4096*q8a*q9a*s11**3*s13*s22**3*s34*s35*s36 -  \
    8192*q8a*q9a*s11**3*s12*s22**2*s23*s34*s35*s36 + 4096*q7a*q9a*s11**3*s22**3*s23*s34*s35*s36 - 2048*q8a*q9a*s11**3*s14*s22**2*s24*s34*s35*s36 + 2048*q8a*q9a*s11**3*s12*s22*s24**2*s34*s35*s36 +  \
    512*q8a*q9a*s11**2*s16**2*s22*s24**2*s34*s35*s36 + 2048*q8a*q9a*s11**2*s12*s15*s22**2*s25*s34*s35*s36 + 2048*q8a*q9a*s11**2*s14*s16*s22**2*s25*s34*s35*s36 - 2048*q7a*q9a*s11**2*s15*s22**3*s25*s34*s35*s36 -  \
    2048*q8a*q9a*s11**2*s12*s16*s22*s24*s25*s34*s35*s36 - 512*q8a*q9a*s11*s16**3*s22*s24*s25*s34*s35*s36 - 2048*q8a*q9a*s11**2*s12**2*s22*s25**2*s34*s35*s36 + 1536*q8a*q9a*s11*s12*s16**2*s22*s25**2*s34*s35*s36 +  \
    2048*q7a*q9a*s11**2*s12*s22**2*s25**2*s34*s35*s36 - 1536*q7a*q9a*s11*s16**2*s22**2*s25**2*s34*s35*s36 - 4096*q8a*q9a*s11**2*s13*s16*s22**2*s26*s34*s35*s36 + 8192*q8a*q9a*s11**2*s12*s16*s22*s23*s26*s34*s35*s36 -  \
    4096*q7a*q9a*s11**2*s16*s22**2*s23*s26*s34*s35*s36 + 1024*q8a*q9a*s11**2*s14*s16*s22*s24*s26*s34*s35*s36 - 512*q8a*q9a*s11*s15*s16**2*s22*s24*s26*s34*s35*s36 - 2048*q8a*q9a*s11**2*s12*s16*s24**2*s26*s34*s35*s36 -  \
    4096*q8a*q9a*s11**2*s12*s14*s22*s25*s26*s34*s35*s36 - 1024*q8a*q9a*s11*s12*s15*s16*s22*s25*s26*s34*s35*s36 - 512*q8a*q9a*s11*s14*s16**2*s22*s25*s26*s34*s35*s36 + 2048*q7a*q9a*s11**2*s14*s22**2*s25*s26*s34*s35*s36 +  \
    1024*q7a*q9a*s11*s15*s16*s22**2*s25*s26*s34*s35*s36 + 4096*q8a*q9a*s11**2*s12**2*s24*s25*s26*s34*s35*s36 + 1024*q8a*q9a*s11*s12*s16**2*s24*s25*s26*s34*s35*s36 - 2048*q7a*q9a*s11**2*s12*s22*s24*s25*s26*s34*s35*s36 +  \
    1536*q7a*q9a*s11*s16**2*s22*s24*s25*s26*s34*s35*s36 - 2048*q8a*q9a*s11*s12**2*s16*s25**2*s26*s34*s35*s36 + 2048*q7a*q9a*s11*s12*s16*s22*s25**2*s26*s34*s35*s36 + 4096*q8a*q9a*s11**2*s12*s13*s22*s26**2*s34*s35*s36 -  \
    1536*q8a*q9a*s11**2*s14**2*s22*s26**2*s34*s35*s36 - 512*q8a*q9a*s11*s12*s15**2*s22*s26**2*s34*s35*s36 + 1536*q8a*q9a*s11*s14*s15*s16*s22*s26**2*s34*s35*s36 + 512*q7a*q9a*s11*s15**2*s22**2*s26**2*s34*s35*s36 -  \
    8192*q8a*q9a*s11**2*s12**2*s23*s26**2*s34*s35*s36 + 4096*q7a*q9a*s11**2*s12*s22*s23*s26**2*s34*s35*s36 + 2048*q8a*q9a*s11**2*s12*s14*s24*s26**2*s34*s35*s36 - 512*q7a*q9a*s11*s15*s16*s22*s24*s26**2*s34*s35*s36 +  \
    2048*q8a*q9a*s11*s12**2*s15*s25*s26**2*s34*s35*s36 - 2048*q7a*q9a*s11*s12*s15*s22*s25*s26**2*s34*s35*s36 - 512*q7a*q9a*s11*s14*s16*s22*s25*s26**2*s34*s35*s36 - 2048*q7a*q9a*s11*s12*s16*s24*s25*s26**2*s34*s35*s36 -  \
    1024*q8a*q9a*s11*s12*s14*s15*s26**3*s34*s35*s36 - 512*q7a*q9a*s11*s14*s15*s22*s26**3*s34*s35*s36 + 1024*q7a*q9a*s11*s12*s15*s24*s26**3*s34*s35*s36 + 1024*q7a*q9a*s11*s12*s14*s25*s26**3*s34*s35*s36 -  \
    4096*q8a**2*s11**3*s12*s22**2*s33*s34*s35*s36 + 4096*q7a*q8a*s11**3*s22**3*s33*s34*s35*s36 + 4096*q8a**2*s11**2*s12*s16*s22*s26*s33*s34*s35*s36 - 4096*q7a*q8a*s11**2*s16*s22**2*s26*s33*s34*s35*s36 -  \
    4096*q8a**2*s11**2*s12**2*s26**2*s33*s34*s35*s36 + 4096*q7a*q8a*s11**2*s12*s22*s26**2*s33*s34*s35*s36 + 2048*q8a**2*s11**3*s14*s22**2*s34**2*s35*s36 - 1024*q8a**2*s11**3*s12*s22*s24*s34**2*s35*s36 -  \
    512*q8a**2*s11**2*s16**2*s22*s24*s34**2*s35*s36 - 1024*q7a*q8a*s11**3*s22**2*s24*s34**2*s35*s36 - 512*q8a**2*s11**2*s12*s16*s22*s25*s34**2*s35*s36 + 512*q8a**2*s11*s16**3*s22*s25*s34**2*s35*s36 +  \
    512*q7a*q8a*s11**2*s16*s22**2*s25*s34**2*s35*s36 + 512*q8a**2*s11**2*s12*s15*s22*s26*s34**2*s35*s36 - 2560*q8a**2*s11**2*s14*s16*s22*s26*s34**2*s35*s36 + 512*q8a**2*s11*s15*s16**2*s22*s26*s34**2*s35*s36 -  \
    512*q7a*q8a*s11**2*s15*s22**2*s26*s34**2*s35*s36 + 1536*q8a**2*s11**2*s12*s16*s24*s26*s34**2*s35*s36 + 2048*q7a*q8a*s11**2*s16*s22*s24*s26*s34**2*s35*s36 - 1024*q8a**2*s11**2*s12**2*s25*s26*s34**2*s35*s36 -  \
    512*q8a**2*s11*s12*s16**2*s25*s26*s34**2*s35*s36 + 2048*q7a*q8a*s11**2*s12*s22*s25*s26*s34**2*s35*s36 - 1536*q7a*q8a*s11*s16**2*s22*s25*s26*s34**2*s35*s36 - 1024*q7a**2*s11**2*s22**2*s25*s26*s34**2*s35*s36 +  \
    1536*q8a**2*s11**2*s12*s14*s26**2*s34**2*s35*s36 - 512*q8a**2*s11*s12*s15*s16*s26**2*s34**2*s35*s36 + 1024*q7a*q8a*s11**2*s14*s22*s26**2*s34**2*s35*s36 - 512*q7a*q8a*s11*s15*s16*s22*s26**2*s34**2*s35*s36 -  \
    3072*q7a*q8a*s11**2*s12*s24*s26**2*s34**2*s35*s36 + 1536*q7a*q8a*s11*s12*s16*s25*s26**2*s34**2*s35*s36 + 1024*q7a**2*s11*s16*s22*s25*s26**2*s34**2*s35*s36 + 512*q7a*q8a*s11*s12*s15*s26**3*s34**2*s35*s36 -  \
    1024*q7a**2*s11*s12*s25*s26**3*s34**2*s35*s36 + 1024*q8a*q9a*s11**2*s14*s15*s22**3*s35**2*s36 - 1024*q8a*q9a*s11**2*s13*s16*s22**3*s35**2*s36 + 2048*q8a*q9a*s11**2*s12*s16*s22**2*s23*s35**2*s36 -  \
    1024*q7a*q9a*s11**2*s16*s22**3*s23*s35**2*s36 - 2560*q8a*q9a*s11**2*s12*s15*s22**2*s24*s35**2*s36 + 1536*q8a*q9a*s11**2*s14*s16*s22**2*s24*s35**2*s36 + 1536*q7a*q9a*s11**2*s15*s22**3*s24*s35**2*s36 -  \
    1024*q8a*q9a*s11**2*s12*s16*s22*s24**2*s35**2*s36 - 512*q7a*q9a*s11**2*s16*s22**2*s24**2*s35**2*s36 - 512*q8a*q9a*s11**2*s12*s14*s22**2*s25*s35**2*s36 + 1024*q8a*q9a*s11*s12*s15*s16*s22**2*s25*s35**2*s36 -  \
    1024*q8a*q9a*s11*s14*s16**2*s22**2*s25*s35**2*s36 - 512*q7a*q9a*s11**2*s14*s22**3*s25*s35**2*s36 - 1024*q7a*q9a*s11*s15*s16*s22**3*s25*s35**2*s36 + 2048*q8a*q9a*s11**2*s12**2*s22*s24*s25*s35**2*s36 +  \
    512*q8a*q9a*s11*s12*s16**2*s22*s24*s25*s35**2*s36 - 1024*q7a*q9a*s11**2*s12*s22**2*s24*s25*s35**2*s36 + 512*q7a*q9a*s11*s16**2*s22**2*s24*s25*s35**2*s36 - 1024*q8a*q9a*s11*s12**2*s16*s22*s25**2*s35**2*s36 +  \
    1024*q7a*q9a*s11*s12*s16*s22**2*s25**2*s35**2*s36 - 1024*q8a*q9a*s11**2*s12*s13*s22**2*s26*s35**2*s36 + 512*q8a*q9a*s11**2*s14**2*s22**2*s26*s35**2*s36 + 1024*q8a*q9a*s11*s12*s15**2*s22**2*s26*s35**2*s36 -  \
    2048*q8a*q9a*s11*s14*s15*s16*s22**2*s26*s35**2*s36 + 1024*q8a*q9a*s11*s13*s16**2*s22**2*s26*s35**2*s36 + 2048*q7a*q9a*s11**2*s13*s22**3*s26*s35**2*s36 - 1024*q7a*q9a*s11*s15**2*s22**3*s26*s35**2*s36 -  \
    2048*q8a*q9a*s11*s12*s16**2*s22*s23*s26*s35**2*s36 - 1024*q7a*q9a*s11**2*s12*s22**2*s23*s26*s35**2*s36 + 1024*q7a*q9a*s11*s16**2*s22**2*s23*s26*s35**2*s36 - 1024*q8a*q9a*s11**2*s12*s14*s22*s24*s26*s35**2*s36 +  \
    2048*q8a*q9a*s11*s12*s15*s16*s22*s24*s26*s35**2*s36 - 512*q8a*q9a*s11*s14*s16**2*s22*s24*s26*s35**2*s36 - 1536*q7a*q9a*s11**2*s14*s22**2*s24*s26*s35**2*s36 + 512*q8a*q9a*s11*s12*s16**2*s24**2*s26*s35**2*s36 +  \
    2048*q7a*q9a*s11**2*s12*s22*s24**2*s26*s35**2*s36 - 3072*q8a*q9a*s11*s12**2*s15*s22*s25*s26*s35**2*s36 + 3072*q8a*q9a*s11*s12*s14*s16*s22*s25*s26*s35**2*s36 + 3072*q7a*q9a*s11*s12*s15*s22**2*s25*s26*s35**2*s36;
    v2_3= \
    1024*q7a*q9a*s11*s14*s16*s22**2*s25*s26*s35**2*s36 - 2048*q8a*q9a*s11*s12**2*s16*s24*s25*s26*s35**2*s36 - 2048*q7a*q9a*s11*s12*s16*s22*s24*s25*s26*s35**2*s36 + 2048*q8a*q9a*s11*s12**3*s25**2*s26*s35**2*s36 -  \
    2048*q7a*q9a*s11*s12**2*s22*s25**2*s26*s35**2*s36 + 512*q8a*q9a*s11*s12*s14*s15*s22*s26**2*s35**2*s36 + 512*q8a*q9a*s11*s14**2*s16*s22*s26**2*s35**2*s36 + 1536*q7a*q9a*s11*s14*s15*s22**2*s26**2*s35**2*s36 -  \
    2048*q7a*q9a*s11*s13*s16*s22**2*s26**2*s35**2*s36 + 2048*q8a*q9a*s11*s12**2*s16*s23*s26**2*s35**2*s36 - 512*q8a*q9a*s11*s12**2*s15*s24*s26**2*s35**2*s36 - 512*q8a*q9a*s11*s12*s14*s16*s24*s26**2*s35**2*s36 -  \
    1536*q7a*q9a*s11*s12*s15*s22*s24*s26**2*s35**2*s36 + 512*q7a*q9a*s11*s14*s16*s22*s24*s26**2*s35**2*s36 - 512*q7a*q9a*s11*s12*s16*s24**2*s26**2*s35**2*s36 - 512*q8a*q9a*s11*s12**2*s14*s25*s26**2*s35**2*s36 -  \
    2560*q7a*q9a*s11*s12*s14*s22*s25*s26**2*s35**2*s36 + 3072*q7a*q9a*s11*s12**2*s24*s25*s26**2*s35**2*s36 - 1024*q8a*q9a*s11*s12**2*s13*s26**3*s35**2*s36 + 2048*q7a*q9a*s11*s12*s13*s22*s26**3*s35**2*s36 -  \
    512*q7a*q9a*s11*s14**2*s22*s26**3*s35**2*s36 - 1024*q7a*q9a*s11*s12**2*s23*s26**3*s35**2*s36 + 512*q7a*q9a*s11*s12*s14*s24*s26**3*s35**2*s36 + 1024*q8a**2*s11**2*s12*s16*s22**2*s33*s35**2*s36 -  \
    1024*q7a*q8a*s11**2*s16*s22**3*s33*s35**2*s36 - 1024*q8a**2*s11*s12*s16**2*s22*s26*s33*s35**2*s36 - 1024*q7a*q8a*s11**2*s12*s22**2*s26*s33*s35**2*s36 + 1024*q7a*q8a*s11*s16**2*s22**2*s26*s33*s35**2*s36 +  \
    1024*q7a**2*s11**2*s22**3*s26*s33*s35**2*s36 + 1024*q8a**2*s11*s12**2*s16*s26**2*s33*s35**2*s36 - 1024*q7a**2*s11*s16*s22**2*s26**2*s33*s35**2*s36 - 1024*q7a*q8a*s11*s12**2*s26**3*s33*s35**2*s36 +  \
    1024*q7a**2*s11*s12*s22*s26**3*s33*s35**2*s36 + 1024*q8a**2*s11**2*s12*s15*s22**2*s34*s35**2*s36 - 1024*q8a**2*s11**2*s14*s16*s22**2*s34*s35**2*s36 - 1024*q7a*q8a*s11**2*s15*s22**3*s34*s35**2*s36 +  \
    1536*q8a**2*s11**2*s12*s16*s22*s24*s34*s35**2*s36 - 512*q7a*q8a*s11**2*s16*s22**2*s24*s34*s35**2*s36 + 1024*q8a**2*s11**2*s12**2*s22*s25*s34*s35**2*s36 - 1024*q8a**2*s11*s12*s16**2*s22*s25*s34*s35**2*s36 -  \
    3072*q7a*q8a*s11**2*s12*s22**2*s25*s34*s35**2*s36 + 1024*q7a*q8a*s11*s16**2*s22**2*s25*s34*s35**2*s36 + 2048*q7a**2*s11**2*s22**3*s25*s34*s35**2*s36 + 512*q8a**2*s11**2*s12*s14*s22*s26*s34*s35**2*s36 -  \
    2048*q8a**2*s11*s12*s15*s16*s22*s26*s34*s35**2*s36 + 1024*q8a**2*s11*s14*s16**2*s22*s26*s34*s35**2*s36 + 512*q7a*q8a*s11**2*s14*s22**2*s26*s34*s35**2*s36 + 2048*q7a*q8a*s11*s15*s16*s22**2*s26*s34*s35**2*s36 -  \
    1024*q8a**2*s11**2*s12**2*s24*s26*s34*s35**2*s36 - 512*q8a**2*s11*s12*s16**2*s24*s26*s34*s35**2*s36 - 512*q7a*q8a*s11*s16**2*s22*s24*s26*s34*s35**2*s36 + 1536*q8a**2*s11*s12**2*s16*s25*s26*s34*s35**2*s36 +  \
    1024*q7a*q8a*s11*s12*s16*s22*s25*s26*s34*s35**2*s36 - 2560*q7a**2*s11*s16*s22**2*s25*s26*s34*s35**2*s36 + 1536*q8a**2*s11*s12**2*s15*s26**2*s34*s35**2*s36 - 512*q8a**2*s11*s12*s14*s16*s26**2*s34*s35**2*s36 -  \
    1024*q7a*q8a*s11*s12*s15*s22*s26**2*s34*s35**2*s36 - 1536*q7a*q8a*s11*s14*s16*s22*s26**2*s34*s35**2*s36 - 512*q7a**2*s11*s15*s22**2*s26**2*s34*s35**2*s36 + 1536*q7a*q8a*s11*s12*s16*s24*s26**2*s34*s35**2*s36 +  \
    512*q7a**2*s11*s16*s22*s24*s26**2*s34*s35**2*s36 - 3072*q7a*q8a*s11*s12**2*s25*s26**2*s34*s35**2*s36 + 3072*q7a**2*s11*s12*s22*s25*s26**2*s34*s35**2*s36 + 512*q7a*q8a*s11*s12*s14*s26**3*s34*s35**2*s36 +  \
    512*q7a**2*s11*s14*s22*s26**3*s34*s35**2*s36 - 1024*q7a**2*s11*s12*s24*s26**3*s34*s35**2*s36 - 1024*q8a**2*s11**2*s12**2*s22*s24*s35**3*s36 + 2048*q7a*q8a*s11**2*s12*s22**2*s24*s35**3*s36 - 1024*q7a**2*s11**2*s22**3*s24*s35**3*s36 +  \
    512*q8a**2*s11*s12**2*s16*s22*s25*s35**3*s36 - 1024*q7a*q8a*s11*s12*s16*s22**2*s25*s35**3*s36 + 512*q7a**2*s11*s16*s22**3*s25*s35**3*s36 + 512*q8a**2*s11*s12**2*s15*s22*s26*s35**3*s36 -  \
    1024*q7a*q8a*s11*s12*s15*s22**2*s26*s35**3*s36 + 512*q7a**2*s11*s15*s22**3*s26*s35**3*s36 + 512*q8a**2*s11*s12**2*s16*s24*s26*s35**3*s36 - 1024*q7a*q8a*s11*s12*s16*s22*s24*s26*s35**3*s36 +  \
    512*q7a**2*s11*s16*s22**2*s24*s26*s35**3*s36 - 1024*q8a**2*s11*s12**3*s25*s26*s35**3*s36 + 2048*q7a*q8a*s11*s12**2*s22*s25*s26*s35**3*s36 - 1024*q7a**2*s11*s12*s22**2*s25*s26*s35**3*s36 -  \
    512*q8a**2*s11*s12**2*s14*s26**2*s35**3*s36 + 1024*q7a*q8a*s11*s12*s14*s22*s26**2*s35**3*s36 - 512*q7a**2*s11*s14*s22**2*s26**2*s35**3*s36 - 1024*q9a**2*s11**3*s13*s22**3*s23*s36**2 + 512*q9a**2*s11**2*s15**2*s22**3*s23*s36**2 +  \
    1024*q9a**2*s11**3*s12*s22**2*s23**2*s36**2 + 512*q9a**2*s11**2*s16**2*s22**2*s23**2*s36**2 + 1536*q9a**2*s11**3*s14*s22**2*s23*s24*s36**2 - 512*q9a**2*s11**2*s15*s16*s22**2*s23*s24*s36**2 + 512*q9a**2*s11**3*s13*s22**2*s24**2*s36**2 -  \
    2048*q9a**2*s11**3*s12*s22*s23*s24**2*s36**2 - 512*q9a**2*s11**3*s14*s22*s24**3*s36**2 + 512*q9a**2*s11**3*s12*s24**4*s36**2 + 1536*q9a**2*s11**2*s13*s15*s22**3*s25*s36**2 - 512*q9a**2*s11*s15**3*s22**3*s25*s36**2 -  \
    2560*q9a**2*s11**2*s12*s15*s22**2*s23*s25*s36**2 + 512*q9a**2*s11**2*s14*s16*s22**2*s23*s25*s36**2 - 512*q9a**2*s11*s15*s16**2*s22**2*s23*s25*s36**2 - 512*q9a**2*s11**2*s14*s15*s22**2*s24*s25*s36**2 -  \
    2048*q9a**2*s11**2*s13*s16*s22**2*s24*s25*s36**2 + 512*q9a**2*s11*s15**2*s16*s22**2*s24*s25*s36**2 + 2048*q9a**2*s11**2*s12*s16*s22*s23*s24*s25*s36**2 + 512*q9a**2*s11**2*s12*s15*s22*s24**2*s25*s36**2 +  \
    1024*q9a**2*s11**2*s14*s16*s22*s24**2*s25*s36**2 - 1024*q9a**2*s11**2*s12*s16*s24**3*s25*s36**2 - 1536*q9a**2*s11**2*s12*s13*s22**2*s25**2*s36**2 + 1536*q9a**2*s11**2*s14**2*s22**2*s25**2*s36**2 +  \
    1536*q9a**2*s11*s12*s15**2*s22**2*s25**2*s36**2 - 1536*q9a**2*s11*s14*s15*s16*s22**2*s25**2*s36**2 + 1536*q9a**2*s11*s13*s16**2*s22**2*s25**2*s36**2 + 2048*q9a**2*s11**2*s12**2*s22*s23*s25**2*s36**2 -  \
    1024*q9a**2*s11*s12*s16**2*s22*s23*s25**2*s36**2 - 2560*q9a**2*s11**2*s12*s14*s22*s24*s25**2*s36**2 + 512*q9a**2*s11*s12*s15*s16*s22*s24*s25**2*s36**2 - 512*q9a**2*s11*s14*s16**2*s22*s24*s25**2*s36**2 +  \
    1024*q9a**2*s11**2*s12**2*s24**2*s25**2*s36**2 + 512*q9a**2*s11*s12*s16**2*s24**2*s25**2*s36**2 - 1536*q9a**2*s11*s12**2*s15*s22*s25**3*s36**2 + 1536*q9a**2*s11*s12*s14*s16*s22*s25**3*s36**2 -  \
    1024*q9a**2*s11*s12**2*s16*s24*s25**3*s36**2 + 512*q9a**2*s11*s12**3*s25**4*s36**2 - 2048*q9a**2*s11**2*s14*s15*s22**2*s23*s26*s36**2 + 2048*q9a**2*s11**2*s13*s16*s22**2*s23*s26*s36**2 -  \
    3072*q9a**2*s11**2*s12*s16*s22*s23**2*s26*s36**2 - 512*q9a**2*s11**2*s13*s15*s22**2*s24*s26*s36**2 + 3072*q9a**2*s11**2*s12*s15*s22*s23*s24*s26*s36**2 - 1024*q9a**2*s11**2*s14*s16*s22*s23*s24*s26*s36**2 +  \
    512*q9a**2*s11**2*s14*s15*s22*s24**2*s26*s36**2 + 1024*q9a**2*s11**2*s12*s16*s23*s24**2*s26*s36**2 - 512*q9a**2*s11**2*s12*s15*s24**3*s26*s36**2 + 512*q9a**2*s11**2*s13*s14*s22**2*s25*s26*s36**2 +  \
    1024*q9a**2*s11*s14*s15**2*s22**2*s25*s26*s36**2 - 1024*q9a**2*s11*s13*s15*s16*s22**2*s25*s26*s36**2 + 1024*q9a**2*s11**2*s12*s14*s22*s23*s25*s26*s36**2 + 2048*q9a**2*s11*s12*s15*s16*s22*s23*s25*s26*s36**2 +  \
    2048*q9a**2*s11**2*s12*s13*s22*s24*s25*s26*s36**2 - 1536*q9a**2*s11**2*s14**2*s22*s24*s25*s26*s36**2 - 1536*q9a**2*s11*s12*s15**2*s22*s24*s25*s26*s36**2 + 512*q9a**2*s11*s14*s15*s16*s22*s24*s25*s26*s36**2 -  \
    4096*q9a**2*s11**2*s12**2*s23*s24*s25*s26*s36**2 + 1536*q9a**2*s11**2*s12*s14*s24**2*s25*s26*s36**2 - 512*q9a**2*s11*s12*s15*s16*s24**2*s25*s26*s36**2 - 512*q9a**2*s11*s12*s14*s15*s22*s25**2*s26*s36**2 -  \
    2048*q9a**2*s11*s12*s13*s16*s22*s25**2*s26*s36**2 + 512*q9a**2*s11*s14**2*s16*s22*s25**2*s26*s36**2 + 1024*q9a**2*s11*s12**2*s16*s23*s25**2*s26*s36**2 + 1536*q9a**2*s11*s12**2*s15*s24*s25**2*s26*s36**2 -  \
    512*q9a**2*s11*s12*s14*s16*s24*s25**2*s26*s36**2 - 512*q9a**2*s11*s12**2*s14*s25**3*s26*s36**2 + 512*q9a**2*s11**2*s13**2*s22**2*s26**2*s36**2 - 3072*q9a**2*s11**2*s12*s13*s22*s23*s26**2*s36**2 +  \
    1536*q9a**2*s11**2*s14**2*s22*s23*s26**2*s36**2 + 3072*q9a**2*s11**2*s12**2*s23**2*s26**2*s36**2 - 512*q9a**2*s11**2*s13*s14*s22*s24*s26**2*s36**2 - 1536*q9a**2*s11**2*s12*s14*s23*s24*s26**2*s36**2 +  \
    512*q9a**2*s11**2*s12*s13*s24**2*s26**2*s36**2 + 1024*q9a**2*s11*s12*s13*s15*s22*s25*s26**2*s36**2 - 512*q9a**2*s11*s14**2*s15*s22*s25*s26**2*s36**2 - 1536*q9a**2*s11*s12**2*s15*s23*s25*s26**2*s36**2 +  \
    512*q9a**2*s11*s12*s14*s15*s24*s25*s26**2*s36**2 + 512*q9a**2*s11*s12**2*s13*s25**2*s26**2*s36**2 - 2048*q8a*q9a*s11**3*s13*s22**3*s33*s36**2 + 1024*q8a*q9a*s11**2*s15**2*s22**3*s33*s36**2 +  \
    4096*q8a*q9a*s11**3*s12*s22**2*s23*s33*s36**2 + 2048*q8a*q9a*s11**2*s16**2*s22**2*s23*s33*s36**2 - 2048*q7a*q9a*s11**3*s22**3*s23*s33*s36**2 + 3072*q8a*q9a*s11**3*s14*s22**2*s24*s33*s36**2 -  \
    1024*q8a*q9a*s11**2*s15*s16*s22**2*s24*s33*s36**2 - 4096*q8a*q9a*s11**3*s12*s22*s24**2*s33*s36**2 + 1024*q7a*q9a*s11**3*s22**2*s24**2*s33*s36**2 - 5120*q8a*q9a*s11**2*s12*s15*s22**2*s25*s33*s36**2 +  \
    1024*q8a*q9a*s11**2*s14*s16*s22**2*s25*s33*s36**2 - 1024*q8a*q9a*s11*s15*s16**2*s22**2*s25*s33*s36**2 + 3072*q7a*q9a*s11**2*s15*s22**3*s25*s33*s36**2 + 4096*q8a*q9a*s11**2*s12*s16*s22*s24*s25*s33*s36**2 -  \
    4096*q7a*q9a*s11**2*s16*s22**2*s24*s25*s33*s36**2 + 4096*q8a*q9a*s11**2*s12**2*s22*s25**2*s33*s36**2 - 2048*q8a*q9a*s11*s12*s16**2*s22*s25**2*s33*s36**2 - 3072*q7a*q9a*s11**2*s12*s22**2*s25**2*s33*s36**2 +  \
    3072*q7a*q9a*s11*s16**2*s22**2*s25**2*s33*s36**2 - 4096*q8a*q9a*s11**2*s14*s15*s22**2*s26*s33*s36**2 + 4096*q8a*q9a*s11**2*s13*s16*s22**2*s26*s33*s36**2 - 12288*q8a*q9a*s11**2*s12*s16*s22*s23*s26*s33*s36**2 +  \
    4096*q7a*q9a*s11**2*s16*s22**2*s23*s26*s33*s36**2 + 6144*q8a*q9a*s11**2*s12*s15*s22*s24*s26*s33*s36**2 - 2048*q8a*q9a*s11**2*s14*s16*s22*s24*s26*s33*s36**2 - 1024*q7a*q9a*s11**2*s15*s22**2*s24*s26*s33*s36**2 +  \
    2048*q8a*q9a*s11**2*s12*s16*s24**2*s26*s33*s36**2 + 2048*q8a*q9a*s11**2*s12*s14*s22*s25*s26*s33*s36**2 + 4096*q8a*q9a*s11*s12*s15*s16*s22*s25*s26*s33*s36**2 + 1024*q7a*q9a*s11**2*s14*s22**2*s25*s26*s33*s36**2 -  \
    2048*q7a*q9a*s11*s15*s16*s22**2*s25*s26*s33*s36**2 - 8192*q8a*q9a*s11**2*s12**2*s24*s25*s26*s33*s36**2 + 4096*q7a*q9a*s11**2*s12*s22*s24*s25*s26*s33*s36**2 + 2048*q8a*q9a*s11*s12**2*s16*s25**2*s26*s33*s36**2 -  \
    4096*q7a*q9a*s11*s12*s16*s22*s25**2*s26*s33*s36**2 - 6144*q8a*q9a*s11**2*s12*s13*s22*s26**2*s33*s36**2 + 3072*q8a*q9a*s11**2*s14**2*s22*s26**2*s33*s36**2 + 2048*q7a*q9a*s11**2*s13*s22**2*s26**2*s33*s36**2 +  \
    12288*q8a*q9a*s11**2*s12**2*s23*s26**2*s33*s36**2 - 6144*q7a*q9a*s11**2*s12*s22*s23*s26**2*s33*s36**2 - 3072*q8a*q9a*s11**2*s12*s14*s24*s26**2*s33*s36**2 - 1024*q7a*q9a*s11**2*s14*s22*s24*s26**2*s33*s36**2 +  \
    1024*q7a*q9a*s11**2*s12*s24**2*s26**2*s33*s36**2 - 3072*q8a*q9a*s11*s12**2*s15*s25*s26**2*s33*s36**2 + 2048*q7a*q9a*s11*s12*s15*s22*s25*s26**2*s33*s36**2 + 1024*q7a*q9a*s11*s12**2*s25**2*s26**2*s33*s36**2 +  \
    1024*q8a**2*s11**3*s12*s22**2*s33**2*s36**2 + 512*q8a**2*s11**2*s16**2*s22**2*s33**2*s36**2 - 1024*q7a*q8a*s11**3*s22**3*s33**2*s36**2 - 3072*q8a**2*s11**2*s12*s16*s22*s26*s33**2*s36**2 + 2048*q7a*q8a*s11**2*s16*s22**2*s26*s33**2*s36**2 +  \
    3072*q8a**2*s11**2*s12**2*s26**2*s33**2*s36**2 - 3072*q7a*q8a*s11**2*s12*s22*s26**2*s33**2*s36**2 + 512*q7a**2*s11**2*s22**2*s26**2*s33**2*s36**2 - 4096*q8a*q9a*s11**3*s14*s22**2*s23*s34*s36**2 -  \
    1024*q8a*q9a*s11**3*s13*s22**2*s24*s34*s36**2 - 512*q8a*q9a*s11**2*s15**2*s22**2*s24*s34*s36**2 + 6144*q8a*q9a*s11**3*s12*s22*s23*s24*s34*s36**2 - 1024*q8a*q9a*s11**2*s16**2*s22*s23*s24*s34*s36**2 -  \
    1024*q7a*q9a*s11**3*s22**2*s23*s24*s34*s36**2 + 1024*q8a*q9a*s11**3*s14*s22*s24**2*s34*s36**2 + 512*q8a*q9a*s11**2*s15*s16*s22*s24**2*s34*s36**2 - 1024*q8a*q9a*s11**3*s12*s24**3*s34*s36**2 +  \
    2560*q8a*q9a*s11**2*s14*s15*s22**2*s25*s34*s36**2 + 512*q8a*q9a*s11**2*s13*s16*s22**2*s25*s34*s36**2 - 1024*q8a*q9a*s11**2*s12*s16*s22*s23*s25*s34*s36**2 + 512*q7a*q9a*s11**2*s16*s22**2*s23*s25*s34*s36**2 -  \
    1024*q8a*q9a*s11**2*s12*s15*s22*s24*s25*s34*s36**2 - 3072*q8a*q9a*s11**2*s14*s16*s22*s24*s25*s34*s36**2 + 512*q8a*q9a*s11*s15*s16**2*s22*s24*s25*s34*s36**2 - 512*q7a*q9a*s11**2*s15*s22**2*s24*s25*s34*s36**2 +  \
    1536*q8a*q9a*s11**2*s12*s16*s24**2*s25*s34*s36**2 + 1024*q7a*q9a*s11**2*s16*s22*s24**2*s25*s34*s36**2 - 512*q8a*q9a*s11*s12*s15*s16*s22*s25**2*s34*s36**2 + 1024*q8a*q9a*s11*s14*s16**2*s22*s25**2*s34*s36**2 -  \
    2560*q7a*q9a*s11**2*s14*s22**2*s25**2*s34*s36**2 + 512*q7a*q9a*s11*s15*s16*s22**2*s25**2*s34*s36**2 - 1024*q8a*q9a*s11**2*s12**2*s24*s25**2*s34*s36**2 - 512*q8a*q9a*s11*s12*s16**2*s24*s25**2*s34*s36**2 +  \
    3072*q7a*q9a*s11**2*s12*s22*s24*s25**2*s34*s36**2 - 1024*q7a*q9a*s11*s16**2*s22*s24*s25**2*s34*s36**2 + 512*q8a*q9a*s11*s12**2*s16*s25**3*s34*s36**2 - 512*q7a*q9a*s11*s12*s16*s22*s25**3*s34*s36**2;
    v2_4 =  \
    1536*q8a*q9a*s11**2*s13*s15*s22**2*s26*s34*s36**2 - 3072*q8a*q9a*s11**2*s12*s15*s22*s23*s26*s34*s36**2 + 5120*q8a*q9a*s11**2*s14*s16*s22*s23*s26*s34*s36**2 + 1536*q7a*q9a*s11**2*s15*s22**2*s23*s26*s34*s36**2 -  \
    1024*q8a*q9a*s11**2*s13*s16*s22*s24*s26*s34*s36**2 - 1024*q8a*q9a*s11**2*s12*s16*s23*s24*s26*s34*s36**2 - 1024*q7a*q9a*s11**2*s16*s22*s23*s24*s26*s34*s36**2 - 512*q8a*q9a*s11**2*s12*s15*s24**2*s26*s34*s36**2 -  \
    512*q8a*q9a*s11**2*s14**2*s22*s25*s26*s34*s36**2 + 512*q8a*q9a*s11*s12*s15**2*s22*s25*s26*s34*s36**2 - 1536*q8a*q9a*s11*s14*s15*s16*s22*s25*s26*s34*s36**2 - 2048*q7a*q9a*s11**2*s13*s22**2*s25*s26*s34*s36**2 -  \
    512*q7a*q9a*s11*s15**2*s22**2*s25*s26*s34*s36**2 + 2048*q8a*q9a*s11**2*s12**2*s23*s25*s26*s34*s36**2 + 2048*q8a*q9a*s11**2*s12*s14*s24*s25*s26*s34*s36**2 + 2048*q7a*q9a*s11**2*s14*s22*s24*s25*s26*s34*s36**2 +  \
    512*q7a*q9a*s11*s15*s16*s22*s24*s25*s26*s34*s36**2 - 3072*q7a*q9a*s11**2*s12*s24**2*s25*s26*s34*s36**2 - 512*q8a*q9a*s11*s12**2*s15*s25**2*s26*s34*s36**2 - 512*q8a*q9a*s11*s12*s14*s16*s25**2*s26*s34*s36**2 +  \
    512*q7a*q9a*s11*s12*s15*s22*s25**2*s26*s34*s36**2 + 1536*q7a*q9a*s11*s12*s16*s24*s25**2*s26*s34*s36**2 - 1024*q8a*q9a*s11**2*s13*s14*s22*s26**2*s34*s36**2 - 3072*q8a*q9a*s11**2*s12*s14*s23*s26**2*s34*s36**2 -  \
    1024*q7a*q9a*s11**2*s14*s22*s23*s26**2*s34*s36**2 + 2048*q8a*q9a*s11**2*s12*s13*s24*s26**2*s34*s36**2 + 2048*q7a*q9a*s11**2*s12*s23*s24*s26**2*s34*s36**2 + 1024*q8a*q9a*s11*s12*s14*s15*s25*s26**2*s34*s36**2 +  \
    512*q7a*q9a*s11*s14*s15*s22*s25*s26**2*s34*s36**2 - 1024*q7a*q9a*s11*s12*s15*s24*s25*s26**2*s34*s36**2 - 512*q7a*q9a*s11*s12*s14*s25**2*s26**2*s34*s36**2 - 2048*q8a**2*s11**3*s14*s22**2*s33*s34*s36**2 +  \
    3072*q8a**2*s11**3*s12*s22*s24*s33*s34*s36**2 - 512*q8a**2*s11**2*s16**2*s22*s24*s33*s34*s36**2 - 1024*q7a*q8a*s11**3*s22**2*s24*s33*s34*s36**2 - 512*q8a**2*s11**2*s12*s16*s22*s25*s33*s34*s36**2 +  \
    512*q7a*q8a*s11**2*s16*s22**2*s25*s33*s34*s36**2 - 1536*q8a**2*s11**2*s12*s15*s22*s26*s33*s34*s36**2 + 2560*q8a**2*s11**2*s14*s16*s22*s26*s33*s34*s36**2 + 1536*q7a*q8a*s11**2*s15*s22**2*s26*s33*s34*s36**2 -  \
    512*q8a**2*s11**2*s12*s16*s24*s26*s33*s34*s36**2 - 1024*q7a*q8a*s11**2*s16*s22*s24*s26*s33*s34*s36**2 + 1024*q8a**2*s11**2*s12**2*s25*s26*s33*s34*s36**2 - 1024*q7a**2*s11**2*s22**2*s25*s26*s33*s34*s36**2 -  \
    1536*q8a**2*s11**2*s12*s14*s26**2*s33*s34*s36**2 - 1024*q7a*q8a*s11**2*s14*s22*s26**2*s33*s34*s36**2 + 2048*q7a*q8a*s11**2*s12*s24*s26**2*s33*s34*s36**2 + 1024*q8a**2*s11**3*s13*s22**2*s34**2*s36**2 +  \
    512*q8a**2*s11**2*s15**2*s22**2*s34**2*s36**2 - 3072*q8a**2*s11**3*s12*s22*s23*s34**2*s36**2 + 1536*q8a**2*s11**2*s16**2*s22*s23*s34**2*s36**2 + 2048*q7a*q8a*s11**3*s22**2*s23*s34**2*s36**2 - 512*q8a**2*s11**3*s14*s22*s24*s34**2*s36**2 -  \
    512*q8a**2*s11**2*s15*s16*s22*s24*s34**2*s36**2 + 512*q8a**2*s11**3*s12*s24**2*s34**2*s36**2 + 1024*q8a**2*s11**2*s14*s16*s22*s25*s34**2*s36**2 - 512*q8a**2*s11*s15*s16**2*s22*s25*s34**2*s36**2 -  \
    1024*q7a*q8a*s11**2*s15*s22**2*s25*s34**2*s36**2 - 512*q8a**2*s11**2*s12*s16*s24*s25*s34**2*s36**2 + 512*q8a**2*s11**2*s12**2*s25**2*s34**2*s36**2 - 1024*q7a*q8a*s11**2*s12*s22*s25**2*s34**2*s36**2 +  \
    512*q7a*q8a*s11*s16**2*s22*s25**2*s34**2*s36**2 + 1024*q7a**2*s11**2*s22**2*s25**2*s34**2*s36**2 - 512*q8a**2*s11**2*s14*s15*s22*s26*s34**2*s36**2 - 512*q8a**2*s11**2*s13*s16*s22*s26*s34**2*s36**2 -  \
    1536*q8a**2*s11**2*s12*s16*s23*s26*s34**2*s36**2 - 1024*q7a*q8a*s11**2*s16*s22*s23*s26*s34**2*s36**2 + 1024*q8a**2*s11**2*s12*s15*s24*s26*s34**2*s36**2 - 512*q8a**2*s11**2*s12*s14*s25*s26*s34**2*s36**2 +  \
    512*q8a**2*s11*s12*s15*s16*s25*s26*s34**2*s36**2 + 512*q7a*q8a*s11*s15*s16*s22*s25*s26*s34**2*s36**2 - 512*q7a*q8a*s11*s12*s16*s25**2*s26*s34**2*s36**2 - 512*q7a**2*s11*s16*s22*s25**2*s26*s34**2*s36**2 +  \
    512*q8a**2*s11**2*s12*s13*s26**2*s34**2*s36**2 + 1024*q7a*q8a*s11**2*s12*s23*s26**2*s34**2*s36**2 - 512*q7a*q8a*s11*s12*s15*s25*s26**2*s34**2*s36**2 + 512*q7a**2*s11*s12*s25**2*s26**2*s34**2*s36**2 -  \
    1024*q8a*q9a*s11**2*s13*s15*s22**3*s35*s36**2 + 2048*q8a*q9a*s11**2*s12*s15*s22**2*s23*s35*s36**2 - 2048*q8a*q9a*s11**2*s14*s16*s22**2*s23*s35*s36**2 - 1024*q7a*q9a*s11**2*s15*s22**3*s23*s35*s36**2 -  \
    512*q8a*q9a*s11**2*s14*s15*s22**2*s24*s35*s36**2 + 1536*q8a*q9a*s11**2*s13*s16*s22**2*s24*s35*s36**2 - 1024*q8a*q9a*s11**2*s12*s16*s22*s23*s24*s35*s36**2 + 1536*q7a*q9a*s11**2*s16*s22**2*s23*s24*s35*s36**2 +  \
    1024*q8a*q9a*s11**2*s12*s15*s22*s24**2*s35*s36**2 - 512*q8a*q9a*s11**2*s14*s16*s22*s24**2*s35*s36**2 - 512*q7a*q9a*s11**2*s15*s22**2*s24**2*s35*s36**2 + 512*q8a*q9a*s11**2*s12*s16*s24**3*s35*s36**2 +  \
    5120*q8a*q9a*s11**2*s12*s13*s22**2*s25*s35*s36**2 - 2560*q8a*q9a*s11**2*s14**2*s22**2*s25*s35*s36**2 - 1024*q8a*q9a*s11*s12*s15**2*s22**2*s25*s35*s36**2 + 2048*q8a*q9a*s11*s14*s15*s16*s22**2*s25*s35*s36**2 -  \
    1024*q8a*q9a*s11*s13*s16**2*s22**2*s25*s35*s36**2 - 4096*q7a*q9a*s11**2*s13*s22**3*s25*s35*s36**2 + 1024*q7a*q9a*s11*s15**2*s22**3*s25*s35*s36**2 - 6144*q8a*q9a*s11**2*s12**2*s22*s23*s25*s35*s36**2 +  \
    2048*q8a*q9a*s11*s12*s16**2*s22*s23*s25*s35*s36**2 + 5120*q7a*q9a*s11**2*s12*s22**2*s23*s25*s35*s36**2 - 1024*q7a*q9a*s11*s16**2*s22**2*s23*s25*s35*s36**2 + 3072*q8a*q9a*s11**2*s12*s14*s22*s24*s25*s35*s36**2 -  \
    2048*q8a*q9a*s11*s12*s15*s16*s22*s24*s25*s35*s36**2 + 512*q8a*q9a*s11*s14*s16**2*s22*s24*s25*s35*s36**2 + 2560*q7a*q9a*s11**2*s14*s22**2*s24*s25*s35*s36**2 - 1024*q8a*q9a*s11**2*s12**2*s24**2*s25*s35*s36**2 -  \
    512*q8a*q9a*s11*s12*s16**2*s24**2*s25*s35*s36**2 - 2048*q7a*q9a*s11**2*s12*s22*s24**2*s25*s35*s36**2 + 2048*q8a*q9a*s11*s12**2*s15*s22*s25**2*s35*s36**2 - 1536*q8a*q9a*s11*s12*s14*s16*s22*s25**2*s35*s36**2 -  \
    2048*q7a*q9a*s11*s12*s15*s22**2*s25**2*s35*s36**2 - 512*q7a*q9a*s11*s14*s16*s22**2*s25**2*s35*s36**2 + 1536*q8a*q9a*s11*s12**2*s16*s24*s25**2*s35*s36**2 + 512*q7a*q9a*s11*s12*s16*s22*s24*s25**2*s35*s36**2 -  \
    1024*q8a*q9a*s11*s12**3*s25**3*s35*s36**2 + 1024*q7a*q9a*s11*s12**2*s22*s25**3*s35*s36**2 + 512*q8a*q9a*s11**2*s13*s14*s22**2*s26*s35*s36**2 + 1024*q8a*q9a*s11*s14*s15**2*s22**2*s26*s35*s36**2 -  \
    1024*q8a*q9a*s11*s13*s15*s16*s22**2*s26*s35*s36**2 + 1024*q8a*q9a*s11**2*s12*s14*s22*s23*s26*s35*s36**2 + 2048*q8a*q9a*s11*s12*s15*s16*s22*s23*s26*s35*s36**2 + 512*q7a*q9a*s11**2*s14*s22**2*s23*s26*s35*s36**2 -  \
    1024*q7a*q9a*s11*s15*s16*s22**2*s23*s26*s35*s36**2 - 2048*q8a*q9a*s11**2*s12*s13*s22*s24*s26*s35*s36**2 + 512*q8a*q9a*s11**2*s14**2*s22*s24*s26*s35*s36**2 - 1536*q8a*q9a*s11*s12*s15**2*s22*s24*s26*s35*s36**2 +  \
    512*q8a*q9a*s11*s14*s15*s16*s22*s24*s26*s35*s36**2 + 512*q7a*q9a*s11*s15**2*s22**2*s24*s26*s35*s36**2 + 2048*q8a*q9a*s11**2*s12**2*s23*s24*s26*s35*s36**2 - 2048*q7a*q9a*s11**2*s12*s22*s23*s24*s26*s35*s36**2 -  \
    512*q8a*q9a*s11**2*s12*s14*s24**2*s26*s35*s36**2 - 512*q8a*q9a*s11*s12*s15*s16*s24**2*s26*s35*s36**2 - 1024*q8a*q9a*s11*s12*s14*s15*s22*s25*s26*s35*s36**2 - 2048*q8a*q9a*s11*s12*s13*s16*s22*s25*s26*s35*s36**2 -  \
    3072*q7a*q9a*s11*s14*s15*s22**2*s25*s26*s35*s36**2 + 5120*q7a*q9a*s11*s13*s16*s22**2*s25*s26*s35*s36**2 - 1024*q8a*q9a*s11*s12**2*s16*s23*s25*s26*s35*s36**2 - 2048*q7a*q9a*s11*s12*s16*s22*s23*s25*s26*s35*s36**2 +  \
    2048*q8a*q9a*s11*s12**2*s15*s24*s25*s26*s35*s36**2 + 2048*q7a*q9a*s11*s12*s15*s22*s24*s25*s26*s35*s36**2 - 1536*q7a*q9a*s11*s14*s16*s22*s24*s25*s26*s35*s36**2 + 1536*q7a*q9a*s11*s12*s16*s24**2*s25*s26*s35*s36**2 -  \
    512*q8a*q9a*s11*s12**2*s14*s25**2*s26*s35*s36**2 + 3584*q7a*q9a*s11*s12*s14*s22*s25**2*s26*s35*s36**2 - 3072*q7a*q9a*s11*s12**2*s24*s25**2*s26*s35*s36**2 + 2048*q8a*q9a*s11*s12*s13*s15*s22*s26**2*s35*s36**2 -  \
    1024*q8a*q9a*s11*s14**2*s15*s22*s26**2*s35*s36**2 - 1024*q7a*q9a*s11*s13*s15*s22**2*s26**2*s35*s36**2 - 3072*q8a*q9a*s11*s12**2*s15*s23*s26**2*s35*s36**2 + 2048*q7a*q9a*s11*s12*s15*s22*s23*s26**2*s35*s36**2 +  \
    1024*q8a*q9a*s11*s12*s14*s15*s24*s26**2*s35*s36**2 + 512*q7a*q9a*s11*s14*s15*s22*s24*s26**2*s35*s36**2 - 512*q7a*q9a*s11*s12*s15*s24**2*s26**2*s35*s36**2 + 2048*q8a*q9a*s11*s12**2*s13*s25*s26**2*s35*s36**2 -  \
    4096*q7a*q9a*s11*s12*s13*s22*s25*s26**2*s35*s36**2 + 1024*q7a*q9a*s11*s14**2*s22*s25*s26**2*s35*s36**2 + 2048*q7a*q9a*s11*s12**2*s23*s25*s26**2*s35*s36**2 - 1024*q7a*q9a*s11*s12*s14*s24*s25*s26**2*s35*s36**2 +  \
    1024*q8a**2*s11**2*s12*s15*s22**2*s33*s35*s36**2 - 1024*q8a**2*s11**2*s14*s16*s22**2*s33*s35*s36**2 - 1024*q7a*q8a*s11**2*s15*s22**3*s33*s35*s36**2 - 512*q8a**2*s11**2*s12*s16*s22*s24*s33*s35*s36**2 +  \
    1536*q7a*q8a*s11**2*s16*s22**2*s24*s33*s35*s36**2 - 3072*q8a**2*s11**2*s12**2*s22*s25*s33*s35*s36**2 + 1024*q8a**2*s11*s12*s16**2*s22*s25*s33*s35*s36**2 + 5120*q7a*q8a*s11**2*s12*s22**2*s25*s33*s35*s36**2 -  \
    1024*q7a*q8a*s11*s16**2*s22**2*s25*s33*s35*s36**2 - 2048*q7a**2*s11**2*s22**3*s25*s33*s35*s36**2 + 512*q8a**2*s11**2*s12*s14*s22*s26*s33*s35*s36**2 + 1024*q8a**2*s11*s12*s15*s16*s22*s26*s33*s35*s36**2 +  \
    512*q7a*q8a*s11**2*s14*s22**2*s26*s33*s35*s36**2 - 1024*q7a*q8a*s11*s15*s16*s22**2*s26*s33*s35*s36**2 + 1024*q8a**2*s11**2*s12**2*s24*s26*s33*s35*s36**2 - 2048*q7a*q8a*s11**2*s12*s22*s24*s26*s33*s35*s36**2 -  \
    512*q8a**2*s11*s12**2*s16*s25*s26*s33*s35*s36**2 - 2048*q7a*q8a*s11*s12*s16*s22*s25*s26*s33*s35*s36**2 + 2560*q7a**2*s11*s16*s22**2*s25*s26*s33*s35*s36**2 - 1536*q8a**2*s11*s12**2*s15*s26**2*s33*s35*s36**2 +  \
    2048*q7a*q8a*s11*s12*s15*s22*s26**2*s33*s35*s36**2 - 512*q7a**2*s11*s15*s22**2*s26**2*s33*s35*s36**2 + 2048*q7a*q8a*s11*s12**2*s25*s26**2*s33*s35*s36**2 - 2048*q7a**2*s11*s12*s22*s25*s26**2*s33*s35*s36**2 -  \
    1024*q8a**2*s11**2*s14*s15*s22**2*s34*s35*s36**2 + 1024*q8a**2*s11**2*s13*s16*s22**2*s34*s35*s36**2 - 3072*q8a**2*s11**2*s12*s16*s22*s23*s34*s35*s36**2 + 2048*q7a*q8a*s11**2*s16*s22**2*s23*s34*s35*s36**2 -  \
    512*q8a**2*s11**2*s12*s15*s22*s24*s34*s35*s36**2 + 1536*q8a**2*s11**2*s14*s16*s22*s24*s34*s35*s36**2 + 1536*q7a*q8a*s11**2*s15*s22**2*s24*s34*s35*s36**2 - 512*q8a**2*s11**2*s12*s16*s24**2*s34*s35*s36**2 -  \
    1024*q7a*q8a*s11**2*s16*s22*s24**2*s34*s35*s36**2 + 512*q8a**2*s11**2*s12*s14*s22*s25*s34*s35*s36**2 + 2048*q8a**2*s11*s12*s15*s16*s22*s25*s34*s35*s36**2 - 1024*q8a**2*s11*s14*s16**2*s22*s25*s34*s35*s36**2 +  \
    512*q7a*q8a*s11**2*s14*s22**2*s25*s34*s35*s36**2 - 2048*q7a*q8a*s11*s15*s16*s22**2*s25*s34*s35*s36**2 + 512*q8a**2*s11*s12*s16**2*s24*s25*s34*s35*s36**2 + 512*q7a*q8a*s11*s16**2*s22*s24*s25*s34*s35*s36**2 -  \
    1024*q7a**2*s11**2*s22**2*s24*s25*s34*s35*s36**2 - 512*q8a**2*s11*s12**2*s16*s25**2*s34*s35*s36**2 - 1024*q7a*q8a*s11*s12*s16*s22*s25**2*s34*s35*s36**2 + 1536*q7a**2*s11*s16*s22**2*s25**2*s34*s35*s36**2 -  \
    3072*q8a**2*s11**2*s12*s13*s22*s26*s34*s35*s36**2 + 1536*q8a**2*s11**2*s14**2*s22*s26*s34*s35*s36**2 + 1024*q8a**2*s11*s12*s15**2*s22*s26*s34*s35*s36**2 - 1024*q8a**2*s11*s14*s15*s16*s22*s26*s34*s35*s36**2 +  \
    2048*q7a*q8a*s11**2*s13*s22**2*s26*s34*s35*s36**2 - 1024*q7a*q8a*s11*s15**2*s22**2*s26*s34*s35*s36**2 + 6144*q8a**2*s11**2*s12**2*s23*s26*s34*s35*s36**2 - 6144*q7a*q8a*s11**2*s12*s22*s23*s26*s34*s35*s36**2 +  \
    1024*q7a**2*s11**2*s22**2*s23*s26*s34*s35*s36**2 - 2560*q8a**2*s11**2*s12*s14*s24*s26*s34*s35*s36**2 + 512*q8a**2*s11*s12*s15*s16*s24*s26*s34*s35*s36**2 - 2048*q7a*q8a*s11**2*s14*s22*s24*s26*s34*s35*s36**2 +  \
    512*q7a*q8a*s11*s15*s16*s22*s24*s26*s34*s35*s36**2 + 3072*q7a*q8a*s11**2*s12*s24**2*s26*s34*s35*s36**2 - 2560*q8a**2*s11*s12**2*s15*s25*s26*s34*s35*s36**2 + 512*q8a**2*s11*s12*s14*s16*s25*s26*s34*s35*s36**2 +  \
    1024*q7a*q8a*s11*s12*s15*s22*s25*s26*s34*s35*s36**2 + 2560*q7a*q8a*s11*s14*s16*s22*s25*s26*s34*s35*s36**2 + 1536*q7a**2*s11*s15*s22**2*s25*s26*s34*s35*s36**2 - 2048*q7a*q8a*s11*s12*s16*s24*s25*s26*s34*s35*s36**2 -  \
    1024*q7a**2*s11*s16*s22*s24*s25*s26*s34*s35*s36**2 + 3072*q7a*q8a*s11*s12**2*s25**2*s26*s34*s35*s36**2 - 3072*q7a**2*s11*s12*s22*s25**2*s26*s34*s35*s36**2 + 512*q8a**2*s11*s12*s14*s15*s26**2*s34*s35*s36**2 +  \
    512*q7a*q8a*s11*s14*s15*s22*s26**2*s34*s35*s36**2 - 1024*q7a*q8a*s11*s12*s15*s24*s26**2*s34*s35*s36**2 - 1024*q7a*q8a*s11*s12*s14*s25*s26**2*s34*s35*s36**2 - 1024*q7a**2*s11*s14*s22*s25*s26**2*s34*s35*s36**2 +  \
    2048*q7a**2*s11*s12*s24*s25*s26**2*s34*s35*s36**2 - 2048*q8a**2*s11**2*s12*s13*s22**2*s35**2*s36**2 + 1024*q8a**2*s11**2*s14**2*s22**2*s35**2*s36**2 + 2048*q7a*q8a*s11**2*s13*s22**3*s35**2*s36**2 +  \
    3072*q8a**2*s11**2*s12**2*s22*s23*s35**2*s36**2 - 4096*q7a*q8a*s11**2*s12*s22**2*s23*s35**2*s36**2 + 1024*q7a**2*s11**2*s22**3*s23*s35**2*s36**2 - 1024*q8a**2*s11**2*s12*s14*s22*s24*s35**2*s36**2 -  \
    1024*q7a*q8a*s11**2*s14*s22**2*s24*s35**2*s36**2 + 512*q8a**2*s11**2*s12**2*s24**2*s35**2*s36**2 + 512*q7a**2*s11**2*s22**2*s24**2*s35**2*s36**2 - 512*q8a**2*s11*s12**2*s15*s22*s25*s35**2*s36**2;
    v2_5= \
    1024*q7a*q8a*s11*s12*s15*s22**2*s25*s35**2*s36**2 - 512*q7a**2*s11*s15*s22**3*s25*s35**2*s36**2 - 512*q8a**2*s11*s12**2*s16*s24*s25*s35**2*s36**2 + 1024*q7a*q8a*s11*s12*s16*s22*s24*s25*s35**2*s36**2 -  \
    512*q7a**2*s11*s16*s22**2*s24*s25*s35**2*s36**2 + 512*q8a**2*s11*s12**3*s25**2*s35**2*s36**2 - 1024*q7a*q8a*s11*s12**2*s22*s25**2*s35**2*s36**2 + 512*q7a**2*s11*s12*s22**2*s25**2*s35**2*s36**2 +  \
    1024*q8a**2*s11*s12*s13*s16*s22*s26*s35**2*s36**2 - 512*q8a**2*s11*s14**2*s16*s22*s26*s35**2*s36**2 - 1024*q7a*q8a*s11*s13*s16*s22**2*s26*s35**2*s36**2 - 1536*q8a**2*s11*s12**2*s16*s23*s26*s35**2*s36**2 +  \
    2048*q7a*q8a*s11*s12*s16*s22*s23*s26*s35**2*s36**2 - 512*q7a**2*s11*s16*s22**2*s23*s26*s35**2*s36**2 - 512*q8a**2*s11*s12**2*s15*s24*s26*s35**2*s36**2 + 512*q8a**2*s11*s12*s14*s16*s24*s26*s35**2*s36**2 +  \
    1024*q7a*q8a*s11*s12*s15*s22*s24*s26*s35**2*s36**2 + 512*q7a*q8a*s11*s14*s16*s22*s24*s26*s35**2*s36**2 - 512*q7a**2*s11*s15*s22**2*s24*s26*s35**2*s36**2 - 512*q7a*q8a*s11*s12*s16*s24**2*s26*s35**2*s36**2 +  \
    1024*q8a**2*s11*s12**2*s14*s25*s26*s35**2*s36**2 - 2048*q7a*q8a*s11*s12*s14*s22*s25*s26*s35**2*s36**2 + 1024*q7a**2*s11*s14*s22**2*s25*s26*s35**2*s36**2 + 512*q8a**2*s11*s12**2*s13*s26**2*s35**2*s36**2 -  \
    2048*q7a*q8a*s11*s12*s13*s22*s26**2*s35**2*s36**2 + 512*q7a*q8a*s11*s14**2*s22*s26**2*s35**2*s36**2 + 1536*q7a**2*s11*s13*s22**2*s26**2*s35**2*s36**2 + 1024*q7a*q8a*s11*s12**2*s23*s26**2*s35**2*s36**2 -  \
    1024*q7a**2*s11*s12*s22*s23*s26**2*s35**2*s36**2 - 512*q7a*q8a*s11*s12*s14*s24*s26**2*s35**2*s36**2 - 512*q7a**2*s11*s14*s22*s24*s26**2*s35**2*s36**2 + 512*q7a**2*s11*s12*s24**2*s26**2*s35**2*s36**2 +  \
    2048*q8a*q9a*s11**2*s14*s15*s22**2*s23*s36**3 - 2048*q8a*q9a*s11**2*s13*s16*s22**2*s23*s36**3 + 3072*q8a*q9a*s11**2*s12*s16*s22*s23**2*s36**3 - 1024*q7a*q9a*s11**2*s16*s22**2*s23**2*s36**3 +  \
    512*q8a*q9a*s11**2*s13*s15*s22**2*s24*s36**3 - 3072*q8a*q9a*s11**2*s12*s15*s22*s23*s24*s36**3 + 1024*q8a*q9a*s11**2*s14*s16*s22*s23*s24*s36**3 + 512*q7a*q9a*s11**2*s15*s22**2*s23*s24*s36**3 -  \
    512*q8a*q9a*s11**2*s14*s15*s22*s24**2*s36**3 - 1024*q8a*q9a*s11**2*s12*s16*s23*s24**2*s36**3 + 512*q8a*q9a*s11**2*s12*s15*s24**3*s36**3 - 512*q8a*q9a*s11**2*s13*s14*s22**2*s25*s36**3 -  \
    1024*q8a*q9a*s11*s14*s15**2*s22**2*s25*s36**3 + 1024*q8a*q9a*s11*s13*s15*s16*s22**2*s25*s36**3 - 1024*q8a*q9a*s11**2*s12*s14*s22*s23*s25*s36**3 - 2048*q8a*q9a*s11*s12*s15*s16*s22*s23*s25*s36**3 -  \
    512*q7a*q9a*s11**2*s14*s22**2*s23*s25*s36**3 + 1024*q7a*q9a*s11*s15*s16*s22**2*s23*s25*s36**3 - 2048*q8a*q9a*s11**2*s12*s13*s22*s24*s25*s36**3 + 1536*q8a*q9a*s11**2*s14**2*s22*s24*s25*s36**3 +  \
    1536*q8a*q9a*s11*s12*s15**2*s22*s24*s25*s36**3 - 512*q8a*q9a*s11*s14*s15*s16*s22*s24*s25*s36**3 + 2048*q7a*q9a*s11**2*s13*s22**2*s24*s25*s36**3 - 512*q7a*q9a*s11*s15**2*s22**2*s24*s25*s36**3 +  \
    4096*q8a*q9a*s11**2*s12**2*s23*s24*s25*s36**3 - 2048*q7a*q9a*s11**2*s12*s22*s23*s24*s25*s36**3 - 1536*q8a*q9a*s11**2*s12*s14*s24**2*s25*s36**3 + 512*q8a*q9a*s11*s12*s15*s16*s24**2*s25*s36**3 -  \
    1024*q7a*q9a*s11**2*s14*s22*s24**2*s25*s36**3 + 1024*q7a*q9a*s11**2*s12*s24**3*s25*s36**3 + 512*q8a*q9a*s11*s12*s14*s15*s22*s25**2*s36**3 + 2048*q8a*q9a*s11*s12*s13*s16*s22*s25**2*s36**3 -  \
    512*q8a*q9a*s11*s14**2*s16*s22*s25**2*s36**3 + 1536*q7a*q9a*s11*s14*s15*s22**2*s25**2*s36**3 - 3072*q7a*q9a*s11*s13*s16*s22**2*s25**2*s36**3 - 1024*q8a*q9a*s11*s12**2*s16*s23*s25**2*s36**3 +  \
    2048*q7a*q9a*s11*s12*s16*s22*s23*s25**2*s36**3 - 1536*q8a*q9a*s11*s12**2*s15*s24*s25**2*s36**3 + 512*q8a*q9a*s11*s12*s14*s16*s24*s25**2*s36**3 - 512*q7a*q9a*s11*s12*s15*s22*s24*s25**2*s36**3 +  \
    1024*q7a*q9a*s11*s14*s16*s22*s24*s25**2*s36**3 - 1024*q7a*q9a*s11*s12*s16*s24**2*s25**2*s36**3 + 512*q8a*q9a*s11*s12**2*s14*s25**3*s36**3 - 1536*q7a*q9a*s11*s12*s14*s22*s25**3*s36**3 + 1024*q7a*q9a*s11*s12**2*s24*s25**3*s36**3 -  \
    1024*q8a*q9a*s11**2*s13**2*s22**2*s26*s36**3 + 6144*q8a*q9a*s11**2*s12*s13*s22*s23*s26*s36**3 - 3072*q8a*q9a*s11**2*s14**2*s22*s23*s26*s36**3 - 2048*q7a*q9a*s11**2*s13*s22**2*s23*s26*s36**3 -  \
    6144*q8a*q9a*s11**2*s12**2*s23**2*s26*s36**3 + 3072*q7a*q9a*s11**2*s12*s22*s23**2*s26*s36**3 + 1024*q8a*q9a*s11**2*s13*s14*s22*s24*s26*s36**3 + 3072*q8a*q9a*s11**2*s12*s14*s23*s24*s26*s36**3 +  \
    1024*q7a*q9a*s11**2*s14*s22*s23*s24*s26*s36**3 - 1024*q8a*q9a*s11**2*s12*s13*s24**2*s26*s36**3 - 1024*q7a*q9a*s11**2*s12*s23*s24**2*s26*s36**3 - 2048*q8a*q9a*s11*s12*s13*s15*s22*s25*s26*s36**3 +  \
    1024*q8a*q9a*s11*s14**2*s15*s22*s25*s26*s36**3 + 1024*q7a*q9a*s11*s13*s15*s22**2*s25*s26*s36**3 + 3072*q8a*q9a*s11*s12**2*s15*s23*s25*s26*s36**3 - 2048*q7a*q9a*s11*s12*s15*s22*s23*s25*s26*s36**3 -  \
    1024*q8a*q9a*s11*s12*s14*s15*s24*s25*s26*s36**3 - 512*q7a*q9a*s11*s14*s15*s22*s24*s25*s26*s36**3 + 512*q7a*q9a*s11*s12*s15*s24**2*s25*s26*s36**3 - 1024*q8a*q9a*s11*s12**2*s13*s25**2*s26*s36**3 +  \
    2048*q7a*q9a*s11*s12*s13*s22*s25**2*s26*s36**3 - 512*q7a*q9a*s11*s14**2*s22*s25**2*s26*s36**3 - 1024*q7a*q9a*s11*s12**2*s23*s25**2*s26*s36**3 + 512*q7a*q9a*s11*s12*s14*s24*s25**2*s26*s36**3 +  \
    1024*q8a**2*s11**2*s14*s15*s22**2*s33*s36**3 - 1024*q8a**2*s11**2*s13*s16*s22**2*s33*s36**3 + 3072*q8a**2*s11**2*s12*s16*s22*s23*s33*s36**3 - 2048*q7a*q8a*s11**2*s16*s22**2*s23*s33*s36**3 -  \
    1536*q8a**2*s11**2*s12*s15*s22*s24*s33*s36**3 + 512*q8a**2*s11**2*s14*s16*s22*s24*s33*s36**3 + 512*q7a*q8a*s11**2*s15*s22**2*s24*s33*s36**3 - 512*q8a**2*s11**2*s12*s16*s24**2*s33*s36**3 -  \
    512*q8a**2*s11**2*s12*s14*s22*s25*s33*s36**3 - 1024*q8a**2*s11*s12*s15*s16*s22*s25*s33*s36**3 - 512*q7a*q8a*s11**2*s14*s22**2*s25*s33*s36**3 + 1024*q7a*q8a*s11*s15*s16*s22**2*s25*s33*s36**3 +  \
    2048*q8a**2*s11**2*s12**2*s24*s25*s33*s36**3 - 2048*q7a*q8a*s11**2*s12*s22*s24*s25*s33*s36**3 + 1024*q7a**2*s11**2*s22**2*s24*s25*s33*s36**3 - 512*q8a**2*s11*s12**2*s16*s25**2*s33*s36**3 +  \
    2048*q7a*q8a*s11*s12*s16*s22*s25**2*s33*s36**3 - 1536*q7a**2*s11*s16*s22**2*s25**2*s33*s36**3 + 3072*q8a**2*s11**2*s12*s13*s22*s26*s33*s36**3 - 1536*q8a**2*s11**2*s14**2*s22*s26*s33*s36**3 -  \
    2048*q7a*q8a*s11**2*s13*s22**2*s26*s33*s36**3 - 6144*q8a**2*s11**2*s12**2*s23*s26*s33*s36**3 + 6144*q7a*q8a*s11**2*s12*s22*s23*s26*s33*s36**3 - 1024*q7a**2*s11**2*s22**2*s23*s26*s33*s36**3 +  \
    1536*q8a**2*s11**2*s12*s14*s24*s26*s33*s36**3 + 1024*q7a*q8a*s11**2*s14*s22*s24*s26*s33*s36**3 - 1024*q7a*q8a*s11**2*s12*s24**2*s26*s33*s36**3 + 1536*q8a**2*s11*s12**2*s15*s25*s26*s33*s36**3 -  \
    2048*q7a*q8a*s11*s12*s15*s22*s25*s26*s33*s36**3 + 512*q7a**2*s11*s15*s22**2*s25*s26*s33*s36**3 - 1024*q7a*q8a*s11*s12**2*s25**2*s26*s33*s36**3 + 1024*q7a**2*s11*s12*s22*s25**2*s26*s33*s36**3 -  \
    1024*q8a**2*s11**2*s13*s15*s22**2*s34*s36**3 + 3072*q8a**2*s11**2*s12*s15*s22*s23*s34*s36**3 - 3072*q8a**2*s11**2*s14*s16*s22*s23*s34*s36**3 - 2048*q7a*q8a*s11**2*s15*s22**2*s23*s34*s36**3 +  \
    512*q8a**2*s11**2*s14*s15*s22*s24*s34*s36**3 + 512*q8a**2*s11**2*s13*s16*s22*s24*s34*s36**3 + 1536*q8a**2*s11**2*s12*s16*s23*s24*s34*s36**3 + 1024*q7a*q8a*s11**2*s16*s22*s23*s24*s34*s36**3 -  \
    512*q8a**2*s11**2*s12*s15*s24**2*s34*s36**3 + 1024*q8a**2*s11**2*s12*s13*s22*s25*s34*s36**3 - 512*q8a**2*s11**2*s14**2*s22*s25*s34*s36**3 - 1024*q8a**2*s11*s12*s15**2*s22*s25*s34*s36**3 +  \
    1024*q8a**2*s11*s14*s15*s16*s22*s25*s34*s36**3 + 1024*q7a*q8a*s11*s15**2*s22**2*s25*s34*s36**3 - 3072*q8a**2*s11**2*s12**2*s23*s25*s34*s36**3 + 2048*q7a*q8a*s11**2*s12*s22*s23*s25*s34*s36**3 +  \
    512*q8a**2*s11**2*s12*s14*s24*s25*s34*s36**3 - 512*q8a**2*s11*s12*s15*s16*s24*s25*s34*s36**3 - 512*q7a*q8a*s11*s15*s16*s22*s24*s25*s34*s36**3 + 1024*q8a**2*s11*s12**2*s15*s25**2*s34*s36**3 -  \
    1024*q7a*q8a*s11*s14*s16*s22*s25**2*s34*s36**3 - 1024*q7a**2*s11*s15*s22**2*s25**2*s34*s36**3 + 512*q7a*q8a*s11*s12*s16*s24*s25**2*s34*s36**3 + 512*q7a**2*s11*s16*s22*s24*s25**2*s34*s36**3 -  \
    1024*q7a*q8a*s11*s12**2*s25**3*s34*s36**3 + 1024*q7a**2*s11*s12*s22*s25**3*s34*s36**3 + 512*q8a**2*s11**2*s13*s14*s22*s26*s34*s36**3 + 1536*q8a**2*s11**2*s12*s14*s23*s26*s34*s36**3 +  \
    1024*q7a*q8a*s11**2*s14*s22*s23*s26*s34*s36**3 - 1024*q8a**2*s11**2*s12*s13*s24*s26*s34*s36**3 - 2048*q7a*q8a*s11**2*s12*s23*s24*s26*s34*s36**3 - 512*q8a**2*s11*s12*s14*s15*s25*s26*s34*s36**3 -  \
    512*q7a*q8a*s11*s14*s15*s22*s25*s26*s34*s36**3 + 1024*q7a*q8a*s11*s12*s15*s24*s25*s26*s34*s36**3 + 512*q7a*q8a*s11*s12*s14*s25**2*s26*s34*s36**3 + 512*q7a**2*s11*s14*s22*s25**2*s26*s34*s36**3 -  \
    1024*q7a**2*s11*s12*s24*s25**2*s26*s34*s36**3 + 2048*q8a**2*s11**2*s12*s13*s22*s24*s35*s36**3 - 1024*q8a**2*s11**2*s14**2*s22*s24*s35*s36**3 - 2048*q7a*q8a*s11**2*s13*s22**2*s24*s35*s36**3 -  \
    3072*q8a**2*s11**2*s12**2*s23*s24*s35*s36**3 + 4096*q7a*q8a*s11**2*s12*s22*s23*s24*s35*s36**3 - 1024*q7a**2*s11**2*s22**2*s23*s24*s35*s36**3 + 1024*q8a**2*s11**2*s12*s14*s24**2*s35*s36**3 +  \
    1024*q7a*q8a*s11**2*s14*s22*s24**2*s35*s36**3 - 1024*q7a*q8a*s11**2*s12*s24**3*s35*s36**3 - 1024*q8a**2*s11*s12*s13*s16*s22*s25*s35*s36**3 + 512*q8a**2*s11*s14**2*s16*s22*s25*s35*s36**3 +  \
    1024*q7a*q8a*s11*s13*s16*s22**2*s25*s35*s36**3 + 1536*q8a**2*s11*s12**2*s16*s23*s25*s35*s36**3 - 2048*q7a*q8a*s11*s12*s16*s22*s23*s25*s35*s36**3 + 512*q7a**2*s11*s16*s22**2*s23*s25*s35*s36**3 +  \
    512*q8a**2*s11*s12**2*s15*s24*s25*s35*s36**3 - 512*q8a**2*s11*s12*s14*s16*s24*s25*s35*s36**3 - 1024*q7a*q8a*s11*s12*s15*s22*s24*s25*s35*s36**3 - 512*q7a*q8a*s11*s14*s16*s22*s24*s25*s35*s36**3 +  \
    512*q7a**2*s11*s15*s22**2*s24*s25*s35*s36**3 + 512*q7a*q8a*s11*s12*s16*s24**2*s25*s35*s36**3 - 512*q8a**2*s11*s12**2*s14*s25**2*s35*s36**3 + 1024*q7a*q8a*s11*s12*s14*s22*s25**2*s35*s36**3 -  \
    512*q7a**2*s11*s14*s22**2*s25**2*s35*s36**3 - 1024*q8a**2*s11*s12*s13*s15*s22*s26*s35*s36**3 + 512*q8a**2*s11*s14**2*s15*s22*s26*s35*s36**3 + 1024*q7a*q8a*s11*s13*s15*s22**2*s26*s35*s36**3 +  \
    1536*q8a**2*s11*s12**2*s15*s23*s26*s35*s36**3 - 2048*q7a*q8a*s11*s12*s15*s22*s23*s26*s35*s36**3 + 512*q7a**2*s11*s15*s22**2*s23*s26*s35*s36**3 - 512*q8a**2*s11*s12*s14*s15*s24*s26*s35*s36**3 -  \
    512*q7a*q8a*s11*s14*s15*s22*s24*s26*s35*s36**3 + 512*q7a*q8a*s11*s12*s15*s24**2*s26*s35*s36**3 - 1024*q8a**2*s11*s12**2*s13*s25*s26*s35*s36**3 + 4096*q7a*q8a*s11*s12*s13*s22*s25*s26*s35*s36**3 -  \
    1024*q7a*q8a*s11*s14**2*s22*s25*s26*s35*s36**3 - 3072*q7a**2*s11*s13*s22**2*s25*s26*s35*s36**3 - 2048*q7a*q8a*s11*s12**2*s23*s25*s26*s35*s36**3 + 2048*q7a**2*s11*s12*s22*s23*s25*s26*s35*s36**3 +  \
    1024*q7a*q8a*s11*s12*s14*s24*s25*s26*s35*s36**3 + 1024*q7a**2*s11*s14*s22*s24*s25*s26*s35*s36**3 - 1024*q7a**2*s11*s12*s24**2*s25*s26*s35*s36**3 + 512*q8a**2*s11**2*s13**2*s22**2*s36**4 -  \
    3072*q8a**2*s11**2*s12*s13*s22*s23*s36**4 + 1536*q8a**2*s11**2*s14**2*s22*s23*s36**4 + 2048*q7a*q8a*s11**2*s13*s22**2*s23*s36**4 + 3072*q8a**2*s11**2*s12**2*s23**2*s36**4 - 3072*q7a*q8a*s11**2*s12*s22*s23**2*s36**4 +  \
    512*q7a**2*s11**2*s22**2*s23**2*s36**4 - 512*q8a**2*s11**2*s13*s14*s22*s24*s36**4 - 1536*q8a**2*s11**2*s12*s14*s23*s24*s36**4 - 1024*q7a*q8a*s11**2*s14*s22*s23*s24*s36**4 + 512*q8a**2*s11**2*s12*s13*s24**2*s36**4 +  \
    1024*q7a*q8a*s11**2*s12*s23*s24**2*s36**4 + 1024*q8a**2*s11*s12*s13*s15*s22*s25*s36**4 - 512*q8a**2*s11*s14**2*s15*s22*s25*s36**4 - 1024*q7a*q8a*s11*s13*s15*s22**2*s25*s36**4 - 1536*q8a**2*s11*s12**2*s15*s23*s25*s36**4 +  \
    2048*q7a*q8a*s11*s12*s15*s22*s23*s25*s36**4 - 512*q7a**2*s11*s15*s22**2*s23*s25*s36**4 + 512*q8a**2*s11*s12*s14*s15*s24*s25*s36**4 + 512*q7a*q8a*s11*s14*s15*s22*s24*s25*s36**4 - 512*q7a*q8a*s11*s12*s15*s24**2*s25*s36**4 +  \
    512*q8a**2*s11*s12**2*s13*s25**2*s36**4 - 2048*q7a*q8a*s11*s12*s13*s22*s25**2*s36**4 + 512*q7a*q8a*s11*s14**2*s22*s25**2*s36**4 + 1536*q7a**2*s11*s13*s22**2*s25**2*s36**4 + 1024*q7a*q8a*s11*s12**2*s23*s25**2*s36**4 -  \
    1024*q7a**2*s11*s12*s22*s23*s25**2*s36**4 - 512*q7a*q8a*s11*s12*s14*s24*s25**2*s36**4 - 512*q7a**2*s11*s14*s22*s24*s25**2*s36**4 + 512*q7a**2*s11*s12*s24**2*s25**2*s36**4;
    v2=v2_0+v2_1+v2_2+v2_3+v2_4+v2_5;

    v3_0=-2048*q9a*s11**4*s22**4*s33**3 + 4096*q9a*s11**3*s16*s22**3*s26*s33**3 - 4096*q9a*s11**3*s12*s22**2*s26**2*s33**3 - 2048*q9a*s11**2*s16**2*s22**2*s26**2*s33**3 + 4096*q9a*s11**2*s12*s16*s22*s26**3*s33**3 -  \
    2048*q9a*s11**2*s12**2*s26**4*s33**3 + 3072*q9a*s11**4*s22**3*s24*s33**2*s34 - 1536*q9a*s11**3*s16*s22**3*s25*s33**2*s34 - 1536*q9a*s11**3*s15*s22**3*s26*s33**2*s34 - 4608*q9a*s11**3*s16*s22**2*s24*s26*s33**2*s34 +  \
    3072*q9a*s11**3*s12*s22**2*s25*s26*s33**2*s34 + 1536*q9a*s11**2*s16**2*s22**2*s25*s26*s33**2*s34 + 1536*q9a*s11**3*s14*s22**2*s26**2*s33**2*s34 + 1536*q9a*s11**2*s15*s16*s22**2*s26**2*s33**2*s34 +  \
    3072*q9a*s11**3*s12*s22*s24*s26**2*s33**2*s34 + 1536*q9a*s11**2*s16**2*s22*s24*s26**2*s33**2*s34 - 4608*q9a*s11**2*s12*s16*s22*s25*s26**2*s33**2*s34 - 1536*q9a*s11**2*s12*s15*s22*s26**3*s33**2*s34 -  \
    1536*q9a*s11**2*s14*s16*s22*s26**3*s33**2*s34 - 1536*q9a*s11**2*s12*s16*s24*s26**3*s33**2*s34 + 3072*q9a*s11**2*s12**2*s25*s26**3*s33**2*s34 + 1536*q9a*s11**2*s12*s14*s26**4*s33**2*s34 - 2048*q9a*s11**4*s22**3*s23*s33*s34**2 -  \
    1024*q9a*s11**4*s22**2*s24**2*s33*s34**2 + 1024*q9a*s11**3*s15*s22**3*s25*s33*s34**2 + 1024*q9a*s11**3*s16*s22**2*s24*s25*s33*s34**2 - 1024*q9a*s11**3*s12*s22**2*s25**2*s33*s34**2 + 3072*q9a*s11**3*s16*s22**2*s23*s26*s33*s34**2 +  \
    1024*q9a*s11**3*s15*s22**2*s24*s26*s33*s34**2 + 1024*q9a*s11**3*s16*s22*s24**2*s26*s33*s34**2 - 2048*q9a*s11**3*s14*s22**2*s25*s26*s33*s34**2 - 1024*q9a*s11**2*s15*s16*s22**2*s25*s26*s33*s34**2 -  \
    1024*q9a*s11**2*s16**2*s22*s24*s25*s26*s33*s34**2 + 1024*q9a*s11**2*s12*s16*s22*s25**2*s26*s33*s34**2 - 1024*q9a*s11**3*s13*s22**2*s26**2*s33*s34**2 - 2048*q9a*s11**3*s12*s22*s23*s26**2*s33*s34**2 -  \
    1024*q9a*s11**2*s16**2*s22*s23*s26**2*s33*s34**2 - 1024*q9a*s11**2*s15*s16*s22*s24*s26**2*s33*s34**2 - 1024*q9a*s11**3*s12*s24**2*s26**2*s33*s34**2 + 1024*q9a*s11**2*s12*s15*s22*s25*s26**2*s33*s34**2 +  \
    2048*q9a*s11**2*s14*s16*s22*s25*s26**2*s33*s34**2 + 1024*q9a*s11**2*s12*s16*s24*s25*s26**2*s33*s34**2 - 1024*q9a*s11**2*s12**2*s25**2*s26**2*s33*s34**2 + 1024*q9a*s11**2*s13*s16*s22*s26**3*s33*s34**2 +  \
    1024*q9a*s11**2*s12*s16*s23*s26**3*s33*s34**2 + 1024*q9a*s11**2*s12*s15*s24*s26**3*s33*s34**2 - 2048*q9a*s11**2*s12*s14*s25*s26**3*s33*s34**2 - 1024*q9a*s11**2*s12*s13*s26**4*s33*s34**2 - 1024*q8a*s11**4*s22**3*s33**2*s34**2 +  \
    1536*q8a*s11**3*s16*s22**2*s26*s33**2*s34**2 - 1024*q8a*s11**3*s12*s22*s26**2*s33**2*s34**2 - 512*q8a*s11**2*s16**2*s22*s26**2*s33**2*s34**2 - 512*q7a*s11**3*s22**2*s26**2*s33**2*s34**2 + 512*q8a*s11**2*s12*s16*s26**3*s33**2*s34**2 +  \
    512*q7a*s11**2*s16*s22*s26**3*s33**2*s34**2 - 512*q7a*s11**2*s12*s26**4*s33**2*s34**2 + 1024*q9a*s11**4*s22**2*s23*s24*s34**3 - 512*q9a*s11**3*s16*s22**2*s23*s25*s34**3 - 512*q9a*s11**3*s15*s22**2*s24*s25*s34**3 +  \
    512*q9a*s11**3*s14*s22**2*s25**2*s34**3 - 512*q9a*s11**3*s15*s22**2*s23*s26*s34**3 - 1024*q9a*s11**3*s16*s22*s23*s24*s26*s34**3 + 1024*q9a*s11**3*s13*s22**2*s25*s26*s34**3 + 512*q9a*s11**2*s16**2*s22*s23*s25*s26*s34**3 +  \
    512*q9a*s11**2*s15*s16*s22*s24*s25*s26*s34**3 - 512*q9a*s11**2*s14*s16*s22*s25**2*s26*s34**3 + 512*q9a*s11**2*s15*s16*s22*s23*s26**2*s34**3 + 1024*q9a*s11**3*s12*s23*s24*s26**2*s34**3 -  \
    1024*q9a*s11**2*s13*s16*s22*s25*s26**2*s34**3 - 512*q9a*s11**2*s12*s16*s23*s25*s26**2*s34**3 - 512*q9a*s11**2*s12*s15*s24*s25*s26**2*s34**3 + 512*q9a*s11**2*s12*s14*s25**2*s26**2*s34**3 - 512*q9a*s11**2*s12*s15*s23*s26**3*s34**3 +  \
    1024*q9a*s11**2*s12*s13*s25*s26**3*s34**3 + 1024*q8a*s11**4*s22**2*s24*s33*s34**3 - 512*q8a*s11**3*s16*s22**2*s25*s33*s34**3 - 512*q8a*s11**3*s15*s22**2*s26*s33*s34**3 - 1024*q8a*s11**3*s16*s22*s24*s26*s33*s34**3 +  \
    512*q8a*s11**2*s16**2*s22*s25*s26*s33*s34**3 + 1024*q7a*s11**3*s22**2*s25*s26*s33*s34**3 + 512*q8a*s11**2*s15*s16*s22*s26**2*s33*s34**3 + 1024*q8a*s11**3*s12*s24*s26**2*s33*s34**3 - 512*q8a*s11**2*s12*s16*s25*s26**2*s33*s34**3 -  \
    1024*q7a*s11**2*s16*s22*s25*s26**2*s33*s34**3 - 512*q8a*s11**2*s12*s15*s26**3*s33*s34**3 + 1024*q7a*s11**2*s12*s25*s26**3*s33*s34**3 - 1024*q8a*s11**4*s22**2*s23*s34**4 + 512*q8a*s11**3*s15*s22**2*s25*s34**4 -  \
    512*q7a*s11**3*s22**2*s25**2*s34**4 + 1024*q8a*s11**3*s16*s22*s23*s26*s34**4 - 512*q8a*s11**2*s15*s16*s22*s25*s26*s34**4 + 512*q7a*s11**2*s16*s22*s25**2*s26*s34**4 - 1024*q8a*s11**3*s12*s23*s26**2*s34**4 +  \
    512*q8a*s11**2*s12*s15*s25*s26**2*s34**4 - 512*q7a*s11**2*s12*s25**2*s26**2*s34**4 + 3072*q9a*s11**3*s15*s22**4*s33**2*s35 - 1536*q9a*s11**3*s16*s22**3*s24*s33**2*s35 - 3072*q9a*s11**3*s12*s22**3*s25*s33**2*s35 +  \
    1536*q9a*s11**2*s16**2*s22**3*s25*s33**2*s35 - 1536*q9a*s11**3*s14*s22**3*s26*s33**2*s35 - 4608*q9a*s11**2*s15*s16*s22**3*s26*s33**2*s35 + 3072*q9a*s11**3*s12*s22**2*s24*s26*s33**2*s35 +  \
    1536*q9a*s11**2*s16**2*s22**2*s24*s26*s33**2*s35 + 1536*q9a*s11**2*s12*s16*s22**2*s25*s26*s33**2*s35 - 1536*q9a*s11*s16**3*s22**2*s25*s26*s33**2*s35 + 4608*q9a*s11**2*s12*s15*s22**2*s26**2*s33**2*s35 +  \
    1536*q9a*s11**2*s14*s16*s22**2*s26**2*s33**2*s35 + 1536*q9a*s11*s15*s16**2*s22**2*s26**2*s33**2*s35 - 4608*q9a*s11**2*s12*s16*s22*s24*s26**2*s33**2*s35 - 3072*q9a*s11**2*s12**2*s22*s25*s26**2*s33**2*s35 +  \
    3072*q9a*s11*s12*s16**2*s22*s25*s26**2*s33**2*s35 - 1536*q9a*s11**2*s12*s14*s22*s26**3*s33**2*s35 - 3072*q9a*s11*s12*s15*s16*s22*s26**3*s33**2*s35 + 3072*q9a*s11**2*s12**2*s24*s26**3*s33**2*s35 -  \
    1536*q9a*s11*s12**2*s16*s25*s26**3*s33**2*s35 + 1536*q9a*s11*s12**2*s15*s26**4*s33**2*s35 + 2048*q9a*s11**3*s16*s22**3*s23*s33*s34*s35 - 3072*q9a*s11**3*s15*s22**3*s24*s33*s34*s35 + 1024*q9a*s11**3*s16*s22**2*s24**2*s33*s34*s35 +  \
    3072*q9a*s11**3*s14*s22**3*s25*s33*s34*s35 - 1024*q9a*s11**2*s15*s16*s22**3*s25*s33*s34*s35 - 1024*q9a*s11**2*s16**2*s22**2*s24*s25*s33*s34*s35 + 1024*q9a*s11**2*s12*s16*s22**2*s25**2*s33*s34*s35 +  \
    2048*q9a*s11**3*s13*s22**3*s26*s33*s34*s35 + 1024*q9a*s11**2*s15**2*s22**3*s26*s33*s34*s35 - 4096*q9a*s11**3*s12*s22**2*s23*s26*s33*s34*s35 - 2048*q9a*s11**2*s16**2*s22**2*s23*s26*s33*s34*s35 -  \
    1024*q9a*s11**3*s14*s22**2*s24*s26*s33*s34*s35 + 4096*q9a*s11**2*s15*s16*s22**2*s24*s26*s33*s34*s35 - 1024*q9a*s11**2*s16**2*s22*s24**2*s26*s33*s34*s35 - 1024*q9a*s11**2*s12*s15*s22**2*s25*s26*s33*s34*s35 -  \
    2048*q9a*s11**2*s14*s16*s22**2*s25*s26*s33*s34*s35 + 1024*q9a*s11*s15*s16**2*s22**2*s25*s26*s33*s34*s35 + 1024*q9a*s11*s16**3*s22*s24*s25*s26*s33*s34*s35 - 1024*q9a*s11*s12*s16**2*s22*s25**2*s26*s33*s34*s35 -  \
    1024*q9a*s11**2*s14*s15*s22**2*s26**2*s33*s34*s35 - 2048*q9a*s11**2*s13*s16*s22**2*s26**2*s33*s34*s35 - 1024*q9a*s11*s15**2*s16*s22**2*s26**2*s33*s34*s35 + 6144*q9a*s11**2*s12*s16*s22*s23*s26**2*s33*s34*s35 -  \
    3072*q9a*s11**2*s12*s15*s22*s24*s26**2*s33*s34*s35 + 1024*q9a*s11**2*s14*s16*s22*s24*s26**2*s33*s34*s35 - 1024*q9a*s11*s15*s16**2*s22*s24*s26**2*s33*s34*s35 + 1024*q9a*s11**2*s12*s16*s24**2*s26**2*s33*s34*s35 +  \
    3072*q9a*s11**2*s12*s14*s22*s25*s26**2*s33*s34*s35 - 1024*q9a*s11*s14*s16**2*s22*s25*s26**2*s33*s34*s35 - 1024*q9a*s11*s12*s16**2*s24*s25*s26**2*s33*s34*s35 + 1024*q9a*s11*s12**2*s16*s25**2*s26**2*s33*s34*s35 +  \
    2048*q9a*s11**2*s12*s13*s22*s26**3*s33*s34*s35 + 1024*q9a*s11*s12*s15**2*s22*s26**3*s33*s34*s35 + 1024*q9a*s11*s14*s15*s16*s22*s26**3*s33*s34*s35 - 4096*q9a*s11**2*s12**2*s23*s26**3*s33*s34*s35 -  \
    1024*q9a*s11**2*s12*s14*s24*s26**3*s33*s34*s35 + 1024*q9a*s11*s12*s15*s16*s24*s26**3*s33*s34*s35 - 1024*q9a*s11*s12**2*s15*s25*s26**3*s33*s34*s35 + 1024*q9a*s11*s12*s14*s16*s25*s26**3*s33*s34*s35 -  \
    1024*q9a*s11*s12*s14*s15*s26**4*s33*s34*s35 + 1024*q8a*s11**3*s16*s22**3*s33**2*s34*s35 - 2048*q8a*s11**3*s12*s22**2*s26*s33**2*s34*s35 - 1024*q8a*s11**2*s16**2*s22**2*s26*s33**2*s34*s35 +  \
    1024*q7a*s11**3*s22**3*s26*s33**2*s34*s35 + 3072*q8a*s11**2*s12*s16*s22*s26**2*s33**2*s34*s35 - 1024*q7a*s11**2*s16*s22**2*s26**2*s33**2*s34*s35 - 2048*q8a*s11**2*s12**2*s26**3*s33**2*s34*s35 +  \
    1024*q7a*s11**2*s12*s22*s26**3*s33**2*s34*s35 + 1024*q9a*s11**3*s15*s22**3*s23*s34**2*s35 - 1536*q9a*s11**3*s16*s22**2*s23*s24*s34**2*s35 + 512*q9a*s11**3*s15*s22**2*s24**2*s34**2*s35 - 2048*q9a*s11**3*s13*s22**3*s25*s34**2*s35 +  \
    1024*q9a*s11**3*s12*s22**2*s23*s25*s34**2*s35 + 512*q9a*s11**2*s16**2*s22**2*s23*s25*s34**2*s35 - 512*q9a*s11**3*s14*s22**2*s24*s25*s34**2*s35 + 512*q9a*s11**2*s15*s16*s22**2*s24*s25*s34**2*s35 -  \
    512*q9a*s11**2*s14*s16*s22**2*s25**2*s34**2*s35 + 1536*q9a*s11**3*s14*s22**2*s23*s26*s34**2*s35 - 1536*q9a*s11**2*s15*s16*s22**2*s23*s26*s34**2*s35 - 512*q9a*s11**2*s15**2*s22**2*s24*s26*s34**2*s35 +  \
    1536*q9a*s11**2*s16**2*s22*s23*s24*s26*s34**2*s35 - 512*q9a*s11**2*s15*s16*s22*s24**2*s26*s34**2*s35 + 512*q9a*s11**2*s14*s15*s22**2*s25*s26*s34**2*s35 + 1536*q9a*s11**2*s13*s16*s22**2*s25*s26*s34**2*s35 -  \
    1024*q9a*s11**2*s12*s16*s22*s23*s25*s26*s34**2*s35 - 512*q9a*s11*s16**3*s22*s23*s25*s26*s34**2*s35 + 512*q9a*s11**2*s14*s16*s22*s24*s25*s26*s34**2*s35 - 512*q9a*s11*s15*s16**2*s22*s24*s25*s26*s34**2*s35 +  \
    512*q9a*s11*s14*s16**2*s22*s25**2*s26*s34**2*s35 + 512*q9a*s11**2*s13*s15*s22**2*s26**2*s34**2*s35 + 1024*q9a*s11**2*s12*s15*s22*s23*s26**2*s34**2*s35 - 1536*q9a*s11**2*s14*s16*s22*s23*s26**2*s34**2*s35 +  \
    512*q9a*s11*s15*s16**2*s22*s23*s26**2*s34**2*s35 + 512*q9a*s11*s15**2*s16*s22*s24*s26**2*s34**2*s35 - 1536*q9a*s11**2*s12*s16*s23*s24*s26**2*s34**2*s35 + 512*q9a*s11**2*s12*s15*s24**2*s26**2*s34**2*s35 -  \
    2048*q9a*s11**2*s12*s13*s22*s25*s26**2*s34**2*s35 - 512*q9a*s11*s14*s15*s16*s22*s25*s26**2*s34**2*s35 + 512*q9a*s11*s13*s16**2*s22*s25*s26**2*s34**2*s35 + 1024*q9a*s11**2*s12**2*s23*s25*s26**2*s34**2*s35 +  \
    512*q9a*s11*s12*s16**2*s23*s25*s26**2*s34**2*s35 - 512*q9a*s11**2*s12*s14*s24*s25*s26**2*s34**2*s35 + 512*q9a*s11*s12*s15*s16*s24*s25*s26**2*s34**2*s35 - 512*q9a*s11*s12*s14*s16*s25**2*s26**2*s34**2*s35 -  \
    512*q9a*s11*s13*s15*s16*s22*s26**3*s34**2*s35 + 1536*q9a*s11**2*s12*s14*s23*s26**3*s34**2*s35 - 512*q9a*s11*s12*s15*s16*s23*s26**3*s34**2*s35 - 512*q9a*s11*s12*s15**2*s24*s26**3*s34**2*s35 +  \
    512*q9a*s11*s12*s14*s15*s25*s26**3*s34**2*s35 - 512*q9a*s11*s12*s13*s16*s25*s26**3*s34**2*s35 + 512*q9a*s11*s12*s13*s15*s26**4*s34**2*s35 + 1024*q8a*s11**3*s15*s22**3*s33*s34**2*s35 -  \
    1536*q8a*s11**3*s16*s22**2*s24*s33*s34**2*s35 + 1024*q8a*s11**3*s12*s22**2*s25*s33*s34**2*s35 + 512*q8a*s11**2*s16**2*s22**2*s25*s33*s34**2*s35 - 2048*q7a*s11**3*s22**3*s25*s33*s34**2*s35 +  \
    1536*q8a*s11**3*s14*s22**2*s26*s33*s34**2*s35 - 1536*q8a*s11**2*s15*s16*s22**2*s26*s33*s34**2*s35 + 1536*q8a*s11**2*s16**2*s22*s24*s26*s33*s34**2*s35 - 1024*q8a*s11**2*s12*s16*s22*s25*s26*s33*s34**2*s35 -  \
    512*q8a*s11*s16**3*s22*s25*s26*s33*s34**2*s35 + 1536*q7a*s11**2*s16*s22**2*s25*s26*s33*s34**2*s35 + 1024*q8a*s11**2*s12*s15*s22*s26**2*s33*s34**2*s35 - 1536*q8a*s11**2*s14*s16*s22*s26**2*s33*s34**2*s35 +  \
    512*q8a*s11*s15*s16**2*s22*s26**2*s33*s34**2*s35 + 512*q7a*s11**2*s15*s22**2*s26**2*s33*s34**2*s35 - 1536*q8a*s11**2*s12*s16*s24*s26**2*s33*s34**2*s35 + 1024*q8a*s11**2*s12**2*s25*s26**2*s33*s34**2*s35 +  \
    512*q8a*s11*s12*s16**2*s25*s26**2*s33*s34**2*s35 - 2048*q7a*s11**2*s12*s22*s25*s26**2*s33*s34**2*s35 + 512*q7a*s11*s16**2*s22*s25*s26**2*s33*s34**2*s35 + 1536*q8a*s11**2*s12*s14*s26**3*s33*s34**2*s35 -  \
    512*q8a*s11*s12*s15*s16*s26**3*s33*s34**2*s35 - 512*q7a*s11*s15*s16*s22*s26**3*s33*s34**2*s35 - 512*q7a*s11*s12*s16*s25*s26**3*s33*s34**2*s35 + 512*q7a*s11*s12*s15*s26**4*s33*s34**2*s35 +  \
    2048*q8a*s11**3*s16*s22**2*s23*s34**3*s35 - 512*q8a*s11**3*s15*s22**2*s24*s34**3*s35 - 512*q8a*s11**3*s14*s22**2*s25*s34**3*s35 - 512*q8a*s11**2*s15*s16*s22**2*s25*s34**3*s35 + 1024*q7a*s11**3*s22**2*s24*s25*s34**3*s35 +  \
    512*q7a*s11**2*s16*s22**2*s25**2*s34**3*s35 - 1024*q8a*s11**3*s13*s22**2*s26*s34**3*s35 + 512*q8a*s11**2*s15**2*s22**2*s26*s34**3*s35 - 2048*q8a*s11**2*s16**2*s22*s23*s26*s34**3*s35 - 1024*q7a*s11**3*s22**2*s23*s26*s34**3*s35;
    v3_1= \
    512*q8a*s11**2*s15*s16*s22*s24*s26*s34**3*s35 + 512*q8a*s11**2*s14*s16*s22*s25*s26*s34**3*s35 + 512*q8a*s11*s15*s16**2*s22*s25*s26*s34**3*s35 - 512*q7a*s11**2*s15*s22**2*s25*s26*s34**3*s35 -  \
    1024*q7a*s11**2*s16*s22*s24*s25*s26*s34**3*s35 - 512*q7a*s11*s16**2*s22*s25**2*s26*s34**3*s35 + 1024*q8a*s11**2*s13*s16*s22*s26**2*s34**3*s35 - 512*q8a*s11*s15**2*s16*s22*s26**2*s34**3*s35 +  \
    2048*q8a*s11**2*s12*s16*s23*s26**2*s34**3*s35 + 1024*q7a*s11**2*s16*s22*s23*s26**2*s34**3*s35 - 512*q8a*s11**2*s12*s15*s24*s26**2*s34**3*s35 - 512*q8a*s11**2*s12*s14*s25*s26**2*s34**3*s35 -  \
    512*q8a*s11*s12*s15*s16*s25*s26**2*s34**3*s35 + 512*q7a*s11*s15*s16*s22*s25*s26**2*s34**3*s35 + 1024*q7a*s11**2*s12*s24*s25*s26**2*s34**3*s35 + 512*q7a*s11*s12*s16*s25**2*s26**2*s34**3*s35 -  \
    1024*q8a*s11**2*s12*s13*s26**3*s34**3*s35 + 512*q8a*s11*s12*s15**2*s26**3*s34**3*s35 - 1024*q7a*s11**2*s12*s23*s26**3*s34**3*s35 - 512*q7a*s11*s12*s15*s25*s26**3*s34**3*s35 - 2048*q9a*s11**3*s13*s22**4*s33*s35**2 -  \
    1024*q9a*s11**2*s15**2*s22**4*s33*s35**2 + 2048*q9a*s11**3*s12*s22**3*s23*s33*s35**2 - 1024*q9a*s11**2*s16**2*s22**3*s23*s33*s35**2 + 1024*q9a*s11**3*s14*s22**3*s24*s33*s35**2 + 1024*q9a*s11**2*s15*s16*s22**3*s24*s33*s35**2 -  \
    1024*q9a*s11**3*s12*s22**2*s24**2*s33*s35**2 + 2048*q9a*s11**2*s12*s15*s22**3*s25*s33*s35**2 - 2048*q9a*s11**2*s14*s16*s22**3*s25*s33*s35**2 + 1024*q9a*s11**2*s12*s16*s22**2*s24*s25*s33*s35**2 -  \
    1024*q9a*s11**2*s12**2*s22**2*s25**2*s33*s35**2 + 1024*q9a*s11**2*s14*s15*s22**3*s26*s33*s35**2 + 3072*q9a*s11**2*s13*s16*s22**3*s26*s33*s35**2 + 1024*q9a*s11*s15**2*s16*s22**3*s26*s33*s35**2 -  \
    1024*q9a*s11**2*s12*s16*s22**2*s23*s26*s33*s35**2 + 1024*q9a*s11*s16**3*s22**2*s23*s26*s33*s35**2 - 2048*q9a*s11**2*s12*s15*s22**2*s24*s26*s33*s35**2 - 1024*q9a*s11**2*s14*s16*s22**2*s24*s26*s33*s35**2 -  \
    1024*q9a*s11*s15*s16**2*s22**2*s24*s26*s33*s35**2 + 1024*q9a*s11**2*s12*s16*s22*s24**2*s26*s33*s35**2 + 1024*q9a*s11**2*s12*s14*s22**2*s25*s26*s33*s35**2 - 2048*q9a*s11*s12*s15*s16*s22**2*s25*s26*s33*s35**2 +  \
    2048*q9a*s11*s14*s16**2*s22**2*s25*s26*s33*s35**2 - 1024*q9a*s11*s12*s16**2*s22*s24*s25*s26*s33*s35**2 + 1024*q9a*s11*s12**2*s16*s22*s25**2*s26*s33*s35**2 - 3072*q9a*s11**2*s12*s13*s22**2*s26**2*s33*s35**2 -  \
    1024*q9a*s11*s12*s15**2*s22**2*s26**2*s33*s35**2 - 1024*q9a*s11*s14*s15*s16*s22**2*s26**2*s33*s35**2 - 1024*q9a*s11*s13*s16**2*s22**2*s26**2*s33*s35**2 + 2048*q9a*s11**2*s12**2*s22*s23*s26**2*s33*s35**2 -  \
    2048*q9a*s11*s12*s16**2*s22*s23*s26**2*s33*s35**2 + 1024*q9a*s11**2*s12*s14*s22*s24*s26**2*s33*s35**2 + 3072*q9a*s11*s12*s15*s16*s22*s24*s26**2*s33*s35**2 - 1024*q9a*s11**2*s12**2*s24**2*s26**2*s33*s35**2 +  \
    2048*q9a*s11*s12**2*s15*s22*s25*s26**2*s33*s35**2 - 3072*q9a*s11*s12*s14*s16*s22*s25*s26**2*s33*s35**2 + 1024*q9a*s11*s12**2*s16*s24*s25*s26**2*s33*s35**2 - 1024*q9a*s11*s12**3*s25**2*s26**2*s33*s35**2 +  \
    1024*q9a*s11*s12*s14*s15*s22*s26**3*s33*s35**2 + 2048*q9a*s11*s12*s13*s16*s22*s26**3*s33*s35**2 + 1024*q9a*s11*s12**2*s16*s23*s26**3*s33*s35**2 - 2048*q9a*s11*s12**2*s15*s24*s26**3*s33*s35**2 +  \
    1024*q9a*s11*s12**2*s14*s25*s26**3*s33*s35**2 - 1024*q9a*s11*s12**2*s13*s26**4*s33*s35**2 + 1024*q8a*s11**3*s12*s22**3*s33**2*s35**2 - 512*q8a*s11**2*s16**2*s22**3*s33**2*s35**2 - 1024*q7a*s11**3*s22**4*s33**2*s35**2 -  \
    512*q8a*s11**2*s12*s16*s22**2*s26*s33**2*s35**2 + 512*q8a*s11*s16**3*s22**2*s26*s33**2*s35**2 + 1536*q7a*s11**2*s16*s22**3*s26*s33**2*s35**2 + 1024*q8a*s11**2*s12**2*s22*s26**2*s33**2*s35**2 -  \
    1024*q8a*s11*s12*s16**2*s22*s26**2*s33**2*s35**2 - 1536*q7a*s11**2*s12*s22**2*s26**2*s33**2*s35**2 - 512*q7a*s11*s16**2*s22**2*s26**2*s33**2*s35**2 + 512*q8a*s11*s12**2*s16*s26**3*s33**2*s35**2 +  \
    1024*q7a*s11*s12*s16*s22*s26**3*s33**2*s35**2 - 512*q7a*s11*s12**2*s26**4*s33**2*s35**2 - 2048*q9a*s11**3*s14*s22**3*s23*s34*s35**2 + 1024*q9a*s11**3*s13*s22**3*s24*s34*s35**2 + 512*q9a*s11**2*s15**2*s22**3*s24*s34*s35**2 +  \
    1024*q9a*s11**3*s12*s22**2*s23*s24*s34*s35**2 + 512*q9a*s11**2*s16**2*s22**2*s23*s24*s34*s35**2 - 512*q9a*s11**2*s15*s16*s22**2*s24**2*s34*s35**2 - 512*q9a*s11**2*s14*s15*s22**3*s25*s34*s35**2 +  \
    1536*q9a*s11**2*s13*s16*s22**3*s25*s34*s35**2 - 1536*q9a*s11**2*s12*s16*s22**2*s23*s25*s34*s35**2 - 512*q9a*s11**2*s12*s15*s22**2*s24*s25*s34*s35**2 + 512*q9a*s11**2*s14*s16*s22**2*s24*s25*s34*s35**2 +  \
    512*q9a*s11**2*s12*s14*s22**2*s25**2*s34*s35**2 - 1536*q9a*s11**2*s13*s15*s22**3*s26*s34*s35**2 + 1536*q9a*s11**2*s12*s15*s22**2*s23*s26*s34*s35**2 + 1536*q9a*s11**2*s14*s16*s22**2*s23*s26*s34*s35**2 +  \
    512*q9a*s11**2*s14*s15*s22**2*s24*s26*s34*s35**2 - 1536*q9a*s11**2*s13*s16*s22**2*s24*s26*s34*s35**2 - 512*q9a*s11*s15**2*s16*s22**2*s24*s26*s34*s35**2 - 1024*q9a*s11**2*s12*s16*s22*s23*s24*s26*s34*s35**2 -  \
    512*q9a*s11*s16**3*s22*s23*s24*s26*s34*s35**2 + 512*q9a*s11*s15*s16**2*s22*s24**2*s26*s34*s35**2 - 512*q9a*s11**2*s14**2*s22**2*s25*s26*s34*s35**2 + 512*q9a*s11*s14*s15*s16*s22**2*s25*s26*s34*s35**2 -  \
    1536*q9a*s11*s13*s16**2*s22**2*s25*s26*s34*s35**2 + 1536*q9a*s11*s12*s16**2*s22*s23*s25*s26*s34*s35**2 + 512*q9a*s11*s12*s15*s16*s22*s24*s25*s26*s34*s35**2 - 512*q9a*s11*s14*s16**2*s22*s24*s25*s26*s34*s35**2 -  \
    512*q9a*s11*s12*s14*s16*s22*s25**2*s26*s34*s35**2 + 512*q9a*s11**2*s13*s14*s22**2*s26**2*s34*s35**2 + 1536*q9a*s11*s13*s15*s16*s22**2*s26**2*s34*s35**2 - 2048*q9a*s11**2*s12*s14*s22*s23*s26**2*s34*s35**2 -  \
    1536*q9a*s11*s12*s15*s16*s22*s23*s26**2*s34*s35**2 + 512*q9a*s11*s14*s16**2*s22*s23*s26**2*s34*s35**2 + 1024*q9a*s11**2*s12*s13*s22*s24*s26**2*s34*s35**2 + 512*q9a*s11*s12*s15**2*s22*s24*s26**2*s34*s35**2 -  \
    512*q9a*s11*s14*s15*s16*s22*s24*s26**2*s34*s35**2 + 512*q9a*s11*s13*s16**2*s22*s24*s26**2*s34*s35**2 + 1024*q9a*s11**2*s12**2*s23*s24*s26**2*s34*s35**2 + 512*q9a*s11*s12*s16**2*s23*s24*s26**2*s34*s35**2 -  \
    512*q9a*s11*s12*s15*s16*s24**2*s26**2*s34*s35**2 - 512*q9a*s11*s12*s14*s15*s22*s25*s26**2*s34*s35**2 + 1536*q9a*s11*s12*s13*s16*s22*s25*s26**2*s34*s35**2 + 512*q9a*s11*s14**2*s16*s22*s25*s26**2*s34*s35**2 -  \
    1536*q9a*s11*s12**2*s16*s23*s25*s26**2*s34*s35**2 - 512*q9a*s11*s12**2*s15*s24*s25*s26**2*s34*s35**2 + 512*q9a*s11*s12*s14*s16*s24*s25*s26**2*s34*s35**2 + 512*q9a*s11*s12**2*s14*s25**2*s26**2*s34*s35**2 -  \
    1536*q9a*s11*s12*s13*s15*s22*s26**3*s34*s35**2 - 512*q9a*s11*s13*s14*s16*s22*s26**3*s34*s35**2 + 1536*q9a*s11*s12**2*s15*s23*s26**3*s34*s35**2 - 512*q9a*s11*s12*s14*s16*s23*s26**3*s34*s35**2 +  \
    512*q9a*s11*s12*s14*s15*s24*s26**3*s34*s35**2 - 512*q9a*s11*s12*s13*s16*s24*s26**3*s34*s35**2 - 512*q9a*s11*s12*s14**2*s25*s26**3*s34*s35**2 + 512*q9a*s11*s12*s13*s14*s26**4*s34*s35**2 -  \
    2048*q8a*s11**3*s14*s22**3*s33*s34*s35**2 + 1024*q8a*s11**3*s12*s22**2*s24*s33*s34*s35**2 + 512*q8a*s11**2*s16**2*s22**2*s24*s33*s34*s35**2 + 1024*q7a*s11**3*s22**3*s24*s33*s34*s35**2 -  \
    1536*q8a*s11**2*s12*s16*s22**2*s25*s33*s34*s35**2 + 1536*q7a*s11**2*s16*s22**3*s25*s33*s34*s35**2 + 1536*q8a*s11**2*s12*s15*s22**2*s26*s33*s34*s35**2 + 1536*q8a*s11**2*s14*s16*s22**2*s26*s33*s34*s35**2 -  \
    1536*q7a*s11**2*s15*s22**3*s26*s33*s34*s35**2 - 1024*q8a*s11**2*s12*s16*s22*s24*s26*s33*s34*s35**2 - 512*q8a*s11*s16**3*s22*s24*s26*s33*s34*s35**2 - 1536*q7a*s11**2*s16*s22**2*s24*s26*s33*s34*s35**2 +  \
    1536*q8a*s11*s12*s16**2*s22*s25*s26*s33*s34*s35**2 - 1536*q7a*s11*s16**2*s22**2*s25*s26*s33*s34*s35**2 - 2048*q8a*s11**2*s12*s14*s22*s26**2*s33*s34*s35**2 - 1536*q8a*s11*s12*s15*s16*s22*s26**2*s33*s34*s35**2 +  \
    512*q8a*s11*s14*s16**2*s22*s26**2*s33*s34*s35**2 + 512*q7a*s11**2*s14*s22**2*s26**2*s33*s34*s35**2 + 1536*q7a*s11*s15*s16*s22**2*s26**2*s33*s34*s35**2 + 1024*q8a*s11**2*s12**2*s24*s26**2*s33*s34*s35**2 +  \
    512*q8a*s11*s12*s16**2*s24*s26**2*s33*s34*s35**2 + 1024*q7a*s11**2*s12*s22*s24*s26**2*s33*s34*s35**2 + 512*q7a*s11*s16**2*s22*s24*s26**2*s33*s34*s35**2 - 1536*q8a*s11*s12**2*s16*s25*s26**2*s33*s34*s35**2 +  \
    1536*q7a*s11*s12*s16*s22*s25*s26**2*s33*s34*s35**2 + 1536*q8a*s11*s12**2*s15*s26**3*s33*s34*s35**2 - 512*q8a*s11*s12*s14*s16*s26**3*s33*s34*s35**2 - 1536*q7a*s11*s12*s15*s22*s26**3*s33*s34*s35**2 -  \
    512*q7a*s11*s14*s16*s22*s26**3*s33*s34*s35**2 - 512*q7a*s11*s12*s16*s24*s26**3*s33*s34*s35**2 + 512*q7a*s11*s12*s14*s26**4*s33*s34*s35**2 + 1024*q8a*s11**3*s13*s22**3*s34**2*s35**2 - 512*q8a*s11**2*s15**2*s22**3*s34**2*s35**2 -  \
    2048*q8a*s11**3*s12*s22**2*s23*s34**2*s35**2 - 1024*q8a*s11**2*s16**2*s22**2*s23*s34**2*s35**2 + 1024*q7a*s11**3*s22**3*s23*s34**2*s35**2 + 512*q8a*s11**3*s14*s22**2*s24*s34**2*s35**2 + 512*q8a*s11**2*s15*s16*s22**2*s24*s34**2*s35**2 -  \
    512*q7a*s11**3*s22**2*s24**2*s34**2*s35**2 + 512*q8a*s11**2*s12*s15*s22**2*s25*s34**2*s35**2 + 512*q8a*s11**2*s14*s16*s22**2*s25*s34**2*s35**2 + 512*q7a*s11**2*s15*s22**3*s25*s34**2*s35**2 -  \
    1024*q7a*s11**2*s16*s22**2*s24*s25*s34**2*s35**2 - 512*q7a*s11**2*s12*s22**2*s25**2*s34**2*s35**2 - 1024*q8a*s11**2*s14*s15*s22**2*s26*s34**2*s35**2 + 512*q8a*s11*s15**2*s16*s22**2*s26*s34**2*s35**2 +  \
    2048*q8a*s11**2*s12*s16*s22*s23*s26*s34**2*s35**2 + 1024*q8a*s11*s16**3*s22*s23*s26*s34**2*s35**2 - 512*q8a*s11**2*s14*s16*s22*s24*s26*s34**2*s35**2 - 512*q8a*s11*s15*s16**2*s22*s24*s26*s34**2*s35**2 +  \
    512*q7a*s11**2*s15*s22**2*s24*s26*s34**2*s35**2 + 512*q7a*s11**2*s16*s22*s24**2*s26*s34**2*s35**2 - 512*q8a*s11*s12*s15*s16*s22*s25*s26*s34**2*s35**2 - 512*q8a*s11*s14*s16**2*s22*s25*s26*s34**2*s35**2 +  \
    512*q7a*s11**2*s14*s22**2*s25*s26*s34**2*s35**2 - 512*q7a*s11*s15*s16*s22**2*s25*s26*s34**2*s35**2 + 1024*q7a*s11*s16**2*s22*s24*s25*s26*s34**2*s35**2 + 512*q7a*s11*s12*s16*s22*s25**2*s26*s34**2*s35**2 +  \
    1024*q8a*s11**2*s12*s13*s22*s26**2*s34**2*s35**2 - 512*q8a*s11*s12*s15**2*s22*s26**2*s34**2*s35**2 + 1024*q8a*s11*s14*s15*s16*s22*s26**2*s34**2*s35**2 - 1024*q8a*s11*s13*s16**2*s22*s26**2*s34**2*s35**2 -  \
    1024*q7a*s11**2*s13*s22**2*s26**2*s34**2*s35**2 - 2048*q8a*s11**2*s12**2*s23*s26**2*s34**2*s35**2 - 1024*q8a*s11*s12*s16**2*s23*s26**2*s34**2*s35**2 + 1024*q7a*s11**2*s12*s22*s23*s26**2*s34**2*s35**2 -  \
    1024*q7a*s11*s16**2*s22*s23*s26**2*s34**2*s35**2 + 512*q8a*s11**2*s12*s14*s24*s26**2*s34**2*s35**2 + 512*q8a*s11*s12*s15*s16*s24*s26**2*s34**2*s35**2 - 512*q7a*s11*s15*s16*s22*s24*s26**2*s34**2*s35**2 -  \
    512*q7a*s11**2*s12*s24**2*s26**2*s34**2*s35**2 + 512*q8a*s11*s12**2*s15*s25*s26**2*s34**2*s35**2 + 512*q8a*s11*s12*s14*s16*s25*s26**2*s34**2*s35**2 + 512*q7a*s11*s12*s15*s22*s25*s26**2*s34**2*s35**2 -  \
    512*q7a*s11*s14*s16*s22*s25*s26**2*s34**2*s35**2 - 1024*q7a*s11*s12*s16*s24*s25*s26**2*s34**2*s35**2 - 512*q7a*s11*s12**2*s25**2*s26**2*s34**2*s35**2 - 1024*q8a*s11*s12*s14*s15*s26**3*s34**2*s35**2 +  \
    1024*q8a*s11*s12*s13*s16*s26**3*s34**2*s35**2 + 1024*q7a*s11*s13*s16*s22*s26**3*s34**2*s35**2 + 1024*q7a*s11*s12*s16*s23*s26**3*s34**2*s35**2 + 512*q7a*s11*s12*s15*s24*s26**3*s34**2*s35**2 +  \
    512*q7a*s11*s12*s14*s25*s26**3*s34**2*s35**2 - 1024*q7a*s11*s12*s13*s26**4*s34**2*s35**2 + 1024*q9a*s11**2*s13*s15*s22**4*s35**3 - 1024*q9a*s11**2*s12*s15*s22**3*s23*s35**3 + 1024*q9a*s11**2*s14*s16*s22**3*s23*s35**3 -  \
    512*q9a*s11**2*s14*s15*s22**3*s24*s35**3 - 512*q9a*s11**2*s13*s16*s22**3*s24*s35**3 - 512*q9a*s11**2*s12*s16*s22**2*s23*s24*s35**3 + 512*q9a*s11**2*s12*s15*s22**2*s24**2*s35**3 - 1024*q9a*s11**2*s12*s13*s22**3*s25*s35**3 +  \
    512*q9a*s11**2*s14**2*s22**3*s25*s35**3 + 1024*q9a*s11**2*s12**2*s22**2*s23*s25*s35**3 - 512*q9a*s11**2*s12*s14*s22**2*s24*s25*s35**3 - 512*q9a*s11**2*s13*s14*s22**3*s26*s35**3 - 1024*q9a*s11*s13*s15*s16*s22**3*s26*s35**3 -  \
    512*q9a*s11**2*s12*s14*s22**2*s23*s26*s35**3 + 1024*q9a*s11*s12*s15*s16*s22**2*s23*s26*s35**3 - 1024*q9a*s11*s14*s16**2*s22**2*s23*s26*s35**3 + 1024*q9a*s11**2*s12*s13*s22**2*s24*s26*s35**3 +  \
    512*q9a*s11*s14*s15*s16*s22**2*s24*s26*s35**3 + 512*q9a*s11*s13*s16**2*s22**2*s24*s26*s35**3 + 512*q9a*s11*s12*s16**2*s22*s23*s24*s26*s35**3 - 512*q9a*s11*s12*s15*s16*s22*s24**2*s26*s35**3 +  \
    1024*q9a*s11*s12*s13*s16*s22**2*s25*s26*s35**3 - 512*q9a*s11*s14**2*s16*s22**2*s25*s26*s35**3 - 1024*q9a*s11*s12**2*s16*s22*s23*s25*s26*s35**3 + 512*q9a*s11*s12*s14*s16*s22*s24*s25*s26*s35**3;
    v3_2= \
    1024*q9a*s11*s12*s13*s15*s22**2*s26**2*s35**3 + 512*q9a*s11*s13*s14*s16*s22**2*s26**2*s35**3 - 1024*q9a*s11*s12**2*s15*s22*s23*s26**2*s35**3 + 1536*q9a*s11*s12*s14*s16*s22*s23*s26**2*s35**3 -  \
    512*q9a*s11*s12*s14*s15*s22*s24*s26**2*s35**3 - 1536*q9a*s11*s12*s13*s16*s22*s24*s26**2*s35**3 - 512*q9a*s11*s12**2*s16*s23*s24*s26**2*s35**3 + 512*q9a*s11*s12**2*s15*s24**2*s26**2*s35**3 -  \
    1024*q9a*s11*s12**2*s13*s22*s25*s26**2*s35**3 + 512*q9a*s11*s12*s14**2*s22*s25*s26**2*s35**3 + 1024*q9a*s11*s12**3*s23*s25*s26**2*s35**3 - 512*q9a*s11*s12**2*s14*s24*s25*s26**2*s35**3 - 512*q9a*s11*s12*s13*s14*s22*s26**3*s35**3 -  \
    512*q9a*s11*s12**2*s14*s23*s26**3*s35**3 + 1024*q9a*s11*s12**2*s13*s24*s26**3*s35**3 - 1024*q8a*s11**2*s12*s15*s22**3*s33*s35**3 + 1024*q8a*s11**2*s14*s16*s22**3*s33*s35**3 + 1024*q7a*s11**2*s15*s22**4*s33*s35**3 -  \
    512*q8a*s11**2*s12*s16*s22**2*s24*s33*s35**3 - 512*q7a*s11**2*s16*s22**3*s24*s33*s35**3 + 1024*q8a*s11**2*s12**2*s22**2*s25*s33*s35**3 - 1024*q7a*s11**2*s12*s22**3*s25*s33*s35**3 - 512*q8a*s11**2*s12*s14*s22**2*s26*s33*s35**3 +  \
    1024*q8a*s11*s12*s15*s16*s22**2*s26*s33*s35**3 - 1024*q8a*s11*s14*s16**2*s22**2*s26*s33*s35**3 - 512*q7a*s11**2*s14*s22**3*s26*s33*s35**3 - 1024*q7a*s11*s15*s16*s22**3*s26*s33*s35**3 +  \
    512*q8a*s11*s12*s16**2*s22*s24*s26*s33*s35**3 + 1024*q7a*s11**2*s12*s22**2*s24*s26*s33*s35**3 + 512*q7a*s11*s16**2*s22**2*s24*s26*s33*s35**3 - 1024*q8a*s11*s12**2*s16*s22*s25*s26*s33*s35**3 +  \
    1024*q7a*s11*s12*s16*s22**2*s25*s26*s33*s35**3 - 1024*q8a*s11*s12**2*s15*s22*s26**2*s33*s35**3 + 1536*q8a*s11*s12*s14*s16*s22*s26**2*s33*s35**3 + 1024*q7a*s11*s12*s15*s22**2*s26**2*s33*s35**3 +  \
    512*q7a*s11*s14*s16*s22**2*s26**2*s33*s35**3 - 512*q8a*s11*s12**2*s16*s24*s26**2*s33*s35**3 - 1536*q7a*s11*s12*s16*s22*s24*s26**2*s33*s35**3 + 1024*q8a*s11*s12**3*s25*s26**2*s33*s35**3 -  \
    1024*q7a*s11*s12**2*s22*s25*s26**2*s33*s35**3 - 512*q8a*s11*s12**2*s14*s26**3*s33*s35**3 - 512*q7a*s11*s12*s14*s22*s26**3*s33*s35**3 + 1024*q7a*s11*s12**2*s24*s26**3*s33*s35**3 + 1024*q8a*s11**2*s14*s15*s22**3*s34*s35**3 -  \
    1024*q8a*s11**2*s13*s16*s22**3*s34*s35**3 + 2048*q8a*s11**2*s12*s16*s22**2*s23*s34*s35**3 - 1024*q7a*s11**2*s16*s22**3*s23*s34*s35**3 - 512*q8a*s11**2*s12*s15*s22**2*s24*s34*s35**3 - 512*q8a*s11**2*s14*s16*s22**2*s24*s34*s35**3 -  \
    512*q7a*s11**2*s15*s22**3*s24*s34*s35**3 + 512*q7a*s11**2*s16*s22**2*s24**2*s34*s35**3 - 512*q8a*s11**2*s12*s14*s22**2*s25*s34*s35**3 - 512*q7a*s11**2*s14*s22**3*s25*s34*s35**3 + 1024*q7a*s11**2*s12*s22**2*s24*s25*s34*s35**3 -  \
    1024*q8a*s11**2*s12*s13*s22**2*s26*s34*s35**3 + 512*q8a*s11**2*s14**2*s22**2*s26*s34*s35**3 - 1024*q8a*s11*s14*s15*s16*s22**2*s26*s34*s35**3 + 1024*q8a*s11*s13*s16**2*s22**2*s26*s34*s35**3 +  \
    2048*q7a*s11**2*s13*s22**3*s26*s34*s35**3 - 2048*q8a*s11*s12*s16**2*s22*s23*s26*s34*s35**3 - 1024*q7a*s11**2*s12*s22**2*s23*s26*s34*s35**3 + 1024*q7a*s11*s16**2*s22**2*s23*s26*s34*s35**3 +  \
    512*q8a*s11*s12*s15*s16*s22*s24*s26*s34*s35**3 + 512*q8a*s11*s14*s16**2*s22*s24*s26*s34*s35**3 - 512*q7a*s11**2*s14*s22**2*s24*s26*s34*s35**3 + 512*q7a*s11*s15*s16*s22**2*s24*s26*s34*s35**3 -  \
    512*q7a*s11*s16**2*s22*s24**2*s26*s34*s35**3 + 512*q8a*s11*s12*s14*s16*s22*s25*s26*s34*s35**3 + 512*q7a*s11*s14*s16*s22**2*s25*s26*s34*s35**3 - 1024*q7a*s11*s12*s16*s22*s24*s25*s26*s34*s35**3 +  \
    1024*q8a*s11*s12*s14*s15*s22*s26**2*s34*s35**3 - 512*q8a*s11*s14**2*s16*s22*s26**2*s34*s35**3 - 2048*q7a*s11*s13*s16*s22**2*s26**2*s34*s35**3 + 2048*q8a*s11*s12**2*s16*s23*s26**2*s34*s35**3 -  \
    512*q8a*s11*s12**2*s15*s24*s26**2*s34*s35**3 - 512*q8a*s11*s12*s14*s16*s24*s26**2*s34*s35**3 - 512*q7a*s11*s12*s15*s22*s24*s26**2*s34*s35**3 + 512*q7a*s11*s14*s16*s22*s24*s26**2*s34*s35**3 +  \
    512*q7a*s11*s12*s16*s24**2*s26**2*s34*s35**3 - 512*q8a*s11*s12**2*s14*s25*s26**2*s34*s35**3 - 512*q7a*s11*s12*s14*s22*s25*s26**2*s34*s35**3 + 1024*q7a*s11*s12**2*s24*s25*s26**2*s34*s35**3 -  \
    1024*q8a*s11*s12**2*s13*s26**3*s34*s35**3 + 512*q8a*s11*s12*s14**2*s26**3*s34*s35**3 + 2048*q7a*s11*s12*s13*s22*s26**3*s34*s35**3 - 1024*q7a*s11*s12**2*s23*s26**3*s34*s35**3 - 512*q7a*s11*s12*s14*s24*s26**3*s34*s35**3 +  \
    1024*q8a*s11**2*s12*s13*s22**3*s35**4 - 512*q8a*s11**2*s14**2*s22**3*s35**4 - 1024*q7a*s11**2*s13*s22**4*s35**4 - 1024*q8a*s11**2*s12**2*s22**2*s23*s35**4 + 1024*q7a*s11**2*s12*s22**3*s23*s35**4 +  \
    512*q8a*s11**2*s12*s14*s22**2*s24*s35**4 + 512*q7a*s11**2*s14*s22**3*s24*s35**4 - 512*q7a*s11**2*s12*s22**2*s24**2*s35**4 - 1024*q8a*s11*s12*s13*s16*s22**2*s26*s35**4 + 512*q8a*s11*s14**2*s16*s22**2*s26*s35**4 +  \
    1024*q7a*s11*s13*s16*s22**3*s26*s35**4 + 1024*q8a*s11*s12**2*s16*s22*s23*s26*s35**4 - 1024*q7a*s11*s12*s16*s22**2*s23*s26*s35**4 - 512*q8a*s11*s12*s14*s16*s22*s24*s26*s35**4 - 512*q7a*s11*s14*s16*s22**2*s24*s26*s35**4 +  \
    512*q7a*s11*s12*s16*s22*s24**2*s26*s35**4 + 1024*q8a*s11*s12**2*s13*s22*s26**2*s35**4 - 512*q8a*s11*s12*s14**2*s22*s26**2*s35**4 - 1024*q7a*s11*s12*s13*s22**2*s26**2*s35**4 - 1024*q8a*s11*s12**3*s23*s26**2*s35**4 +  \
    1024*q7a*s11*s12**2*s22*s23*s26**2*s35**4 + 512*q8a*s11*s12**2*s14*s24*s26**2*s35**4 + 512*q7a*s11*s12*s14*s22*s24*s26**2*s35**4 - 512*q7a*s11*s12**2*s24**2*s26**2*s35**4 - 3072*q9a*s11**3*s16*s22**3*s23*s33**2*s36 -  \
    1536*q9a*s11**3*s15*s22**3*s24*s33**2*s36 + 1536*q9a*s11**3*s16*s22**2*s24**2*s33**2*s36 - 4608*q9a*s11**3*s14*s22**3*s25*s33**2*s36 + 4608*q9a*s11**2*s15*s16*s22**3*s25*s33**2*s36 + 6144*q9a*s11**3*s12*s22**2*s24*s25*s33**2*s36 -  \
    3072*q9a*s11**2*s16**2*s22**2*s24*s25*s33**2*s36 - 4608*q9a*s11**2*s12*s16*s22**2*s25**2*s33**2*s36 + 1536*q9a*s11*s16**3*s22**2*s25**2*s33**2*s36 - 3072*q9a*s11**3*s13*s22**3*s26*s33**2*s36 +  \
    1536*q9a*s11**2*s15**2*s22**3*s26*s33**2*s36 + 6144*q9a*s11**3*s12*s22**2*s23*s26*s33**2*s36 + 3072*q9a*s11**2*s16**2*s22**2*s23*s26*s33**2*s36 + 4608*q9a*s11**3*s14*s22**2*s24*s26*s33**2*s36 -  \
    1536*q9a*s11**2*s15*s16*s22**2*s24*s26*s33**2*s36 - 6144*q9a*s11**3*s12*s22*s24**2*s26*s33**2*s36 - 7680*q9a*s11**2*s12*s15*s22**2*s25*s26*s33**2*s36 + 1536*q9a*s11**2*s14*s16*s22**2*s25*s26*s33**2*s36 -  \
    1536*q9a*s11*s15*s16**2*s22**2*s25*s26*s33**2*s36 + 6144*q9a*s11**2*s12*s16*s22*s24*s25*s26*s33**2*s36 + 6144*q9a*s11**2*s12**2*s22*s25**2*s26*s33**2*s36 - 3072*q9a*s11*s12*s16**2*s22*s25**2*s26*s33**2*s36 -  \
    3072*q9a*s11**2*s14*s15*s22**2*s26**2*s33**2*s36 + 3072*q9a*s11**2*s13*s16*s22**2*s26**2*s33**2*s36 - 9216*q9a*s11**2*s12*s16*s22*s23*s26**2*s33**2*s36 + 4608*q9a*s11**2*s12*s15*s22*s24*s26**2*s33**2*s36 -  \
    1536*q9a*s11**2*s14*s16*s22*s24*s26**2*s33**2*s36 + 1536*q9a*s11**2*s12*s16*s24**2*s26**2*s33**2*s36 + 1536*q9a*s11**2*s12*s14*s22*s25*s26**2*s33**2*s36 + 3072*q9a*s11*s12*s15*s16*s22*s25*s26**2*s33**2*s36 -  \
    6144*q9a*s11**2*s12**2*s24*s25*s26**2*s33**2*s36 + 1536*q9a*s11*s12**2*s16*s25**2*s26**2*s33**2*s36 - 3072*q9a*s11**2*s12*s13*s22*s26**3*s33**2*s36 + 1536*q9a*s11**2*s14**2*s22*s26**3*s33**2*s36 +  \
    6144*q9a*s11**2*s12**2*s23*s26**3*s33**2*s36 - 1536*q9a*s11**2*s12*s14*s24*s26**3*s33**2*s36 - 1536*q9a*s11*s12**2*s15*s25*s26**3*s33**2*s36 - 1024*q8a*s11**3*s16*s22**3*s33**3*s36 + 2048*q8a*s11**3*s12*s22**2*s26*s33**3*s36 +  \
    1024*q8a*s11**2*s16**2*s22**2*s26*s33**3*s36 - 1024*q7a*s11**3*s22**3*s26*s33**3*s36 - 3072*q8a*s11**2*s12*s16*s22*s26**2*s33**3*s36 + 1024*q7a*s11**2*s16*s22**2*s26**2*s33**3*s36 + 2048*q8a*s11**2*s12**2*s26**3*s33**3*s36 -  \
    1024*q7a*s11**2*s12*s22*s26**3*s33**3*s36 + 2048*q9a*s11**3*s15*s22**3*s23*s33*s34*s36 + 1024*q9a*s11**3*s16*s22**2*s23*s24*s33*s34*s36 + 1024*q9a*s11**3*s15*s22**2*s24**2*s33*s34*s36 -  \
    1024*q9a*s11**3*s16*s22*s24**3*s33*s34*s36 + 4096*q9a*s11**3*s13*s22**3*s25*s33*s34*s36 - 2048*q9a*s11**2*s15**2*s22**3*s25*s33*s34*s36 - 6144*q9a*s11**3*s12*s22**2*s23*s25*s33*s34*s36 +  \
    1024*q9a*s11**2*s16**2*s22**2*s23*s25*s33*s34*s36 + 1024*q9a*s11**3*s14*s22**2*s24*s25*s33*s34*s36 - 2048*q9a*s11**2*s15*s16*s22**2*s24*s25*s33*s34*s36 - 2048*q9a*s11**3*s12*s22*s24**2*s25*s33*s34*s36 +  \
    2048*q9a*s11**2*s16**2*s22*s24**2*s25*s33*s34*s36 + 4096*q9a*s11**2*s12*s15*s22**2*s25**2*s33*s34*s36 + 1024*q9a*s11**2*s14*s16*s22**2*s25**2*s33*s34*s36 - 1024*q9a*s11*s15*s16**2*s22**2*s25**2*s33*s34*s36 +  \
    1024*q9a*s11**2*s12*s16*s22*s24*s25**2*s33*s34*s36 - 1024*q9a*s11*s16**3*s22*s24*s25**2*s33*s34*s36 - 2048*q9a*s11**2*s12**2*s22*s25**3*s33*s34*s36 + 1024*q9a*s11*s12*s16**2*s22*s25**3*s33*s34*s36 -  \
    5120*q9a*s11**3*s14*s22**2*s23*s26*s33*s34*s36 - 1024*q9a*s11**2*s15*s16*s22**2*s23*s26*s33*s34*s36 - 1024*q9a*s11**2*s15**2*s22**2*s24*s26*s33*s34*s36 + 4096*q9a*s11**3*s12*s22*s23*s24*s26*s33*s34*s36 -  \
    2048*q9a*s11**2*s16**2*s22*s23*s24*s26*s33*s34*s36 - 1024*q9a*s11**3*s14*s22*s24**2*s26*s33*s34*s36 + 1024*q9a*s11**2*s15*s16*s22*s24**2*s26*s33*s34*s36 + 2048*q9a*s11**3*s12*s24**3*s26*s33*s34*s36 +  \
    4096*q9a*s11**2*s14*s15*s22**2*s25*s26*s33*s34*s36 - 3072*q9a*s11**2*s13*s16*s22**2*s25*s26*s33*s34*s36 + 1024*q9a*s11*s15**2*s16*s22**2*s25*s26*s33*s34*s36 + 2048*q9a*s11**2*s12*s16*s22*s23*s25*s26*s33*s34*s36 -  \
    2048*q9a*s11**2*s14*s16*s22*s24*s25*s26*s33*s34*s36 + 1024*q9a*s11*s15*s16**2*s22*s24*s25*s26*s33*s34*s36 - 3072*q9a*s11**2*s12*s16*s24**2*s25*s26*s33*s34*s36 - 5120*q9a*s11**2*s12*s14*s22*s25**2*s26*s33*s34*s36 +  \
    1024*q9a*s11*s14*s16**2*s22*s25**2*s26*s33*s34*s36 + 2048*q9a*s11**2*s12**2*s24*s25**2*s26*s33*s34*s36 + 1024*q9a*s11*s12*s16**2*s24*s25**2*s26*s33*s34*s36 - 1024*q9a*s11*s12**2*s16*s25**3*s26*s33*s34*s36 +  \
    1024*q9a*s11**2*s13*s15*s22**2*s26**2*s33*s34*s36 + 4096*q9a*s11**2*s14*s16*s22*s23*s26**2*s33*s34*s36 + 1024*q9a*s11**2*s14*s15*s22*s24*s26**2*s33*s34*s36 - 1024*q9a*s11**2*s13*s16*s22*s24*s26**2*s33*s34*s36 +  \
    1024*q9a*s11**2*s12*s16*s23*s24*s26**2*s33*s34*s36 - 2048*q9a*s11**2*s12*s15*s24**2*s26**2*s33*s34*s36 + 2048*q9a*s11**2*s12*s13*s22*s25*s26**2*s33*s34*s36 - 2048*q9a*s11**2*s14**2*s22*s25*s26**2*s33*s34*s36 -  \
    1024*q9a*s11*s12*s15**2*s22*s25*s26**2*s33*s34*s36 - 1024*q9a*s11*s14*s15*s16*s22*s25*s26**2*s33*s34*s36 - 2048*q9a*s11**2*s12**2*s23*s25*s26**2*s33*s34*s36 + 5120*q9a*s11**2*s12*s14*s24*s25*s26**2*s33*s34*s36 -  \
    1024*q9a*s11*s12*s15*s16*s24*s25*s26**2*s33*s34*s36 + 1024*q9a*s11*s12**2*s15*s25**2*s26**2*s33*s34*s36 - 1024*q9a*s11*s12*s14*s16*s25**2*s26**2*s33*s34*s36 - 1024*q9a*s11**2*s13*s14*s22*s26**3*s33*s34*s36 -  \
    3072*q9a*s11**2*s12*s14*s23*s26**3*s33*s34*s36 + 2048*q9a*s11**2*s12*s13*s24*s26**3*s33*s34*s36 + 1024*q9a*s11*s12*s14*s15*s25*s26**3*s33*s34*s36 + 1024*q8a*s11**3*s15*s22**3*s33**2*s34*s36 +  \
    512*q8a*s11**3*s16*s22**2*s24*s33**2*s34*s36 - 3072*q8a*s11**3*s12*s22**2*s25*s33**2*s34*s36 + 512*q8a*s11**2*s16**2*s22**2*s25*s33**2*s34*s36 + 2048*q7a*s11**3*s22**3*s25*s33**2*s34*s36 -  \
    2560*q8a*s11**3*s14*s22**2*s26*s33**2*s34*s36 - 512*q8a*s11**2*s15*s16*s22**2*s26*s33**2*s34*s36 + 2048*q8a*s11**3*s12*s22*s24*s26*s33**2*s34*s36 - 1024*q8a*s11**2*s16**2*s22*s24*s26*s33**2*s34*s36 +  \
    1024*q8a*s11**2*s12*s16*s22*s25*s26*s33**2*s34*s36 - 1536*q7a*s11**2*s16*s22**2*s25*s26*s33**2*s34*s36 + 2048*q8a*s11**2*s14*s16*s22*s26**2*s33**2*s34*s36 + 512*q7a*s11**2*s15*s22**2*s26**2*s33**2*s34*s36;
    v3_3= \
    512*q8a*s11**2*s12*s16*s24*s26**2*s33**2*s34*s36 - 512*q7a*s11**2*s16*s22*s24*s26**2*s33**2*s34*s36 - 1024*q8a*s11**2*s12**2*s25*s26**2*s33**2*s34*s36 + 1024*q7a*s11**2*s12*s22*s25*s26**2*s33**2*s34*s36 -  \
    1536*q8a*s11**2*s12*s14*s26**3*s33**2*s34*s36 - 512*q7a*s11**2*s14*s22*s26**3*s33**2*s34*s36 + 1024*q7a*s11**2*s12*s24*s26**3*s33**2*s34*s36 - 1024*q9a*s11**3*s16*s22**2*s23**2*s34**2*s36 -  \
    1536*q9a*s11**3*s15*s22**2*s23*s24*s34**2*s36 + 1024*q9a*s11**3*s16*s22*s23*s24**2*s34**2*s36 + 512*q9a*s11**3*s14*s22**2*s23*s25*s34**2*s36 + 1536*q9a*s11**2*s15*s16*s22**2*s23*s25*s34**2*s36 -  \
    1024*q9a*s11**3*s13*s22**2*s24*s25*s34**2*s36 + 1024*q9a*s11**2*s15**2*s22**2*s24*s25*s34**2*s36 + 2048*q9a*s11**3*s12*s22*s23*s24*s25*s34**2*s36 - 1536*q9a*s11**2*s16**2*s22*s23*s24*s25*s34**2*s36 -  \
    512*q9a*s11**2*s15*s16*s22*s24**2*s25*s34**2*s36 - 1024*q9a*s11**2*s14*s15*s22**2*s25**2*s34**2*s36 - 512*q9a*s11**2*s13*s16*s22**2*s25**2*s34**2*s36 - 1024*q9a*s11**2*s12*s16*s22*s23*s25**2*s34**2*s36 +  \
    512*q9a*s11*s16**3*s22*s23*s25**2*s34**2*s36 - 1024*q9a*s11**2*s12*s15*s22*s24*s25**2*s34**2*s36 + 512*q9a*s11**2*s14*s16*s22*s24*s25**2*s34**2*s36 + 512*q9a*s11*s15*s16**2*s22*s24*s25**2*s34**2*s36 +  \
    1024*q9a*s11**2*s12*s14*s22*s25**3*s34**2*s36 - 512*q9a*s11*s14*s16**2*s22*s25**3*s34**2*s36 + 1024*q9a*s11**3*s13*s22**2*s23*s26*s34**2*s36 + 512*q9a*s11**2*s15**2*s22**2*s23*s26*s34**2*s36 +  \
    1024*q9a*s11**2*s16**2*s22*s23**2*s26*s34**2*s36 + 1024*q9a*s11**3*s14*s22*s23*s24*s26*s34**2*s36 - 2048*q9a*s11**3*s12*s23*s24**2*s26*s34**2*s36 - 1536*q9a*s11**2*s13*s15*s22**2*s25*s26*s34**2*s36 -  \
    1024*q9a*s11**2*s12*s15*s22*s23*s25*s26*s34**2*s36 - 1024*q9a*s11**2*s14*s16*s22*s23*s25*s26*s34**2*s36 - 512*q9a*s11*s15*s16**2*s22*s23*s25*s26*s34**2*s36 - 512*q9a*s11**2*s14*s15*s22*s24*s25*s26*s34**2*s36 +  \
    2048*q9a*s11**2*s13*s16*s22*s24*s25*s26*s34**2*s36 - 512*q9a*s11*s15**2*s16*s22*s24*s25*s26*s34**2*s36 + 2048*q9a*s11**2*s12*s16*s23*s24*s25*s26*s34**2*s36 + 1024*q9a*s11**2*s12*s15*s24**2*s25*s26*s34**2*s36 +  \
    2048*q9a*s11**2*s12*s13*s22*s25**2*s26*s34**2*s36 + 512*q9a*s11**2*s14**2*s22*s25**2*s26*s34**2*s36 + 512*q9a*s11*s14*s15*s16*s22*s25**2*s26*s34**2*s36 - 512*q9a*s11*s13*s16**2*s22*s25**2*s26*s34**2*s36 -  \
    512*q9a*s11*s12*s16**2*s23*s25**2*s26*s34**2*s36 - 1024*q9a*s11**2*s12*s14*s24*s25**2*s26*s34**2*s36 - 512*q9a*s11*s12*s15*s16*s24*s25**2*s26*s34**2*s36 + 512*q9a*s11*s12*s14*s16*s25**3*s26*s34**2*s36 -  \
    512*q9a*s11**2*s14*s15*s22*s23*s26**2*s34**2*s36 - 1024*q9a*s11**2*s13*s16*s22*s23*s26**2*s34**2*s36 - 1024*q9a*s11**2*s12*s16*s23**2*s26**2*s34**2*s36 + 512*q9a*s11**2*s12*s15*s23*s24*s26**2*s34**2*s36 +  \
    1024*q9a*s11**2*s13*s14*s22*s25*s26**2*s34**2*s36 + 512*q9a*s11*s13*s15*s16*s22*s25*s26**2*s34**2*s36 + 512*q9a*s11**2*s12*s14*s23*s25*s26**2*s34**2*s36 + 512*q9a*s11*s12*s15*s16*s23*s25*s26**2*s34**2*s36 -  \
    3072*q9a*s11**2*s12*s13*s24*s25*s26**2*s34**2*s36 + 512*q9a*s11*s12*s15**2*s24*s25*s26**2*s34**2*s36 - 512*q9a*s11*s12*s14*s15*s25**2*s26**2*s34**2*s36 + 512*q9a*s11*s12*s13*s16*s25**2*s26**2*s34**2*s36 +  \
    1024*q9a*s11**2*s12*s13*s23*s26**3*s34**2*s36 - 512*q9a*s11*s12*s13*s15*s25*s26**3*s34**2*s36 - 2048*q8a*s11**3*s16*s22**2*s23*s33*s34**2*s36 - 1536*q8a*s11**3*s15*s22**2*s24*s33*s34**2*s36 +  \
    1024*q8a*s11**3*s16*s22*s24**2*s33*s34**2*s36 + 512*q8a*s11**3*s14*s22**2*s25*s33*s34**2*s36 + 1536*q8a*s11**2*s15*s16*s22**2*s25*s33*s34**2*s36 + 2048*q8a*s11**3*s12*s22*s24*s25*s33*s34**2*s36 -  \
    1536*q8a*s11**2*s16**2*s22*s24*s25*s33*s34**2*s36 - 1024*q7a*s11**3*s22**2*s24*s25*s33*s34**2*s36 - 1024*q8a*s11**2*s12*s16*s22*s25**2*s33*s34**2*s36 + 512*q8a*s11*s16**3*s22*s25**2*s33*s34**2*s36 -  \
    512*q7a*s11**2*s16*s22**2*s25**2*s33*s34**2*s36 + 1024*q8a*s11**3*s13*s22**2*s26*s33*s34**2*s36 + 512*q8a*s11**2*s15**2*s22**2*s26*s33*s34**2*s36 + 2048*q8a*s11**2*s16**2*s22*s23*s26*s33*s34**2*s36 +  \
    1024*q7a*s11**3*s22**2*s23*s26*s33*s34**2*s36 + 1024*q8a*s11**3*s14*s22*s24*s26*s33*s34**2*s36 - 2048*q8a*s11**3*s12*s24**2*s26*s33*s34**2*s36 - 1024*q8a*s11**2*s12*s15*s22*s25*s26*s33*s34**2*s36 -  \
    1024*q8a*s11**2*s14*s16*s22*s25*s26*s33*s34**2*s36 - 512*q8a*s11*s15*s16**2*s22*s25*s26*s33*s34**2*s36 - 1536*q7a*s11**2*s15*s22**2*s25*s26*s33*s34**2*s36 + 2048*q8a*s11**2*s12*s16*s24*s25*s26*s33*s34**2*s36 +  \
    2048*q7a*s11**2*s16*s22*s24*s25*s26*s33*s34**2*s36 - 512*q8a*s11*s12*s16**2*s25**2*s26*s33*s34**2*s36 + 2048*q7a*s11**2*s12*s22*s25**2*s26*s33*s34**2*s36 - 512*q7a*s11*s16**2*s22*s25**2*s26*s33*s34**2*s36 -  \
    512*q8a*s11**2*s14*s15*s22*s26**2*s33*s34**2*s36 - 1024*q8a*s11**2*s13*s16*s22*s26**2*s33*s34**2*s36 - 2048*q8a*s11**2*s12*s16*s23*s26**2*s33*s34**2*s36 - 1024*q7a*s11**2*s16*s22*s23*s26**2*s33*s34**2*s36 +  \
    512*q8a*s11**2*s12*s15*s24*s26**2*s33*s34**2*s36 + 512*q8a*s11**2*s12*s14*s25*s26**2*s33*s34**2*s36 + 512*q8a*s11*s12*s15*s16*s25*s26**2*s33*s34**2*s36 + 1024*q7a*s11**2*s14*s22*s25*s26**2*s33*s34**2*s36 +  \
    512*q7a*s11*s15*s16*s22*s25*s26**2*s33*s34**2*s36 - 3072*q7a*s11**2*s12*s24*s25*s26**2*s33*s34**2*s36 + 512*q7a*s11*s12*s16*s25**2*s26**2*s33*s34**2*s36 + 1024*q8a*s11**2*s12*s13*s26**3*s33*s34**2*s36 +  \
    1024*q7a*s11**2*s12*s23*s26**3*s33*s34**2*s36 - 512*q7a*s11*s12*s15*s25*s26**3*s33*s34**2*s36 + 2048*q8a*s11**3*s15*s22**2*s23*s34**3*s36 - 1024*q8a*s11**3*s16*s22*s23*s24*s34**3*s36 -  \
    1024*q8a*s11**2*s15**2*s22**2*s25*s34**3*s36 - 2048*q8a*s11**3*s12*s22*s23*s25*s34**3*s36 + 1024*q8a*s11**2*s16**2*s22*s23*s25*s34**3*s36 + 512*q8a*s11**2*s15*s16*s22*s24*s25*s34**3*s36 +  \
    1024*q8a*s11**2*s12*s15*s22*s25**2*s34**3*s36 - 512*q8a*s11*s15*s16**2*s22*s25**2*s34**3*s36 + 1024*q7a*s11**2*s15*s22**2*s25**2*s34**3*s36 - 512*q7a*s11**2*s16*s22*s24*s25**2*s34**3*s36 - 1024*q7a*s11**2*s12*s22*s25**3*s34**3*s36 +  \
    512*q7a*s11*s16**2*s22*s25**3*s34**3*s36 - 1024*q8a*s11**3*s14*s22*s23*s26*s34**3*s36 - 1024*q8a*s11**2*s15*s16*s22*s23*s26*s34**3*s36 + 2048*q8a*s11**3*s12*s23*s24*s26*s34**3*s36 +  \
    512*q8a*s11**2*s14*s15*s22*s25*s26*s34**3*s36 + 512*q8a*s11*s15**2*s16*s22*s25*s26*s34**3*s36 - 1024*q8a*s11**2*s12*s16*s23*s25*s26*s34**3*s36 - 1024*q8a*s11**2*s12*s15*s24*s25*s26*s34**3*s36 +  \
    512*q8a*s11*s12*s15*s16*s25**2*s26*s34**3*s36 - 512*q7a*s11**2*s14*s22*s25**2*s26*s34**3*s36 - 512*q7a*s11*s15*s16*s22*s25**2*s26*s34**3*s36 + 1024*q7a*s11**2*s12*s24*s25**2*s26*s34**3*s36 -  \
    512*q7a*s11*s12*s16*s25**3*s26*s34**3*s36 + 1024*q8a*s11**2*s12*s15*s23*s26**2*s34**3*s36 - 512*q8a*s11*s12*s15**2*s25*s26**2*s34**3*s36 + 512*q7a*s11*s12*s15*s25**2*s26**2*s34**3*s36 +  \
    4096*q9a*s11**3*s14*s22**3*s23*s33*s35*s36 + 2048*q9a*s11**3*s13*s22**3*s24*s33*s35*s36 + 1024*q9a*s11**2*s15**2*s22**3*s24*s33*s35*s36 - 6144*q9a*s11**3*s12*s22**2*s23*s24*s33*s35*s36 +  \
    1024*q9a*s11**2*s16**2*s22**2*s23*s24*s33*s35*s36 - 2048*q9a*s11**3*s14*s22**2*s24**2*s33*s35*s36 - 1024*q9a*s11**2*s15*s16*s22**2*s24**2*s33*s35*s36 + 2048*q9a*s11**3*s12*s22*s24**3*s33*s35*s36 +  \
    1024*q9a*s11**2*s14*s15*s22**3*s25*s33*s35*s36 - 5120*q9a*s11**2*s13*s16*s22**3*s25*s33*s35*s36 - 1024*q9a*s11*s15**2*s16*s22**3*s25*s33*s35*s36 + 5120*q9a*s11**2*s12*s16*s22**2*s23*s25*s33*s35*s36 -  \
    1024*q9a*s11*s16**3*s22**2*s23*s25*s33*s35*s36 - 3072*q9a*s11**2*s12*s15*s22**2*s24*s25*s33*s35*s36 + 4096*q9a*s11**2*s14*s16*s22**2*s24*s25*s33*s35*s36 + 1024*q9a*s11*s15*s16**2*s22**2*s24*s25*s33*s35*s36 -  \
    3072*q9a*s11**2*s12*s16*s22*s24**2*s25*s33*s35*s36 - 1024*q9a*s11**2*s12*s14*s22**2*s25**2*s33*s35*s36 + 2048*q9a*s11*s12*s15*s16*s22**2*s25**2*s33*s35*s36 - 2048*q9a*s11*s14*s16**2*s22**2*s25**2*s33*s35*s36 +  \
    2048*q9a*s11**2*s12**2*s22*s24*s25**2*s33*s35*s36 + 1024*q9a*s11*s12*s16**2*s22*s24*s25**2*s33*s35*s36 - 1024*q9a*s11*s12**2*s16*s22*s25**3*s33*s35*s36 + 1024*q9a*s11**2*s13*s15*s22**3*s26*s33*s35*s36 -  \
    1024*q9a*s11*s15**3*s22**3*s26*s33*s35*s36 - 1024*q9a*s11**2*s12*s15*s22**2*s23*s26*s33*s35*s36 - 3072*q9a*s11**2*s14*s16*s22**2*s23*s26*s33*s35*s36 - 1024*q9a*s11*s15*s16**2*s22**2*s23*s26*s33*s35*s36 -  \
    2048*q9a*s11**2*s14*s15*s22**2*s24*s26*s33*s35*s36 - 1024*q9a*s11**2*s13*s16*s22**2*s24*s26*s33*s35*s36 + 1024*q9a*s11*s15**2*s16*s22**2*s24*s26*s33*s35*s36 + 2048*q9a*s11**2*s12*s16*s22*s23*s24*s26*s33*s35*s36 +  \
    3072*q9a*s11**2*s12*s15*s22*s24**2*s26*s33*s35*s36 + 1024*q9a*s11**2*s14*s16*s22*s24**2*s26*s33*s35*s36 - 1024*q9a*s11**2*s12*s16*s24**3*s26*s33*s35*s36 + 4096*q9a*s11**2*s12*s13*s22**2*s25*s26*s33*s35*s36 +  \
    1024*q9a*s11**2*s14**2*s22**2*s25*s26*s33*s35*s36 + 4096*q9a*s11*s12*s15**2*s22**2*s25*s26*s33*s35*s36 - 2048*q9a*s11*s14*s15*s16*s22**2*s25*s26*s33*s35*s36 + 4096*q9a*s11*s13*s16**2*s22**2*s25*s26*s33*s35*s36 -  \
    4096*q9a*s11**2*s12**2*s22*s23*s25*s26*s33*s35*s36 - 4096*q9a*s11**2*s12*s14*s22*s24*s25*s26*s33*s35*s36 - 2048*q9a*s11*s12*s15*s16*s22*s24*s25*s26*s33*s35*s36 - 1024*q9a*s11*s14*s16**2*s22*s24*s25*s26*s33*s35*s36 +  \
    2048*q9a*s11**2*s12**2*s24**2*s25*s26*s33*s35*s36 + 1024*q9a*s11*s12*s16**2*s24**2*s25*s26*s33*s35*s36 - 5120*q9a*s11*s12**2*s15*s22*s25**2*s26*s33*s35*s36 + 6144*q9a*s11*s12*s14*s16*s22*s25**2*s26*s33*s35*s36 -  \
    3072*q9a*s11*s12**2*s16*s24*s25**2*s26*s33*s35*s36 + 2048*q9a*s11*s12**3*s25**3*s26*s33*s35*s36 + 1024*q9a*s11**2*s13*s14*s22**2*s26**2*s33*s35*s36 + 2048*q9a*s11*s14*s15**2*s22**2*s26**2*s33*s35*s36 -  \
    2048*q9a*s11*s13*s15*s16*s22**2*s26**2*s33*s35*s36 + 2048*q9a*s11**2*s12*s14*s22*s23*s26**2*s33*s35*s36 + 4096*q9a*s11*s12*s15*s16*s22*s23*s26**2*s33*s35*s36 - 1024*q9a*s11**2*s14**2*s22*s24*s26**2*s33*s35*s36 -  \
    3072*q9a*s11*s12*s15**2*s22*s24*s26**2*s33*s35*s36 + 1024*q9a*s11*s14*s15*s16*s22*s24*s26**2*s33*s35*s36 - 2048*q9a*s11**2*s12**2*s23*s24*s26**2*s33*s35*s36 + 1024*q9a*s11**2*s12*s14*s24**2*s26**2*s33*s35*s36 -  \
    1024*q9a*s11*s12*s15*s16*s24**2*s26**2*s33*s35*s36 - 2048*q9a*s11*s12*s14*s15*s22*s25*s26**2*s33*s35*s36 - 6144*q9a*s11*s12*s13*s16*s22*s25*s26**2*s33*s35*s36 + 1024*q9a*s11*s14**2*s16*s22*s25*s26**2*s33*s35*s36 +  \
    1024*q9a*s11*s12**2*s16*s23*s25*s26**2*s33*s35*s36 + 5120*q9a*s11*s12**2*s15*s24*s25*s26**2*s33*s35*s36 - 1024*q9a*s11*s12*s14*s16*s24*s25*s26**2*s33*s35*s36 - 2048*q9a*s11*s12**2*s14*s25**2*s26**2*s33*s35*s36 +  \
    2048*q9a*s11*s12*s13*s15*s22*s26**3*s33*s35*s36 - 1024*q9a*s11*s14**2*s15*s22*s26**3*s33*s35*s36 - 3072*q9a*s11*s12**2*s15*s23*s26**3*s33*s35*s36 + 1024*q9a*s11*s12*s14*s15*s24*s26**3*s33*s35*s36 +  \
    2048*q9a*s11*s12**2*s13*s25*s26**3*s33*s35*s36 + 2048*q8a*s11**3*s14*s22**3*s33**2*s35*s36 - 3072*q8a*s11**3*s12*s22**2*s24*s33**2*s35*s36 + 512*q8a*s11**2*s16**2*s22**2*s24*s33**2*s35*s36 +  \
    1024*q7a*s11**3*s22**3*s24*s33**2*s35*s36 + 2560*q8a*s11**2*s12*s16*s22**2*s25*s33**2*s35*s36 - 512*q8a*s11*s16**3*s22**2*s25*s33**2*s35*s36 - 2560*q7a*s11**2*s16*s22**3*s25*s33**2*s35*s36 -  \
    512*q8a*s11**2*s12*s15*s22**2*s26*s33**2*s35*s36 - 1536*q8a*s11**2*s14*s16*s22**2*s26*s33**2*s35*s36 - 512*q8a*s11*s15*s16**2*s22**2*s26*s33**2*s35*s36 + 512*q7a*s11**2*s15*s22**3*s26*s33**2*s35*s36 +  \
    1024*q8a*s11**2*s12*s16*s22*s24*s26*s33**2*s35*s36 - 512*q7a*s11**2*s16*s22**2*s24*s26*s33**2*s35*s36 - 2048*q8a*s11**2*s12**2*s22*s25*s26*s33**2*s35*s36 + 2048*q7a*s11**2*s12*s22**2*s25*s26*s33**2*s35*s36 +  \
    2048*q7a*s11*s16**2*s22**2*s25*s26*s33**2*s35*s36 + 1024*q8a*s11**2*s12*s14*s22*s26**2*s33**2*s35*s36 + 2048*q8a*s11*s12*s15*s16*s22*s26**2*s33**2*s35*s36 + 512*q7a*s11**2*s14*s22**2*s26**2*s33**2*s35*s36 -  \
    1024*q7a*s11*s15*s16*s22**2*s26**2*s33**2*s35*s36 - 1024*q8a*s11**2*s12**2*s24*s26**2*s33**2*s35*s36 + 512*q8a*s11*s12**2*s16*s25*s26**2*s33**2*s35*s36 - 3072*q7a*s11*s12*s16*s22*s25*s26**2*s33**2*s35*s36 -  \
    1536*q8a*s11*s12**2*s15*s26**3*s33**2*s35*s36 + 1024*q7a*s11*s12*s15*s22*s26**3*s33**2*s35*s36 + 1024*q7a*s11*s12**2*s25*s26**3*s33**2*s35*s36 - 4096*q9a*s11**3*s13*s22**3*s23*s34*s35*s36 +  \
    4096*q9a*s11**3*s12*s22**2*s23**2*s34*s35*s36 + 2048*q9a*s11**3*s14*s22**2*s23*s24*s34*s35*s36 - 512*q9a*s11**2*s15**2*s22**2*s24**2*s34*s35*s36 - 2048*q9a*s11**3*s12*s22*s23*s24**2*s34*s35*s36 -  \
    512*q9a*s11**2*s16**2*s22*s23*s24**2*s34*s35*s36 + 512*q9a*s11**2*s15*s16*s22*s24**3*s34*s35*s36 + 2048*q9a*s11**2*s13*s15*s22**3*s25*s34*s35*s36 - 2048*q9a*s11**2*s12*s15*s22**2*s23*s25*s34*s35*s36 -  \
    2048*q9a*s11**2*s14*s16*s22**2*s23*s25*s34*s35*s36 - 1024*q9a*s11**2*s14*s15*s22**2*s24*s25*s34*s35*s36 + 512*q9a*s11*s15**2*s16*s22**2*s24*s25*s34*s35*s36 + 2048*q9a*s11**2*s12*s16*s22*s23*s24*s25*s34*s35*s36;
    v3_4= \
    512*q9a*s11*s16**3*s22*s23*s24*s25*s34*s35*s36 + 2048*q9a*s11**2*s12*s15*s22*s24**2*s25*s34*s35*s36 - 512*q9a*s11**2*s14*s16*s22*s24**2*s25*s34*s35*s36 - 512*q9a*s11*s15*s16**2*s22*s24**2*s25*s34*s35*s36 -  \
    2048*q9a*s11**2*s12*s13*s22**2*s25**2*s34*s35*s36 + 1536*q9a*s11**2*s14**2*s22**2*s25**2*s34*s35*s36 - 512*q9a*s11*s14*s15*s16*s22**2*s25**2*s34*s35*s36 + 1536*q9a*s11*s13*s16**2*s22**2*s25**2*s34*s35*s36 +  \
    2048*q9a*s11**2*s12**2*s22*s23*s25**2*s34*s35*s36 - 1536*q9a*s11*s12*s16**2*s22*s23*s25**2*s34*s35*s36 - 2048*q9a*s11**2*s12*s14*s22*s24*s25**2*s34*s35*s36 - 512*q9a*s11*s12*s15*s16*s22*s24*s25**2*s34*s35*s36 +  \
    512*q9a*s11*s14*s16**2*s22*s24*s25**2*s34*s35*s36 + 512*q9a*s11*s12*s14*s16*s22*s25**3*s34*s35*s36 + 4096*q9a*s11**2*s13*s16*s22**2*s23*s26*s34*s35*s36 - 4096*q9a*s11**2*s12*s16*s22*s23**2*s26*s34*s35*s36 +  \
    512*q9a*s11*s15**3*s22**2*s24*s26*s34*s35*s36 - 1024*q9a*s11**2*s14*s16*s22*s23*s24*s26*s34*s35*s36 + 512*q9a*s11*s15*s16**2*s22*s23*s24*s26*s34*s35*s36 + 512*q9a*s11**2*s14*s15*s22*s24**2*s26*s34*s35*s36 -  \
    512*q9a*s11*s15**2*s16*s22*s24**2*s26*s34*s35*s36 + 2048*q9a*s11**2*s12*s16*s23*s24**2*s26*s34*s35*s36 - 1024*q9a*s11**2*s12*s15*s24**3*s26*s34*s35*s36 - 2048*q9a*s11**2*s13*s14*s22**2*s25*s26*s34*s35*s36 -  \
    512*q9a*s11*s14*s15**2*s22**2*s25*s26*s34*s35*s36 - 1024*q9a*s11*s13*s15*s16*s22**2*s25*s26*s34*s35*s36 + 4096*q9a*s11**2*s12*s14*s22*s23*s25*s26*s34*s35*s36 + 1024*q9a*s11*s12*s15*s16*s22*s23*s25*s26*s34*s35*s36 +  \
    512*q9a*s11*s14*s16**2*s22*s23*s25*s26*s34*s35*s36 + 2048*q9a*s11**2*s12*s13*s22*s24*s25*s26*s34*s35*s36 - 512*q9a*s11**2*s14**2*s22*s24*s25*s26*s34*s35*s36 - 1536*q9a*s11*s12*s15**2*s22*s24*s25*s26*s34*s35*s36 +  \
    2048*q9a*s11*s14*s15*s16*s22*s24*s25*s26*s34*s35*s36 - 1536*q9a*s11*s13*s16**2*s22*s24*s25*s26*s34*s35*s36 - 4096*q9a*s11**2*s12**2*s23*s24*s25*s26*s34*s35*s36 - 1024*q9a*s11*s12*s16**2*s23*s24*s25*s26*s34*s35*s36 +  \
    1024*q9a*s11**2*s12*s14*s24**2*s25*s26*s34*s35*s36 + 1536*q9a*s11*s12*s14*s15*s22*s25**2*s26*s34*s35*s36 - 2048*q9a*s11*s12*s13*s16*s22*s25**2*s26*s34*s35*s36 - 1536*q9a*s11*s14**2*s16*s22*s25**2*s26*s34*s35*s36 +  \
    2048*q9a*s11*s12**2*s16*s23*s25**2*s26*s34*s35*s36 + 1024*q9a*s11*s12**2*s15*s24*s25**2*s26*s34*s35*s36 - 1024*q9a*s11*s12**2*s14*s25**3*s26*s34*s35*s36 - 512*q9a*s11*s13*s15**2*s22**2*s26**2*s34*s35*s36 -  \
    4096*q9a*s11**2*s12*s13*s22*s23*s26**2*s34*s35*s36 + 1536*q9a*s11**2*s14**2*s22*s23*s26**2*s34*s35*s36 + 512*q9a*s11*s12*s15**2*s22*s23*s26**2*s34*s35*s36 - 1536*q9a*s11*s14*s15*s16*s22*s23*s26**2*s34*s35*s36 +  \
    4096*q9a*s11**2*s12**2*s23**2*s26**2*s34*s35*s36 - 512*q9a*s11*s14*s15**2*s22*s24*s26**2*s34*s35*s36 + 512*q9a*s11*s13*s15*s16*s22*s24*s26**2*s34*s35*s36 - 2048*q9a*s11**2*s12*s14*s23*s24*s26**2*s34*s35*s36 +  \
    1024*q9a*s11*s12*s15**2*s24**2*s26**2*s34*s35*s36 + 2048*q9a*s11*s12*s13*s15*s22*s25*s26**2*s34*s35*s36 + 512*q9a*s11*s14**2*s15*s22*s25*s26**2*s34*s35*s36 + 512*q9a*s11*s13*s14*s16*s22*s25*s26**2*s34*s35*s36 -  \
    2048*q9a*s11*s12**2*s15*s23*s25*s26**2*s34*s35*s36 - 2048*q9a*s11*s12*s14*s15*s24*s25*s26**2*s34*s35*s36 + 2048*q9a*s11*s12*s13*s16*s24*s25*s26**2*s34*s35*s36 + 1024*q9a*s11*s12*s14**2*s25**2*s26**2*s34*s35*s36 +  \
    512*q9a*s11*s13*s14*s15*s22*s26**3*s34*s35*s36 + 1024*q9a*s11*s12*s14*s15*s23*s26**3*s34*s35*s36 - 1024*q9a*s11*s12*s13*s15*s24*s26**3*s34*s35*s36 - 1024*q9a*s11*s12*s13*s14*s25*s26**3*s34*s35*s36 -  \
    4096*q8a*s11**3*s13*s22**3*s33*s34*s35*s36 + 8192*q8a*s11**3*s12*s22**2*s23*s33*s34*s35*s36 - 4096*q7a*s11**3*s22**3*s23*s33*s34*s35*s36 + 2048*q8a*s11**3*s14*s22**2*s24*s33*s34*s35*s36 -  \
    2048*q8a*s11**3*s12*s22*s24**2*s33*s34*s35*s36 - 512*q8a*s11**2*s16**2*s22*s24**2*s33*s34*s35*s36 - 2048*q8a*s11**2*s12*s15*s22**2*s25*s33*s34*s35*s36 - 2048*q8a*s11**2*s14*s16*s22**2*s25*s33*s34*s35*s36 +  \
    2048*q7a*s11**2*s15*s22**3*s25*s33*s34*s35*s36 + 2048*q8a*s11**2*s12*s16*s22*s24*s25*s33*s34*s35*s36 + 512*q8a*s11*s16**3*s22*s24*s25*s33*s34*s35*s36 + 2048*q8a*s11**2*s12**2*s22*s25**2*s33*s34*s35*s36 -  \
    1536*q8a*s11*s12*s16**2*s22*s25**2*s33*s34*s35*s36 - 2048*q7a*s11**2*s12*s22**2*s25**2*s33*s34*s35*s36 + 1536*q7a*s11*s16**2*s22**2*s25**2*s33*s34*s35*s36 + 4096*q8a*s11**2*s13*s16*s22**2*s26*s33*s34*s35*s36 -  \
    8192*q8a*s11**2*s12*s16*s22*s23*s26*s33*s34*s35*s36 + 4096*q7a*s11**2*s16*s22**2*s23*s26*s33*s34*s35*s36 - 1024*q8a*s11**2*s14*s16*s22*s24*s26*s33*s34*s35*s36 + 512*q8a*s11*s15*s16**2*s22*s24*s26*s33*s34*s35*s36 +  \
    2048*q8a*s11**2*s12*s16*s24**2*s26*s33*s34*s35*s36 + 4096*q8a*s11**2*s12*s14*s22*s25*s26*s33*s34*s35*s36 + 1024*q8a*s11*s12*s15*s16*s22*s25*s26*s33*s34*s35*s36 + 512*q8a*s11*s14*s16**2*s22*s25*s26*s33*s34*s35*s36 -  \
    2048*q7a*s11**2*s14*s22**2*s25*s26*s33*s34*s35*s36 - 1024*q7a*s11*s15*s16*s22**2*s25*s26*s33*s34*s35*s36 - 4096*q8a*s11**2*s12**2*s24*s25*s26*s33*s34*s35*s36 - 1024*q8a*s11*s12*s16**2*s24*s25*s26*s33*s34*s35*s36 +  \
    2048*q7a*s11**2*s12*s22*s24*s25*s26*s33*s34*s35*s36 - 1536*q7a*s11*s16**2*s22*s24*s25*s26*s33*s34*s35*s36 + 2048*q8a*s11*s12**2*s16*s25**2*s26*s33*s34*s35*s36 - 2048*q7a*s11*s12*s16*s22*s25**2*s26*s33*s34*s35*s36 -  \
    4096*q8a*s11**2*s12*s13*s22*s26**2*s33*s34*s35*s36 + 1536*q8a*s11**2*s14**2*s22*s26**2*s33*s34*s35*s36 + 512*q8a*s11*s12*s15**2*s22*s26**2*s33*s34*s35*s36 - 1536*q8a*s11*s14*s15*s16*s22*s26**2*s33*s34*s35*s36 -  \
    512*q7a*s11*s15**2*s22**2*s26**2*s33*s34*s35*s36 + 8192*q8a*s11**2*s12**2*s23*s26**2*s33*s34*s35*s36 - 4096*q7a*s11**2*s12*s22*s23*s26**2*s33*s34*s35*s36 - 2048*q8a*s11**2*s12*s14*s24*s26**2*s33*s34*s35*s36 +  \
    512*q7a*s11*s15*s16*s22*s24*s26**2*s33*s34*s35*s36 - 2048*q8a*s11*s12**2*s15*s25*s26**2*s33*s34*s35*s36 + 2048*q7a*s11*s12*s15*s22*s25*s26**2*s33*s34*s35*s36 + 512*q7a*s11*s14*s16*s22*s25*s26**2*s33*s34*s35*s36 +  \
    2048*q7a*s11*s12*s16*s24*s25*s26**2*s33*s34*s35*s36 + 1024*q8a*s11*s12*s14*s15*s26**3*s33*s34*s35*s36 + 512*q7a*s11*s14*s15*s22*s26**3*s33*s34*s35*s36 - 1024*q7a*s11*s12*s15*s24*s26**3*s33*s34*s35*s36 -  \
    1024*q7a*s11*s12*s14*s25*s26**3*s33*s34*s35*s36 - 4096*q8a*s11**3*s14*s22**2*s23*s34**2*s35*s36 + 1024*q8a*s11**3*s13*s22**2*s24*s34**2*s35*s36 + 512*q8a*s11**2*s15**2*s22**2*s24*s34**2*s35*s36 +  \
    2048*q8a*s11**3*s12*s22*s23*s24*s34**2*s35*s36 + 1024*q8a*s11**2*s16**2*s22*s23*s24*s34**2*s35*s36 + 1024*q7a*s11**3*s22**2*s23*s24*s34**2*s35*s36 - 512*q8a*s11**2*s15*s16*s22*s24**2*s34**2*s35*s36 +  \
    2560*q8a*s11**2*s14*s15*s22**2*s25*s34**2*s35*s36 - 512*q8a*s11**2*s13*s16*s22**2*s25*s34**2*s35*s36 - 512*q8a*s11*s15**2*s16*s22**2*s25*s34**2*s35*s36 + 1024*q8a*s11**2*s12*s16*s22*s23*s25*s34**2*s35*s36 -  \
    1024*q8a*s11*s16**3*s22*s23*s25*s34**2*s35*s36 - 512*q7a*s11**2*s16*s22**2*s23*s25*s34**2*s35*s36 - 2048*q8a*s11**2*s12*s15*s22*s24*s25*s34**2*s35*s36 - 512*q8a*s11**2*s14*s16*s22*s24*s25*s34**2*s35*s36 +  \
    512*q8a*s11*s15*s16**2*s22*s24*s25*s34**2*s35*s36 - 1536*q7a*s11**2*s15*s22**2*s24*s25*s34**2*s35*s36 + 1024*q7a*s11**2*s16*s22*s24**2*s25*s34**2*s35*s36 - 1024*q8a*s11**2*s12*s14*s22*s25**2*s34**2*s35*s36 +  \
    512*q8a*s11*s12*s15*s16*s22*s25**2*s34**2*s35*s36 + 512*q8a*s11*s14*s16**2*s22*s25**2*s34**2*s35*s36 - 1536*q7a*s11**2*s14*s22**2*s25**2*s34**2*s35*s36 + 512*q7a*s11*s15*s16*s22**2*s25**2*s34**2*s35*s36 +  \
    3072*q7a*s11**2*s12*s22*s24*s25**2*s34**2*s35*s36 - 1024*q7a*s11*s16**2*s22*s24*s25**2*s34**2*s35*s36 - 512*q7a*s11*s12*s16*s22*s25**3*s34**2*s35*s36 + 512*q8a*s11**2*s13*s15*s22**2*s26*s34**2*s35*s36 -  \
    512*q8a*s11*s15**3*s22**2*s26*s34**2*s35*s36 - 1024*q8a*s11**2*s12*s15*s22*s23*s26*s34**2*s35*s36 + 5120*q8a*s11**2*s14*s16*s22*s23*s26*s34**2*s35*s36 - 1024*q8a*s11*s15*s16**2*s22*s23*s26*s34**2*s35*s36 +  \
    512*q7a*s11**2*s15*s22**2*s23*s26*s34**2*s35*s36 - 512*q8a*s11**2*s14*s15*s22*s24*s26*s34**2*s35*s36 - 2048*q8a*s11**2*s13*s16*s22*s24*s26*s34**2*s35*s36 + 512*q8a*s11*s15**2*s16*s22*s24*s26*s34**2*s35*s36 -  \
    3072*q8a*s11**2*s12*s16*s23*s24*s26*s34**2*s35*s36 - 2048*q7a*s11**2*s16*s22*s23*s24*s26*s34**2*s35*s36 + 1024*q8a*s11**2*s12*s15*s24**2*s26*s34**2*s35*s36 - 2048*q8a*s11**2*s12*s13*s22*s25*s26*s34**2*s35*s36 -  \
    512*q8a*s11**2*s14**2*s22*s25*s26*s34**2*s35*s36 + 1536*q8a*s11*s12*s15**2*s22*s25*s26*s34**2*s35*s36 - 2048*q8a*s11*s14*s15*s16*s22*s25*s26*s34**2*s35*s36 + 1536*q8a*s11*s13*s16**2*s22*s25*s26*s34**2*s35*s36 +  \
    2048*q7a*s11**2*s13*s22**2*s25*s26*s34**2*s35*s36 + 512*q7a*s11*s15**2*s22**2*s25*s26*s34**2*s35*s36 + 2048*q8a*s11**2*s12**2*s23*s25*s26*s34**2*s35*s36 + 1024*q8a*s11*s12*s16**2*s23*s25*s26*s34**2*s35*s36 -  \
    2048*q7a*s11**2*s12*s22*s23*s25*s26*s34**2*s35*s36 + 1536*q7a*s11*s16**2*s22*s23*s25*s26*s34**2*s35*s36 + 1024*q8a*s11**2*s12*s14*s24*s25*s26*s34**2*s35*s36 + 1024*q7a*s11**2*s14*s22*s24*s25*s26*s34**2*s35*s36 -  \
    2048*q7a*s11**2*s12*s24**2*s25*s26*s34**2*s35*s36 - 1024*q8a*s11*s12**2*s15*s25**2*s26*s34**2*s35*s36 - 512*q8a*s11*s12*s14*s16*s25**2*s26*s34**2*s35*s36 - 1536*q7a*s11*s12*s15*s22*s25**2*s26*s34**2*s35*s36 +  \
    1536*q7a*s11*s14*s16*s22*s25**2*s26*s34**2*s35*s36 + 512*q7a*s11*s12*s16*s24*s25**2*s26*s34**2*s35*s36 + 1024*q7a*s11*s12**2*s25**3*s26*s34**2*s35*s36 - 1024*q8a*s11**2*s13*s14*s22*s26**2*s34**2*s35*s36 +  \
    512*q8a*s11*s14*s15**2*s22*s26**2*s34**2*s35*s36 + 512*q8a*s11*s13*s15*s16*s22*s26**2*s34**2*s35*s36 - 3072*q8a*s11**2*s12*s14*s23*s26**2*s34**2*s35*s36 + 1024*q8a*s11*s12*s15*s16*s23*s26**2*s34**2*s35*s36 -  \
    1024*q7a*s11**2*s14*s22*s23*s26**2*s34**2*s35*s36 + 512*q7a*s11*s15*s16*s22*s23*s26**2*s34**2*s35*s36 + 3072*q8a*s11**2*s12*s13*s24*s26**2*s34**2*s35*s36 - 1024*q8a*s11*s12*s15**2*s24*s26**2*s34**2*s35*s36 +  \
    3072*q7a*s11**2*s12*s23*s24*s26**2*s34**2*s35*s36 + 1536*q8a*s11*s12*s14*s15*s25*s26**2*s34**2*s35*s36 - 1536*q8a*s11*s12*s13*s16*s25*s26**2*s34**2*s35*s36 - 512*q7a*s11*s14*s15*s22*s25*s26**2*s34**2*s35*s36 -  \
    2048*q7a*s11*s13*s16*s22*s25*s26**2*s34**2*s35*s36 - 1536*q7a*s11*s12*s16*s23*s25*s26**2*s34**2*s35*s36 + 512*q7a*s11*s12*s15*s24*s25*s26**2*s34**2*s35*s36 - 1024*q7a*s11*s12*s14*s25**2*s26**2*s34**2*s35*s36 -  \
    512*q8a*s11*s12*s13*s15*s26**3*s34**2*s35*s36 - 512*q7a*s11*s12*s15*s23*s26**3*s34**2*s35*s36 + 2048*q7a*s11*s12*s13*s25*s26**3*s34**2*s35*s36 - 1024*q9a*s11**2*s14*s15*s22**3*s23*s35**2*s36 +  \
    1024*q9a*s11**2*s13*s16*s22**3*s23*s35**2*s36 - 1024*q9a*s11**2*s12*s16*s22**2*s23**2*s35**2*s36 - 1536*q9a*s11**2*s13*s15*s22**3*s24*s35**2*s36 + 2560*q9a*s11**2*s12*s15*s22**2*s23*s24*s35**2*s36 -  \
    1536*q9a*s11**2*s14*s16*s22**2*s23*s24*s35**2*s36 + 1024*q9a*s11**2*s14*s15*s22**2*s24**2*s35**2*s36 + 512*q9a*s11**2*s13*s16*s22**2*s24**2*s35**2*s36 + 1024*q9a*s11**2*s12*s16*s22*s23*s24**2*s35**2*s36 -  \
    1024*q9a*s11**2*s12*s15*s22*s24**3*s35**2*s36 + 512*q9a*s11**2*s13*s14*s22**3*s25*s35**2*s36 + 1024*q9a*s11*s13*s15*s16*s22**3*s25*s35**2*s36 + 512*q9a*s11**2*s12*s14*s22**2*s23*s25*s35**2*s36 -  \
    1024*q9a*s11*s12*s15*s16*s22**2*s23*s25*s35**2*s36 + 1024*q9a*s11*s14*s16**2*s22**2*s23*s25*s35**2*s36 + 1024*q9a*s11**2*s12*s13*s22**2*s24*s25*s35**2*s36 - 1024*q9a*s11**2*s14**2*s22**2*s24*s25*s35**2*s36 -  \
    512*q9a*s11*s14*s15*s16*s22**2*s24*s25*s35**2*s36 - 512*q9a*s11*s13*s16**2*s22**2*s24*s25*s35**2*s36 - 2048*q9a*s11**2*s12**2*s22*s23*s24*s25*s35**2*s36 - 512*q9a*s11*s12*s16**2*s22*s23*s24*s25*s35**2*s36 +  \
    1024*q9a*s11**2*s12*s14*s22*s24**2*s25*s35**2*s36 + 512*q9a*s11*s12*s15*s16*s22*s24**2*s25*s35**2*s36 - 1024*q9a*s11*s12*s13*s16*s22**2*s25**2*s35**2*s36 + 512*q9a*s11*s14**2*s16*s22**2*s25**2*s35**2*s36 +  \
    1024*q9a*s11*s12**2*s16*s22*s23*s25**2*s35**2*s36 - 512*q9a*s11*s12*s14*s16*s22*s24*s25**2*s35**2*s36 - 1024*q9a*s11**2*s13**2*s22**3*s26*s35**2*s36 + 1024*q9a*s11*s13*s15**2*s22**3*s26*s35**2*s36 +  \
    1024*q9a*s11**2*s12*s13*s22**2*s23*s26*s35**2*s36 - 512*q9a*s11**2*s14**2*s22**2*s23*s26*s35**2*s36 - 1024*q9a*s11*s12*s15**2*s22**2*s23*s26*s35**2*s36 + 2048*q9a*s11*s14*s15*s16*s22**2*s23*s26*s35**2*s36 -  \
    1024*q9a*s11*s13*s16**2*s22**2*s23*s26*s35**2*s36 + 1024*q9a*s11*s12*s16**2*s22*s23**2*s26*s35**2*s36 + 1536*q9a*s11**2*s13*s14*s22**2*s24*s26*s35**2*s36 - 512*q9a*s11*s14*s15**2*s22**2*s24*s26*s35**2*s36 +  \
    1024*q9a*s11**2*s12*s14*s22*s23*s24*s26*s35**2*s36 - 2048*q9a*s11*s12*s15*s16*s22*s23*s24*s26*s35**2*s36 + 512*q9a*s11*s14*s16**2*s22*s23*s24*s26*s35**2*s36 - 2048*q9a*s11**2*s12*s13*s22*s24**2*s26*s35**2*s36;
    v3_5= \
    512*q9a*s11*s12*s15**2*s22*s24**2*s26*s35**2*s36 - 512*q9a*s11*s14*s15*s16*s22*s24**2*s26*s35**2*s36 - 512*q9a*s11*s12*s16**2*s23*s24**2*s26*s35**2*s36 + 512*q9a*s11*s12*s15*s16*s24**3*s26*s35**2*s36 -  \
    3072*q9a*s11*s12*s13*s15*s22**2*s25*s26*s35**2*s36 + 512*q9a*s11*s14**2*s15*s22**2*s25*s26*s35**2*s36 - 1024*q9a*s11*s13*s14*s16*s22**2*s25*s26*s35**2*s36 + 3072*q9a*s11*s12**2*s15*s22*s23*s25*s26*s35**2*s36 -  \
    3072*q9a*s11*s12*s14*s16*s22*s23*s25*s26*s35**2*s36 + 512*q9a*s11*s12*s14*s15*s22*s24*s25*s26*s35**2*s36 + 2048*q9a*s11*s12*s13*s16*s22*s24*s25*s26*s35**2*s36 + 512*q9a*s11*s14**2*s16*s22*s24*s25*s26*s35**2*s36 +  \
    2048*q9a*s11*s12**2*s16*s23*s24*s25*s26*s35**2*s36 - 1024*q9a*s11*s12**2*s15*s24**2*s25*s26*s35**2*s36 - 512*q9a*s11*s12*s14*s16*s24**2*s25*s26*s35**2*s36 + 2048*q9a*s11*s12**2*s13*s22*s25**2*s26*s35**2*s36 -  \
    1024*q9a*s11*s12*s14**2*s22*s25**2*s26*s35**2*s36 - 2048*q9a*s11*s12**3*s23*s25**2*s26*s35**2*s36 + 1024*q9a*s11*s12**2*s14*s24*s25**2*s26*s35**2*s36 - 1536*q9a*s11*s13*s14*s15*s22**2*s26**2*s35**2*s36 +  \
    1024*q9a*s11*s13**2*s16*s22**2*s26**2*s35**2*s36 - 512*q9a*s11*s12*s14*s15*s22*s23*s26**2*s35**2*s36 - 512*q9a*s11*s14**2*s16*s22*s23*s26**2*s35**2*s36 - 1024*q9a*s11*s12**2*s16*s23**2*s26**2*s35**2*s36 +  \
    1536*q9a*s11*s12*s13*s15*s22*s24*s26**2*s35**2*s36 + 512*q9a*s11*s14**2*s15*s22*s24*s26**2*s35**2*s36 - 512*q9a*s11*s13*s14*s16*s22*s24*s26**2*s35**2*s36 + 512*q9a*s11*s12**2*s15*s23*s24*s26**2*s35**2*s36 +  \
    512*q9a*s11*s12*s14*s16*s23*s24*s26**2*s35**2*s36 - 512*q9a*s11*s12*s14*s15*s24**2*s26**2*s35**2*s36 + 512*q9a*s11*s12*s13*s16*s24**2*s26**2*s35**2*s36 + 2560*q9a*s11*s12*s13*s14*s22*s25*s26**2*s35**2*s36 -  \
    512*q9a*s11*s14**3*s22*s25*s26**2*s35**2*s36 + 512*q9a*s11*s12**2*s14*s23*s25*s26**2*s35**2*s36 - 3072*q9a*s11*s12**2*s13*s24*s25*s26**2*s35**2*s36 + 512*q9a*s11*s12*s14**2*s24*s25*s26**2*s35**2*s36 -  \
    1024*q9a*s11*s12*s13**2*s22*s26**3*s35**2*s36 + 512*q9a*s11*s13*s14**2*s22*s26**3*s35**2*s36 + 1024*q9a*s11*s12**2*s13*s23*s26**3*s35**2*s36 - 512*q9a*s11*s12*s13*s14*s24*s26**3*s35**2*s36 -  \
    1024*q8a*s11**2*s14*s15*s22**3*s33*s35**2*s36 + 1024*q8a*s11**2*s13*s16*s22**3*s33*s35**2*s36 - 2048*q8a*s11**2*s12*s16*s22**2*s23*s33*s35**2*s36 + 1024*q7a*s11**2*s16*s22**3*s23*s33*s35**2*s36 +  \
    2560*q8a*s11**2*s12*s15*s22**2*s24*s33*s35**2*s36 - 1536*q8a*s11**2*s14*s16*s22**2*s24*s33*s35**2*s36 - 1536*q7a*s11**2*s15*s22**3*s24*s33*s35**2*s36 + 1024*q8a*s11**2*s12*s16*s22*s24**2*s33*s35**2*s36 +  \
    512*q7a*s11**2*s16*s22**2*s24**2*s33*s35**2*s36 + 512*q8a*s11**2*s12*s14*s22**2*s25*s33*s35**2*s36 - 1024*q8a*s11*s12*s15*s16*s22**2*s25*s33*s35**2*s36 + 1024*q8a*s11*s14*s16**2*s22**2*s25*s33*s35**2*s36 +  \
    512*q7a*s11**2*s14*s22**3*s25*s33*s35**2*s36 + 1024*q7a*s11*s15*s16*s22**3*s25*s33*s35**2*s36 - 2048*q8a*s11**2*s12**2*s22*s24*s25*s33*s35**2*s36 - 512*q8a*s11*s12*s16**2*s22*s24*s25*s33*s35**2*s36 +  \
    1024*q7a*s11**2*s12*s22**2*s24*s25*s33*s35**2*s36 - 512*q7a*s11*s16**2*s22**2*s24*s25*s33*s35**2*s36 + 1024*q8a*s11*s12**2*s16*s22*s25**2*s33*s35**2*s36 - 1024*q7a*s11*s12*s16*s22**2*s25**2*s33*s35**2*s36 +  \
    1024*q8a*s11**2*s12*s13*s22**2*s26*s33*s35**2*s36 - 512*q8a*s11**2*s14**2*s22**2*s26*s33*s35**2*s36 - 1024*q8a*s11*s12*s15**2*s22**2*s26*s33*s35**2*s36 + 2048*q8a*s11*s14*s15*s16*s22**2*s26*s33*s35**2*s36 -  \
    1024*q8a*s11*s13*s16**2*s22**2*s26*s33*s35**2*s36 - 2048*q7a*s11**2*s13*s22**3*s26*s33*s35**2*s36 + 1024*q7a*s11*s15**2*s22**3*s26*s33*s35**2*s36 + 2048*q8a*s11*s12*s16**2*s22*s23*s26*s33*s35**2*s36 +  \
    1024*q7a*s11**2*s12*s22**2*s23*s26*s33*s35**2*s36 - 1024*q7a*s11*s16**2*s22**2*s23*s26*s33*s35**2*s36 + 1024*q8a*s11**2*s12*s14*s22*s24*s26*s33*s35**2*s36 - 2048*q8a*s11*s12*s15*s16*s22*s24*s26*s33*s35**2*s36 +  \
    512*q8a*s11*s14*s16**2*s22*s24*s26*s33*s35**2*s36 + 1536*q7a*s11**2*s14*s22**2*s24*s26*s33*s35**2*s36 - 512*q8a*s11*s12*s16**2*s24**2*s26*s33*s35**2*s36 - 2048*q7a*s11**2*s12*s22*s24**2*s26*s33*s35**2*s36 +  \
    3072*q8a*s11*s12**2*s15*s22*s25*s26*s33*s35**2*s36 - 3072*q8a*s11*s12*s14*s16*s22*s25*s26*s33*s35**2*s36 - 3072*q7a*s11*s12*s15*s22**2*s25*s26*s33*s35**2*s36 - 1024*q7a*s11*s14*s16*s22**2*s25*s26*s33*s35**2*s36 +  \
    2048*q8a*s11*s12**2*s16*s24*s25*s26*s33*s35**2*s36 + 2048*q7a*s11*s12*s16*s22*s24*s25*s26*s33*s35**2*s36 - 2048*q8a*s11*s12**3*s25**2*s26*s33*s35**2*s36 + 2048*q7a*s11*s12**2*s22*s25**2*s26*s33*s35**2*s36 -  \
    512*q8a*s11*s12*s14*s15*s22*s26**2*s33*s35**2*s36 - 512*q8a*s11*s14**2*s16*s22*s26**2*s33*s35**2*s36 - 1536*q7a*s11*s14*s15*s22**2*s26**2*s33*s35**2*s36 + 2048*q7a*s11*s13*s16*s22**2*s26**2*s33*s35**2*s36 -  \
    2048*q8a*s11*s12**2*s16*s23*s26**2*s33*s35**2*s36 + 512*q8a*s11*s12**2*s15*s24*s26**2*s33*s35**2*s36 + 512*q8a*s11*s12*s14*s16*s24*s26**2*s33*s35**2*s36 + 1536*q7a*s11*s12*s15*s22*s24*s26**2*s33*s35**2*s36 -  \
    512*q7a*s11*s14*s16*s22*s24*s26**2*s33*s35**2*s36 + 512*q7a*s11*s12*s16*s24**2*s26**2*s33*s35**2*s36 + 512*q8a*s11*s12**2*s14*s25*s26**2*s33*s35**2*s36 + 2560*q7a*s11*s12*s14*s22*s25*s26**2*s33*s35**2*s36 -  \
    3072*q7a*s11*s12**2*s24*s25*s26**2*s33*s35**2*s36 + 1024*q8a*s11*s12**2*s13*s26**3*s33*s35**2*s36 - 2048*q7a*s11*s12*s13*s22*s26**3*s33*s35**2*s36 + 512*q7a*s11*s14**2*s22*s26**3*s33*s35**2*s36 +  \
    1024*q7a*s11*s12**2*s23*s26**3*s33*s35**2*s36 - 512*q7a*s11*s12*s14*s24*s26**3*s33*s35**2*s36 + 1024*q8a*s11**2*s13*s15*s22**3*s34*s35**2*s36 - 2048*q8a*s11**2*s12*s15*s22**2*s23*s34*s35**2*s36 +  \
    2048*q8a*s11**2*s14*s16*s22**2*s23*s34*s35**2*s36 + 1024*q7a*s11**2*s15*s22**3*s23*s34*s35**2*s36 - 1536*q8a*s11**2*s14*s15*s22**2*s24*s34*s35**2*s36 + 512*q8a*s11**2*s13*s16*s22**2*s24*s34*s35**2*s36 -  \
    3072*q8a*s11**2*s12*s16*s22*s23*s24*s34*s35**2*s36 + 512*q7a*s11**2*s16*s22**2*s23*s24*s34*s35**2*s36 + 1024*q8a*s11**2*s12*s15*s22*s24**2*s34*s35**2*s36 + 512*q8a*s11**2*s14*s16*s22*s24**2*s34*s35**2*s36 +  \
    512*q7a*s11**2*s15*s22**2*s24**2*s34*s35**2*s36 - 512*q7a*s11**2*s16*s22*s24**3*s34*s35**2*s36 + 3072*q8a*s11**2*s12*s13*s22**2*s25*s34*s35**2*s36 - 1536*q8a*s11**2*s14**2*s22**2*s25*s34*s35**2*s36 +  \
    1024*q8a*s11*s14*s15*s16*s22**2*s25*s34*s35**2*s36 - 1024*q8a*s11*s13*s16**2*s22**2*s25*s34*s35**2*s36 - 4096*q7a*s11**2*s13*s22**3*s25*s34*s35**2*s36 - 2048*q8a*s11**2*s12**2*s22*s23*s25*s34*s35**2*s36 +  \
    2048*q8a*s11*s12*s16**2*s22*s23*s25*s34*s35**2*s36 + 3072*q7a*s11**2*s12*s22**2*s23*s25*s34*s35**2*s36 - 1024*q7a*s11*s16**2*s22**2*s23*s25*s34*s35**2*s36 + 2048*q8a*s11**2*s12*s14*s22*s24*s25*s34*s35**2*s36 -  \
    512*q8a*s11*s12*s15*s16*s22*s24*s25*s34*s35**2*s36 - 512*q8a*s11*s14*s16**2*s22*s24*s25*s34*s35**2*s36 + 2560*q7a*s11**2*s14*s22**2*s24*s25*s34*s35**2*s36 - 512*q7a*s11*s15*s16*s22**2*s24*s25*s34*s35**2*s36 -  \
    3072*q7a*s11**2*s12*s22*s24**2*s25*s34*s35**2*s36 + 512*q7a*s11*s16**2*s22*s24**2*s25*s34*s35**2*s36 - 512*q8a*s11*s12*s14*s16*s22*s25**2*s34*s35**2*s36 - 512*q7a*s11*s14*s16*s22**2*s25**2*s34*s35**2*s36 +  \
    1024*q7a*s11*s12*s16*s22*s24*s25**2*s34*s35**2*s36 - 512*q8a*s11**2*s13*s14*s22**2*s26*s34*s35**2*s36 + 1024*q8a*s11*s14*s15**2*s22**2*s26*s34*s35**2*s36 - 2048*q8a*s11*s13*s15*s16*s22**2*s26*s34*s35**2*s36 -  \
    1024*q8a*s11**2*s12*s14*s22*s23*s26*s34*s35**2*s36 + 4096*q8a*s11*s12*s15*s16*s22*s23*s26*s34*s35**2*s36 - 2048*q8a*s11*s14*s16**2*s22*s23*s26*s34*s35**2*s36 - 512*q7a*s11**2*s14*s22**2*s23*s26*s34*s35**2*s36 -  \
    2048*q7a*s11*s15*s16*s22**2*s23*s26*s34*s35**2*s36 + 512*q8a*s11**2*s14**2*s22*s24*s26*s34*s35**2*s36 - 512*q8a*s11*s12*s15**2*s22*s24*s26*s34*s35**2*s36 + 512*q8a*s11*s13*s16**2*s22*s24*s26*s34*s35**2*s36 -  \
    512*q7a*s11*s15**2*s22**2*s24*s26*s34*s35**2*s36 + 2048*q8a*s11**2*s12**2*s23*s24*s26*s34*s35**2*s36 + 1024*q8a*s11*s12*s16**2*s23*s24*s26*s34*s35**2*s36 + 512*q7a*s11*s16**2*s22*s23*s24*s26*s34*s35**2*s36 -  \
    1024*q8a*s11**2*s12*s14*s24**2*s26*s34*s35**2*s36 - 512*q8a*s11*s12*s15*s16*s24**2*s26*s34*s35**2*s36 - 512*q7a*s11**2*s14*s22*s24**2*s26*s34*s35**2*s36 + 512*q7a*s11*s15*s16*s22*s24**2*s26*s34*s35**2*s36 +  \
    1024*q7a*s11**2*s12*s24**3*s26*s34*s35**2*s36 - 2560*q8a*s11*s12*s14*s15*s22*s25*s26*s34*s35**2*s36 - 1024*q8a*s11*s12*s13*s16*s22*s25*s26*s34*s35**2*s36 + 1536*q8a*s11*s14**2*s16*s22*s25*s26*s34*s35**2*s36 -  \
    512*q7a*s11*s14*s15*s22**2*s25*s26*s34*s35**2*s36 + 5120*q7a*s11*s13*s16*s22**2*s25*s26*s34*s35**2*s36 - 3072*q8a*s11*s12**2*s16*s23*s25*s26*s34*s35**2*s36 - 1024*q7a*s11*s12*s16*s22*s23*s25*s26*s34*s35**2*s36 +  \
    1024*q8a*s11*s12**2*s15*s24*s25*s26*s34*s35**2*s36 + 2048*q7a*s11*s12*s15*s22*s24*s25*s26*s34*s35**2*s36 - 2048*q7a*s11*s14*s16*s22*s24*s25*s26*s34*s35**2*s36 + 512*q7a*s11*s12*s16*s24**2*s25*s26*s34*s35**2*s36 +  \
    1024*q8a*s11*s12**2*s14*s25**2*s26*s34*s35**2*s36 + 1024*q7a*s11*s12*s14*s22*s25**2*s26*s34*s35**2*s36 - 2048*q7a*s11*s12**2*s24*s25**2*s26*s34*s35**2*s36 + 1024*q8a*s11*s12*s13*s15*s22*s26**2*s34*s35**2*s36 -  \
    1024*q8a*s11*s14**2*s15*s22*s26**2*s34*s35**2*s36 + 1536*q8a*s11*s13*s14*s16*s22*s26**2*s34*s35**2*s36 + 1024*q7a*s11*s13*s15*s22**2*s26**2*s34*s35**2*s36 - 3072*q8a*s11*s12**2*s15*s23*s26**2*s34*s35**2*s36 +  \
    1024*q8a*s11*s12*s14*s16*s23*s26**2*s34*s35**2*s36 + 1024*q7a*s11*s12*s15*s22*s23*s26**2*s34*s35**2*s36 + 1536*q7a*s11*s14*s16*s22*s23*s26**2*s34*s35**2*s36 + 1536*q8a*s11*s12*s14*s15*s24*s26**2*s34*s35**2*s36 -  \
    1536*q8a*s11*s12*s13*s16*s24*s26**2*s34*s35**2*s36 + 512*q7a*s11*s14*s15*s22*s24*s26**2*s34*s35**2*s36 - 1024*q7a*s11*s13*s16*s22*s24*s26**2*s34*s35**2*s36 - 1536*q7a*s11*s12*s16*s23*s24*s26**2*s34*s35**2*s36 -  \
    1024*q7a*s11*s12*s15*s24**2*s26**2*s34*s35**2*s36 + 3072*q8a*s11*s12**2*s13*s25*s26**2*s34*s35**2*s36 - 1024*q8a*s11*s12*s14**2*s25*s26**2*s34*s35**2*s36 - 6144*q7a*s11*s12*s13*s22*s25*s26**2*s34*s35**2*s36 +  \
    512*q7a*s11*s14**2*s22*s25*s26**2*s34*s35**2*s36 + 3072*q7a*s11*s12**2*s23*s25*s26**2*s34*s35**2*s36 + 512*q7a*s11*s12*s14*s24*s25*s26**2*s34*s35**2*s36 - 512*q8a*s11*s12*s13*s14*s26**3*s34*s35**2*s36 -  \
    1024*q7a*s11*s13*s14*s22*s26**3*s34*s35**2*s36 - 512*q7a*s11*s12*s14*s23*s26**3*s34*s35**2*s36 + 2048*q7a*s11*s12*s13*s24*s26**3*s34*s35**2*s36 - 2048*q8a*s11**2*s12*s13*s22**2*s24*s35**3*s36 +  \
    1024*q8a*s11**2*s14**2*s22**2*s24*s35**3*s36 + 2048*q7a*s11**2*s13*s22**3*s24*s35**3*s36 + 2048*q8a*s11**2*s12**2*s22*s23*s24*s35**3*s36 - 2048*q7a*s11**2*s12*s22**2*s23*s24*s35**3*s36 -  \
    1024*q8a*s11**2*s12*s14*s22*s24**2*s35**3*s36 - 1024*q7a*s11**2*s14*s22**2*s24**2*s35**3*s36 + 1024*q7a*s11**2*s12*s22*s24**3*s35**3*s36 + 1024*q8a*s11*s12*s13*s16*s22**2*s25*s35**3*s36 -  \
    512*q8a*s11*s14**2*s16*s22**2*s25*s35**3*s36 - 1024*q7a*s11*s13*s16*s22**3*s25*s35**3*s36 - 1024*q8a*s11*s12**2*s16*s22*s23*s25*s35**3*s36 + 1024*q7a*s11*s12*s16*s22**2*s23*s25*s35**3*s36 +  \
    512*q8a*s11*s12*s14*s16*s22*s24*s25*s35**3*s36 + 512*q7a*s11*s14*s16*s22**2*s24*s25*s35**3*s36 - 512*q7a*s11*s12*s16*s22*s24**2*s25*s35**3*s36 + 1024*q8a*s11*s12*s13*s15*s22**2*s26*s35**3*s36 -  \
    512*q8a*s11*s14**2*s15*s22**2*s26*s35**3*s36 - 1024*q7a*s11*s13*s15*s22**3*s26*s35**3*s36 - 1024*q8a*s11*s12**2*s15*s22*s23*s26*s35**3*s36 + 1024*q7a*s11*s12*s15*s22**2*s23*s26*s35**3*s36 +  \
    512*q8a*s11*s12*s14*s15*s22*s24*s26*s35**3*s36 + 1024*q8a*s11*s12*s13*s16*s22*s24*s26*s35**3*s36 - 512*q8a*s11*s14**2*s16*s22*s24*s26*s35**3*s36 + 512*q7a*s11*s14*s15*s22**2*s24*s26*s35**3*s36 -  \
    1024*q7a*s11*s13*s16*s22**2*s24*s26*s35**3*s36 - 1024*q8a*s11*s12**2*s16*s23*s24*s26*s35**3*s36 + 1024*q7a*s11*s12*s16*s22*s23*s24*s26*s35**3*s36 + 512*q8a*s11*s12*s14*s16*s24**2*s26*s35**3*s36 -  \
    512*q7a*s11*s12*s15*s22*s24**2*s26*s35**3*s36 + 512*q7a*s11*s14*s16*s22*s24**2*s26*s35**3*s36 - 512*q7a*s11*s12*s16*s24**3*s26*s35**3*s36 - 2048*q8a*s11*s12**2*s13*s22*s25*s26*s35**3*s36 +  \
    1024*q8a*s11*s12*s14**2*s22*s25*s26*s35**3*s36 + 2048*q7a*s11*s12*s13*s22**2*s25*s26*s35**3*s36 + 2048*q8a*s11*s12**3*s23*s25*s26*s35**3*s36 - 2048*q7a*s11*s12**2*s22*s23*s25*s26*s35**3*s36 -  \
    1024*q8a*s11*s12**2*s14*s24*s25*s26*s35**3*s36 - 1024*q7a*s11*s12*s14*s22*s24*s25*s26*s35**3*s36 + 1024*q7a*s11*s12**2*s24**2*s25*s26*s35**3*s36 - 1024*q8a*s11*s12*s13*s14*s22*s26**2*s35**3*s36;
    v3_6= \
    512*q8a*s11*s14**3*s22*s26**2*s35**3*s36 + 1024*q7a*s11*s13*s14*s22**2*s26**2*s35**3*s36 + 1024*q8a*s11*s12**2*s14*s23*s26**2*s35**3*s36 - 1024*q7a*s11*s12*s14*s22*s23*s26**2*s35**3*s36 -  \
    512*q8a*s11*s12*s14**2*s24*s26**2*s35**3*s36 - 512*q7a*s11*s14**2*s22*s24*s26**2*s35**3*s36 + 512*q7a*s11*s12*s14*s24**2*s26**2*s35**3*s36 + 2048*q9a*s11**3*s13*s22**3*s23*s33*s36**2 - 1024*q9a*s11**2*s15**2*s22**3*s23*s33*s36**2 -  \
    2048*q9a*s11**3*s12*s22**2*s23**2*s33*s36**2 - 1024*q9a*s11**2*s16**2*s22**2*s23**2*s33*s36**2 - 3072*q9a*s11**3*s14*s22**2*s23*s24*s33*s36**2 + 1024*q9a*s11**2*s15*s16*s22**2*s23*s24*s33*s36**2 -  \
    1024*q9a*s11**3*s13*s22**2*s24**2*s33*s36**2 + 4096*q9a*s11**3*s12*s22*s23*s24**2*s33*s36**2 + 1024*q9a*s11**3*s14*s22*s24**3*s33*s36**2 - 1024*q9a*s11**3*s12*s24**4*s33*s36**2 - 3072*q9a*s11**2*s13*s15*s22**3*s25*s33*s36**2 +  \
    1024*q9a*s11*s15**3*s22**3*s25*s33*s36**2 + 5120*q9a*s11**2*s12*s15*s22**2*s23*s25*s33*s36**2 - 1024*q9a*s11**2*s14*s16*s22**2*s23*s25*s33*s36**2 + 1024*q9a*s11*s15*s16**2*s22**2*s23*s25*s33*s36**2 +  \
    1024*q9a*s11**2*s14*s15*s22**2*s24*s25*s33*s36**2 + 4096*q9a*s11**2*s13*s16*s22**2*s24*s25*s33*s36**2 - 1024*q9a*s11*s15**2*s16*s22**2*s24*s25*s33*s36**2 - 4096*q9a*s11**2*s12*s16*s22*s23*s24*s25*s33*s36**2 -  \
    1024*q9a*s11**2*s12*s15*s22*s24**2*s25*s33*s36**2 - 2048*q9a*s11**2*s14*s16*s22*s24**2*s25*s33*s36**2 + 2048*q9a*s11**2*s12*s16*s24**3*s25*s33*s36**2 + 3072*q9a*s11**2*s12*s13*s22**2*s25**2*s33*s36**2 -  \
    3072*q9a*s11**2*s14**2*s22**2*s25**2*s33*s36**2 - 3072*q9a*s11*s12*s15**2*s22**2*s25**2*s33*s36**2 + 3072*q9a*s11*s14*s15*s16*s22**2*s25**2*s33*s36**2 - 3072*q9a*s11*s13*s16**2*s22**2*s25**2*s33*s36**2 -  \
    4096*q9a*s11**2*s12**2*s22*s23*s25**2*s33*s36**2 + 2048*q9a*s11*s12*s16**2*s22*s23*s25**2*s33*s36**2 + 5120*q9a*s11**2*s12*s14*s22*s24*s25**2*s33*s36**2 - 1024*q9a*s11*s12*s15*s16*s22*s24*s25**2*s33*s36**2 +  \
    1024*q9a*s11*s14*s16**2*s22*s24*s25**2*s33*s36**2 - 2048*q9a*s11**2*s12**2*s24**2*s25**2*s33*s36**2 - 1024*q9a*s11*s12*s16**2*s24**2*s25**2*s33*s36**2 + 3072*q9a*s11*s12**2*s15*s22*s25**3*s33*s36**2 -  \
    3072*q9a*s11*s12*s14*s16*s22*s25**3*s33*s36**2 + 2048*q9a*s11*s12**2*s16*s24*s25**3*s33*s36**2 - 1024*q9a*s11*s12**3*s25**4*s33*s36**2 + 4096*q9a*s11**2*s14*s15*s22**2*s23*s26*s33*s36**2 -  \
    4096*q9a*s11**2*s13*s16*s22**2*s23*s26*s33*s36**2 + 6144*q9a*s11**2*s12*s16*s22*s23**2*s26*s33*s36**2 + 1024*q9a*s11**2*s13*s15*s22**2*s24*s26*s33*s36**2 - 6144*q9a*s11**2*s12*s15*s22*s23*s24*s26*s33*s36**2 +  \
    2048*q9a*s11**2*s14*s16*s22*s23*s24*s26*s33*s36**2 - 1024*q9a*s11**2*s14*s15*s22*s24**2*s26*s33*s36**2 - 2048*q9a*s11**2*s12*s16*s23*s24**2*s26*s33*s36**2 + 1024*q9a*s11**2*s12*s15*s24**3*s26*s33*s36**2 -  \
    1024*q9a*s11**2*s13*s14*s22**2*s25*s26*s33*s36**2 - 2048*q9a*s11*s14*s15**2*s22**2*s25*s26*s33*s36**2 + 2048*q9a*s11*s13*s15*s16*s22**2*s25*s26*s33*s36**2 - 2048*q9a*s11**2*s12*s14*s22*s23*s25*s26*s33*s36**2 -  \
    4096*q9a*s11*s12*s15*s16*s22*s23*s25*s26*s33*s36**2 - 4096*q9a*s11**2*s12*s13*s22*s24*s25*s26*s33*s36**2 + 3072*q9a*s11**2*s14**2*s22*s24*s25*s26*s33*s36**2 + 3072*q9a*s11*s12*s15**2*s22*s24*s25*s26*s33*s36**2 -  \
    1024*q9a*s11*s14*s15*s16*s22*s24*s25*s26*s33*s36**2 + 8192*q9a*s11**2*s12**2*s23*s24*s25*s26*s33*s36**2 - 3072*q9a*s11**2*s12*s14*s24**2*s25*s26*s33*s36**2 + 1024*q9a*s11*s12*s15*s16*s24**2*s25*s26*s33*s36**2 +  \
    1024*q9a*s11*s12*s14*s15*s22*s25**2*s26*s33*s36**2 + 4096*q9a*s11*s12*s13*s16*s22*s25**2*s26*s33*s36**2 - 1024*q9a*s11*s14**2*s16*s22*s25**2*s26*s33*s36**2 - 2048*q9a*s11*s12**2*s16*s23*s25**2*s26*s33*s36**2 -  \
    3072*q9a*s11*s12**2*s15*s24*s25**2*s26*s33*s36**2 + 1024*q9a*s11*s12*s14*s16*s24*s25**2*s26*s33*s36**2 + 1024*q9a*s11*s12**2*s14*s25**3*s26*s33*s36**2 - 1024*q9a*s11**2*s13**2*s22**2*s26**2*s33*s36**2 +  \
    6144*q9a*s11**2*s12*s13*s22*s23*s26**2*s33*s36**2 - 3072*q9a*s11**2*s14**2*s22*s23*s26**2*s33*s36**2 - 6144*q9a*s11**2*s12**2*s23**2*s26**2*s33*s36**2 + 1024*q9a*s11**2*s13*s14*s22*s24*s26**2*s33*s36**2 +  \
    3072*q9a*s11**2*s12*s14*s23*s24*s26**2*s33*s36**2 - 1024*q9a*s11**2*s12*s13*s24**2*s26**2*s33*s36**2 - 2048*q9a*s11*s12*s13*s15*s22*s25*s26**2*s33*s36**2 + 1024*q9a*s11*s14**2*s15*s22*s25*s26**2*s33*s36**2 +  \
    3072*q9a*s11*s12**2*s15*s23*s25*s26**2*s33*s36**2 - 1024*q9a*s11*s12*s14*s15*s24*s25*s26**2*s33*s36**2 - 1024*q9a*s11*s12**2*s13*s25**2*s26**2*s33*s36**2 + 1024*q8a*s11**3*s13*s22**3*s33**2*s36**2 -  \
    512*q8a*s11**2*s15**2*s22**3*s33**2*s36**2 - 2048*q8a*s11**3*s12*s22**2*s23*s33**2*s36**2 - 1024*q8a*s11**2*s16**2*s22**2*s23*s33**2*s36**2 + 1024*q7a*s11**3*s22**3*s23*s33**2*s36**2 - 1536*q8a*s11**3*s14*s22**2*s24*s33**2*s36**2 +  \
    512*q8a*s11**2*s15*s16*s22**2*s24*s33**2*s36**2 + 2048*q8a*s11**3*s12*s22*s24**2*s33**2*s36**2 - 512*q7a*s11**3*s22**2*s24**2*s33**2*s36**2 + 2560*q8a*s11**2*s12*s15*s22**2*s25*s33**2*s36**2 -  \
    512*q8a*s11**2*s14*s16*s22**2*s25*s33**2*s36**2 + 512*q8a*s11*s15*s16**2*s22**2*s25*s33**2*s36**2 - 1536*q7a*s11**2*s15*s22**3*s25*s33**2*s36**2 - 2048*q8a*s11**2*s12*s16*s22*s24*s25*s33**2*s36**2 +  \
    2048*q7a*s11**2*s16*s22**2*s24*s25*s33**2*s36**2 - 2048*q8a*s11**2*s12**2*s22*s25**2*s33**2*s36**2 + 1024*q8a*s11*s12*s16**2*s22*s25**2*s33**2*s36**2 + 1536*q7a*s11**2*s12*s22**2*s25**2*s33**2*s36**2 -  \
    1536*q7a*s11*s16**2*s22**2*s25**2*s33**2*s36**2 + 2048*q8a*s11**2*s14*s15*s22**2*s26*s33**2*s36**2 - 2048*q8a*s11**2*s13*s16*s22**2*s26*s33**2*s36**2 + 6144*q8a*s11**2*s12*s16*s22*s23*s26*s33**2*s36**2 -  \
    2048*q7a*s11**2*s16*s22**2*s23*s26*s33**2*s36**2 - 3072*q8a*s11**2*s12*s15*s22*s24*s26*s33**2*s36**2 + 1024*q8a*s11**2*s14*s16*s22*s24*s26*s33**2*s36**2 + 512*q7a*s11**2*s15*s22**2*s24*s26*s33**2*s36**2 -  \
    1024*q8a*s11**2*s12*s16*s24**2*s26*s33**2*s36**2 - 1024*q8a*s11**2*s12*s14*s22*s25*s26*s33**2*s36**2 - 2048*q8a*s11*s12*s15*s16*s22*s25*s26*s33**2*s36**2 - 512*q7a*s11**2*s14*s22**2*s25*s26*s33**2*s36**2 +  \
    1024*q7a*s11*s15*s16*s22**2*s25*s26*s33**2*s36**2 + 4096*q8a*s11**2*s12**2*s24*s25*s26*s33**2*s36**2 - 2048*q7a*s11**2*s12*s22*s24*s25*s26*s33**2*s36**2 - 1024*q8a*s11*s12**2*s16*s25**2*s26*s33**2*s36**2 +  \
    2048*q7a*s11*s12*s16*s22*s25**2*s26*s33**2*s36**2 + 3072*q8a*s11**2*s12*s13*s22*s26**2*s33**2*s36**2 - 1536*q8a*s11**2*s14**2*s22*s26**2*s33**2*s36**2 - 1024*q7a*s11**2*s13*s22**2*s26**2*s33**2*s36**2 -  \
    6144*q8a*s11**2*s12**2*s23*s26**2*s33**2*s36**2 + 3072*q7a*s11**2*s12*s22*s23*s26**2*s33**2*s36**2 + 1536*q8a*s11**2*s12*s14*s24*s26**2*s33**2*s36**2 + 512*q7a*s11**2*s14*s22*s24*s26**2*s33**2*s36**2 -  \
    512*q7a*s11**2*s12*s24**2*s26**2*s33**2*s36**2 + 1536*q8a*s11*s12**2*s15*s25*s26**2*s33**2*s36**2 - 1024*q7a*s11*s12*s15*s22*s25*s26**2*s33**2*s36**2 - 512*q7a*s11*s12**2*s25**2*s26**2*s33**2*s36**2 +  \
    2048*q9a*s11**3*s14*s22**2*s23**2*s34*s36**2 + 1024*q9a*s11**3*s13*s22**2*s23*s24*s34*s36**2 + 512*q9a*s11**2*s15**2*s22**2*s23*s24*s34*s36**2 - 3072*q9a*s11**3*s12*s22*s23**2*s24*s34*s36**2 +  \
    512*q9a*s11**2*s16**2*s22*s23**2*s24*s34*s36**2 - 1024*q9a*s11**3*s14*s22*s23*s24**2*s34*s36**2 - 512*q9a*s11**2*s15*s16*s22*s23*s24**2*s34*s36**2 + 1024*q9a*s11**3*s12*s23*s24**3*s34*s36**2 -  \
    2560*q9a*s11**2*s14*s15*s22**2*s23*s25*s34*s36**2 - 512*q9a*s11**2*s13*s16*s22**2*s23*s25*s34*s36**2 + 512*q9a*s11**2*s12*s16*s22*s23**2*s25*s34*s36**2 + 512*q9a*s11**2*s13*s15*s22**2*s24*s25*s34*s36**2 -  \
    512*q9a*s11*s15**3*s22**2*s24*s25*s34*s36**2 + 1024*q9a*s11**2*s12*s15*s22*s23*s24*s25*s34*s36**2 + 3072*q9a*s11**2*s14*s16*s22*s23*s24*s25*s34*s36**2 - 512*q9a*s11*s15*s16**2*s22*s23*s24*s25*s34*s36**2 +  \
    512*q9a*s11**2*s14*s15*s22*s24**2*s25*s34*s36**2 - 1024*q9a*s11**2*s13*s16*s22*s24**2*s25*s34*s36**2 + 512*q9a*s11*s15**2*s16*s22*s24**2*s25*s34*s36**2 - 1536*q9a*s11**2*s12*s16*s23*s24**2*s25*s34*s36**2 -  \
    512*q9a*s11**2*s12*s15*s24**3*s25*s34*s36**2 + 2560*q9a*s11**2*s13*s14*s22**2*s25**2*s34*s36**2 + 512*q9a*s11*s14*s15**2*s22**2*s25**2*s34*s36**2 - 512*q9a*s11*s13*s15*s16*s22**2*s25**2*s34*s36**2 +  \
    512*q9a*s11*s12*s15*s16*s22*s23*s25**2*s34*s36**2 - 1024*q9a*s11*s14*s16**2*s22*s23*s25**2*s34*s36**2 - 3072*q9a*s11**2*s12*s13*s22*s24*s25**2*s34*s36**2 - 512*q9a*s11**2*s14**2*s22*s24*s25**2*s34*s36**2 +  \
    1024*q9a*s11*s12*s15**2*s22*s24*s25**2*s34*s36**2 - 1536*q9a*s11*s14*s15*s16*s22*s24*s25**2*s34*s36**2 + 1024*q9a*s11*s13*s16**2*s22*s24*s25**2*s34*s36**2 + 1024*q9a*s11**2*s12**2*s23*s24*s25**2*s34*s36**2 +  \
    512*q9a*s11*s12*s16**2*s23*s24*s25**2*s34*s36**2 + 512*q9a*s11**2*s12*s14*s24**2*s25**2*s34*s36**2 + 512*q9a*s11*s12*s15*s16*s24**2*s25**2*s34*s36**2 - 1024*q9a*s11*s12*s14*s15*s22*s25**3*s34*s36**2 +  \
    512*q9a*s11*s12*s13*s16*s22*s25**3*s34*s36**2 + 1024*q9a*s11*s14**2*s16*s22*s25**3*s34*s36**2 - 512*q9a*s11*s12**2*s16*s23*s25**3*s34*s36**2 - 512*q9a*s11*s12**2*s15*s24*s25**3*s34*s36**2 -  \
    512*q9a*s11*s12*s14*s16*s24*s25**3*s34*s36**2 + 512*q9a*s11*s12**2*s14*s25**4*s34*s36**2 - 1536*q9a*s11**2*s13*s15*s22**2*s23*s26*s34*s36**2 + 1536*q9a*s11**2*s12*s15*s22*s23**2*s26*s34*s36**2 -  \
    2560*q9a*s11**2*s14*s16*s22*s23**2*s26*s34*s36**2 + 1024*q9a*s11**2*s13*s16*s22*s23*s24*s26*s34*s36**2 + 512*q9a*s11**2*s12*s16*s23**2*s24*s26*s34*s36**2 + 512*q9a*s11**2*s12*s15*s23*s24**2*s26*s34*s36**2 +  \
    1024*q9a*s11**2*s13**2*s22**2*s25*s26*s34*s36**2 + 512*q9a*s11*s13*s15**2*s22**2*s25*s26*s34*s36**2 + 512*q9a*s11**2*s14**2*s22*s23*s25*s26*s34*s36**2 - 512*q9a*s11*s12*s15**2*s22*s23*s25*s26*s34*s36**2 +  \
    1536*q9a*s11*s14*s15*s16*s22*s23*s25*s26*s34*s36**2 - 1024*q9a*s11**2*s12**2*s23**2*s25*s26*s34*s36**2 - 2048*q9a*s11**2*s13*s14*s22*s24*s25*s26*s34*s36**2 + 512*q9a*s11*s14*s15**2*s22*s24*s25*s26*s34*s36**2 -  \
    512*q9a*s11*s13*s15*s16*s22*s24*s25*s26*s34*s36**2 - 2048*q9a*s11**2*s12*s14*s23*s24*s25*s26*s34*s36**2 + 3072*q9a*s11**2*s12*s13*s24**2*s25*s26*s34*s36**2 - 1024*q9a*s11*s12*s15**2*s24**2*s25*s26*s34*s36**2 -  \
    512*q9a*s11*s12*s13*s15*s22*s25**2*s26*s34*s36**2 - 512*q9a*s11*s14**2*s15*s22*s25**2*s26*s34*s36**2 + 512*q9a*s11*s12**2*s15*s23*s25**2*s26*s34*s36**2 + 512*q9a*s11*s12*s14*s16*s23*s25**2*s26*s34*s36**2 +  \
    1536*q9a*s11*s12*s14*s15*s24*s25**2*s26*s34*s36**2 - 1536*q9a*s11*s12*s13*s16*s24*s25**2*s26*s34*s36**2 - 512*q9a*s11*s12*s14**2*s25**3*s26*s34*s36**2 + 1024*q9a*s11**2*s13*s14*s22*s23*s26**2*s34*s36**2 +  \
    1536*q9a*s11**2*s12*s14*s23**2*s26**2*s34*s36**2 - 2048*q9a*s11**2*s12*s13*s23*s24*s26**2*s34*s36**2 - 512*q9a*s11*s13*s14*s15*s22*s25*s26**2*s34*s36**2 - 1024*q9a*s11*s12*s14*s15*s23*s25*s26**2*s34*s36**2 +  \
    1024*q9a*s11*s12*s13*s15*s24*s25*s26**2*s34*s36**2 + 512*q9a*s11*s12*s13*s14*s25**2*s26**2*s34*s36**2 + 4096*q8a*s11**3*s14*s22**2*s23*s33*s34*s36**2 + 1024*q8a*s11**3*s13*s22**2*s24*s33*s34*s36**2 +  \
    512*q8a*s11**2*s15**2*s22**2*s24*s33*s34*s36**2 - 6144*q8a*s11**3*s12*s22*s23*s24*s33*s34*s36**2 + 1024*q8a*s11**2*s16**2*s22*s23*s24*s33*s34*s36**2 + 1024*q7a*s11**3*s22**2*s23*s24*s33*s34*s36**2 -  \
    1024*q8a*s11**3*s14*s22*s24**2*s33*s34*s36**2 - 512*q8a*s11**2*s15*s16*s22*s24**2*s33*s34*s36**2 + 1024*q8a*s11**3*s12*s24**3*s33*s34*s36**2 - 2560*q8a*s11**2*s14*s15*s22**2*s25*s33*s34*s36**2 -  \
    512*q8a*s11**2*s13*s16*s22**2*s25*s33*s34*s36**2 + 1024*q8a*s11**2*s12*s16*s22*s23*s25*s33*s34*s36**2 - 512*q7a*s11**2*s16*s22**2*s23*s25*s33*s34*s36**2 + 1024*q8a*s11**2*s12*s15*s22*s24*s25*s33*s34*s36**2 +  \
    3072*q8a*s11**2*s14*s16*s22*s24*s25*s33*s34*s36**2 - 512*q8a*s11*s15*s16**2*s22*s24*s25*s33*s34*s36**2 + 512*q7a*s11**2*s15*s22**2*s24*s25*s33*s34*s36**2 - 1536*q8a*s11**2*s12*s16*s24**2*s25*s33*s34*s36**2 -  \
    1024*q7a*s11**2*s16*s22*s24**2*s25*s33*s34*s36**2 + 512*q8a*s11*s12*s15*s16*s22*s25**2*s33*s34*s36**2 - 1024*q8a*s11*s14*s16**2*s22*s25**2*s33*s34*s36**2 + 2560*q7a*s11**2*s14*s22**2*s25**2*s33*s34*s36**2 -  \
    512*q7a*s11*s15*s16*s22**2*s25**2*s33*s34*s36**2 + 1024*q8a*s11**2*s12**2*s24*s25**2*s33*s34*s36**2 + 512*q8a*s11*s12*s16**2*s24*s25**2*s33*s34*s36**2 - 3072*q7a*s11**2*s12*s22*s24*s25**2*s33*s34*s36**2 +  \
    1024*q7a*s11*s16**2*s22*s24*s25**2*s33*s34*s36**2 - 512*q8a*s11*s12**2*s16*s25**3*s33*s34*s36**2 + 512*q7a*s11*s12*s16*s22*s25**3*s33*s34*s36**2 - 1536*q8a*s11**2*s13*s15*s22**2*s26*s33*s34*s36**2 +  \
    3072*q8a*s11**2*s12*s15*s22*s23*s26*s33*s34*s36**2 - 5120*q8a*s11**2*s14*s16*s22*s23*s26*s33*s34*s36**2 - 1536*q7a*s11**2*s15*s22**2*s23*s26*s33*s34*s36**2 + 1024*q8a*s11**2*s13*s16*s22*s24*s26*s33*s34*s36**2;
    v3_7= \
    1024*q8a*s11**2*s12*s16*s23*s24*s26*s33*s34*s36**2 + 1024*q7a*s11**2*s16*s22*s23*s24*s26*s33*s34*s36**2 + 512*q8a*s11**2*s12*s15*s24**2*s26*s33*s34*s36**2 + 512*q8a*s11**2*s14**2*s22*s25*s26*s33*s34*s36**2 -  \
    512*q8a*s11*s12*s15**2*s22*s25*s26*s33*s34*s36**2 + 1536*q8a*s11*s14*s15*s16*s22*s25*s26*s33*s34*s36**2 + 2048*q7a*s11**2*s13*s22**2*s25*s26*s33*s34*s36**2 + 512*q7a*s11*s15**2*s22**2*s25*s26*s33*s34*s36**2 -  \
    2048*q8a*s11**2*s12**2*s23*s25*s26*s33*s34*s36**2 - 2048*q8a*s11**2*s12*s14*s24*s25*s26*s33*s34*s36**2 - 2048*q7a*s11**2*s14*s22*s24*s25*s26*s33*s34*s36**2 - 512*q7a*s11*s15*s16*s22*s24*s25*s26*s33*s34*s36**2 +  \
    3072*q7a*s11**2*s12*s24**2*s25*s26*s33*s34*s36**2 + 512*q8a*s11*s12**2*s15*s25**2*s26*s33*s34*s36**2 + 512*q8a*s11*s12*s14*s16*s25**2*s26*s33*s34*s36**2 - 512*q7a*s11*s12*s15*s22*s25**2*s26*s33*s34*s36**2 -  \
    1536*q7a*s11*s12*s16*s24*s25**2*s26*s33*s34*s36**2 + 1024*q8a*s11**2*s13*s14*s22*s26**2*s33*s34*s36**2 + 3072*q8a*s11**2*s12*s14*s23*s26**2*s33*s34*s36**2 + 1024*q7a*s11**2*s14*s22*s23*s26**2*s33*s34*s36**2 -  \
    2048*q8a*s11**2*s12*s13*s24*s26**2*s33*s34*s36**2 - 2048*q7a*s11**2*s12*s23*s24*s26**2*s33*s34*s36**2 - 1024*q8a*s11*s12*s14*s15*s25*s26**2*s33*s34*s36**2 - 512*q7a*s11*s14*s15*s22*s25*s26**2*s33*s34*s36**2 +  \
    1024*q7a*s11*s12*s15*s24*s25*s26**2*s33*s34*s36**2 + 512*q7a*s11*s12*s14*s25**2*s26**2*s33*s34*s36**2 - 2048*q8a*s11**3*s13*s22**2*s23*s34**2*s36**2 - 1024*q8a*s11**2*s15**2*s22**2*s23*s34**2*s36**2 +  \
    3072*q8a*s11**3*s12*s22*s23**2*s34**2*s36**2 - 1536*q8a*s11**2*s16**2*s22*s23**2*s34**2*s36**2 - 1024*q7a*s11**3*s22**2*s23**2*s34**2*s36**2 + 1024*q8a*s11**3*s14*s22*s23*s24*s34**2*s36**2 +  \
    1024*q8a*s11**2*s15*s16*s22*s23*s24*s34**2*s36**2 - 1024*q8a*s11**3*s12*s23*s24**2*s34**2*s36**2 + 1024*q8a*s11**2*s13*s15*s22**2*s25*s34**2*s36**2 + 512*q8a*s11*s15**3*s22**2*s25*s34**2*s36**2 -  \
    2048*q8a*s11**2*s14*s16*s22*s23*s25*s34**2*s36**2 + 1024*q8a*s11*s15*s16**2*s22*s23*s25*s34**2*s36**2 + 1024*q7a*s11**2*s15*s22**2*s23*s25*s34**2*s36**2 - 512*q8a*s11**2*s14*s15*s22*s24*s25*s34**2*s36**2 -  \
    512*q8a*s11*s15**2*s16*s22*s24*s25*s34**2*s36**2 + 1024*q8a*s11**2*s12*s16*s23*s24*s25*s34**2*s36**2 + 512*q8a*s11**2*s12*s15*s24**2*s25*s34**2*s36**2 + 1024*q8a*s11**2*s12*s13*s22*s25**2*s34**2*s36**2 -  \
    1024*q8a*s11*s12*s15**2*s22*s25**2*s34**2*s36**2 + 1024*q8a*s11*s14*s15*s16*s22*s25**2*s34**2*s36**2 - 512*q8a*s11*s13*s16**2*s22*s25**2*s34**2*s36**2 - 2048*q7a*s11**2*s13*s22**2*s25**2*s34**2*s36**2 -  \
    512*q7a*s11*s15**2*s22**2*s25**2*s34**2*s36**2 - 1024*q8a*s11**2*s12**2*s23*s25**2*s34**2*s36**2 + 1024*q7a*s11**2*s12*s22*s23*s25**2*s34**2*s36**2 - 512*q7a*s11*s16**2*s22*s23*s25**2*s34**2*s36**2 -  \
    512*q8a*s11*s12*s15*s16*s24*s25**2*s34**2*s36**2 + 512*q7a*s11**2*s14*s22*s24*s25**2*s34**2*s36**2 + 512*q7a*s11*s15*s16*s22*s24*s25**2*s34**2*s36**2 - 512*q7a*s11**2*s12*s24**2*s25**2*s34**2*s36**2 +  \
    512*q8a*s11*s12**2*s15*s25**3*s34**2*s36**2 + 1024*q7a*s11*s12*s15*s22*s25**3*s34**2*s36**2 - 1024*q7a*s11*s14*s16*s22*s25**3*s34**2*s36**2 + 512*q7a*s11*s12*s16*s24*s25**3*s34**2*s36**2 - 512*q7a*s11*s12**2*s25**4*s34**2*s36**2 +  \
    1024*q8a*s11**2*s14*s15*s22*s23*s26*s34**2*s36**2 + 1024*q8a*s11**2*s13*s16*s22*s23*s26*s34**2*s36**2 + 1536*q8a*s11**2*s12*s16*s23**2*s26*s34**2*s36**2 + 512*q7a*s11**2*s16*s22*s23**2*s26*s34**2*s36**2 -  \
    2048*q8a*s11**2*s12*s15*s23*s24*s26*s34**2*s36**2 - 512*q8a*s11*s14*s15**2*s22*s25*s26*s34**2*s36**2 - 512*q8a*s11*s13*s15*s16*s22*s25*s26*s34**2*s36**2 + 1024*q8a*s11**2*s12*s14*s23*s25*s26*s34**2*s36**2 -  \
    1024*q8a*s11*s12*s15*s16*s23*s25*s26*s34**2*s36**2 - 512*q7a*s11*s15*s16*s22*s23*s25*s26*s34**2*s36**2 + 1024*q8a*s11*s12*s15**2*s24*s25*s26*s34**2*s36**2 - 512*q8a*s11*s12*s14*s15*s25**2*s26*s34**2*s36**2 +  \
    512*q8a*s11*s12*s13*s16*s25**2*s26*s34**2*s36**2 + 512*q7a*s11*s14*s15*s22*s25**2*s26*s34**2*s36**2 + 1024*q7a*s11*s13*s16*s22*s25**2*s26*s34**2*s36**2 + 512*q7a*s11*s12*s16*s23*s25**2*s26*s34**2*s36**2 -  \
    1024*q7a*s11*s12*s15*s24*s25**2*s26*s34**2*s36**2 + 512*q7a*s11*s12*s14*s25**3*s26*s34**2*s36**2 - 1024*q8a*s11**2*s12*s13*s23*s26**2*s34**2*s36**2 - 512*q7a*s11**2*s12*s23**2*s26**2*s34**2*s36**2 +  \
    512*q8a*s11*s12*s13*s15*s25*s26**2*s34**2*s36**2 + 512*q7a*s11*s12*s15*s23*s25*s26**2*s34**2*s36**2 - 1024*q7a*s11*s12*s13*s25**2*s26**2*s34**2*s36**2 + 1024*q9a*s11**2*s13*s15*s22**3*s23*s35*s36**2 -  \
    1024*q9a*s11**2*s12*s15*s22**2*s23**2*s35*s36**2 + 1024*q9a*s11**2*s14*s16*s22**2*s23**2*s35*s36**2 + 512*q9a*s11**2*s14*s15*s22**2*s23*s24*s35*s36**2 - 1536*q9a*s11**2*s13*s16*s22**2*s23*s24*s35*s36**2 +  \
    512*q9a*s11**2*s12*s16*s22*s23**2*s24*s35*s36**2 + 512*q9a*s11**2*s13*s15*s22**2*s24**2*s35*s36**2 - 1024*q9a*s11**2*s12*s15*s22*s23*s24**2*s35*s36**2 + 512*q9a*s11**2*s14*s16*s22*s23*s24**2*s35*s36**2 -  \
    512*q9a*s11**2*s14*s15*s22*s24**3*s35*s36**2 - 512*q9a*s11**2*s12*s16*s23*s24**3*s35*s36**2 + 512*q9a*s11**2*s12*s15*s24**4*s35*s36**2 + 2048*q9a*s11**2*s13**2*s22**3*s25*s35*s36**2 - 1024*q9a*s11*s13*s15**2*s22**3*s25*s35*s36**2 -  \
    5120*q9a*s11**2*s12*s13*s22**2*s23*s25*s35*s36**2 + 2560*q9a*s11**2*s14**2*s22**2*s23*s25*s35*s36**2 + 1024*q9a*s11*s12*s15**2*s22**2*s23*s25*s35*s36**2 - 2048*q9a*s11*s14*s15*s16*s22**2*s23*s25*s35*s36**2 +  \
    1024*q9a*s11*s13*s16**2*s22**2*s23*s25*s35*s36**2 + 3072*q9a*s11**2*s12**2*s22*s23**2*s25*s35*s36**2 - 1024*q9a*s11*s12*s16**2*s22*s23**2*s25*s35*s36**2 - 2560*q9a*s11**2*s13*s14*s22**2*s24*s25*s35*s36**2 +  \
    512*q9a*s11*s14*s15**2*s22**2*s24*s25*s35*s36**2 - 3072*q9a*s11**2*s12*s14*s22*s23*s24*s25*s35*s36**2 + 2048*q9a*s11*s12*s15*s16*s22*s23*s24*s25*s35*s36**2 - 512*q9a*s11*s14*s16**2*s22*s23*s24*s25*s35*s36**2 +  \
    2048*q9a*s11**2*s12*s13*s22*s24**2*s25*s35*s36**2 + 512*q9a*s11**2*s14**2*s22*s24**2*s25*s35*s36**2 - 512*q9a*s11*s12*s15**2*s22*s24**2*s25*s35*s36**2 + 512*q9a*s11*s14*s15*s16*s22*s24**2*s25*s35*s36**2 +  \
    1024*q9a*s11**2*s12**2*s23*s24**2*s25*s35*s36**2 + 512*q9a*s11*s12*s16**2*s23*s24**2*s25*s35*s36**2 - 512*q9a*s11**2*s12*s14*s24**3*s25*s35*s36**2 - 512*q9a*s11*s12*s15*s16*s24**3*s25*s35*s36**2 +  \
    2048*q9a*s11*s12*s13*s15*s22**2*s25**2*s35*s36**2 - 512*q9a*s11*s14**2*s15*s22**2*s25**2*s35*s36**2 + 512*q9a*s11*s13*s14*s16*s22**2*s25**2*s35*s36**2 - 2048*q9a*s11*s12**2*s15*s22*s23*s25**2*s35*s36**2 +  \
    1536*q9a*s11*s12*s14*s16*s22*s23*s25**2*s35*s36**2 - 512*q9a*s11*s12*s13*s16*s22*s24*s25**2*s35*s36**2 - 512*q9a*s11*s14**2*s16*s22*s24*s25**2*s35*s36**2 - 1536*q9a*s11*s12**2*s16*s23*s24*s25**2*s35*s36**2 +  \
    512*q9a*s11*s12**2*s15*s24**2*s25**2*s35*s36**2 + 512*q9a*s11*s12*s14*s16*s24**2*s25**2*s35*s36**2 - 1024*q9a*s11*s12**2*s13*s22*s25**3*s35*s36**2 + 512*q9a*s11*s12*s14**2*s22*s25**3*s35*s36**2 +  \
    1024*q9a*s11*s12**3*s23*s25**3*s35*s36**2 - 512*q9a*s11*s12**2*s14*s24*s25**3*s35*s36**2 - 512*q9a*s11**2*s13*s14*s22**2*s23*s26*s35*s36**2 - 1024*q9a*s11*s14*s15**2*s22**2*s23*s26*s35*s36**2 +  \
    1024*q9a*s11*s13*s15*s16*s22**2*s23*s26*s35*s36**2 - 512*q9a*s11**2*s12*s14*s22*s23**2*s26*s35*s36**2 - 1024*q9a*s11*s12*s15*s16*s22*s23**2*s26*s35*s36**2 - 512*q9a*s11*s13*s15**2*s22**2*s24*s26*s35*s36**2 +  \
    2048*q9a*s11**2*s12*s13*s22*s23*s24*s26*s35*s36**2 - 512*q9a*s11**2*s14**2*s22*s23*s24*s26*s35*s36**2 + 1536*q9a*s11*s12*s15**2*s22*s23*s24*s26*s35*s36**2 - 512*q9a*s11*s14*s15*s16*s22*s23*s24*s26*s35*s36**2 -  \
    1024*q9a*s11**2*s12**2*s23**2*s24*s26*s35*s36**2 + 512*q9a*s11*s14*s15**2*s22*s24**2*s26*s35*s36**2 + 512*q9a*s11**2*s12*s14*s23*s24**2*s26*s35*s36**2 + 512*q9a*s11*s12*s15*s16*s23*s24**2*s26*s35*s36**2 -  \
    512*q9a*s11*s12*s15**2*s24**3*s26*s35*s36**2 + 3072*q9a*s11*s13*s14*s15*s22**2*s25*s26*s35*s36**2 - 2560*q9a*s11*s13**2*s16*s22**2*s25*s26*s35*s36**2 + 1024*q9a*s11*s12*s14*s15*s22*s23*s25*s26*s35*s36**2 +  \
    2048*q9a*s11*s12*s13*s16*s22*s23*s25*s26*s35*s36**2 + 512*q9a*s11*s12**2*s16*s23**2*s25*s26*s35*s36**2 - 2048*q9a*s11*s12*s13*s15*s22*s24*s25*s26*s35*s36**2 - 1536*q9a*s11*s14**2*s15*s22*s24*s25*s26*s35*s36**2 +  \
    1536*q9a*s11*s13*s14*s16*s22*s24*s25*s26*s35*s36**2 - 2048*q9a*s11*s12**2*s15*s23*s24*s25*s26*s35*s36**2 + 1536*q9a*s11*s12*s14*s15*s24**2*s25*s26*s35*s36**2 - 1536*q9a*s11*s12*s13*s16*s24**2*s25*s26*s35*s36**2 -  \
    3584*q9a*s11*s12*s13*s14*s22*s25**2*s26*s35*s36**2 + 1024*q9a*s11*s14**3*s22*s25**2*s26*s35*s36**2 + 512*q9a*s11*s12**2*s14*s23*s25**2*s26*s35*s36**2 + 3072*q9a*s11*s12**2*s13*s24*s25**2*s26*s35*s36**2 -  \
    1024*q9a*s11*s12*s14**2*s24*s25**2*s26*s35*s36**2 + 512*q9a*s11*s13**2*s15*s22**2*s26**2*s35*s36**2 - 2048*q9a*s11*s12*s13*s15*s22*s23*s26**2*s35*s36**2 + 1024*q9a*s11*s14**2*s15*s22*s23*s26**2*s35*s36**2 +  \
    1536*q9a*s11*s12**2*s15*s23**2*s26**2*s35*s36**2 - 512*q9a*s11*s13*s14*s15*s22*s24*s26**2*s35*s36**2 - 1024*q9a*s11*s12*s14*s15*s23*s24*s26**2*s35*s36**2 + 512*q9a*s11*s12*s13*s15*s24**2*s26**2*s35*s36**2 +  \
    2048*q9a*s11*s12*s13**2*s22*s25*s26**2*s35*s36**2 - 1024*q9a*s11*s13*s14**2*s22*s25*s26**2*s35*s36**2 - 2048*q9a*s11*s12**2*s13*s23*s25*s26**2*s35*s36**2 + 1024*q9a*s11*s12*s13*s14*s24*s25*s26**2*s35*s36**2 +  \
    1024*q8a*s11**2*s13*s15*s22**3*s33*s35*s36**2 - 2048*q8a*s11**2*s12*s15*s22**2*s23*s33*s35*s36**2 + 2048*q8a*s11**2*s14*s16*s22**2*s23*s33*s35*s36**2 + 1024*q7a*s11**2*s15*s22**3*s23*s33*s35*s36**2 +  \
    512*q8a*s11**2*s14*s15*s22**2*s24*s33*s35*s36**2 - 1536*q8a*s11**2*s13*s16*s22**2*s24*s33*s35*s36**2 + 1024*q8a*s11**2*s12*s16*s22*s23*s24*s33*s35*s36**2 - 1536*q7a*s11**2*s16*s22**2*s23*s24*s33*s35*s36**2 -  \
    1024*q8a*s11**2*s12*s15*s22*s24**2*s33*s35*s36**2 + 512*q8a*s11**2*s14*s16*s22*s24**2*s33*s35*s36**2 + 512*q7a*s11**2*s15*s22**2*s24**2*s33*s35*s36**2 - 512*q8a*s11**2*s12*s16*s24**3*s33*s35*s36**2 -  \
    5120*q8a*s11**2*s12*s13*s22**2*s25*s33*s35*s36**2 + 2560*q8a*s11**2*s14**2*s22**2*s25*s33*s35*s36**2 + 1024*q8a*s11*s12*s15**2*s22**2*s25*s33*s35*s36**2 - 2048*q8a*s11*s14*s15*s16*s22**2*s25*s33*s35*s36**2 +  \
    1024*q8a*s11*s13*s16**2*s22**2*s25*s33*s35*s36**2 + 4096*q7a*s11**2*s13*s22**3*s25*s33*s35*s36**2 - 1024*q7a*s11*s15**2*s22**3*s25*s33*s35*s36**2 + 6144*q8a*s11**2*s12**2*s22*s23*s25*s33*s35*s36**2 -  \
    2048*q8a*s11*s12*s16**2*s22*s23*s25*s33*s35*s36**2 - 5120*q7a*s11**2*s12*s22**2*s23*s25*s33*s35*s36**2 + 1024*q7a*s11*s16**2*s22**2*s23*s25*s33*s35*s36**2 - 3072*q8a*s11**2*s12*s14*s22*s24*s25*s33*s35*s36**2 +  \
    2048*q8a*s11*s12*s15*s16*s22*s24*s25*s33*s35*s36**2 - 512*q8a*s11*s14*s16**2*s22*s24*s25*s33*s35*s36**2 - 2560*q7a*s11**2*s14*s22**2*s24*s25*s33*s35*s36**2 + 1024*q8a*s11**2*s12**2*s24**2*s25*s33*s35*s36**2 +  \
    512*q8a*s11*s12*s16**2*s24**2*s25*s33*s35*s36**2 + 2048*q7a*s11**2*s12*s22*s24**2*s25*s33*s35*s36**2 - 2048*q8a*s11*s12**2*s15*s22*s25**2*s33*s35*s36**2 + 1536*q8a*s11*s12*s14*s16*s22*s25**2*s33*s35*s36**2 +  \
    2048*q7a*s11*s12*s15*s22**2*s25**2*s33*s35*s36**2 + 512*q7a*s11*s14*s16*s22**2*s25**2*s33*s35*s36**2 - 1536*q8a*s11*s12**2*s16*s24*s25**2*s33*s35*s36**2 - 512*q7a*s11*s12*s16*s22*s24*s25**2*s33*s35*s36**2 +  \
    1024*q8a*s11*s12**3*s25**3*s33*s35*s36**2 - 1024*q7a*s11*s12**2*s22*s25**3*s33*s35*s36**2 - 512*q8a*s11**2*s13*s14*s22**2*s26*s33*s35*s36**2 - 1024*q8a*s11*s14*s15**2*s22**2*s26*s33*s35*s36**2 +  \
    1024*q8a*s11*s13*s15*s16*s22**2*s26*s33*s35*s36**2 - 1024*q8a*s11**2*s12*s14*s22*s23*s26*s33*s35*s36**2 - 2048*q8a*s11*s12*s15*s16*s22*s23*s26*s33*s35*s36**2 - 512*q7a*s11**2*s14*s22**2*s23*s26*s33*s35*s36**2 +  \
    1024*q7a*s11*s15*s16*s22**2*s23*s26*s33*s35*s36**2 + 2048*q8a*s11**2*s12*s13*s22*s24*s26*s33*s35*s36**2 - 512*q8a*s11**2*s14**2*s22*s24*s26*s33*s35*s36**2 + 1536*q8a*s11*s12*s15**2*s22*s24*s26*s33*s35*s36**2 -  \
    512*q8a*s11*s14*s15*s16*s22*s24*s26*s33*s35*s36**2 - 512*q7a*s11*s15**2*s22**2*s24*s26*s33*s35*s36**2 - 2048*q8a*s11**2*s12**2*s23*s24*s26*s33*s35*s36**2 + 2048*q7a*s11**2*s12*s22*s23*s24*s26*s33*s35*s36**2 +  \
    512*q8a*s11**2*s12*s14*s24**2*s26*s33*s35*s36**2 + 512*q8a*s11*s12*s15*s16*s24**2*s26*s33*s35*s36**2 + 1024*q8a*s11*s12*s14*s15*s22*s25*s26*s33*s35*s36**2 + 2048*q8a*s11*s12*s13*s16*s22*s25*s26*s33*s35*s36**2 +  \
    3072*q7a*s11*s14*s15*s22**2*s25*s26*s33*s35*s36**2 - 5120*q7a*s11*s13*s16*s22**2*s25*s26*s33*s35*s36**2 + 1024*q8a*s11*s12**2*s16*s23*s25*s26*s33*s35*s36**2 + 2048*q7a*s11*s12*s16*s22*s23*s25*s26*s33*s35*s36**2 -  \
    2048*q8a*s11*s12**2*s15*s24*s25*s26*s33*s35*s36**2 - 2048*q7a*s11*s12*s15*s22*s24*s25*s26*s33*s35*s36**2 + 1536*q7a*s11*s14*s16*s22*s24*s25*s26*s33*s35*s36**2 - 1536*q7a*s11*s12*s16*s24**2*s25*s26*s33*s35*s36**2 +  \
    512*q8a*s11*s12**2*s14*s25**2*s26*s33*s35*s36**2 - 3584*q7a*s11*s12*s14*s22*s25**2*s26*s33*s35*s36**2 + 3072*q7a*s11*s12**2*s24*s25**2*s26*s33*s35*s36**2 - 2048*q8a*s11*s12*s13*s15*s22*s26**2*s33*s35*s36**2;
    v3_8= \
    1024*q8a*s11*s14**2*s15*s22*s26**2*s33*s35*s36**2 + 1024*q7a*s11*s13*s15*s22**2*s26**2*s33*s35*s36**2 + 3072*q8a*s11*s12**2*s15*s23*s26**2*s33*s35*s36**2 - 2048*q7a*s11*s12*s15*s22*s23*s26**2*s33*s35*s36**2 -  \
    1024*q8a*s11*s12*s14*s15*s24*s26**2*s33*s35*s36**2 - 512*q7a*s11*s14*s15*s22*s24*s26**2*s33*s35*s36**2 + 512*q7a*s11*s12*s15*s24**2*s26**2*s33*s35*s36**2 - 2048*q8a*s11*s12**2*s13*s25*s26**2*s33*s35*s36**2 +  \
    4096*q7a*s11*s12*s13*s22*s25*s26**2*s33*s35*s36**2 - 1024*q7a*s11*s14**2*s22*s25*s26**2*s33*s35*s36**2 - 2048*q7a*s11*s12**2*s23*s25*s26**2*s33*s35*s36**2 + 1024*q7a*s11*s12*s14*s24*s25*s26**2*s33*s35*s36**2 +  \
    2048*q8a*s11**2*s14*s15*s22**2*s23*s34*s35*s36**2 - 2048*q8a*s11**2*s13*s16*s22**2*s23*s34*s35*s36**2 + 3072*q8a*s11**2*s12*s16*s22*s23**2*s34*s35*s36**2 - 1024*q7a*s11**2*s16*s22**2*s23**2*s34*s35*s36**2 -  \
    1536*q8a*s11**2*s13*s15*s22**2*s24*s34*s35*s36**2 + 1024*q8a*s11**2*s12*s15*s22*s23*s24*s34*s35*s36**2 - 3072*q8a*s11**2*s14*s16*s22*s23*s24*s34*s35*s36**2 - 1536*q7a*s11**2*s15*s22**2*s23*s24*s34*s35*s36**2 +  \
    512*q8a*s11**2*s14*s15*s22*s24**2*s34*s35*s36**2 + 1024*q8a*s11**2*s13*s16*s22*s24**2*s34*s35*s36**2 + 1024*q8a*s11**2*s12*s16*s23*s24**2*s34*s35*s36**2 + 1024*q7a*s11**2*s16*s22*s23*s24**2*s34*s35*s36**2 -  \
    512*q8a*s11**2*s12*s15*s24**3*s34*s35*s36**2 - 512*q8a*s11**2*s13*s14*s22**2*s25*s34*s35*s36**2 - 1024*q8a*s11*s14*s15**2*s22**2*s25*s34*s35*s36**2 + 2048*q8a*s11*s13*s15*s16*s22**2*s25*s34*s35*s36**2 -  \
    1024*q8a*s11**2*s12*s14*s22*s23*s25*s34*s35*s36**2 - 4096*q8a*s11*s12*s15*s16*s22*s23*s25*s34*s35*s36**2 + 2048*q8a*s11*s14*s16**2*s22*s23*s25*s34*s35*s36**2 - 512*q7a*s11**2*s14*s22**2*s23*s25*s34*s35*s36**2 +  \
    2048*q7a*s11*s15*s16*s22**2*s23*s25*s34*s35*s36**2 + 512*q8a*s11**2*s14**2*s22*s24*s25*s34*s35*s36**2 + 512*q8a*s11*s12*s15**2*s22*s24*s25*s34*s35*s36**2 - 512*q8a*s11*s13*s16**2*s22*s24*s25*s34*s35*s36**2 +  \
    2048*q7a*s11**2*s13*s22**2*s24*s25*s34*s35*s36**2 + 512*q7a*s11*s15**2*s22**2*s24*s25*s34*s35*s36**2 - 1024*q8a*s11*s12*s16**2*s23*s24*s25*s34*s35*s36**2 - 512*q7a*s11*s16**2*s22*s23*s24*s25*s34*s35*s36**2 -  \
    512*q8a*s11**2*s12*s14*s24**2*s25*s34*s35*s36**2 + 512*q8a*s11*s12*s15*s16*s24**2*s25*s34*s35*s36**2 - 1024*q7a*s11**2*s14*s22*s24**2*s25*s34*s35*s36**2 - 512*q7a*s11*s15*s16*s22*s24**2*s25*s34*s35*s36**2 +  \
    1024*q7a*s11**2*s12*s24**3*s25*s34*s35*s36**2 + 1536*q8a*s11*s12*s14*s15*s22*s25**2*s34*s35*s36**2 + 1024*q8a*s11*s12*s13*s16*s22*s25**2*s34*s35*s36**2 - 1024*q8a*s11*s14**2*s16*s22*s25**2*s34*s35*s36**2 +  \
    512*q7a*s11*s14*s15*s22**2*s25**2*s34*s35*s36**2 - 3072*q7a*s11*s13*s16*s22**2*s25**2*s34*s35*s36**2 + 1024*q8a*s11*s12**2*s16*s23*s25**2*s34*s35*s36**2 + 1024*q7a*s11*s12*s16*s22*s23*s25**2*s34*s35*s36**2 -  \
    512*q8a*s11*s12**2*s15*s24*s25**2*s34*s35*s36**2 + 512*q8a*s11*s12*s14*s16*s24*s25**2*s34*s35*s36**2 - 1536*q7a*s11*s12*s15*s22*s24*s25**2*s34*s35*s36**2 + 1536*q7a*s11*s14*s16*s22*s24*s25**2*s34*s35*s36**2 -  \
    1024*q7a*s11*s12*s16*s24**2*s25**2*s34*s35*s36**2 - 512*q8a*s11*s12**2*s14*s25**3*s34*s35*s36**2 - 512*q7a*s11*s12*s14*s22*s25**3*s34*s35*s36**2 + 1024*q7a*s11*s12**2*s24*s25**3*s34*s35*s36**2 -  \
    1024*q8a*s11**2*s13**2*s22**2*s26*s34*s35*s36**2 + 1024*q8a*s11*s13*s15**2*s22**2*s26*s34*s35*s36**2 + 6144*q8a*s11**2*s12*s13*s22*s23*s26*s34*s35*s36**2 - 3072*q8a*s11**2*s14**2*s22*s23*s26*s34*s35*s36**2 -  \
    2048*q8a*s11*s12*s15**2*s22*s23*s26*s34*s35*s36**2 + 2048*q8a*s11*s14*s15*s16*s22*s23*s26*s34*s35*s36**2 - 2048*q7a*s11**2*s13*s22**2*s23*s26*s34*s35*s36**2 + 1024*q7a*s11*s15**2*s22**2*s23*s26*s34*s35*s36**2 -  \
    6144*q8a*s11**2*s12**2*s23**2*s26*s34*s35*s36**2 + 3072*q7a*s11**2*s12*s22*s23**2*s26*s34*s35*s36**2 + 2048*q8a*s11**2*s13*s14*s22*s24*s26*s34*s35*s36**2 - 512*q8a*s11*s14*s15**2*s22*s24*s26*s34*s35*s36**2 -  \
    512*q8a*s11*s13*s15*s16*s22*s24*s26*s34*s35*s36**2 + 5120*q8a*s11**2*s12*s14*s23*s24*s26*s34*s35*s36**2 - 1024*q8a*s11*s12*s15*s16*s23*s24*s26*s34*s35*s36**2 + 2048*q7a*s11**2*s14*s22*s23*s24*s26*s34*s35*s36**2 -  \
    512*q7a*s11*s15*s16*s22*s23*s24*s26*s34*s35*s36**2 - 3072*q8a*s11**2*s12*s13*s24**2*s26*s34*s35*s36**2 + 512*q8a*s11*s12*s15**2*s24**2*s26*s34*s35*s36**2 - 3072*q7a*s11**2*s12*s23*s24**2*s26*s34*s35*s36**2 -  \
    1024*q8a*s11*s12*s13*s15*s22*s25*s26*s34*s35*s36**2 + 1536*q8a*s11*s14**2*s15*s22*s25*s26*s34*s35*s36**2 - 2560*q8a*s11*s13*s14*s16*s22*s25*s26*s34*s35*s36**2 - 3072*q7a*s11*s13*s15*s22**2*s25*s26*s34*s35*s36**2 +  \
    5120*q8a*s11*s12**2*s15*s23*s25*s26*s34*s35*s36**2 - 1024*q8a*s11*s12*s14*s16*s23*s25*s26*s34*s35*s36**2 - 1024*q7a*s11*s12*s15*s22*s23*s25*s26*s34*s35*s36**2 - 2560*q7a*s11*s14*s16*s22*s23*s25*s26*s34*s35*s36**2 -  \
    2048*q8a*s11*s12*s14*s15*s24*s25*s26*s34*s35*s36**2 + 2048*q8a*s11*s12*s13*s16*s24*s25*s26*s34*s35*s36**2 + 2048*q7a*s11*s13*s16*s22*s24*s25*s26*s34*s35*s36**2 + 2048*q7a*s11*s12*s16*s23*s24*s25*s26*s34*s35*s36**2 +  \
    512*q7a*s11*s12*s15*s24**2*s25*s26*s34*s35*s36**2 - 3072*q8a*s11*s12**2*s13*s25**2*s26*s34*s35*s36**2 + 512*q8a*s11*s12*s14**2*s25**2*s26*s34*s35*s36**2 + 6144*q7a*s11*s12*s13*s22*s25**2*s26*s34*s35*s36**2 -  \
    1024*q7a*s11*s14**2*s22*s25**2*s26*s34*s35*s36**2 - 3072*q7a*s11*s12**2*s23*s25**2*s26*s34*s35*s36**2 + 512*q7a*s11*s12*s14*s24*s25**2*s26*s34*s35*s36**2 - 512*q8a*s11*s13*s14*s15*s22*s26**2*s34*s35*s36**2 -  \
    1024*q8a*s11*s12*s14*s15*s23*s26**2*s34*s35*s36**2 - 512*q7a*s11*s14*s15*s22*s23*s26**2*s34*s35*s36**2 + 1024*q8a*s11*s12*s13*s15*s24*s26**2*s34*s35*s36**2 + 1024*q7a*s11*s12*s15*s23*s24*s26**2*s34*s35*s36**2 +  \
    1024*q8a*s11*s12*s13*s14*s25*s26**2*s34*s35*s36**2 + 2048*q7a*s11*s13*s14*s22*s25*s26**2*s34*s35*s36**2 + 1024*q7a*s11*s12*s14*s23*s25*s26**2*s34*s35*s36**2 - 4096*q7a*s11*s12*s13*s24*s25*s26**2*s34*s35*s36**2 -  \
    1024*q8a*s11**2*s13**2*s22**3*s35**2*s36**2 + 4096*q8a*s11**2*s12*s13*s22**2*s23*s35**2*s36**2 - 2048*q8a*s11**2*s14**2*s22**2*s23*s35**2*s36**2 - 2048*q7a*s11**2*s13*s22**3*s23*s35**2*s36**2 -  \
    3072*q8a*s11**2*s12**2*s22*s23**2*s35**2*s36**2 + 2048*q7a*s11**2*s12*s22**2*s23**2*s35**2*s36**2 + 1024*q8a*s11**2*s13*s14*s22**2*s24*s35**2*s36**2 + 2048*q8a*s11**2*s12*s14*s22*s23*s24*s35**2*s36**2 +  \
    1024*q7a*s11**2*s14*s22**2*s23*s24*s35**2*s36**2 - 512*q8a*s11**2*s14**2*s22*s24**2*s35**2*s36**2 - 1024*q7a*s11**2*s13*s22**2*s24**2*s35**2*s36**2 - 1024*q8a*s11**2*s12**2*s23*s24**2*s35**2*s36**2 +  \
    512*q8a*s11**2*s12*s14*s24**3*s35**2*s36**2 + 512*q7a*s11**2*s14*s22*s24**3*s35**2*s36**2 - 512*q7a*s11**2*s12*s24**4*s35**2*s36**2 - 1024*q8a*s11*s12*s13*s15*s22**2*s25*s35**2*s36**2 + 512*q8a*s11*s14**2*s15*s22**2*s25*s35**2*s36**2 +  \
    1024*q7a*s11*s13*s15*s22**3*s25*s35**2*s36**2 + 1024*q8a*s11*s12**2*s15*s22*s23*s25*s35**2*s36**2 - 1024*q7a*s11*s12*s15*s22**2*s23*s25*s35**2*s36**2 - 512*q8a*s11*s12*s14*s15*s22*s24*s25*s35**2*s36**2 -  \
    1024*q8a*s11*s12*s13*s16*s22*s24*s25*s35**2*s36**2 + 512*q8a*s11*s14**2*s16*s22*s24*s25*s35**2*s36**2 - 512*q7a*s11*s14*s15*s22**2*s24*s25*s35**2*s36**2 + 1024*q7a*s11*s13*s16*s22**2*s24*s25*s35**2*s36**2 +  \
    1024*q8a*s11*s12**2*s16*s23*s24*s25*s35**2*s36**2 - 1024*q7a*s11*s12*s16*s22*s23*s24*s25*s35**2*s36**2 - 512*q8a*s11*s12*s14*s16*s24**2*s25*s35**2*s36**2 + 512*q7a*s11*s12*s15*s22*s24**2*s25*s35**2*s36**2 -  \
    512*q7a*s11*s14*s16*s22*s24**2*s25*s35**2*s36**2 + 512*q7a*s11*s12*s16*s24**3*s25*s35**2*s36**2 + 1024*q8a*s11*s12**2*s13*s22*s25**2*s35**2*s36**2 - 512*q8a*s11*s12*s14**2*s22*s25**2*s35**2*s36**2 -  \
    1024*q7a*s11*s12*s13*s22**2*s25**2*s35**2*s36**2 - 1024*q8a*s11*s12**3*s23*s25**2*s35**2*s36**2 + 1024*q7a*s11*s12**2*s22*s23*s25**2*s35**2*s36**2 + 512*q8a*s11*s12**2*s14*s24*s25**2*s35**2*s36**2 +  \
    512*q7a*s11*s12*s14*s22*s24*s25**2*s35**2*s36**2 - 512*q7a*s11*s12**2*s24**2*s25**2*s35**2*s36**2 + 512*q8a*s11*s13**2*s16*s22**2*s26*s35**2*s36**2 - 2048*q8a*s11*s12*s13*s16*s22*s23*s26*s35**2*s36**2 +  \
    1024*q8a*s11*s14**2*s16*s22*s23*s26*s35**2*s36**2 + 1024*q7a*s11*s13*s16*s22**2*s23*s26*s35**2*s36**2 + 1536*q8a*s11*s12**2*s16*s23**2*s26*s35**2*s36**2 - 1024*q7a*s11*s12*s16*s22*s23**2*s26*s35**2*s36**2 -  \
    1024*q8a*s11*s12*s13*s15*s22*s24*s26*s35**2*s36**2 + 512*q8a*s11*s14**2*s15*s22*s24*s26*s35**2*s36**2 - 512*q8a*s11*s13*s14*s16*s22*s24*s26*s35**2*s36**2 + 1024*q7a*s11*s13*s15*s22**2*s24*s26*s35**2*s36**2 +  \
    1024*q8a*s11*s12**2*s15*s23*s24*s26*s35**2*s36**2 - 1024*q8a*s11*s12*s14*s16*s23*s24*s26*s35**2*s36**2 - 1024*q7a*s11*s12*s15*s22*s23*s24*s26*s35**2*s36**2 - 512*q7a*s11*s14*s16*s22*s23*s24*s26*s35**2*s36**2 -  \
    512*q8a*s11*s12*s14*s15*s24**2*s26*s35**2*s36**2 + 512*q8a*s11*s12*s13*s16*s24**2*s26*s35**2*s36**2 - 512*q7a*s11*s14*s15*s22*s24**2*s26*s35**2*s36**2 + 512*q7a*s11*s12*s16*s23*s24**2*s26*s35**2*s36**2 +  \
    512*q7a*s11*s12*s15*s24**3*s26*s35**2*s36**2 + 2048*q8a*s11*s12*s13*s14*s22*s25*s26*s35**2*s36**2 - 1024*q8a*s11*s14**3*s22*s25*s26*s35**2*s36**2 - 2048*q7a*s11*s13*s14*s22**2*s25*s26*s35**2*s36**2 -  \
    2048*q8a*s11*s12**2*s14*s23*s25*s26*s35**2*s36**2 + 2048*q7a*s11*s12*s14*s22*s23*s25*s26*s35**2*s36**2 + 1024*q8a*s11*s12*s14**2*s24*s25*s26*s35**2*s36**2 + 1024*q7a*s11*s14**2*s22*s24*s25*s26*s35**2*s36**2 -  \
    1024*q7a*s11*s12*s14*s24**2*s25*s26*s35**2*s36**2 + 1024*q8a*s11*s12*s13**2*s22*s26**2*s35**2*s36**2 - 512*q8a*s11*s13*s14**2*s22*s26**2*s35**2*s36**2 - 1536*q7a*s11*s13**2*s22**2*s26**2*s35**2*s36**2 -  \
    1024*q8a*s11*s12**2*s13*s23*s26**2*s35**2*s36**2 + 2048*q7a*s11*s12*s13*s22*s23*s26**2*s35**2*s36**2 - 512*q7a*s11*s14**2*s22*s23*s26**2*s35**2*s36**2 - 512*q7a*s11*s12**2*s23**2*s26**2*s35**2*s36**2 +  \
    512*q8a*s11*s12*s13*s14*s24*s26**2*s35**2*s36**2 + 1024*q7a*s11*s13*s14*s22*s24*s26**2*s35**2*s36**2 + 512*q7a*s11*s12*s14*s23*s24*s26**2*s35**2*s36**2 - 1024*q7a*s11*s12*s13*s24**2*s26**2*s35**2*s36**2 -  \
    1024*q9a*s11**2*s14*s15*s22**2*s23**2*s36**3 + 1024*q9a*s11**2*s13*s16*s22**2*s23**2*s36**3 - 1024*q9a*s11**2*s12*s16*s22*s23**3*s36**3 - 512*q9a*s11**2*s13*s15*s22**2*s23*s24*s36**3 + 1536*q9a*s11**2*s12*s15*s22*s23**2*s24*s36**3 -  \
    512*q9a*s11**2*s14*s16*s22*s23**2*s24*s36**3 + 512*q9a*s11**2*s14*s15*s22*s23*s24**2*s36**3 + 512*q9a*s11**2*s12*s16*s23**2*s24**2*s36**3 - 512*q9a*s11**2*s12*s15*s23*s24**3*s36**3 + 512*q9a*s11**2*s13*s14*s22**2*s23*s25*s36**3 +  \
    1024*q9a*s11*s14*s15**2*s22**2*s23*s25*s36**3 - 1024*q9a*s11*s13*s15*s16*s22**2*s23*s25*s36**3 + 512*q9a*s11**2*s12*s14*s22*s23**2*s25*s36**3 + 1024*q9a*s11*s12*s15*s16*s22*s23**2*s25*s36**3 -  \
    1024*q9a*s11**2*s13**2*s22**2*s24*s25*s36**3 + 512*q9a*s11*s13*s15**2*s22**2*s24*s25*s36**3 + 2048*q9a*s11**2*s12*s13*s22*s23*s24*s25*s36**3 - 1536*q9a*s11**2*s14**2*s22*s23*s24*s25*s36**3 -  \
    1536*q9a*s11*s12*s15**2*s22*s23*s24*s25*s36**3 + 512*q9a*s11*s14*s15*s16*s22*s23*s24*s25*s36**3 - 2048*q9a*s11**2*s12**2*s23**2*s24*s25*s36**3 + 1024*q9a*s11**2*s13*s14*s22*s24**2*s25*s36**3 -  \
    512*q9a*s11*s14*s15**2*s22*s24**2*s25*s36**3 + 1536*q9a*s11**2*s12*s14*s23*s24**2*s25*s36**3 - 512*q9a*s11*s12*s15*s16*s23*s24**2*s25*s36**3 - 1024*q9a*s11**2*s12*s13*s24**3*s25*s36**3 + 512*q9a*s11*s12*s15**2*s24**3*s25*s36**3 -  \
    1536*q9a*s11*s13*s14*s15*s22**2*s25**2*s36**3 + 1536*q9a*s11*s13**2*s16*s22**2*s25**2*s36**3 - 512*q9a*s11*s12*s14*s15*s22*s23*s25**2*s36**3 - 2048*q9a*s11*s12*s13*s16*s22*s23*s25**2*s36**3 +  \
    512*q9a*s11*s14**2*s16*s22*s23*s25**2*s36**3 + 512*q9a*s11*s12**2*s16*s23**2*s25**2*s36**3 + 512*q9a*s11*s12*s13*s15*s22*s24*s25**2*s36**3 + 1024*q9a*s11*s14**2*s15*s22*s24*s25**2*s36**3 -  \
    1024*q9a*s11*s13*s14*s16*s22*s24*s25**2*s36**3 + 1536*q9a*s11*s12**2*s15*s23*s24*s25**2*s36**3 - 512*q9a*s11*s12*s14*s16*s23*s24*s25**2*s36**3 - 1024*q9a*s11*s12*s14*s15*s24**2*s25**2*s36**3 +  \
    1024*q9a*s11*s12*s13*s16*s24**2*s25**2*s36**3 + 1536*q9a*s11*s12*s13*s14*s22*s25**3*s36**3 - 512*q9a*s11*s14**3*s22*s25**3*s36**3 - 512*q9a*s11*s12**2*s14*s23*s25**3*s36**3 - 1024*q9a*s11*s12**2*s13*s24*s25**3*s36**3 +  \
    512*q9a*s11*s12*s14**2*s24*s25**3*s36**3 + 1024*q9a*s11**2*s13**2*s22**2*s23*s26*s36**3 - 3072*q9a*s11**2*s12*s13*s22*s23**2*s26*s36**3 + 1536*q9a*s11**2*s14**2*s22*s23**2*s26*s36**3 + 2048*q9a*s11**2*s12**2*s23**3*s26*s36**3 -  \
    1024*q9a*s11**2*s13*s14*s22*s23*s24*s26*s36**3 - 1536*q9a*s11**2*s12*s14*s23**2*s24*s26*s36**3 + 1024*q9a*s11**2*s12*s13*s23*s24**2*s26*s36**3 - 512*q9a*s11*s13**2*s15*s22**2*s25*s26*s36**3 +  \
    2048*q9a*s11*s12*s13*s15*s22*s23*s25*s26*s36**3 - 1024*q9a*s11*s14**2*s15*s22*s23*s25*s26*s36**3 - 1536*q9a*s11*s12**2*s15*s23**2*s25*s26*s36**3 + 512*q9a*s11*s13*s14*s15*s22*s24*s25*s26*s36**3;
    v3_9= \
    1024*q9a*s11*s12*s14*s15*s23*s24*s25*s26*s36**3 - 512*q9a*s11*s12*s13*s15*s24**2*s25*s26*s36**3 - 1024*q9a*s11*s12*s13**2*s22*s25**2*s26*s36**3 + 512*q9a*s11*s13*s14**2*s22*s25**2*s26*s36**3 +  \
    1024*q9a*s11*s12**2*s13*s23*s25**2*s26*s36**3 - 512*q9a*s11*s12*s13*s14*s24*s25**2*s26*s36**3 - 2048*q8a*s11**2*s14*s15*s22**2*s23*s33*s36**3 + 2048*q8a*s11**2*s13*s16*s22**2*s23*s33*s36**3 -  \
    3072*q8a*s11**2*s12*s16*s22*s23**2*s33*s36**3 + 1024*q7a*s11**2*s16*s22**2*s23**2*s33*s36**3 - 512*q8a*s11**2*s13*s15*s22**2*s24*s33*s36**3 + 3072*q8a*s11**2*s12*s15*s22*s23*s24*s33*s36**3 -  \
    1024*q8a*s11**2*s14*s16*s22*s23*s24*s33*s36**3 - 512*q7a*s11**2*s15*s22**2*s23*s24*s33*s36**3 + 512*q8a*s11**2*s14*s15*s22*s24**2*s33*s36**3 + 1024*q8a*s11**2*s12*s16*s23*s24**2*s33*s36**3 -  \
    512*q8a*s11**2*s12*s15*s24**3*s33*s36**3 + 512*q8a*s11**2*s13*s14*s22**2*s25*s33*s36**3 + 1024*q8a*s11*s14*s15**2*s22**2*s25*s33*s36**3 - 1024*q8a*s11*s13*s15*s16*s22**2*s25*s33*s36**3 +  \
    1024*q8a*s11**2*s12*s14*s22*s23*s25*s33*s36**3 + 2048*q8a*s11*s12*s15*s16*s22*s23*s25*s33*s36**3 + 512*q7a*s11**2*s14*s22**2*s23*s25*s33*s36**3 - 1024*q7a*s11*s15*s16*s22**2*s23*s25*s33*s36**3 +  \
    2048*q8a*s11**2*s12*s13*s22*s24*s25*s33*s36**3 - 1536*q8a*s11**2*s14**2*s22*s24*s25*s33*s36**3 - 1536*q8a*s11*s12*s15**2*s22*s24*s25*s33*s36**3 + 512*q8a*s11*s14*s15*s16*s22*s24*s25*s33*s36**3 -  \
    2048*q7a*s11**2*s13*s22**2*s24*s25*s33*s36**3 + 512*q7a*s11*s15**2*s22**2*s24*s25*s33*s36**3 - 4096*q8a*s11**2*s12**2*s23*s24*s25*s33*s36**3 + 2048*q7a*s11**2*s12*s22*s23*s24*s25*s33*s36**3 +  \
    1536*q8a*s11**2*s12*s14*s24**2*s25*s33*s36**3 - 512*q8a*s11*s12*s15*s16*s24**2*s25*s33*s36**3 + 1024*q7a*s11**2*s14*s22*s24**2*s25*s33*s36**3 - 1024*q7a*s11**2*s12*s24**3*s25*s33*s36**3 -  \
    512*q8a*s11*s12*s14*s15*s22*s25**2*s33*s36**3 - 2048*q8a*s11*s12*s13*s16*s22*s25**2*s33*s36**3 + 512*q8a*s11*s14**2*s16*s22*s25**2*s33*s36**3 - 1536*q7a*s11*s14*s15*s22**2*s25**2*s33*s36**3 +  \
    3072*q7a*s11*s13*s16*s22**2*s25**2*s33*s36**3 + 1024*q8a*s11*s12**2*s16*s23*s25**2*s33*s36**3 - 2048*q7a*s11*s12*s16*s22*s23*s25**2*s33*s36**3 + 1536*q8a*s11*s12**2*s15*s24*s25**2*s33*s36**3 -  \
    512*q8a*s11*s12*s14*s16*s24*s25**2*s33*s36**3 + 512*q7a*s11*s12*s15*s22*s24*s25**2*s33*s36**3 - 1024*q7a*s11*s14*s16*s22*s24*s25**2*s33*s36**3 + 1024*q7a*s11*s12*s16*s24**2*s25**2*s33*s36**3 -  \
    512*q8a*s11*s12**2*s14*s25**3*s33*s36**3 + 1536*q7a*s11*s12*s14*s22*s25**3*s33*s36**3 - 1024*q7a*s11*s12**2*s24*s25**3*s33*s36**3 + 1024*q8a*s11**2*s13**2*s22**2*s26*s33*s36**3 - 6144*q8a*s11**2*s12*s13*s22*s23*s26*s33*s36**3 +  \
    3072*q8a*s11**2*s14**2*s22*s23*s26*s33*s36**3 + 2048*q7a*s11**2*s13*s22**2*s23*s26*s33*s36**3 + 6144*q8a*s11**2*s12**2*s23**2*s26*s33*s36**3 - 3072*q7a*s11**2*s12*s22*s23**2*s26*s33*s36**3 -  \
    1024*q8a*s11**2*s13*s14*s22*s24*s26*s33*s36**3 - 3072*q8a*s11**2*s12*s14*s23*s24*s26*s33*s36**3 - 1024*q7a*s11**2*s14*s22*s23*s24*s26*s33*s36**3 + 1024*q8a*s11**2*s12*s13*s24**2*s26*s33*s36**3 +  \
    1024*q7a*s11**2*s12*s23*s24**2*s26*s33*s36**3 + 2048*q8a*s11*s12*s13*s15*s22*s25*s26*s33*s36**3 - 1024*q8a*s11*s14**2*s15*s22*s25*s26*s33*s36**3 - 1024*q7a*s11*s13*s15*s22**2*s25*s26*s33*s36**3 -  \
    3072*q8a*s11*s12**2*s15*s23*s25*s26*s33*s36**3 + 2048*q7a*s11*s12*s15*s22*s23*s25*s26*s33*s36**3 + 1024*q8a*s11*s12*s14*s15*s24*s25*s26*s33*s36**3 + 512*q7a*s11*s14*s15*s22*s24*s25*s26*s33*s36**3 -  \
    512*q7a*s11*s12*s15*s24**2*s25*s26*s33*s36**3 + 1024*q8a*s11*s12**2*s13*s25**2*s26*s33*s36**3 - 2048*q7a*s11*s12*s13*s22*s25**2*s26*s33*s36**3 + 512*q7a*s11*s14**2*s22*s25**2*s26*s33*s36**3 +  \
    1024*q7a*s11*s12**2*s23*s25**2*s26*s33*s36**3 - 512*q7a*s11*s12*s14*s24*s25**2*s26*s33*s36**3 + 2048*q8a*s11**2*s13*s15*s22**2*s23*s34*s36**3 - 3072*q8a*s11**2*s12*s15*s22*s23**2*s34*s36**3 +  \
    3072*q8a*s11**2*s14*s16*s22*s23**2*s34*s36**3 + 1024*q7a*s11**2*s15*s22**2*s23**2*s34*s36**3 - 1024*q8a*s11**2*s14*s15*s22*s23*s24*s34*s36**3 - 1024*q8a*s11**2*s13*s16*s22*s23*s24*s34*s36**3 -  \
    1536*q8a*s11**2*s12*s16*s23**2*s24*s34*s36**3 - 512*q7a*s11**2*s16*s22*s23**2*s24*s34*s36**3 + 1024*q8a*s11**2*s12*s15*s23*s24**2*s34*s36**3 - 1024*q8a*s11*s13*s15**2*s22**2*s25*s34*s36**3 -  \
    2048*q8a*s11**2*s12*s13*s22*s23*s25*s34*s36**3 + 1024*q8a*s11**2*s14**2*s22*s23*s25*s34*s36**3 + 2048*q8a*s11*s12*s15**2*s22*s23*s25*s34*s36**3 - 2048*q8a*s11*s14*s15*s16*s22*s23*s25*s34*s36**3 -  \
    1024*q7a*s11*s15**2*s22**2*s23*s25*s34*s36**3 + 3072*q8a*s11**2*s12**2*s23**2*s25*s34*s36**3 - 1024*q7a*s11**2*s12*s22*s23**2*s25*s34*s36**3 + 512*q8a*s11*s14*s15**2*s22*s24*s25*s34*s36**3 +  \
    512*q8a*s11*s13*s15*s16*s22*s24*s25*s34*s36**3 - 1024*q8a*s11**2*s12*s14*s23*s24*s25*s34*s36**3 + 1024*q8a*s11*s12*s15*s16*s23*s24*s25*s34*s36**3 + 512*q7a*s11*s15*s16*s22*s23*s24*s25*s34*s36**3 -  \
    512*q8a*s11*s12*s15**2*s24**2*s25*s34*s36**3 - 512*q8a*s11*s14**2*s15*s22*s25**2*s34*s36**3 + 1024*q8a*s11*s13*s14*s16*s22*s25**2*s34*s36**3 + 2048*q7a*s11*s13*s15*s22**2*s25**2*s34*s36**3 -  \
    2048*q8a*s11*s12**2*s15*s23*s25**2*s34*s36**3 + 1024*q7a*s11*s14*s16*s22*s23*s25**2*s34*s36**3 + 512*q8a*s11*s12*s14*s15*s24*s25**2*s34*s36**3 - 512*q8a*s11*s12*s13*s16*s24*s25**2*s34*s36**3 -  \
    512*q7a*s11*s14*s15*s22*s24*s25**2*s34*s36**3 - 1024*q7a*s11*s13*s16*s22*s24*s25**2*s34*s36**3 - 512*q7a*s11*s12*s16*s23*s24*s25**2*s34*s36**3 + 512*q7a*s11*s12*s15*s24**2*s25**2*s34*s36**3 +  \
    1024*q8a*s11*s12**2*s13*s25**3*s34*s36**3 - 2048*q7a*s11*s12*s13*s22*s25**3*s34*s36**3 + 512*q7a*s11*s14**2*s22*s25**3*s34*s36**3 + 1024*q7a*s11*s12**2*s23*s25**3*s34*s36**3 - 512*q7a*s11*s12*s14*s24*s25**3*s34*s36**3 -  \
    1024*q8a*s11**2*s13*s14*s22*s23*s26*s34*s36**3 - 1536*q8a*s11**2*s12*s14*s23**2*s26*s34*s36**3 - 512*q7a*s11**2*s14*s22*s23**2*s26*s34*s36**3 + 2048*q8a*s11**2*s12*s13*s23*s24*s26*s34*s36**3 +  \
    1024*q7a*s11**2*s12*s23**2*s24*s26*s34*s36**3 + 512*q8a*s11*s13*s14*s15*s22*s25*s26*s34*s36**3 + 1024*q8a*s11*s12*s14*s15*s23*s25*s26*s34*s36**3 + 512*q7a*s11*s14*s15*s22*s23*s25*s26*s34*s36**3 -  \
    1024*q8a*s11*s12*s13*s15*s24*s25*s26*s34*s36**3 - 1024*q7a*s11*s12*s15*s23*s24*s25*s26*s34*s36**3 - 512*q8a*s11*s12*s13*s14*s25**2*s26*s34*s36**3 - 1024*q7a*s11*s13*s14*s22*s25**2*s26*s34*s36**3 -  \
    512*q7a*s11*s12*s14*s23*s25**2*s26*s34*s36**3 + 2048*q7a*s11*s12*s13*s24*s25**2*s26*s34*s36**3 + 1024*q8a*s11**2*s13**2*s22**2*s24*s35*s36**3 - 4096*q8a*s11**2*s12*s13*s22*s23*s24*s35*s36**3 +  \
    2048*q8a*s11**2*s14**2*s22*s23*s24*s35*s36**3 + 2048*q7a*s11**2*s13*s22**2*s23*s24*s35*s36**3 + 3072*q8a*s11**2*s12**2*s23**2*s24*s35*s36**3 - 2048*q7a*s11**2*s12*s22*s23**2*s24*s35*s36**3 -  \
    1024*q8a*s11**2*s13*s14*s22*s24**2*s35*s36**3 - 2048*q8a*s11**2*s12*s14*s23*s24**2*s35*s36**3 - 1024*q7a*s11**2*s14*s22*s23*s24**2*s35*s36**3 + 1024*q8a*s11**2*s12*s13*s24**3*s35*s36**3 + 1024*q7a*s11**2*s12*s23*s24**3*s35*s36**3 -  \
    512*q8a*s11*s13**2*s16*s22**2*s25*s35*s36**3 + 2048*q8a*s11*s12*s13*s16*s22*s23*s25*s35*s36**3 - 1024*q8a*s11*s14**2*s16*s22*s23*s25*s35*s36**3 - 1024*q7a*s11*s13*s16*s22**2*s23*s25*s35*s36**3 -  \
    1536*q8a*s11*s12**2*s16*s23**2*s25*s35*s36**3 + 1024*q7a*s11*s12*s16*s22*s23**2*s25*s35*s36**3 + 1024*q8a*s11*s12*s13*s15*s22*s24*s25*s35*s36**3 - 512*q8a*s11*s14**2*s15*s22*s24*s25*s35*s36**3 +  \
    512*q8a*s11*s13*s14*s16*s22*s24*s25*s35*s36**3 - 1024*q7a*s11*s13*s15*s22**2*s24*s25*s35*s36**3 - 1024*q8a*s11*s12**2*s15*s23*s24*s25*s35*s36**3 + 1024*q8a*s11*s12*s14*s16*s23*s24*s25*s35*s36**3 +  \
    1024*q7a*s11*s12*s15*s22*s23*s24*s25*s35*s36**3 + 512*q7a*s11*s14*s16*s22*s23*s24*s25*s35*s36**3 + 512*q8a*s11*s12*s14*s15*s24**2*s25*s35*s36**3 - 512*q8a*s11*s12*s13*s16*s24**2*s25*s35*s36**3 +  \
    512*q7a*s11*s14*s15*s22*s24**2*s25*s35*s36**3 - 512*q7a*s11*s12*s16*s23*s24**2*s25*s35*s36**3 - 512*q7a*s11*s12*s15*s24**3*s25*s35*s36**3 - 1024*q8a*s11*s12*s13*s14*s22*s25**2*s35*s36**3 +  \
    512*q8a*s11*s14**3*s22*s25**2*s35*s36**3 + 1024*q7a*s11*s13*s14*s22**2*s25**2*s35*s36**3 + 1024*q8a*s11*s12**2*s14*s23*s25**2*s35*s36**3 - 1024*q7a*s11*s12*s14*s22*s23*s25**2*s35*s36**3 -  \
    512*q8a*s11*s12*s14**2*s24*s25**2*s35*s36**3 - 512*q7a*s11*s14**2*s22*s24*s25**2*s35*s36**3 + 512*q7a*s11*s12*s14*s24**2*s25**2*s35*s36**3 - 512*q8a*s11*s13**2*s15*s22**2*s26*s35*s36**3 +  \
    2048*q8a*s11*s12*s13*s15*s22*s23*s26*s35*s36**3 - 1024*q8a*s11*s14**2*s15*s22*s23*s26*s35*s36**3 - 1024*q7a*s11*s13*s15*s22**2*s23*s26*s35*s36**3 - 1536*q8a*s11*s12**2*s15*s23**2*s26*s35*s36**3 +  \
    1024*q7a*s11*s12*s15*s22*s23**2*s26*s35*s36**3 + 512*q8a*s11*s13*s14*s15*s22*s24*s26*s35*s36**3 + 1024*q8a*s11*s12*s14*s15*s23*s24*s26*s35*s36**3 + 512*q7a*s11*s14*s15*s22*s23*s24*s26*s35*s36**3 -  \
    512*q8a*s11*s12*s13*s15*s24**2*s26*s35*s36**3 - 512*q7a*s11*s12*s15*s23*s24**2*s26*s35*s36**3 - 2048*q8a*s11*s12*s13**2*s22*s25*s26*s35*s36**3 + 1024*q8a*s11*s13*s14**2*s22*s25*s26*s35*s36**3 +  \
    3072*q7a*s11*s13**2*s22**2*s25*s26*s35*s36**3 + 2048*q8a*s11*s12**2*s13*s23*s25*s26*s35*s36**3 - 4096*q7a*s11*s12*s13*s22*s23*s25*s26*s35*s36**3 + 1024*q7a*s11*s14**2*s22*s23*s25*s26*s35*s36**3 +  \
    1024*q7a*s11*s12**2*s23**2*s25*s26*s35*s36**3 - 1024*q8a*s11*s12*s13*s14*s24*s25*s26*s35*s36**3 - 2048*q7a*s11*s13*s14*s22*s24*s25*s26*s35*s36**3 - 1024*q7a*s11*s12*s14*s23*s24*s25*s26*s35*s36**3 +  \
    2048*q7a*s11*s12*s13*s24**2*s25*s26*s35*s36**3 - 1024*q8a*s11**2*s13**2*s22**2*s23*s36**4 + 3072*q8a*s11**2*s12*s13*s22*s23**2*s36**4 - 1536*q8a*s11**2*s14**2*s22*s23**2*s36**4 - 1024*q7a*s11**2*s13*s22**2*s23**2*s36**4 -  \
    2048*q8a*s11**2*s12**2*s23**3*s36**4 + 1024*q7a*s11**2*s12*s22*s23**3*s36**4 + 1024*q8a*s11**2*s13*s14*s22*s23*s24*s36**4 + 1536*q8a*s11**2*s12*s14*s23**2*s24*s36**4 + 512*q7a*s11**2*s14*s22*s23**2*s24*s36**4 -  \
    1024*q8a*s11**2*s12*s13*s23*s24**2*s36**4 - 512*q7a*s11**2*s12*s23**2*s24**2*s36**4 + 512*q8a*s11*s13**2*s15*s22**2*s25*s36**4 - 2048*q8a*s11*s12*s13*s15*s22*s23*s25*s36**4 + 1024*q8a*s11*s14**2*s15*s22*s23*s25*s36**4 +  \
    1024*q7a*s11*s13*s15*s22**2*s23*s25*s36**4 + 1536*q8a*s11*s12**2*s15*s23**2*s25*s36**4 - 1024*q7a*s11*s12*s15*s22*s23**2*s25*s36**4 - 512*q8a*s11*s13*s14*s15*s22*s24*s25*s36**4 - 1024*q8a*s11*s12*s14*s15*s23*s24*s25*s36**4 -  \
    512*q7a*s11*s14*s15*s22*s23*s24*s25*s36**4 + 512*q8a*s11*s12*s13*s15*s24**2*s25*s36**4 + 512*q7a*s11*s12*s15*s23*s24**2*s25*s36**4 + 1024*q8a*s11*s12*s13**2*s22*s25**2*s36**4 - 512*q8a*s11*s13*s14**2*s22*s25**2*s36**4 -  \
    1536*q7a*s11*s13**2*s22**2*s25**2*s36**4 - 1024*q8a*s11*s12**2*s13*s23*s25**2*s36**4 + 2048*q7a*s11*s12*s13*s22*s23*s25**2*s36**4 - 512*q7a*s11*s14**2*s22*s23*s25**2*s36**4 - 512*q7a*s11*s12**2*s23**2*s25**2*s36**4 +  \
    512*q8a*s11*s12*s13*s14*s24*s25**2*s36**4 + 1024*q7a*s11*s13*s14*s22*s24*s25**2*s36**4 + 512*q7a*s11*s12*s14*s23*s24*s25**2*s36**4 - 1024*q7a*s11*s12*s13*s24**2*s25**2*s36**4;
    v3=v3_0+v3_1+v3_2+v3_3+v3_4+v3_5+v3_6+v3_7+v3_8+v3_9

    v4_0= \
    512*s11**4*s22**4*s33**4 - 1024*s11**3*s16*s22**3*s26*s33**4 + 1024*s11**3*s12*s22**2*s26**2*s33**4 + 512*s11**2*s16**2*s22**2*s26**2*s33**4 - 1024*s11**2*s12*s16*s22*s26**3*s33**4 + 512*s11**2*s12**2*s26**4*s33**4 -  \
    1024*s11**4*s22**3*s24*s33**3*s34 + 512*s11**3*s16*s22**3*s25*s33**3*s34 + 512*s11**3*s15*s22**3*s26*s33**3*s34 + 1536*s11**3*s16*s22**2*s24*s26*s33**3*s34 - 1024*s11**3*s12*s22**2*s25*s26*s33**3*s34 -  \
    512*s11**2*s16**2*s22**2*s25*s26*s33**3*s34 - 512*s11**3*s14*s22**2*s26**2*s33**3*s34 - 512*s11**2*s15*s16*s22**2*s26**2*s33**3*s34 - 1024*s11**3*s12*s22*s24*s26**2*s33**3*s34 - 512*s11**2*s16**2*s22*s24*s26**2*s33**3*s34 +  \
    1536*s11**2*s12*s16*s22*s25*s26**2*s33**3*s34 + 512*s11**2*s12*s15*s22*s26**3*s33**3*s34 + 512*s11**2*s14*s16*s22*s26**3*s33**3*s34 + 512*s11**2*s12*s16*s24*s26**3*s33**3*s34 - 1024*s11**2*s12**2*s25*s26**3*s33**3*s34 -  \
    512*s11**2*s12*s14*s26**4*s33**3*s34 + 1024*s11**4*s22**3*s23*s33**2*s34**2 + 512*s11**4*s22**2*s24**2*s33**2*s34**2 - 512*s11**3*s15*s22**3*s25*s33**2*s34**2 - 512*s11**3*s16*s22**2*s24*s25*s33**2*s34**2 +  \
    512*s11**3*s12*s22**2*s25**2*s33**2*s34**2 - 1536*s11**3*s16*s22**2*s23*s26*s33**2*s34**2 - 512*s11**3*s15*s22**2*s24*s26*s33**2*s34**2 - 512*s11**3*s16*s22*s24**2*s26*s33**2*s34**2 + 1024*s11**3*s14*s22**2*s25*s26*s33**2*s34**2 +  \
    512*s11**2*s15*s16*s22**2*s25*s26*s33**2*s34**2 + 512*s11**2*s16**2*s22*s24*s25*s26*s33**2*s34**2 - 512*s11**2*s12*s16*s22*s25**2*s26*s33**2*s34**2 + 512*s11**3*s13*s22**2*s26**2*s33**2*s34**2 +  \
    1024*s11**3*s12*s22*s23*s26**2*s33**2*s34**2 + 512*s11**2*s16**2*s22*s23*s26**2*s33**2*s34**2 + 512*s11**2*s15*s16*s22*s24*s26**2*s33**2*s34**2 + 512*s11**3*s12*s24**2*s26**2*s33**2*s34**2 -  \
    512*s11**2*s12*s15*s22*s25*s26**2*s33**2*s34**2 - 1024*s11**2*s14*s16*s22*s25*s26**2*s33**2*s34**2 - 512*s11**2*s12*s16*s24*s25*s26**2*s33**2*s34**2 + 512*s11**2*s12**2*s25**2*s26**2*s33**2*s34**2 -  \
    512*s11**2*s13*s16*s22*s26**3*s33**2*s34**2 - 512*s11**2*s12*s16*s23*s26**3*s33**2*s34**2 - 512*s11**2*s12*s15*s24*s26**3*s33**2*s34**2 + 1024*s11**2*s12*s14*s25*s26**3*s33**2*s34**2 + 512*s11**2*s12*s13*s26**4*s33**2*s34**2 -  \
    1024*s11**4*s22**2*s23*s24*s33*s34**3 + 512*s11**3*s16*s22**2*s23*s25*s33*s34**3 + 512*s11**3*s15*s22**2*s24*s25*s33*s34**3 - 512*s11**3*s14*s22**2*s25**2*s33*s34**3 + 512*s11**3*s15*s22**2*s23*s26*s33*s34**3 +  \
    1024*s11**3*s16*s22*s23*s24*s26*s33*s34**3 - 1024*s11**3*s13*s22**2*s25*s26*s33*s34**3 - 512*s11**2*s16**2*s22*s23*s25*s26*s33*s34**3 - 512*s11**2*s15*s16*s22*s24*s25*s26*s33*s34**3 +  \
    512*s11**2*s14*s16*s22*s25**2*s26*s33*s34**3 - 512*s11**2*s15*s16*s22*s23*s26**2*s33*s34**3 - 1024*s11**3*s12*s23*s24*s26**2*s33*s34**3 + 1024*s11**2*s13*s16*s22*s25*s26**2*s33*s34**3;
    v4_1= \
    512*s11**2*s12*s16*s23*s25*s26**2*s33*s34**3 + 512*s11**2*s12*s15*s24*s25*s26**2*s33*s34**3 - 512*s11**2*s12*s14*s25**2*s26**2*s33*s34**3 + 512*s11**2*s12*s15*s23*s26**3*s33*s34**3 - 1024*s11**2*s12*s13*s25*s26**3*s33*s34**3 +  \
    512*s11**4*s22**2*s23**2*s34**4 - 512*s11**3*s15*s22**2*s23*s25*s34**4 + 512*s11**3*s13*s22**2*s25**2*s34**4 - 512*s11**3*s16*s22*s23**2*s26*s34**4 + 512*s11**2*s15*s16*s22*s23*s25*s26*s34**4 -  \
    512*s11**2*s13*s16*s22*s25**2*s26*s34**4 + 512*s11**3*s12*s23**2*s26**2*s34**4 - 512*s11**2*s12*s15*s23*s25*s26**2*s34**4 + 512*s11**2*s12*s13*s25**2*s26**2*s34**4 - 1024*s11**3*s15*s22**4*s33**3*s35 +  \
    512*s11**3*s16*s22**3*s24*s33**3*s35 + 1024*s11**3*s12*s22**3*s25*s33**3*s35 - 512*s11**2*s16**2*s22**3*s25*s33**3*s35 + 512*s11**3*s14*s22**3*s26*s33**3*s35 + 1536*s11**2*s15*s16*s22**3*s26*s33**3*s35 -  \
    1024*s11**3*s12*s22**2*s24*s26*s33**3*s35 - 512*s11**2*s16**2*s22**2*s24*s26*s33**3*s35 - 512*s11**2*s12*s16*s22**2*s25*s26*s33**3*s35 + 512*s11*s16**3*s22**2*s25*s26*s33**3*s35 - 1536*s11**2*s12*s15*s22**2*s26**2*s33**3*s35 -  \
    512*s11**2*s14*s16*s22**2*s26**2*s33**3*s35 - 512*s11*s15*s16**2*s22**2*s26**2*s33**3*s35 + 1536*s11**2*s12*s16*s22*s24*s26**2*s33**3*s35 + 1024*s11**2*s12**2*s22*s25*s26**2*s33**3*s35 - 1024*s11*s12*s16**2*s22*s25*s26**2*s33**3*s35 +  \
    512*s11**2*s12*s14*s22*s26**3*s33**3*s35 + 1024*s11*s12*s15*s16*s22*s26**3*s33**3*s35 - 1024*s11**2*s12**2*s24*s26**3*s33**3*s35 + 512*s11*s12**2*s16*s25*s26**3*s33**3*s35 - 512*s11*s12**2*s15*s26**4*s33**3*s35 -  \
    1024*s11**3*s16*s22**3*s23*s33**2*s34*s35 + 1536*s11**3*s15*s22**3*s24*s33**2*s34*s35 - 512*s11**3*s16*s22**2*s24**2*s33**2*s34*s35 - 1536*s11**3*s14*s22**3*s25*s33**2*s34*s35 + 512*s11**2*s15*s16*s22**3*s25*s33**2*s34*s35 +  \
    512*s11**2*s16**2*s22**2*s24*s25*s33**2*s34*s35 - 512*s11**2*s12*s16*s22**2*s25**2*s33**2*s34*s35 - 1024*s11**3*s13*s22**3*s26*s33**2*s34*s35 - 512*s11**2*s15**2*s22**3*s26*s33**2*s34*s35 +  \
    2048*s11**3*s12*s22**2*s23*s26*s33**2*s34*s35 + 1024*s11**2*s16**2*s22**2*s23*s26*s33**2*s34*s35 + 512*s11**3*s14*s22**2*s24*s26*s33**2*s34*s35 - 2048*s11**2*s15*s16*s22**2*s24*s26*s33**2*s34*s35 +  \
    512*s11**2*s16**2*s22*s24**2*s26*s33**2*s34*s35 + 512*s11**2*s12*s15*s22**2*s25*s26*s33**2*s34*s35 + 1024*s11**2*s14*s16*s22**2*s25*s26*s33**2*s34*s35 - 512*s11*s15*s16**2*s22**2*s25*s26*s33**2*s34*s35 -  \
    512*s11*s16**3*s22*s24*s25*s26*s33**2*s34*s35 + 512*s11*s12*s16**2*s22*s25**2*s26*s33**2*s34*s35 + 512*s11**2*s14*s15*s22**2*s26**2*s33**2*s34*s35 + 1024*s11**2*s13*s16*s22**2*s26**2*s33**2*s34*s35 +  \
    512*s11*s15**2*s16*s22**2*s26**2*s33**2*s34*s35 - 3072*s11**2*s12*s16*s22*s23*s26**2*s33**2*s34*s35 + 1536*s11**2*s12*s15*s22*s24*s26**2*s33**2*s34*s35 - 512*s11**2*s14*s16*s22*s24*s26**2*s33**2*s34*s35 +  \
    512*s11*s15*s16**2*s22*s24*s26**2*s33**2*s34*s35 - 512*s11**2*s12*s16*s24**2*s26**2*s33**2*s34*s35 - 1536*s11**2*s12*s14*s22*s25*s26**2*s33**2*s34*s35 + 512*s11*s14*s16**2*s22*s25*s26**2*s33**2*s34*s35 +  \
    512*s11*s12*s16**2*s24*s25*s26**2*s33**2*s34*s35 - 512*s11*s12**2*s16*s25**2*s26**2*s33**2*s34*s35 - 1024*s11**2*s12*s13*s22*s26**3*s33**2*s34*s35 - 512*s11*s12*s15**2*s22*s26**3*s33**2*s34*s35 -  \
    512*s11*s14*s15*s16*s22*s26**3*s33**2*s34*s35 + 2048*s11**2*s12**2*s23*s26**3*s33**2*s34*s35 + 512*s11**2*s12*s14*s24*s26**3*s33**2*s34*s35 - 512*s11*s12*s15*s16*s24*s26**3*s33**2*s34*s35 +  \
    512*s11*s12**2*s15*s25*s26**3*s33**2*s34*s35 - 512*s11*s12*s14*s16*s25*s26**3*s33**2*s34*s35 + 512*s11*s12*s14*s15*s26**4*s33**2*s34*s35 - 1024*s11**3*s15*s22**3*s23*s33*s34**2*s35 +  \
    1536*s11**3*s16*s22**2*s23*s24*s33*s34**2*s35 - 512*s11**3*s15*s22**2*s24**2*s33*s34**2*s35 + 2048*s11**3*s13*s22**3*s25*s33*s34**2*s35 - 1024*s11**3*s12*s22**2*s23*s25*s33*s34**2*s35 -  \
    512*s11**2*s16**2*s22**2*s23*s25*s33*s34**2*s35 + 512*s11**3*s14*s22**2*s24*s25*s33*s34**2*s35 - 512*s11**2*s15*s16*s22**2*s24*s25*s33*s34**2*s35 + 512*s11**2*s14*s16*s22**2*s25**2*s33*s34**2*s35 -  \
    1536*s11**3*s14*s22**2*s23*s26*s33*s34**2*s35 + 1536*s11**2*s15*s16*s22**2*s23*s26*s33*s34**2*s35 + 512*s11**2*s15**2*s22**2*s24*s26*s33*s34**2*s35 - 1536*s11**2*s16**2*s22*s23*s24*s26*s33*s34**2*s35 +  \
    512*s11**2*s15*s16*s22*s24**2*s26*s33*s34**2*s35 - 512*s11**2*s14*s15*s22**2*s25*s26*s33*s34**2*s35 - 1536*s11**2*s13*s16*s22**2*s25*s26*s33*s34**2*s35 + 1024*s11**2*s12*s16*s22*s23*s25*s26*s33*s34**2*s35 +  \
    512*s11*s16**3*s22*s23*s25*s26*s33*s34**2*s35 - 512*s11**2*s14*s16*s22*s24*s25*s26*s33*s34**2*s35 + 512*s11*s15*s16**2*s22*s24*s25*s26*s33*s34**2*s35 - 512*s11*s14*s16**2*s22*s25**2*s26*s33*s34**2*s35 -  \
    512*s11**2*s13*s15*s22**2*s26**2*s33*s34**2*s35 - 1024*s11**2*s12*s15*s22*s23*s26**2*s33*s34**2*s35 + 1536*s11**2*s14*s16*s22*s23*s26**2*s33*s34**2*s35 - 512*s11*s15*s16**2*s22*s23*s26**2*s33*s34**2*s35 -  \
    512*s11*s15**2*s16*s22*s24*s26**2*s33*s34**2*s35 + 1536*s11**2*s12*s16*s23*s24*s26**2*s33*s34**2*s35 - 512*s11**2*s12*s15*s24**2*s26**2*s33*s34**2*s35 + 2048*s11**2*s12*s13*s22*s25*s26**2*s33*s34**2*s35 +  \
    512*s11*s14*s15*s16*s22*s25*s26**2*s33*s34**2*s35 - 512*s11*s13*s16**2*s22*s25*s26**2*s33*s34**2*s35 - 1024*s11**2*s12**2*s23*s25*s26**2*s33*s34**2*s35 - 512*s11*s12*s16**2*s23*s25*s26**2*s33*s34**2*s35 +  \
    512*s11**2*s12*s14*s24*s25*s26**2*s33*s34**2*s35 - 512*s11*s12*s15*s16*s24*s25*s26**2*s33*s34**2*s35 + 512*s11*s12*s14*s16*s25**2*s26**2*s33*s34**2*s35 + 512*s11*s13*s15*s16*s22*s26**3*s33*s34**2*s35 -  \
    1536*s11**2*s12*s14*s23*s26**3*s33*s34**2*s35 + 512*s11*s12*s15*s16*s23*s26**3*s33*s34**2*s35 + 512*s11*s12*s15**2*s24*s26**3*s33*s34**2*s35 - 512*s11*s12*s14*s15*s25*s26**3*s33*s34**2*s35 +  \
    512*s11*s12*s13*s16*s25*s26**3*s33*s34**2*s35 - 512*s11*s12*s13*s15*s26**4*s33*s34**2*s35 - 1024*s11**3*s16*s22**2*s23**2*s34**3*s35 + 512*s11**3*s15*s22**2*s23*s24*s34**3*s35 + 512*s11**3*s14*s22**2*s23*s25*s34**3*s35 +  \
    512*s11**2*s15*s16*s22**2*s23*s25*s34**3*s35 - 1024*s11**3*s13*s22**2*s24*s25*s34**3*s35 - 512*s11**2*s13*s16*s22**2*s25**2*s34**3*s35 + 1024*s11**3*s13*s22**2*s23*s26*s34**3*s35 - 512*s11**2*s15**2*s22**2*s23*s26*s34**3*s35 +  \
    1024*s11**2*s16**2*s22*s23**2*s26*s34**3*s35 - 512*s11**2*s15*s16*s22*s23*s24*s26*s34**3*s35 + 512*s11**2*s13*s15*s22**2*s25*s26*s34**3*s35 - 512*s11**2*s14*s16*s22*s23*s25*s26*s34**3*s35 -  \
    512*s11*s15*s16**2*s22*s23*s25*s26*s34**3*s35 + 1024*s11**2*s13*s16*s22*s24*s25*s26*s34**3*s35 + 512*s11*s13*s16**2*s22*s25**2*s26*s34**3*s35 - 1024*s11**2*s13*s16*s22*s23*s26**2*s34**3*s35 +  \
    512*s11*s15**2*s16*s22*s23*s26**2*s34**3*s35 - 1024*s11**2*s12*s16*s23**2*s26**2*s34**3*s35 + 512*s11**2*s12*s15*s23*s24*s26**2*s34**3*s35 - 512*s11*s13*s15*s16*s22*s25*s26**2*s34**3*s35 +  \
    512*s11**2*s12*s14*s23*s25*s26**2*s34**3*s35 + 512*s11*s12*s15*s16*s23*s25*s26**2*s34**3*s35 - 1024*s11**2*s12*s13*s24*s25*s26**2*s34**3*s35 - 512*s11*s12*s13*s16*s25**2*s26**2*s34**3*s35 +  \
    1024*s11**2*s12*s13*s23*s26**3*s34**3*s35 - 512*s11*s12*s15**2*s23*s26**3*s34**3*s35 + 512*s11*s12*s13*s15*s25*s26**3*s34**3*s35 + 1024*s11**3*s13*s22**4*s33**2*s35**2 + 512*s11**2*s15**2*s22**4*s33**2*s35**2 -  \
    1024*s11**3*s12*s22**3*s23*s33**2*s35**2 + 512*s11**2*s16**2*s22**3*s23*s33**2*s35**2 - 512*s11**3*s14*s22**3*s24*s33**2*s35**2 - 512*s11**2*s15*s16*s22**3*s24*s33**2*s35**2 + 512*s11**3*s12*s22**2*s24**2*s33**2*s35**2 -  \
    1024*s11**2*s12*s15*s22**3*s25*s33**2*s35**2 + 1024*s11**2*s14*s16*s22**3*s25*s33**2*s35**2 - 512*s11**2*s12*s16*s22**2*s24*s25*s33**2*s35**2 + 512*s11**2*s12**2*s22**2*s25**2*s33**2*s35**2 - 512*s11**2*s14*s15*s22**3*s26*s33**2*s35**2 -  \
    1536*s11**2*s13*s16*s22**3*s26*s33**2*s35**2 - 512*s11*s15**2*s16*s22**3*s26*s33**2*s35**2 + 512*s11**2*s12*s16*s22**2*s23*s26*s33**2*s35**2 - 512*s11*s16**3*s22**2*s23*s26*s33**2*s35**2;
    v4_2= \
    1024*s11**2*s12*s15*s22**2*s24*s26*s33**2*s35**2 + 512*s11**2*s14*s16*s22**2*s24*s26*s33**2*s35**2 + 512*s11*s15*s16**2*s22**2*s24*s26*s33**2*s35**2 - 512*s11**2*s12*s16*s22*s24**2*s26*s33**2*s35**2 -  \
    512*s11**2*s12*s14*s22**2*s25*s26*s33**2*s35**2 + 1024*s11*s12*s15*s16*s22**2*s25*s26*s33**2*s35**2 - 1024*s11*s14*s16**2*s22**2*s25*s26*s33**2*s35**2 + 512*s11*s12*s16**2*s22*s24*s25*s26*s33**2*s35**2 -  \
    512*s11*s12**2*s16*s22*s25**2*s26*s33**2*s35**2 + 1536*s11**2*s12*s13*s22**2*s26**2*s33**2*s35**2 + 512*s11*s12*s15**2*s22**2*s26**2*s33**2*s35**2 + 512*s11*s14*s15*s16*s22**2*s26**2*s33**2*s35**2 +  \
    512*s11*s13*s16**2*s22**2*s26**2*s33**2*s35**2 - 1024*s11**2*s12**2*s22*s23*s26**2*s33**2*s35**2 + 1024*s11*s12*s16**2*s22*s23*s26**2*s33**2*s35**2 - 512*s11**2*s12*s14*s22*s24*s26**2*s33**2*s35**2 -  \
    1536*s11*s12*s15*s16*s22*s24*s26**2*s33**2*s35**2 + 512*s11**2*s12**2*s24**2*s26**2*s33**2*s35**2 - 1024*s11*s12**2*s15*s22*s25*s26**2*s33**2*s35**2 + 1536*s11*s12*s14*s16*s22*s25*s26**2*s33**2*s35**2 -  \
    512*s11*s12**2*s16*s24*s25*s26**2*s33**2*s35**2 + 512*s11*s12**3*s25**2*s26**2*s33**2*s35**2 - 512*s11*s12*s14*s15*s22*s26**3*s33**2*s35**2 - 1024*s11*s12*s13*s16*s22*s26**3*s33**2*s35**2 - 512*s11*s12**2*s16*s23*s26**3*s33**2*s35**2 +  \
    1024*s11*s12**2*s15*s24*s26**3*s33**2*s35**2 - 512*s11*s12**2*s14*s25*s26**3*s33**2*s35**2 + 512*s11*s12**2*s13*s26**4*s33**2*s35**2 + 2048*s11**3*s14*s22**3*s23*s33*s34*s35**2 - 1024*s11**3*s13*s22**3*s24*s33*s34*s35**2 -  \
    512*s11**2*s15**2*s22**3*s24*s33*s34*s35**2 - 1024*s11**3*s12*s22**2*s23*s24*s33*s34*s35**2 - 512*s11**2*s16**2*s22**2*s23*s24*s33*s34*s35**2 + 512*s11**2*s15*s16*s22**2*s24**2*s33*s34*s35**2 +  \
    512*s11**2*s14*s15*s22**3*s25*s33*s34*s35**2 - 1536*s11**2*s13*s16*s22**3*s25*s33*s34*s35**2 + 1536*s11**2*s12*s16*s22**2*s23*s25*s33*s34*s35**2 + 512*s11**2*s12*s15*s22**2*s24*s25*s33*s34*s35**2 -  \
    512*s11**2*s14*s16*s22**2*s24*s25*s33*s34*s35**2 - 512*s11**2*s12*s14*s22**2*s25**2*s33*s34*s35**2 + 1536*s11**2*s13*s15*s22**3*s26*s33*s34*s35**2 - 1536*s11**2*s12*s15*s22**2*s23*s26*s33*s34*s35**2 -  \
    1536*s11**2*s14*s16*s22**2*s23*s26*s33*s34*s35**2 - 512*s11**2*s14*s15*s22**2*s24*s26*s33*s34*s35**2 + 1536*s11**2*s13*s16*s22**2*s24*s26*s33*s34*s35**2 + 512*s11*s15**2*s16*s22**2*s24*s26*s33*s34*s35**2 +  \
    1024*s11**2*s12*s16*s22*s23*s24*s26*s33*s34*s35**2 + 512*s11*s16**3*s22*s23*s24*s26*s33*s34*s35**2 - 512*s11*s15*s16**2*s22*s24**2*s26*s33*s34*s35**2 + 512*s11**2*s14**2*s22**2*s25*s26*s33*s34*s35**2 -  \
    512*s11*s14*s15*s16*s22**2*s25*s26*s33*s34*s35**2 + 1536*s11*s13*s16**2*s22**2*s25*s26*s33*s34*s35**2 - 1536*s11*s12*s16**2*s22*s23*s25*s26*s33*s34*s35**2 - 512*s11*s12*s15*s16*s22*s24*s25*s26*s33*s34*s35**2 +  \
    512*s11*s14*s16**2*s22*s24*s25*s26*s33*s34*s35**2 + 512*s11*s12*s14*s16*s22*s25**2*s26*s33*s34*s35**2 - 512*s11**2*s13*s14*s22**2*s26**2*s33*s34*s35**2 - 1536*s11*s13*s15*s16*s22**2*s26**2*s33*s34*s35**2 +  \
    2048*s11**2*s12*s14*s22*s23*s26**2*s33*s34*s35**2 + 1536*s11*s12*s15*s16*s22*s23*s26**2*s33*s34*s35**2 - 512*s11*s14*s16**2*s22*s23*s26**2*s33*s34*s35**2 - 1024*s11**2*s12*s13*s22*s24*s26**2*s33*s34*s35**2 -  \
    512*s11*s12*s15**2*s22*s24*s26**2*s33*s34*s35**2 + 512*s11*s14*s15*s16*s22*s24*s26**2*s33*s34*s35**2 - 512*s11*s13*s16**2*s22*s24*s26**2*s33*s34*s35**2 - 1024*s11**2*s12**2*s23*s24*s26**2*s33*s34*s35**2 -  \
    512*s11*s12*s16**2*s23*s24*s26**2*s33*s34*s35**2 + 512*s11*s12*s15*s16*s24**2*s26**2*s33*s34*s35**2 + 512*s11*s12*s14*s15*s22*s25*s26**2*s33*s34*s35**2 - 1536*s11*s12*s13*s16*s22*s25*s26**2*s33*s34*s35**2 -  \
    512*s11*s14**2*s16*s22*s25*s26**2*s33*s34*s35**2 + 1536*s11*s12**2*s16*s23*s25*s26**2*s33*s34*s35**2 + 512*s11*s12**2*s15*s24*s25*s26**2*s33*s34*s35**2 - 512*s11*s12*s14*s16*s24*s25*s26**2*s33*s34*s35**2 -  \
    512*s11*s12**2*s14*s25**2*s26**2*s33*s34*s35**2 + 1536*s11*s12*s13*s15*s22*s26**3*s33*s34*s35**2 + 512*s11*s13*s14*s16*s22*s26**3*s33*s34*s35**2 - 1536*s11*s12**2*s15*s23*s26**3*s33*s34*s35**2 +  \
    512*s11*s12*s14*s16*s23*s26**3*s33*s34*s35**2 - 512*s11*s12*s14*s15*s24*s26**3*s33*s34*s35**2 + 512*s11*s12*s13*s16*s24*s26**3*s33*s34*s35**2 + 512*s11*s12*s14**2*s25*s26**3*s33*s34*s35**2 -  \
    512*s11*s12*s13*s14*s26**4*s33*s34*s35**2 - 1024*s11**3*s13*s22**3*s23*s34**2*s35**2 + 512*s11**2*s15**2*s22**3*s23*s34**2*s35**2 + 1024*s11**3*s12*s22**2*s23**2*s34**2*s35**2 + 512*s11**2*s16**2*s22**2*s23**2*s34**2*s35**2 -  \
    512*s11**3*s14*s22**2*s23*s24*s34**2*s35**2 - 512*s11**2*s15*s16*s22**2*s23*s24*s34**2*s35**2 + 512*s11**3*s13*s22**2*s24**2*s34**2*s35**2 - 512*s11**2*s13*s15*s22**3*s25*s34**2*s35**2 - 512*s11**2*s12*s15*s22**2*s23*s25*s34**2*s35**2 -  \
    512*s11**2*s14*s16*s22**2*s23*s25*s34**2*s35**2 + 1024*s11**2*s13*s16*s22**2*s24*s25*s34**2*s35**2 + 512*s11**2*s12*s13*s22**2*s25**2*s34**2*s35**2 + 1024*s11**2*s14*s15*s22**2*s23*s26*s34**2*s35**2 -  \
    512*s11*s15**2*s16*s22**2*s23*s26*s34**2*s35**2 - 1024*s11**2*s12*s16*s22*s23**2*s26*s34**2*s35**2 - 512*s11*s16**3*s22*s23**2*s26*s34**2*s35**2 - 512*s11**2*s13*s15*s22**2*s24*s26*s34**2*s35**2 +  \
    512*s11**2*s14*s16*s22*s23*s24*s26*s34**2*s35**2 + 512*s11*s15*s16**2*s22*s23*s24*s26*s34**2*s35**2 - 512*s11**2*s13*s16*s22*s24**2*s26*s34**2*s35**2 - 512*s11**2*s13*s14*s22**2*s25*s26*s34**2*s35**2 +  \
    512*s11*s13*s15*s16*s22**2*s25*s26*s34**2*s35**2 + 512*s11*s12*s15*s16*s22*s23*s25*s26*s34**2*s35**2 + 512*s11*s14*s16**2*s22*s23*s25*s26*s34**2*s35**2 - 1024*s11*s13*s16**2*s22*s24*s25*s26*s34**2*s35**2 -  \
    512*s11*s12*s13*s16*s22*s25**2*s26*s34**2*s35**2 + 512*s11**2*s13**2*s22**2*s26**2*s34**2*s35**2 - 1024*s11**2*s12*s13*s22*s23*s26**2*s34**2*s35**2 + 512*s11*s12*s15**2*s22*s23*s26**2*s34**2*s35**2 -  \
    1024*s11*s14*s15*s16*s22*s23*s26**2*s34**2*s35**2 + 1024*s11*s13*s16**2*s22*s23*s26**2*s34**2*s35**2 + 1024*s11**2*s12**2*s23**2*s26**2*s34**2*s35**2 + 512*s11*s12*s16**2*s23**2*s26**2*s34**2*s35**2 +  \
    512*s11*s13*s15*s16*s22*s24*s26**2*s34**2*s35**2 - 512*s11**2*s12*s14*s23*s24*s26**2*s34**2*s35**2 - 512*s11*s12*s15*s16*s23*s24*s26**2*s34**2*s35**2 + 512*s11**2*s12*s13*s24**2*s26**2*s34**2*s35**2 -  \
    512*s11*s12*s13*s15*s22*s25*s26**2*s34**2*s35**2 + 512*s11*s13*s14*s16*s22*s25*s26**2*s34**2*s35**2 - 512*s11*s12**2*s15*s23*s25*s26**2*s34**2*s35**2 - 512*s11*s12*s14*s16*s23*s25*s26**2*s34**2*s35**2 +  \
    1024*s11*s12*s13*s16*s24*s25*s26**2*s34**2*s35**2 + 512*s11*s12**2*s13*s25**2*s26**2*s34**2*s35**2 - 512*s11*s13**2*s16*s22*s26**3*s34**2*s35**2 + 1024*s11*s12*s14*s15*s23*s26**3*s34**2*s35**2 -  \
    1024*s11*s12*s13*s16*s23*s26**3*s34**2*s35**2 - 512*s11*s12*s13*s15*s24*s26**3*s34**2*s35**2 - 512*s11*s12*s13*s14*s25*s26**3*s34**2*s35**2 + 512*s11*s12*s13**2*s26**4*s34**2*s35**2 - 1024*s11**2*s13*s15*s22**4*s33*s35**3 +  \
    1024*s11**2*s12*s15*s22**3*s23*s33*s35**3 - 1024*s11**2*s14*s16*s22**3*s23*s33*s35**3 + 512*s11**2*s14*s15*s22**3*s24*s33*s35**3 + 512*s11**2*s13*s16*s22**3*s24*s33*s35**3 + 512*s11**2*s12*s16*s22**2*s23*s24*s33*s35**3 -  \
    512*s11**2*s12*s15*s22**2*s24**2*s33*s35**3 + 1024*s11**2*s12*s13*s22**3*s25*s33*s35**3 - 512*s11**2*s14**2*s22**3*s25*s33*s35**3 - 1024*s11**2*s12**2*s22**2*s23*s25*s33*s35**3 + 512*s11**2*s12*s14*s22**2*s24*s25*s33*s35**3 +  \
    512*s11**2*s13*s14*s22**3*s26*s33*s35**3 + 1024*s11*s13*s15*s16*s22**3*s26*s33*s35**3 + 512*s11**2*s12*s14*s22**2*s23*s26*s33*s35**3 - 1024*s11*s12*s15*s16*s22**2*s23*s26*s33*s35**3 +  \
    1024*s11*s14*s16**2*s22**2*s23*s26*s33*s35**3 - 1024*s11**2*s12*s13*s22**2*s24*s26*s33*s35**3 - 512*s11*s14*s15*s16*s22**2*s24*s26*s33*s35**3 - 512*s11*s13*s16**2*s22**2*s24*s26*s33*s35**3 -  \
    512*s11*s12*s16**2*s22*s23*s24*s26*s33*s35**3 + 512*s11*s12*s15*s16*s22*s24**2*s26*s33*s35**3 - 1024*s11*s12*s13*s16*s22**2*s25*s26*s33*s35**3 + 512*s11*s14**2*s16*s22**2*s25*s26*s33*s35**3 +  \
    1024*s11*s12**2*s16*s22*s23*s25*s26*s33*s35**3 - 512*s11*s12*s14*s16*s22*s24*s25*s26*s33*s35**3 - 1024*s11*s12*s13*s15*s22**2*s26**2*s33*s35**3 - 512*s11*s13*s14*s16*s22**2*s26**2*s33*s35**3 +  \
    1024*s11*s12**2*s15*s22*s23*s26**2*s33*s35**3 - 1536*s11*s12*s14*s16*s22*s23*s26**2*s33*s35**3 + 512*s11*s12*s14*s15*s22*s24*s26**2*s33*s35**3 + 1536*s11*s12*s13*s16*s22*s24*s26**2*s33*s35**3 +  \
    512*s11*s12**2*s16*s23*s24*s26**2*s33*s35**3 - 512*s11*s12**2*s15*s24**2*s26**2*s33*s35**3 + 1024*s11*s12**2*s13*s22*s25*s26**2*s33*s35**3 - 512*s11*s12*s14**2*s22*s25*s26**2*s33*s35**3 - 1024*s11*s12**3*s23*s25*s26**2*s33*s35**3 +  \
    512*s11*s12**2*s14*s24*s25*s26**2*s33*s35**3 + 512*s11*s12*s13*s14*s22*s26**3*s33*s35**3 + 512*s11*s12**2*s14*s23*s26**3*s33*s35**3 - 1024*s11*s12**2*s13*s24*s26**3*s33*s35**3 - 1024*s11**2*s14*s15*s22**3*s23*s34*s35**3 +  \
    1024*s11**2*s13*s16*s22**3*s23*s34*s35**3 - 1024*s11**2*s12*s16*s22**2*s23**2*s34*s35**3 + 512*s11**2*s13*s15*s22**3*s24*s34*s35**3 + 512*s11**2*s12*s15*s22**2*s23*s24*s34*s35**3 + 512*s11**2*s14*s16*s22**2*s23*s24*s34*s35**3 -  \
    512*s11**2*s13*s16*s22**2*s24**2*s34*s35**3 + 512*s11**2*s13*s14*s22**3*s25*s34*s35**3 + 512*s11**2*s12*s14*s22**2*s23*s25*s34*s35**3 - 1024*s11**2*s12*s13*s22**2*s24*s25*s34*s35**3 - 1024*s11**2*s13**2*s22**3*s26*s34*s35**3 +  \
    1024*s11**2*s12*s13*s22**2*s23*s26*s34*s35**3 - 512*s11**2*s14**2*s22**2*s23*s26*s34*s35**3 + 1024*s11*s14*s15*s16*s22**2*s23*s26*s34*s35**3 - 1024*s11*s13*s16**2*s22**2*s23*s26*s34*s35**3 +  \
    1024*s11*s12*s16**2*s22*s23**2*s26*s34*s35**3 + 512*s11**2*s13*s14*s22**2*s24*s26*s34*s35**3 - 512*s11*s13*s15*s16*s22**2*s24*s26*s34*s35**3 - 512*s11*s12*s15*s16*s22*s23*s24*s26*s34*s35**3 -  \
    512*s11*s14*s16**2*s22*s23*s24*s26*s34*s35**3 + 512*s11*s13*s16**2*s22*s24**2*s26*s34*s35**3 - 512*s11*s13*s14*s16*s22**2*s25*s26*s34*s35**3 - 512*s11*s12*s14*s16*s22*s23*s25*s26*s34*s35**3;
    v4_3= \
    1024*s11*s12*s13*s16*s22*s24*s25*s26*s34*s35**3 + 1024*s11*s13**2*s16*s22**2*s26**2*s34*s35**3 - 1024*s11*s12*s14*s15*s22*s23*s26**2*s34*s35**3 + 512*s11*s14**2*s16*s22*s23*s26**2*s34*s35**3 -  \
    1024*s11*s12**2*s16*s23**2*s26**2*s34*s35**3 + 512*s11*s12*s13*s15*s22*s24*s26**2*s34*s35**3 - 512*s11*s13*s14*s16*s22*s24*s26**2*s34*s35**3 + 512*s11*s12**2*s15*s23*s24*s26**2*s34*s35**3 +  \
    512*s11*s12*s14*s16*s23*s24*s26**2*s34*s35**3 - 512*s11*s12*s13*s16*s24**2*s26**2*s34*s35**3 + 512*s11*s12*s13*s14*s22*s25*s26**2*s34*s35**3 + 512*s11*s12**2*s14*s23*s25*s26**2*s34*s35**3 -  \
    1024*s11*s12**2*s13*s24*s25*s26**2*s34*s35**3 - 1024*s11*s12*s13**2*s22*s26**3*s34*s35**3 + 1024*s11*s12**2*s13*s23*s26**3*s34*s35**3 - 512*s11*s12*s14**2*s23*s26**3*s34*s35**3 + 512*s11*s12*s13*s14*s24*s26**3*s34*s35**3 +  \
    512*s11**2*s13**2*s22**4*s35**4 - 1024*s11**2*s12*s13*s22**3*s23*s35**4 + 512*s11**2*s14**2*s22**3*s23*s35**4 + 512*s11**2*s12**2*s22**2*s23**2*s35**4 - 512*s11**2*s13*s14*s22**3*s24*s35**4 - 512*s11**2*s12*s14*s22**2*s23*s24*s35**4 +  \
    512*s11**2*s12*s13*s22**2*s24**2*s35**4 - 512*s11*s13**2*s16*s22**3*s26*s35**4 + 1024*s11*s12*s13*s16*s22**2*s23*s26*s35**4 - 512*s11*s14**2*s16*s22**2*s23*s26*s35**4 - 512*s11*s12**2*s16*s22*s23**2*s26*s35**4 +  \
    512*s11*s13*s14*s16*s22**2*s24*s26*s35**4 + 512*s11*s12*s14*s16*s22*s23*s24*s26*s35**4 - 512*s11*s12*s13*s16*s22*s24**2*s26*s35**4 + 512*s11*s12*s13**2*s22**2*s26**2*s35**4 - 1024*s11*s12**2*s13*s22*s23*s26**2*s35**4 +  \
    512*s11*s12*s14**2*s22*s23*s26**2*s35**4 + 512*s11*s12**3*s23**2*s26**2*s35**4 - 512*s11*s12*s13*s14*s22*s24*s26**2*s35**4 - 512*s11*s12**2*s14*s23*s24*s26**2*s35**4 + 512*s11*s12**2*s13*s24**2*s26**2*s35**4 +  \
    1024*s11**3*s16*s22**3*s23*s33**3*s36 + 512*s11**3*s15*s22**3*s24*s33**3*s36 - 512*s11**3*s16*s22**2*s24**2*s33**3*s36 + 1536*s11**3*s14*s22**3*s25*s33**3*s36 - 1536*s11**2*s15*s16*s22**3*s25*s33**3*s36 -  \
    2048*s11**3*s12*s22**2*s24*s25*s33**3*s36 + 1024*s11**2*s16**2*s22**2*s24*s25*s33**3*s36 + 1536*s11**2*s12*s16*s22**2*s25**2*s33**3*s36 - 512*s11*s16**3*s22**2*s25**2*s33**3*s36 + 1024*s11**3*s13*s22**3*s26*s33**3*s36 -  \
    512*s11**2*s15**2*s22**3*s26*s33**3*s36 - 2048*s11**3*s12*s22**2*s23*s26*s33**3*s36 - 1024*s11**2*s16**2*s22**2*s23*s26*s33**3*s36 - 1536*s11**3*s14*s22**2*s24*s26*s33**3*s36 + 512*s11**2*s15*s16*s22**2*s24*s26*s33**3*s36 +  \
    2048*s11**3*s12*s22*s24**2*s26*s33**3*s36 + 2560*s11**2*s12*s15*s22**2*s25*s26*s33**3*s36 - 512*s11**2*s14*s16*s22**2*s25*s26*s33**3*s36 + 512*s11*s15*s16**2*s22**2*s25*s26*s33**3*s36 -  \
    2048*s11**2*s12*s16*s22*s24*s25*s26*s33**3*s36 - 2048*s11**2*s12**2*s22*s25**2*s26*s33**3*s36 + 1024*s11*s12*s16**2*s22*s25**2*s26*s33**3*s36 + 1024*s11**2*s14*s15*s22**2*s26**2*s33**3*s36 -  \
    1024*s11**2*s13*s16*s22**2*s26**2*s33**3*s36 + 3072*s11**2*s12*s16*s22*s23*s26**2*s33**3*s36 - 1536*s11**2*s12*s15*s22*s24*s26**2*s33**3*s36 + 512*s11**2*s14*s16*s22*s24*s26**2*s33**3*s36 -  \
    512*s11**2*s12*s16*s24**2*s26**2*s33**3*s36 - 512*s11**2*s12*s14*s22*s25*s26**2*s33**3*s36 - 1024*s11*s12*s15*s16*s22*s25*s26**2*s33**3*s36 + 2048*s11**2*s12**2*s24*s25*s26**2*s33**3*s36 -  \
    512*s11*s12**2*s16*s25**2*s26**2*s33**3*s36 + 1024*s11**2*s12*s13*s22*s26**3*s33**3*s36 - 512*s11**2*s14**2*s22*s26**3*s33**3*s36 - 2048*s11**2*s12**2*s23*s26**3*s33**3*s36 + 512*s11**2*s12*s14*s24*s26**3*s33**3*s36 +  \
    512*s11*s12**2*s15*s25*s26**3*s33**3*s36 - 1024*s11**3*s15*s22**3*s23*s33**2*s34*s36 - 512*s11**3*s16*s22**2*s23*s24*s33**2*s34*s36 - 512*s11**3*s15*s22**2*s24**2*s33**2*s34*s36 + 512*s11**3*s16*s22*s24**3*s33**2*s34*s36 -  \
    2048*s11**3*s13*s22**3*s25*s33**2*s34*s36 + 1024*s11**2*s15**2*s22**3*s25*s33**2*s34*s36 + 3072*s11**3*s12*s22**2*s23*s25*s33**2*s34*s36 - 512*s11**2*s16**2*s22**2*s23*s25*s33**2*s34*s36 -  \
    512*s11**3*s14*s22**2*s24*s25*s33**2*s34*s36 + 1024*s11**2*s15*s16*s22**2*s24*s25*s33**2*s34*s36 + 1024*s11**3*s12*s22*s24**2*s25*s33**2*s34*s36 - 1024*s11**2*s16**2*s22*s24**2*s25*s33**2*s34*s36 -  \
    2048*s11**2*s12*s15*s22**2*s25**2*s33**2*s34*s36 - 512*s11**2*s14*s16*s22**2*s25**2*s33**2*s34*s36 + 512*s11*s15*s16**2*s22**2*s25**2*s33**2*s34*s36 - 512*s11**2*s12*s16*s22*s24*s25**2*s33**2*s34*s36 +  \
    512*s11*s16**3*s22*s24*s25**2*s33**2*s34*s36 + 1024*s11**2*s12**2*s22*s25**3*s33**2*s34*s36 - 512*s11*s12*s16**2*s22*s25**3*s33**2*s34*s36 + 2560*s11**3*s14*s22**2*s23*s26*s33**2*s34*s36 +  \
    512*s11**2*s15*s16*s22**2*s23*s26*s33**2*s34*s36 + 512*s11**2*s15**2*s22**2*s24*s26*s33**2*s34*s36 - 2048*s11**3*s12*s22*s23*s24*s26*s33**2*s34*s36 + 1024*s11**2*s16**2*s22*s23*s24*s26*s33**2*s34*s36 +  \
    512*s11**3*s14*s22*s24**2*s26*s33**2*s34*s36 - 512*s11**2*s15*s16*s22*s24**2*s26*s33**2*s34*s36 - 1024*s11**3*s12*s24**3*s26*s33**2*s34*s36 - 2048*s11**2*s14*s15*s22**2*s25*s26*s33**2*s34*s36 +  \
    1536*s11**2*s13*s16*s22**2*s25*s26*s33**2*s34*s36 - 512*s11*s15**2*s16*s22**2*s25*s26*s33**2*s34*s36 - 1024*s11**2*s12*s16*s22*s23*s25*s26*s33**2*s34*s36 + 1024*s11**2*s14*s16*s22*s24*s25*s26*s33**2*s34*s36 -  \
    512*s11*s15*s16**2*s22*s24*s25*s26*s33**2*s34*s36 + 1536*s11**2*s12*s16*s24**2*s25*s26*s33**2*s34*s36 + 2560*s11**2*s12*s14*s22*s25**2*s26*s33**2*s34*s36 - 512*s11*s14*s16**2*s22*s25**2*s26*s33**2*s34*s36 -  \
    1024*s11**2*s12**2*s24*s25**2*s26*s33**2*s34*s36 - 512*s11*s12*s16**2*s24*s25**2*s26*s33**2*s34*s36 + 512*s11*s12**2*s16*s25**3*s26*s33**2*s34*s36 - 512*s11**2*s13*s15*s22**2*s26**2*s33**2*s34*s36 -  \
    2048*s11**2*s14*s16*s22*s23*s26**2*s33**2*s34*s36 - 512*s11**2*s14*s15*s22*s24*s26**2*s33**2*s34*s36 + 512*s11**2*s13*s16*s22*s24*s26**2*s33**2*s34*s36 - 512*s11**2*s12*s16*s23*s24*s26**2*s33**2*s34*s36 +  \
    1024*s11**2*s12*s15*s24**2*s26**2*s33**2*s34*s36 - 1024*s11**2*s12*s13*s22*s25*s26**2*s33**2*s34*s36 + 1024*s11**2*s14**2*s22*s25*s26**2*s33**2*s34*s36 + 512*s11*s12*s15**2*s22*s25*s26**2*s33**2*s34*s36 +  \
    512*s11*s14*s15*s16*s22*s25*s26**2*s33**2*s34*s36 + 1024*s11**2*s12**2*s23*s25*s26**2*s33**2*s34*s36 - 2560*s11**2*s12*s14*s24*s25*s26**2*s33**2*s34*s36 + 512*s11*s12*s15*s16*s24*s25*s26**2*s33**2*s34*s36 -  \
    512*s11*s12**2*s15*s25**2*s26**2*s33**2*s34*s36 + 512*s11*s12*s14*s16*s25**2*s26**2*s33**2*s34*s36 + 512*s11**2*s13*s14*s22*s26**3*s33**2*s34*s36 + 1536*s11**2*s12*s14*s23*s26**3*s33**2*s34*s36 -  \
    1024*s11**2*s12*s13*s24*s26**3*s33**2*s34*s36 - 512*s11*s12*s14*s15*s25*s26**3*s33**2*s34*s36 + 1024*s11**3*s16*s22**2*s23**2*s33*s34**2*s36 + 1536*s11**3*s15*s22**2*s23*s24*s33*s34**2*s36 -  \
    1024*s11**3*s16*s22*s23*s24**2*s33*s34**2*s36 - 512*s11**3*s14*s22**2*s23*s25*s33*s34**2*s36 - 1536*s11**2*s15*s16*s22**2*s23*s25*s33*s34**2*s36 + 1024*s11**3*s13*s22**2*s24*s25*s33*s34**2*s36 -  \
    1024*s11**2*s15**2*s22**2*s24*s25*s33*s34**2*s36 - 2048*s11**3*s12*s22*s23*s24*s25*s33*s34**2*s36 + 1536*s11**2*s16**2*s22*s23*s24*s25*s33*s34**2*s36 + 512*s11**2*s15*s16*s22*s24**2*s25*s33*s34**2*s36 +  \
    1024*s11**2*s14*s15*s22**2*s25**2*s33*s34**2*s36 + 512*s11**2*s13*s16*s22**2*s25**2*s33*s34**2*s36 + 1024*s11**2*s12*s16*s22*s23*s25**2*s33*s34**2*s36 - 512*s11*s16**3*s22*s23*s25**2*s33*s34**2*s36 +  \
    1024*s11**2*s12*s15*s22*s24*s25**2*s33*s34**2*s36 - 512*s11**2*s14*s16*s22*s24*s25**2*s33*s34**2*s36 - 512*s11*s15*s16**2*s22*s24*s25**2*s33*s34**2*s36 - 1024*s11**2*s12*s14*s22*s25**3*s33*s34**2*s36 +  \
    512*s11*s14*s16**2*s22*s25**3*s33*s34**2*s36 - 1024*s11**3*s13*s22**2*s23*s26*s33*s34**2*s36 - 512*s11**2*s15**2*s22**2*s23*s26*s33*s34**2*s36 - 1024*s11**2*s16**2*s22*s23**2*s26*s33*s34**2*s36 -  \
    1024*s11**3*s14*s22*s23*s24*s26*s33*s34**2*s36 + 2048*s11**3*s12*s23*s24**2*s26*s33*s34**2*s36 + 1536*s11**2*s13*s15*s22**2*s25*s26*s33*s34**2*s36 + 1024*s11**2*s12*s15*s22*s23*s25*s26*s33*s34**2*s36 +  \
    1024*s11**2*s14*s16*s22*s23*s25*s26*s33*s34**2*s36 + 512*s11*s15*s16**2*s22*s23*s25*s26*s33*s34**2*s36 + 512*s11**2*s14*s15*s22*s24*s25*s26*s33*s34**2*s36 - 2048*s11**2*s13*s16*s22*s24*s25*s26*s33*s34**2*s36 +  \
    512*s11*s15**2*s16*s22*s24*s25*s26*s33*s34**2*s36 - 2048*s11**2*s12*s16*s23*s24*s25*s26*s33*s34**2*s36 - 1024*s11**2*s12*s15*s24**2*s25*s26*s33*s34**2*s36 - 2048*s11**2*s12*s13*s22*s25**2*s26*s33*s34**2*s36 -  \
    512*s11**2*s14**2*s22*s25**2*s26*s33*s34**2*s36 - 512*s11*s14*s15*s16*s22*s25**2*s26*s33*s34**2*s36 + 512*s11*s13*s16**2*s22*s25**2*s26*s33*s34**2*s36 + 512*s11*s12*s16**2*s23*s25**2*s26*s33*s34**2*s36 +  \
    1024*s11**2*s12*s14*s24*s25**2*s26*s33*s34**2*s36 + 512*s11*s12*s15*s16*s24*s25**2*s26*s33*s34**2*s36 - 512*s11*s12*s14*s16*s25**3*s26*s33*s34**2*s36 + 512*s11**2*s14*s15*s22*s23*s26**2*s33*s34**2*s36 +  \
    1024*s11**2*s13*s16*s22*s23*s26**2*s33*s34**2*s36 + 1024*s11**2*s12*s16*s23**2*s26**2*s33*s34**2*s36 - 512*s11**2*s12*s15*s23*s24*s26**2*s33*s34**2*s36 - 1024*s11**2*s13*s14*s22*s25*s26**2*s33*s34**2*s36 -  \
    512*s11*s13*s15*s16*s22*s25*s26**2*s33*s34**2*s36 - 512*s11**2*s12*s14*s23*s25*s26**2*s33*s34**2*s36 - 512*s11*s12*s15*s16*s23*s25*s26**2*s33*s34**2*s36 + 3072*s11**2*s12*s13*s24*s25*s26**2*s33*s34**2*s36 -  \
    512*s11*s12*s15**2*s24*s25*s26**2*s33*s34**2*s36 + 512*s11*s12*s14*s15*s25**2*s26**2*s33*s34**2*s36 - 512*s11*s12*s13*s16*s25**2*s26**2*s33*s34**2*s36 - 1024*s11**2*s12*s13*s23*s26**3*s33*s34**2*s36 +  \
    512*s11*s12*s13*s15*s25*s26**3*s33*s34**2*s36 - 1024*s11**3*s15*s22**2*s23**2*s34**3*s36 + 512*s11**3*s16*s22*s23**2*s24*s34**3*s36 + 1024*s11**2*s15**2*s22**2*s23*s25*s34**3*s36 + 1024*s11**3*s12*s22*s23**2*s25*s34**3*s36 -  \
    512*s11**2*s16**2*s22*s23**2*s25*s34**3*s36 - 512*s11**2*s15*s16*s22*s23*s24*s25*s34**3*s36 - 1024*s11**2*s13*s15*s22**2*s25**2*s34**3*s36 - 1024*s11**2*s12*s15*s22*s23*s25**2*s34**3*s36 +  \
    512*s11*s15*s16**2*s22*s23*s25**2*s34**3*s36 + 512*s11**2*s13*s16*s22*s24*s25**2*s34**3*s36 + 1024*s11**2*s12*s13*s22*s25**3*s34**3*s36 - 512*s11*s13*s16**2*s22*s25**3*s34**3*s36 + 512*s11**3*s14*s22*s23**2*s26*s34**3*s36 +  \
    512*s11**2*s15*s16*s22*s23**2*s26*s34**3*s36 - 1024*s11**3*s12*s23**2*s24*s26*s34**3*s36 - 512*s11**2*s14*s15*s22*s23*s25*s26*s34**3*s36 - 512*s11*s15**2*s16*s22*s23*s25*s26*s34**3*s36 +  \
    512*s11**2*s12*s16*s23**2*s25*s26*s34**3*s36 + 1024*s11**2*s12*s15*s23*s24*s25*s26*s34**3*s36 + 512*s11**2*s13*s14*s22*s25**2*s26*s34**3*s36 + 512*s11*s13*s15*s16*s22*s25**2*s26*s34**3*s36 -  \
    512*s11*s12*s15*s16*s23*s25**2*s26*s34**3*s36 - 1024*s11**2*s12*s13*s24*s25**2*s26*s34**3*s36 + 512*s11*s12*s13*s16*s25**3*s26*s34**3*s36 - 512*s11**2*s12*s15*s23**2*s26**2*s34**3*s36 +  \
    512*s11*s12*s15**2*s23*s25*s26**2*s34**3*s36 - 512*s11*s12*s13*s15*s25**2*s26**2*s34**3*s36 - 2048*s11**3*s14*s22**3*s23*s33**2*s35*s36 - 1024*s11**3*s13*s22**3*s24*s33**2*s35*s36 - 512*s11**2*s15**2*s22**3*s24*s33**2*s35*s36 +  \
    3072*s11**3*s12*s22**2*s23*s24*s33**2*s35*s36 - 512*s11**2*s16**2*s22**2*s23*s24*s33**2*s35*s36 + 1024*s11**3*s14*s22**2*s24**2*s33**2*s35*s36 + 512*s11**2*s15*s16*s22**2*s24**2*s33**2*s35*s36 -  \
    1024*s11**3*s12*s22*s24**3*s33**2*s35*s36 - 512*s11**2*s14*s15*s22**3*s25*s33**2*s35*s36 + 2560*s11**2*s13*s16*s22**3*s25*s33**2*s35*s36 + 512*s11*s15**2*s16*s22**3*s25*s33**2*s35*s36 -  \
    2560*s11**2*s12*s16*s22**2*s23*s25*s33**2*s35*s36 + 512*s11*s16**3*s22**2*s23*s25*s33**2*s35*s36 + 1536*s11**2*s12*s15*s22**2*s24*s25*s33**2*s35*s36 - 2048*s11**2*s14*s16*s22**2*s24*s25*s33**2*s35*s36 -  \
    512*s11*s15*s16**2*s22**2*s24*s25*s33**2*s35*s36 + 1536*s11**2*s12*s16*s22*s24**2*s25*s33**2*s35*s36 + 512*s11**2*s12*s14*s22**2*s25**2*s33**2*s35*s36 - 1024*s11*s12*s15*s16*s22**2*s25**2*s33**2*s35*s36 +  \
    1024*s11*s14*s16**2*s22**2*s25**2*s33**2*s35*s36 - 1024*s11**2*s12**2*s22*s24*s25**2*s33**2*s35*s36 - 512*s11*s12*s16**2*s22*s24*s25**2*s33**2*s35*s36 + 512*s11*s12**2*s16*s22*s25**3*s33**2*s35*s36 -  \
    512*s11**2*s13*s15*s22**3*s26*s33**2*s35*s36 + 512*s11*s15**3*s22**3*s26*s33**2*s35*s36 + 512*s11**2*s12*s15*s22**2*s23*s26*s33**2*s35*s36 + 1536*s11**2*s14*s16*s22**2*s23*s26*s33**2*s35*s36;
    v4_4= \
    512*s11*s15*s16**2*s22**2*s23*s26*s33**2*s35*s36 + 1024*s11**2*s14*s15*s22**2*s24*s26*s33**2*s35*s36 + 512*s11**2*s13*s16*s22**2*s24*s26*s33**2*s35*s36 - 512*s11*s15**2*s16*s22**2*s24*s26*s33**2*s35*s36 -  \
    1024*s11**2*s12*s16*s22*s23*s24*s26*s33**2*s35*s36 - 1536*s11**2*s12*s15*s22*s24**2*s26*s33**2*s35*s36 - 512*s11**2*s14*s16*s22*s24**2*s26*s33**2*s35*s36 + 512*s11**2*s12*s16*s24**3*s26*s33**2*s35*s36 -  \
    2048*s11**2*s12*s13*s22**2*s25*s26*s33**2*s35*s36 - 512*s11**2*s14**2*s22**2*s25*s26*s33**2*s35*s36 - 2048*s11*s12*s15**2*s22**2*s25*s26*s33**2*s35*s36 + 1024*s11*s14*s15*s16*s22**2*s25*s26*s33**2*s35*s36 -  \
    2048*s11*s13*s16**2*s22**2*s25*s26*s33**2*s35*s36 + 2048*s11**2*s12**2*s22*s23*s25*s26*s33**2*s35*s36 + 2048*s11**2*s12*s14*s22*s24*s25*s26*s33**2*s35*s36 + 1024*s11*s12*s15*s16*s22*s24*s25*s26*s33**2*s35*s36 +  \
    512*s11*s14*s16**2*s22*s24*s25*s26*s33**2*s35*s36 - 1024*s11**2*s12**2*s24**2*s25*s26*s33**2*s35*s36 - 512*s11*s12*s16**2*s24**2*s25*s26*s33**2*s35*s36 + 2560*s11*s12**2*s15*s22*s25**2*s26*s33**2*s35*s36 -  \
    3072*s11*s12*s14*s16*s22*s25**2*s26*s33**2*s35*s36 + 1536*s11*s12**2*s16*s24*s25**2*s26*s33**2*s35*s36 - 1024*s11*s12**3*s25**3*s26*s33**2*s35*s36 - 512*s11**2*s13*s14*s22**2*s26**2*s33**2*s35*s36 -  \
    1024*s11*s14*s15**2*s22**2*s26**2*s33**2*s35*s36 + 1024*s11*s13*s15*s16*s22**2*s26**2*s33**2*s35*s36 - 1024*s11**2*s12*s14*s22*s23*s26**2*s33**2*s35*s36 - 2048*s11*s12*s15*s16*s22*s23*s26**2*s33**2*s35*s36 +  \
    512*s11**2*s14**2*s22*s24*s26**2*s33**2*s35*s36 + 1536*s11*s12*s15**2*s22*s24*s26**2*s33**2*s35*s36 - 512*s11*s14*s15*s16*s22*s24*s26**2*s33**2*s35*s36 + 1024*s11**2*s12**2*s23*s24*s26**2*s33**2*s35*s36 -  \
    512*s11**2*s12*s14*s24**2*s26**2*s33**2*s35*s36 + 512*s11*s12*s15*s16*s24**2*s26**2*s33**2*s35*s36 + 1024*s11*s12*s14*s15*s22*s25*s26**2*s33**2*s35*s36 + 3072*s11*s12*s13*s16*s22*s25*s26**2*s33**2*s35*s36 -  \
    512*s11*s14**2*s16*s22*s25*s26**2*s33**2*s35*s36 - 512*s11*s12**2*s16*s23*s25*s26**2*s33**2*s35*s36 - 2560*s11*s12**2*s15*s24*s25*s26**2*s33**2*s35*s36 + 512*s11*s12*s14*s16*s24*s25*s26**2*s33**2*s35*s36 +  \
    1024*s11*s12**2*s14*s25**2*s26**2*s33**2*s35*s36 - 1024*s11*s12*s13*s15*s22*s26**3*s33**2*s35*s36 + 512*s11*s14**2*s15*s22*s26**3*s33**2*s35*s36 + 1536*s11*s12**2*s15*s23*s26**3*s33**2*s35*s36 -  \
    512*s11*s12*s14*s15*s24*s26**3*s33**2*s35*s36 - 1024*s11*s12**2*s13*s25*s26**3*s33**2*s35*s36 + 4096*s11**3*s13*s22**3*s23*s33*s34*s35*s36 - 4096*s11**3*s12*s22**2*s23**2*s33*s34*s35*s36 -  \
    2048*s11**3*s14*s22**2*s23*s24*s33*s34*s35*s36 + 512*s11**2*s15**2*s22**2*s24**2*s33*s34*s35*s36 + 2048*s11**3*s12*s22*s23*s24**2*s33*s34*s35*s36 + 512*s11**2*s16**2*s22*s23*s24**2*s33*s34*s35*s36 -  \
    512*s11**2*s15*s16*s22*s24**3*s33*s34*s35*s36 - 2048*s11**2*s13*s15*s22**3*s25*s33*s34*s35*s36 + 2048*s11**2*s12*s15*s22**2*s23*s25*s33*s34*s35*s36 + 2048*s11**2*s14*s16*s22**2*s23*s25*s33*s34*s35*s36 +  \
    1024*s11**2*s14*s15*s22**2*s24*s25*s33*s34*s35*s36 - 512*s11*s15**2*s16*s22**2*s24*s25*s33*s34*s35*s36 - 2048*s11**2*s12*s16*s22*s23*s24*s25*s33*s34*s35*s36 - 512*s11*s16**3*s22*s23*s24*s25*s33*s34*s35*s36 -  \
    2048*s11**2*s12*s15*s22*s24**2*s25*s33*s34*s35*s36 + 512*s11**2*s14*s16*s22*s24**2*s25*s33*s34*s35*s36 + 512*s11*s15*s16**2*s22*s24**2*s25*s33*s34*s35*s36 + 2048*s11**2*s12*s13*s22**2*s25**2*s33*s34*s35*s36 -  \
    1536*s11**2*s14**2*s22**2*s25**2*s33*s34*s35*s36 + 512*s11*s14*s15*s16*s22**2*s25**2*s33*s34*s35*s36 - 1536*s11*s13*s16**2*s22**2*s25**2*s33*s34*s35*s36 - 2048*s11**2*s12**2*s22*s23*s25**2*s33*s34*s35*s36 +  \
    1536*s11*s12*s16**2*s22*s23*s25**2*s33*s34*s35*s36 + 2048*s11**2*s12*s14*s22*s24*s25**2*s33*s34*s35*s36 + 512*s11*s12*s15*s16*s22*s24*s25**2*s33*s34*s35*s36 - 512*s11*s14*s16**2*s22*s24*s25**2*s33*s34*s35*s36 -  \
    512*s11*s12*s14*s16*s22*s25**3*s33*s34*s35*s36 - 4096*s11**2*s13*s16*s22**2*s23*s26*s33*s34*s35*s36 + 4096*s11**2*s12*s16*s22*s23**2*s26*s33*s34*s35*s36 - 512*s11*s15**3*s22**2*s24*s26*s33*s34*s35*s36 +  \
    1024*s11**2*s14*s16*s22*s23*s24*s26*s33*s34*s35*s36 - 512*s11*s15*s16**2*s22*s23*s24*s26*s33*s34*s35*s36 - 512*s11**2*s14*s15*s22*s24**2*s26*s33*s34*s35*s36 + 512*s11*s15**2*s16*s22*s24**2*s26*s33*s34*s35*s36 -  \
    2048*s11**2*s12*s16*s23*s24**2*s26*s33*s34*s35*s36 + 1024*s11**2*s12*s15*s24**3*s26*s33*s34*s35*s36 + 2048*s11**2*s13*s14*s22**2*s25*s26*s33*s34*s35*s36 + 512*s11*s14*s15**2*s22**2*s25*s26*s33*s34*s35*s36 +  \
    1024*s11*s13*s15*s16*s22**2*s25*s26*s33*s34*s35*s36 - 4096*s11**2*s12*s14*s22*s23*s25*s26*s33*s34*s35*s36 - 1024*s11*s12*s15*s16*s22*s23*s25*s26*s33*s34*s35*s36 - 512*s11*s14*s16**2*s22*s23*s25*s26*s33*s34*s35*s36 -  \
    2048*s11**2*s12*s13*s22*s24*s25*s26*s33*s34*s35*s36 + 512*s11**2*s14**2*s22*s24*s25*s26*s33*s34*s35*s36 + 1536*s11*s12*s15**2*s22*s24*s25*s26*s33*s34*s35*s36 - 2048*s11*s14*s15*s16*s22*s24*s25*s26*s33*s34*s35*s36 +  \
    1536*s11*s13*s16**2*s22*s24*s25*s26*s33*s34*s35*s36 + 4096*s11**2*s12**2*s23*s24*s25*s26*s33*s34*s35*s36 + 1024*s11*s12*s16**2*s23*s24*s25*s26*s33*s34*s35*s36 - 1024*s11**2*s12*s14*s24**2*s25*s26*s33*s34*s35*s36 -  \
    1536*s11*s12*s14*s15*s22*s25**2*s26*s33*s34*s35*s36 + 2048*s11*s12*s13*s16*s22*s25**2*s26*s33*s34*s35*s36 + 1536*s11*s14**2*s16*s22*s25**2*s26*s33*s34*s35*s36 - 2048*s11*s12**2*s16*s23*s25**2*s26*s33*s34*s35*s36 -  \
    1024*s11*s12**2*s15*s24*s25**2*s26*s33*s34*s35*s36 + 1024*s11*s12**2*s14*s25**3*s26*s33*s34*s35*s36 + 512*s11*s13*s15**2*s22**2*s26**2*s33*s34*s35*s36 + 4096*s11**2*s12*s13*s22*s23*s26**2*s33*s34*s35*s36 -  \
    1536*s11**2*s14**2*s22*s23*s26**2*s33*s34*s35*s36 - 512*s11*s12*s15**2*s22*s23*s26**2*s33*s34*s35*s36 + 1536*s11*s14*s15*s16*s22*s23*s26**2*s33*s34*s35*s36 - 4096*s11**2*s12**2*s23**2*s26**2*s33*s34*s35*s36 +  \
    512*s11*s14*s15**2*s22*s24*s26**2*s33*s34*s35*s36 - 512*s11*s13*s15*s16*s22*s24*s26**2*s33*s34*s35*s36 + 2048*s11**2*s12*s14*s23*s24*s26**2*s33*s34*s35*s36 - 1024*s11*s12*s15**2*s24**2*s26**2*s33*s34*s35*s36 -  \
    2048*s11*s12*s13*s15*s22*s25*s26**2*s33*s34*s35*s36 - 512*s11*s14**2*s15*s22*s25*s26**2*s33*s34*s35*s36 - 512*s11*s13*s14*s16*s22*s25*s26**2*s33*s34*s35*s36 + 2048*s11*s12**2*s15*s23*s25*s26**2*s33*s34*s35*s36 +  \
    2048*s11*s12*s14*s15*s24*s25*s26**2*s33*s34*s35*s36 - 2048*s11*s12*s13*s16*s24*s25*s26**2*s33*s34*s35*s36 - 1024*s11*s12*s14**2*s25**2*s26**2*s33*s34*s35*s36 - 512*s11*s13*s14*s15*s22*s26**3*s33*s34*s35*s36 -  \
    1024*s11*s12*s14*s15*s23*s26**3*s33*s34*s35*s36 + 1024*s11*s12*s13*s15*s24*s26**3*s33*s34*s35*s36 + 1024*s11*s12*s13*s14*s25*s26**3*s33*s34*s35*s36 + 2048*s11**3*s14*s22**2*s23**2*s34**2*s35*s36 -  \
    1024*s11**3*s13*s22**2*s23*s24*s34**2*s35*s36 - 512*s11**2*s15**2*s22**2*s23*s24*s34**2*s35*s36 - 1024*s11**3*s12*s22*s23**2*s24*s34**2*s35*s36 - 512*s11**2*s16**2*s22*s23**2*s24*s34**2*s35*s36 +  \
    512*s11**2*s15*s16*s22*s23*s24**2*s34**2*s35*s36 - 2560*s11**2*s14*s15*s22**2*s23*s25*s34**2*s35*s36 + 512*s11**2*s13*s16*s22**2*s23*s25*s34**2*s35*s36 + 512*s11*s15**2*s16*s22**2*s23*s25*s34**2*s35*s36 -  \
    512*s11**2*s12*s16*s22*s23**2*s25*s34**2*s35*s36 + 512*s11*s16**3*s22*s23**2*s25*s34**2*s35*s36 + 1536*s11**2*s13*s15*s22**2*s24*s25*s34**2*s35*s36 + 2048*s11**2*s12*s15*s22*s23*s24*s25*s34**2*s35*s36 +  \
    512*s11**2*s14*s16*s22*s23*s24*s25*s34**2*s35*s36 - 512*s11*s15*s16**2*s22*s23*s24*s25*s34**2*s35*s36 - 1024*s11**2*s13*s16*s22*s24**2*s25*s34**2*s35*s36 + 1536*s11**2*s13*s14*s22**2*s25**2*s34**2*s35*s36 -  \
    512*s11*s13*s15*s16*s22**2*s25**2*s34**2*s35*s36 + 1024*s11**2*s12*s14*s22*s23*s25**2*s34**2*s35*s36 - 512*s11*s12*s15*s16*s22*s23*s25**2*s34**2*s35*s36 - 512*s11*s14*s16**2*s22*s23*s25**2*s34**2*s35*s36 -  \
    3072*s11**2*s12*s13*s22*s24*s25**2*s34**2*s35*s36 + 1024*s11*s13*s16**2*s22*s24*s25**2*s34**2*s35*s36 + 512*s11*s12*s13*s16*s22*s25**3*s34**2*s35*s36 - 512*s11**2*s13*s15*s22**2*s23*s26*s34**2*s35*s36 +  \
    512*s11*s15**3*s22**2*s23*s26*s34**2*s35*s36 + 512*s11**2*s12*s15*s22*s23**2*s26*s34**2*s35*s36 - 2560*s11**2*s14*s16*s22*s23**2*s26*s34**2*s35*s36 + 512*s11*s15*s16**2*s22*s23**2*s26*s34**2*s35*s36 +  \
    512*s11**2*s14*s15*s22*s23*s24*s26*s34**2*s35*s36 + 2048*s11**2*s13*s16*s22*s23*s24*s26*s34**2*s35*s36 - 512*s11*s15**2*s16*s22*s23*s24*s26*s34**2*s35*s36 + 1536*s11**2*s12*s16*s23**2*s24*s26*s34**2*s35*s36 -  \
    1024*s11**2*s12*s15*s23*s24**2*s26*s34**2*s35*s36 - 1024*s11**2*s13**2*s22**2*s25*s26*s34**2*s35*s36 - 512*s11*s13*s15**2*s22**2*s25*s26*s34**2*s35*s36 + 2048*s11**2*s12*s13*s22*s23*s25*s26*s34**2*s35*s36 +  \
    512*s11**2*s14**2*s22*s23*s25*s26*s34**2*s35*s36 - 1536*s11*s12*s15**2*s22*s23*s25*s26*s34**2*s35*s36 + 2048*s11*s14*s15*s16*s22*s23*s25*s26*s34**2*s35*s36 - 1536*s11*s13*s16**2*s22*s23*s25*s26*s34**2*s35*s36 -  \
    1024*s11**2*s12**2*s23**2*s25*s26*s34**2*s35*s36 - 512*s11*s12*s16**2*s23**2*s25*s26*s34**2*s35*s36 - 1024*s11**2*s13*s14*s22*s24*s25*s26*s34**2*s35*s36 - 1024*s11**2*s12*s14*s23*s24*s25*s26*s34**2*s35*s36 +  \
    2048*s11**2*s12*s13*s24**2*s25*s26*s34**2*s35*s36 + 1536*s11*s12*s13*s15*s22*s25**2*s26*s34**2*s35*s36 - 1536*s11*s13*s14*s16*s22*s25**2*s26*s34**2*s35*s36 + 1024*s11*s12**2*s15*s23*s25**2*s26*s34**2*s35*s36 +  \
    512*s11*s12*s14*s16*s23*s25**2*s26*s34**2*s35*s36 - 512*s11*s12*s13*s16*s24*s25**2*s26*s34**2*s35*s36 - 1024*s11*s12**2*s13*s25**3*s26*s34**2*s35*s36 + 1024*s11**2*s13*s14*s22*s23*s26**2*s34**2*s35*s36 -  \
    512*s11*s14*s15**2*s22*s23*s26**2*s34**2*s35*s36 - 512*s11*s13*s15*s16*s22*s23*s26**2*s34**2*s35*s36 + 1536*s11**2*s12*s14*s23**2*s26**2*s34**2*s35*s36 - 512*s11*s12*s15*s16*s23**2*s26**2*s34**2*s35*s36 -  \
    3072*s11**2*s12*s13*s23*s24*s26**2*s34**2*s35*s36 + 1024*s11*s12*s15**2*s23*s24*s26**2*s34**2*s35*s36 + 512*s11*s13*s14*s15*s22*s25*s26**2*s34**2*s35*s36 + 1024*s11*s13**2*s16*s22*s25*s26**2*s34**2*s35*s36 -  \
    1536*s11*s12*s14*s15*s23*s25*s26**2*s34**2*s35*s36 + 1536*s11*s12*s13*s16*s23*s25*s26**2*s34**2*s35*s36 - 512*s11*s12*s13*s15*s24*s25*s26**2*s34**2*s35*s36 + 1024*s11*s12*s13*s14*s25**2*s26**2*s34**2*s35*s36 +  \
    512*s11*s12*s13*s15*s23*s26**3*s34**2*s35*s36 - 1024*s11*s12*s13**2*s25*s26**3*s34**2*s35*s36 + 1024*s11**2*s14*s15*s22**3*s23*s33*s35**2*s36 - 1024*s11**2*s13*s16*s22**3*s23*s33*s35**2*s36 +  \
    1024*s11**2*s12*s16*s22**2*s23**2*s33*s35**2*s36 + 1536*s11**2*s13*s15*s22**3*s24*s33*s35**2*s36 - 2560*s11**2*s12*s15*s22**2*s23*s24*s33*s35**2*s36 + 1536*s11**2*s14*s16*s22**2*s23*s24*s33*s35**2*s36 -  \
    1024*s11**2*s14*s15*s22**2*s24**2*s33*s35**2*s36 - 512*s11**2*s13*s16*s22**2*s24**2*s33*s35**2*s36 - 1024*s11**2*s12*s16*s22*s23*s24**2*s33*s35**2*s36 + 1024*s11**2*s12*s15*s22*s24**3*s33*s35**2*s36 -  \
    512*s11**2*s13*s14*s22**3*s25*s33*s35**2*s36 - 1024*s11*s13*s15*s16*s22**3*s25*s33*s35**2*s36 - 512*s11**2*s12*s14*s22**2*s23*s25*s33*s35**2*s36 + 1024*s11*s12*s15*s16*s22**2*s23*s25*s33*s35**2*s36 -  \
    1024*s11*s14*s16**2*s22**2*s23*s25*s33*s35**2*s36 - 1024*s11**2*s12*s13*s22**2*s24*s25*s33*s35**2*s36 + 1024*s11**2*s14**2*s22**2*s24*s25*s33*s35**2*s36 + 512*s11*s14*s15*s16*s22**2*s24*s25*s33*s35**2*s36 +  \
    512*s11*s13*s16**2*s22**2*s24*s25*s33*s35**2*s36 + 2048*s11**2*s12**2*s22*s23*s24*s25*s33*s35**2*s36 + 512*s11*s12*s16**2*s22*s23*s24*s25*s33*s35**2*s36 - 1024*s11**2*s12*s14*s22*s24**2*s25*s33*s35**2*s36 -  \
    512*s11*s12*s15*s16*s22*s24**2*s25*s33*s35**2*s36 + 1024*s11*s12*s13*s16*s22**2*s25**2*s33*s35**2*s36 - 512*s11*s14**2*s16*s22**2*s25**2*s33*s35**2*s36 - 1024*s11*s12**2*s16*s22*s23*s25**2*s33*s35**2*s36 +  \
    512*s11*s12*s14*s16*s22*s24*s25**2*s33*s35**2*s36 + 1024*s11**2*s13**2*s22**3*s26*s33*s35**2*s36 - 1024*s11*s13*s15**2*s22**3*s26*s33*s35**2*s36 - 1024*s11**2*s12*s13*s22**2*s23*s26*s33*s35**2*s36 +  \
    512*s11**2*s14**2*s22**2*s23*s26*s33*s35**2*s36 + 1024*s11*s12*s15**2*s22**2*s23*s26*s33*s35**2*s36 - 2048*s11*s14*s15*s16*s22**2*s23*s26*s33*s35**2*s36 + 1024*s11*s13*s16**2*s22**2*s23*s26*s33*s35**2*s36 -  \
    1024*s11*s12*s16**2*s22*s23**2*s26*s33*s35**2*s36 - 1536*s11**2*s13*s14*s22**2*s24*s26*s33*s35**2*s36 + 512*s11*s14*s15**2*s22**2*s24*s26*s33*s35**2*s36 - 1024*s11**2*s12*s14*s22*s23*s24*s26*s33*s35**2*s36 +  \
    2048*s11*s12*s15*s16*s22*s23*s24*s26*s33*s35**2*s36 - 512*s11*s14*s16**2*s22*s23*s24*s26*s33*s35**2*s36 + 2048*s11**2*s12*s13*s22*s24**2*s26*s33*s35**2*s36 - 512*s11*s12*s15**2*s22*s24**2*s26*s33*s35**2*s36;
    v4_5= \
    512*s11*s14*s15*s16*s22*s24**2*s26*s33*s35**2*s36 + 512*s11*s12*s16**2*s23*s24**2*s26*s33*s35**2*s36 - 512*s11*s12*s15*s16*s24**3*s26*s33*s35**2*s36 + 3072*s11*s12*s13*s15*s22**2*s25*s26*s33*s35**2*s36 -  \
    512*s11*s14**2*s15*s22**2*s25*s26*s33*s35**2*s36 + 1024*s11*s13*s14*s16*s22**2*s25*s26*s33*s35**2*s36 - 3072*s11*s12**2*s15*s22*s23*s25*s26*s33*s35**2*s36 + 3072*s11*s12*s14*s16*s22*s23*s25*s26*s33*s35**2*s36 -  \
    512*s11*s12*s14*s15*s22*s24*s25*s26*s33*s35**2*s36 - 2048*s11*s12*s13*s16*s22*s24*s25*s26*s33*s35**2*s36 - 512*s11*s14**2*s16*s22*s24*s25*s26*s33*s35**2*s36 - 2048*s11*s12**2*s16*s23*s24*s25*s26*s33*s35**2*s36 +  \
    1024*s11*s12**2*s15*s24**2*s25*s26*s33*s35**2*s36 + 512*s11*s12*s14*s16*s24**2*s25*s26*s33*s35**2*s36 - 2048*s11*s12**2*s13*s22*s25**2*s26*s33*s35**2*s36 + 1024*s11*s12*s14**2*s22*s25**2*s26*s33*s35**2*s36 +  \
    2048*s11*s12**3*s23*s25**2*s26*s33*s35**2*s36 - 1024*s11*s12**2*s14*s24*s25**2*s26*s33*s35**2*s36 + 1536*s11*s13*s14*s15*s22**2*s26**2*s33*s35**2*s36 - 1024*s11*s13**2*s16*s22**2*s26**2*s33*s35**2*s36 +  \
    512*s11*s12*s14*s15*s22*s23*s26**2*s33*s35**2*s36 + 512*s11*s14**2*s16*s22*s23*s26**2*s33*s35**2*s36 + 1024*s11*s12**2*s16*s23**2*s26**2*s33*s35**2*s36 - 1536*s11*s12*s13*s15*s22*s24*s26**2*s33*s35**2*s36 -  \
    512*s11*s14**2*s15*s22*s24*s26**2*s33*s35**2*s36 + 512*s11*s13*s14*s16*s22*s24*s26**2*s33*s35**2*s36 - 512*s11*s12**2*s15*s23*s24*s26**2*s33*s35**2*s36 - 512*s11*s12*s14*s16*s23*s24*s26**2*s33*s35**2*s36 +  \
    512*s11*s12*s14*s15*s24**2*s26**2*s33*s35**2*s36 - 512*s11*s12*s13*s16*s24**2*s26**2*s33*s35**2*s36 - 2560*s11*s12*s13*s14*s22*s25*s26**2*s33*s35**2*s36 + 512*s11*s14**3*s22*s25*s26**2*s33*s35**2*s36 -  \
    512*s11*s12**2*s14*s23*s25*s26**2*s33*s35**2*s36 + 3072*s11*s12**2*s13*s24*s25*s26**2*s33*s35**2*s36 - 512*s11*s12*s14**2*s24*s25*s26**2*s33*s35**2*s36 + 1024*s11*s12*s13**2*s22*s26**3*s33*s35**2*s36 -  \
    512*s11*s13*s14**2*s22*s26**3*s33*s35**2*s36 - 1024*s11*s12**2*s13*s23*s26**3*s33*s35**2*s36 + 512*s11*s12*s13*s14*s24*s26**3*s33*s35**2*s36 - 1024*s11**2*s13*s15*s22**3*s23*s34*s35**2*s36 +  \
    1024*s11**2*s12*s15*s22**2*s23**2*s34*s35**2*s36 - 1024*s11**2*s14*s16*s22**2*s23**2*s34*s35**2*s36 + 1536*s11**2*s14*s15*s22**2*s23*s24*s34*s35**2*s36 - 512*s11**2*s13*s16*s22**2*s23*s24*s34*s35**2*s36 +  \
    1536*s11**2*s12*s16*s22*s23**2*s24*s34*s35**2*s36 - 512*s11**2*s13*s15*s22**2*s24**2*s34*s35**2*s36 - 1024*s11**2*s12*s15*s22*s23*s24**2*s34*s35**2*s36 - 512*s11**2*s14*s16*s22*s23*s24**2*s34*s35**2*s36 +  \
    512*s11**2*s13*s16*s22*s24**3*s34*s35**2*s36 + 2048*s11**2*s13**2*s22**3*s25*s34*s35**2*s36 - 3072*s11**2*s12*s13*s22**2*s23*s25*s34*s35**2*s36 + 1536*s11**2*s14**2*s22**2*s23*s25*s34*s35**2*s36 -  \
    1024*s11*s14*s15*s16*s22**2*s23*s25*s34*s35**2*s36 + 1024*s11*s13*s16**2*s22**2*s23*s25*s34*s35**2*s36 + 1024*s11**2*s12**2*s22*s23**2*s25*s34*s35**2*s36 - 1024*s11*s12*s16**2*s22*s23**2*s25*s34*s35**2*s36 -  \
    2560*s11**2*s13*s14*s22**2*s24*s25*s34*s35**2*s36 + 512*s11*s13*s15*s16*s22**2*s24*s25*s34*s35**2*s36 - 2048*s11**2*s12*s14*s22*s23*s24*s25*s34*s35**2*s36 + 512*s11*s12*s15*s16*s22*s23*s24*s25*s34*s35**2*s36 +  \
    512*s11*s14*s16**2*s22*s23*s24*s25*s34*s35**2*s36 + 3072*s11**2*s12*s13*s22*s24**2*s25*s34*s35**2*s36 - 512*s11*s13*s16**2*s22*s24**2*s25*s34*s35**2*s36 + 512*s11*s13*s14*s16*s22**2*s25**2*s34*s35**2*s36 +  \
    512*s11*s12*s14*s16*s22*s23*s25**2*s34*s35**2*s36 - 1024*s11*s12*s13*s16*s22*s24*s25**2*s34*s35**2*s36 + 512*s11**2*s13*s14*s22**2*s23*s26*s34*s35**2*s36 - 1024*s11*s14*s15**2*s22**2*s23*s26*s34*s35**2*s36 +  \
    2048*s11*s13*s15*s16*s22**2*s23*s26*s34*s35**2*s36 + 512*s11**2*s12*s14*s22*s23**2*s26*s34*s35**2*s36 - 2048*s11*s12*s15*s16*s22*s23**2*s26*s34*s35**2*s36 + 1024*s11*s14*s16**2*s22*s23**2*s26*s34*s35**2*s36 +  \
    512*s11*s13*s15**2*s22**2*s24*s26*s34*s35**2*s36 - 512*s11**2*s14**2*s22*s23*s24*s26*s34*s35**2*s36 + 512*s11*s12*s15**2*s22*s23*s24*s26*s34*s35**2*s36 - 512*s11*s13*s16**2*s22*s23*s24*s26*s34*s35**2*s36 -  \
    1024*s11**2*s12**2*s23**2*s24*s26*s34*s35**2*s36 - 512*s11*s12*s16**2*s23**2*s24*s26*s34*s35**2*s36 + 512*s11**2*s13*s14*s22*s24**2*s26*s34*s35**2*s36 - 512*s11*s13*s15*s16*s22*s24**2*s26*s34*s35**2*s36 +  \
    1024*s11**2*s12*s14*s23*s24**2*s26*s34*s35**2*s36 + 512*s11*s12*s15*s16*s23*s24**2*s26*s34*s35**2*s36 - 1024*s11**2*s12*s13*s24**3*s26*s34*s35**2*s36 + 512*s11*s13*s14*s15*s22**2*s25*s26*s34*s35**2*s36 -  \
    2560*s11*s13**2*s16*s22**2*s25*s26*s34*s35**2*s36 + 2560*s11*s12*s14*s15*s22*s23*s25*s26*s34*s35**2*s36 + 1024*s11*s12*s13*s16*s22*s23*s25*s26*s34*s35**2*s36 - 1536*s11*s14**2*s16*s22*s23*s25*s26*s34*s35**2*s36 +  \
    1536*s11*s12**2*s16*s23**2*s25*s26*s34*s35**2*s36 - 2048*s11*s12*s13*s15*s22*s24*s25*s26*s34*s35**2*s36 + 2048*s11*s13*s14*s16*s22*s24*s25*s26*s34*s35**2*s36 - 1024*s11*s12**2*s15*s23*s24*s25*s26*s34*s35**2*s36 -  \
    512*s11*s12*s13*s16*s24**2*s25*s26*s34*s35**2*s36 - 1024*s11*s12*s13*s14*s22*s25**2*s26*s34*s35**2*s36 - 1024*s11*s12**2*s14*s23*s25**2*s26*s34*s35**2*s36 + 2048*s11*s12**2*s13*s24*s25**2*s26*s34*s35**2*s36 -  \
    512*s11*s13**2*s15*s22**2*s26**2*s34*s35**2*s36 - 1024*s11*s12*s13*s15*s22*s23*s26**2*s34*s35**2*s36 + 1024*s11*s14**2*s15*s22*s23*s26**2*s34*s35**2*s36 - 1536*s11*s13*s14*s16*s22*s23*s26**2*s34*s35**2*s36 +  \
    1536*s11*s12**2*s15*s23**2*s26**2*s34*s35**2*s36 - 512*s11*s12*s14*s16*s23**2*s26**2*s34*s35**2*s36 - 512*s11*s13*s14*s15*s22*s24*s26**2*s34*s35**2*s36 + 512*s11*s13**2*s16*s22*s24*s26**2*s34*s35**2*s36 -  \
    1536*s11*s12*s14*s15*s23*s24*s26**2*s34*s35**2*s36 + 1536*s11*s12*s13*s16*s23*s24*s26**2*s34*s35**2*s36 + 1024*s11*s12*s13*s15*s24**2*s26**2*s34*s35**2*s36 + 3072*s11*s12*s13**2*s22*s25*s26**2*s34*s35**2*s36 -  \
    512*s11*s13*s14**2*s22*s25*s26**2*s34*s35**2*s36 - 3072*s11*s12**2*s13*s23*s25*s26**2*s34*s35**2*s36 + 1024*s11*s12*s14**2*s23*s25*s26**2*s34*s35**2*s36 - 512*s11*s12*s13*s14*s24*s25*s26**2*s34*s35**2*s36 +  \
    512*s11*s13**2*s14*s22*s26**3*s34*s35**2*s36 + 512*s11*s12*s13*s14*s23*s26**3*s34*s35**2*s36 - 1024*s11*s12*s13**2*s24*s26**3*s34*s35**2*s36 - 1024*s11**2*s13**2*s22**3*s24*s35**3*s36 +  \
    2048*s11**2*s12*s13*s22**2*s23*s24*s35**3*s36 - 1024*s11**2*s14**2*s22**2*s23*s24*s35**3*s36 - 1024*s11**2*s12**2*s22*s23**2*s24*s35**3*s36 + 1024*s11**2*s13*s14*s22**2*s24**2*s35**3*s36 +  \
    1024*s11**2*s12*s14*s22*s23*s24**2*s35**3*s36 - 1024*s11**2*s12*s13*s22*s24**3*s35**3*s36 + 512*s11*s13**2*s16*s22**3*s25*s35**3*s36 - 1024*s11*s12*s13*s16*s22**2*s23*s25*s35**3*s36 + 512*s11*s14**2*s16*s22**2*s23*s25*s35**3*s36 +  \
    512*s11*s12**2*s16*s22*s23**2*s25*s35**3*s36 - 512*s11*s13*s14*s16*s22**2*s24*s25*s35**3*s36 - 512*s11*s12*s14*s16*s22*s23*s24*s25*s35**3*s36 + 512*s11*s12*s13*s16*s22*s24**2*s25*s35**3*s36 +  \
    512*s11*s13**2*s15*s22**3*s26*s35**3*s36 - 1024*s11*s12*s13*s15*s22**2*s23*s26*s35**3*s36 + 512*s11*s14**2*s15*s22**2*s23*s26*s35**3*s36 + 512*s11*s12**2*s15*s22*s23**2*s26*s35**3*s36 -  \
    512*s11*s13*s14*s15*s22**2*s24*s26*s35**3*s36 + 512*s11*s13**2*s16*s22**2*s24*s26*s35**3*s36 - 512*s11*s12*s14*s15*s22*s23*s24*s26*s35**3*s36 - 1024*s11*s12*s13*s16*s22*s23*s24*s26*s35**3*s36 +  \
    512*s11*s14**2*s16*s22*s23*s24*s26*s35**3*s36 + 512*s11*s12**2*s16*s23**2*s24*s26*s35**3*s36 + 512*s11*s12*s13*s15*s22*s24**2*s26*s35**3*s36 - 512*s11*s13*s14*s16*s22*s24**2*s26*s35**3*s36 -  \
    512*s11*s12*s14*s16*s23*s24**2*s26*s35**3*s36 + 512*s11*s12*s13*s16*s24**3*s26*s35**3*s36 - 1024*s11*s12*s13**2*s22**2*s25*s26*s35**3*s36 + 2048*s11*s12**2*s13*s22*s23*s25*s26*s35**3*s36 -  \
    1024*s11*s12*s14**2*s22*s23*s25*s26*s35**3*s36 - 1024*s11*s12**3*s23**2*s25*s26*s35**3*s36 + 1024*s11*s12*s13*s14*s22*s24*s25*s26*s35**3*s36 + 1024*s11*s12**2*s14*s23*s24*s25*s26*s35**3*s36 -  \
    1024*s11*s12**2*s13*s24**2*s25*s26*s35**3*s36 - 512*s11*s13**2*s14*s22**2*s26**2*s35**3*s36 + 1024*s11*s12*s13*s14*s22*s23*s26**2*s35**3*s36 - 512*s11*s14**3*s22*s23*s26**2*s35**3*s36 - 512*s11*s12**2*s14*s23**2*s26**2*s35**3*s36 +  \
    512*s11*s13*s14**2*s22*s24*s26**2*s35**3*s36 + 512*s11*s12*s14**2*s23*s24*s26**2*s35**3*s36 - 512*s11*s12*s13*s14*s24**2*s26**2*s35**3*s36 - 1024*s11**3*s13*s22**3*s23*s33**2*s36**2 + 512*s11**2*s15**2*s22**3*s23*s33**2*s36**2 +  \
    1024*s11**3*s12*s22**2*s23**2*s33**2*s36**2 + 512*s11**2*s16**2*s22**2*s23**2*s33**2*s36**2 + 1536*s11**3*s14*s22**2*s23*s24*s33**2*s36**2 - 512*s11**2*s15*s16*s22**2*s23*s24*s33**2*s36**2 + 512*s11**3*s13*s22**2*s24**2*s33**2*s36**2 -  \
    2048*s11**3*s12*s22*s23*s24**2*s33**2*s36**2 - 512*s11**3*s14*s22*s24**3*s33**2*s36**2 + 512*s11**3*s12*s24**4*s33**2*s36**2 + 1536*s11**2*s13*s15*s22**3*s25*s33**2*s36**2 - 512*s11*s15**3*s22**3*s25*s33**2*s36**2 -  \
    2560*s11**2*s12*s15*s22**2*s23*s25*s33**2*s36**2 + 512*s11**2*s14*s16*s22**2*s23*s25*s33**2*s36**2 - 512*s11*s15*s16**2*s22**2*s23*s25*s33**2*s36**2 - 512*s11**2*s14*s15*s22**2*s24*s25*s33**2*s36**2 -  \
    2048*s11**2*s13*s16*s22**2*s24*s25*s33**2*s36**2 + 512*s11*s15**2*s16*s22**2*s24*s25*s33**2*s36**2 + 2048*s11**2*s12*s16*s22*s23*s24*s25*s33**2*s36**2 + 512*s11**2*s12*s15*s22*s24**2*s25*s33**2*s36**2 +  \
    1024*s11**2*s14*s16*s22*s24**2*s25*s33**2*s36**2 - 1024*s11**2*s12*s16*s24**3*s25*s33**2*s36**2 - 1536*s11**2*s12*s13*s22**2*s25**2*s33**2*s36**2 + 1536*s11**2*s14**2*s22**2*s25**2*s33**2*s36**2 +  \
    1536*s11*s12*s15**2*s22**2*s25**2*s33**2*s36**2 - 1536*s11*s14*s15*s16*s22**2*s25**2*s33**2*s36**2 + 1536*s11*s13*s16**2*s22**2*s25**2*s33**2*s36**2 + 2048*s11**2*s12**2*s22*s23*s25**2*s33**2*s36**2 -  \
    1024*s11*s12*s16**2*s22*s23*s25**2*s33**2*s36**2 - 2560*s11**2*s12*s14*s22*s24*s25**2*s33**2*s36**2 + 512*s11*s12*s15*s16*s22*s24*s25**2*s33**2*s36**2 - 512*s11*s14*s16**2*s22*s24*s25**2*s33**2*s36**2 +  \
    1024*s11**2*s12**2*s24**2*s25**2*s33**2*s36**2 + 512*s11*s12*s16**2*s24**2*s25**2*s33**2*s36**2 - 1536*s11*s12**2*s15*s22*s25**3*s33**2*s36**2 + 1536*s11*s12*s14*s16*s22*s25**3*s33**2*s36**2 -  \
    1024*s11*s12**2*s16*s24*s25**3*s33**2*s36**2 + 512*s11*s12**3*s25**4*s33**2*s36**2 - 2048*s11**2*s14*s15*s22**2*s23*s26*s33**2*s36**2 + 2048*s11**2*s13*s16*s22**2*s23*s26*s33**2*s36**2 -  \
    3072*s11**2*s12*s16*s22*s23**2*s26*s33**2*s36**2 - 512*s11**2*s13*s15*s22**2*s24*s26*s33**2*s36**2 + 3072*s11**2*s12*s15*s22*s23*s24*s26*s33**2*s36**2 - 1024*s11**2*s14*s16*s22*s23*s24*s26*s33**2*s36**2 +  \
    512*s11**2*s14*s15*s22*s24**2*s26*s33**2*s36**2 + 1024*s11**2*s12*s16*s23*s24**2*s26*s33**2*s36**2 - 512*s11**2*s12*s15*s24**3*s26*s33**2*s36**2 + 512*s11**2*s13*s14*s22**2*s25*s26*s33**2*s36**2 +  \
    1024*s11*s14*s15**2*s22**2*s25*s26*s33**2*s36**2 - 1024*s11*s13*s15*s16*s22**2*s25*s26*s33**2*s36**2 + 1024*s11**2*s12*s14*s22*s23*s25*s26*s33**2*s36**2 + 2048*s11*s12*s15*s16*s22*s23*s25*s26*s33**2*s36**2 +  \
    2048*s11**2*s12*s13*s22*s24*s25*s26*s33**2*s36**2 - 1536*s11**2*s14**2*s22*s24*s25*s26*s33**2*s36**2 - 1536*s11*s12*s15**2*s22*s24*s25*s26*s33**2*s36**2 + 512*s11*s14*s15*s16*s22*s24*s25*s26*s33**2*s36**2 -  \
    4096*s11**2*s12**2*s23*s24*s25*s26*s33**2*s36**2 + 1536*s11**2*s12*s14*s24**2*s25*s26*s33**2*s36**2 - 512*s11*s12*s15*s16*s24**2*s25*s26*s33**2*s36**2 - 512*s11*s12*s14*s15*s22*s25**2*s26*s33**2*s36**2 -  \
    2048*s11*s12*s13*s16*s22*s25**2*s26*s33**2*s36**2 + 512*s11*s14**2*s16*s22*s25**2*s26*s33**2*s36**2 + 1024*s11*s12**2*s16*s23*s25**2*s26*s33**2*s36**2 + 1536*s11*s12**2*s15*s24*s25**2*s26*s33**2*s36**2 -  \
    512*s11*s12*s14*s16*s24*s25**2*s26*s33**2*s36**2 - 512*s11*s12**2*s14*s25**3*s26*s33**2*s36**2 + 512*s11**2*s13**2*s22**2*s26**2*s33**2*s36**2 - 3072*s11**2*s12*s13*s22*s23*s26**2*s33**2*s36**2 +  \
    1536*s11**2*s14**2*s22*s23*s26**2*s33**2*s36**2 + 3072*s11**2*s12**2*s23**2*s26**2*s33**2*s36**2 - 512*s11**2*s13*s14*s22*s24*s26**2*s33**2*s36**2 - 1536*s11**2*s12*s14*s23*s24*s26**2*s33**2*s36**2 +  \
    512*s11**2*s12*s13*s24**2*s26**2*s33**2*s36**2 + 1024*s11*s12*s13*s15*s22*s25*s26**2*s33**2*s36**2 - 512*s11*s14**2*s15*s22*s25*s26**2*s33**2*s36**2 - 1536*s11*s12**2*s15*s23*s25*s26**2*s33**2*s36**2 +  \
    512*s11*s12*s14*s15*s24*s25*s26**2*s33**2*s36**2 + 512*s11*s12**2*s13*s25**2*s26**2*s33**2*s36**2 - 2048*s11**3*s14*s22**2*s23**2*s33*s34*s36**2 - 1024*s11**3*s13*s22**2*s23*s24*s33*s34*s36**2 -  \
    512*s11**2*s15**2*s22**2*s23*s24*s33*s34*s36**2 + 3072*s11**3*s12*s22*s23**2*s24*s33*s34*s36**2 - 512*s11**2*s16**2*s22*s23**2*s24*s33*s34*s36**2 + 1024*s11**3*s14*s22*s23*s24**2*s33*s34*s36**2
    v4_6= \
    512*s11**2*s15*s16*s22*s23*s24**2*s33*s34*s36**2 - 1024*s11**3*s12*s23*s24**3*s33*s34*s36**2 + 2560*s11**2*s14*s15*s22**2*s23*s25*s33*s34*s36**2 + 512*s11**2*s13*s16*s22**2*s23*s25*s33*s34*s36**2 -  \
    512*s11**2*s12*s16*s22*s23**2*s25*s33*s34*s36**2 - 512*s11**2*s13*s15*s22**2*s24*s25*s33*s34*s36**2 + 512*s11*s15**3*s22**2*s24*s25*s33*s34*s36**2 - 1024*s11**2*s12*s15*s22*s23*s24*s25*s33*s34*s36**2 -  \
    3072*s11**2*s14*s16*s22*s23*s24*s25*s33*s34*s36**2 + 512*s11*s15*s16**2*s22*s23*s24*s25*s33*s34*s36**2 - 512*s11**2*s14*s15*s22*s24**2*s25*s33*s34*s36**2 + 1024*s11**2*s13*s16*s22*s24**2*s25*s33*s34*s36**2 -  \
    512*s11*s15**2*s16*s22*s24**2*s25*s33*s34*s36**2 + 1536*s11**2*s12*s16*s23*s24**2*s25*s33*s34*s36**2 + 512*s11**2*s12*s15*s24**3*s25*s33*s34*s36**2 - 2560*s11**2*s13*s14*s22**2*s25**2*s33*s34*s36**2 -  \
    512*s11*s14*s15**2*s22**2*s25**2*s33*s34*s36**2 + 512*s11*s13*s15*s16*s22**2*s25**2*s33*s34*s36**2 - 512*s11*s12*s15*s16*s22*s23*s25**2*s33*s34*s36**2 + 1024*s11*s14*s16**2*s22*s23*s25**2*s33*s34*s36**2 +  \
    3072*s11**2*s12*s13*s22*s24*s25**2*s33*s34*s36**2 + 512*s11**2*s14**2*s22*s24*s25**2*s33*s34*s36**2 - 1024*s11*s12*s15**2*s22*s24*s25**2*s33*s34*s36**2 + 1536*s11*s14*s15*s16*s22*s24*s25**2*s33*s34*s36**2 -  \
    1024*s11*s13*s16**2*s22*s24*s25**2*s33*s34*s36**2 - 1024*s11**2*s12**2*s23*s24*s25**2*s33*s34*s36**2 - 512*s11*s12*s16**2*s23*s24*s25**2*s33*s34*s36**2 - 512*s11**2*s12*s14*s24**2*s25**2*s33*s34*s36**2 -  \
    512*s11*s12*s15*s16*s24**2*s25**2*s33*s34*s36**2 + 1024*s11*s12*s14*s15*s22*s25**3*s33*s34*s36**2 - 512*s11*s12*s13*s16*s22*s25**3*s33*s34*s36**2 - 1024*s11*s14**2*s16*s22*s25**3*s33*s34*s36**2 +  \
    512*s11*s12**2*s16*s23*s25**3*s33*s34*s36**2 + 512*s11*s12**2*s15*s24*s25**3*s33*s34*s36**2 + 512*s11*s12*s14*s16*s24*s25**3*s33*s34*s36**2 - 512*s11*s12**2*s14*s25**4*s33*s34*s36**2 +  \
    1536*s11**2*s13*s15*s22**2*s23*s26*s33*s34*s36**2 - 1536*s11**2*s12*s15*s22*s23**2*s26*s33*s34*s36**2 + 2560*s11**2*s14*s16*s22*s23**2*s26*s33*s34*s36**2 - 1024*s11**2*s13*s16*s22*s23*s24*s26*s33*s34*s36**2 -  \
    512*s11**2*s12*s16*s23**2*s24*s26*s33*s34*s36**2 - 512*s11**2*s12*s15*s23*s24**2*s26*s33*s34*s36**2 - 1024*s11**2*s13**2*s22**2*s25*s26*s33*s34*s36**2 - 512*s11*s13*s15**2*s22**2*s25*s26*s33*s34*s36**2 -  \
    512*s11**2*s14**2*s22*s23*s25*s26*s33*s34*s36**2 + 512*s11*s12*s15**2*s22*s23*s25*s26*s33*s34*s36**2 - 1536*s11*s14*s15*s16*s22*s23*s25*s26*s33*s34*s36**2 + 1024*s11**2*s12**2*s23**2*s25*s26*s33*s34*s36**2 +  \
    2048*s11**2*s13*s14*s22*s24*s25*s26*s33*s34*s36**2 - 512*s11*s14*s15**2*s22*s24*s25*s26*s33*s34*s36**2 + 512*s11*s13*s15*s16*s22*s24*s25*s26*s33*s34*s36**2 + 2048*s11**2*s12*s14*s23*s24*s25*s26*s33*s34*s36**2 -  \
    3072*s11**2*s12*s13*s24**2*s25*s26*s33*s34*s36**2 + 1024*s11*s12*s15**2*s24**2*s25*s26*s33*s34*s36**2 + 512*s11*s12*s13*s15*s22*s25**2*s26*s33*s34*s36**2 + 512*s11*s14**2*s15*s22*s25**2*s26*s33*s34*s36**2 -  \
    512*s11*s12**2*s15*s23*s25**2*s26*s33*s34*s36**2 - 512*s11*s12*s14*s16*s23*s25**2*s26*s33*s34*s36**2 - 1536*s11*s12*s14*s15*s24*s25**2*s26*s33*s34*s36**2 + 1536*s11*s12*s13*s16*s24*s25**2*s26*s33*s34*s36**2 +  \
    512*s11*s12*s14**2*s25**3*s26*s33*s34*s36**2 - 1024*s11**2*s13*s14*s22*s23*s26**2*s33*s34*s36**2 - 1536*s11**2*s12*s14*s23**2*s26**2*s33*s34*s36**2 + 2048*s11**2*s12*s13*s23*s24*s26**2*s33*s34*s36**2 +  \
    512*s11*s13*s14*s15*s22*s25*s26**2*s33*s34*s36**2 + 1024*s11*s12*s14*s15*s23*s25*s26**2*s33*s34*s36**2 - 1024*s11*s12*s13*s15*s24*s25*s26**2*s33*s34*s36**2 - 512*s11*s12*s13*s14*s25**2*s26**2*s33*s34*s36**2 +  \
    1024*s11**3*s13*s22**2*s23**2*s34**2*s36**2 + 512*s11**2*s15**2*s22**2*s23**2*s34**2*s36**2 - 1024*s11**3*s12*s22*s23**3*s34**2*s36**2 + 512*s11**2*s16**2*s22*s23**3*s34**2*s36**2 - 512*s11**3*s14*s22*s23**2*s24*s34**2*s36**2 -  \
    512*s11**2*s15*s16*s22*s23**2*s24*s34**2*s36**2 + 512*s11**3*s12*s23**2*s24**2*s34**2*s36**2 - 1024*s11**2*s13*s15*s22**2*s23*s25*s34**2*s36**2 - 512*s11*s15**3*s22**2*s23*s25*s34**2*s36**2 +  \
    1024*s11**2*s14*s16*s22*s23**2*s25*s34**2*s36**2 - 512*s11*s15*s16**2*s22*s23**2*s25*s34**2*s36**2 + 512*s11**2*s14*s15*s22*s23*s24*s25*s34**2*s36**2 + 512*s11*s15**2*s16*s22*s23*s24*s25*s34**2*s36**2 -  \
    512*s11**2*s12*s16*s23**2*s24*s25*s34**2*s36**2 - 512*s11**2*s12*s15*s23*s24**2*s25*s34**2*s36**2 + 1024*s11**2*s13**2*s22**2*s25**2*s34**2*s36**2 + 512*s11*s13*s15**2*s22**2*s25**2*s34**2*s36**2 -  \
    1024*s11**2*s12*s13*s22*s23*s25**2*s34**2*s36**2 + 1024*s11*s12*s15**2*s22*s23*s25**2*s34**2*s36**2 - 1024*s11*s14*s15*s16*s22*s23*s25**2*s34**2*s36**2 + 512*s11*s13*s16**2*s22*s23*s25**2*s34**2*s36**2 +  \
    512*s11**2*s12**2*s23**2*s25**2*s34**2*s36**2 - 512*s11**2*s13*s14*s22*s24*s25**2*s34**2*s36**2 - 512*s11*s13*s15*s16*s22*s24*s25**2*s34**2*s36**2 + 512*s11*s12*s15*s16*s23*s24*s25**2*s34**2*s36**2 +  \
    512*s11**2*s12*s13*s24**2*s25**2*s34**2*s36**2 - 1024*s11*s12*s13*s15*s22*s25**3*s34**2*s36**2 + 1024*s11*s13*s14*s16*s22*s25**3*s34**2*s36**2 - 512*s11*s12**2*s15*s23*s25**3*s34**2*s36**2 -  \
    512*s11*s12*s13*s16*s24*s25**3*s34**2*s36**2 + 512*s11*s12**2*s13*s25**4*s34**2*s36**2 - 512*s11**2*s14*s15*s22*s23**2*s26*s34**2*s36**2 - 512*s11**2*s13*s16*s22*s23**2*s26*s34**2*s36**2 - 512*s11**2*s12*s16*s23**3*s26*s34**2*s36**2 +  \
    1024*s11**2*s12*s15*s23**2*s24*s26*s34**2*s36**2 + 512*s11*s14*s15**2*s22*s23*s25*s26*s34**2*s36**2 + 512*s11*s13*s15*s16*s22*s23*s25*s26*s34**2*s36**2 - 512*s11**2*s12*s14*s23**2*s25*s26*s34**2*s36**2 +  \
    512*s11*s12*s15*s16*s23**2*s25*s26*s34**2*s36**2 - 1024*s11*s12*s15**2*s23*s24*s25*s26*s34**2*s36**2 - 512*s11*s13*s14*s15*s22*s25**2*s26*s34**2*s36**2 - 512*s11*s13**2*s16*s22*s25**2*s26*s34**2*s36**2 +  \
    512*s11*s12*s14*s15*s23*s25**2*s26*s34**2*s36**2 - 512*s11*s12*s13*s16*s23*s25**2*s26*s34**2*s36**2 + 1024*s11*s12*s13*s15*s24*s25**2*s26*s34**2*s36**2 - 512*s11*s12*s13*s14*s25**3*s26*s34**2*s36**2 +  \
    512*s11**2*s12*s13*s23**2*s26**2*s34**2*s36**2 - 512*s11*s12*s13*s15*s23*s25*s26**2*s34**2*s36**2 + 512*s11*s12*s13**2*s25**2*s26**2*s34**2*s36**2 - 1024*s11**2*s13*s15*s22**3*s23*s33*s35*s36**2 +  \
    1024*s11**2*s12*s15*s22**2*s23**2*s33*s35*s36**2 - 1024*s11**2*s14*s16*s22**2*s23**2*s33*s35*s36**2 - 512*s11**2*s14*s15*s22**2*s23*s24*s33*s35*s36**2 + 1536*s11**2*s13*s16*s22**2*s23*s24*s33*s35*s36**2 -  \
    512*s11**2*s12*s16*s22*s23**2*s24*s33*s35*s36**2 - 512*s11**2*s13*s15*s22**2*s24**2*s33*s35*s36**2 + 1024*s11**2*s12*s15*s22*s23*s24**2*s33*s35*s36**2 - 512*s11**2*s14*s16*s22*s23*s24**2*s33*s35*s36**2 +  \
    512*s11**2*s14*s15*s22*s24**3*s33*s35*s36**2 + 512*s11**2*s12*s16*s23*s24**3*s33*s35*s36**2 - 512*s11**2*s12*s15*s24**4*s33*s35*s36**2 - 2048*s11**2*s13**2*s22**3*s25*s33*s35*s36**2 + 1024*s11*s13*s15**2*s22**3*s25*s33*s35*s36**2 +  \
    5120*s11**2*s12*s13*s22**2*s23*s25*s33*s35*s36**2 - 2560*s11**2*s14**2*s22**2*s23*s25*s33*s35*s36**2 - 1024*s11*s12*s15**2*s22**2*s23*s25*s33*s35*s36**2 + 2048*s11*s14*s15*s16*s22**2*s23*s25*s33*s35*s36**2 -  \
    1024*s11*s13*s16**2*s22**2*s23*s25*s33*s35*s36**2 - 3072*s11**2*s12**2*s22*s23**2*s25*s33*s35*s36**2 + 1024*s11*s12*s16**2*s22*s23**2*s25*s33*s35*s36**2 + 2560*s11**2*s13*s14*s22**2*s24*s25*s33*s35*s36**2 -  \
    512*s11*s14*s15**2*s22**2*s24*s25*s33*s35*s36**2 + 3072*s11**2*s12*s14*s22*s23*s24*s25*s33*s35*s36**2 - 2048*s11*s12*s15*s16*s22*s23*s24*s25*s33*s35*s36**2 + 512*s11*s14*s16**2*s22*s23*s24*s25*s33*s35*s36**2 -  \
    2048*s11**2*s12*s13*s22*s24**2*s25*s33*s35*s36**2 - 512*s11**2*s14**2*s22*s24**2*s25*s33*s35*s36**2 + 512*s11*s12*s15**2*s22*s24**2*s25*s33*s35*s36**2 - 512*s11*s14*s15*s16*s22*s24**2*s25*s33*s35*s36**2 -  \
    1024*s11**2*s12**2*s23*s24**2*s25*s33*s35*s36**2 - 512*s11*s12*s16**2*s23*s24**2*s25*s33*s35*s36**2 + 512*s11**2*s12*s14*s24**3*s25*s33*s35*s36**2 + 512*s11*s12*s15*s16*s24**3*s25*s33*s35*s36**2 -  \
    2048*s11*s12*s13*s15*s22**2*s25**2*s33*s35*s36**2 + 512*s11*s14**2*s15*s22**2*s25**2*s33*s35*s36**2 - 512*s11*s13*s14*s16*s22**2*s25**2*s33*s35*s36**2 + 2048*s11*s12**2*s15*s22*s23*s25**2*s33*s35*s36**2 -  \
    1536*s11*s12*s14*s16*s22*s23*s25**2*s33*s35*s36**2 + 512*s11*s12*s13*s16*s22*s24*s25**2*s33*s35*s36**2 + 512*s11*s14**2*s16*s22*s24*s25**2*s33*s35*s36**2 + 1536*s11*s12**2*s16*s23*s24*s25**2*s33*s35*s36**2 -  \
    512*s11*s12**2*s15*s24**2*s25**2*s33*s35*s36**2 - 512*s11*s12*s14*s16*s24**2*s25**2*s33*s35*s36**2 + 1024*s11*s12**2*s13*s22*s25**3*s33*s35*s36**2 - 512*s11*s12*s14**2*s22*s25**3*s33*s35*s36**2 -  \
    1024*s11*s12**3*s23*s25**3*s33*s35*s36**2 + 512*s11*s12**2*s14*s24*s25**3*s33*s35*s36**2 + 512*s11**2*s13*s14*s22**2*s23*s26*s33*s35*s36**2 + 1024*s11*s14*s15**2*s22**2*s23*s26*s33*s35*s36**2 -  \
    1024*s11*s13*s15*s16*s22**2*s23*s26*s33*s35*s36**2 + 512*s11**2*s12*s14*s22*s23**2*s26*s33*s35*s36**2 + 1024*s11*s12*s15*s16*s22*s23**2*s26*s33*s35*s36**2 + 512*s11*s13*s15**2*s22**2*s24*s26*s33*s35*s36**2 -  \
    2048*s11**2*s12*s13*s22*s23*s24*s26*s33*s35*s36**2 + 512*s11**2*s14**2*s22*s23*s24*s26*s33*s35*s36**2 - 1536*s11*s12*s15**2*s22*s23*s24*s26*s33*s35*s36**2 + 512*s11*s14*s15*s16*s22*s23*s24*s26*s33*s35*s36**2 +  \
    1024*s11**2*s12**2*s23**2*s24*s26*s33*s35*s36**2 - 512*s11*s14*s15**2*s22*s24**2*s26*s33*s35*s36**2 - 512*s11**2*s12*s14*s23*s24**2*s26*s33*s35*s36**2 - 512*s11*s12*s15*s16*s23*s24**2*s26*s33*s35*s36**2 +  \
    512*s11*s12*s15**2*s24**3*s26*s33*s35*s36**2 - 3072*s11*s13*s14*s15*s22**2*s25*s26*s33*s35*s36**2 + 2560*s11*s13**2*s16*s22**2*s25*s26*s33*s35*s36**2 - 1024*s11*s12*s14*s15*s22*s23*s25*s26*s33*s35*s36**2 -  \
    2048*s11*s12*s13*s16*s22*s23*s25*s26*s33*s35*s36**2 - 512*s11*s12**2*s16*s23**2*s25*s26*s33*s35*s36**2 + 2048*s11*s12*s13*s15*s22*s24*s25*s26*s33*s35*s36**2 + 1536*s11*s14**2*s15*s22*s24*s25*s26*s33*s35*s36**2 -  \
    1536*s11*s13*s14*s16*s22*s24*s25*s26*s33*s35*s36**2 + 2048*s11*s12**2*s15*s23*s24*s25*s26*s33*s35*s36**2 - 1536*s11*s12*s14*s15*s24**2*s25*s26*s33*s35*s36**2 + 1536*s11*s12*s13*s16*s24**2*s25*s26*s33*s35*s36**2 +  \
    3584*s11*s12*s13*s14*s22*s25**2*s26*s33*s35*s36**2 - 1024*s11*s14**3*s22*s25**2*s26*s33*s35*s36**2 - 512*s11*s12**2*s14*s23*s25**2*s26*s33*s35*s36**2 - 3072*s11*s12**2*s13*s24*s25**2*s26*s33*s35*s36**2 +  \
    1024*s11*s12*s14**2*s24*s25**2*s26*s33*s35*s36**2 - 512*s11*s13**2*s15*s22**2*s26**2*s33*s35*s36**2 + 2048*s11*s12*s13*s15*s22*s23*s26**2*s33*s35*s36**2 - 1024*s11*s14**2*s15*s22*s23*s26**2*s33*s35*s36**2 -  \
    1536*s11*s12**2*s15*s23**2*s26**2*s33*s35*s36**2 + 512*s11*s13*s14*s15*s22*s24*s26**2*s33*s35*s36**2 + 1024*s11*s12*s14*s15*s23*s24*s26**2*s33*s35*s36**2 - 512*s11*s12*s13*s15*s24**2*s26**2*s33*s35*s36**2 -  \
    2048*s11*s12*s13**2*s22*s25*s26**2*s33*s35*s36**2 + 1024*s11*s13*s14**2*s22*s25*s26**2*s33*s35*s36**2 + 2048*s11*s12**2*s13*s23*s25*s26**2*s33*s35*s36**2 - 1024*s11*s12*s13*s14*s24*s25*s26**2*s33*s35*s36**2 -  \
    1024*s11**2*s14*s15*s22**2*s23**2*s34*s35*s36**2 + 1024*s11**2*s13*s16*s22**2*s23**2*s34*s35*s36**2 - 1024*s11**2*s12*s16*s22*s23**3*s34*s35*s36**2 + 1536*s11**2*s13*s15*s22**2*s23*s24*s34*s35*s36**2 -  \
    512*s11**2*s12*s15*s22*s23**2*s24*s34*s35*s36**2 + 1536*s11**2*s14*s16*s22*s23**2*s24*s34*s35*s36**2 - 512*s11**2*s14*s15*s22*s23*s24**2*s34*s35*s36**2 - 1024*s11**2*s13*s16*s22*s23*s24**2*s34*s35*s36**2 -  \
    512*s11**2*s12*s16*s23**2*s24**2*s34*s35*s36**2 + 512*s11**2*s12*s15*s23*s24**3*s34*s35*s36**2 + 512*s11**2*s13*s14*s22**2*s23*s25*s34*s35*s36**2 + 1024*s11*s14*s15**2*s22**2*s23*s25*s34*s35*s36**2 -  \
    2048*s11*s13*s15*s16*s22**2*s23*s25*s34*s35*s36**2 + 512*s11**2*s12*s14*s22*s23**2*s25*s34*s35*s36**2 + 2048*s11*s12*s15*s16*s22*s23**2*s25*s34*s35*s36**2 - 1024*s11*s14*s16**2*s22*s23**2*s25*s34*s35*s36**2 -  \
    1024*s11**2*s13**2*s22**2*s24*s25*s34*s35*s36**2 - 512*s11*s13*s15**2*s22**2*s24*s25*s34*s35*s36**2 - 512*s11**2*s14**2*s22*s23*s24*s25*s34*s35*s36**2 - 512*s11*s12*s15**2*s22*s23*s24*s25*s34*s35*s36**2 +  \
    512*s11*s13*s16**2*s22*s23*s24*s25*s34*s35*s36**2 + 512*s11*s12*s16**2*s23**2*s24*s25*s34*s35*s36**2 + 1024*s11**2*s13*s14*s22*s24**2*s25*s34*s35*s36**2 + 512*s11*s13*s15*s16*s22*s24**2*s25*s34*s35*s36**2 +  \
    512*s11**2*s12*s14*s23*s24**2*s25*s34*s35*s36**2 - 512*s11*s12*s15*s16*s23*s24**2*s25*s34*s35*s36**2 - 1024*s11**2*s12*s13*s24**3*s25*s34*s35*s36**2 - 512*s11*s13*s14*s15*s22**2*s25**2*s34*s35*s36**2
    v4_7= \
    1536*s11*s13**2*s16*s22**2*s25**2*s34*s35*s36**2 - 1536*s11*s12*s14*s15*s22*s23*s25**2*s34*s35*s36**2 - 1024*s11*s12*s13*s16*s22*s23*s25**2*s34*s35*s36**2 + 1024*s11*s14**2*s16*s22*s23*s25**2*s34*s35*s36**2 -  \
    512*s11*s12**2*s16*s23**2*s25**2*s34*s35*s36**2 + 1536*s11*s12*s13*s15*s22*s24*s25**2*s34*s35*s36**2 - 1536*s11*s13*s14*s16*s22*s24*s25**2*s34*s35*s36**2 + 512*s11*s12**2*s15*s23*s24*s25**2*s34*s35*s36**2 -  \
    512*s11*s12*s14*s16*s23*s24*s25**2*s34*s35*s36**2 + 1024*s11*s12*s13*s16*s24**2*s25**2*s34*s35*s36**2 + 512*s11*s12*s13*s14*s22*s25**3*s34*s35*s36**2 + 512*s11*s12**2*s14*s23*s25**3*s34*s35*s36**2 -  \
    1024*s11*s12**2*s13*s24*s25**3*s34*s35*s36**2 + 1024*s11**2*s13**2*s22**2*s23*s26*s34*s35*s36**2 - 1024*s11*s13*s15**2*s22**2*s23*s26*s34*s35*s36**2 - 3072*s11**2*s12*s13*s22*s23**2*s26*s34*s35*s36**2 +  \
    1536*s11**2*s14**2*s22*s23**2*s26*s34*s35*s36**2 + 1024*s11*s12*s15**2*s22*s23**2*s26*s34*s35*s36**2 - 1024*s11*s14*s15*s16*s22*s23**2*s26*s34*s35*s36**2 + 2048*s11**2*s12**2*s23**3*s26*s34*s35*s36**2 -  \
    2048*s11**2*s13*s14*s22*s23*s24*s26*s34*s35*s36**2 + 512*s11*s14*s15**2*s22*s23*s24*s26*s34*s35*s36**2 + 512*s11*s13*s15*s16*s22*s23*s24*s26*s34*s35*s36**2 - 2560*s11**2*s12*s14*s23**2*s24*s26*s34*s35*s36**2 +  \
    512*s11*s12*s15*s16*s23**2*s24*s26*s34*s35*s36**2 + 3072*s11**2*s12*s13*s23*s24**2*s26*s34*s35*s36**2 - 512*s11*s12*s15**2*s23*s24**2*s26*s34*s35*s36**2 + 1536*s11*s13**2*s15*s22**2*s25*s26*s34*s35*s36**2 +  \
    1024*s11*s12*s13*s15*s22*s23*s25*s26*s34*s35*s36**2 - 1536*s11*s14**2*s15*s22*s23*s25*s26*s34*s35*s36**2 + 2560*s11*s13*s14*s16*s22*s23*s25*s26*s34*s35*s36**2 - 2560*s11*s12**2*s15*s23**2*s25*s26*s34*s35*s36**2 +  \
    512*s11*s12*s14*s16*s23**2*s25*s26*s34*s35*s36**2 - 1024*s11*s13**2*s16*s22*s24*s25*s26*s34*s35*s36**2 + 2048*s11*s12*s14*s15*s23*s24*s25*s26*s34*s35*s36**2 - 2048*s11*s12*s13*s16*s23*s24*s25*s26*s34*s35*s36**2 -  \
    512*s11*s12*s13*s15*s24**2*s25*s26*s34*s35*s36**2 - 3072*s11*s12*s13**2*s22*s25**2*s26*s34*s35*s36**2 + 1024*s11*s13*s14**2*s22*s25**2*s26*s34*s35*s36**2 + 3072*s11*s12**2*s13*s23*s25**2*s26*s34*s35*s36**2 -  \
    512*s11*s12*s14**2*s23*s25**2*s26*s34*s35*s36**2 - 512*s11*s12*s13*s14*s24*s25**2*s26*s34*s35*s36**2 + 512*s11*s13*s14*s15*s22*s23*s26**2*s34*s35*s36**2 + 512*s11*s12*s14*s15*s23**2*s26**2*s34*s35*s36**2 -  \
    1024*s11*s12*s13*s15*s23*s24*s26**2*s34*s35*s36**2 - 1024*s11*s13**2*s14*s22*s25*s26**2*s34*s35*s36**2 - 1024*s11*s12*s13*s14*s23*s25*s26**2*s34*s35*s36**2 + 2048*s11*s12*s13**2*s24*s25*s26**2*s34*s35*s36**2 +  \
    1024*s11**2*s13**2*s22**3*s23*s35**2*s36**2 - 2048*s11**2*s12*s13*s22**2*s23**2*s35**2*s36**2 + 1024*s11**2*s14**2*s22**2*s23**2*s35**2*s36**2 + 1024*s11**2*s12**2*s22*s23**3*s35**2*s36**2 - 1024*s11**2*s13*s14*s22**2*s23*s24*s35**2*s36**2 -  \
    1024*s11**2*s12*s14*s22*s23**2*s24*s35**2*s36**2 + 512*s11**2*s13**2*s22**2*s24**2*s35**2*s36**2 + 512*s11**2*s14**2*s22*s23*s24**2*s35**2*s36**2 + 512*s11**2*s12**2*s23**2*s24**2*s35**2*s36**2 - 512*s11**2*s13*s14*s22*s24**3*s35**2*s36**2 -  \
    512*s11**2*s12*s14*s23*s24**3*s35**2*s36**2 + 512*s11**2*s12*s13*s24**4*s35**2*s36**2 - 512*s11*s13**2*s15*s22**3*s25*s35**2*s36**2 + 1024*s11*s12*s13*s15*s22**2*s23*s25*s35**2*s36**2 - 512*s11*s14**2*s15*s22**2*s23*s25*s35**2*s36**2 -  \
    512*s11*s12**2*s15*s22*s23**2*s25*s35**2*s36**2 + 512*s11*s13*s14*s15*s22**2*s24*s25*s35**2*s36**2 - 512*s11*s13**2*s16*s22**2*s24*s25*s35**2*s36**2 + 512*s11*s12*s14*s15*s22*s23*s24*s25*s35**2*s36**2 +  \
    1024*s11*s12*s13*s16*s22*s23*s24*s25*s35**2*s36**2 - 512*s11*s14**2*s16*s22*s23*s24*s25*s35**2*s36**2 - 512*s11*s12**2*s16*s23**2*s24*s25*s35**2*s36**2 - 512*s11*s12*s13*s15*s22*s24**2*s25*s35**2*s36**2 +  \
    512*s11*s13*s14*s16*s22*s24**2*s25*s35**2*s36**2 + 512*s11*s12*s14*s16*s23*s24**2*s25*s35**2*s36**2 - 512*s11*s12*s13*s16*s24**3*s25*s35**2*s36**2 + 512*s11*s12*s13**2*s22**2*s25**2*s35**2*s36**2 -  \
    1024*s11*s12**2*s13*s22*s23*s25**2*s35**2*s36**2 + 512*s11*s12*s14**2*s22*s23*s25**2*s35**2*s36**2 + 512*s11*s12**3*s23**2*s25**2*s35**2*s36**2 - 512*s11*s12*s13*s14*s22*s24*s25**2*s35**2*s36**2 -  \
    512*s11*s12**2*s14*s23*s24*s25**2*s35**2*s36**2 + 512*s11*s12**2*s13*s24**2*s25**2*s35**2*s36**2 - 512*s11*s13**2*s16*s22**2*s23*s26*s35**2*s36**2 + 1024*s11*s12*s13*s16*s22*s23**2*s26*s35**2*s36**2 -  \
    512*s11*s14**2*s16*s22*s23**2*s26*s35**2*s36**2 - 512*s11*s12**2*s16*s23**3*s26*s35**2*s36**2 - 512*s11*s13**2*s15*s22**2*s24*s26*s35**2*s36**2 + 1024*s11*s12*s13*s15*s22*s23*s24*s26*s35**2*s36**2 -  \
    512*s11*s14**2*s15*s22*s23*s24*s26*s35**2*s36**2 + 512*s11*s13*s14*s16*s22*s23*s24*s26*s35**2*s36**2 - 512*s11*s12**2*s15*s23**2*s24*s26*s35**2*s36**2 + 512*s11*s12*s14*s16*s23**2*s24*s26*s35**2*s36**2 +  \
    512*s11*s13*s14*s15*s22*s24**2*s26*s35**2*s36**2 + 512*s11*s12*s14*s15*s23*s24**2*s26*s35**2*s36**2 - 512*s11*s12*s13*s16*s23*s24**2*s26*s35**2*s36**2 - 512*s11*s12*s13*s15*s24**3*s26*s35**2*s36**2 +  \
    1024*s11*s13**2*s14*s22**2*s25*s26*s35**2*s36**2 - 2048*s11*s12*s13*s14*s22*s23*s25*s26*s35**2*s36**2 + 1024*s11*s14**3*s22*s23*s25*s26*s35**2*s36**2 + 1024*s11*s12**2*s14*s23**2*s25*s26*s35**2*s36**2 -  \
    1024*s11*s13*s14**2*s22*s24*s25*s26*s35**2*s36**2 - 1024*s11*s12*s14**2*s23*s24*s25*s26*s35**2*s36**2 + 1024*s11*s12*s13*s14*s24**2*s25*s26*s35**2*s36**2 + 512*s11*s13**3*s22**2*s26**2*s35**2*s36**2 -  \
    1024*s11*s12*s13**2*s22*s23*s26**2*s35**2*s36**2 + 512*s11*s13*s14**2*s22*s23*s26**2*s35**2*s36**2 + 512*s11*s12**2*s13*s23**2*s26**2*s35**2*s36**2 - 512*s11*s13**2*s14*s22*s24*s26**2*s35**2*s36**2 -  \
    512*s11*s12*s13*s14*s23*s24*s26**2*s35**2*s36**2 + 512*s11*s12*s13**2*s24**2*s26**2*s35**2*s36**2 + 1024*s11**2*s14*s15*s22**2*s23**2*s33*s36**3 - 1024*s11**2*s13*s16*s22**2*s23**2*s33*s36**3 +  \
    1024*s11**2*s12*s16*s22*s23**3*s33*s36**3 + 512*s11**2*s13*s15*s22**2*s23*s24*s33*s36**3 - 1536*s11**2*s12*s15*s22*s23**2*s24*s33*s36**3 + 512*s11**2*s14*s16*s22*s23**2*s24*s33*s36**3 -  \
    512*s11**2*s14*s15*s22*s23*s24**2*s33*s36**3 - 512*s11**2*s12*s16*s23**2*s24**2*s33*s36**3 + 512*s11**2*s12*s15*s23*s24**3*s33*s36**3 - 512*s11**2*s13*s14*s22**2*s23*s25*s33*s36**3 - 1024*s11*s14*s15**2*s22**2*s23*s25*s33*s36**3 +  \
    1024*s11*s13*s15*s16*s22**2*s23*s25*s33*s36**3 - 512*s11**2*s12*s14*s22*s23**2*s25*s33*s36**3 - 1024*s11*s12*s15*s16*s22*s23**2*s25*s33*s36**3 + 1024*s11**2*s13**2*s22**2*s24*s25*s33*s36**3 -  \
    512*s11*s13*s15**2*s22**2*s24*s25*s33*s36**3 - 2048*s11**2*s12*s13*s22*s23*s24*s25*s33*s36**3 + 1536*s11**2*s14**2*s22*s23*s24*s25*s33*s36**3 + 1536*s11*s12*s15**2*s22*s23*s24*s25*s33*s36**3 -  \
    512*s11*s14*s15*s16*s22*s23*s24*s25*s33*s36**3 + 2048*s11**2*s12**2*s23**2*s24*s25*s33*s36**3 - 1024*s11**2*s13*s14*s22*s24**2*s25*s33*s36**3 + 512*s11*s14*s15**2*s22*s24**2*s25*s33*s36**3 -  \
    1536*s11**2*s12*s14*s23*s24**2*s25*s33*s36**3 + 512*s11*s12*s15*s16*s23*s24**2*s25*s33*s36**3 + 1024*s11**2*s12*s13*s24**3*s25*s33*s36**3 - 512*s11*s12*s15**2*s24**3*s25*s33*s36**3 + 1536*s11*s13*s14*s15*s22**2*s25**2*s33*s36**3 -  \
    1536*s11*s13**2*s16*s22**2*s25**2*s33*s36**3 + 512*s11*s12*s14*s15*s22*s23*s25**2*s33*s36**3 + 2048*s11*s12*s13*s16*s22*s23*s25**2*s33*s36**3 - 512*s11*s14**2*s16*s22*s23*s25**2*s33*s36**3 -  \
    512*s11*s12**2*s16*s23**2*s25**2*s33*s36**3 - 512*s11*s12*s13*s15*s22*s24*s25**2*s33*s36**3 - 1024*s11*s14**2*s15*s22*s24*s25**2*s33*s36**3 + 1024*s11*s13*s14*s16*s22*s24*s25**2*s33*s36**3 -  \
    1536*s11*s12**2*s15*s23*s24*s25**2*s33*s36**3 + 512*s11*s12*s14*s16*s23*s24*s25**2*s33*s36**3 + 1024*s11*s12*s14*s15*s24**2*s25**2*s33*s36**3 - 1024*s11*s12*s13*s16*s24**2*s25**2*s33*s36**3 -  \
    1536*s11*s12*s13*s14*s22*s25**3*s33*s36**3 + 512*s11*s14**3*s22*s25**3*s33*s36**3 + 512*s11*s12**2*s14*s23*s25**3*s33*s36**3 + 1024*s11*s12**2*s13*s24*s25**3*s33*s36**3 - 512*s11*s12*s14**2*s24*s25**3*s33*s36**3 -  \
    1024*s11**2*s13**2*s22**2*s23*s26*s33*s36**3 + 3072*s11**2*s12*s13*s22*s23**2*s26*s33*s36**3 - 1536*s11**2*s14**2*s22*s23**2*s26*s33*s36**3 - 2048*s11**2*s12**2*s23**3*s26*s33*s36**3 + 1024*s11**2*s13*s14*s22*s23*s24*s26*s33*s36**3 +  \
    1536*s11**2*s12*s14*s23**2*s24*s26*s33*s36**3 - 1024*s11**2*s12*s13*s23*s24**2*s26*s33*s36**3 + 512*s11*s13**2*s15*s22**2*s25*s26*s33*s36**3 - 2048*s11*s12*s13*s15*s22*s23*s25*s26*s33*s36**3 +  \
    1024*s11*s14**2*s15*s22*s23*s25*s26*s33*s36**3 + 1536*s11*s12**2*s15*s23**2*s25*s26*s33*s36**3 - 512*s11*s13*s14*s15*s22*s24*s25*s26*s33*s36**3 - 1024*s11*s12*s14*s15*s23*s24*s25*s26*s33*s36**3 +  \
    512*s11*s12*s13*s15*s24**2*s25*s26*s33*s36**3 + 1024*s11*s12*s13**2*s22*s25**2*s26*s33*s36**3 - 512*s11*s13*s14**2*s22*s25**2*s26*s33*s36**3 - 1024*s11*s12**2*s13*s23*s25**2*s26*s33*s36**3 +  \
    512*s11*s12*s13*s14*s24*s25**2*s26*s33*s36**3 - 1024*s11**2*s13*s15*s22**2*s23**2*s34*s36**3 + 1024*s11**2*s12*s15*s22*s23**3*s34*s36**3 - 1024*s11**2*s14*s16*s22*s23**3*s34*s36**3 + 512*s11**2*s14*s15*s22*s23**2*s24*s34*s36**3 +  \
    512*s11**2*s13*s16*s22*s23**2*s24*s34*s36**3 + 512*s11**2*s12*s16*s23**3*s24*s34*s36**3 - 512*s11**2*s12*s15*s23**2*s24**2*s34*s36**3 + 1024*s11*s13*s15**2*s22**2*s23*s25*s34*s36**3 + 1024*s11**2*s12*s13*s22*s23**2*s25*s34*s36**3 -  \
    512*s11**2*s14**2*s22*s23**2*s25*s34*s36**3 - 1024*s11*s12*s15**2*s22*s23**2*s25*s34*s36**3 + 1024*s11*s14*s15*s16*s22*s23**2*s25*s34*s36**3 - 1024*s11**2*s12**2*s23**3*s25*s34*s36**3 -  \
    512*s11*s14*s15**2*s22*s23*s24*s25*s34*s36**3 - 512*s11*s13*s15*s16*s22*s23*s24*s25*s34*s36**3 + 512*s11**2*s12*s14*s23**2*s24*s25*s34*s36**3 - 512*s11*s12*s15*s16*s23**2*s24*s25*s34*s36**3 +  \
    512*s11*s12*s15**2*s23*s24**2*s25*s34*s36**3 - 1024*s11*s13**2*s15*s22**2*s25**2*s34*s36**3 + 512*s11*s14**2*s15*s22*s23*s25**2*s34*s36**3 - 1024*s11*s13*s14*s16*s22*s23*s25**2*s34*s36**3 +  \
    1024*s11*s12**2*s15*s23**2*s25**2*s34*s36**3 + 512*s11*s13*s14*s15*s22*s24*s25**2*s34*s36**3 + 512*s11*s13**2*s16*s22*s24*s25**2*s34*s36**3 - 512*s11*s12*s14*s15*s23*s24*s25**2*s34*s36**3 +  \
    512*s11*s12*s13*s16*s23*s24*s25**2*s34*s36**3 - 512*s11*s12*s13*s15*s24**2*s25**2*s34*s36**3 + 1024*s11*s12*s13**2*s22*s25**3*s34*s36**3 - 512*s11*s13*s14**2*s22*s25**3*s34*s36**3 - 1024*s11*s12**2*s13*s23*s25**3*s34*s36**3 +  \
    512*s11*s12*s13*s14*s24*s25**3*s34*s36**3 + 512*s11**2*s13*s14*s22*s23**2*s26*s34*s36**3 + 512*s11**2*s12*s14*s23**3*s26*s34*s36**3 - 1024*s11**2*s12*s13*s23**2*s24*s26*s34*s36**3 -  \
    512*s11*s13*s14*s15*s22*s23*s25*s26*s34*s36**3 - 512*s11*s12*s14*s15*s23**2*s25*s26*s34*s36**3 + 1024*s11*s12*s13*s15*s23*s24*s25*s26*s34*s36**3 + 512*s11*s13**2*s14*s22*s25**2*s26*s34*s36**3 +  \
    512*s11*s12*s13*s14*s23*s25**2*s26*s34*s36**3 - 1024*s11*s12*s13**2*s24*s25**2*s26*s34*s36**3 - 1024*s11**2*s13**2*s22**2*s23*s24*s35*s36**3 + 2048*s11**2*s12*s13*s22*s23**2*s24*s35*s36**3 -  \
    1024*s11**2*s14**2*s22*s23**2*s24*s35*s36**3 - 1024*s11**2*s12**2*s23**3*s24*s35*s36**3 + 1024*s11**2*s13*s14*s22*s23*s24**2*s35*s36**3 + 1024*s11**2*s12*s14*s23**2*s24**2*s35*s36**3 - 1024*s11**2*s12*s13*s23*s24**3*s35*s36**3 +  \
    512*s11*s13**2*s16*s22**2*s23*s25*s35*s36**3 - 1024*s11*s12*s13*s16*s22*s23**2*s25*s35*s36**3 + 512*s11*s14**2*s16*s22*s23**2*s25*s35*s36**3 + 512*s11*s12**2*s16*s23**3*s25*s35*s36**3 +  \
    512*s11*s13**2*s15*s22**2*s24*s25*s35*s36**3 - 1024*s11*s12*s13*s15*s22*s23*s24*s25*s35*s36**3 + 512*s11*s14**2*s15*s22*s23*s24*s25*s35*s36**3 - 512*s11*s13*s14*s16*s22*s23*s24*s25*s35*s36**3 +  \
    512*s11*s12**2*s15*s23**2*s24*s25*s35*s36**3 - 512*s11*s12*s14*s16*s23**2*s24*s25*s35*s36**3 - 512*s11*s13*s14*s15*s22*s24**2*s25*s35*s36**3 - 512*s11*s12*s14*s15*s23*s24**2*s25*s35*s36**3 +  \
    512*s11*s12*s13*s16*s23*s24**2*s25*s35*s36**3 + 512*s11*s12*s13*s15*s24**3*s25*s35*s36**3 - 512*s11*s13**2*s14*s22**2*s25**2*s35*s36**3 + 1024*s11*s12*s13*s14*s22*s23*s25**2*s35*s36**3 - 512*s11*s14**3*s22*s23*s25**2*s35*s36**3 -  \
    512*s11*s12**2*s14*s23**2*s25**2*s35*s36**3 + 512*s11*s13*s14**2*s22*s24*s25**2*s35*s36**3 + 512*s11*s12*s14**2*s23*s24*s25**2*s35*s36**3 - 512*s11*s12*s13*s14*s24**2*s25**2*s35*s36**3
    v4_8= \
    512*s11*s13**2*s15*s22**2*s23*s26*s35*s36**3 - 1024*s11*s12*s13*s15*s22*s23**2*s26*s35*s36**3 + 512*s11*s14**2*s15*s22*s23**2*s26*s35*s36**3 + 512*s11*s12**2*s15*s23**3*s26*s35*s36**3 -  \
    512*s11*s13*s14*s15*s22*s23*s24*s26*s35*s36**3 - 512*s11*s12*s14*s15*s23**2*s24*s26*s35*s36**3 + 512*s11*s12*s13*s15*s23*s24**2*s26*s35*s36**3 - 1024*s11*s13**3*s22**2*s25*s26*s35*s36**3 +  \
    2048*s11*s12*s13**2*s22*s23*s25*s26*s35*s36**3 - 1024*s11*s13*s14**2*s22*s23*s25*s26*s35*s36**3 - 1024*s11*s12**2*s13*s23**2*s25*s26*s35*s36**3 + 1024*s11*s13**2*s14*s22*s24*s25*s26*s35*s36**3 +  \
    1024*s11*s12*s13*s14*s23*s24*s25*s26*s35*s36**3 - 1024*s11*s12*s13**2*s24**2*s25*s26*s35*s36**3 + 512*s11**2*s13**2*s22**2*s23**2*s36**4 - 1024*s11**2*s12*s13*s22*s23**3*s36**4 + 512*s11**2*s14**2*s22*s23**3*s36**4 +  \
    512*s11**2*s12**2*s23**4*s36**4 - 512*s11**2*s13*s14*s22*s23**2*s24*s36**4 - 512*s11**2*s12*s14*s23**3*s24*s36**4 + 512*s11**2*s12*s13*s23**2*s24**2*s36**4 - 512*s11*s13**2*s15*s22**2*s23*s25*s36**4 +  \
    1024*s11*s12*s13*s15*s22*s23**2*s25*s36**4 - 512*s11*s14**2*s15*s22*s23**2*s25*s36**4 - 512*s11*s12**2*s15*s23**3*s25*s36**4 + 512*s11*s13*s14*s15*s22*s23*s24*s25*s36**4 + 512*s11*s12*s14*s15*s23**2*s24*s25*s36**4 -  \
    512*s11*s12*s13*s15*s23*s24**2*s25*s36**4 + 512*s11*s13**3*s22**2*s25**2*s36**4 - 1024*s11*s12*s13**2*s22*s23*s25**2*s36**4 + 512*s11*s13*s14**2*s22*s23*s25**2*s36**4 + 512*s11*s12**2*s13*s23**2*s25**2*s36**4 -  \
    512*s11*s13**2*s14*s22*s24*s25**2*s36**4 - 512*s11*s12*s13*s14*s23*s24*s25**2*s36**4 + 512*s11*s12*s13**2*s24**2*s25**2*s36**4
    v4=v4_0+v4_1+v4_2+v4_3+v4_4+v4_5+v4_6+v4_7+v4_8

    return v0, v1, v2, v3, v4

def R_bn_spherical(phi,theta, Rb_om=None, Rn_om=None):
    # R_BN_SPHERICAL from navigation frame to body frame

    assert(len(phi) == len(theta), "error")
    N = len(phi)
    R = np.zeros((3,3, N))

    for n in range(N):
        t = theta[n]
        p = phi[n]
        if Rb_om is None:
            R[:,:,n] = np.array([[np.cos(t)*np.sin(p), np.sin(t)*np.sin(p), np.cos(p)],
                                 [np.cos(t)*np.cos(p), np.sin(t)*np.cos(p), -np.sin(p)],
                                 [-np.sin(t), np.cos(t), 0]])
        else:
            Rom_n=Rn_om(t).T
            R[:,:,n]=Rb_om(p)*Rom_n

    return R


def getVcoeffs(s11,s12,s13,s14,s15,s16,s22,s23,s24,s25,s26,s33,s34,s35,s36,q7a,q8a,q9a):
    v0=512*s11**2*(q9a**2*(s11*s22**2 + s26*((-s16)*s22 + s12*s26)) +  \
        q9a*(q8a*s16*s22 - 2*q8a*s12*s26 + q7a*s22*s26)*s36 +  \
        q8a*(q8a*s12 - q7a*s22)*s36**2)**2

    v1_0=-2048*q9a**3*s11**4*s22**4*s33 + 4096*q9a**3*s11**3*s16*s22**3*s26*s33 - 4096*q9a**3*s11**3*s12*s22**2*s26**2*s33 - 2048*q9a**3*s11**2*s16**2*s22**2*s26**2*s33 + 4096*q9a**3*s11**2*s12*s16*s22*s26**3*s33 -  \
    2048*q9a**3*s11**2*s12**2*s26**4*s33 + 1024*q9a**3*s11**4*s22**3*s24*s34 - 512*q9a**3*s11**3*s16*s22**3*s25*s34 - 512*q9a**3*s11**3*s15*s22**3*s26*s34 - 1536*q9a**3*s11**3*s16*s22**2*s24*s26*s34 +  \
    1024*q9a**3*s11**3*s12*s22**2*s25*s26*s34 + 512*q9a**3*s11**2*s16**2*s22**2*s25*s26*s34 + 512*q9a**3*s11**3*s14*s22**2*s26**2*s34 + 512*q9a**3*s11**2*s15*s16*s22**2*s26**2*s34 + 1024*q9a**3*s11**3*s12*s22*s24*s26**2*s34 +  \
    512*q9a**3*s11**2*s16**2*s22*s24*s26**2*s34 - 1536*q9a**3*s11**2*s12*s16*s22*s25*s26**2*s34 - 512*q9a**3*s11**2*s12*s15*s22*s26**3*s34 - 512*q9a**3*s11**2*s14*s16*s22*s26**3*s34 - 512*q9a**3*s11**2*s12*s16*s24*s26**3*s34 +  \
    1024*q9a**3*s11**2*s12**2*s25*s26**3*s34 + 512*q9a**3*s11**2*s12*s14*s26**4*s34 - 1024*q8a*q9a**2*s11**4*s22**3*s34**2 + 1536*q8a*q9a**2*s11**3*s16*s22**2*s26*s34**2 - 1024*q8a*q9a**2*s11**3*s12*s22*s26**2*s34**2 -  \
    512*q8a*q9a**2*s11**2*s16**2*s22*s26**2*s34**2 - 512*q7a*q9a**2*s11**3*s22**2*s26**2*s34**2 + 512*q8a*q9a**2*s11**2*s12*s16*s26**3*s34**2 + 512*q7a*q9a**2*s11**2*s16*s22*s26**3*s34**2 - 512*q7a*q9a**2*s11**2*s12*s26**4*s34**2 +  \
    1024*q9a**3*s11**3*s15*s22**4*s35 - 512*q9a**3*s11**3*s16*s22**3*s24*s35 - 1024*q9a**3*s11**3*s12*s22**3*s25*s35 + 512*q9a**3*s11**2*s16**2*s22**3*s25*s35 - 512*q9a**3*s11**3*s14*s22**3*s26*s35 -  \
    1536*q9a**3*s11**2*s15*s16*s22**3*s26*s35 + 1024*q9a**3*s11**3*s12*s22**2*s24*s26*s35 + 512*q9a**3*s11**2*s16**2*s22**2*s24*s26*s35 + 512*q9a**3*s11**2*s12*s16*s22**2*s25*s26*s35 - 512*q9a**3*s11*s16**3*s22**2*s25*s26*s35 +  \
    1536*q9a**3*s11**2*s12*s15*s22**2*s26**2*s35 + 512*q9a**3*s11**2*s14*s16*s22**2*s26**2*s35 + 512*q9a**3*s11*s15*s16**2*s22**2*s26**2*s35 - 1536*q9a**3*s11**2*s12*s16*s22*s24*s26**2*s35 - 1024*q9a**3*s11**2*s12**2*s22*s25*s26**2*s35 +  \
    1024*q9a**3*s11*s12*s16**2*s22*s25*s26**2*s35 - 512*q9a**3*s11**2*s12*s14*s22*s26**3*s35 - 1024*q9a**3*s11*s12*s15*s16*s22*s26**3*s35 + 1024*q9a**3*s11**2*s12**2*s24*s26**3*s35 - 512*q9a**3*s11*s12**2*s16*s25*s26**3*s35 +  \
    512*q9a**3*s11*s12**2*s15*s26**4*s35 + 1024*q8a*q9a**2*s11**3*s16*s22**3*s34*s35 - 2048*q8a*q9a**2*s11**3*s12*s22**2*s26*s34*s35 - 1024*q8a*q9a**2*s11**2*s16**2*s22**2*s26*s34*s35 + 1024*q7a*q9a**2*s11**3*s22**3*s26*s34*s35 +  \
    3072*q8a*q9a**2*s11**2*s12*s16*s22*s26**2*s34*s35 - 1024*q7a*q9a**2*s11**2*s16*s22**2*s26**2*s34*s35 - 2048*q8a*q9a**2*s11**2*s12**2*s26**3*s34*s35 + 1024*q7a*q9a**2*s11**2*s12*s22*s26**3*s34*s35 +  \
    1024*q8a*q9a**2*s11**3*s12*s22**3*s35**2 - 512*q8a*q9a**2*s11**2*s16**2*s22**3*s35**2 - 1024*q7a*q9a**2*s11**3*s22**4*s35**2 - 512*q8a*q9a**2*s11**2*s12*s16*s22**2*s26*s35**2 + 512*q8a*q9a**2*s11*s16**3*s22**2*s26*s35**2 +  \
    1536*q7a*q9a**2*s11**2*s16*s22**3*s26*s35**2 + 1024*q8a*q9a**2*s11**2*s12**2*s22*s26**2*s35**2 - 1024*q8a*q9a**2*s11*s12*s16**2*s22*s26**2*s35**2 - 1536*q7a*q9a**2*s11**2*s12*s22**2*s26**2*s35**2 -  \
    512*q7a*q9a**2*s11*s16**2*s22**2*s26**2*s35**2 + 512*q8a*q9a**2*s11*s12**2*s16*s26**3*s35**2 + 1024*q7a*q9a**2*s11*s12*s16*s22*s26**3*s35**2 - 512*q7a*q9a**2*s11*s12**2*s26**4*s35**2 - 1024*q9a**3*s11**3*s16*s22**3*s23*s36 -  \
    512*q9a**3*s11**3*s15*s22**3*s24*s36 + 512*q9a**3*s11**3*s16*s22**2*s24**2*s36 - 1536*q9a**3*s11**3*s14*s22**3*s25*s36 + 1536*q9a**3*s11**2*s15*s16*s22**3*s25*s36 + 2048*q9a**3*s11**3*s12*s22**2*s24*s25*s36 -  \
    1024*q9a**3*s11**2*s16**2*s22**2*s24*s25*s36 - 1536*q9a**3*s11**2*s12*s16*s22**2*s25**2*s36 + 512*q9a**3*s11*s16**3*s22**2*s25**2*s36 - 1024*q9a**3*s11**3*s13*s22**3*s26*s36 + 512*q9a**3*s11**2*s15**2*s22**3*s26*s36 +  \
    2048*q9a**3*s11**3*s12*s22**2*s23*s26*s36 + 1024*q9a**3*s11**2*s16**2*s22**2*s23*s26*s36 + 1536*q9a**3*s11**3*s14*s22**2*s24*s26*s36 - 512*q9a**3*s11**2*s15*s16*s22**2*s24*s26*s36 - 2048*q9a**3*s11**3*s12*s22*s24**2*s26*s36 -  \
    2560*q9a**3*s11**2*s12*s15*s22**2*s25*s26*s36 + 512*q9a**3*s11**2*s14*s16*s22**2*s25*s26*s36 - 512*q9a**3*s11*s15*s16**2*s22**2*s25*s26*s36 + 2048*q9a**3*s11**2*s12*s16*s22*s24*s25*s26*s36 +  \
    2048*q9a**3*s11**2*s12**2*s22*s25**2*s26*s36 - 1024*q9a**3*s11*s12*s16**2*s22*s25**2*s26*s36 - 1024*q9a**3*s11**2*s14*s15*s22**2*s26**2*s36 + 1024*q9a**3*s11**2*s13*s16*s22**2*s26**2*s36 -  \
    3072*q9a**3*s11**2*s12*s16*s22*s23*s26**2*s36 + 1536*q9a**3*s11**2*s12*s15*s22*s24*s26**2*s36 - 512*q9a**3*s11**2*s14*s16*s22*s24*s26**2*s36 + 512*q9a**3*s11**2*s12*s16*s24**2*s26**2*s36 +  \
    512*q9a**3*s11**2*s12*s14*s22*s25*s26**2*s36 + 1024*q9a**3*s11*s12*s15*s16*s22*s25*s26**2*s36 - 2048*q9a**3*s11**2*s12**2*s24*s25*s26**2*s36 + 512*q9a**3*s11*s12**2*s16*s25**2*s26**2*s36 - 1024*q9a**3*s11**2*s12*s13*s22*s26**3*s36 +  \
    512*q9a**3*s11**2*s14**2*s22*s26**3*s36 + 2048*q9a**3*s11**2*s12**2*s23*s26**3*s36 - 512*q9a**3*s11**2*s12*s14*s24*s26**3*s36 - 512*q9a**3*s11*s12**2*s15*s25*s26**3*s36 - 3072*q8a*q9a**2*s11**3*s16*s22**3*s33*s36 +  \
    6144*q8a*q9a**2*s11**3*s12*s22**2*s26*s33*s36 + 3072*q8a*q9a**2*s11**2*s16**2*s22**2*s26*s33*s36 - 3072*q7a*q9a**2*s11**3*s22**3*s26*s33*s36 - 9216*q8a*q9a**2*s11**2*s12*s16*s22*s26**2*s33*s36 +  \
    3072*q7a*q9a**2*s11**2*s16*s22**2*s26**2*s33*s36 + 6144*q8a*q9a**2*s11**2*s12**2*s26**3*s33*s36 - 3072*q7a*q9a**2*s11**2*s12*s22*s26**3*s33*s36 + 1024*q8a*q9a**2*s11**3*s15*s22**3*s34*s36 +  \
    512*q8a*q9a**2*s11**3*s16*s22**2*s24*s34*s36 - 3072*q8a*q9a**2*s11**3*s12*s22**2*s25*s34*s36 + 512*q8a*q9a**2*s11**2*s16**2*s22**2*s25*s34*s36 + 2048*q7a*q9a**2*s11**3*s22**3*s25*s34*s36 -  \
    2560*q8a*q9a**2*s11**3*s14*s22**2*s26*s34*s36 - 512*q8a*q9a**2*s11**2*s15*s16*s22**2*s26*s34*s36 + 2048*q8a*q9a**2*s11**3*s12*s22*s24*s26*s34*s36 - 1024*q8a*q9a**2*s11**2*s16**2*s22*s24*s26*s34*s36 +  \
    1024*q8a*q9a**2*s11**2*s12*s16*s22*s25*s26*s34*s36 - 1536*q7a*q9a**2*s11**2*s16*s22**2*s25*s26*s34*s36 + 2048*q8a*q9a**2*s11**2*s14*s16*s22*s26**2*s34*s36 + 512*q7a*q9a**2*s11**2*s15*s22**2*s26**2*s34*s36 +  \
    512*q8a*q9a**2*s11**2*s12*s16*s24*s26**2*s34*s36 - 512*q7a*q9a**2*s11**2*s16*s22*s24*s26**2*s34*s36 - 1024*q8a*q9a**2*s11**2*s12**2*s25*s26**2*s34*s36 + 1024*q7a*q9a**2*s11**2*s12*s22*s25*s26**2*s34*s36 -  \
    1536*q8a*q9a**2*s11**2*s12*s14*s26**3*s34*s36 - 512*q7a*q9a**2*s11**2*s14*s22*s26**3*s34*s36 + 1024*q7a*q9a**2*s11**2*s12*s24*s26**3*s34*s36 - 1024*q8a**2*q9a*s11**3*s16*s22**2*s34**2*s36 +  \
    1024*q8a**2*q9a*s11**2*s16**2*s22*s26*s34**2*s36 + 1024*q7a*q8a*q9a*s11**3*s22**2*s26*s34**2*s36 - 1024*q8a**2*q9a*s11**2*s12*s16*s26**2*s34**2*s36 - 1024*q7a*q8a*q9a*s11**2*s16*s22*s26**2*s34**2*s36 +  \
    1024*q7a*q8a*q9a*s11**2*s12*s26**3*s34**2*s36 + 2048*q8a*q9a**2*s11**3*s14*s22**3*s35*s36 - 3072*q8a*q9a**2*s11**3*s12*s22**2*s24*s35*s36 + 512*q8a*q9a**2*s11**2*s16**2*s22**2*s24*s35*s36 +  \
    1024*q7a*q9a**2*s11**3*s22**3*s24*s35*s36 + 2560*q8a*q9a**2*s11**2*s12*s16*s22**2*s25*s35*s36 - 512*q8a*q9a**2*s11*s16**3*s22**2*s25*s35*s36 - 2560*q7a*q9a**2*s11**2*s16*s22**3*s25*s35*s36 -  \
    512*q8a*q9a**2*s11**2*s12*s15*s22**2*s26*s35*s36 - 1536*q8a*q9a**2*s11**2*s14*s16*s22**2*s26*s35*s36 - 512*q8a*q9a**2*s11*s15*s16**2*s22**2*s26*s35*s36 + 512*q7a*q9a**2*s11**2*s15*s22**3*s26*s35*s36 +  \
    1024*q8a*q9a**2*s11**2*s12*s16*s22*s24*s26*s35*s36 - 512*q7a*q9a**2*s11**2*s16*s22**2*s24*s26*s35*s36 - 2048*q8a*q9a**2*s11**2*s12**2*s22*s25*s26*s35*s36 + 2048*q7a*q9a**2*s11**2*s12*s22**2*s25*s26*s35*s36 +  \
    2048*q7a*q9a**2*s11*s16**2*s22**2*s25*s26*s35*s36 + 1024*q8a*q9a**2*s11**2*s12*s14*s22*s26**2*s35*s36 + 2048*q8a*q9a**2*s11*s12*s15*s16*s22*s26**2*s35*s36 + 512*q7a*q9a**2*s11**2*s14*s22**2*s26**2*s35*s36 -  \
    1024*q7a*q9a**2*s11*s15*s16*s22**2*s26**2*s35*s36 - 1024*q8a*q9a**2*s11**2*s12**2*s24*s26**2*s35*s36 + 512*q8a*q9a**2*s11*s12**2*s16*s25*s26**2*s35*s36 - 3072*q7a*q9a**2*s11*s12*s16*s22*s25*s26**2*s35*s36 -  \
    1536*q8a*q9a**2*s11*s12**2*s15*s26**3*s35*s36 + 1024*q7a*q9a**2*s11*s12*s15*s22*s26**3*s35*s36 + 1024*q7a*q9a**2*s11*s12**2*s25*s26**3*s35*s36 + 4096*q8a**2*q9a*s11**3*s12*s22**2*s34*s35*s36 -  \
    4096*q7a*q8a*q9a*s11**3*s22**3*s34*s35*s36 - 4096*q8a**2*q9a*s11**2*s12*s16*s22*s26*s34*s35*s36 + 4096*q7a*q8a*q9a*s11**2*s16*s22**2*s26*s34*s35*s36 + 4096*q8a**2*q9a*s11**2*s12**2*s26**2*s34*s35*s36 -  \
    4096*q7a*q8a*q9a*s11**2*s12*s22*s26**2*s34*s35*s36 - 1024*q8a**2*q9a*s11**2*s12*s16*s22**2*s35**2*s36 + 1024*q7a*q8a*q9a*s11**2*s16*s22**3*s35**2*s36 + 1024*q8a**2*q9a*s11*s12*s16**2*s22*s26*s35**2*s36 +  \
    1024*q7a*q8a*q9a*s11**2*s12*s22**2*s26*s35**2*s36 - 1024*q7a*q8a*q9a*s11*s16**2*s22**2*s26*s35**2*s36 - 1024*q7a**2*q9a*s11**2*s22**3*s26*s35**2*s36 - 1024*q8a**2*q9a*s11*s12**2*s16*s26**2*s35**2*s36 +  \
    1024*q7a**2*q9a*s11*s16*s22**2*s26**2*s35**2*s36 + 1024*q7a*q8a*q9a*s11*s12**2*s26**3*s35**2*s36 - 1024*q7a**2*q9a*s11*s12*s22*s26**3*s35**2*s36 + 1024*q8a*q9a**2*s11**3*s13*s22**3*s36**2 - 512*q8a*q9a**2*s11**2*s15**2*s22**3*s36**2 -  \
    2048*q8a*q9a**2*s11**3*s12*s22**2*s23*s36**2 - 1024*q8a*q9a**2*s11**2*s16**2*s22**2*s23*s36**2 + 1024*q7a*q9a**2*s11**3*s22**3*s23*s36**2 - 1536*q8a*q9a**2*s11**3*s14*s22**2*s24*s36**2 + 512*q8a*q9a**2*s11**2*s15*s16*s22**2*s24*s36**2 +  \
    2048*q8a*q9a**2*s11**3*s12*s22*s24**2*s36**2 - 512*q7a*q9a**2*s11**3*s22**2*s24**2*s36**2 + 2560*q8a*q9a**2*s11**2*s12*s15*s22**2*s25*s36**2 - 512*q8a*q9a**2*s11**2*s14*s16*s22**2*s25*s36**2 +  \
    512*q8a*q9a**2*s11*s15*s16**2*s22**2*s25*s36**2 - 1536*q7a*q9a**2*s11**2*s15*s22**3*s25*s36**2 - 2048*q8a*q9a**2*s11**2*s12*s16*s22*s24*s25*s36**2 + 2048*q7a*q9a**2*s11**2*s16*s22**2*s24*s25*s36**2 -  \
    2048*q8a*q9a**2*s11**2*s12**2*s22*s25**2*s36**2 + 1024*q8a*q9a**2*s11*s12*s16**2*s22*s25**2*s36**2 + 1536*q7a*q9a**2*s11**2*s12*s22**2*s25**2*s36**2 - 1536*q7a*q9a**2*s11*s16**2*s22**2*s25**2*s36**2 +  \
    2048*q8a*q9a**2*s11**2*s14*s15*s22**2*s26*s36**2 - 2048*q8a*q9a**2*s11**2*s13*s16*s22**2*s26*s36**2 + 6144*q8a*q9a**2*s11**2*s12*s16*s22*s23*s26*s36**2 - 2048*q7a*q9a**2*s11**2*s16*s22**2*s23*s26*s36**2 -  \
    3072*q8a*q9a**2*s11**2*s12*s15*s22*s24*s26*s36**2 + 1024*q8a*q9a**2*s11**2*s14*s16*s22*s24*s26*s36**2 + 512*q7a*q9a**2*s11**2*s15*s22**2*s24*s26*s36**2 - 1024*q8a*q9a**2*s11**2*s12*s16*s24**2*s26*s36**2 -  \
    1024*q8a*q9a**2*s11**2*s12*s14*s22*s25*s26*s36**2 - 2048*q8a*q9a**2*s11*s12*s15*s16*s22*s25*s26*s36**2 - 512*q7a*q9a**2*s11**2*s14*s22**2*s25*s26*s36**2 + 1024*q7a*q9a**2*s11*s15*s16*s22**2*s25*s26*s36**2 +  \
    4096*q8a*q9a**2*s11**2*s12**2*s24*s25*s26*s36**2 - 2048*q7a*q9a**2*s11**2*s12*s22*s24*s25*s26*s36**2 - 1024*q8a*q9a**2*s11*s12**2*s16*s25**2*s26*s36**2 + 2048*q7a*q9a**2*s11*s12*s16*s22*s25**2*s26*s36**2 +  \
    3072*q8a*q9a**2*s11**2*s12*s13*s22*s26**2*s36**2 - 1536*q8a*q9a**2*s11**2*s14**2*s22*s26**2*s36**2 - 1024*q7a*q9a**2*s11**2*s13*s22**2*s26**2*s36**2 - 6144*q8a*q9a**2*s11**2*s12**2*s23*s26**2*s36**2 +  \
    3072*q7a*q9a**2*s11**2*s12*s22*s23*s26**2*s36**2 + 1536*q8a*q9a**2*s11**2*s12*s14*s24*s26**2*s36**2 + 512*q7a*q9a**2*s11**2*s14*s22*s24*s26**2*s36**2 - 512*q7a*q9a**2*s11**2*s12*s24**2*s26**2*s36**2 +  \
    1536*q8a*q9a**2*s11*s12**2*s15*s25*s26**2*s36**2 - 1024*q7a*q9a**2*s11*s12*s15*s22*s25*s26**2*s36**2 - 512*q7a*q9a**2*s11*s12**2*s25**2*s26**2*s36**2 - 2048*q8a**2*q9a*s11**3*s12*s22**2*s33*s36**2 -  \
    1024*q8a**2*q9a*s11**2*s16**2*s22**2*s33*s36**2 + 2048*q7a*q8a*q9a*s11**3*s22**3*s33*s36**2 + 6144*q8a**2*q9a*s11**2*s12*s16*s22*s26*s33*s36**2 - 4096*q7a*q8a*q9a*s11**2*s16*s22**2*s26*s33*s36**2 -  \
    6144*q8a**2*q9a*s11**2*s12**2*s26**2*s33*s36**2 + 6144*q7a*q8a*q9a*s11**2*s12*s22*s26**2*s33*s36**2 - 1024*q7a**2*q9a*s11**2*s22**2*s26**2*s33*s36**2 + 2048*q8a**2*q9a*s11**3*s14*s22**2*s34*s36**2 -  \
    3072*q8a**2*q9a*s11**3*s12*s22*s24*s34*s36**2 + 512*q8a**2*q9a*s11**2*s16**2*s22*s24*s34*s36**2 + 1024*q7a*q8a*q9a*s11**3*s22**2*s24*s34*s36**2 + 512*q8a**2*q9a*s11**2*s12*s16*s22*s25*s34*s36**2 -  \
    512*q7a*q8a*q9a*s11**2*s16*s22**2*s25*s34*s36**2 + 1536*q8a**2*q9a*s11**2*s12*s15*s22*s26*s34*s36**2 - 2560*q8a**2*q9a*s11**2*s14*s16*s22*s26*s34*s36**2 - 1536*q7a*q8a*q9a*s11**2*s15*s22**2*s26*s34*s36**2

v1_1=  \
    512*q8a**2*q9a*s11**2*s12*s16*s24*s26*s34*s36**2 + 1024*q7a*q8a*q9a*s11**2*s16*s22*s24*s26*s34*s36**2 - 1024*q8a**2*q9a*s11**2*s12**2*s25*s26*s34*s36**2 + 1024*q7a**2*q9a*s11**2*s22**2*s25*s26*s34*s36**2 +  \
    1536*q8a**2*q9a*s11**2*s12*s14*s26**2*s34*s36**2 + 1024*q7a*q8a*q9a*s11**2*s14*s22*s26**2*s34*s36**2 - 2048*q7a*q8a*q9a*s11**2*s12*s24*s26**2*s34*s36**2 + 1024*q8a**3*s11**3*s12*s22*s34**2*s36**2 -  \
    512*q8a**3*s11**2*s16**2*s22*s34**2*s36**2 - 1024*q7a*q8a**2*s11**3*s22**2*s34**2*s36**2 + 512*q8a**3*s11**2*s12*s16*s26*s34**2*s36**2 + 512*q7a*q8a**2*s11**2*s16*s22*s26*s34**2*s36**2 - 512*q7a*q8a**2*s11**2*s12*s26**2*s34**2*s36**2 -  \
    1024*q8a**2*q9a*s11**2*s12*s15*s22**2*s35*s36**2 + 1024*q8a**2*q9a*s11**2*s14*s16*s22**2*s35*s36**2 + 1024*q7a*q8a*q9a*s11**2*s15*s22**3*s35*s36**2 + 512*q8a**2*q9a*s11**2*s12*s16*s22*s24*s35*s36**2 -  \
    1536*q7a*q8a*q9a*s11**2*s16*s22**2*s24*s35*s36**2 + 3072*q8a**2*q9a*s11**2*s12**2*s22*s25*s35*s36**2 - 1024*q8a**2*q9a*s11*s12*s16**2*s22*s25*s35*s36**2 - 5120*q7a*q8a*q9a*s11**2*s12*s22**2*s25*s35*s36**2 +  \
    1024*q7a*q8a*q9a*s11*s16**2*s22**2*s25*s35*s36**2 + 2048*q7a**2*q9a*s11**2*s22**3*s25*s35*s36**2 - 512*q8a**2*q9a*s11**2*s12*s14*s22*s26*s35*s36**2 - 1024*q8a**2*q9a*s11*s12*s15*s16*s22*s26*s35*s36**2 -  \
    512*q7a*q8a*q9a*s11**2*s14*s22**2*s26*s35*s36**2 + 1024*q7a*q8a*q9a*s11*s15*s16*s22**2*s26*s35*s36**2 - 1024*q8a**2*q9a*s11**2*s12**2*s24*s26*s35*s36**2 + 2048*q7a*q8a*q9a*s11**2*s12*s22*s24*s26*s35*s36**2 +  \
    512*q8a**2*q9a*s11*s12**2*s16*s25*s26*s35*s36**2 + 2048*q7a*q8a*q9a*s11*s12*s16*s22*s25*s26*s35*s36**2 - 2560*q7a**2*q9a*s11*s16*s22**2*s25*s26*s35*s36**2 + 1536*q8a**2*q9a*s11*s12**2*s15*s26**2*s35*s36**2 -  \
    2048*q7a*q8a*q9a*s11*s12*s15*s22*s26**2*s35*s36**2 + 512*q7a**2*q9a*s11*s15*s22**2*s26**2*s35*s36**2 - 2048*q7a*q8a*q9a*s11*s12**2*s25*s26**2*s35*s36**2 + 2048*q7a**2*q9a*s11*s12*s22*s25*s26**2*s35*s36**2 +  \
    1024*q8a**3*s11**2*s12*s16*s22*s34*s35*s36**2 - 1024*q7a*q8a**2*s11**2*s16*s22**2*s34*s35*s36**2 - 2048*q8a**3*s11**2*s12**2*s26*s34*s35*s36**2 + 3072*q7a*q8a**2*s11**2*s12*s22*s26*s34*s35*s36**2 -  \
    1024*q7a**2*q8a*s11**2*s22**2*s26*s34*s35*s36**2 - 1024*q8a**3*s11**2*s12**2*s22*s35**2*s36**2 + 2048*q7a*q8a**2*s11**2*s12*s22**2*s35**2*s36**2 - 1024*q7a**2*q8a*s11**2*s22**3*s35**2*s36**2 + 512*q8a**3*s11*s12**2*s16*s26*s35**2*s36**2 -  \
    1024*q7a*q8a**2*s11*s12*s16*s22*s26*s35**2*s36**2 + 512*q7a**2*q8a*s11*s16*s22**2*s26*s35**2*s36**2 - 512*q7a*q8a**2*s11*s12**2*s26**2*s35**2*s36**2 + 1024*q7a**2*q8a*s11*s12*s22*s26**2*s35**2*s36**2 -  \
    512*q7a**3*s11*s22**2*s26**2*s35**2*s36**2 - 1024*q8a**2*q9a*s11**2*s14*s15*s22**2*s36**3 + 1024*q8a**2*q9a*s11**2*s13*s16*s22**2*s36**3 - 3072*q8a**2*q9a*s11**2*s12*s16*s22*s23*s36**3 + 2048*q7a*q8a*q9a*s11**2*s16*s22**2*s23*s36**3 +  \
    1536*q8a**2*q9a*s11**2*s12*s15*s22*s24*s36**3 - 512*q8a**2*q9a*s11**2*s14*s16*s22*s24*s36**3 - 512*q7a*q8a*q9a*s11**2*s15*s22**2*s24*s36**3 + 512*q8a**2*q9a*s11**2*s12*s16*s24**2*s36**3 +  \
    512*q8a**2*q9a*s11**2*s12*s14*s22*s25*s36**3 + 1024*q8a**2*q9a*s11*s12*s15*s16*s22*s25*s36**3 + 512*q7a*q8a*q9a*s11**2*s14*s22**2*s25*s36**3 - 1024*q7a*q8a*q9a*s11*s15*s16*s22**2*s25*s36**3 -  \
    2048*q8a**2*q9a*s11**2*s12**2*s24*s25*s36**3 + 2048*q7a*q8a*q9a*s11**2*s12*s22*s24*s25*s36**3 - 1024*q7a**2*q9a*s11**2*s22**2*s24*s25*s36**3 + 512*q8a**2*q9a*s11*s12**2*s16*s25**2*s36**3 -  \
    2048*q7a*q8a*q9a*s11*s12*s16*s22*s25**2*s36**3 + 1536*q7a**2*q9a*s11*s16*s22**2*s25**2*s36**3 - 3072*q8a**2*q9a*s11**2*s12*s13*s22*s26*s36**3 + 1536*q8a**2*q9a*s11**2*s14**2*s22*s26*s36**3 +  \
    2048*q7a*q8a*q9a*s11**2*s13*s22**2*s26*s36**3 + 6144*q8a**2*q9a*s11**2*s12**2*s23*s26*s36**3 - 6144*q7a*q8a*q9a*s11**2*s12*s22*s23*s26*s36**3 + 1024*q7a**2*q9a*s11**2*s22**2*s23*s26*s36**3 -  \
    1536*q8a**2*q9a*s11**2*s12*s14*s24*s26*s36**3 - 1024*q7a*q8a*q9a*s11**2*s14*s22*s24*s26*s36**3 + 1024*q7a*q8a*q9a*s11**2*s12*s24**2*s26*s36**3 - 1536*q8a**2*q9a*s11*s12**2*s15*s25*s26*s36**3 +  \
    2048*q7a*q8a*q9a*s11*s12*s15*s22*s25*s26*s36**3 - 512*q7a**2*q9a*s11*s15*s22**2*s25*s26*s36**3 + 1024*q7a*q8a*q9a*s11*s12**2*s25**2*s26*s36**3 - 1024*q7a**2*q9a*s11*s12*s22*s25**2*s26*s36**3 -  \
    1024*q8a**3*s11**2*s12*s16*s22*s33*s36**3 + 1024*q7a*q8a**2*s11**2*s16*s22**2*s33*s36**3 + 2048*q8a**3*s11**2*s12**2*s26*s33*s36**3 - 3072*q7a*q8a**2*s11**2*s12*s22*s26*s33*s36**3 + 1024*q7a**2*q8a*s11**2*s22**2*s26*s33*s36**3 -  \
    1024*q8a**3*s11**2*s12*s15*s22*s34*s36**3 + 1024*q8a**3*s11**2*s14*s16*s22*s34*s36**3 + 1024*q7a*q8a**2*s11**2*s15*s22**2*s34*s36**3 - 512*q8a**3*s11**2*s12*s16*s24*s34*s36**3 - 512*q7a*q8a**2*s11**2*s16*s22*s24*s34*s36**3 +  \
    1024*q8a**3*s11**2*s12**2*s25*s34*s36**3 - 1024*q7a*q8a**2*s11**2*s12*s22*s25*s34*s36**3 - 512*q8a**3*s11**2*s12*s14*s26*s34*s36**3 - 512*q7a*q8a**2*s11**2*s14*s22*s26*s34*s36**3 + 1024*q7a*q8a**2*s11**2*s12*s24*s26*s34*s36**3 +  \
    1024*q8a**3*s11**2*s12**2*s24*s35*s36**3 - 2048*q7a*q8a**2*s11**2*s12*s22*s24*s35*s36**3 + 1024*q7a**2*q8a*s11**2*s22**2*s24*s35*s36**3 - 512*q8a**3*s11*s12**2*s16*s25*s35*s36**3 + 1024*q7a*q8a**2*s11*s12*s16*s22*s25*s35*s36**3 -  \
    512*q7a**2*q8a*s11*s16*s22**2*s25*s35*s36**3 - 512*q8a**3*s11*s12**2*s15*s26*s35*s36**3 + 1024*q7a*q8a**2*s11*s12*s15*s22*s26*s35*s36**3 - 512*q7a**2*q8a*s11*s15*s22**2*s26*s35*s36**3 +  \
    1024*q7a*q8a**2*s11*s12**2*s25*s26*s35*s36**3 - 2048*q7a**2*q8a*s11*s12*s22*s25*s26*s35*s36**3 + 1024*q7a**3*s11*s22**2*s25*s26*s35*s36**3 + 1024*q8a**3*s11**2*s12*s13*s22*s36**4 - 512*q8a**3*s11**2*s14**2*s22*s36**4 -  \
    1024*q7a*q8a**2*s11**2*s13*s22**2*s36**4 - 2048*q8a**3*s11**2*s12**2*s23*s36**4 + 3072*q7a*q8a**2*s11**2*s12*s22*s23*s36**4 - 1024*q7a**2*q8a*s11**2*s22**2*s23*s36**4 + 512*q8a**3*s11**2*s12*s14*s24*s36**4 +  \
    512*q7a*q8a**2*s11**2*s14*s22*s24*s36**4 - 512*q7a*q8a**2*s11**2*s12*s24**2*s36**4 + 512*q8a**3*s11*s12**2*s15*s25*s36**4 - 1024*q7a*q8a**2*s11*s12*s15*s22*s25*s36**4 + 512*q7a**2*q8a*s11*s15*s22**2*s25*s36**4 -  \
    512*q7a*q8a**2*s11*s12**2*s25**2*s36**4 + 1024*q7a**2*q8a*s11*s12*s22*s25**2*s36**4 - 512*q7a**3*s11*s22**2*s25**2*s36**4
v1=v1_0+v1_1

v2_0=3072*q9a**2*s11**4*s22**4*s33**2 - 6144*q9a**2*s11**3*s16*s22**3*s26*s33**2 + 6144*q9a**2*s11**3*s12*s22**2*s26**2*s33**2 + 3072*q9a**2*s11**2*s16**2*s22**2*s26**2*s33**2 - 6144*q9a**2*s11**2*s12*s16*s22*s26**3*s33**2 +  \
    3072*q9a**2*s11**2*s12**2*s26**4*s33**2 - 3072*q9a**2*s11**4*s22**3*s24*s33*s34 + 1536*q9a**2*s11**3*s16*s22**3*s25*s33*s34 + 1536*q9a**2*s11**3*s15*s22**3*s26*s33*s34 + 4608*q9a**2*s11**3*s16*s22**2*s24*s26*s33*s34 -  \
    3072*q9a**2*s11**3*s12*s22**2*s25*s26*s33*s34 - 1536*q9a**2*s11**2*s16**2*s22**2*s25*s26*s33*s34 - 1536*q9a**2*s11**3*s14*s22**2*s26**2*s33*s34 - 1536*q9a**2*s11**2*s15*s16*s22**2*s26**2*s33*s34 -  \
    3072*q9a**2*s11**3*s12*s22*s24*s26**2*s33*s34 - 1536*q9a**2*s11**2*s16**2*s22*s24*s26**2*s33*s34 + 4608*q9a**2*s11**2*s12*s16*s22*s25*s26**2*s33*s34 + 1536*q9a**2*s11**2*s12*s15*s22*s26**3*s33*s34 +  \
    1536*q9a**2*s11**2*s14*s16*s22*s26**3*s33*s34 + 1536*q9a**2*s11**2*s12*s16*s24*s26**3*s33*s34 - 3072*q9a**2*s11**2*s12**2*s25*s26**3*s33*s34 - 1536*q9a**2*s11**2*s12*s14*s26**4*s33*s34 + 1024*q9a**2*s11**4*s22**3*s23*s34**2 +  \
    512*q9a**2*s11**4*s22**2*s24**2*s34**2 - 512*q9a**2*s11**3*s15*s22**3*s25*s34**2 - 512*q9a**2*s11**3*s16*s22**2*s24*s25*s34**2 + 512*q9a**2*s11**3*s12*s22**2*s25**2*s34**2 - 1536*q9a**2*s11**3*s16*s22**2*s23*s26*s34**2 -  \
    512*q9a**2*s11**3*s15*s22**2*s24*s26*s34**2 - 512*q9a**2*s11**3*s16*s22*s24**2*s26*s34**2 + 1024*q9a**2*s11**3*s14*s22**2*s25*s26*s34**2 + 512*q9a**2*s11**2*s15*s16*s22**2*s25*s26*s34**2 +  \
    512*q9a**2*s11**2*s16**2*s22*s24*s25*s26*s34**2 - 512*q9a**2*s11**2*s12*s16*s22*s25**2*s26*s34**2 + 512*q9a**2*s11**3*s13*s22**2*s26**2*s34**2 + 1024*q9a**2*s11**3*s12*s22*s23*s26**2*s34**2 +  \
    512*q9a**2*s11**2*s16**2*s22*s23*s26**2*s34**2 + 512*q9a**2*s11**2*s15*s16*s22*s24*s26**2*s34**2 + 512*q9a**2*s11**3*s12*s24**2*s26**2*s34**2 - 512*q9a**2*s11**2*s12*s15*s22*s25*s26**2*s34**2 -  \
    1024*q9a**2*s11**2*s14*s16*s22*s25*s26**2*s34**2 - 512*q9a**2*s11**2*s12*s16*s24*s25*s26**2*s34**2 + 512*q9a**2*s11**2*s12**2*s25**2*s26**2*s34**2 - 512*q9a**2*s11**2*s13*s16*s22*s26**3*s34**2 -  \
    512*q9a**2*s11**2*s12*s16*s23*s26**3*s34**2 - 512*q9a**2*s11**2*s12*s15*s24*s26**3*s34**2 + 1024*q9a**2*s11**2*s12*s14*s25*s26**3*s34**2 + 512*q9a**2*s11**2*s12*s13*s26**4*s34**2 + 2048*q8a*q9a*s11**4*s22**3*s33*s34**2 -  \
    3072*q8a*q9a*s11**3*s16*s22**2*s26*s33*s34**2 + 2048*q8a*q9a*s11**3*s12*s22*s26**2*s33*s34**2 + 1024*q8a*q9a*s11**2*s16**2*s22*s26**2*s33*s34**2 + 1024*q7a*q9a*s11**3*s22**2*s26**2*s33*s34**2 -  \
    1024*q8a*q9a*s11**2*s12*s16*s26**3*s33*s34**2 - 1024*q7a*q9a*s11**2*s16*s22*s26**3*s33*s34**2 + 1024*q7a*q9a*s11**2*s12*s26**4*s33*s34**2 - 1024*q8a*q9a*s11**4*s22**2*s24*s34**3 + 512*q8a*q9a*s11**3*s16*s22**2*s25*s34**3 +  \
    512*q8a*q9a*s11**3*s15*s22**2*s26*s34**3 + 1024*q8a*q9a*s11**3*s16*s22*s24*s26*s34**3 - 512*q8a*q9a*s11**2*s16**2*s22*s25*s26*s34**3 - 1024*q7a*q9a*s11**3*s22**2*s25*s26*s34**3 - 512*q8a*q9a*s11**2*s15*s16*s22*s26**2*s34**3 -  \
    1024*q8a*q9a*s11**3*s12*s24*s26**2*s34**3 + 512*q8a*q9a*s11**2*s12*s16*s25*s26**2*s34**3 + 1024*q7a*q9a*s11**2*s16*s22*s25*s26**2*s34**3 + 512*q8a*q9a*s11**2*s12*s15*s26**3*s34**3 - 1024*q7a*q9a*s11**2*s12*s25*s26**3*s34**3 +  \
    512*q8a**2*s11**4*s22**2*s34**4 - 512*q8a**2*s11**3*s16*s22*s26*s34**4 + 512*q8a**2*s11**3*s12*s26**2*s34**4 - 3072*q9a**2*s11**3*s15*s22**4*s33*s35 + 1536*q9a**2*s11**3*s16*s22**3*s24*s33*s35 +  \
    3072*q9a**2*s11**3*s12*s22**3*s25*s33*s35 - 1536*q9a**2*s11**2*s16**2*s22**3*s25*s33*s35 + 1536*q9a**2*s11**3*s14*s22**3*s26*s33*s35 + 4608*q9a**2*s11**2*s15*s16*s22**3*s26*s33*s35 - 3072*q9a**2*s11**3*s12*s22**2*s24*s26*s33*s35 -  \
    1536*q9a**2*s11**2*s16**2*s22**2*s24*s26*s33*s35 - 1536*q9a**2*s11**2*s12*s16*s22**2*s25*s26*s33*s35 + 1536*q9a**2*s11*s16**3*s22**2*s25*s26*s33*s35 - 4608*q9a**2*s11**2*s12*s15*s22**2*s26**2*s33*s35 -  \
    1536*q9a**2*s11**2*s14*s16*s22**2*s26**2*s33*s35 - 1536*q9a**2*s11*s15*s16**2*s22**2*s26**2*s33*s35 + 4608*q9a**2*s11**2*s12*s16*s22*s24*s26**2*s33*s35 + 3072*q9a**2*s11**2*s12**2*s22*s25*s26**2*s33*s35 -  \
    3072*q9a**2*s11*s12*s16**2*s22*s25*s26**2*s33*s35 + 1536*q9a**2*s11**2*s12*s14*s22*s26**3*s33*s35 + 3072*q9a**2*s11*s12*s15*s16*s22*s26**3*s33*s35 - 3072*q9a**2*s11**2*s12**2*s24*s26**3*s33*s35 +  \
    1536*q9a**2*s11*s12**2*s16*s25*s26**3*s33*s35 - 1536*q9a**2*s11*s12**2*s15*s26**4*s33*s35 - 1024*q9a**2*s11**3*s16*s22**3*s23*s34*s35 + 1536*q9a**2*s11**3*s15*s22**3*s24*s34*s35 - 512*q9a**2*s11**3*s16*s22**2*s24**2*s34*s35 -  \
    1536*q9a**2*s11**3*s14*s22**3*s25*s34*s35 + 512*q9a**2*s11**2*s15*s16*s22**3*s25*s34*s35 + 512*q9a**2*s11**2*s16**2*s22**2*s24*s25*s34*s35 - 512*q9a**2*s11**2*s12*s16*s22**2*s25**2*s34*s35 -  \
    1024*q9a**2*s11**3*s13*s22**3*s26*s34*s35 - 512*q9a**2*s11**2*s15**2*s22**3*s26*s34*s35 + 2048*q9a**2*s11**3*s12*s22**2*s23*s26*s34*s35 + 1024*q9a**2*s11**2*s16**2*s22**2*s23*s26*s34*s35 +  \
    512*q9a**2*s11**3*s14*s22**2*s24*s26*s34*s35 - 2048*q9a**2*s11**2*s15*s16*s22**2*s24*s26*s34*s35 + 512*q9a**2*s11**2*s16**2*s22*s24**2*s26*s34*s35 + 512*q9a**2*s11**2*s12*s15*s22**2*s25*s26*s34*s35 +  \
    1024*q9a**2*s11**2*s14*s16*s22**2*s25*s26*s34*s35 - 512*q9a**2*s11*s15*s16**2*s22**2*s25*s26*s34*s35 - 512*q9a**2*s11*s16**3*s22*s24*s25*s26*s34*s35 + 512*q9a**2*s11*s12*s16**2*s22*s25**2*s26*s34*s35 +  \
    512*q9a**2*s11**2*s14*s15*s22**2*s26**2*s34*s35 + 1024*q9a**2*s11**2*s13*s16*s22**2*s26**2*s34*s35 + 512*q9a**2*s11*s15**2*s16*s22**2*s26**2*s34*s35 - 3072*q9a**2*s11**2*s12*s16*s22*s23*s26**2*s34*s35 +  \
    1536*q9a**2*s11**2*s12*s15*s22*s24*s26**2*s34*s35 - 512*q9a**2*s11**2*s14*s16*s22*s24*s26**2*s34*s35 + 512*q9a**2*s11*s15*s16**2*s22*s24*s26**2*s34*s35 - 512*q9a**2*s11**2*s12*s16*s24**2*s26**2*s34*s35 -  \
    1536*q9a**2*s11**2*s12*s14*s22*s25*s26**2*s34*s35 + 512*q9a**2*s11*s14*s16**2*s22*s25*s26**2*s34*s35 + 512*q9a**2*s11*s12*s16**2*s24*s25*s26**2*s34*s35 - 512*q9a**2*s11*s12**2*s16*s25**2*s26**2*s34*s35 -  \
    1024*q9a**2*s11**2*s12*s13*s22*s26**3*s34*s35 - 512*q9a**2*s11*s12*s15**2*s22*s26**3*s34*s35 - 512*q9a**2*s11*s14*s15*s16*s22*s26**3*s34*s35 + 2048*q9a**2*s11**2*s12**2*s23*s26**3*s34*s35 +  \
    512*q9a**2*s11**2*s12*s14*s24*s26**3*s34*s35 - 512*q9a**2*s11*s12*s15*s16*s24*s26**3*s34*s35 + 512*q9a**2*s11*s12**2*s15*s25*s26**3*s34*s35 - 512*q9a**2*s11*s12*s14*s16*s25*s26**3*s34*s35 +  \
    512*q9a**2*s11*s12*s14*s15*s26**4*s34*s35 - 2048*q8a*q9a*s11**3*s16*s22**3*s33*s34*s35 + 4096*q8a*q9a*s11**3*s12*s22**2*s26*s33*s34*s35 + 2048*q8a*q9a*s11**2*s16**2*s22**2*s26*s33*s34*s35 -  \
    2048*q7a*q9a*s11**3*s22**3*s26*s33*s34*s35 - 6144*q8a*q9a*s11**2*s12*s16*s22*s26**2*s33*s34*s35 + 2048*q7a*q9a*s11**2*s16*s22**2*s26**2*s33*s34*s35 + 4096*q8a*q9a*s11**2*s12**2*s26**3*s33*s34*s35 -  \
    2048*q7a*q9a*s11**2*s12*s22*s26**3*s33*s34*s35 - 1024*q8a*q9a*s11**3*s15*s22**3*s34**2*s35 + 1536*q8a*q9a*s11**3*s16*s22**2*s24*s34**2*s35 - 1024*q8a*q9a*s11**3*s12*s22**2*s25*s34**2*s35 -  \
    512*q8a*q9a*s11**2*s16**2*s22**2*s25*s34**2*s35 + 2048*q7a*q9a*s11**3*s22**3*s25*s34**2*s35 - 1536*q8a*q9a*s11**3*s14*s22**2*s26*s34**2*s35 + 1536*q8a*q9a*s11**2*s15*s16*s22**2*s26*s34**2*s35 -  \
    1536*q8a*q9a*s11**2*s16**2*s22*s24*s26*s34**2*s35 + 1024*q8a*q9a*s11**2*s12*s16*s22*s25*s26*s34**2*s35 + 512*q8a*q9a*s11*s16**3*s22*s25*s26*s34**2*s35 - 1536*q7a*q9a*s11**2*s16*s22**2*s25*s26*s34**2*s35 -  \
    1024*q8a*q9a*s11**2*s12*s15*s22*s26**2*s34**2*s35 + 1536*q8a*q9a*s11**2*s14*s16*s22*s26**2*s34**2*s35 - 512*q8a*q9a*s11*s15*s16**2*s22*s26**2*s34**2*s35 - 512*q7a*q9a*s11**2*s15*s22**2*s26**2*s34**2*s35 +  \
    1536*q8a*q9a*s11**2*s12*s16*s24*s26**2*s34**2*s35 - 1024*q8a*q9a*s11**2*s12**2*s25*s26**2*s34**2*s35 - 512*q8a*q9a*s11*s12*s16**2*s25*s26**2*s34**2*s35 + 2048*q7a*q9a*s11**2*s12*s22*s25*s26**2*s34**2*s35 -  \
    512*q7a*q9a*s11*s16**2*s22*s25*s26**2*s34**2*s35 - 1536*q8a*q9a*s11**2*s12*s14*s26**3*s34**2*s35 + 512*q8a*q9a*s11*s12*s15*s16*s26**3*s34**2*s35 + 512*q7a*q9a*s11*s15*s16*s22*s26**3*s34**2*s35 +  \
    512*q7a*q9a*s11*s12*s16*s25*s26**3*s34**2*s35 - 512*q7a*q9a*s11*s12*s15*s26**4*s34**2*s35 - 1024*q8a**2*s11**3*s16*s22**2*s34**3*s35 + 1024*q8a**2*s11**2*s16**2*s22*s26*s34**3*s35 + 1024*q7a*q8a*s11**3*s22**2*s26*s34**3*s35 -  \
    1024*q8a**2*s11**2*s12*s16*s26**2*s34**3*s35 - 1024*q7a*q8a*s11**2*s16*s22*s26**2*s34**3*s35 + 1024*q7a*q8a*s11**2*s12*s26**3*s34**3*s35 + 1024*q9a**2*s11**3*s13*s22**4*s35**2 + 512*q9a**2*s11**2*s15**2*s22**4*s35**2 -  \
    1024*q9a**2*s11**3*s12*s22**3*s23*s35**2 + 512*q9a**2*s11**2*s16**2*s22**3*s23*s35**2 - 512*q9a**2*s11**3*s14*s22**3*s24*s35**2 - 512*q9a**2*s11**2*s15*s16*s22**3*s24*s35**2 + 512*q9a**2*s11**3*s12*s22**2*s24**2*s35**2 -  \
    1024*q9a**2*s11**2*s12*s15*s22**3*s25*s35**2 + 1024*q9a**2*s11**2*s14*s16*s22**3*s25*s35**2 - 512*q9a**2*s11**2*s12*s16*s22**2*s24*s25*s35**2 + 512*q9a**2*s11**2*s12**2*s22**2*s25**2*s35**2 - 512*q9a**2*s11**2*s14*s15*s22**3*s26*s35**2 -  \
    1536*q9a**2*s11**2*s13*s16*s22**3*s26*s35**2 - 512*q9a**2*s11*s15**2*s16*s22**3*s26*s35**2 + 512*q9a**2*s11**2*s12*s16*s22**2*s23*s26*s35**2 - 512*q9a**2*s11*s16**3*s22**2*s23*s26*s35**2 +  \
    1024*q9a**2*s11**2*s12*s15*s22**2*s24*s26*s35**2 + 512*q9a**2*s11**2*s14*s16*s22**2*s24*s26*s35**2 + 512*q9a**2*s11*s15*s16**2*s22**2*s24*s26*s35**2 - 512*q9a**2*s11**2*s12*s16*s22*s24**2*s26*s35**2 -  \
    512*q9a**2*s11**2*s12*s14*s22**2*s25*s26*s35**2 + 1024*q9a**2*s11*s12*s15*s16*s22**2*s25*s26*s35**2 - 1024*q9a**2*s11*s14*s16**2*s22**2*s25*s26*s35**2 + 512*q9a**2*s11*s12*s16**2*s22*s24*s25*s26*s35**2 -  \
    512*q9a**2*s11*s12**2*s16*s22*s25**2*s26*s35**2 + 1536*q9a**2*s11**2*s12*s13*s22**2*s26**2*s35**2 + 512*q9a**2*s11*s12*s15**2*s22**2*s26**2*s35**2 + 512*q9a**2*s11*s14*s15*s16*s22**2*s26**2*s35**2 +  \
    512*q9a**2*s11*s13*s16**2*s22**2*s26**2*s35**2 - 1024*q9a**2*s11**2*s12**2*s22*s23*s26**2*s35**2 + 1024*q9a**2*s11*s12*s16**2*s22*s23*s26**2*s35**2 - 512*q9a**2*s11**2*s12*s14*s22*s24*s26**2*s35**2 -  \
    1536*q9a**2*s11*s12*s15*s16*s22*s24*s26**2*s35**2 + 512*q9a**2*s11**2*s12**2*s24**2*s26**2*s35**2 - 1024*q9a**2*s11*s12**2*s15*s22*s25*s26**2*s35**2 + 1536*q9a**2*s11*s12*s14*s16*s22*s25*s26**2*s35**2 -  \
    512*q9a**2*s11*s12**2*s16*s24*s25*s26**2*s35**2 + 512*q9a**2*s11*s12**3*s25**2*s26**2*s35**2 - 512*q9a**2*s11*s12*s14*s15*s22*s26**3*s35**2 - 1024*q9a**2*s11*s12*s13*s16*s22*s26**3*s35**2 - 512*q9a**2*s11*s12**2*s16*s23*s26**3*s35**2 +  \
    1024*q9a**2*s11*s12**2*s15*s24*s26**3*s35**2 - 512*q9a**2*s11*s12**2*s14*s25*s26**3*s35**2 + 512*q9a**2*s11*s12**2*s13*s26**4*s35**2 - 2048*q8a*q9a*s11**3*s12*s22**3*s33*s35**2 + 1024*q8a*q9a*s11**2*s16**2*s22**3*s33*s35**2 +  \
    2048*q7a*q9a*s11**3*s22**4*s33*s35**2 + 1024*q8a*q9a*s11**2*s12*s16*s22**2*s26*s33*s35**2 - 1024*q8a*q9a*s11*s16**3*s22**2*s26*s33*s35**2 - 3072*q7a*q9a*s11**2*s16*s22**3*s26*s33*s35**2 -  \
    2048*q8a*q9a*s11**2*s12**2*s22*s26**2*s33*s35**2 + 2048*q8a*q9a*s11*s12*s16**2*s22*s26**2*s33*s35**2 + 3072*q7a*q9a*s11**2*s12*s22**2*s26**2*s33*s35**2 + 1024*q7a*q9a*s11*s16**2*s22**2*s26**2*s33*s35**2 -  \
    1024*q8a*q9a*s11*s12**2*s16*s26**3*s33*s35**2 - 2048*q7a*q9a*s11*s12*s16*s22*s26**3*s33*s35**2 + 1024*q7a*q9a*s11*s12**2*s26**4*s33*s35**2 + 2048*q8a*q9a*s11**3*s14*s22**3*s34*s35**2 -  \
    1024*q8a*q9a*s11**3*s12*s22**2*s24*s34*s35**2 - 512*q8a*q9a*s11**2*s16**2*s22**2*s24*s34*s35**2 - 1024*q7a*q9a*s11**3*s22**3*s24*s34*s35**2 + 1536*q8a*q9a*s11**2*s12*s16*s22**2*s25*s34*s35**2 -  \
    1536*q7a*q9a*s11**2*s16*s22**3*s25*s34*s35**2 - 1536*q8a*q9a*s11**2*s12*s15*s22**2*s26*s34*s35**2 - 1536*q8a*q9a*s11**2*s14*s16*s22**2*s26*s34*s35**2 + 1536*q7a*q9a*s11**2*s15*s22**3*s26*s34*s35**2 +  \
    1024*q8a*q9a*s11**2*s12*s16*s22*s24*s26*s34*s35**2 + 512*q8a*q9a*s11*s16**3*s22*s24*s26*s34*s35**2 + 1536*q7a*q9a*s11**2*s16*s22**2*s24*s26*s34*s35**2 - 1536*q8a*q9a*s11*s12*s16**2*s22*s25*s26*s34*s35**2
v2_1= \
    1536*q7a*q9a*s11*s16**2*s22**2*s25*s26*s34*s35**2 + 2048*q8a*q9a*s11**2*s12*s14*s22*s26**2*s34*s35**2 + 1536*q8a*q9a*s11*s12*s15*s16*s22*s26**2*s34*s35**2 - 512*q8a*q9a*s11*s14*s16**2*s22*s26**2*s34*s35**2 -  \
    512*q7a*q9a*s11**2*s14*s22**2*s26**2*s34*s35**2 - 1536*q7a*q9a*s11*s15*s16*s22**2*s26**2*s34*s35**2 - 1024*q8a*q9a*s11**2*s12**2*s24*s26**2*s34*s35**2 - 512*q8a*q9a*s11*s12*s16**2*s24*s26**2*s34*s35**2 -  \
    1024*q7a*q9a*s11**2*s12*s22*s24*s26**2*s34*s35**2 - 512*q7a*q9a*s11*s16**2*s22*s24*s26**2*s34*s35**2 + 1536*q8a*q9a*s11*s12**2*s16*s25*s26**2*s34*s35**2 - 1536*q7a*q9a*s11*s12*s16*s22*s25*s26**2*s34*s35**2 -  \
    1536*q8a*q9a*s11*s12**2*s15*s26**3*s34*s35**2 + 512*q8a*q9a*s11*s12*s14*s16*s26**3*s34*s35**2 + 1536*q7a*q9a*s11*s12*s15*s22*s26**3*s34*s35**2 + 512*q7a*q9a*s11*s14*s16*s22*s26**3*s34*s35**2 +  \
    512*q7a*q9a*s11*s12*s16*s24*s26**3*s34*s35**2 - 512*q7a*q9a*s11*s12*s14*s26**4*s34*s35**2 + 1024*q8a**2*s11**3*s12*s22**2*s34**2*s35**2 + 512*q8a**2*s11**2*s16**2*s22**2*s34**2*s35**2 - 1024*q7a*q8a*s11**3*s22**3*s34**2*s35**2 -  \
    1024*q8a**2*s11**2*s12*s16*s22*s26*s34**2*s35**2 - 512*q8a**2*s11*s16**3*s22*s26*s34**2*s35**2 + 1024*q8a**2*s11**2*s12**2*s26**2*s34**2*s35**2 + 512*q8a**2*s11*s12*s16**2*s26**2*s34**2*s35**2 -  \
    1024*q7a*q8a*s11**2*s12*s22*s26**2*s34**2*s35**2 + 1024*q7a*q8a*s11*s16**2*s22*s26**2*s34**2*s35**2 + 512*q7a**2*s11**2*s22**2*s26**2*s34**2*s35**2 - 1024*q7a*q8a*s11*s12*s16*s26**3*s34**2*s35**2 -  \
    512*q7a**2*s11*s16*s22*s26**3*s34**2*s35**2 + 512*q7a**2*s11*s12*s26**4*s34**2*s35**2 + 1024*q8a*q9a*s11**2*s12*s15*s22**3*s35**3 - 1024*q8a*q9a*s11**2*s14*s16*s22**3*s35**3 - 1024*q7a*q9a*s11**2*s15*s22**4*s35**3 +  \
    512*q8a*q9a*s11**2*s12*s16*s22**2*s24*s35**3 + 512*q7a*q9a*s11**2*s16*s22**3*s24*s35**3 - 1024*q8a*q9a*s11**2*s12**2*s22**2*s25*s35**3 + 1024*q7a*q9a*s11**2*s12*s22**3*s25*s35**3 + 512*q8a*q9a*s11**2*s12*s14*s22**2*s26*s35**3 -  \
    1024*q8a*q9a*s11*s12*s15*s16*s22**2*s26*s35**3 + 1024*q8a*q9a*s11*s14*s16**2*s22**2*s26*s35**3 + 512*q7a*q9a*s11**2*s14*s22**3*s26*s35**3 + 1024*q7a*q9a*s11*s15*s16*s22**3*s26*s35**3 -  \
    512*q8a*q9a*s11*s12*s16**2*s22*s24*s26*s35**3 - 1024*q7a*q9a*s11**2*s12*s22**2*s24*s26*s35**3 - 512*q7a*q9a*s11*s16**2*s22**2*s24*s26*s35**3 + 1024*q8a*q9a*s11*s12**2*s16*s22*s25*s26*s35**3 -  \
    1024*q7a*q9a*s11*s12*s16*s22**2*s25*s26*s35**3 + 1024*q8a*q9a*s11*s12**2*s15*s22*s26**2*s35**3 - 1536*q8a*q9a*s11*s12*s14*s16*s22*s26**2*s35**3 - 1024*q7a*q9a*s11*s12*s15*s22**2*s26**2*s35**3 -  \
    512*q7a*q9a*s11*s14*s16*s22**2*s26**2*s35**3 + 512*q8a*q9a*s11*s12**2*s16*s24*s26**2*s35**3 + 1536*q7a*q9a*s11*s12*s16*s22*s24*s26**2*s35**3 - 1024*q8a*q9a*s11*s12**3*s25*s26**2*s35**3 +  \
    1024*q7a*q9a*s11*s12**2*s22*s25*s26**2*s35**3 + 512*q8a*q9a*s11*s12**2*s14*s26**3*s35**3 + 512*q7a*q9a*s11*s12*s14*s22*s26**3*s35**3 - 1024*q7a*q9a*s11*s12**2*s24*s26**3*s35**3 - 1024*q8a**2*s11**2*s12*s16*s22**2*s34*s35**3 +  \
    1024*q7a*q8a*s11**2*s16*s22**3*s34*s35**3 + 1024*q8a**2*s11*s12*s16**2*s22*s26*s34*s35**3 + 1024*q7a*q8a*s11**2*s12*s22**2*s26*s34*s35**3 - 1024*q7a*q8a*s11*s16**2*s22**2*s26*s34*s35**3 - 1024*q7a**2*s11**2*s22**3*s26*s34*s35**3 -  \
    1024*q8a**2*s11*s12**2*s16*s26**2*s34*s35**3 + 1024*q7a**2*s11*s16*s22**2*s26**2*s34*s35**3 + 1024*q7a*q8a*s11*s12**2*s26**3*s34*s35**3 - 1024*q7a**2*s11*s12*s22*s26**3*s34*s35**3 + 512*q8a**2*s11**2*s12**2*s22**2*s35**4 -  \
    1024*q7a*q8a*s11**2*s12*s22**3*s35**4 + 512*q7a**2*s11**2*s22**4*s35**4 - 512*q8a**2*s11*s12**2*s16*s22*s26*s35**4 + 1024*q7a*q8a*s11*s12*s16*s22**2*s26*s35**4 - 512*q7a**2*s11*s16*s22**3*s26*s35**4 +  \
    512*q8a**2*s11*s12**3*s26**2*s35**4 - 1024*q7a*q8a*s11*s12**2*s22*s26**2*s35**4 + 512*q7a**2*s11*s12*s22**2*s26**2*s35**4 + 3072*q9a**2*s11**3*s16*s22**3*s23*s33*s36 + 1536*q9a**2*s11**3*s15*s22**3*s24*s33*s36 -  \
    1536*q9a**2*s11**3*s16*s22**2*s24**2*s33*s36 + 4608*q9a**2*s11**3*s14*s22**3*s25*s33*s36 - 4608*q9a**2*s11**2*s15*s16*s22**3*s25*s33*s36 - 6144*q9a**2*s11**3*s12*s22**2*s24*s25*s33*s36 +  \
    3072*q9a**2*s11**2*s16**2*s22**2*s24*s25*s33*s36 + 4608*q9a**2*s11**2*s12*s16*s22**2*s25**2*s33*s36 - 1536*q9a**2*s11*s16**3*s22**2*s25**2*s33*s36 + 3072*q9a**2*s11**3*s13*s22**3*s26*s33*s36 -  \
    1536*q9a**2*s11**2*s15**2*s22**3*s26*s33*s36 - 6144*q9a**2*s11**3*s12*s22**2*s23*s26*s33*s36 - 3072*q9a**2*s11**2*s16**2*s22**2*s23*s26*s33*s36 - 4608*q9a**2*s11**3*s14*s22**2*s24*s26*s33*s36 +  \
    1536*q9a**2*s11**2*s15*s16*s22**2*s24*s26*s33*s36 + 6144*q9a**2*s11**3*s12*s22*s24**2*s26*s33*s36 + 7680*q9a**2*s11**2*s12*s15*s22**2*s25*s26*s33*s36 - 1536*q9a**2*s11**2*s14*s16*s22**2*s25*s26*s33*s36 +  \
    1536*q9a**2*s11*s15*s16**2*s22**2*s25*s26*s33*s36 - 6144*q9a**2*s11**2*s12*s16*s22*s24*s25*s26*s33*s36 - 6144*q9a**2*s11**2*s12**2*s22*s25**2*s26*s33*s36 + 3072*q9a**2*s11*s12*s16**2*s22*s25**2*s26*s33*s36 +  \
    3072*q9a**2*s11**2*s14*s15*s22**2*s26**2*s33*s36 - 3072*q9a**2*s11**2*s13*s16*s22**2*s26**2*s33*s36 + 9216*q9a**2*s11**2*s12*s16*s22*s23*s26**2*s33*s36 - 4608*q9a**2*s11**2*s12*s15*s22*s24*s26**2*s33*s36 +  \
    1536*q9a**2*s11**2*s14*s16*s22*s24*s26**2*s33*s36 - 1536*q9a**2*s11**2*s12*s16*s24**2*s26**2*s33*s36 - 1536*q9a**2*s11**2*s12*s14*s22*s25*s26**2*s33*s36 - 3072*q9a**2*s11*s12*s15*s16*s22*s25*s26**2*s33*s36 +  \
    6144*q9a**2*s11**2*s12**2*s24*s25*s26**2*s33*s36 - 1536*q9a**2*s11*s12**2*s16*s25**2*s26**2*s33*s36 + 3072*q9a**2*s11**2*s12*s13*s22*s26**3*s33*s36 - 1536*q9a**2*s11**2*s14**2*s22*s26**3*s33*s36 -  \
    6144*q9a**2*s11**2*s12**2*s23*s26**3*s33*s36 + 1536*q9a**2*s11**2*s12*s14*s24*s26**3*s33*s36 + 1536*q9a**2*s11*s12**2*s15*s25*s26**3*s33*s36 + 3072*q8a*q9a*s11**3*s16*s22**3*s33**2*s36 -  \
    6144*q8a*q9a*s11**3*s12*s22**2*s26*s33**2*s36 - 3072*q8a*q9a*s11**2*s16**2*s22**2*s26*s33**2*s36 + 3072*q7a*q9a*s11**3*s22**3*s26*s33**2*s36 + 9216*q8a*q9a*s11**2*s12*s16*s22*s26**2*s33**2*s36 -  \
    3072*q7a*q9a*s11**2*s16*s22**2*s26**2*s33**2*s36 - 6144*q8a*q9a*s11**2*s12**2*s26**3*s33**2*s36 + 3072*q7a*q9a*s11**2*s12*s22*s26**3*s33**2*s36 - 1024*q9a**2*s11**3*s15*s22**3*s23*s34*s36 -  \
    512*q9a**2*s11**3*s16*s22**2*s23*s24*s34*s36 - 512*q9a**2*s11**3*s15*s22**2*s24**2*s34*s36 + 512*q9a**2*s11**3*s16*s22*s24**3*s34*s36 - 2048*q9a**2*s11**3*s13*s22**3*s25*s34*s36 + 1024*q9a**2*s11**2*s15**2*s22**3*s25*s34*s36 +  \
    3072*q9a**2*s11**3*s12*s22**2*s23*s25*s34*s36 - 512*q9a**2*s11**2*s16**2*s22**2*s23*s25*s34*s36 - 512*q9a**2*s11**3*s14*s22**2*s24*s25*s34*s36 + 1024*q9a**2*s11**2*s15*s16*s22**2*s24*s25*s34*s36 +  \
    1024*q9a**2*s11**3*s12*s22*s24**2*s25*s34*s36 - 1024*q9a**2*s11**2*s16**2*s22*s24**2*s25*s34*s36 - 2048*q9a**2*s11**2*s12*s15*s22**2*s25**2*s34*s36 - 512*q9a**2*s11**2*s14*s16*s22**2*s25**2*s34*s36 +  \
    512*q9a**2*s11*s15*s16**2*s22**2*s25**2*s34*s36 - 512*q9a**2*s11**2*s12*s16*s22*s24*s25**2*s34*s36 + 512*q9a**2*s11*s16**3*s22*s24*s25**2*s34*s36 + 1024*q9a**2*s11**2*s12**2*s22*s25**3*s34*s36 -  \
    512*q9a**2*s11*s12*s16**2*s22*s25**3*s34*s36 + 2560*q9a**2*s11**3*s14*s22**2*s23*s26*s34*s36 + 512*q9a**2*s11**2*s15*s16*s22**2*s23*s26*s34*s36 + 512*q9a**2*s11**2*s15**2*s22**2*s24*s26*s34*s36 -  \
    2048*q9a**2*s11**3*s12*s22*s23*s24*s26*s34*s36 + 1024*q9a**2*s11**2*s16**2*s22*s23*s24*s26*s34*s36 + 512*q9a**2*s11**3*s14*s22*s24**2*s26*s34*s36 - 512*q9a**2*s11**2*s15*s16*s22*s24**2*s26*s34*s36 -  \
    1024*q9a**2*s11**3*s12*s24**3*s26*s34*s36 - 2048*q9a**2*s11**2*s14*s15*s22**2*s25*s26*s34*s36 + 1536*q9a**2*s11**2*s13*s16*s22**2*s25*s26*s34*s36 - 512*q9a**2*s11*s15**2*s16*s22**2*s25*s26*s34*s36 -  \
    1024*q9a**2*s11**2*s12*s16*s22*s23*s25*s26*s34*s36 + 1024*q9a**2*s11**2*s14*s16*s22*s24*s25*s26*s34*s36 - 512*q9a**2*s11*s15*s16**2*s22*s24*s25*s26*s34*s36 + 1536*q9a**2*s11**2*s12*s16*s24**2*s25*s26*s34*s36 +  \
    2560*q9a**2*s11**2*s12*s14*s22*s25**2*s26*s34*s36 - 512*q9a**2*s11*s14*s16**2*s22*s25**2*s26*s34*s36 - 1024*q9a**2*s11**2*s12**2*s24*s25**2*s26*s34*s36 - 512*q9a**2*s11*s12*s16**2*s24*s25**2*s26*s34*s36 +  \
    512*q9a**2*s11*s12**2*s16*s25**3*s26*s34*s36 - 512*q9a**2*s11**2*s13*s15*s22**2*s26**2*s34*s36 - 2048*q9a**2*s11**2*s14*s16*s22*s23*s26**2*s34*s36 - 512*q9a**2*s11**2*s14*s15*s22*s24*s26**2*s34*s36 +  \
    512*q9a**2*s11**2*s13*s16*s22*s24*s26**2*s34*s36 - 512*q9a**2*s11**2*s12*s16*s23*s24*s26**2*s34*s36 + 1024*q9a**2*s11**2*s12*s15*s24**2*s26**2*s34*s36 - 1024*q9a**2*s11**2*s12*s13*s22*s25*s26**2*s34*s36 +  \
    1024*q9a**2*s11**2*s14**2*s22*s25*s26**2*s34*s36 + 512*q9a**2*s11*s12*s15**2*s22*s25*s26**2*s34*s36 + 512*q9a**2*s11*s14*s15*s16*s22*s25*s26**2*s34*s36 + 1024*q9a**2*s11**2*s12**2*s23*s25*s26**2*s34*s36 -  \
    2560*q9a**2*s11**2*s12*s14*s24*s25*s26**2*s34*s36 + 512*q9a**2*s11*s12*s15*s16*s24*s25*s26**2*s34*s36 - 512*q9a**2*s11*s12**2*s15*s25**2*s26**2*s34*s36 + 512*q9a**2*s11*s12*s14*s16*s25**2*s26**2*s34*s36 +  \
    512*q9a**2*s11**2*s13*s14*s22*s26**3*s34*s36 + 1536*q9a**2*s11**2*s12*s14*s23*s26**3*s34*s36 - 1024*q9a**2*s11**2*s12*s13*s24*s26**3*s34*s36 - 512*q9a**2*s11*s12*s14*s15*s25*s26**3*s34*s36 -  \
    2048*q8a*q9a*s11**3*s15*s22**3*s33*s34*s36 - 1024*q8a*q9a*s11**3*s16*s22**2*s24*s33*s34*s36 + 6144*q8a*q9a*s11**3*s12*s22**2*s25*s33*s34*s36 - 1024*q8a*q9a*s11**2*s16**2*s22**2*s25*s33*s34*s36 -  \
    4096*q7a*q9a*s11**3*s22**3*s25*s33*s34*s36 + 5120*q8a*q9a*s11**3*s14*s22**2*s26*s33*s34*s36 + 1024*q8a*q9a*s11**2*s15*s16*s22**2*s26*s33*s34*s36 - 4096*q8a*q9a*s11**3*s12*s22*s24*s26*s33*s34*s36 +  \
    2048*q8a*q9a*s11**2*s16**2*s22*s24*s26*s33*s34*s36 - 2048*q8a*q9a*s11**2*s12*s16*s22*s25*s26*s33*s34*s36 + 3072*q7a*q9a*s11**2*s16*s22**2*s25*s26*s33*s34*s36 - 4096*q8a*q9a*s11**2*s14*s16*s22*s26**2*s33*s34*s36 -  \
    1024*q7a*q9a*s11**2*s15*s22**2*s26**2*s33*s34*s36 - 1024*q8a*q9a*s11**2*s12*s16*s24*s26**2*s33*s34*s36 + 1024*q7a*q9a*s11**2*s16*s22*s24*s26**2*s33*s34*s36 + 2048*q8a*q9a*s11**2*s12**2*s25*s26**2*s33*s34*s36 -  \
    2048*q7a*q9a*s11**2*s12*s22*s25*s26**2*s33*s34*s36 + 3072*q8a*q9a*s11**2*s12*s14*s26**3*s33*s34*s36 + 1024*q7a*q9a*s11**2*s14*s22*s26**3*s33*s34*s36 - 2048*q7a*q9a*s11**2*s12*s24*s26**3*s33*s34*s36 +  \
    2048*q8a*q9a*s11**3*s16*s22**2*s23*s34**2*s36 + 1536*q8a*q9a*s11**3*s15*s22**2*s24*s34**2*s36 - 1024*q8a*q9a*s11**3*s16*s22*s24**2*s34**2*s36 - 512*q8a*q9a*s11**3*s14*s22**2*s25*s34**2*s36 -  \
    1536*q8a*q9a*s11**2*s15*s16*s22**2*s25*s34**2*s36 - 2048*q8a*q9a*s11**3*s12*s22*s24*s25*s34**2*s36 + 1536*q8a*q9a*s11**2*s16**2*s22*s24*s25*s34**2*s36 + 1024*q7a*q9a*s11**3*s22**2*s24*s25*s34**2*s36 +  \
    1024*q8a*q9a*s11**2*s12*s16*s22*s25**2*s34**2*s36 - 512*q8a*q9a*s11*s16**3*s22*s25**2*s34**2*s36 + 512*q7a*q9a*s11**2*s16*s22**2*s25**2*s34**2*s36 - 1024*q8a*q9a*s11**3*s13*s22**2*s26*s34**2*s36 -  \
    512*q8a*q9a*s11**2*s15**2*s22**2*s26*s34**2*s36 - 2048*q8a*q9a*s11**2*s16**2*s22*s23*s26*s34**2*s36 - 1024*q7a*q9a*s11**3*s22**2*s23*s26*s34**2*s36 - 1024*q8a*q9a*s11**3*s14*s22*s24*s26*s34**2*s36 +  \
    2048*q8a*q9a*s11**3*s12*s24**2*s26*s34**2*s36 + 1024*q8a*q9a*s11**2*s12*s15*s22*s25*s26*s34**2*s36 + 1024*q8a*q9a*s11**2*s14*s16*s22*s25*s26*s34**2*s36 + 512*q8a*q9a*s11*s15*s16**2*s22*s25*s26*s34**2*s36 +  \
    1536*q7a*q9a*s11**2*s15*s22**2*s25*s26*s34**2*s36 - 2048*q8a*q9a*s11**2*s12*s16*s24*s25*s26*s34**2*s36 - 2048*q7a*q9a*s11**2*s16*s22*s24*s25*s26*s34**2*s36 + 512*q8a*q9a*s11*s12*s16**2*s25**2*s26*s34**2*s36 -  \
    2048*q7a*q9a*s11**2*s12*s22*s25**2*s26*s34**2*s36 + 512*q7a*q9a*s11*s16**2*s22*s25**2*s26*s34**2*s36 + 512*q8a*q9a*s11**2*s14*s15*s22*s26**2*s34**2*s36 + 1024*q8a*q9a*s11**2*s13*s16*s22*s26**2*s34**2*s36 +  \
    2048*q8a*q9a*s11**2*s12*s16*s23*s26**2*s34**2*s36 + 1024*q7a*q9a*s11**2*s16*s22*s23*s26**2*s34**2*s36 - 512*q8a*q9a*s11**2*s12*s15*s24*s26**2*s34**2*s36 - 512*q8a*q9a*s11**2*s12*s14*s25*s26**2*s34**2*s36 -  \
    512*q8a*q9a*s11*s12*s15*s16*s25*s26**2*s34**2*s36 - 1024*q7a*q9a*s11**2*s14*s22*s25*s26**2*s34**2*s36 - 512*q7a*q9a*s11*s15*s16*s22*s25*s26**2*s34**2*s36 + 3072*q7a*q9a*s11**2*s12*s24*s25*s26**2*s34**2*s36 -  \
    512*q7a*q9a*s11*s12*s16*s25**2*s26**2*s34**2*s36 - 1024*q8a*q9a*s11**2*s12*s13*s26**3*s34**2*s36 - 1024*q7a*q9a*s11**2*s12*s23*s26**3*s34**2*s36 + 512*q7a*q9a*s11*s12*s15*s25*s26**3*s34**2*s36 +  \
    1024*q8a**2*s11**3*s16*s22**2*s33*s34**2*s36 - 1024*q8a**2*s11**2*s16**2*s22*s26*s33*s34**2*s36 - 1024*q7a*q8a*s11**3*s22**2*s26*s33*s34**2*s36 + 1024*q8a**2*s11**2*s12*s16*s26**2*s33*s34**2*s36
v2_2= \
    1024*q7a*q8a*s11**2*s16*s22*s26**2*s33*s34**2*s36 - 1024*q7a*q8a*s11**2*s12*s26**3*s33*s34**2*s36 - 1024*q8a**2*s11**3*s15*s22**2*s34**3*s36 + 512*q8a**2*s11**3*s16*s22*s24*s34**3*s36 + 1024*q8a**2*s11**3*s12*s22*s25*s34**3*s36 -  \
    512*q8a**2*s11**2*s16**2*s22*s25*s34**3*s36 + 512*q8a**2*s11**3*s14*s22*s26*s34**3*s36 + 512*q8a**2*s11**2*s15*s16*s22*s26*s34**3*s36 - 1024*q8a**2*s11**3*s12*s24*s26*s34**3*s36 + 512*q8a**2*s11**2*s12*s16*s25*s26*s34**3*s36 -  \
    512*q8a**2*s11**2*s12*s15*s26**2*s34**3*s36 - 2048*q9a**2*s11**3*s14*s22**3*s23*s35*s36 - 1024*q9a**2*s11**3*s13*s22**3*s24*s35*s36 - 512*q9a**2*s11**2*s15**2*s22**3*s24*s35*s36 + 3072*q9a**2*s11**3*s12*s22**2*s23*s24*s35*s36 -  \
    512*q9a**2*s11**2*s16**2*s22**2*s23*s24*s35*s36 + 1024*q9a**2*s11**3*s14*s22**2*s24**2*s35*s36 + 512*q9a**2*s11**2*s15*s16*s22**2*s24**2*s35*s36 - 1024*q9a**2*s11**3*s12*s22*s24**3*s35*s36 -  \
    512*q9a**2*s11**2*s14*s15*s22**3*s25*s35*s36 + 2560*q9a**2*s11**2*s13*s16*s22**3*s25*s35*s36 + 512*q9a**2*s11*s15**2*s16*s22**3*s25*s35*s36 - 2560*q9a**2*s11**2*s12*s16*s22**2*s23*s25*s35*s36 +  \
    512*q9a**2*s11*s16**3*s22**2*s23*s25*s35*s36 + 1536*q9a**2*s11**2*s12*s15*s22**2*s24*s25*s35*s36 - 2048*q9a**2*s11**2*s14*s16*s22**2*s24*s25*s35*s36 - 512*q9a**2*s11*s15*s16**2*s22**2*s24*s25*s35*s36 +  \
    1536*q9a**2*s11**2*s12*s16*s22*s24**2*s25*s35*s36 + 512*q9a**2*s11**2*s12*s14*s22**2*s25**2*s35*s36 - 1024*q9a**2*s11*s12*s15*s16*s22**2*s25**2*s35*s36 + 1024*q9a**2*s11*s14*s16**2*s22**2*s25**2*s35*s36 -  \
    1024*q9a**2*s11**2*s12**2*s22*s24*s25**2*s35*s36 - 512*q9a**2*s11*s12*s16**2*s22*s24*s25**2*s35*s36 + 512*q9a**2*s11*s12**2*s16*s22*s25**3*s35*s36 - 512*q9a**2*s11**2*s13*s15*s22**3*s26*s35*s36 +  \
    512*q9a**2*s11*s15**3*s22**3*s26*s35*s36 + 512*q9a**2*s11**2*s12*s15*s22**2*s23*s26*s35*s36 + 1536*q9a**2*s11**2*s14*s16*s22**2*s23*s26*s35*s36 + 512*q9a**2*s11*s15*s16**2*s22**2*s23*s26*s35*s36 +  \
    1024*q9a**2*s11**2*s14*s15*s22**2*s24*s26*s35*s36 + 512*q9a**2*s11**2*s13*s16*s22**2*s24*s26*s35*s36 - 512*q9a**2*s11*s15**2*s16*s22**2*s24*s26*s35*s36 - 1024*q9a**2*s11**2*s12*s16*s22*s23*s24*s26*s35*s36 -  \
    1536*q9a**2*s11**2*s12*s15*s22*s24**2*s26*s35*s36 - 512*q9a**2*s11**2*s14*s16*s22*s24**2*s26*s35*s36 + 512*q9a**2*s11**2*s12*s16*s24**3*s26*s35*s36 - 2048*q9a**2*s11**2*s12*s13*s22**2*s25*s26*s35*s36 -  \
    512*q9a**2*s11**2*s14**2*s22**2*s25*s26*s35*s36 - 2048*q9a**2*s11*s12*s15**2*s22**2*s25*s26*s35*s36 + 1024*q9a**2*s11*s14*s15*s16*s22**2*s25*s26*s35*s36 - 2048*q9a**2*s11*s13*s16**2*s22**2*s25*s26*s35*s36 +  \
    2048*q9a**2*s11**2*s12**2*s22*s23*s25*s26*s35*s36 + 2048*q9a**2*s11**2*s12*s14*s22*s24*s25*s26*s35*s36 + 1024*q9a**2*s11*s12*s15*s16*s22*s24*s25*s26*s35*s36 + 512*q9a**2*s11*s14*s16**2*s22*s24*s25*s26*s35*s36 -  \
    1024*q9a**2*s11**2*s12**2*s24**2*s25*s26*s35*s36 - 512*q9a**2*s11*s12*s16**2*s24**2*s25*s26*s35*s36 + 2560*q9a**2*s11*s12**2*s15*s22*s25**2*s26*s35*s36 - 3072*q9a**2*s11*s12*s14*s16*s22*s25**2*s26*s35*s36 +  \
    1536*q9a**2*s11*s12**2*s16*s24*s25**2*s26*s35*s36 - 1024*q9a**2*s11*s12**3*s25**3*s26*s35*s36 - 512*q9a**2*s11**2*s13*s14*s22**2*s26**2*s35*s36 - 1024*q9a**2*s11*s14*s15**2*s22**2*s26**2*s35*s36 +  \
    1024*q9a**2*s11*s13*s15*s16*s22**2*s26**2*s35*s36 - 1024*q9a**2*s11**2*s12*s14*s22*s23*s26**2*s35*s36 - 2048*q9a**2*s11*s12*s15*s16*s22*s23*s26**2*s35*s36 + 512*q9a**2*s11**2*s14**2*s22*s24*s26**2*s35*s36 +  \
    1536*q9a**2*s11*s12*s15**2*s22*s24*s26**2*s35*s36 - 512*q9a**2*s11*s14*s15*s16*s22*s24*s26**2*s35*s36 + 1024*q9a**2*s11**2*s12**2*s23*s24*s26**2*s35*s36 - 512*q9a**2*s11**2*s12*s14*s24**2*s26**2*s35*s36 +  \
    512*q9a**2*s11*s12*s15*s16*s24**2*s26**2*s35*s36 + 1024*q9a**2*s11*s12*s14*s15*s22*s25*s26**2*s35*s36 + 3072*q9a**2*s11*s12*s13*s16*s22*s25*s26**2*s35*s36 - 512*q9a**2*s11*s14**2*s16*s22*s25*s26**2*s35*s36 -  \
    512*q9a**2*s11*s12**2*s16*s23*s25*s26**2*s35*s36 - 2560*q9a**2*s11*s12**2*s15*s24*s25*s26**2*s35*s36 + 512*q9a**2*s11*s12*s14*s16*s24*s25*s26**2*s35*s36 + 1024*q9a**2*s11*s12**2*s14*s25**2*s26**2*s35*s36 -  \
    1024*q9a**2*s11*s12*s13*s15*s22*s26**3*s35*s36 + 512*q9a**2*s11*s14**2*s15*s22*s26**3*s35*s36 + 1536*q9a**2*s11*s12**2*s15*s23*s26**3*s35*s36 - 512*q9a**2*s11*s12*s14*s15*s24*s26**3*s35*s36 -  \
    1024*q9a**2*s11*s12**2*s13*s25*s26**3*s35*s36 - 4096*q8a*q9a*s11**3*s14*s22**3*s33*s35*s36 + 6144*q8a*q9a*s11**3*s12*s22**2*s24*s33*s35*s36 - 1024*q8a*q9a*s11**2*s16**2*s22**2*s24*s33*s35*s36 -  \
    2048*q7a*q9a*s11**3*s22**3*s24*s33*s35*s36 - 5120*q8a*q9a*s11**2*s12*s16*s22**2*s25*s33*s35*s36 + 1024*q8a*q9a*s11*s16**3*s22**2*s25*s33*s35*s36 + 5120*q7a*q9a*s11**2*s16*s22**3*s25*s33*s35*s36 +  \
    1024*q8a*q9a*s11**2*s12*s15*s22**2*s26*s33*s35*s36 + 3072*q8a*q9a*s11**2*s14*s16*s22**2*s26*s33*s35*s36 + 1024*q8a*q9a*s11*s15*s16**2*s22**2*s26*s33*s35*s36 - 1024*q7a*q9a*s11**2*s15*s22**3*s26*s33*s35*s36 -  \
    2048*q8a*q9a*s11**2*s12*s16*s22*s24*s26*s33*s35*s36 + 1024*q7a*q9a*s11**2*s16*s22**2*s24*s26*s33*s35*s36 + 4096*q8a*q9a*s11**2*s12**2*s22*s25*s26*s33*s35*s36 - 4096*q7a*q9a*s11**2*s12*s22**2*s25*s26*s33*s35*s36 -  \
    4096*q7a*q9a*s11*s16**2*s22**2*s25*s26*s33*s35*s36 - 2048*q8a*q9a*s11**2*s12*s14*s22*s26**2*s33*s35*s36 - 4096*q8a*q9a*s11*s12*s15*s16*s22*s26**2*s33*s35*s36 - 1024*q7a*q9a*s11**2*s14*s22**2*s26**2*s33*s35*s36 +  \
    2048*q7a*q9a*s11*s15*s16*s22**2*s26**2*s33*s35*s36 + 2048*q8a*q9a*s11**2*s12**2*s24*s26**2*s33*s35*s36 - 1024*q8a*q9a*s11*s12**2*s16*s25*s26**2*s33*s35*s36 + 6144*q7a*q9a*s11*s12*s16*s22*s25*s26**2*s33*s35*s36 +  \
    3072*q8a*q9a*s11*s12**2*s15*s26**3*s33*s35*s36 - 2048*q7a*q9a*s11*s12*s15*s22*s26**3*s33*s35*s36 - 2048*q7a*q9a*s11*s12**2*s25*s26**3*s33*s35*s36 + 4096*q8a*q9a*s11**3*s13*s22**3*s34*s35*s36 -  \
    8192*q8a*q9a*s11**3*s12*s22**2*s23*s34*s35*s36 + 4096*q7a*q9a*s11**3*s22**3*s23*s34*s35*s36 - 2048*q8a*q9a*s11**3*s14*s22**2*s24*s34*s35*s36 + 2048*q8a*q9a*s11**3*s12*s22*s24**2*s34*s35*s36 +  \
    512*q8a*q9a*s11**2*s16**2*s22*s24**2*s34*s35*s36 + 2048*q8a*q9a*s11**2*s12*s15*s22**2*s25*s34*s35*s36 + 2048*q8a*q9a*s11**2*s14*s16*s22**2*s25*s34*s35*s36 - 2048*q7a*q9a*s11**2*s15*s22**3*s25*s34*s35*s36 -  \
    2048*q8a*q9a*s11**2*s12*s16*s22*s24*s25*s34*s35*s36 - 512*q8a*q9a*s11*s16**3*s22*s24*s25*s34*s35*s36 - 2048*q8a*q9a*s11**2*s12**2*s22*s25**2*s34*s35*s36 + 1536*q8a*q9a*s11*s12*s16**2*s22*s25**2*s34*s35*s36 +  \
    2048*q7a*q9a*s11**2*s12*s22**2*s25**2*s34*s35*s36 - 1536*q7a*q9a*s11*s16**2*s22**2*s25**2*s34*s35*s36 - 4096*q8a*q9a*s11**2*s13*s16*s22**2*s26*s34*s35*s36 + 8192*q8a*q9a*s11**2*s12*s16*s22*s23*s26*s34*s35*s36 -  \
    4096*q7a*q9a*s11**2*s16*s22**2*s23*s26*s34*s35*s36 + 1024*q8a*q9a*s11**2*s14*s16*s22*s24*s26*s34*s35*s36 - 512*q8a*q9a*s11*s15*s16**2*s22*s24*s26*s34*s35*s36 - 2048*q8a*q9a*s11**2*s12*s16*s24**2*s26*s34*s35*s36 -  \
    4096*q8a*q9a*s11**2*s12*s14*s22*s25*s26*s34*s35*s36 - 1024*q8a*q9a*s11*s12*s15*s16*s22*s25*s26*s34*s35*s36 - 512*q8a*q9a*s11*s14*s16**2*s22*s25*s26*s34*s35*s36 + 2048*q7a*q9a*s11**2*s14*s22**2*s25*s26*s34*s35*s36 +  \
    1024*q7a*q9a*s11*s15*s16*s22**2*s25*s26*s34*s35*s36 + 4096*q8a*q9a*s11**2*s12**2*s24*s25*s26*s34*s35*s36 + 1024*q8a*q9a*s11*s12*s16**2*s24*s25*s26*s34*s35*s36 - 2048*q7a*q9a*s11**2*s12*s22*s24*s25*s26*s34*s35*s36 +  \
    1536*q7a*q9a*s11*s16**2*s22*s24*s25*s26*s34*s35*s36 - 2048*q8a*q9a*s11*s12**2*s16*s25**2*s26*s34*s35*s36 + 2048*q7a*q9a*s11*s12*s16*s22*s25**2*s26*s34*s35*s36 + 4096*q8a*q9a*s11**2*s12*s13*s22*s26**2*s34*s35*s36 -  \
    1536*q8a*q9a*s11**2*s14**2*s22*s26**2*s34*s35*s36 - 512*q8a*q9a*s11*s12*s15**2*s22*s26**2*s34*s35*s36 + 1536*q8a*q9a*s11*s14*s15*s16*s22*s26**2*s34*s35*s36 + 512*q7a*q9a*s11*s15**2*s22**2*s26**2*s34*s35*s36 -  \
    8192*q8a*q9a*s11**2*s12**2*s23*s26**2*s34*s35*s36 + 4096*q7a*q9a*s11**2*s12*s22*s23*s26**2*s34*s35*s36 + 2048*q8a*q9a*s11**2*s12*s14*s24*s26**2*s34*s35*s36 - 512*q7a*q9a*s11*s15*s16*s22*s24*s26**2*s34*s35*s36 +  \
    2048*q8a*q9a*s11*s12**2*s15*s25*s26**2*s34*s35*s36 - 2048*q7a*q9a*s11*s12*s15*s22*s25*s26**2*s34*s35*s36 - 512*q7a*q9a*s11*s14*s16*s22*s25*s26**2*s34*s35*s36 - 2048*q7a*q9a*s11*s12*s16*s24*s25*s26**2*s34*s35*s36 -  \
    1024*q8a*q9a*s11*s12*s14*s15*s26**3*s34*s35*s36 - 512*q7a*q9a*s11*s14*s15*s22*s26**3*s34*s35*s36 + 1024*q7a*q9a*s11*s12*s15*s24*s26**3*s34*s35*s36 + 1024*q7a*q9a*s11*s12*s14*s25*s26**3*s34*s35*s36 -  \
    4096*q8a**2*s11**3*s12*s22**2*s33*s34*s35*s36 + 4096*q7a*q8a*s11**3*s22**3*s33*s34*s35*s36 + 4096*q8a**2*s11**2*s12*s16*s22*s26*s33*s34*s35*s36 - 4096*q7a*q8a*s11**2*s16*s22**2*s26*s33*s34*s35*s36 -  \
    4096*q8a**2*s11**2*s12**2*s26**2*s33*s34*s35*s36 + 4096*q7a*q8a*s11**2*s12*s22*s26**2*s33*s34*s35*s36 + 2048*q8a**2*s11**3*s14*s22**2*s34**2*s35*s36 - 1024*q8a**2*s11**3*s12*s22*s24*s34**2*s35*s36 -  \
    512*q8a**2*s11**2*s16**2*s22*s24*s34**2*s35*s36 - 1024*q7a*q8a*s11**3*s22**2*s24*s34**2*s35*s36 - 512*q8a**2*s11**2*s12*s16*s22*s25*s34**2*s35*s36 + 512*q8a**2*s11*s16**3*s22*s25*s34**2*s35*s36 +  \
    512*q7a*q8a*s11**2*s16*s22**2*s25*s34**2*s35*s36 + 512*q8a**2*s11**2*s12*s15*s22*s26*s34**2*s35*s36 - 2560*q8a**2*s11**2*s14*s16*s22*s26*s34**2*s35*s36 + 512*q8a**2*s11*s15*s16**2*s22*s26*s34**2*s35*s36 -  \
    512*q7a*q8a*s11**2*s15*s22**2*s26*s34**2*s35*s36 + 1536*q8a**2*s11**2*s12*s16*s24*s26*s34**2*s35*s36 + 2048*q7a*q8a*s11**2*s16*s22*s24*s26*s34**2*s35*s36 - 1024*q8a**2*s11**2*s12**2*s25*s26*s34**2*s35*s36 -  \
    512*q8a**2*s11*s12*s16**2*s25*s26*s34**2*s35*s36 + 2048*q7a*q8a*s11**2*s12*s22*s25*s26*s34**2*s35*s36 - 1536*q7a*q8a*s11*s16**2*s22*s25*s26*s34**2*s35*s36 - 1024*q7a**2*s11**2*s22**2*s25*s26*s34**2*s35*s36 +  \
    1536*q8a**2*s11**2*s12*s14*s26**2*s34**2*s35*s36 - 512*q8a**2*s11*s12*s15*s16*s26**2*s34**2*s35*s36 + 1024*q7a*q8a*s11**2*s14*s22*s26**2*s34**2*s35*s36 - 512*q7a*q8a*s11*s15*s16*s22*s26**2*s34**2*s35*s36 -  \
    3072*q7a*q8a*s11**2*s12*s24*s26**2*s34**2*s35*s36 + 1536*q7a*q8a*s11*s12*s16*s25*s26**2*s34**2*s35*s36 + 1024*q7a**2*s11*s16*s22*s25*s26**2*s34**2*s35*s36 + 512*q7a*q8a*s11*s12*s15*s26**3*s34**2*s35*s36 -  \
    1024*q7a**2*s11*s12*s25*s26**3*s34**2*s35*s36 + 1024*q8a*q9a*s11**2*s14*s15*s22**3*s35**2*s36 - 1024*q8a*q9a*s11**2*s13*s16*s22**3*s35**2*s36 + 2048*q8a*q9a*s11**2*s12*s16*s22**2*s23*s35**2*s36 -  \
    1024*q7a*q9a*s11**2*s16*s22**3*s23*s35**2*s36 - 2560*q8a*q9a*s11**2*s12*s15*s22**2*s24*s35**2*s36 + 1536*q8a*q9a*s11**2*s14*s16*s22**2*s24*s35**2*s36 + 1536*q7a*q9a*s11**2*s15*s22**3*s24*s35**2*s36 -  \
    1024*q8a*q9a*s11**2*s12*s16*s22*s24**2*s35**2*s36 - 512*q7a*q9a*s11**2*s16*s22**2*s24**2*s35**2*s36 - 512*q8a*q9a*s11**2*s12*s14*s22**2*s25*s35**2*s36 + 1024*q8a*q9a*s11*s12*s15*s16*s22**2*s25*s35**2*s36 -  \
    1024*q8a*q9a*s11*s14*s16**2*s22**2*s25*s35**2*s36 - 512*q7a*q9a*s11**2*s14*s22**3*s25*s35**2*s36 - 1024*q7a*q9a*s11*s15*s16*s22**3*s25*s35**2*s36 + 2048*q8a*q9a*s11**2*s12**2*s22*s24*s25*s35**2*s36 +  \
    512*q8a*q9a*s11*s12*s16**2*s22*s24*s25*s35**2*s36 - 1024*q7a*q9a*s11**2*s12*s22**2*s24*s25*s35**2*s36 + 512*q7a*q9a*s11*s16**2*s22**2*s24*s25*s35**2*s36 - 1024*q8a*q9a*s11*s12**2*s16*s22*s25**2*s35**2*s36 +  \
    1024*q7a*q9a*s11*s12*s16*s22**2*s25**2*s35**2*s36 - 1024*q8a*q9a*s11**2*s12*s13*s22**2*s26*s35**2*s36 + 512*q8a*q9a*s11**2*s14**2*s22**2*s26*s35**2*s36 + 1024*q8a*q9a*s11*s12*s15**2*s22**2*s26*s35**2*s36 -  \
    2048*q8a*q9a*s11*s14*s15*s16*s22**2*s26*s35**2*s36 + 1024*q8a*q9a*s11*s13*s16**2*s22**2*s26*s35**2*s36 + 2048*q7a*q9a*s11**2*s13*s22**3*s26*s35**2*s36 - 1024*q7a*q9a*s11*s15**2*s22**3*s26*s35**2*s36 -  \
    2048*q8a*q9a*s11*s12*s16**2*s22*s23*s26*s35**2*s36 - 1024*q7a*q9a*s11**2*s12*s22**2*s23*s26*s35**2*s36 + 1024*q7a*q9a*s11*s16**2*s22**2*s23*s26*s35**2*s36 - 1024*q8a*q9a*s11**2*s12*s14*s22*s24*s26*s35**2*s36 +  \
    2048*q8a*q9a*s11*s12*s15*s16*s22*s24*s26*s35**2*s36 - 512*q8a*q9a*s11*s14*s16**2*s22*s24*s26*s35**2*s36 - 1536*q7a*q9a*s11**2*s14*s22**2*s24*s26*s35**2*s36 + 512*q8a*q9a*s11*s12*s16**2*s24**2*s26*s35**2*s36 +  \
    2048*q7a*q9a*s11**2*s12*s22*s24**2*s26*s35**2*s36 - 3072*q8a*q9a*s11*s12**2*s15*s22*s25*s26*s35**2*s36 + 3072*q8a*q9a*s11*s12*s14*s16*s22*s25*s26*s35**2*s36 + 3072*q7a*q9a*s11*s12*s15*s22**2*s25*s26*s35**2*s36
v2_3= \
    1024*q7a*q9a*s11*s14*s16*s22**2*s25*s26*s35**2*s36 - 2048*q8a*q9a*s11*s12**2*s16*s24*s25*s26*s35**2*s36 - 2048*q7a*q9a*s11*s12*s16*s22*s24*s25*s26*s35**2*s36 + 2048*q8a*q9a*s11*s12**3*s25**2*s26*s35**2*s36 -  \
    2048*q7a*q9a*s11*s12**2*s22*s25**2*s26*s35**2*s36 + 512*q8a*q9a*s11*s12*s14*s15*s22*s26**2*s35**2*s36 + 512*q8a*q9a*s11*s14**2*s16*s22*s26**2*s35**2*s36 + 1536*q7a*q9a*s11*s14*s15*s22**2*s26**2*s35**2*s36 -  \
    2048*q7a*q9a*s11*s13*s16*s22**2*s26**2*s35**2*s36 + 2048*q8a*q9a*s11*s12**2*s16*s23*s26**2*s35**2*s36 - 512*q8a*q9a*s11*s12**2*s15*s24*s26**2*s35**2*s36 - 512*q8a*q9a*s11*s12*s14*s16*s24*s26**2*s35**2*s36 -  \
    1536*q7a*q9a*s11*s12*s15*s22*s24*s26**2*s35**2*s36 + 512*q7a*q9a*s11*s14*s16*s22*s24*s26**2*s35**2*s36 - 512*q7a*q9a*s11*s12*s16*s24**2*s26**2*s35**2*s36 - 512*q8a*q9a*s11*s12**2*s14*s25*s26**2*s35**2*s36 -  \
    2560*q7a*q9a*s11*s12*s14*s22*s25*s26**2*s35**2*s36 + 3072*q7a*q9a*s11*s12**2*s24*s25*s26**2*s35**2*s36 - 1024*q8a*q9a*s11*s12**2*s13*s26**3*s35**2*s36 + 2048*q7a*q9a*s11*s12*s13*s22*s26**3*s35**2*s36 -  \
    512*q7a*q9a*s11*s14**2*s22*s26**3*s35**2*s36 - 1024*q7a*q9a*s11*s12**2*s23*s26**3*s35**2*s36 + 512*q7a*q9a*s11*s12*s14*s24*s26**3*s35**2*s36 + 1024*q8a**2*s11**2*s12*s16*s22**2*s33*s35**2*s36 -  \
    1024*q7a*q8a*s11**2*s16*s22**3*s33*s35**2*s36 - 1024*q8a**2*s11*s12*s16**2*s22*s26*s33*s35**2*s36 - 1024*q7a*q8a*s11**2*s12*s22**2*s26*s33*s35**2*s36 + 1024*q7a*q8a*s11*s16**2*s22**2*s26*s33*s35**2*s36 +  \
    1024*q7a**2*s11**2*s22**3*s26*s33*s35**2*s36 + 1024*q8a**2*s11*s12**2*s16*s26**2*s33*s35**2*s36 - 1024*q7a**2*s11*s16*s22**2*s26**2*s33*s35**2*s36 - 1024*q7a*q8a*s11*s12**2*s26**3*s33*s35**2*s36 +  \
    1024*q7a**2*s11*s12*s22*s26**3*s33*s35**2*s36 + 1024*q8a**2*s11**2*s12*s15*s22**2*s34*s35**2*s36 - 1024*q8a**2*s11**2*s14*s16*s22**2*s34*s35**2*s36 - 1024*q7a*q8a*s11**2*s15*s22**3*s34*s35**2*s36 +  \
    1536*q8a**2*s11**2*s12*s16*s22*s24*s34*s35**2*s36 - 512*q7a*q8a*s11**2*s16*s22**2*s24*s34*s35**2*s36 + 1024*q8a**2*s11**2*s12**2*s22*s25*s34*s35**2*s36 - 1024*q8a**2*s11*s12*s16**2*s22*s25*s34*s35**2*s36 -  \
    3072*q7a*q8a*s11**2*s12*s22**2*s25*s34*s35**2*s36 + 1024*q7a*q8a*s11*s16**2*s22**2*s25*s34*s35**2*s36 + 2048*q7a**2*s11**2*s22**3*s25*s34*s35**2*s36 + 512*q8a**2*s11**2*s12*s14*s22*s26*s34*s35**2*s36 -  \
    2048*q8a**2*s11*s12*s15*s16*s22*s26*s34*s35**2*s36 + 1024*q8a**2*s11*s14*s16**2*s22*s26*s34*s35**2*s36 + 512*q7a*q8a*s11**2*s14*s22**2*s26*s34*s35**2*s36 + 2048*q7a*q8a*s11*s15*s16*s22**2*s26*s34*s35**2*s36 -  \
    1024*q8a**2*s11**2*s12**2*s24*s26*s34*s35**2*s36 - 512*q8a**2*s11*s12*s16**2*s24*s26*s34*s35**2*s36 - 512*q7a*q8a*s11*s16**2*s22*s24*s26*s34*s35**2*s36 + 1536*q8a**2*s11*s12**2*s16*s25*s26*s34*s35**2*s36 +  \
    1024*q7a*q8a*s11*s12*s16*s22*s25*s26*s34*s35**2*s36 - 2560*q7a**2*s11*s16*s22**2*s25*s26*s34*s35**2*s36 + 1536*q8a**2*s11*s12**2*s15*s26**2*s34*s35**2*s36 - 512*q8a**2*s11*s12*s14*s16*s26**2*s34*s35**2*s36 -  \
    1024*q7a*q8a*s11*s12*s15*s22*s26**2*s34*s35**2*s36 - 1536*q7a*q8a*s11*s14*s16*s22*s26**2*s34*s35**2*s36 - 512*q7a**2*s11*s15*s22**2*s26**2*s34*s35**2*s36 + 1536*q7a*q8a*s11*s12*s16*s24*s26**2*s34*s35**2*s36 +  \
    512*q7a**2*s11*s16*s22*s24*s26**2*s34*s35**2*s36 - 3072*q7a*q8a*s11*s12**2*s25*s26**2*s34*s35**2*s36 + 3072*q7a**2*s11*s12*s22*s25*s26**2*s34*s35**2*s36 + 512*q7a*q8a*s11*s12*s14*s26**3*s34*s35**2*s36 +  \
    512*q7a**2*s11*s14*s22*s26**3*s34*s35**2*s36 - 1024*q7a**2*s11*s12*s24*s26**3*s34*s35**2*s36 - 1024*q8a**2*s11**2*s12**2*s22*s24*s35**3*s36 + 2048*q7a*q8a*s11**2*s12*s22**2*s24*s35**3*s36 - 1024*q7a**2*s11**2*s22**3*s24*s35**3*s36 +  \
    512*q8a**2*s11*s12**2*s16*s22*s25*s35**3*s36 - 1024*q7a*q8a*s11*s12*s16*s22**2*s25*s35**3*s36 + 512*q7a**2*s11*s16*s22**3*s25*s35**3*s36 + 512*q8a**2*s11*s12**2*s15*s22*s26*s35**3*s36 -  \
    1024*q7a*q8a*s11*s12*s15*s22**2*s26*s35**3*s36 + 512*q7a**2*s11*s15*s22**3*s26*s35**3*s36 + 512*q8a**2*s11*s12**2*s16*s24*s26*s35**3*s36 - 1024*q7a*q8a*s11*s12*s16*s22*s24*s26*s35**3*s36 +  \
    512*q7a**2*s11*s16*s22**2*s24*s26*s35**3*s36 - 1024*q8a**2*s11*s12**3*s25*s26*s35**3*s36 + 2048*q7a*q8a*s11*s12**2*s22*s25*s26*s35**3*s36 - 1024*q7a**2*s11*s12*s22**2*s25*s26*s35**3*s36 -  \
    512*q8a**2*s11*s12**2*s14*s26**2*s35**3*s36 + 1024*q7a*q8a*s11*s12*s14*s22*s26**2*s35**3*s36 - 512*q7a**2*s11*s14*s22**2*s26**2*s35**3*s36 - 1024*q9a**2*s11**3*s13*s22**3*s23*s36**2 + 512*q9a**2*s11**2*s15**2*s22**3*s23*s36**2 +  \
    1024*q9a**2*s11**3*s12*s22**2*s23**2*s36**2 + 512*q9a**2*s11**2*s16**2*s22**2*s23**2*s36**2 + 1536*q9a**2*s11**3*s14*s22**2*s23*s24*s36**2 - 512*q9a**2*s11**2*s15*s16*s22**2*s23*s24*s36**2 + 512*q9a**2*s11**3*s13*s22**2*s24**2*s36**2 -  \
    2048*q9a**2*s11**3*s12*s22*s23*s24**2*s36**2 - 512*q9a**2*s11**3*s14*s22*s24**3*s36**2 + 512*q9a**2*s11**3*s12*s24**4*s36**2 + 1536*q9a**2*s11**2*s13*s15*s22**3*s25*s36**2 - 512*q9a**2*s11*s15**3*s22**3*s25*s36**2 -  \
    2560*q9a**2*s11**2*s12*s15*s22**2*s23*s25*s36**2 + 512*q9a**2*s11**2*s14*s16*s22**2*s23*s25*s36**2 - 512*q9a**2*s11*s15*s16**2*s22**2*s23*s25*s36**2 - 512*q9a**2*s11**2*s14*s15*s22**2*s24*s25*s36**2 -  \
    2048*q9a**2*s11**2*s13*s16*s22**2*s24*s25*s36**2 + 512*q9a**2*s11*s15**2*s16*s22**2*s24*s25*s36**2 + 2048*q9a**2*s11**2*s12*s16*s22*s23*s24*s25*s36**2 + 512*q9a**2*s11**2*s12*s15*s22*s24**2*s25*s36**2 +  \
    1024*q9a**2*s11**2*s14*s16*s22*s24**2*s25*s36**2 - 1024*q9a**2*s11**2*s12*s16*s24**3*s25*s36**2 - 1536*q9a**2*s11**2*s12*s13*s22**2*s25**2*s36**2 + 1536*q9a**2*s11**2*s14**2*s22**2*s25**2*s36**2 +  \
    1536*q9a**2*s11*s12*s15**2*s22**2*s25**2*s36**2 - 1536*q9a**2*s11*s14*s15*s16*s22**2*s25**2*s36**2 + 1536*q9a**2*s11*s13*s16**2*s22**2*s25**2*s36**2 + 2048*q9a**2*s11**2*s12**2*s22*s23*s25**2*s36**2 -  \
    1024*q9a**2*s11*s12*s16**2*s22*s23*s25**2*s36**2 - 2560*q9a**2*s11**2*s12*s14*s22*s24*s25**2*s36**2 + 512*q9a**2*s11*s12*s15*s16*s22*s24*s25**2*s36**2 - 512*q9a**2*s11*s14*s16**2*s22*s24*s25**2*s36**2 +  \
    1024*q9a**2*s11**2*s12**2*s24**2*s25**2*s36**2 + 512*q9a**2*s11*s12*s16**2*s24**2*s25**2*s36**2 - 1536*q9a**2*s11*s12**2*s15*s22*s25**3*s36**2 + 1536*q9a**2*s11*s12*s14*s16*s22*s25**3*s36**2 -  \
    1024*q9a**2*s11*s12**2*s16*s24*s25**3*s36**2 + 512*q9a**2*s11*s12**3*s25**4*s36**2 - 2048*q9a**2*s11**2*s14*s15*s22**2*s23*s26*s36**2 + 2048*q9a**2*s11**2*s13*s16*s22**2*s23*s26*s36**2 -  \
    3072*q9a**2*s11**2*s12*s16*s22*s23**2*s26*s36**2 - 512*q9a**2*s11**2*s13*s15*s22**2*s24*s26*s36**2 + 3072*q9a**2*s11**2*s12*s15*s22*s23*s24*s26*s36**2 - 1024*q9a**2*s11**2*s14*s16*s22*s23*s24*s26*s36**2 +  \
    512*q9a**2*s11**2*s14*s15*s22*s24**2*s26*s36**2 + 1024*q9a**2*s11**2*s12*s16*s23*s24**2*s26*s36**2 - 512*q9a**2*s11**2*s12*s15*s24**3*s26*s36**2 + 512*q9a**2*s11**2*s13*s14*s22**2*s25*s26*s36**2 +  \
    1024*q9a**2*s11*s14*s15**2*s22**2*s25*s26*s36**2 - 1024*q9a**2*s11*s13*s15*s16*s22**2*s25*s26*s36**2 + 1024*q9a**2*s11**2*s12*s14*s22*s23*s25*s26*s36**2 + 2048*q9a**2*s11*s12*s15*s16*s22*s23*s25*s26*s36**2 +  \
    2048*q9a**2*s11**2*s12*s13*s22*s24*s25*s26*s36**2 - 1536*q9a**2*s11**2*s14**2*s22*s24*s25*s26*s36**2 - 1536*q9a**2*s11*s12*s15**2*s22*s24*s25*s26*s36**2 + 512*q9a**2*s11*s14*s15*s16*s22*s24*s25*s26*s36**2 -  \
    4096*q9a**2*s11**2*s12**2*s23*s24*s25*s26*s36**2 + 1536*q9a**2*s11**2*s12*s14*s24**2*s25*s26*s36**2 - 512*q9a**2*s11*s12*s15*s16*s24**2*s25*s26*s36**2 - 512*q9a**2*s11*s12*s14*s15*s22*s25**2*s26*s36**2 -  \
    2048*q9a**2*s11*s12*s13*s16*s22*s25**2*s26*s36**2 + 512*q9a**2*s11*s14**2*s16*s22*s25**2*s26*s36**2 + 1024*q9a**2*s11*s12**2*s16*s23*s25**2*s26*s36**2 + 1536*q9a**2*s11*s12**2*s15*s24*s25**2*s26*s36**2 -  \
    512*q9a**2*s11*s12*s14*s16*s24*s25**2*s26*s36**2 - 512*q9a**2*s11*s12**2*s14*s25**3*s26*s36**2 + 512*q9a**2*s11**2*s13**2*s22**2*s26**2*s36**2 - 3072*q9a**2*s11**2*s12*s13*s22*s23*s26**2*s36**2 +  \
    1536*q9a**2*s11**2*s14**2*s22*s23*s26**2*s36**2 + 3072*q9a**2*s11**2*s12**2*s23**2*s26**2*s36**2 - 512*q9a**2*s11**2*s13*s14*s22*s24*s26**2*s36**2 - 1536*q9a**2*s11**2*s12*s14*s23*s24*s26**2*s36**2 +  \
    512*q9a**2*s11**2*s12*s13*s24**2*s26**2*s36**2 + 1024*q9a**2*s11*s12*s13*s15*s22*s25*s26**2*s36**2 - 512*q9a**2*s11*s14**2*s15*s22*s25*s26**2*s36**2 - 1536*q9a**2*s11*s12**2*s15*s23*s25*s26**2*s36**2 +  \
    512*q9a**2*s11*s12*s14*s15*s24*s25*s26**2*s36**2 + 512*q9a**2*s11*s12**2*s13*s25**2*s26**2*s36**2 - 2048*q8a*q9a*s11**3*s13*s22**3*s33*s36**2 + 1024*q8a*q9a*s11**2*s15**2*s22**3*s33*s36**2 +  \
    4096*q8a*q9a*s11**3*s12*s22**2*s23*s33*s36**2 + 2048*q8a*q9a*s11**2*s16**2*s22**2*s23*s33*s36**2 - 2048*q7a*q9a*s11**3*s22**3*s23*s33*s36**2 + 3072*q8a*q9a*s11**3*s14*s22**2*s24*s33*s36**2 -  \
    1024*q8a*q9a*s11**2*s15*s16*s22**2*s24*s33*s36**2 - 4096*q8a*q9a*s11**3*s12*s22*s24**2*s33*s36**2 + 1024*q7a*q9a*s11**3*s22**2*s24**2*s33*s36**2 - 5120*q8a*q9a*s11**2*s12*s15*s22**2*s25*s33*s36**2 +  \
    1024*q8a*q9a*s11**2*s14*s16*s22**2*s25*s33*s36**2 - 1024*q8a*q9a*s11*s15*s16**2*s22**2*s25*s33*s36**2 + 3072*q7a*q9a*s11**2*s15*s22**3*s25*s33*s36**2 + 4096*q8a*q9a*s11**2*s12*s16*s22*s24*s25*s33*s36**2 -  \
    4096*q7a*q9a*s11**2*s16*s22**2*s24*s25*s33*s36**2 + 4096*q8a*q9a*s11**2*s12**2*s22*s25**2*s33*s36**2 - 2048*q8a*q9a*s11*s12*s16**2*s22*s25**2*s33*s36**2 - 3072*q7a*q9a*s11**2*s12*s22**2*s25**2*s33*s36**2 +  \
    3072*q7a*q9a*s11*s16**2*s22**2*s25**2*s33*s36**2 - 4096*q8a*q9a*s11**2*s14*s15*s22**2*s26*s33*s36**2 + 4096*q8a*q9a*s11**2*s13*s16*s22**2*s26*s33*s36**2 - 12288*q8a*q9a*s11**2*s12*s16*s22*s23*s26*s33*s36**2 +  \
    4096*q7a*q9a*s11**2*s16*s22**2*s23*s26*s33*s36**2 + 6144*q8a*q9a*s11**2*s12*s15*s22*s24*s26*s33*s36**2 - 2048*q8a*q9a*s11**2*s14*s16*s22*s24*s26*s33*s36**2 - 1024*q7a*q9a*s11**2*s15*s22**2*s24*s26*s33*s36**2 +  \
    2048*q8a*q9a*s11**2*s12*s16*s24**2*s26*s33*s36**2 + 2048*q8a*q9a*s11**2*s12*s14*s22*s25*s26*s33*s36**2 + 4096*q8a*q9a*s11*s12*s15*s16*s22*s25*s26*s33*s36**2 + 1024*q7a*q9a*s11**2*s14*s22**2*s25*s26*s33*s36**2 -  \
    2048*q7a*q9a*s11*s15*s16*s22**2*s25*s26*s33*s36**2 - 8192*q8a*q9a*s11**2*s12**2*s24*s25*s26*s33*s36**2 + 4096*q7a*q9a*s11**2*s12*s22*s24*s25*s26*s33*s36**2 + 2048*q8a*q9a*s11*s12**2*s16*s25**2*s26*s33*s36**2 -  \
    4096*q7a*q9a*s11*s12*s16*s22*s25**2*s26*s33*s36**2 - 6144*q8a*q9a*s11**2*s12*s13*s22*s26**2*s33*s36**2 + 3072*q8a*q9a*s11**2*s14**2*s22*s26**2*s33*s36**2 + 2048*q7a*q9a*s11**2*s13*s22**2*s26**2*s33*s36**2 +  \
    12288*q8a*q9a*s11**2*s12**2*s23*s26**2*s33*s36**2 - 6144*q7a*q9a*s11**2*s12*s22*s23*s26**2*s33*s36**2 - 3072*q8a*q9a*s11**2*s12*s14*s24*s26**2*s33*s36**2 - 1024*q7a*q9a*s11**2*s14*s22*s24*s26**2*s33*s36**2 +  \
    1024*q7a*q9a*s11**2*s12*s24**2*s26**2*s33*s36**2 - 3072*q8a*q9a*s11*s12**2*s15*s25*s26**2*s33*s36**2 + 2048*q7a*q9a*s11*s12*s15*s22*s25*s26**2*s33*s36**2 + 1024*q7a*q9a*s11*s12**2*s25**2*s26**2*s33*s36**2 +  \
    1024*q8a**2*s11**3*s12*s22**2*s33**2*s36**2 + 512*q8a**2*s11**2*s16**2*s22**2*s33**2*s36**2 - 1024*q7a*q8a*s11**3*s22**3*s33**2*s36**2 - 3072*q8a**2*s11**2*s12*s16*s22*s26*s33**2*s36**2 + 2048*q7a*q8a*s11**2*s16*s22**2*s26*s33**2*s36**2 +  \
    3072*q8a**2*s11**2*s12**2*s26**2*s33**2*s36**2 - 3072*q7a*q8a*s11**2*s12*s22*s26**2*s33**2*s36**2 + 512*q7a**2*s11**2*s22**2*s26**2*s33**2*s36**2 - 4096*q8a*q9a*s11**3*s14*s22**2*s23*s34*s36**2 -  \
    1024*q8a*q9a*s11**3*s13*s22**2*s24*s34*s36**2 - 512*q8a*q9a*s11**2*s15**2*s22**2*s24*s34*s36**2 + 6144*q8a*q9a*s11**3*s12*s22*s23*s24*s34*s36**2 - 1024*q8a*q9a*s11**2*s16**2*s22*s23*s24*s34*s36**2 -  \
    1024*q7a*q9a*s11**3*s22**2*s23*s24*s34*s36**2 + 1024*q8a*q9a*s11**3*s14*s22*s24**2*s34*s36**2 + 512*q8a*q9a*s11**2*s15*s16*s22*s24**2*s34*s36**2 - 1024*q8a*q9a*s11**3*s12*s24**3*s34*s36**2 +  \
    2560*q8a*q9a*s11**2*s14*s15*s22**2*s25*s34*s36**2 + 512*q8a*q9a*s11**2*s13*s16*s22**2*s25*s34*s36**2 - 1024*q8a*q9a*s11**2*s12*s16*s22*s23*s25*s34*s36**2 + 512*q7a*q9a*s11**2*s16*s22**2*s23*s25*s34*s36**2 -  \
    1024*q8a*q9a*s11**2*s12*s15*s22*s24*s25*s34*s36**2 - 3072*q8a*q9a*s11**2*s14*s16*s22*s24*s25*s34*s36**2 + 512*q8a*q9a*s11*s15*s16**2*s22*s24*s25*s34*s36**2 - 512*q7a*q9a*s11**2*s15*s22**2*s24*s25*s34*s36**2 +  \
    1536*q8a*q9a*s11**2*s12*s16*s24**2*s25*s34*s36**2 + 1024*q7a*q9a*s11**2*s16*s22*s24**2*s25*s34*s36**2 - 512*q8a*q9a*s11*s12*s15*s16*s22*s25**2*s34*s36**2 + 1024*q8a*q9a*s11*s14*s16**2*s22*s25**2*s34*s36**2 -  \
    2560*q7a*q9a*s11**2*s14*s22**2*s25**2*s34*s36**2 + 512*q7a*q9a*s11*s15*s16*s22**2*s25**2*s34*s36**2 - 1024*q8a*q9a*s11**2*s12**2*s24*s25**2*s34*s36**2 - 512*q8a*q9a*s11*s12*s16**2*s24*s25**2*s34*s36**2 +  \
    3072*q7a*q9a*s11**2*s12*s22*s24*s25**2*s34*s36**2 - 1024*q7a*q9a*s11*s16**2*s22*s24*s25**2*s34*s36**2 + 512*q8a*q9a*s11*s12**2*s16*s25**3*s34*s36**2 - 512*q7a*q9a*s11*s12*s16*s22*s25**3*s34*s36**2
v2_4 =  \
    1536*q8a*q9a*s11**2*s13*s15*s22**2*s26*s34*s36**2 - 3072*q8a*q9a*s11**2*s12*s15*s22*s23*s26*s34*s36**2 + 5120*q8a*q9a*s11**2*s14*s16*s22*s23*s26*s34*s36**2 + 1536*q7a*q9a*s11**2*s15*s22**2*s23*s26*s34*s36**2 -  \
    1024*q8a*q9a*s11**2*s13*s16*s22*s24*s26*s34*s36**2 - 1024*q8a*q9a*s11**2*s12*s16*s23*s24*s26*s34*s36**2 - 1024*q7a*q9a*s11**2*s16*s22*s23*s24*s26*s34*s36**2 - 512*q8a*q9a*s11**2*s12*s15*s24**2*s26*s34*s36**2 -  \
    512*q8a*q9a*s11**2*s14**2*s22*s25*s26*s34*s36**2 + 512*q8a*q9a*s11*s12*s15**2*s22*s25*s26*s34*s36**2 - 1536*q8a*q9a*s11*s14*s15*s16*s22*s25*s26*s34*s36**2 - 2048*q7a*q9a*s11**2*s13*s22**2*s25*s26*s34*s36**2 -  \
    512*q7a*q9a*s11*s15**2*s22**2*s25*s26*s34*s36**2 + 2048*q8a*q9a*s11**2*s12**2*s23*s25*s26*s34*s36**2 + 2048*q8a*q9a*s11**2*s12*s14*s24*s25*s26*s34*s36**2 + 2048*q7a*q9a*s11**2*s14*s22*s24*s25*s26*s34*s36**2 +  \
    512*q7a*q9a*s11*s15*s16*s22*s24*s25*s26*s34*s36**2 - 3072*q7a*q9a*s11**2*s12*s24**2*s25*s26*s34*s36**2 - 512*q8a*q9a*s11*s12**2*s15*s25**2*s26*s34*s36**2 - 512*q8a*q9a*s11*s12*s14*s16*s25**2*s26*s34*s36**2 +  \
    512*q7a*q9a*s11*s12*s15*s22*s25**2*s26*s34*s36**2 + 1536*q7a*q9a*s11*s12*s16*s24*s25**2*s26*s34*s36**2 - 1024*q8a*q9a*s11**2*s13*s14*s22*s26**2*s34*s36**2 - 3072*q8a*q9a*s11**2*s12*s14*s23*s26**2*s34*s36**2 -  \
    1024*q7a*q9a*s11**2*s14*s22*s23*s26**2*s34*s36**2 + 2048*q8a*q9a*s11**2*s12*s13*s24*s26**2*s34*s36**2 + 2048*q7a*q9a*s11**2*s12*s23*s24*s26**2*s34*s36**2 + 1024*q8a*q9a*s11*s12*s14*s15*s25*s26**2*s34*s36**2 +  \
    512*q7a*q9a*s11*s14*s15*s22*s25*s26**2*s34*s36**2 - 1024*q7a*q9a*s11*s12*s15*s24*s25*s26**2*s34*s36**2 - 512*q7a*q9a*s11*s12*s14*s25**2*s26**2*s34*s36**2 - 2048*q8a**2*s11**3*s14*s22**2*s33*s34*s36**2 +  \
    3072*q8a**2*s11**3*s12*s22*s24*s33*s34*s36**2 - 512*q8a**2*s11**2*s16**2*s22*s24*s33*s34*s36**2 - 1024*q7a*q8a*s11**3*s22**2*s24*s33*s34*s36**2 - 512*q8a**2*s11**2*s12*s16*s22*s25*s33*s34*s36**2 +  \
    512*q7a*q8a*s11**2*s16*s22**2*s25*s33*s34*s36**2 - 1536*q8a**2*s11**2*s12*s15*s22*s26*s33*s34*s36**2 + 2560*q8a**2*s11**2*s14*s16*s22*s26*s33*s34*s36**2 + 1536*q7a*q8a*s11**2*s15*s22**2*s26*s33*s34*s36**2 -  \
    512*q8a**2*s11**2*s12*s16*s24*s26*s33*s34*s36**2 - 1024*q7a*q8a*s11**2*s16*s22*s24*s26*s33*s34*s36**2 + 1024*q8a**2*s11**2*s12**2*s25*s26*s33*s34*s36**2 - 1024*q7a**2*s11**2*s22**2*s25*s26*s33*s34*s36**2 -  \
    1536*q8a**2*s11**2*s12*s14*s26**2*s33*s34*s36**2 - 1024*q7a*q8a*s11**2*s14*s22*s26**2*s33*s34*s36**2 + 2048*q7a*q8a*s11**2*s12*s24*s26**2*s33*s34*s36**2 + 1024*q8a**2*s11**3*s13*s22**2*s34**2*s36**2 +  \
    512*q8a**2*s11**2*s15**2*s22**2*s34**2*s36**2 - 3072*q8a**2*s11**3*s12*s22*s23*s34**2*s36**2 + 1536*q8a**2*s11**2*s16**2*s22*s23*s34**2*s36**2 + 2048*q7a*q8a*s11**3*s22**2*s23*s34**2*s36**2 - 512*q8a**2*s11**3*s14*s22*s24*s34**2*s36**2 -  \
    512*q8a**2*s11**2*s15*s16*s22*s24*s34**2*s36**2 + 512*q8a**2*s11**3*s12*s24**2*s34**2*s36**2 + 1024*q8a**2*s11**2*s14*s16*s22*s25*s34**2*s36**2 - 512*q8a**2*s11*s15*s16**2*s22*s25*s34**2*s36**2 -  \
    1024*q7a*q8a*s11**2*s15*s22**2*s25*s34**2*s36**2 - 512*q8a**2*s11**2*s12*s16*s24*s25*s34**2*s36**2 + 512*q8a**2*s11**2*s12**2*s25**2*s34**2*s36**2 - 1024*q7a*q8a*s11**2*s12*s22*s25**2*s34**2*s36**2 +  \
    512*q7a*q8a*s11*s16**2*s22*s25**2*s34**2*s36**2 + 1024*q7a**2*s11**2*s22**2*s25**2*s34**2*s36**2 - 512*q8a**2*s11**2*s14*s15*s22*s26*s34**2*s36**2 - 512*q8a**2*s11**2*s13*s16*s22*s26*s34**2*s36**2 -  \
    1536*q8a**2*s11**2*s12*s16*s23*s26*s34**2*s36**2 - 1024*q7a*q8a*s11**2*s16*s22*s23*s26*s34**2*s36**2 + 1024*q8a**2*s11**2*s12*s15*s24*s26*s34**2*s36**2 - 512*q8a**2*s11**2*s12*s14*s25*s26*s34**2*s36**2 +  \
    512*q8a**2*s11*s12*s15*s16*s25*s26*s34**2*s36**2 + 512*q7a*q8a*s11*s15*s16*s22*s25*s26*s34**2*s36**2 - 512*q7a*q8a*s11*s12*s16*s25**2*s26*s34**2*s36**2 - 512*q7a**2*s11*s16*s22*s25**2*s26*s34**2*s36**2 +  \
    512*q8a**2*s11**2*s12*s13*s26**2*s34**2*s36**2 + 1024*q7a*q8a*s11**2*s12*s23*s26**2*s34**2*s36**2 - 512*q7a*q8a*s11*s12*s15*s25*s26**2*s34**2*s36**2 + 512*q7a**2*s11*s12*s25**2*s26**2*s34**2*s36**2 -  \
    1024*q8a*q9a*s11**2*s13*s15*s22**3*s35*s36**2 + 2048*q8a*q9a*s11**2*s12*s15*s22**2*s23*s35*s36**2 - 2048*q8a*q9a*s11**2*s14*s16*s22**2*s23*s35*s36**2 - 1024*q7a*q9a*s11**2*s15*s22**3*s23*s35*s36**2 -  \
    512*q8a*q9a*s11**2*s14*s15*s22**2*s24*s35*s36**2 + 1536*q8a*q9a*s11**2*s13*s16*s22**2*s24*s35*s36**2 - 1024*q8a*q9a*s11**2*s12*s16*s22*s23*s24*s35*s36**2 + 1536*q7a*q9a*s11**2*s16*s22**2*s23*s24*s35*s36**2 +  \
    1024*q8a*q9a*s11**2*s12*s15*s22*s24**2*s35*s36**2 - 512*q8a*q9a*s11**2*s14*s16*s22*s24**2*s35*s36**2 - 512*q7a*q9a*s11**2*s15*s22**2*s24**2*s35*s36**2 + 512*q8a*q9a*s11**2*s12*s16*s24**3*s35*s36**2 +  \
    5120*q8a*q9a*s11**2*s12*s13*s22**2*s25*s35*s36**2 - 2560*q8a*q9a*s11**2*s14**2*s22**2*s25*s35*s36**2 - 1024*q8a*q9a*s11*s12*s15**2*s22**2*s25*s35*s36**2 + 2048*q8a*q9a*s11*s14*s15*s16*s22**2*s25*s35*s36**2 -  \
    1024*q8a*q9a*s11*s13*s16**2*s22**2*s25*s35*s36**2 - 4096*q7a*q9a*s11**2*s13*s22**3*s25*s35*s36**2 + 1024*q7a*q9a*s11*s15**2*s22**3*s25*s35*s36**2 - 6144*q8a*q9a*s11**2*s12**2*s22*s23*s25*s35*s36**2 +  \
    2048*q8a*q9a*s11*s12*s16**2*s22*s23*s25*s35*s36**2 + 5120*q7a*q9a*s11**2*s12*s22**2*s23*s25*s35*s36**2 - 1024*q7a*q9a*s11*s16**2*s22**2*s23*s25*s35*s36**2 + 3072*q8a*q9a*s11**2*s12*s14*s22*s24*s25*s35*s36**2 -  \
    2048*q8a*q9a*s11*s12*s15*s16*s22*s24*s25*s35*s36**2 + 512*q8a*q9a*s11*s14*s16**2*s22*s24*s25*s35*s36**2 + 2560*q7a*q9a*s11**2*s14*s22**2*s24*s25*s35*s36**2 - 1024*q8a*q9a*s11**2*s12**2*s24**2*s25*s35*s36**2 -  \
    512*q8a*q9a*s11*s12*s16**2*s24**2*s25*s35*s36**2 - 2048*q7a*q9a*s11**2*s12*s22*s24**2*s25*s35*s36**2 + 2048*q8a*q9a*s11*s12**2*s15*s22*s25**2*s35*s36**2 - 1536*q8a*q9a*s11*s12*s14*s16*s22*s25**2*s35*s36**2 -  \
    2048*q7a*q9a*s11*s12*s15*s22**2*s25**2*s35*s36**2 - 512*q7a*q9a*s11*s14*s16*s22**2*s25**2*s35*s36**2 + 1536*q8a*q9a*s11*s12**2*s16*s24*s25**2*s35*s36**2 + 512*q7a*q9a*s11*s12*s16*s22*s24*s25**2*s35*s36**2 -  \
    1024*q8a*q9a*s11*s12**3*s25**3*s35*s36**2 + 1024*q7a*q9a*s11*s12**2*s22*s25**3*s35*s36**2 + 512*q8a*q9a*s11**2*s13*s14*s22**2*s26*s35*s36**2 + 1024*q8a*q9a*s11*s14*s15**2*s22**2*s26*s35*s36**2 -  \
    1024*q8a*q9a*s11*s13*s15*s16*s22**2*s26*s35*s36**2 + 1024*q8a*q9a*s11**2*s12*s14*s22*s23*s26*s35*s36**2 + 2048*q8a*q9a*s11*s12*s15*s16*s22*s23*s26*s35*s36**2 + 512*q7a*q9a*s11**2*s14*s22**2*s23*s26*s35*s36**2 -  \
    1024*q7a*q9a*s11*s15*s16*s22**2*s23*s26*s35*s36**2 - 2048*q8a*q9a*s11**2*s12*s13*s22*s24*s26*s35*s36**2 + 512*q8a*q9a*s11**2*s14**2*s22*s24*s26*s35*s36**2 - 1536*q8a*q9a*s11*s12*s15**2*s22*s24*s26*s35*s36**2 +  \
    512*q8a*q9a*s11*s14*s15*s16*s22*s24*s26*s35*s36**2 + 512*q7a*q9a*s11*s15**2*s22**2*s24*s26*s35*s36**2 + 2048*q8a*q9a*s11**2*s12**2*s23*s24*s26*s35*s36**2 - 2048*q7a*q9a*s11**2*s12*s22*s23*s24*s26*s35*s36**2 -  \
    512*q8a*q9a*s11**2*s12*s14*s24**2*s26*s35*s36**2 - 512*q8a*q9a*s11*s12*s15*s16*s24**2*s26*s35*s36**2 - 1024*q8a*q9a*s11*s12*s14*s15*s22*s25*s26*s35*s36**2 - 2048*q8a*q9a*s11*s12*s13*s16*s22*s25*s26*s35*s36**2 -  \
    3072*q7a*q9a*s11*s14*s15*s22**2*s25*s26*s35*s36**2 + 5120*q7a*q9a*s11*s13*s16*s22**2*s25*s26*s35*s36**2 - 1024*q8a*q9a*s11*s12**2*s16*s23*s25*s26*s35*s36**2 - 2048*q7a*q9a*s11*s12*s16*s22*s23*s25*s26*s35*s36**2 +  \
    2048*q8a*q9a*s11*s12**2*s15*s24*s25*s26*s35*s36**2 + 2048*q7a*q9a*s11*s12*s15*s22*s24*s25*s26*s35*s36**2 - 1536*q7a*q9a*s11*s14*s16*s22*s24*s25*s26*s35*s36**2 + 1536*q7a*q9a*s11*s12*s16*s24**2*s25*s26*s35*s36**2 -  \
    512*q8a*q9a*s11*s12**2*s14*s25**2*s26*s35*s36**2 + 3584*q7a*q9a*s11*s12*s14*s22*s25**2*s26*s35*s36**2 - 3072*q7a*q9a*s11*s12**2*s24*s25**2*s26*s35*s36**2 + 2048*q8a*q9a*s11*s12*s13*s15*s22*s26**2*s35*s36**2 -  \
    1024*q8a*q9a*s11*s14**2*s15*s22*s26**2*s35*s36**2 - 1024*q7a*q9a*s11*s13*s15*s22**2*s26**2*s35*s36**2 - 3072*q8a*q9a*s11*s12**2*s15*s23*s26**2*s35*s36**2 + 2048*q7a*q9a*s11*s12*s15*s22*s23*s26**2*s35*s36**2 +  \
    1024*q8a*q9a*s11*s12*s14*s15*s24*s26**2*s35*s36**2 + 512*q7a*q9a*s11*s14*s15*s22*s24*s26**2*s35*s36**2 - 512*q7a*q9a*s11*s12*s15*s24**2*s26**2*s35*s36**2 + 2048*q8a*q9a*s11*s12**2*s13*s25*s26**2*s35*s36**2 -  \
    4096*q7a*q9a*s11*s12*s13*s22*s25*s26**2*s35*s36**2 + 1024*q7a*q9a*s11*s14**2*s22*s25*s26**2*s35*s36**2 + 2048*q7a*q9a*s11*s12**2*s23*s25*s26**2*s35*s36**2 - 1024*q7a*q9a*s11*s12*s14*s24*s25*s26**2*s35*s36**2 +  \
    1024*q8a**2*s11**2*s12*s15*s22**2*s33*s35*s36**2 - 1024*q8a**2*s11**2*s14*s16*s22**2*s33*s35*s36**2 - 1024*q7a*q8a*s11**2*s15*s22**3*s33*s35*s36**2 - 512*q8a**2*s11**2*s12*s16*s22*s24*s33*s35*s36**2 +  \
    1536*q7a*q8a*s11**2*s16*s22**2*s24*s33*s35*s36**2 - 3072*q8a**2*s11**2*s12**2*s22*s25*s33*s35*s36**2 + 1024*q8a**2*s11*s12*s16**2*s22*s25*s33*s35*s36**2 + 5120*q7a*q8a*s11**2*s12*s22**2*s25*s33*s35*s36**2 -  \
    1024*q7a*q8a*s11*s16**2*s22**2*s25*s33*s35*s36**2 - 2048*q7a**2*s11**2*s22**3*s25*s33*s35*s36**2 + 512*q8a**2*s11**2*s12*s14*s22*s26*s33*s35*s36**2 + 1024*q8a**2*s11*s12*s15*s16*s22*s26*s33*s35*s36**2 +  \
    512*q7a*q8a*s11**2*s14*s22**2*s26*s33*s35*s36**2 - 1024*q7a*q8a*s11*s15*s16*s22**2*s26*s33*s35*s36**2 + 1024*q8a**2*s11**2*s12**2*s24*s26*s33*s35*s36**2 - 2048*q7a*q8a*s11**2*s12*s22*s24*s26*s33*s35*s36**2 -  \
    512*q8a**2*s11*s12**2*s16*s25*s26*s33*s35*s36**2 - 2048*q7a*q8a*s11*s12*s16*s22*s25*s26*s33*s35*s36**2 + 2560*q7a**2*s11*s16*s22**2*s25*s26*s33*s35*s36**2 - 1536*q8a**2*s11*s12**2*s15*s26**2*s33*s35*s36**2 +  \
    2048*q7a*q8a*s11*s12*s15*s22*s26**2*s33*s35*s36**2 - 512*q7a**2*s11*s15*s22**2*s26**2*s33*s35*s36**2 + 2048*q7a*q8a*s11*s12**2*s25*s26**2*s33*s35*s36**2 - 2048*q7a**2*s11*s12*s22*s25*s26**2*s33*s35*s36**2 -  \
    1024*q8a**2*s11**2*s14*s15*s22**2*s34*s35*s36**2 + 1024*q8a**2*s11**2*s13*s16*s22**2*s34*s35*s36**2 - 3072*q8a**2*s11**2*s12*s16*s22*s23*s34*s35*s36**2 + 2048*q7a*q8a*s11**2*s16*s22**2*s23*s34*s35*s36**2 -  \
    512*q8a**2*s11**2*s12*s15*s22*s24*s34*s35*s36**2 + 1536*q8a**2*s11**2*s14*s16*s22*s24*s34*s35*s36**2 + 1536*q7a*q8a*s11**2*s15*s22**2*s24*s34*s35*s36**2 - 512*q8a**2*s11**2*s12*s16*s24**2*s34*s35*s36**2 -  \
    1024*q7a*q8a*s11**2*s16*s22*s24**2*s34*s35*s36**2 + 512*q8a**2*s11**2*s12*s14*s22*s25*s34*s35*s36**2 + 2048*q8a**2*s11*s12*s15*s16*s22*s25*s34*s35*s36**2 - 1024*q8a**2*s11*s14*s16**2*s22*s25*s34*s35*s36**2 +  \
    512*q7a*q8a*s11**2*s14*s22**2*s25*s34*s35*s36**2 - 2048*q7a*q8a*s11*s15*s16*s22**2*s25*s34*s35*s36**2 + 512*q8a**2*s11*s12*s16**2*s24*s25*s34*s35*s36**2 + 512*q7a*q8a*s11*s16**2*s22*s24*s25*s34*s35*s36**2 -  \
    1024*q7a**2*s11**2*s22**2*s24*s25*s34*s35*s36**2 - 512*q8a**2*s11*s12**2*s16*s25**2*s34*s35*s36**2 - 1024*q7a*q8a*s11*s12*s16*s22*s25**2*s34*s35*s36**2 + 1536*q7a**2*s11*s16*s22**2*s25**2*s34*s35*s36**2 -  \
    3072*q8a**2*s11**2*s12*s13*s22*s26*s34*s35*s36**2 + 1536*q8a**2*s11**2*s14**2*s22*s26*s34*s35*s36**2 + 1024*q8a**2*s11*s12*s15**2*s22*s26*s34*s35*s36**2 - 1024*q8a**2*s11*s14*s15*s16*s22*s26*s34*s35*s36**2 +  \
    2048*q7a*q8a*s11**2*s13*s22**2*s26*s34*s35*s36**2 - 1024*q7a*q8a*s11*s15**2*s22**2*s26*s34*s35*s36**2 + 6144*q8a**2*s11**2*s12**2*s23*s26*s34*s35*s36**2 - 6144*q7a*q8a*s11**2*s12*s22*s23*s26*s34*s35*s36**2 +  \
    1024*q7a**2*s11**2*s22**2*s23*s26*s34*s35*s36**2 - 2560*q8a**2*s11**2*s12*s14*s24*s26*s34*s35*s36**2 + 512*q8a**2*s11*s12*s15*s16*s24*s26*s34*s35*s36**2 - 2048*q7a*q8a*s11**2*s14*s22*s24*s26*s34*s35*s36**2 +  \
    512*q7a*q8a*s11*s15*s16*s22*s24*s26*s34*s35*s36**2 + 3072*q7a*q8a*s11**2*s12*s24**2*s26*s34*s35*s36**2 - 2560*q8a**2*s11*s12**2*s15*s25*s26*s34*s35*s36**2 + 512*q8a**2*s11*s12*s14*s16*s25*s26*s34*s35*s36**2 +  \
    1024*q7a*q8a*s11*s12*s15*s22*s25*s26*s34*s35*s36**2 + 2560*q7a*q8a*s11*s14*s16*s22*s25*s26*s34*s35*s36**2 + 1536*q7a**2*s11*s15*s22**2*s25*s26*s34*s35*s36**2 - 2048*q7a*q8a*s11*s12*s16*s24*s25*s26*s34*s35*s36**2 -  \
    1024*q7a**2*s11*s16*s22*s24*s25*s26*s34*s35*s36**2 + 3072*q7a*q8a*s11*s12**2*s25**2*s26*s34*s35*s36**2 - 3072*q7a**2*s11*s12*s22*s25**2*s26*s34*s35*s36**2 + 512*q8a**2*s11*s12*s14*s15*s26**2*s34*s35*s36**2 +  \
    512*q7a*q8a*s11*s14*s15*s22*s26**2*s34*s35*s36**2 - 1024*q7a*q8a*s11*s12*s15*s24*s26**2*s34*s35*s36**2 - 1024*q7a*q8a*s11*s12*s14*s25*s26**2*s34*s35*s36**2 - 1024*q7a**2*s11*s14*s22*s25*s26**2*s34*s35*s36**2 +  \
    2048*q7a**2*s11*s12*s24*s25*s26**2*s34*s35*s36**2 - 2048*q8a**2*s11**2*s12*s13*s22**2*s35**2*s36**2 + 1024*q8a**2*s11**2*s14**2*s22**2*s35**2*s36**2 + 2048*q7a*q8a*s11**2*s13*s22**3*s35**2*s36**2 +  \
    3072*q8a**2*s11**2*s12**2*s22*s23*s35**2*s36**2 - 4096*q7a*q8a*s11**2*s12*s22**2*s23*s35**2*s36**2 + 1024*q7a**2*s11**2*s22**3*s23*s35**2*s36**2 - 1024*q8a**2*s11**2*s12*s14*s22*s24*s35**2*s36**2 -  \
    1024*q7a*q8a*s11**2*s14*s22**2*s24*s35**2*s36**2 + 512*q8a**2*s11**2*s12**2*s24**2*s35**2*s36**2 + 512*q7a**2*s11**2*s22**2*s24**2*s35**2*s36**2 - 512*q8a**2*s11*s12**2*s15*s22*s25*s35**2*s36**2
v2_5= \
    1024*q7a*q8a*s11*s12*s15*s22**2*s25*s35**2*s36**2 - 512*q7a**2*s11*s15*s22**3*s25*s35**2*s36**2 - 512*q8a**2*s11*s12**2*s16*s24*s25*s35**2*s36**2 + 1024*q7a*q8a*s11*s12*s16*s22*s24*s25*s35**2*s36**2 -  \
    512*q7a**2*s11*s16*s22**2*s24*s25*s35**2*s36**2 + 512*q8a**2*s11*s12**3*s25**2*s35**2*s36**2 - 1024*q7a*q8a*s11*s12**2*s22*s25**2*s35**2*s36**2 + 512*q7a**2*s11*s12*s22**2*s25**2*s35**2*s36**2 +  \
    1024*q8a**2*s11*s12*s13*s16*s22*s26*s35**2*s36**2 - 512*q8a**2*s11*s14**2*s16*s22*s26*s35**2*s36**2 - 1024*q7a*q8a*s11*s13*s16*s22**2*s26*s35**2*s36**2 - 1536*q8a**2*s11*s12**2*s16*s23*s26*s35**2*s36**2 +  \
    2048*q7a*q8a*s11*s12*s16*s22*s23*s26*s35**2*s36**2 - 512*q7a**2*s11*s16*s22**2*s23*s26*s35**2*s36**2 - 512*q8a**2*s11*s12**2*s15*s24*s26*s35**2*s36**2 + 512*q8a**2*s11*s12*s14*s16*s24*s26*s35**2*s36**2 +  \
    1024*q7a*q8a*s11*s12*s15*s22*s24*s26*s35**2*s36**2 + 512*q7a*q8a*s11*s14*s16*s22*s24*s26*s35**2*s36**2 - 512*q7a**2*s11*s15*s22**2*s24*s26*s35**2*s36**2 - 512*q7a*q8a*s11*s12*s16*s24**2*s26*s35**2*s36**2 +  \
    1024*q8a**2*s11*s12**2*s14*s25*s26*s35**2*s36**2 - 2048*q7a*q8a*s11*s12*s14*s22*s25*s26*s35**2*s36**2 + 1024*q7a**2*s11*s14*s22**2*s25*s26*s35**2*s36**2 + 512*q8a**2*s11*s12**2*s13*s26**2*s35**2*s36**2 -  \
    2048*q7a*q8a*s11*s12*s13*s22*s26**2*s35**2*s36**2 + 512*q7a*q8a*s11*s14**2*s22*s26**2*s35**2*s36**2 + 1536*q7a**2*s11*s13*s22**2*s26**2*s35**2*s36**2 + 1024*q7a*q8a*s11*s12**2*s23*s26**2*s35**2*s36**2 -  \
    1024*q7a**2*s11*s12*s22*s23*s26**2*s35**2*s36**2 - 512*q7a*q8a*s11*s12*s14*s24*s26**2*s35**2*s36**2 - 512*q7a**2*s11*s14*s22*s24*s26**2*s35**2*s36**2 + 512*q7a**2*s11*s12*s24**2*s26**2*s35**2*s36**2 +  \
    2048*q8a*q9a*s11**2*s14*s15*s22**2*s23*s36**3 - 2048*q8a*q9a*s11**2*s13*s16*s22**2*s23*s36**3 + 3072*q8a*q9a*s11**2*s12*s16*s22*s23**2*s36**3 - 1024*q7a*q9a*s11**2*s16*s22**2*s23**2*s36**3 +  \
    512*q8a*q9a*s11**2*s13*s15*s22**2*s24*s36**3 - 3072*q8a*q9a*s11**2*s12*s15*s22*s23*s24*s36**3 + 1024*q8a*q9a*s11**2*s14*s16*s22*s23*s24*s36**3 + 512*q7a*q9a*s11**2*s15*s22**2*s23*s24*s36**3 -  \
    512*q8a*q9a*s11**2*s14*s15*s22*s24**2*s36**3 - 1024*q8a*q9a*s11**2*s12*s16*s23*s24**2*s36**3 + 512*q8a*q9a*s11**2*s12*s15*s24**3*s36**3 - 512*q8a*q9a*s11**2*s13*s14*s22**2*s25*s36**3 -  \
    1024*q8a*q9a*s11*s14*s15**2*s22**2*s25*s36**3 + 1024*q8a*q9a*s11*s13*s15*s16*s22**2*s25*s36**3 - 1024*q8a*q9a*s11**2*s12*s14*s22*s23*s25*s36**3 - 2048*q8a*q9a*s11*s12*s15*s16*s22*s23*s25*s36**3 -  \
    512*q7a*q9a*s11**2*s14*s22**2*s23*s25*s36**3 + 1024*q7a*q9a*s11*s15*s16*s22**2*s23*s25*s36**3 - 2048*q8a*q9a*s11**2*s12*s13*s22*s24*s25*s36**3 + 1536*q8a*q9a*s11**2*s14**2*s22*s24*s25*s36**3 +  \
    1536*q8a*q9a*s11*s12*s15**2*s22*s24*s25*s36**3 - 512*q8a*q9a*s11*s14*s15*s16*s22*s24*s25*s36**3 + 2048*q7a*q9a*s11**2*s13*s22**2*s24*s25*s36**3 - 512*q7a*q9a*s11*s15**2*s22**2*s24*s25*s36**3 +  \
    4096*q8a*q9a*s11**2*s12**2*s23*s24*s25*s36**3 - 2048*q7a*q9a*s11**2*s12*s22*s23*s24*s25*s36**3 - 1536*q8a*q9a*s11**2*s12*s14*s24**2*s25*s36**3 + 512*q8a*q9a*s11*s12*s15*s16*s24**2*s25*s36**3 -  \
    1024*q7a*q9a*s11**2*s14*s22*s24**2*s25*s36**3 + 1024*q7a*q9a*s11**2*s12*s24**3*s25*s36**3 + 512*q8a*q9a*s11*s12*s14*s15*s22*s25**2*s36**3 + 2048*q8a*q9a*s11*s12*s13*s16*s22*s25**2*s36**3 -  \
    512*q8a*q9a*s11*s14**2*s16*s22*s25**2*s36**3 + 1536*q7a*q9a*s11*s14*s15*s22**2*s25**2*s36**3 - 3072*q7a*q9a*s11*s13*s16*s22**2*s25**2*s36**3 - 1024*q8a*q9a*s11*s12**2*s16*s23*s25**2*s36**3 +  \
    2048*q7a*q9a*s11*s12*s16*s22*s23*s25**2*s36**3 - 1536*q8a*q9a*s11*s12**2*s15*s24*s25**2*s36**3 + 512*q8a*q9a*s11*s12*s14*s16*s24*s25**2*s36**3 - 512*q7a*q9a*s11*s12*s15*s22*s24*s25**2*s36**3 +  \
    1024*q7a*q9a*s11*s14*s16*s22*s24*s25**2*s36**3 - 1024*q7a*q9a*s11*s12*s16*s24**2*s25**2*s36**3 + 512*q8a*q9a*s11*s12**2*s14*s25**3*s36**3 - 1536*q7a*q9a*s11*s12*s14*s22*s25**3*s36**3 + 1024*q7a*q9a*s11*s12**2*s24*s25**3*s36**3 -  \
    1024*q8a*q9a*s11**2*s13**2*s22**2*s26*s36**3 + 6144*q8a*q9a*s11**2*s12*s13*s22*s23*s26*s36**3 - 3072*q8a*q9a*s11**2*s14**2*s22*s23*s26*s36**3 - 2048*q7a*q9a*s11**2*s13*s22**2*s23*s26*s36**3 -  \
    6144*q8a*q9a*s11**2*s12**2*s23**2*s26*s36**3 + 3072*q7a*q9a*s11**2*s12*s22*s23**2*s26*s36**3 + 1024*q8a*q9a*s11**2*s13*s14*s22*s24*s26*s36**3 + 3072*q8a*q9a*s11**2*s12*s14*s23*s24*s26*s36**3 +  \
    1024*q7a*q9a*s11**2*s14*s22*s23*s24*s26*s36**3 - 1024*q8a*q9a*s11**2*s12*s13*s24**2*s26*s36**3 - 1024*q7a*q9a*s11**2*s12*s23*s24**2*s26*s36**3 - 2048*q8a*q9a*s11*s12*s13*s15*s22*s25*s26*s36**3 +  \
    1024*q8a*q9a*s11*s14**2*s15*s22*s25*s26*s36**3 + 1024*q7a*q9a*s11*s13*s15*s22**2*s25*s26*s36**3 + 3072*q8a*q9a*s11*s12**2*s15*s23*s25*s26*s36**3 - 2048*q7a*q9a*s11*s12*s15*s22*s23*s25*s26*s36**3 -  \
    1024*q8a*q9a*s11*s12*s14*s15*s24*s25*s26*s36**3 - 512*q7a*q9a*s11*s14*s15*s22*s24*s25*s26*s36**3 + 512*q7a*q9a*s11*s12*s15*s24**2*s25*s26*s36**3 - 1024*q8a*q9a*s11*s12**2*s13*s25**2*s26*s36**3 +  \
    2048*q7a*q9a*s11*s12*s13*s22*s25**2*s26*s36**3 - 512*q7a*q9a*s11*s14**2*s22*s25**2*s26*s36**3 - 1024*q7a*q9a*s11*s12**2*s23*s25**2*s26*s36**3 + 512*q7a*q9a*s11*s12*s14*s24*s25**2*s26*s36**3 +  \
    1024*q8a**2*s11**2*s14*s15*s22**2*s33*s36**3 - 1024*q8a**2*s11**2*s13*s16*s22**2*s33*s36**3 + 3072*q8a**2*s11**2*s12*s16*s22*s23*s33*s36**3 - 2048*q7a*q8a*s11**2*s16*s22**2*s23*s33*s36**3 -  \
    1536*q8a**2*s11**2*s12*s15*s22*s24*s33*s36**3 + 512*q8a**2*s11**2*s14*s16*s22*s24*s33*s36**3 + 512*q7a*q8a*s11**2*s15*s22**2*s24*s33*s36**3 - 512*q8a**2*s11**2*s12*s16*s24**2*s33*s36**3 -  \
    512*q8a**2*s11**2*s12*s14*s22*s25*s33*s36**3 - 1024*q8a**2*s11*s12*s15*s16*s22*s25*s33*s36**3 - 512*q7a*q8a*s11**2*s14*s22**2*s25*s33*s36**3 + 1024*q7a*q8a*s11*s15*s16*s22**2*s25*s33*s36**3 +  \
    2048*q8a**2*s11**2*s12**2*s24*s25*s33*s36**3 - 2048*q7a*q8a*s11**2*s12*s22*s24*s25*s33*s36**3 + 1024*q7a**2*s11**2*s22**2*s24*s25*s33*s36**3 - 512*q8a**2*s11*s12**2*s16*s25**2*s33*s36**3 +  \
    2048*q7a*q8a*s11*s12*s16*s22*s25**2*s33*s36**3 - 1536*q7a**2*s11*s16*s22**2*s25**2*s33*s36**3 + 3072*q8a**2*s11**2*s12*s13*s22*s26*s33*s36**3 - 1536*q8a**2*s11**2*s14**2*s22*s26*s33*s36**3 -  \
    2048*q7a*q8a*s11**2*s13*s22**2*s26*s33*s36**3 - 6144*q8a**2*s11**2*s12**2*s23*s26*s33*s36**3 + 6144*q7a*q8a*s11**2*s12*s22*s23*s26*s33*s36**3 - 1024*q7a**2*s11**2*s22**2*s23*s26*s33*s36**3 +  \
    1536*q8a**2*s11**2*s12*s14*s24*s26*s33*s36**3 + 1024*q7a*q8a*s11**2*s14*s22*s24*s26*s33*s36**3 - 1024*q7a*q8a*s11**2*s12*s24**2*s26*s33*s36**3 + 1536*q8a**2*s11*s12**2*s15*s25*s26*s33*s36**3 -  \
    2048*q7a*q8a*s11*s12*s15*s22*s25*s26*s33*s36**3 + 512*q7a**2*s11*s15*s22**2*s25*s26*s33*s36**3 - 1024*q7a*q8a*s11*s12**2*s25**2*s26*s33*s36**3 + 1024*q7a**2*s11*s12*s22*s25**2*s26*s33*s36**3 -  \
    1024*q8a**2*s11**2*s13*s15*s22**2*s34*s36**3 + 3072*q8a**2*s11**2*s12*s15*s22*s23*s34*s36**3 - 3072*q8a**2*s11**2*s14*s16*s22*s23*s34*s36**3 - 2048*q7a*q8a*s11**2*s15*s22**2*s23*s34*s36**3 +  \
    512*q8a**2*s11**2*s14*s15*s22*s24*s34*s36**3 + 512*q8a**2*s11**2*s13*s16*s22*s24*s34*s36**3 + 1536*q8a**2*s11**2*s12*s16*s23*s24*s34*s36**3 + 1024*q7a*q8a*s11**2*s16*s22*s23*s24*s34*s36**3 -  \
    512*q8a**2*s11**2*s12*s15*s24**2*s34*s36**3 + 1024*q8a**2*s11**2*s12*s13*s22*s25*s34*s36**3 - 512*q8a**2*s11**2*s14**2*s22*s25*s34*s36**3 - 1024*q8a**2*s11*s12*s15**2*s22*s25*s34*s36**3 +  \
    1024*q8a**2*s11*s14*s15*s16*s22*s25*s34*s36**3 + 1024*q7a*q8a*s11*s15**2*s22**2*s25*s34*s36**3 - 3072*q8a**2*s11**2*s12**2*s23*s25*s34*s36**3 + 2048*q7a*q8a*s11**2*s12*s22*s23*s25*s34*s36**3 +  \
    512*q8a**2*s11**2*s12*s14*s24*s25*s34*s36**3 - 512*q8a**2*s11*s12*s15*s16*s24*s25*s34*s36**3 - 512*q7a*q8a*s11*s15*s16*s22*s24*s25*s34*s36**3 + 1024*q8a**2*s11*s12**2*s15*s25**2*s34*s36**3 -  \
    1024*q7a*q8a*s11*s14*s16*s22*s25**2*s34*s36**3 - 1024*q7a**2*s11*s15*s22**2*s25**2*s34*s36**3 + 512*q7a*q8a*s11*s12*s16*s24*s25**2*s34*s36**3 + 512*q7a**2*s11*s16*s22*s24*s25**2*s34*s36**3 -  \
    1024*q7a*q8a*s11*s12**2*s25**3*s34*s36**3 + 1024*q7a**2*s11*s12*s22*s25**3*s34*s36**3 + 512*q8a**2*s11**2*s13*s14*s22*s26*s34*s36**3 + 1536*q8a**2*s11**2*s12*s14*s23*s26*s34*s36**3 +  \
    1024*q7a*q8a*s11**2*s14*s22*s23*s26*s34*s36**3 - 1024*q8a**2*s11**2*s12*s13*s24*s26*s34*s36**3 - 2048*q7a*q8a*s11**2*s12*s23*s24*s26*s34*s36**3 - 512*q8a**2*s11*s12*s14*s15*s25*s26*s34*s36**3 -  \
    512*q7a*q8a*s11*s14*s15*s22*s25*s26*s34*s36**3 + 1024*q7a*q8a*s11*s12*s15*s24*s25*s26*s34*s36**3 + 512*q7a*q8a*s11*s12*s14*s25**2*s26*s34*s36**3 + 512*q7a**2*s11*s14*s22*s25**2*s26*s34*s36**3 -  \
    1024*q7a**2*s11*s12*s24*s25**2*s26*s34*s36**3 + 2048*q8a**2*s11**2*s12*s13*s22*s24*s35*s36**3 - 1024*q8a**2*s11**2*s14**2*s22*s24*s35*s36**3 - 2048*q7a*q8a*s11**2*s13*s22**2*s24*s35*s36**3 -  \
    3072*q8a**2*s11**2*s12**2*s23*s24*s35*s36**3 + 4096*q7a*q8a*s11**2*s12*s22*s23*s24*s35*s36**3 - 1024*q7a**2*s11**2*s22**2*s23*s24*s35*s36**3 + 1024*q8a**2*s11**2*s12*s14*s24**2*s35*s36**3 +  \
    1024*q7a*q8a*s11**2*s14*s22*s24**2*s35*s36**3 - 1024*q7a*q8a*s11**2*s12*s24**3*s35*s36**3 - 1024*q8a**2*s11*s12*s13*s16*s22*s25*s35*s36**3 + 512*q8a**2*s11*s14**2*s16*s22*s25*s35*s36**3 +  \
    1024*q7a*q8a*s11*s13*s16*s22**2*s25*s35*s36**3 + 1536*q8a**2*s11*s12**2*s16*s23*s25*s35*s36**3 - 2048*q7a*q8a*s11*s12*s16*s22*s23*s25*s35*s36**3 + 512*q7a**2*s11*s16*s22**2*s23*s25*s35*s36**3 +  \
    512*q8a**2*s11*s12**2*s15*s24*s25*s35*s36**3 - 512*q8a**2*s11*s12*s14*s16*s24*s25*s35*s36**3 - 1024*q7a*q8a*s11*s12*s15*s22*s24*s25*s35*s36**3 - 512*q7a*q8a*s11*s14*s16*s22*s24*s25*s35*s36**3 +  \
    512*q7a**2*s11*s15*s22**2*s24*s25*s35*s36**3 + 512*q7a*q8a*s11*s12*s16*s24**2*s25*s35*s36**3 - 512*q8a**2*s11*s12**2*s14*s25**2*s35*s36**3 + 1024*q7a*q8a*s11*s12*s14*s22*s25**2*s35*s36**3 -  \
    512*q7a**2*s11*s14*s22**2*s25**2*s35*s36**3 - 1024*q8a**2*s11*s12*s13*s15*s22*s26*s35*s36**3 + 512*q8a**2*s11*s14**2*s15*s22*s26*s35*s36**3 + 1024*q7a*q8a*s11*s13*s15*s22**2*s26*s35*s36**3 +  \
    1536*q8a**2*s11*s12**2*s15*s23*s26*s35*s36**3 - 2048*q7a*q8a*s11*s12*s15*s22*s23*s26*s35*s36**3 + 512*q7a**2*s11*s15*s22**2*s23*s26*s35*s36**3 - 512*q8a**2*s11*s12*s14*s15*s24*s26*s35*s36**3 -  \
    512*q7a*q8a*s11*s14*s15*s22*s24*s26*s35*s36**3 + 512*q7a*q8a*s11*s12*s15*s24**2*s26*s35*s36**3 - 1024*q8a**2*s11*s12**2*s13*s25*s26*s35*s36**3 + 4096*q7a*q8a*s11*s12*s13*s22*s25*s26*s35*s36**3 -  \
    1024*q7a*q8a*s11*s14**2*s22*s25*s26*s35*s36**3 - 3072*q7a**2*s11*s13*s22**2*s25*s26*s35*s36**3 - 2048*q7a*q8a*s11*s12**2*s23*s25*s26*s35*s36**3 + 2048*q7a**2*s11*s12*s22*s23*s25*s26*s35*s36**3 +  \
    1024*q7a*q8a*s11*s12*s14*s24*s25*s26*s35*s36**3 + 1024*q7a**2*s11*s14*s22*s24*s25*s26*s35*s36**3 - 1024*q7a**2*s11*s12*s24**2*s25*s26*s35*s36**3 + 512*q8a**2*s11**2*s13**2*s22**2*s36**4 -  \
    3072*q8a**2*s11**2*s12*s13*s22*s23*s36**4 + 1536*q8a**2*s11**2*s14**2*s22*s23*s36**4 + 2048*q7a*q8a*s11**2*s13*s22**2*s23*s36**4 + 3072*q8a**2*s11**2*s12**2*s23**2*s36**4 - 3072*q7a*q8a*s11**2*s12*s22*s23**2*s36**4 +  \
    512*q7a**2*s11**2*s22**2*s23**2*s36**4 - 512*q8a**2*s11**2*s13*s14*s22*s24*s36**4 - 1536*q8a**2*s11**2*s12*s14*s23*s24*s36**4 - 1024*q7a*q8a*s11**2*s14*s22*s23*s24*s36**4 + 512*q8a**2*s11**2*s12*s13*s24**2*s36**4 +  \
    1024*q7a*q8a*s11**2*s12*s23*s24**2*s36**4 + 1024*q8a**2*s11*s12*s13*s15*s22*s25*s36**4 - 512*q8a**2*s11*s14**2*s15*s22*s25*s36**4 - 1024*q7a*q8a*s11*s13*s15*s22**2*s25*s36**4 - 1536*q8a**2*s11*s12**2*s15*s23*s25*s36**4 +  \
    2048*q7a*q8a*s11*s12*s15*s22*s23*s25*s36**4 - 512*q7a**2*s11*s15*s22**2*s23*s25*s36**4 + 512*q8a**2*s11*s12*s14*s15*s24*s25*s36**4 + 512*q7a*q8a*s11*s14*s15*s22*s24*s25*s36**4 - 512*q7a*q8a*s11*s12*s15*s24**2*s25*s36**4 +  \
    512*q8a**2*s11*s12**2*s13*s25**2*s36**4 - 2048*q7a*q8a*s11*s12*s13*s22*s25**2*s36**4 + 512*q7a*q8a*s11*s14**2*s22*s25**2*s36**4 + 1536*q7a**2*s11*s13*s22**2*s25**2*s36**4 + 1024*q7a*q8a*s11*s12**2*s23*s25**2*s36**4 -  \
    1024*q7a**2*s11*s12*s22*s23*s25**2*s36**4 - 512*q7a*q8a*s11*s12*s14*s24*s25**2*s36**4 - 512*q7a**2*s11*s14*s22*s24*s25**2*s36**4 + 512*q7a**2*s11*s12*s24**2*s25**2*s36**4
v2=v2_0+v2_1+v2_2+v2_3+v2_4+v2_5

v3_0=-2048*q9a*s11**4*s22**4*s33**3 + 4096*q9a*s11**3*s16*s22**3*s26*s33**3 - 4096*q9a*s11**3*s12*s22**2*s26**2*s33**3 - 2048*q9a*s11**2*s16**2*s22**2*s26**2*s33**3 + 4096*q9a*s11**2*s12*s16*s22*s26**3*s33**3 -  \
    2048*q9a*s11**2*s12**2*s26**4*s33**3 + 3072*q9a*s11**4*s22**3*s24*s33**2*s34 - 1536*q9a*s11**3*s16*s22**3*s25*s33**2*s34 - 1536*q9a*s11**3*s15*s22**3*s26*s33**2*s34 - 4608*q9a*s11**3*s16*s22**2*s24*s26*s33**2*s34 +  \
    3072*q9a*s11**3*s12*s22**2*s25*s26*s33**2*s34 + 1536*q9a*s11**2*s16**2*s22**2*s25*s26*s33**2*s34 + 1536*q9a*s11**3*s14*s22**2*s26**2*s33**2*s34 + 1536*q9a*s11**2*s15*s16*s22**2*s26**2*s33**2*s34 +  \
    3072*q9a*s11**3*s12*s22*s24*s26**2*s33**2*s34 + 1536*q9a*s11**2*s16**2*s22*s24*s26**2*s33**2*s34 - 4608*q9a*s11**2*s12*s16*s22*s25*s26**2*s33**2*s34 - 1536*q9a*s11**2*s12*s15*s22*s26**3*s33**2*s34 -  \
    1536*q9a*s11**2*s14*s16*s22*s26**3*s33**2*s34 - 1536*q9a*s11**2*s12*s16*s24*s26**3*s33**2*s34 + 3072*q9a*s11**2*s12**2*s25*s26**3*s33**2*s34 + 1536*q9a*s11**2*s12*s14*s26**4*s33**2*s34 - 2048*q9a*s11**4*s22**3*s23*s33*s34**2 -  \
    1024*q9a*s11**4*s22**2*s24**2*s33*s34**2 + 1024*q9a*s11**3*s15*s22**3*s25*s33*s34**2 + 1024*q9a*s11**3*s16*s22**2*s24*s25*s33*s34**2 - 1024*q9a*s11**3*s12*s22**2*s25**2*s33*s34**2 + 3072*q9a*s11**3*s16*s22**2*s23*s26*s33*s34**2 +  \
    1024*q9a*s11**3*s15*s22**2*s24*s26*s33*s34**2 + 1024*q9a*s11**3*s16*s22*s24**2*s26*s33*s34**2 - 2048*q9a*s11**3*s14*s22**2*s25*s26*s33*s34**2 - 1024*q9a*s11**2*s15*s16*s22**2*s25*s26*s33*s34**2 -  \
    1024*q9a*s11**2*s16**2*s22*s24*s25*s26*s33*s34**2 + 1024*q9a*s11**2*s12*s16*s22*s25**2*s26*s33*s34**2 - 1024*q9a*s11**3*s13*s22**2*s26**2*s33*s34**2 - 2048*q9a*s11**3*s12*s22*s23*s26**2*s33*s34**2 -  \
    1024*q9a*s11**2*s16**2*s22*s23*s26**2*s33*s34**2 - 1024*q9a*s11**2*s15*s16*s22*s24*s26**2*s33*s34**2 - 1024*q9a*s11**3*s12*s24**2*s26**2*s33*s34**2 + 1024*q9a*s11**2*s12*s15*s22*s25*s26**2*s33*s34**2 +  \
    2048*q9a*s11**2*s14*s16*s22*s25*s26**2*s33*s34**2 + 1024*q9a*s11**2*s12*s16*s24*s25*s26**2*s33*s34**2 - 1024*q9a*s11**2*s12**2*s25**2*s26**2*s33*s34**2 + 1024*q9a*s11**2*s13*s16*s22*s26**3*s33*s34**2 +  \
    1024*q9a*s11**2*s12*s16*s23*s26**3*s33*s34**2 + 1024*q9a*s11**2*s12*s15*s24*s26**3*s33*s34**2 - 2048*q9a*s11**2*s12*s14*s25*s26**3*s33*s34**2 - 1024*q9a*s11**2*s12*s13*s26**4*s33*s34**2 - 1024*q8a*s11**4*s22**3*s33**2*s34**2 +  \
    1536*q8a*s11**3*s16*s22**2*s26*s33**2*s34**2 - 1024*q8a*s11**3*s12*s22*s26**2*s33**2*s34**2 - 512*q8a*s11**2*s16**2*s22*s26**2*s33**2*s34**2 - 512*q7a*s11**3*s22**2*s26**2*s33**2*s34**2 + 512*q8a*s11**2*s12*s16*s26**3*s33**2*s34**2 +  \
    512*q7a*s11**2*s16*s22*s26**3*s33**2*s34**2 - 512*q7a*s11**2*s12*s26**4*s33**2*s34**2 + 1024*q9a*s11**4*s22**2*s23*s24*s34**3 - 512*q9a*s11**3*s16*s22**2*s23*s25*s34**3 - 512*q9a*s11**3*s15*s22**2*s24*s25*s34**3 +  \
    512*q9a*s11**3*s14*s22**2*s25**2*s34**3 - 512*q9a*s11**3*s15*s22**2*s23*s26*s34**3 - 1024*q9a*s11**3*s16*s22*s23*s24*s26*s34**3 + 1024*q9a*s11**3*s13*s22**2*s25*s26*s34**3 + 512*q9a*s11**2*s16**2*s22*s23*s25*s26*s34**3 +  \
    512*q9a*s11**2*s15*s16*s22*s24*s25*s26*s34**3 - 512*q9a*s11**2*s14*s16*s22*s25**2*s26*s34**3 + 512*q9a*s11**2*s15*s16*s22*s23*s26**2*s34**3 + 1024*q9a*s11**3*s12*s23*s24*s26**2*s34**3 -  \
    1024*q9a*s11**2*s13*s16*s22*s25*s26**2*s34**3 - 512*q9a*s11**2*s12*s16*s23*s25*s26**2*s34**3 - 512*q9a*s11**2*s12*s15*s24*s25*s26**2*s34**3 + 512*q9a*s11**2*s12*s14*s25**2*s26**2*s34**3 - 512*q9a*s11**2*s12*s15*s23*s26**3*s34**3 +  \
    1024*q9a*s11**2*s12*s13*s25*s26**3*s34**3 + 1024*q8a*s11**4*s22**2*s24*s33*s34**3 - 512*q8a*s11**3*s16*s22**2*s25*s33*s34**3 - 512*q8a*s11**3*s15*s22**2*s26*s33*s34**3 - 1024*q8a*s11**3*s16*s22*s24*s26*s33*s34**3 +  \
    512*q8a*s11**2*s16**2*s22*s25*s26*s33*s34**3 + 1024*q7a*s11**3*s22**2*s25*s26*s33*s34**3 + 512*q8a*s11**2*s15*s16*s22*s26**2*s33*s34**3 + 1024*q8a*s11**3*s12*s24*s26**2*s33*s34**3 - 512*q8a*s11**2*s12*s16*s25*s26**2*s33*s34**3 -  \
    1024*q7a*s11**2*s16*s22*s25*s26**2*s33*s34**3 - 512*q8a*s11**2*s12*s15*s26**3*s33*s34**3 + 1024*q7a*s11**2*s12*s25*s26**3*s33*s34**3 - 1024*q8a*s11**4*s22**2*s23*s34**4 + 512*q8a*s11**3*s15*s22**2*s25*s34**4 -  \
    512*q7a*s11**3*s22**2*s25**2*s34**4 + 1024*q8a*s11**3*s16*s22*s23*s26*s34**4 - 512*q8a*s11**2*s15*s16*s22*s25*s26*s34**4 + 512*q7a*s11**2*s16*s22*s25**2*s26*s34**4 - 1024*q8a*s11**3*s12*s23*s26**2*s34**4 +  \
    512*q8a*s11**2*s12*s15*s25*s26**2*s34**4 - 512*q7a*s11**2*s12*s25**2*s26**2*s34**4 + 3072*q9a*s11**3*s15*s22**4*s33**2*s35 - 1536*q9a*s11**3*s16*s22**3*s24*s33**2*s35 - 3072*q9a*s11**3*s12*s22**3*s25*s33**2*s35 +  \
    1536*q9a*s11**2*s16**2*s22**3*s25*s33**2*s35 - 1536*q9a*s11**3*s14*s22**3*s26*s33**2*s35 - 4608*q9a*s11**2*s15*s16*s22**3*s26*s33**2*s35 + 3072*q9a*s11**3*s12*s22**2*s24*s26*s33**2*s35 +  \
    1536*q9a*s11**2*s16**2*s22**2*s24*s26*s33**2*s35 + 1536*q9a*s11**2*s12*s16*s22**2*s25*s26*s33**2*s35 - 1536*q9a*s11*s16**3*s22**2*s25*s26*s33**2*s35 + 4608*q9a*s11**2*s12*s15*s22**2*s26**2*s33**2*s35 +  \
    1536*q9a*s11**2*s14*s16*s22**2*s26**2*s33**2*s35 + 1536*q9a*s11*s15*s16**2*s22**2*s26**2*s33**2*s35 - 4608*q9a*s11**2*s12*s16*s22*s24*s26**2*s33**2*s35 - 3072*q9a*s11**2*s12**2*s22*s25*s26**2*s33**2*s35 +  \
    3072*q9a*s11*s12*s16**2*s22*s25*s26**2*s33**2*s35 - 1536*q9a*s11**2*s12*s14*s22*s26**3*s33**2*s35 - 3072*q9a*s11*s12*s15*s16*s22*s26**3*s33**2*s35 + 3072*q9a*s11**2*s12**2*s24*s26**3*s33**2*s35 -  \
    1536*q9a*s11*s12**2*s16*s25*s26**3*s33**2*s35 + 1536*q9a*s11*s12**2*s15*s26**4*s33**2*s35 + 2048*q9a*s11**3*s16*s22**3*s23*s33*s34*s35 - 3072*q9a*s11**3*s15*s22**3*s24*s33*s34*s35 + 1024*q9a*s11**3*s16*s22**2*s24**2*s33*s34*s35 +  \
    3072*q9a*s11**3*s14*s22**3*s25*s33*s34*s35 - 1024*q9a*s11**2*s15*s16*s22**3*s25*s33*s34*s35 - 1024*q9a*s11**2*s16**2*s22**2*s24*s25*s33*s34*s35 + 1024*q9a*s11**2*s12*s16*s22**2*s25**2*s33*s34*s35 +  \
    2048*q9a*s11**3*s13*s22**3*s26*s33*s34*s35 + 1024*q9a*s11**2*s15**2*s22**3*s26*s33*s34*s35 - 4096*q9a*s11**3*s12*s22**2*s23*s26*s33*s34*s35 - 2048*q9a*s11**2*s16**2*s22**2*s23*s26*s33*s34*s35 -  \
    1024*q9a*s11**3*s14*s22**2*s24*s26*s33*s34*s35 + 4096*q9a*s11**2*s15*s16*s22**2*s24*s26*s33*s34*s35 - 1024*q9a*s11**2*s16**2*s22*s24**2*s26*s33*s34*s35 - 1024*q9a*s11**2*s12*s15*s22**2*s25*s26*s33*s34*s35 -  \
    2048*q9a*s11**2*s14*s16*s22**2*s25*s26*s33*s34*s35 + 1024*q9a*s11*s15*s16**2*s22**2*s25*s26*s33*s34*s35 + 1024*q9a*s11*s16**3*s22*s24*s25*s26*s33*s34*s35 - 1024*q9a*s11*s12*s16**2*s22*s25**2*s26*s33*s34*s35 -  \
    1024*q9a*s11**2*s14*s15*s22**2*s26**2*s33*s34*s35 - 2048*q9a*s11**2*s13*s16*s22**2*s26**2*s33*s34*s35 - 1024*q9a*s11*s15**2*s16*s22**2*s26**2*s33*s34*s35 + 6144*q9a*s11**2*s12*s16*s22*s23*s26**2*s33*s34*s35 -  \
    3072*q9a*s11**2*s12*s15*s22*s24*s26**2*s33*s34*s35 + 1024*q9a*s11**2*s14*s16*s22*s24*s26**2*s33*s34*s35 - 1024*q9a*s11*s15*s16**2*s22*s24*s26**2*s33*s34*s35 + 1024*q9a*s11**2*s12*s16*s24**2*s26**2*s33*s34*s35 +  \
    3072*q9a*s11**2*s12*s14*s22*s25*s26**2*s33*s34*s35 - 1024*q9a*s11*s14*s16**2*s22*s25*s26**2*s33*s34*s35 - 1024*q9a*s11*s12*s16**2*s24*s25*s26**2*s33*s34*s35 + 1024*q9a*s11*s12**2*s16*s25**2*s26**2*s33*s34*s35 +  \
    2048*q9a*s11**2*s12*s13*s22*s26**3*s33*s34*s35 + 1024*q9a*s11*s12*s15**2*s22*s26**3*s33*s34*s35 + 1024*q9a*s11*s14*s15*s16*s22*s26**3*s33*s34*s35 - 4096*q9a*s11**2*s12**2*s23*s26**3*s33*s34*s35 -  \
    1024*q9a*s11**2*s12*s14*s24*s26**3*s33*s34*s35 + 1024*q9a*s11*s12*s15*s16*s24*s26**3*s33*s34*s35 - 1024*q9a*s11*s12**2*s15*s25*s26**3*s33*s34*s35 + 1024*q9a*s11*s12*s14*s16*s25*s26**3*s33*s34*s35 -  \
    1024*q9a*s11*s12*s14*s15*s26**4*s33*s34*s35 + 1024*q8a*s11**3*s16*s22**3*s33**2*s34*s35 - 2048*q8a*s11**3*s12*s22**2*s26*s33**2*s34*s35 - 1024*q8a*s11**2*s16**2*s22**2*s26*s33**2*s34*s35 +  \
    1024*q7a*s11**3*s22**3*s26*s33**2*s34*s35 + 3072*q8a*s11**2*s12*s16*s22*s26**2*s33**2*s34*s35 - 1024*q7a*s11**2*s16*s22**2*s26**2*s33**2*s34*s35 - 2048*q8a*s11**2*s12**2*s26**3*s33**2*s34*s35 +  \
    1024*q7a*s11**2*s12*s22*s26**3*s33**2*s34*s35 + 1024*q9a*s11**3*s15*s22**3*s23*s34**2*s35 - 1536*q9a*s11**3*s16*s22**2*s23*s24*s34**2*s35 + 512*q9a*s11**3*s15*s22**2*s24**2*s34**2*s35 - 2048*q9a*s11**3*s13*s22**3*s25*s34**2*s35 +  \
    1024*q9a*s11**3*s12*s22**2*s23*s25*s34**2*s35 + 512*q9a*s11**2*s16**2*s22**2*s23*s25*s34**2*s35 - 512*q9a*s11**3*s14*s22**2*s24*s25*s34**2*s35 + 512*q9a*s11**2*s15*s16*s22**2*s24*s25*s34**2*s35 -  \
    512*q9a*s11**2*s14*s16*s22**2*s25**2*s34**2*s35 + 1536*q9a*s11**3*s14*s22**2*s23*s26*s34**2*s35 - 1536*q9a*s11**2*s15*s16*s22**2*s23*s26*s34**2*s35 - 512*q9a*s11**2*s15**2*s22**2*s24*s26*s34**2*s35 +  \
    1536*q9a*s11**2*s16**2*s22*s23*s24*s26*s34**2*s35 - 512*q9a*s11**2*s15*s16*s22*s24**2*s26*s34**2*s35 + 512*q9a*s11**2*s14*s15*s22**2*s25*s26*s34**2*s35 + 1536*q9a*s11**2*s13*s16*s22**2*s25*s26*s34**2*s35 -  \
    1024*q9a*s11**2*s12*s16*s22*s23*s25*s26*s34**2*s35 - 512*q9a*s11*s16**3*s22*s23*s25*s26*s34**2*s35 + 512*q9a*s11**2*s14*s16*s22*s24*s25*s26*s34**2*s35 - 512*q9a*s11*s15*s16**2*s22*s24*s25*s26*s34**2*s35 +  \
    512*q9a*s11*s14*s16**2*s22*s25**2*s26*s34**2*s35 + 512*q9a*s11**2*s13*s15*s22**2*s26**2*s34**2*s35 + 1024*q9a*s11**2*s12*s15*s22*s23*s26**2*s34**2*s35 - 1536*q9a*s11**2*s14*s16*s22*s23*s26**2*s34**2*s35 +  \
    512*q9a*s11*s15*s16**2*s22*s23*s26**2*s34**2*s35 + 512*q9a*s11*s15**2*s16*s22*s24*s26**2*s34**2*s35 - 1536*q9a*s11**2*s12*s16*s23*s24*s26**2*s34**2*s35 + 512*q9a*s11**2*s12*s15*s24**2*s26**2*s34**2*s35 -  \
    2048*q9a*s11**2*s12*s13*s22*s25*s26**2*s34**2*s35 - 512*q9a*s11*s14*s15*s16*s22*s25*s26**2*s34**2*s35 + 512*q9a*s11*s13*s16**2*s22*s25*s26**2*s34**2*s35 + 1024*q9a*s11**2*s12**2*s23*s25*s26**2*s34**2*s35 +  \
    512*q9a*s11*s12*s16**2*s23*s25*s26**2*s34**2*s35 - 512*q9a*s11**2*s12*s14*s24*s25*s26**2*s34**2*s35 + 512*q9a*s11*s12*s15*s16*s24*s25*s26**2*s34**2*s35 - 512*q9a*s11*s12*s14*s16*s25**2*s26**2*s34**2*s35 -  \
    512*q9a*s11*s13*s15*s16*s22*s26**3*s34**2*s35 + 1536*q9a*s11**2*s12*s14*s23*s26**3*s34**2*s35 - 512*q9a*s11*s12*s15*s16*s23*s26**3*s34**2*s35 - 512*q9a*s11*s12*s15**2*s24*s26**3*s34**2*s35 +  \
    512*q9a*s11*s12*s14*s15*s25*s26**3*s34**2*s35 - 512*q9a*s11*s12*s13*s16*s25*s26**3*s34**2*s35 + 512*q9a*s11*s12*s13*s15*s26**4*s34**2*s35 + 1024*q8a*s11**3*s15*s22**3*s33*s34**2*s35 -  \
    1536*q8a*s11**3*s16*s22**2*s24*s33*s34**2*s35 + 1024*q8a*s11**3*s12*s22**2*s25*s33*s34**2*s35 + 512*q8a*s11**2*s16**2*s22**2*s25*s33*s34**2*s35 - 2048*q7a*s11**3*s22**3*s25*s33*s34**2*s35 +  \
    1536*q8a*s11**3*s14*s22**2*s26*s33*s34**2*s35 - 1536*q8a*s11**2*s15*s16*s22**2*s26*s33*s34**2*s35 + 1536*q8a*s11**2*s16**2*s22*s24*s26*s33*s34**2*s35 - 1024*q8a*s11**2*s12*s16*s22*s25*s26*s33*s34**2*s35 -  \
    512*q8a*s11*s16**3*s22*s25*s26*s33*s34**2*s35 + 1536*q7a*s11**2*s16*s22**2*s25*s26*s33*s34**2*s35 + 1024*q8a*s11**2*s12*s15*s22*s26**2*s33*s34**2*s35 - 1536*q8a*s11**2*s14*s16*s22*s26**2*s33*s34**2*s35 +  \
    512*q8a*s11*s15*s16**2*s22*s26**2*s33*s34**2*s35 + 512*q7a*s11**2*s15*s22**2*s26**2*s33*s34**2*s35 - 1536*q8a*s11**2*s12*s16*s24*s26**2*s33*s34**2*s35 + 1024*q8a*s11**2*s12**2*s25*s26**2*s33*s34**2*s35 +  \
    512*q8a*s11*s12*s16**2*s25*s26**2*s33*s34**2*s35 - 2048*q7a*s11**2*s12*s22*s25*s26**2*s33*s34**2*s35 + 512*q7a*s11*s16**2*s22*s25*s26**2*s33*s34**2*s35 + 1536*q8a*s11**2*s12*s14*s26**3*s33*s34**2*s35 -  \
    512*q8a*s11*s12*s15*s16*s26**3*s33*s34**2*s35 - 512*q7a*s11*s15*s16*s22*s26**3*s33*s34**2*s35 - 512*q7a*s11*s12*s16*s25*s26**3*s33*s34**2*s35 + 512*q7a*s11*s12*s15*s26**4*s33*s34**2*s35 +  \
    2048*q8a*s11**3*s16*s22**2*s23*s34**3*s35 - 512*q8a*s11**3*s15*s22**2*s24*s34**3*s35 - 512*q8a*s11**3*s14*s22**2*s25*s34**3*s35 - 512*q8a*s11**2*s15*s16*s22**2*s25*s34**3*s35 + 1024*q7a*s11**3*s22**2*s24*s25*s34**3*s35 +  \
    512*q7a*s11**2*s16*s22**2*s25**2*s34**3*s35 - 1024*q8a*s11**3*s13*s22**2*s26*s34**3*s35 + 512*q8a*s11**2*s15**2*s22**2*s26*s34**3*s35 - 2048*q8a*s11**2*s16**2*s22*s23*s26*s34**3*s35 - 1024*q7a*s11**3*s22**2*s23*s26*s34**3*s35
v3_1= \
    512*q8a*s11**2*s15*s16*s22*s24*s26*s34**3*s35 + 512*q8a*s11**2*s14*s16*s22*s25*s26*s34**3*s35 + 512*q8a*s11*s15*s16**2*s22*s25*s26*s34**3*s35 - 512*q7a*s11**2*s15*s22**2*s25*s26*s34**3*s35 -  \
    1024*q7a*s11**2*s16*s22*s24*s25*s26*s34**3*s35 - 512*q7a*s11*s16**2*s22*s25**2*s26*s34**3*s35 + 1024*q8a*s11**2*s13*s16*s22*s26**2*s34**3*s35 - 512*q8a*s11*s15**2*s16*s22*s26**2*s34**3*s35 +  \
    2048*q8a*s11**2*s12*s16*s23*s26**2*s34**3*s35 + 1024*q7a*s11**2*s16*s22*s23*s26**2*s34**3*s35 - 512*q8a*s11**2*s12*s15*s24*s26**2*s34**3*s35 - 512*q8a*s11**2*s12*s14*s25*s26**2*s34**3*s35 -  \
    512*q8a*s11*s12*s15*s16*s25*s26**2*s34**3*s35 + 512*q7a*s11*s15*s16*s22*s25*s26**2*s34**3*s35 + 1024*q7a*s11**2*s12*s24*s25*s26**2*s34**3*s35 + 512*q7a*s11*s12*s16*s25**2*s26**2*s34**3*s35 -  \
    1024*q8a*s11**2*s12*s13*s26**3*s34**3*s35 + 512*q8a*s11*s12*s15**2*s26**3*s34**3*s35 - 1024*q7a*s11**2*s12*s23*s26**3*s34**3*s35 - 512*q7a*s11*s12*s15*s25*s26**3*s34**3*s35 - 2048*q9a*s11**3*s13*s22**4*s33*s35**2 -  \
    1024*q9a*s11**2*s15**2*s22**4*s33*s35**2 + 2048*q9a*s11**3*s12*s22**3*s23*s33*s35**2 - 1024*q9a*s11**2*s16**2*s22**3*s23*s33*s35**2 + 1024*q9a*s11**3*s14*s22**3*s24*s33*s35**2 + 1024*q9a*s11**2*s15*s16*s22**3*s24*s33*s35**2 -  \
    1024*q9a*s11**3*s12*s22**2*s24**2*s33*s35**2 + 2048*q9a*s11**2*s12*s15*s22**3*s25*s33*s35**2 - 2048*q9a*s11**2*s14*s16*s22**3*s25*s33*s35**2 + 1024*q9a*s11**2*s12*s16*s22**2*s24*s25*s33*s35**2 -  \
    1024*q9a*s11**2*s12**2*s22**2*s25**2*s33*s35**2 + 1024*q9a*s11**2*s14*s15*s22**3*s26*s33*s35**2 + 3072*q9a*s11**2*s13*s16*s22**3*s26*s33*s35**2 + 1024*q9a*s11*s15**2*s16*s22**3*s26*s33*s35**2 -  \
    1024*q9a*s11**2*s12*s16*s22**2*s23*s26*s33*s35**2 + 1024*q9a*s11*s16**3*s22**2*s23*s26*s33*s35**2 - 2048*q9a*s11**2*s12*s15*s22**2*s24*s26*s33*s35**2 - 1024*q9a*s11**2*s14*s16*s22**2*s24*s26*s33*s35**2 -  \
    1024*q9a*s11*s15*s16**2*s22**2*s24*s26*s33*s35**2 + 1024*q9a*s11**2*s12*s16*s22*s24**2*s26*s33*s35**2 + 1024*q9a*s11**2*s12*s14*s22**2*s25*s26*s33*s35**2 - 2048*q9a*s11*s12*s15*s16*s22**2*s25*s26*s33*s35**2 +  \
    2048*q9a*s11*s14*s16**2*s22**2*s25*s26*s33*s35**2 - 1024*q9a*s11*s12*s16**2*s22*s24*s25*s26*s33*s35**2 + 1024*q9a*s11*s12**2*s16*s22*s25**2*s26*s33*s35**2 - 3072*q9a*s11**2*s12*s13*s22**2*s26**2*s33*s35**2 -  \
    1024*q9a*s11*s12*s15**2*s22**2*s26**2*s33*s35**2 - 1024*q9a*s11*s14*s15*s16*s22**2*s26**2*s33*s35**2 - 1024*q9a*s11*s13*s16**2*s22**2*s26**2*s33*s35**2 + 2048*q9a*s11**2*s12**2*s22*s23*s26**2*s33*s35**2 -  \
    2048*q9a*s11*s12*s16**2*s22*s23*s26**2*s33*s35**2 + 1024*q9a*s11**2*s12*s14*s22*s24*s26**2*s33*s35**2 + 3072*q9a*s11*s12*s15*s16*s22*s24*s26**2*s33*s35**2 - 1024*q9a*s11**2*s12**2*s24**2*s26**2*s33*s35**2 +  \
    2048*q9a*s11*s12**2*s15*s22*s25*s26**2*s33*s35**2 - 3072*q9a*s11*s12*s14*s16*s22*s25*s26**2*s33*s35**2 + 1024*q9a*s11*s12**2*s16*s24*s25*s26**2*s33*s35**2 - 1024*q9a*s11*s12**3*s25**2*s26**2*s33*s35**2 +  \
    1024*q9a*s11*s12*s14*s15*s22*s26**3*s33*s35**2 + 2048*q9a*s11*s12*s13*s16*s22*s26**3*s33*s35**2 + 1024*q9a*s11*s12**2*s16*s23*s26**3*s33*s35**2 - 2048*q9a*s11*s12**2*s15*s24*s26**3*s33*s35**2 +  \
    1024*q9a*s11*s12**2*s14*s25*s26**3*s33*s35**2 - 1024*q9a*s11*s12**2*s13*s26**4*s33*s35**2 + 1024*q8a*s11**3*s12*s22**3*s33**2*s35**2 - 512*q8a*s11**2*s16**2*s22**3*s33**2*s35**2 - 1024*q7a*s11**3*s22**4*s33**2*s35**2 -  \
    512*q8a*s11**2*s12*s16*s22**2*s26*s33**2*s35**2 + 512*q8a*s11*s16**3*s22**2*s26*s33**2*s35**2 + 1536*q7a*s11**2*s16*s22**3*s26*s33**2*s35**2 + 1024*q8a*s11**2*s12**2*s22*s26**2*s33**2*s35**2 -  \
    1024*q8a*s11*s12*s16**2*s22*s26**2*s33**2*s35**2 - 1536*q7a*s11**2*s12*s22**2*s26**2*s33**2*s35**2 - 512*q7a*s11*s16**2*s22**2*s26**2*s33**2*s35**2 + 512*q8a*s11*s12**2*s16*s26**3*s33**2*s35**2 +  \
    1024*q7a*s11*s12*s16*s22*s26**3*s33**2*s35**2 - 512*q7a*s11*s12**2*s26**4*s33**2*s35**2 - 2048*q9a*s11**3*s14*s22**3*s23*s34*s35**2 + 1024*q9a*s11**3*s13*s22**3*s24*s34*s35**2 + 512*q9a*s11**2*s15**2*s22**3*s24*s34*s35**2 +  \
    1024*q9a*s11**3*s12*s22**2*s23*s24*s34*s35**2 + 512*q9a*s11**2*s16**2*s22**2*s23*s24*s34*s35**2 - 512*q9a*s11**2*s15*s16*s22**2*s24**2*s34*s35**2 - 512*q9a*s11**2*s14*s15*s22**3*s25*s34*s35**2 +  \
    1536*q9a*s11**2*s13*s16*s22**3*s25*s34*s35**2 - 1536*q9a*s11**2*s12*s16*s22**2*s23*s25*s34*s35**2 - 512*q9a*s11**2*s12*s15*s22**2*s24*s25*s34*s35**2 + 512*q9a*s11**2*s14*s16*s22**2*s24*s25*s34*s35**2 +  \
    512*q9a*s11**2*s12*s14*s22**2*s25**2*s34*s35**2 - 1536*q9a*s11**2*s13*s15*s22**3*s26*s34*s35**2 + 1536*q9a*s11**2*s12*s15*s22**2*s23*s26*s34*s35**2 + 1536*q9a*s11**2*s14*s16*s22**2*s23*s26*s34*s35**2 +  \
    512*q9a*s11**2*s14*s15*s22**2*s24*s26*s34*s35**2 - 1536*q9a*s11**2*s13*s16*s22**2*s24*s26*s34*s35**2 - 512*q9a*s11*s15**2*s16*s22**2*s24*s26*s34*s35**2 - 1024*q9a*s11**2*s12*s16*s22*s23*s24*s26*s34*s35**2 -  \
    512*q9a*s11*s16**3*s22*s23*s24*s26*s34*s35**2 + 512*q9a*s11*s15*s16**2*s22*s24**2*s26*s34*s35**2 - 512*q9a*s11**2*s14**2*s22**2*s25*s26*s34*s35**2 + 512*q9a*s11*s14*s15*s16*s22**2*s25*s26*s34*s35**2 -  \
    1536*q9a*s11*s13*s16**2*s22**2*s25*s26*s34*s35**2 + 1536*q9a*s11*s12*s16**2*s22*s23*s25*s26*s34*s35**2 + 512*q9a*s11*s12*s15*s16*s22*s24*s25*s26*s34*s35**2 - 512*q9a*s11*s14*s16**2*s22*s24*s25*s26*s34*s35**2 -  \
    512*q9a*s11*s12*s14*s16*s22*s25**2*s26*s34*s35**2 + 512*q9a*s11**2*s13*s14*s22**2*s26**2*s34*s35**2 + 1536*q9a*s11*s13*s15*s16*s22**2*s26**2*s34*s35**2 - 2048*q9a*s11**2*s12*s14*s22*s23*s26**2*s34*s35**2 -  \
    1536*q9a*s11*s12*s15*s16*s22*s23*s26**2*s34*s35**2 + 512*q9a*s11*s14*s16**2*s22*s23*s26**2*s34*s35**2 + 1024*q9a*s11**2*s12*s13*s22*s24*s26**2*s34*s35**2 + 512*q9a*s11*s12*s15**2*s22*s24*s26**2*s34*s35**2 -  \
    512*q9a*s11*s14*s15*s16*s22*s24*s26**2*s34*s35**2 + 512*q9a*s11*s13*s16**2*s22*s24*s26**2*s34*s35**2 + 1024*q9a*s11**2*s12**2*s23*s24*s26**2*s34*s35**2 + 512*q9a*s11*s12*s16**2*s23*s24*s26**2*s34*s35**2 -  \
    512*q9a*s11*s12*s15*s16*s24**2*s26**2*s34*s35**2 - 512*q9a*s11*s12*s14*s15*s22*s25*s26**2*s34*s35**2 + 1536*q9a*s11*s12*s13*s16*s22*s25*s26**2*s34*s35**2 + 512*q9a*s11*s14**2*s16*s22*s25*s26**2*s34*s35**2 -  \
    1536*q9a*s11*s12**2*s16*s23*s25*s26**2*s34*s35**2 - 512*q9a*s11*s12**2*s15*s24*s25*s26**2*s34*s35**2 + 512*q9a*s11*s12*s14*s16*s24*s25*s26**2*s34*s35**2 + 512*q9a*s11*s12**2*s14*s25**2*s26**2*s34*s35**2 -  \
    1536*q9a*s11*s12*s13*s15*s22*s26**3*s34*s35**2 - 512*q9a*s11*s13*s14*s16*s22*s26**3*s34*s35**2 + 1536*q9a*s11*s12**2*s15*s23*s26**3*s34*s35**2 - 512*q9a*s11*s12*s14*s16*s23*s26**3*s34*s35**2 +  \
    512*q9a*s11*s12*s14*s15*s24*s26**3*s34*s35**2 - 512*q9a*s11*s12*s13*s16*s24*s26**3*s34*s35**2 - 512*q9a*s11*s12*s14**2*s25*s26**3*s34*s35**2 + 512*q9a*s11*s12*s13*s14*s26**4*s34*s35**2 -  \
    2048*q8a*s11**3*s14*s22**3*s33*s34*s35**2 + 1024*q8a*s11**3*s12*s22**2*s24*s33*s34*s35**2 + 512*q8a*s11**2*s16**2*s22**2*s24*s33*s34*s35**2 + 1024*q7a*s11**3*s22**3*s24*s33*s34*s35**2 -  \
    1536*q8a*s11**2*s12*s16*s22**2*s25*s33*s34*s35**2 + 1536*q7a*s11**2*s16*s22**3*s25*s33*s34*s35**2 + 1536*q8a*s11**2*s12*s15*s22**2*s26*s33*s34*s35**2 + 1536*q8a*s11**2*s14*s16*s22**2*s26*s33*s34*s35**2 -  \
    1536*q7a*s11**2*s15*s22**3*s26*s33*s34*s35**2 - 1024*q8a*s11**2*s12*s16*s22*s24*s26*s33*s34*s35**2 - 512*q8a*s11*s16**3*s22*s24*s26*s33*s34*s35**2 - 1536*q7a*s11**2*s16*s22**2*s24*s26*s33*s34*s35**2 +  \
    1536*q8a*s11*s12*s16**2*s22*s25*s26*s33*s34*s35**2 - 1536*q7a*s11*s16**2*s22**2*s25*s26*s33*s34*s35**2 - 2048*q8a*s11**2*s12*s14*s22*s26**2*s33*s34*s35**2 - 1536*q8a*s11*s12*s15*s16*s22*s26**2*s33*s34*s35**2 +  \
    512*q8a*s11*s14*s16**2*s22*s26**2*s33*s34*s35**2 + 512*q7a*s11**2*s14*s22**2*s26**2*s33*s34*s35**2 + 1536*q7a*s11*s15*s16*s22**2*s26**2*s33*s34*s35**2 + 1024*q8a*s11**2*s12**2*s24*s26**2*s33*s34*s35**2 +  \
    512*q8a*s11*s12*s16**2*s24*s26**2*s33*s34*s35**2 + 1024*q7a*s11**2*s12*s22*s24*s26**2*s33*s34*s35**2 + 512*q7a*s11*s16**2*s22*s24*s26**2*s33*s34*s35**2 - 1536*q8a*s11*s12**2*s16*s25*s26**2*s33*s34*s35**2 +  \
    1536*q7a*s11*s12*s16*s22*s25*s26**2*s33*s34*s35**2 + 1536*q8a*s11*s12**2*s15*s26**3*s33*s34*s35**2 - 512*q8a*s11*s12*s14*s16*s26**3*s33*s34*s35**2 - 1536*q7a*s11*s12*s15*s22*s26**3*s33*s34*s35**2 -  \
    512*q7a*s11*s14*s16*s22*s26**3*s33*s34*s35**2 - 512*q7a*s11*s12*s16*s24*s26**3*s33*s34*s35**2 + 512*q7a*s11*s12*s14*s26**4*s33*s34*s35**2 + 1024*q8a*s11**3*s13*s22**3*s34**2*s35**2 - 512*q8a*s11**2*s15**2*s22**3*s34**2*s35**2 -  \
    2048*q8a*s11**3*s12*s22**2*s23*s34**2*s35**2 - 1024*q8a*s11**2*s16**2*s22**2*s23*s34**2*s35**2 + 1024*q7a*s11**3*s22**3*s23*s34**2*s35**2 + 512*q8a*s11**3*s14*s22**2*s24*s34**2*s35**2 + 512*q8a*s11**2*s15*s16*s22**2*s24*s34**2*s35**2 -  \
    512*q7a*s11**3*s22**2*s24**2*s34**2*s35**2 + 512*q8a*s11**2*s12*s15*s22**2*s25*s34**2*s35**2 + 512*q8a*s11**2*s14*s16*s22**2*s25*s34**2*s35**2 + 512*q7a*s11**2*s15*s22**3*s25*s34**2*s35**2 -  \
    1024*q7a*s11**2*s16*s22**2*s24*s25*s34**2*s35**2 - 512*q7a*s11**2*s12*s22**2*s25**2*s34**2*s35**2 - 1024*q8a*s11**2*s14*s15*s22**2*s26*s34**2*s35**2 + 512*q8a*s11*s15**2*s16*s22**2*s26*s34**2*s35**2 +  \
    2048*q8a*s11**2*s12*s16*s22*s23*s26*s34**2*s35**2 + 1024*q8a*s11*s16**3*s22*s23*s26*s34**2*s35**2 - 512*q8a*s11**2*s14*s16*s22*s24*s26*s34**2*s35**2 - 512*q8a*s11*s15*s16**2*s22*s24*s26*s34**2*s35**2 +  \
    512*q7a*s11**2*s15*s22**2*s24*s26*s34**2*s35**2 + 512*q7a*s11**2*s16*s22*s24**2*s26*s34**2*s35**2 - 512*q8a*s11*s12*s15*s16*s22*s25*s26*s34**2*s35**2 - 512*q8a*s11*s14*s16**2*s22*s25*s26*s34**2*s35**2 +  \
    512*q7a*s11**2*s14*s22**2*s25*s26*s34**2*s35**2 - 512*q7a*s11*s15*s16*s22**2*s25*s26*s34**2*s35**2 + 1024*q7a*s11*s16**2*s22*s24*s25*s26*s34**2*s35**2 + 512*q7a*s11*s12*s16*s22*s25**2*s26*s34**2*s35**2 +  \
    1024*q8a*s11**2*s12*s13*s22*s26**2*s34**2*s35**2 - 512*q8a*s11*s12*s15**2*s22*s26**2*s34**2*s35**2 + 1024*q8a*s11*s14*s15*s16*s22*s26**2*s34**2*s35**2 - 1024*q8a*s11*s13*s16**2*s22*s26**2*s34**2*s35**2 -  \
    1024*q7a*s11**2*s13*s22**2*s26**2*s34**2*s35**2 - 2048*q8a*s11**2*s12**2*s23*s26**2*s34**2*s35**2 - 1024*q8a*s11*s12*s16**2*s23*s26**2*s34**2*s35**2 + 1024*q7a*s11**2*s12*s22*s23*s26**2*s34**2*s35**2 -  \
    1024*q7a*s11*s16**2*s22*s23*s26**2*s34**2*s35**2 + 512*q8a*s11**2*s12*s14*s24*s26**2*s34**2*s35**2 + 512*q8a*s11*s12*s15*s16*s24*s26**2*s34**2*s35**2 - 512*q7a*s11*s15*s16*s22*s24*s26**2*s34**2*s35**2 -  \
    512*q7a*s11**2*s12*s24**2*s26**2*s34**2*s35**2 + 512*q8a*s11*s12**2*s15*s25*s26**2*s34**2*s35**2 + 512*q8a*s11*s12*s14*s16*s25*s26**2*s34**2*s35**2 + 512*q7a*s11*s12*s15*s22*s25*s26**2*s34**2*s35**2 -  \
    512*q7a*s11*s14*s16*s22*s25*s26**2*s34**2*s35**2 - 1024*q7a*s11*s12*s16*s24*s25*s26**2*s34**2*s35**2 - 512*q7a*s11*s12**2*s25**2*s26**2*s34**2*s35**2 - 1024*q8a*s11*s12*s14*s15*s26**3*s34**2*s35**2 +  \
    1024*q8a*s11*s12*s13*s16*s26**3*s34**2*s35**2 + 1024*q7a*s11*s13*s16*s22*s26**3*s34**2*s35**2 + 1024*q7a*s11*s12*s16*s23*s26**3*s34**2*s35**2 + 512*q7a*s11*s12*s15*s24*s26**3*s34**2*s35**2 +  \
    512*q7a*s11*s12*s14*s25*s26**3*s34**2*s35**2 - 1024*q7a*s11*s12*s13*s26**4*s34**2*s35**2 + 1024*q9a*s11**2*s13*s15*s22**4*s35**3 - 1024*q9a*s11**2*s12*s15*s22**3*s23*s35**3 + 1024*q9a*s11**2*s14*s16*s22**3*s23*s35**3 -  \
    512*q9a*s11**2*s14*s15*s22**3*s24*s35**3 - 512*q9a*s11**2*s13*s16*s22**3*s24*s35**3 - 512*q9a*s11**2*s12*s16*s22**2*s23*s24*s35**3 + 512*q9a*s11**2*s12*s15*s22**2*s24**2*s35**3 - 1024*q9a*s11**2*s12*s13*s22**3*s25*s35**3 +  \
    512*q9a*s11**2*s14**2*s22**3*s25*s35**3 + 1024*q9a*s11**2*s12**2*s22**2*s23*s25*s35**3 - 512*q9a*s11**2*s12*s14*s22**2*s24*s25*s35**3 - 512*q9a*s11**2*s13*s14*s22**3*s26*s35**3 - 1024*q9a*s11*s13*s15*s16*s22**3*s26*s35**3 -  \
    512*q9a*s11**2*s12*s14*s22**2*s23*s26*s35**3 + 1024*q9a*s11*s12*s15*s16*s22**2*s23*s26*s35**3 - 1024*q9a*s11*s14*s16**2*s22**2*s23*s26*s35**3 + 1024*q9a*s11**2*s12*s13*s22**2*s24*s26*s35**3 +  \
    512*q9a*s11*s14*s15*s16*s22**2*s24*s26*s35**3 + 512*q9a*s11*s13*s16**2*s22**2*s24*s26*s35**3 + 512*q9a*s11*s12*s16**2*s22*s23*s24*s26*s35**3 - 512*q9a*s11*s12*s15*s16*s22*s24**2*s26*s35**3 +  \
    1024*q9a*s11*s12*s13*s16*s22**2*s25*s26*s35**3 - 512*q9a*s11*s14**2*s16*s22**2*s25*s26*s35**3 - 1024*q9a*s11*s12**2*s16*s22*s23*s25*s26*s35**3 + 512*q9a*s11*s12*s14*s16*s22*s24*s25*s26*s35**3
v3_2= \
    1024*q9a*s11*s12*s13*s15*s22**2*s26**2*s35**3 + 512*q9a*s11*s13*s14*s16*s22**2*s26**2*s35**3 - 1024*q9a*s11*s12**2*s15*s22*s23*s26**2*s35**3 + 1536*q9a*s11*s12*s14*s16*s22*s23*s26**2*s35**3 -  \
    512*q9a*s11*s12*s14*s15*s22*s24*s26**2*s35**3 - 1536*q9a*s11*s12*s13*s16*s22*s24*s26**2*s35**3 - 512*q9a*s11*s12**2*s16*s23*s24*s26**2*s35**3 + 512*q9a*s11*s12**2*s15*s24**2*s26**2*s35**3 -  \
    1024*q9a*s11*s12**2*s13*s22*s25*s26**2*s35**3 + 512*q9a*s11*s12*s14**2*s22*s25*s26**2*s35**3 + 1024*q9a*s11*s12**3*s23*s25*s26**2*s35**3 - 512*q9a*s11*s12**2*s14*s24*s25*s26**2*s35**3 - 512*q9a*s11*s12*s13*s14*s22*s26**3*s35**3 -  \
    512*q9a*s11*s12**2*s14*s23*s26**3*s35**3 + 1024*q9a*s11*s12**2*s13*s24*s26**3*s35**3 - 1024*q8a*s11**2*s12*s15*s22**3*s33*s35**3 + 1024*q8a*s11**2*s14*s16*s22**3*s33*s35**3 + 1024*q7a*s11**2*s15*s22**4*s33*s35**3 -  \
    512*q8a*s11**2*s12*s16*s22**2*s24*s33*s35**3 - 512*q7a*s11**2*s16*s22**3*s24*s33*s35**3 + 1024*q8a*s11**2*s12**2*s22**2*s25*s33*s35**3 - 1024*q7a*s11**2*s12*s22**3*s25*s33*s35**3 - 512*q8a*s11**2*s12*s14*s22**2*s26*s33*s35**3 +  \
    1024*q8a*s11*s12*s15*s16*s22**2*s26*s33*s35**3 - 1024*q8a*s11*s14*s16**2*s22**2*s26*s33*s35**3 - 512*q7a*s11**2*s14*s22**3*s26*s33*s35**3 - 1024*q7a*s11*s15*s16*s22**3*s26*s33*s35**3 +  \
    512*q8a*s11*s12*s16**2*s22*s24*s26*s33*s35**3 + 1024*q7a*s11**2*s12*s22**2*s24*s26*s33*s35**3 + 512*q7a*s11*s16**2*s22**2*s24*s26*s33*s35**3 - 1024*q8a*s11*s12**2*s16*s22*s25*s26*s33*s35**3 +  \
    1024*q7a*s11*s12*s16*s22**2*s25*s26*s33*s35**3 - 1024*q8a*s11*s12**2*s15*s22*s26**2*s33*s35**3 + 1536*q8a*s11*s12*s14*s16*s22*s26**2*s33*s35**3 + 1024*q7a*s11*s12*s15*s22**2*s26**2*s33*s35**3 +  \
    512*q7a*s11*s14*s16*s22**2*s26**2*s33*s35**3 - 512*q8a*s11*s12**2*s16*s24*s26**2*s33*s35**3 - 1536*q7a*s11*s12*s16*s22*s24*s26**2*s33*s35**3 + 1024*q8a*s11*s12**3*s25*s26**2*s33*s35**3 -  \
    1024*q7a*s11*s12**2*s22*s25*s26**2*s33*s35**3 - 512*q8a*s11*s12**2*s14*s26**3*s33*s35**3 - 512*q7a*s11*s12*s14*s22*s26**3*s33*s35**3 + 1024*q7a*s11*s12**2*s24*s26**3*s33*s35**3 + 1024*q8a*s11**2*s14*s15*s22**3*s34*s35**3 -  \
    1024*q8a*s11**2*s13*s16*s22**3*s34*s35**3 + 2048*q8a*s11**2*s12*s16*s22**2*s23*s34*s35**3 - 1024*q7a*s11**2*s16*s22**3*s23*s34*s35**3 - 512*q8a*s11**2*s12*s15*s22**2*s24*s34*s35**3 - 512*q8a*s11**2*s14*s16*s22**2*s24*s34*s35**3 -  \
    512*q7a*s11**2*s15*s22**3*s24*s34*s35**3 + 512*q7a*s11**2*s16*s22**2*s24**2*s34*s35**3 - 512*q8a*s11**2*s12*s14*s22**2*s25*s34*s35**3 - 512*q7a*s11**2*s14*s22**3*s25*s34*s35**3 + 1024*q7a*s11**2*s12*s22**2*s24*s25*s34*s35**3 -  \
    1024*q8a*s11**2*s12*s13*s22**2*s26*s34*s35**3 + 512*q8a*s11**2*s14**2*s22**2*s26*s34*s35**3 - 1024*q8a*s11*s14*s15*s16*s22**2*s26*s34*s35**3 + 1024*q8a*s11*s13*s16**2*s22**2*s26*s34*s35**3 +  \
    2048*q7a*s11**2*s13*s22**3*s26*s34*s35**3 - 2048*q8a*s11*s12*s16**2*s22*s23*s26*s34*s35**3 - 1024*q7a*s11**2*s12*s22**2*s23*s26*s34*s35**3 + 1024*q7a*s11*s16**2*s22**2*s23*s26*s34*s35**3 +  \
    512*q8a*s11*s12*s15*s16*s22*s24*s26*s34*s35**3 + 512*q8a*s11*s14*s16**2*s22*s24*s26*s34*s35**3 - 512*q7a*s11**2*s14*s22**2*s24*s26*s34*s35**3 + 512*q7a*s11*s15*s16*s22**2*s24*s26*s34*s35**3 -  \
    512*q7a*s11*s16**2*s22*s24**2*s26*s34*s35**3 + 512*q8a*s11*s12*s14*s16*s22*s25*s26*s34*s35**3 + 512*q7a*s11*s14*s16*s22**2*s25*s26*s34*s35**3 - 1024*q7a*s11*s12*s16*s22*s24*s25*s26*s34*s35**3 +  \
    1024*q8a*s11*s12*s14*s15*s22*s26**2*s34*s35**3 - 512*q8a*s11*s14**2*s16*s22*s26**2*s34*s35**3 - 2048*q7a*s11*s13*s16*s22**2*s26**2*s34*s35**3 + 2048*q8a*s11*s12**2*s16*s23*s26**2*s34*s35**3 -  \
    512*q8a*s11*s12**2*s15*s24*s26**2*s34*s35**3 - 512*q8a*s11*s12*s14*s16*s24*s26**2*s34*s35**3 - 512*q7a*s11*s12*s15*s22*s24*s26**2*s34*s35**3 + 512*q7a*s11*s14*s16*s22*s24*s26**2*s34*s35**3 +  \
    512*q7a*s11*s12*s16*s24**2*s26**2*s34*s35**3 - 512*q8a*s11*s12**2*s14*s25*s26**2*s34*s35**3 - 512*q7a*s11*s12*s14*s22*s25*s26**2*s34*s35**3 + 1024*q7a*s11*s12**2*s24*s25*s26**2*s34*s35**3 -  \
    1024*q8a*s11*s12**2*s13*s26**3*s34*s35**3 + 512*q8a*s11*s12*s14**2*s26**3*s34*s35**3 + 2048*q7a*s11*s12*s13*s22*s26**3*s34*s35**3 - 1024*q7a*s11*s12**2*s23*s26**3*s34*s35**3 - 512*q7a*s11*s12*s14*s24*s26**3*s34*s35**3 +  \
    1024*q8a*s11**2*s12*s13*s22**3*s35**4 - 512*q8a*s11**2*s14**2*s22**3*s35**4 - 1024*q7a*s11**2*s13*s22**4*s35**4 - 1024*q8a*s11**2*s12**2*s22**2*s23*s35**4 + 1024*q7a*s11**2*s12*s22**3*s23*s35**4 +  \
    512*q8a*s11**2*s12*s14*s22**2*s24*s35**4 + 512*q7a*s11**2*s14*s22**3*s24*s35**4 - 512*q7a*s11**2*s12*s22**2*s24**2*s35**4 - 1024*q8a*s11*s12*s13*s16*s22**2*s26*s35**4 + 512*q8a*s11*s14**2*s16*s22**2*s26*s35**4 +  \
    1024*q7a*s11*s13*s16*s22**3*s26*s35**4 + 1024*q8a*s11*s12**2*s16*s22*s23*s26*s35**4 - 1024*q7a*s11*s12*s16*s22**2*s23*s26*s35**4 - 512*q8a*s11*s12*s14*s16*s22*s24*s26*s35**4 - 512*q7a*s11*s14*s16*s22**2*s24*s26*s35**4 +  \
    512*q7a*s11*s12*s16*s22*s24**2*s26*s35**4 + 1024*q8a*s11*s12**2*s13*s22*s26**2*s35**4 - 512*q8a*s11*s12*s14**2*s22*s26**2*s35**4 - 1024*q7a*s11*s12*s13*s22**2*s26**2*s35**4 - 1024*q8a*s11*s12**3*s23*s26**2*s35**4 +  \
    1024*q7a*s11*s12**2*s22*s23*s26**2*s35**4 + 512*q8a*s11*s12**2*s14*s24*s26**2*s35**4 + 512*q7a*s11*s12*s14*s22*s24*s26**2*s35**4 - 512*q7a*s11*s12**2*s24**2*s26**2*s35**4 - 3072*q9a*s11**3*s16*s22**3*s23*s33**2*s36 -  \
    1536*q9a*s11**3*s15*s22**3*s24*s33**2*s36 + 1536*q9a*s11**3*s16*s22**2*s24**2*s33**2*s36 - 4608*q9a*s11**3*s14*s22**3*s25*s33**2*s36 + 4608*q9a*s11**2*s15*s16*s22**3*s25*s33**2*s36 + 6144*q9a*s11**3*s12*s22**2*s24*s25*s33**2*s36 -  \
    3072*q9a*s11**2*s16**2*s22**2*s24*s25*s33**2*s36 - 4608*q9a*s11**2*s12*s16*s22**2*s25**2*s33**2*s36 + 1536*q9a*s11*s16**3*s22**2*s25**2*s33**2*s36 - 3072*q9a*s11**3*s13*s22**3*s26*s33**2*s36 +  \
    1536*q9a*s11**2*s15**2*s22**3*s26*s33**2*s36 + 6144*q9a*s11**3*s12*s22**2*s23*s26*s33**2*s36 + 3072*q9a*s11**2*s16**2*s22**2*s23*s26*s33**2*s36 + 4608*q9a*s11**3*s14*s22**2*s24*s26*s33**2*s36 -  \
    1536*q9a*s11**2*s15*s16*s22**2*s24*s26*s33**2*s36 - 6144*q9a*s11**3*s12*s22*s24**2*s26*s33**2*s36 - 7680*q9a*s11**2*s12*s15*s22**2*s25*s26*s33**2*s36 + 1536*q9a*s11**2*s14*s16*s22**2*s25*s26*s33**2*s36 -  \
    1536*q9a*s11*s15*s16**2*s22**2*s25*s26*s33**2*s36 + 6144*q9a*s11**2*s12*s16*s22*s24*s25*s26*s33**2*s36 + 6144*q9a*s11**2*s12**2*s22*s25**2*s26*s33**2*s36 - 3072*q9a*s11*s12*s16**2*s22*s25**2*s26*s33**2*s36 -  \
    3072*q9a*s11**2*s14*s15*s22**2*s26**2*s33**2*s36 + 3072*q9a*s11**2*s13*s16*s22**2*s26**2*s33**2*s36 - 9216*q9a*s11**2*s12*s16*s22*s23*s26**2*s33**2*s36 + 4608*q9a*s11**2*s12*s15*s22*s24*s26**2*s33**2*s36 -  \
    1536*q9a*s11**2*s14*s16*s22*s24*s26**2*s33**2*s36 + 1536*q9a*s11**2*s12*s16*s24**2*s26**2*s33**2*s36 + 1536*q9a*s11**2*s12*s14*s22*s25*s26**2*s33**2*s36 + 3072*q9a*s11*s12*s15*s16*s22*s25*s26**2*s33**2*s36 -  \
    6144*q9a*s11**2*s12**2*s24*s25*s26**2*s33**2*s36 + 1536*q9a*s11*s12**2*s16*s25**2*s26**2*s33**2*s36 - 3072*q9a*s11**2*s12*s13*s22*s26**3*s33**2*s36 + 1536*q9a*s11**2*s14**2*s22*s26**3*s33**2*s36 +  \
    6144*q9a*s11**2*s12**2*s23*s26**3*s33**2*s36 - 1536*q9a*s11**2*s12*s14*s24*s26**3*s33**2*s36 - 1536*q9a*s11*s12**2*s15*s25*s26**3*s33**2*s36 - 1024*q8a*s11**3*s16*s22**3*s33**3*s36 + 2048*q8a*s11**3*s12*s22**2*s26*s33**3*s36 +  \
    1024*q8a*s11**2*s16**2*s22**2*s26*s33**3*s36 - 1024*q7a*s11**3*s22**3*s26*s33**3*s36 - 3072*q8a*s11**2*s12*s16*s22*s26**2*s33**3*s36 + 1024*q7a*s11**2*s16*s22**2*s26**2*s33**3*s36 + 2048*q8a*s11**2*s12**2*s26**3*s33**3*s36 -  \
    1024*q7a*s11**2*s12*s22*s26**3*s33**3*s36 + 2048*q9a*s11**3*s15*s22**3*s23*s33*s34*s36 + 1024*q9a*s11**3*s16*s22**2*s23*s24*s33*s34*s36 + 1024*q9a*s11**3*s15*s22**2*s24**2*s33*s34*s36 -  \
    1024*q9a*s11**3*s16*s22*s24**3*s33*s34*s36 + 4096*q9a*s11**3*s13*s22**3*s25*s33*s34*s36 - 2048*q9a*s11**2*s15**2*s22**3*s25*s33*s34*s36 - 6144*q9a*s11**3*s12*s22**2*s23*s25*s33*s34*s36 +  \
    1024*q9a*s11**2*s16**2*s22**2*s23*s25*s33*s34*s36 + 1024*q9a*s11**3*s14*s22**2*s24*s25*s33*s34*s36 - 2048*q9a*s11**2*s15*s16*s22**2*s24*s25*s33*s34*s36 - 2048*q9a*s11**3*s12*s22*s24**2*s25*s33*s34*s36 +  \
    2048*q9a*s11**2*s16**2*s22*s24**2*s25*s33*s34*s36 + 4096*q9a*s11**2*s12*s15*s22**2*s25**2*s33*s34*s36 + 1024*q9a*s11**2*s14*s16*s22**2*s25**2*s33*s34*s36 - 1024*q9a*s11*s15*s16**2*s22**2*s25**2*s33*s34*s36 +  \
    1024*q9a*s11**2*s12*s16*s22*s24*s25**2*s33*s34*s36 - 1024*q9a*s11*s16**3*s22*s24*s25**2*s33*s34*s36 - 2048*q9a*s11**2*s12**2*s22*s25**3*s33*s34*s36 + 1024*q9a*s11*s12*s16**2*s22*s25**3*s33*s34*s36 -  \
    5120*q9a*s11**3*s14*s22**2*s23*s26*s33*s34*s36 - 1024*q9a*s11**2*s15*s16*s22**2*s23*s26*s33*s34*s36 - 1024*q9a*s11**2*s15**2*s22**2*s24*s26*s33*s34*s36 + 4096*q9a*s11**3*s12*s22*s23*s24*s26*s33*s34*s36 -  \
    2048*q9a*s11**2*s16**2*s22*s23*s24*s26*s33*s34*s36 - 1024*q9a*s11**3*s14*s22*s24**2*s26*s33*s34*s36 + 1024*q9a*s11**2*s15*s16*s22*s24**2*s26*s33*s34*s36 + 2048*q9a*s11**3*s12*s24**3*s26*s33*s34*s36 +  \
    4096*q9a*s11**2*s14*s15*s22**2*s25*s26*s33*s34*s36 - 3072*q9a*s11**2*s13*s16*s22**2*s25*s26*s33*s34*s36 + 1024*q9a*s11*s15**2*s16*s22**2*s25*s26*s33*s34*s36 + 2048*q9a*s11**2*s12*s16*s22*s23*s25*s26*s33*s34*s36 -  \
    2048*q9a*s11**2*s14*s16*s22*s24*s25*s26*s33*s34*s36 + 1024*q9a*s11*s15*s16**2*s22*s24*s25*s26*s33*s34*s36 - 3072*q9a*s11**2*s12*s16*s24**2*s25*s26*s33*s34*s36 - 5120*q9a*s11**2*s12*s14*s22*s25**2*s26*s33*s34*s36 +  \
    1024*q9a*s11*s14*s16**2*s22*s25**2*s26*s33*s34*s36 + 2048*q9a*s11**2*s12**2*s24*s25**2*s26*s33*s34*s36 + 1024*q9a*s11*s12*s16**2*s24*s25**2*s26*s33*s34*s36 - 1024*q9a*s11*s12**2*s16*s25**3*s26*s33*s34*s36 +  \
    1024*q9a*s11**2*s13*s15*s22**2*s26**2*s33*s34*s36 + 4096*q9a*s11**2*s14*s16*s22*s23*s26**2*s33*s34*s36 + 1024*q9a*s11**2*s14*s15*s22*s24*s26**2*s33*s34*s36 - 1024*q9a*s11**2*s13*s16*s22*s24*s26**2*s33*s34*s36 +  \
    1024*q9a*s11**2*s12*s16*s23*s24*s26**2*s33*s34*s36 - 2048*q9a*s11**2*s12*s15*s24**2*s26**2*s33*s34*s36 + 2048*q9a*s11**2*s12*s13*s22*s25*s26**2*s33*s34*s36 - 2048*q9a*s11**2*s14**2*s22*s25*s26**2*s33*s34*s36 -  \
    1024*q9a*s11*s12*s15**2*s22*s25*s26**2*s33*s34*s36 - 1024*q9a*s11*s14*s15*s16*s22*s25*s26**2*s33*s34*s36 - 2048*q9a*s11**2*s12**2*s23*s25*s26**2*s33*s34*s36 + 5120*q9a*s11**2*s12*s14*s24*s25*s26**2*s33*s34*s36 -  \
    1024*q9a*s11*s12*s15*s16*s24*s25*s26**2*s33*s34*s36 + 1024*q9a*s11*s12**2*s15*s25**2*s26**2*s33*s34*s36 - 1024*q9a*s11*s12*s14*s16*s25**2*s26**2*s33*s34*s36 - 1024*q9a*s11**2*s13*s14*s22*s26**3*s33*s34*s36 -  \
    3072*q9a*s11**2*s12*s14*s23*s26**3*s33*s34*s36 + 2048*q9a*s11**2*s12*s13*s24*s26**3*s33*s34*s36 + 1024*q9a*s11*s12*s14*s15*s25*s26**3*s33*s34*s36 + 1024*q8a*s11**3*s15*s22**3*s33**2*s34*s36 +  \
    512*q8a*s11**3*s16*s22**2*s24*s33**2*s34*s36 - 3072*q8a*s11**3*s12*s22**2*s25*s33**2*s34*s36 + 512*q8a*s11**2*s16**2*s22**2*s25*s33**2*s34*s36 + 2048*q7a*s11**3*s22**3*s25*s33**2*s34*s36 -  \
    2560*q8a*s11**3*s14*s22**2*s26*s33**2*s34*s36 - 512*q8a*s11**2*s15*s16*s22**2*s26*s33**2*s34*s36 + 2048*q8a*s11**3*s12*s22*s24*s26*s33**2*s34*s36 - 1024*q8a*s11**2*s16**2*s22*s24*s26*s33**2*s34*s36 +  \
    1024*q8a*s11**2*s12*s16*s22*s25*s26*s33**2*s34*s36 - 1536*q7a*s11**2*s16*s22**2*s25*s26*s33**2*s34*s36 + 2048*q8a*s11**2*s14*s16*s22*s26**2*s33**2*s34*s36 + 512*q7a*s11**2*s15*s22**2*s26**2*s33**2*s34*s36
v3_3= \
    512*q8a*s11**2*s12*s16*s24*s26**2*s33**2*s34*s36 - 512*q7a*s11**2*s16*s22*s24*s26**2*s33**2*s34*s36 - 1024*q8a*s11**2*s12**2*s25*s26**2*s33**2*s34*s36 + 1024*q7a*s11**2*s12*s22*s25*s26**2*s33**2*s34*s36 -  \
    1536*q8a*s11**2*s12*s14*s26**3*s33**2*s34*s36 - 512*q7a*s11**2*s14*s22*s26**3*s33**2*s34*s36 + 1024*q7a*s11**2*s12*s24*s26**3*s33**2*s34*s36 - 1024*q9a*s11**3*s16*s22**2*s23**2*s34**2*s36 -  \
    1536*q9a*s11**3*s15*s22**2*s23*s24*s34**2*s36 + 1024*q9a*s11**3*s16*s22*s23*s24**2*s34**2*s36 + 512*q9a*s11**3*s14*s22**2*s23*s25*s34**2*s36 + 1536*q9a*s11**2*s15*s16*s22**2*s23*s25*s34**2*s36 -  \
    1024*q9a*s11**3*s13*s22**2*s24*s25*s34**2*s36 + 1024*q9a*s11**2*s15**2*s22**2*s24*s25*s34**2*s36 + 2048*q9a*s11**3*s12*s22*s23*s24*s25*s34**2*s36 - 1536*q9a*s11**2*s16**2*s22*s23*s24*s25*s34**2*s36 -  \
    512*q9a*s11**2*s15*s16*s22*s24**2*s25*s34**2*s36 - 1024*q9a*s11**2*s14*s15*s22**2*s25**2*s34**2*s36 - 512*q9a*s11**2*s13*s16*s22**2*s25**2*s34**2*s36 - 1024*q9a*s11**2*s12*s16*s22*s23*s25**2*s34**2*s36 +  \
    512*q9a*s11*s16**3*s22*s23*s25**2*s34**2*s36 - 1024*q9a*s11**2*s12*s15*s22*s24*s25**2*s34**2*s36 + 512*q9a*s11**2*s14*s16*s22*s24*s25**2*s34**2*s36 + 512*q9a*s11*s15*s16**2*s22*s24*s25**2*s34**2*s36 +  \
    1024*q9a*s11**2*s12*s14*s22*s25**3*s34**2*s36 - 512*q9a*s11*s14*s16**2*s22*s25**3*s34**2*s36 + 1024*q9a*s11**3*s13*s22**2*s23*s26*s34**2*s36 + 512*q9a*s11**2*s15**2*s22**2*s23*s26*s34**2*s36 +  \
    1024*q9a*s11**2*s16**2*s22*s23**2*s26*s34**2*s36 + 1024*q9a*s11**3*s14*s22*s23*s24*s26*s34**2*s36 - 2048*q9a*s11**3*s12*s23*s24**2*s26*s34**2*s36 - 1536*q9a*s11**2*s13*s15*s22**2*s25*s26*s34**2*s36 -  \
    1024*q9a*s11**2*s12*s15*s22*s23*s25*s26*s34**2*s36 - 1024*q9a*s11**2*s14*s16*s22*s23*s25*s26*s34**2*s36 - 512*q9a*s11*s15*s16**2*s22*s23*s25*s26*s34**2*s36 - 512*q9a*s11**2*s14*s15*s22*s24*s25*s26*s34**2*s36 +  \
    2048*q9a*s11**2*s13*s16*s22*s24*s25*s26*s34**2*s36 - 512*q9a*s11*s15**2*s16*s22*s24*s25*s26*s34**2*s36 + 2048*q9a*s11**2*s12*s16*s23*s24*s25*s26*s34**2*s36 + 1024*q9a*s11**2*s12*s15*s24**2*s25*s26*s34**2*s36 +  \
    2048*q9a*s11**2*s12*s13*s22*s25**2*s26*s34**2*s36 + 512*q9a*s11**2*s14**2*s22*s25**2*s26*s34**2*s36 + 512*q9a*s11*s14*s15*s16*s22*s25**2*s26*s34**2*s36 - 512*q9a*s11*s13*s16**2*s22*s25**2*s26*s34**2*s36 -  \
    512*q9a*s11*s12*s16**2*s23*s25**2*s26*s34**2*s36 - 1024*q9a*s11**2*s12*s14*s24*s25**2*s26*s34**2*s36 - 512*q9a*s11*s12*s15*s16*s24*s25**2*s26*s34**2*s36 + 512*q9a*s11*s12*s14*s16*s25**3*s26*s34**2*s36 -  \
    512*q9a*s11**2*s14*s15*s22*s23*s26**2*s34**2*s36 - 1024*q9a*s11**2*s13*s16*s22*s23*s26**2*s34**2*s36 - 1024*q9a*s11**2*s12*s16*s23**2*s26**2*s34**2*s36 + 512*q9a*s11**2*s12*s15*s23*s24*s26**2*s34**2*s36 +  \
    1024*q9a*s11**2*s13*s14*s22*s25*s26**2*s34**2*s36 + 512*q9a*s11*s13*s15*s16*s22*s25*s26**2*s34**2*s36 + 512*q9a*s11**2*s12*s14*s23*s25*s26**2*s34**2*s36 + 512*q9a*s11*s12*s15*s16*s23*s25*s26**2*s34**2*s36 -  \
    3072*q9a*s11**2*s12*s13*s24*s25*s26**2*s34**2*s36 + 512*q9a*s11*s12*s15**2*s24*s25*s26**2*s34**2*s36 - 512*q9a*s11*s12*s14*s15*s25**2*s26**2*s34**2*s36 + 512*q9a*s11*s12*s13*s16*s25**2*s26**2*s34**2*s36 +  \
    1024*q9a*s11**2*s12*s13*s23*s26**3*s34**2*s36 - 512*q9a*s11*s12*s13*s15*s25*s26**3*s34**2*s36 - 2048*q8a*s11**3*s16*s22**2*s23*s33*s34**2*s36 - 1536*q8a*s11**3*s15*s22**2*s24*s33*s34**2*s36 +  \
    1024*q8a*s11**3*s16*s22*s24**2*s33*s34**2*s36 + 512*q8a*s11**3*s14*s22**2*s25*s33*s34**2*s36 + 1536*q8a*s11**2*s15*s16*s22**2*s25*s33*s34**2*s36 + 2048*q8a*s11**3*s12*s22*s24*s25*s33*s34**2*s36 -  \
    1536*q8a*s11**2*s16**2*s22*s24*s25*s33*s34**2*s36 - 1024*q7a*s11**3*s22**2*s24*s25*s33*s34**2*s36 - 1024*q8a*s11**2*s12*s16*s22*s25**2*s33*s34**2*s36 + 512*q8a*s11*s16**3*s22*s25**2*s33*s34**2*s36 -  \
    512*q7a*s11**2*s16*s22**2*s25**2*s33*s34**2*s36 + 1024*q8a*s11**3*s13*s22**2*s26*s33*s34**2*s36 + 512*q8a*s11**2*s15**2*s22**2*s26*s33*s34**2*s36 + 2048*q8a*s11**2*s16**2*s22*s23*s26*s33*s34**2*s36 +  \
    1024*q7a*s11**3*s22**2*s23*s26*s33*s34**2*s36 + 1024*q8a*s11**3*s14*s22*s24*s26*s33*s34**2*s36 - 2048*q8a*s11**3*s12*s24**2*s26*s33*s34**2*s36 - 1024*q8a*s11**2*s12*s15*s22*s25*s26*s33*s34**2*s36 -  \
    1024*q8a*s11**2*s14*s16*s22*s25*s26*s33*s34**2*s36 - 512*q8a*s11*s15*s16**2*s22*s25*s26*s33*s34**2*s36 - 1536*q7a*s11**2*s15*s22**2*s25*s26*s33*s34**2*s36 + 2048*q8a*s11**2*s12*s16*s24*s25*s26*s33*s34**2*s36 +  \
    2048*q7a*s11**2*s16*s22*s24*s25*s26*s33*s34**2*s36 - 512*q8a*s11*s12*s16**2*s25**2*s26*s33*s34**2*s36 + 2048*q7a*s11**2*s12*s22*s25**2*s26*s33*s34**2*s36 - 512*q7a*s11*s16**2*s22*s25**2*s26*s33*s34**2*s36 -  \
    512*q8a*s11**2*s14*s15*s22*s26**2*s33*s34**2*s36 - 1024*q8a*s11**2*s13*s16*s22*s26**2*s33*s34**2*s36 - 2048*q8a*s11**2*s12*s16*s23*s26**2*s33*s34**2*s36 - 1024*q7a*s11**2*s16*s22*s23*s26**2*s33*s34**2*s36 +  \
    512*q8a*s11**2*s12*s15*s24*s26**2*s33*s34**2*s36 + 512*q8a*s11**2*s12*s14*s25*s26**2*s33*s34**2*s36 + 512*q8a*s11*s12*s15*s16*s25*s26**2*s33*s34**2*s36 + 1024*q7a*s11**2*s14*s22*s25*s26**2*s33*s34**2*s36 +  \
    512*q7a*s11*s15*s16*s22*s25*s26**2*s33*s34**2*s36 - 3072*q7a*s11**2*s12*s24*s25*s26**2*s33*s34**2*s36 + 512*q7a*s11*s12*s16*s25**2*s26**2*s33*s34**2*s36 + 1024*q8a*s11**2*s12*s13*s26**3*s33*s34**2*s36 +  \
    1024*q7a*s11**2*s12*s23*s26**3*s33*s34**2*s36 - 512*q7a*s11*s12*s15*s25*s26**3*s33*s34**2*s36 + 2048*q8a*s11**3*s15*s22**2*s23*s34**3*s36 - 1024*q8a*s11**3*s16*s22*s23*s24*s34**3*s36 -  \
    1024*q8a*s11**2*s15**2*s22**2*s25*s34**3*s36 - 2048*q8a*s11**3*s12*s22*s23*s25*s34**3*s36 + 1024*q8a*s11**2*s16**2*s22*s23*s25*s34**3*s36 + 512*q8a*s11**2*s15*s16*s22*s24*s25*s34**3*s36 +  \
    1024*q8a*s11**2*s12*s15*s22*s25**2*s34**3*s36 - 512*q8a*s11*s15*s16**2*s22*s25**2*s34**3*s36 + 1024*q7a*s11**2*s15*s22**2*s25**2*s34**3*s36 - 512*q7a*s11**2*s16*s22*s24*s25**2*s34**3*s36 - 1024*q7a*s11**2*s12*s22*s25**3*s34**3*s36 +  \
    512*q7a*s11*s16**2*s22*s25**3*s34**3*s36 - 1024*q8a*s11**3*s14*s22*s23*s26*s34**3*s36 - 1024*q8a*s11**2*s15*s16*s22*s23*s26*s34**3*s36 + 2048*q8a*s11**3*s12*s23*s24*s26*s34**3*s36 +  \
    512*q8a*s11**2*s14*s15*s22*s25*s26*s34**3*s36 + 512*q8a*s11*s15**2*s16*s22*s25*s26*s34**3*s36 - 1024*q8a*s11**2*s12*s16*s23*s25*s26*s34**3*s36 - 1024*q8a*s11**2*s12*s15*s24*s25*s26*s34**3*s36 +  \
    512*q8a*s11*s12*s15*s16*s25**2*s26*s34**3*s36 - 512*q7a*s11**2*s14*s22*s25**2*s26*s34**3*s36 - 512*q7a*s11*s15*s16*s22*s25**2*s26*s34**3*s36 + 1024*q7a*s11**2*s12*s24*s25**2*s26*s34**3*s36 -  \
    512*q7a*s11*s12*s16*s25**3*s26*s34**3*s36 + 1024*q8a*s11**2*s12*s15*s23*s26**2*s34**3*s36 - 512*q8a*s11*s12*s15**2*s25*s26**2*s34**3*s36 + 512*q7a*s11*s12*s15*s25**2*s26**2*s34**3*s36 +  \
    4096*q9a*s11**3*s14*s22**3*s23*s33*s35*s36 + 2048*q9a*s11**3*s13*s22**3*s24*s33*s35*s36 + 1024*q9a*s11**2*s15**2*s22**3*s24*s33*s35*s36 - 6144*q9a*s11**3*s12*s22**2*s23*s24*s33*s35*s36 +  \
    1024*q9a*s11**2*s16**2*s22**2*s23*s24*s33*s35*s36 - 2048*q9a*s11**3*s14*s22**2*s24**2*s33*s35*s36 - 1024*q9a*s11**2*s15*s16*s22**2*s24**2*s33*s35*s36 + 2048*q9a*s11**3*s12*s22*s24**3*s33*s35*s36 +  \
    1024*q9a*s11**2*s14*s15*s22**3*s25*s33*s35*s36 - 5120*q9a*s11**2*s13*s16*s22**3*s25*s33*s35*s36 - 1024*q9a*s11*s15**2*s16*s22**3*s25*s33*s35*s36 + 5120*q9a*s11**2*s12*s16*s22**2*s23*s25*s33*s35*s36 -  \
    1024*q9a*s11*s16**3*s22**2*s23*s25*s33*s35*s36 - 3072*q9a*s11**2*s12*s15*s22**2*s24*s25*s33*s35*s36 + 4096*q9a*s11**2*s14*s16*s22**2*s24*s25*s33*s35*s36 + 1024*q9a*s11*s15*s16**2*s22**2*s24*s25*s33*s35*s36 -  \
    3072*q9a*s11**2*s12*s16*s22*s24**2*s25*s33*s35*s36 - 1024*q9a*s11**2*s12*s14*s22**2*s25**2*s33*s35*s36 + 2048*q9a*s11*s12*s15*s16*s22**2*s25**2*s33*s35*s36 - 2048*q9a*s11*s14*s16**2*s22**2*s25**2*s33*s35*s36 +  \
    2048*q9a*s11**2*s12**2*s22*s24*s25**2*s33*s35*s36 + 1024*q9a*s11*s12*s16**2*s22*s24*s25**2*s33*s35*s36 - 1024*q9a*s11*s12**2*s16*s22*s25**3*s33*s35*s36 + 1024*q9a*s11**2*s13*s15*s22**3*s26*s33*s35*s36 -  \
    1024*q9a*s11*s15**3*s22**3*s26*s33*s35*s36 - 1024*q9a*s11**2*s12*s15*s22**2*s23*s26*s33*s35*s36 - 3072*q9a*s11**2*s14*s16*s22**2*s23*s26*s33*s35*s36 - 1024*q9a*s11*s15*s16**2*s22**2*s23*s26*s33*s35*s36 -  \
    2048*q9a*s11**2*s14*s15*s22**2*s24*s26*s33*s35*s36 - 1024*q9a*s11**2*s13*s16*s22**2*s24*s26*s33*s35*s36 + 1024*q9a*s11*s15**2*s16*s22**2*s24*s26*s33*s35*s36 + 2048*q9a*s11**2*s12*s16*s22*s23*s24*s26*s33*s35*s36 +  \
    3072*q9a*s11**2*s12*s15*s22*s24**2*s26*s33*s35*s36 + 1024*q9a*s11**2*s14*s16*s22*s24**2*s26*s33*s35*s36 - 1024*q9a*s11**2*s12*s16*s24**3*s26*s33*s35*s36 + 4096*q9a*s11**2*s12*s13*s22**2*s25*s26*s33*s35*s36 +  \
    1024*q9a*s11**2*s14**2*s22**2*s25*s26*s33*s35*s36 + 4096*q9a*s11*s12*s15**2*s22**2*s25*s26*s33*s35*s36 - 2048*q9a*s11*s14*s15*s16*s22**2*s25*s26*s33*s35*s36 + 4096*q9a*s11*s13*s16**2*s22**2*s25*s26*s33*s35*s36 -  \
    4096*q9a*s11**2*s12**2*s22*s23*s25*s26*s33*s35*s36 - 4096*q9a*s11**2*s12*s14*s22*s24*s25*s26*s33*s35*s36 - 2048*q9a*s11*s12*s15*s16*s22*s24*s25*s26*s33*s35*s36 - 1024*q9a*s11*s14*s16**2*s22*s24*s25*s26*s33*s35*s36 +  \
    2048*q9a*s11**2*s12**2*s24**2*s25*s26*s33*s35*s36 + 1024*q9a*s11*s12*s16**2*s24**2*s25*s26*s33*s35*s36 - 5120*q9a*s11*s12**2*s15*s22*s25**2*s26*s33*s35*s36 + 6144*q9a*s11*s12*s14*s16*s22*s25**2*s26*s33*s35*s36 -  \
    3072*q9a*s11*s12**2*s16*s24*s25**2*s26*s33*s35*s36 + 2048*q9a*s11*s12**3*s25**3*s26*s33*s35*s36 + 1024*q9a*s11**2*s13*s14*s22**2*s26**2*s33*s35*s36 + 2048*q9a*s11*s14*s15**2*s22**2*s26**2*s33*s35*s36 -  \
    2048*q9a*s11*s13*s15*s16*s22**2*s26**2*s33*s35*s36 + 2048*q9a*s11**2*s12*s14*s22*s23*s26**2*s33*s35*s36 + 4096*q9a*s11*s12*s15*s16*s22*s23*s26**2*s33*s35*s36 - 1024*q9a*s11**2*s14**2*s22*s24*s26**2*s33*s35*s36 -  \
    3072*q9a*s11*s12*s15**2*s22*s24*s26**2*s33*s35*s36 + 1024*q9a*s11*s14*s15*s16*s22*s24*s26**2*s33*s35*s36 - 2048*q9a*s11**2*s12**2*s23*s24*s26**2*s33*s35*s36 + 1024*q9a*s11**2*s12*s14*s24**2*s26**2*s33*s35*s36 -  \
    1024*q9a*s11*s12*s15*s16*s24**2*s26**2*s33*s35*s36 - 2048*q9a*s11*s12*s14*s15*s22*s25*s26**2*s33*s35*s36 - 6144*q9a*s11*s12*s13*s16*s22*s25*s26**2*s33*s35*s36 + 1024*q9a*s11*s14**2*s16*s22*s25*s26**2*s33*s35*s36 +  \
    1024*q9a*s11*s12**2*s16*s23*s25*s26**2*s33*s35*s36 + 5120*q9a*s11*s12**2*s15*s24*s25*s26**2*s33*s35*s36 - 1024*q9a*s11*s12*s14*s16*s24*s25*s26**2*s33*s35*s36 - 2048*q9a*s11*s12**2*s14*s25**2*s26**2*s33*s35*s36 +  \
    2048*q9a*s11*s12*s13*s15*s22*s26**3*s33*s35*s36 - 1024*q9a*s11*s14**2*s15*s22*s26**3*s33*s35*s36 - 3072*q9a*s11*s12**2*s15*s23*s26**3*s33*s35*s36 + 1024*q9a*s11*s12*s14*s15*s24*s26**3*s33*s35*s36 +  \
    2048*q9a*s11*s12**2*s13*s25*s26**3*s33*s35*s36 + 2048*q8a*s11**3*s14*s22**3*s33**2*s35*s36 - 3072*q8a*s11**3*s12*s22**2*s24*s33**2*s35*s36 + 512*q8a*s11**2*s16**2*s22**2*s24*s33**2*s35*s36 +  \
    1024*q7a*s11**3*s22**3*s24*s33**2*s35*s36 + 2560*q8a*s11**2*s12*s16*s22**2*s25*s33**2*s35*s36 - 512*q8a*s11*s16**3*s22**2*s25*s33**2*s35*s36 - 2560*q7a*s11**2*s16*s22**3*s25*s33**2*s35*s36 -  \
    512*q8a*s11**2*s12*s15*s22**2*s26*s33**2*s35*s36 - 1536*q8a*s11**2*s14*s16*s22**2*s26*s33**2*s35*s36 - 512*q8a*s11*s15*s16**2*s22**2*s26*s33**2*s35*s36 + 512*q7a*s11**2*s15*s22**3*s26*s33**2*s35*s36 +  \
    1024*q8a*s11**2*s12*s16*s22*s24*s26*s33**2*s35*s36 - 512*q7a*s11**2*s16*s22**2*s24*s26*s33**2*s35*s36 - 2048*q8a*s11**2*s12**2*s22*s25*s26*s33**2*s35*s36 + 2048*q7a*s11**2*s12*s22**2*s25*s26*s33**2*s35*s36 +  \
    2048*q7a*s11*s16**2*s22**2*s25*s26*s33**2*s35*s36 + 1024*q8a*s11**2*s12*s14*s22*s26**2*s33**2*s35*s36 + 2048*q8a*s11*s12*s15*s16*s22*s26**2*s33**2*s35*s36 + 512*q7a*s11**2*s14*s22**2*s26**2*s33**2*s35*s36 -  \
    1024*q7a*s11*s15*s16*s22**2*s26**2*s33**2*s35*s36 - 1024*q8a*s11**2*s12**2*s24*s26**2*s33**2*s35*s36 + 512*q8a*s11*s12**2*s16*s25*s26**2*s33**2*s35*s36 - 3072*q7a*s11*s12*s16*s22*s25*s26**2*s33**2*s35*s36 -  \
    1536*q8a*s11*s12**2*s15*s26**3*s33**2*s35*s36 + 1024*q7a*s11*s12*s15*s22*s26**3*s33**2*s35*s36 + 1024*q7a*s11*s12**2*s25*s26**3*s33**2*s35*s36 - 4096*q9a*s11**3*s13*s22**3*s23*s34*s35*s36 +  \
    4096*q9a*s11**3*s12*s22**2*s23**2*s34*s35*s36 + 2048*q9a*s11**3*s14*s22**2*s23*s24*s34*s35*s36 - 512*q9a*s11**2*s15**2*s22**2*s24**2*s34*s35*s36 - 2048*q9a*s11**3*s12*s22*s23*s24**2*s34*s35*s36 -  \
    512*q9a*s11**2*s16**2*s22*s23*s24**2*s34*s35*s36 + 512*q9a*s11**2*s15*s16*s22*s24**3*s34*s35*s36 + 2048*q9a*s11**2*s13*s15*s22**3*s25*s34*s35*s36 - 2048*q9a*s11**2*s12*s15*s22**2*s23*s25*s34*s35*s36 -  \
    2048*q9a*s11**2*s14*s16*s22**2*s23*s25*s34*s35*s36 - 1024*q9a*s11**2*s14*s15*s22**2*s24*s25*s34*s35*s36 + 512*q9a*s11*s15**2*s16*s22**2*s24*s25*s34*s35*s36 + 2048*q9a*s11**2*s12*s16*s22*s23*s24*s25*s34*s35*s36
v3_4= \
    512*q9a*s11*s16**3*s22*s23*s24*s25*s34*s35*s36 + 2048*q9a*s11**2*s12*s15*s22*s24**2*s25*s34*s35*s36 - 512*q9a*s11**2*s14*s16*s22*s24**2*s25*s34*s35*s36 - 512*q9a*s11*s15*s16**2*s22*s24**2*s25*s34*s35*s36 -  \
    2048*q9a*s11**2*s12*s13*s22**2*s25**2*s34*s35*s36 + 1536*q9a*s11**2*s14**2*s22**2*s25**2*s34*s35*s36 - 512*q9a*s11*s14*s15*s16*s22**2*s25**2*s34*s35*s36 + 1536*q9a*s11*s13*s16**2*s22**2*s25**2*s34*s35*s36 +  \
    2048*q9a*s11**2*s12**2*s22*s23*s25**2*s34*s35*s36 - 1536*q9a*s11*s12*s16**2*s22*s23*s25**2*s34*s35*s36 - 2048*q9a*s11**2*s12*s14*s22*s24*s25**2*s34*s35*s36 - 512*q9a*s11*s12*s15*s16*s22*s24*s25**2*s34*s35*s36 +  \
    512*q9a*s11*s14*s16**2*s22*s24*s25**2*s34*s35*s36 + 512*q9a*s11*s12*s14*s16*s22*s25**3*s34*s35*s36 + 4096*q9a*s11**2*s13*s16*s22**2*s23*s26*s34*s35*s36 - 4096*q9a*s11**2*s12*s16*s22*s23**2*s26*s34*s35*s36 +  \
    512*q9a*s11*s15**3*s22**2*s24*s26*s34*s35*s36 - 1024*q9a*s11**2*s14*s16*s22*s23*s24*s26*s34*s35*s36 + 512*q9a*s11*s15*s16**2*s22*s23*s24*s26*s34*s35*s36 + 512*q9a*s11**2*s14*s15*s22*s24**2*s26*s34*s35*s36 -  \
    512*q9a*s11*s15**2*s16*s22*s24**2*s26*s34*s35*s36 + 2048*q9a*s11**2*s12*s16*s23*s24**2*s26*s34*s35*s36 - 1024*q9a*s11**2*s12*s15*s24**3*s26*s34*s35*s36 - 2048*q9a*s11**2*s13*s14*s22**2*s25*s26*s34*s35*s36 -  \
    512*q9a*s11*s14*s15**2*s22**2*s25*s26*s34*s35*s36 - 1024*q9a*s11*s13*s15*s16*s22**2*s25*s26*s34*s35*s36 + 4096*q9a*s11**2*s12*s14*s22*s23*s25*s26*s34*s35*s36 + 1024*q9a*s11*s12*s15*s16*s22*s23*s25*s26*s34*s35*s36 +  \
    512*q9a*s11*s14*s16**2*s22*s23*s25*s26*s34*s35*s36 + 2048*q9a*s11**2*s12*s13*s22*s24*s25*s26*s34*s35*s36 - 512*q9a*s11**2*s14**2*s22*s24*s25*s26*s34*s35*s36 - 1536*q9a*s11*s12*s15**2*s22*s24*s25*s26*s34*s35*s36 +  \
    2048*q9a*s11*s14*s15*s16*s22*s24*s25*s26*s34*s35*s36 - 1536*q9a*s11*s13*s16**2*s22*s24*s25*s26*s34*s35*s36 - 4096*q9a*s11**2*s12**2*s23*s24*s25*s26*s34*s35*s36 - 1024*q9a*s11*s12*s16**2*s23*s24*s25*s26*s34*s35*s36 +  \
    1024*q9a*s11**2*s12*s14*s24**2*s25*s26*s34*s35*s36 + 1536*q9a*s11*s12*s14*s15*s22*s25**2*s26*s34*s35*s36 - 2048*q9a*s11*s12*s13*s16*s22*s25**2*s26*s34*s35*s36 - 1536*q9a*s11*s14**2*s16*s22*s25**2*s26*s34*s35*s36 +  \
    2048*q9a*s11*s12**2*s16*s23*s25**2*s26*s34*s35*s36 + 1024*q9a*s11*s12**2*s15*s24*s25**2*s26*s34*s35*s36 - 1024*q9a*s11*s12**2*s14*s25**3*s26*s34*s35*s36 - 512*q9a*s11*s13*s15**2*s22**2*s26**2*s34*s35*s36 -  \
    4096*q9a*s11**2*s12*s13*s22*s23*s26**2*s34*s35*s36 + 1536*q9a*s11**2*s14**2*s22*s23*s26**2*s34*s35*s36 + 512*q9a*s11*s12*s15**2*s22*s23*s26**2*s34*s35*s36 - 1536*q9a*s11*s14*s15*s16*s22*s23*s26**2*s34*s35*s36 +  \
    4096*q9a*s11**2*s12**2*s23**2*s26**2*s34*s35*s36 - 512*q9a*s11*s14*s15**2*s22*s24*s26**2*s34*s35*s36 + 512*q9a*s11*s13*s15*s16*s22*s24*s26**2*s34*s35*s36 - 2048*q9a*s11**2*s12*s14*s23*s24*s26**2*s34*s35*s36 +  \
    1024*q9a*s11*s12*s15**2*s24**2*s26**2*s34*s35*s36 + 2048*q9a*s11*s12*s13*s15*s22*s25*s26**2*s34*s35*s36 + 512*q9a*s11*s14**2*s15*s22*s25*s26**2*s34*s35*s36 + 512*q9a*s11*s13*s14*s16*s22*s25*s26**2*s34*s35*s36 -  \
    2048*q9a*s11*s12**2*s15*s23*s25*s26**2*s34*s35*s36 - 2048*q9a*s11*s12*s14*s15*s24*s25*s26**2*s34*s35*s36 + 2048*q9a*s11*s12*s13*s16*s24*s25*s26**2*s34*s35*s36 + 1024*q9a*s11*s12*s14**2*s25**2*s26**2*s34*s35*s36 +  \
    512*q9a*s11*s13*s14*s15*s22*s26**3*s34*s35*s36 + 1024*q9a*s11*s12*s14*s15*s23*s26**3*s34*s35*s36 - 1024*q9a*s11*s12*s13*s15*s24*s26**3*s34*s35*s36 - 1024*q9a*s11*s12*s13*s14*s25*s26**3*s34*s35*s36 -  \
    4096*q8a*s11**3*s13*s22**3*s33*s34*s35*s36 + 8192*q8a*s11**3*s12*s22**2*s23*s33*s34*s35*s36 - 4096*q7a*s11**3*s22**3*s23*s33*s34*s35*s36 + 2048*q8a*s11**3*s14*s22**2*s24*s33*s34*s35*s36 -  \
    2048*q8a*s11**3*s12*s22*s24**2*s33*s34*s35*s36 - 512*q8a*s11**2*s16**2*s22*s24**2*s33*s34*s35*s36 - 2048*q8a*s11**2*s12*s15*s22**2*s25*s33*s34*s35*s36 - 2048*q8a*s11**2*s14*s16*s22**2*s25*s33*s34*s35*s36 +  \
    2048*q7a*s11**2*s15*s22**3*s25*s33*s34*s35*s36 + 2048*q8a*s11**2*s12*s16*s22*s24*s25*s33*s34*s35*s36 + 512*q8a*s11*s16**3*s22*s24*s25*s33*s34*s35*s36 + 2048*q8a*s11**2*s12**2*s22*s25**2*s33*s34*s35*s36 -  \
    1536*q8a*s11*s12*s16**2*s22*s25**2*s33*s34*s35*s36 - 2048*q7a*s11**2*s12*s22**2*s25**2*s33*s34*s35*s36 + 1536*q7a*s11*s16**2*s22**2*s25**2*s33*s34*s35*s36 + 4096*q8a*s11**2*s13*s16*s22**2*s26*s33*s34*s35*s36 -  \
    8192*q8a*s11**2*s12*s16*s22*s23*s26*s33*s34*s35*s36 + 4096*q7a*s11**2*s16*s22**2*s23*s26*s33*s34*s35*s36 - 1024*q8a*s11**2*s14*s16*s22*s24*s26*s33*s34*s35*s36 + 512*q8a*s11*s15*s16**2*s22*s24*s26*s33*s34*s35*s36 +  \
    2048*q8a*s11**2*s12*s16*s24**2*s26*s33*s34*s35*s36 + 4096*q8a*s11**2*s12*s14*s22*s25*s26*s33*s34*s35*s36 + 1024*q8a*s11*s12*s15*s16*s22*s25*s26*s33*s34*s35*s36 + 512*q8a*s11*s14*s16**2*s22*s25*s26*s33*s34*s35*s36 -  \
    2048*q7a*s11**2*s14*s22**2*s25*s26*s33*s34*s35*s36 - 1024*q7a*s11*s15*s16*s22**2*s25*s26*s33*s34*s35*s36 - 4096*q8a*s11**2*s12**2*s24*s25*s26*s33*s34*s35*s36 - 1024*q8a*s11*s12*s16**2*s24*s25*s26*s33*s34*s35*s36 +  \
    2048*q7a*s11**2*s12*s22*s24*s25*s26*s33*s34*s35*s36 - 1536*q7a*s11*s16**2*s22*s24*s25*s26*s33*s34*s35*s36 + 2048*q8a*s11*s12**2*s16*s25**2*s26*s33*s34*s35*s36 - 2048*q7a*s11*s12*s16*s22*s25**2*s26*s33*s34*s35*s36 -  \
    4096*q8a*s11**2*s12*s13*s22*s26**2*s33*s34*s35*s36 + 1536*q8a*s11**2*s14**2*s22*s26**2*s33*s34*s35*s36 + 512*q8a*s11*s12*s15**2*s22*s26**2*s33*s34*s35*s36 - 1536*q8a*s11*s14*s15*s16*s22*s26**2*s33*s34*s35*s36 -  \
    512*q7a*s11*s15**2*s22**2*s26**2*s33*s34*s35*s36 + 8192*q8a*s11**2*s12**2*s23*s26**2*s33*s34*s35*s36 - 4096*q7a*s11**2*s12*s22*s23*s26**2*s33*s34*s35*s36 - 2048*q8a*s11**2*s12*s14*s24*s26**2*s33*s34*s35*s36 +  \
    512*q7a*s11*s15*s16*s22*s24*s26**2*s33*s34*s35*s36 - 2048*q8a*s11*s12**2*s15*s25*s26**2*s33*s34*s35*s36 + 2048*q7a*s11*s12*s15*s22*s25*s26**2*s33*s34*s35*s36 + 512*q7a*s11*s14*s16*s22*s25*s26**2*s33*s34*s35*s36 +  \
    2048*q7a*s11*s12*s16*s24*s25*s26**2*s33*s34*s35*s36 + 1024*q8a*s11*s12*s14*s15*s26**3*s33*s34*s35*s36 + 512*q7a*s11*s14*s15*s22*s26**3*s33*s34*s35*s36 - 1024*q7a*s11*s12*s15*s24*s26**3*s33*s34*s35*s36 -  \
    1024*q7a*s11*s12*s14*s25*s26**3*s33*s34*s35*s36 - 4096*q8a*s11**3*s14*s22**2*s23*s34**2*s35*s36 + 1024*q8a*s11**3*s13*s22**2*s24*s34**2*s35*s36 + 512*q8a*s11**2*s15**2*s22**2*s24*s34**2*s35*s36 +  \
    2048*q8a*s11**3*s12*s22*s23*s24*s34**2*s35*s36 + 1024*q8a*s11**2*s16**2*s22*s23*s24*s34**2*s35*s36 + 1024*q7a*s11**3*s22**2*s23*s24*s34**2*s35*s36 - 512*q8a*s11**2*s15*s16*s22*s24**2*s34**2*s35*s36 +  \
    2560*q8a*s11**2*s14*s15*s22**2*s25*s34**2*s35*s36 - 512*q8a*s11**2*s13*s16*s22**2*s25*s34**2*s35*s36 - 512*q8a*s11*s15**2*s16*s22**2*s25*s34**2*s35*s36 + 1024*q8a*s11**2*s12*s16*s22*s23*s25*s34**2*s35*s36 -  \
    1024*q8a*s11*s16**3*s22*s23*s25*s34**2*s35*s36 - 512*q7a*s11**2*s16*s22**2*s23*s25*s34**2*s35*s36 - 2048*q8a*s11**2*s12*s15*s22*s24*s25*s34**2*s35*s36 - 512*q8a*s11**2*s14*s16*s22*s24*s25*s34**2*s35*s36 +  \
    512*q8a*s11*s15*s16**2*s22*s24*s25*s34**2*s35*s36 - 1536*q7a*s11**2*s15*s22**2*s24*s25*s34**2*s35*s36 + 1024*q7a*s11**2*s16*s22*s24**2*s25*s34**2*s35*s36 - 1024*q8a*s11**2*s12*s14*s22*s25**2*s34**2*s35*s36 +  \
    512*q8a*s11*s12*s15*s16*s22*s25**2*s34**2*s35*s36 + 512*q8a*s11*s14*s16**2*s22*s25**2*s34**2*s35*s36 - 1536*q7a*s11**2*s14*s22**2*s25**2*s34**2*s35*s36 + 512*q7a*s11*s15*s16*s22**2*s25**2*s34**2*s35*s36 +  \
    3072*q7a*s11**2*s12*s22*s24*s25**2*s34**2*s35*s36 - 1024*q7a*s11*s16**2*s22*s24*s25**2*s34**2*s35*s36 - 512*q7a*s11*s12*s16*s22*s25**3*s34**2*s35*s36 + 512*q8a*s11**2*s13*s15*s22**2*s26*s34**2*s35*s36 -  \
    512*q8a*s11*s15**3*s22**2*s26*s34**2*s35*s36 - 1024*q8a*s11**2*s12*s15*s22*s23*s26*s34**2*s35*s36 + 5120*q8a*s11**2*s14*s16*s22*s23*s26*s34**2*s35*s36 - 1024*q8a*s11*s15*s16**2*s22*s23*s26*s34**2*s35*s36 +  \
    512*q7a*s11**2*s15*s22**2*s23*s26*s34**2*s35*s36 - 512*q8a*s11**2*s14*s15*s22*s24*s26*s34**2*s35*s36 - 2048*q8a*s11**2*s13*s16*s22*s24*s26*s34**2*s35*s36 + 512*q8a*s11*s15**2*s16*s22*s24*s26*s34**2*s35*s36 -  \
    3072*q8a*s11**2*s12*s16*s23*s24*s26*s34**2*s35*s36 - 2048*q7a*s11**2*s16*s22*s23*s24*s26*s34**2*s35*s36 + 1024*q8a*s11**2*s12*s15*s24**2*s26*s34**2*s35*s36 - 2048*q8a*s11**2*s12*s13*s22*s25*s26*s34**2*s35*s36 -  \
    512*q8a*s11**2*s14**2*s22*s25*s26*s34**2*s35*s36 + 1536*q8a*s11*s12*s15**2*s22*s25*s26*s34**2*s35*s36 - 2048*q8a*s11*s14*s15*s16*s22*s25*s26*s34**2*s35*s36 + 1536*q8a*s11*s13*s16**2*s22*s25*s26*s34**2*s35*s36 +  \
    2048*q7a*s11**2*s13*s22**2*s25*s26*s34**2*s35*s36 + 512*q7a*s11*s15**2*s22**2*s25*s26*s34**2*s35*s36 + 2048*q8a*s11**2*s12**2*s23*s25*s26*s34**2*s35*s36 + 1024*q8a*s11*s12*s16**2*s23*s25*s26*s34**2*s35*s36 -  \
    2048*q7a*s11**2*s12*s22*s23*s25*s26*s34**2*s35*s36 + 1536*q7a*s11*s16**2*s22*s23*s25*s26*s34**2*s35*s36 + 1024*q8a*s11**2*s12*s14*s24*s25*s26*s34**2*s35*s36 + 1024*q7a*s11**2*s14*s22*s24*s25*s26*s34**2*s35*s36 -  \
    2048*q7a*s11**2*s12*s24**2*s25*s26*s34**2*s35*s36 - 1024*q8a*s11*s12**2*s15*s25**2*s26*s34**2*s35*s36 - 512*q8a*s11*s12*s14*s16*s25**2*s26*s34**2*s35*s36 - 1536*q7a*s11*s12*s15*s22*s25**2*s26*s34**2*s35*s36 +  \
    1536*q7a*s11*s14*s16*s22*s25**2*s26*s34**2*s35*s36 + 512*q7a*s11*s12*s16*s24*s25**2*s26*s34**2*s35*s36 + 1024*q7a*s11*s12**2*s25**3*s26*s34**2*s35*s36 - 1024*q8a*s11**2*s13*s14*s22*s26**2*s34**2*s35*s36 +  \
    512*q8a*s11*s14*s15**2*s22*s26**2*s34**2*s35*s36 + 512*q8a*s11*s13*s15*s16*s22*s26**2*s34**2*s35*s36 - 3072*q8a*s11**2*s12*s14*s23*s26**2*s34**2*s35*s36 + 1024*q8a*s11*s12*s15*s16*s23*s26**2*s34**2*s35*s36 -  \
    1024*q7a*s11**2*s14*s22*s23*s26**2*s34**2*s35*s36 + 512*q7a*s11*s15*s16*s22*s23*s26**2*s34**2*s35*s36 + 3072*q8a*s11**2*s12*s13*s24*s26**2*s34**2*s35*s36 - 1024*q8a*s11*s12*s15**2*s24*s26**2*s34**2*s35*s36 +  \
    3072*q7a*s11**2*s12*s23*s24*s26**2*s34**2*s35*s36 + 1536*q8a*s11*s12*s14*s15*s25*s26**2*s34**2*s35*s36 - 1536*q8a*s11*s12*s13*s16*s25*s26**2*s34**2*s35*s36 - 512*q7a*s11*s14*s15*s22*s25*s26**2*s34**2*s35*s36 -  \
    2048*q7a*s11*s13*s16*s22*s25*s26**2*s34**2*s35*s36 - 1536*q7a*s11*s12*s16*s23*s25*s26**2*s34**2*s35*s36 + 512*q7a*s11*s12*s15*s24*s25*s26**2*s34**2*s35*s36 - 1024*q7a*s11*s12*s14*s25**2*s26**2*s34**2*s35*s36 -  \
    512*q8a*s11*s12*s13*s15*s26**3*s34**2*s35*s36 - 512*q7a*s11*s12*s15*s23*s26**3*s34**2*s35*s36 + 2048*q7a*s11*s12*s13*s25*s26**3*s34**2*s35*s36 - 1024*q9a*s11**2*s14*s15*s22**3*s23*s35**2*s36 +  \
    1024*q9a*s11**2*s13*s16*s22**3*s23*s35**2*s36 - 1024*q9a*s11**2*s12*s16*s22**2*s23**2*s35**2*s36 - 1536*q9a*s11**2*s13*s15*s22**3*s24*s35**2*s36 + 2560*q9a*s11**2*s12*s15*s22**2*s23*s24*s35**2*s36 -  \
    1536*q9a*s11**2*s14*s16*s22**2*s23*s24*s35**2*s36 + 1024*q9a*s11**2*s14*s15*s22**2*s24**2*s35**2*s36 + 512*q9a*s11**2*s13*s16*s22**2*s24**2*s35**2*s36 + 1024*q9a*s11**2*s12*s16*s22*s23*s24**2*s35**2*s36 -  \
    1024*q9a*s11**2*s12*s15*s22*s24**3*s35**2*s36 + 512*q9a*s11**2*s13*s14*s22**3*s25*s35**2*s36 + 1024*q9a*s11*s13*s15*s16*s22**3*s25*s35**2*s36 + 512*q9a*s11**2*s12*s14*s22**2*s23*s25*s35**2*s36 -  \
    1024*q9a*s11*s12*s15*s16*s22**2*s23*s25*s35**2*s36 + 1024*q9a*s11*s14*s16**2*s22**2*s23*s25*s35**2*s36 + 1024*q9a*s11**2*s12*s13*s22**2*s24*s25*s35**2*s36 - 1024*q9a*s11**2*s14**2*s22**2*s24*s25*s35**2*s36 -  \
    512*q9a*s11*s14*s15*s16*s22**2*s24*s25*s35**2*s36 - 512*q9a*s11*s13*s16**2*s22**2*s24*s25*s35**2*s36 - 2048*q9a*s11**2*s12**2*s22*s23*s24*s25*s35**2*s36 - 512*q9a*s11*s12*s16**2*s22*s23*s24*s25*s35**2*s36 +  \
    1024*q9a*s11**2*s12*s14*s22*s24**2*s25*s35**2*s36 + 512*q9a*s11*s12*s15*s16*s22*s24**2*s25*s35**2*s36 - 1024*q9a*s11*s12*s13*s16*s22**2*s25**2*s35**2*s36 + 512*q9a*s11*s14**2*s16*s22**2*s25**2*s35**2*s36 +  \
    1024*q9a*s11*s12**2*s16*s22*s23*s25**2*s35**2*s36 - 512*q9a*s11*s12*s14*s16*s22*s24*s25**2*s35**2*s36 - 1024*q9a*s11**2*s13**2*s22**3*s26*s35**2*s36 + 1024*q9a*s11*s13*s15**2*s22**3*s26*s35**2*s36 +  \
    1024*q9a*s11**2*s12*s13*s22**2*s23*s26*s35**2*s36 - 512*q9a*s11**2*s14**2*s22**2*s23*s26*s35**2*s36 - 1024*q9a*s11*s12*s15**2*s22**2*s23*s26*s35**2*s36 + 2048*q9a*s11*s14*s15*s16*s22**2*s23*s26*s35**2*s36 -  \
    1024*q9a*s11*s13*s16**2*s22**2*s23*s26*s35**2*s36 + 1024*q9a*s11*s12*s16**2*s22*s23**2*s26*s35**2*s36 + 1536*q9a*s11**2*s13*s14*s22**2*s24*s26*s35**2*s36 - 512*q9a*s11*s14*s15**2*s22**2*s24*s26*s35**2*s36 +  \
    1024*q9a*s11**2*s12*s14*s22*s23*s24*s26*s35**2*s36 - 2048*q9a*s11*s12*s15*s16*s22*s23*s24*s26*s35**2*s36 + 512*q9a*s11*s14*s16**2*s22*s23*s24*s26*s35**2*s36 - 2048*q9a*s11**2*s12*s13*s22*s24**2*s26*s35**2*s36
v3_5= \
    512*q9a*s11*s12*s15**2*s22*s24**2*s26*s35**2*s36 - 512*q9a*s11*s14*s15*s16*s22*s24**2*s26*s35**2*s36 - 512*q9a*s11*s12*s16**2*s23*s24**2*s26*s35**2*s36 + 512*q9a*s11*s12*s15*s16*s24**3*s26*s35**2*s36 -  \
    3072*q9a*s11*s12*s13*s15*s22**2*s25*s26*s35**2*s36 + 512*q9a*s11*s14**2*s15*s22**2*s25*s26*s35**2*s36 - 1024*q9a*s11*s13*s14*s16*s22**2*s25*s26*s35**2*s36 + 3072*q9a*s11*s12**2*s15*s22*s23*s25*s26*s35**2*s36 -  \
    3072*q9a*s11*s12*s14*s16*s22*s23*s25*s26*s35**2*s36 + 512*q9a*s11*s12*s14*s15*s22*s24*s25*s26*s35**2*s36 + 2048*q9a*s11*s12*s13*s16*s22*s24*s25*s26*s35**2*s36 + 512*q9a*s11*s14**2*s16*s22*s24*s25*s26*s35**2*s36 +  \
    2048*q9a*s11*s12**2*s16*s23*s24*s25*s26*s35**2*s36 - 1024*q9a*s11*s12**2*s15*s24**2*s25*s26*s35**2*s36 - 512*q9a*s11*s12*s14*s16*s24**2*s25*s26*s35**2*s36 + 2048*q9a*s11*s12**2*s13*s22*s25**2*s26*s35**2*s36 -  \
    1024*q9a*s11*s12*s14**2*s22*s25**2*s26*s35**2*s36 - 2048*q9a*s11*s12**3*s23*s25**2*s26*s35**2*s36 + 1024*q9a*s11*s12**2*s14*s24*s25**2*s26*s35**2*s36 - 1536*q9a*s11*s13*s14*s15*s22**2*s26**2*s35**2*s36 +  \
    1024*q9a*s11*s13**2*s16*s22**2*s26**2*s35**2*s36 - 512*q9a*s11*s12*s14*s15*s22*s23*s26**2*s35**2*s36 - 512*q9a*s11*s14**2*s16*s22*s23*s26**2*s35**2*s36 - 1024*q9a*s11*s12**2*s16*s23**2*s26**2*s35**2*s36 +  \
    1536*q9a*s11*s12*s13*s15*s22*s24*s26**2*s35**2*s36 + 512*q9a*s11*s14**2*s15*s22*s24*s26**2*s35**2*s36 - 512*q9a*s11*s13*s14*s16*s22*s24*s26**2*s35**2*s36 + 512*q9a*s11*s12**2*s15*s23*s24*s26**2*s35**2*s36 +  \
    512*q9a*s11*s12*s14*s16*s23*s24*s26**2*s35**2*s36 - 512*q9a*s11*s12*s14*s15*s24**2*s26**2*s35**2*s36 + 512*q9a*s11*s12*s13*s16*s24**2*s26**2*s35**2*s36 + 2560*q9a*s11*s12*s13*s14*s22*s25*s26**2*s35**2*s36 -  \
    512*q9a*s11*s14**3*s22*s25*s26**2*s35**2*s36 + 512*q9a*s11*s12**2*s14*s23*s25*s26**2*s35**2*s36 - 3072*q9a*s11*s12**2*s13*s24*s25*s26**2*s35**2*s36 + 512*q9a*s11*s12*s14**2*s24*s25*s26**2*s35**2*s36 -  \
    1024*q9a*s11*s12*s13**2*s22*s26**3*s35**2*s36 + 512*q9a*s11*s13*s14**2*s22*s26**3*s35**2*s36 + 1024*q9a*s11*s12**2*s13*s23*s26**3*s35**2*s36 - 512*q9a*s11*s12*s13*s14*s24*s26**3*s35**2*s36 -  \
    1024*q8a*s11**2*s14*s15*s22**3*s33*s35**2*s36 + 1024*q8a*s11**2*s13*s16*s22**3*s33*s35**2*s36 - 2048*q8a*s11**2*s12*s16*s22**2*s23*s33*s35**2*s36 + 1024*q7a*s11**2*s16*s22**3*s23*s33*s35**2*s36 +  \
    2560*q8a*s11**2*s12*s15*s22**2*s24*s33*s35**2*s36 - 1536*q8a*s11**2*s14*s16*s22**2*s24*s33*s35**2*s36 - 1536*q7a*s11**2*s15*s22**3*s24*s33*s35**2*s36 + 1024*q8a*s11**2*s12*s16*s22*s24**2*s33*s35**2*s36 +  \
    512*q7a*s11**2*s16*s22**2*s24**2*s33*s35**2*s36 + 512*q8a*s11**2*s12*s14*s22**2*s25*s33*s35**2*s36 - 1024*q8a*s11*s12*s15*s16*s22**2*s25*s33*s35**2*s36 + 1024*q8a*s11*s14*s16**2*s22**2*s25*s33*s35**2*s36 +  \
    512*q7a*s11**2*s14*s22**3*s25*s33*s35**2*s36 + 1024*q7a*s11*s15*s16*s22**3*s25*s33*s35**2*s36 - 2048*q8a*s11**2*s12**2*s22*s24*s25*s33*s35**2*s36 - 512*q8a*s11*s12*s16**2*s22*s24*s25*s33*s35**2*s36 +  \
    1024*q7a*s11**2*s12*s22**2*s24*s25*s33*s35**2*s36 - 512*q7a*s11*s16**2*s22**2*s24*s25*s33*s35**2*s36 + 1024*q8a*s11*s12**2*s16*s22*s25**2*s33*s35**2*s36 - 1024*q7a*s11*s12*s16*s22**2*s25**2*s33*s35**2*s36 +  \
    1024*q8a*s11**2*s12*s13*s22**2*s26*s33*s35**2*s36 - 512*q8a*s11**2*s14**2*s22**2*s26*s33*s35**2*s36 - 1024*q8a*s11*s12*s15**2*s22**2*s26*s33*s35**2*s36 + 2048*q8a*s11*s14*s15*s16*s22**2*s26*s33*s35**2*s36 -  \
    1024*q8a*s11*s13*s16**2*s22**2*s26*s33*s35**2*s36 - 2048*q7a*s11**2*s13*s22**3*s26*s33*s35**2*s36 + 1024*q7a*s11*s15**2*s22**3*s26*s33*s35**2*s36 + 2048*q8a*s11*s12*s16**2*s22*s23*s26*s33*s35**2*s36 +  \
    1024*q7a*s11**2*s12*s22**2*s23*s26*s33*s35**2*s36 - 1024*q7a*s11*s16**2*s22**2*s23*s26*s33*s35**2*s36 + 1024*q8a*s11**2*s12*s14*s22*s24*s26*s33*s35**2*s36 - 2048*q8a*s11*s12*s15*s16*s22*s24*s26*s33*s35**2*s36 +  \
    512*q8a*s11*s14*s16**2*s22*s24*s26*s33*s35**2*s36 + 1536*q7a*s11**2*s14*s22**2*s24*s26*s33*s35**2*s36 - 512*q8a*s11*s12*s16**2*s24**2*s26*s33*s35**2*s36 - 2048*q7a*s11**2*s12*s22*s24**2*s26*s33*s35**2*s36 +  \
    3072*q8a*s11*s12**2*s15*s22*s25*s26*s33*s35**2*s36 - 3072*q8a*s11*s12*s14*s16*s22*s25*s26*s33*s35**2*s36 - 3072*q7a*s11*s12*s15*s22**2*s25*s26*s33*s35**2*s36 - 1024*q7a*s11*s14*s16*s22**2*s25*s26*s33*s35**2*s36 +  \
    2048*q8a*s11*s12**2*s16*s24*s25*s26*s33*s35**2*s36 + 2048*q7a*s11*s12*s16*s22*s24*s25*s26*s33*s35**2*s36 - 2048*q8a*s11*s12**3*s25**2*s26*s33*s35**2*s36 + 2048*q7a*s11*s12**2*s22*s25**2*s26*s33*s35**2*s36 -  \
    512*q8a*s11*s12*s14*s15*s22*s26**2*s33*s35**2*s36 - 512*q8a*s11*s14**2*s16*s22*s26**2*s33*s35**2*s36 - 1536*q7a*s11*s14*s15*s22**2*s26**2*s33*s35**2*s36 + 2048*q7a*s11*s13*s16*s22**2*s26**2*s33*s35**2*s36 -  \
    2048*q8a*s11*s12**2*s16*s23*s26**2*s33*s35**2*s36 + 512*q8a*s11*s12**2*s15*s24*s26**2*s33*s35**2*s36 + 512*q8a*s11*s12*s14*s16*s24*s26**2*s33*s35**2*s36 + 1536*q7a*s11*s12*s15*s22*s24*s26**2*s33*s35**2*s36 -  \
    512*q7a*s11*s14*s16*s22*s24*s26**2*s33*s35**2*s36 + 512*q7a*s11*s12*s16*s24**2*s26**2*s33*s35**2*s36 + 512*q8a*s11*s12**2*s14*s25*s26**2*s33*s35**2*s36 + 2560*q7a*s11*s12*s14*s22*s25*s26**2*s33*s35**2*s36 -  \
    3072*q7a*s11*s12**2*s24*s25*s26**2*s33*s35**2*s36 + 1024*q8a*s11*s12**2*s13*s26**3*s33*s35**2*s36 - 2048*q7a*s11*s12*s13*s22*s26**3*s33*s35**2*s36 + 512*q7a*s11*s14**2*s22*s26**3*s33*s35**2*s36 +  \
    1024*q7a*s11*s12**2*s23*s26**3*s33*s35**2*s36 - 512*q7a*s11*s12*s14*s24*s26**3*s33*s35**2*s36 + 1024*q8a*s11**2*s13*s15*s22**3*s34*s35**2*s36 - 2048*q8a*s11**2*s12*s15*s22**2*s23*s34*s35**2*s36 +  \
    2048*q8a*s11**2*s14*s16*s22**2*s23*s34*s35**2*s36 + 1024*q7a*s11**2*s15*s22**3*s23*s34*s35**2*s36 - 1536*q8a*s11**2*s14*s15*s22**2*s24*s34*s35**2*s36 + 512*q8a*s11**2*s13*s16*s22**2*s24*s34*s35**2*s36 -  \
    3072*q8a*s11**2*s12*s16*s22*s23*s24*s34*s35**2*s36 + 512*q7a*s11**2*s16*s22**2*s23*s24*s34*s35**2*s36 + 1024*q8a*s11**2*s12*s15*s22*s24**2*s34*s35**2*s36 + 512*q8a*s11**2*s14*s16*s22*s24**2*s34*s35**2*s36 +  \
    512*q7a*s11**2*s15*s22**2*s24**2*s34*s35**2*s36 - 512*q7a*s11**2*s16*s22*s24**3*s34*s35**2*s36 + 3072*q8a*s11**2*s12*s13*s22**2*s25*s34*s35**2*s36 - 1536*q8a*s11**2*s14**2*s22**2*s25*s34*s35**2*s36 +  \
    1024*q8a*s11*s14*s15*s16*s22**2*s25*s34*s35**2*s36 - 1024*q8a*s11*s13*s16**2*s22**2*s25*s34*s35**2*s36 - 4096*q7a*s11**2*s13*s22**3*s25*s34*s35**2*s36 - 2048*q8a*s11**2*s12**2*s22*s23*s25*s34*s35**2*s36 +  \
    2048*q8a*s11*s12*s16**2*s22*s23*s25*s34*s35**2*s36 + 3072*q7a*s11**2*s12*s22**2*s23*s25*s34*s35**2*s36 - 1024*q7a*s11*s16**2*s22**2*s23*s25*s34*s35**2*s36 + 2048*q8a*s11**2*s12*s14*s22*s24*s25*s34*s35**2*s36 -  \
    512*q8a*s11*s12*s15*s16*s22*s24*s25*s34*s35**2*s36 - 512*q8a*s11*s14*s16**2*s22*s24*s25*s34*s35**2*s36 + 2560*q7a*s11**2*s14*s22**2*s24*s25*s34*s35**2*s36 - 512*q7a*s11*s15*s16*s22**2*s24*s25*s34*s35**2*s36 -  \
    3072*q7a*s11**2*s12*s22*s24**2*s25*s34*s35**2*s36 + 512*q7a*s11*s16**2*s22*s24**2*s25*s34*s35**2*s36 - 512*q8a*s11*s12*s14*s16*s22*s25**2*s34*s35**2*s36 - 512*q7a*s11*s14*s16*s22**2*s25**2*s34*s35**2*s36 +  \
    1024*q7a*s11*s12*s16*s22*s24*s25**2*s34*s35**2*s36 - 512*q8a*s11**2*s13*s14*s22**2*s26*s34*s35**2*s36 + 1024*q8a*s11*s14*s15**2*s22**2*s26*s34*s35**2*s36 - 2048*q8a*s11*s13*s15*s16*s22**2*s26*s34*s35**2*s36 -  \
    1024*q8a*s11**2*s12*s14*s22*s23*s26*s34*s35**2*s36 + 4096*q8a*s11*s12*s15*s16*s22*s23*s26*s34*s35**2*s36 - 2048*q8a*s11*s14*s16**2*s22*s23*s26*s34*s35**2*s36 - 512*q7a*s11**2*s14*s22**2*s23*s26*s34*s35**2*s36 -  \
    2048*q7a*s11*s15*s16*s22**2*s23*s26*s34*s35**2*s36 + 512*q8a*s11**2*s14**2*s22*s24*s26*s34*s35**2*s36 - 512*q8a*s11*s12*s15**2*s22*s24*s26*s34*s35**2*s36 + 512*q8a*s11*s13*s16**2*s22*s24*s26*s34*s35**2*s36 -  \
    512*q7a*s11*s15**2*s22**2*s24*s26*s34*s35**2*s36 + 2048*q8a*s11**2*s12**2*s23*s24*s26*s34*s35**2*s36 + 1024*q8a*s11*s12*s16**2*s23*s24*s26*s34*s35**2*s36 + 512*q7a*s11*s16**2*s22*s23*s24*s26*s34*s35**2*s36 -  \
    1024*q8a*s11**2*s12*s14*s24**2*s26*s34*s35**2*s36 - 512*q8a*s11*s12*s15*s16*s24**2*s26*s34*s35**2*s36 - 512*q7a*s11**2*s14*s22*s24**2*s26*s34*s35**2*s36 + 512*q7a*s11*s15*s16*s22*s24**2*s26*s34*s35**2*s36 +  \
    1024*q7a*s11**2*s12*s24**3*s26*s34*s35**2*s36 - 2560*q8a*s11*s12*s14*s15*s22*s25*s26*s34*s35**2*s36 - 1024*q8a*s11*s12*s13*s16*s22*s25*s26*s34*s35**2*s36 + 1536*q8a*s11*s14**2*s16*s22*s25*s26*s34*s35**2*s36 -  \
    512*q7a*s11*s14*s15*s22**2*s25*s26*s34*s35**2*s36 + 5120*q7a*s11*s13*s16*s22**2*s25*s26*s34*s35**2*s36 - 3072*q8a*s11*s12**2*s16*s23*s25*s26*s34*s35**2*s36 - 1024*q7a*s11*s12*s16*s22*s23*s25*s26*s34*s35**2*s36 +  \
    1024*q8a*s11*s12**2*s15*s24*s25*s26*s34*s35**2*s36 + 2048*q7a*s11*s12*s15*s22*s24*s25*s26*s34*s35**2*s36 - 2048*q7a*s11*s14*s16*s22*s24*s25*s26*s34*s35**2*s36 + 512*q7a*s11*s12*s16*s24**2*s25*s26*s34*s35**2*s36 +  \
    1024*q8a*s11*s12**2*s14*s25**2*s26*s34*s35**2*s36 + 1024*q7a*s11*s12*s14*s22*s25**2*s26*s34*s35**2*s36 - 2048*q7a*s11*s12**2*s24*s25**2*s26*s34*s35**2*s36 + 1024*q8a*s11*s12*s13*s15*s22*s26**2*s34*s35**2*s36 -  \
    1024*q8a*s11*s14**2*s15*s22*s26**2*s34*s35**2*s36 + 1536*q8a*s11*s13*s14*s16*s22*s26**2*s34*s35**2*s36 + 1024*q7a*s11*s13*s15*s22**2*s26**2*s34*s35**2*s36 - 3072*q8a*s11*s12**2*s15*s23*s26**2*s34*s35**2*s36 +  \
    1024*q8a*s11*s12*s14*s16*s23*s26**2*s34*s35**2*s36 + 1024*q7a*s11*s12*s15*s22*s23*s26**2*s34*s35**2*s36 + 1536*q7a*s11*s14*s16*s22*s23*s26**2*s34*s35**2*s36 + 1536*q8a*s11*s12*s14*s15*s24*s26**2*s34*s35**2*s36 -  \
    1536*q8a*s11*s12*s13*s16*s24*s26**2*s34*s35**2*s36 + 512*q7a*s11*s14*s15*s22*s24*s26**2*s34*s35**2*s36 - 1024*q7a*s11*s13*s16*s22*s24*s26**2*s34*s35**2*s36 - 1536*q7a*s11*s12*s16*s23*s24*s26**2*s34*s35**2*s36 -  \
    1024*q7a*s11*s12*s15*s24**2*s26**2*s34*s35**2*s36 + 3072*q8a*s11*s12**2*s13*s25*s26**2*s34*s35**2*s36 - 1024*q8a*s11*s12*s14**2*s25*s26**2*s34*s35**2*s36 - 6144*q7a*s11*s12*s13*s22*s25*s26**2*s34*s35**2*s36 +  \
    512*q7a*s11*s14**2*s22*s25*s26**2*s34*s35**2*s36 + 3072*q7a*s11*s12**2*s23*s25*s26**2*s34*s35**2*s36 + 512*q7a*s11*s12*s14*s24*s25*s26**2*s34*s35**2*s36 - 512*q8a*s11*s12*s13*s14*s26**3*s34*s35**2*s36 -  \
    1024*q7a*s11*s13*s14*s22*s26**3*s34*s35**2*s36 - 512*q7a*s11*s12*s14*s23*s26**3*s34*s35**2*s36 + 2048*q7a*s11*s12*s13*s24*s26**3*s34*s35**2*s36 - 2048*q8a*s11**2*s12*s13*s22**2*s24*s35**3*s36 +  \
    1024*q8a*s11**2*s14**2*s22**2*s24*s35**3*s36 + 2048*q7a*s11**2*s13*s22**3*s24*s35**3*s36 + 2048*q8a*s11**2*s12**2*s22*s23*s24*s35**3*s36 - 2048*q7a*s11**2*s12*s22**2*s23*s24*s35**3*s36 -  \
    1024*q8a*s11**2*s12*s14*s22*s24**2*s35**3*s36 - 1024*q7a*s11**2*s14*s22**2*s24**2*s35**3*s36 + 1024*q7a*s11**2*s12*s22*s24**3*s35**3*s36 + 1024*q8a*s11*s12*s13*s16*s22**2*s25*s35**3*s36 -  \
    512*q8a*s11*s14**2*s16*s22**2*s25*s35**3*s36 - 1024*q7a*s11*s13*s16*s22**3*s25*s35**3*s36 - 1024*q8a*s11*s12**2*s16*s22*s23*s25*s35**3*s36 + 1024*q7a*s11*s12*s16*s22**2*s23*s25*s35**3*s36 +  \
    512*q8a*s11*s12*s14*s16*s22*s24*s25*s35**3*s36 + 512*q7a*s11*s14*s16*s22**2*s24*s25*s35**3*s36 - 512*q7a*s11*s12*s16*s22*s24**2*s25*s35**3*s36 + 1024*q8a*s11*s12*s13*s15*s22**2*s26*s35**3*s36 -  \
    512*q8a*s11*s14**2*s15*s22**2*s26*s35**3*s36 - 1024*q7a*s11*s13*s15*s22**3*s26*s35**3*s36 - 1024*q8a*s11*s12**2*s15*s22*s23*s26*s35**3*s36 + 1024*q7a*s11*s12*s15*s22**2*s23*s26*s35**3*s36 +  \
    512*q8a*s11*s12*s14*s15*s22*s24*s26*s35**3*s36 + 1024*q8a*s11*s12*s13*s16*s22*s24*s26*s35**3*s36 - 512*q8a*s11*s14**2*s16*s22*s24*s26*s35**3*s36 + 512*q7a*s11*s14*s15*s22**2*s24*s26*s35**3*s36 -  \
    1024*q7a*s11*s13*s16*s22**2*s24*s26*s35**3*s36 - 1024*q8a*s11*s12**2*s16*s23*s24*s26*s35**3*s36 + 1024*q7a*s11*s12*s16*s22*s23*s24*s26*s35**3*s36 + 512*q8a*s11*s12*s14*s16*s24**2*s26*s35**3*s36 -  \
    512*q7a*s11*s12*s15*s22*s24**2*s26*s35**3*s36 + 512*q7a*s11*s14*s16*s22*s24**2*s26*s35**3*s36 - 512*q7a*s11*s12*s16*s24**3*s26*s35**3*s36 - 2048*q8a*s11*s12**2*s13*s22*s25*s26*s35**3*s36 +  \
    1024*q8a*s11*s12*s14**2*s22*s25*s26*s35**3*s36 + 2048*q7a*s11*s12*s13*s22**2*s25*s26*s35**3*s36 + 2048*q8a*s11*s12**3*s23*s25*s26*s35**3*s36 - 2048*q7a*s11*s12**2*s22*s23*s25*s26*s35**3*s36 -  \
    1024*q8a*s11*s12**2*s14*s24*s25*s26*s35**3*s36 - 1024*q7a*s11*s12*s14*s22*s24*s25*s26*s35**3*s36 + 1024*q7a*s11*s12**2*s24**2*s25*s26*s35**3*s36 - 1024*q8a*s11*s12*s13*s14*s22*s26**2*s35**3*s36
v3_6= \
    512*q8a*s11*s14**3*s22*s26**2*s35**3*s36 + 1024*q7a*s11*s13*s14*s22**2*s26**2*s35**3*s36 + 1024*q8a*s11*s12**2*s14*s23*s26**2*s35**3*s36 - 1024*q7a*s11*s12*s14*s22*s23*s26**2*s35**3*s36 -  \
    512*q8a*s11*s12*s14**2*s24*s26**2*s35**3*s36 - 512*q7a*s11*s14**2*s22*s24*s26**2*s35**3*s36 + 512*q7a*s11*s12*s14*s24**2*s26**2*s35**3*s36 + 2048*q9a*s11**3*s13*s22**3*s23*s33*s36**2 - 1024*q9a*s11**2*s15**2*s22**3*s23*s33*s36**2 -  \
    2048*q9a*s11**3*s12*s22**2*s23**2*s33*s36**2 - 1024*q9a*s11**2*s16**2*s22**2*s23**2*s33*s36**2 - 3072*q9a*s11**3*s14*s22**2*s23*s24*s33*s36**2 + 1024*q9a*s11**2*s15*s16*s22**2*s23*s24*s33*s36**2 -  \
    1024*q9a*s11**3*s13*s22**2*s24**2*s33*s36**2 + 4096*q9a*s11**3*s12*s22*s23*s24**2*s33*s36**2 + 1024*q9a*s11**3*s14*s22*s24**3*s33*s36**2 - 1024*q9a*s11**3*s12*s24**4*s33*s36**2 - 3072*q9a*s11**2*s13*s15*s22**3*s25*s33*s36**2 +  \
    1024*q9a*s11*s15**3*s22**3*s25*s33*s36**2 + 5120*q9a*s11**2*s12*s15*s22**2*s23*s25*s33*s36**2 - 1024*q9a*s11**2*s14*s16*s22**2*s23*s25*s33*s36**2 + 1024*q9a*s11*s15*s16**2*s22**2*s23*s25*s33*s36**2 +  \
    1024*q9a*s11**2*s14*s15*s22**2*s24*s25*s33*s36**2 + 4096*q9a*s11**2*s13*s16*s22**2*s24*s25*s33*s36**2 - 1024*q9a*s11*s15**2*s16*s22**2*s24*s25*s33*s36**2 - 4096*q9a*s11**2*s12*s16*s22*s23*s24*s25*s33*s36**2 -  \
    1024*q9a*s11**2*s12*s15*s22*s24**2*s25*s33*s36**2 - 2048*q9a*s11**2*s14*s16*s22*s24**2*s25*s33*s36**2 + 2048*q9a*s11**2*s12*s16*s24**3*s25*s33*s36**2 + 3072*q9a*s11**2*s12*s13*s22**2*s25**2*s33*s36**2 -  \
    3072*q9a*s11**2*s14**2*s22**2*s25**2*s33*s36**2 - 3072*q9a*s11*s12*s15**2*s22**2*s25**2*s33*s36**2 + 3072*q9a*s11*s14*s15*s16*s22**2*s25**2*s33*s36**2 - 3072*q9a*s11*s13*s16**2*s22**2*s25**2*s33*s36**2 -  \
    4096*q9a*s11**2*s12**2*s22*s23*s25**2*s33*s36**2 + 2048*q9a*s11*s12*s16**2*s22*s23*s25**2*s33*s36**2 + 5120*q9a*s11**2*s12*s14*s22*s24*s25**2*s33*s36**2 - 1024*q9a*s11*s12*s15*s16*s22*s24*s25**2*s33*s36**2 +  \
    1024*q9a*s11*s14*s16**2*s22*s24*s25**2*s33*s36**2 - 2048*q9a*s11**2*s12**2*s24**2*s25**2*s33*s36**2 - 1024*q9a*s11*s12*s16**2*s24**2*s25**2*s33*s36**2 + 3072*q9a*s11*s12**2*s15*s22*s25**3*s33*s36**2 -  \
    3072*q9a*s11*s12*s14*s16*s22*s25**3*s33*s36**2 + 2048*q9a*s11*s12**2*s16*s24*s25**3*s33*s36**2 - 1024*q9a*s11*s12**3*s25**4*s33*s36**2 + 4096*q9a*s11**2*s14*s15*s22**2*s23*s26*s33*s36**2 -  \
    4096*q9a*s11**2*s13*s16*s22**2*s23*s26*s33*s36**2 + 6144*q9a*s11**2*s12*s16*s22*s23**2*s26*s33*s36**2 + 1024*q9a*s11**2*s13*s15*s22**2*s24*s26*s33*s36**2 - 6144*q9a*s11**2*s12*s15*s22*s23*s24*s26*s33*s36**2 +  \
    2048*q9a*s11**2*s14*s16*s22*s23*s24*s26*s33*s36**2 - 1024*q9a*s11**2*s14*s15*s22*s24**2*s26*s33*s36**2 - 2048*q9a*s11**2*s12*s16*s23*s24**2*s26*s33*s36**2 + 1024*q9a*s11**2*s12*s15*s24**3*s26*s33*s36**2 -  \
    1024*q9a*s11**2*s13*s14*s22**2*s25*s26*s33*s36**2 - 2048*q9a*s11*s14*s15**2*s22**2*s25*s26*s33*s36**2 + 2048*q9a*s11*s13*s15*s16*s22**2*s25*s26*s33*s36**2 - 2048*q9a*s11**2*s12*s14*s22*s23*s25*s26*s33*s36**2 -  \
    4096*q9a*s11*s12*s15*s16*s22*s23*s25*s26*s33*s36**2 - 4096*q9a*s11**2*s12*s13*s22*s24*s25*s26*s33*s36**2 + 3072*q9a*s11**2*s14**2*s22*s24*s25*s26*s33*s36**2 + 3072*q9a*s11*s12*s15**2*s22*s24*s25*s26*s33*s36**2 -  \
    1024*q9a*s11*s14*s15*s16*s22*s24*s25*s26*s33*s36**2 + 8192*q9a*s11**2*s12**2*s23*s24*s25*s26*s33*s36**2 - 3072*q9a*s11**2*s12*s14*s24**2*s25*s26*s33*s36**2 + 1024*q9a*s11*s12*s15*s16*s24**2*s25*s26*s33*s36**2 +  \
    1024*q9a*s11*s12*s14*s15*s22*s25**2*s26*s33*s36**2 + 4096*q9a*s11*s12*s13*s16*s22*s25**2*s26*s33*s36**2 - 1024*q9a*s11*s14**2*s16*s22*s25**2*s26*s33*s36**2 - 2048*q9a*s11*s12**2*s16*s23*s25**2*s26*s33*s36**2 -  \
    3072*q9a*s11*s12**2*s15*s24*s25**2*s26*s33*s36**2 + 1024*q9a*s11*s12*s14*s16*s24*s25**2*s26*s33*s36**2 + 1024*q9a*s11*s12**2*s14*s25**3*s26*s33*s36**2 - 1024*q9a*s11**2*s13**2*s22**2*s26**2*s33*s36**2 +  \
    6144*q9a*s11**2*s12*s13*s22*s23*s26**2*s33*s36**2 - 3072*q9a*s11**2*s14**2*s22*s23*s26**2*s33*s36**2 - 6144*q9a*s11**2*s12**2*s23**2*s26**2*s33*s36**2 + 1024*q9a*s11**2*s13*s14*s22*s24*s26**2*s33*s36**2 +  \
    3072*q9a*s11**2*s12*s14*s23*s24*s26**2*s33*s36**2 - 1024*q9a*s11**2*s12*s13*s24**2*s26**2*s33*s36**2 - 2048*q9a*s11*s12*s13*s15*s22*s25*s26**2*s33*s36**2 + 1024*q9a*s11*s14**2*s15*s22*s25*s26**2*s33*s36**2 +  \
    3072*q9a*s11*s12**2*s15*s23*s25*s26**2*s33*s36**2 - 1024*q9a*s11*s12*s14*s15*s24*s25*s26**2*s33*s36**2 - 1024*q9a*s11*s12**2*s13*s25**2*s26**2*s33*s36**2 + 1024*q8a*s11**3*s13*s22**3*s33**2*s36**2 -  \
    512*q8a*s11**2*s15**2*s22**3*s33**2*s36**2 - 2048*q8a*s11**3*s12*s22**2*s23*s33**2*s36**2 - 1024*q8a*s11**2*s16**2*s22**2*s23*s33**2*s36**2 + 1024*q7a*s11**3*s22**3*s23*s33**2*s36**2 - 1536*q8a*s11**3*s14*s22**2*s24*s33**2*s36**2 +  \
    512*q8a*s11**2*s15*s16*s22**2*s24*s33**2*s36**2 + 2048*q8a*s11**3*s12*s22*s24**2*s33**2*s36**2 - 512*q7a*s11**3*s22**2*s24**2*s33**2*s36**2 + 2560*q8a*s11**2*s12*s15*s22**2*s25*s33**2*s36**2 -  \
    512*q8a*s11**2*s14*s16*s22**2*s25*s33**2*s36**2 + 512*q8a*s11*s15*s16**2*s22**2*s25*s33**2*s36**2 - 1536*q7a*s11**2*s15*s22**3*s25*s33**2*s36**2 - 2048*q8a*s11**2*s12*s16*s22*s24*s25*s33**2*s36**2 +  \
    2048*q7a*s11**2*s16*s22**2*s24*s25*s33**2*s36**2 - 2048*q8a*s11**2*s12**2*s22*s25**2*s33**2*s36**2 + 1024*q8a*s11*s12*s16**2*s22*s25**2*s33**2*s36**2 + 1536*q7a*s11**2*s12*s22**2*s25**2*s33**2*s36**2 -  \
    1536*q7a*s11*s16**2*s22**2*s25**2*s33**2*s36**2 + 2048*q8a*s11**2*s14*s15*s22**2*s26*s33**2*s36**2 - 2048*q8a*s11**2*s13*s16*s22**2*s26*s33**2*s36**2 + 6144*q8a*s11**2*s12*s16*s22*s23*s26*s33**2*s36**2 -  \
    2048*q7a*s11**2*s16*s22**2*s23*s26*s33**2*s36**2 - 3072*q8a*s11**2*s12*s15*s22*s24*s26*s33**2*s36**2 + 1024*q8a*s11**2*s14*s16*s22*s24*s26*s33**2*s36**2 + 512*q7a*s11**2*s15*s22**2*s24*s26*s33**2*s36**2 -  \
    1024*q8a*s11**2*s12*s16*s24**2*s26*s33**2*s36**2 - 1024*q8a*s11**2*s12*s14*s22*s25*s26*s33**2*s36**2 - 2048*q8a*s11*s12*s15*s16*s22*s25*s26*s33**2*s36**2 - 512*q7a*s11**2*s14*s22**2*s25*s26*s33**2*s36**2 +  \
    1024*q7a*s11*s15*s16*s22**2*s25*s26*s33**2*s36**2 + 4096*q8a*s11**2*s12**2*s24*s25*s26*s33**2*s36**2 - 2048*q7a*s11**2*s12*s22*s24*s25*s26*s33**2*s36**2 - 1024*q8a*s11*s12**2*s16*s25**2*s26*s33**2*s36**2 +  \
    2048*q7a*s11*s12*s16*s22*s25**2*s26*s33**2*s36**2 + 3072*q8a*s11**2*s12*s13*s22*s26**2*s33**2*s36**2 - 1536*q8a*s11**2*s14**2*s22*s26**2*s33**2*s36**2 - 1024*q7a*s11**2*s13*s22**2*s26**2*s33**2*s36**2 -  \
    6144*q8a*s11**2*s12**2*s23*s26**2*s33**2*s36**2 + 3072*q7a*s11**2*s12*s22*s23*s26**2*s33**2*s36**2 + 1536*q8a*s11**2*s12*s14*s24*s26**2*s33**2*s36**2 + 512*q7a*s11**2*s14*s22*s24*s26**2*s33**2*s36**2 -  \
    512*q7a*s11**2*s12*s24**2*s26**2*s33**2*s36**2 + 1536*q8a*s11*s12**2*s15*s25*s26**2*s33**2*s36**2 - 1024*q7a*s11*s12*s15*s22*s25*s26**2*s33**2*s36**2 - 512*q7a*s11*s12**2*s25**2*s26**2*s33**2*s36**2 +  \
    2048*q9a*s11**3*s14*s22**2*s23**2*s34*s36**2 + 1024*q9a*s11**3*s13*s22**2*s23*s24*s34*s36**2 + 512*q9a*s11**2*s15**2*s22**2*s23*s24*s34*s36**2 - 3072*q9a*s11**3*s12*s22*s23**2*s24*s34*s36**2 +  \
    512*q9a*s11**2*s16**2*s22*s23**2*s24*s34*s36**2 - 1024*q9a*s11**3*s14*s22*s23*s24**2*s34*s36**2 - 512*q9a*s11**2*s15*s16*s22*s23*s24**2*s34*s36**2 + 1024*q9a*s11**3*s12*s23*s24**3*s34*s36**2 -  \
    2560*q9a*s11**2*s14*s15*s22**2*s23*s25*s34*s36**2 - 512*q9a*s11**2*s13*s16*s22**2*s23*s25*s34*s36**2 + 512*q9a*s11**2*s12*s16*s22*s23**2*s25*s34*s36**2 + 512*q9a*s11**2*s13*s15*s22**2*s24*s25*s34*s36**2 -  \
    512*q9a*s11*s15**3*s22**2*s24*s25*s34*s36**2 + 1024*q9a*s11**2*s12*s15*s22*s23*s24*s25*s34*s36**2 + 3072*q9a*s11**2*s14*s16*s22*s23*s24*s25*s34*s36**2 - 512*q9a*s11*s15*s16**2*s22*s23*s24*s25*s34*s36**2 +  \
    512*q9a*s11**2*s14*s15*s22*s24**2*s25*s34*s36**2 - 1024*q9a*s11**2*s13*s16*s22*s24**2*s25*s34*s36**2 + 512*q9a*s11*s15**2*s16*s22*s24**2*s25*s34*s36**2 - 1536*q9a*s11**2*s12*s16*s23*s24**2*s25*s34*s36**2 -  \
    512*q9a*s11**2*s12*s15*s24**3*s25*s34*s36**2 + 2560*q9a*s11**2*s13*s14*s22**2*s25**2*s34*s36**2 + 512*q9a*s11*s14*s15**2*s22**2*s25**2*s34*s36**2 - 512*q9a*s11*s13*s15*s16*s22**2*s25**2*s34*s36**2 +  \
    512*q9a*s11*s12*s15*s16*s22*s23*s25**2*s34*s36**2 - 1024*q9a*s11*s14*s16**2*s22*s23*s25**2*s34*s36**2 - 3072*q9a*s11**2*s12*s13*s22*s24*s25**2*s34*s36**2 - 512*q9a*s11**2*s14**2*s22*s24*s25**2*s34*s36**2 +  \
    1024*q9a*s11*s12*s15**2*s22*s24*s25**2*s34*s36**2 - 1536*q9a*s11*s14*s15*s16*s22*s24*s25**2*s34*s36**2 + 1024*q9a*s11*s13*s16**2*s22*s24*s25**2*s34*s36**2 + 1024*q9a*s11**2*s12**2*s23*s24*s25**2*s34*s36**2 +  \
    512*q9a*s11*s12*s16**2*s23*s24*s25**2*s34*s36**2 + 512*q9a*s11**2*s12*s14*s24**2*s25**2*s34*s36**2 + 512*q9a*s11*s12*s15*s16*s24**2*s25**2*s34*s36**2 - 1024*q9a*s11*s12*s14*s15*s22*s25**3*s34*s36**2 +  \
    512*q9a*s11*s12*s13*s16*s22*s25**3*s34*s36**2 + 1024*q9a*s11*s14**2*s16*s22*s25**3*s34*s36**2 - 512*q9a*s11*s12**2*s16*s23*s25**3*s34*s36**2 - 512*q9a*s11*s12**2*s15*s24*s25**3*s34*s36**2 -  \
    512*q9a*s11*s12*s14*s16*s24*s25**3*s34*s36**2 + 512*q9a*s11*s12**2*s14*s25**4*s34*s36**2 - 1536*q9a*s11**2*s13*s15*s22**2*s23*s26*s34*s36**2 + 1536*q9a*s11**2*s12*s15*s22*s23**2*s26*s34*s36**2 -  \
    2560*q9a*s11**2*s14*s16*s22*s23**2*s26*s34*s36**2 + 1024*q9a*s11**2*s13*s16*s22*s23*s24*s26*s34*s36**2 + 512*q9a*s11**2*s12*s16*s23**2*s24*s26*s34*s36**2 + 512*q9a*s11**2*s12*s15*s23*s24**2*s26*s34*s36**2 +  \
    1024*q9a*s11**2*s13**2*s22**2*s25*s26*s34*s36**2 + 512*q9a*s11*s13*s15**2*s22**2*s25*s26*s34*s36**2 + 512*q9a*s11**2*s14**2*s22*s23*s25*s26*s34*s36**2 - 512*q9a*s11*s12*s15**2*s22*s23*s25*s26*s34*s36**2 +  \
    1536*q9a*s11*s14*s15*s16*s22*s23*s25*s26*s34*s36**2 - 1024*q9a*s11**2*s12**2*s23**2*s25*s26*s34*s36**2 - 2048*q9a*s11**2*s13*s14*s22*s24*s25*s26*s34*s36**2 + 512*q9a*s11*s14*s15**2*s22*s24*s25*s26*s34*s36**2 -  \
    512*q9a*s11*s13*s15*s16*s22*s24*s25*s26*s34*s36**2 - 2048*q9a*s11**2*s12*s14*s23*s24*s25*s26*s34*s36**2 + 3072*q9a*s11**2*s12*s13*s24**2*s25*s26*s34*s36**2 - 1024*q9a*s11*s12*s15**2*s24**2*s25*s26*s34*s36**2 -  \
    512*q9a*s11*s12*s13*s15*s22*s25**2*s26*s34*s36**2 - 512*q9a*s11*s14**2*s15*s22*s25**2*s26*s34*s36**2 + 512*q9a*s11*s12**2*s15*s23*s25**2*s26*s34*s36**2 + 512*q9a*s11*s12*s14*s16*s23*s25**2*s26*s34*s36**2 +  \
    1536*q9a*s11*s12*s14*s15*s24*s25**2*s26*s34*s36**2 - 1536*q9a*s11*s12*s13*s16*s24*s25**2*s26*s34*s36**2 - 512*q9a*s11*s12*s14**2*s25**3*s26*s34*s36**2 + 1024*q9a*s11**2*s13*s14*s22*s23*s26**2*s34*s36**2 +  \
    1536*q9a*s11**2*s12*s14*s23**2*s26**2*s34*s36**2 - 2048*q9a*s11**2*s12*s13*s23*s24*s26**2*s34*s36**2 - 512*q9a*s11*s13*s14*s15*s22*s25*s26**2*s34*s36**2 - 1024*q9a*s11*s12*s14*s15*s23*s25*s26**2*s34*s36**2 +  \
    1024*q9a*s11*s12*s13*s15*s24*s25*s26**2*s34*s36**2 + 512*q9a*s11*s12*s13*s14*s25**2*s26**2*s34*s36**2 + 4096*q8a*s11**3*s14*s22**2*s23*s33*s34*s36**2 + 1024*q8a*s11**3*s13*s22**2*s24*s33*s34*s36**2 +  \
    512*q8a*s11**2*s15**2*s22**2*s24*s33*s34*s36**2 - 6144*q8a*s11**3*s12*s22*s23*s24*s33*s34*s36**2 + 1024*q8a*s11**2*s16**2*s22*s23*s24*s33*s34*s36**2 + 1024*q7a*s11**3*s22**2*s23*s24*s33*s34*s36**2 -  \
    1024*q8a*s11**3*s14*s22*s24**2*s33*s34*s36**2 - 512*q8a*s11**2*s15*s16*s22*s24**2*s33*s34*s36**2 + 1024*q8a*s11**3*s12*s24**3*s33*s34*s36**2 - 2560*q8a*s11**2*s14*s15*s22**2*s25*s33*s34*s36**2 -  \
    512*q8a*s11**2*s13*s16*s22**2*s25*s33*s34*s36**2 + 1024*q8a*s11**2*s12*s16*s22*s23*s25*s33*s34*s36**2 - 512*q7a*s11**2*s16*s22**2*s23*s25*s33*s34*s36**2 + 1024*q8a*s11**2*s12*s15*s22*s24*s25*s33*s34*s36**2 +  \
    3072*q8a*s11**2*s14*s16*s22*s24*s25*s33*s34*s36**2 - 512*q8a*s11*s15*s16**2*s22*s24*s25*s33*s34*s36**2 + 512*q7a*s11**2*s15*s22**2*s24*s25*s33*s34*s36**2 - 1536*q8a*s11**2*s12*s16*s24**2*s25*s33*s34*s36**2 -  \
    1024*q7a*s11**2*s16*s22*s24**2*s25*s33*s34*s36**2 + 512*q8a*s11*s12*s15*s16*s22*s25**2*s33*s34*s36**2 - 1024*q8a*s11*s14*s16**2*s22*s25**2*s33*s34*s36**2 + 2560*q7a*s11**2*s14*s22**2*s25**2*s33*s34*s36**2 -  \
    512*q7a*s11*s15*s16*s22**2*s25**2*s33*s34*s36**2 + 1024*q8a*s11**2*s12**2*s24*s25**2*s33*s34*s36**2 + 512*q8a*s11*s12*s16**2*s24*s25**2*s33*s34*s36**2 - 3072*q7a*s11**2*s12*s22*s24*s25**2*s33*s34*s36**2 +  \
    1024*q7a*s11*s16**2*s22*s24*s25**2*s33*s34*s36**2 - 512*q8a*s11*s12**2*s16*s25**3*s33*s34*s36**2 + 512*q7a*s11*s12*s16*s22*s25**3*s33*s34*s36**2 - 1536*q8a*s11**2*s13*s15*s22**2*s26*s33*s34*s36**2 +  \
    3072*q8a*s11**2*s12*s15*s22*s23*s26*s33*s34*s36**2 - 5120*q8a*s11**2*s14*s16*s22*s23*s26*s33*s34*s36**2 - 1536*q7a*s11**2*s15*s22**2*s23*s26*s33*s34*s36**2 + 1024*q8a*s11**2*s13*s16*s22*s24*s26*s33*s34*s36**2
v3_7= \
    1024*q8a*s11**2*s12*s16*s23*s24*s26*s33*s34*s36**2 + 1024*q7a*s11**2*s16*s22*s23*s24*s26*s33*s34*s36**2 + 512*q8a*s11**2*s12*s15*s24**2*s26*s33*s34*s36**2 + 512*q8a*s11**2*s14**2*s22*s25*s26*s33*s34*s36**2 -  \
    512*q8a*s11*s12*s15**2*s22*s25*s26*s33*s34*s36**2 + 1536*q8a*s11*s14*s15*s16*s22*s25*s26*s33*s34*s36**2 + 2048*q7a*s11**2*s13*s22**2*s25*s26*s33*s34*s36**2 + 512*q7a*s11*s15**2*s22**2*s25*s26*s33*s34*s36**2 -  \
    2048*q8a*s11**2*s12**2*s23*s25*s26*s33*s34*s36**2 - 2048*q8a*s11**2*s12*s14*s24*s25*s26*s33*s34*s36**2 - 2048*q7a*s11**2*s14*s22*s24*s25*s26*s33*s34*s36**2 - 512*q7a*s11*s15*s16*s22*s24*s25*s26*s33*s34*s36**2 +  \
    3072*q7a*s11**2*s12*s24**2*s25*s26*s33*s34*s36**2 + 512*q8a*s11*s12**2*s15*s25**2*s26*s33*s34*s36**2 + 512*q8a*s11*s12*s14*s16*s25**2*s26*s33*s34*s36**2 - 512*q7a*s11*s12*s15*s22*s25**2*s26*s33*s34*s36**2 -  \
    1536*q7a*s11*s12*s16*s24*s25**2*s26*s33*s34*s36**2 + 1024*q8a*s11**2*s13*s14*s22*s26**2*s33*s34*s36**2 + 3072*q8a*s11**2*s12*s14*s23*s26**2*s33*s34*s36**2 + 1024*q7a*s11**2*s14*s22*s23*s26**2*s33*s34*s36**2 -  \
    2048*q8a*s11**2*s12*s13*s24*s26**2*s33*s34*s36**2 - 2048*q7a*s11**2*s12*s23*s24*s26**2*s33*s34*s36**2 - 1024*q8a*s11*s12*s14*s15*s25*s26**2*s33*s34*s36**2 - 512*q7a*s11*s14*s15*s22*s25*s26**2*s33*s34*s36**2 +  \
    1024*q7a*s11*s12*s15*s24*s25*s26**2*s33*s34*s36**2 + 512*q7a*s11*s12*s14*s25**2*s26**2*s33*s34*s36**2 - 2048*q8a*s11**3*s13*s22**2*s23*s34**2*s36**2 - 1024*q8a*s11**2*s15**2*s22**2*s23*s34**2*s36**2 +  \
    3072*q8a*s11**3*s12*s22*s23**2*s34**2*s36**2 - 1536*q8a*s11**2*s16**2*s22*s23**2*s34**2*s36**2 - 1024*q7a*s11**3*s22**2*s23**2*s34**2*s36**2 + 1024*q8a*s11**3*s14*s22*s23*s24*s34**2*s36**2 +  \
    1024*q8a*s11**2*s15*s16*s22*s23*s24*s34**2*s36**2 - 1024*q8a*s11**3*s12*s23*s24**2*s34**2*s36**2 + 1024*q8a*s11**2*s13*s15*s22**2*s25*s34**2*s36**2 + 512*q8a*s11*s15**3*s22**2*s25*s34**2*s36**2 -  \
    2048*q8a*s11**2*s14*s16*s22*s23*s25*s34**2*s36**2 + 1024*q8a*s11*s15*s16**2*s22*s23*s25*s34**2*s36**2 + 1024*q7a*s11**2*s15*s22**2*s23*s25*s34**2*s36**2 - 512*q8a*s11**2*s14*s15*s22*s24*s25*s34**2*s36**2 -  \
    512*q8a*s11*s15**2*s16*s22*s24*s25*s34**2*s36**2 + 1024*q8a*s11**2*s12*s16*s23*s24*s25*s34**2*s36**2 + 512*q8a*s11**2*s12*s15*s24**2*s25*s34**2*s36**2 + 1024*q8a*s11**2*s12*s13*s22*s25**2*s34**2*s36**2 -  \
    1024*q8a*s11*s12*s15**2*s22*s25**2*s34**2*s36**2 + 1024*q8a*s11*s14*s15*s16*s22*s25**2*s34**2*s36**2 - 512*q8a*s11*s13*s16**2*s22*s25**2*s34**2*s36**2 - 2048*q7a*s11**2*s13*s22**2*s25**2*s34**2*s36**2 -  \
    512*q7a*s11*s15**2*s22**2*s25**2*s34**2*s36**2 - 1024*q8a*s11**2*s12**2*s23*s25**2*s34**2*s36**2 + 1024*q7a*s11**2*s12*s22*s23*s25**2*s34**2*s36**2 - 512*q7a*s11*s16**2*s22*s23*s25**2*s34**2*s36**2 -  \
    512*q8a*s11*s12*s15*s16*s24*s25**2*s34**2*s36**2 + 512*q7a*s11**2*s14*s22*s24*s25**2*s34**2*s36**2 + 512*q7a*s11*s15*s16*s22*s24*s25**2*s34**2*s36**2 - 512*q7a*s11**2*s12*s24**2*s25**2*s34**2*s36**2 +  \
    512*q8a*s11*s12**2*s15*s25**3*s34**2*s36**2 + 1024*q7a*s11*s12*s15*s22*s25**3*s34**2*s36**2 - 1024*q7a*s11*s14*s16*s22*s25**3*s34**2*s36**2 + 512*q7a*s11*s12*s16*s24*s25**3*s34**2*s36**2 - 512*q7a*s11*s12**2*s25**4*s34**2*s36**2 +  \
    1024*q8a*s11**2*s14*s15*s22*s23*s26*s34**2*s36**2 + 1024*q8a*s11**2*s13*s16*s22*s23*s26*s34**2*s36**2 + 1536*q8a*s11**2*s12*s16*s23**2*s26*s34**2*s36**2 + 512*q7a*s11**2*s16*s22*s23**2*s26*s34**2*s36**2 -  \
    2048*q8a*s11**2*s12*s15*s23*s24*s26*s34**2*s36**2 - 512*q8a*s11*s14*s15**2*s22*s25*s26*s34**2*s36**2 - 512*q8a*s11*s13*s15*s16*s22*s25*s26*s34**2*s36**2 + 1024*q8a*s11**2*s12*s14*s23*s25*s26*s34**2*s36**2 -  \
    1024*q8a*s11*s12*s15*s16*s23*s25*s26*s34**2*s36**2 - 512*q7a*s11*s15*s16*s22*s23*s25*s26*s34**2*s36**2 + 1024*q8a*s11*s12*s15**2*s24*s25*s26*s34**2*s36**2 - 512*q8a*s11*s12*s14*s15*s25**2*s26*s34**2*s36**2 +  \
    512*q8a*s11*s12*s13*s16*s25**2*s26*s34**2*s36**2 + 512*q7a*s11*s14*s15*s22*s25**2*s26*s34**2*s36**2 + 1024*q7a*s11*s13*s16*s22*s25**2*s26*s34**2*s36**2 + 512*q7a*s11*s12*s16*s23*s25**2*s26*s34**2*s36**2 -  \
    1024*q7a*s11*s12*s15*s24*s25**2*s26*s34**2*s36**2 + 512*q7a*s11*s12*s14*s25**3*s26*s34**2*s36**2 - 1024*q8a*s11**2*s12*s13*s23*s26**2*s34**2*s36**2 - 512*q7a*s11**2*s12*s23**2*s26**2*s34**2*s36**2 +  \
    512*q8a*s11*s12*s13*s15*s25*s26**2*s34**2*s36**2 + 512*q7a*s11*s12*s15*s23*s25*s26**2*s34**2*s36**2 - 1024*q7a*s11*s12*s13*s25**2*s26**2*s34**2*s36**2 + 1024*q9a*s11**2*s13*s15*s22**3*s23*s35*s36**2 -  \
    1024*q9a*s11**2*s12*s15*s22**2*s23**2*s35*s36**2 + 1024*q9a*s11**2*s14*s16*s22**2*s23**2*s35*s36**2 + 512*q9a*s11**2*s14*s15*s22**2*s23*s24*s35*s36**2 - 1536*q9a*s11**2*s13*s16*s22**2*s23*s24*s35*s36**2 +  \
    512*q9a*s11**2*s12*s16*s22*s23**2*s24*s35*s36**2 + 512*q9a*s11**2*s13*s15*s22**2*s24**2*s35*s36**2 - 1024*q9a*s11**2*s12*s15*s22*s23*s24**2*s35*s36**2 + 512*q9a*s11**2*s14*s16*s22*s23*s24**2*s35*s36**2 -  \
    512*q9a*s11**2*s14*s15*s22*s24**3*s35*s36**2 - 512*q9a*s11**2*s12*s16*s23*s24**3*s35*s36**2 + 512*q9a*s11**2*s12*s15*s24**4*s35*s36**2 + 2048*q9a*s11**2*s13**2*s22**3*s25*s35*s36**2 - 1024*q9a*s11*s13*s15**2*s22**3*s25*s35*s36**2 -  \
    5120*q9a*s11**2*s12*s13*s22**2*s23*s25*s35*s36**2 + 2560*q9a*s11**2*s14**2*s22**2*s23*s25*s35*s36**2 + 1024*q9a*s11*s12*s15**2*s22**2*s23*s25*s35*s36**2 - 2048*q9a*s11*s14*s15*s16*s22**2*s23*s25*s35*s36**2 +  \
    1024*q9a*s11*s13*s16**2*s22**2*s23*s25*s35*s36**2 + 3072*q9a*s11**2*s12**2*s22*s23**2*s25*s35*s36**2 - 1024*q9a*s11*s12*s16**2*s22*s23**2*s25*s35*s36**2 - 2560*q9a*s11**2*s13*s14*s22**2*s24*s25*s35*s36**2 +  \
    512*q9a*s11*s14*s15**2*s22**2*s24*s25*s35*s36**2 - 3072*q9a*s11**2*s12*s14*s22*s23*s24*s25*s35*s36**2 + 2048*q9a*s11*s12*s15*s16*s22*s23*s24*s25*s35*s36**2 - 512*q9a*s11*s14*s16**2*s22*s23*s24*s25*s35*s36**2 +  \
    2048*q9a*s11**2*s12*s13*s22*s24**2*s25*s35*s36**2 + 512*q9a*s11**2*s14**2*s22*s24**2*s25*s35*s36**2 - 512*q9a*s11*s12*s15**2*s22*s24**2*s25*s35*s36**2 + 512*q9a*s11*s14*s15*s16*s22*s24**2*s25*s35*s36**2 +  \
    1024*q9a*s11**2*s12**2*s23*s24**2*s25*s35*s36**2 + 512*q9a*s11*s12*s16**2*s23*s24**2*s25*s35*s36**2 - 512*q9a*s11**2*s12*s14*s24**3*s25*s35*s36**2 - 512*q9a*s11*s12*s15*s16*s24**3*s25*s35*s36**2 +  \
    2048*q9a*s11*s12*s13*s15*s22**2*s25**2*s35*s36**2 - 512*q9a*s11*s14**2*s15*s22**2*s25**2*s35*s36**2 + 512*q9a*s11*s13*s14*s16*s22**2*s25**2*s35*s36**2 - 2048*q9a*s11*s12**2*s15*s22*s23*s25**2*s35*s36**2 +  \
    1536*q9a*s11*s12*s14*s16*s22*s23*s25**2*s35*s36**2 - 512*q9a*s11*s12*s13*s16*s22*s24*s25**2*s35*s36**2 - 512*q9a*s11*s14**2*s16*s22*s24*s25**2*s35*s36**2 - 1536*q9a*s11*s12**2*s16*s23*s24*s25**2*s35*s36**2 +  \
    512*q9a*s11*s12**2*s15*s24**2*s25**2*s35*s36**2 + 512*q9a*s11*s12*s14*s16*s24**2*s25**2*s35*s36**2 - 1024*q9a*s11*s12**2*s13*s22*s25**3*s35*s36**2 + 512*q9a*s11*s12*s14**2*s22*s25**3*s35*s36**2 +  \
    1024*q9a*s11*s12**3*s23*s25**3*s35*s36**2 - 512*q9a*s11*s12**2*s14*s24*s25**3*s35*s36**2 - 512*q9a*s11**2*s13*s14*s22**2*s23*s26*s35*s36**2 - 1024*q9a*s11*s14*s15**2*s22**2*s23*s26*s35*s36**2 +  \
    1024*q9a*s11*s13*s15*s16*s22**2*s23*s26*s35*s36**2 - 512*q9a*s11**2*s12*s14*s22*s23**2*s26*s35*s36**2 - 1024*q9a*s11*s12*s15*s16*s22*s23**2*s26*s35*s36**2 - 512*q9a*s11*s13*s15**2*s22**2*s24*s26*s35*s36**2 +  \
    2048*q9a*s11**2*s12*s13*s22*s23*s24*s26*s35*s36**2 - 512*q9a*s11**2*s14**2*s22*s23*s24*s26*s35*s36**2 + 1536*q9a*s11*s12*s15**2*s22*s23*s24*s26*s35*s36**2 - 512*q9a*s11*s14*s15*s16*s22*s23*s24*s26*s35*s36**2 -  \
    1024*q9a*s11**2*s12**2*s23**2*s24*s26*s35*s36**2 + 512*q9a*s11*s14*s15**2*s22*s24**2*s26*s35*s36**2 + 512*q9a*s11**2*s12*s14*s23*s24**2*s26*s35*s36**2 + 512*q9a*s11*s12*s15*s16*s23*s24**2*s26*s35*s36**2 -  \
    512*q9a*s11*s12*s15**2*s24**3*s26*s35*s36**2 + 3072*q9a*s11*s13*s14*s15*s22**2*s25*s26*s35*s36**2 - 2560*q9a*s11*s13**2*s16*s22**2*s25*s26*s35*s36**2 + 1024*q9a*s11*s12*s14*s15*s22*s23*s25*s26*s35*s36**2 +  \
    2048*q9a*s11*s12*s13*s16*s22*s23*s25*s26*s35*s36**2 + 512*q9a*s11*s12**2*s16*s23**2*s25*s26*s35*s36**2 - 2048*q9a*s11*s12*s13*s15*s22*s24*s25*s26*s35*s36**2 - 1536*q9a*s11*s14**2*s15*s22*s24*s25*s26*s35*s36**2 +  \
    1536*q9a*s11*s13*s14*s16*s22*s24*s25*s26*s35*s36**2 - 2048*q9a*s11*s12**2*s15*s23*s24*s25*s26*s35*s36**2 + 1536*q9a*s11*s12*s14*s15*s24**2*s25*s26*s35*s36**2 - 1536*q9a*s11*s12*s13*s16*s24**2*s25*s26*s35*s36**2 -  \
    3584*q9a*s11*s12*s13*s14*s22*s25**2*s26*s35*s36**2 + 1024*q9a*s11*s14**3*s22*s25**2*s26*s35*s36**2 + 512*q9a*s11*s12**2*s14*s23*s25**2*s26*s35*s36**2 + 3072*q9a*s11*s12**2*s13*s24*s25**2*s26*s35*s36**2 -  \
    1024*q9a*s11*s12*s14**2*s24*s25**2*s26*s35*s36**2 + 512*q9a*s11*s13**2*s15*s22**2*s26**2*s35*s36**2 - 2048*q9a*s11*s12*s13*s15*s22*s23*s26**2*s35*s36**2 + 1024*q9a*s11*s14**2*s15*s22*s23*s26**2*s35*s36**2 +  \
    1536*q9a*s11*s12**2*s15*s23**2*s26**2*s35*s36**2 - 512*q9a*s11*s13*s14*s15*s22*s24*s26**2*s35*s36**2 - 1024*q9a*s11*s12*s14*s15*s23*s24*s26**2*s35*s36**2 + 512*q9a*s11*s12*s13*s15*s24**2*s26**2*s35*s36**2 +  \
    2048*q9a*s11*s12*s13**2*s22*s25*s26**2*s35*s36**2 - 1024*q9a*s11*s13*s14**2*s22*s25*s26**2*s35*s36**2 - 2048*q9a*s11*s12**2*s13*s23*s25*s26**2*s35*s36**2 + 1024*q9a*s11*s12*s13*s14*s24*s25*s26**2*s35*s36**2 +  \
    1024*q8a*s11**2*s13*s15*s22**3*s33*s35*s36**2 - 2048*q8a*s11**2*s12*s15*s22**2*s23*s33*s35*s36**2 + 2048*q8a*s11**2*s14*s16*s22**2*s23*s33*s35*s36**2 + 1024*q7a*s11**2*s15*s22**3*s23*s33*s35*s36**2 +  \
    512*q8a*s11**2*s14*s15*s22**2*s24*s33*s35*s36**2 - 1536*q8a*s11**2*s13*s16*s22**2*s24*s33*s35*s36**2 + 1024*q8a*s11**2*s12*s16*s22*s23*s24*s33*s35*s36**2 - 1536*q7a*s11**2*s16*s22**2*s23*s24*s33*s35*s36**2 -  \
    1024*q8a*s11**2*s12*s15*s22*s24**2*s33*s35*s36**2 + 512*q8a*s11**2*s14*s16*s22*s24**2*s33*s35*s36**2 + 512*q7a*s11**2*s15*s22**2*s24**2*s33*s35*s36**2 - 512*q8a*s11**2*s12*s16*s24**3*s33*s35*s36**2 -  \
    5120*q8a*s11**2*s12*s13*s22**2*s25*s33*s35*s36**2 + 2560*q8a*s11**2*s14**2*s22**2*s25*s33*s35*s36**2 + 1024*q8a*s11*s12*s15**2*s22**2*s25*s33*s35*s36**2 - 2048*q8a*s11*s14*s15*s16*s22**2*s25*s33*s35*s36**2 +  \
    1024*q8a*s11*s13*s16**2*s22**2*s25*s33*s35*s36**2 + 4096*q7a*s11**2*s13*s22**3*s25*s33*s35*s36**2 - 1024*q7a*s11*s15**2*s22**3*s25*s33*s35*s36**2 + 6144*q8a*s11**2*s12**2*s22*s23*s25*s33*s35*s36**2 -  \
    2048*q8a*s11*s12*s16**2*s22*s23*s25*s33*s35*s36**2 - 5120*q7a*s11**2*s12*s22**2*s23*s25*s33*s35*s36**2 + 1024*q7a*s11*s16**2*s22**2*s23*s25*s33*s35*s36**2 - 3072*q8a*s11**2*s12*s14*s22*s24*s25*s33*s35*s36**2 +  \
    2048*q8a*s11*s12*s15*s16*s22*s24*s25*s33*s35*s36**2 - 512*q8a*s11*s14*s16**2*s22*s24*s25*s33*s35*s36**2 - 2560*q7a*s11**2*s14*s22**2*s24*s25*s33*s35*s36**2 + 1024*q8a*s11**2*s12**2*s24**2*s25*s33*s35*s36**2 +  \
    512*q8a*s11*s12*s16**2*s24**2*s25*s33*s35*s36**2 + 2048*q7a*s11**2*s12*s22*s24**2*s25*s33*s35*s36**2 - 2048*q8a*s11*s12**2*s15*s22*s25**2*s33*s35*s36**2 + 1536*q8a*s11*s12*s14*s16*s22*s25**2*s33*s35*s36**2 +  \
    2048*q7a*s11*s12*s15*s22**2*s25**2*s33*s35*s36**2 + 512*q7a*s11*s14*s16*s22**2*s25**2*s33*s35*s36**2 - 1536*q8a*s11*s12**2*s16*s24*s25**2*s33*s35*s36**2 - 512*q7a*s11*s12*s16*s22*s24*s25**2*s33*s35*s36**2 +  \
    1024*q8a*s11*s12**3*s25**3*s33*s35*s36**2 - 1024*q7a*s11*s12**2*s22*s25**3*s33*s35*s36**2 - 512*q8a*s11**2*s13*s14*s22**2*s26*s33*s35*s36**2 - 1024*q8a*s11*s14*s15**2*s22**2*s26*s33*s35*s36**2 +  \
    1024*q8a*s11*s13*s15*s16*s22**2*s26*s33*s35*s36**2 - 1024*q8a*s11**2*s12*s14*s22*s23*s26*s33*s35*s36**2 - 2048*q8a*s11*s12*s15*s16*s22*s23*s26*s33*s35*s36**2 - 512*q7a*s11**2*s14*s22**2*s23*s26*s33*s35*s36**2 +  \
    1024*q7a*s11*s15*s16*s22**2*s23*s26*s33*s35*s36**2 + 2048*q8a*s11**2*s12*s13*s22*s24*s26*s33*s35*s36**2 - 512*q8a*s11**2*s14**2*s22*s24*s26*s33*s35*s36**2 + 1536*q8a*s11*s12*s15**2*s22*s24*s26*s33*s35*s36**2 -  \
    512*q8a*s11*s14*s15*s16*s22*s24*s26*s33*s35*s36**2 - 512*q7a*s11*s15**2*s22**2*s24*s26*s33*s35*s36**2 - 2048*q8a*s11**2*s12**2*s23*s24*s26*s33*s35*s36**2 + 2048*q7a*s11**2*s12*s22*s23*s24*s26*s33*s35*s36**2 +  \
    512*q8a*s11**2*s12*s14*s24**2*s26*s33*s35*s36**2 + 512*q8a*s11*s12*s15*s16*s24**2*s26*s33*s35*s36**2 + 1024*q8a*s11*s12*s14*s15*s22*s25*s26*s33*s35*s36**2 + 2048*q8a*s11*s12*s13*s16*s22*s25*s26*s33*s35*s36**2 +  \
    3072*q7a*s11*s14*s15*s22**2*s25*s26*s33*s35*s36**2 - 5120*q7a*s11*s13*s16*s22**2*s25*s26*s33*s35*s36**2 + 1024*q8a*s11*s12**2*s16*s23*s25*s26*s33*s35*s36**2 + 2048*q7a*s11*s12*s16*s22*s23*s25*s26*s33*s35*s36**2 -  \
    2048*q8a*s11*s12**2*s15*s24*s25*s26*s33*s35*s36**2 - 2048*q7a*s11*s12*s15*s22*s24*s25*s26*s33*s35*s36**2 + 1536*q7a*s11*s14*s16*s22*s24*s25*s26*s33*s35*s36**2 - 1536*q7a*s11*s12*s16*s24**2*s25*s26*s33*s35*s36**2 +  \
    512*q8a*s11*s12**2*s14*s25**2*s26*s33*s35*s36**2 - 3584*q7a*s11*s12*s14*s22*s25**2*s26*s33*s35*s36**2 + 3072*q7a*s11*s12**2*s24*s25**2*s26*s33*s35*s36**2 - 2048*q8a*s11*s12*s13*s15*s22*s26**2*s33*s35*s36**2
v3_8= \
    1024*q8a*s11*s14**2*s15*s22*s26**2*s33*s35*s36**2 + 1024*q7a*s11*s13*s15*s22**2*s26**2*s33*s35*s36**2 + 3072*q8a*s11*s12**2*s15*s23*s26**2*s33*s35*s36**2 - 2048*q7a*s11*s12*s15*s22*s23*s26**2*s33*s35*s36**2 -  \
    1024*q8a*s11*s12*s14*s15*s24*s26**2*s33*s35*s36**2 - 512*q7a*s11*s14*s15*s22*s24*s26**2*s33*s35*s36**2 + 512*q7a*s11*s12*s15*s24**2*s26**2*s33*s35*s36**2 - 2048*q8a*s11*s12**2*s13*s25*s26**2*s33*s35*s36**2 +  \
    4096*q7a*s11*s12*s13*s22*s25*s26**2*s33*s35*s36**2 - 1024*q7a*s11*s14**2*s22*s25*s26**2*s33*s35*s36**2 - 2048*q7a*s11*s12**2*s23*s25*s26**2*s33*s35*s36**2 + 1024*q7a*s11*s12*s14*s24*s25*s26**2*s33*s35*s36**2 +  \
    2048*q8a*s11**2*s14*s15*s22**2*s23*s34*s35*s36**2 - 2048*q8a*s11**2*s13*s16*s22**2*s23*s34*s35*s36**2 + 3072*q8a*s11**2*s12*s16*s22*s23**2*s34*s35*s36**2 - 1024*q7a*s11**2*s16*s22**2*s23**2*s34*s35*s36**2 -  \
    1536*q8a*s11**2*s13*s15*s22**2*s24*s34*s35*s36**2 + 1024*q8a*s11**2*s12*s15*s22*s23*s24*s34*s35*s36**2 - 3072*q8a*s11**2*s14*s16*s22*s23*s24*s34*s35*s36**2 - 1536*q7a*s11**2*s15*s22**2*s23*s24*s34*s35*s36**2 +  \
    512*q8a*s11**2*s14*s15*s22*s24**2*s34*s35*s36**2 + 1024*q8a*s11**2*s13*s16*s22*s24**2*s34*s35*s36**2 + 1024*q8a*s11**2*s12*s16*s23*s24**2*s34*s35*s36**2 + 1024*q7a*s11**2*s16*s22*s23*s24**2*s34*s35*s36**2 -  \
    512*q8a*s11**2*s12*s15*s24**3*s34*s35*s36**2 - 512*q8a*s11**2*s13*s14*s22**2*s25*s34*s35*s36**2 - 1024*q8a*s11*s14*s15**2*s22**2*s25*s34*s35*s36**2 + 2048*q8a*s11*s13*s15*s16*s22**2*s25*s34*s35*s36**2 -  \
    1024*q8a*s11**2*s12*s14*s22*s23*s25*s34*s35*s36**2 - 4096*q8a*s11*s12*s15*s16*s22*s23*s25*s34*s35*s36**2 + 2048*q8a*s11*s14*s16**2*s22*s23*s25*s34*s35*s36**2 - 512*q7a*s11**2*s14*s22**2*s23*s25*s34*s35*s36**2 +  \
    2048*q7a*s11*s15*s16*s22**2*s23*s25*s34*s35*s36**2 + 512*q8a*s11**2*s14**2*s22*s24*s25*s34*s35*s36**2 + 512*q8a*s11*s12*s15**2*s22*s24*s25*s34*s35*s36**2 - 512*q8a*s11*s13*s16**2*s22*s24*s25*s34*s35*s36**2 +  \
    2048*q7a*s11**2*s13*s22**2*s24*s25*s34*s35*s36**2 + 512*q7a*s11*s15**2*s22**2*s24*s25*s34*s35*s36**2 - 1024*q8a*s11*s12*s16**2*s23*s24*s25*s34*s35*s36**2 - 512*q7a*s11*s16**2*s22*s23*s24*s25*s34*s35*s36**2 -  \
    512*q8a*s11**2*s12*s14*s24**2*s25*s34*s35*s36**2 + 512*q8a*s11*s12*s15*s16*s24**2*s25*s34*s35*s36**2 - 1024*q7a*s11**2*s14*s22*s24**2*s25*s34*s35*s36**2 - 512*q7a*s11*s15*s16*s22*s24**2*s25*s34*s35*s36**2 +  \
    1024*q7a*s11**2*s12*s24**3*s25*s34*s35*s36**2 + 1536*q8a*s11*s12*s14*s15*s22*s25**2*s34*s35*s36**2 + 1024*q8a*s11*s12*s13*s16*s22*s25**2*s34*s35*s36**2 - 1024*q8a*s11*s14**2*s16*s22*s25**2*s34*s35*s36**2 +  \
    512*q7a*s11*s14*s15*s22**2*s25**2*s34*s35*s36**2 - 3072*q7a*s11*s13*s16*s22**2*s25**2*s34*s35*s36**2 + 1024*q8a*s11*s12**2*s16*s23*s25**2*s34*s35*s36**2 + 1024*q7a*s11*s12*s16*s22*s23*s25**2*s34*s35*s36**2 -  \
    512*q8a*s11*s12**2*s15*s24*s25**2*s34*s35*s36**2 + 512*q8a*s11*s12*s14*s16*s24*s25**2*s34*s35*s36**2 - 1536*q7a*s11*s12*s15*s22*s24*s25**2*s34*s35*s36**2 + 1536*q7a*s11*s14*s16*s22*s24*s25**2*s34*s35*s36**2 -  \
    1024*q7a*s11*s12*s16*s24**2*s25**2*s34*s35*s36**2 - 512*q8a*s11*s12**2*s14*s25**3*s34*s35*s36**2 - 512*q7a*s11*s12*s14*s22*s25**3*s34*s35*s36**2 + 1024*q7a*s11*s12**2*s24*s25**3*s34*s35*s36**2 -  \
    1024*q8a*s11**2*s13**2*s22**2*s26*s34*s35*s36**2 + 1024*q8a*s11*s13*s15**2*s22**2*s26*s34*s35*s36**2 + 6144*q8a*s11**2*s12*s13*s22*s23*s26*s34*s35*s36**2 - 3072*q8a*s11**2*s14**2*s22*s23*s26*s34*s35*s36**2 -  \
    2048*q8a*s11*s12*s15**2*s22*s23*s26*s34*s35*s36**2 + 2048*q8a*s11*s14*s15*s16*s22*s23*s26*s34*s35*s36**2 - 2048*q7a*s11**2*s13*s22**2*s23*s26*s34*s35*s36**2 + 1024*q7a*s11*s15**2*s22**2*s23*s26*s34*s35*s36**2 -  \
    6144*q8a*s11**2*s12**2*s23**2*s26*s34*s35*s36**2 + 3072*q7a*s11**2*s12*s22*s23**2*s26*s34*s35*s36**2 + 2048*q8a*s11**2*s13*s14*s22*s24*s26*s34*s35*s36**2 - 512*q8a*s11*s14*s15**2*s22*s24*s26*s34*s35*s36**2 -  \
    512*q8a*s11*s13*s15*s16*s22*s24*s26*s34*s35*s36**2 + 5120*q8a*s11**2*s12*s14*s23*s24*s26*s34*s35*s36**2 - 1024*q8a*s11*s12*s15*s16*s23*s24*s26*s34*s35*s36**2 + 2048*q7a*s11**2*s14*s22*s23*s24*s26*s34*s35*s36**2 -  \
    512*q7a*s11*s15*s16*s22*s23*s24*s26*s34*s35*s36**2 - 3072*q8a*s11**2*s12*s13*s24**2*s26*s34*s35*s36**2 + 512*q8a*s11*s12*s15**2*s24**2*s26*s34*s35*s36**2 - 3072*q7a*s11**2*s12*s23*s24**2*s26*s34*s35*s36**2 -  \
    1024*q8a*s11*s12*s13*s15*s22*s25*s26*s34*s35*s36**2 + 1536*q8a*s11*s14**2*s15*s22*s25*s26*s34*s35*s36**2 - 2560*q8a*s11*s13*s14*s16*s22*s25*s26*s34*s35*s36**2 - 3072*q7a*s11*s13*s15*s22**2*s25*s26*s34*s35*s36**2 +  \
    5120*q8a*s11*s12**2*s15*s23*s25*s26*s34*s35*s36**2 - 1024*q8a*s11*s12*s14*s16*s23*s25*s26*s34*s35*s36**2 - 1024*q7a*s11*s12*s15*s22*s23*s25*s26*s34*s35*s36**2 - 2560*q7a*s11*s14*s16*s22*s23*s25*s26*s34*s35*s36**2 -  \
    2048*q8a*s11*s12*s14*s15*s24*s25*s26*s34*s35*s36**2 + 2048*q8a*s11*s12*s13*s16*s24*s25*s26*s34*s35*s36**2 + 2048*q7a*s11*s13*s16*s22*s24*s25*s26*s34*s35*s36**2 + 2048*q7a*s11*s12*s16*s23*s24*s25*s26*s34*s35*s36**2 +  \
    512*q7a*s11*s12*s15*s24**2*s25*s26*s34*s35*s36**2 - 3072*q8a*s11*s12**2*s13*s25**2*s26*s34*s35*s36**2 + 512*q8a*s11*s12*s14**2*s25**2*s26*s34*s35*s36**2 + 6144*q7a*s11*s12*s13*s22*s25**2*s26*s34*s35*s36**2 -  \
    1024*q7a*s11*s14**2*s22*s25**2*s26*s34*s35*s36**2 - 3072*q7a*s11*s12**2*s23*s25**2*s26*s34*s35*s36**2 + 512*q7a*s11*s12*s14*s24*s25**2*s26*s34*s35*s36**2 - 512*q8a*s11*s13*s14*s15*s22*s26**2*s34*s35*s36**2 -  \
    1024*q8a*s11*s12*s14*s15*s23*s26**2*s34*s35*s36**2 - 512*q7a*s11*s14*s15*s22*s23*s26**2*s34*s35*s36**2 + 1024*q8a*s11*s12*s13*s15*s24*s26**2*s34*s35*s36**2 + 1024*q7a*s11*s12*s15*s23*s24*s26**2*s34*s35*s36**2 +  \
    1024*q8a*s11*s12*s13*s14*s25*s26**2*s34*s35*s36**2 + 2048*q7a*s11*s13*s14*s22*s25*s26**2*s34*s35*s36**2 + 1024*q7a*s11*s12*s14*s23*s25*s26**2*s34*s35*s36**2 - 4096*q7a*s11*s12*s13*s24*s25*s26**2*s34*s35*s36**2 -  \
    1024*q8a*s11**2*s13**2*s22**3*s35**2*s36**2 + 4096*q8a*s11**2*s12*s13*s22**2*s23*s35**2*s36**2 - 2048*q8a*s11**2*s14**2*s22**2*s23*s35**2*s36**2 - 2048*q7a*s11**2*s13*s22**3*s23*s35**2*s36**2 -  \
    3072*q8a*s11**2*s12**2*s22*s23**2*s35**2*s36**2 + 2048*q7a*s11**2*s12*s22**2*s23**2*s35**2*s36**2 + 1024*q8a*s11**2*s13*s14*s22**2*s24*s35**2*s36**2 + 2048*q8a*s11**2*s12*s14*s22*s23*s24*s35**2*s36**2 +  \
    1024*q7a*s11**2*s14*s22**2*s23*s24*s35**2*s36**2 - 512*q8a*s11**2*s14**2*s22*s24**2*s35**2*s36**2 - 1024*q7a*s11**2*s13*s22**2*s24**2*s35**2*s36**2 - 1024*q8a*s11**2*s12**2*s23*s24**2*s35**2*s36**2 +  \
    512*q8a*s11**2*s12*s14*s24**3*s35**2*s36**2 + 512*q7a*s11**2*s14*s22*s24**3*s35**2*s36**2 - 512*q7a*s11**2*s12*s24**4*s35**2*s36**2 - 1024*q8a*s11*s12*s13*s15*s22**2*s25*s35**2*s36**2 + 512*q8a*s11*s14**2*s15*s22**2*s25*s35**2*s36**2 +  \
    1024*q7a*s11*s13*s15*s22**3*s25*s35**2*s36**2 + 1024*q8a*s11*s12**2*s15*s22*s23*s25*s35**2*s36**2 - 1024*q7a*s11*s12*s15*s22**2*s23*s25*s35**2*s36**2 - 512*q8a*s11*s12*s14*s15*s22*s24*s25*s35**2*s36**2 -  \
    1024*q8a*s11*s12*s13*s16*s22*s24*s25*s35**2*s36**2 + 512*q8a*s11*s14**2*s16*s22*s24*s25*s35**2*s36**2 - 512*q7a*s11*s14*s15*s22**2*s24*s25*s35**2*s36**2 + 1024*q7a*s11*s13*s16*s22**2*s24*s25*s35**2*s36**2 +  \
    1024*q8a*s11*s12**2*s16*s23*s24*s25*s35**2*s36**2 - 1024*q7a*s11*s12*s16*s22*s23*s24*s25*s35**2*s36**2 - 512*q8a*s11*s12*s14*s16*s24**2*s25*s35**2*s36**2 + 512*q7a*s11*s12*s15*s22*s24**2*s25*s35**2*s36**2 -  \
    512*q7a*s11*s14*s16*s22*s24**2*s25*s35**2*s36**2 + 512*q7a*s11*s12*s16*s24**3*s25*s35**2*s36**2 + 1024*q8a*s11*s12**2*s13*s22*s25**2*s35**2*s36**2 - 512*q8a*s11*s12*s14**2*s22*s25**2*s35**2*s36**2 -  \
    1024*q7a*s11*s12*s13*s22**2*s25**2*s35**2*s36**2 - 1024*q8a*s11*s12**3*s23*s25**2*s35**2*s36**2 + 1024*q7a*s11*s12**2*s22*s23*s25**2*s35**2*s36**2 + 512*q8a*s11*s12**2*s14*s24*s25**2*s35**2*s36**2 +  \
    512*q7a*s11*s12*s14*s22*s24*s25**2*s35**2*s36**2 - 512*q7a*s11*s12**2*s24**2*s25**2*s35**2*s36**2 + 512*q8a*s11*s13**2*s16*s22**2*s26*s35**2*s36**2 - 2048*q8a*s11*s12*s13*s16*s22*s23*s26*s35**2*s36**2 +  \
    1024*q8a*s11*s14**2*s16*s22*s23*s26*s35**2*s36**2 + 1024*q7a*s11*s13*s16*s22**2*s23*s26*s35**2*s36**2 + 1536*q8a*s11*s12**2*s16*s23**2*s26*s35**2*s36**2 - 1024*q7a*s11*s12*s16*s22*s23**2*s26*s35**2*s36**2 -  \
    1024*q8a*s11*s12*s13*s15*s22*s24*s26*s35**2*s36**2 + 512*q8a*s11*s14**2*s15*s22*s24*s26*s35**2*s36**2 - 512*q8a*s11*s13*s14*s16*s22*s24*s26*s35**2*s36**2 + 1024*q7a*s11*s13*s15*s22**2*s24*s26*s35**2*s36**2 +  \
    1024*q8a*s11*s12**2*s15*s23*s24*s26*s35**2*s36**2 - 1024*q8a*s11*s12*s14*s16*s23*s24*s26*s35**2*s36**2 - 1024*q7a*s11*s12*s15*s22*s23*s24*s26*s35**2*s36**2 - 512*q7a*s11*s14*s16*s22*s23*s24*s26*s35**2*s36**2 -  \
    512*q8a*s11*s12*s14*s15*s24**2*s26*s35**2*s36**2 + 512*q8a*s11*s12*s13*s16*s24**2*s26*s35**2*s36**2 - 512*q7a*s11*s14*s15*s22*s24**2*s26*s35**2*s36**2 + 512*q7a*s11*s12*s16*s23*s24**2*s26*s35**2*s36**2 +  \
    512*q7a*s11*s12*s15*s24**3*s26*s35**2*s36**2 + 2048*q8a*s11*s12*s13*s14*s22*s25*s26*s35**2*s36**2 - 1024*q8a*s11*s14**3*s22*s25*s26*s35**2*s36**2 - 2048*q7a*s11*s13*s14*s22**2*s25*s26*s35**2*s36**2 -  \
    2048*q8a*s11*s12**2*s14*s23*s25*s26*s35**2*s36**2 + 2048*q7a*s11*s12*s14*s22*s23*s25*s26*s35**2*s36**2 + 1024*q8a*s11*s12*s14**2*s24*s25*s26*s35**2*s36**2 + 1024*q7a*s11*s14**2*s22*s24*s25*s26*s35**2*s36**2 -  \
    1024*q7a*s11*s12*s14*s24**2*s25*s26*s35**2*s36**2 + 1024*q8a*s11*s12*s13**2*s22*s26**2*s35**2*s36**2 - 512*q8a*s11*s13*s14**2*s22*s26**2*s35**2*s36**2 - 1536*q7a*s11*s13**2*s22**2*s26**2*s35**2*s36**2 -  \
    1024*q8a*s11*s12**2*s13*s23*s26**2*s35**2*s36**2 + 2048*q7a*s11*s12*s13*s22*s23*s26**2*s35**2*s36**2 - 512*q7a*s11*s14**2*s22*s23*s26**2*s35**2*s36**2 - 512*q7a*s11*s12**2*s23**2*s26**2*s35**2*s36**2 +  \
    512*q8a*s11*s12*s13*s14*s24*s26**2*s35**2*s36**2 + 1024*q7a*s11*s13*s14*s22*s24*s26**2*s35**2*s36**2 + 512*q7a*s11*s12*s14*s23*s24*s26**2*s35**2*s36**2 - 1024*q7a*s11*s12*s13*s24**2*s26**2*s35**2*s36**2 -  \
    1024*q9a*s11**2*s14*s15*s22**2*s23**2*s36**3 + 1024*q9a*s11**2*s13*s16*s22**2*s23**2*s36**3 - 1024*q9a*s11**2*s12*s16*s22*s23**3*s36**3 - 512*q9a*s11**2*s13*s15*s22**2*s23*s24*s36**3 + 1536*q9a*s11**2*s12*s15*s22*s23**2*s24*s36**3 -  \
    512*q9a*s11**2*s14*s16*s22*s23**2*s24*s36**3 + 512*q9a*s11**2*s14*s15*s22*s23*s24**2*s36**3 + 512*q9a*s11**2*s12*s16*s23**2*s24**2*s36**3 - 512*q9a*s11**2*s12*s15*s23*s24**3*s36**3 + 512*q9a*s11**2*s13*s14*s22**2*s23*s25*s36**3 +  \
    1024*q9a*s11*s14*s15**2*s22**2*s23*s25*s36**3 - 1024*q9a*s11*s13*s15*s16*s22**2*s23*s25*s36**3 + 512*q9a*s11**2*s12*s14*s22*s23**2*s25*s36**3 + 1024*q9a*s11*s12*s15*s16*s22*s23**2*s25*s36**3 -  \
    1024*q9a*s11**2*s13**2*s22**2*s24*s25*s36**3 + 512*q9a*s11*s13*s15**2*s22**2*s24*s25*s36**3 + 2048*q9a*s11**2*s12*s13*s22*s23*s24*s25*s36**3 - 1536*q9a*s11**2*s14**2*s22*s23*s24*s25*s36**3 -  \
    1536*q9a*s11*s12*s15**2*s22*s23*s24*s25*s36**3 + 512*q9a*s11*s14*s15*s16*s22*s23*s24*s25*s36**3 - 2048*q9a*s11**2*s12**2*s23**2*s24*s25*s36**3 + 1024*q9a*s11**2*s13*s14*s22*s24**2*s25*s36**3 -  \
    512*q9a*s11*s14*s15**2*s22*s24**2*s25*s36**3 + 1536*q9a*s11**2*s12*s14*s23*s24**2*s25*s36**3 - 512*q9a*s11*s12*s15*s16*s23*s24**2*s25*s36**3 - 1024*q9a*s11**2*s12*s13*s24**3*s25*s36**3 + 512*q9a*s11*s12*s15**2*s24**3*s25*s36**3 -  \
    1536*q9a*s11*s13*s14*s15*s22**2*s25**2*s36**3 + 1536*q9a*s11*s13**2*s16*s22**2*s25**2*s36**3 - 512*q9a*s11*s12*s14*s15*s22*s23*s25**2*s36**3 - 2048*q9a*s11*s12*s13*s16*s22*s23*s25**2*s36**3 +  \
    512*q9a*s11*s14**2*s16*s22*s23*s25**2*s36**3 + 512*q9a*s11*s12**2*s16*s23**2*s25**2*s36**3 + 512*q9a*s11*s12*s13*s15*s22*s24*s25**2*s36**3 + 1024*q9a*s11*s14**2*s15*s22*s24*s25**2*s36**3 -  \
    1024*q9a*s11*s13*s14*s16*s22*s24*s25**2*s36**3 + 1536*q9a*s11*s12**2*s15*s23*s24*s25**2*s36**3 - 512*q9a*s11*s12*s14*s16*s23*s24*s25**2*s36**3 - 1024*q9a*s11*s12*s14*s15*s24**2*s25**2*s36**3 +  \
    1024*q9a*s11*s12*s13*s16*s24**2*s25**2*s36**3 + 1536*q9a*s11*s12*s13*s14*s22*s25**3*s36**3 - 512*q9a*s11*s14**3*s22*s25**3*s36**3 - 512*q9a*s11*s12**2*s14*s23*s25**3*s36**3 - 1024*q9a*s11*s12**2*s13*s24*s25**3*s36**3 +  \
    512*q9a*s11*s12*s14**2*s24*s25**3*s36**3 + 1024*q9a*s11**2*s13**2*s22**2*s23*s26*s36**3 - 3072*q9a*s11**2*s12*s13*s22*s23**2*s26*s36**3 + 1536*q9a*s11**2*s14**2*s22*s23**2*s26*s36**3 + 2048*q9a*s11**2*s12**2*s23**3*s26*s36**3 -  \
    1024*q9a*s11**2*s13*s14*s22*s23*s24*s26*s36**3 - 1536*q9a*s11**2*s12*s14*s23**2*s24*s26*s36**3 + 1024*q9a*s11**2*s12*s13*s23*s24**2*s26*s36**3 - 512*q9a*s11*s13**2*s15*s22**2*s25*s26*s36**3 +  \
    2048*q9a*s11*s12*s13*s15*s22*s23*s25*s26*s36**3 - 1024*q9a*s11*s14**2*s15*s22*s23*s25*s26*s36**3 - 1536*q9a*s11*s12**2*s15*s23**2*s25*s26*s36**3 + 512*q9a*s11*s13*s14*s15*s22*s24*s25*s26*s36**3
v3_9= \
    1024*q9a*s11*s12*s14*s15*s23*s24*s25*s26*s36**3 - 512*q9a*s11*s12*s13*s15*s24**2*s25*s26*s36**3 - 1024*q9a*s11*s12*s13**2*s22*s25**2*s26*s36**3 + 512*q9a*s11*s13*s14**2*s22*s25**2*s26*s36**3 +  \
    1024*q9a*s11*s12**2*s13*s23*s25**2*s26*s36**3 - 512*q9a*s11*s12*s13*s14*s24*s25**2*s26*s36**3 - 2048*q8a*s11**2*s14*s15*s22**2*s23*s33*s36**3 + 2048*q8a*s11**2*s13*s16*s22**2*s23*s33*s36**3 -  \
    3072*q8a*s11**2*s12*s16*s22*s23**2*s33*s36**3 + 1024*q7a*s11**2*s16*s22**2*s23**2*s33*s36**3 - 512*q8a*s11**2*s13*s15*s22**2*s24*s33*s36**3 + 3072*q8a*s11**2*s12*s15*s22*s23*s24*s33*s36**3 -  \
    1024*q8a*s11**2*s14*s16*s22*s23*s24*s33*s36**3 - 512*q7a*s11**2*s15*s22**2*s23*s24*s33*s36**3 + 512*q8a*s11**2*s14*s15*s22*s24**2*s33*s36**3 + 1024*q8a*s11**2*s12*s16*s23*s24**2*s33*s36**3 -  \
    512*q8a*s11**2*s12*s15*s24**3*s33*s36**3 + 512*q8a*s11**2*s13*s14*s22**2*s25*s33*s36**3 + 1024*q8a*s11*s14*s15**2*s22**2*s25*s33*s36**3 - 1024*q8a*s11*s13*s15*s16*s22**2*s25*s33*s36**3 +  \
    1024*q8a*s11**2*s12*s14*s22*s23*s25*s33*s36**3 + 2048*q8a*s11*s12*s15*s16*s22*s23*s25*s33*s36**3 + 512*q7a*s11**2*s14*s22**2*s23*s25*s33*s36**3 - 1024*q7a*s11*s15*s16*s22**2*s23*s25*s33*s36**3 +  \
    2048*q8a*s11**2*s12*s13*s22*s24*s25*s33*s36**3 - 1536*q8a*s11**2*s14**2*s22*s24*s25*s33*s36**3 - 1536*q8a*s11*s12*s15**2*s22*s24*s25*s33*s36**3 + 512*q8a*s11*s14*s15*s16*s22*s24*s25*s33*s36**3 -  \
    2048*q7a*s11**2*s13*s22**2*s24*s25*s33*s36**3 + 512*q7a*s11*s15**2*s22**2*s24*s25*s33*s36**3 - 4096*q8a*s11**2*s12**2*s23*s24*s25*s33*s36**3 + 2048*q7a*s11**2*s12*s22*s23*s24*s25*s33*s36**3 +  \
    1536*q8a*s11**2*s12*s14*s24**2*s25*s33*s36**3 - 512*q8a*s11*s12*s15*s16*s24**2*s25*s33*s36**3 + 1024*q7a*s11**2*s14*s22*s24**2*s25*s33*s36**3 - 1024*q7a*s11**2*s12*s24**3*s25*s33*s36**3 -  \
    512*q8a*s11*s12*s14*s15*s22*s25**2*s33*s36**3 - 2048*q8a*s11*s12*s13*s16*s22*s25**2*s33*s36**3 + 512*q8a*s11*s14**2*s16*s22*s25**2*s33*s36**3 - 1536*q7a*s11*s14*s15*s22**2*s25**2*s33*s36**3 +  \
    3072*q7a*s11*s13*s16*s22**2*s25**2*s33*s36**3 + 1024*q8a*s11*s12**2*s16*s23*s25**2*s33*s36**3 - 2048*q7a*s11*s12*s16*s22*s23*s25**2*s33*s36**3 + 1536*q8a*s11*s12**2*s15*s24*s25**2*s33*s36**3 -  \
    512*q8a*s11*s12*s14*s16*s24*s25**2*s33*s36**3 + 512*q7a*s11*s12*s15*s22*s24*s25**2*s33*s36**3 - 1024*q7a*s11*s14*s16*s22*s24*s25**2*s33*s36**3 + 1024*q7a*s11*s12*s16*s24**2*s25**2*s33*s36**3 -  \
    512*q8a*s11*s12**2*s14*s25**3*s33*s36**3 + 1536*q7a*s11*s12*s14*s22*s25**3*s33*s36**3 - 1024*q7a*s11*s12**2*s24*s25**3*s33*s36**3 + 1024*q8a*s11**2*s13**2*s22**2*s26*s33*s36**3 - 6144*q8a*s11**2*s12*s13*s22*s23*s26*s33*s36**3 +  \
    3072*q8a*s11**2*s14**2*s22*s23*s26*s33*s36**3 + 2048*q7a*s11**2*s13*s22**2*s23*s26*s33*s36**3 + 6144*q8a*s11**2*s12**2*s23**2*s26*s33*s36**3 - 3072*q7a*s11**2*s12*s22*s23**2*s26*s33*s36**3 -  \
    1024*q8a*s11**2*s13*s14*s22*s24*s26*s33*s36**3 - 3072*q8a*s11**2*s12*s14*s23*s24*s26*s33*s36**3 - 1024*q7a*s11**2*s14*s22*s23*s24*s26*s33*s36**3 + 1024*q8a*s11**2*s12*s13*s24**2*s26*s33*s36**3 +  \
    1024*q7a*s11**2*s12*s23*s24**2*s26*s33*s36**3 + 2048*q8a*s11*s12*s13*s15*s22*s25*s26*s33*s36**3 - 1024*q8a*s11*s14**2*s15*s22*s25*s26*s33*s36**3 - 1024*q7a*s11*s13*s15*s22**2*s25*s26*s33*s36**3 -  \
    3072*q8a*s11*s12**2*s15*s23*s25*s26*s33*s36**3 + 2048*q7a*s11*s12*s15*s22*s23*s25*s26*s33*s36**3 + 1024*q8a*s11*s12*s14*s15*s24*s25*s26*s33*s36**3 + 512*q7a*s11*s14*s15*s22*s24*s25*s26*s33*s36**3 -  \
    512*q7a*s11*s12*s15*s24**2*s25*s26*s33*s36**3 + 1024*q8a*s11*s12**2*s13*s25**2*s26*s33*s36**3 - 2048*q7a*s11*s12*s13*s22*s25**2*s26*s33*s36**3 + 512*q7a*s11*s14**2*s22*s25**2*s26*s33*s36**3 +  \
    1024*q7a*s11*s12**2*s23*s25**2*s26*s33*s36**3 - 512*q7a*s11*s12*s14*s24*s25**2*s26*s33*s36**3 + 2048*q8a*s11**2*s13*s15*s22**2*s23*s34*s36**3 - 3072*q8a*s11**2*s12*s15*s22*s23**2*s34*s36**3 +  \
    3072*q8a*s11**2*s14*s16*s22*s23**2*s34*s36**3 + 1024*q7a*s11**2*s15*s22**2*s23**2*s34*s36**3 - 1024*q8a*s11**2*s14*s15*s22*s23*s24*s34*s36**3 - 1024*q8a*s11**2*s13*s16*s22*s23*s24*s34*s36**3 -  \
    1536*q8a*s11**2*s12*s16*s23**2*s24*s34*s36**3 - 512*q7a*s11**2*s16*s22*s23**2*s24*s34*s36**3 + 1024*q8a*s11**2*s12*s15*s23*s24**2*s34*s36**3 - 1024*q8a*s11*s13*s15**2*s22**2*s25*s34*s36**3 -  \
    2048*q8a*s11**2*s12*s13*s22*s23*s25*s34*s36**3 + 1024*q8a*s11**2*s14**2*s22*s23*s25*s34*s36**3 + 2048*q8a*s11*s12*s15**2*s22*s23*s25*s34*s36**3 - 2048*q8a*s11*s14*s15*s16*s22*s23*s25*s34*s36**3 -  \
    1024*q7a*s11*s15**2*s22**2*s23*s25*s34*s36**3 + 3072*q8a*s11**2*s12**2*s23**2*s25*s34*s36**3 - 1024*q7a*s11**2*s12*s22*s23**2*s25*s34*s36**3 + 512*q8a*s11*s14*s15**2*s22*s24*s25*s34*s36**3 +  \
    512*q8a*s11*s13*s15*s16*s22*s24*s25*s34*s36**3 - 1024*q8a*s11**2*s12*s14*s23*s24*s25*s34*s36**3 + 1024*q8a*s11*s12*s15*s16*s23*s24*s25*s34*s36**3 + 512*q7a*s11*s15*s16*s22*s23*s24*s25*s34*s36**3 -  \
    512*q8a*s11*s12*s15**2*s24**2*s25*s34*s36**3 - 512*q8a*s11*s14**2*s15*s22*s25**2*s34*s36**3 + 1024*q8a*s11*s13*s14*s16*s22*s25**2*s34*s36**3 + 2048*q7a*s11*s13*s15*s22**2*s25**2*s34*s36**3 -  \
    2048*q8a*s11*s12**2*s15*s23*s25**2*s34*s36**3 + 1024*q7a*s11*s14*s16*s22*s23*s25**2*s34*s36**3 + 512*q8a*s11*s12*s14*s15*s24*s25**2*s34*s36**3 - 512*q8a*s11*s12*s13*s16*s24*s25**2*s34*s36**3 -  \
    512*q7a*s11*s14*s15*s22*s24*s25**2*s34*s36**3 - 1024*q7a*s11*s13*s16*s22*s24*s25**2*s34*s36**3 - 512*q7a*s11*s12*s16*s23*s24*s25**2*s34*s36**3 + 512*q7a*s11*s12*s15*s24**2*s25**2*s34*s36**3 +  \
    1024*q8a*s11*s12**2*s13*s25**3*s34*s36**3 - 2048*q7a*s11*s12*s13*s22*s25**3*s34*s36**3 + 512*q7a*s11*s14**2*s22*s25**3*s34*s36**3 + 1024*q7a*s11*s12**2*s23*s25**3*s34*s36**3 - 512*q7a*s11*s12*s14*s24*s25**3*s34*s36**3 -  \
    1024*q8a*s11**2*s13*s14*s22*s23*s26*s34*s36**3 - 1536*q8a*s11**2*s12*s14*s23**2*s26*s34*s36**3 - 512*q7a*s11**2*s14*s22*s23**2*s26*s34*s36**3 + 2048*q8a*s11**2*s12*s13*s23*s24*s26*s34*s36**3 +  \
    1024*q7a*s11**2*s12*s23**2*s24*s26*s34*s36**3 + 512*q8a*s11*s13*s14*s15*s22*s25*s26*s34*s36**3 + 1024*q8a*s11*s12*s14*s15*s23*s25*s26*s34*s36**3 + 512*q7a*s11*s14*s15*s22*s23*s25*s26*s34*s36**3 -  \
    1024*q8a*s11*s12*s13*s15*s24*s25*s26*s34*s36**3 - 1024*q7a*s11*s12*s15*s23*s24*s25*s26*s34*s36**3 - 512*q8a*s11*s12*s13*s14*s25**2*s26*s34*s36**3 - 1024*q7a*s11*s13*s14*s22*s25**2*s26*s34*s36**3 -  \
    512*q7a*s11*s12*s14*s23*s25**2*s26*s34*s36**3 + 2048*q7a*s11*s12*s13*s24*s25**2*s26*s34*s36**3 + 1024*q8a*s11**2*s13**2*s22**2*s24*s35*s36**3 - 4096*q8a*s11**2*s12*s13*s22*s23*s24*s35*s36**3 +  \
    2048*q8a*s11**2*s14**2*s22*s23*s24*s35*s36**3 + 2048*q7a*s11**2*s13*s22**2*s23*s24*s35*s36**3 + 3072*q8a*s11**2*s12**2*s23**2*s24*s35*s36**3 - 2048*q7a*s11**2*s12*s22*s23**2*s24*s35*s36**3 -  \
    1024*q8a*s11**2*s13*s14*s22*s24**2*s35*s36**3 - 2048*q8a*s11**2*s12*s14*s23*s24**2*s35*s36**3 - 1024*q7a*s11**2*s14*s22*s23*s24**2*s35*s36**3 + 1024*q8a*s11**2*s12*s13*s24**3*s35*s36**3 + 1024*q7a*s11**2*s12*s23*s24**3*s35*s36**3 -  \
    512*q8a*s11*s13**2*s16*s22**2*s25*s35*s36**3 + 2048*q8a*s11*s12*s13*s16*s22*s23*s25*s35*s36**3 - 1024*q8a*s11*s14**2*s16*s22*s23*s25*s35*s36**3 - 1024*q7a*s11*s13*s16*s22**2*s23*s25*s35*s36**3 -  \
    1536*q8a*s11*s12**2*s16*s23**2*s25*s35*s36**3 + 1024*q7a*s11*s12*s16*s22*s23**2*s25*s35*s36**3 + 1024*q8a*s11*s12*s13*s15*s22*s24*s25*s35*s36**3 - 512*q8a*s11*s14**2*s15*s22*s24*s25*s35*s36**3 +  \
    512*q8a*s11*s13*s14*s16*s22*s24*s25*s35*s36**3 - 1024*q7a*s11*s13*s15*s22**2*s24*s25*s35*s36**3 - 1024*q8a*s11*s12**2*s15*s23*s24*s25*s35*s36**3 + 1024*q8a*s11*s12*s14*s16*s23*s24*s25*s35*s36**3 +  \
    1024*q7a*s11*s12*s15*s22*s23*s24*s25*s35*s36**3 + 512*q7a*s11*s14*s16*s22*s23*s24*s25*s35*s36**3 + 512*q8a*s11*s12*s14*s15*s24**2*s25*s35*s36**3 - 512*q8a*s11*s12*s13*s16*s24**2*s25*s35*s36**3 +  \
    512*q7a*s11*s14*s15*s22*s24**2*s25*s35*s36**3 - 512*q7a*s11*s12*s16*s23*s24**2*s25*s35*s36**3 - 512*q7a*s11*s12*s15*s24**3*s25*s35*s36**3 - 1024*q8a*s11*s12*s13*s14*s22*s25**2*s35*s36**3 +  \
    512*q8a*s11*s14**3*s22*s25**2*s35*s36**3 + 1024*q7a*s11*s13*s14*s22**2*s25**2*s35*s36**3 + 1024*q8a*s11*s12**2*s14*s23*s25**2*s35*s36**3 - 1024*q7a*s11*s12*s14*s22*s23*s25**2*s35*s36**3 -  \
    512*q8a*s11*s12*s14**2*s24*s25**2*s35*s36**3 - 512*q7a*s11*s14**2*s22*s24*s25**2*s35*s36**3 + 512*q7a*s11*s12*s14*s24**2*s25**2*s35*s36**3 - 512*q8a*s11*s13**2*s15*s22**2*s26*s35*s36**3 +  \
    2048*q8a*s11*s12*s13*s15*s22*s23*s26*s35*s36**3 - 1024*q8a*s11*s14**2*s15*s22*s23*s26*s35*s36**3 - 1024*q7a*s11*s13*s15*s22**2*s23*s26*s35*s36**3 - 1536*q8a*s11*s12**2*s15*s23**2*s26*s35*s36**3 +  \
    1024*q7a*s11*s12*s15*s22*s23**2*s26*s35*s36**3 + 512*q8a*s11*s13*s14*s15*s22*s24*s26*s35*s36**3 + 1024*q8a*s11*s12*s14*s15*s23*s24*s26*s35*s36**3 + 512*q7a*s11*s14*s15*s22*s23*s24*s26*s35*s36**3 -  \
    512*q8a*s11*s12*s13*s15*s24**2*s26*s35*s36**3 - 512*q7a*s11*s12*s15*s23*s24**2*s26*s35*s36**3 - 2048*q8a*s11*s12*s13**2*s22*s25*s26*s35*s36**3 + 1024*q8a*s11*s13*s14**2*s22*s25*s26*s35*s36**3 +  \
    3072*q7a*s11*s13**2*s22**2*s25*s26*s35*s36**3 + 2048*q8a*s11*s12**2*s13*s23*s25*s26*s35*s36**3 - 4096*q7a*s11*s12*s13*s22*s23*s25*s26*s35*s36**3 + 1024*q7a*s11*s14**2*s22*s23*s25*s26*s35*s36**3 +  \
    1024*q7a*s11*s12**2*s23**2*s25*s26*s35*s36**3 - 1024*q8a*s11*s12*s13*s14*s24*s25*s26*s35*s36**3 - 2048*q7a*s11*s13*s14*s22*s24*s25*s26*s35*s36**3 - 1024*q7a*s11*s12*s14*s23*s24*s25*s26*s35*s36**3 +  \
    2048*q7a*s11*s12*s13*s24**2*s25*s26*s35*s36**3 - 1024*q8a*s11**2*s13**2*s22**2*s23*s36**4 + 3072*q8a*s11**2*s12*s13*s22*s23**2*s36**4 - 1536*q8a*s11**2*s14**2*s22*s23**2*s36**4 - 1024*q7a*s11**2*s13*s22**2*s23**2*s36**4 -  \
    2048*q8a*s11**2*s12**2*s23**3*s36**4 + 1024*q7a*s11**2*s12*s22*s23**3*s36**4 + 1024*q8a*s11**2*s13*s14*s22*s23*s24*s36**4 + 1536*q8a*s11**2*s12*s14*s23**2*s24*s36**4 + 512*q7a*s11**2*s14*s22*s23**2*s24*s36**4 -  \
    1024*q8a*s11**2*s12*s13*s23*s24**2*s36**4 - 512*q7a*s11**2*s12*s23**2*s24**2*s36**4 + 512*q8a*s11*s13**2*s15*s22**2*s25*s36**4 - 2048*q8a*s11*s12*s13*s15*s22*s23*s25*s36**4 + 1024*q8a*s11*s14**2*s15*s22*s23*s25*s36**4 +  \
    1024*q7a*s11*s13*s15*s22**2*s23*s25*s36**4 + 1536*q8a*s11*s12**2*s15*s23**2*s25*s36**4 - 1024*q7a*s11*s12*s15*s22*s23**2*s25*s36**4 - 512*q8a*s11*s13*s14*s15*s22*s24*s25*s36**4 - 1024*q8a*s11*s12*s14*s15*s23*s24*s25*s36**4 -  \
    512*q7a*s11*s14*s15*s22*s23*s24*s25*s36**4 + 512*q8a*s11*s12*s13*s15*s24**2*s25*s36**4 + 512*q7a*s11*s12*s15*s23*s24**2*s25*s36**4 + 1024*q8a*s11*s12*s13**2*s22*s25**2*s36**4 - 512*q8a*s11*s13*s14**2*s22*s25**2*s36**4 -  \
    1536*q7a*s11*s13**2*s22**2*s25**2*s36**4 - 1024*q8a*s11*s12**2*s13*s23*s25**2*s36**4 + 2048*q7a*s11*s12*s13*s22*s23*s25**2*s36**4 - 512*q7a*s11*s14**2*s22*s23*s25**2*s36**4 - 512*q7a*s11*s12**2*s23**2*s25**2*s36**4 +  \
    512*q8a*s11*s12*s13*s14*s24*s25**2*s36**4 + 1024*q7a*s11*s13*s14*s22*s24*s25**2*s36**4 + 512*q7a*s11*s12*s14*s23*s24*s25**2*s36**4 - 1024*q7a*s11*s12*s13*s24**2*s25**2*s36**4
v3=v3_0+v3_1+v3_2+v3_3+v3_4+v3_5+v3_6+v3_7+v3_8+v3_9

v4_0= \
    512*s11**4*s22**4*s33**4 - 1024*s11**3*s16*s22**3*s26*s33**4 + 1024*s11**3*s12*s22**2*s26**2*s33**4 + 512*s11**2*s16**2*s22**2*s26**2*s33**4 - 1024*s11**2*s12*s16*s22*s26**3*s33**4 + 512*s11**2*s12**2*s26**4*s33**4 -  \
    1024*s11**4*s22**3*s24*s33**3*s34 + 512*s11**3*s16*s22**3*s25*s33**3*s34 + 512*s11**3*s15*s22**3*s26*s33**3*s34 + 1536*s11**3*s16*s22**2*s24*s26*s33**3*s34 - 1024*s11**3*s12*s22**2*s25*s26*s33**3*s34 -  \
    512*s11**2*s16**2*s22**2*s25*s26*s33**3*s34 - 512*s11**3*s14*s22**2*s26**2*s33**3*s34 - 512*s11**2*s15*s16*s22**2*s26**2*s33**3*s34 - 1024*s11**3*s12*s22*s24*s26**2*s33**3*s34 - 512*s11**2*s16**2*s22*s24*s26**2*s33**3*s34 +  \
    1536*s11**2*s12*s16*s22*s25*s26**2*s33**3*s34 + 512*s11**2*s12*s15*s22*s26**3*s33**3*s34 + 512*s11**2*s14*s16*s22*s26**3*s33**3*s34 + 512*s11**2*s12*s16*s24*s26**3*s33**3*s34 - 1024*s11**2*s12**2*s25*s26**3*s33**3*s34 -  \
    512*s11**2*s12*s14*s26**4*s33**3*s34 + 1024*s11**4*s22**3*s23*s33**2*s34**2 + 512*s11**4*s22**2*s24**2*s33**2*s34**2 - 512*s11**3*s15*s22**3*s25*s33**2*s34**2 - 512*s11**3*s16*s22**2*s24*s25*s33**2*s34**2 +  \
    512*s11**3*s12*s22**2*s25**2*s33**2*s34**2 - 1536*s11**3*s16*s22**2*s23*s26*s33**2*s34**2 - 512*s11**3*s15*s22**2*s24*s26*s33**2*s34**2 - 512*s11**3*s16*s22*s24**2*s26*s33**2*s34**2 + 1024*s11**3*s14*s22**2*s25*s26*s33**2*s34**2 +  \
    512*s11**2*s15*s16*s22**2*s25*s26*s33**2*s34**2 + 512*s11**2*s16**2*s22*s24*s25*s26*s33**2*s34**2 - 512*s11**2*s12*s16*s22*s25**2*s26*s33**2*s34**2 + 512*s11**3*s13*s22**2*s26**2*s33**2*s34**2 +  \
    1024*s11**3*s12*s22*s23*s26**2*s33**2*s34**2 + 512*s11**2*s16**2*s22*s23*s26**2*s33**2*s34**2 + 512*s11**2*s15*s16*s22*s24*s26**2*s33**2*s34**2 + 512*s11**3*s12*s24**2*s26**2*s33**2*s34**2 -  \
    512*s11**2*s12*s15*s22*s25*s26**2*s33**2*s34**2 - 1024*s11**2*s14*s16*s22*s25*s26**2*s33**2*s34**2 - 512*s11**2*s12*s16*s24*s25*s26**2*s33**2*s34**2 + 512*s11**2*s12**2*s25**2*s26**2*s33**2*s34**2 -  \
    512*s11**2*s13*s16*s22*s26**3*s33**2*s34**2 - 512*s11**2*s12*s16*s23*s26**3*s33**2*s34**2 - 512*s11**2*s12*s15*s24*s26**3*s33**2*s34**2 + 1024*s11**2*s12*s14*s25*s26**3*s33**2*s34**2 + 512*s11**2*s12*s13*s26**4*s33**2*s34**2 -  \
    1024*s11**4*s22**2*s23*s24*s33*s34**3 + 512*s11**3*s16*s22**2*s23*s25*s33*s34**3 + 512*s11**3*s15*s22**2*s24*s25*s33*s34**3 - 512*s11**3*s14*s22**2*s25**2*s33*s34**3 + 512*s11**3*s15*s22**2*s23*s26*s33*s34**3 +  \
    1024*s11**3*s16*s22*s23*s24*s26*s33*s34**3 - 1024*s11**3*s13*s22**2*s25*s26*s33*s34**3 - 512*s11**2*s16**2*s22*s23*s25*s26*s33*s34**3 - 512*s11**2*s15*s16*s22*s24*s25*s26*s33*s34**3 +  \
    512*s11**2*s14*s16*s22*s25**2*s26*s33*s34**3 - 512*s11**2*s15*s16*s22*s23*s26**2*s33*s34**3 - 1024*s11**3*s12*s23*s24*s26**2*s33*s34**3 + 1024*s11**2*s13*s16*s22*s25*s26**2*s33*s34**3
v4_1= \
    512*s11**2*s12*s16*s23*s25*s26**2*s33*s34**3 + 512*s11**2*s12*s15*s24*s25*s26**2*s33*s34**3 - 512*s11**2*s12*s14*s25**2*s26**2*s33*s34**3 + 512*s11**2*s12*s15*s23*s26**3*s33*s34**3 - 1024*s11**2*s12*s13*s25*s26**3*s33*s34**3 +  \
    512*s11**4*s22**2*s23**2*s34**4 - 512*s11**3*s15*s22**2*s23*s25*s34**4 + 512*s11**3*s13*s22**2*s25**2*s34**4 - 512*s11**3*s16*s22*s23**2*s26*s34**4 + 512*s11**2*s15*s16*s22*s23*s25*s26*s34**4 -  \
    512*s11**2*s13*s16*s22*s25**2*s26*s34**4 + 512*s11**3*s12*s23**2*s26**2*s34**4 - 512*s11**2*s12*s15*s23*s25*s26**2*s34**4 + 512*s11**2*s12*s13*s25**2*s26**2*s34**4 - 1024*s11**3*s15*s22**4*s33**3*s35 +  \
    512*s11**3*s16*s22**3*s24*s33**3*s35 + 1024*s11**3*s12*s22**3*s25*s33**3*s35 - 512*s11**2*s16**2*s22**3*s25*s33**3*s35 + 512*s11**3*s14*s22**3*s26*s33**3*s35 + 1536*s11**2*s15*s16*s22**3*s26*s33**3*s35 -  \
    1024*s11**3*s12*s22**2*s24*s26*s33**3*s35 - 512*s11**2*s16**2*s22**2*s24*s26*s33**3*s35 - 512*s11**2*s12*s16*s22**2*s25*s26*s33**3*s35 + 512*s11*s16**3*s22**2*s25*s26*s33**3*s35 - 1536*s11**2*s12*s15*s22**2*s26**2*s33**3*s35 -  \
    512*s11**2*s14*s16*s22**2*s26**2*s33**3*s35 - 512*s11*s15*s16**2*s22**2*s26**2*s33**3*s35 + 1536*s11**2*s12*s16*s22*s24*s26**2*s33**3*s35 + 1024*s11**2*s12**2*s22*s25*s26**2*s33**3*s35 - 1024*s11*s12*s16**2*s22*s25*s26**2*s33**3*s35 +  \
    512*s11**2*s12*s14*s22*s26**3*s33**3*s35 + 1024*s11*s12*s15*s16*s22*s26**3*s33**3*s35 - 1024*s11**2*s12**2*s24*s26**3*s33**3*s35 + 512*s11*s12**2*s16*s25*s26**3*s33**3*s35 - 512*s11*s12**2*s15*s26**4*s33**3*s35 -  \
    1024*s11**3*s16*s22**3*s23*s33**2*s34*s35 + 1536*s11**3*s15*s22**3*s24*s33**2*s34*s35 - 512*s11**3*s16*s22**2*s24**2*s33**2*s34*s35 - 1536*s11**3*s14*s22**3*s25*s33**2*s34*s35 + 512*s11**2*s15*s16*s22**3*s25*s33**2*s34*s35 +  \
    512*s11**2*s16**2*s22**2*s24*s25*s33**2*s34*s35 - 512*s11**2*s12*s16*s22**2*s25**2*s33**2*s34*s35 - 1024*s11**3*s13*s22**3*s26*s33**2*s34*s35 - 512*s11**2*s15**2*s22**3*s26*s33**2*s34*s35 +  \
    2048*s11**3*s12*s22**2*s23*s26*s33**2*s34*s35 + 1024*s11**2*s16**2*s22**2*s23*s26*s33**2*s34*s35 + 512*s11**3*s14*s22**2*s24*s26*s33**2*s34*s35 - 2048*s11**2*s15*s16*s22**2*s24*s26*s33**2*s34*s35 +  \
    512*s11**2*s16**2*s22*s24**2*s26*s33**2*s34*s35 + 512*s11**2*s12*s15*s22**2*s25*s26*s33**2*s34*s35 + 1024*s11**2*s14*s16*s22**2*s25*s26*s33**2*s34*s35 - 512*s11*s15*s16**2*s22**2*s25*s26*s33**2*s34*s35 -  \
    512*s11*s16**3*s22*s24*s25*s26*s33**2*s34*s35 + 512*s11*s12*s16**2*s22*s25**2*s26*s33**2*s34*s35 + 512*s11**2*s14*s15*s22**2*s26**2*s33**2*s34*s35 + 1024*s11**2*s13*s16*s22**2*s26**2*s33**2*s34*s35 +  \
    512*s11*s15**2*s16*s22**2*s26**2*s33**2*s34*s35 - 3072*s11**2*s12*s16*s22*s23*s26**2*s33**2*s34*s35 + 1536*s11**2*s12*s15*s22*s24*s26**2*s33**2*s34*s35 - 512*s11**2*s14*s16*s22*s24*s26**2*s33**2*s34*s35 +  \
    512*s11*s15*s16**2*s22*s24*s26**2*s33**2*s34*s35 - 512*s11**2*s12*s16*s24**2*s26**2*s33**2*s34*s35 - 1536*s11**2*s12*s14*s22*s25*s26**2*s33**2*s34*s35 + 512*s11*s14*s16**2*s22*s25*s26**2*s33**2*s34*s35 +  \
    512*s11*s12*s16**2*s24*s25*s26**2*s33**2*s34*s35 - 512*s11*s12**2*s16*s25**2*s26**2*s33**2*s34*s35 - 1024*s11**2*s12*s13*s22*s26**3*s33**2*s34*s35 - 512*s11*s12*s15**2*s22*s26**3*s33**2*s34*s35 -  \
    512*s11*s14*s15*s16*s22*s26**3*s33**2*s34*s35 + 2048*s11**2*s12**2*s23*s26**3*s33**2*s34*s35 + 512*s11**2*s12*s14*s24*s26**3*s33**2*s34*s35 - 512*s11*s12*s15*s16*s24*s26**3*s33**2*s34*s35 +  \
    512*s11*s12**2*s15*s25*s26**3*s33**2*s34*s35 - 512*s11*s12*s14*s16*s25*s26**3*s33**2*s34*s35 + 512*s11*s12*s14*s15*s26**4*s33**2*s34*s35 - 1024*s11**3*s15*s22**3*s23*s33*s34**2*s35 +  \
    1536*s11**3*s16*s22**2*s23*s24*s33*s34**2*s35 - 512*s11**3*s15*s22**2*s24**2*s33*s34**2*s35 + 2048*s11**3*s13*s22**3*s25*s33*s34**2*s35 - 1024*s11**3*s12*s22**2*s23*s25*s33*s34**2*s35 -  \
    512*s11**2*s16**2*s22**2*s23*s25*s33*s34**2*s35 + 512*s11**3*s14*s22**2*s24*s25*s33*s34**2*s35 - 512*s11**2*s15*s16*s22**2*s24*s25*s33*s34**2*s35 + 512*s11**2*s14*s16*s22**2*s25**2*s33*s34**2*s35 -  \
    1536*s11**3*s14*s22**2*s23*s26*s33*s34**2*s35 + 1536*s11**2*s15*s16*s22**2*s23*s26*s33*s34**2*s35 + 512*s11**2*s15**2*s22**2*s24*s26*s33*s34**2*s35 - 1536*s11**2*s16**2*s22*s23*s24*s26*s33*s34**2*s35 +  \
    512*s11**2*s15*s16*s22*s24**2*s26*s33*s34**2*s35 - 512*s11**2*s14*s15*s22**2*s25*s26*s33*s34**2*s35 - 1536*s11**2*s13*s16*s22**2*s25*s26*s33*s34**2*s35 + 1024*s11**2*s12*s16*s22*s23*s25*s26*s33*s34**2*s35 +  \
    512*s11*s16**3*s22*s23*s25*s26*s33*s34**2*s35 - 512*s11**2*s14*s16*s22*s24*s25*s26*s33*s34**2*s35 + 512*s11*s15*s16**2*s22*s24*s25*s26*s33*s34**2*s35 - 512*s11*s14*s16**2*s22*s25**2*s26*s33*s34**2*s35 -  \
    512*s11**2*s13*s15*s22**2*s26**2*s33*s34**2*s35 - 1024*s11**2*s12*s15*s22*s23*s26**2*s33*s34**2*s35 + 1536*s11**2*s14*s16*s22*s23*s26**2*s33*s34**2*s35 - 512*s11*s15*s16**2*s22*s23*s26**2*s33*s34**2*s35 -  \
    512*s11*s15**2*s16*s22*s24*s26**2*s33*s34**2*s35 + 1536*s11**2*s12*s16*s23*s24*s26**2*s33*s34**2*s35 - 512*s11**2*s12*s15*s24**2*s26**2*s33*s34**2*s35 + 2048*s11**2*s12*s13*s22*s25*s26**2*s33*s34**2*s35 +  \
    512*s11*s14*s15*s16*s22*s25*s26**2*s33*s34**2*s35 - 512*s11*s13*s16**2*s22*s25*s26**2*s33*s34**2*s35 - 1024*s11**2*s12**2*s23*s25*s26**2*s33*s34**2*s35 - 512*s11*s12*s16**2*s23*s25*s26**2*s33*s34**2*s35 +  \
    512*s11**2*s12*s14*s24*s25*s26**2*s33*s34**2*s35 - 512*s11*s12*s15*s16*s24*s25*s26**2*s33*s34**2*s35 + 512*s11*s12*s14*s16*s25**2*s26**2*s33*s34**2*s35 + 512*s11*s13*s15*s16*s22*s26**3*s33*s34**2*s35 -  \
    1536*s11**2*s12*s14*s23*s26**3*s33*s34**2*s35 + 512*s11*s12*s15*s16*s23*s26**3*s33*s34**2*s35 + 512*s11*s12*s15**2*s24*s26**3*s33*s34**2*s35 - 512*s11*s12*s14*s15*s25*s26**3*s33*s34**2*s35 +  \
    512*s11*s12*s13*s16*s25*s26**3*s33*s34**2*s35 - 512*s11*s12*s13*s15*s26**4*s33*s34**2*s35 - 1024*s11**3*s16*s22**2*s23**2*s34**3*s35 + 512*s11**3*s15*s22**2*s23*s24*s34**3*s35 + 512*s11**3*s14*s22**2*s23*s25*s34**3*s35 +  \
    512*s11**2*s15*s16*s22**2*s23*s25*s34**3*s35 - 1024*s11**3*s13*s22**2*s24*s25*s34**3*s35 - 512*s11**2*s13*s16*s22**2*s25**2*s34**3*s35 + 1024*s11**3*s13*s22**2*s23*s26*s34**3*s35 - 512*s11**2*s15**2*s22**2*s23*s26*s34**3*s35 +  \
    1024*s11**2*s16**2*s22*s23**2*s26*s34**3*s35 - 512*s11**2*s15*s16*s22*s23*s24*s26*s34**3*s35 + 512*s11**2*s13*s15*s22**2*s25*s26*s34**3*s35 - 512*s11**2*s14*s16*s22*s23*s25*s26*s34**3*s35 -  \
    512*s11*s15*s16**2*s22*s23*s25*s26*s34**3*s35 + 1024*s11**2*s13*s16*s22*s24*s25*s26*s34**3*s35 + 512*s11*s13*s16**2*s22*s25**2*s26*s34**3*s35 - 1024*s11**2*s13*s16*s22*s23*s26**2*s34**3*s35 +  \
    512*s11*s15**2*s16*s22*s23*s26**2*s34**3*s35 - 1024*s11**2*s12*s16*s23**2*s26**2*s34**3*s35 + 512*s11**2*s12*s15*s23*s24*s26**2*s34**3*s35 - 512*s11*s13*s15*s16*s22*s25*s26**2*s34**3*s35 +  \
    512*s11**2*s12*s14*s23*s25*s26**2*s34**3*s35 + 512*s11*s12*s15*s16*s23*s25*s26**2*s34**3*s35 - 1024*s11**2*s12*s13*s24*s25*s26**2*s34**3*s35 - 512*s11*s12*s13*s16*s25**2*s26**2*s34**3*s35 +  \
    1024*s11**2*s12*s13*s23*s26**3*s34**3*s35 - 512*s11*s12*s15**2*s23*s26**3*s34**3*s35 + 512*s11*s12*s13*s15*s25*s26**3*s34**3*s35 + 1024*s11**3*s13*s22**4*s33**2*s35**2 + 512*s11**2*s15**2*s22**4*s33**2*s35**2 -  \
    1024*s11**3*s12*s22**3*s23*s33**2*s35**2 + 512*s11**2*s16**2*s22**3*s23*s33**2*s35**2 - 512*s11**3*s14*s22**3*s24*s33**2*s35**2 - 512*s11**2*s15*s16*s22**3*s24*s33**2*s35**2 + 512*s11**3*s12*s22**2*s24**2*s33**2*s35**2 -  \
    1024*s11**2*s12*s15*s22**3*s25*s33**2*s35**2 + 1024*s11**2*s14*s16*s22**3*s25*s33**2*s35**2 - 512*s11**2*s12*s16*s22**2*s24*s25*s33**2*s35**2 + 512*s11**2*s12**2*s22**2*s25**2*s33**2*s35**2 - 512*s11**2*s14*s15*s22**3*s26*s33**2*s35**2 -  \
    1536*s11**2*s13*s16*s22**3*s26*s33**2*s35**2 - 512*s11*s15**2*s16*s22**3*s26*s33**2*s35**2 + 512*s11**2*s12*s16*s22**2*s23*s26*s33**2*s35**2 - 512*s11*s16**3*s22**2*s23*s26*s33**2*s35**2
v4_2= \
    1024*s11**2*s12*s15*s22**2*s24*s26*s33**2*s35**2 + 512*s11**2*s14*s16*s22**2*s24*s26*s33**2*s35**2 + 512*s11*s15*s16**2*s22**2*s24*s26*s33**2*s35**2 - 512*s11**2*s12*s16*s22*s24**2*s26*s33**2*s35**2 -  \
    512*s11**2*s12*s14*s22**2*s25*s26*s33**2*s35**2 + 1024*s11*s12*s15*s16*s22**2*s25*s26*s33**2*s35**2 - 1024*s11*s14*s16**2*s22**2*s25*s26*s33**2*s35**2 + 512*s11*s12*s16**2*s22*s24*s25*s26*s33**2*s35**2 -  \
    512*s11*s12**2*s16*s22*s25**2*s26*s33**2*s35**2 + 1536*s11**2*s12*s13*s22**2*s26**2*s33**2*s35**2 + 512*s11*s12*s15**2*s22**2*s26**2*s33**2*s35**2 + 512*s11*s14*s15*s16*s22**2*s26**2*s33**2*s35**2 +  \
    512*s11*s13*s16**2*s22**2*s26**2*s33**2*s35**2 - 1024*s11**2*s12**2*s22*s23*s26**2*s33**2*s35**2 + 1024*s11*s12*s16**2*s22*s23*s26**2*s33**2*s35**2 - 512*s11**2*s12*s14*s22*s24*s26**2*s33**2*s35**2 -  \
    1536*s11*s12*s15*s16*s22*s24*s26**2*s33**2*s35**2 + 512*s11**2*s12**2*s24**2*s26**2*s33**2*s35**2 - 1024*s11*s12**2*s15*s22*s25*s26**2*s33**2*s35**2 + 1536*s11*s12*s14*s16*s22*s25*s26**2*s33**2*s35**2 -  \
    512*s11*s12**2*s16*s24*s25*s26**2*s33**2*s35**2 + 512*s11*s12**3*s25**2*s26**2*s33**2*s35**2 - 512*s11*s12*s14*s15*s22*s26**3*s33**2*s35**2 - 1024*s11*s12*s13*s16*s22*s26**3*s33**2*s35**2 - 512*s11*s12**2*s16*s23*s26**3*s33**2*s35**2 +  \
    1024*s11*s12**2*s15*s24*s26**3*s33**2*s35**2 - 512*s11*s12**2*s14*s25*s26**3*s33**2*s35**2 + 512*s11*s12**2*s13*s26**4*s33**2*s35**2 + 2048*s11**3*s14*s22**3*s23*s33*s34*s35**2 - 1024*s11**3*s13*s22**3*s24*s33*s34*s35**2 -  \
    512*s11**2*s15**2*s22**3*s24*s33*s34*s35**2 - 1024*s11**3*s12*s22**2*s23*s24*s33*s34*s35**2 - 512*s11**2*s16**2*s22**2*s23*s24*s33*s34*s35**2 + 512*s11**2*s15*s16*s22**2*s24**2*s33*s34*s35**2 +  \
    512*s11**2*s14*s15*s22**3*s25*s33*s34*s35**2 - 1536*s11**2*s13*s16*s22**3*s25*s33*s34*s35**2 + 1536*s11**2*s12*s16*s22**2*s23*s25*s33*s34*s35**2 + 512*s11**2*s12*s15*s22**2*s24*s25*s33*s34*s35**2 -  \
    512*s11**2*s14*s16*s22**2*s24*s25*s33*s34*s35**2 - 512*s11**2*s12*s14*s22**2*s25**2*s33*s34*s35**2 + 1536*s11**2*s13*s15*s22**3*s26*s33*s34*s35**2 - 1536*s11**2*s12*s15*s22**2*s23*s26*s33*s34*s35**2 -  \
    1536*s11**2*s14*s16*s22**2*s23*s26*s33*s34*s35**2 - 512*s11**2*s14*s15*s22**2*s24*s26*s33*s34*s35**2 + 1536*s11**2*s13*s16*s22**2*s24*s26*s33*s34*s35**2 + 512*s11*s15**2*s16*s22**2*s24*s26*s33*s34*s35**2 +  \
    1024*s11**2*s12*s16*s22*s23*s24*s26*s33*s34*s35**2 + 512*s11*s16**3*s22*s23*s24*s26*s33*s34*s35**2 - 512*s11*s15*s16**2*s22*s24**2*s26*s33*s34*s35**2 + 512*s11**2*s14**2*s22**2*s25*s26*s33*s34*s35**2 -  \
    512*s11*s14*s15*s16*s22**2*s25*s26*s33*s34*s35**2 + 1536*s11*s13*s16**2*s22**2*s25*s26*s33*s34*s35**2 - 1536*s11*s12*s16**2*s22*s23*s25*s26*s33*s34*s35**2 - 512*s11*s12*s15*s16*s22*s24*s25*s26*s33*s34*s35**2 +  \
    512*s11*s14*s16**2*s22*s24*s25*s26*s33*s34*s35**2 + 512*s11*s12*s14*s16*s22*s25**2*s26*s33*s34*s35**2 - 512*s11**2*s13*s14*s22**2*s26**2*s33*s34*s35**2 - 1536*s11*s13*s15*s16*s22**2*s26**2*s33*s34*s35**2 +  \
    2048*s11**2*s12*s14*s22*s23*s26**2*s33*s34*s35**2 + 1536*s11*s12*s15*s16*s22*s23*s26**2*s33*s34*s35**2 - 512*s11*s14*s16**2*s22*s23*s26**2*s33*s34*s35**2 - 1024*s11**2*s12*s13*s22*s24*s26**2*s33*s34*s35**2 -  \
    512*s11*s12*s15**2*s22*s24*s26**2*s33*s34*s35**2 + 512*s11*s14*s15*s16*s22*s24*s26**2*s33*s34*s35**2 - 512*s11*s13*s16**2*s22*s24*s26**2*s33*s34*s35**2 - 1024*s11**2*s12**2*s23*s24*s26**2*s33*s34*s35**2 -  \
    512*s11*s12*s16**2*s23*s24*s26**2*s33*s34*s35**2 + 512*s11*s12*s15*s16*s24**2*s26**2*s33*s34*s35**2 + 512*s11*s12*s14*s15*s22*s25*s26**2*s33*s34*s35**2 - 1536*s11*s12*s13*s16*s22*s25*s26**2*s33*s34*s35**2 -  \
    512*s11*s14**2*s16*s22*s25*s26**2*s33*s34*s35**2 + 1536*s11*s12**2*s16*s23*s25*s26**2*s33*s34*s35**2 + 512*s11*s12**2*s15*s24*s25*s26**2*s33*s34*s35**2 - 512*s11*s12*s14*s16*s24*s25*s26**2*s33*s34*s35**2 -  \
    512*s11*s12**2*s14*s25**2*s26**2*s33*s34*s35**2 + 1536*s11*s12*s13*s15*s22*s26**3*s33*s34*s35**2 + 512*s11*s13*s14*s16*s22*s26**3*s33*s34*s35**2 - 1536*s11*s12**2*s15*s23*s26**3*s33*s34*s35**2 +  \
    512*s11*s12*s14*s16*s23*s26**3*s33*s34*s35**2 - 512*s11*s12*s14*s15*s24*s26**3*s33*s34*s35**2 + 512*s11*s12*s13*s16*s24*s26**3*s33*s34*s35**2 + 512*s11*s12*s14**2*s25*s26**3*s33*s34*s35**2 -  \
    512*s11*s12*s13*s14*s26**4*s33*s34*s35**2 - 1024*s11**3*s13*s22**3*s23*s34**2*s35**2 + 512*s11**2*s15**2*s22**3*s23*s34**2*s35**2 + 1024*s11**3*s12*s22**2*s23**2*s34**2*s35**2 + 512*s11**2*s16**2*s22**2*s23**2*s34**2*s35**2 -  \
    512*s11**3*s14*s22**2*s23*s24*s34**2*s35**2 - 512*s11**2*s15*s16*s22**2*s23*s24*s34**2*s35**2 + 512*s11**3*s13*s22**2*s24**2*s34**2*s35**2 - 512*s11**2*s13*s15*s22**3*s25*s34**2*s35**2 - 512*s11**2*s12*s15*s22**2*s23*s25*s34**2*s35**2 -  \
    512*s11**2*s14*s16*s22**2*s23*s25*s34**2*s35**2 + 1024*s11**2*s13*s16*s22**2*s24*s25*s34**2*s35**2 + 512*s11**2*s12*s13*s22**2*s25**2*s34**2*s35**2 + 1024*s11**2*s14*s15*s22**2*s23*s26*s34**2*s35**2 -  \
    512*s11*s15**2*s16*s22**2*s23*s26*s34**2*s35**2 - 1024*s11**2*s12*s16*s22*s23**2*s26*s34**2*s35**2 - 512*s11*s16**3*s22*s23**2*s26*s34**2*s35**2 - 512*s11**2*s13*s15*s22**2*s24*s26*s34**2*s35**2 +  \
    512*s11**2*s14*s16*s22*s23*s24*s26*s34**2*s35**2 + 512*s11*s15*s16**2*s22*s23*s24*s26*s34**2*s35**2 - 512*s11**2*s13*s16*s22*s24**2*s26*s34**2*s35**2 - 512*s11**2*s13*s14*s22**2*s25*s26*s34**2*s35**2 +  \
    512*s11*s13*s15*s16*s22**2*s25*s26*s34**2*s35**2 + 512*s11*s12*s15*s16*s22*s23*s25*s26*s34**2*s35**2 + 512*s11*s14*s16**2*s22*s23*s25*s26*s34**2*s35**2 - 1024*s11*s13*s16**2*s22*s24*s25*s26*s34**2*s35**2 -  \
    512*s11*s12*s13*s16*s22*s25**2*s26*s34**2*s35**2 + 512*s11**2*s13**2*s22**2*s26**2*s34**2*s35**2 - 1024*s11**2*s12*s13*s22*s23*s26**2*s34**2*s35**2 + 512*s11*s12*s15**2*s22*s23*s26**2*s34**2*s35**2 -  \
    1024*s11*s14*s15*s16*s22*s23*s26**2*s34**2*s35**2 + 1024*s11*s13*s16**2*s22*s23*s26**2*s34**2*s35**2 + 1024*s11**2*s12**2*s23**2*s26**2*s34**2*s35**2 + 512*s11*s12*s16**2*s23**2*s26**2*s34**2*s35**2 +  \
    512*s11*s13*s15*s16*s22*s24*s26**2*s34**2*s35**2 - 512*s11**2*s12*s14*s23*s24*s26**2*s34**2*s35**2 - 512*s11*s12*s15*s16*s23*s24*s26**2*s34**2*s35**2 + 512*s11**2*s12*s13*s24**2*s26**2*s34**2*s35**2 -  \
    512*s11*s12*s13*s15*s22*s25*s26**2*s34**2*s35**2 + 512*s11*s13*s14*s16*s22*s25*s26**2*s34**2*s35**2 - 512*s11*s12**2*s15*s23*s25*s26**2*s34**2*s35**2 - 512*s11*s12*s14*s16*s23*s25*s26**2*s34**2*s35**2 +  \
    1024*s11*s12*s13*s16*s24*s25*s26**2*s34**2*s35**2 + 512*s11*s12**2*s13*s25**2*s26**2*s34**2*s35**2 - 512*s11*s13**2*s16*s22*s26**3*s34**2*s35**2 + 1024*s11*s12*s14*s15*s23*s26**3*s34**2*s35**2 -  \
    1024*s11*s12*s13*s16*s23*s26**3*s34**2*s35**2 - 512*s11*s12*s13*s15*s24*s26**3*s34**2*s35**2 - 512*s11*s12*s13*s14*s25*s26**3*s34**2*s35**2 + 512*s11*s12*s13**2*s26**4*s34**2*s35**2 - 1024*s11**2*s13*s15*s22**4*s33*s35**3 +  \
    1024*s11**2*s12*s15*s22**3*s23*s33*s35**3 - 1024*s11**2*s14*s16*s22**3*s23*s33*s35**3 + 512*s11**2*s14*s15*s22**3*s24*s33*s35**3 + 512*s11**2*s13*s16*s22**3*s24*s33*s35**3 + 512*s11**2*s12*s16*s22**2*s23*s24*s33*s35**3 -  \
    512*s11**2*s12*s15*s22**2*s24**2*s33*s35**3 + 1024*s11**2*s12*s13*s22**3*s25*s33*s35**3 - 512*s11**2*s14**2*s22**3*s25*s33*s35**3 - 1024*s11**2*s12**2*s22**2*s23*s25*s33*s35**3 + 512*s11**2*s12*s14*s22**2*s24*s25*s33*s35**3 +  \
    512*s11**2*s13*s14*s22**3*s26*s33*s35**3 + 1024*s11*s13*s15*s16*s22**3*s26*s33*s35**3 + 512*s11**2*s12*s14*s22**2*s23*s26*s33*s35**3 - 1024*s11*s12*s15*s16*s22**2*s23*s26*s33*s35**3 +  \
    1024*s11*s14*s16**2*s22**2*s23*s26*s33*s35**3 - 1024*s11**2*s12*s13*s22**2*s24*s26*s33*s35**3 - 512*s11*s14*s15*s16*s22**2*s24*s26*s33*s35**3 - 512*s11*s13*s16**2*s22**2*s24*s26*s33*s35**3 -  \
    512*s11*s12*s16**2*s22*s23*s24*s26*s33*s35**3 + 512*s11*s12*s15*s16*s22*s24**2*s26*s33*s35**3 - 1024*s11*s12*s13*s16*s22**2*s25*s26*s33*s35**3 + 512*s11*s14**2*s16*s22**2*s25*s26*s33*s35**3 +  \
    1024*s11*s12**2*s16*s22*s23*s25*s26*s33*s35**3 - 512*s11*s12*s14*s16*s22*s24*s25*s26*s33*s35**3 - 1024*s11*s12*s13*s15*s22**2*s26**2*s33*s35**3 - 512*s11*s13*s14*s16*s22**2*s26**2*s33*s35**3 +  \
    1024*s11*s12**2*s15*s22*s23*s26**2*s33*s35**3 - 1536*s11*s12*s14*s16*s22*s23*s26**2*s33*s35**3 + 512*s11*s12*s14*s15*s22*s24*s26**2*s33*s35**3 + 1536*s11*s12*s13*s16*s22*s24*s26**2*s33*s35**3 +  \
    512*s11*s12**2*s16*s23*s24*s26**2*s33*s35**3 - 512*s11*s12**2*s15*s24**2*s26**2*s33*s35**3 + 1024*s11*s12**2*s13*s22*s25*s26**2*s33*s35**3 - 512*s11*s12*s14**2*s22*s25*s26**2*s33*s35**3 - 1024*s11*s12**3*s23*s25*s26**2*s33*s35**3 +  \
    512*s11*s12**2*s14*s24*s25*s26**2*s33*s35**3 + 512*s11*s12*s13*s14*s22*s26**3*s33*s35**3 + 512*s11*s12**2*s14*s23*s26**3*s33*s35**3 - 1024*s11*s12**2*s13*s24*s26**3*s33*s35**3 - 1024*s11**2*s14*s15*s22**3*s23*s34*s35**3 +  \
    1024*s11**2*s13*s16*s22**3*s23*s34*s35**3 - 1024*s11**2*s12*s16*s22**2*s23**2*s34*s35**3 + 512*s11**2*s13*s15*s22**3*s24*s34*s35**3 + 512*s11**2*s12*s15*s22**2*s23*s24*s34*s35**3 + 512*s11**2*s14*s16*s22**2*s23*s24*s34*s35**3 -  \
    512*s11**2*s13*s16*s22**2*s24**2*s34*s35**3 + 512*s11**2*s13*s14*s22**3*s25*s34*s35**3 + 512*s11**2*s12*s14*s22**2*s23*s25*s34*s35**3 - 1024*s11**2*s12*s13*s22**2*s24*s25*s34*s35**3 - 1024*s11**2*s13**2*s22**3*s26*s34*s35**3 +  \
    1024*s11**2*s12*s13*s22**2*s23*s26*s34*s35**3 - 512*s11**2*s14**2*s22**2*s23*s26*s34*s35**3 + 1024*s11*s14*s15*s16*s22**2*s23*s26*s34*s35**3 - 1024*s11*s13*s16**2*s22**2*s23*s26*s34*s35**3 +  \
    1024*s11*s12*s16**2*s22*s23**2*s26*s34*s35**3 + 512*s11**2*s13*s14*s22**2*s24*s26*s34*s35**3 - 512*s11*s13*s15*s16*s22**2*s24*s26*s34*s35**3 - 512*s11*s12*s15*s16*s22*s23*s24*s26*s34*s35**3 -  \
    512*s11*s14*s16**2*s22*s23*s24*s26*s34*s35**3 + 512*s11*s13*s16**2*s22*s24**2*s26*s34*s35**3 - 512*s11*s13*s14*s16*s22**2*s25*s26*s34*s35**3 - 512*s11*s12*s14*s16*s22*s23*s25*s26*s34*s35**3
v4_3= \
    1024*s11*s12*s13*s16*s22*s24*s25*s26*s34*s35**3 + 1024*s11*s13**2*s16*s22**2*s26**2*s34*s35**3 - 1024*s11*s12*s14*s15*s22*s23*s26**2*s34*s35**3 + 512*s11*s14**2*s16*s22*s23*s26**2*s34*s35**3 -  \
    1024*s11*s12**2*s16*s23**2*s26**2*s34*s35**3 + 512*s11*s12*s13*s15*s22*s24*s26**2*s34*s35**3 - 512*s11*s13*s14*s16*s22*s24*s26**2*s34*s35**3 + 512*s11*s12**2*s15*s23*s24*s26**2*s34*s35**3 +  \
    512*s11*s12*s14*s16*s23*s24*s26**2*s34*s35**3 - 512*s11*s12*s13*s16*s24**2*s26**2*s34*s35**3 + 512*s11*s12*s13*s14*s22*s25*s26**2*s34*s35**3 + 512*s11*s12**2*s14*s23*s25*s26**2*s34*s35**3 -  \
    1024*s11*s12**2*s13*s24*s25*s26**2*s34*s35**3 - 1024*s11*s12*s13**2*s22*s26**3*s34*s35**3 + 1024*s11*s12**2*s13*s23*s26**3*s34*s35**3 - 512*s11*s12*s14**2*s23*s26**3*s34*s35**3 + 512*s11*s12*s13*s14*s24*s26**3*s34*s35**3 +  \
    512*s11**2*s13**2*s22**4*s35**4 - 1024*s11**2*s12*s13*s22**3*s23*s35**4 + 512*s11**2*s14**2*s22**3*s23*s35**4 + 512*s11**2*s12**2*s22**2*s23**2*s35**4 - 512*s11**2*s13*s14*s22**3*s24*s35**4 - 512*s11**2*s12*s14*s22**2*s23*s24*s35**4 +  \
    512*s11**2*s12*s13*s22**2*s24**2*s35**4 - 512*s11*s13**2*s16*s22**3*s26*s35**4 + 1024*s11*s12*s13*s16*s22**2*s23*s26*s35**4 - 512*s11*s14**2*s16*s22**2*s23*s26*s35**4 - 512*s11*s12**2*s16*s22*s23**2*s26*s35**4 +  \
    512*s11*s13*s14*s16*s22**2*s24*s26*s35**4 + 512*s11*s12*s14*s16*s22*s23*s24*s26*s35**4 - 512*s11*s12*s13*s16*s22*s24**2*s26*s35**4 + 512*s11*s12*s13**2*s22**2*s26**2*s35**4 - 1024*s11*s12**2*s13*s22*s23*s26**2*s35**4 +  \
    512*s11*s12*s14**2*s22*s23*s26**2*s35**4 + 512*s11*s12**3*s23**2*s26**2*s35**4 - 512*s11*s12*s13*s14*s22*s24*s26**2*s35**4 - 512*s11*s12**2*s14*s23*s24*s26**2*s35**4 + 512*s11*s12**2*s13*s24**2*s26**2*s35**4 +  \
    1024*s11**3*s16*s22**3*s23*s33**3*s36 + 512*s11**3*s15*s22**3*s24*s33**3*s36 - 512*s11**3*s16*s22**2*s24**2*s33**3*s36 + 1536*s11**3*s14*s22**3*s25*s33**3*s36 - 1536*s11**2*s15*s16*s22**3*s25*s33**3*s36 -  \
    2048*s11**3*s12*s22**2*s24*s25*s33**3*s36 + 1024*s11**2*s16**2*s22**2*s24*s25*s33**3*s36 + 1536*s11**2*s12*s16*s22**2*s25**2*s33**3*s36 - 512*s11*s16**3*s22**2*s25**2*s33**3*s36 + 1024*s11**3*s13*s22**3*s26*s33**3*s36 -  \
    512*s11**2*s15**2*s22**3*s26*s33**3*s36 - 2048*s11**3*s12*s22**2*s23*s26*s33**3*s36 - 1024*s11**2*s16**2*s22**2*s23*s26*s33**3*s36 - 1536*s11**3*s14*s22**2*s24*s26*s33**3*s36 + 512*s11**2*s15*s16*s22**2*s24*s26*s33**3*s36 +  \
    2048*s11**3*s12*s22*s24**2*s26*s33**3*s36 + 2560*s11**2*s12*s15*s22**2*s25*s26*s33**3*s36 - 512*s11**2*s14*s16*s22**2*s25*s26*s33**3*s36 + 512*s11*s15*s16**2*s22**2*s25*s26*s33**3*s36 -  \
    2048*s11**2*s12*s16*s22*s24*s25*s26*s33**3*s36 - 2048*s11**2*s12**2*s22*s25**2*s26*s33**3*s36 + 1024*s11*s12*s16**2*s22*s25**2*s26*s33**3*s36 + 1024*s11**2*s14*s15*s22**2*s26**2*s33**3*s36 -  \
    1024*s11**2*s13*s16*s22**2*s26**2*s33**3*s36 + 3072*s11**2*s12*s16*s22*s23*s26**2*s33**3*s36 - 1536*s11**2*s12*s15*s22*s24*s26**2*s33**3*s36 + 512*s11**2*s14*s16*s22*s24*s26**2*s33**3*s36 -  \
    512*s11**2*s12*s16*s24**2*s26**2*s33**3*s36 - 512*s11**2*s12*s14*s22*s25*s26**2*s33**3*s36 - 1024*s11*s12*s15*s16*s22*s25*s26**2*s33**3*s36 + 2048*s11**2*s12**2*s24*s25*s26**2*s33**3*s36 -  \
    512*s11*s12**2*s16*s25**2*s26**2*s33**3*s36 + 1024*s11**2*s12*s13*s22*s26**3*s33**3*s36 - 512*s11**2*s14**2*s22*s26**3*s33**3*s36 - 2048*s11**2*s12**2*s23*s26**3*s33**3*s36 + 512*s11**2*s12*s14*s24*s26**3*s33**3*s36 +  \
    512*s11*s12**2*s15*s25*s26**3*s33**3*s36 - 1024*s11**3*s15*s22**3*s23*s33**2*s34*s36 - 512*s11**3*s16*s22**2*s23*s24*s33**2*s34*s36 - 512*s11**3*s15*s22**2*s24**2*s33**2*s34*s36 + 512*s11**3*s16*s22*s24**3*s33**2*s34*s36 -  \
    2048*s11**3*s13*s22**3*s25*s33**2*s34*s36 + 1024*s11**2*s15**2*s22**3*s25*s33**2*s34*s36 + 3072*s11**3*s12*s22**2*s23*s25*s33**2*s34*s36 - 512*s11**2*s16**2*s22**2*s23*s25*s33**2*s34*s36 -  \
    512*s11**3*s14*s22**2*s24*s25*s33**2*s34*s36 + 1024*s11**2*s15*s16*s22**2*s24*s25*s33**2*s34*s36 + 1024*s11**3*s12*s22*s24**2*s25*s33**2*s34*s36 - 1024*s11**2*s16**2*s22*s24**2*s25*s33**2*s34*s36 -  \
    2048*s11**2*s12*s15*s22**2*s25**2*s33**2*s34*s36 - 512*s11**2*s14*s16*s22**2*s25**2*s33**2*s34*s36 + 512*s11*s15*s16**2*s22**2*s25**2*s33**2*s34*s36 - 512*s11**2*s12*s16*s22*s24*s25**2*s33**2*s34*s36 +  \
    512*s11*s16**3*s22*s24*s25**2*s33**2*s34*s36 + 1024*s11**2*s12**2*s22*s25**3*s33**2*s34*s36 - 512*s11*s12*s16**2*s22*s25**3*s33**2*s34*s36 + 2560*s11**3*s14*s22**2*s23*s26*s33**2*s34*s36 +  \
    512*s11**2*s15*s16*s22**2*s23*s26*s33**2*s34*s36 + 512*s11**2*s15**2*s22**2*s24*s26*s33**2*s34*s36 - 2048*s11**3*s12*s22*s23*s24*s26*s33**2*s34*s36 + 1024*s11**2*s16**2*s22*s23*s24*s26*s33**2*s34*s36 +  \
    512*s11**3*s14*s22*s24**2*s26*s33**2*s34*s36 - 512*s11**2*s15*s16*s22*s24**2*s26*s33**2*s34*s36 - 1024*s11**3*s12*s24**3*s26*s33**2*s34*s36 - 2048*s11**2*s14*s15*s22**2*s25*s26*s33**2*s34*s36 +  \
    1536*s11**2*s13*s16*s22**2*s25*s26*s33**2*s34*s36 - 512*s11*s15**2*s16*s22**2*s25*s26*s33**2*s34*s36 - 1024*s11**2*s12*s16*s22*s23*s25*s26*s33**2*s34*s36 + 1024*s11**2*s14*s16*s22*s24*s25*s26*s33**2*s34*s36 -  \
    512*s11*s15*s16**2*s22*s24*s25*s26*s33**2*s34*s36 + 1536*s11**2*s12*s16*s24**2*s25*s26*s33**2*s34*s36 + 2560*s11**2*s12*s14*s22*s25**2*s26*s33**2*s34*s36 - 512*s11*s14*s16**2*s22*s25**2*s26*s33**2*s34*s36 -  \
    1024*s11**2*s12**2*s24*s25**2*s26*s33**2*s34*s36 - 512*s11*s12*s16**2*s24*s25**2*s26*s33**2*s34*s36 + 512*s11*s12**2*s16*s25**3*s26*s33**2*s34*s36 - 512*s11**2*s13*s15*s22**2*s26**2*s33**2*s34*s36 -  \
    2048*s11**2*s14*s16*s22*s23*s26**2*s33**2*s34*s36 - 512*s11**2*s14*s15*s22*s24*s26**2*s33**2*s34*s36 + 512*s11**2*s13*s16*s22*s24*s26**2*s33**2*s34*s36 - 512*s11**2*s12*s16*s23*s24*s26**2*s33**2*s34*s36 +  \
    1024*s11**2*s12*s15*s24**2*s26**2*s33**2*s34*s36 - 1024*s11**2*s12*s13*s22*s25*s26**2*s33**2*s34*s36 + 1024*s11**2*s14**2*s22*s25*s26**2*s33**2*s34*s36 + 512*s11*s12*s15**2*s22*s25*s26**2*s33**2*s34*s36 +  \
    512*s11*s14*s15*s16*s22*s25*s26**2*s33**2*s34*s36 + 1024*s11**2*s12**2*s23*s25*s26**2*s33**2*s34*s36 - 2560*s11**2*s12*s14*s24*s25*s26**2*s33**2*s34*s36 + 512*s11*s12*s15*s16*s24*s25*s26**2*s33**2*s34*s36 -  \
    512*s11*s12**2*s15*s25**2*s26**2*s33**2*s34*s36 + 512*s11*s12*s14*s16*s25**2*s26**2*s33**2*s34*s36 + 512*s11**2*s13*s14*s22*s26**3*s33**2*s34*s36 + 1536*s11**2*s12*s14*s23*s26**3*s33**2*s34*s36 -  \
    1024*s11**2*s12*s13*s24*s26**3*s33**2*s34*s36 - 512*s11*s12*s14*s15*s25*s26**3*s33**2*s34*s36 + 1024*s11**3*s16*s22**2*s23**2*s33*s34**2*s36 + 1536*s11**3*s15*s22**2*s23*s24*s33*s34**2*s36 -  \
    1024*s11**3*s16*s22*s23*s24**2*s33*s34**2*s36 - 512*s11**3*s14*s22**2*s23*s25*s33*s34**2*s36 - 1536*s11**2*s15*s16*s22**2*s23*s25*s33*s34**2*s36 + 1024*s11**3*s13*s22**2*s24*s25*s33*s34**2*s36 -  \
    1024*s11**2*s15**2*s22**2*s24*s25*s33*s34**2*s36 - 2048*s11**3*s12*s22*s23*s24*s25*s33*s34**2*s36 + 1536*s11**2*s16**2*s22*s23*s24*s25*s33*s34**2*s36 + 512*s11**2*s15*s16*s22*s24**2*s25*s33*s34**2*s36 +  \
    1024*s11**2*s14*s15*s22**2*s25**2*s33*s34**2*s36 + 512*s11**2*s13*s16*s22**2*s25**2*s33*s34**2*s36 + 1024*s11**2*s12*s16*s22*s23*s25**2*s33*s34**2*s36 - 512*s11*s16**3*s22*s23*s25**2*s33*s34**2*s36 +  \
    1024*s11**2*s12*s15*s22*s24*s25**2*s33*s34**2*s36 - 512*s11**2*s14*s16*s22*s24*s25**2*s33*s34**2*s36 - 512*s11*s15*s16**2*s22*s24*s25**2*s33*s34**2*s36 - 1024*s11**2*s12*s14*s22*s25**3*s33*s34**2*s36 +  \
    512*s11*s14*s16**2*s22*s25**3*s33*s34**2*s36 - 1024*s11**3*s13*s22**2*s23*s26*s33*s34**2*s36 - 512*s11**2*s15**2*s22**2*s23*s26*s33*s34**2*s36 - 1024*s11**2*s16**2*s22*s23**2*s26*s33*s34**2*s36 -  \
    1024*s11**3*s14*s22*s23*s24*s26*s33*s34**2*s36 + 2048*s11**3*s12*s23*s24**2*s26*s33*s34**2*s36 + 1536*s11**2*s13*s15*s22**2*s25*s26*s33*s34**2*s36 + 1024*s11**2*s12*s15*s22*s23*s25*s26*s33*s34**2*s36 +  \
    1024*s11**2*s14*s16*s22*s23*s25*s26*s33*s34**2*s36 + 512*s11*s15*s16**2*s22*s23*s25*s26*s33*s34**2*s36 + 512*s11**2*s14*s15*s22*s24*s25*s26*s33*s34**2*s36 - 2048*s11**2*s13*s16*s22*s24*s25*s26*s33*s34**2*s36 +  \
    512*s11*s15**2*s16*s22*s24*s25*s26*s33*s34**2*s36 - 2048*s11**2*s12*s16*s23*s24*s25*s26*s33*s34**2*s36 - 1024*s11**2*s12*s15*s24**2*s25*s26*s33*s34**2*s36 - 2048*s11**2*s12*s13*s22*s25**2*s26*s33*s34**2*s36 -  \
    512*s11**2*s14**2*s22*s25**2*s26*s33*s34**2*s36 - 512*s11*s14*s15*s16*s22*s25**2*s26*s33*s34**2*s36 + 512*s11*s13*s16**2*s22*s25**2*s26*s33*s34**2*s36 + 512*s11*s12*s16**2*s23*s25**2*s26*s33*s34**2*s36 +  \
    1024*s11**2*s12*s14*s24*s25**2*s26*s33*s34**2*s36 + 512*s11*s12*s15*s16*s24*s25**2*s26*s33*s34**2*s36 - 512*s11*s12*s14*s16*s25**3*s26*s33*s34**2*s36 + 512*s11**2*s14*s15*s22*s23*s26**2*s33*s34**2*s36 +  \
    1024*s11**2*s13*s16*s22*s23*s26**2*s33*s34**2*s36 + 1024*s11**2*s12*s16*s23**2*s26**2*s33*s34**2*s36 - 512*s11**2*s12*s15*s23*s24*s26**2*s33*s34**2*s36 - 1024*s11**2*s13*s14*s22*s25*s26**2*s33*s34**2*s36 -  \
    512*s11*s13*s15*s16*s22*s25*s26**2*s33*s34**2*s36 - 512*s11**2*s12*s14*s23*s25*s26**2*s33*s34**2*s36 - 512*s11*s12*s15*s16*s23*s25*s26**2*s33*s34**2*s36 + 3072*s11**2*s12*s13*s24*s25*s26**2*s33*s34**2*s36 -  \
    512*s11*s12*s15**2*s24*s25*s26**2*s33*s34**2*s36 + 512*s11*s12*s14*s15*s25**2*s26**2*s33*s34**2*s36 - 512*s11*s12*s13*s16*s25**2*s26**2*s33*s34**2*s36 - 1024*s11**2*s12*s13*s23*s26**3*s33*s34**2*s36 +  \
    512*s11*s12*s13*s15*s25*s26**3*s33*s34**2*s36 - 1024*s11**3*s15*s22**2*s23**2*s34**3*s36 + 512*s11**3*s16*s22*s23**2*s24*s34**3*s36 + 1024*s11**2*s15**2*s22**2*s23*s25*s34**3*s36 + 1024*s11**3*s12*s22*s23**2*s25*s34**3*s36 -  \
    512*s11**2*s16**2*s22*s23**2*s25*s34**3*s36 - 512*s11**2*s15*s16*s22*s23*s24*s25*s34**3*s36 - 1024*s11**2*s13*s15*s22**2*s25**2*s34**3*s36 - 1024*s11**2*s12*s15*s22*s23*s25**2*s34**3*s36 +  \
    512*s11*s15*s16**2*s22*s23*s25**2*s34**3*s36 + 512*s11**2*s13*s16*s22*s24*s25**2*s34**3*s36 + 1024*s11**2*s12*s13*s22*s25**3*s34**3*s36 - 512*s11*s13*s16**2*s22*s25**3*s34**3*s36 + 512*s11**3*s14*s22*s23**2*s26*s34**3*s36 +  \
    512*s11**2*s15*s16*s22*s23**2*s26*s34**3*s36 - 1024*s11**3*s12*s23**2*s24*s26*s34**3*s36 - 512*s11**2*s14*s15*s22*s23*s25*s26*s34**3*s36 - 512*s11*s15**2*s16*s22*s23*s25*s26*s34**3*s36 +  \
    512*s11**2*s12*s16*s23**2*s25*s26*s34**3*s36 + 1024*s11**2*s12*s15*s23*s24*s25*s26*s34**3*s36 + 512*s11**2*s13*s14*s22*s25**2*s26*s34**3*s36 + 512*s11*s13*s15*s16*s22*s25**2*s26*s34**3*s36 -  \
    512*s11*s12*s15*s16*s23*s25**2*s26*s34**3*s36 - 1024*s11**2*s12*s13*s24*s25**2*s26*s34**3*s36 + 512*s11*s12*s13*s16*s25**3*s26*s34**3*s36 - 512*s11**2*s12*s15*s23**2*s26**2*s34**3*s36 +  \
    512*s11*s12*s15**2*s23*s25*s26**2*s34**3*s36 - 512*s11*s12*s13*s15*s25**2*s26**2*s34**3*s36 - 2048*s11**3*s14*s22**3*s23*s33**2*s35*s36 - 1024*s11**3*s13*s22**3*s24*s33**2*s35*s36 - 512*s11**2*s15**2*s22**3*s24*s33**2*s35*s36 +  \
    3072*s11**3*s12*s22**2*s23*s24*s33**2*s35*s36 - 512*s11**2*s16**2*s22**2*s23*s24*s33**2*s35*s36 + 1024*s11**3*s14*s22**2*s24**2*s33**2*s35*s36 + 512*s11**2*s15*s16*s22**2*s24**2*s33**2*s35*s36 -  \
    1024*s11**3*s12*s22*s24**3*s33**2*s35*s36 - 512*s11**2*s14*s15*s22**3*s25*s33**2*s35*s36 + 2560*s11**2*s13*s16*s22**3*s25*s33**2*s35*s36 + 512*s11*s15**2*s16*s22**3*s25*s33**2*s35*s36 -  \
    2560*s11**2*s12*s16*s22**2*s23*s25*s33**2*s35*s36 + 512*s11*s16**3*s22**2*s23*s25*s33**2*s35*s36 + 1536*s11**2*s12*s15*s22**2*s24*s25*s33**2*s35*s36 - 2048*s11**2*s14*s16*s22**2*s24*s25*s33**2*s35*s36 -  \
    512*s11*s15*s16**2*s22**2*s24*s25*s33**2*s35*s36 + 1536*s11**2*s12*s16*s22*s24**2*s25*s33**2*s35*s36 + 512*s11**2*s12*s14*s22**2*s25**2*s33**2*s35*s36 - 1024*s11*s12*s15*s16*s22**2*s25**2*s33**2*s35*s36 +  \
    1024*s11*s14*s16**2*s22**2*s25**2*s33**2*s35*s36 - 1024*s11**2*s12**2*s22*s24*s25**2*s33**2*s35*s36 - 512*s11*s12*s16**2*s22*s24*s25**2*s33**2*s35*s36 + 512*s11*s12**2*s16*s22*s25**3*s33**2*s35*s36 -  \
    512*s11**2*s13*s15*s22**3*s26*s33**2*s35*s36 + 512*s11*s15**3*s22**3*s26*s33**2*s35*s36 + 512*s11**2*s12*s15*s22**2*s23*s26*s33**2*s35*s36 + 1536*s11**2*s14*s16*s22**2*s23*s26*s33**2*s35*s36
v4_4= \
    512*s11*s15*s16**2*s22**2*s23*s26*s33**2*s35*s36 + 1024*s11**2*s14*s15*s22**2*s24*s26*s33**2*s35*s36 + 512*s11**2*s13*s16*s22**2*s24*s26*s33**2*s35*s36 - 512*s11*s15**2*s16*s22**2*s24*s26*s33**2*s35*s36 -  \
    1024*s11**2*s12*s16*s22*s23*s24*s26*s33**2*s35*s36 - 1536*s11**2*s12*s15*s22*s24**2*s26*s33**2*s35*s36 - 512*s11**2*s14*s16*s22*s24**2*s26*s33**2*s35*s36 + 512*s11**2*s12*s16*s24**3*s26*s33**2*s35*s36 -  \
    2048*s11**2*s12*s13*s22**2*s25*s26*s33**2*s35*s36 - 512*s11**2*s14**2*s22**2*s25*s26*s33**2*s35*s36 - 2048*s11*s12*s15**2*s22**2*s25*s26*s33**2*s35*s36 + 1024*s11*s14*s15*s16*s22**2*s25*s26*s33**2*s35*s36 -  \
    2048*s11*s13*s16**2*s22**2*s25*s26*s33**2*s35*s36 + 2048*s11**2*s12**2*s22*s23*s25*s26*s33**2*s35*s36 + 2048*s11**2*s12*s14*s22*s24*s25*s26*s33**2*s35*s36 + 1024*s11*s12*s15*s16*s22*s24*s25*s26*s33**2*s35*s36 +  \
    512*s11*s14*s16**2*s22*s24*s25*s26*s33**2*s35*s36 - 1024*s11**2*s12**2*s24**2*s25*s26*s33**2*s35*s36 - 512*s11*s12*s16**2*s24**2*s25*s26*s33**2*s35*s36 + 2560*s11*s12**2*s15*s22*s25**2*s26*s33**2*s35*s36 -  \
    3072*s11*s12*s14*s16*s22*s25**2*s26*s33**2*s35*s36 + 1536*s11*s12**2*s16*s24*s25**2*s26*s33**2*s35*s36 - 1024*s11*s12**3*s25**3*s26*s33**2*s35*s36 - 512*s11**2*s13*s14*s22**2*s26**2*s33**2*s35*s36 -  \
    1024*s11*s14*s15**2*s22**2*s26**2*s33**2*s35*s36 + 1024*s11*s13*s15*s16*s22**2*s26**2*s33**2*s35*s36 - 1024*s11**2*s12*s14*s22*s23*s26**2*s33**2*s35*s36 - 2048*s11*s12*s15*s16*s22*s23*s26**2*s33**2*s35*s36 +  \
    512*s11**2*s14**2*s22*s24*s26**2*s33**2*s35*s36 + 1536*s11*s12*s15**2*s22*s24*s26**2*s33**2*s35*s36 - 512*s11*s14*s15*s16*s22*s24*s26**2*s33**2*s35*s36 + 1024*s11**2*s12**2*s23*s24*s26**2*s33**2*s35*s36 -  \
    512*s11**2*s12*s14*s24**2*s26**2*s33**2*s35*s36 + 512*s11*s12*s15*s16*s24**2*s26**2*s33**2*s35*s36 + 1024*s11*s12*s14*s15*s22*s25*s26**2*s33**2*s35*s36 + 3072*s11*s12*s13*s16*s22*s25*s26**2*s33**2*s35*s36 -  \
    512*s11*s14**2*s16*s22*s25*s26**2*s33**2*s35*s36 - 512*s11*s12**2*s16*s23*s25*s26**2*s33**2*s35*s36 - 2560*s11*s12**2*s15*s24*s25*s26**2*s33**2*s35*s36 + 512*s11*s12*s14*s16*s24*s25*s26**2*s33**2*s35*s36 +  \
    1024*s11*s12**2*s14*s25**2*s26**2*s33**2*s35*s36 - 1024*s11*s12*s13*s15*s22*s26**3*s33**2*s35*s36 + 512*s11*s14**2*s15*s22*s26**3*s33**2*s35*s36 + 1536*s11*s12**2*s15*s23*s26**3*s33**2*s35*s36 -  \
    512*s11*s12*s14*s15*s24*s26**3*s33**2*s35*s36 - 1024*s11*s12**2*s13*s25*s26**3*s33**2*s35*s36 + 4096*s11**3*s13*s22**3*s23*s33*s34*s35*s36 - 4096*s11**3*s12*s22**2*s23**2*s33*s34*s35*s36 -  \
    2048*s11**3*s14*s22**2*s23*s24*s33*s34*s35*s36 + 512*s11**2*s15**2*s22**2*s24**2*s33*s34*s35*s36 + 2048*s11**3*s12*s22*s23*s24**2*s33*s34*s35*s36 + 512*s11**2*s16**2*s22*s23*s24**2*s33*s34*s35*s36 -  \
    512*s11**2*s15*s16*s22*s24**3*s33*s34*s35*s36 - 2048*s11**2*s13*s15*s22**3*s25*s33*s34*s35*s36 + 2048*s11**2*s12*s15*s22**2*s23*s25*s33*s34*s35*s36 + 2048*s11**2*s14*s16*s22**2*s23*s25*s33*s34*s35*s36 +  \
    1024*s11**2*s14*s15*s22**2*s24*s25*s33*s34*s35*s36 - 512*s11*s15**2*s16*s22**2*s24*s25*s33*s34*s35*s36 - 2048*s11**2*s12*s16*s22*s23*s24*s25*s33*s34*s35*s36 - 512*s11*s16**3*s22*s23*s24*s25*s33*s34*s35*s36 -  \
    2048*s11**2*s12*s15*s22*s24**2*s25*s33*s34*s35*s36 + 512*s11**2*s14*s16*s22*s24**2*s25*s33*s34*s35*s36 + 512*s11*s15*s16**2*s22*s24**2*s25*s33*s34*s35*s36 + 2048*s11**2*s12*s13*s22**2*s25**2*s33*s34*s35*s36 -  \
    1536*s11**2*s14**2*s22**2*s25**2*s33*s34*s35*s36 + 512*s11*s14*s15*s16*s22**2*s25**2*s33*s34*s35*s36 - 1536*s11*s13*s16**2*s22**2*s25**2*s33*s34*s35*s36 - 2048*s11**2*s12**2*s22*s23*s25**2*s33*s34*s35*s36 +  \
    1536*s11*s12*s16**2*s22*s23*s25**2*s33*s34*s35*s36 + 2048*s11**2*s12*s14*s22*s24*s25**2*s33*s34*s35*s36 + 512*s11*s12*s15*s16*s22*s24*s25**2*s33*s34*s35*s36 - 512*s11*s14*s16**2*s22*s24*s25**2*s33*s34*s35*s36 -  \
    512*s11*s12*s14*s16*s22*s25**3*s33*s34*s35*s36 - 4096*s11**2*s13*s16*s22**2*s23*s26*s33*s34*s35*s36 + 4096*s11**2*s12*s16*s22*s23**2*s26*s33*s34*s35*s36 - 512*s11*s15**3*s22**2*s24*s26*s33*s34*s35*s36 +  \
    1024*s11**2*s14*s16*s22*s23*s24*s26*s33*s34*s35*s36 - 512*s11*s15*s16**2*s22*s23*s24*s26*s33*s34*s35*s36 - 512*s11**2*s14*s15*s22*s24**2*s26*s33*s34*s35*s36 + 512*s11*s15**2*s16*s22*s24**2*s26*s33*s34*s35*s36 -  \
    2048*s11**2*s12*s16*s23*s24**2*s26*s33*s34*s35*s36 + 1024*s11**2*s12*s15*s24**3*s26*s33*s34*s35*s36 + 2048*s11**2*s13*s14*s22**2*s25*s26*s33*s34*s35*s36 + 512*s11*s14*s15**2*s22**2*s25*s26*s33*s34*s35*s36 +  \
    1024*s11*s13*s15*s16*s22**2*s25*s26*s33*s34*s35*s36 - 4096*s11**2*s12*s14*s22*s23*s25*s26*s33*s34*s35*s36 - 1024*s11*s12*s15*s16*s22*s23*s25*s26*s33*s34*s35*s36 - 512*s11*s14*s16**2*s22*s23*s25*s26*s33*s34*s35*s36 -  \
    2048*s11**2*s12*s13*s22*s24*s25*s26*s33*s34*s35*s36 + 512*s11**2*s14**2*s22*s24*s25*s26*s33*s34*s35*s36 + 1536*s11*s12*s15**2*s22*s24*s25*s26*s33*s34*s35*s36 - 2048*s11*s14*s15*s16*s22*s24*s25*s26*s33*s34*s35*s36 +  \
    1536*s11*s13*s16**2*s22*s24*s25*s26*s33*s34*s35*s36 + 4096*s11**2*s12**2*s23*s24*s25*s26*s33*s34*s35*s36 + 1024*s11*s12*s16**2*s23*s24*s25*s26*s33*s34*s35*s36 - 1024*s11**2*s12*s14*s24**2*s25*s26*s33*s34*s35*s36 -  \
    1536*s11*s12*s14*s15*s22*s25**2*s26*s33*s34*s35*s36 + 2048*s11*s12*s13*s16*s22*s25**2*s26*s33*s34*s35*s36 + 1536*s11*s14**2*s16*s22*s25**2*s26*s33*s34*s35*s36 - 2048*s11*s12**2*s16*s23*s25**2*s26*s33*s34*s35*s36 -  \
    1024*s11*s12**2*s15*s24*s25**2*s26*s33*s34*s35*s36 + 1024*s11*s12**2*s14*s25**3*s26*s33*s34*s35*s36 + 512*s11*s13*s15**2*s22**2*s26**2*s33*s34*s35*s36 + 4096*s11**2*s12*s13*s22*s23*s26**2*s33*s34*s35*s36 -  \
    1536*s11**2*s14**2*s22*s23*s26**2*s33*s34*s35*s36 - 512*s11*s12*s15**2*s22*s23*s26**2*s33*s34*s35*s36 + 1536*s11*s14*s15*s16*s22*s23*s26**2*s33*s34*s35*s36 - 4096*s11**2*s12**2*s23**2*s26**2*s33*s34*s35*s36 +  \
    512*s11*s14*s15**2*s22*s24*s26**2*s33*s34*s35*s36 - 512*s11*s13*s15*s16*s22*s24*s26**2*s33*s34*s35*s36 + 2048*s11**2*s12*s14*s23*s24*s26**2*s33*s34*s35*s36 - 1024*s11*s12*s15**2*s24**2*s26**2*s33*s34*s35*s36 -  \
    2048*s11*s12*s13*s15*s22*s25*s26**2*s33*s34*s35*s36 - 512*s11*s14**2*s15*s22*s25*s26**2*s33*s34*s35*s36 - 512*s11*s13*s14*s16*s22*s25*s26**2*s33*s34*s35*s36 + 2048*s11*s12**2*s15*s23*s25*s26**2*s33*s34*s35*s36 +  \
    2048*s11*s12*s14*s15*s24*s25*s26**2*s33*s34*s35*s36 - 2048*s11*s12*s13*s16*s24*s25*s26**2*s33*s34*s35*s36 - 1024*s11*s12*s14**2*s25**2*s26**2*s33*s34*s35*s36 - 512*s11*s13*s14*s15*s22*s26**3*s33*s34*s35*s36 -  \
    1024*s11*s12*s14*s15*s23*s26**3*s33*s34*s35*s36 + 1024*s11*s12*s13*s15*s24*s26**3*s33*s34*s35*s36 + 1024*s11*s12*s13*s14*s25*s26**3*s33*s34*s35*s36 + 2048*s11**3*s14*s22**2*s23**2*s34**2*s35*s36 -  \
    1024*s11**3*s13*s22**2*s23*s24*s34**2*s35*s36 - 512*s11**2*s15**2*s22**2*s23*s24*s34**2*s35*s36 - 1024*s11**3*s12*s22*s23**2*s24*s34**2*s35*s36 - 512*s11**2*s16**2*s22*s23**2*s24*s34**2*s35*s36 +  \
    512*s11**2*s15*s16*s22*s23*s24**2*s34**2*s35*s36 - 2560*s11**2*s14*s15*s22**2*s23*s25*s34**2*s35*s36 + 512*s11**2*s13*s16*s22**2*s23*s25*s34**2*s35*s36 + 512*s11*s15**2*s16*s22**2*s23*s25*s34**2*s35*s36 -  \
    512*s11**2*s12*s16*s22*s23**2*s25*s34**2*s35*s36 + 512*s11*s16**3*s22*s23**2*s25*s34**2*s35*s36 + 1536*s11**2*s13*s15*s22**2*s24*s25*s34**2*s35*s36 + 2048*s11**2*s12*s15*s22*s23*s24*s25*s34**2*s35*s36 +  \
    512*s11**2*s14*s16*s22*s23*s24*s25*s34**2*s35*s36 - 512*s11*s15*s16**2*s22*s23*s24*s25*s34**2*s35*s36 - 1024*s11**2*s13*s16*s22*s24**2*s25*s34**2*s35*s36 + 1536*s11**2*s13*s14*s22**2*s25**2*s34**2*s35*s36 -  \
    512*s11*s13*s15*s16*s22**2*s25**2*s34**2*s35*s36 + 1024*s11**2*s12*s14*s22*s23*s25**2*s34**2*s35*s36 - 512*s11*s12*s15*s16*s22*s23*s25**2*s34**2*s35*s36 - 512*s11*s14*s16**2*s22*s23*s25**2*s34**2*s35*s36 -  \
    3072*s11**2*s12*s13*s22*s24*s25**2*s34**2*s35*s36 + 1024*s11*s13*s16**2*s22*s24*s25**2*s34**2*s35*s36 + 512*s11*s12*s13*s16*s22*s25**3*s34**2*s35*s36 - 512*s11**2*s13*s15*s22**2*s23*s26*s34**2*s35*s36 +  \
    512*s11*s15**3*s22**2*s23*s26*s34**2*s35*s36 + 512*s11**2*s12*s15*s22*s23**2*s26*s34**2*s35*s36 - 2560*s11**2*s14*s16*s22*s23**2*s26*s34**2*s35*s36 + 512*s11*s15*s16**2*s22*s23**2*s26*s34**2*s35*s36 +  \
    512*s11**2*s14*s15*s22*s23*s24*s26*s34**2*s35*s36 + 2048*s11**2*s13*s16*s22*s23*s24*s26*s34**2*s35*s36 - 512*s11*s15**2*s16*s22*s23*s24*s26*s34**2*s35*s36 + 1536*s11**2*s12*s16*s23**2*s24*s26*s34**2*s35*s36 -  \
    1024*s11**2*s12*s15*s23*s24**2*s26*s34**2*s35*s36 - 1024*s11**2*s13**2*s22**2*s25*s26*s34**2*s35*s36 - 512*s11*s13*s15**2*s22**2*s25*s26*s34**2*s35*s36 + 2048*s11**2*s12*s13*s22*s23*s25*s26*s34**2*s35*s36 +  \
    512*s11**2*s14**2*s22*s23*s25*s26*s34**2*s35*s36 - 1536*s11*s12*s15**2*s22*s23*s25*s26*s34**2*s35*s36 + 2048*s11*s14*s15*s16*s22*s23*s25*s26*s34**2*s35*s36 - 1536*s11*s13*s16**2*s22*s23*s25*s26*s34**2*s35*s36 -  \
    1024*s11**2*s12**2*s23**2*s25*s26*s34**2*s35*s36 - 512*s11*s12*s16**2*s23**2*s25*s26*s34**2*s35*s36 - 1024*s11**2*s13*s14*s22*s24*s25*s26*s34**2*s35*s36 - 1024*s11**2*s12*s14*s23*s24*s25*s26*s34**2*s35*s36 +  \
    2048*s11**2*s12*s13*s24**2*s25*s26*s34**2*s35*s36 + 1536*s11*s12*s13*s15*s22*s25**2*s26*s34**2*s35*s36 - 1536*s11*s13*s14*s16*s22*s25**2*s26*s34**2*s35*s36 + 1024*s11*s12**2*s15*s23*s25**2*s26*s34**2*s35*s36 +  \
    512*s11*s12*s14*s16*s23*s25**2*s26*s34**2*s35*s36 - 512*s11*s12*s13*s16*s24*s25**2*s26*s34**2*s35*s36 - 1024*s11*s12**2*s13*s25**3*s26*s34**2*s35*s36 + 1024*s11**2*s13*s14*s22*s23*s26**2*s34**2*s35*s36 -  \
    512*s11*s14*s15**2*s22*s23*s26**2*s34**2*s35*s36 - 512*s11*s13*s15*s16*s22*s23*s26**2*s34**2*s35*s36 + 1536*s11**2*s12*s14*s23**2*s26**2*s34**2*s35*s36 - 512*s11*s12*s15*s16*s23**2*s26**2*s34**2*s35*s36 -  \
    3072*s11**2*s12*s13*s23*s24*s26**2*s34**2*s35*s36 + 1024*s11*s12*s15**2*s23*s24*s26**2*s34**2*s35*s36 + 512*s11*s13*s14*s15*s22*s25*s26**2*s34**2*s35*s36 + 1024*s11*s13**2*s16*s22*s25*s26**2*s34**2*s35*s36 -  \
    1536*s11*s12*s14*s15*s23*s25*s26**2*s34**2*s35*s36 + 1536*s11*s12*s13*s16*s23*s25*s26**2*s34**2*s35*s36 - 512*s11*s12*s13*s15*s24*s25*s26**2*s34**2*s35*s36 + 1024*s11*s12*s13*s14*s25**2*s26**2*s34**2*s35*s36 +  \
    512*s11*s12*s13*s15*s23*s26**3*s34**2*s35*s36 - 1024*s11*s12*s13**2*s25*s26**3*s34**2*s35*s36 + 1024*s11**2*s14*s15*s22**3*s23*s33*s35**2*s36 - 1024*s11**2*s13*s16*s22**3*s23*s33*s35**2*s36 +  \
    1024*s11**2*s12*s16*s22**2*s23**2*s33*s35**2*s36 + 1536*s11**2*s13*s15*s22**3*s24*s33*s35**2*s36 - 2560*s11**2*s12*s15*s22**2*s23*s24*s33*s35**2*s36 + 1536*s11**2*s14*s16*s22**2*s23*s24*s33*s35**2*s36 -  \
    1024*s11**2*s14*s15*s22**2*s24**2*s33*s35**2*s36 - 512*s11**2*s13*s16*s22**2*s24**2*s33*s35**2*s36 - 1024*s11**2*s12*s16*s22*s23*s24**2*s33*s35**2*s36 + 1024*s11**2*s12*s15*s22*s24**3*s33*s35**2*s36 -  \
    512*s11**2*s13*s14*s22**3*s25*s33*s35**2*s36 - 1024*s11*s13*s15*s16*s22**3*s25*s33*s35**2*s36 - 512*s11**2*s12*s14*s22**2*s23*s25*s33*s35**2*s36 + 1024*s11*s12*s15*s16*s22**2*s23*s25*s33*s35**2*s36 -  \
    1024*s11*s14*s16**2*s22**2*s23*s25*s33*s35**2*s36 - 1024*s11**2*s12*s13*s22**2*s24*s25*s33*s35**2*s36 + 1024*s11**2*s14**2*s22**2*s24*s25*s33*s35**2*s36 + 512*s11*s14*s15*s16*s22**2*s24*s25*s33*s35**2*s36 +  \
    512*s11*s13*s16**2*s22**2*s24*s25*s33*s35**2*s36 + 2048*s11**2*s12**2*s22*s23*s24*s25*s33*s35**2*s36 + 512*s11*s12*s16**2*s22*s23*s24*s25*s33*s35**2*s36 - 1024*s11**2*s12*s14*s22*s24**2*s25*s33*s35**2*s36 -  \
    512*s11*s12*s15*s16*s22*s24**2*s25*s33*s35**2*s36 + 1024*s11*s12*s13*s16*s22**2*s25**2*s33*s35**2*s36 - 512*s11*s14**2*s16*s22**2*s25**2*s33*s35**2*s36 - 1024*s11*s12**2*s16*s22*s23*s25**2*s33*s35**2*s36 +  \
    512*s11*s12*s14*s16*s22*s24*s25**2*s33*s35**2*s36 + 1024*s11**2*s13**2*s22**3*s26*s33*s35**2*s36 - 1024*s11*s13*s15**2*s22**3*s26*s33*s35**2*s36 - 1024*s11**2*s12*s13*s22**2*s23*s26*s33*s35**2*s36 +  \
    512*s11**2*s14**2*s22**2*s23*s26*s33*s35**2*s36 + 1024*s11*s12*s15**2*s22**2*s23*s26*s33*s35**2*s36 - 2048*s11*s14*s15*s16*s22**2*s23*s26*s33*s35**2*s36 + 1024*s11*s13*s16**2*s22**2*s23*s26*s33*s35**2*s36 -  \
    1024*s11*s12*s16**2*s22*s23**2*s26*s33*s35**2*s36 - 1536*s11**2*s13*s14*s22**2*s24*s26*s33*s35**2*s36 + 512*s11*s14*s15**2*s22**2*s24*s26*s33*s35**2*s36 - 1024*s11**2*s12*s14*s22*s23*s24*s26*s33*s35**2*s36 +  \
    2048*s11*s12*s15*s16*s22*s23*s24*s26*s33*s35**2*s36 - 512*s11*s14*s16**2*s22*s23*s24*s26*s33*s35**2*s36 + 2048*s11**2*s12*s13*s22*s24**2*s26*s33*s35**2*s36 - 512*s11*s12*s15**2*s22*s24**2*s26*s33*s35**2*s36
v4_5= \
    512*s11*s14*s15*s16*s22*s24**2*s26*s33*s35**2*s36 + 512*s11*s12*s16**2*s23*s24**2*s26*s33*s35**2*s36 - 512*s11*s12*s15*s16*s24**3*s26*s33*s35**2*s36 + 3072*s11*s12*s13*s15*s22**2*s25*s26*s33*s35**2*s36 -  \
    512*s11*s14**2*s15*s22**2*s25*s26*s33*s35**2*s36 + 1024*s11*s13*s14*s16*s22**2*s25*s26*s33*s35**2*s36 - 3072*s11*s12**2*s15*s22*s23*s25*s26*s33*s35**2*s36 + 3072*s11*s12*s14*s16*s22*s23*s25*s26*s33*s35**2*s36 -  \
    512*s11*s12*s14*s15*s22*s24*s25*s26*s33*s35**2*s36 - 2048*s11*s12*s13*s16*s22*s24*s25*s26*s33*s35**2*s36 - 512*s11*s14**2*s16*s22*s24*s25*s26*s33*s35**2*s36 - 2048*s11*s12**2*s16*s23*s24*s25*s26*s33*s35**2*s36 +  \
    1024*s11*s12**2*s15*s24**2*s25*s26*s33*s35**2*s36 + 512*s11*s12*s14*s16*s24**2*s25*s26*s33*s35**2*s36 - 2048*s11*s12**2*s13*s22*s25**2*s26*s33*s35**2*s36 + 1024*s11*s12*s14**2*s22*s25**2*s26*s33*s35**2*s36 +  \
    2048*s11*s12**3*s23*s25**2*s26*s33*s35**2*s36 - 1024*s11*s12**2*s14*s24*s25**2*s26*s33*s35**2*s36 + 1536*s11*s13*s14*s15*s22**2*s26**2*s33*s35**2*s36 - 1024*s11*s13**2*s16*s22**2*s26**2*s33*s35**2*s36 +  \
    512*s11*s12*s14*s15*s22*s23*s26**2*s33*s35**2*s36 + 512*s11*s14**2*s16*s22*s23*s26**2*s33*s35**2*s36 + 1024*s11*s12**2*s16*s23**2*s26**2*s33*s35**2*s36 - 1536*s11*s12*s13*s15*s22*s24*s26**2*s33*s35**2*s36 -  \
    512*s11*s14**2*s15*s22*s24*s26**2*s33*s35**2*s36 + 512*s11*s13*s14*s16*s22*s24*s26**2*s33*s35**2*s36 - 512*s11*s12**2*s15*s23*s24*s26**2*s33*s35**2*s36 - 512*s11*s12*s14*s16*s23*s24*s26**2*s33*s35**2*s36 +  \
    512*s11*s12*s14*s15*s24**2*s26**2*s33*s35**2*s36 - 512*s11*s12*s13*s16*s24**2*s26**2*s33*s35**2*s36 - 2560*s11*s12*s13*s14*s22*s25*s26**2*s33*s35**2*s36 + 512*s11*s14**3*s22*s25*s26**2*s33*s35**2*s36 -  \
    512*s11*s12**2*s14*s23*s25*s26**2*s33*s35**2*s36 + 3072*s11*s12**2*s13*s24*s25*s26**2*s33*s35**2*s36 - 512*s11*s12*s14**2*s24*s25*s26**2*s33*s35**2*s36 + 1024*s11*s12*s13**2*s22*s26**3*s33*s35**2*s36 -  \
    512*s11*s13*s14**2*s22*s26**3*s33*s35**2*s36 - 1024*s11*s12**2*s13*s23*s26**3*s33*s35**2*s36 + 512*s11*s12*s13*s14*s24*s26**3*s33*s35**2*s36 - 1024*s11**2*s13*s15*s22**3*s23*s34*s35**2*s36 +  \
    1024*s11**2*s12*s15*s22**2*s23**2*s34*s35**2*s36 - 1024*s11**2*s14*s16*s22**2*s23**2*s34*s35**2*s36 + 1536*s11**2*s14*s15*s22**2*s23*s24*s34*s35**2*s36 - 512*s11**2*s13*s16*s22**2*s23*s24*s34*s35**2*s36 +  \
    1536*s11**2*s12*s16*s22*s23**2*s24*s34*s35**2*s36 - 512*s11**2*s13*s15*s22**2*s24**2*s34*s35**2*s36 - 1024*s11**2*s12*s15*s22*s23*s24**2*s34*s35**2*s36 - 512*s11**2*s14*s16*s22*s23*s24**2*s34*s35**2*s36 +  \
    512*s11**2*s13*s16*s22*s24**3*s34*s35**2*s36 + 2048*s11**2*s13**2*s22**3*s25*s34*s35**2*s36 - 3072*s11**2*s12*s13*s22**2*s23*s25*s34*s35**2*s36 + 1536*s11**2*s14**2*s22**2*s23*s25*s34*s35**2*s36 -  \
    1024*s11*s14*s15*s16*s22**2*s23*s25*s34*s35**2*s36 + 1024*s11*s13*s16**2*s22**2*s23*s25*s34*s35**2*s36 + 1024*s11**2*s12**2*s22*s23**2*s25*s34*s35**2*s36 - 1024*s11*s12*s16**2*s22*s23**2*s25*s34*s35**2*s36 -  \
    2560*s11**2*s13*s14*s22**2*s24*s25*s34*s35**2*s36 + 512*s11*s13*s15*s16*s22**2*s24*s25*s34*s35**2*s36 - 2048*s11**2*s12*s14*s22*s23*s24*s25*s34*s35**2*s36 + 512*s11*s12*s15*s16*s22*s23*s24*s25*s34*s35**2*s36 +  \
    512*s11*s14*s16**2*s22*s23*s24*s25*s34*s35**2*s36 + 3072*s11**2*s12*s13*s22*s24**2*s25*s34*s35**2*s36 - 512*s11*s13*s16**2*s22*s24**2*s25*s34*s35**2*s36 + 512*s11*s13*s14*s16*s22**2*s25**2*s34*s35**2*s36 +  \
    512*s11*s12*s14*s16*s22*s23*s25**2*s34*s35**2*s36 - 1024*s11*s12*s13*s16*s22*s24*s25**2*s34*s35**2*s36 + 512*s11**2*s13*s14*s22**2*s23*s26*s34*s35**2*s36 - 1024*s11*s14*s15**2*s22**2*s23*s26*s34*s35**2*s36 +  \
    2048*s11*s13*s15*s16*s22**2*s23*s26*s34*s35**2*s36 + 512*s11**2*s12*s14*s22*s23**2*s26*s34*s35**2*s36 - 2048*s11*s12*s15*s16*s22*s23**2*s26*s34*s35**2*s36 + 1024*s11*s14*s16**2*s22*s23**2*s26*s34*s35**2*s36 +  \
    512*s11*s13*s15**2*s22**2*s24*s26*s34*s35**2*s36 - 512*s11**2*s14**2*s22*s23*s24*s26*s34*s35**2*s36 + 512*s11*s12*s15**2*s22*s23*s24*s26*s34*s35**2*s36 - 512*s11*s13*s16**2*s22*s23*s24*s26*s34*s35**2*s36 -  \
    1024*s11**2*s12**2*s23**2*s24*s26*s34*s35**2*s36 - 512*s11*s12*s16**2*s23**2*s24*s26*s34*s35**2*s36 + 512*s11**2*s13*s14*s22*s24**2*s26*s34*s35**2*s36 - 512*s11*s13*s15*s16*s22*s24**2*s26*s34*s35**2*s36 +  \
    1024*s11**2*s12*s14*s23*s24**2*s26*s34*s35**2*s36 + 512*s11*s12*s15*s16*s23*s24**2*s26*s34*s35**2*s36 - 1024*s11**2*s12*s13*s24**3*s26*s34*s35**2*s36 + 512*s11*s13*s14*s15*s22**2*s25*s26*s34*s35**2*s36 -  \
    2560*s11*s13**2*s16*s22**2*s25*s26*s34*s35**2*s36 + 2560*s11*s12*s14*s15*s22*s23*s25*s26*s34*s35**2*s36 + 1024*s11*s12*s13*s16*s22*s23*s25*s26*s34*s35**2*s36 - 1536*s11*s14**2*s16*s22*s23*s25*s26*s34*s35**2*s36 +  \
    1536*s11*s12**2*s16*s23**2*s25*s26*s34*s35**2*s36 - 2048*s11*s12*s13*s15*s22*s24*s25*s26*s34*s35**2*s36 + 2048*s11*s13*s14*s16*s22*s24*s25*s26*s34*s35**2*s36 - 1024*s11*s12**2*s15*s23*s24*s25*s26*s34*s35**2*s36 -  \
    512*s11*s12*s13*s16*s24**2*s25*s26*s34*s35**2*s36 - 1024*s11*s12*s13*s14*s22*s25**2*s26*s34*s35**2*s36 - 1024*s11*s12**2*s14*s23*s25**2*s26*s34*s35**2*s36 + 2048*s11*s12**2*s13*s24*s25**2*s26*s34*s35**2*s36 -  \
    512*s11*s13**2*s15*s22**2*s26**2*s34*s35**2*s36 - 1024*s11*s12*s13*s15*s22*s23*s26**2*s34*s35**2*s36 + 1024*s11*s14**2*s15*s22*s23*s26**2*s34*s35**2*s36 - 1536*s11*s13*s14*s16*s22*s23*s26**2*s34*s35**2*s36 +  \
    1536*s11*s12**2*s15*s23**2*s26**2*s34*s35**2*s36 - 512*s11*s12*s14*s16*s23**2*s26**2*s34*s35**2*s36 - 512*s11*s13*s14*s15*s22*s24*s26**2*s34*s35**2*s36 + 512*s11*s13**2*s16*s22*s24*s26**2*s34*s35**2*s36 -  \
    1536*s11*s12*s14*s15*s23*s24*s26**2*s34*s35**2*s36 + 1536*s11*s12*s13*s16*s23*s24*s26**2*s34*s35**2*s36 + 1024*s11*s12*s13*s15*s24**2*s26**2*s34*s35**2*s36 + 3072*s11*s12*s13**2*s22*s25*s26**2*s34*s35**2*s36 -  \
    512*s11*s13*s14**2*s22*s25*s26**2*s34*s35**2*s36 - 3072*s11*s12**2*s13*s23*s25*s26**2*s34*s35**2*s36 + 1024*s11*s12*s14**2*s23*s25*s26**2*s34*s35**2*s36 - 512*s11*s12*s13*s14*s24*s25*s26**2*s34*s35**2*s36 +  \
    512*s11*s13**2*s14*s22*s26**3*s34*s35**2*s36 + 512*s11*s12*s13*s14*s23*s26**3*s34*s35**2*s36 - 1024*s11*s12*s13**2*s24*s26**3*s34*s35**2*s36 - 1024*s11**2*s13**2*s22**3*s24*s35**3*s36 +  \
    2048*s11**2*s12*s13*s22**2*s23*s24*s35**3*s36 - 1024*s11**2*s14**2*s22**2*s23*s24*s35**3*s36 - 1024*s11**2*s12**2*s22*s23**2*s24*s35**3*s36 + 1024*s11**2*s13*s14*s22**2*s24**2*s35**3*s36 +  \
    1024*s11**2*s12*s14*s22*s23*s24**2*s35**3*s36 - 1024*s11**2*s12*s13*s22*s24**3*s35**3*s36 + 512*s11*s13**2*s16*s22**3*s25*s35**3*s36 - 1024*s11*s12*s13*s16*s22**2*s23*s25*s35**3*s36 + 512*s11*s14**2*s16*s22**2*s23*s25*s35**3*s36 +  \
    512*s11*s12**2*s16*s22*s23**2*s25*s35**3*s36 - 512*s11*s13*s14*s16*s22**2*s24*s25*s35**3*s36 - 512*s11*s12*s14*s16*s22*s23*s24*s25*s35**3*s36 + 512*s11*s12*s13*s16*s22*s24**2*s25*s35**3*s36 +  \
    512*s11*s13**2*s15*s22**3*s26*s35**3*s36 - 1024*s11*s12*s13*s15*s22**2*s23*s26*s35**3*s36 + 512*s11*s14**2*s15*s22**2*s23*s26*s35**3*s36 + 512*s11*s12**2*s15*s22*s23**2*s26*s35**3*s36 -  \
    512*s11*s13*s14*s15*s22**2*s24*s26*s35**3*s36 + 512*s11*s13**2*s16*s22**2*s24*s26*s35**3*s36 - 512*s11*s12*s14*s15*s22*s23*s24*s26*s35**3*s36 - 1024*s11*s12*s13*s16*s22*s23*s24*s26*s35**3*s36 +  \
    512*s11*s14**2*s16*s22*s23*s24*s26*s35**3*s36 + 512*s11*s12**2*s16*s23**2*s24*s26*s35**3*s36 + 512*s11*s12*s13*s15*s22*s24**2*s26*s35**3*s36 - 512*s11*s13*s14*s16*s22*s24**2*s26*s35**3*s36 -  \
    512*s11*s12*s14*s16*s23*s24**2*s26*s35**3*s36 + 512*s11*s12*s13*s16*s24**3*s26*s35**3*s36 - 1024*s11*s12*s13**2*s22**2*s25*s26*s35**3*s36 + 2048*s11*s12**2*s13*s22*s23*s25*s26*s35**3*s36 -  \
    1024*s11*s12*s14**2*s22*s23*s25*s26*s35**3*s36 - 1024*s11*s12**3*s23**2*s25*s26*s35**3*s36 + 1024*s11*s12*s13*s14*s22*s24*s25*s26*s35**3*s36 + 1024*s11*s12**2*s14*s23*s24*s25*s26*s35**3*s36 -  \
    1024*s11*s12**2*s13*s24**2*s25*s26*s35**3*s36 - 512*s11*s13**2*s14*s22**2*s26**2*s35**3*s36 + 1024*s11*s12*s13*s14*s22*s23*s26**2*s35**3*s36 - 512*s11*s14**3*s22*s23*s26**2*s35**3*s36 - 512*s11*s12**2*s14*s23**2*s26**2*s35**3*s36 +  \
    512*s11*s13*s14**2*s22*s24*s26**2*s35**3*s36 + 512*s11*s12*s14**2*s23*s24*s26**2*s35**3*s36 - 512*s11*s12*s13*s14*s24**2*s26**2*s35**3*s36 - 1024*s11**3*s13*s22**3*s23*s33**2*s36**2 + 512*s11**2*s15**2*s22**3*s23*s33**2*s36**2 +  \
    1024*s11**3*s12*s22**2*s23**2*s33**2*s36**2 + 512*s11**2*s16**2*s22**2*s23**2*s33**2*s36**2 + 1536*s11**3*s14*s22**2*s23*s24*s33**2*s36**2 - 512*s11**2*s15*s16*s22**2*s23*s24*s33**2*s36**2 + 512*s11**3*s13*s22**2*s24**2*s33**2*s36**2 -  \
    2048*s11**3*s12*s22*s23*s24**2*s33**2*s36**2 - 512*s11**3*s14*s22*s24**3*s33**2*s36**2 + 512*s11**3*s12*s24**4*s33**2*s36**2 + 1536*s11**2*s13*s15*s22**3*s25*s33**2*s36**2 - 512*s11*s15**3*s22**3*s25*s33**2*s36**2 -  \
    2560*s11**2*s12*s15*s22**2*s23*s25*s33**2*s36**2 + 512*s11**2*s14*s16*s22**2*s23*s25*s33**2*s36**2 - 512*s11*s15*s16**2*s22**2*s23*s25*s33**2*s36**2 - 512*s11**2*s14*s15*s22**2*s24*s25*s33**2*s36**2 -  \
    2048*s11**2*s13*s16*s22**2*s24*s25*s33**2*s36**2 + 512*s11*s15**2*s16*s22**2*s24*s25*s33**2*s36**2 + 2048*s11**2*s12*s16*s22*s23*s24*s25*s33**2*s36**2 + 512*s11**2*s12*s15*s22*s24**2*s25*s33**2*s36**2 +  \
    1024*s11**2*s14*s16*s22*s24**2*s25*s33**2*s36**2 - 1024*s11**2*s12*s16*s24**3*s25*s33**2*s36**2 - 1536*s11**2*s12*s13*s22**2*s25**2*s33**2*s36**2 + 1536*s11**2*s14**2*s22**2*s25**2*s33**2*s36**2 +  \
    1536*s11*s12*s15**2*s22**2*s25**2*s33**2*s36**2 - 1536*s11*s14*s15*s16*s22**2*s25**2*s33**2*s36**2 + 1536*s11*s13*s16**2*s22**2*s25**2*s33**2*s36**2 + 2048*s11**2*s12**2*s22*s23*s25**2*s33**2*s36**2 -  \
    1024*s11*s12*s16**2*s22*s23*s25**2*s33**2*s36**2 - 2560*s11**2*s12*s14*s22*s24*s25**2*s33**2*s36**2 + 512*s11*s12*s15*s16*s22*s24*s25**2*s33**2*s36**2 - 512*s11*s14*s16**2*s22*s24*s25**2*s33**2*s36**2 +  \
    1024*s11**2*s12**2*s24**2*s25**2*s33**2*s36**2 + 512*s11*s12*s16**2*s24**2*s25**2*s33**2*s36**2 - 1536*s11*s12**2*s15*s22*s25**3*s33**2*s36**2 + 1536*s11*s12*s14*s16*s22*s25**3*s33**2*s36**2 -  \
    1024*s11*s12**2*s16*s24*s25**3*s33**2*s36**2 + 512*s11*s12**3*s25**4*s33**2*s36**2 - 2048*s11**2*s14*s15*s22**2*s23*s26*s33**2*s36**2 + 2048*s11**2*s13*s16*s22**2*s23*s26*s33**2*s36**2 -  \
    3072*s11**2*s12*s16*s22*s23**2*s26*s33**2*s36**2 - 512*s11**2*s13*s15*s22**2*s24*s26*s33**2*s36**2 + 3072*s11**2*s12*s15*s22*s23*s24*s26*s33**2*s36**2 - 1024*s11**2*s14*s16*s22*s23*s24*s26*s33**2*s36**2 +  \
    512*s11**2*s14*s15*s22*s24**2*s26*s33**2*s36**2 + 1024*s11**2*s12*s16*s23*s24**2*s26*s33**2*s36**2 - 512*s11**2*s12*s15*s24**3*s26*s33**2*s36**2 + 512*s11**2*s13*s14*s22**2*s25*s26*s33**2*s36**2 +  \
    1024*s11*s14*s15**2*s22**2*s25*s26*s33**2*s36**2 - 1024*s11*s13*s15*s16*s22**2*s25*s26*s33**2*s36**2 + 1024*s11**2*s12*s14*s22*s23*s25*s26*s33**2*s36**2 + 2048*s11*s12*s15*s16*s22*s23*s25*s26*s33**2*s36**2 +  \
    2048*s11**2*s12*s13*s22*s24*s25*s26*s33**2*s36**2 - 1536*s11**2*s14**2*s22*s24*s25*s26*s33**2*s36**2 - 1536*s11*s12*s15**2*s22*s24*s25*s26*s33**2*s36**2 + 512*s11*s14*s15*s16*s22*s24*s25*s26*s33**2*s36**2 -  \
    4096*s11**2*s12**2*s23*s24*s25*s26*s33**2*s36**2 + 1536*s11**2*s12*s14*s24**2*s25*s26*s33**2*s36**2 - 512*s11*s12*s15*s16*s24**2*s25*s26*s33**2*s36**2 - 512*s11*s12*s14*s15*s22*s25**2*s26*s33**2*s36**2 -  \
    2048*s11*s12*s13*s16*s22*s25**2*s26*s33**2*s36**2 + 512*s11*s14**2*s16*s22*s25**2*s26*s33**2*s36**2 + 1024*s11*s12**2*s16*s23*s25**2*s26*s33**2*s36**2 + 1536*s11*s12**2*s15*s24*s25**2*s26*s33**2*s36**2 -  \
    512*s11*s12*s14*s16*s24*s25**2*s26*s33**2*s36**2 - 512*s11*s12**2*s14*s25**3*s26*s33**2*s36**2 + 512*s11**2*s13**2*s22**2*s26**2*s33**2*s36**2 - 3072*s11**2*s12*s13*s22*s23*s26**2*s33**2*s36**2 +  \
    1536*s11**2*s14**2*s22*s23*s26**2*s33**2*s36**2 + 3072*s11**2*s12**2*s23**2*s26**2*s33**2*s36**2 - 512*s11**2*s13*s14*s22*s24*s26**2*s33**2*s36**2 - 1536*s11**2*s12*s14*s23*s24*s26**2*s33**2*s36**2 +  \
    512*s11**2*s12*s13*s24**2*s26**2*s33**2*s36**2 + 1024*s11*s12*s13*s15*s22*s25*s26**2*s33**2*s36**2 - 512*s11*s14**2*s15*s22*s25*s26**2*s33**2*s36**2 - 1536*s11*s12**2*s15*s23*s25*s26**2*s33**2*s36**2 +  \
    512*s11*s12*s14*s15*s24*s25*s26**2*s33**2*s36**2 + 512*s11*s12**2*s13*s25**2*s26**2*s33**2*s36**2 - 2048*s11**3*s14*s22**2*s23**2*s33*s34*s36**2 - 1024*s11**3*s13*s22**2*s23*s24*s33*s34*s36**2 -  \
    512*s11**2*s15**2*s22**2*s23*s24*s33*s34*s36**2 + 3072*s11**3*s12*s22*s23**2*s24*s33*s34*s36**2 - 512*s11**2*s16**2*s22*s23**2*s24*s33*s34*s36**2 + 1024*s11**3*s14*s22*s23*s24**2*s33*s34*s36**2
v4_6= \
    512*s11**2*s15*s16*s22*s23*s24**2*s33*s34*s36**2 - 1024*s11**3*s12*s23*s24**3*s33*s34*s36**2 + 2560*s11**2*s14*s15*s22**2*s23*s25*s33*s34*s36**2 + 512*s11**2*s13*s16*s22**2*s23*s25*s33*s34*s36**2 -  \
    512*s11**2*s12*s16*s22*s23**2*s25*s33*s34*s36**2 - 512*s11**2*s13*s15*s22**2*s24*s25*s33*s34*s36**2 + 512*s11*s15**3*s22**2*s24*s25*s33*s34*s36**2 - 1024*s11**2*s12*s15*s22*s23*s24*s25*s33*s34*s36**2 -  \
    3072*s11**2*s14*s16*s22*s23*s24*s25*s33*s34*s36**2 + 512*s11*s15*s16**2*s22*s23*s24*s25*s33*s34*s36**2 - 512*s11**2*s14*s15*s22*s24**2*s25*s33*s34*s36**2 + 1024*s11**2*s13*s16*s22*s24**2*s25*s33*s34*s36**2 -  \
    512*s11*s15**2*s16*s22*s24**2*s25*s33*s34*s36**2 + 1536*s11**2*s12*s16*s23*s24**2*s25*s33*s34*s36**2 + 512*s11**2*s12*s15*s24**3*s25*s33*s34*s36**2 - 2560*s11**2*s13*s14*s22**2*s25**2*s33*s34*s36**2 -  \
    512*s11*s14*s15**2*s22**2*s25**2*s33*s34*s36**2 + 512*s11*s13*s15*s16*s22**2*s25**2*s33*s34*s36**2 - 512*s11*s12*s15*s16*s22*s23*s25**2*s33*s34*s36**2 + 1024*s11*s14*s16**2*s22*s23*s25**2*s33*s34*s36**2 +  \
    3072*s11**2*s12*s13*s22*s24*s25**2*s33*s34*s36**2 + 512*s11**2*s14**2*s22*s24*s25**2*s33*s34*s36**2 - 1024*s11*s12*s15**2*s22*s24*s25**2*s33*s34*s36**2 + 1536*s11*s14*s15*s16*s22*s24*s25**2*s33*s34*s36**2 -  \
    1024*s11*s13*s16**2*s22*s24*s25**2*s33*s34*s36**2 - 1024*s11**2*s12**2*s23*s24*s25**2*s33*s34*s36**2 - 512*s11*s12*s16**2*s23*s24*s25**2*s33*s34*s36**2 - 512*s11**2*s12*s14*s24**2*s25**2*s33*s34*s36**2 -  \
    512*s11*s12*s15*s16*s24**2*s25**2*s33*s34*s36**2 + 1024*s11*s12*s14*s15*s22*s25**3*s33*s34*s36**2 - 512*s11*s12*s13*s16*s22*s25**3*s33*s34*s36**2 - 1024*s11*s14**2*s16*s22*s25**3*s33*s34*s36**2 +  \
    512*s11*s12**2*s16*s23*s25**3*s33*s34*s36**2 + 512*s11*s12**2*s15*s24*s25**3*s33*s34*s36**2 + 512*s11*s12*s14*s16*s24*s25**3*s33*s34*s36**2 - 512*s11*s12**2*s14*s25**4*s33*s34*s36**2 +  \
    1536*s11**2*s13*s15*s22**2*s23*s26*s33*s34*s36**2 - 1536*s11**2*s12*s15*s22*s23**2*s26*s33*s34*s36**2 + 2560*s11**2*s14*s16*s22*s23**2*s26*s33*s34*s36**2 - 1024*s11**2*s13*s16*s22*s23*s24*s26*s33*s34*s36**2 -  \
    512*s11**2*s12*s16*s23**2*s24*s26*s33*s34*s36**2 - 512*s11**2*s12*s15*s23*s24**2*s26*s33*s34*s36**2 - 1024*s11**2*s13**2*s22**2*s25*s26*s33*s34*s36**2 - 512*s11*s13*s15**2*s22**2*s25*s26*s33*s34*s36**2 -  \
    512*s11**2*s14**2*s22*s23*s25*s26*s33*s34*s36**2 + 512*s11*s12*s15**2*s22*s23*s25*s26*s33*s34*s36**2 - 1536*s11*s14*s15*s16*s22*s23*s25*s26*s33*s34*s36**2 + 1024*s11**2*s12**2*s23**2*s25*s26*s33*s34*s36**2 +  \
    2048*s11**2*s13*s14*s22*s24*s25*s26*s33*s34*s36**2 - 512*s11*s14*s15**2*s22*s24*s25*s26*s33*s34*s36**2 + 512*s11*s13*s15*s16*s22*s24*s25*s26*s33*s34*s36**2 + 2048*s11**2*s12*s14*s23*s24*s25*s26*s33*s34*s36**2 -  \
    3072*s11**2*s12*s13*s24**2*s25*s26*s33*s34*s36**2 + 1024*s11*s12*s15**2*s24**2*s25*s26*s33*s34*s36**2 + 512*s11*s12*s13*s15*s22*s25**2*s26*s33*s34*s36**2 + 512*s11*s14**2*s15*s22*s25**2*s26*s33*s34*s36**2 -  \
    512*s11*s12**2*s15*s23*s25**2*s26*s33*s34*s36**2 - 512*s11*s12*s14*s16*s23*s25**2*s26*s33*s34*s36**2 - 1536*s11*s12*s14*s15*s24*s25**2*s26*s33*s34*s36**2 + 1536*s11*s12*s13*s16*s24*s25**2*s26*s33*s34*s36**2 +  \
    512*s11*s12*s14**2*s25**3*s26*s33*s34*s36**2 - 1024*s11**2*s13*s14*s22*s23*s26**2*s33*s34*s36**2 - 1536*s11**2*s12*s14*s23**2*s26**2*s33*s34*s36**2 + 2048*s11**2*s12*s13*s23*s24*s26**2*s33*s34*s36**2 +  \
    512*s11*s13*s14*s15*s22*s25*s26**2*s33*s34*s36**2 + 1024*s11*s12*s14*s15*s23*s25*s26**2*s33*s34*s36**2 - 1024*s11*s12*s13*s15*s24*s25*s26**2*s33*s34*s36**2 - 512*s11*s12*s13*s14*s25**2*s26**2*s33*s34*s36**2 +  \
    1024*s11**3*s13*s22**2*s23**2*s34**2*s36**2 + 512*s11**2*s15**2*s22**2*s23**2*s34**2*s36**2 - 1024*s11**3*s12*s22*s23**3*s34**2*s36**2 + 512*s11**2*s16**2*s22*s23**3*s34**2*s36**2 - 512*s11**3*s14*s22*s23**2*s24*s34**2*s36**2 -  \
    512*s11**2*s15*s16*s22*s23**2*s24*s34**2*s36**2 + 512*s11**3*s12*s23**2*s24**2*s34**2*s36**2 - 1024*s11**2*s13*s15*s22**2*s23*s25*s34**2*s36**2 - 512*s11*s15**3*s22**2*s23*s25*s34**2*s36**2 +  \
    1024*s11**2*s14*s16*s22*s23**2*s25*s34**2*s36**2 - 512*s11*s15*s16**2*s22*s23**2*s25*s34**2*s36**2 + 512*s11**2*s14*s15*s22*s23*s24*s25*s34**2*s36**2 + 512*s11*s15**2*s16*s22*s23*s24*s25*s34**2*s36**2 -  \
    512*s11**2*s12*s16*s23**2*s24*s25*s34**2*s36**2 - 512*s11**2*s12*s15*s23*s24**2*s25*s34**2*s36**2 + 1024*s11**2*s13**2*s22**2*s25**2*s34**2*s36**2 + 512*s11*s13*s15**2*s22**2*s25**2*s34**2*s36**2 -  \
    1024*s11**2*s12*s13*s22*s23*s25**2*s34**2*s36**2 + 1024*s11*s12*s15**2*s22*s23*s25**2*s34**2*s36**2 - 1024*s11*s14*s15*s16*s22*s23*s25**2*s34**2*s36**2 + 512*s11*s13*s16**2*s22*s23*s25**2*s34**2*s36**2 +  \
    512*s11**2*s12**2*s23**2*s25**2*s34**2*s36**2 - 512*s11**2*s13*s14*s22*s24*s25**2*s34**2*s36**2 - 512*s11*s13*s15*s16*s22*s24*s25**2*s34**2*s36**2 + 512*s11*s12*s15*s16*s23*s24*s25**2*s34**2*s36**2 +  \
    512*s11**2*s12*s13*s24**2*s25**2*s34**2*s36**2 - 1024*s11*s12*s13*s15*s22*s25**3*s34**2*s36**2 + 1024*s11*s13*s14*s16*s22*s25**3*s34**2*s36**2 - 512*s11*s12**2*s15*s23*s25**3*s34**2*s36**2 -  \
    512*s11*s12*s13*s16*s24*s25**3*s34**2*s36**2 + 512*s11*s12**2*s13*s25**4*s34**2*s36**2 - 512*s11**2*s14*s15*s22*s23**2*s26*s34**2*s36**2 - 512*s11**2*s13*s16*s22*s23**2*s26*s34**2*s36**2 - 512*s11**2*s12*s16*s23**3*s26*s34**2*s36**2 +  \
    1024*s11**2*s12*s15*s23**2*s24*s26*s34**2*s36**2 + 512*s11*s14*s15**2*s22*s23*s25*s26*s34**2*s36**2 + 512*s11*s13*s15*s16*s22*s23*s25*s26*s34**2*s36**2 - 512*s11**2*s12*s14*s23**2*s25*s26*s34**2*s36**2 +  \
    512*s11*s12*s15*s16*s23**2*s25*s26*s34**2*s36**2 - 1024*s11*s12*s15**2*s23*s24*s25*s26*s34**2*s36**2 - 512*s11*s13*s14*s15*s22*s25**2*s26*s34**2*s36**2 - 512*s11*s13**2*s16*s22*s25**2*s26*s34**2*s36**2 +  \
    512*s11*s12*s14*s15*s23*s25**2*s26*s34**2*s36**2 - 512*s11*s12*s13*s16*s23*s25**2*s26*s34**2*s36**2 + 1024*s11*s12*s13*s15*s24*s25**2*s26*s34**2*s36**2 - 512*s11*s12*s13*s14*s25**3*s26*s34**2*s36**2 +  \
    512*s11**2*s12*s13*s23**2*s26**2*s34**2*s36**2 - 512*s11*s12*s13*s15*s23*s25*s26**2*s34**2*s36**2 + 512*s11*s12*s13**2*s25**2*s26**2*s34**2*s36**2 - 1024*s11**2*s13*s15*s22**3*s23*s33*s35*s36**2 +  \
    1024*s11**2*s12*s15*s22**2*s23**2*s33*s35*s36**2 - 1024*s11**2*s14*s16*s22**2*s23**2*s33*s35*s36**2 - 512*s11**2*s14*s15*s22**2*s23*s24*s33*s35*s36**2 + 1536*s11**2*s13*s16*s22**2*s23*s24*s33*s35*s36**2 -  \
    512*s11**2*s12*s16*s22*s23**2*s24*s33*s35*s36**2 - 512*s11**2*s13*s15*s22**2*s24**2*s33*s35*s36**2 + 1024*s11**2*s12*s15*s22*s23*s24**2*s33*s35*s36**2 - 512*s11**2*s14*s16*s22*s23*s24**2*s33*s35*s36**2 +  \
    512*s11**2*s14*s15*s22*s24**3*s33*s35*s36**2 + 512*s11**2*s12*s16*s23*s24**3*s33*s35*s36**2 - 512*s11**2*s12*s15*s24**4*s33*s35*s36**2 - 2048*s11**2*s13**2*s22**3*s25*s33*s35*s36**2 + 1024*s11*s13*s15**2*s22**3*s25*s33*s35*s36**2 +  \
    5120*s11**2*s12*s13*s22**2*s23*s25*s33*s35*s36**2 - 2560*s11**2*s14**2*s22**2*s23*s25*s33*s35*s36**2 - 1024*s11*s12*s15**2*s22**2*s23*s25*s33*s35*s36**2 + 2048*s11*s14*s15*s16*s22**2*s23*s25*s33*s35*s36**2 -  \
    1024*s11*s13*s16**2*s22**2*s23*s25*s33*s35*s36**2 - 3072*s11**2*s12**2*s22*s23**2*s25*s33*s35*s36**2 + 1024*s11*s12*s16**2*s22*s23**2*s25*s33*s35*s36**2 + 2560*s11**2*s13*s14*s22**2*s24*s25*s33*s35*s36**2 -  \
    512*s11*s14*s15**2*s22**2*s24*s25*s33*s35*s36**2 + 3072*s11**2*s12*s14*s22*s23*s24*s25*s33*s35*s36**2 - 2048*s11*s12*s15*s16*s22*s23*s24*s25*s33*s35*s36**2 + 512*s11*s14*s16**2*s22*s23*s24*s25*s33*s35*s36**2 -  \
    2048*s11**2*s12*s13*s22*s24**2*s25*s33*s35*s36**2 - 512*s11**2*s14**2*s22*s24**2*s25*s33*s35*s36**2 + 512*s11*s12*s15**2*s22*s24**2*s25*s33*s35*s36**2 - 512*s11*s14*s15*s16*s22*s24**2*s25*s33*s35*s36**2 -  \
    1024*s11**2*s12**2*s23*s24**2*s25*s33*s35*s36**2 - 512*s11*s12*s16**2*s23*s24**2*s25*s33*s35*s36**2 + 512*s11**2*s12*s14*s24**3*s25*s33*s35*s36**2 + 512*s11*s12*s15*s16*s24**3*s25*s33*s35*s36**2 -  \
    2048*s11*s12*s13*s15*s22**2*s25**2*s33*s35*s36**2 + 512*s11*s14**2*s15*s22**2*s25**2*s33*s35*s36**2 - 512*s11*s13*s14*s16*s22**2*s25**2*s33*s35*s36**2 + 2048*s11*s12**2*s15*s22*s23*s25**2*s33*s35*s36**2 -  \
    1536*s11*s12*s14*s16*s22*s23*s25**2*s33*s35*s36**2 + 512*s11*s12*s13*s16*s22*s24*s25**2*s33*s35*s36**2 + 512*s11*s14**2*s16*s22*s24*s25**2*s33*s35*s36**2 + 1536*s11*s12**2*s16*s23*s24*s25**2*s33*s35*s36**2 -  \
    512*s11*s12**2*s15*s24**2*s25**2*s33*s35*s36**2 - 512*s11*s12*s14*s16*s24**2*s25**2*s33*s35*s36**2 + 1024*s11*s12**2*s13*s22*s25**3*s33*s35*s36**2 - 512*s11*s12*s14**2*s22*s25**3*s33*s35*s36**2 -  \
    1024*s11*s12**3*s23*s25**3*s33*s35*s36**2 + 512*s11*s12**2*s14*s24*s25**3*s33*s35*s36**2 + 512*s11**2*s13*s14*s22**2*s23*s26*s33*s35*s36**2 + 1024*s11*s14*s15**2*s22**2*s23*s26*s33*s35*s36**2 -  \
    1024*s11*s13*s15*s16*s22**2*s23*s26*s33*s35*s36**2 + 512*s11**2*s12*s14*s22*s23**2*s26*s33*s35*s36**2 + 1024*s11*s12*s15*s16*s22*s23**2*s26*s33*s35*s36**2 + 512*s11*s13*s15**2*s22**2*s24*s26*s33*s35*s36**2 -  \
    2048*s11**2*s12*s13*s22*s23*s24*s26*s33*s35*s36**2 + 512*s11**2*s14**2*s22*s23*s24*s26*s33*s35*s36**2 - 1536*s11*s12*s15**2*s22*s23*s24*s26*s33*s35*s36**2 + 512*s11*s14*s15*s16*s22*s23*s24*s26*s33*s35*s36**2 +  \
    1024*s11**2*s12**2*s23**2*s24*s26*s33*s35*s36**2 - 512*s11*s14*s15**2*s22*s24**2*s26*s33*s35*s36**2 - 512*s11**2*s12*s14*s23*s24**2*s26*s33*s35*s36**2 - 512*s11*s12*s15*s16*s23*s24**2*s26*s33*s35*s36**2 +  \
    512*s11*s12*s15**2*s24**3*s26*s33*s35*s36**2 - 3072*s11*s13*s14*s15*s22**2*s25*s26*s33*s35*s36**2 + 2560*s11*s13**2*s16*s22**2*s25*s26*s33*s35*s36**2 - 1024*s11*s12*s14*s15*s22*s23*s25*s26*s33*s35*s36**2 -  \
    2048*s11*s12*s13*s16*s22*s23*s25*s26*s33*s35*s36**2 - 512*s11*s12**2*s16*s23**2*s25*s26*s33*s35*s36**2 + 2048*s11*s12*s13*s15*s22*s24*s25*s26*s33*s35*s36**2 + 1536*s11*s14**2*s15*s22*s24*s25*s26*s33*s35*s36**2 -  \
    1536*s11*s13*s14*s16*s22*s24*s25*s26*s33*s35*s36**2 + 2048*s11*s12**2*s15*s23*s24*s25*s26*s33*s35*s36**2 - 1536*s11*s12*s14*s15*s24**2*s25*s26*s33*s35*s36**2 + 1536*s11*s12*s13*s16*s24**2*s25*s26*s33*s35*s36**2 +  \
    3584*s11*s12*s13*s14*s22*s25**2*s26*s33*s35*s36**2 - 1024*s11*s14**3*s22*s25**2*s26*s33*s35*s36**2 - 512*s11*s12**2*s14*s23*s25**2*s26*s33*s35*s36**2 - 3072*s11*s12**2*s13*s24*s25**2*s26*s33*s35*s36**2 +  \
    1024*s11*s12*s14**2*s24*s25**2*s26*s33*s35*s36**2 - 512*s11*s13**2*s15*s22**2*s26**2*s33*s35*s36**2 + 2048*s11*s12*s13*s15*s22*s23*s26**2*s33*s35*s36**2 - 1024*s11*s14**2*s15*s22*s23*s26**2*s33*s35*s36**2 -  \
    1536*s11*s12**2*s15*s23**2*s26**2*s33*s35*s36**2 + 512*s11*s13*s14*s15*s22*s24*s26**2*s33*s35*s36**2 + 1024*s11*s12*s14*s15*s23*s24*s26**2*s33*s35*s36**2 - 512*s11*s12*s13*s15*s24**2*s26**2*s33*s35*s36**2 -  \
    2048*s11*s12*s13**2*s22*s25*s26**2*s33*s35*s36**2 + 1024*s11*s13*s14**2*s22*s25*s26**2*s33*s35*s36**2 + 2048*s11*s12**2*s13*s23*s25*s26**2*s33*s35*s36**2 - 1024*s11*s12*s13*s14*s24*s25*s26**2*s33*s35*s36**2 -  \
    1024*s11**2*s14*s15*s22**2*s23**2*s34*s35*s36**2 + 1024*s11**2*s13*s16*s22**2*s23**2*s34*s35*s36**2 - 1024*s11**2*s12*s16*s22*s23**3*s34*s35*s36**2 + 1536*s11**2*s13*s15*s22**2*s23*s24*s34*s35*s36**2 -  \
    512*s11**2*s12*s15*s22*s23**2*s24*s34*s35*s36**2 + 1536*s11**2*s14*s16*s22*s23**2*s24*s34*s35*s36**2 - 512*s11**2*s14*s15*s22*s23*s24**2*s34*s35*s36**2 - 1024*s11**2*s13*s16*s22*s23*s24**2*s34*s35*s36**2 -  \
    512*s11**2*s12*s16*s23**2*s24**2*s34*s35*s36**2 + 512*s11**2*s12*s15*s23*s24**3*s34*s35*s36**2 + 512*s11**2*s13*s14*s22**2*s23*s25*s34*s35*s36**2 + 1024*s11*s14*s15**2*s22**2*s23*s25*s34*s35*s36**2 -  \
    2048*s11*s13*s15*s16*s22**2*s23*s25*s34*s35*s36**2 + 512*s11**2*s12*s14*s22*s23**2*s25*s34*s35*s36**2 + 2048*s11*s12*s15*s16*s22*s23**2*s25*s34*s35*s36**2 - 1024*s11*s14*s16**2*s22*s23**2*s25*s34*s35*s36**2 -  \
    1024*s11**2*s13**2*s22**2*s24*s25*s34*s35*s36**2 - 512*s11*s13*s15**2*s22**2*s24*s25*s34*s35*s36**2 - 512*s11**2*s14**2*s22*s23*s24*s25*s34*s35*s36**2 - 512*s11*s12*s15**2*s22*s23*s24*s25*s34*s35*s36**2 +  \
    512*s11*s13*s16**2*s22*s23*s24*s25*s34*s35*s36**2 + 512*s11*s12*s16**2*s23**2*s24*s25*s34*s35*s36**2 + 1024*s11**2*s13*s14*s22*s24**2*s25*s34*s35*s36**2 + 512*s11*s13*s15*s16*s22*s24**2*s25*s34*s35*s36**2 +  \
    512*s11**2*s12*s14*s23*s24**2*s25*s34*s35*s36**2 - 512*s11*s12*s15*s16*s23*s24**2*s25*s34*s35*s36**2 - 1024*s11**2*s12*s13*s24**3*s25*s34*s35*s36**2 - 512*s11*s13*s14*s15*s22**2*s25**2*s34*s35*s36**2
v4_7= \
    1536*s11*s13**2*s16*s22**2*s25**2*s34*s35*s36**2 - 1536*s11*s12*s14*s15*s22*s23*s25**2*s34*s35*s36**2 - 1024*s11*s12*s13*s16*s22*s23*s25**2*s34*s35*s36**2 + 1024*s11*s14**2*s16*s22*s23*s25**2*s34*s35*s36**2 -  \
    512*s11*s12**2*s16*s23**2*s25**2*s34*s35*s36**2 + 1536*s11*s12*s13*s15*s22*s24*s25**2*s34*s35*s36**2 - 1536*s11*s13*s14*s16*s22*s24*s25**2*s34*s35*s36**2 + 512*s11*s12**2*s15*s23*s24*s25**2*s34*s35*s36**2 -  \
    512*s11*s12*s14*s16*s23*s24*s25**2*s34*s35*s36**2 + 1024*s11*s12*s13*s16*s24**2*s25**2*s34*s35*s36**2 + 512*s11*s12*s13*s14*s22*s25**3*s34*s35*s36**2 + 512*s11*s12**2*s14*s23*s25**3*s34*s35*s36**2 -  \
    1024*s11*s12**2*s13*s24*s25**3*s34*s35*s36**2 + 1024*s11**2*s13**2*s22**2*s23*s26*s34*s35*s36**2 - 1024*s11*s13*s15**2*s22**2*s23*s26*s34*s35*s36**2 - 3072*s11**2*s12*s13*s22*s23**2*s26*s34*s35*s36**2 +  \
    1536*s11**2*s14**2*s22*s23**2*s26*s34*s35*s36**2 + 1024*s11*s12*s15**2*s22*s23**2*s26*s34*s35*s36**2 - 1024*s11*s14*s15*s16*s22*s23**2*s26*s34*s35*s36**2 + 2048*s11**2*s12**2*s23**3*s26*s34*s35*s36**2 -  \
    2048*s11**2*s13*s14*s22*s23*s24*s26*s34*s35*s36**2 + 512*s11*s14*s15**2*s22*s23*s24*s26*s34*s35*s36**2 + 512*s11*s13*s15*s16*s22*s23*s24*s26*s34*s35*s36**2 - 2560*s11**2*s12*s14*s23**2*s24*s26*s34*s35*s36**2 +  \
    512*s11*s12*s15*s16*s23**2*s24*s26*s34*s35*s36**2 + 3072*s11**2*s12*s13*s23*s24**2*s26*s34*s35*s36**2 - 512*s11*s12*s15**2*s23*s24**2*s26*s34*s35*s36**2 + 1536*s11*s13**2*s15*s22**2*s25*s26*s34*s35*s36**2 +  \
    1024*s11*s12*s13*s15*s22*s23*s25*s26*s34*s35*s36**2 - 1536*s11*s14**2*s15*s22*s23*s25*s26*s34*s35*s36**2 + 2560*s11*s13*s14*s16*s22*s23*s25*s26*s34*s35*s36**2 - 2560*s11*s12**2*s15*s23**2*s25*s26*s34*s35*s36**2 +  \
    512*s11*s12*s14*s16*s23**2*s25*s26*s34*s35*s36**2 - 1024*s11*s13**2*s16*s22*s24*s25*s26*s34*s35*s36**2 + 2048*s11*s12*s14*s15*s23*s24*s25*s26*s34*s35*s36**2 - 2048*s11*s12*s13*s16*s23*s24*s25*s26*s34*s35*s36**2 -  \
    512*s11*s12*s13*s15*s24**2*s25*s26*s34*s35*s36**2 - 3072*s11*s12*s13**2*s22*s25**2*s26*s34*s35*s36**2 + 1024*s11*s13*s14**2*s22*s25**2*s26*s34*s35*s36**2 + 3072*s11*s12**2*s13*s23*s25**2*s26*s34*s35*s36**2 -  \
    512*s11*s12*s14**2*s23*s25**2*s26*s34*s35*s36**2 - 512*s11*s12*s13*s14*s24*s25**2*s26*s34*s35*s36**2 + 512*s11*s13*s14*s15*s22*s23*s26**2*s34*s35*s36**2 + 512*s11*s12*s14*s15*s23**2*s26**2*s34*s35*s36**2 -  \
    1024*s11*s12*s13*s15*s23*s24*s26**2*s34*s35*s36**2 - 1024*s11*s13**2*s14*s22*s25*s26**2*s34*s35*s36**2 - 1024*s11*s12*s13*s14*s23*s25*s26**2*s34*s35*s36**2 + 2048*s11*s12*s13**2*s24*s25*s26**2*s34*s35*s36**2 +  \
    1024*s11**2*s13**2*s22**3*s23*s35**2*s36**2 - 2048*s11**2*s12*s13*s22**2*s23**2*s35**2*s36**2 + 1024*s11**2*s14**2*s22**2*s23**2*s35**2*s36**2 + 1024*s11**2*s12**2*s22*s23**3*s35**2*s36**2 - 1024*s11**2*s13*s14*s22**2*s23*s24*s35**2*s36**2 -  \
    1024*s11**2*s12*s14*s22*s23**2*s24*s35**2*s36**2 + 512*s11**2*s13**2*s22**2*s24**2*s35**2*s36**2 + 512*s11**2*s14**2*s22*s23*s24**2*s35**2*s36**2 + 512*s11**2*s12**2*s23**2*s24**2*s35**2*s36**2 - 512*s11**2*s13*s14*s22*s24**3*s35**2*s36**2 -  \
    512*s11**2*s12*s14*s23*s24**3*s35**2*s36**2 + 512*s11**2*s12*s13*s24**4*s35**2*s36**2 - 512*s11*s13**2*s15*s22**3*s25*s35**2*s36**2 + 1024*s11*s12*s13*s15*s22**2*s23*s25*s35**2*s36**2 - 512*s11*s14**2*s15*s22**2*s23*s25*s35**2*s36**2 -  \
    512*s11*s12**2*s15*s22*s23**2*s25*s35**2*s36**2 + 512*s11*s13*s14*s15*s22**2*s24*s25*s35**2*s36**2 - 512*s11*s13**2*s16*s22**2*s24*s25*s35**2*s36**2 + 512*s11*s12*s14*s15*s22*s23*s24*s25*s35**2*s36**2 +  \
    1024*s11*s12*s13*s16*s22*s23*s24*s25*s35**2*s36**2 - 512*s11*s14**2*s16*s22*s23*s24*s25*s35**2*s36**2 - 512*s11*s12**2*s16*s23**2*s24*s25*s35**2*s36**2 - 512*s11*s12*s13*s15*s22*s24**2*s25*s35**2*s36**2 +  \
    512*s11*s13*s14*s16*s22*s24**2*s25*s35**2*s36**2 + 512*s11*s12*s14*s16*s23*s24**2*s25*s35**2*s36**2 - 512*s11*s12*s13*s16*s24**3*s25*s35**2*s36**2 + 512*s11*s12*s13**2*s22**2*s25**2*s35**2*s36**2 -  \
    1024*s11*s12**2*s13*s22*s23*s25**2*s35**2*s36**2 + 512*s11*s12*s14**2*s22*s23*s25**2*s35**2*s36**2 + 512*s11*s12**3*s23**2*s25**2*s35**2*s36**2 - 512*s11*s12*s13*s14*s22*s24*s25**2*s35**2*s36**2 -  \
    512*s11*s12**2*s14*s23*s24*s25**2*s35**2*s36**2 + 512*s11*s12**2*s13*s24**2*s25**2*s35**2*s36**2 - 512*s11*s13**2*s16*s22**2*s23*s26*s35**2*s36**2 + 1024*s11*s12*s13*s16*s22*s23**2*s26*s35**2*s36**2 -  \
    512*s11*s14**2*s16*s22*s23**2*s26*s35**2*s36**2 - 512*s11*s12**2*s16*s23**3*s26*s35**2*s36**2 - 512*s11*s13**2*s15*s22**2*s24*s26*s35**2*s36**2 + 1024*s11*s12*s13*s15*s22*s23*s24*s26*s35**2*s36**2 -  \
    512*s11*s14**2*s15*s22*s23*s24*s26*s35**2*s36**2 + 512*s11*s13*s14*s16*s22*s23*s24*s26*s35**2*s36**2 - 512*s11*s12**2*s15*s23**2*s24*s26*s35**2*s36**2 + 512*s11*s12*s14*s16*s23**2*s24*s26*s35**2*s36**2 +  \
    512*s11*s13*s14*s15*s22*s24**2*s26*s35**2*s36**2 + 512*s11*s12*s14*s15*s23*s24**2*s26*s35**2*s36**2 - 512*s11*s12*s13*s16*s23*s24**2*s26*s35**2*s36**2 - 512*s11*s12*s13*s15*s24**3*s26*s35**2*s36**2 +  \
    1024*s11*s13**2*s14*s22**2*s25*s26*s35**2*s36**2 - 2048*s11*s12*s13*s14*s22*s23*s25*s26*s35**2*s36**2 + 1024*s11*s14**3*s22*s23*s25*s26*s35**2*s36**2 + 1024*s11*s12**2*s14*s23**2*s25*s26*s35**2*s36**2 -  \
    1024*s11*s13*s14**2*s22*s24*s25*s26*s35**2*s36**2 - 1024*s11*s12*s14**2*s23*s24*s25*s26*s35**2*s36**2 + 1024*s11*s12*s13*s14*s24**2*s25*s26*s35**2*s36**2 + 512*s11*s13**3*s22**2*s26**2*s35**2*s36**2 -  \
    1024*s11*s12*s13**2*s22*s23*s26**2*s35**2*s36**2 + 512*s11*s13*s14**2*s22*s23*s26**2*s35**2*s36**2 + 512*s11*s12**2*s13*s23**2*s26**2*s35**2*s36**2 - 512*s11*s13**2*s14*s22*s24*s26**2*s35**2*s36**2 -  \
    512*s11*s12*s13*s14*s23*s24*s26**2*s35**2*s36**2 + 512*s11*s12*s13**2*s24**2*s26**2*s35**2*s36**2 + 1024*s11**2*s14*s15*s22**2*s23**2*s33*s36**3 - 1024*s11**2*s13*s16*s22**2*s23**2*s33*s36**3 +  \
    1024*s11**2*s12*s16*s22*s23**3*s33*s36**3 + 512*s11**2*s13*s15*s22**2*s23*s24*s33*s36**3 - 1536*s11**2*s12*s15*s22*s23**2*s24*s33*s36**3 + 512*s11**2*s14*s16*s22*s23**2*s24*s33*s36**3 -  \
    512*s11**2*s14*s15*s22*s23*s24**2*s33*s36**3 - 512*s11**2*s12*s16*s23**2*s24**2*s33*s36**3 + 512*s11**2*s12*s15*s23*s24**3*s33*s36**3 - 512*s11**2*s13*s14*s22**2*s23*s25*s33*s36**3 - 1024*s11*s14*s15**2*s22**2*s23*s25*s33*s36**3 +  \
    1024*s11*s13*s15*s16*s22**2*s23*s25*s33*s36**3 - 512*s11**2*s12*s14*s22*s23**2*s25*s33*s36**3 - 1024*s11*s12*s15*s16*s22*s23**2*s25*s33*s36**3 + 1024*s11**2*s13**2*s22**2*s24*s25*s33*s36**3 -  \
    512*s11*s13*s15**2*s22**2*s24*s25*s33*s36**3 - 2048*s11**2*s12*s13*s22*s23*s24*s25*s33*s36**3 + 1536*s11**2*s14**2*s22*s23*s24*s25*s33*s36**3 + 1536*s11*s12*s15**2*s22*s23*s24*s25*s33*s36**3 -  \
    512*s11*s14*s15*s16*s22*s23*s24*s25*s33*s36**3 + 2048*s11**2*s12**2*s23**2*s24*s25*s33*s36**3 - 1024*s11**2*s13*s14*s22*s24**2*s25*s33*s36**3 + 512*s11*s14*s15**2*s22*s24**2*s25*s33*s36**3 -  \
    1536*s11**2*s12*s14*s23*s24**2*s25*s33*s36**3 + 512*s11*s12*s15*s16*s23*s24**2*s25*s33*s36**3 + 1024*s11**2*s12*s13*s24**3*s25*s33*s36**3 - 512*s11*s12*s15**2*s24**3*s25*s33*s36**3 + 1536*s11*s13*s14*s15*s22**2*s25**2*s33*s36**3 -  \
    1536*s11*s13**2*s16*s22**2*s25**2*s33*s36**3 + 512*s11*s12*s14*s15*s22*s23*s25**2*s33*s36**3 + 2048*s11*s12*s13*s16*s22*s23*s25**2*s33*s36**3 - 512*s11*s14**2*s16*s22*s23*s25**2*s33*s36**3 -  \
    512*s11*s12**2*s16*s23**2*s25**2*s33*s36**3 - 512*s11*s12*s13*s15*s22*s24*s25**2*s33*s36**3 - 1024*s11*s14**2*s15*s22*s24*s25**2*s33*s36**3 + 1024*s11*s13*s14*s16*s22*s24*s25**2*s33*s36**3 -  \
    1536*s11*s12**2*s15*s23*s24*s25**2*s33*s36**3 + 512*s11*s12*s14*s16*s23*s24*s25**2*s33*s36**3 + 1024*s11*s12*s14*s15*s24**2*s25**2*s33*s36**3 - 1024*s11*s12*s13*s16*s24**2*s25**2*s33*s36**3 -  \
    1536*s11*s12*s13*s14*s22*s25**3*s33*s36**3 + 512*s11*s14**3*s22*s25**3*s33*s36**3 + 512*s11*s12**2*s14*s23*s25**3*s33*s36**3 + 1024*s11*s12**2*s13*s24*s25**3*s33*s36**3 - 512*s11*s12*s14**2*s24*s25**3*s33*s36**3 -  \
    1024*s11**2*s13**2*s22**2*s23*s26*s33*s36**3 + 3072*s11**2*s12*s13*s22*s23**2*s26*s33*s36**3 - 1536*s11**2*s14**2*s22*s23**2*s26*s33*s36**3 - 2048*s11**2*s12**2*s23**3*s26*s33*s36**3 + 1024*s11**2*s13*s14*s22*s23*s24*s26*s33*s36**3 +  \
    1536*s11**2*s12*s14*s23**2*s24*s26*s33*s36**3 - 1024*s11**2*s12*s13*s23*s24**2*s26*s33*s36**3 + 512*s11*s13**2*s15*s22**2*s25*s26*s33*s36**3 - 2048*s11*s12*s13*s15*s22*s23*s25*s26*s33*s36**3 +  \
    1024*s11*s14**2*s15*s22*s23*s25*s26*s33*s36**3 + 1536*s11*s12**2*s15*s23**2*s25*s26*s33*s36**3 - 512*s11*s13*s14*s15*s22*s24*s25*s26*s33*s36**3 - 1024*s11*s12*s14*s15*s23*s24*s25*s26*s33*s36**3 +  \
    512*s11*s12*s13*s15*s24**2*s25*s26*s33*s36**3 + 1024*s11*s12*s13**2*s22*s25**2*s26*s33*s36**3 - 512*s11*s13*s14**2*s22*s25**2*s26*s33*s36**3 - 1024*s11*s12**2*s13*s23*s25**2*s26*s33*s36**3 +  \
    512*s11*s12*s13*s14*s24*s25**2*s26*s33*s36**3 - 1024*s11**2*s13*s15*s22**2*s23**2*s34*s36**3 + 1024*s11**2*s12*s15*s22*s23**3*s34*s36**3 - 1024*s11**2*s14*s16*s22*s23**3*s34*s36**3 + 512*s11**2*s14*s15*s22*s23**2*s24*s34*s36**3 +  \
    512*s11**2*s13*s16*s22*s23**2*s24*s34*s36**3 + 512*s11**2*s12*s16*s23**3*s24*s34*s36**3 - 512*s11**2*s12*s15*s23**2*s24**2*s34*s36**3 + 1024*s11*s13*s15**2*s22**2*s23*s25*s34*s36**3 + 1024*s11**2*s12*s13*s22*s23**2*s25*s34*s36**3 -  \
    512*s11**2*s14**2*s22*s23**2*s25*s34*s36**3 - 1024*s11*s12*s15**2*s22*s23**2*s25*s34*s36**3 + 1024*s11*s14*s15*s16*s22*s23**2*s25*s34*s36**3 - 1024*s11**2*s12**2*s23**3*s25*s34*s36**3 -  \
    512*s11*s14*s15**2*s22*s23*s24*s25*s34*s36**3 - 512*s11*s13*s15*s16*s22*s23*s24*s25*s34*s36**3 + 512*s11**2*s12*s14*s23**2*s24*s25*s34*s36**3 - 512*s11*s12*s15*s16*s23**2*s24*s25*s34*s36**3 +  \
    512*s11*s12*s15**2*s23*s24**2*s25*s34*s36**3 - 1024*s11*s13**2*s15*s22**2*s25**2*s34*s36**3 + 512*s11*s14**2*s15*s22*s23*s25**2*s34*s36**3 - 1024*s11*s13*s14*s16*s22*s23*s25**2*s34*s36**3 +  \
    1024*s11*s12**2*s15*s23**2*s25**2*s34*s36**3 + 512*s11*s13*s14*s15*s22*s24*s25**2*s34*s36**3 + 512*s11*s13**2*s16*s22*s24*s25**2*s34*s36**3 - 512*s11*s12*s14*s15*s23*s24*s25**2*s34*s36**3 +  \
    512*s11*s12*s13*s16*s23*s24*s25**2*s34*s36**3 - 512*s11*s12*s13*s15*s24**2*s25**2*s34*s36**3 + 1024*s11*s12*s13**2*s22*s25**3*s34*s36**3 - 512*s11*s13*s14**2*s22*s25**3*s34*s36**3 - 1024*s11*s12**2*s13*s23*s25**3*s34*s36**3 +  \
    512*s11*s12*s13*s14*s24*s25**3*s34*s36**3 + 512*s11**2*s13*s14*s22*s23**2*s26*s34*s36**3 + 512*s11**2*s12*s14*s23**3*s26*s34*s36**3 - 1024*s11**2*s12*s13*s23**2*s24*s26*s34*s36**3 -  \
    512*s11*s13*s14*s15*s22*s23*s25*s26*s34*s36**3 - 512*s11*s12*s14*s15*s23**2*s25*s26*s34*s36**3 + 1024*s11*s12*s13*s15*s23*s24*s25*s26*s34*s36**3 + 512*s11*s13**2*s14*s22*s25**2*s26*s34*s36**3 +  \
    512*s11*s12*s13*s14*s23*s25**2*s26*s34*s36**3 - 1024*s11*s12*s13**2*s24*s25**2*s26*s34*s36**3 - 1024*s11**2*s13**2*s22**2*s23*s24*s35*s36**3 + 2048*s11**2*s12*s13*s22*s23**2*s24*s35*s36**3 -  \
    1024*s11**2*s14**2*s22*s23**2*s24*s35*s36**3 - 1024*s11**2*s12**2*s23**3*s24*s35*s36**3 + 1024*s11**2*s13*s14*s22*s23*s24**2*s35*s36**3 + 1024*s11**2*s12*s14*s23**2*s24**2*s35*s36**3 - 1024*s11**2*s12*s13*s23*s24**3*s35*s36**3 +  \
    512*s11*s13**2*s16*s22**2*s23*s25*s35*s36**3 - 1024*s11*s12*s13*s16*s22*s23**2*s25*s35*s36**3 + 512*s11*s14**2*s16*s22*s23**2*s25*s35*s36**3 + 512*s11*s12**2*s16*s23**3*s25*s35*s36**3 +  \
    512*s11*s13**2*s15*s22**2*s24*s25*s35*s36**3 - 1024*s11*s12*s13*s15*s22*s23*s24*s25*s35*s36**3 + 512*s11*s14**2*s15*s22*s23*s24*s25*s35*s36**3 - 512*s11*s13*s14*s16*s22*s23*s24*s25*s35*s36**3 +  \
    512*s11*s12**2*s15*s23**2*s24*s25*s35*s36**3 - 512*s11*s12*s14*s16*s23**2*s24*s25*s35*s36**3 - 512*s11*s13*s14*s15*s22*s24**2*s25*s35*s36**3 - 512*s11*s12*s14*s15*s23*s24**2*s25*s35*s36**3 +  \
    512*s11*s12*s13*s16*s23*s24**2*s25*s35*s36**3 + 512*s11*s12*s13*s15*s24**3*s25*s35*s36**3 - 512*s11*s13**2*s14*s22**2*s25**2*s35*s36**3 + 1024*s11*s12*s13*s14*s22*s23*s25**2*s35*s36**3 - 512*s11*s14**3*s22*s23*s25**2*s35*s36**3 -  \
    512*s11*s12**2*s14*s23**2*s25**2*s35*s36**3 + 512*s11*s13*s14**2*s22*s24*s25**2*s35*s36**3 + 512*s11*s12*s14**2*s23*s24*s25**2*s35*s36**3 - 512*s11*s12*s13*s14*s24**2*s25**2*s35*s36**3
v4_8= \
    512*s11*s13**2*s15*s22**2*s23*s26*s35*s36**3 - 1024*s11*s12*s13*s15*s22*s23**2*s26*s35*s36**3 + 512*s11*s14**2*s15*s22*s23**2*s26*s35*s36**3 + 512*s11*s12**2*s15*s23**3*s26*s35*s36**3 -  \
    512*s11*s13*s14*s15*s22*s23*s24*s26*s35*s36**3 - 512*s11*s12*s14*s15*s23**2*s24*s26*s35*s36**3 + 512*s11*s12*s13*s15*s23*s24**2*s26*s35*s36**3 - 1024*s11*s13**3*s22**2*s25*s26*s35*s36**3 +  \
    2048*s11*s12*s13**2*s22*s23*s25*s26*s35*s36**3 - 1024*s11*s13*s14**2*s22*s23*s25*s26*s35*s36**3 - 1024*s11*s12**2*s13*s23**2*s25*s26*s35*s36**3 + 1024*s11*s13**2*s14*s22*s24*s25*s26*s35*s36**3 +  \
    1024*s11*s12*s13*s14*s23*s24*s25*s26*s35*s36**3 - 1024*s11*s12*s13**2*s24**2*s25*s26*s35*s36**3 + 512*s11**2*s13**2*s22**2*s23**2*s36**4 - 1024*s11**2*s12*s13*s22*s23**3*s36**4 + 512*s11**2*s14**2*s22*s23**3*s36**4 +  \
    512*s11**2*s12**2*s23**4*s36**4 - 512*s11**2*s13*s14*s22*s23**2*s24*s36**4 - 512*s11**2*s12*s14*s23**3*s24*s36**4 + 512*s11**2*s12*s13*s23**2*s24**2*s36**4 - 512*s11*s13**2*s15*s22**2*s23*s25*s36**4 +  \
    1024*s11*s12*s13*s15*s22*s23**2*s25*s36**4 - 512*s11*s14**2*s15*s22*s23**2*s25*s36**4 - 512*s11*s12**2*s15*s23**3*s25*s36**4 + 512*s11*s13*s14*s15*s22*s23*s24*s25*s36**4 + 512*s11*s12*s14*s15*s23**2*s24*s25*s36**4 -  \
    512*s11*s12*s13*s15*s23*s24**2*s25*s36**4 + 512*s11*s13**3*s22**2*s25**2*s36**4 - 1024*s11*s12*s13**2*s22*s23*s25**2*s36**4 + 512*s11*s13*s14**2*s22*s23*s25**2*s36**4 + 512*s11*s12**2*s13*s23**2*s25**2*s36**4 -  \
    512*s11*s13**2*s14*s22*s24*s25**2*s36**4 - 512*s11*s12*s13*s14*s23*s24*s25**2*s36**4 + 512*s11*s12*s13**2*s24**2*s25**2*s36**4
v4=v4_0+v4_1+v4_2+v4_3+v4_4+v4_5+v4_6+v4_7+v4_8

def R_bn_spherical(phi,theta, Rb_om=None, Rn_om=None):
    # R_BN_SPHERICAL from navigation frame to body frame

    assert(len(phi) == len(theta), "error")
    N = len(phi)
    R = np.zeros((3,3, N))

    for n in range(N):
        t = theta[n]
        p = phi[n]
        if Rb_om is None:
            R[:,:,n] = np.array([[np.cos(t)*np.sin(p), np.sin(t)*np.sin(p), np.cos(p)],
                                 [np.cos(t)*np.cos(p), np.sin(t)*np.cos(p), -np.sin(p)],
                                 [-np.sin(t), np.cos(t), 0]])
        else:
            Rom_n=Rn_om(t).T
            R[:,:,n]=Rb_om(p)*Rom_n
    return R


end

