import numpy as np
from scipy.optimize import minimize, fmin_bfgs
from scipy.interpolate import UnivariateSpline, splev, splrep
from scipy.misc import derivative
import multiprocessing as mp
from functools import partial
import warnings

def myfminunc(objfun, x0, options=None):
    """
    MATLAB fminunc equivalent using scipy.optimize.minimize
    """
    if options is None:
        options = {}
    
    # Convert MATLAB options to scipy options
    scipy_options = {}
    if 'MaxIter' in options:
        scipy_options['maxiter'] = options['MaxIter']
    if 'TolFun' in options:
        scipy_options['ftol'] = options['TolFun']
    if 'TolX' in options:
        scipy_options['xtol'] = options['TolX']
    
    # Use BFGS method (closest to MATLAB's fminunc default)
    result = minimize(objfun, x0, method='BFGS', options=scipy_options, 
                     jac=True if 'GradObj' in options and options['GradObj'] == 'on' else False)
    
    # Create output structure similar to MATLAB
    out = {
        'x_opt': result.x,
        'fval': result.fun,
        'exitflag': 1 if result.success else 0,
        'output': {
            'iterations': result.nit,
            'funcCount': result.nfev,
            'message': result.message
        },
        'grad': result.jac if hasattr(result, 'jac') else None,
        'hessian': result.hess_inv if hasattr(result, 'hess_inv') else None
    }
    
    return out


def numerical_hessian_forward_diff(f, x):
    """
    Compute numerical Hessian using forward differences
    """
    epsilon = 1e-5
    epsilon_inv = 1 / epsilon
    
    nx = len(x)  # Dimension of the input x
    f0 = f(x)    # calculate f0, when no perturbation happens
    
    f_e = np.zeros(nx)
    g = np.zeros(nx)
    
    # Do perturbation
    for i in range(nx):
        x_ = x.copy()
        x_[i] = x[i] + epsilon
        f_e[i] = f(x_)
        g[i] = (f_e[i] - f0) * epsilon_inv
    
    # Create boolean index for upper triangular matrix
    boolind = np.triu(np.ones((nx, nx), dtype=bool))
    boolind_flat = boolind.flatten()
    
    n_2e = nx * (nx + 1) // 2
    f_2e_vec = np.zeros(n_2e)
    
    # Get indices for upper triangular matrix
    inds = np.triu(np.arange(nx**2).reshape(nx, nx))
    inds_vec = inds[boolind]
    
    for n in range(n_2e):
        ind = inds_vec[n]
        i, j = np.unravel_index(ind, (nx, nx))
        
        x_ = x.copy()
        if i == j:
            x_[i] = x[i] + 2 * epsilon
        else:
            x_[i] = x[i] + epsilon
            x_[j] = x[j] + epsilon
        f_2e_vec[n] = f(x_)
    
    f_2e = np.zeros((nx, nx))
    f_2e[boolind] = f_2e_vec
    
    # Make symmetric
    f_2e = f_2e + np.triu(f_2e, 1).T
    
    h = np.zeros((nx, nx))
    for i in range(nx):
        for j in range(nx):
            h[i, j] = (f_2e[i, j] - f_e[i] - f_e[j] + f0) * epsilon_inv**2
    
    return h, g


def numeric_jacobian(f, x):
    """
    Calculate Jacobian of function f at given x
    """
    epsilon = 1e-6
    epsilon_inv = 1 / epsilon
    nx = len(x)  # Dimension of the input x
    f0 = f(x)    # calculate f0, when no perturbation happens
    nf = len(f0) if hasattr(f0, '__len__') else 1
    if nf == 1:
        nf = 1
        f0 = np.array([f0])
    
    jac = np.zeros((nf, nx))
    
    # Do perturbation
    for i in range(nx):
        x_ = x.copy()
        x_[i] = x[i] + epsilon
        f_pert = f(x_)
        if not hasattr(f_pert, '__len__'):
            f_pert = np.array([f_pert])
        jac[:, i] = (f_pert - f0) * epsilon_inv
    
    return jac


def _hessian_worker(args):
    """Worker function for parallel hessian computation"""
    f, x, epsilon, n, inds_vec, nx = args
    ind = inds_vec[n]
    i, j = np.unravel_index(ind, (nx, nx))
    
    x_ = x.copy()
    if i == j:
        x_[i] = x[i] + 2 * epsilon
    else:
        x_[i] = x[i] + epsilon
        x_[j] = x[j] + epsilon
    return f(x_)


def pnumerical_hessian_forward_diff(f, x):
    """
    Parallel version of numerical Hessian computation using forward differences
    """
    epsilon = 1e-5
    epsilon_inv = 1 / epsilon
    
    nx = len(x)  # Dimension of the input x
    f0 = f(x)    # calculate f0, when no perturbation happens
    
    f_e = np.zeros(nx)
    g = np.zeros(nx)
    
    # Parallel perturbation for gradient
    with mp.Pool() as pool:
        def grad_worker(i):
            x_ = x.copy()
            x_[i] = x[i] + epsilon
            return f(x_)
        
        f_e = pool.map(grad_worker, range(nx))
        f_e = np.array(f_e)
    
    for i in range(nx):
        g[i] = (f_e[i] - f0) * epsilon_inv
    
    # Create boolean index for upper triangular matrix
    boolind = np.triu(np.ones((nx, nx), dtype=bool))
    
    n_2e = nx * (nx + 1) // 2
    
    # Get indices for upper triangular matrix
    inds = np.triu(np.arange(nx**2).reshape(nx, nx))
    inds_vec = inds[boolind]
    
    # Parallel computation of second derivatives
    with mp.Pool() as pool:
        worker_args = [(f, x, epsilon, n, inds_vec, nx) for n in range(n_2e)]
        f_2e_vec = pool.map(_hessian_worker, worker_args)
    
    f_2e_vec = np.array(f_2e_vec)
    f_2e = np.zeros((nx, nx))
    f_2e[boolind] = f_2e_vec
    
    # Make symmetric
    f_2e = f_2e + np.triu(f_2e, 1).T
    
    h = np.zeros((nx, nx))
    for i in range(nx):
        for j in range(nx):
            h[i, j] = (f_2e[i, j] - f_e[i] - f_e[j] + f0) * epsilon_inv**2
    
    return h, g


def runfmincon(objfun, x0, options=None, **kwargs):
    """
    Run constrained optimization with output function
    """
    opt = {}
    if kwargs:
        for key, value in kwargs.items():
            opt[key] = value
    
    # Set up shared variables with outfun
    history = {'x': [], 'fval': []}
    searchdir = []
    x0 = np.array(x0).flatten()
    
    # Call optimization
    callback_data = {'history': history, 'searchdir': searchdir}
    
    def callback_func(xk, *args):
        """Callback function to store optimization history"""
        # This is called after each iteration
        fval = objfun(xk)
        callback_data['history']['x'].append(xk.copy())
        callback_data['history']['fval'].append(fval)
        return False  # Don't stop optimization
    
    if options is None:
        options = {}
    
    scipy_options = {}
    if 'MaxIter' in options:
        scipy_options['maxiter'] = options['MaxIter']
    
    # Check if we need hessian and gradient
    compute_hess_grad = opt.get("hessian_and_grad", False)
    if compute_hess_grad:
        print("Save final hessian and gradient")
    
    # Use scipy.optimize.minimize for constrained optimization
    result = minimize(objfun, x0, method='SLSQP', options=scipy_options, callback=callback_func)
    
    # Convert history to proper format
    if callback_data['history']['x']:
        history['x'] = np.array(callback_data['history']['x']).T
    else:
        history['x'] = np.array([]).reshape(len(x0), 0)
    
    out = {
        'x_opt': result.x,
        'fval': result.fun,
        'exitflag': 1 if result.success else 0,
        'output': {
            'iterations': result.nit,
            'funcCount': result.nfev,
            'message': result.message
        },
        'history': history,
        'searchdir': searchdir
    }
    
    if compute_hess_grad:
        # Compute numerical gradient and hessian at optimal point
        h, g = numerical_hessian_forward_diff(objfun, result.x)
        out['grad'] = g
        out['hessian'] = h
    
    return out


def solveHessian(test_function, a):
    """
    Objective: Generates Hessian of a function at some point
    -----------------------------------------------------------------------
    hf = solveHessian(test_function, a)
    where a = input vector
          test_function = objective function
    -----------------------------------------------------------------------
    Output: hf = Hessian matrix
    -----------------------------------------------------------------------
    
    Code by:
    Salil Sharma
    May 3, 2017
    -----------------------------------------------------------------------
    """
    
    l = len(a)  # Hessian would be l x l matrix
    ep = 0.0001  # step size for numerical differentiation
    valf = test_function(a)  # value of obj function at a
    ep2 = ep * ep
    ep3 = 4 * ep * ep
    hf = np.zeros((l, l))
    
    for i in range(l):
        x1 = a.copy()
        x1[i] = a[i] - ep  # Change ith element in x1
        x2 = a.copy()
        x2[i] = a[i] + ep  # Change ith element in x2
        hf[i, i] = (test_function(x2) - 2*valf + test_function(x1)) / ep2  # diagonal entries
        
        j = i + 1
        while j < l:  # Loop computes the rest of the elements of the Hessian matrix
            x1[j] = a[j] - ep  # Lower the value of step size
            x2[j] = a[j] + ep  # Increment the value of step size
            v4 = test_function(x1)  # compute the respective values
            v1 = test_function(x2)  # compute the respective values
            x1[j] = x1[j] + 2*ep
            x2[j] = x2[j] - 2*ep
            v2 = test_function(x1)
            v3 = test_function(x2)
            hf[i, j] = (v1 + v4 - v2 - v3) / ep3
            hf[j, i] = hf[i, j]  # d2f/dxdy is same as that of d2f/dydx
            x1[j] = a[j]
            x2[j] = a[j]
            j += 1
    
    return hf


def spline_derivative(t, w, j=1):
    """
    Compute spline derivatives
    """
    if j < 1:
        j = 1
    
    w = np.array(w)
    d = np.zeros_like(w)
    
    for i in range(3):
        # Create spline representation
        tck = splrep(t, w[i, :], s=0)  # s=0 for interpolation
        
        # Compute derivative
        tck_deriv = tck
        for _ in range(j):
            # Compute derivative of spline
            knots, coeffs, degree = tck_deriv
            if degree == 0:
                d[i, :] = np.zeros_like(t)
                break
            # Derivative reduces degree by 1
            new_coeffs = []
            for k in range(len(coeffs)):
                if k < len(coeffs) - 1:
                    new_coeffs.append(coeffs[k] * degree)
            tck_deriv = (knots[1:-1], new_coeffs, degree - 1)
        
        if degree > 0:
            d[i, :] = splev(t, tck_deriv)
    
    return d
