#ifndef SO3_H
#define SO3_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * SO(3) Library - C Implementation
 * Converted from MATLAB code for Special Orthogonal Group operations
 */

// Helper function declarations
double matrix_norm(double matrix[3][3]);
double vector_norm(double v[3]);
double matrix_trace(double matrix[3][3]);
void matrix_multiply(double A[3][3], double B[3][3], double result[3][3]);
void set_identity(double matrix[3][3]);
void copy_matrix(double src[3][3], double dest[3][3]);

/**
 * Adjoint representation of SO(3)
 * @param w: 3D vector
 * @param adw: Output 3x3 matrix (adjoint map)
 */
void adjSO3(double w[3], double adw[3][3]);

/**
 * Adjoint action of SO(3)
 * @param R: 3x3 rotation matrix
 * @param AdR: Output 3x3 matrix (Adjoint map)
 */
void AdSO3(double R[3][3], double AdR[3][3]);

/**
 * Dual hat operator for SO(3)
 * w_hat*x = x_dual_hat*w
 * @param x: 3D vector
 * @param x_dual_hat: Output 3x3 skew-symmetric matrix
 */
void dualHatSO3(double x[3], double x_dual_hat[3][3]);

/**
 * Calculate rotation matrix errors
 * UNTITLED Summary of this function goes here
 * Need to be same length
 * @param Rhat: Array of estimated rotation matrices
 * @param Rtrue: True rotation matrix (or array)
 * @param err: Output error vectors
 * @param N: Number of matrices
 * @param Rtrue_is_array: Flag indicating if Rtrue is an array
 */
void errorSO3(double Rhat[][3][3], double Rtrue[3][3], double err[][3], int N, int Rtrue_is_array);

/**
 * Exponential map from so(3) to SO(3)
 * @param w: 3D axis-angle vector
 * @param R: Output 3x3 rotation matrix
 */
void expSO3(double w[3], double R[3][3]);

/**
 * Hat operator - converts vector to skew-symmetric matrix
 * @param w: 3D vector
 * @param w_hat: Output 3x3 skew-symmetric matrix
 */
void HatSO3(double w[3], double w_hat[3][3]);

/**
 * Inverse of rotation matrix (transpose)
 * @param R: Input 3x3 rotation matrix
 * @param Rinv: Output 3x3 inverse matrix
 * @param errorFlag: Optional error flag pointer
 */
void invSO3(double R[3][3], double Rinv[3][3], int* errorFlag);

/**
 * Logarithm map from SO(3) to so(3)
 * @param R: 3x3 rotation matrix
 * @param w: Output 3D axis-angle vector
 * @param errorflag: Optional error flag pointer
 */
void logSO3(double R[3][3], double w[3], int* errorflag);

/**
 * Normalize rotation matrix using simplified SVD approach
 * [u,s,v] = svd(R);
 * Rnorm = u*v';
 * @param R: Input 3x3 matrix
 * @param Rnorm: Output normalized 3x3 rotation matrix
 */
void normalizeSO3(double R[3][3], double Rnorm[3][3]);

/**
 * Left-Jacobian to SO(3)
 * sum_k 1/(k + 1)! ad(w)^k
 * @param w: 3D vector
 * @param Phiw: Output 3x3 Jacobian matrix
 */
void PhiSO3(double w[3], double Phiw[3][3]);

/**
 * Vector extraction from skew-symmetric matrix
 * @param w_hat: 3x3 skew-symmetric matrix
 * @param w: Output 3D vector
 */
void VecSO3(double w_hat[3][3], double w[3]);

#endif // SO3_H 