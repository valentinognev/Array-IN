#ifndef ROTATION_UTILS_H
#define ROTATION_UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
 * average_rotation Average rotation from multiple rotation matrices
 * 
 * Parameters:
 *   Rs: Array of rotation matrices [3][3][N], N number of rotations
 *   number_rotations: N number of rotations
 *   nb_it_max: max number of iterations, default 20 
 *   tol_r: tolerance for deviations, default 1e-10
 *   R_mean: Average rotation (output)
 *   list_r: residuals in lie algebra (output) [3][N]
 * 
 * Returns: number of iterations used, or -1 if max iterations reached
 */
int average_rotation(double Rs[][3][3], int number_rotations, int nb_it_max, double tol_r, 
                     double R_mean[3][3], double** list_r);

/**
 * INTERPOLATE_POS_AND_ROTATION Interpolate IMU pos and rotation estimates 
 * 
 * Parameters:
 *   t: time points (input data)
 *   R: rotation matrices [3][3][n_points]
 *   n_points: number of input data points
 *   t_inter: time points where interpolation should occur
 *   n_inter: number of interpolation points
 *   R_inter: interpolated rotation matrices [3][3][n_inter] (output)
 * 
 * t and R are data points 
 * t_inter is the time points where interpolation should occur
 * R_inter is the interpolated rotation matrix
 */
void interpolate_rotation(double* t, double R[][3][3], int n_points,
                         double* t_inter, int n_inter, double R_inter[][3][3]);

/**
 * MY_ROTM2EUL Rotation matrix to euler angles [roll, pitch, yaw]
 *
 * Parameters:
 *   R: rotation matrix [3][3][N]
 *   N: number of rotation matrices
 *   E: euler angles [3][N] (output)
 *
 *   Roll: around x-axis
 *   Pitch: around y-axis
 *   Yaw: around z-axis (heading)
 *   R (3,3,N) ->  E (3, N) 
 */
void my_rotm2eul(double R[][3][3], int N, double E[][3]);

/**
 * R2W_CENTRAL_DIFF Rotation matrix 2 angular velocity using central 
 * difference 
 * 
 * Parameters:
 *   R: rotation matrices [3][3][n_points]
 *   t: time vector
 *   n_points: number of time points
 *   w: angular velocity [3][n_points] (output)
 *
 *   Based on:
 *   R_{t+1} = R_{t} exp_SO3(w*t)
 *   w in body frame
 */
void R2w_central_diff(double R[][3][3], double* t, int n_points, double w[][3]);

/**
 * rotationMatrixFromTwoUnitVectors Find rotation matrix from a to b
 * 
 * Parameters:
 *   a: first unit vector [3]
 *   b: second unit vector [3]
 *   R: rotation matrix [3][3] (output)
 */
void rotationMatrixFromTwoUnitVectors(double a[3], double b[3], double R[3][3]);

// Helper functions
double matrix_norm(double* vec, int size);
void compute_mean_columns(double** matrix, int rows, int cols, double* result);
void matrix_multiply_3x3(double A[3][3], double B[3][3], double result[3][3]);
void transpose_3x3(double A[3][3], double result[3][3]);
void copy_matrix_3x3(double src[3][3], double dest[3][3]);

// Forward declarations for utility functions (need to be implemented elsewhere)
void logSO3(double R_input[3][3], double result[3]);
void expSO3(double omega[3], double R_result[3][3]);
void invSO3(double R[3][3], double R_inv[3][3]);
void skew_sym(double v[3], double skew[3][3]);
void rotm2eul(double R[3][3], char* sequence, double euler[3]);
double interp1_linear(double* x, double* v, int n, double xi);

#endif // ROTATION_UTILS_H 