#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

// Forward declarations for utility functions (assumed to be implemented elsewhere)
void logSO3(double R_input[3][3], double result[3]);
void expSO3(double omega[3], double R_result[3][3]);
void invSO3(double R[3][3], double R_inv[3][3]);
void skew_sym(double v[3], double skew[3][3]);
void rotm2eul(double R[3][3], char* sequence, double euler[3]);
double interp1_linear(double* x, double* v, int n, double xi);

// Helper function to compute matrix norm
double matrix_norm(double* vec, int size) {
    double sum = 0.0;
    for (int i = 0; i < size; i++) {
        sum += vec[i] * vec[i];
    }
    return sqrt(sum);
}

// Helper function to compute mean of columns
void compute_mean_columns(double** matrix, int rows, int cols, double* result) {
    for (int i = 0; i < rows; i++) {
        result[i] = 0.0;
        for (int j = 0; j < cols; j++) {
            result[i] += matrix[i][j];
        }
        result[i] /= cols;
    }
}

// Helper function to matrix multiply 3x3 matrices
void matrix_multiply_3x3(double A[3][3], double B[3][3], double result[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Helper function to transpose 3x3 matrix
void transpose_3x3(double A[3][3], double result[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[j][i] = A[i][j];
        }
    }
}

// Helper function to copy 3x3 matrix
void copy_matrix_3x3(double src[3][3], double dest[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

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
                     double R_mean[3][3], double** list_r) {
    
    // Default parameters
    if (nb_it_max <= 0) nb_it_max = 20;
    if (tol_r <= 0) tol_r = 1e-10;   // [1]
    
    // First approx of R [1]
    copy_matrix_3x3(Rs[0], R_mean);
    
    int nb_it;
    for (nb_it = 1; nb_it <= nb_it_max; nb_it++) { // [2]
        // Initialize list_r with NaN equivalent (using a large value) [3]
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < number_rotations; j++) {
                list_r[i][j] = NAN;
            }
        }
        
        for (int i = 0; i < number_rotations; i++) {
            double R_mean_transpose[3][3];
            double temp_product[3][3];
            double temp_log[3];
            
            transpose_3x3(R_mean, R_mean_transpose);
            matrix_multiply_3x3(R_mean_transpose, Rs[i], temp_product);
            logSO3(temp_product, temp_log);
            
            for (int j = 0; j < 3; j++) {
                list_r[j][i] = temp_log[j];
            }
        }
        
        double r[3];
        compute_mean_columns(list_r, 3, number_rotations, r);
        
        printf("%d/%d: tol: %.3e / %.3e\n", nb_it, nb_it_max, matrix_norm(r, 3), tol_r);
        
        if (matrix_norm(r, 3) < tol_r) { // [4]
            break;
        }
        
        // Update [7]
        double exp_r[3][3];
        double new_R_mean[3][3];
        expSO3(r, exp_r);
        matrix_multiply_3x3(R_mean, exp_r, new_R_mean);
        copy_matrix_3x3(new_R_mean, R_mean);
        
    } // [8]
    
    if (nb_it > nb_it_max) {
        fprintf(stderr, "Error: the maximum number of iteration were reached\n");
        return -1;
    }
    
    return nb_it;
}

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
                         double* t_inter, int n_inter, double R_inter[][3][3]) {
    
    // Find the fractional indices using linear interpolation 
    double* inds_imu_time = (double*)malloc(n_inter * sizeof(double));
    
    for (int i = 0; i < n_inter; i++) {
        inds_imu_time[i] = interp1_linear(t, NULL, n_points, t_inter[i]);
    }
    
    // Initialize R_inter with NaN
    for (int n = 0; n < n_inter; n++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                R_inter[n][i][j] = NAN;
            }
        }
    }
    
    for (int n = 0; n < n_inter; n++) {
        // Extrapolation set to NaN
        if (isnan(inds_imu_time[n])) {
            continue;
        }
        
        double frac = inds_imu_time[n] - floor(inds_imu_time[n]);
        if (frac > 1) {
            printf("Warning: fraction larger than 1\n");
        } else if (frac == 0) {
            // Same point in time
            int round_ind = (int)round(inds_imu_time[n]);
            copy_matrix_3x3(R[round_ind], R_inter[n]);
        } else {
            // Calculate the rotation vector and scale it
            int left_ind = (int)floor(inds_imu_time[n]);
            int right_ind = left_ind + 1;
            
            double R_left_inv[3][3];
            double temp_product[3][3];
            double theta[3];
            double scaled_theta[3];
            double exp_scaled_theta[3][3];
            
            invSO3(R[left_ind], R_left_inv);
            matrix_multiply_3x3(R_left_inv, R[right_ind], temp_product);
            logSO3(temp_product, theta);
            
            for (int i = 0; i < 3; i++) {
                scaled_theta[i] = frac * theta[i];
            }
            
            expSO3(scaled_theta, exp_scaled_theta);
            matrix_multiply_3x3(R[left_ind], exp_scaled_theta, R_inter[n]);
        }
    }
    
    free(inds_imu_time);
}

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
void my_rotm2eul(double R[][3][3], int N, double E[][3]) {
    // rotm2eul gives [yaw, pitch, roll] intrinsic rotation
    // R = R_z(yaw)*R_y(pitch)*R_z(roll)
    // unwrap: adds 2pi when wrapping 
    // flipud to get in order [roll, pitch, yaw]
    
    for (int n = 0; n < N; n++) {
        double euler_temp[3];
        rotm2eul(R[n], "ZYX", euler_temp);
        
        // Flip order to get [roll, pitch, yaw]
        E[n][0] = euler_temp[2];  // roll
        E[n][1] = euler_temp[1];  // pitch
        E[n][2] = euler_temp[0];  // yaw
    }
    
    // Note: unwrap functionality would need to be implemented separately
    // as it requires processing the entire sequence for phase unwrapping
}

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
void R2w_central_diff(double R[][3][3], double* t, int n_points, double w[][3]) {
    
    // Initialize w with NaN
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < n_points; j++) {
            w[i][j] = NAN;
        }
    }
    
    for (int n = 1; n < n_points - 1; n++) {
        double dt = t[n+1] - t[n-1];
        
        double R_prev_transpose[3][3];
        double temp_product[3][3];
        double log_result[3];
        
        transpose_3x3(R[n-1], R_prev_transpose);
        matrix_multiply_3x3(R_prev_transpose, R[n+1], temp_product);
        logSO3(temp_product, log_result);
        
        for (int i = 0; i < 3; i++) {
            w[i][n] = log_result[i] / dt;
        }
    }
}

/**
 * rotationMatrixFromTwoUnitVectors Find rotation matrix from a to b
 * 
 * Parameters:
 *   a: first unit vector [3]
 *   b: second unit vector [3]
 *   R: rotation matrix [3][3] (output)
 */
void rotationMatrixFromTwoUnitVectors(double a[3], double b[3], double R[3][3]) {
    // Normalize vectors
    double norm_a = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    double norm_b = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
    
    double a_norm[3] = {a[0]/norm_a, a[1]/norm_a, a[2]/norm_a};
    double b_norm[3] = {b[0]/norm_b, b[1]/norm_b, b[2]/norm_b};
    
    // Cross product v = a × b
    double v[3];
    v[0] = a_norm[1]*b_norm[2] - a_norm[2]*b_norm[1];
    v[1] = a_norm[2]*b_norm[0] - a_norm[0]*b_norm[2];
    v[2] = a_norm[0]*b_norm[1] - a_norm[1]*b_norm[0];
    
    double s = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    double c = a_norm[0]*b_norm[0] + a_norm[1]*b_norm[1] + a_norm[2]*b_norm[2];
    
    // Create skew symmetric matrices
    double skew_v[3][3];
    double skew_v_squared[3][3];
    
    skew_sym(v, skew_v);
    matrix_multiply_3x3(skew_v, skew_v, skew_v_squared);
    
    // R = I + skew(v) + skew(v)^2 * (1-c)/s^2
    double identity[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
    double factor = (1-c)/(s*s);
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R[i][j] = identity[i][j] + skew_v[i][j] + factor * skew_v_squared[i][j];
        }
    }
} 