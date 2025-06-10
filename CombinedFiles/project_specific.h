#ifndef PROJECT_SPECIFIC_H
#define PROJECT_SPECIFIC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

// Data structure definitions
typedef struct {
    double **data;
    int rows;
    int cols;
} Matrix2D;

typedef struct {
    double ***data;
    int rows;
    int cols;
    int depth;
} Matrix3D;

typedef struct {
    double *data;
    int length;
} Vector;

// State structure - equivalent to MATLAB struct S
typedef struct {
    Matrix3D *R;        // Rotation matrices
    Matrix2D *w;        // Angular velocity
    Matrix2D *p;        // Position
    Matrix2D *v;        // Velocity
    Matrix2D *b_a;      // Accelerometer bias
    Matrix2D *b_g;      // Gyroscope bias
    Matrix2D *omega_dot; // Angular acceleration
    Matrix2D *v_dot;    // Linear acceleration
    Matrix2D *s;        // Scale factors
    Matrix2D *b_s;      // Scale bias
    Matrix3D *T_a;      // Transformation matrices
    Matrix2D *b_omega_dot; // Angular acceleration bias
    
    // For interpolation results
    Matrix2D *p_rig_time;
    Matrix3D *R_rig_time;
} State;

// Error structure
typedef struct {
    Vector *R;          // Rotation error
    Vector *R_deg;      // Rotation error in degrees
    Matrix2D *w;        // Angular velocity error
    Matrix2D *w_deg;    // Angular velocity error in degrees
    Matrix2D *p;        // Position error
    Matrix2D *v;        // Velocity error
    Matrix2D *b_a;      // Accelerometer bias error
    Matrix2D *b_g;      // Gyroscope bias error
    Matrix2D *b_g_deg;  // Gyroscope bias error in degrees
    Matrix2D *omega_dot; // Angular acceleration error
    Matrix2D *v_dot;    // Linear acceleration error
    Matrix2D *s;        // Scale factor error
    Matrix2D *b_s;      // Scale bias error
    Matrix3D *T_a;      // Transformation matrix error
    Matrix2D *b_omega_dot; // Angular acceleration bias error
} ErrorStruct;

// Measurement structure
typedef struct {
    Matrix2D *y;        // Measurements
    Matrix2D *Q;        // Covariance matrix
    Matrix2D *Q_inv;    // Inverse covariance matrix
} Measurements;

// Filter results structure
typedef struct {
    State *filt;        // Filtered state
    State *pred;        // Predicted state
} FilterResults;

// Complete result structure
typedef struct {
    FilterResults *results;
    struct {
        ErrorStruct *filt;
        ErrorStruct *pred;
    } err;
} CompleteResults;

// Release indices structure
typedef struct {
    Matrix2D *inds_growth;
    Vector *IN_time_array;
} ReleaseIndices;

// Function declarations

// Main trajectory error calculation function
ErrorStruct* calculate_trajectory_error(State *S_hat, State *S_true, int *section, int section_size);

// Covariance compensation function
Matrix2D* compensate_covariance(Matrix2D *Q_y, Matrix2D *T);

// Measurement compensation function
Vector* compensate_measurements(Vector *y, Matrix2D *T, Vector *b);

// General error computation function
ErrorStruct* compute_error(State *S, State *S_ref);

// T and b estimation function
void estimate_T_and_b(Matrix2D *y, Matrix2D *u, Matrix2D *Q, Matrix2D **T_out, Vector **b_out);

// Release indices calculation
ReleaseIndices* get_release_inds(Vector *time, Vector *release_times, double IN_time, double T);

// Position and rotation interpolation
void interpolate_pos_and_rotation(State *S, Vector *imu_time, Vector *rig_time);

// Least squares triad functions
void lsq_triad(Vector *y, Matrix2D *Q, Vector **u_out, Matrix2D **Qu_out);
void lsq_triad_naive(Vector *y, Matrix2D *Q, Vector **u_out, Matrix2D **Qu_out);

// Measurement rotation function
Measurements* rotate_measurements(Measurements *S_in, Matrix2D *R);

// Filter running functions
CompleteResults* run_filter(Measurements *sensorData, State *initData, void *my_settings, State *S_ref, void *myFilter);
CompleteResults* run_filter_w_error(Measurements *sensorData, State *initData, void *my_settings, State *S_ref, void *myFilter);

// Utility functions for matrix operations (these would need to be implemented or use external library)
double errorSO3(Matrix2D *R1, Matrix2D *R2);
Vector* errorSO3_vector(Matrix3D *R_hat, Matrix3D *R_true);
double rad2deg(double rad);
Vector* rad2deg_vector(Vector *rad_vec);
Matrix2D* rad2deg_matrix(Matrix2D *rad_mat);
Matrix2D* chol_lower(Matrix2D *A);
Matrix2D* matrix_multiply(Matrix2D *A, Matrix2D *B);
Matrix2D* matrix_transpose(Matrix2D *A);
Matrix2D* matrix_inverse(Matrix2D *A);
Vector* matrix_vector_multiply(Matrix2D *A, Vector *x);
Matrix2D* matrix3d2blkdiag(Matrix3D *T);
Matrix2D* logSO3(Matrix2D *R);
Matrix2D* expSO3(Vector *omega);
Matrix2D* invSO3(Matrix2D *R);

// Memory management functions
Matrix2D* create_matrix2d(int rows, int cols);
Matrix3D* create_matrix3d(int rows, int cols, int depth);
Vector* create_vector(int length);
State* create_state(void);
ErrorStruct* create_error_struct(void);
Measurements* create_measurements(void);
void free_matrix2d(Matrix2D *mat);
void free_matrix3d(Matrix3D *mat);
void free_vector(Vector *vec);
void free_state(State *state);
void free_error_struct(ErrorStruct *err);
void free_measurements(Measurements *meas);

#endif // PROJECT_SPECIFIC_H 