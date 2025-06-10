#include "project_specific.h"

#define M_PI 3.14159265358979323846
/**
 * CALCULATE_TRAJECTORY_ERROR Summary of this function goes here
 * Detailed explanation goes here
 */
ErrorStruct* calculate_trajectory_error(State *S_hat, State *S_true, int *section, int section_size) {
    ErrorStruct *err = create_error_struct();
    
    // If no section specified, use all indices
    int *actual_section;
    int actual_section_size;
    
    if (section == NULL) {
        actual_section_size = S_true->R->depth;
        actual_section = malloc(actual_section_size * sizeof(int));
        for (int i = 0; i < actual_section_size; i++) {
            actual_section[i] = i;
        }
    } else {
        actual_section = section;
        actual_section_size = section_size;
    }
    
    // Calculate rotation error
    err->R = errorSO3_vector(S_hat->R, S_true->R);
    err->R_deg = rad2deg_vector(err->R);
    
    // Calculate angular velocity error if available
    if (S_hat->w != NULL) {
        err->w = create_matrix2d(S_hat->w->rows, actual_section_size);
        for (int i = 0; i < S_hat->w->rows; i++) {
            for (int j = 0; j < actual_section_size; j++) {
                int idx = actual_section[j];
                err->w->data[i][j] = S_hat->w->data[i][idx] - S_true->w->data[i][idx];
            }
        }
        err->w_deg = rad2deg_matrix(err->w);
    }
    
    // Calculate position error
    err->p = create_matrix2d(S_hat->p->rows, actual_section_size);
    for (int i = 0; i < S_hat->p->rows; i++) {
        for (int j = 0; j < actual_section_size; j++) {
            int idx = actual_section[j];
            err->p->data[i][j] = S_hat->p->data[i][idx] - S_true->p->data[i][idx];
        }
    }
    
    // Calculate velocity error
    err->v = create_matrix2d(S_hat->v->rows, actual_section_size);
    for (int i = 0; i < S_hat->v->rows; i++) {
        for (int j = 0; j < actual_section_size; j++) {
            int idx = actual_section[j];
            err->v->data[i][j] = S_hat->v->data[i][idx] - S_true->v->data[i][idx];
        }
    }
    
    // Calculate accelerometer bias error if available
    if (S_hat->b_a != NULL) {
        int inds_bias = S_hat->b_a->rows;
        err->b_a = create_matrix2d(inds_bias, actual_section_size);
        for (int i = 0; i < inds_bias; i++) {
            for (int j = 0; j < actual_section_size; j++) {
                int idx = actual_section[j];
                err->b_a->data[i][j] = S_hat->b_a->data[i][idx] - S_true->b_a->data[i][idx];
            }
        }
    }
    
    // Calculate gyroscope bias error if available
    if (S_hat->b_g != NULL) {
        err->b_g = create_matrix2d(S_hat->b_g->rows, actual_section_size);
        for (int i = 0; i < S_hat->b_g->rows; i++) {
            for (int j = 0; j < actual_section_size; j++) {
                int idx = actual_section[j];
                err->b_g->data[i][j] = S_hat->b_g->data[i][idx] - S_true->b_g->data[i][idx];
            }
        }
        err->b_g_deg = rad2deg_matrix(err->b_g);
    }
    
    if (section == NULL) {
        free(actual_section);
    }
    
    return err;
}

/**
 * COMPENSATE_COVARIANCE Summary of this function goes here
 * Detailed explanation goes here
 */
Matrix2D* compensate_covariance(Matrix2D *Q_y, Matrix2D *T) {
    Matrix2D *L = chol_lower(Q_y);
    Matrix2D *T_inv = matrix_inverse(T);
    Matrix2D *q = matrix_multiply(T_inv, L);
    Matrix2D *q_t = matrix_transpose(q);
    Matrix2D *Q_u = matrix_multiply(q, q_t);
    
    free_matrix2d(L);
    free_matrix2d(T_inv);
    free_matrix2d(q);
    free_matrix2d(q_t);
    
    return Q_u;
}

/**
 * COMPENSATE_MEASUREMENTS Summary of this function goes here
 * Detailed explanation goes here
 */
Vector* compensate_measurements(Vector *y, Matrix2D *T, Vector *b) {
    Matrix2D *T_diag;
    
    // Check if T is 3D (converted to block diagonal) or already 2D
    // For simplicity, assuming T is already in the correct format
    T_diag = T;
    
    // Create vector from b (reshape to column vector)
    Vector *b_col = create_vector(b->length);
    memcpy(b_col->data, b->data, b->length * sizeof(double));
    
    // Calculate y - b
    Vector *y_minus_b = create_vector(y->length);
    for (int i = 0; i < y->length; i++) {
        y_minus_b->data[i] = y->data[i] - b_col->data[i];
    }
    
    // Solve T_diag \ (y - b)
    Matrix2D *T_inv = matrix_inverse(T_diag);
    Vector *u = matrix_vector_multiply(T_inv, y_minus_b);
    
    free_vector(b_col);
    free_vector(y_minus_b);
    free_matrix2d(T_inv);
    
    return u;
}

// Memory management functions
Matrix2D* create_matrix2d(int rows, int cols) {
    Matrix2D *mat = malloc(sizeof(Matrix2D));
    mat->rows = rows;
    mat->cols = cols;
    mat->data = malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        mat->data[i] = calloc(cols, sizeof(double));
    }
    return mat;
}

Matrix3D* create_matrix3d(int rows, int cols, int depth) {
    Matrix3D *mat = malloc(sizeof(Matrix3D));
    mat->rows = rows;
    mat->cols = cols;
    mat->depth = depth;
    mat->data = malloc(rows * sizeof(double**));
    for (int i = 0; i < rows; i++) {
        mat->data[i] = malloc(cols * sizeof(double*));
        for (int j = 0; j < cols; j++) {
            mat->data[i][j] = calloc(depth, sizeof(double));
        }
    }
    return mat;
}

Vector* create_vector(int length) {
    Vector *vec = malloc(sizeof(Vector));
    vec->length = length;
    vec->data = calloc(length, sizeof(double));
    return vec;
}

State* create_state(void) {
    State *state = malloc(sizeof(State));
    memset(state, 0, sizeof(State));
    return state;
}

ErrorStruct* create_error_struct(void) {
    ErrorStruct *err = malloc(sizeof(ErrorStruct));
    memset(err, 0, sizeof(ErrorStruct));
    return err;
}

Measurements* create_measurements(void) {
    Measurements *meas = malloc(sizeof(Measurements));
    memset(meas, 0, sizeof(Measurements));
    return meas;
}

void free_matrix2d(Matrix2D *mat) {
    if (mat) {
        for (int i = 0; i < mat->rows; i++) {
            free(mat->data[i]);
        }
        free(mat->data);
        free(mat);
    }
}

void free_matrix3d(Matrix3D *mat) {
    if (mat) {
        for (int i = 0; i < mat->rows; i++) {
            for (int j = 0; j < mat->cols; j++) {
                free(mat->data[i][j]);
            }
            free(mat->data[i]);
        }
        free(mat->data);
        free(mat);
    }
}

void free_vector(Vector *vec) {
    if (vec) {
        free(vec->data);
        free(vec);
    }
}

void free_state(State *state) {
    if (state) {
        free_matrix3d(state->R);
        free_matrix2d(state->w);
        free_matrix2d(state->p);
        free_matrix2d(state->v);
        free_matrix2d(state->b_a);
        free_matrix2d(state->b_g);
        free_matrix2d(state->omega_dot);
        free_matrix2d(state->v_dot);
        free_matrix2d(state->s);
        free_matrix2d(state->b_s);
        free_matrix3d(state->T_a);
        free_matrix2d(state->b_omega_dot);
        free_matrix2d(state->p_rig_time);
        free_matrix3d(state->R_rig_time);
        free(state);
    }
}

void free_error_struct(ErrorStruct *err) {
    if (err) {
        free_vector(err->R);
        free_vector(err->R_deg);
        free_matrix2d(err->w);
        free_matrix2d(err->w_deg);
        free_matrix2d(err->p);
        free_matrix2d(err->v);
        free_matrix2d(err->b_a);
        free_matrix2d(err->b_g);
        free_matrix2d(err->b_g_deg);
        free_matrix2d(err->omega_dot);
        free_matrix2d(err->v_dot);
        free_matrix2d(err->s);
        free_matrix2d(err->b_s);
        free_matrix3d(err->T_a);
        free_matrix2d(err->b_omega_dot);
        free(err);
    }
}

void free_measurements(Measurements *meas) {
    if (meas) {
        free_matrix2d(meas->y);
        free_matrix2d(meas->Q);
        free_matrix2d(meas->Q_inv);
        free(meas);
    }
}

// Utility function implementations
double rad2deg(double rad) {
    return rad * 180.0 / M_PI;
}

Vector* rad2deg_vector(Vector *rad_vec) {
    Vector *deg_vec = create_vector(rad_vec->length);
    for (int i = 0; i < rad_vec->length; i++) {
        deg_vec->data[i] = rad2deg(rad_vec->data[i]);
    }
    return deg_vec;
}

Matrix2D* rad2deg_matrix(Matrix2D *rad_mat) {
    Matrix2D *deg_mat = create_matrix2d(rad_mat->rows, rad_mat->cols);
    for (int i = 0; i < rad_mat->rows; i++) {
        for (int j = 0; j < rad_mat->cols; j++) {
            deg_mat->data[i][j] = rad2deg(rad_mat->data[i][j]);
        }
    }
    return deg_mat;
}

// Placeholder implementations for specialized matrix functions
double errorSO3(Matrix2D *R1, Matrix2D *R2) {
    // Placeholder - implement proper SO(3) error calculation
    return 0.0;
}

Vector* errorSO3_vector(Matrix3D *R_hat, Matrix3D *R_true) {
    // Placeholder - implement proper SO(3) error calculation for vector of rotations
    Vector *err = create_vector(R_hat->depth);
    return err;
}

Matrix2D* chol_lower(Matrix2D *A) {
    // Placeholder - implement Cholesky decomposition
    return create_matrix2d(A->rows, A->cols);
}

Matrix2D* matrix_multiply(Matrix2D *A, Matrix2D *B) {
    // Placeholder - implement matrix multiplication
    Matrix2D *C = create_matrix2d(A->rows, B->cols);
    return C;
}

Matrix2D* matrix_transpose(Matrix2D *A) {
    Matrix2D *At = create_matrix2d(A->cols, A->rows);
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            At->data[j][i] = A->data[i][j];
        }
    }
    return At;
}

Matrix2D* matrix_inverse(Matrix2D *A) {
    // Placeholder - implement matrix inversion
    return create_matrix2d(A->rows, A->cols);
}

Vector* matrix_vector_multiply(Matrix2D *A, Vector *x) {
    Vector *y = create_vector(A->rows);
    for (int i = 0; i < A->rows; i++) {
        y->data[i] = 0.0;
        for (int j = 0; j < A->cols; j++) {
            y->data[i] += A->data[i][j] * x->data[j];
        }
    }
    return y;
}

/**
 * INTERPOLATE_POS_AND_ROTATION Interpolate IMU pos and rotation estimates 
 * to rig time
 */
void interpolate_pos_and_rotation(State *S, Vector *imu_time, Vector *rig_time) {
    S->p_rig_time = create_matrix2d(3, rig_time->length);
    
    // Interpolate position for each component
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < rig_time->length; j++) {
            // Simple linear interpolation (PCHIP would be more complex to implement)
            double t = rig_time->data[j];
            int idx = 0;
            
            // Find interpolation indices
            while (idx < imu_time->length - 1 && imu_time->data[idx + 1] < t) {
                idx++;
            }
            
            if (idx < imu_time->length - 1) {
                double t1 = imu_time->data[idx];
                double t2 = imu_time->data[idx + 1];
                double p1 = S->p->data[i][idx];
                double p2 = S->p->data[i][idx + 1];
                double alpha = (t - t1) / (t2 - t1);
                S->p_rig_time->data[i][j] = p1 + alpha * (p2 - p1);
            } else {
                S->p_rig_time->data[i][j] = S->p->data[i][idx];
            }
        }
    }
    
    // Find the fractional indices using linear interpolation 
    Vector *inds_time = create_vector(rig_time->length);
    for (int j = 0; j < rig_time->length; j++) {
        double t = rig_time->data[j];
        int idx = 0;
        
        while (idx < imu_time->length - 1 && imu_time->data[idx + 1] < t) {
            idx++;
        }
        
        if (idx < imu_time->length - 1) {
            double t1 = imu_time->data[idx];
            double t2 = imu_time->data[idx + 1];
            double alpha = (t - t1) / (t2 - t1);
            inds_time->data[j] = idx + alpha;
        } else {
            inds_time->data[j] = idx;
        }
    }
    
    S->R_rig_time = create_matrix3d(3, 3, rig_time->length);
    
    for (int n = 0; n < rig_time->length; n++) {
        double frac = inds_time->data[n] - floor(inds_time->data[n]);
        
        if (frac > 1) {
            printf("Warning: fraction larger than 1\n");
        } else if (frac == 0) {
            // Same point in time
            int idx = (int)round(inds_time->data[n]);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    S->R_rig_time->data[i][j][n] = S->R->data[i][j][idx];
                }
            }
        } else {
            // Calculate the rotation vector and scale it
            int left_ind = (int)floor(inds_time->data[n]);
            int right_ind = left_ind + 1;
            
            // Extract rotation matrices
            Matrix2D *R_left = create_matrix2d(3, 3);
            Matrix2D *R_right = create_matrix2d(3, 3);
            
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    R_left->data[i][j] = S->R->data[i][j][left_ind];
                    R_right->data[i][j] = S->R->data[i][j][right_ind];
                }
            }
            
            Matrix2D *R_left_inv = invSO3(R_left);
            Matrix2D *R_diff = matrix_multiply(R_left_inv, R_right);
            Matrix2D *theta_mat = logSO3(R_diff);
            
            // Convert matrix to vector (simplified)
            Vector *theta = create_vector(3);
            theta->data[0] = theta_mat->data[2][1];
            theta->data[1] = theta_mat->data[0][2];
            theta->data[2] = theta_mat->data[1][0];
            
            // Scale by fraction
            Vector *scaled_theta = create_vector(3);
            for (int i = 0; i < 3; i++) {
                scaled_theta->data[i] = frac * theta->data[i];
            }
            
            Matrix2D *exp_scaled = expSO3(scaled_theta);
            Matrix2D *R_result = matrix_multiply(R_left, exp_scaled);
            
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    S->R_rig_time->data[i][j][n] = R_result->data[i][j];
                }
            }
            
            free_matrix2d(R_left);
            free_matrix2d(R_right);
            free_matrix2d(R_left_inv);
            free_matrix2d(R_diff);
            free_matrix2d(theta_mat);
            free_vector(theta);
            free_vector(scaled_theta);
            free_matrix2d(exp_scaled);
            free_matrix2d(R_result);
        }
    }
    
    free_vector(inds_time);
}

/**
 * LSQ_TRIAD Weighted Mean of triad
 * y = Hu + e , e ~ N(0,Q)
 * u = (H'*Q^{-1}*H)^{-1}(H'*Q^{-1}*y)
 * Where u is triad
 */
void lsq_triad(Vector *y, Matrix2D *Q, Vector **u_out, Matrix2D **Qu_out) {
    assert(y->length % 3 == 0);
    assert(Q->rows == Q->cols);
    
    int N = y->length / 3;
    Matrix2D *H = create_matrix2d(y->length, 3);
    
    // H = repmat(eye(3), N, 1)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                H->data[i*3 + j][k] = (j == k) ? 1.0 : 0.0;
            }
        }
    }
    
    Matrix2D *L = chol_lower(Q);
    Matrix2D *H_t = matrix_transpose(H);
    Matrix2D *L_inv = matrix_inverse(L);
    Matrix2D *t1 = matrix_multiply(H_t, L_inv);
    
    Vector *L_inv_y = matrix_vector_multiply(L_inv, y);
    Matrix2D *t1_t1t = matrix_multiply(t1, matrix_transpose(t1));
    Matrix2D *t1_t1t_inv = matrix_inverse(t1_t1t);
    Vector *t1_t2 = matrix_vector_multiply(t1, L_inv_y);
    
    *u_out = matrix_vector_multiply(t1_t1t_inv, t1_t2);
    *Qu_out = t1_t1t_inv; // Don't free this as it's returned
    
    free_matrix2d(H);
    free_matrix2d(L);
    free_matrix2d(H_t);
    free_matrix2d(L_inv);
    free_matrix2d(t1);
    free_vector(L_inv_y);
    free_matrix2d(t1_t1t);
    free_vector(t1_t2);
}

/**
 * LSQ_TRIAD Summary of this function goes here
 * y = Hu + e , e ~ N(0,Q)
 * u = (H'*Q^{-1}*H)^{-1}(H'*Q^{-1}*y)
 * Where u is triad
 */
void lsq_triad_naive(Vector *y, Matrix2D *Q, Vector **u_out, Matrix2D **Qu_out) {
    assert(y->length % 3 == 0);
    assert(Q->rows == Q->cols);
    
    int N = y->length / 3;
    Matrix2D *H = create_matrix2d(y->length, 3);
    
    // H = kron(ones(N,1),eye(3))
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                H->data[i*3 + j][k] = (j == k) ? 1.0 : 0.0;
            }
        }
    }
    
    Matrix2D *H_t = matrix_transpose(H);
    Matrix2D *Q_inv = matrix_inverse(Q);
    Matrix2D *HT_Q_inv = matrix_multiply(H_t, Q_inv);
    Matrix2D *HT_Q_inv_H = matrix_multiply(HT_Q_inv, H);
    Matrix2D *HT_Q_inv_H_inv = matrix_inverse(HT_Q_inv_H);
    Vector *HT_Q_inv_y = matrix_vector_multiply(HT_Q_inv, y);
    
    *u_out = matrix_vector_multiply(HT_Q_inv_H_inv, HT_Q_inv_y);
    *Qu_out = HT_Q_inv_H_inv; // Don't free this as it's returned
    
    free_matrix2d(H);
    free_matrix2d(H_t);
    free_matrix2d(Q_inv);
    free_matrix2d(HT_Q_inv);
    free_matrix2d(HT_Q_inv_H);
    free_vector(HT_Q_inv_y);
}

/**
 * ROTATE_MEASUREMENTS Summary of this function goes here
 * Detailed explanation goes here
 */
Measurements* rotate_measurements(Measurements *S_in, Matrix2D *R) {
    Measurements *S_out = create_measurements();
    
    S_out->y = matrix_multiply(R, S_in->y);
    Matrix2D *R_t = matrix_transpose(R);
    Matrix2D *temp = matrix_multiply(R, S_in->Q);
    S_out->Q = matrix_multiply(temp, R_t);
    S_out->Q_inv = matrix_inverse(S_out->Q);
    
    free_matrix2d(R_t);
    free_matrix2d(temp);
    
    return S_out;
}

/**
 * RUN_FILTER Run filter and calculate error 
 * Detailed explanation goes here
 */
CompleteResults* run_filter(Measurements *sensorData, State *initData, void *my_settings, State *S_ref, void *myFilter) {
    CompleteResults *res = malloc(sizeof(CompleteResults));
    res->results = malloc(sizeof(FilterResults));
    
    // Note: This is a placeholder since myFilter is a function pointer
    // In actual implementation, you would call: myFilter(sensorData, initData, my_settings)
    // For now, we'll assume the filter results are populated elsewhere
    res->results->filt = create_state();
    res->results->pred = create_state();
    
    ErrorStruct *err = create_error_struct();
    
    // Calculate rotation error if available
    if (S_ref->R != NULL) {
        // Try to calculate rotation error
        // In case of error, print warning (equivalent to MATLAB's try/catch)
        double error_val = errorSO3(res->results->filt->R->data[0], S_ref->R->data[0]);
        if (isnan(error_val) || isinf(error_val)) {
            printf("Warning: Angle Error is too high.\n");
        } else {
            err->R = create_vector(1);
            err->R->data[0] = error_val;
        }
    }
    
    // Calculate velocity error if available
    if (S_ref->v != NULL) {
        err->v = create_matrix2d(res->results->filt->v->rows, res->results->filt->v->cols);
        for (int i = 0; i < res->results->filt->v->rows; i++) {
            for (int j = 0; j < res->results->filt->v->cols; j++) {
                err->v->data[i][j] = res->results->filt->v->data[i][j] - S_ref->v->data[i][j];
            }
        }
    }
    
    // Calculate position error if available
    if (S_ref->p != NULL) {
        err->p = create_matrix2d(res->results->filt->p->rows, res->results->filt->p->cols);
        for (int i = 0; i < res->results->filt->p->rows; i++) {
            for (int j = 0; j < res->results->filt->p->cols; j++) {
                err->p->data[i][j] = res->results->filt->p->data[i][j] - S_ref->p->data[i][j];
            }
        }
    }
    
    // Calculate angular velocity error if available
    if (S_ref->w != NULL && res->results->filt->w != NULL) {
        err->w = create_matrix2d(res->results->filt->w->rows, res->results->filt->w->cols);
        for (int i = 0; i < res->results->filt->w->rows; i++) {
            for (int j = 0; j < res->results->filt->w->cols; j++) {
                err->w->data[i][j] = res->results->filt->w->data[i][j] - S_ref->w->data[i][j];
            }
        }
    }
    
    // Calculate angular acceleration error if available
    if (S_ref->omega_dot != NULL) {
        err->omega_dot = create_matrix2d(res->results->pred->omega_dot->rows, res->results->pred->omega_dot->cols);
        for (int i = 0; i < res->results->pred->omega_dot->rows; i++) {
            for (int j = 0; j < res->results->pred->omega_dot->cols; j++) {
                err->omega_dot->data[i][j] = res->results->pred->omega_dot->data[i][j] - S_ref->omega_dot->data[i][j];
            }
        }
    }
    
    // Calculate linear acceleration error if available
    if (S_ref->v_dot != NULL) {
        err->v_dot = create_matrix2d(res->results->pred->v_dot->rows, res->results->pred->v_dot->cols);
        for (int i = 0; i < res->results->pred->v_dot->rows; i++) {
            for (int j = 0; j < res->results->pred->v_dot->cols; j++) {
                err->v_dot->data[i][j] = res->results->pred->v_dot->data[i][j] - S_ref->v_dot->data[i][j];
            }
        }
    }
    
    // Store single error structure (not nested like in run_filter_w_error)
    res->err.filt = err;
    res->err.pred = NULL;
    
    return res;
}

/**
 * RUN_FILTER Run filter and calculate error 
 * Detailed explanation goes here
 */
CompleteResults* run_filter_w_error(Measurements *sensorData, State *initData, void *my_settings, State *S_ref, void *myFilter) {
    CompleteResults *res = malloc(sizeof(CompleteResults));
    res->results = malloc(sizeof(FilterResults));
    
    // Note: This is a placeholder since myFilter is a function pointer
    // In actual implementation, you would call: myFilter(sensorData, initData, my_settings)
    res->results->filt = create_state();
    res->results->pred = create_state();
    
    // Calculate errors using compute_error function
    res->err.filt = compute_error(res->results->filt, S_ref);
    res->err.pred = compute_error(res->results->pred, S_ref);
    
    return res;
}

Vector* errorSO3_vector(Matrix3D *R_hat, Matrix3D *R_true) {
    // Placeholder - implement proper SO(3) error calculation for vector of rotations
    Vector *err = create_vector(R_hat->depth);
    return err;
}

Matrix2D* chol_lower(Matrix2D *A) {
    // Placeholder - implement Cholesky decomposition
    return create_matrix2d(A->rows, A->cols);
}

Matrix2D* matrix_multiply(Matrix2D *A, Matrix2D *B) {
    // Placeholder - implement matrix multiplication
    Matrix2D *C = create_matrix2d(A->rows, B->cols);
    return C;
}

Matrix2D* matrix_transpose(Matrix2D *A) {
    Matrix2D *At = create_matrix2d(A->cols, A->rows);
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            At->data[j][i] = A->data[i][j];
        }
    }
    return At;
}

Matrix2D* matrix_inverse(Matrix2D *A) {
    // Placeholder - implement matrix inversion
    return create_matrix2d(A->rows, A->cols);
}

Vector* matrix_vector_multiply(Matrix2D *A, Vector *x) {
    Vector *y = create_vector(A->rows);
    for (int i = 0; i < A->rows; i++) {
        y->data[i] = 0.0;
        for (int j = 0; j < A->cols; j++) {
            y->data[i] += A->data[i][j] * x->data[j];
        }
    }
    return y;
}