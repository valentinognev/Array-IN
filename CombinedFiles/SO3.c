#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define M_PI 3.14159265358979323846

// Helper function to compute matrix norm (Frobenius norm)
double matrix_norm(double matrix[3][3]) {
    double sum = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            sum += matrix[i][j] * matrix[i][j];
        }
    }
    return sqrt(sum);
}

// Helper function to compute vector norm
double vector_norm(double v[3]) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// Helper function to compute matrix trace
double matrix_trace(double matrix[3][3]) {
    return matrix[0][0] + matrix[1][1] + matrix[2][2];
}

// Helper function for matrix multiplication
void matrix_multiply(double A[3][3], double B[3][3], double result[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i][j] = 0.0;
            for (int k = 0; k < 3; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Helper function to set identity matrix
void set_identity(double matrix[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            matrix[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

// Helper function to copy matrix
void copy_matrix(double src[3][3], double dest[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

void adjSO3(double w[3], double adw[3][3]) {
    HatSO3(w, adw);
}

void AdSO3(double R[3][3], double AdR[3][3]) {
    copy_matrix(R, AdR);
}

void dualHatSO3(double x[3], double x_dual_hat[3][3]) {
    // w_hat*x = x_dual_hat*w
    double x_hat[3][3];
    HatSO3(x, x_hat);
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            x_dual_hat[i][j] = -x_hat[i][j];
        }
    }
}

void errorSO3(double Rhat[][3][3], double Rtrue[3][3], double err[][3], int N, int Rtrue_is_array) {
    // UNTITLED Summary of this function goes here
    // Need to be same length
    
    // Calculate errors
    for (int n = 0; n < N; n++) {
        double Rhat_inv[3][3];
        double temp_product[3][3];
        double log_result[3];
        
        invSO3(Rhat[n], Rhat_inv, NULL);
        
        if (Rtrue_is_array) {
            matrix_multiply(Rhat_inv, ((double(*)[3][3])Rtrue)[n], temp_product);
        } else {
            matrix_multiply(Rhat_inv, Rtrue, temp_product);
        }
        
        logSO3(temp_product, log_result, NULL);
        
        for (int i = 0; i < 3; i++) {
            err[i][n] = log_result[i];
        }
    }
}

void expSO3(double w[3], double R[3][3]) {
    double normw = vector_norm(w);
    
    if (normw == 0) {
        set_identity(R);
        return;
    }
    
    double w_hat[3][3];
    HatSO3(w, w_hat);
    
    double w_hat_squared[3][3];
    matrix_multiply(w_hat, w_hat, w_hat_squared);
    
    // R = eye(3) + sin(normw)*w_hat/normw + (1-cos(normw))*w_hat*w_hat/(normw^2);
    set_identity(R);
    
    double sin_normw = sin(normw);
    double cos_normw = cos(normw);
    
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            R[i][j] += (sin_normw / normw) * w_hat[i][j] + 
                      ((1 - cos_normw) / (normw * normw)) * w_hat_squared[i][j];
        }
    }
}

void HatSO3(double w[3], double w_hat[3][3]) {
    w_hat[0][0] = 0;      w_hat[0][1] = -w[2];   w_hat[0][2] = w[1];
    w_hat[1][0] = w[2];   w_hat[1][1] = 0;       w_hat[1][2] = -w[0];
    w_hat[2][0] = -w[1];  w_hat[2][1] = w[0];    w_hat[2][2] = 0;
    
    // Alternative implementation (commented out in original):
    // w_hat = zeros(3,3);
    // w_hat(2,1) = w(3);
    // w_hat(3,1) = -w(2);
    // 
    // w_hat(1,2) = -w(3);
    // w_hat(3,2) = w(1);
    // 
    // w_hat(1,3) = w(2);
    // w_hat(2,3) = -w(1);
}

void invSO3(double R[3][3], double Rinv[3][3], int* errorFlag) {
    if (errorFlag) *errorFlag = 0;
    
    // Rinv = R';  (transpose)
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            Rinv[i][j] = R[j][i];
        }
    }
}

void logSO3(double R[3][3], double w[3], int* errorflag) {
    double phy = acos((matrix_trace(R) - 1) / 2);
    
    if (fabs(phy) > M_PI) {
        fprintf(stderr, "Error: angle supérieur à pi\n");
        if (errorflag) *errorflag = 1;
        return;
    }
    
    if (phy == 0) {
        w[0] = w[1] = w[2] = 0.0;
    } else if (fabs(phy) == M_PI) {
        double A[3][3];
        double eye[3][3];
        set_identity(eye);
        
        // A = (R-eye(3))/2;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                A[i][j] = (R[i][j] - eye[i][j]) / 2;
            }
        }
        
        double w1 = sqrt(-((A[1][1] + A[2][2] - A[0][0]) / 2));
        double w2 = sqrt(-((A[0][0] + A[2][2] - A[1][1]) / 2));
        double w3 = sqrt(-((A[0][0] + A[1][1] - A[2][2]) / 2));
        
        if (w1 != 0) {
            if (A[0][1] < 0) {
                w2 = -w2;
            }
            if (A[0][2] < 0) {
                w3 = -w3;
            }
        } else if (w2 != 0) {
            if (A[1][2] < 0) {
                w3 = -w3;
            }
        }
        
        w[0] = w1 * phy;
        w[1] = w2 * phy;
        w[2] = w3 * phy;
    } else {
        double w_hat[3][3];
        double R_transpose[3][3];
        
        // Transpose R
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                R_transpose[i][j] = R[j][i];
            }
        }
        
        // w_hat = (R-R.')/(2*sin(phy))*phy; on remultiplie par phy pour retrouver le vecteur avec sa norme originale
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                w_hat[i][j] = (R[i][j] - R_transpose[i][j]) / (2 * sin(phy)) * phy;
            }
        }
        
        VecSO3(w_hat, w);
    }
    
    if (errorflag) *errorflag = 0;
}

// Simple SVD implementation would be complex, using a simplified normalization
void normalizeSO3(double R[3][3], double Rnorm[3][3]) {
    // [u,s,v] = svd(R);
    // Rnorm = u*v';
    // 
    // For simplicity, we'll use Gram-Schmidt orthogonalization
    // This is a simplified approach - for production use, implement proper SVD
    
    double u1[3], u2[3], u3[3];
    
    // First column
    double norm1 = sqrt(R[0][0]*R[0][0] + R[1][0]*R[1][0] + R[2][0]*R[2][0]);
    u1[0] = R[0][0] / norm1;
    u1[1] = R[1][0] / norm1;
    u1[2] = R[2][0] / norm1;
    
    // Second column (orthogonalize)
    double dot = u1[0]*R[0][1] + u1[1]*R[1][1] + u1[2]*R[2][1];
    u2[0] = R[0][1] - dot * u1[0];
    u2[1] = R[1][1] - dot * u1[1];
    u2[2] = R[2][1] - dot * u1[2];
    
    double norm2 = sqrt(u2[0]*u2[0] + u2[1]*u2[1] + u2[2]*u2[2]);
    u2[0] /= norm2;
    u2[1] /= norm2;
    u2[2] /= norm2;
    
    // Third column (cross product)
    u3[0] = u1[1]*u2[2] - u1[2]*u2[1];
    u3[1] = u1[2]*u2[0] - u1[0]*u2[2];
    u3[2] = u1[0]*u2[1] - u1[1]*u2[0];
    
    Rnorm[0][0] = u1[0]; Rnorm[0][1] = u2[0]; Rnorm[0][2] = u3[0];
    Rnorm[1][0] = u1[1]; Rnorm[1][1] = u2[1]; Rnorm[1][2] = u3[1];
    Rnorm[2][0] = u1[2]; Rnorm[2][1] = u2[2]; Rnorm[2][2] = u3[2];
}

void PhiSO3(double w[3], double Phiw[3][3]) {
    // Left-Jacobian to SO(3)
    // sum_k 1/(k + 1)! ad(w)^k
    double normw = vector_norm(w);
    
    // if(normw > pi/2)
    //     error('formula not sure')
    // end
    
    if (normw > 0) {
        double adw[3][3];
        adjSO3(w, adw);
        
        double adw2[3][3], adw3[3][3], adw4[3][3];
        matrix_multiply(adw, adw, adw2);
        matrix_multiply(adw2, adw, adw3);
        matrix_multiply(adw3, adw, adw4);
        
        set_identity(Phiw);
        
        double coeff1 = (1.0 / (2 * normw * normw)) * (4 - normw * sin(normw) - 4 * cos(normw));
        double coeff2 = (1.0 / (2 * normw * normw * normw)) * (4 * normw - 5 * sin(normw) + normw * cos(normw));
        double coeff3 = (1.0 / (2 * normw * normw * normw * normw)) * (2 - normw * sin(normw) - 2 * cos(normw));
        double coeff4 = (1.0 / (2 * normw * normw * normw * normw * normw)) * (2 * normw - 3 * sin(normw) + normw * cos(normw));
        
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                Phiw[i][j] += coeff1 * adw[i][j] + coeff2 * adw2[i][j] + 
                              coeff3 * adw3[i][j] + coeff4 * adw4[i][j];
            }
        }
    } else {
        set_identity(Phiw);
    }
}

void VecSO3(double w_hat[3][3], double w[3]) {
    w[0] = w_hat[2][1];
    w[1] = w_hat[0][2];
    w[2] = w_hat[1][0];
} 