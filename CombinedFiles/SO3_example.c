#include "SO3.h"
#include <stdio.h>

// Helper function to print a 3x3 matrix
void print_matrix(const char* name, double matrix[3][3]) {
    printf("%s:\n", name);
    for (int i = 0; i < 3; i++) {
        printf("  [");
        for (int j = 0; j < 3; j++) {
            printf("%8.4f", matrix[i][j]);
            if (j < 2) printf(", ");
        }
        printf("]\n");
    }
    printf("\n");
}

// Helper function to print a 3D vector
void print_vector(const char* name, double vector[3]) {
    printf("%s: [%8.4f, %8.4f, %8.4f]\n\n", name, vector[0], vector[1], vector[2]);
}

int main() {
    printf("SO(3) Library Example\n");
    printf("=====================\n\n");
    
    // Example 1: Hat and Vec operators
    printf("Example 1: Hat and Vec operators\n");
    printf("--------------------------------\n");
    double w[3] = {0.1, 0.2, 0.3};
    double w_hat[3][3];
    double w_recovered[3];
    
    print_vector("Original vector w", w);
    
    HatSO3(w, w_hat);
    print_matrix("Hat(w) - skew-symmetric matrix", w_hat);
    
    VecSO3(w_hat, w_recovered);
    print_vector("Vec(Hat(w)) - recovered vector", w_recovered);
    
    // Example 2: Exponential and Logarithm maps
    printf("Example 2: Exponential and Logarithm maps\n");
    printf("-----------------------------------------\n");
    double R[3][3];
    double w_log[3];
    
    expSO3(w, R);
    print_matrix("R = exp(w) - rotation matrix", R);
    
    logSO3(R, w_log, NULL);
    print_vector("log(R) - recovered axis-angle", w_log);
    
    // Example 3: Matrix inverse
    printf("Example 3: Matrix inverse\n");
    printf("------------------------\n");
    double R_inv[3][3];
    double identity_check[3][3];
    
    invSO3(R, R_inv, NULL);
    print_matrix("R^(-1) - inverse matrix", R_inv);
    
    matrix_multiply(R, R_inv, identity_check);
    print_matrix("R * R^(-1) - should be identity", identity_check);
    
    // Example 4: Adjoint representation
    printf("Example 4: Adjoint representation\n");
    printf("---------------------------------\n");
    double adw[3][3];
    
    adjSO3(w, adw);
    print_matrix("adj(w) - adjoint representation", adw);
    
    // Example 5: Left Jacobian
    printf("Example 5: Left Jacobian\n");
    printf("------------------------\n");
    double Phi[3][3];
    
    PhiSO3(w, Phi);
    print_matrix("Phi(w) - Left Jacobian", Phi);
    
    printf("All examples completed successfully!\n");
    
    return 0;
} 