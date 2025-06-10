#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

// Define constants
#define MAX_STRING_LEN 256
#define MAX_CASES 50
#define MAX_COMPONENTS 3
#define MAX_TIME_POINTS 10000

// Data structures to replace MATLAB structs
typedef struct {
    char ls_c[16];      // Line style character
    char m_c[16];       // Marker character  
    int m_space;        // Marker spacing
    int m_offset;       // Marker offset
} LineSpecs;

typedef struct {
    char variable_name[MAX_STRING_LEN];
    char title_label[MAX_STRING_LEN];
    char y_unit[MAX_STRING_LEN];
    double (*transformation)(double);  // Function pointer for transformations
} FigurePlotOpts;

typedef struct {
    double p[3][MAX_TIME_POINTS];      // Position data [x,y,z][time]
    double R[3][MAX_TIME_POINTS];      // Rotation data [x,y,z][time]  
    double v[3][MAX_TIME_POINTS];      // Velocity data [x,y,z][time]
    double w[3][MAX_TIME_POINTS];      // Angular velocity [x,y,z][time]
    char label[MAX_STRING_LEN];        // Data label
    int num_points;                    // Number of time points
} NavigationData;

typedef struct {
    NavigationData* data;
    LineSpecs line_specs;
    char case_name[MAX_STRING_LEN];
} CaseData;

// Function prototypes
LineSpecs get_ls_specs(bool has_ls, const char* ls, bool has_m, const char* m, 
                      bool has_m_space, int m_space, bool has_m_offset, int m_offset);
double identity_transform(double x);
double rad2deg_transform(double x);
int plot_error_general(NavigationData* data, double* time_array, int num_time_points,
                      CaseData* cases_plot, int num_cases, FigurePlotOpts* opts);

/**
 * Get line specification structure with default values
 * Equivalent to MATLAB function get_ls_specs(opt)
 */
LineSpecs get_ls_specs(bool has_ls, const char* ls, bool has_m, const char* m, 
                      bool has_m_space, int m_space, bool has_m_offset, int m_offset) {
    LineSpecs res;
    
    // Set line style
    if (has_ls && ls != NULL) {
        strncpy(res.ls_c, ls, sizeof(res.ls_c) - 1);
        res.ls_c[sizeof(res.ls_c) - 1] = '\0';
    } else {
        strcpy(res.ls_c, "-");  // Default line style
    }
    
    // Set marker
    if (has_m && m != NULL) {
        strncpy(res.m_c, m, sizeof(res.m_c) - 1);
        res.m_c[sizeof(res.m_c) - 1] = '\0';
    } else {
        strcpy(res.m_c, "none");  // Default marker
    }
    
    // Set marker spacing
    if (has_m_space) {
        res.m_space = m_space;
    } else {
        res.m_space = 1;  // Default marker spacing
    }
    
    // Set marker offset
    if (has_m_offset) {
        res.m_offset = m_offset;
    } else {
        res.m_offset = 1;  // Default marker offset
    }
    
    return res;
}

/**
 * Identity transformation function
 */
double identity_transform(double x) {
    return x;
}

/**
 * Radians to degrees transformation
 */
double rad2deg_transform(double x) {
    return x * 180.0 / M_PI;
}

/**
 * General error plotting function
 * Equivalent to MATLAB function IN_error_general
 * Note: This function prepares data for plotting but doesn't actually plot
 * since C doesn't have built-in plotting capabilities
 */
int plot_error_general(NavigationData* data, double* time_array, int num_time_points,
                      CaseData* cases_plot, int num_cases, FigurePlotOpts* opts) {
    
    printf("=== %s ===\n", opts->title_label);
    printf("Time points: %d, Cases: %d\n", num_time_points, num_cases);
    
    // Process each component (x, y, z)
    const char* directions[] = {"x", "y", "z"};
    
    for (int component = 0; component < 3; component++) {
        printf("\n--- %s component ---\n", directions[component]);
        
        for (int case_idx = 0; case_idx < num_cases; case_idx++) {
            NavigationData* case_data = &data[case_idx];
            
            printf("Case: %s\n", case_data->label);
            
            // Get appropriate data based on variable name
            double* component_data = NULL;
            if (strcmp(opts->variable_name, "p") == 0) {
                component_data = case_data->p[component];
            } else if (strcmp(opts->variable_name, "R") == 0) {
                component_data = case_data->R[component];
            } else if (strcmp(opts->variable_name, "v") == 0) {
                component_data = case_data->v[component];
            } else if (strcmp(opts->variable_name, "w") == 0) {
                component_data = case_data->w[component];
            }
            
            if (component_data == NULL) {
                printf("Warning: Variable %s not found\n", opts->variable_name);
                continue;
            }
            
            // Apply transformation and print sample data
            printf("Sample transformed data: ");
            for (int i = 0; i < (num_time_points < 5 ? num_time_points : 5); i++) {
                double transformed_val = opts->transformation(component_data[i]);
                printf("%.3f ", transformed_val);
            }
            printf("... %s\n", opts->y_unit);
            
            // Here you would call external plotting library
            // For example, using gnuplot:
            /*
            FILE *gnuplot = popen("gnuplot -persistent", "w");
            if (gnuplot) {
                fprintf(gnuplot, "set title '%s %s component'\n", opts->title_label, directions[component]);
                fprintf(gnuplot, "set xlabel 'Time [s]'\n");
                fprintf(gnuplot, "set ylabel '%s %s'\n", directions[component], opts->y_unit);
                fprintf(gnuplot, "plot '-' with lines title '%s'\n", case_data->label);
                
                for (int i = 0; i < num_time_points; i++) {
                    fprintf(gnuplot, "%f %f\n", time_array[i], 
                           opts->transformation(component_data[i]));
                }
                fprintf(gnuplot, "e\n");
                pclose(gnuplot);
            }
            */
        }
    }
    
    return 0;  // Success
}

/**
 * Position error plotting function
 * Equivalent to MATLAB function IN_error_position
 */
int plot_error_position(NavigationData* data, double* time_array, int num_time_points,
                       CaseData* cases_plot, int num_cases) {
    
    FigurePlotOpts opts;
    strcpy(opts.variable_name, "p");
    opts.transformation = identity_transform;
    strcpy(opts.title_label, "Position Error");
    strcpy(opts.y_unit, "[m]");
    
    return plot_error_general(data, time_array, num_time_points, cases_plot, num_cases, &opts);
}

/**
 * Rotation error plotting function  
 * Equivalent to MATLAB function IN_error_rotation
 */
int plot_error_rotation(NavigationData* data, double* time_array, int num_time_points,
                       CaseData* cases_plot, int num_cases) {
    
    FigurePlotOpts opts;
    strcpy(opts.variable_name, "R");
    opts.transformation = rad2deg_transform;
    strcpy(opts.title_label, "Rotation Error");
    strcpy(opts.y_unit, "[deg]");
    
    return plot_error_general(data, time_array, num_time_points, cases_plot, num_cases, &opts);
}

/**
 * Calculate RMSE for all components
 * Equivalent to MATLAB function rmse_all_components_general
 */
int plot_rmse_all_components(NavigationData* data, double* time_array, int num_time_points,
                            CaseData* cases_plot, int num_cases, FigurePlotOpts* opts) {
    
    printf("=== %s (RMSE All Components) ===\n", opts->title_label);
    
    for (int case_idx = 0; case_idx < num_cases; case_idx++) {
        NavigationData* case_data = &data[case_idx];
        printf("Case: %s\n", case_data->label);
        
        // Calculate RMSE over all components for each time point
        for (int t = 0; t < num_time_points; t++) {
            double sum_squares = 0.0;
            
            // Get data based on variable name
            double* component_data[3];
            if (strcmp(opts->variable_name, "p") == 0) {
                for (int i = 0; i < 3; i++) {
                    component_data[i] = case_data->p[i];
                }
            } else if (strcmp(opts->variable_name, "R") == 0) {
                for (int i = 0; i < 3; i++) {
                    component_data[i] = case_data->R[i];
                }
            }
            
            // Calculate sum of squares for all components
            for (int comp = 0; comp < 3; comp++) {
                double transformed_val = opts->transformation(component_data[comp][t]);
                sum_squares += transformed_val * transformed_val;
            }
            
            // Calculate RMSE
            double rmse = sqrt(sum_squares / 3.0);
            
            // Print sample RMSE values (first few time points)
            if (t < 5) {
                printf("t=%.2f: RMSE=%.6f %s\n", time_array[t], rmse, opts->y_unit);
            }
        }
        printf("\n");
    }
    
    return 0;
}

/**
 * RMSE position plotting function
 * Equivalent to MATLAB function rmse_all_components_position
 */
int plot_rmse_position(NavigationData* data, double* time_array, int num_time_points,
                      CaseData* cases_plot, int num_cases) {
    
    FigurePlotOpts opts;
    strcpy(opts.variable_name, "p");
    opts.transformation = identity_transform;
    strcpy(opts.title_label, "RMSE Position");
    strcpy(opts.y_unit, "[m]");
    
    return plot_rmse_all_components(data, time_array, num_time_points, cases_plot, num_cases, &opts);
}

/**
 * RMSE rotation plotting function
 * Equivalent to MATLAB function rmse_all_components_rotation  
 */
int plot_rmse_rotation(NavigationData* data, double* time_array, int num_time_points,
                      CaseData* cases_plot, int num_cases) {
    
    FigurePlotOpts opts;
    strcpy(opts.variable_name, "R");
    opts.transformation = rad2deg_transform;
    strcpy(opts.title_label, "RMSE Rotation");
    strcpy(opts.y_unit, "[deg]");
    
    return plot_rmse_all_components(data, time_array, num_time_points, cases_plot, num_cases, &opts);
}

/**
 * Utility function to calculate norm over time
 * Equivalent to MATLAB's norm function applied over time dimension
 */
double calculate_norm_at_time(double** data, int num_components, int time_index) {
    double sum_squares = 0.0;
    for (int i = 0; i < num_components; i++) {
        sum_squares += data[i][time_index] * data[i][time_index];
    }
    return sqrt(sum_squares);
}

/**
 * Memory allocation helper for navigation data
 */
NavigationData* allocate_navigation_data(int num_cases) {
    NavigationData* data = (NavigationData*)malloc(num_cases * sizeof(NavigationData));
    if (!data) {
        fprintf(stderr, "Error: Failed to allocate memory for navigation data\n");
        return NULL;
    }
    
    // Initialize all data to zero
    for (int i = 0; i < num_cases; i++) {
        memset(&data[i], 0, sizeof(NavigationData));
        data[i].num_points = 0;
    }
    
    return data;
}

/**
 * Free navigation data memory
 */
void free_navigation_data(NavigationData* data) {
    if (data) {
        free(data);
    }
}

/**
 * Example usage and test function
 */
int main() {
    printf("C Navigation Data Plotting Functions\n");
    printf("=====================================\n");
    
    // Example: Create some test data
    const int num_cases = 2;
    const int num_time_points = 100;
    
    NavigationData* test_data = allocate_navigation_data(num_cases);
    if (!test_data) {
        return -1;
    }
    
    double time_array[100];
    
    // Initialize test data
    for (int case_idx = 0; case_idx < num_cases; case_idx++) {
        sprintf(test_data[case_idx].label, "Test Case %d", case_idx + 1);
        test_data[case_idx].num_points = num_time_points;
        
        for (int t = 0; t < num_time_points; t++) {
            time_array[t] = t * 0.1;  // 0.1 second intervals
            
            // Generate some test sinusoidal data
            for (int comp = 0; comp < 3; comp++) {
                double phase = (case_idx + 1) * (comp + 1) * 0.1;
                test_data[case_idx].p[comp][t] = sin(time_array[t] + phase) * (comp + 1);
                test_data[case_idx].R[comp][t] = cos(time_array[t] + phase) * 0.1;
            }
        }
    }
    
    CaseData cases[2];
    for (int i = 0; i < num_cases; i++) {
        sprintf(cases[i].case_name, "case_%d", i);
        cases[i].line_specs = get_ls_specs(false, NULL, false, NULL, false, 0, false, 0);
    }
    
    // Test the plotting functions
    printf("\nTesting position error plotting:\n");
    plot_error_position(test_data, time_array, num_time_points, cases, num_cases);
    
    printf("\nTesting rotation error plotting:\n");
    plot_error_rotation(test_data, time_array, num_time_points, cases, num_cases);
    
    printf("\nTesting RMSE position plotting:\n");
    plot_rmse_position(test_data, time_array, num_time_points, cases, num_cases);
    
    // Clean up
    free_navigation_data(test_data);
    
    printf("\nNote: This C version provides data processing and preparation.\n");
    printf("For actual plotting, integrate with libraries like:\n");
    printf("- gnuplot (pipe data to gnuplot commands)\n");
    printf("- Cairo/GTK+ (for GUI applications)\n");
    printf("- OpenGL (for 3D visualization)\n");
    printf("- Write CSV files for external plotting tools\n");
    
    return 0;
} 