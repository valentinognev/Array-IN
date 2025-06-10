#ifndef PROJECT_SPECIFIC_PLOT_FUNCTIONS_H
#define PROJECT_SPECIFIC_PLOT_FUNCTIONS_H

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

/**
 * Line specification structure
 * Equivalent to MATLAB line specification options
 */
typedef struct {
    char ls_c[16];      // Line style character ("-", "--", ":", etc.)
    char m_c[16];       // Marker character ("o", "s", "^", "none", etc.)
    int m_space;        // Marker spacing (every Nth point)
    int m_offset;       // Marker offset (starting point)
} LineSpecs;

/**
 * Figure plotting options structure
 * Contains plotting configuration parameters
 */
typedef struct {
    char variable_name[MAX_STRING_LEN];     // Variable name to plot ("p", "R", "v", "w")
    char title_label[MAX_STRING_LEN];       // Plot title
    char y_unit[MAX_STRING_LEN];           // Y-axis unit label
    double (*transformation)(double);       // Function pointer for data transformation
} FigurePlotOpts;

/**
 * Navigation data structure
 * Contains time series data for inertial navigation analysis
 */
typedef struct {
    double p[3][MAX_TIME_POINTS];      // Position data [x,y,z][time]
    double R[3][MAX_TIME_POINTS];      // Rotation data [x,y,z][time]  
    double v[3][MAX_TIME_POINTS];      // Velocity data [x,y,z][time]
    double w[3][MAX_TIME_POINTS];      // Angular velocity [x,y,z][time]
    double b_a[3][MAX_TIME_POINTS];    // Accelerometer bias [x,y,z][time]
    double b_g[3][MAX_TIME_POINTS];    // Gyroscope bias [x,y,z][time]
    double s[3][MAX_TIME_POINTS];      // Specific force [x,y,z][time]
    char label[MAX_STRING_LEN];        // Data label for plotting
    int num_points;                    // Number of time points
} NavigationData;

/**
 * Case data structure
 * Combines navigation data with plotting specifications
 */
typedef struct {
    NavigationData* data;              // Pointer to navigation data
    LineSpecs line_specs;              // Line specification for plotting
    char case_name[MAX_STRING_LEN];    // Case identifier
} CaseData;

// Function prototypes

/**
 * Get line specification structure with default values
 * @param has_ls Whether line style is specified
 * @param ls Line style string
 * @param has_m Whether marker is specified
 * @param m Marker string
 * @param has_m_space Whether marker spacing is specified
 * @param m_space Marker spacing value
 * @param has_m_offset Whether marker offset is specified
 * @param m_offset Marker offset value
 * @return LineSpecs structure with specified or default values
 */
LineSpecs get_ls_specs(bool has_ls, const char* ls, bool has_m, const char* m, 
                      bool has_m_space, int m_space, bool has_m_offset, int m_offset);

/**
 * Identity transformation function
 * @param x Input value
 * @return Same value (no transformation)
 */
double identity_transform(double x);

/**
 * Radians to degrees transformation
 * @param x Input value in radians
 * @return Value converted to degrees
 */
double rad2deg_transform(double x);

/**
 * General error plotting function
 * Equivalent to MATLAB function IN_error_general
 * @param data Array of navigation data
 * @param time_array Time vector
 * @param num_time_points Number of time points
 * @param cases_plot Array of case data
 * @param num_cases Number of cases
 * @param opts Plotting options
 * @return 0 on success, -1 on error
 */
int plot_error_general(NavigationData* data, double* time_array, int num_time_points,
                      CaseData* cases_plot, int num_cases, FigurePlotOpts* opts);

/**
 * Position error plotting function
 * Equivalent to MATLAB function IN_error_position
 * @param data Array of navigation data
 * @param time_array Time vector
 * @param num_time_points Number of time points
 * @param cases_plot Array of case data
 * @param num_cases Number of cases
 * @return 0 on success, -1 on error
 */
int plot_error_position(NavigationData* data, double* time_array, int num_time_points,
                       CaseData* cases_plot, int num_cases);

/**
 * Rotation error plotting function  
 * Equivalent to MATLAB function IN_error_rotation
 * @param data Array of navigation data
 * @param time_array Time vector
 * @param num_time_points Number of time points
 * @param cases_plot Array of case data
 * @param num_cases Number of cases
 * @return 0 on success, -1 on error
 */
int plot_error_rotation(NavigationData* data, double* time_array, int num_time_points,
                       CaseData* cases_plot, int num_cases);

/**
 * Calculate RMSE for all components
 * Equivalent to MATLAB function rmse_all_components_general
 * @param data Array of navigation data
 * @param time_array Time vector
 * @param num_time_points Number of time points
 * @param cases_plot Array of case data
 * @param num_cases Number of cases
 * @param opts Plotting options
 * @return 0 on success, -1 on error
 */
int plot_rmse_all_components(NavigationData* data, double* time_array, int num_time_points,
                            CaseData* cases_plot, int num_cases, FigurePlotOpts* opts);

/**
 * RMSE position plotting function
 * Equivalent to MATLAB function rmse_all_components_position
 * @param data Array of navigation data
 * @param time_array Time vector
 * @param num_time_points Number of time points
 * @param cases_plot Array of case data
 * @param num_cases Number of cases
 * @return 0 on success, -1 on error
 */
int plot_rmse_position(NavigationData* data, double* time_array, int num_time_points,
                      CaseData* cases_plot, int num_cases);

/**
 * RMSE rotation plotting function
 * Equivalent to MATLAB function rmse_all_components_rotation  
 * @param data Array of navigation data
 * @param time_array Time vector
 * @param num_time_points Number of time points
 * @param cases_plot Array of case data
 * @param num_cases Number of cases
 * @return 0 on success, -1 on error
 */
int plot_rmse_rotation(NavigationData* data, double* time_array, int num_time_points,
                      CaseData* cases_plot, int num_cases);

/**
 * Utility function to calculate norm over time
 * Equivalent to MATLAB's norm function applied over time dimension
 * @param data 2D array of component data
 * @param num_components Number of components (typically 3)
 * @param time_index Time index to calculate norm at
 * @return Euclidean norm of components at specified time
 */
double calculate_norm_at_time(double** data, int num_components, int time_index);

/**
 * Memory allocation helper for navigation data
 * @param num_cases Number of cases to allocate
 * @return Pointer to allocated NavigationData array, NULL on failure
 */
NavigationData* allocate_navigation_data(int num_cases);

/**
 * Free navigation data memory
 * @param data Pointer to navigation data to free
 */
void free_navigation_data(NavigationData* data);

// Additional plotting functions that could be implemented
// (These would require more complex data structures and external plotting libraries)

/**
 * Plot angular velocity data
 * Equivalent to MATLAB function plot_angular_velocity
 */
int plot_angular_velocity(NavigationData* data, double* time_array, int num_time_points,
                         CaseData* cases_plot, int num_cases, const char* plot_type);

/**
 * Plot navigation position data
 * Equivalent to MATLAB function plot_navigation_position
 */
int plot_navigation_position(NavigationData* data, double* time_array, int num_time_points,
                            CaseData* cases_plot, int num_cases, const char* plot_type);

/**
 * Plot navigation velocity data
 * Equivalent to MATLAB function plot_navigation_velocity
 */
int plot_navigation_velocity(NavigationData* data, double* time_array, int num_time_points,
                            CaseData* cases_plot, int num_cases, const char* plot_type);

/**
 * Plot accelerometer bias data
 * Equivalent to MATLAB function plot_bias_accelerometer
 */
int plot_bias_accelerometer(NavigationData* data, double* time_array, int num_time_points,
                           CaseData* cases_plot, int num_cases, const char* plot_type);

/**
 * Plot gyroscope bias data
 * Equivalent to MATLAB function plot_bias_gyroscopes
 */
int plot_bias_gyroscopes(NavigationData* data, double* time_array, int num_time_points,
                        CaseData* cases_plot, int num_cases, const char* plot_type);

/**
 * Plot specific force data
 * Equivalent to MATLAB function plot_specific_force
 */
int plot_specific_force(NavigationData* data, double* time_array, int num_time_points,
                       CaseData* cases_plot, int num_cases, const char* plot_type);

// Utility macros for common operations
#define DEG_TO_RAD(x) ((x) * M_PI / 180.0)
#define RAD_TO_DEG(x) ((x) * 180.0 / M_PI)
#define SQUARE(x) ((x) * (x))
#define MAX_VAL(a, b) ((a) > (b) ? (a) : (b))
#define MIN_VAL(a, b) ((a) < (b) ? (a) : (b))

// Error codes
#define PLOT_SUCCESS 0
#define PLOT_ERROR_NULL_POINTER -1
#define PLOT_ERROR_INVALID_PARAMS -2
#define PLOT_ERROR_MEMORY_ALLOCATION -3
#define PLOT_ERROR_FILE_IO -4

#endif // PROJECT_SPECIFIC_PLOT_FUNCTIONS_H 