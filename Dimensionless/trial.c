#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "trial.h"
#include <math.h>
#include <time.h>

int main() {
    int M; // Total number of cycles (inactive + active)
    int active_cycle; // Number of active cycles 
    int inactive_cycle; // Number of inactive cycles
    char particle_name[100]; // Name of the particle
    double scaling_factor; // Scaling factor 

    neutron *parent_fission_bank = NULL; // Pointer to the parent fission bank (stack)
    neutron *child_fission_bank = NULL; // Pointer to the child fission bank (stack)
    neutron **parent_fission_bank_double_pointer = &parent_fission_bank; // Double pointer to the parent fission bank 
    neutron **child_fission_bank_double_pointer = &child_fission_bank; // Double pointer to the child fission bank
    
    double k1; // Theoretical eigenvalue for the first harmonic mode 
    double k2; // Theoretical eigenvalue for the second harmonic mode
    
    const char *input = "input.txt"; // File for reading the inputs
    const char *filename1 = "flux_data.txt"; // File for storing neutron flux data
    const char *filename2 = "text_data.txt"; // File for storing data of each cycle
    const char *filename3 = "position_data.txt"; // File for storing neutron position data
    const char *filename4 = "reweight_bank.txt"; // File for storing reweight data of the child fission bank 
    const char *filename5 = "k_data.txt"; // File for storing k data of each cycle
    const char *filename6 = "E_data.txt"; // File for storing E data of each cycle
    const char *filename7 = "fission_matrix.txt"; // File for storing the final fission matrix data

    // Register cleanup function
    atexit(cleanup);

    // Open text files
    FILE *input_file = fopen(input, "r");
    FILE *file1 = fopen(filename1, "w");
    FILE *file2 = fopen(filename2, "w");
    FILE *file3 = fopen(filename3, "w");
    FILE *file4 = fopen(filename4, "w");
    FILE *file5 = fopen(filename5, "w");
    FILE *file6 = fopen(filename6, "w");
    FILE *file7 = fopen(filename7, "w");
    if(input_file == NULL || file1 == NULL || file2 == NULL || file3 == NULL || file4 == NULL || file5 == NULL || file6 == NULL || file7 == NULL) {
        printf("Unable to open files.\n");
        exit(EXIT_FAILURE); // Terminate program immediatedly
    }
    register_pointer(file1, TYPE_FILE);
    register_pointer(file2, TYPE_FILE);
    register_pointer(file3, TYPE_FILE);
    register_pointer(file4, TYPE_FILE);
    register_pointer(file5, TYPE_FILE);
    register_pointer(file6, TYPE_FILE);
    register_pointer(file7, TYPE_FILE);

    // Read input_file
    char buffer[100]; // Buffer to store each line read from the file 
    fgets(buffer, sizeof(buffer), input_file);
    sscanf(buffer, "%*[^0-9]%d", &boundary_condition); // %*[^0-9] means to skip all non-numeric numbers
    fgets(buffer, sizeof(buffer), input_file);
    sscanf(buffer, "%*[^0-9]%d", &N);
    fgets(buffer, sizeof(buffer), input_file);
    sscanf(buffer, "%*[^0-9]%d", &inactive_cycle);
    fgets(buffer, sizeof(buffer), input_file);
    sscanf(buffer, "%*[^0-9]%d", &active_cycle);
    fgets(buffer, sizeof(buffer), input_file);
    sscanf(buffer, "%*[^:]: %[^\n]", particle_name);  
    fgets(buffer, sizeof(buffer), input_file); 
    sscanf(buffer, "%*[^:]: %lf", &particle_mass);
    fgets(buffer, sizeof(buffer), input_file);
    sscanf(buffer, "%*[^:]: %lf", &potential);
    fgets(buffer, sizeof(buffer), input_file);
    sscanf(buffer, "%*[^0-9]%lf", &scaling_factor);

    // Close input_file
    fclose(input_file); // close input_file

    // Check the boundary condition to simulate
    while(boundary_condition != 0 && boundary_condition != 1) {
        printf("Enter the boundary condition (0 = vacuum, 1 = reflective): ");
        if (scanf("%d", &boundary_condition) != 1) { 
            printf("Invalid input. Enter 0 (vacuum) or 1 (reflective).\n");
            while (getchar() != '\n'); // Clear input buffer by reading and removing characters until a newline is encountered
        } else if (boundary_condition != 0 && boundary_condition != 1) {
            printf("Invalid value. Boundary condition must be 0 (vacuum) or 1 (reflective).\n");
        } else {
            break; // Valid input, exit loop
        }
    }

    // Check the initial number of neutrons to simulate
    while (N <= 0) {
        printf("Enter the number of neutrons (must be greater than 0): ");
        if (scanf("%d", &N) != 1 || N <= 0) { 
            printf("Invalid input. The number of neutrons must be greater than 0.\n");
            while (getchar() != '\n'); 
        } else {
            break; 
        }
    }

    // Check the number of cycles to simulate
    M = inactive_cycle + active_cycle;
    while (inactive_cycle < 0 || active_cycle < 0 || active_cycle < inactive_cycle) {
        printf("Enter the number of cycles (inactive and active): ");
        if (scanf("%d %d", &inactive_cycle, &active_cycle) != 2 || inactive_cycle < 0 || active_cycle < 0 || active_cycle < inactive_cycle) {
            printf("Invalid input. Active cycles must be greater than or equal to inactive cycles, and both must be non-negative.\n");
            while (getchar() != '\n'); 
        } else {
            M = inactive_cycle + active_cycle;
            break; 
        } 
    }

    // Check the particle mass to simulate
    while (particle_mass <= 0) {
        printf("Enter the particle mass (must be greater than 0 and in scientific notation and SI unit): ");
        if (scanf("%lf", &particle_mass) != 1 || particle_mass <= 0) { 
            printf("Invalid input. The particle mass must be greater than 0.\n");
            while (getchar() != '\n'); 
        } else {
            break; 
        }
    }
    
    // Check the potential to simulate
    while (potential < 0) {
        printf("Enter the potential (must be greater than or equal to 0 and in scientific notation and SI unit): ");
        if (scanf("%lf", &potential) != 1 || potential < 0) { 
            printf("Invalid input. The potential must be greater than or equal to 0.\n");
            while (getchar() != '\n'); 
        } else {
            break; 
        }
    }

    // Change units to g and cm
    particle_mass *= 1000;
    potential *= 10000000;

    // Initialize RNG instance
    initialize_rng();

    // Initialize cross sections
    sigma_t = (2/(3*width))*scaling_factor;

    double r;
    if (potential == 0) {
        do {
            r = random_number_generator(); // r = [0,1)
        } while (r == 0); // avoid being extremely close to 0 and 1
        sigma_a = r * sigma_t;
    }
    else {
        sigma_a = ((particle_mass*width/(h_bar*h_bar))*potential)*(1.0/scaling_factor);
    }
    sigma_s = sigma_t - sigma_a;

    do {
        r = random_number_generator(); 
    } while (r == 0); 
    sigma_c = r * sigma_a;

    sigma_f = sigma_a - sigma_c;
    while (sigma_f < 0.02) {
        //printf("sigma_f is smaller than 0.02: %f\n", sigma_f);
        do {
            r = random_number_generator(); 
        } while (r == 0); 

        if (potential == 0) {
            sigma_a = r * sigma_t;
        }
        else {
            printf("sigma_f is smaller than 0.02: %f\n", sigma_f);
            exit(EXIT_FAILURE); // Terminate program if sigma_f is smaller than 0.02
        }

        do {
        r = random_number_generator(); 
        } while (r == 0); 
        sigma_c = r * sigma_a;
        sigma_f = sigma_a - sigma_c;
        sigma_s = sigma_t - sigma_a;
    }

    // Print input data 
    printf("\n");
    printf("===========================================================\n");
    printf("Boundary condition (0 = vacuum, 1 = reflective): %d\n", boundary_condition);
    printf("Number of neutrons: %d\n", N);
    printf("Total number of cycles: %d\n", M);
    printf("Number of inactive cycles: %d\n", inactive_cycle);
    printf("Number of active cycles: %d\n", active_cycle);
    printf("Particle: %s\n", particle_name);
    printf("Particle mass (g): %.5e\n", particle_mass);
    printf("Potential (g*cm^2/s) : %.5e\n", potential);
    printf("Scaling factor : %.5e\n", scaling_factor);
    printf("===========================================================\n");
    printf("sigma_t = %f\n", sigma_t);
    printf("sigma_a = %f\n", sigma_a);
    printf("sigma_s = %f\n", sigma_s);
    printf("sigma_c = %f\n", sigma_c);
    printf("sigma_f = %f\n", sigma_f);
    printf("===========================================================\n");

    // Calculate the theoretical eigenvalues for the first and second harmonic modes
    if (boundary_condition == 0) { // Vacuum (flux zero)
        k1 = (nu*sigma_f)/(((1/(3*(sigma_s+sigma_c+sigma_f)))*(M_PI/width)*(M_PI/width))+(sigma_c+sigma_f));
        k2 = (nu*sigma_f)/(((1/(3*(sigma_s+sigma_c+sigma_f)))*(M_PI/(width/2))*(M_PI/(width/2)))+(sigma_c+sigma_f));
    }
    else { // Reflective
        k1 = (nu*sigma_f)/(sigma_c+sigma_f);
        k2 = (nu*sigma_f)/(((1/(3*(sigma_s+sigma_c+sigma_f)))*(M_PI/width)*(M_PI/width))+(sigma_c+sigma_f));
    }
    printf("k1 = %f\n", k1);
    printf("k2 = %f\n", k2);

    // Allocate memory for fission banks
    parent_fission_bank = (neutron *) calloc(2 * N, sizeof(neutron));
    child_fission_bank = (neutron *) calloc(2 * N, sizeof(neutron)); 
    if (parent_fission_bank == NULL || child_fission_bank == NULL) {
        printf("Memory allocation failed for one or more fission banks.\n");
        exit(EXIT_FAILURE); 
    }
    register_pointer(parent_fission_bank, TYPE_CALLOC);
    register_pointer(child_fission_bank, TYPE_CALLOC);

    // Initialize the base (bottom) of the parent and child fission bank stacks
    parent_fission_bank[parent_neutron].position = -1.0;
    parent_fission_bank[parent_neutron].weight = 0.0;
    child_fission_bank[child_neutron].position = -1.0;
    child_fission_bank[child_neutron].weight = 0.0;

    // Initialize the parent fission bank
    initialize_fission_bank(parent_fission_bank_double_pointer, parent_pointer);

    // Allocate memory to store initial position of parent neutrons
    gsl_vector *initial_neutron_flux = gsl_vector_calloc(num_bins); 
    if (initial_neutron_flux == NULL) {
        printf("Memory allocation failed for initial neutron positions.\n");
        exit(EXIT_FAILURE);
    }

    // Tally the positions of parent neutrons into initial_neutron_flux
    position_tally(initial_neutron_flux, parent_fission_bank, parent_neutron);

    // Debug Start: Logging initial neutron flux and parent neutron positions
    // Log initial neutron flux data
    double n_avg = (double)N/num_bins; // Average number of neutrons in each bin
    fprintf(file1, "----------- initial neutron flux -----------\n");
    for (size_t i = 0; i < initial_neutron_flux->size; i++) {
        fprintf(file1, "%f %g\n", (width/(2*num_bins))+i*(width/num_bins), gsl_vector_get(initial_neutron_flux, i)/n_avg); // Write the bin index and corresponding neutron flux value
    }

    // Free initial_neutron_flux
    gsl_vector_free(initial_neutron_flux);

    // Log initial parent neutron positions
    fprintf(file2, "initial parent neutron numbers : %d\n", parent_neutron);
    for(int i = 0; i < N+1; i++) {
        fprintf(file2, "%dth neutron: %f \n", i+1, parent_fission_bank[i].position);
    }

    // Log the start of the simulation
    fprintf(file2, "\n----------- starting simulation -----------\n");
    // Debug End

    // Simulate over M cycles
    double k = 1.0; // Initial k value for scaling => nu*simga_f/sigma_a 로 시작해보기 
    double k_sum = 0.0; // Accumulated sum of k values
    double k_squared_sum = 0.0; // Accumulated sum of the squares of k values

    double k_average; // Average of k values over all cycles
    double k_squared_average; // Average of squared k values over all cycles
    double k_sample_standard_deviation; // Sample standard deviation of k values over all cycles

    double E_sum = 0.0; // Accumulated sum of E values
    double E_squared_sum = 0.0; // Accumulated sum of the squares of E values

    double E_average; // Average of E values over all cycles
    double E_squared_average; // Average of squared E values over all cycles
    double E_sample_standard_deviation; // Sample standard deviation of E values over all cycles

    gsl_vector *neutron_flux_sum = gsl_vector_calloc(num_bins); // Vector accumulating the sum of neutron flux values over all cycles
    gsl_vector *neutron_flux_squared_sum = gsl_vector_calloc(num_bins); // Vector accumulating the sum of the square of neutron flux values over all cycles
    if (neutron_flux_sum == NULL || neutron_flux_squared_sum == NULL) {
        printf("Memory allocation failed for one or more vectors.\n");
        exit(EXIT_FAILURE); 
    }
    register_pointer(neutron_flux_sum, TYPE_GSL_VECTOR);
    register_pointer(neutron_flux_squared_sum, TYPE_GSL_VECTOR);

    gsl_vector *initial_position_vector = gsl_vector_calloc(num_bins); // Vector tracking the initial positions of parent neutrons across cycles
    if (initial_position_vector == NULL) {
        printf("Memory allocation failed for the initial position vector.\n");
        exit(EXIT_FAILURE); 
    }
    register_pointer(initial_position_vector, TYPE_GSL_VECTOR);

    gsl_matrix *A_sum = gsl_matrix_calloc(num_bins, num_bins); // Matrix to store the sum of all A matrices over all cycles
    if (A_sum == NULL) {
        printf("Memory allocation failed for A_sum matrix.\n");
        exit(EXIT_FAILURE);
    }
    register_pointer(A_sum, TYPE_GSL_MATRIX);

    for(int i = 0; i < M; i++) {
        double k_cycle = 0.0; // k value for a single cycle
        double *k_pointer = &k_cycle; // Pointer to k_cycle 
        int size_of_fb = parent_neutron; // Number of neutrons in parent_fission_bank
        double total_weight = 0.0; // Total weight of parent neutrons for a single cycle 

        gsl_matrix *A = gsl_matrix_calloc(num_bins, num_bins); // matrix A storing data for each cycle
        if (A == NULL) {
            printf("Memory allocation failed for matrix A.\n");
            exit(EXIT_FAILURE); 
        }
        register_pointer(A, TYPE_GSL_MATRIX);

        sum_of_weights(parent_fission_bank, &total_weight, parent_neutron); // Total sum of neutron weights in parent_fission_bank
        position_tally(initial_position_vector, parent_fission_bank, parent_neutron); // Tally positions of parent neutrons into initial_position_vector

        // Log total weight of parent neutrons before simulation
        fprintf(file2, "----------- cycle# = %d, size of fb = %d -----------\n", i+1, size_of_fb);
        fprintf(file2, "total weight = %f\n", total_weight);

        // Simulate parent fission bank
        gsl_vector *neutron_flux = simulate_neutron_diffusion(file2, parent_fission_bank, child_fission_bank_double_pointer, k_pointer, k, A);
        if (neutron_flux == NULL) {
            printf("Memory allocation failed for neutron flux.\n");
            exit(EXIT_FAILURE);
        }
        register_pointer(neutron_flux, TYPE_GSL_VECTOR);

        // Scale k value and normalize neutron flux per neutron per bin width
        k_cycle /= total_weight; // k value for a single cycle
        gsl_vector_scale(neutron_flux, 1.0/(size_of_fb*(width/(double)num_bins))); // Neutron flux for one cycle per neutron per width of bin

        // Accumulate matrix A over all cycles
        gsl_matrix_add(A_sum, A); // A_sum += A

        // Debug Start: Logging neutron flux, k values, matrices, and fission bank data after the simulation 
        // Log neutron flux values
        fprintf(file1, "----------- cycle# = %d -----------\n", i+1);
        for (size_t i = 0; i < neutron_flux_squared_sum->size; i++) {
            fprintf(file1, "%f %g\n", (width/(2.0*num_bins))+i*(width/num_bins), gsl_vector_get(neutron_flux, i));
        }

        // Log parent and child fission bank size
        fprintf(file2, "++++++++++++++++++++++++++++++++++++++\n");
        fprintf(file2, "size of fb, before / after = %d, %d\n", size_of_fb, child_neutron);
        fprintf(file2, "k_cycle = %10.5f\n", k_cycle);
        fprintf(file2, "++++++++++++++++++++++++++++++++++++++\n\n");
        // Debug End

        // Accumulate k value only for active cycles
        if (i >= inactive_cycle) {
            k_sum += k_cycle;
            k_squared_sum += k_cycle * k_cycle;

            if (potential == 0) {
                double temp = ((nu*sigma_f/k_cycle)-sigma_a)*(h_bar*h_bar/(particle_mass * width))*scaling_factor;
                E_sum += temp;
                E_squared_sum += temp * temp;
            } else {
                double temp = ((nu * sigma_f * h_bar * h_bar)/(particle_mass * width * k_cycle))*scaling_factor;
                E_sum += temp;
                E_squared_sum += temp * temp;
            }

            gsl_vector_add(neutron_flux_sum, neutron_flux); // neutron_flux_sum += neutron_flux
            gsl_vector_mul(neutron_flux, neutron_flux); // neutron_flux *= neutron_flux
            gsl_vector_add(neutron_flux_squared_sum, neutron_flux); // neutron_flux_squared_sum += neutron_flux
        }

        // Terminate program if child fission bank is empty (no neutrons left)
        if (child_neutron == 0) {
            fprintf(file2, "Fission bank is empty\n");
            printf("Fission bank is empty\n");
            exit(EXIT_FAILURE);
        } 

        // Reweight neutrons to conserve the total weight of every cycle
        reweight_neutrons(boundary_condition, child_fission_bank_double_pointer);

        // Swap parent and child fission banks for the next cycle
        neutron *temp = parent_fission_bank; 
        parent_fission_bank = child_fission_bank; 
        child_fission_bank = temp; 

        // Swap parent and child neutron counts for the next cycle
        int tempp = parent_neutron;
        parent_neutron = child_neutron;
        child_neutron = tempp; 

        // Update k value
        k = k_cycle; 

        // Free neutron flux and matrix A
        unregister_pointer(neutron_flux);
        unregister_pointer(A);
    }

    // Calculate the mean and standard deviation of the neutron flux
    gsl_vector_scale(neutron_flux_sum, 1.0/active_cycle); // Compute mean neutron flux
    gsl_vector *mean_neutron_flux = gsl_vector_calloc(neutron_flux_sum->size);
    if (mean_neutron_flux == NULL) {
        printf("Memory allocation failed for mean neutron flux.\n");
        exit(EXIT_FAILURE);
    }
    register_pointer(mean_neutron_flux, TYPE_GSL_VECTOR);

    gsl_vector_memcpy(mean_neutron_flux, neutron_flux_sum);
    gsl_vector_scale(neutron_flux_squared_sum, 1.0/active_cycle); // Compute mean squared neutron flux
    
    gsl_vector_mul(neutron_flux_sum, neutron_flux_sum); // Square mean neutron flux
    gsl_vector_sub(neutron_flux_squared_sum, neutron_flux_sum); // Compute variance
    
    vector_square_root(neutron_flux_squared_sum); // Compute standard deviation 
    
    // Scale by 1/sqrt(active_cycle - 1) for sample standard deviation (Bessel's correction)
    if (active_cycle > 1) {
        gsl_vector_scale(neutron_flux_squared_sum, 1.0/sqrt(active_cycle-1));
    } else {
        printf("Cannot calculate the sample standard deviation of neutron flux. The number of active cycles must be greater than 1.\n");
        gsl_vector_set_all(neutron_flux_squared_sum, 0.0);
    }

    // Calculate the mean and standard deviation of the k value and energy
    k_average = k_sum/active_cycle;
    k_squared_average = k_squared_sum/active_cycle;
    E_average = E_sum/active_cycle;
    E_squared_average = E_squared_sum/active_cycle;

    // Compute sample standard deviation using Bessel's correction
    if (active_cycle > 1) {
        k_sample_standard_deviation = sqrt((k_squared_average - k_average * k_average) / (active_cycle - 1));
        E_sample_standard_deviation = sqrt((E_squared_average - E_average * E_average) / (active_cycle - 1));
    } else {
        printf("Cannot calculate the sample standard deviation of k and E. The number of active cycles must be greater than 1.\n");
        k_sample_standard_deviation = 0.0; 
        E_sample_standard_deviation = 0.0; 
    }

    // Calculate the fission matrix
    scale_columns(A_sum, initial_position_vector);

    // Calculate the theoretical and simulated ground state energy of the particle
    //double E_theoretical = (h_bar*h_bar*PI*PI)/(2*particle_mass*width*width); // infinite potential well 
    double E_theoretical = (h_bar*h_bar*PI*PI)/(2*particle_mass*width*width) + potential;

    // Debug Start: Logging final neutron flux, k values, and fission matrix
    // Log mean and standard deviation of neutron flux
    fprintf(file1, "----------- mean -----------\n");
    for (size_t i = 0; i < mean_neutron_flux->size; i++) {
        fprintf(file1, "%f %g\n", (width/(2.0*num_bins))+i*(width/num_bins), gsl_vector_get(mean_neutron_flux, i)); // write the index and value
    }
    fprintf(file1, "----------- standard deviation -----------\n");
    for (size_t i = 0; i < neutron_flux_squared_sum->size; i++) {
        fprintf(file1, "%f %g\n", (width/(2.0*num_bins))+i*(width/num_bins), gsl_vector_get(neutron_flux_squared_sum, i)); // write the index and value
    }
    
    // Log mean and standard deviation of k, theoretical k values for the first and second harmonic modes
    fprintf(file2, "k_avg=%f, k_sa=%f, std=%f\n", k_average, k_squared_average, k_sample_standard_deviation);
    fprintf(file2, "k1=%f, k2=%f, k2/k1=%f\n", k1, k2, k2/k1);
    fprintf(file2, "E_avg=%f, E_sa=%f, std=%f\n", E_average, E_squared_average, E_sample_standard_deviation);
    fprintf(file2, "E_theoretical=%f\n", E_theoretical);

    printf("===========================================================\n");
    printf("Theoretical values: k1=%f, k2=%f, k2/k1=%f\n", k1, k2, k2/k1);
    printf("Simulated values: k_avg=%f, k_sa=%f, std=%f\n", k_average, k_squared_average, k_sample_standard_deviation);
    printf("Theoretical ground state enegry: %e\n", E_theoretical);
    printf("Simulated ground state enegry: E_avg=%0.5e, E_sa=%f, std=%f\n", E_average, E_squared_average, E_sample_standard_deviation);
    printf("===========================================================\n");

    // Save simulated ground state energy to a file
    FILE *energy_file = fopen("electron_infinite.txt", "a");
    if (energy_file == NULL) {
        printf("Unable to open electron_infinite.txt for writing.\n");
        exit(EXIT_FAILURE);
    }
    fprintf(energy_file, "sigma_s: %f, sigma_c: %f, simga_f: %f, k1: %f, k_avg: %f, std: %f, E_avg: %0.5e, E_sa: %0.5e, std: %0.5e\n", 
        sigma_s, sigma_c, sigma_f, k1, k_average, k_sample_standard_deviation, E_average, E_squared_average, E_sample_standard_deviation);
    fclose(energy_file);
    // Debug End

    // Log fission matrix
    for (size_t i = 0; i < A_sum ->size1; i++) {
        for (size_t j = 0; j < A_sum ->size2; j++) {
            fprintf(file7, "%.5f ", gsl_matrix_get(A_sum, i, j));
            printf("%.5f ", gsl_matrix_get(A_sum, i, j));
        }
        fprintf(file7, "\n");
        printf("\n");
    }
   // Debug End

    // Cleanup
    return 0;
}

