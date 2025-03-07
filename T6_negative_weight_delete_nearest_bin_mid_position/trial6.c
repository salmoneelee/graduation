#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include "trial6.h"
#include <math.h>
#include <time.h>

int main() {
    int M; // Total number of cycles (inactive + active)
    int active_cycle; // Number of active cycles 
    int inactive_cycle; // Number of inactive cycles

    neutron *parent_fission_bank = NULL; // Pointer to the parent fission bank (stack)
    neutron *child_fission_bank = NULL; // Pointer to the child fission bank (stack)
    neutron **parent_fission_bank_double_pointer = &parent_fission_bank; // Double pointer to the parent fission bank 
    neutron **child_fission_bank_double_pointer = &child_fission_bank; // Double pointer to the child fission bank
    
    double k1; // Theoretical eigenvalue for the first harmonic mode 
    double k2; // Theoretical eigenvalue for the second harmonic mode
    
    const char *input = "input_6.txt"; // File for reading the inputs
    const char *filename1 = "flux_data_6.txt"; // File for storing neutron flux data
    const char *filename2 = "text_data_6.txt"; // File for storing data of each cycle
    const char *filename3 = "position_data_6.txt"; // File for storing neutron position data
    const char *filename4 = "remove_data_6.txt"; // File for storing neutron data related to removal

    // Register cleanup function
    atexit(cleanup); 

    // Open text files
    FILE *input_file = fopen(input, "r");
    FILE *file1 = fopen(filename1, "w");
    FILE *file2 = fopen(filename2, "w");
    FILE *file3 = fopen(filename3, "w");
    FILE *file4 = fopen(filename4, "w");
    if(input_file == NULL || file1 == NULL || file2 == NULL || file3 == NULL || file4 == NULL) {
        printf("Unable to open files.\n");
        exit(EXIT_FAILURE); // Terminate program immediatedly
    }
    register_pointer(file1, TYPE_FILE);
    register_pointer(file2, TYPE_FILE);
    register_pointer(file3, TYPE_FILE);
    register_pointer(file4, TYPE_FILE);

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

    // Close input_file
    fclose(input_file); // close input_file

    // Calculate the theoretical eigenvalues for the first and second harmonic modes
    if (boundary_condition == 0) { // Vacuum (flux zero)
        k1 = (nu*sigma_f)/(((1/(3*(sigma_s+sigma_c+sigma_f)))*(M_PI/width)*(M_PI/width))+(sigma_c+sigma_f));
        k2 = (nu*sigma_f)/(((1/(3*(sigma_s+sigma_c+sigma_f)))*(M_PI/(width/2))*(M_PI/(width/2)))+(sigma_c+sigma_f));
    }
    else { // Reflective
        k1 = (nu*sigma_f)/(sigma_c+sigma_f);
        k2 = (nu*sigma_f)/(((1/(3*(sigma_s+sigma_c+sigma_f)))*(M_PI/width)*(M_PI/width))+(sigma_c+sigma_f));
    }

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
            while (getchar() != '\n'); // Clear input buffer by reading and removing characters until a newline is encountered
        } else {
            break; // Valid input, exit loop
        }
    }

    // Check the number of cycles to simulate
    M = inactive_cycle + active_cycle;
    while (inactive_cycle < 0 || active_cycle < 0 || active_cycle < inactive_cycle) {
        printf("Enter the number of cycles (inactive and active): ");
        if (scanf("%d %d", &inactive_cycle, &active_cycle) != 2 || inactive_cycle < 0 || active_cycle < 0 || active_cycle < inactive_cycle) {
            printf("Invalid input: Active cycles must be greater than or equal to inactive cycles, and both must be non-negative.\n");
            while (getchar() != '\n'); // Clear input buffer by reading and removing characters until a newline is encountered
        } else {
            M = inactive_cycle + active_cycle;
            break; // Valid input, exit loop
        } 
    }

    // Print input data 
    printf("\n");
    printf("Boundary condition (0 = vacuum, 1 = reflective): %d\n", boundary_condition);
    printf("Number of neutrons: %d\n", N);
    printf("Total number of cycles: %d\n", M);
    printf("Number of inactive cycles: %d\n", inactive_cycle);
    printf("Number of active cycles: %d\n", active_cycle);

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
    
    // Initialize RNG instance
    initialize_rng();

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
    double k = 1.0; // Initial k value for scaling
    double k_sum = 0.0; // Accumulated sum of k values
    double k_squared_sum = 0.0; // Accumulated sum of the squares of k values

    double k_average; // Average of k values over all cycles
    double k_squared_average; // Average of squared k values over all cycles
    double k_sample_standard_deviation; // Sample standard deviation of k values over all cycles

    gsl_vector *neutron_flux_sum = gsl_vector_calloc(num_bins); // Vector accumulating the sum of neutron flux values over all cycles
    gsl_vector *neutron_flux_squared_sum = gsl_vector_calloc(num_bins); // Vector accumulating the sum of the square of neutron flux values over all cycles
    if (neutron_flux_sum == NULL || neutron_flux_squared_sum == NULL) {
        printf("Memory allocation failed for one or more vectors.\n");
        exit(EXIT_FAILURE); 
    }
    register_pointer(neutron_flux_sum, TYPE_GSL_VECTOR);
    register_pointer(neutron_flux_squared_sum, TYPE_GSL_VECTOR);

    for(int i = 0; i < M; i++) {
        double k_cycle = 0.0; // k value for a single cycle
        double *k_pointer = &k_cycle; // Pointer to k_cycle 
        int size_of_fb = parent_neutron; // Number of neutrons in parent_fission_bank
        double total_weight = 0.0; // Total weight of parent neutrons for a single cycle 

        sum_of_weights(parent_fission_bank, &total_weight, parent_neutron); // Total sum of neutron weights in parent_fission_bank

        // Debug Start: Logging parent neutron, k values, and initial position vector before simulation
        // Log parent neutrons to simulate
        fprintf(file3, "----------- cycle# = %d -----------\n", i+1);
        fprintf(file3, "----------- parent -----------\n");
        for (int i = 1; i < parent_neutron + 1; i++) {
            fprintf(file3, "%d %0.10f %0.10f\n", i, parent_fission_bank[i].position, parent_fission_bank[i].weight); 
        }

        // Log total weight of parent neutrons before simulation
        fprintf(file2, "----------- cycle# = %d, size of fb = %d -----------\n", i+1, size_of_fb);
        fprintf(file2, "total weight = %f\n", total_weight);
        // Debug End
        
        // Simulate parent fission bank
        gsl_vector *neutron_flux = simulate_neutron_diffusion(file2, parent_fission_bank, child_fission_bank_double_pointer, k_pointer, k);
        if (neutron_flux == NULL) {
            printf("Memory allocation failed for neutron flux.\n");
            exit(EXIT_FAILURE);
        }
        register_pointer(neutron_flux, TYPE_GSL_VECTOR);

        // Scale k value and normalize neutron flux per neutron per bin width
        k_cycle /= total_weight; // k value for a single cycle
        gsl_vector_scale(neutron_flux, 1.0/(size_of_fb*(width/(double)num_bins))); // Neutron flux for one cycle per neutron per width of bin

        // Debug Start: Logging neutron flux, k values, and matrices after the simulation 
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

        if (boundary_condition == 0) {
            // Log the position and weight of child neutrons before removal
            fprintf(file4, "----------- child before removal -----------\n");
            for (int i = 1; i < child_neutron + 1; i++) {
                fprintf(file4, "%d %0.10f %0.10f\n", i, child_fission_bank[i].position, child_fission_bank[i].weight); 
            }
        
            // Remove positive weight neutrons using negative weight neutrons within the child fission bank
            remove_and_redistribute_neutrons(file4, child_fission_bank_double_pointer);

            // Log the position and weight of child neutrons after removal
            fprintf(file4, "----------- child after removal-----------\n");
            for (int i = 1; i < child_neutron + 1; i++) {
                fprintf(file4, "%d %0.10f %0.10f\n", i, child_fission_bank[i].position, child_fission_bank[i].weight); 
            }
        }

        // Debug Start: Logging fission bank data after the simulation 
        // Log child fission bank data 
        fprintf(file3, "----------- child -----------\n");
        for (int i = 1; i < child_neutron + 1; i++) {
            fprintf(file3, "%d %0.10f %0.10f\n", i, child_fission_bank[i].position, child_fission_bank[i].weight); 
        }
        // Debug End

        // Adjust the size of the child_fission_bank to N
        fission_bank_size_adjustment(file3, child_fission_bank_double_pointer); 

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

    // Calculate the mean and standard deviation of the k value
    k_average = k_sum/active_cycle;
    k_squared_average = k_squared_sum/active_cycle;

    // Compute sample standard deviation using Bessel's correction
    if (active_cycle > 1) {
        k_sample_standard_deviation = sqrt((k_squared_average - k_average * k_average) / (active_cycle - 1));
    } else {
        printf("Cannot calculate the sample standard deviation of k. The number of active cycles must be greater than 1.\n");
        k_sample_standard_deviation = 0.0; 
    }

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
    fprintf(file2, "k_a=%f, k_sa=%f, std=%f\n", k_average, k_squared_average, k_sample_standard_deviation);
    fprintf(file2, "k1=%f, k2=%f, k2/k1=%f\n", k1, k2, k2/k1);
    printf("k_a=%f, k_sa=%f, std=%f\n", k_average, k_squared_average, k_sample_standard_deviation);
    printf("k1=%f, k2=%f, k2/k1=%f\n\n", k1, k2, k2/k1);
   // Debug End

    // Cleanup
    return 0;
}

