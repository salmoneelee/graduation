#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include "trial1.h"
#include <math.h>
#include <string.h>

int N; // Total number of neutrons
int boundary_condition; // Boundary condition (0 = vacuum, 1 = reflective)

const double width = 20.0; // Width of the reactor
const int num_bins = 20; // Number of bins in the reactor

const double sigma_s = 0.1; // Scattering cross section
const double sigma_c = 0.07; // Capture cross section
const double sigma_f = 0.06; // Fission cross section 
const double nu = 2.43; // Average number of neutrons produced per fission

int parent_neutron = 0; // Number of neutrons in the parent fission bank
int child_neutron = 0; // Number of neutrons in the child fission bank
int *parent_pointer = &parent_neutron; // Pointer to parent_neutron
int *child_pointer = &child_neutron; // Pointer to child_neutron

cleanup_entry chained_pointers[MAX_POINTERS]; // Array of pointers that need to be freed or closed
int chain_count = 0; // Number of pointers that need to be freed or closed

const gsl_rng_type *T = NULL; // Type of GSL RNG (random number generator)
gsl_rng *r = NULL; // Pointer to GSL RNG instance

// Add pointer to the cleanup list
void register_pointer(void *ptr, int type) {
    if (ptr == NULL || chain_count >= MAX_POINTERS) return;

    for (int i = 0; i < chain_count; i++) {
        if (chained_pointers[i].data == ptr) {
            printf("Warning: Pointer %p is already registered! (Type: %d)\n", ptr, chained_pointers[i].type);
            return;  // Prevent duplicate registration
        }
    }

    chained_pointers[chain_count].data = ptr;
    chained_pointers[chain_count].type = type;
    chain_count++;
}

// Remove pointer from cleanup list and free memory
void unregister_pointer(void *ptr) {
    if (ptr == NULL) return;  // Ignore NULL pointers

    for (int i = 0; i < chain_count; i++) {
        if (chained_pointers[i].data == ptr) {
            if (chained_pointers[i].data == NULL) {
                printf("Warning: Attempted to free already freed pointer: %p\n", ptr);
                return;  // Prevent double-free
            }

            // Free the pointer based on its type
            switch (chained_pointers[i].type) {
                case TYPE_FILE:
                    fclose((FILE *)ptr);
                    break;
                case TYPE_CALLOC:
                    free(ptr);
                    break;
                case TYPE_GSL_VECTOR:
                    gsl_vector_free((gsl_vector *)ptr);
                    break;
                case TYPE_RNG:
                    free_rng();
                    break;
            }

            // Set pointer to NULL to prevent double-free
            chained_pointers[i].data = NULL;

            // Shift remaining elements to fill the gap
            for (int j = i; j < chain_count - 1; j++) {
                chained_pointers[j] = chained_pointers[j + 1];
            }

            chain_count--;  // Reduce count
            return;  // Exit after removing
        }
    }
}

// Cleanup function for exit()
void cleanup() {
    for (int i = 0; i < chain_count; i++) {
        if (chained_pointers[i].data == NULL) continue;

        switch (chained_pointers[i].type) {
            case TYPE_FILE:
                fclose((FILE *)chained_pointers[i].data);
                break;
            case TYPE_CALLOC:
                free(chained_pointers[i].data);
                break;
            case TYPE_GSL_VECTOR:
                gsl_vector_free((gsl_vector *)chained_pointers[i].data);
                break;
            case TYPE_RNG:
                free_rng();
                break;
        }

        chained_pointers[i].data = NULL; // Avoid double free
    }
    chain_count = 0; // Reset count after cleanup
}


// Initialize RNG instance
void initialize_rng(void) {
    T = gsl_rng_default; // Assign the default RNG algorithm (Mersenne Twister = gsl_rng_mt19937)
    r = gsl_rng_alloc(T); // Allocate the RNG instance based on T
    if (r == NULL) {
        printf("RNG allocation failed.\n");
        exit(EXIT_FAILURE);
    }
    register_pointer(r, TYPE_RNG);
    //gsl_rng_set(r, time(NULL)); // Set seed to current time
    gsl_rng_set(r, 10); // Set a fixed seed
}

// Generate a uniformly distributed random number in the range [0, 1)
double random_number_generator(void) {
    if (r == NULL) {
        printf("RNG not initialized. Call initialize_rng() first.\n");
        exit(EXIT_FAILURE);
    }
    return gsl_rng_uniform(r); 
}

// Free the RNG instance 
void free_rng(void) {
    if (r != NULL) {
        gsl_rng_free(r);
        r = NULL; // Prevent accidental reuse
        T = NULL; // Reset RNG type pointer 
    }
}

// Initialize fission bank
void initialize_fission_bank(neutron **fission_bank_double_pointer, int *index_pointer) {
    if (fission_bank_double_pointer == NULL || *fission_bank_double_pointer == NULL || index_pointer == NULL) {
        printf("Null pointer in initialize_fission_bank function.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; i++) {
        (*index_pointer)++;
        (*fission_bank_double_pointer)[*index_pointer].position = width * random_number_generator();
        (*fission_bank_double_pointer)[*index_pointer].direction = determine_direction();
    }
}

// Tally neutron weights into position bins in the vector
void position_tally(gsl_vector *vector, neutron *fission_bank_pointer, int neutron) {
    if (vector == NULL || fission_bank_pointer == NULL) {
        printf("Null pointer in position_tally function.\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i < neutron+1; i++) {
        int bin_number = (int)floor(fission_bank_pointer[i].position/(width/(double)num_bins));

        if (bin_number < 0) {
            printf("Neutron is out of bin range.\n");
            bin_number = 0;
        }
        if (bin_number == num_bins) bin_number = num_bins - 1;

        double current_value = gsl_vector_get(vector, bin_number);
        gsl_vector_set(vector, bin_number, current_value + 1.0);
    }
}

// Determine direction of neutron
double determine_direction(void) {
    double random_number= random_number_generator();
    if (random_number < 0.5) {
        return -1;
    } 
    else {
        return 1;
    }
}

// Determine distance of neutron
double determine_distance(void) {
    double random_number= random_number_generator();
    return (-log(random_number)/(sigma_s + sigma_c + sigma_f));
}

// Simulate neutron diffusion for one cycle
gsl_vector *simulate_neutron_diffusion(FILE *file2, neutron *parent_fission_bank, neutron **child_fission_bank_double_pointer, double *k_pointer) {
    gsl_vector *flux = gsl_vector_calloc(num_bins); 
    if (flux == NULL) {
        printf("Memory allocation for flux failed.\n");
        exit(EXIT_FAILURE);
    }

    while (parent_fission_bank[parent_neutron].position != -1.0) {
        double initial_position = parent_fission_bank[parent_neutron].position;
        double distance = determine_distance(); // Distance to the next collision

        // Debug Start: Logging simulation step
        fprintf(file2, "---- simulating neutron# %d ----\n ", parent_neutron);
        fprintf(file2, "pos: %f, dis: %f, direction: %d ",initial_position, distance, parent_fission_bank[parent_neutron].direction);
        // Debug End

        int fission_neutrons = 0; // Number of neutrons produced by fission
        double new_position;

        if (parent_fission_bank[parent_neutron].direction > 0) {
            new_position = parent_fission_bank[parent_neutron].position + distance;
            parent_fission_bank[parent_neutron].position = new_position;

            if (new_position < width && new_position > 0) {
                fprintf(file2, "newpos: %f\n", new_position);

                sample_collision_type(file2, parent_fission_bank, &fission_neutrons);
                add_fission_bank(child_fission_bank_double_pointer, fission_neutrons, new_position);
                
                k_tally(k_pointer, &fission_neutrons);
                flux_tally(flux, new_position);
            }
                
            else {
                if (boundary_condition == 0) { // vacuum
                    fprintf(file2, "\n");

                    parent_neutron--;
                }
                else { // reflective
                    if (((int)(new_position/width))%2 != 0) {
                        new_position = width - fmod(new_position, width);
                    }
                    else {
                        new_position = fmod(new_position, width);
                    }
                
                    parent_fission_bank[parent_neutron].position = new_position;
                    fprintf(file2, "newpos: %f\n", new_position);

                    sample_collision_type(file2, parent_fission_bank, &fission_neutrons);
                    add_fission_bank(child_fission_bank_double_pointer, fission_neutrons, new_position);
                    
                    k_tally(k_pointer, &fission_neutrons);
                    flux_tally(flux, new_position);
                }
            }
        }
        else {
            new_position = parent_fission_bank[parent_neutron].position - distance;
            parent_fission_bank[parent_neutron].position = new_position;

            if (new_position < width && new_position > 0){
                fprintf(file2, "newpos: %f\n", new_position);

                sample_collision_type(file2, parent_fission_bank, &fission_neutrons);
                add_fission_bank(child_fission_bank_double_pointer, fission_neutrons, new_position);
                
                k_tally(k_pointer, &fission_neutrons);
                flux_tally(flux, new_position);
            }
            else {
                if (boundary_condition == 0) {
                    fprintf(file2, "\n");

                    parent_neutron--;
                }
                else {
                    if (((int)(new_position/width))%2 != 0) {
                        new_position = width + fmod(new_position, width);
                    }
                    else {
                        new_position = -fmod(new_position, width);
                    }

                    parent_fission_bank[parent_neutron].position = new_position;
                    fprintf(file2, "newpos: %f\n", new_position);

                    sample_collision_type(file2, parent_fission_bank, &fission_neutrons);
                    add_fission_bank(child_fission_bank_double_pointer, fission_neutrons, new_position);
                    
                    k_tally(k_pointer, &fission_neutrons);
                    flux_tally(flux, new_position);
                }
            }
        }
    }   
    return flux;
}

// Sample collision type
void sample_collision_type(FILE *file2, neutron* fission_bank, int *neutrons) {
    double random_number_1 = random_number_generator(); 

    if (random_number_1 <= sigma_f/(sigma_s + sigma_c + sigma_f)) { 
        fission(file2, neutrons);
    } 
    else if (random_number_1 <= (sigma_f + sigma_s)/(sigma_s + sigma_c + sigma_f)) {
        scattering(fission_bank);
    }
    else {
        capture(file2);
    }
}

// Simulate fission
void fission(FILE *file2, int *neutrons) {
    double random_number = random_number_generator();

    if (random_number <= nu - floor(nu)) { // 3 new neutrons are produced
        (*neutrons) = floor(nu)+1.0;

        for (int i = 0; i < (*neutrons); i++) {
            fprintf(file2, "FISSION!!\n");
        }
    }
    else { 
        (*neutrons) = floor(nu);

        for (int i = 0; i < (*neutrons); i++) {
            fprintf(file2, "FISSION!!\n");
        }
    }
    parent_neutron--;
}

// Simulate scattering
void scattering(neutron* fission_bank) {
    double new_position = fission_bank[parent_neutron].position;
    fission_bank[parent_neutron].direction = determine_direction();
    if (boundary_condition == 1){ // reflective
        if (new_position == 0.0 && fission_bank[parent_neutron].direction < 0) fission_bank[parent_neutron].direction = 1; // re-enter the well
        else if (new_position == width && fission_bank[parent_neutron].direction > 0) fission_bank[parent_neutron].direction = -1; // re-enter the well
    }
}

// Simulate capture
void capture(FILE *file2) {
    fprintf(file2, "CAPTURE!!\n");
    (*parent_pointer)--; 
}

// Add fission-produced neutrons to the fission bank
void add_fission_bank(neutron **fission_bank_double_pointer, int n, double position) {
    neutron *new_fission_bank_pointer = NULL;

    if ((*child_pointer) + n >= 2 * N - 1) {
        size_t new_size = ((*child_pointer) == 1) ? 10 * (*child_pointer) : 2 * (*child_pointer);

        new_fission_bank_pointer = realloc((*fission_bank_double_pointer), new_size * sizeof(neutron));
        if (new_fission_bank_pointer == NULL) {
            printf("Memory reallocation for fission bank failed.\n");
            exit(EXIT_FAILURE); // 
        }

        size_t new_bytes = (new_size - ((*child_pointer) + 1)) * sizeof(neutron); // Memory that is not initialized 
        memset((char *)new_fission_bank_pointer + (((*child_pointer) + 1) * sizeof(neutron)), 0, new_bytes); // Initialize memory with 0x00 = initialize to 0.0 and 0

        for (int i = 0; i < chain_count; i++) {
            if (chained_pointers[i].data == *fission_bank_double_pointer) {
                chained_pointers[i].data = new_fission_bank_pointer; // Update chained pointers
                break;
            }
        }

        (*fission_bank_double_pointer) = new_fission_bank_pointer; // Old fission bank memory is automatically freed
    }

    for (int i = 0; i < n; i++) {
        (*child_pointer)++;
        (*fission_bank_double_pointer)[*child_pointer].position = position;
        (*fission_bank_double_pointer)[*child_pointer].direction = determine_direction();
    }
}

// Update the k value for collision
void k_tally(double *k_pointer, int *neutrons) {
    *k_pointer += (*neutrons);
}

// Update the neutron flux for collision
void flux_tally(gsl_vector *flux, double position) {
    int bin_number = (int)floor(position/(width/num_bins)); // Bin number of the position
    if (bin_number < 0) {
        printf("Neutron is out of bin range.\n");
        bin_number = 0;
    }
    if (bin_number == num_bins) {
        bin_number = num_bins - 1; 
    }

    double current_value = gsl_vector_get(flux, bin_number); // Current value of the flux vector at bin_number element
    gsl_vector_set(flux, bin_number, current_value + 1/(sigma_s + sigma_c + sigma_f)); // Set element at index bin_number of the flux vector to current_value +1/(sigma_s + sigma_c + sigma_f)
}

// Adjust the size of the child_fission_bank to N
void fission_bank_size_adjustment(neutron **child_fission_bank_double_pointer) {
    neutron *new_fission_bank = calloc(2 * N, sizeof(neutron)); // for the size adjusted child_fission_bank
    if (new_fission_bank == NULL) {
        printf("Memory allocation for size adjustment failed.\n");
        exit(EXIT_FAILURE);
    }
    int index; 

    new_fission_bank[0].position = -1.0;
    new_fission_bank[0].direction = 0;
    for (int i = 1; i < N+1; i++)
    { 
        index = ((int)(random_number_generator() * child_neutron)) ; // [0, child_neutron -1]

        new_fission_bank[i].position = (*child_fission_bank_double_pointer)[index+1].position;
        new_fission_bank[i].direction = (*child_fission_bank_double_pointer)[index+1].direction;
    }

    child_neutron = N;

    for (int i = 0; i < chain_count; i++) {
        if (chained_pointers[i].data == *child_fission_bank_double_pointer) {
            chained_pointers[i].data = new_fission_bank;  // update chained_pointers
            break;
        }
    }

    free((*child_fission_bank_double_pointer));
    (*child_fission_bank_double_pointer) = new_fission_bank;
}

// Calculate the square root of the vector 
void vector_square_root(gsl_vector *vector) {
    if (vector == NULL) {
        printf("Null pointer passed to vector_square_root function.\n");
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < vector->size; i++) {
        double value = gsl_vector_get(vector, i); 
        if (value < 0) value = 0.0; // Prevent NaN
        gsl_vector_set(vector, i, sqrt(value));
    }
}
