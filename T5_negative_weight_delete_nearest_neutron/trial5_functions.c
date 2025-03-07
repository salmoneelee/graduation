#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include "trial5.h"
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
    gsl_rng_set(r, time(NULL)); // Set seed to current time
    //gsl_rng_set(r, 10); // Set a fixed seed
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
        (*fission_bank_double_pointer)[*index_pointer].weight = 1.0;
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
        gsl_vector_set(vector, bin_number, current_value + fission_bank_pointer[i].weight);
    }
}

// Sum of weights of neutrons in a fission bank
void sum_of_weights(neutron *fission_bank, double *total_weight_pointer, int neutron) {
    for (int i = 1; i < neutron + 1; i++) {
        (*total_weight_pointer) += fission_bank[i].weight;
    }
}

// Simulate neutron diffusion for one cycle
gsl_vector *simulate_neutron_diffusion(FILE *file2, neutron *parent_fission_bank, neutron **child_fission_bank_double_pointer, double *k_pointer, double k) {
    gsl_vector *flux = gsl_vector_calloc(num_bins); 
    if (flux == NULL) {
        printf("Memory allocation for flux failed.\n");
        exit(EXIT_FAILURE);
    }

    while (parent_fission_bank[parent_neutron].position != -1.0) {
        double initial_position = parent_fission_bank[parent_neutron].position;
        double distance = determine_direction_and_distance(); // Distance to the next collision
        double new_position = parent_fission_bank[parent_neutron].position + distance;  // New position for the neutron
        parent_fission_bank[parent_neutron].position = new_position;

        // Debug Start: Logging simulation step
        fprintf(file2, "---- simulating neutron# %d ----\n ", parent_neutron);
        fprintf(file2, "weight: %f, pos: %f, dis: %f, ", parent_fission_bank[parent_neutron].weight, initial_position, distance);
        // Debug End

        int fission_neutrons = 0; // Number of neutrons produced by fission
        double child_weight = 1.0; // Weight of child neutron
        if (new_position < width && new_position > 0) {
            fprintf(file2, "newpos: %f\n", new_position);

            fission(file2, parent_fission_bank, &fission_neutrons, k); // Number of neutrons produced by fission
            if (parent_fission_bank[parent_neutron].weight < 0) {
                child_weight = -1.0; 
            }
            add_fission_bank(child_weight, fission_neutrons, child_fission_bank_double_pointer, new_position);
            
            k_tally(k_pointer, parent_fission_bank);
            flux_tally(flux, parent_fission_bank, new_position);
        
            force_scattering(parent_fission_bank);
            russian_roulette(parent_fission_bank);
        }
            
        else {
            if (new_position > 0) {
                if (((int)(new_position/width))%2 != 0) {
                    new_position = width - fmod(new_position, width);
                }
                else {
                    new_position = fmod(new_position, width);
                }
            }
            else {
                if (((int)(new_position/width))%2 != 0) {
                    new_position = width + fmod(new_position, width);
                }
                else {
                    new_position = -fmod(new_position, width);
                }
            }
    
                parent_fission_bank[parent_neutron].position = new_position;
                fprintf(file2, "newpos: %f\n", new_position);

                if (boundary_condition == 0) { // Vacuum
                    parent_fission_bank[parent_neutron].weight *= -1; // Negative weight neutron
                    if (parent_fission_bank[parent_neutron].weight < 0) {
                        child_weight = -1.0;
                    }
                }
    
                fission(file2, parent_fission_bank, &fission_neutrons, k); 
                add_fission_bank(child_weight, fission_neutrons, child_fission_bank_double_pointer, new_position);
                
                k_tally(k_pointer, parent_fission_bank);
                flux_tally(flux, parent_fission_bank, new_position);
            
                force_scattering(parent_fission_bank);
                russian_roulette(parent_fission_bank);
        }
    }
    return flux;
} 

// Determine direction and distance of neutron
double determine_direction_and_distance(void) {
    double random_number_1 = random_number_generator();
    if (random_number_1 == 0) {
        random_number_1 = 1e-10; // Small positive number to prevent log(0)
    }
    double distance_3D = (-log(random_number_1))/(sigma_s + sigma_c + sigma_f); // 3D distance to the next collision
    
    double random_number_2 = random_number_generator();
    double cos_theta = 2 * random_number_2 - 1; // Cosine of the angle between neutron's path and x-axis, in the range [-1, 1)

    return distance_3D * cos_theta;
}

// Simulate fission
void fission(FILE *file2, neutron *parent_fission_bank, int *fission_neutrons, double k) {
    double init_weight = parent_fission_bank[parent_neutron].weight; // Weight of the parent neutron before collision
    double abs_init_weight = fabs(init_weight); // Absolute value of init_weight
    double R = abs_init_weight * nu * (sigma_f/(sigma_s + sigma_c + sigma_f)) * (1/k); 
    int n = (int)floor(R); 
    double random_number = random_number_generator(); 

    if (random_number < R - n) { // n+1 new neutrons are produced
        (*fission_neutrons) = n+1;

        for (int i = 0; i < (*fission_neutrons); i++) {
            fprintf(file2, "FISSION!!\n");
        }
    }
    else { // n new neutrons are produced
        (*fission_neutrons) = n;

        for (int i = 0; i < (*fission_neutrons); i++) {
            fprintf(file2, "FISSION!!\n");
        }
    }

    fprintf(file2, "init_weight: %f, abs_init_weight: %f, k: %f R: %f, n: %d, rand: %f, fission_neutrons: %d\n\n", init_weight, abs_init_weight, k, R, n, random_number, (*fission_neutrons));
}

// Add fission-produced neutrons to the fission bank
void add_fission_bank(double child_weight, int n, neutron **fission_bank_double_pointer, double position) {
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
        (*fission_bank_double_pointer)[*child_pointer].weight = child_weight;
    }
}

// Update the k value for collision
void k_tally(double *k_pointer, neutron *parent_fission_bank) {
    double init_weight = parent_fission_bank[parent_neutron].weight; // Weight of the parent neutron before collision
    
    *k_pointer += init_weight * nu * (sigma_f/(sigma_s + sigma_c + sigma_f));
}

// Update the neutron flux for collision
void flux_tally(gsl_vector *flux, neutron *fission_bank_pointer, double position) {
    int bin_number = (int)floor(position/(width/num_bins)); // Bin number of the position
    if (bin_number < 0) {
        printf("Neutron is out of bin range.\n");
        bin_number = 0;
    }
    if (bin_number == num_bins) {
        bin_number = num_bins - 1; 
    }

    double current_value = gsl_vector_get(flux, bin_number); // Current value of the flux vector at bin_number element
    gsl_vector_set(flux, bin_number, current_value + (1/(sigma_s + sigma_c + sigma_f)) * fission_bank_pointer[parent_neutron].weight); // Set element at index bin_number of the flux vector to current_value +1/(sigma_s + sigma_c + sigma_f)
}

// Force scattering
void force_scattering(neutron *parent_fission_bank) {
    parent_fission_bank[parent_neutron].weight *= (sigma_s/(sigma_s + sigma_c + sigma_f));
}

// Adjust the weight of the neutron
void russian_roulette(neutron *parent_fission_bank) {
    double fin_weight = parent_fission_bank[parent_neutron].weight; // weight of the parent neutron after collision

    if (fin_weight < 0.1 && fin_weight > -0.1) {
        double random_number = random_number_generator();
        if (random_number < fin_weight) {
            if (fin_weight > 0) {
                parent_fission_bank[parent_neutron].weight = 1.0;
            }
            else {
                parent_fission_bank[parent_neutron].weight = -1.0;
            }
        }
        else {
            parent_neutron--;
        }
    }
}


// Quick sort the child fission bank
void quicksort(neutron *fission_bank, int start_index, int end_index) {
    if (start_index < end_index) {
        int pivot = random_partition(fission_bank, start_index, end_index);
        quicksort(fission_bank, start_index, pivot-1);
        quicksort(fission_bank, pivot+1, end_index);
    }
}

// In-place partitioning
int partition(neutron *fission_bank, int start_index, int end_index) {
    double pivot_value = fission_bank[start_index].position;
    int left = start_index + 1;
    int right = end_index;

    while (left <= right) {
        while (left <= end_index && fission_bank[left].position <= pivot_value) {
            left++;
        }
        while (right > start_index && fission_bank[right].position > pivot_value) {
            right--;
        }
        if (left < right) {
            swap(&fission_bank[left], &fission_bank[right]);
        }
    }

    swap(&fission_bank[start_index], &fission_bank[right]);

    return right;
}

int random_partition(neutron *fission_bank, int start_index, int end_index) {
    int random_pivot = start_index + rand() % (end_index - start_index+ 1);
    swap(&fission_bank[random_pivot], &fission_bank[start_index]);
    return partition(fission_bank, start_index, end_index);
}

// Swap two elements
void swap(neutron *element1, neutron *element2) {
    neutron temp = *element1;
    *element1 = *element2;
    *element2 = temp;
}

// Remove positive weight neutrons using negative weight neutrons within the child fission bank
void remove_negative_weight(neutron **child_fission_bank_double_pointer) {
    neutron *new_fission_bank = calloc(child_neutron + 5, sizeof(neutron)); // New child_fission_bank
    if (new_fission_bank == NULL) {
        printf("Memory allocation for size adjustment failed.\n");
        exit(EXIT_FAILURE);
    }

    new_fission_bank[0].position = -1.0;
    new_fission_bank[0].weight = 0.0;

    int negative_neutrons = 0; // number of negative weight neutrons
    int index = 1; // index of new_fission_bank

    for (int i = 1; i < child_neutron + 1; i++) {
        if ((*child_fission_bank_double_pointer)[i].weight == 10) {
            continue;
        }

        if ((*child_fission_bank_double_pointer)[i].weight < 0) {
            int left_positive_index = 0; // Index of the nearest neutron with a positive weight that is to the left of the negative weight neutron
            int right_positive_index = 0; // Index of the nearest neutron with a positive weight that is to the right of the negative weight neutron
            
            for (int j = i-1; j > 0; j--) {
                if ((*child_fission_bank_double_pointer)[j].weight > 0 && (*child_fission_bank_double_pointer)[j].weight != 10.0) {
                    left_positive_index = j;
                    break;
                }
            }

            for (int j = i+1; j < child_neutron + 1; j++) {
                if ((*child_fission_bank_double_pointer)[j].weight > 0 && (*child_fission_bank_double_pointer)[j].weight != 10.0) {
                    right_positive_index = j;
                    break;
                }
            }

            double diff_left = fabs((*child_fission_bank_double_pointer)[i].position - (*child_fission_bank_double_pointer)[left_positive_index].position);
            double diff_right = fabs((*child_fission_bank_double_pointer)[right_positive_index].position - (*child_fission_bank_double_pointer)[i].position);

            if (diff_left < diff_right) {
                if (left_positive_index == 0) {
                    (*child_fission_bank_double_pointer)[right_positive_index].weight = 10.0; // Remove the positive weight neutron by changing the weight to 10.0
                }
                else {
                      (*child_fission_bank_double_pointer)[left_positive_index].weight = 10.0; // Remove the positive weight neutron by changing the weight to 10.0
                    index--;
                    new_fission_bank[index].weight = 0.0; // Remove the most recent neutron by changing the weight to 0.0
                    new_fission_bank[index].position = 0.0; // Remove the most recent neutron by changing the weight to 0.0
                }
            }
            else {
                (*child_fission_bank_double_pointer)[right_positive_index].weight = 10.0; // Remove the positive weight neutron by changing the weight to 10.0
            }

            (*child_fission_bank_double_pointer)[i].weight = 10.0; // Remove the negative weight neutron by changing the weight to 10.0
        }

        else if ((*child_fission_bank_double_pointer)[i].weight > 0) {
            new_fission_bank[index] = (*child_fission_bank_double_pointer)[i];
            index++;
        }
    }

    while (negative_neutrons > 0) {
        index--;
        new_fission_bank[index].position = 0.0;
        new_fission_bank[index].weight = 0.0;
        negative_neutrons--;
    }

    child_neutron = index - 1; // Update the number of child neutrons

    for (int i = 0; i < chain_count; i++) {
        if (chained_pointers[i].data == *child_fission_bank_double_pointer) {
            chained_pointers[i].data = new_fission_bank; 
            break;
        }
    }

    free((*child_fission_bank_double_pointer));
    (*child_fission_bank_double_pointer) = new_fission_bank;
}

// Adjust the size of the child_fission_bank to N
void fission_bank_size_adjustment(FILE *file3, neutron **child_fission_bank_double_pointer) {
    neutron *new_fission_bank = calloc(2 * N, sizeof(neutron)); // for the size adjusted child_fission_bank
    if (new_fission_bank == NULL) {
        printf("Memory allocation for size adjustment failed.\n");
        exit(EXIT_FAILURE);
    }

    // represents the base of the fission bank
    new_fission_bank[0].position = -1.0;
    new_fission_bank[0].weight = 0.0;

    for (int i = 1; i < N+1; i++)
    {   
        double random_number = random_number_generator();
        int index = ((int)(random_number * (child_neutron))); //  randomly chosen index = [0, child_neutron)

        ////----- debug -----////
        fprintf(file3, "child_neutron: %d, random: %f, index: %d\n", child_neutron, random_number, index+1); 
        ////-------------------////

        new_fission_bank[i].position = (*child_fission_bank_double_pointer)[index+1].position; // index = [1, child_neutron]
        new_fission_bank[i].weight = (*child_fission_bank_double_pointer)[index+1].weight;
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

// Reweight child neutrons to conserve the total weight of each cycles 
void reweight_neutrons(FILE *file4, neutron **child_fission_bank_double_pointer) {
    if (file4 == NULL || child_fission_bank_double_pointer == NULL || *child_fission_bank_double_pointer == NULL) {
        printf("Null pointer passed to reweight_neutrons function.\n");
        exit(EXIT_FAILURE);
    }

    double total_weight = 0.0;
    sum_of_weights((*child_fission_bank_double_pointer), &total_weight, child_neutron);

    double renormalize = (double)N/child_neutron;

    fprintf(file4, "----------- child before reweight -----------\n");
    fprintf(file4, "total_weight: %f\n", total_weight); 
    for (int i = 1; i < child_neutron + 1; i++) {
        fprintf(file4, "%d %0.10f %0.10f\n", i, (*child_fission_bank_double_pointer)[i].position, (*child_fission_bank_double_pointer)[i].weight); 
    }
    fprintf(file4, "\n");

   for (int i = 1; i < child_neutron + 1; i++) {
        (*child_fission_bank_double_pointer)[i].weight = renormalize;
    }

    fprintf(file4, "----------- child after reweight -----------\n");
    fprintf(file4, "total_weight: %f\n", total_weight); 
    for (int i = 1; i < child_neutron + 1; i++) {
        fprintf(file4, "%d %0.10f %0.10f\n", i, (*child_fission_bank_double_pointer)[i].position, (*child_fission_bank_double_pointer)[i].weight); 
    }
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