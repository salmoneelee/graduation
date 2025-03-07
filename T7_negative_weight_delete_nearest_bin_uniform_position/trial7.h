#ifndef trial7_h
#define trial7_h

#include <gsl/gsl_rng.h>

#define MAX_POINTERS 30 // Limit number of dynamically allocated pointers

extern int N; // Total number of neutrons
extern int boundary_condition; // Boundary condition (0 = vacuum, 1 = reflective)

extern const double width; // Width of the reactor
extern const int num_bins; // Number of bins in the reactor

extern const double sigma_s; // Scattering cross section
extern const double sigma_c; // Capture cross section
extern const double sigma_f; // Fission cross section
extern const double nu; // Average number of neutrons produced per fission

extern int parent_neutron; // Number of neutrons in the parent fission bank
extern int child_neutron; // Number of neutrons in the child fission bank
extern int *parent_pointer; // Pointer to parent_neutron
extern int *child_pointer; // Pointer to child_neutron

extern const gsl_rng_type *T; // Type of GSL RNG 
extern gsl_rng *r; // Pointer to GSL RNG instance

// Struct to store dynamically allocated pointer and type
typedef struct {
    void *data;
    int type;
} cleanup_entry;

extern cleanup_entry chained_pointers[MAX_POINTERS]; // Array of pointers with metadata that need to be freed of closed
extern int chain_count; // Number of pointers that need to be freed or closed

typedef struct {
    double position; // Position of neutron 
    double weight; // Weight of neutron 
} neutron;

// Define types for better readability
typedef enum {
    TYPE_FILE,
    TYPE_CALLOC,
    TYPE_GSL_VECTOR,
    TYPE_RNG
} pointertype;

void register_pointer(void *ptr, int type);
void unregister_pointer(void *ptr);
void cleanup(void);

void initialize_rng(void);
double random_number_generator(void);
void free_rng(void);

void initialize_fission_bank(neutron **, int *);

gsl_vector *simulate_neutron_diffusion(FILE *,neutron *, neutron **, double *, double);
double determine_direction_and_distance(void);
void fission(FILE *, neutron *, int *, double);
void add_fission_bank(double, int, neutron **, double);
void k_tally(double *, neutron *);
void flux_tally(gsl_vector *, neutron *, double);
void force_scattering(neutron *);
void russian_roulette(neutron *);

void remove_and_redistribute_neutrons(FILE *, neutron **);

void fission_bank_size_adjustment(FILE *, neutron **);

void position_tally(gsl_vector *, neutron *, int);
void sum_of_weights(neutron *, double *, int);
void vector_square_root(gsl_vector *);

#endif