#ifndef APTASWITCH_CONNECTED_DP_HPP
#define APTASWITCH_CONNECTED_DP_HPP

/* Standalone C ABI. No external thermodynamic dependency is linked by the fixture
 * build. Buffers remain owned by the caller; create copies all input arrays.
 * Indices and lengths are nucleotide coordinates, preserving input order.
 * Calls on one handle are serialized internally. The caller must not destroy
 * a handle concurrently with its use; caller-owned buffers must not overlap. */
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define APTADP_ABI_VERSION 1u
enum aptadp_status {
    APTADP_OK = 0, APTADP_BAD_INPUT = 1, APTADP_NUMERICAL_ERROR = 2,
    APTADP_EMPTY_ENSEMBLE = 3, APTADP_INTERNAL_ERROR = 4, APTADP_ALLOCATION_ERROR = 5
};
enum aptadp_fixture_provider { APTADP_FIXTURE_UNIT = 1, APTADP_FIXTURE_ASYMMETRIC = 2 };

typedef struct aptadp_input_v1 {
    uint32_t struct_size;
    uint32_t length;
    const char *sequence;               /* exactly length ASCII bases */
    const uint8_t *materials;           /* length entries: DNA=1, RNA=2 */
    const uint8_t *nicks_before;        /* length+1 entries; endpoints zero */
    const uint8_t *pair_allowed;        /* length^2, symmetric binary, diagonal zero */
    const int32_t *constraints;         /* length entries, or NULL for unrestricted */
    int32_t min_hairpin_size;
    int32_t max_interior_size;          /* 0..30; explicit ensemble convention */
    uint64_t memory_limit_bytes;        /* zero selects 512 MiB; no silent truncation */
} aptadp_input_v1;

typedef struct aptadp_model_v1 {
    uint32_t struct_size;
    uint32_t fixture_provider;          /* synthetic fixtures only, never physical */
    double thermal_energy_kcal_mol;
    double log_global_factor;          /* association/symmetry once; -inf permits zero */
} aptadp_model_v1;

typedef struct aptadp_handle aptadp_handle;

uint32_t aptadp_abi_version(void);
int aptadp_create_v1(const aptadp_input_v1 *, const aptadp_model_v1 *, aptadp_handle **,
                     char *error, size_t error_capacity);
void aptadp_destroy(aptadp_handle *);
int aptadp_log_partition(aptadp_handle *, double *log_z, char *error, size_t error_capacity);
int aptadp_pair_probabilities(aptadp_handle *, double *matrix, size_t matrix_cells,
                              double *maximum_endpoint_adjustment, char *error, size_t error_capacity);
int aptadp_minimum(aptadp_handle *, double *energy_kcal_mol, char *structure,
                   size_t structure_capacity, char *error, size_t error_capacity);
/* Minimum total labelled-state energy subject to containing each pair.
 * Matrix is symmetric; diagonal and impossible pairs are +infinity. The
 * distinct output minimum is the unrestricted minimum within input constraints.
 * Global association/uniform cyclic-symmetry energy is included once in every
 * finite value. These are not structure-dependent physical symmetry minima.
 * An empty connected ensemble returns APTADP_EMPTY_ENSEMBLE. */
int aptadp_minimum_pair_energies(aptadp_handle *, double *matrix, size_t matrix_cells,
                                double *minimum_total_energy, char *error, size_t error_capacity);
/* Deterministic enumeration in the global energy window [minimum,minimum+atol].
 * First output equals aptadp_minimum's traceback. Other outputs are not sorted.
 * Structures have stride length+1 including NUL; energies are their actual
 * total labelled-state energies. truncated=1 proves another admissible structure
 * exists beyond count_limit; no completeness claim is made then. count_limit=0
 * is allowed and returns count=0,truncated=1 for a nonempty ensemble.
 * The hybrid window is tested before adding the shared global constant, avoiding
 * cancellation; floating-point energies are never rounded or clamped. */
int aptadp_minimum_structures(aptadp_handle *, double atol_kcal, uint64_t count_limit,
                              char *structures, size_t structure_bytes, double *energies, size_t energy_cells,
                              uint64_t *count, uint8_t *truncated, double *minimum_total_energy,
                              char *error, size_t error_capacity);
/* Each sample is length dot-bracket bytes followed by NUL, with stride length+1.
 * Fixture/hybrid factories use mt19937_64. Other explicit factories document
 * their own stream (the official RNA factory uses RNAlib's generator). */
int aptadp_sample(aptadp_handle *, uint64_t seed, uint64_t count, char *structures,
                  size_t structures_capacity, char *error, size_t error_capacity);
/* Peak DP workspace estimate (cached PF inside + max inside + max outside +
 * output matrix and auxiliaries), excluding provider tables and process RSS.
 * Hybrid enumeration frontier allocations are also checked against the input
 * limit. See each alternative factory's contract for algorithm-specific bounds. */
int aptadp_memory_estimate(aptadp_handle *, uint64_t *peak_bytes, char *error, size_t error_capacity);

#ifdef __cplusplus
}
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace aptaswitch { namespace thermodynamics {
enum class FaceMode : uint8_t { Multiloop, Exterior };

/* Internal dispatch only: the fixture build remains independent of RNAlib.
 * An explicit factory identifies which algorithm implements each handle. */
struct EmptyEnsemble : std::runtime_error { using std::runtime_error::runtime_error; };
class Calculation {
public:
    virtual ~Calculation() = default;
    virtual int length() const = 0;
    virtual uint64_t memory() const = 0;
    virtual double log_partition() = 0;
    virtual void pairs(double*, double&) = 0;
    virtual std::pair<double, std::string> minimum() = 0;
    virtual double minimum_pair_energies(double*) const = 0;
    virtual double minimum_structures(double, uint64_t, char*, double*, uint64_t&, uint8_t&) const = 0;
    virtual void sample(uint64_t, uint64_t, char*) = 0;
};
aptadp_handle *make_calculation_handle_v1(std::unique_ptr<Calculation> calculation);

/* Native extension point for physical providers. Values are natural-log
 * Boltzmann weights, using one immutable parameter snapshot and one RT.
 * The provider owns sequence/material/condition context as needed. No Python
 * callbacks occur inside a recurrence. Branch i<j; closing includes its own
 * correctly oriented terminal penalty, plus multiloop initiation if relevant.
 * A family is the completed face's material mask, not the complex's mask. */
class EnergyProvider {
public:
    virtual ~EnergyProvider() = default;
    virtual double log_hairpin(int i, int j) const = 0;
    virtual double log_interior(int i, int j, int k, int l) const = 0;
    virtual double log_unpaired(FaceMode mode, int family, int i) const = 0;
    virtual double log_branch(FaceMode mode, int family, int i, int j) const = 0;
    virtual double log_closing(FaceMode mode, int family, int i, int j) const = 0;
};
/* Native-only injection point. Input arrays are validated/copied before the
 * provider is queried. Throws C++ exceptions on failure; a C entrypoint must
 * catch and translate them. On success the caller owns the returned handle
 * and destroys it with aptadp_destroy. A NULL provider is rejected. */
aptadp_handle *make_handle_v1(const aptadp_input_v1& input,
                            double thermal_energy_kcal_mol, double log_global_factor,
                            std::unique_ptr<EnergyProvider> provider);
}}
#endif
#endif
