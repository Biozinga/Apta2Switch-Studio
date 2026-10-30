#ifndef APTASWITCH_HYBRID_API_HPP
#define APTASWITCH_HYBRID_API_HPP
#include "connected_dp.hpp"
#ifdef __cplusplus
extern "C" {
#endif
typedef struct aptadp_hybrid_model_v1 {
    uint32_t struct_size;
    uint32_t salt_policy; /* 1: Na=1 M reference; 2: published pair correction */
    double celsius;
    double sodium_molar;
    double magnesium_molar; /* must be zero */
} aptadp_hybrid_model_v1;
/* This entrypoint enforces the complete chemical pair mask and grammar of the
 * identified hybrid model; conditional constraints remain explicit. */
int aptadp_create_hybrid_v1(const aptadp_input_v1 *, const aptadp_hybrid_model_v1 *,
                          aptadp_handle **, double *thermal_energy,
                          char *error, size_t error_capacity);
/* Official ViennaRNA single-RNA PF/BPP/sampling and rounded-energy MFE.
 * Same chemical parameters and exact PF salt correction as the hybrid model.
 * Rejects multistrand/DNA input. Pair-conditioned minima use repeated official
 * constrained folds (O(n^5) worst-case). Minimum enumeration bounds retained
 * output; RNAlib still traverses the complete energy window and its internal
 * enumeration frontier is not covered by aptadp_memory_estimate.
 * MFE salt rounds to nearest centikcal (ties away from zero); PF salt does not.
 * Sampling uses RNAlib's stream seeded with uint32(seed ^ (seed >> 32)). */
int aptadp_create_rna_v1(const aptadp_input_v1 *, const aptadp_hybrid_model_v1 *,
                       aptadp_handle **, double *thermal_energy,
                       char *error, size_t error_capacity);
#ifdef __cplusplus
}
#endif
#endif
