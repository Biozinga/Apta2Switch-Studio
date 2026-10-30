#include "hybrid_api.hpp"
#include "hybrid_provider.hpp"
#include <algorithm>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <new>

extern "C" int aptadp_create_hybrid_v1(const aptadp_input_v1 *input,
        const aptadp_hybrid_model_v1 *model, aptadp_handle **output, double *thermal,
        char *error, size_t error_capacity) {
    using namespace aptaswitch::thermodynamics;
    if (output) *output = nullptr;
    auto fail = [&](int status, const char *message) {
        if (error && error_capacity) std::snprintf(error, error_capacity, "%s", message);
        return status;
    };
    try {
        if (!input || !model || !output || !thermal
            || input->struct_size != sizeof(aptadp_input_v1)
            || model->struct_size != sizeof(aptadp_hybrid_model_v1))
            throw std::invalid_argument("Invalid hybrid API input or ABI structure size");
        if (!input->length || input->length > 100000 || !input->sequence
            || !input->materials || !input->nicks_before || !input->pair_allowed)
            throw std::invalid_argument("Invalid hybrid input buffers");
        if (input->min_hairpin_size != 3 || input->max_interior_size != 30)
            throw std::invalid_argument("Hybrid model requires min_hairpin=3, max_interior=30");
        if (model->magnesium_molar != 0 || (model->salt_policy != 1 && model->salt_policy != 2))
            throw std::invalid_argument("Hybrid model requires Mg=0 and an explicit salt policy");
        const size_t n = input->length;
        const uint64_t limit = input->memory_limit_bytes ? input->memory_limit_bytes : 512ULL*1024*1024;
        if (n*n > limit / sizeof(double)) throw std::bad_alloc();
        auto provider = std::make_unique<HybridEnergyProvider>(
            std::string(input->sequence, n), std::vector<uint8_t>(input->materials,input->materials+n),
            std::vector<uint8_t>(input->nicks_before,input->nicks_before+n+1), model->celsius,
            model->sodium_molar, model->salt_policy == 1 ? HybridEnergyProvider::SaltCorrection::Reference
                                                      : HybridEnergyProvider::SaltCorrection::PublishedPair);
        for (size_t i=0; i<n; ++i)
            for (size_t j=0; j<n; ++j) {
                const uint8_t expected = i!=j && provider->pair_allowed(static_cast<int>(i),static_cast<int>(j));
                if (input->pair_allowed[i*n+j] != expected)
                    throw std::invalid_argument("Pair mask differs from the hybrid model chemistry");
            }
        const double rt=provider->thermal_energy(), global=provider->log_global_factor();
        *output=make_handle_v1(*input,rt,global,std::move(provider));
        *thermal=rt;
        if (error && error_capacity) error[0]='\0';
        return APTADP_OK;
    } catch (const std::bad_alloc& e) { return fail(APTADP_ALLOCATION_ERROR,e.what()); }
      catch (const std::invalid_argument& e) { return fail(APTADP_BAD_INPUT,e.what()); }
      catch (const std::range_error& e) { return fail(APTADP_NUMERICAL_ERROR,e.what()); }
      catch (const std::exception& e) { return fail(APTADP_INTERNAL_ERROR,e.what()); }
      catch (...) { return fail(APTADP_INTERNAL_ERROR,"Unknown native hybrid error"); }
}

extern "C" int aptadp_create_rna_v1(const aptadp_input_v1 *input,
        const aptadp_hybrid_model_v1 *model, aptadp_handle **output, double *thermal,
        char *error, size_t error_capacity) {
    using namespace aptaswitch::thermodynamics;
    if (output) *output = nullptr;
    auto fail = [&](int status, const char *message) {
        if (error && error_capacity) std::snprintf(error, error_capacity, "%s", message);
        return status;
    };
    try {
        if (!input || !model || !output || !thermal
            || input->struct_size != sizeof(aptadp_input_v1)
            || model->struct_size != sizeof(aptadp_hybrid_model_v1))
            throw std::invalid_argument("Invalid RNA API input or ABI structure size");
        if (model->magnesium_molar != 0 || (model->salt_policy != 1 && model->salt_policy != 2))
            throw std::invalid_argument("Single-RNA model requires Mg=0 and an explicit salt policy");
        double calculated_thermal = 0.;
        auto calculation = make_vienna_rna_calculation(*input, model->celsius, model->sodium_molar,
            model->salt_policy == 1 ? HybridEnergyProvider::SaltCorrection::Reference
                                    : HybridEnergyProvider::SaltCorrection::PublishedPair,
            calculated_thermal);
        *output = make_calculation_handle_v1(std::move(calculation));
        *thermal = calculated_thermal;
        if (error && error_capacity) error[0] = '\0';
        return APTADP_OK;
    } catch (const std::bad_alloc& e) { return fail(APTADP_ALLOCATION_ERROR, e.what()); }
      catch (const std::invalid_argument& e) { return fail(APTADP_BAD_INPUT, e.what()); }
      catch (const std::range_error& e) { return fail(APTADP_NUMERICAL_ERROR, e.what()); }
      catch (const std::exception& e) { return fail(APTADP_INTERNAL_ERROR, e.what()); }
      catch (...) { return fail(APTADP_INTERNAL_ERROR, "Unknown native RNA error"); }
}
