/* Standalone provider-injection and C boundary smoke test. No physical tables.
 * build_connected_native.py --output /tmp/fixture.dylib --check-factory --sanitize
 */
#include "connected_dp.hpp"
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

using aptaswitch::thermodynamics::EnergyProvider;
using aptaswitch::thermodynamics::FaceMode;
using aptaswitch::thermodynamics::make_handle_v1;

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
class CountingProvider final : public EnergyProvider {
    int& calls_;
    bool invalid_;
public:
    explicit CountingProvider(int& calls, bool invalid = false) : calls_(calls), invalid_(invalid) {}
    double log_hairpin(int, int) const override { ++calls_; return std::log(7.0); }
    double log_interior(int, int, int, int) const override { ++calls_; return 0; }
    double log_unpaired(FaceMode, int, int) const override {
        ++calls_; return invalid_ ? std::numeric_limits<double>::quiet_NaN() : 0;
    }
    double log_branch(FaceMode, int, int, int) const override { ++calls_; return std::log(3.0); }
    double log_closing(FaceMode, int, int, int) const override { ++calls_; return std::log(2.0); }
};
class SmallPenaltyProvider final : public EnergyProvider {
public:
    double log_hairpin(int, int) const override { return -.06; }
    double log_interior(int, int, int, int) const override { return 0; }
    double log_unpaired(FaceMode, int, int) const override { return 0; }
    double log_branch(FaceMode, int, int, int) const override { return 0; }
    double log_closing(FaceMode, int, int, int) const override { return 0; }
};
void success(int status, const char* error) {
    if (status != APTADP_OK) throw std::runtime_error(error);
}
}

int main() {
    try {
        char error[512]{};
        const uint8_t material[2]{2, 2}, cuts[3]{0, 1, 0}, allowed[4]{0, 1, 1, 0};
        aptadp_input_v1 input{sizeof(aptadp_input_v1), 2, "AA", material, cuts, allowed, nullptr, 3, 30, 0};
        int calls = 0;
        // The only connected structure is (). Its two terminal contributions
        // are 2 and 3; the distinct global factor is applied exactly once.
        std::unique_ptr<aptadp_handle, decltype(&aptadp_destroy)> handle(
            make_handle_v1(input, .6, std::log(5.0), std::make_unique<CountingProvider>(calls)), aptadp_destroy);
        double log_z = 0, matrix[4]{}, adjustment = -1, minimum = 0;
        char structure[3]{}, samples[3 * 17]{};
        success(aptadp_log_partition(handle.get(), &log_z, error, sizeof(error)), error);
        require(std::abs(log_z - std::log(30.0)) < 1e-14, "Injected provider partition mismatch");
        success(aptadp_pair_probabilities(handle.get(), matrix, 4, &adjustment, error, sizeof(error)), error);
        require(matrix[0] == 0 && matrix[1] == 1 && matrix[2] == 1 && matrix[3] == 0,
                "Injected provider probabilities mismatch");
        success(aptadp_minimum(handle.get(), &minimum, structure, sizeof(structure), error, sizeof(error)), error);
        require(std::abs(minimum + .6 * std::log(30.0)) < 1e-14 && std::string(structure) == "()",
                "Injected provider minimum mismatch");
        double pair_minimum = 0;
        success(aptadp_minimum_pair_energies(handle.get(), matrix, 4, &pair_minimum, error, sizeof(error)), error);
        require(pair_minimum == minimum && matrix[1] == minimum && matrix[2] == minimum
                && matrix[0] == std::numeric_limits<double>::infinity()
                && matrix[3] == std::numeric_limits<double>::infinity(), "Injected provider pair minima mismatch");
        success(aptadp_sample(handle.get(), 19, 17, samples, sizeof(samples), error, sizeof(error)), error);
        for (int i = 0; i < 17; ++i) require(std::string(samples + 3 * i) == "()", "Injected traceback mismatch");
        require(calls > 0, "Native provider was not queried");

        // Validation must precede any energy callback, including construction.
        calls = 0;
        auto invalid_input = input;
        invalid_input.struct_size = 0;
        bool rejected = false;
        try {
            aptadp_destroy(make_handle_v1(invalid_input, .6, 0, std::make_unique<CountingProvider>(calls)));
        } catch (const std::invalid_argument&) { rejected = true; }
        require(rejected && calls == 0, "Invalid input reached the provider");
        rejected = false;
        try { aptadp_destroy(make_handle_v1(input, .6, 0, nullptr)); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "NULL native provider accepted");
        rejected = false;
        try { aptadp_destroy(make_handle_v1(input, .6, 0, std::make_unique<CountingProvider>(calls, true))); }
        catch (const std::range_error&) { rejected = true; }
        require(rejected, "Nonfinite provider weight accepted");

        // All errors below are returned through the C ABI without exceptions.
        require(aptadp_sample(handle.get(), 0, 17, samples, sizeof(samples) - 1, error, sizeof(error))
                == APTADP_BAD_INPUT, "Short sample buffer accepted");
        require(aptadp_pair_probabilities(handle.get(), matrix, 3, &adjustment, error, sizeof(error))
                == APTADP_BAD_INPUT, "Short probability buffer accepted");
        require(aptadp_minimum_pair_energies(handle.get(), matrix, 3, &pair_minimum, error, sizeof(error))
                == APTADP_BAD_INPUT, "Short pair-energy buffer accepted");
        require(aptadp_log_partition(nullptr, &log_z, error, sizeof(error)) == APTADP_BAD_INPUT,
                "NULL handle accepted");
        // Two independent .06-kcal branches: accepting local errors of .1
        // would incorrectly admit the double-branch .12-kcal structure.
        const uint8_t chemistry4[4]{2, 2, 2, 2}, cuts4[5]{0, 0, 0, 0, 0};
        const uint8_t allowed4[16]{0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0};
        const aptadp_input_v1 input4{sizeof(aptadp_input_v1), 4, "AAAA", chemistry4, cuts4, allowed4,
                                    nullptr, 0, 30, 0};
        std::unique_ptr<aptadp_handle, decltype(&aptadp_destroy)> tolerance_handle(
            make_handle_v1(input4, 1, 123, std::make_unique<SmallPenaltyProvider>()), aptadp_destroy);
        char alternatives[50]{}; double alternative_energies[10]{};
        uint64_t count = 0; uint8_t truncated = 0;
        success(aptadp_minimum_structures(tolerance_handle.get(), .1, 10, alternatives, sizeof(alternatives),
            alternative_energies, 10, &count, &truncated, &minimum, error, sizeof(error)), error);
        require(count == 3 && truncated == 0 && minimum == -123 && std::string(alternatives) == "....",
                "Global minimum enumeration window or first structure mismatch");
        for (uint64_t i = 0; i < count; ++i) {
            require(std::string(alternatives + 5 * i) != "()()", "Local tolerances accumulated across branches");
            require(alternative_energies[i] <= minimum + .1, "Global energy window was exceeded");
        }
        success(aptadp_minimum_structures(tolerance_handle.get(), .1, 3, alternatives, sizeof(alternatives),
            alternative_energies, 10, &count, &truncated, &minimum, error, sizeof(error)), error);
        require(count == 3 && !truncated, "Exact output limit incorrectly marked truncated");
        success(aptadp_minimum_structures(tolerance_handle.get(), .1, 1, alternatives, sizeof(alternatives),
            alternative_energies, 10, &count, &truncated, &minimum, error, sizeof(error)), error);
        require(count == 1 && truncated, "Omitted minimum alternative was not reported");
        success(aptadp_minimum_structures(tolerance_handle.get(), .1, 0, nullptr, 0, nullptr, 0,
            &count, &truncated, &minimum, error, sizeof(error)), error);
        require(count == 0 && truncated, "Zero-limit minimum query mismatch");
        std::cout << "Provider injection and C boundary checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
