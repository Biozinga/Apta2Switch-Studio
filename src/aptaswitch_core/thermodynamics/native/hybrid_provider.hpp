#ifndef APTASWITCH_HYBRID_PROVIDER_HPP
#define APTASWITCH_HYBRID_PROVIDER_HPP

#include "connected_dp.hpp"
#include <memory>
#include <string>
#include <vector>

namespace aptaswitch { namespace thermodynamics {

/* Native counterpart of HybridMotifModel. This is the identified experimental
 * Sugimoto95 / Vienna04 loop-mean model, not an empirical mixed-loop fit.
 * Snapshot construction is serialized; evaluating a snapshot changes no global
 * Vienna parameter state. Individual input strands must be homogeneous.
 */
class HybridEnergyProvider final : public EnergyProvider {
public:
    enum class SaltCorrection { Reference, PublishedPair };
    HybridEnergyProvider(std::string sequence, std::vector<uint8_t> materials,
                         std::vector<uint8_t> nicks_before, double celsius,
                         double sodium, SaltCorrection salt);
    ~HybridEnergyProvider() override;
    HybridEnergyProvider(const HybridEnergyProvider&) = delete;
    HybridEnergyProvider& operator=(const HybridEnergyProvider&) = delete;
    double thermal_energy() const;
    double log_global_factor() const;
    bool pair_allowed(int i, int j) const;
    double log_hairpin(int i, int j) const override;
    double log_interior(int i, int j, int k, int l) const override;
    double log_unpaired(FaceMode mode, int family, int i) const override;
    double log_branch(FaceMode mode, int family, int i, int j) const override;
    double log_closing(FaceMode mode, int family, int i, int j) const override;
private:
    struct Impl;
    std::unique_ptr<Impl> p_;
};

/* Official RNAlib single-RNA algorithm. Its PF uses the same model and exact
 * continuous pair-salt correction as HybridEnergyProvider. The MFE uses
 * RNAlib's integer energy convention; it is a distinct observable. */
std::unique_ptr<Calculation> make_vienna_rna_calculation(
    const aptadp_input_v1& input, double celsius, double sodium,
    HybridEnergyProvider::SaltCorrection salt, double& thermal);
}}
#endif
