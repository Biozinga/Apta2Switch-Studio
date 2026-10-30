/* Model definitions and the attributed Sugimoto95 H/S data source are recorded
 * in ../hybrid_motifs.py. The 16-value table is passed through a separate header
 * generated from that attributed source during validation builds. Homogeneous
 * loop motifs are supplied by the declared RNAlib provider below.
 */
#include "hybrid_provider.hpp"
#include "hybrid_parameters.hpp"
extern "C" {
#include <ViennaRNA/model.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/params/io.h>
#include <ViennaRNA/eval/hairpin.h>
#include <ViennaRNA/eval/internal.h>
#include <ViennaRNA/eval/exterior.h>
#include <ViennaRNA/eval/multibranch.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/constraints/hard.h>
#include <ViennaRNA/constraints/soft.h>
#include <ViennaRNA/datastructures/dp_matrices.h>
#include <ViennaRNA/mfe/global.h>
#include <ViennaRNA/partfunc/global.h>
#include <ViennaRNA/sampling/basic.h>
#include <ViennaRNA/subopt/wuchty.h>
#include <ViennaRNA/utils/basic.h>
}
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <utility>

namespace aptaswitch { namespace thermodynamics {
namespace {
constexpr double negative_infinity = -std::numeric_limits<double>::infinity();
std::mutex parameter_mutex;
double log_weight(double weight) {
    if (!std::isfinite(weight) || weight < 0)
        throw std::range_error("Nonfinite or negative Vienna motif weight");
    return weight == 0 ? negative_infinity : std::log(weight);
}
int base_code(char base) {
    switch (base) { case 'A': return 1; case 'C': return 2; case 'G': return 3;
        case 'T': case 'U': return 4; default: throw std::invalid_argument("Invalid base"); }
}
int pair_code(char a, char b) {
    const int x = base_code(a), y = base_code(b);
    if (x == 2 && y == 3) return 1;
    if (x == 3 && y == 2) return 2;
    if (x == 3 && y == 4) return 3;
    if (x == 4 && y == 3) return 4;
    if (x == 1 && y == 4) return 5;
    if (x == 4 && y == 1) return 6;
    return 0;
}
void family_valid(int family) {
    if (family < 1 || family > 3) throw std::invalid_argument("Invalid face material family");
}
}

struct HybridEnergyProvider::Impl {
    std::string sequence, rna_projection;
    std::vector<uint8_t> materials, nicks;
    std::vector<int> codes, cut_prefix;
    using Parameters = std::unique_ptr<vrna_exp_param_t, decltype(&std::free)>;
    Parameters dna{nullptr, &std::free}, rna{nullptr, &std::free};
    double kelvin, thermal, salt_log, global_log;
    int n;

    Impl(std::string seq, std::vector<uint8_t> mat, std::vector<uint8_t> cut,
         double celsius, double sodium, SaltCorrection salt)
        : sequence(std::move(seq)), rna_projection(sequence), materials(std::move(mat)),
          nicks(std::move(cut)), kelvin(celsius + 273.15), thermal(0), salt_log(0), global_log(0),
          n(static_cast<int>(sequence.size())) {
        if (sequence.empty() || sequence.size() > static_cast<size_t>(std::numeric_limits<int>::max())
            || materials.size() != sequence.size() || nicks.size() != sequence.size() + 1)
            throw std::invalid_argument("Invalid native hybrid sequence dimensions");
        if (!std::isfinite(celsius) || celsius < 25 || celsius > 55 || !std::isfinite(sodium))
            throw std::invalid_argument("Experimental model requires 25..55 Celsius and finite sodium");
        if ((salt == SaltCorrection::Reference && sodium != 1.)
            || (salt == SaltCorrection::PublishedPair && (sodium < .12 || sodium > 1.))
            || (salt != SaltCorrection::Reference && salt != SaltCorrection::PublishedPair))
            throw std::invalid_argument("Unsupported experimental salt model or concentration");
        if (nicks.front() || nicks.back()) throw std::invalid_argument("Unexpected endpoint nick");
        cut_prefix.resize(n + 1, 0);
        std::vector<std::pair<int, std::string>> strands;
        int start = 0, rna_count = 0;
        for (int i = 0; i < n; ++i) {
            if (materials[i] != 1 && materials[i] != 2)
                throw std::invalid_argument("Material must be DNA=1 or RNA=2");
            if (nicks[i] > 1) throw std::invalid_argument("Nicks must be binary");
            if ((materials[i] == 1 && sequence[i] == 'U') || (materials[i] == 2 && sequence[i] == 'T'))
                throw std::invalid_argument("Base alphabet disagrees with strand chemistry");
            codes.push_back(base_code(sequence[i]));
            if (sequence[i] == 'T') rna_projection[i] = 'U';
            if (i && !nicks[i] && materials[i] != materials[i-1])
                throw std::invalid_argument("Chimeric backbones are outside the model");
            if (nicks[i]) { strands.emplace_back(materials[start], sequence.substr(start, i-start)); start = i; }
            cut_prefix[i+1] = cut_prefix[i] + nicks[i];
            rna_count += materials[i] == 2;
        }
        strands.emplace_back(materials[start], sequence.substr(start));
        {
            std::lock_guard<std::mutex> lock(parameter_mutex);
            // Vienna's parameter cache compares model details, not the table
            // contents. Construct details AFTER each table load so its updated
            // RNA/DNA geometry cannot accidentally reuse the other snapshot.
            auto details = [&]() {
                vrna_md_t md;
                vrna_md_set_default(&md);
                md.temperature = celsius; md.salt = VRNA_MODEL_DEFAULT_SALT; md.dangles = 0;
                md.noGU = 0; md.noGUclosure = 0; md.noLP = 0; md.uniq_ML = 1; md.min_loop_size = 3;
                md.special_hp = 1; md.circ = 0; md.gquad = 0; md.logML = 0; md.betaScale = 1.;
                md.pf_smooth = 1;
                return md;
            };
            if (!vrna_params_load_DNA_Mathews2004()) throw std::runtime_error("Cannot load DNA Mathews2004");
            auto md = details();
            dna.reset(vrna_exp_params(&md));
            if (!vrna_params_load_RNA_Turner2004()) throw std::runtime_error("Cannot load RNA Turner2004");
            md = details();
            rna.reset(vrna_exp_params(&md));
        }
        if (!dna || !rna) throw std::bad_alloc();
        thermal = rna->kT / 1000.;
        if (dna->kT != rna->kT || !std::isfinite(thermal) || thermal <= 0)
            throw std::range_error("Inconsistent thermal energies");
        salt_log = salt == SaltCorrection::PublishedPair ? .114 * std::log(sodium) * kelvin / 310.15 / thermal : 0.;
        const double rna_fraction = static_cast<double>(rna_count) / n;
        const bool rna_majority = rna_fraction >= .5;
        const double majority_fraction = rna_majority ? rna_fraction : 1.-rna_fraction;
        const double mono_init = -thermal * log_weight((rna_majority ? rna : dna)->expDuplexInit);
        const double hybrid_init = hybrid_hs_init[0] - kelvin * hybrid_hs_init[1] / 1000.;
        const double association = (strands.size()-1) * (hybrid_init + (2*majority_fraction-1)*(mono_init-hybrid_init));
        size_t symmetry = 0;
        for (size_t offset = 0; offset < strands.size(); ++offset) {
            bool same = true;
            for (size_t j = 0; j < strands.size(); ++j)
                same = same && strands[j] == strands[(j+offset) % strands.size()];
            symmetry += same;
        }
        global_log = -association / thermal - std::log(static_cast<double>(symmetry));
    }
    void indices(int i, int j) const {
        if (i < 0 || j < 0 || i >= n || j >= n || i == j)
            throw std::invalid_argument("Invalid motif pair coordinates");
    }
    bool allowed(int i, int j) const {
        indices(i, j);
        const int type = pair_code(sequence[i], sequence[j]);
        return type && ((type != 3 && type != 4) || (materials[i] == 2 && materials[j] == 2));
    }
    bool unbroken(int start, int end) const { return cut_prefix[end+1] == cut_prefix[start+1]; }
    vrna_exp_param_t *params(int material) const { return material == 1 ? dna.get() : rna.get(); }
    double terminal(int i, int j) const {
        if (!allowed(i,j)) return negative_infinity;
        if (materials[i] != materials[j]) return 0.;
        return log_weight(vrna_exp_E_exterior_stem(pair_code(sequence[i], sequence[j]), -1, -1, params(materials[i])));
    }
    double stack(int i, int j, int k, int l) const {
        const char a = materials[i] == 2 ? sequence[i] : sequence[l];
        const char b = materials[i] == 2 ? sequence[k] : sequence[j];
        const int index = (base_code(a)-1)*4 + base_code(b)-1;
        return -(hybrid_hs_stacks[index][0] - kelvin*hybrid_hs_stacks[index][1]/1000.) / thermal;
    }
    double interior(int i, int j, int k, int l, int material) const {
        return log_weight(vrna_exp_E_internal(k-i-1, j-l-1,
            pair_code(sequence[i],sequence[j]), pair_code(sequence[l],sequence[k]),
            codes[i+1],codes[j-1],codes[k-1],codes[l+1],params(material)));
    }
};

HybridEnergyProvider::HybridEnergyProvider(std::string sequence, std::vector<uint8_t> materials,
        std::vector<uint8_t> nicks, double celsius, double sodium, SaltCorrection salt)
    : p_(std::make_unique<Impl>(std::move(sequence),std::move(materials),std::move(nicks),celsius,sodium,salt)) {}
HybridEnergyProvider::~HybridEnergyProvider() = default;
double HybridEnergyProvider::thermal_energy() const { return p_->thermal; }
double HybridEnergyProvider::log_global_factor() const { return p_->global_log; }
bool HybridEnergyProvider::pair_allowed(int i, int j) const { return p_->allowed(i,j); }

double HybridEnergyProvider::log_hairpin(int i, int j) const {
    if (!p_->allowed(i,j) || i >= j || j-i-1 < 3 || !p_->unbroken(i,j)) return negative_infinity;
    return log_weight(vrna_exp_E_hairpin(j-i-1, pair_code(p_->sequence[i],p_->sequence[j]),
        p_->codes[i+1],p_->codes[j-1],p_->rna_projection.c_str()+i,p_->params(p_->materials[i]))) + p_->salt_log;
}
double HybridEnergyProvider::log_interior(int i, int j, int k, int l) const {
    if (!p_->allowed(i,j) || !p_->allowed(k,l)) return negative_infinity;
    if (!(i<k && k<l && l<j)) throw std::invalid_argument("Invalid interior nesting");
    const int left=k-i-1, right=j-l-1, size=left+right;
    if (size > 30 || !p_->unbroken(i,k) || !p_->unbroken(l,j)) return negative_infinity;
    if (p_->materials[i] == p_->materials[j])
        return p_->interior(i,j,k,l,p_->materials[i]) + p_->salt_log;
    double result;
    if (!size) result=p_->stack(i,j,k,l);
    else if (left && right) result=(p_->interior(i,j,k,l,1)+p_->interior(i,j,k,l,2))/2;
    else {
        result=(log_weight(p_->dna->expbulge[size])+log_weight(p_->rna->expbulge[size]))/2;
        result+= size==1 ? p_->stack(i,j,k,l) : p_->terminal(i,j)+p_->terminal(k,l);
    }
    return result+p_->salt_log;
}
double HybridEnergyProvider::log_unpaired(FaceMode mode, int family, int i) const {
    family_valid(family);
    if (i < 0 || i >= p_->n) throw std::invalid_argument("Invalid free-base coordinate");
    if (mode == FaceMode::Exterior) return 0.;
    return family == 3 ? (log_weight(p_->dna->expMLbase)+log_weight(p_->rna->expMLbase))/2
                       : log_weight(p_->params(family)->expMLbase);
}
double HybridEnergyProvider::log_branch(FaceMode mode, int family, int i, int j) const {
    family_valid(family);
    if (!p_->allowed(i,j)) return negative_infinity;
    if (mode == FaceMode::Exterior) return p_->terminal(i,j);
    if (family == 3)
        return (log_weight(p_->dna->expMLintern[1])+log_weight(p_->rna->expMLintern[1]))/2 + p_->terminal(i,j);
    return log_weight(vrna_exp_E_multibranch_stem(pair_code(p_->sequence[i],p_->sequence[j]), -1,-1,p_->params(family)));
}
double HybridEnergyProvider::log_closing(FaceMode mode, int family, int i, int j) const {
    double result=log_branch(mode,family,j,i)+p_->salt_log;
    if (mode == FaceMode::Multiloop)
        result+=family==3 ? (log_weight(p_->dna->expMLclosing)+log_weight(p_->rna->expMLclosing))/2
                          : log_weight(p_->params(family)->expMLclosing);
    return result;
}

namespace {
/* RNAlib parameter loading and its RNG both use process-global state. This
 * mutex is shared with HybridEnergyProvider's parameter snapshot constructor.
 * Fold compounds and their DP arrays remain private to individual handles. */
class ViennaRnaCalculation final : public Calculation {
    using Compound = std::unique_ptr<vrna_fold_compound_t, decltype(&vrna_fold_compound_free)>;
    Compound fold_{nullptr, &vrna_fold_compound_free};
    int n_;
    uint64_t memory_;
    double thermal_ = 0., pair_factor_ = 1.;
    int pair_mfe_ = 0;
    bool empty_ = false, pf_ready_ = false, mfe_ready_ = false;
    double log_z_ = negative_infinity, minimum_ = 0.;
    std::string structure_;
    std::string sequence_;
    std::vector<int32_t> required_;
    std::vector<uint8_t> pair_mask_;
    double celsius_, sodium_;
    HybridEnergyProvider::SaltCorrection salt_;
    uint64_t memory_limit_;

    struct Enumeration {
        int length;
        uint64_t limit;
        char* output;
        double* energies;
        uint64_t* count;
        uint8_t* truncated;
        const std::string* first;
        int minimum_centikcal, delta;
        bool invalid = false;
    };
    static void suboptimal(const char* structure, float energy, void* data) noexcept {
        if (!structure) return;
        auto& result = *static_cast<Enumeration*>(data);
        if (!std::isfinite(energy) || std::strlen(structure) != static_cast<size_t>(result.length)) {
            result.invalid = true;
            return;
        }
        const long centikcal = std::lround(static_cast<double>(energy) * 100.);
        if (centikcal < result.minimum_centikcal) { result.invalid = true; return; }
        if (centikcal > static_cast<long>(result.minimum_centikcal) + result.delta
            || *result.first == structure) return;
        if (*result.count == result.limit) { *result.truncated = 1; return; }
        // uniq_ML and dangles=0 give each secondary structure one traceback.
        std::memcpy(result.output + static_cast<size_t>(*result.count) * (result.length + 1),
                    structure, result.length + 1);
        result.energies[*result.count] = centikcal / 100.;
        ++*result.count;
    }

    static bool closes_pair(unsigned char decomposition) {
        return decomposition == VRNA_DECOMP_PAIR_HP || decomposition == VRNA_DECOMP_PAIR_IL
            || decomposition == VRNA_DECOMP_PAIR_ML;
    }
    static FLT_OR_DBL salt_pf(int, int, int, int, unsigned char decomposition, void* data) {
        return closes_pair(decomposition) ? static_cast<ViennaRnaCalculation*>(data)->pair_factor_ : 1.;
    }
    static int salt_mfe(int, int, int, int, unsigned char decomposition, void* data) {
        return closes_pair(decomposition) ? static_cast<ViennaRnaCalculation*>(data)->pair_mfe_ : 0;
    }
    void ensure_mfe() {
        if (mfe_ready_) return;
        if (empty_) throw EmptyEnsemble("The single-RNA ensemble is empty");
        std::vector<char> result(static_cast<size_t>(n_) + 1, '\0');
        const double reported = vrna_mfe(fold_.get(), result.data());
        if (!std::isfinite(reported)) throw std::range_error("Nonfinite RNAlib minimum");
        if (reported >= INF / 100.) {
            empty_ = true;
            throw EmptyEnsemble("The single-RNA ensemble is empty");
        }
        // The public return type is float. Read the official integer DP result
        // to avoid adding a float conversion error to its centikcal convention.
        minimum_ = fold_->matrices->f5[n_] / 100.;
        structure_.assign(result.data());
        if (structure_.size() != static_cast<size_t>(n_))
            throw std::range_error("RNAlib minimum traceback has the wrong length");
        mfe_ready_ = true;
    }
    void ensure_pf() {
        if (pf_ready_) return;
        if (empty_) { pf_ready_ = true; return; }
        try { ensure_mfe(); }
        catch (const EmptyEnsemble&) { pf_ready_ = true; return; }
        double scaling_minimum = minimum_;
        vrna_exp_params_rescale(fold_.get(), &scaling_minimum);
        vrna_pf(fold_.get(), nullptr);
        const double scaled_q = fold_->exp_matrices->q[fold_->iindx[1] - n_];
        if (!std::isfinite(scaled_q) || scaled_q < 0.)
            throw std::range_error("Invalid RNAlib partition function");
        if (scaled_q == 0.) {
            // A successful constrained MFE proves the ensemble is nonempty.
            // Do not disguise numerical underflow as an impossible structure.
            throw std::range_error("RNAlib partition function underflowed for a nonempty ensemble");
        }
        log_z_ = std::log(scaled_q) + n_ * std::log(fold_->exp_params->pf_scale);
        if (!std::isfinite(log_z_)) throw std::range_error("Nonfinite RNAlib log partition");
        pf_ready_ = true;
    }
public:
    ViennaRnaCalculation(const aptadp_input_v1& input, double celsius, double sodium,
                        HybridEnergyProvider::SaltCorrection salt) : n_(input.length), memory_(0),
            celsius_(celsius), sodium_(sodium), salt_(salt), memory_limit_(input.memory_limit_bytes) {
        if (input.struct_size != sizeof(input) || !input.length || input.length > 32767
            || !input.sequence || !input.materials || !input.nicks_before || !input.pair_allowed)
            throw std::invalid_argument("Invalid single-RNA input buffers or length (1..32767 required)");
        if (input.min_hairpin_size != 3 || input.max_interior_size != 30)
            throw std::invalid_argument("Single-RNA model requires min_hairpin=3, max_interior=30");
        if (!std::isfinite(celsius) || celsius < 25. || celsius > 55. || !std::isfinite(sodium)
            || (salt == HybridEnergyProvider::SaltCorrection::Reference && sodium != 1.)
            || (salt == HybridEnergyProvider::SaltCorrection::PublishedPair && (sodium < .12 || sodium > 1.))
            || (salt != HybridEnergyProvider::SaltCorrection::Reference
                && salt != HybridEnergyProvider::SaltCorrection::PublishedPair))
            throw std::invalid_argument("Unsupported single-RNA temperature or salt conditions");
        // Conservative bound for both RNAlib DP families, hard constraints,
        // probabilities, returned matrix, parameter snapshots and auxiliaries.
        memory_ = 128ULL * input.length * input.length + 2048ULL * input.length + 4ULL * 1024 * 1024;
        const uint64_t limit = input.memory_limit_bytes ? input.memory_limit_bytes : 512ULL * 1024 * 1024;
        if (memory_ > limit) throw std::bad_alloc();
        std::string sequence(input.sequence, input.length);
        std::vector<int> required(n_, -1);
        for (int i = 0; i <= n_; ++i) {
            if (input.nicks_before[i]) throw std::invalid_argument("Official RNA factory requires exactly one strand");
            if (i == n_) break;
            if (input.materials[i] != 2 || std::string("ACGU").find(sequence[i]) == std::string::npos)
                throw std::invalid_argument("Official RNA factory requires a homogeneous RNA strand");
            required[i] = input.constraints ? input.constraints[i] : -1;
            if (required[i] < -1 || required[i] >= n_)
                throw std::invalid_argument("Constraint index is out of range");
        }
        for (int i = 0; i < n_; ++i) {
            for (int j = 0; j < n_; ++j) {
                const uint8_t expected = i != j && pair_code(sequence[i], sequence[j]);
                if (input.pair_allowed[static_cast<size_t>(i) * n_ + j] != expected)
                    throw std::invalid_argument("Pair mask differs from the RNA model chemistry");
            }
            const int j = required[i];
            if (j >= 0 && required[j] != i) throw std::invalid_argument("Constraints are not reciprocal");
            if (j > i) {
                if (j - i <= 3 || !pair_code(sequence[i], sequence[j])) empty_ = true;
                for (int k = i + 1; k < j; ++k)
                    if (required[k] > j) throw std::invalid_argument("Pseudoknotted constraints are unsupported");
            }
        }
        sequence_ = sequence;
        required_.assign(required.begin(), required.end());
        pair_mask_.assign(input.pair_allowed, input.pair_allowed + static_cast<size_t>(n_) * n_);
        {
            std::lock_guard<std::mutex> lock(parameter_mutex);
            if (!vrna_params_load_RNA_Turner2004()) throw std::runtime_error("Cannot load RNA Turner2004");
            vrna_md_t md;
            vrna_md_set_default(&md);
            md.temperature = celsius; md.salt = VRNA_MODEL_DEFAULT_SALT; md.dangles = 0;
            md.noGU = 0; md.noGUclosure = 0; md.noLP = 0; md.uniq_ML = 1; md.min_loop_size = 3;
            md.special_hp = 1; md.circ = 0; md.gquad = 0; md.logML = 0; md.betaScale = 1.;
            md.pf_smooth = 1; md.compute_bpp = 1;
            fold_.reset(vrna_fold_compound(sequence.c_str(), &md, VRNA_OPTION_MFE | VRNA_OPTION_PF));
        }
        if (!fold_ || !fold_->exp_params || !fold_->params) throw std::bad_alloc();
        thermal_ = fold_->exp_params->kT / 1000.;
        if (!std::isfinite(thermal_) || thermal_ <= 0.) throw std::range_error("Invalid RNAlib thermal energy");
        const double salt_energy = salt == HybridEnergyProvider::SaltCorrection::PublishedPair
            ? -.114 * std::log(sodium) * (celsius + 273.15) / 310.15 : 0.;
        pair_factor_ = std::exp(-salt_energy / thermal_);
        // MFE callbacks take integer centikcal. Ties round away from zero,
        // matching C++ lround; PF intentionally retains the continuous value.
        pair_mfe_ = static_cast<int>(std::lround(100. * salt_energy));
        vrna_sc_add_data(fold_.get(), this, nullptr);
        vrna_sc_add_exp_f(fold_.get(), salt_pf);
        vrna_sc_add_f(fold_.get(), salt_mfe);
        if (!empty_) {
            for (int i = 0; i < n_; ++i) {
                if (required[i] == i)
                    vrna_hc_add_up(fold_.get(), i + 1, VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS
                        | VRNA_CONSTRAINT_CONTEXT_ENFORCE);
                else if (required[i] > i)
                    vrna_hc_add_bp(fold_.get(), i + 1, required[i] + 1,
                        VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS | VRNA_CONSTRAINT_CONTEXT_ENFORCE);
            }
        }
    }
    int length() const override { return n_; }
    uint64_t memory() const override { return memory_; }
    double thermal_energy() const { return thermal_; }
    double log_partition() override { ensure_pf(); return log_z_; }
    void pairs(double* matrix, double& adjustment) override {
        ensure_pf();
        if (empty_) throw EmptyEnsemble("The single-RNA ensemble is empty");
        std::fill(matrix, matrix + static_cast<size_t>(n_) * n_, 0.);
        adjustment = 0.;
        for (int i = 0; i < n_; ++i) for (int j = i + 1; j < n_; ++j) {
            const double value = fold_->exp_matrices->probs[fold_->iindx[i + 1] - (j + 1)];
            if (!std::isfinite(value) || value < -1e-12 || value > 1. + 1e-12)
                throw std::range_error("RNAlib probability is outside [0,1]");
            const double bounded = std::max(0., std::min(1., value));
            adjustment = std::max(adjustment, std::abs(value - bounded));
            matrix[static_cast<size_t>(i) * n_ + j] = matrix[static_cast<size_t>(j) * n_ + i] = bounded;
        }
        for (int i = 0; i < n_; ++i) {
            long double sum = 0.;
            for (int j = 0; j < n_; ++j) sum += matrix[static_cast<size_t>(i) * n_ + j];
            const double unpaired = static_cast<double>(1. - sum);
            if (unpaired < -1e-12) throw std::range_error("RNAlib probability row sums to more than one");
            matrix[static_cast<size_t>(i) * n_ + i] = std::max(0., unpaired);
            adjustment = std::max(adjustment, std::max(0., -unpaired));
        }
    }
    std::pair<double, std::string> minimum() override { ensure_mfe(); return {minimum_, structure_}; }
    double minimum_pair_energies(double* matrix) const override {
        auto& self = *const_cast<ViennaRnaCalculation*>(this);
        self.ensure_mfe();
        std::fill(matrix, matrix + static_cast<size_t>(n_) * n_, std::numeric_limits<double>::infinity());
        std::vector<uint8_t> materials(n_, 2), nicks(static_cast<size_t>(n_) + 1, 0);
        auto required = required_;
        aptadp_input_v1 input{sizeof(aptadp_input_v1), static_cast<uint32_t>(n_), sequence_.c_str(),
            materials.data(), nicks.data(), pair_mask_.data(), required.data(), 3, 30, memory_limit_};
        // RNAlib exposes no all-pairs MFE outside pass. Repeated constrained
        // official folds preserve the exact observable, at O(n^5) worst-case
        // time and bounded O(n^2) workspace. This is not used for ordinary PF.
        for (int i = 0; i < n_; ++i) for (int j = i + 4; j < n_; ++j) {
            if (!pair_mask_[static_cast<size_t>(i) * n_ + j]
                || (required_[i] != -1 && required_[i] != j)
                || (required_[j] != -1 && required_[j] != i)) continue;
            bool crossing = false;
            for (int k = 0; k < n_; ++k) {
                const int l = required_[k];
                if (l > k && ((k < i && i < l && l < j) || (i < k && k < j && j < l))) {
                    crossing = true;
                    break;
                }
            }
            if (crossing) continue;
            if (required_[i] == j) {
                matrix[static_cast<size_t>(i) * n_ + j] = matrix[static_cast<size_t>(j) * n_ + i] = minimum_;
                continue;
            }
            required[i] = j; required[j] = i;
            try {
                ViennaRnaCalculation conditioned(input, celsius_, sodium_, salt_);
                const double energy = conditioned.minimum().first;
                matrix[static_cast<size_t>(i) * n_ + j] = matrix[static_cast<size_t>(j) * n_ + i] = energy;
            } catch (const EmptyEnsemble&) {
                // A chemically possible forced pair can still conflict with
                // loop-size bounds imposed by the remaining constraints.
            }
            required[i] = required_[i]; required[j] = required_[j];
        }
        return minimum_;
    }
    double minimum_structures(double tolerance, uint64_t limit, char* output, double* energies,
                              uint64_t& count, uint8_t& truncated) const override {
        if (!std::isfinite(tolerance) || tolerance < 0. || tolerance > (INF / 200.))
            throw std::invalid_argument("RNAlib minimum tolerance is outside the supported energy range");
        auto& self = *const_cast<ViennaRnaCalculation*>(this);
        self.ensure_mfe();
        count = 0; truncated = 0;
        if (!limit) { truncated = 1; return minimum_; }
        std::memcpy(output, structure_.c_str(), n_ + 1);
        energies[0] = minimum_; count = 1;
        // RNAlib's integer energy band is exact for this rounded-energy MFE
        // observable. The tiny epsilon only compensates decimal-to-double
        // representation at a centikcal boundary.
        const int delta = static_cast<int>(std::floor(tolerance * 100. + 1e-9));
        Enumeration result{n_, limit, output, energies, &count, &truncated, &structure_,
            fold_->matrices->f5[n_], delta};
        {
            // Wuchty's public API updates global density-of-states counters.
            std::lock_guard<std::mutex> lock(parameter_mutex);
            std::fill(density_of_states, density_of_states + MAXDOS + 1, 0);
            vrna_subopt_cb(fold_.get(), delta, suboptimal, &result);
        }
        if (result.invalid) throw std::range_error("RNAlib alternatives violate the minimum energy window");
        // The callback bounds retained output, not RNAlib's internal search.
        // RNAlib completes traversal even after another solution proves that
        // output was truncated; it has no public early-stop callback.
        return minimum_;
    }
    void sample(uint64_t seed, uint64_t count, char* output) override {
        ensure_pf();
        if (empty_) throw EmptyEnsemble("The single-RNA ensemble is empty");
        std::lock_guard<std::mutex> lock(parameter_mutex);
        // RNAlib accepts a 32-bit seed; record this deterministic 64-bit fold
        // in the factory contract instead of pretending to share native RNG.
        vrna_init_rand_seed(static_cast<unsigned int>(seed ^ (seed >> 32)));
        for (uint64_t k = 0; k < count; ++k) {
            std::unique_ptr<char, decltype(&std::free)> result(vrna_pbacktrack(fold_.get()), &std::free);
            if (!result || std::strlen(result.get()) != static_cast<size_t>(n_))
                throw std::range_error("RNAlib stochastic traceback failed");
            std::memcpy(output + static_cast<size_t>(k) * (n_ + 1), result.get(), n_ + 1);
        }
    }
};
}

std::unique_ptr<Calculation> make_vienna_rna_calculation(const aptadp_input_v1& input,
        double celsius, double sodium, HybridEnergyProvider::SaltCorrection salt, double& thermal) {
    auto calculation = std::make_unique<ViennaRnaCalculation>(input, celsius, sodium, salt);
    thermal = calculation->thermal_energy();
    return calculation;
}
}}
