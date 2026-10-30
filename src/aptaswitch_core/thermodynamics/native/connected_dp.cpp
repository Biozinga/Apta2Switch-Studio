#include "connected_dp.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <new>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using aptaswitch::thermodynamics::EnergyProvider;
using aptaswitch::thermodynamics::FaceMode;
constexpr double NEG = -std::numeric_limits<double>::infinity();
/* One paired field and twenty forest fields, all triangular. A multiloop
 * forest has no exposed nick, one or >=2 branches, and a material mask. An
 * exterior forest has zero or one exposed nick; its branch count is irrelevant.
 * The zero-branch multiloop path is an O(1) prefix query. Families 1/2/3 denote
 * the completed face's DNA/RNA/mixed mask, not its enclosed child interiors.
 * Same productions serve log-sum, max-log, outside adjoints and traceback.
 * Time O(n^3 + n^2 L^2), L<=30; memory O(n^2), without a production tape or
 * precomputed four-index interior-energy tensor. */
constexpr int BF = 20, ROOT = 21, NF = 21;
constexpr std::array<int, 5> FAMILY{{1, 2, 3, 3, 3}};
constexpr std::array<int, 5> MASK{{1, 2, 1, 2, 3}};

using Empty = aptaswitch::thermodynamics::EmptyEnsemble;
uint64_t mul(uint64_t a, uint64_t b) {
    if (b && a > std::numeric_limits<uint64_t>::max() / b)
        throw std::invalid_argument("Native workspace size overflows uint64");
    return a * b;
}
uint64_t add(uint64_t a, uint64_t b) {
    if (a > std::numeric_limits<uint64_t>::max() - b)
        throw std::invalid_argument("Native workspace size overflows uint64");
    return a + b;
}
double logadd(double a, double b) {
    if (a == NEG) return b;
    if (b == NEG) return a;
    const double hi = std::max(a, b), lo = std::min(a, b);
    const double result = hi + std::log1p(std::exp(lo - hi));
    if (!std::isfinite(result)) throw std::range_error("Native log-sum exceeds floating-point range");
    return result;
}
double logweight(double value) {
    if (value != NEG && !std::isfinite(value))
        throw std::range_error("Energy provider returned an invalid log weight");
    return value;
}
int combo(int family, int mask) {
    if (family == 1 && mask == 1) return 0;
    if (family == 2 && mask == 2) return 1;
    if (family == 3 && mask >= 1 && mask <= 3) return mask + 1;
    return -1;
}
int field(FaceMode mode, int family, int mask, int count) {
    const int c = combo(family, mask);
    if (c < 0) return -1;
    return mode == FaceMode::Multiloop ? c * 2 + count - 1 : 10 + c * 2 + count;
}
struct Cell { int f, a, b; };
struct Production {
    double coefficient;
    std::array<Cell, 2> children{};
    int count = 0;
};

/* Clearly synthetic providers. These formulas mirror the independent Python
 * fixture, and are deliberately not exposed as physical RNA/DNA models. */
class FixtureProvider final : public EnergyProvider {
    bool unit_;
    double thermal_;
    std::vector<uint8_t> materials_;
    std::vector<int> dna_prefix_;
    static double bias(int family) { return family == 1 ? .23 : family == 2 ? .79 : 1.67; }
    double weight(double energy) const {
        const double result = -energy / thermal_;
        if (!std::isfinite(result)) throw std::range_error("Fixture log weight exceeds floating-point range");
        return result;
    }
    int mask(int first, int stop) const {
        const int dna = dna_prefix_[stop] - dna_prefix_[first];
        return (dna ? 1 : 0) | (stop - first - dna ? 2 : 0);
    }
    double branch_energy(FaceMode mode, int family, int i, int j) const {
        return .09 + .011 * (i + 1) + .017 * (j + 1) + .13 * (materials_[i] | materials_[j])
            + (mode == FaceMode::Multiloop ? .29 : .11) * bias(family);
    }
public:
    FixtureProvider(bool unit, double thermal, const std::vector<uint8_t>& materials)
        : unit_(unit), thermal_(thermal), materials_(materials), dna_prefix_(materials.size() + 1, 0) {
        for (size_t i = 0; i < materials.size(); ++i) dna_prefix_[i + 1] = dna_prefix_[i] + (materials[i] == 1);
    }
    double log_hairpin(int i, int j) const override {
        return unit_ ? 0 : weight(.19 + .07 * (j - i - 1) + .011 * i + .017 * j + bias(mask(i, j + 1)));
    }
    double log_interior(int i, int j, int k, int l) const override {
        return unit_ ? 0 : weight(.23 + .03 * (k - i - 1) + .07 * (j - l - 1)
            + .009 * (i + 2 * k + 3 * l + 5 * j) + bias(mask(i, k + 1) | mask(l, j + 1)));
    }
    double log_unpaired(FaceMode mode, int family, int i) const override {
        return unit_ ? 0 : weight(.013 * (i + 1) + .021 * materials_[i]
            + (mode == FaceMode::Multiloop ? .041 : .019) * bias(family));
    }
    double log_branch(FaceMode mode, int family, int i, int j) const override {
        return unit_ ? 0 : weight(branch_energy(mode, family, i, j));
    }
    double log_closing(FaceMode mode, int family, int i, int j) const override {
        return unit_ ? 0 : weight(branch_energy(mode, family, j, i)
            + (mode == FaceMode::Multiloop ? .57 + .19 * bias(family) : 0));
    }
};

class DP final : public aptaswitch::thermodynamics::Calculation {
    int n_, min_hairpin_, max_interior_;
    double thermal_, global_;
    size_t triangle_, cells_;
    uint64_t memory_, memory_limit_;
    std::string sequence_;
    std::vector<uint8_t> materials_, cuts_, allowed_;
    std::vector<int32_t> required_;
    std::vector<int> cut_prefix_, blocked_prefix_, dna_prefix_;
    std::array<std::vector<double>, 6> unpaired_, free_prefix_;
    std::array<std::vector<int>, 6> zero_prefix_;
    std::unique_ptr<EnergyProvider> provider_;
    std::vector<double> inside_;
    bool ready_ = false;

    size_t index(Cell cell) const {
        const size_t a = static_cast<size_t>(cell.a + 1), b = static_cast<size_t>(cell.b + 1);
        const size_t vertices = static_cast<size_t>(n_) + 2;
        return static_cast<size_t>(cell.f) * triangle_ + a * (2 * vertices - a - 1) / 2 + b - a - 1;
    }
    double value(const std::vector<double>& table, Cell cell) const { return table[index(cell)]; }
    bool free_run(int first, int stop) const { return blocked_prefix_[first] == blocked_prefix_[stop]; }
    int nicks(int a, int b) const { return cut_prefix_[b] - cut_prefix_[a]; }
    int material_mask(int first, int stop) const {
        const int dna = dna_prefix_[stop] - dna_prefix_[first];
        return (dna ? 1 : 0) | (stop - first - dna ? 2 : 0);
    }
    static int ui(FaceMode mode, int family) { return (mode == FaceMode::Multiloop ? 0 : 3) + family - 1; }
    double free_multiloop(int family, int a, int b) const {
        if (nicks(a, b) || !free_run(a + 1, b)) return NEG;
        const int u = ui(FaceMode::Multiloop, family);
        if (zero_prefix_[u][b] != zero_prefix_[u][a + 1]) return NEG;
        return free_prefix_[u][b] - free_prefix_[u][a + 1];
    }
    template<class Emit> void forest(Cell cell, const std::vector<double>& table, Emit&& emit) const {
        const bool multi = cell.f < 10;
        const int c = (cell.f % 10) / 2, family = FAMILY[c], target_mask = MASK[c];
        const int count = multi ? cell.f % 2 + 1 : cell.f % 2;
        const FaceMode mode = multi ? FaceMode::Multiloop : FaceMode::Exterior;
        const int a = cell.a, b = cell.b, t = a + 1;
        if (b == a + 1) return; // empty walks have mask=0 and are handled directly
        const int edge = cuts_[t];
        if ((multi && edge) || (!multi && edge > count)) return;
        auto tail = [&](double coefficient, Cell stem, bool has_stem, int start, int added_mask, int branches) {
            if (coefficient == NEG) return;
            if (multi) {
                if (branches == 1 && count == 1) {
                    const int mask = material_mask(start + 1, b);
                    const double empty = free_multiloop(family, start, b);
                    if (empty != NEG && (mask | added_mask) == target_mask)
                        emit(Production{coefficient + empty, {stem, {}}, 1});
                    return;
                }
                const int first_h = branches ? 1 : count;
                const int last_h = branches ? 2 : count;
                for (int h = first_h; h <= last_h; ++h) {
                    for (int mask = 1; mask <= 3; ++mask) {
                        if ((mask | added_mask) != target_mask) continue;
                        const int f = field(mode, family, mask, h);
                        if (f < 0) continue;
                        const Cell suffix{f, start, b};
                        if (value(table, suffix) == NEG) continue;
                        emit(has_stem ? Production{coefficient, {stem, suffix}, 2}
                                      : Production{coefficient, {suffix, {}}, 1});
                    }
                }
            } else {
                const int remaining = count - edge;
                if (b == start + 1) {
                    if (remaining == cuts_[b] && added_mask == target_mask)
                        emit(has_stem ? Production{coefficient, {stem, {}}, 1}
                                      : Production{coefficient, {}, 0});
                    return;
                }
                for (int mask = 1; mask <= 3; ++mask) {
                    if ((mask | added_mask) != target_mask) continue;
                    const int f = field(mode, family, mask, remaining);
                    if (f < 0) continue;
                    const Cell suffix{f, start, b};
                    if (value(table, suffix) == NEG) continue;
                    emit(has_stem ? Production{coefficient, {stem, suffix}, 2}
                                  : Production{coefficient, {suffix, {}}, 1});
                }
            }
        };
        if (free_run(t, t + 1))
            tail(unpaired_[ui(mode, family)][t], {}, false, t, materials_[t], 0);
        for (int k = t + 1; k < b; ++k) {
            const Cell stem{BF, t, k};
            if (value(table, stem) == NEG) continue;
            tail(logweight(provider_->log_branch(mode, family, t, k)), stem, true, k,
                 materials_[t] | materials_[k], 1);
        }
    }
    template<class Emit> void paired(Cell cell, const std::vector<double>& table, Emit&& emit) const {
        const int i = cell.a, j = cell.b;
        if (!allowed_[static_cast<size_t>(i) * n_ + j]) return;
        if (!nicks(i, j) && free_run(i + 1, j) && j - i - 1 >= min_hairpin_)
            emit(Production{logweight(provider_->log_hairpin(i, j)), {}, 0});
        for (int left = 0; left <= std::min(max_interior_, j - i - 3); ++left) {
            const int k = i + left + 1;
            for (int right = 0; right <= std::min(max_interior_ - left, j - i - 3 - left); ++right) {
                const int l = j - right - 1;
                const Cell child{BF, k, l};
                if (nicks(i, k) + nicks(l, j) || !free_run(i + 1, k) || !free_run(l + 1, j)
                    || value(table, child) == NEG) continue;
                emit(Production{logweight(provider_->log_interior(i, j, k, l)), {child, {}}, 1});
            }
        }
        const int boundary = materials_[i] | materials_[j];
        if (j == i + 1 && cuts_[j])
            emit(Production{logweight(provider_->log_closing(FaceMode::Exterior, boundary, i, j)), {}, 0});
        for (int family = 1; family <= 3; ++family) {
            for (int mask = 1; mask <= 3; ++mask) {
                if ((mask | boundary) != family) continue;
                for (FaceMode mode : {FaceMode::Multiloop, FaceMode::Exterior}) {
                    const int f = field(mode, family, mask, mode == FaceMode::Multiloop ? 2 : 1);
                    if (f < 0) continue;
                    const Cell child{f, i, j};
                    if (value(table, child) != NEG)
                        emit(Production{logweight(provider_->log_closing(mode, family, i, j)), {child, {}}, 1});
                }
            }
        }
    }
    template<class Emit> void productions(Cell cell, const std::vector<double>& table, Emit&& emit) const {
        if (cell.f == ROOT) {
            for (int family = 1; family <= 3; ++family) {
                const Cell child{field(FaceMode::Exterior, family, family, 0), -1, n_};
                // The global constant cancels from every conditional law.
                // Apply it only to reported Z/energy, never a traceback draw.
                if (value(table, child) != NEG) emit(Production{0, {child, {}}, 1});
            }
        } else if (cell.f == BF) paired(cell, table, emit);
        else forest(cell, table, emit);
    }
    double score(const Production& p, const std::vector<double>& table) const {
        if (p.coefficient == NEG) return NEG;
        double result = p.coefficient;
        for (int k = 0; k < p.count; ++k) {
            const double v = value(table, p.children[k]);
            if (v == NEG) return NEG;
            result += v;
        }
        if (!std::isfinite(result)) throw std::range_error("Native log-product exceeds floating-point range");
        return result;
    }
    double evaluate(Cell cell, const std::vector<double>& table, bool maximum) const {
        double result = NEG;
        productions(cell, table, [&](const Production& p) {
            const double term = score(p, table);
            result = maximum ? std::max(result, term) : logadd(result, term);
        });
        return result;
    }
    void fill(std::vector<double>& table, bool maximum) const {
        for (int span = 1; span <= n_ + 1; ++span) {
            for (int a = -1; a + span <= n_; ++a) {
                const int b = a + span;
                for (int f = 0; f < BF; ++f) {
                    const Cell cell{f, a, b};
                    table[index(cell)] = evaluate(cell, table, maximum);
                }
            }
            for (int i = 0; i + span < n_; ++i) {
                const Cell cell{BF, i, i + span};
                table[index(cell)] = evaluate(cell, table, maximum);
            }
        }
    }
    void ensure_inside() {
        if (!ready_) {
            // Minimum-only requests use their own max-log table. Do not keep
            // an unused log-sum table alongside it (one quadratic allocation).
            inside_.assign(cells_, NEG);
            fill(inside_, false);
            ready_ = true;
        }
    }
    std::vector<double> outside_values(const std::vector<double>& table, bool maximum) const {
        std::vector<double> outside(cells_, NEG);
        auto propagate = [&](Cell cell, double adjoint) {
            if (adjoint == NEG) return;
            productions(cell, table, [&](const Production& p) {
                if (score(p, table) == NEG) return;
                for (int k = 0; k < p.count; ++k) {
                    double contribution = adjoint + p.coefficient;
                    for (int other = 0; other < p.count; ++other)
                        if (other != k) contribution += value(table, p.children[other]);
                    if (!std::isfinite(contribution))
                        throw std::range_error("Outside log weight exceeds floating-point range");
                    const size_t at = index(p.children[k]);
                    outside[at] = maximum ? std::max(outside[at], contribution) : logadd(outside[at], contribution);
                }
            });
        };
        // Both semirings differentiate the same root and production graph.
        // The global constant is excluded here and added to energies once.
        propagate({ROOT, -1, n_}, 0);
        for (int span = n_ + 1; span >= 1; --span) {
            for (int i = 0; i + span < n_; ++i) {
                const Cell cell{BF, i, i + span}; propagate(cell, value(outside, cell));
            }
            for (int a = -1; a + span <= n_; ++a) for (int f = 0; f < BF; ++f) {
                const Cell cell{f, a, a + span}; propagate(cell, value(outside, cell));
            }
        }
        return outside;
    }
    std::string traceback(const std::vector<double>& table, std::mt19937_64* random) const {
        std::string structure(static_cast<size_t>(n_), '.');
        std::vector<Cell> pending{{ROOT, -1, n_}};
        while (!pending.empty()) {
            const Cell cell = pending.back(); pending.pop_back();
            if (cell.f == BF) { structure[cell.a] = '('; structure[cell.b] = ')'; }
            std::vector<std::pair<Production, double>> choices;
            double maximum = NEG;
            productions(cell, table, [&](const Production& p) {
                const double s = score(p, table);
                if (s != NEG) { choices.emplace_back(p, s); maximum = std::max(maximum, s); }
            });
            if (choices.empty()) throw std::logic_error("Traceback has no admissible production");
            size_t selected = 0;
            if (random) {
                double total = 0;
                for (auto& choice : choices) { choice.second = std::exp(choice.second - maximum); total += choice.second; }
                const double u = static_cast<double>((*random)() >> 11) * 0x1p-53;
                const double cutoff = u * total;
                double cumulative = 0;
                selected = choices.size() - 1;
                for (size_t i = 0; i < choices.size(); ++i) {
                    cumulative += choices[i].second;
                    if (cutoff < cumulative) { selected = i; break; }
                }
            } else {
                for (size_t i = 0; i < choices.size(); ++i)
                    if (choices[i].second == maximum) { selected = i; break; }
            }
            const auto& chosen = choices[selected].first;
            for (int k = chosen.count - 1; k >= 0; --k) pending.push_back(chosen.children[k]);
        }
        validate_traceback(structure);
        return structure;
    }
    void validate_traceback(const std::string& structure) const {
        std::vector<int> partners(n_), parents(n_), stack;
        for (int i = 0; i < n_; ++i) partners[i] = parents[i] = i;
        auto find = [&](int i) { while (parents[i] != i) i = parents[i]; return i; };
        for (int i = 0; i < n_; ++i) {
            if (i && !cuts_[i]) parents[find(i - 1)] = find(i);
            if (structure[i] == '(') stack.push_back(i);
            else if (structure[i] == ')') {
                if (stack.empty()) throw std::logic_error("Unbalanced native traceback");
                const int j = stack.back(); stack.pop_back();
                partners[i] = j; partners[j] = i; parents[find(i)] = find(j);
            }
        }
        if (!stack.empty()) throw std::logic_error("Unbalanced native traceback");
        for (int i = 0; i < n_; ++i)
            if (find(i) != find(0) || (required_[i] >= 0 && required_[i] != partners[i]))
                throw std::logic_error("Native traceback violated connectivity or constraints");
    }
public:
    DP(const aptadp_input_v1& input, const aptadp_model_v1& model,
       std::unique_ptr<EnergyProvider> provider = nullptr) {
        if (input.struct_size != sizeof(input) || model.struct_size != sizeof(model))
            throw std::invalid_argument("C ABI struct_size mismatch");
        if (!input.length || input.length > static_cast<uint32_t>(std::numeric_limits<int>::max() - 2))
            throw std::invalid_argument("Invalid native nucleotide count");
        if (!input.sequence || !input.materials || !input.nicks_before || !input.pair_allowed)
            throw std::invalid_argument("Required input array is NULL");
        if (input.min_hairpin_size < 0 || input.max_interior_size < 0 || input.max_interior_size > 30)
            throw std::invalid_argument("Invalid hairpin or maximum interior size (must be <=30)");
        if (!std::isfinite(model.thermal_energy_kcal_mol) || model.thermal_energy_kcal_mol <= 0
            || (model.log_global_factor != NEG && !std::isfinite(model.log_global_factor)))
            throw std::invalid_argument("Invalid thermal energy or log_global_factor");
        if (!provider && model.fixture_provider != APTADP_FIXTURE_UNIT && model.fixture_provider != APTADP_FIXTURE_ASYMMETRIC)
            throw std::invalid_argument("This build supports explicitly synthetic fixture providers only");
        n_ = static_cast<int>(input.length); min_hairpin_ = input.min_hairpin_size; max_interior_ = input.max_interior_size;
        thermal_ = model.thermal_energy_kcal_mol; global_ = model.log_global_factor;
        const uint64_t vertices = static_cast<uint64_t>(n_) + 2;
        const uint64_t triangular = mul(vertices, vertices - 1) / 2;
        const uint64_t cells = mul(triangular, NF);
        memory_ = add(add(mul(cells, 3 * sizeof(double)), mul(mul(n_, n_), 9)), mul(n_ + 1, 256));
        const uint64_t limit = input.memory_limit_bytes ? input.memory_limit_bytes : 512ull * 1024 * 1024;
        memory_limit_ = limit;
        if (memory_ > limit || cells > std::numeric_limits<size_t>::max() / sizeof(double))
            throw std::invalid_argument("Native workspace exceeds the explicit memory limit");
        triangle_ = static_cast<size_t>(triangular); cells_ = static_cast<size_t>(cells);
        sequence_.assign(input.sequence, input.sequence + n_);
        materials_.assign(input.materials, input.materials + n_);
        cuts_.assign(input.nicks_before, input.nicks_before + n_ + 1);
        if (cuts_[0] || cuts_[n_]) throw std::invalid_argument("Endpoint cuts must be zero; root nick is implicit");
        required_.resize(n_);
        cut_prefix_.assign(n_ + 1, 0); blocked_prefix_.assign(n_ + 1, 0); dna_prefix_.assign(n_ + 1, 0);
        for (int i = 0; i < n_; ++i) {
            if (materials_[i] != 1 && materials_[i] != 2) throw std::invalid_argument("Material must be DNA=1 or RNA=2");
            if (cuts_[i] > 1 || cuts_[i + 1] > 1) throw std::invalid_argument("Cuts must be binary");
            if (i && !cuts_[i] && materials_[i] != materials_[i - 1])
                throw std::invalid_argument("A strand must have one explicit material");
            const std::string alphabet = materials_[i] == 1 ? "ACGT" : "ACGU";
            if (alphabet.find(sequence_[i]) == std::string::npos)
                throw std::invalid_argument("Sequence alphabet disagrees with explicit material");
            required_[i] = input.constraints ? input.constraints[i] : -1;
            if (required_[i] < -1 || required_[i] >= n_) throw std::invalid_argument("Constraint index is out of range");
            cut_prefix_[i + 1] = cut_prefix_[i] + cuts_[i + 1];
            blocked_prefix_[i + 1] = blocked_prefix_[i] + (required_[i] != -1 && required_[i] != i);
            dna_prefix_[i + 1] = dna_prefix_[i] + (materials_[i] == 1);
        }
        for (int i = 0; i < n_; ++i) {
            const int j = required_[i];
            if (j >= 0 && required_[j] != i) throw std::invalid_argument("Constraints are not reciprocal");
            if (j > i) for (int k = i + 1; k < j; ++k)
                if (required_[k] > j) throw std::invalid_argument("Pseudoknotted constraints are unsupported");
        }
        allowed_.assign(static_cast<size_t>(n_) * n_, 0);
        for (int i = 0; i < n_; ++i) for (int j = 0; j < n_; ++j) {
            const size_t at = static_cast<size_t>(i) * n_ + j;
            const uint8_t permit = input.pair_allowed[at];
            if (permit > 1 || permit != input.pair_allowed[static_cast<size_t>(j) * n_ + i] || (i == j && permit))
                throw std::invalid_argument("Pair mask must be symmetric binary with zero diagonal");
            if (j > i && permit && (required_[i] == -1 || required_[i] == j)
                && (required_[j] == -1 || required_[j] == i)
                && (nicks(i, j) || j - i - 1 >= min_hairpin_)) allowed_[at] = 1;
        }
        if (provider) provider_ = std::move(provider);
        else provider_ = std::make_unique<FixtureProvider>(model.fixture_provider == APTADP_FIXTURE_UNIT, thermal_, materials_);
        for (FaceMode mode : {FaceMode::Multiloop, FaceMode::Exterior}) for (int family = 1; family <= 3; ++family) {
            const int u = ui(mode, family);
            unpaired_[u].resize(n_); free_prefix_[u].assign(n_ + 1, 0); zero_prefix_[u].assign(n_ + 1, 0);
            for (int i = 0; i < n_; ++i) {
                const double weight = logweight(provider_->log_unpaired(mode, family, i));
                unpaired_[u][i] = weight;
                free_prefix_[u][i + 1] = free_prefix_[u][i] + (weight == NEG ? 0 : weight);
                zero_prefix_[u][i + 1] = zero_prefix_[u][i] + (weight == NEG);
            }
        }
    }
    int length() const { return n_; }
    uint64_t memory() const { return memory_; }
    double log_partition() {
        ensure_inside();
        const double raw = evaluate({ROOT, -1, n_}, inside_, false);
        if (raw == NEG || global_ == NEG) return NEG;
        const double result = raw + global_;
        if (!std::isfinite(result)) throw std::range_error("Global log partition exceeds floating-point range");
        return result;
    }
    void pairs(double* matrix, double& adjustment) {
        if (log_partition() == NEG) throw Empty("The connected ensemble is empty");
        const double log_z = evaluate({ROOT, -1, n_}, inside_, false);
        const auto outside = outside_values(inside_, false);
        std::fill(matrix, matrix + static_cast<size_t>(n_) * n_, 0.0); adjustment = 0;
        for (int i = 0; i < n_; ++i) for (int j = i + 1; j < n_; ++j) {
            const Cell cell{BF, i, j};
            if (value(inside_, cell) == NEG || value(outside, cell) == NEG) continue;
            const double p = std::exp(value(inside_, cell) + value(outside, cell) - log_z);
            if (!std::isfinite(p) || p > 1 + 1e-12) throw std::range_error("Outside probability exceeds one");
            const double corrected = std::min(1.0, p);
            adjustment = std::max(adjustment, std::abs(p - corrected));
            matrix[static_cast<size_t>(i) * n_ + j] = matrix[static_cast<size_t>(j) * n_ + i] = corrected;
        }
        for (int i = 0; i < n_; ++i) {
            long double sum = 0;
            for (int j = 0; j < n_; ++j) sum += matrix[static_cast<size_t>(i) * n_ + j];
            const double unpaired = static_cast<double>(1 - sum);
            if (unpaired < -1e-12) throw std::range_error("Outside row sums to more than one");
            matrix[static_cast<size_t>(i) * n_ + i] = std::max(0.0, unpaired);
            adjustment = std::max(adjustment, std::max(0.0, -unpaired));
        }
    }
    std::pair<double, std::string> minimum() {
        std::vector<double> best(cells_, NEG); fill(best, true);
        const double weight = evaluate({ROOT, -1, n_}, best, true);
        if (weight == NEG || global_ == NEG) throw Empty("The connected ensemble is empty");
        const double energy = -thermal_ * (weight + global_);
        if (!std::isfinite(energy)) throw std::range_error("Minimum total energy exceeds floating-point range");
        return {energy, traceback(best, nullptr)};
    }
    double minimum_pair_energies(double* matrix) const {
        std::vector<double> best(cells_, NEG); fill(best, true);
        const double weight = evaluate({ROOT, -1, n_}, best, true);
        if (weight == NEG || global_ == NEG) throw Empty("The connected ensemble is empty");
        const double minimum = -thermal_ * (weight + global_);
        if (!std::isfinite(minimum)) throw std::range_error("Minimum total energy exceeds floating-point range");
        const auto outside = outside_values(best, true);
        std::fill(matrix, matrix + static_cast<size_t>(n_) * n_, std::numeric_limits<double>::infinity());
        for (int i = 0; i < n_; ++i) for (int j = i + 1; j < n_; ++j) {
            const Cell cell{BF, i, j};
            if (value(best, cell) == NEG || value(outside, cell) == NEG) continue;
            const double paired_weight = value(best, cell) + value(outside, cell);
            if (!std::isfinite(paired_weight)) throw std::range_error("Conditioned maximum log weight exceeds floating-point range");
            // Only permit accumulation roundoff, measured before adding the
            // potentially large global constant. Never clamp energies or gaps.
            const double tolerance = 64 * std::numeric_limits<double>::epsilon()
                * std::max({1.0, std::abs(weight), std::abs(paired_weight)});
            if (paired_weight > weight + tolerance)
                throw std::range_error("Pair-conditioned minimum is below the ensemble minimum");
            const double energy = -thermal_ * (paired_weight + global_);
            if (!std::isfinite(energy)) throw std::range_error("Pair-conditioned energy exceeds floating-point range");
            matrix[static_cast<size_t>(i) * n_ + j] = matrix[static_cast<size_t>(j) * n_ + i] = energy;
        }
        return minimum;
    }
    double minimum_structures(double atol, uint64_t limit, char* output, double* energies,
                              uint64_t& count, uint8_t& truncated) const {
        if (!std::isfinite(atol) || atol < 0 || !std::isfinite(atol / thermal_))
            throw std::invalid_argument("Minimum-structure tolerance and atol/RT must be finite and nonnegative");
        std::vector<double> best(cells_, NEG); fill(best, true);
        const double weight = evaluate({ROOT, -1, n_}, best, true);
        if (weight == NEG || global_ == NEG) throw Empty("The connected ensemble is empty");
        const double minimum = -thermal_ * (weight + global_);
        if (!std::isfinite(minimum)) throw std::range_error("Minimum total energy exceeds floating-point range");
        count = 0; truncated = 0;
        if (!limit) { truncated = 1; return minimum; }
        const std::string first = traceback(best, nullptr);
        std::memcpy(output, first.c_str(), n_ + 1); energies[0] = minimum; count = 1;

        constexpr size_t NONE = std::numeric_limits<size_t>::max();
        struct Work { Cell cell; size_t next; };
        struct Frame {
            Cell cell;
            size_t rest, pool_size;
            long double prefix, remaining;
            std::vector<Production> choices;
            size_t next = 0, active = 0;
        };
        // Persistent pending-list indices let backtracking restore an earlier
        // continuation without copying whole strings/stacks per alternative.
        // No C++ call-stack recursion is used, including final rescoring.
        const size_t maximum_depth = static_cast<size_t>(n_) * 2 + 3;
        std::vector<Frame> frames; frames.reserve(maximum_depth);
        std::vector<Work> pool; pool.reserve(maximum_depth);
        std::vector<double> evaluated; evaluated.reserve(maximum_depth);
        pool.push_back({{ROOT, -1, n_}, NONE});
        std::string structure(static_cast<size_t>(n_), '.');
        size_t head = 0;
        long double prefix = 0, remaining = weight;
        uint64_t choice_bytes = 0;
        const uint64_t fixed_bytes = add(memory_ - mul(cells_, sizeof(double)),
            add(add(mul(frames.capacity(), sizeof(Frame)), mul(pool.capacity(), sizeof(Work))),
                add(mul(evaluated.capacity(), sizeof(double)), mul(n_ + 1, 2))));
        if (fixed_bytes > memory_limit_)
            throw std::invalid_argument("Minimum enumeration workspace exceeds the explicit memory limit");
        const long double allowed_gap = static_cast<long double>(atol) / thermal_;
        // This guard relaxes ONLY pruning. Every completed derivation is
        // rescored in the original double DP order and tested strictly below.
        const long double pruning_guard = 64 * std::numeric_limits<double>::epsilon()
            * (n_ + 1) * std::max(1.0, std::abs(weight));
        auto cell_best = [&](Cell cell) { return cell.f == ROOT ? weight : value(best, cell); };
        for (;;) {
            if (head != NONE) {
                if (frames.size() >= maximum_depth) throw std::logic_error("Enumeration depth invariant failed");
                const auto work = pool[head];
                Frame frame{work.cell, work.next, pool.size(), prefix, remaining, {}};
                productions(work.cell, best, [&](const Production& p) {
                    if (score(p, best) != NEG) frame.choices.push_back(p);
                });
                const uint64_t bytes = mul(frame.choices.capacity(), sizeof(Production));
                if (add(fixed_bytes, add(choice_bytes, bytes)) > memory_limit_)
                    throw std::invalid_argument("Minimum enumeration frontier exceeds the explicit memory limit");
                choice_bytes += bytes;
                frames.push_back(std::move(frame));
                if (work.cell.f == BF) { structure[work.cell.a] = '('; structure[work.cell.b] = ')'; }
            } else {
                // Frames are a preorder derivation. Reverse-order reduction
                // reproduces score(p,best)'s child order exactly; a flat sum
                // could spuriously lose or gain tied minima through rounding.
                evaluated.clear();
                for (auto frame = frames.rbegin(); frame != frames.rend(); ++frame) {
                    const auto& p = frame->choices[frame->active];
                    double total = p.coefficient;
                    for (int k = 0; k < p.count; ++k) {
                        if (evaluated.empty()) throw std::logic_error("Incomplete enumeration derivation");
                        total += evaluated.back(); evaluated.pop_back();
                    }
                    if (!std::isfinite(total)) throw std::range_error("Enumerated log weight exceeds floating-point range");
                    evaluated.push_back(total);
                }
                if (evaluated.size() != 1) throw std::logic_error("Invalid enumeration reduction");
                const double actual_weight = evaluated.back();
                if (actual_weight > weight) throw std::range_error("Enumerated energy is below the DP minimum");
                const long double gap = static_cast<long double>(weight) - actual_weight;
                if (gap <= allowed_gap && structure != first) {
                    validate_traceback(structure);
                    if (count == limit) { truncated = 1; return minimum; }
                    const double energy = -thermal_ * (actual_weight + global_);
                    if (!std::isfinite(energy)) throw std::range_error("Enumerated energy exceeds floating-point range");
                    std::memcpy(output + static_cast<size_t>(count) * (n_ + 1), structure.c_str(), n_ + 1);
                    energies[count++] = energy;
                }
            }
            bool advanced = false;
            while (!frames.empty()) {
                auto& frame = frames.back();
                pool.resize(frame.pool_size);
                while (frame.next < frame.choices.size()) {
                    const auto& p = frame.choices[frame.next++];
                    const long double next_prefix = frame.prefix + p.coefficient;
                    long double next_remaining = frame.remaining - cell_best(frame.cell);
                    for (int k = 0; k < p.count; ++k) next_remaining += value(best, p.children[k]);
                    if (next_prefix + next_remaining < static_cast<long double>(weight) - allowed_gap - pruning_guard)
                        continue;
                    frame.active = frame.next - 1;
                    prefix = next_prefix; remaining = next_remaining; head = frame.rest;
                    for (int k = p.count - 1; k >= 0; --k) {
                        if (pool.size() >= maximum_depth) throw std::logic_error("Enumeration continuation invariant failed");
                        pool.push_back({p.children[k], head}); head = pool.size() - 1;
                    }
                    advanced = true; break;
                }
                if (advanced) break;
                if (frame.cell.f == BF) { structure[frame.cell.a] = '.'; structure[frame.cell.b] = '.'; }
                choice_bytes -= mul(frame.choices.capacity(), sizeof(Production));
                frames.pop_back();
            }
            if (!advanced) return minimum;
        }
    }
    void sample(uint64_t seed, uint64_t count, char* output) {
        if (!count) return;
        if (log_partition() == NEG) throw Empty("The connected ensemble is empty");
        std::mt19937_64 rng(seed);
        for (uint64_t k = 0; k < count; ++k) {
            const std::string result = traceback(inside_, &rng);
            std::memcpy(output + static_cast<size_t>(k) * (n_ + 1), result.c_str(), n_ + 1);
        }
    }
};

void error_text(char* output, size_t capacity, const char* message) {
    if (output && capacity) {
        const size_t count = std::min(capacity - 1, std::strlen(message));
        std::memcpy(output, message, count); output[count] = '\0';
    }
}
template<class Work> int boundary(char* error, size_t capacity, Work&& work) noexcept {
    try { work(); error_text(error, capacity, ""); return APTADP_OK; }
    catch (const Empty& e) { error_text(error, capacity, e.what()); return APTADP_EMPTY_ENSEMBLE; }
    catch (const std::bad_alloc& e) { error_text(error, capacity, e.what()); return APTADP_ALLOCATION_ERROR; }
    catch (const std::invalid_argument& e) { error_text(error, capacity, e.what()); return APTADP_BAD_INPUT; }
    catch (const std::range_error& e) { error_text(error, capacity, e.what()); return APTADP_NUMERICAL_ERROR; }
    catch (const std::exception& e) { error_text(error, capacity, e.what()); return APTADP_INTERNAL_ERROR; }
    catch (...) { error_text(error, capacity, "Unknown native exception"); return APTADP_INTERNAL_ERROR; }
}
}

struct aptadp_handle { std::unique_ptr<aptaswitch::thermodynamics::Calculation> dp; std::mutex mutex; };
namespace aptaswitch { namespace thermodynamics {
aptadp_handle* make_calculation_handle_v1(std::unique_ptr<Calculation> calculation) {
    if (!calculation) throw std::invalid_argument("A native calculation is required");
    auto handle = std::make_unique<aptadp_handle>();
    handle->dp = std::move(calculation);
    return handle.release();
}
aptadp_handle* make_handle_v1(const aptadp_input_v1& input, double thermal, double global,
                            std::unique_ptr<EnergyProvider> provider) {
    if (!provider) throw std::invalid_argument("A native energy provider is required");
    const aptadp_model_v1 model{sizeof(aptadp_model_v1), 0, thermal, global};
    auto handle = std::make_unique<aptadp_handle>();
    handle->dp = std::make_unique<DP>(input, model, std::move(provider));
    return handle.release();
}
}}
extern "C" {
uint32_t aptadp_abi_version(void) { return APTADP_ABI_VERSION; }
int aptadp_create_v1(const aptadp_input_v1* input, const aptadp_model_v1* model, aptadp_handle** output,
                     char* error, size_t capacity) {
    if (output) *output = nullptr;
    return boundary(error, capacity, [&] {
        if (!input || !model || !output) throw std::invalid_argument("NULL constructor argument");
        auto handle = std::make_unique<aptadp_handle>();
        handle->dp = std::make_unique<DP>(*input, *model); *output = handle.release();
    });
}
void aptadp_destroy(aptadp_handle* handle) { delete handle; }
int aptadp_log_partition(aptadp_handle* handle, double* output, char* error, size_t capacity) {
    return boundary(error, capacity, [&] {
        if (!handle || !output) throw std::invalid_argument("NULL partition argument");
        std::lock_guard<std::mutex> lock(handle->mutex);
        *output = handle->dp->log_partition();
    });
}
int aptadp_pair_probabilities(aptadp_handle* handle, double* matrix, size_t cells, double* adjustment,
                              char* error, size_t capacity) {
    return boundary(error, capacity, [&] {
        if (!handle || !matrix || !adjustment) throw std::invalid_argument("NULL probability argument");
        std::lock_guard<std::mutex> lock(handle->mutex);
        const size_t n = static_cast<size_t>(handle->dp->length());
        if (cells < n * n) throw std::invalid_argument("Probability buffer is too small");
        handle->dp->pairs(matrix, *adjustment);
    });
}
int aptadp_minimum(aptadp_handle* handle, double* energy, char* structure, size_t size,
                   char* error, size_t capacity) {
    return boundary(error, capacity, [&] {
        if (!handle || !energy || !structure) throw std::invalid_argument("NULL minimum argument");
        std::lock_guard<std::mutex> lock(handle->mutex);
        if (size < static_cast<size_t>(handle->dp->length()) + 1) throw std::invalid_argument("Structure buffer is too small");
        const auto result = handle->dp->minimum(); *energy = result.first;
        std::memcpy(structure, result.second.c_str(), result.second.size() + 1);
    });
}
int aptadp_minimum_pair_energies(aptadp_handle* handle, double* matrix, size_t cells, double* minimum,
                                char* error, size_t capacity) {
    return boundary(error, capacity, [&] {
        if (!handle || !matrix || !minimum) throw std::invalid_argument("NULL conditioned-minimum argument");
        std::lock_guard<std::mutex> lock(handle->mutex);
        const size_t n = static_cast<size_t>(handle->dp->length());
        if (cells < n * n) throw std::invalid_argument("Pair-energy buffer is too small");
        *minimum = handle->dp->minimum_pair_energies(matrix);
    });
}
int aptadp_minimum_structures(aptadp_handle* handle, double atol, uint64_t limit,
                              char* structures, size_t bytes, double* energies, size_t energy_cells,
                              uint64_t* count, uint8_t* truncated, double* minimum, char* error, size_t capacity) {
    if (count) *count = 0;
    if (truncated) *truncated = 0;
    return boundary(error, capacity, [&] {
        if (!handle || !count || !truncated || !minimum || (limit && (!structures || !energies)))
            throw std::invalid_argument("NULL minimum-enumeration argument");
        std::lock_guard<std::mutex> lock(handle->mutex);
        if (mul(limit, static_cast<uint64_t>(handle->dp->length()) + 1) > bytes || limit > energy_cells)
            throw std::invalid_argument("Minimum-enumeration buffer is too small");
        *minimum = handle->dp->minimum_structures(atol, limit, structures, energies, *count, *truncated);
    });
}
int aptadp_sample(aptadp_handle* handle, uint64_t seed, uint64_t count, char* structures, size_t size,
                  char* error, size_t capacity) {
    return boundary(error, capacity, [&] {
        if (!handle || (count && !structures)) throw std::invalid_argument("NULL sample argument");
        std::lock_guard<std::mutex> lock(handle->mutex);
        const uint64_t bytes = mul(count, static_cast<uint64_t>(handle->dp->length()) + 1);
        if (bytes > size) throw std::invalid_argument("Sample buffer is too small");
        handle->dp->sample(seed, count, structures);
    });
}
int aptadp_memory_estimate(aptadp_handle* handle, uint64_t* bytes, char* error, size_t capacity) {
    return boundary(error, capacity, [&] {
        if (!handle || !bytes) throw std::invalid_argument("NULL memory-estimate argument");
        std::lock_guard<std::mutex> lock(handle->mutex);
        *bytes = handle->dp->memory();
    });
}
}
