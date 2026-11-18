#include "gmp_factor.hpp"
#include <charconv>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string_view>

namespace {

std::string quote(const std::string& text) {
    std::ostringstream out; out << '"';
    for (unsigned char c : text) {
        switch (c) {
            case '"': out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << unsigned(c) << std::dec;
                else out << char(c);
        }
    }
    out << '"'; return out.str();
}

std::string integer(const vf::Z& value) { return quote(value.get_str()); }
const char* boolean(bool value) { return value ? "true" : "false"; }

void write_result(std::ostream& out, const vf::Result& r) {
    out << std::setprecision(17)
        << "{\"algorithm\":\"fused_gmp\",\"n\":" << integer(r.n)
        << ",\"factor\":" << (r.factor == 0 ? "null" : integer(r.factor))
        << ",\"cofactor\":" << (r.factor == 0 ? "null" : integer(vf::Z(r.n/r.factor)))
        << ",\"status\":" << quote(r.status)
        << ",\"elapsed_seconds\":" << r.elapsed_seconds << ",\"counters\":{";
    bool first = true;
    for (const auto& item : r.counters.counts) {
        if (!first) out << ',';
        first = false;
        out << quote(item.first) << ':' << item.second;
    }
    for (const auto& item : r.counters.seconds) {
        if (!first) out << ',';
        first = false;
        out << quote(item.first) << ':' << item.second;
    }
    for (const auto& item : r.counters.integers) {
        if (!first) out << ',';
        first = false;
        out << quote(item.first) << ':' << integer(item.second);
    }
    const auto& opt = r.options;
    out << "},\"details\":{\"phase\":" << quote(r.phase)
        << ",\"lo\":" << integer(r.lo) << ",\"hi\":" << integer(r.hi)
        << ",\"fermat_next\":" << (r.fermat_next ? integer(*r.fermat_next) : "null")
        << ",\"fermat_covered\":" << (r.fermat_covered ? integer(*r.fermat_covered) : "null")
        << ",\"fermat_backend\":" << quote(r.fermat_backend)
        << ",\"polynomial_backend\":\"gmp_kronecker_monic_tree\""
        << ",\"shared_bounds\":" << boolean(opt.feedback)
        << ",\"shared_polynomial\":" << boolean(opt.shared_polynomial)
        << ",\"vf_probes\":" << boolean(opt.vf_probes)
        << ",\"fermat_enabled\":" << boolean(opt.fermat_enabled)
        << ",\"probable_prime_check\":" << boolean(opt.probable_prime_check)
        << ",\"probable_prime_reps\":25"
        << ",\"initial_degree\":" << opt.initial_degree
        << ",\"max_degree\":" << opt.max_degree
        << ",\"fermat_ratio\":" << opt.fermat_ratio
        << ",\"fermat_period\":160,\"vf_levels\":[";
    for (std::size_t i=0; i<r.vf_levels.size(); ++i) {
        if (i) out << ',';
        out << r.vf_levels[i];
    }
    out << ']';
    if (!r.reason.empty()) out << ",\"reason\":" << quote(r.reason);
    if (opt.trace) {
        out << ",\"trace\":[";
        for (std::size_t i=0; i<r.trace.size(); ++i) {
            if (i) out << ',';
            out << "{\"event\":" << quote(r.trace[i].event);
            for (const auto& item : r.trace[i].fields)
                out << ',' << quote(item.first) << ':' << (item.second ? integer(*item.second) : "null");
            out << '}';
        }
        out << ']';
    }
    out << "}}";
}

void write_integers(std::ostream& out, const std::vector<vf::Z>& values) {
    out << '[';
    for (std::size_t i=0; i<values.size(); ++i) { if (i) out << ','; out << integer(values[i]); }
    out << ']';
}

void write_full(std::ostream& out, const vf::FullResult& r) {
    out << std::setprecision(17) << "{\"input\":" << integer(r.input)
        << ",\"status\":" << quote(!r.unresolved.empty() ? "partial" :
             !r.probable_prime_factors.empty() ? "probable_complete" : "complete")
        << ",\"prime_factors\":"; write_integers(out,r.prime_factors);
    out << ",\"probable_prime_factors\":"; write_integers(out,r.probable_prime_factors);
    out << ",\"unresolved\":"; write_integers(out,r.unresolved);
    out << ",\"elapsed_seconds\":" << r.elapsed_seconds << ",\"work_units\":" << r.work_units
        << ",\"runs\":[";
    for (std::size_t i=0; i<r.runs.size(); ++i) { if (i) out << ','; write_result(out,r.runs[i]); }
    out << "]}";
}

vf::Z parse_integer(const std::string& text) {
    std::size_t start = !text.empty() && text[0] == '+' ? 1 : 0;
    if (start == text.size()) throw std::invalid_argument("N must be a decimal integer >= 2");
    for (std::size_t i=start; i<text.size(); ++i)
        if (text[i] < '0' || text[i] > '9') throw std::invalid_argument("N must be a decimal integer >= 2");
    vf::Z value;
    if (mpz_set_str(value.get_mpz_t(),text.c_str()+start,10) != 0 || value < 2)
        throw std::invalid_argument("N must be a decimal integer >= 2");
    return value;
}

std::uint64_t parse_u64(const std::string& text, const std::string& name) {
    std::uint64_t value = 0;
    const auto parsed = std::from_chars(text.data(),text.data()+text.size(),value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data()+text.size())
        throw std::invalid_argument(name+" requires an unsigned 64-bit integer");
    return value;
}

double parse_double(const std::string& text, const std::string& name) {
    std::size_t end = 0;
    double value;
    try { value = std::stod(text,&end); }
    catch (...) { throw std::invalid_argument(name+" requires a finite number"); }
    if (end != text.size() || !std::isfinite(value))
        throw std::invalid_argument(name+" requires a finite number");
    return value;
}

void metadata(std::ostream& out) {
    out << "{\"program\":\"vf_factor_gmp\",\"version\":\"0.3.0\",\"gmp_version\":" << quote(gmp_version)
        << ",\"compiler\":" << quote(__VERSION__) << ",\"cplusplus\":" << __cplusplus
        << ",\"gmp_limb_bits\":" << GMP_LIMB_BITS
        << ",\"integer_backend\":\"GMP mpz_class\""
        << ",\"polynomial_backend\":\"GMP Kronecker multiplication and monic product/remainder trees\""
        << ",\"fermat_backends\":[\"native64\",\"gmp\"]"
        << ",\"fermat_period\":160,\"probable_prime_reps\":25"
        << ",\"big_integer_json\":\"decimal strings\",\"runtime_dependencies\":[\"libgmpxx\",\"libgmp\"]}"
        << '\n';
}

const char* help = R"HELP(Unified Fermat / VF / product factorization with GMP (C++17)

Usage: vf_factor_gmp N [options]
       vf_factor_gmp --batch [options] < decimal-inputs.txt
       vf_factor_gmp --metadata

  --timeout SECONDS          One cooperative deadline (default 10)
  --max-steps U64            Overall represented-work allowance
  --full                     Recursively split with one shared deadline/budget
  --trace                    Include exact shared-bound execution trace
  --gmp-fermat                Use GMP Fermat arithmetic even below 2^64
  --no-vf                    Disable auxiliary VF progression probes
  --no-feedback              Disable reciprocal bound feedback
  --no-sharing               Build separate block polynomials by step
  --no-fermat                Run the product side alone
  --initial-degree INT       Initial polynomial degree (default 16)
  --max-degree INT           Degree cap (default 512; maximum 65536)
  --fermat-ratio NUMBER       Time allocation ratio (default 2)
  --initial-fermat-tests U64  Initial admissible square tests (default 64)
  --no-probable-prime-check   Continue exact search on unchecked >64-bit primes
  --metadata                 Compiler, GMP and backend information
  --help                     Show this help

Factors, inputs and arbitrary-precision bounds are decimal JSON strings.
Above 64 bits, GMP probable primes are explicitly marked probable_prime;
--full keeps them separate from proven prime_factors. A budget result is
unresolved. Batch mode resets budgets per line and emits one JSON per line,
including an error object for an invalid line, then continues.
)HELP";

} // anonymous namespace

int main(int argc, char** argv) {
    try {
        vf::Options options;
        std::string input;
        bool batch = false, full = false, show_metadata = false;
        for (int i=1; i<argc; ++i) {
            const std::string arg = argv[i];
            auto value = [&]() -> std::string {
                if (i+1 >= argc) throw std::invalid_argument("missing value for "+arg);
                return argv[++i];
            };
            if (arg == "--help" || arg == "-h") { std::cout << help; return 0; }
            else if (arg == "--metadata" || arg == "--version") show_metadata = true;
            else if (arg == "--batch") batch = true;
            else if (arg == "--full") full = true;
            else if (arg == "--trace") options.trace = true;
            else if (arg == "--gmp-fermat") options.native64 = false;
            else if (arg == "--no-vf") options.vf_probes = false;
            else if (arg == "--no-feedback") options.feedback = false;
            else if (arg == "--no-sharing") options.shared_polynomial = false;
            else if (arg == "--no-fermat") options.fermat_enabled = false;
            else if (arg == "--no-probable-prime-check") options.probable_prime_check = false;
            else if (arg == "--timeout") options.timeout = parse_double(value(),arg);
            else if (arg == "--max-steps") options.max_steps = parse_u64(value(),arg);
            else if (arg == "--initial-degree") options.initial_degree = parse_u64(value(),arg);
            else if (arg == "--max-degree") options.max_degree = parse_u64(value(),arg);
            else if (arg == "--fermat-ratio") options.fermat_ratio = parse_double(value(),arg);
            else if (arg == "--initial-fermat-tests") options.initial_fermat_tests = parse_u64(value(),arg);
            else if (!arg.empty() && arg[0] == '-') throw std::invalid_argument("unknown option: "+arg);
            else if (input.empty()) input = arg;
            else throw std::invalid_argument("supply only one N, or use --batch");
        }
        options.validate();
        if (show_metadata) { metadata(std::cout); return 0; }
        auto run = [&](const std::string& text) {
            const auto n = parse_integer(text);
            if (full) write_full(std::cout,vf::factor_full(n,options));
            else write_result(std::cout,vf::factor(n,options));
            std::cout << '\n';
        };
        if (batch) {
            if (!input.empty()) throw std::invalid_argument("--batch takes N values from stdin, not an argument");
            std::string line;
            int exit_code = 0;
            while (std::getline(std::cin,line)) {
                if (!line.empty() && line.back() == '\r') line.pop_back();
                try { run(line); }
                catch (const std::exception& error) {
                    std::cout << "{\"status\":\"error\",\"error\":" << quote(error.what()) << "}\n";
                    exit_code = 2;
                }
            }
            return exit_code;
        }
        if (input.empty()) throw std::invalid_argument("supply N, --batch, or --help");
        run(input);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "{\"status\":\"error\",\"error\":" << quote(error.what()) << "}\n";
        return 2;
    }
}
