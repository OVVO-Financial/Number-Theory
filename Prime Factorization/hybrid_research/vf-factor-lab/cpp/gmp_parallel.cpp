#include "gmp_parallel.hpp"
#include <charconv>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace {

std::string quoted(const std::string& value) {
    std::ostringstream out; out << '"';
    for (unsigned char c : value) {
        if (c=='"') out << "\\\"";
        else if (c=='\\') out << "\\\\";
        else if (c=='\n') out << "\\n";
        else if (c=='\r') out << "\\r";
        else if (c=='\t') out << "\\t";
        else if (c<32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << unsigned(c) << std::dec;
        else out << char(c);
    }
    out << '"'; return out.str();
}
std::string integer(const vf::Z& value) { return quoted(value.get_str()); }
std::string optional_integer(const std::optional<vf::Z>& value) { return value ? integer(*value) : "null"; }

void write_result(std::ostream& out,const vf::ParallelResult& r,bool trace) {
    out << std::setprecision(17) << "{\"algorithm\":\"regional_gmp_parallel\",\"n\":" << integer(r.n)
        << ",\"factor\":" << (r.factor==0 ? "null" : integer(r.factor))
        << ",\"cofactor\":" << (r.factor==0 ? "null" : integer(vf::Z(r.n/r.factor)))
        << ",\"status\":" << quoted(r.status) << ",\"elapsed_seconds\":" << r.elapsed_seconds
        << ",\"process_cpu_seconds\":" << r.process_cpu_seconds << ",\"counters\":{";
    bool first=true;
    for (const auto& [key,value] : r.counters.counts) {
        if (!first) out << ',';
        first=false; out << quoted(key) << ':' << value;
    }
    for (const auto& [key,value] : r.counters.seconds) {
        if (!first) out << ',';
        first=false; out << quoted(key) << ':' << value;
    }
    out << "},\"details\":{\"threads\":" << r.threads << ",\"phase\":" << quoted(r.phase)
        << ",\"reason\":" << quoted(r.reason)
        << ",\"search_lo\":" << integer(r.search_lo) << ",\"search_hi\":" << integer(r.search_hi)
        << ",\"lo\":" << optional_integer(r.lo) << ",\"hi\":" << optional_integer(r.hi)
        << ",\"outstanding_reservations\":" << r.outstanding_reservations
        << ",\"seed_work_units\":" << r.seed_work_units
        << ",\"coverage_intervals\":[";
    for (std::size_t i=0;i<r.certificates.size();++i) {
        if (i) out << ',';
        out << '[' << integer(r.certificates[i].lo) << ',' << integer(r.certificates[i].hi) << ']';
    }
    out << ']';
    if (trace) {
        out << ",\"trace\":[";
        for (std::size_t i=0;i<r.trace.size();++i) {
            if (i) out << ',';
            const auto& t=r.trace[i];
            out << "{\"event\":" << quoted(t.event) << ",\"kind\":" << quoted(t.kind)
                << ",\"status\":" << quoted(t.status) << ",\"job_id\":" << t.job_id
                << ",\"worker_id\":" << t.worker_id << ",\"job_work\":" << t.job_work
                << ",\"work_units\":" << t.work_units
                << ",\"requested_lo\":" << integer(t.requested_lo)
                << ",\"requested_hi\":" << integer(t.requested_hi)
                << ",\"certified_lo\":" << optional_integer(t.certified_lo)
                << ",\"certified_hi\":" << optional_integer(t.certified_hi)
                << ",\"lo\":" << optional_integer(t.lo) << ",\"hi\":" << optional_integer(t.hi) << '}';
        }
        out << ']';
    }
    out << "}}";
}

void write_integers(std::ostream& out,const std::vector<vf::Z>& values) {
    out << '[';
    for (std::size_t i=0;i<values.size();++i) { if (i) out << ','; out << integer(values[i]); }
    out << ']';
}

void write_full(std::ostream& out,const vf::ParallelFullResult& r,bool trace) {
    out << std::setprecision(17) << "{\"input\":" << integer(r.input)
        << ",\"status\":" << quoted(!r.unresolved.empty() ? "partial" : !r.probable_prime_factors.empty() ? "probable_complete" : "complete")
        << ",\"elapsed_seconds\":" << r.elapsed_seconds << ",\"process_cpu_seconds\":" << r.process_cpu_seconds
        << ",\"work_units\":" << r.work_units << ",\"prime_factors\":";
    write_integers(out,r.prime_factors);
    out << ",\"probable_prime_factors\":"; write_integers(out,r.probable_prime_factors);
    out << ",\"unresolved\":"; write_integers(out,r.unresolved);
    out << ",\"runs\":[";
    for (std::size_t i=0;i<r.runs.size();++i) { if (i) out << ','; write_result(out,r.runs[i],trace); }
    out << "]}";
}

vf::Z parse_integer(const std::string& text) {
    const std::size_t start=!text.empty() && text[0]=='+' ? 1 : 0;
    if (start==text.size()) throw std::invalid_argument("expected decimal integer >= 2");
    for (std::size_t i=start;i<text.size();++i)
        if (text[i]<'0' || text[i]>'9') throw std::invalid_argument("expected decimal integer >= 2");
    vf::Z n;
    if (mpz_set_str(n.get_mpz_t(),text.c_str()+start,10) || n<2)
        throw std::invalid_argument("expected decimal integer >= 2");
    return n;
}

std::uint64_t parse_u64(const std::string& text,const std::string& flag) {
    std::uint64_t value=0;
    const auto parsed=std::from_chars(text.data(),text.data()+text.size(),value);
    if (parsed.ec!=std::errc{} || parsed.ptr!=text.data()+text.size())
        throw std::invalid_argument(flag+" requires an unsigned 64-bit integer");
    return value;
}

double parse_double(const std::string& text,const std::string& flag) {
    std::size_t end=0; double value;
    try { value=std::stod(text,&end); }
    catch (...) { throw std::invalid_argument(flag+" requires a finite number"); }
    if (end!=text.size() || !std::isfinite(value)) throw std::invalid_argument(flag+" requires a finite number");
    return value;
}

void metadata() {
    std::cout << "{\"program\":\"vf_factor_gmp_parallel\",\"version\":\"0.4.0\",\"gmp_version\":" << quoted(gmp_version)
              << ",\"compiler\":" << quoted(__VERSION__) << ",\"hardware_concurrency\":" << std::thread::hardware_concurrency()
              << ",\"parallelism\":\"C++17 worker pool with disjoint regional jobs\""
              << ",\"policy\":\"adaptive two-to-one Fermat/product job-time target including in-flight estimates; worker-local caches; exact shared interval ledger\""
              << ",\"default_degree\":128,\"default_fermat_span\":655360,\"default_job_seconds\":0.025"
              << ",\"network_transport\":false,\"explicit_factor_bands\":true,\"integer_json\":\"decimal strings\"}\n";
}

const char* help=R"HELP(Parallel regional factorization with C++17 and GMP

Usage: vf_factor_gmp_parallel N [options]
       vf_factor_gmp_parallel --batch [options] < inputs.txt

  --threads N                 Concurrent workers (default 4; 1..256)
  --timeout SECONDS           Overall cooperative deadline (default 10)
  --max-steps U64             Overall reserved/consumed work allowance
  --full                      Recursive factorization under the same budget
  --trace                     Job results and exact exclusion certificates
  --factor-lo INTEGER         Restrict smaller-factor band (inclusive)
  --factor-hi INTEGER         Restrict smaller-factor band (inclusive)
  --degree N                  Product block degree (default 128; 1..65536)
  --fermat-span U64           Raw midpoint positions per tile (default 655360)
  --job-seconds SECONDS        Maximum time allocation per job (default .025)
  --fermat-ratio NUMBER        Target Fermat/product job-time ratio (default 2)
  --initial-fermat-tests U64   Serial cheap close-factor pass (default 64)
  --no-vf                     Disable auxiliary VF probes
  --no-fermat                 Use divisor-product regions only
  --no-sharing                Separate polynomial constructions by step
  --gmp-fermat                Force GMP cursor below 2^64 too
  --no-probable-prime-check   Continue exact search on unchecked large primes
  --metadata                  Compiler, GMP and regional policy
  --help                      This help

Factor bands can be assigned manually to separate machines. This executable
implements local parallelism and explicit range jobs, not network transport.
A certified partial band returns region_exhausted, never a global prime claim.
--full cannot be combined with restricted factor bands. Probable primes remain
explicitly labeled. Arithmetic values are exact decimal strings in JSON.
)HELP";

} // anonymous namespace

int main(int argc,char** argv) {
    try {
        vf::ParallelOptions options;
        std::string input;
        bool batch=false,full=false,show_metadata=false;
        for (int i=1;i<argc;++i) {
            const std::string arg=argv[i];
            auto value=[&]() -> std::string {
                if (i+1>=argc) throw std::invalid_argument("missing value for "+arg);
                return argv[++i];
            };
            if (arg=="--help" || arg=="-h") { std::cout << help; return 0; }
            else if (arg=="--metadata" || arg=="--version") show_metadata=true;
            else if (arg=="--batch") batch=true;
            else if (arg=="--full") full=true;
            else if (arg=="--trace") options.algorithm.trace=true;
            else if (arg=="--no-vf") options.algorithm.vf_probes=false;
            else if (arg=="--no-fermat") options.algorithm.fermat_enabled=false;
            else if (arg=="--no-sharing") options.algorithm.shared_polynomial=false;
            else if (arg=="--gmp-fermat") options.algorithm.native64=false;
            else if (arg=="--no-probable-prime-check") options.algorithm.probable_prime_check=false;
            else if (arg=="--threads") options.threads=parse_u64(value(),arg);
            else if (arg=="--timeout") options.algorithm.timeout=parse_double(value(),arg);
            else if (arg=="--max-steps") options.algorithm.max_steps=parse_u64(value(),arg);
            else if (arg=="--degree") options.degree=parse_u64(value(),arg);
            else if (arg=="--fermat-span") options.fermat_span=parse_u64(value(),arg);
            else if (arg=="--job-seconds") options.job_seconds=parse_double(value(),arg);
            else if (arg=="--fermat-ratio") options.algorithm.fermat_ratio=parse_double(value(),arg);
            else if (arg=="--initial-fermat-tests") options.algorithm.initial_fermat_tests=parse_u64(value(),arg);
            else if (arg=="--factor-lo") options.factor_lo=parse_integer(value());
            else if (arg=="--factor-hi") options.factor_hi=parse_integer(value());
            else if (!arg.empty() && arg[0]=='-') throw std::invalid_argument("unknown option: "+arg);
            else if (input.empty()) input=arg;
            else throw std::invalid_argument("supply one N or use --batch");
        }
        options.validate();
        if (show_metadata) { metadata(); return 0; }
        if (full && (options.factor_lo || options.factor_hi))
            throw std::invalid_argument("--full cannot be combined with a restricted factor band");
        auto run=[&](const std::string& text) {
            const auto n=parse_integer(text);
            if (full) write_full(std::cout,vf::factor_parallel_full(n,options),options.algorithm.trace);
            else write_result(std::cout,vf::factor_parallel(n,options),options.algorithm.trace);
            std::cout << '\n';
        };
        if (batch) {
            if (!input.empty()) throw std::invalid_argument("--batch accepts inputs on stdin");
            std::string line; int code=0;
            while (std::getline(std::cin,line)) {
                if (!line.empty() && line.back()=='\r') line.pop_back();
                try { run(line); }
                catch (const std::exception& error) {
                    std::cout << "{\"status\":\"error\",\"error\":" << quoted(error.what()) << "}\n"; code=2;
                }
            }
            return code;
        }
        if (input.empty()) throw std::invalid_argument("supply N, --batch, or --help");
        run(input); return 0;
    } catch (const std::exception& error) {
        std::cerr << "{\"status\":\"error\",\"error\":" << quoted(error.what()) << "}\n"; return 2;
    }
}
