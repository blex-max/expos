#include <functional>
#include <vector>
#include <concepts>
#include <cstdint>
#include <format>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <stdexcept>
#include <string>

template <std::signed_integral T>
uint64_t
constexpr inline as_uint (const T &a)
{
    if (a < 0) {
        throw std::runtime_error (
            "cannot convert negative value to uint"
        );
    }
    return static_cast<uint64_t> (a);
}

template <class T>
std::string opt_to_str (
    std::optional<T> opt,
    std::string_view sentinel,
    std::function<std::string (T)> conv = [] (const T &a) {
        return std::to_string (a);
    }
) {
    return opt ? conv (*opt) : std::string (sentinel);
}

inline std::string rdbl2 (const double &a) {
    return std::format ("{:.2f}", a);
}
inline std::string rdbl4 (const double &a) {
    return std::format ("{:.4f}", a);
}


