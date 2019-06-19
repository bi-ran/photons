#ifndef INTEGRAL_ANGLES_H
#define INTEGRAL_ANGLES_H

#include <limits>

constexpr int64_t convert_pis(long double angle) {
    return static_cast<int64_t>(angle * std::numeric_limits<int64_t>::max());
}

constexpr int64_t convert_radians(long double angle) {
    return static_cast<int64_t>((angle / M_PI)
        * std::numeric_limits<int64_t>::max());
}

constexpr int64_t convert_degrees(long double angle) {
    return static_cast<int64_t>((angle / 180.)
        * std::numeric_limits<int64_t>::max());
}

constexpr int64_t operator ""_pi(long double angle) {
    return convert_pis(angle);
}

constexpr int64_t operator ""_rad(long double angle) {
    return convert_radians(angle);
}

constexpr int64_t operator ""_deg(long double angle) {
    return convert_degrees(angle);
}

double revert_radians(int64_t angle) {
    return static_cast<double>(angle) / std::numeric_limits<int64_t>::max();
}

#endif /* INTEGRAL_ANGLES_H */
