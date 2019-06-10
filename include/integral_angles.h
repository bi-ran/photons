#ifndef INTEGRAL_ANGLES_H
#define INTEGRAL_ANGLES_H

int64_t convert_radians(long double angle) {
    return (int64_t)(angle * std::numeric_limits<int64_t>::max() / M_PI);
}

int64_t convert_degrees(long double angle) {
    return (int64_t)(angle * std::numeric_limits<int64_t>::max() / 180.);
}

constexpr int64_t operator ""_pi(long double angle) {
    return (int64_t)(angle * std::numeric_limits<int64_t>::max());
}

constexpr int64_t operator ""_rad(long double angle) {
    return (int64_t)(angle * std::numeric_limits<int64_t>::max() / M_PI);
}

constexpr int64_t operator ""_deg(long double angle) {
    return (int64_t)(angle * std::numeric_limits<int64_t>::max() / 180.);
}

#endif /* INTEGRAL_ANGLES_H */
