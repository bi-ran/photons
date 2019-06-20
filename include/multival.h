#ifndef MULTIVAL_H
#define MULTIVAL_H

#include <iterator>
#include <vector>

#include "interval.h"

class multival {
  public:
    template <typename... T>
    multival(T const&... intervals)
            : _dims(sizeof...(T)) {
        extract(intervals...);
        _size = std::accumulate(std::begin(_shape), std::end(_shape),
                                1, std::multiplies<int64_t>());
    }

    multival(multival const& other) = default;
    multival& operator=(multival const& other) = default;
    ~multival() = default;

    template <template <typename...> class T>
    std::vector<int64_t> indices_for(T<double> const& values) const {
        std::vector<double> x(std::begin(values), std::end(values));
        std::vector<int64_t> indices(_dims);
        for (int64_t i = 0; i < _dims; ++i)
            indices[i] = _intervals[i].index_for(x[i]);

        return indices;
    }

    template <template <typename...> class T, typename U>
    typename std::enable_if<std::is_integral<U>::value, int64_t>::type
    index_for(T<U> const& indices) const {
        std::vector<int64_t> x(std::begin(indices), std::end(indices));
        int64_t index = 0;
        for (int64_t i = 0, block = 1; i < _dims; ++i) {
            index = index + x[i] * block;
            block = block * _shape[i];
        }

        return index;
    }

    template <template <typename...> class T, typename U>
    typename std::enable_if<!std::is_integral<U>::value, int64_t>::type
    index_for(T<U> const& values) const {
        return index_for(indices_for(values)); }

    std::vector<int64_t> const& shape() const { return _shape; }
    int64_t dims() const { return _dims; }
    int64_t size() const { return _size; }

  private:
    template <typename... T>
    void extract(T const&... args) {
        (void) (int [sizeof...(T)]) { (_intervals.emplace_back(args), 0)... };
        for (auto const& i : _intervals) { _shape.push_back(i.size()); }
    }

    std::vector<int64_t> _shape;
    int64_t _dims;
    int64_t _size;

    std::vector<interval> _intervals;
};

#endif /* MULTIVAL_H */
