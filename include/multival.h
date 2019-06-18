#ifndef MULTIVAL_H
#define MULTIVAL_H

#include <iterator>
#include <vector>

#include "interval.h"

class multival {
  public:
    template <typename... T>
    multival(T const&... intervals) : _dims(sizeof...(T)) {
        extract(intervals...);
        _size = std::accumulate(std::begin(_shape), std::end(_shape),
                                1, std::multiplies<int64_t>());
    }

    multival(multival const& other) = default;

    multival& operator=(multival const& other) = default;

    ~multival() = default;

    std::vector<int64_t> indices_for(std::vector<double> const& values) const {
        std::vector<int64_t> indices(_shape);
        for (int64_t i = 0; i < _dims; ++i) {
            indices[i] = _intervals[i].index_for(values[i]);
#ifdef DEBUG
            if (indices[i] < 0 || indices[i] == _shape[i]) {
                printf("value: %f [%i] exceeds boundaries!\n", values[i], i);
                exit(1);
            }
#endif /* DEBUG */
        }

        return indices;
    }

    int64_t dims() const { return _dims; }
    int64_t size() const { return _size; }
    std::vector<int64_t> const& shape() const { return _shape; }

  private:
    template <typename... T>
    void extract(T const&... args) {
        (void) (int [sizeof...(T)]) { (_intervals.emplace_back(args), 0)... };
        for (auto const& i : _intervals) { _shape.push_back(i.size()); }
    }

    int64_t _dims;
    int64_t _size;
    std::vector<int64_t> _shape;

    std::vector<interval> _intervals;
};

#endif /* MULTIVAL_H */
