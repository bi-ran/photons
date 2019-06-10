#ifndef INTERVAL_H
#define INTERVAL_H

#include <vector>

class interval {
  public:
    template <template <typename...> class T>
    interval(T<float> const& edges)
        : _size(edges.size() - 1),
          _edges(std::begin(edges), std::end(edges)) { }

    interval(int64_t number, double min, double max)
            : _size(number),
              _edges(std::vector<double>(number + 1)) {
        double accumulator = min;
        double interval = (max - min) / number;
        for (auto& edge : _edges) {
            edge = accumulator;
            accumulator = accumulator + interval;
        }
    }

    template <typename... T>
    interval(T const&... edges)
        : _size(sizeof...(T) - 1),
          _edges({static_cast<double>(edges)...}) { }

    interval(interval const& other) = default;

    interval& operator=(interval const& other) = default;

    ~interval() = default;

    int64_t index_for(double value) const {
        int64_t index = _size;
        for (std::size_t j = 0; j < _edges.size(); ++j)
            if (value < _edges[j])
                --index;

        return index;
    }

    int64_t operator[](int64_t index) const { return _edges[index]; }

    int64_t size() const { return _size; }
    double const* raw() const { return _edges.data(); }

  private:
    int64_t _size;
    std::vector<double> _edges;
};

#endif /* INTERVAL_H */
