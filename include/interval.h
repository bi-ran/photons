#ifndef INTERVAL_H
#define INTERVAL_H

#include <iterator>
#include <numeric>
#include <vector>

class interval {
  public:
    interval(int64_t number, double min, double max,
             std::string const& abscissa)
            : _size(number),
              _edges(std::vector<double>(number + 1)),
              _abscissa(abscissa) {
        std::iota(std::begin(_edges), std::end(_edges), 0);
        double interval = (max - min) / number;
        for (auto& edge : _edges) { edge = min + edge * interval; }
    }

    interval(int64_t number, double min, double max)
            : interval(number, min, max, std::string()) { }

    template <template <typename...> class T>
    interval(T<float> const& edges, std::string const& abscissa)
        : _size(edges.size() - 1),
          _edges(std::begin(edges), std::end(edges)),
          _abscissa(abscissa) { }

    template <template <typename...> class T>
    interval(T<float> const& edges)
        : interval(edges, std::string()) { }

    template <typename... T>
    interval(T const&... edges, std::string const& abscissa)
        : _size(sizeof...(T) - 1),
          _edges({static_cast<double>(edges)...}),
          _abscissa(abscissa) { }

    template <typename... T>
    interval(T const&... edges)
        : _size(sizeof...(T) - 1),
          _edges({static_cast<double>(edges)...}) { }

    interval(interval const& other) = default;
    interval& operator=(interval const& other) = default;
    ~interval() = default;

    int64_t index_for(double value) const {
        int64_t index = _size;
        for (auto edge : _edges)
            if (value < edge)
                --index;

        return index;
    }

    int64_t operator[](int64_t index) const { return _edges[index]; }

    std::string const abscissa() const { return _abscissa; }
    int64_t size() const { return _size; }

    template <typename T>
    T* book(std::string const& name, std::string const& title) {
        return new T(name.data(), title.data(), _size, _edges.data()); }

  private:
    int64_t _size;
    std::vector<double> _edges;

    std::string const _abscissa;
};

#endif /* INTERVAL_H */
