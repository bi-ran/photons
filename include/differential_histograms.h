#ifndef DIFFERENTIAL_HISTOGRAMS_H
#define DIFFERENTIAL_HISTOGRAMS_H

#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "interval.h"
#include "multival.h"

using x = std::initializer_list<int64_t>;
using v = std::initializer_list<double>;

class differential_histograms {
  public:
    differential_histograms(std::string const& tag,
                            std::shared_ptr<interval> const& bins,
                            std::shared_ptr<multival> const& intervals)
            : _tag(tag),
              _dims(intervals->dims()),
              _size(intervals->size()),
              _shape(intervals->shape()),
              bins(bins),
              intervals(intervals) {
        allocate_histograms();
    }

    template <template <typename...> class T>
    differential_histograms(std::string const& tag,
                            std::shared_ptr<interval> const& bins,
                            T<int64_t> const& shape)
            : _tag(tag),
              _dims(shape.size()),
              _size(std::accumulate(std::begin(shape), std::end(shape), 1,
                                    std::multiplies<int64_t>())),
              _shape(std::vector<int64_t>(std::begin(shape), std::end(shape))),
              bins(bins) {
        allocate_histograms();
    }

    template <typename... T>
    differential_histograms(std::string const& tag,
                            std::shared_ptr<interval> const& bins,
                            T const&... dimensions)
            : _tag(tag),
              _dims(sizeof...(T)),
              _size(size_of(dimensions...)),
              bins(bins) {
        auto array = shape_of(dimensions...);
        _shape = std::vector<int64_t>(std::begin(array), std::end(array));

        allocate_histograms();
    }

    template <template <typename...> class T>
    differential_histograms(std::string const& tag,
                            T<float> const& edges,
                            std::shared_ptr<multival> const& intervals)
        : differential_histograms(tag, std::make_shared<interval>(edges),
                                  intervals) {
    }

    template <template <typename...> class T, template <typename...> class U>
    differential_histograms(std::string const& tag,
                            T<float> const& edges,
                            U<int64_t> const& shape)
        : differential_histograms(tag, std::make_shared<interval>(edges),
                                  shape) {
    }

    template <template <typename...> class T, typename... U>
    differential_histograms(std::string const& tag,
                            T<float> const& edges,
                            U const&... dimensions)
        : differential_histograms(tag, std::make_shared<interval>(edges),
                                  dimensions...) {
    }

    differential_histograms(differential_histograms const&) = delete;

    differential_histograms& operator=(differential_histograms const&) = delete;

    ~differential_histograms() = default;

    template <template <typename...> class T>
    int64_t index_for(T<int64_t> const& indices) {
        std::vector<int64_t> vx(std::begin(indices), std::end(indices));
        int64_t index = 0;
        for (int64_t i = 0, block = 1; i < _dims; ++i) {
            index = index + vx[i] * block;
            block = block * _shape[i];
        }

        return index;
    }

    std::vector<int64_t> indices_for(int64_t index) {
        std::vector<int64_t> indices(_dims);
        for (int64_t i = 0; i < _dims; ++i) {
            indices[i] = index % _shape[i];
            index = index / _shape[i];
        }

        return indices;
    }

    void add(differential_histograms const& other, double c1) {
        for (int64_t i = 0; i < _size; ++i)
            histograms[i]->Add(other[i], c1);
    }

    void operator+=(differential_histograms const& other) {
        this->add(other, 1);
    }

    void operator-=(differential_histograms const& other) {
        this->add(other, -1);
    }

    TH1F*& operator[](int64_t index) {
        return histograms[index];
    }

    TH1F* const& operator[](int64_t index) const {
        return histograms[index];
    }

    template <template <typename...> class T>
    TH1F*& operator[](T<int64_t> const& indices) {
        return histograms[index_for(indices)];
    }

    template <template <typename...> class T>
    TH1F* const& operator[](T<int64_t> const& indices) const {
        return histograms[index_for(indices)];
    }

    TH1F* sum(std::vector<int64_t> const& indices, int64_t axis) {
        using namespace std::literals::string_literals;

        std::vector<int64_t> output = indices;
        output.erase(std::next(std::begin(output), axis));

        std::string full_tag = _tag + "_sum"s + std::to_string(axis);
        for (auto const& index : output)
            full_tag = full_tag + "_"s + std::to_string(index);

        return sum_impl(full_tag, indices, axis, 0, _shape[axis]);
    }

    TH1F* sum(std::vector<int64_t> const& indices, int64_t axis,
              int64_t start, int64_t end) {
        using namespace std::literals::string_literals;

        std::vector<int64_t> output = indices;
        output.erase(std::next(std::begin(output), axis));

        std::string full_tag = _tag + "_sum"s + std::to_string(axis)
            + "f"s + std::to_string(start) + "t"s + std::to_string(end);
        for (auto const& index : output)
            full_tag = full_tag + "_"s + std::to_string(index);

        return sum_impl(full_tag, indices, axis, start, end);
    }

    std::unique_ptr<differential_histograms> sum(int64_t axis) {
        using namespace std::literals::string_literals;

        std::vector<int64_t> output = _shape;
        output.erase(std::next(std::begin(output), axis));

        auto result = std::make_unique<differential_histograms>(
            _tag + "_sum"s + std::to_string(axis), bins, output);

        for (int64_t i = 0; i < result->size(); ++i) {
            auto indices = result->indices_for(i);
            indices.insert(std::next(std::begin(indices), axis), 0);
            (*result)[i] = this->sum(indices, axis);
        }

        return result;
    }

    template <typename T, typename... U>
    T operator()(int64_t index, T (TH1::* function)(U...), U... args) {
        return forward(index, function, args...); }

    template <typename T, typename... U>
    T operator()(int64_t index, T (TH1::* function)(U...) const, U... args) {
        return forward(index, function, args...); }

    template <typename T, template <typename...> class U, typename... V>
    T operator()(U<int64_t> const& indices, T (TH1::* function)(V...),
                 V... args) {
        return forward(index_for(indices), function, args...); }

    template <typename T, template <typename...> class U, typename... V>
    T operator()(U<int64_t> const& indices, T (TH1::* function)(V...) const,
                 V... args) {
        return forward(index_for(indices), function, args...); }

    template <typename T, template <typename...> class U, typename... V>
    T operator()(U<double> const& values, T (TH1::* function)(V...),
                 V... args) {
        auto indices = intervals->indices_for(values);
        return forward(index_for(indices), function, args...);
    }

    template <typename T, template <typename...> class U, typename... V>
    T operator()(U<double> const& values, T (TH1::* function)(V...) const,
                 V... args) {
        auto indices = intervals->indices_for(values);
        return forward(index_for(indices), function, args...);
    }

    int64_t dims() const { return _dims; }
    int64_t size() const { return _size; }
    std::vector<int64_t> const& shape() const { return _shape; }

  private:
    TH1F* sum_impl(std::string const& name, std::vector<int64_t> indices,
                   int64_t axis, int64_t start, int64_t end) {
        auto sum = static_cast<TH1F*>((*this)[indices]->Clone(name.data()));
        sum->Reset("MICES");
        for (int64_t i = start; i < end; ++i) {
            indices[axis] = i;
            sum->Add((*this)[indices]);
        }

        return sum;
    }

    template <typename T, typename U, typename... V>
    T forward(int64_t index, T (U::* function)(V...), V... args) {
        return ((*histograms[index]).*function)(std::forward<V>(args)...); }

    template <typename T, typename U, typename... V>
    T forward(int64_t index, T (U::* function)(V...) const, V... args) {
        return ((*histograms[index]).*function)(std::forward<V>(args)...); }

    void allocate_histograms() {
        using namespace std::literals::string_literals;

        histograms = std::vector<TH1F*>(_size, 0);
        for (int64_t i = 0; i < _size; ++i) {
            std::string index_string;
            for (auto const& index : indices_for(i))
                index_string = index_string + "_"s + std::to_string(index);

            histograms[i] = new TH1F(
                (_tag + index_string).data(), "",
                bins->size(), bins->raw()
            );
        }
    }

    template <typename... T, int64_t N = sizeof...(T)>
    std::array<int64_t, N> shape_of(T const&... dimensions) {
        return std::array<int64_t, N>(dimensions...); }

    template <typename T>
    int64_t size_of(T const& last) {
        return last; }

    template <typename T, typename... U>
    int64_t size_of(T const& first, U const&... rest) {
        return first * size_of(rest...); }

    std::string _tag;

    int64_t _dims;
    int64_t _size;
    std::vector<int64_t> _shape;

    std::shared_ptr<interval> bins;
    std::shared_ptr<multival> intervals;
    std::vector<TH1F*> histograms;
};

#endif /* DIFFERENTIAL_HISTOGRAMS_H */
