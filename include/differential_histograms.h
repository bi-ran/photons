#ifndef DIFFERENTIAL_HISTOGRAMS_H
#define DIFFERENTIAL_HISTOGRAMS_H

#include <algorithm>
#include <array>
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
                            std::string const& ordinate,
                            std::shared_ptr<interval> const& bins,
                            std::shared_ptr<multival> const& intervals)
            : _tag(tag),
              _ordinate(ordinate),
              _dims(intervals->dims()),
              _size(intervals->size()),
              _shape(intervals->shape()),
              bins(bins),
              intervals(intervals) {
        allocate_histograms();
    }

    template <template <typename...> class T>
    differential_histograms(std::string const& tag,
                            std::string const& ordinate,
                            std::shared_ptr<interval> const& bins,
                            T<int64_t> const& shape)
            : _tag(tag),
              _ordinate(ordinate),
              _dims(shape.size()),
              _size(std::accumulate(std::begin(shape), std::end(shape), 1,
                                    std::multiplies<int64_t>())),
              _shape(std::vector<int64_t>(std::begin(shape), std::end(shape))),
              bins(bins) {
        allocate_histograms();
    }

    template <typename... T>
    differential_histograms(std::string const& tag,
                            std::string const& ordinate,
                            std::shared_ptr<interval> const& bins,
                            T const&... dimensions)
            : _tag(tag),
              _ordinate(ordinate),
              _dims(sizeof...(T)),
              _size(size_of(dimensions...)),
              bins(bins) {
        auto shape = std::array<int64_t, sizeof...(T)>(dimensions...);
        _shape = std::vector<int64_t>(std::begin(shape), std::end(shape));

        allocate_histograms();
    }

    template <template <typename...> class T>
    differential_histograms(std::string const& tag,
                            std::string const& ordinate,
                            T<float> const& edges,
                            std::shared_ptr<multival> const& intervals)
        : differential_histograms(tag, ordinate,
                                  std::make_shared<interval>(edges),
                                  intervals) {
    }

    template <template <typename...> class T, template <typename...> class U>
    differential_histograms(std::string const& tag,
                            std::string const& ordinate,
                            T<float> const& edges,
                            U<int64_t> const& shape)
        : differential_histograms(tag, ordinate,
                                  std::make_shared<interval>(edges),
                                  shape) {
    }

    template <template <typename...> class T, typename... U>
    differential_histograms(std::string const& tag,
                            std::string const& ordinate,
                            T<float> const& edges,
                            U const&... dimensions)
        : differential_histograms(tag, ordinate,
                                  std::make_shared<interval>(edges),
                                  dimensions...) {
    }

    template <template <typename...> class T>
    differential_histograms(std::string const& tag,
                            std::string const& ordinate,
                            std::string const& abscissa,
                            T<float> const& edges,
                            std::shared_ptr<multival> const& intervals)
        : differential_histograms(tag, ordinate,
                                  std::make_shared<interval>(edges, abscissa),
                                  intervals) {
    }

    template <template <typename...> class T, template <typename...> class U>
    differential_histograms(std::string const& tag,
                            std::string const& ordinate,
                            std::string const& abscissa,
                            T<float> const& edges,
                            U<int64_t> const& shape)
        : differential_histograms(tag, ordinate,
                                  std::make_shared<interval>(edges, abscissa),
                                  shape) {
    }

    template <template <typename...> class T, typename... U>
    differential_histograms(std::string const& tag,
                            std::string const& ordinate,
                            std::string const& abscissa,
                            T<float> const& edges,
                            U const&... dimensions)
        : differential_histograms(tag, ordinate,
                                  std::make_shared<interval>(edges, abscissa),
                                  dimensions...) {
    }

    differential_histograms(differential_histograms const&) = delete;
    differential_histograms& operator=(differential_histograms const&) = delete;
    ~differential_histograms() = default;

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
        return intervals->index_for(values); }

    std::vector<int64_t> indices_for(int64_t index) const {
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
        this->add(other, 1); }

    void operator-=(differential_histograms const& other) {
        this->add(other, -1); }

    void scale(double c1) {
        for (auto const& hist : histograms)
            hist->Scale(c1);
    }

    void operator*=(double c1) { this->scale(c1); }

    void operator/=(double c1) { this->scale(1. / c1); }

    void multiply(differential_histograms const& other) {
        /* assume self, other have equal shapes */
        for (int64_t j = 0; j < _size; ++j) {
            auto count = other[j]->GetBinContent(1);
            histograms[j]->Scale(count);
        }
    }

    void divide(differential_histograms const& other) {
        /* assume self, other have equal shapes */
        for (int64_t j = 0; j < _size; ++j) {
            auto count = other[j]->GetBinContent(1);
            auto scale = count != 0 ? 1. / count : 0;
            histograms[j]->Scale(scale);
        }
    }

    /* divide histograms duplicated along an axis. useful when all histograms
     * on an axis are filled once per event. */
    void divide(differential_histograms const& other, int64_t axis) {
        /* assume self, other have equal shapes after integrating out
         * duplicated axis */
        for (int64_t j = 0; j < other.size(); ++j) {
            auto count = other[j]->GetBinContent(1);
            auto scale = count != 0 ? 1. / count : 0;
            auto indices = other.indices_for(j);
            indices.insert(std::next(std::begin(indices), axis), 0);
            for (int64_t k = 0; k < _shape[axis]; ++k) {
                indices[axis] = k;
                (*this)[indices]->Scale(scale);
            }
        }
    }

    void operator*=(differential_histograms const& other) {
        this->multiply(other); }

    void operator/=(differential_histograms const& other) {
        this->divide(other); }

    void divide(TH1* const other) {
        for (auto const& hist : histograms)
            hist->Divide(other);
    }

    TH1F*& operator[](int64_t index) { return histograms[index]; }

    TH1F* const& operator[](int64_t index) const { return histograms[index]; }

    template <template <typename...> class T>
    TH1F*& operator[](T<int64_t> const& indices) {
        return histograms[index_for(indices)]; }

    template <template <typename...> class T>
    TH1F* const& operator[](T<int64_t> const& indices) const {
        return histograms[index_for(indices)]; }

    TH1F* sum(std::vector<int64_t> const& indices, int64_t axis) const {
        using namespace std::literals::string_literals;

        std::vector<int64_t> output = indices;
        output.erase(std::next(std::begin(output), axis));

        std::string full_tag = _tag + "_sum"s + std::to_string(axis);
        for (auto const& index : output)
            full_tag = full_tag + "_"s + std::to_string(index);

        return sum_impl(full_tag, indices, axis, 0, _shape[axis]);
    }

    TH1F* sum(std::vector<int64_t> const& indices, int64_t axis,
              int64_t start, int64_t end) const {
        using namespace std::literals::string_literals;

        std::vector<int64_t> output = indices;
        output.erase(std::next(std::begin(output), axis));

        std::string full_tag = _tag + "_sum"s + std::to_string(axis)
            + "f"s + std::to_string(start) + "t"s + std::to_string(end);
        for (auto const& index : output)
            full_tag = full_tag + "_"s + std::to_string(index);

        return sum_impl(full_tag, indices, axis, start, end);
    }

    std::unique_ptr<differential_histograms> sum(int64_t axis) const {
        using namespace std::literals::string_literals;

        std::vector<int64_t> output = _shape;
        output.erase(std::next(std::begin(output), axis));

        auto result = std::make_unique<differential_histograms>(
            _tag + "_sum"s + std::to_string(axis), _ordinate, bins, output);

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

    template <typename T, template <typename...> class U, typename V,
              typename... W>
    T operator()(U<V> const& indices, T (TH1::* function)(W...),
                 W... args) {
        return forward(index_for(indices), function, args...); }

    template <typename T, template <typename...> class U, typename V,
              typename... W>
    T operator()(U<V> const& indices, T (TH1::* function)(W...) const,
                 W... args) {
        return forward(index_for(indices), function, args...); }

    void apply(std::function<void(TH1*)> f) {
        for (auto hist : histograms) { f(hist); } }

    int64_t dims() const { return _dims; }
    int64_t size() const { return _size; }
    std::vector<int64_t> const& shape() const { return _shape; }

  private:
    TH1F* sum_impl(std::string const& name, std::vector<int64_t> indices,
                   int64_t axis, int64_t start, int64_t end) const {
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
    T forward(int64_t index, T (U::* function)(V...) const, V... args) const {
        return ((*histograms[index]).*function)(std::forward<V>(args)...); }

    void allocate_histograms() {
        using namespace std::literals::string_literals;

        histograms = std::vector<TH1F*>(_size, 0);
        for (int64_t i = 0; i < _size; ++i) {
            std::string index_string;
            for (auto const& index : indices_for(i))
                index_string = index_string + "_"s + std::to_string(index);

            histograms[i] = bins->book<TH1F>(_tag + index_string,
                ";"s + bins->abscissa() + ";"s + _ordinate);
        }
    }

    template <typename T>
    int64_t size_of(T const& last) const {
        return last; }

    template <typename T, typename... U>
    int64_t size_of(T const& first, U const&... rest) const {
        return first * size_of(rest...); }

    std::string _tag;
    std::string _ordinate;

    int64_t _dims;
    int64_t _size;
    std::vector<int64_t> _shape;

    std::shared_ptr<interval> bins;
    std::shared_ptr<multival> intervals;
    std::vector<TH1F*> histograms;
};

#endif /* DIFFERENTIAL_HISTOGRAMS_H */
