#ifndef PIGMENT_H
#define PIGMENT_H

#include <iterator>
#include <numeric>

#include "TColor.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TObject.h"

static const std::vector<int32_t> colours = {
    TColor::GetColor("#515151"),
    TColor::GetColor("#f2777a"),
    TColor::GetColor("#f99157"),
    TColor::GetColor("#ffcc66"),
    TColor::GetColor("#99cc99"),
    TColor::GetColor("#6699cc"),
    TColor::GetColor("#9999cc"),
};

static const std::vector<int32_t> solid = {
    20, 21, 22, 23, 29, 41 };

static const std::vector<int32_t> open = {
    24, 25, 26, 32, 30, 40 };

class pigment {
  public:
    pigment() : binary(-1) { }

    pigment(pigment const&) = delete;

    pigment& operator=(pigment const&) = delete;

    ~pigment() = default;

    template <typename T>
    void operator()(T* const obj, std::vector<int64_t> const& attrs) {
        int64_t colour_index = 0;
        int64_t marker_index = 0;
        int64_t marker_type = 0;

        int64_t dims = static_cast<int64_t>(features.size());
        for (int64_t i = 0; i < dims; ++i) {
            switch (features[i]) {
                case -1:
                    marker_type = attrs[i];
                    break;
                case 0:
                    colour_index = attrs[i];
                    break;
                case 1:
                    marker_index = attrs[i];
                    break;
            }
        }

        int32_t colour = colours[colour_index];
        int32_t marker = (*(marker_type ? &open : &solid))[marker_index];

        apply<TH1>(obj, colour, marker);
        apply<TGraph>(obj, colour, marker);
    }

    void set_binary(int64_t index) { binary = index; }

    void set_features(int64_t dims) {
        features = std::vector<int64_t>(dims);
        std::iota(std::begin(features), std::end(features), 0);

        if (binary != -1) {
            features[binary] = -1;
            for (int64_t i = binary + 1; i < dims; ++i)
                --features[i];
        }
    }

  private:
    template <typename T>
    void apply(TObject* const obj, int32_t colour, int32_t marker) {
        if (obj->InheritsFrom(T::Class())) {
            auto cast = static_cast<T*>(obj);
            cast->SetLineColor(colour);
            cast->SetMarkerColor(colour);
            cast->SetMarkerStyle(marker);
        }
    }

    int64_t binary;

    std::vector<int64_t> features;
};

#endif /* PIGMENT_H */
