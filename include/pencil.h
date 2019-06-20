#ifndef PENCIL_H
#define PENCIL_H

#include <array>
#include <iterator>
#include <numeric>
#include <map>
#include <memory>
#include <string>
#include <vector>

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

static const std::vector<int32_t> solid = { 20, 21, 22, 23, 29, 41 };
static const std::vector<int32_t> open = { 24, 25, 26, 32, 30, 40 };

class pencil {
  public:
    pencil() : binary(-1) { }

    pencil(pencil const&) = delete;
    pencil& operator=(pencil const&) = delete;
    ~pencil() = default;

    template <typename... T>
    void category(std::string const& label, T const&... items) {
        (void) (int [sizeof...(T)]) { (categorise(label, items), 0)... }; }

    template <typename... T>
    void describe(TObject* const object, T const&... adjectives) {
        (void) (int [sizeof...(T)]) { (mark(object, adjectives), 0)... }; }

    void set_binary(std::string const& label) {
        binary = categories[label][0]; }

    template <typename T, template <typename...> class U>
    void operator()(T* const obj, U<int64_t> const& attrs) const {
        int64_t colour_index = 0;
        int64_t marker_index = 0;
        int64_t marker_type = 0;

        std::vector<int64_t> attributes(std::begin(attrs), std::end(attrs));
        int64_t dims = static_cast<int64_t>(features.size());
        for (int64_t i = 0; i < dims; ++i) {
            switch (features[i]) {
                case -1: marker_type = attributes[i]; break;
                case 0: colour_index = attributes[i]; break;
                case 1: marker_index = attributes[i]; break;
            }
        }

        int32_t colour = colours[colour_index];
        int32_t marker = (*(marker_type ? &open : &solid))[marker_index];

        apply<TH1>(obj, colour, marker);
        apply<TGraph>(obj, colour, marker);
    }

    template <typename T, template <typename...> class U>
    void operator()(T* const obj, U<std::string> const& adjectives) const {
        std::vector<int64_t> attrs(adjectives.size());
        for (auto const& adj : adjectives) {
            auto attr = attributes[adj];
            attrs[attr[0]] = attr[1];
        }

        (*this)(obj, attrs);
    }

    void sketch() {
        set_features(categories.size());
        for (auto const& obj : objects)
            (*this)(obj.first, obj.second);
    }

    void alias(std::string const& label, std::string const& formal) {
        aliases[label] = formal; }

    auto description() const {
        using namespace std::literals::string_literals;

        std::map<TObject* const, std::string> desc;

        /* build reverse map (intended to be called only once!) */
        std::map<std::array<int64_t, 2>, std::string> reverse;
        for (auto const& attr : attributes)
            reverse[attr.second] = attr.first;

        for (auto const& obj : objects) {
            std::string descriptive_string;
            auto const& indices = obj.second;
            int64_t count = static_cast<int64_t>(indices.size());
            for (int64_t i = 0; i < count; ++i) {
                auto attr = reverse[{ i, indices[i] }];
                auto it = aliases.find(attr);
                if (it != std::end(aliases))
                    attr = it->second;
                descriptive_string += attr + ", "s;
            }

            descriptive_string.pop_back();
            descriptive_string.pop_back();
            desc[obj.first] = descriptive_string;
        }

        return desc;
    }

  private:
    void categorise(std::string const& label, std::string const& item) {
        if (categories.find(label) == categories.end())
            categories[label] = { static_cast<int>(categories.size()) - 1, 0 };

        if (attributes.find(item) == attributes.end()) {
            attributes[item] = categories[label];
            ++categories[label][1];
        }
    }

    void mark(TObject* const object, std::string const& adjective) {
        if (objects.find(object) == objects.end())
            objects[object] = std::vector<int64_t>(categories.size());

        auto attr = attributes[adjective];
        objects[object][attr[0]] = attr[1];
    }

    void set_features(int64_t dims) {
        features = std::vector<int64_t>(dims);
        std::iota(std::begin(features), std::end(features), 0);

        if (binary != -1) {
            features[binary] = -1;
            for (int64_t i = binary + 1; i < dims; ++i)
                --features[i];
        }
    }

    template <typename T>
    void apply(TObject* const obj, int32_t colour, int32_t marker) const {
        if (obj->InheritsFrom(T::Class())) {
            auto cast = static_cast<T*>(obj);
            cast->SetLineColor(colour);
            cast->SetMarkerColor(colour);
            cast->SetMarkerStyle(marker);
        }
    }

    std::map<TObject* const, std::vector<int64_t>> objects;

    std::map<std::string, std::array<int64_t, 2>> categories;
    std::map<std::string, std::array<int64_t, 2>> attributes;

    std::map<std::string, std::string> aliases;

    std::vector<int64_t> features;
    int64_t binary;
};

#endif /* PENCIL_H */
