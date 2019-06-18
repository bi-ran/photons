#ifndef PAPER_H
#define PAPER_H

#include <array>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TObject.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLegend.h"

class paper {
  public:
    paper(std::string const& tag, int64_t size, int64_t cols, int64_t rows)
        : _tag(tag),
          _size(size),
          _cols(cols),
          _rows(rows),
          _desc(nullptr),
          canvas(nullptr) { }

    paper(std::string const& tag)
        : paper(tag, 0, 0, 0) { }

    paper(std::string const& tag, int64_t size)
        : paper(tag, size, 0, 0) { }

    paper(std::string const& tag, int64_t cols, int64_t rows)
        : paper(tag, 0, cols, rows) { }

    paper(paper const&) = delete;

    paper& operator=(paper const&) = delete;

    ~paper() = default;

    void add() { ++_size; }

    void stack(int64_t index, TObject* const object) {
        objects.push_back(object);
        indices.push_back(index);
    }

    void add(TObject* const object) { add(); stack(_size, object); }

    void stack(TObject* const object) { stack(_size, object); }

    void divide(int64_t cols, int64_t rows) { _cols = cols; _rows = rows; }

    void decorate(std::function<void()> d) { _d = d; }

    void format(std::function<void(TH1*)> f) { _f = f; }

    void format(std::function<void(TGraph*)> g) { _g = g; }

    void legend(std::function<std::array<float, 4>()> l) { _l = l; }

    void style(std::function<void(TLegend*)> s) { _s = s; }

    void describe(auto& desc) { _desc = &desc; }

    void draw(char const* ext) {
        using namespace std::literals::string_literals;

        if (_cols * _rows == 0) { split(); }

        if (canvas == NULL) {
            canvas = new TCanvas(("paper_"s + _tag).data(), "",
                             400 * _cols, 400 * _rows);
            canvas->Divide(_cols, _rows);

            int64_t count = static_cast<int64_t>(objects.size());
            for (int64_t i = 0; i < count; ++i) {
                canvas->cd(indices[i]);
                apply(objects[i], _f);
                apply(objects[i], _g);
                objects[i]->Draw("same p e");
                apply(_d);
            }

            legends();
        }

        save(ext);
    }

  private:
    template <typename T>
    void apply(TObject* const obj, std::function<void(T*)> f) {
        if (f && obj->InheritsFrom(T::Class()))
            f(static_cast<T*>(obj));
    }

    template <typename T>
    void apply(std::function<T> f) { if (f) { f(); } }

    void split() {
        float rows = std::ceil(std::sqrt(_size));
        float cols = std::ceil(_size / rows);

        _cols = cols;
        _rows = rows;
    }

    void legends() {
        if (_desc == nullptr) { return; }

        for (int64_t i = 1; i <= _size; ++i) {
            canvas->cd(i);

            std::vector<TObject*> associates;
            int64_t count = static_cast<int64_t>(objects.size());
            for (int64_t j = 0; j < count; ++j)
                if (indices[j] == i)
                    associates.push_back(objects[j]);

            auto point = _l ? _l() : std::array<float, 4>{
                0.5, 0.9, 0.87, 0.04 };
            point[3] = point[2] - associates.size() * point[3];

            TLegend* l = new TLegend(point[0], point[3], point[1], point[2]);
            apply(l, _s);

            for (auto const& obj : associates) {
                auto desc = _desc->find(obj) != std::end(*_desc) ?
                    (*_desc)[obj] : std::string(obj->GetName());
                l->AddEntry(obj, desc.data(), "pl");
            }

            l->Draw();
        }
    }

    void save(char const* ext) {
        using namespace std::literals::string_literals;

        canvas->SaveAs((_tag + "."s + ext).data());
    }

    std::string const _tag;
    int64_t _size;
    int64_t _cols;
    int64_t _rows;

    std::vector<TObject*> objects;
    std::vector<int64_t> indices;

    std::function<void()> _d;
    std::function<void(TH1*)> _f;
    std::function<void(TGraph*)> _g;
    std::function<std::array<float, 4>()> _l;
    std::function<void(TLegend*)> _s;

    std::map<TObject* const, std::string>* _desc;

    TCanvas* canvas;
};

#endif /* PAPER_H */
