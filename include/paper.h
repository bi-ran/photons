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

#include "pencil.h"

class paper {
  public:
    paper(std::string const& tag, pencil* p, int64_t cols, int64_t rows)
        : _tag(tag),
          _size(0),
          _cols(cols),
          _rows(rows),
          _pencil(p),
          canvas(nullptr) { }

    paper(std::string const& tag)
        : paper(tag, 0, 0, 0) { }

    paper(std::string const& tag, pencil* p)
        : paper(tag, p, 0, 0) { }

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

    template <typename... T>
    void stack(int64_t index, TObject* const object, T const&... adjectives) {
        stack(index, object); _pencil->describe(object, adjectives...); }

    void add(TObject* const object) { add(); stack(_size, object); }

    template <typename... T>
    void add(TObject* const object, T const&... adjectives) {
        add(object); _pencil->describe(object, adjectives...); }

    void stack(TObject* const object) { stack(_size, object); }

    template <typename... T>
    void stack(TObject* const object, T const&... adjectives) {
        stack(object); _pencil->describe(object, adjectives...); }

    void divide(int64_t cols, int64_t rows) { _cols = cols; _rows = rows; }

    void decorate(std::function<void()> d) { _d = d; }
    void format(std::function<void(TH1*)> f) { _f = f; }
    void format(std::function<void(TGraph*)> g) { _g = g; }
    void legend(std::function<std::array<float, 4>()> l) { _l = l; }
    void style(std::function<void(TLegend*)> s) { _s = s; }

    void accessory(std::function<void(int64_t)> a) { _a.push_back(a); }

    void link(pencil* pencil) { _pencil = pencil; }

    void draw(std::string const& ext) {
        using namespace std::literals::string_literals;

        if (canvas == nullptr) {
            auto description = _pencil->description();

            layout();

            canvas = new TCanvas(("paper_"s + _tag).data(), "",
                                 400 * _cols, 400 * _rows);
            canvas->Divide(_cols, _rows);

            /* loop pads */
            for (int64_t i = 1; i <= _size; ++i) {
                canvas->cd(i);
                draw_pad(description, i);
            }
        }

        canvas->SaveAs((_tag + "."s + ext).data());
    }

    void split(std::string const& ext) const {
        using namespace std::literals::string_literals;

        auto description = _pencil->description();

        for (int64_t i = 1; i <= _size; ++i) {
            TCanvas* c = new TCanvas(("tmp_"s + std::to_string(i)).data(), "",
                                     400, 400);

            draw_pad(description, i);

            c->SaveAs((_tag + "_p"s + std::to_string(i) + "."s + ext).data());
            c->Delete();
        }
    }

  private:
    template <typename T>
    void apply(TObject* const obj, std::function<void(T*)> f) const {
        if (f && obj->InheritsFrom(T::Class()))
            f(static_cast<T*>(obj));
    }

    template <typename T>
    void apply(std::function<T> f) const { if (f) { f(); } }

    template <typename T>
    void apply(std::function<T> f, int64_t index) const { if (f) f(index); }

    void layout() {
        if (!_cols || !_rows) {
            float rows = std::ceil(std::sqrt(_size));
            float cols = std::ceil(_size / rows);

            _cols = cols;
            _rows = rows;
        } else if (_rows < 0) {
            _rows = std::ceil(_size / _cols);
        } else if (_cols < 0) {
            _cols = std::ceil(_size / _rows);
        }
    }

    std::vector<TObject*> associated(int64_t index) const {
        std::vector<TObject*> associates;
        int64_t count = static_cast<int64_t>(objects.size());
        for (int64_t j = 0; j < count; ++j)
            if (indices[j] == index)
                associates.push_back(objects[j]);

        return associates;
    }

    void draw_pad(auto const& description, int64_t index) const {
        auto associates = associated(index);
        for (auto const& obj : associates) {
            apply(obj, _f);
            apply(obj, _g);
            obj->Draw("same p e");
            apply(_d);

            for (auto const& a : _a)
                apply(a, index);
        }

        draw_legend(description, associates);
    }

    void draw_legend(auto const& description, auto const& associates) const {
        auto xy = _l ? _l() : std::array<float, 4>{ 0.5, 0.9, 0.87, 0.04 };
        xy[3] = xy[2] - associates.size() * xy[3];

        TLegend* l = new TLegend(xy[0], xy[3], xy[1], xy[2]);
        apply(l, _s);

        for (auto const& obj : associates) {
            auto it = description.find(obj);
            auto desc = it != std::end(description) ? it->second
                : std::string(obj->GetName());
            l->AddEntry(obj, desc.data(), "pl");
        }

        l->Draw();
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

    std::vector<std::function<void(int64_t)>> _a;

    pencil* _pencil;

    TCanvas* canvas;
};

#endif /* PAPER_H */
