#ifndef PAPER_H
#define PAPER_H

#include <cmath>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TObject.h"

class paper {
  public:
    paper(std::string const& tag, int64_t size, int64_t cols, int64_t rows)
        : _tag(tag),
          _size(size),
          _cols(cols),
          _rows(rows),
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
                objects[i]->Draw("same p e");
            }
        }

        save(ext);
    }

  private:
    void split() {
        float rows = std::ceil(std::sqrt(_size));
        float cols = std::ceil(_size / rows);

        _cols = cols;
        _rows = rows;
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

    TCanvas* canvas;
};

#endif /* PAPER_H */
