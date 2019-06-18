#ifndef PENCIL_H
#define PENCIL_H

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "pigment.h"

class pencil {
  public:
    pencil() : core(std::make_unique<pigment>()) { }

    pencil(pencil const&) = delete;

    pencil& operator=(pencil const&) = delete;

    ~pencil() = default;

    void category(std::string const& label, std::string const& item) {
        if (categories.find(label) == categories.end()) {
            categories[label] = {
                static_cast<int64_t>(categories.size()) - 1, 0 };
        }

        if (attributes.find(item) == attributes.end()) {
            attributes[item] = categories[label];
            ++categories[label][1];
        }
    }

    template <typename... T>
    void category(std::string const& label, std::string const& first,
                  T const&... items) {
        category(label, first);
        category(label, items...);
    }

    void describe(TObject* const object, std::string const& description) {
        if (objects.find(object) == objects.end())
            objects[object] = std::vector<int64_t>(categories.size());

        auto attr = attributes[description];
        objects[object][attr[0]] = attr[1];
    }

    template <typename... T>
    void describe(TObject* const object, std::string const& first,
                  T const&... items) {
        describe(object, first);
        describe(object, items...);
    }

    void set_binary(std::string const& label) {
        core->set_binary(categories[label][0]);
    }

    void sketch() {
        core->set_features(categories.size());
        for (auto const& obj : objects)
            (*core)(obj.first, obj.second);
    }

  private:
    std::map<TObject* const, std::vector<int64_t>> objects;

    std::map<std::string, std::array<int64_t, 2>> categories;
    std::map<std::string, std::array<int64_t, 2>> attributes;

    std::unique_ptr<pigment> core;
};

#endif /* PENCIL_H */
