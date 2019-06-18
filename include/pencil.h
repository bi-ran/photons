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

    template <typename... T>
    void category(std::string const& label, T const&... items) {
        (void) (int [sizeof...(T)]) { (category(label, items), 0)... };
    }

    template <typename... T>
    void describe(TObject* const object, T const&... adjectives) {
        (void) (int [sizeof...(T)]) { (describe(object, adjectives), 0)... };
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
    void category(std::string const& label, std::string const& item) {
        if (categories.find(label) == categories.end())
            categories[label] = { static_cast<int>(categories.size()) - 1, 0 };

        if (attributes.find(item) == attributes.end()) {
            attributes[item] = categories[label];
            ++categories[label][1];
        }
    }

    void describe(TObject* const object, std::string const& adjective) {
        if (objects.find(object) == objects.end())
            objects[object] = std::vector<int64_t>(categories.size());

        auto attr = attributes[adjective];
        objects[object][attr[0]] = attr[1];
    }

    std::map<TObject* const, std::vector<int64_t>> objects;

    std::map<std::string, std::array<int64_t, 2>> categories;
    std::map<std::string, std::array<int64_t, 2>> attributes;

    std::unique_ptr<pigment> core;
};

#endif /* PENCIL_H */
