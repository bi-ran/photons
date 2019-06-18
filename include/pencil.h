#ifndef PENCIL_H
#define PENCIL_H

#include <array>
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
        (void) (int [sizeof...(T)]) { (categorise(label, items), 0)... }; }

    template <typename... T>
    void describe(TObject* const object, T const&... adjectives) {
        (void) (int [sizeof...(T)]) { (mark(object, adjectives), 0)... }; }

    void set_binary(std::string const& label) {
        core->set_binary(categories[label][0]); }

    void sketch() {
        core->set_features(categories.size());
        for (auto const& obj : objects)
            (*core)(obj.first, obj.second);
    }

    void alias(std::string const& label, std::string const& formal) {
        aliases[label] = formal; }

    auto description() {
        using namespace std::literals::string_literals;

        std::map<TObject* const, std::string> desc;

        /* build reverse map (function intended to be called only once!) */
        std::map<std::array<int64_t, 2>, std::string> reverse;
        for (auto const& attr : attributes)
            reverse[attr.second] = attr.first;

        for (auto const& obj : objects) {
            std::string descriptive_string;
            auto const& indices = obj.second;
            int64_t count = static_cast<int64_t>(indices.size());
            for (int64_t i = 0; i < count; ++i) {
                auto attr = reverse[{ i, indices[i] }];
                if (aliases.find(attr) != std::end(aliases))
                    attr = aliases[attr];
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

    std::map<TObject* const, std::vector<int64_t>> objects;

    std::map<std::string, std::array<int64_t, 2>> categories;
    std::map<std::string, std::array<int64_t, 2>> attributes;

    std::map<std::string, std::string> aliases;

    std::unique_ptr<pigment> core;
};

#endif /* PENCIL_H */
