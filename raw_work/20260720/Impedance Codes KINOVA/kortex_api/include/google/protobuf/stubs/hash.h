#ifndef GOOGLE_PROTOBUF_STUBS_HASH_H__
#define GOOGLE_PROTOBUF_STUBS_HASH_H__

#include <functional>
#include <string>
#include <utility>
#include <unordered_map>
#include <unordered_set>

namespace google {
  namespace protobuf {

    // Generic hash template for all types
    template <typename Key>
    struct hash {
      size_t operator()(const Key& key) const {
        return std::hash<Key>()(key);
      }
      bool operator()(const Key& a, const Key& b) const { return a < b; }
    };

    // Specialization for std::string
    template <>
    struct hash<std::string> {
      size_t operator()(const std::string& key) const {
        return std::hash<std::string>()(key);
      }
    };

    // Specialization for std::pair
    template <typename First, typename Second>
    struct hash<std::pair<First, Second>> {
      size_t operator()(const std::pair<First, Second>& p) const {
        hash<First> h1;
        hash<Second> h2;
        return h1(p.first) ^ (h2(p.second) << 1);
      }
    };

    // Aliases for hash_map and hash_set
    template <typename Key, typename Value, typename Hash = hash<Key>>
    using hash_map = std::unordered_map<Key, Value, Hash>;

    template <typename Key, typename Hash = hash<Key>>
    using hash_set = std::unordered_set<Key, Hash>;

  }  // namespace protobuf
}  // namespace google

#endif  // GOOGLE_PROTOBUF_STUBS_HASH_H__
