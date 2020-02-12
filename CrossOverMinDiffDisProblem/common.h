#pragma once
#ifndef MIN_DIFF_DP_COMMON_H
#define MIN_DIFF_DP_COMMON_H
constexpr auto DISTANCE_MAX = 1.79769e+308;
#include <string>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>

namespace min_diff_dp {

template<typename T,typename N>
using Pair = std::pair<T, N>;

template<typename T>
using List = std::vector<T>;

template<typename T>
using Set = std::set<T>;

template<typename Key, typename Val>
using Map = std::map<Key, Val>;

template<typename Key, typename Val>
using HashMap = std::unordered_map<Key, Val>;

using String = std::string;

using Distance = double;

}

#endif