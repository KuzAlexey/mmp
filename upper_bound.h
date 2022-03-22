// =====================================================================================
//
//       Filename:  upper_bound.h
//
//    Description:
//
//        Version:  1.0
//        Created:  2018年 07月 03日 星期二 11:07:01 CST
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Jinkun Lin, jkunlin@gmail.com
//   Organization:  School of EECS, Peking University
//
// =====================================================================================
#ifndef UPPER_BOUND_H_NMRGIRQI
#define UPPER_BOUND_H_NMRGIRQI

#include "basic.h"

// reduction
vector<int> indicator;
vector<int> vertex_to_removed;
vector<int> working_vertex;
vector<int> next_working_vertex;
vector<bool> is_pending;
vector<long> branch_ub;

// color
vector<int> vertex_color;
vector<int> color_max_weight;
vector<int> one_conflict_vertex;
vector<int> recolor_indicator;
vector<int> to_color_vertex;
vector<vector<int>> color_adj_list;
vector<vector<int>> color_indicator_list;
vector<long> color_ub;
vector<int> color_elem_cnt;

// require: vertex_color of all remaining vertex (not only adjacency_list[v])
// are -1
// corresponding indicator and one_conflict_vertex will be initialized here when
// new color used
bool comp(int v1, int v2) {
  if (vertex_weight[v1] > vertex_weight[v2]) {
    return true;
  } else if (vertex_weight[v1] == vertex_weight[v2]) {
    return color_adj_list[v1].size() > color_adj_list[v2].size();
  } else {
    return false;
  }
}

long branch_bound(int v) {
  if (branch_ub[v] != 0) {
    return branch_ub[v];
  }
  if (adjacency_list[v].empty()) {
    branch_ub[v] = vertex_weight[v];
    return vertex_weight[v];
  }
  int largest_weight_neighbor;
  int largest_weight = 0;
  for (auto u : adjacency_list[v]) {
    if (vertex_weight[u] > largest_weight) {
      largest_weight = vertex_weight[u];
      largest_weight_neighbor = u;
    }
  }

  int v1, v2;
  if (adjacency_list[v].size() <
      adjacency_list[largest_weight_neighbor].size()) {
    v1 = v;
    v2 = largest_weight_neighbor;
  } else {
    v1 = largest_weight_neighbor;
    v2 = v;
  }

  for (auto i : adjacency_list[v1]) {
    indicator[i] = false;
  }
  for (auto i : adjacency_list[v2]) {
    indicator[i] = true;
  }

  long common_neigbor_sum_weight = 0;
  for (auto i : adjacency_list[v1]) {
    if (indicator[i]) {
      common_neigbor_sum_weight += vertex_weight[i];
    }
  }

  long value_with_largest = vertex_weight[v] +
                            vertex_weight[largest_weight_neighbor] +
                            common_neigbor_sum_weight;
  long value_without_largest = vertex_weight[v] + vertex_neighbor_weight[v] -
                               vertex_weight[largest_weight_neighbor];

  branch_ub[v] = value_with_largest > value_without_largest
                     ? value_with_largest
                     : value_without_largest;
  return branch_ub[v];
}

long color_bound2(vector<int> &to_color_vertex) {
  // hack: use vertex_color as status
  for (auto u : to_color_vertex) {
    vertex_color[u] = 0;
    color_adj_list[u].clear();
    color_indicator_list[u].clear();
  }
  for (auto u : to_color_vertex) {
    for (auto w : adjacency_list[u]) {
      if (vertex_color[w] != -1) {
        color_adj_list[u].push_back(w);
      }
    }
    color_indicator_list[u].resize(color_adj_list[u].size(), 0);
  }
  std::sort(to_color_vertex.begin(), to_color_vertex.end(), comp);
  for (auto u : to_color_vertex) {
    vertex_color[u] = -1;
  }

  int color_num = 0;
  for (auto u : to_color_vertex) {
    for (int c = 0; c < color_num; ++c) {
      if (c >= (int)color_indicator_list[u].size() ||
          color_indicator_list[u][c] == 0) {
        vertex_color[u] = c;
        break;
      }
    }
    if (vertex_color[u] == -1) {
      vertex_color[u] = color_num++;
    }
    for (auto w : color_adj_list[u]) {
      if (vertex_color[u] < (int)color_indicator_list[w].size()) {
        color_indicator_list[w][vertex_color[u]]++;
      }
    }
  }

  color_max_weight.clear();
  color_max_weight.resize(color_num, 0);
  for (auto u : to_color_vertex) {
    int color = vertex_color[u];
    if (color_max_weight[color] < vertex_weight[u]) {
      color_max_weight[color] = vertex_weight[u];
    }
  }
  long ub = 0;
  for (int i = 0; i < color_num; ++i) {
    ub += color_max_weight[i];
  }
  // unset vertex_color for next invoke
  for (auto u : to_color_vertex) {
    vertex_color[u] = -1;
  }
  return ub;
}

long branch_color_bound(int v) {
  if (color_ub[v] != 0) {
    return color_ub[v];
  }
  if (adjacency_list[v].empty()) {
    color_ub[v] = vertex_weight[v];
    return color_ub[v];
  }
  // split
  int largest_weight_neighbor;
  int largest_weight = 0;
  for (auto u : adjacency_list[v]) {
    if (vertex_weight[u] > largest_weight) {
      largest_weight = vertex_weight[u];
      largest_weight_neighbor = u;
    }
  }

  int v1, v2;
  if (adjacency_list[v].size() <
      adjacency_list[largest_weight_neighbor].size()) {
    v1 = v;
    v2 = largest_weight_neighbor;
  } else {
    v1 = largest_weight_neighbor;
    v2 = v;
  }

  to_color_vertex.clear();
  for (auto i : adjacency_list[v1]) {
    indicator[i] = false;
  }
  for (auto i : adjacency_list[v2]) {
    indicator[i] = true;
  }

  for (auto i : adjacency_list[v1]) {
    if (indicator[i]) {
      to_color_vertex.push_back(i);
    }
  }

  long ub1 = vertex_weight[v] + vertex_weight[largest_weight_neighbor] +
             color_bound2(to_color_vertex);

  to_color_vertex.clear();
  for (auto u : adjacency_list[v]) {
    if (u != largest_weight_neighbor) {
      to_color_vertex.push_back(u);
    }
  }
  long ub2 = vertex_weight[v] + color_bound2(to_color_vertex);

  long ub = ub1;
  if (ub1 < ub2) {
    ub = ub2;
  }

  color_ub[v] = ub;
  return ub;
}

void rm_vertex(int v) {
  for (vector<int>::size_type i = 0; i < adjacency_list[v].size(); ++i) {
    int u = adjacency_list[v][i];
    vector<int>::size_type j = adjacency_pos_list[v][i];
    adjacency_list[u][j] = *adjacency_list[u].rbegin();
    adjacency_list[u].pop_back();
    // update adjacency_pos_list
    adjacency_pos_list[u][j] = *adjacency_pos_list[u].rbegin();
    adjacency_pos_list[u].pop_back();
    int w = adjacency_list[u][j];
    vector<int>::size_type k = adjacency_pos_list[u][j];
    adjacency_pos_list[w][k] = j;

    vertex_neighbor_weight[u] -= vertex_weight[v];
    color_ub[u] = 0;
    branch_ub[u] = 0;
    is_computed[u] = 0;
  }
  adjacency_list[v].clear();
  adjacency_pos_list[v].clear();
  remaining_vertex.remove(v);
}

bool simplify_branch_color() {
  auto old_size = remaining_vertex.size();

  long threshold_weight_degree = best_solution_weight;
  working_vertex.clear();

  for (auto v : remaining_vertex) {
    if (color_ub[v] == 0 || color_ub[v] <= threshold_weight_degree) {
      working_vertex.push_back(v);
      is_pending[v] = true;
    }
  }
  next_working_vertex.clear();

  while (!working_vertex.empty()) {
    for (vector<int>::size_type i = 0; i < working_vertex.size(); ++i) {
      auto v = working_vertex[i];
      if (vertex_weight[v] + vertex_neighbor_weight[v] <=
              threshold_weight_degree ||
          branch_color_bound(v) <= threshold_weight_degree) {
        for (auto u : adjacency_list[v]) {
          if (!is_pending[u]) {
            next_working_vertex.push_back(u);
            is_pending[u] = true;
          }
        }
        rm_vertex(v);
      }
      is_pending[v] = false;
    }
    working_vertex.clear();
    working_vertex.swap(next_working_vertex);
  }
  return old_size != remaining_vertex.size();
}

bool simplify_branch() {
  auto old_size = remaining_vertex.size();

  long threshold_weight_degree = best_solution_weight;
  working_vertex.clear();
  next_working_vertex.clear();

  for (auto v : remaining_vertex) {
    if (branch_ub[v] == 0 || branch_ub[v] <= threshold_weight_degree) {
      working_vertex.push_back(v);
      // true if v in working_vertex or next_working_vertex
      is_pending[v] = true;
    }
  }

  while (!working_vertex.empty()) {
    for (vector<int>::size_type i = 0; i < working_vertex.size(); ++i) {
      auto v = working_vertex[i];
      if (vertex_weight[v] + vertex_neighbor_weight[v] <=
              threshold_weight_degree ||
          branch_bound(v) <= threshold_weight_degree) {
        for (auto u : adjacency_list[v]) {
          if (!is_pending[u]) {
            next_working_vertex.push_back(u);
            is_pending[u] = true;
          }
        }
        rm_vertex(v);
      }
      is_pending[v] = false;
    }
    working_vertex.swap(next_working_vertex);
    next_working_vertex.clear();
  }

  return old_size != remaining_vertex.size();
}

bool simpilfy() {
  auto old_size = remaining_vertex.size();
  // simplify by degree
  long threshold_weight_degree = best_solution_weight;

  for (auto v : remaining_vertex) {
    if (vertex_weight[v] + vertex_neighbor_weight[v] <=
        threshold_weight_degree) {
      vertex_to_removed.push_back(v);
    }
  }

  while (!vertex_to_removed.empty()) {
    int v = *vertex_to_removed.rbegin();
    vertex_to_removed.pop_back();
    for (auto u : adjacency_list[v]) {
      if (vertex_weight[u] + vertex_neighbor_weight[u] >
              threshold_weight_degree &&
          vertex_weight[u] + vertex_neighbor_weight[u] - vertex_weight[v] <=
              threshold_weight_degree) {
        vertex_to_removed.push_back(u);
      }
    }
    rm_vertex(v);
  }

  return old_size != remaining_vertex.size();
}

#endif /* end of include guard: UPPER_BOUND_H_NMRGIRQI */
