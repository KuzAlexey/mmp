#include "basic.h"
#include "upper_bound.h"

// build adjacency_list
void build(string file_name) {
  ifstream in_file(file_name);
  if (!in_file.is_open()) {
    cout << "in_file error" << endl;
    exit(1);
  }

  int vertex_count;

  // get vertex_count
  string line;
  istringstream is;
  string p, tmp;
  do {
    getline(in_file, line);
    is.clear();
    is.str(line);
    is >> p >> tmp >> vertex_count;
  } while (p != "p");

  // reading vertex weight
  vertex_weight.resize(vertex_count + 1);
  int v, w;
  for (vector<vector<int>>::size_type i = 1; i < vertex_weight.size(); ++i) {
    in_file >> tmp >> v >> w;
    if (tmp != "v")
      break;
    vertex_weight[v] = w;

    // v=i;
    // vertex_weight[v]=v%200+1;
  }

  adjacency_list.resize(vertex_count + 1);
  adjacency_pos_list.resize(vertex_count + 1);
  vertex_neighbor_weight.resize(vertex_count + 1, 0);

  int v1, v2;
  while (in_file >> tmp >> v1 >> v2) {
    adjacency_list[v1].push_back(v2);
    adjacency_list[v2].push_back(v1);
    adjacency_pos_list[v1].push_back(adjacency_list[v2].size() - 1);
    adjacency_pos_list[v2].push_back(adjacency_list[v1].size() - 1);
    vertex_neighbor_weight[v1] += vertex_weight[v2];
    vertex_neighbor_weight[v2] += vertex_weight[v1];

    // neighbor_hash_sum[v1+v2]=1;
  }
  in_file.close();
}

// build adjacency_list
void fast_build(string file_name) {
  ifstream in_file(file_name);
  if (!in_file) {
    cout << "in_file error" << endl;
    exit(1);
  }
  in_file.seekg(0, in_file.end);
  size_t file_len = in_file.tellg();
  in_file.seekg(0, in_file.beg);
  char *data = new char[file_len];

  in_file.read(data, file_len);
  in_file.close();

  // skip comments
  char *pos = data;
  while (*pos == '%') {
    while (*(pos++) != '\n')
      ;
  }

  // read vertex_count
  int vertex_count = 0, edge_count = 0;
  while (*pos < '0' || *pos > '9') {
    ++pos;
  }
  while (*pos != ' ') {
    vertex_count = vertex_count * 10 + *pos - '0';
    ++pos;
  }
  // read edge_count
  while (*pos < '0' || *pos > '9') {
    ++pos;
  }
  while (*pos >= '0' && *pos <= '9') {
    edge_count = edge_count * 10 + *pos - '0';
    ++pos;
  }

  // read vertex_weight
  vertex_weight.resize(vertex_count + 1);
  int v, w;
  for (vector<vector<int>>::size_type i = 1; i < vertex_weight.size(); ++i) {
    v = w = 0;
    // read vertex
    while (*pos < '0' || *pos > '9') {
      ++pos;
    }
    while (*pos != ' ') {
      v = v * 10 + *pos - '0';
      ++pos;
    }
    // read weight
    while (*pos < '0' || *pos > '9') {
      ++pos;
    }
    while (*pos >= '0' && *pos <= '9') {
      w = w * 10 + *pos - '0';
      ++pos;
    }
    vertex_weight[v] = w;
  }

  // read adjacency_list
  adjacency_list.resize(vertex_count + 1);
  adjacency_pos_list.resize(vertex_count + 1);
  vertex_neighbor_weight.resize(vertex_count + 1, 0);

  vector<int> vertex_degree(vertex_count + 1, 0);
  char *stash_pos = pos;
  int v1, v2;
  for (int i = 0; i < edge_count; ++i) {
    v1 = v2 = 0;
    // read v1
    while (*pos < '0' || *pos > '9') {
      ++pos;
    }
    while (*pos != ' ') {
      v1 = v1 * 10 + *pos - '0';
      ++pos;
    }
    // read v2
    while (*pos < '0' || *pos > '9') {
      ++pos;
    }
    while (*pos >= '0' && *pos <= '9') {
      v2 = v2 * 10 + *pos - '0';
      ++pos;
    }
    vertex_degree[v1]++;
    vertex_degree[v2]++;
  }
  for (size_t v = 1; v < adjacency_list.size(); ++v) {
    adjacency_list[v].resize(vertex_degree[v]);
    adjacency_pos_list[v].resize(vertex_degree[v]);
  }
  pos = stash_pos;
  for (int i = 0; i < edge_count; ++i) {
    v1 = v2 = 0;
    // read v1
    while (*pos < '0' || *pos > '9') {
      ++pos;
    }
    while (*pos != ' ') {
      v1 = v1 * 10 + *pos - '0';
      ++pos;
    }
    // read weight
    while (*pos < '0' || *pos > '9') {
      ++pos;
    }
    while (*pos >= '0' && *pos <= '9') {
      v2 = v2 * 10 + *pos - '0';
      ++pos;
    }

    adjacency_list[v1][--vertex_degree[v1]] = v2;
    adjacency_list[v2][--vertex_degree[v2]] = v1;
    adjacency_pos_list[v1][vertex_degree[v1]] = vertex_degree[v2];
    adjacency_pos_list[v2][vertex_degree[v2]] = vertex_degree[v1];

    vertex_neighbor_weight[v1] += vertex_weight[v2];
    vertex_neighbor_weight[v2] += vertex_weight[v1];
  }
  delete[] data;
}

void output_graph_size() {
  int edge_count = 0;
  for (auto v : remaining_vertex) {
    edge_count += adjacency_list[v].size();
  }
  edge_count /= 2;
  cout << "p edge " << remaining_vertex.size() << ' ' << edge_count << " "
       << get_time() << endl;
}

void update_best_solution() {
  best_solution.clear();
  for (auto v : solution) {
    best_solution.push_back(v);
  }
  best_solution_weight = solution_weight;
  best_solution_time = get_time();
  best_solution_try = tries;
}

vector<int> start_vertices;
int untest_pointer;

vector<int> candidates;
vector<long> cand_neighbor_weight;
vector<bool> is_in_candidates;

vector<bool> is_addv_neighbor; // indicates whether a candidate vertex is
                               // adjacent to the add_v

bool is_new_graph = true;

void build_init_clique(const vector<int> &vertex, const vector<int> &tmp_degree,
                       const vector<size_t> &pos) {
  best_solution_weight = 0;
  best_solution.clear();

  int max_degree = tmp_degree[*vertex.rbegin()];
  auto pos_head = pos[max_degree];
  for (auto cur_pos = vertex.size(); cur_pos-- > pos_head;) {
    // for (auto cur_pos = vertex.size() - 1; cur_pos >= pos_head; --cur_pos) {
    solution.clear();
    solution_weight = 0;
    candidates.clear();
    long sum_candidate_weight = 0;

    int v = vertex[cur_pos];
    candidates = adjacency_list[v];
    for (auto u : candidates) {
      is_in_candidates[u] = 1;
      sum_candidate_weight += vertex_weight[u];
    }
    solution.push_back(v);
    solution_weight += vertex_weight[v];

    while (!candidates.empty()) {
      if (solution_weight + sum_candidate_weight <= best_solution_weight) {
        for (auto v : candidates) {
          is_in_candidates[v] = 0;
        }
        break;
      }
      int best_addv = candidates[0];
      for (vector<int>::size_type i = 1; i < candidates.size(); ++i) {
        int u = candidates[i];
        if (tmp_degree[u] > tmp_degree[best_addv] ||
            (tmp_degree[u] == tmp_degree[best_addv] &&
             vertex_weight[u] > vertex_weight[best_addv])) {
          best_addv = u;
        }
      }
      solution.push_back(best_addv);
      solution_weight += vertex_weight[best_addv];

      // update candidates
      for (auto u : adjacency_list[best_addv]) {
        if (is_in_candidates[u]) {
          is_addv_neighbor[u] = 1;
        }
      }
      for (vector<int>::size_type i = 0; i < candidates.size();) {
        int u = candidates[i];
        if (is_addv_neighbor[u] == 0) {
          candidates[i] = *candidates.rbegin();
          candidates.pop_back();
          sum_candidate_weight -= vertex_weight[u];
          is_in_candidates[u] = 0; // reset is_in_candidates
        } else {
          is_addv_neighbor[u] = 0; // reset is_addv_neighbor
          ++i;
        }
      }
    }
    if (solution_weight > best_solution_weight) {
      best_solution_weight = solution_weight;
      best_solution.swap(solution);
    }
  }
}

// TODO: test
void kcore_init() {
  vertex.clear();
  vertex.reserve(adjacency_list.size() - 1);
  index.resize(adjacency_list.size());

  for (vector<vector<int>>::size_type v = 1; v < adjacency_list.size(); ++v) {
    vertex.push_back(v);
    tmp_degree[v] = adjacency_list[v].size();
  }

  int max_degree = 0;
  for (auto v : vertex) {
    if (max_degree < tmp_degree[v]) {
      max_degree = tmp_degree[v];
    }
  }
  vector<size_t> pos(max_degree + 1, 0);
  for (auto v : vertex) {
    pos[tmp_degree[v]]++;
  }
  // end position
  for (size_t i = 1; i < pos.size(); ++i) {
    pos[i] += pos[i - 1];
  }
  // begin position
  for (vector<int>::size_type v = 1; v <= vertex.size(); ++v) {
    auto degree = tmp_degree[v];
    vertex[--pos[degree]] = v;
    // note: destroy consisdent
    index[v] = pos[degree];
  }
  // core decomposition
  for (vector<int>::size_type i = 0; i < vertex.size(); ++i) {
    auto v = vertex[i];
    // is it okay to break here?
    if (tmp_degree[v] == max_degree ||
        pos[tmp_degree[v] + 1] == vertex.size()) {
      break;
    }
    for (auto u : adjacency_list[v]) {
      if (tmp_degree[u] > tmp_degree[v]) {
        // exchange u and begin position
        auto index_w = pos[tmp_degree[u]], index_u = index[u];
        auto w = vertex[index_w];
        vertex[index_u] = w;
        vertex[index_w] = u;
        index[u] = index_w;
        index[w] = index_u;
        ++pos[tmp_degree[u]];
        --tmp_degree[u];
      }
    }
  }
  build_init_clique(vertex, tmp_degree, pos);
  remaining_vertex.init(vertex, index);
}

void init() {
  solution.reserve(adjacency_list.size());
  best_solution.reserve(adjacency_list.size());

  indicator.resize(adjacency_list.size());
  working_vertex.reserve(adjacency_list.size());
  next_working_vertex.reserve(adjacency_list.size());
  is_pending.resize(adjacency_list.size());

  // remaining_vertex.init(adjacency_list.size());

  vertex_color.resize(adjacency_list.size(), -1);
  recolor_indicator.resize(adjacency_list.size());
  color_adj_list.resize(adjacency_list.size());
  color_indicator_list.resize(adjacency_list.size());

  cand_neighbor_weight.resize(adjacency_list.size(), 0);
  is_in_candidates.resize(adjacency_list.size(), 0);
  is_computed.resize(adjacency_list.size(), 0);
  is_addv_neighbor.resize(adjacency_list.size(), 0);
  adjacency_cand_neighbor_weight.resize(adjacency_list.size());
  for (vector<int>::size_type v = 1; v < adjacency_list.size(); ++v) {
    adjacency_cand_neighbor_weight[v].resize(adjacency_list[v].size());
  }

  real_bms_count = min_bms_count;

  color_ub.resize(adjacency_list.size(), 0);
  branch_ub.resize(adjacency_list.size(), 0);

  // using kcore to calculate first clique
  tmp_degree.resize(adjacency_list.size());
  // kcore_init();
  remaining_vertex.init(adjacency_list.size() - 1);
  vertex.reserve(adjacency_list.size() - 1);
  index.resize(adjacency_list.size());
}

void branch(vector<int> &candidates, vector<int> &cur_clq, long cur_weight,
            vector<int> &best_clq, long &best_weight) {
  if (candidates.empty()) {
    if (best_weight < cur_weight) {
      best_weight = cur_weight;
      best_clq = cur_clq;
    }
    return;
  }

  for (size_t i = 0; i < candidates.size(); ++i) {
    int v = candidates[i];
    cur_clq.push_back(v);

    vector<int> tmp_candidates;
    for (auto u : adjacency_list[v]) {
      indicator[u] = 0;
    }
    for (size_t j = i + 1; j < candidates.size(); ++j) {
      int u = candidates[j];
      indicator[u] = 1;
    }
    for (auto u : adjacency_list[v]) {
      if (indicator[u] == 1) {
        tmp_candidates.push_back(u);
      }
    }

    if (cur_weight + vertex_weight[v] + color_bound2(tmp_candidates) >=
        best_weight) {
      branch(tmp_candidates, cur_clq, cur_weight + vertex_weight[v], best_clq,
             best_weight);
    }
    cur_clq.pop_back();
  }
}

void hill_climb() {
  const int BRANCH_AND_BOUND_THREDSHOLD = 50;
  if (solution.size() == 1 &&
      remaining_vertex.size() < BRANCH_AND_BOUND_THREDSHOLD) {
    vector<int> best_clq, cur_clq;

    best_clq = solution;
    long best_weight = solution_weight;

    candidates = remaining_vertex.vertex;
    branch(candidates, cur_clq, 0, best_clq, best_weight);

    // update best_solution
    if (best_weight > solution_weight) {
      solution_weight = best_weight;
      solution = best_clq;
    }
    return;
  }

  while (true) {
    bool is_improve = false;
    size_t end_index = solution.size(); // avoid scanning twice
    for (size_t rm_index = 0; rm_index < solution.size(); ++rm_index) {
      int rm_v = solution[rm_index];

      // TODO: cache indicator to speed up
      // build candidates
      candidates.clear();
      int first_v = -1;
      for (auto v : solution) {
        if (v == rm_v) {
          continue;
        }
        if (first_v == -1) {
          first_v = v;
          for (auto u : adjacency_list[v]) {
            if (u != rm_v) {
              indicator[u] = 1;
            }
          }
        } else {
          for (auto u : adjacency_list[v]) {
            if (u != rm_v) {
              indicator[u]++;
            }
          }
        }
      }
      for (auto u : adjacency_list[first_v]) {
        if (indicator[u] == (int)solution.size() - 1) {
          candidates.push_back(u);
        }
      }

      // build clique candidates
      vector<int> best_clq, cur_clq;
      best_clq.push_back(rm_v);
      long best_weight = vertex_weight[rm_v];

      // skip large branch case
      if (candidates.size() < BRANCH_AND_BOUND_THREDSHOLD) {
        branch(candidates, cur_clq, 0, best_clq, best_weight);
      }

      // update best_solution
      // TODO: equal
      if (best_weight > vertex_weight[rm_v]) {
        solution_weight = solution_weight - vertex_weight[rm_v] + best_weight;
        solution[rm_index] = best_clq[0];
        for (size_t i = 1; i < best_clq.size(); ++i) {
          solution.push_back(best_clq[i]);
        }
        is_improve = true;

        end_index = rm_index;
      } else {
        if (rm_index == end_index) {
          break;
        }
      }
    }
    if (!is_improve) {
      break;
    }
  }
}

void improve_solution() {
  for (auto remove_index = solution.size(); remove_index-- > 0;) {
    if (solution.size() == 1) {
      candidates = remaining_vertex.vertex;
      for (auto v : candidates) {
        is_in_candidates[v] = true;
      }
    } else {
      size_t base_index;
      if (remove_index == 0) {
        base_index = 1;
      } else {
        base_index = 0;
      }
      for (auto v : adjacency_list[solution[base_index]]) {
        candidates.push_back(v);
        is_in_candidates[v] = true;
      }
    }

    for (size_t i = 0; i < solution.size(); ++i) {
      if (i == remove_index) {
        continue;
      }
      int v = solution[i];
      for (auto u : adjacency_list[v]) {
        if (is_in_candidates[u]) {
          is_addv_neighbor[u] = true;
        }
      }
      for (size_t j = 0; j < candidates.size();) {
        int u = candidates[j];
        if (!is_addv_neighbor[u]) {
          candidates[j] = *candidates.rbegin();
          candidates.pop_back();
          is_in_candidates[u] = false;
        } else {
          ++j;
          is_addv_neighbor[u] = false;
        }
      }
    }
    // candidates.size() > 0, because removed vertex is in candidates
    int add_v = candidates[0];
    long best_score = vertex_weight[add_v];
    for (size_t i = 1; i < candidates.size(); ++i) {
      int v = candidates[i];
      if (best_score < vertex_weight[v]) {
        add_v = v;
        best_score = vertex_weight[v];
      }
    }
    if (best_score > vertex_weight[solution[remove_index]]) {
      solution_weight -= vertex_weight[solution[remove_index]];
      solution_weight += best_score;
      solution[remove_index] = add_v;

      while (true) {
        for (auto v : adjacency_list[add_v]) {
          if (is_in_candidates[v]) {
            is_addv_neighbor[v] = true;
          }
        }
        for (size_t j = 0; j < candidates.size();) {
          int u = candidates[j];
          if (!is_addv_neighbor[u]) {
            candidates[j] = *candidates.rbegin();
            candidates.pop_back();
            is_in_candidates[u] = false;
          } else {
            ++j;
            is_addv_neighbor[u] = false;
          }
        }
        if (!candidates.empty()) {
          add_v = candidates[0];
          best_score = vertex_weight[add_v];
          for (size_t i = 1; i < candidates.size(); ++i) {
            int v = candidates[i];
            if (best_score < vertex_weight[v]) {
              add_v = v;
              best_score = vertex_weight[v];
            }
          }
          solution.push_back(add_v);
          solution_weight += vertex_weight[add_v];
        } else {
          break;
        }
      }
    }
    for (auto v : candidates) {
      is_in_candidates[v] = false;
    }
    candidates.clear();
  }
}

int construct2() {
  if (is_new_graph) {
    start_vertices = remaining_vertex.vertex;
    untest_pointer = 0;
  }

  if (untest_pointer == (int)start_vertices.size()) {
    untest_pointer = 0;
    real_bms_count = real_bms_count * 2;
    if (real_bms_count > max_bms_count) {
      real_bms_count = min_bms_count;
    }
  }

  // choose start_v
  int start_index;
  long best_score = -1;
  for (int i = 0; i < start_bms_count; ++i) {
    int index =
        untest_pointer + rand() % (start_vertices.size() - untest_pointer);
    int v = start_vertices[index];
    long score = vertex_weight[v] + vertex_neighbor_weight[v] / t;

    if (best_score < score) {
      start_index = index;
      best_score = score;
    }
  }
  int start_v = start_vertices[start_index];
  swap(start_vertices[untest_pointer++], start_vertices[start_index]);

  solution.clear();
  solution.push_back(start_v);
  solution_weight = vertex_weight[start_v];

  // init candidates
  candidates.clear();
  for (auto v : adjacency_list[start_v]) {
    candidates.push_back(v);
    is_in_candidates[v] = true;
  }

  // init cand_neighbor_weight
  if (is_computed[start_v]) {
    for (size_t i = 0; i < adjacency_list[start_v].size(); ++i) {
      int v = adjacency_list[start_v][i];
      cand_neighbor_weight[v] = adjacency_cand_neighbor_weight[start_v][i];
    }
  } else {
    for (size_t i = 0; i < adjacency_list[start_v].size(); ++i) {
      int v = adjacency_list[start_v][i];
      cand_neighbor_weight[v] = 0;
      for (auto u : adjacency_list[v]) {
        if (is_in_candidates[u]) {
          cand_neighbor_weight[v] += vertex_weight[u];
        }
      }
      adjacency_cand_neighbor_weight[start_v][i] = cand_neighbor_weight[v];
    }
    is_computed[start_v] = 1;
  }

  bool force_hill_climb = false;
  if (rand() % 10000 < 8) {
    force_hill_climb = true;
  }

  // construct clique
  while (!candidates.empty()) {
    int add_index;
    long best_score = -1;
    if ((int)candidates.size() < real_bms_count) {
      for (size_t index = 0; index < candidates.size(); ++index) {
        int v = candidates[index];
        long score = vertex_weight[v] + cand_neighbor_weight[v] / t;
        if (best_score < score) {
          add_index = index;
          best_score = score;
        }
      }
    } else {
      for (int i = 0; i < real_bms_count; ++i) {
        int index = rand() % candidates.size();
        int v = candidates[index];
        long score = vertex_weight[v] + cand_neighbor_weight[v] / t;
        if (best_score < score) {
          add_index = index;
          best_score = score;
        }
      }
    }

    int add_v = candidates[add_index];
    // upper bound prune
    if (!force_hill_climb &&
        solution_weight + cand_neighbor_weight[add_v] + vertex_weight[add_v] <=
            best_solution_weight) {
      for (auto v : candidates) {
        is_in_candidates[v] = false;
      }
      return false;
    }

    solution.push_back(add_v);
    solution_weight += vertex_weight[add_v];

    // update candidates, true is false, false is true
    for (auto v : adjacency_list[add_v]) {
      is_in_candidates[v] = false;
    }
    for (size_t i = 0; i < candidates.size();) {
      int v = candidates[i];
      if (is_in_candidates[v]) {
        candidates[i] = *candidates.rbegin();
        candidates.pop_back();
        is_in_candidates[v] = false;

        // update cand_neighbor_weight
        for (auto u : adjacency_list[v]) {
          cand_neighbor_weight[u] -= vertex_weight[v];
        }
      } else {
        is_in_candidates[v] = true;
        ++i;
      }
    }
  }

  bool is_improve = false;
  if (solution_weight > best_solution_weight) {
    update_best_solution();
    is_improve = true;
  }

  hill_climb();

  // bool is_improve = false;
  if (solution_weight > best_solution_weight) {
    update_best_solution();
    is_improve = true;
  }

  return is_improve;
}

bool construct() {
  solution.clear();
  candidates.clear();
  solution_weight = 0;

  int startv;
  vector<int>::size_type index, tmp_index;

  if (is_new_graph) {
    start_vertices = remaining_vertex.vertex;
    untest_pointer = 0;
  }

  if (untest_pointer == (int)start_vertices.size()) {
    untest_pointer = 0;

    real_bms_count = real_bms_count * 2;
    if (real_bms_count > max_bms_count)
      real_bms_count = min_bms_count;
  }
  long best_score;
  long tmp_score;

  // pick the starting vertex
  index = untest_pointer + rand() % (start_vertices.size() - untest_pointer);
  startv = start_vertices[index];
  best_score = vertex_weight[startv] + vertex_neighbor_weight[startv] / t;
  for (int i = 1; i <= start_bms_count; ++i) {
    tmp_index =
        untest_pointer + rand() % (start_vertices.size() - untest_pointer);
    tmp_score = vertex_weight[start_vertices[tmp_index]] +
                vertex_neighbor_weight[start_vertices[tmp_index]] / t;
    if (tmp_score > best_score) {
      best_score = tmp_score;
      index = tmp_index;
    }
  }
  startv = start_vertices[index];
  std::swap(start_vertices[index], start_vertices[untest_pointer++]);

  solution.push_back(startv);
  solution_weight += vertex_weight[startv];

  // initialize the set of candidate vertices S = N(startv)
  for (auto u : adjacency_list[startv]) {
    candidates.push_back(u);
    is_in_candidates[u] = 1;
  }

  // initialize the cand_neighbor_weight value of vertices in candidate set S
  int i = 0;
  if (is_computed[startv] == 0) {
    for (auto v : candidates) {
      cand_neighbor_weight[v] = 0;
      for (auto n : adjacency_list[v]) {
        if (is_in_candidates[n] == 1)
          cand_neighbor_weight[v] += vertex_weight[n];
      }
      // record the beginning cand_neighbor_weight value of vertices in the
      // beginning candidate set S
      adjacency_cand_neighbor_weight[startv][i++] = cand_neighbor_weight[v];
    }

    is_computed[startv] = 1;
  } else {
    i = 0;
    for (auto v : adjacency_list[startv]) {
      cand_neighbor_weight[v] = adjacency_cand_neighbor_weight[startv][i++];
    }
  }

  int add_v;
  long max_score;

  while (!candidates.empty()) {
    // pick add_v
    if ((int)candidates.size() < real_bms_count) {
      max_score = 0;
      index = 0;
      for (vector<int>::size_type i = 0; i < candidates.size(); ++i) {
        tmp_score = cand_neighbor_weight[candidates[i]] / t +
                    vertex_weight[candidates[i]];
        if (tmp_score > max_score) {
          max_score = tmp_score;
          index = i;
        }
      }

    } else {
      index = rand() % candidates.size();
      max_score = cand_neighbor_weight[candidates[index]] / t +
                  vertex_weight[candidates[index]];

      for (int i = 1; i < real_bms_count; ++i) {
        tmp_index = rand() % candidates.size();
        tmp_score = cand_neighbor_weight[candidates[tmp_index]] / t +
                    vertex_weight[candidates[tmp_index]];
        if (tmp_score > max_score) {
          max_score = tmp_score;
          index = tmp_index;
        }
      }
    }
    add_v = candidates[index];

    // upper bound prune
    if (solution_weight + cand_neighbor_weight[add_v] + vertex_weight[add_v] <=
        best_solution_weight) {
      for (auto v : candidates) {
        is_in_candidates[v] = 0;
      }
      return false;
    }

    solution.push_back(add_v);
    solution_weight += vertex_weight[add_v];

    // remove add_v and update its neighbors information
    for (auto u : adjacency_list[add_v]) {
      if (is_in_candidates[u] == 1)
        cand_neighbor_weight[u] -= vertex_weight[add_v];
    }
    is_in_candidates[add_v] = 0;
    candidates[index] = *(candidates.end() - 1);
    candidates.pop_back();

    for (auto v : adjacency_list[add_v]) {
      if (is_in_candidates[v] == 1)
        is_addv_neighbor[v] = 1;
    }

    for (vector<int>::size_type i = 0; i < candidates.size();) {
      // remove the vertice doesn't belong to the neighborhood of add_v

      int cur_v = candidates[i];

      if (is_addv_neighbor[cur_v] == 0) {
        // update
        for (auto u : adjacency_list[cur_v]) {
          if (is_in_candidates[u] == 1)
            cand_neighbor_weight[u] -= vertex_weight[cur_v];
        }
        is_in_candidates[cur_v] = 0;
        candidates[i] = *(candidates.end() - 1);
        candidates.pop_back();
      } else {
        is_addv_neighbor[cur_v] = 0; // erase the value here, making this
                                     // array clean after the step.
        ++i;
      }
    }
  }

  // improve_solution();
  hill_climb();

  bool is_improve;
  if (solution_weight > best_solution_weight) {
    update_best_solution();
    is_improve = true;
  } else {
    is_improve = false;
  }

  return is_improve;
}

bool verify_simple(string file_name) {
  ifstream in_file(file_name);
  if (!in_file.is_open()) {
    cout << "in_file error" << endl;
    exit(1);
  }

  int vertex_count;

  // get vertex_count
  string line;
  istringstream is;
  string p, tmp;
  do {
    getline(in_file, line);
    is.clear();
    is.str(line);
    is >> p >> tmp >> vertex_count;
  } while (p != "p");

  // reading vertex weight
  vector<int> vw;
  vw.resize(vertex_count + 1);
  int v, w;
  for (vector<vector<int>>::size_type i = 1; i < vw.size(); ++i) {
    in_file >> tmp >> v >> w;
    vw[v] = w;
  }

  // check solution weight
  long sw = 0;
  for (auto v : best_solution) {
    sw += vw[v];
  }
  if (sw != best_solution_weight) {
    cout << "wrong weight!" << endl;
    return false;
  }

  // adjacency list
  vector<vector<int>> adj_list;
  adj_list.resize(vertex_count + 1);
  int v1, v2;
  while (in_file >> tmp >> v1 >> v2) {
    adj_list[v1].push_back(v2);
    adj_list[v2].push_back(v1);
  }
  in_file.close();

  // check clique
  for (vector<int>::size_type i = 0; i < best_solution.size(); ++i) {
    int v1 = best_solution[i];
    for (vector<int>::size_type j = i + 1; j < best_solution.size(); ++j) {
      int v2 = best_solution[j];
      vector<int>::size_type k = 0;
      for (; k < adj_list[v1].size(); ++k) {
        if (v2 == adj_list[v1][k]) {
          break;
        }
      }
      if (k >= adj_list[v1].size()) {
        cerr << "wrong anser: " << v1 << " is not neighbor of " << v2 << "!"
             << endl;
        return false;
      }
    }
  }
  return true;
}

bool verify_fast(string file_name) {
  // release memory
  adjacency_list.clear();
  adjacency_pos_list.clear();
  adjacency_cand_neighbor_weight.clear();

  ifstream in_file(file_name);
  if (!in_file) {
    cout << "in_file error" << endl;
    exit(1);
  }
  in_file.seekg(0, in_file.end);
  size_t file_len = in_file.tellg();
  in_file.seekg(0, in_file.beg);
  char *data = new char[file_len];

  in_file.read(data, file_len);
  in_file.close();

  // skip comments
  char *pos = data;
  while (*pos == '%') {
    while (*(pos++) != '\n')
      ;
  }

  // read vertex_count
  int vertex_count = 0, edge_count = 0;
  while (*pos < '0' || *pos > '9') {
    ++pos;
  }
  while (*pos != ' ') {
    vertex_count = vertex_count * 10 + *pos - '0';
    ++pos;
  }
  // read edge_count
  while (*pos < '0' || *pos > '9') {
    ++pos;
  }
  while (*pos >= '0' && *pos <= '9') {
    edge_count = edge_count * 10 + *pos - '0';
    ++pos;
  }

  // read vertex_weight
  vector<int> vw;
  vw.resize(vertex_count + 1);
  int v, w;
  for (vector<vector<int>>::size_type i = 1; i < vw.size(); ++i) {
    v = w = 0;
    // read vertex
    while (*pos < '0' || *pos > '9') {
      ++pos;
    }
    while (*pos != ' ') {
      v = v * 10 + *pos - '0';
      ++pos;
    }
    // read weight
    while (*pos < '0' || *pos > '9') {
      ++pos;
    }
    while (*pos >= '0' && *pos <= '9') {
      w = w * 10 + *pos - '0';
      ++pos;
    }
    vw[v] = w;
  }

  // check solution weight
  long sw = 0;
  for (auto v : best_solution) {
    sw += vw[v];
  }
  if (sw != best_solution_weight) {
    cout << file_name << " wrong weight!" << endl;
    return false;
  }

  // read adjacency_list
  vector<vector<int>> adj_list;
  adj_list.resize(vertex_count + 1);

  vector<int> vertex_degree(vertex_count + 1, 0);
  char *stash_pos = pos;
  int v1, v2;
  for (int i = 0; i < edge_count; ++i) {
    v1 = v2 = 0;
    // read v1
    while (*pos < '0' || *pos > '9') {
      ++pos;
    }
    while (*pos != ' ') {
      v1 = v1 * 10 + *pos - '0';
      ++pos;
    }
    // read v2
    while (*pos < '0' || *pos > '9') {
      ++pos;
    }
    while (*pos >= '0' && *pos <= '9') {
      v2 = v2 * 10 + *pos - '0';
      ++pos;
    }
    vertex_degree[v1]++;
    vertex_degree[v2]++;
  }
  for (size_t v = 1; v < adj_list.size(); ++v) {
    adj_list[v].reserve(vertex_degree[v]);
  }
  pos = stash_pos;
  for (int i = 0; i < edge_count; ++i) {
    v1 = v2 = 0;
    // read v1
    while (*pos < '0' || *pos > '9') {
      ++pos;
    }
    while (*pos != ' ') {
      v1 = v1 * 10 + *pos - '0';
      ++pos;
    }
    // read weight
    while (*pos < '0' || *pos > '9') {
      ++pos;
    }
    while (*pos >= '0' && *pos <= '9') {
      v2 = v2 * 10 + *pos - '0';
      ++pos;
    }
    adj_list[v1].push_back(v2);
    adj_list[v2].push_back(v1);
  }
  delete[] data;

  // check clique
  for (vector<int>::size_type i = 0; i < best_solution.size(); ++i) {
    int v1 = best_solution[i];
    for (vector<int>::size_type j = i + 1; j < best_solution.size(); ++j) {
      int v2 = best_solution[j];
      vector<int>::size_type k = 0;
      for (; k < adj_list[v1].size(); ++k) {
        if (v2 == adj_list[v1][k]) {
          break;
        }
      }
      if (k >= adj_list[v1].size()) {
        cerr << file_name << " wrong anser: " << v1 << " is not neighbor of "
             << v2 << "!" << endl;
        return false;
      }
    }
  }
  return true;
}

long simp_count = 0;

bool better_since_simp = true;
// bool is_colored = false;
int main(int argc, char const *argv[]) {
  if (argc != 6) {
    std::cout << "usage: ./fast-wclq <instance> <seed> <cutoff> "
                 "<min_bms_count> <max_bms_count>"
              << std::endl;
    return 0;
  }

  int seed;
  long maxtries = 2000000000;
  int cutoff_time, i = 2;
  bool exact = false;

  sscanf(argv[i++], "%d", &seed);
  sscanf(argv[i++], "%d", &cutoff_time);
  sscanf(argv[i++], "%d", &min_bms_count);
  sscanf(argv[i++], "%d", &max_bms_count);

  srand(seed);
  t = 2;
  size_threshold = 300000;

  // clock_t start = clock();
  times(&start);

  fast_build(argv[1]);

  init();

  output_graph_size();

  best_solution_try = 0;
  best_solution_time = get_time();

  if (remaining_vertex.empty()) {
    exact = true;
    maxtries = 0;
  }

  if (best_solution_weight != 0) {
    cout << seed << ' ' << best_solution_weight << ' ' << best_solution_time
         << ' ' << best_solution_try << " h ";
    cout << tries << " " << simp_count << endl;
  }
  int level = 0;
  int reuse = false;
  for (tries = 1; tries <= maxtries; tries++) {

    // if (tries % 1000000 == 0)
    //   srand(seed++);

    if (get_time() > cutoff_time)
      break;

    int is_improve = construct2();

    if (is_improve == 1) {
      cout << seed << ' ' << best_solution_weight << ' ' << best_solution_time
           << ' ' << best_solution_try << " h ";
      cout << tries << " " << simp_count << endl;

      better_since_simp = true;
      level = 0;
      reuse = false;
    }

    is_new_graph = false;

    long branch_tries_thredshold = 100000;
    long color_tries_thredshold = 100000;

    switch (level) {
    case 0:
      if (simpilfy()) {
        is_new_graph = true;
      }
      level++;
      break;

    case 1:
      if (tries < branch_tries_thredshold) {
        break;
      }
      // if (remaining_vertex.size() > (unsigned long)size_threshold) {
      //   level++;
      //   break;
      // }

      if (simplify_branch()) {
        is_new_graph = true;
      }
      level++;
      break;

    case 2:
      if (tries < color_tries_thredshold) {
        break;
      }
      // if (remaining_vertex.size() > (unsigned long)size_threshold) {
      //   level++;
      //   break;
      // }
      if (simplify_branch_color()) {
        is_new_graph = true;
      } else
        level++;
      reuse = true;
      break;

    default:
      level++;
    }

    if (is_new_graph) {
      // is_colored = false;
      better_since_simp = false;
      simp_count++;

      if (remaining_vertex.empty()) {
        exact = true;
        break;
      }
      output_graph_size();
    }
  }

  cout << seed << ' ' << best_solution_weight << ' ' << best_solution_time
       << ' ' << best_solution_try;
  if (exact)
    cout << " x";
  else
    cout << " h";

  cout << " " << tries << " " << simp_count;
  cout << endl;

  verify_fast(argv[1]);

  return 0;
}
