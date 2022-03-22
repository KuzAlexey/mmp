// =====================================================================================
//
//       Filename:  basic.h
//
//    Description:  bas
//
//        Version:  1.0
//        Created:  2018年 07月 03日 星期二 11:12:39 CST
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Jinkun Lin, jkunlin@gmail.com
//   Organization:  School of EECS, Peking University
//
// =====================================================================================
#ifndef BASIC_H_KGHMXNZL
#define BASIC_H_KGHMXNZL

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <sys/times.h>
#include <unistd.h>
#include <vector>
using namespace std;

tms start, finish;
inline double get_time() {
  times(&finish);
  double res = double(finish.tms_utime - start.tms_utime + finish.tms_stime -
                      start.tms_stime) /
               sysconf(_SC_CLK_TCK);
  res = round(res * 100) / 100;
  return res;
}

vector<vector<int>> adjacency_list;
vector<vector<int>> adjacency_pos_list;
vector<long> vertex_neighbor_weight;
vector<int> vertex_weight;

// kcore
vector<vector<int>::size_type> index;
vector<int> vertex;
vector<int> tmp_degree;

// solutions
vector<int> solution;
vector<int> best_solution;
long best_solution_weight = 0;
long solution_weight;

long tries;
double best_solution_time;
long best_solution_try;

// input parameter
int size_threshold;
int t;

// bms
int start_bms_count = 1;
int min_bms_count;
int max_bms_count;
int real_bms_count;

struct Remaining_vertex {
  vector<int> vertex;
  vector<vector<int>::size_type> index;
  vector<int> tmp_vertex;
  vector<int> tmp_degree;

  vector<int>::iterator begin() { return vertex.begin(); }
  vector<int>::iterator end() { return vertex.end(); }
  void init(vector<int>::size_type vertex_size) {
    vertex.reserve(vertex_size);
    index.resize(vertex_size + 1);
    for (vector<int>::size_type i = 1; i <= vertex_size; ++i) {
      vertex.push_back(i);
      index[i] = i - 1;
    }
  }
  void init(vector<int> &vec, vector<vector<int>::size_type> &ind) {
    vertex.swap(vec);
    index.swap(ind);
  }
  void remove(int v) {
    index[*vertex.rbegin()] = index[v];
    vertex[index[v]] = *vertex.rbegin();
    vertex.pop_back();
  }

  vector<int>::size_type size() { return vertex.size(); }

  bool empty() { return vertex.empty(); }

  int &operator[](size_t i) { return vertex[i]; }
  void exchange(size_t i, size_t j) {
    swap(index[vertex[i]], index[vertex[j]]);
    swap(vertex[i], vertex[j]);
  }
};
Remaining_vertex remaining_vertex;

// construct
vector<vector<double>> adjacency_cand_neighbor_weight;
vector<bool> is_computed;

#endif /* end of include guard: BASIC_H_KGHMXNZL */
