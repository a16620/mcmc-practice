#include <Rcpp.h>
#include <queue>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame C_summary_tree(List x) {
  CharacterVector variable_vec;
  IntegerVector   depth_vec;
  
  size_t depth = 1;
  std::queue<List> queue;
  std::queue<List> queue_next;
  queue.push(x);

  while (!queue.empty()) {
    while (!queue.empty()) {
      List& node = queue.front();

      if (node.containsElementNamed("split.rule")) {
        NumericVector&& rule = as<NumericVector>(node["split.rule"]);
        variable_vec.push_back(as<CharacterVector>(rule.names()).at(0));
        depth_vec.push_back(depth);

        queue_next.push(node["left"]);
        queue_next.push(node["right"]);
      }
      queue.pop();
    }
    queue.swap(queue_next);
    depth++;
  }

  return DataFrame::create(Named("depth")=depth_vec, Named("variable")=variable_vec);
}

