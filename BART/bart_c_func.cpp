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

// [[Rcpp::export]]
NumericVector C_predict_single_set(DataFrame x, List tree_set) {
  const size_t X_len = x.nrow(), tree_len = tree_set.length();
  NumericVector pred_output(X_len);
  List ldf = as<List>(x);

  for (size_t i = 0; i < X_len; i++) {
    double y_pred = 0;
    for (size_t k = 0; k < tree_len; k++) {
      List cursor = tree_set.at(k);
      while (cursor.containsElementNamed("split.rule")) {
        NumericVector&& rule = as<NumericVector>(cursor["split.rule"]);
        String rcol = as<CharacterVector>(rule.names())[0];
        auto&& X_val = as<NumericVector>(ldf[rcol]).at(i);
        if (X_val <= rule.at(0)) {
          cursor = cursor["left"];
        } else {
          cursor = cursor["right"];
        }
      }
      y_pred += as<NumericVector>(cursor["node.param"]).at(0);
    }
    pred_output[i] = y_pred;
  }

  return pred_output;
}