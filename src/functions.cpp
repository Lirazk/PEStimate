#include <Rcpp.h>
#include <boost/math/special_functions/owens_t.hpp>

using namespace Rcpp;

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(Rcpp)]]
//' @export
// [[Rcpp::export]]
NumericVector inner_integral_vec_rcpp(NumericVector &x, NumericVector &a, double b)
{
  NumericVector result = NumericVector(x.size());
  double sqrt1pb2 = sqrt(1 + b*b);
  NumericVector t1 = a / (x * sqrt1pb2);
  NumericVector t2 = (x * sqrt1pb2) / a;
  NumericVector t3 = (a + b * x) / x;
  NumericVector t4 = (a * b + x * (1 + b*b)) / a;
  NumericVector o = a / sqrt1pb2;
  //result[0] = boost::math::owens_t(x[0], t1[0]);
  //return result;
  for (int i = 0; i < result.size(); i++) {
    result[i] = boost::math::owens_t(x[i], t1[i]) + boost::math::owens_t(o[i], t2[i]) - boost::math::owens_t(o[i], t4[i]) + R::pnorm(x[i], 0, 1, 1, 0) * R::pnorm(o[i], 0, 1, 1, 0);
    if (!std::isinf(x[i])) {
      result[i] -= boost::math::owens_t(x[i], t3[i]);
    }
  }
  return result;
}

// [[Rcpp::depends(Rcpp)]]
//' @export
// [[Rcpp::export]]
NumericVector selected_score_rcpp(NumericMatrix embryo_scores, int iter, double prs_threshold) {
  NumericVector selected_score(iter);
  
  for (int i = 0; i < iter; ++i) {
    NumericVector row_scores_with_na = embryo_scores(i, _);
    NumericVector row_scores = row_scores_with_na[!is_na(row_scores_with_na)];
    NumericVector eligible_scores;
    
    for (int j = 0; j < row_scores.size(); ++j) {
      if (row_scores[j] < prs_threshold) {
        eligible_scores.push_back(row_scores[j]);
      }
    }
    
    if (eligible_scores.size() > 0) {
      selected_score[i] = Rcpp::sample(eligible_scores, 1)[0];
    } else {
      selected_score[i] = Rcpp::sample(row_scores, 1)[0];
    }
  }
  
  return selected_score;
}
