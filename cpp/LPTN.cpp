#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>
#include <cmath>
#pragma omp declare reduction( + : arma::mat : omp_out += omp_in ) \
  initializer( omp_priv = omp_orig )

#pragma omp declare reduction( + : arma::vec : omp_out += omp_in ) \
  initializer( omp_priv = omp_orig )


#pragma omp declare reduction( + : arma::rowvec : omp_out += omp_in ) \
  initializer( omp_priv = omp_orig )  


// [[Rcpp::plugins(openmp)]]

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/gamma.hpp>


using namespace arma;


#include <sstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>


// [[Rcpp::export]]
arma::mat outer_product(arma::mat & X){

  //Dimensions of X
  int n = X.n_rows;
  int m = X.n_cols;

  arma::mat out(m,m, fill::zeros);

  for(int i = 0; i < n; i++){

    out += X.row(i).t()  * X.row(i);

  }

  return out;

}

// [[Rcpp::export]]
arma::vec row_sums(arma::mat & X){

  int m = X.n_cols;

  arma::vec out(m, fill::zeros);

  for(int j = 0; j < m; j++){

    out(j) = arma::accu(X.col(j));

  }

  return out;

}



// [[Rcpp::export]]
double cdf_chisq(double x, double df){

  return boost::math::gamma_p(0.5 * df, 0.5 * x);

}



//The chi-squared distributed quantities needed to determine if y is an outlier
// [[Rcpp::export]]
arma::colvec chi_squared_sums(arma::mat & X, arma::vec & mu, arma::mat & precision_matrix, int ncores){

  int n = X.n_rows;
  int m = X.n_cols;

  //Dummy matrix for outer products
  arma::mat dummy(m, m, fill::zeros);

  //New mu-vectors of appropriate dimentions for the matrix multiplications
  arma::rowvec mu_row(m, fill::zeros);
  arma::colvec mu_col(m, fill::zeros);

  for(int j = 0; j < m; j++){

    mu_row(j) += mu(j);
    mu_col(j) += mu(j);

  }

  arma::colvec out(n);
  omp_set_num_threads(ncores);
  #pragma omp parallel for 
  for(int i = 0; i < n; i++){

    out(i) = arma::as_scalar((X.row(i) - mu_row) * precision_matrix * (X.row(i).t() - mu_col));

  } 

  return out;  

}



// [[Rcpp::export]]
Rcpp::List eigen_sym(arma::mat & W, arma::mat & sd_error){

  arma::mat vcov_mat = W * W.t() + sd_error * sd_error;
  arma::mat precision_matrix = arma::inv(vcov_mat);

  arma::vec eigval;
  arma::mat eigvec;

  arma::eig_sym(eigval, eigvec, precision_matrix);

  return Rcpp::List::create(Rcpp::Named("values") = eigval,
                            Rcpp::Named("vectors") = eigvec,
                            Rcpp::Named("vcov_matrix") = vcov_mat,
                            Rcpp::Named("precision_matrix") = precision_matrix);

}



// [[Rcpp::export]]
arma::rowvec weighted_mu(arma::mat & X, arma::vec & weights){

  int n = X.n_rows;
  int m = X.n_cols;

  arma::rowvec out(m, fill::zeros);

  for(int i = 0; i < n; i++){

    out += arma::as_scalar(weights(i)) * X.row(i);

  }

  return(out);

}



// [[Rcpp::export]]
arma::mat weighted_vcov(arma::mat & X, arma::vec mu, arma::vec & weights){

  int n = X.n_rows;
  int m = X.n_cols;

  arma::rowvec mu_row(m, fill::zeros);

  for(int j = 0; j < m; j++){

    mu_row(j) += mu(j);

  }

  arma::mat out(m, m, fill::zeros);

  for(int i = 0; i < n; i++){

    out += arma::as_scalar(weights(i)) * (X.row(i) - mu_row).t()  * (X.row(i) - mu_row);

  }

  return(out);

}





//Function to compute the Gaussian LPTN log-likelihood
// [[Rcpp::export]]
double Gaussian_LPTN_log_likelihood_2(arma::mat & X, arma::vec & mu, arma::mat & covariance_matrix, double tau, int ncores){


  //---------------------------------------------------------------------------------------
  //Necessary matrices, determinant and constants that need to be computed wrt to the arguments

  arma::mat precision_matrix = covariance_matrix.i();
  double det_covariance_matrix = 1 / arma::det(precision_matrix);
  double F_chisq = cdf_chisq(std::pow(tau, 2), X.n_cols);

  //---------------------------------------------------------------------------------------
  //Variables to determine if the row of X is an outlier
  arma::colvec chi_squared = chi_squared_sums(X, mu, precision_matrix, ncores);

  //The log-likelihood
  double LL = 0;

  //Dimensions of X
  int n = X.n_rows;
  int m = X.n_cols;

  //The scaling constant of the distribution
  double tau_sq = std::pow(tau, 2);
  double pi = datum::pi;

  //beta normalizing constant
  double beta = m * std::exp(-0.5 * tau_sq) * std::pow(tau, m) * std::log(tau) / (std::pow(2, 0.5*m) * std::tgamma(0.5*m + 1) * (1 - F_chisq));

  //double c0 = 2 * F_chisq - 1;
  //double c1 = std::pow(pi_pow_m * det_covariance_matrix, -0.5) * std::exp(-0.5 * tau_sq);
  //double overall_scaling_constant = std::log(beta - 1) - std::log(c0 * (beta - 1) + 2 * tau * std::log(tau) * c1);

  //The scaling constant of the Gaussian
  double gaussian_scaling_constant = -0.5 * (m * std::log(2*pi) + std::log(det_covariance_matrix));

  //Every pdf includes both the Gaussian the Overall scaling constants
  //double scaling_term_sum = n * (overall_scaling_constant + gaussian_scaling_constant);
  double scaling_term_sum = n * gaussian_scaling_constant;

  //When an observation is an outlier, one needs to multiply by -0.5 * tau^2
  double outlier_log_exp_term = -0.5 * tau_sq;

  //Accumulate the values
  double x;
  omp_set_num_threads(ncores);
  #pragma omp parallel for reduction(+:LL) private(x)
  for(int i = 0; i < n; i++){

    if(chi_squared(i) <= tau_sq){

      LL += -0.5 *chi_squared(i);

    } else {

      x = std::sqrt(chi_squared(i));
      LL += outlier_log_exp_term + m * std::log(tau) - m * std::log(x) + (beta + 1) * (std::log(std::log(tau)) - std::log(std::log(x)));

    }

  }

  return scaling_term_sum + LL;

}





//Function to compute the Gradient of the LPTN log-likelihood wrt W, mu and tau_inv 
//Format: c(W, mu, tau_inv)
// [[Rcpp::export]]
Rcpp::List grad_Gaussian_LPTN_log_likelihood(arma::mat & X, arma::vec & mu, arma::mat & vcov_mat,
                                            double tau, int ncores){

  //---------------------------------------------------------------------------------------
  //Necessary matrices, determinant and constants that need to be computed wrt to the arguments

  arma::mat precision_matrix = vcov_mat.i();
  double F_chisq = cdf_chisq(std::pow(tau, 2), X.n_cols);


  //---------------------------------------------------------------------------------------
  //Placeholder for gradients 
  arma::mat grad_wrt_vcov_mat(vcov_mat.n_rows, vcov_mat.n_cols, fill::zeros);
  arma::vec grad_wrt_mu(X.n_cols, fill::zeros);


  //---------------------------------------------------------------------------------------
  //Variables to determine if the row of X is an outlier
  arma::colvec chi_squared = chi_squared_sums(X, mu, precision_matrix, ncores);

  //Dimensions of X
  int n = X.n_rows;
  int m = X.n_cols;

  //Dummy matrix for outer products
  arma::mat outer_product_matrix(m, m, fill::zeros);
  arma::mat outer_product_matrix_outliers(m, m, fill::zeros);

  //New mu-vectors of appropriate dimentions for the matrix multiplications
  arma::rowvec mu_row(m, fill::zeros);
  arma::colvec mu_col(m, fill::zeros);

  for(int j = 0; j < m; j++){

    mu_row(j) += mu(j);
    mu_col(j) += mu(j);

  }


  //---------------------------------------------------------------------------------------
  //Commonly used constants
  //General
  double tau_sq = std::pow(tau, 2);

  //beta normalizing constant
  double beta = m * std::exp(-0.5 * tau_sq) * std::pow(tau, m) * std::log(tau) / (std::pow(2, 0.5*m) * std::tgamma(0.5*m + 1) * (1 - F_chisq));

  //Distribution scaling factor

  //---------------------------------------------------------------------------------------
  //GRADIENTS 
  //---------------------------------------------------------------------------------------

  //Gradient for the non-outliers
  omp_set_num_threads(ncores);
  #pragma omp parallel for reduction(+:outer_product_matrix, outer_product_matrix_outliers, grad_wrt_mu)
  for(int i = 0; i < n; i++){

    double x;
    double derivative_wrt_tau;

    if(chi_squared(i) <= tau_sq){

      //Add to outer product matrix
      outer_product_matrix += (X.row(i) - mu_row).t()  * (X.row(i) - mu_row);

      //Add the quantity needed for the gradient wrt mu
      //grad_wrt_mu += precision_matrix * (X.row(i).t() - mu);
      grad_wrt_mu += (X.row(i).t() - mu);

    } else {

      //Add to outer product matrix
      //Multiply by the derivative of the function wrt (x - mu) * precision * (x - mu) so to use the chain rule
      x = std::sqrt(chi_squared(i));
      derivative_wrt_tau = (-m * std::pow(x, -1) -(beta + 1.0) / (x * std::log(x))) * (0.5) * std::pow(x, -1);
      //outer_product_matrix_outliers += derivative_wrt_tau * (-2.0) * (X.row(i) - mu_row).t()  * (X.row(i) - mu_row);
      outer_product_matrix_outliers += derivative_wrt_tau * (X.row(i) - mu_row).t()  * (X.row(i) - mu_row);

      //Add the quantity needed for the gradient wrt mu
      grad_wrt_mu += (-2) * derivative_wrt_tau * (X.row(i).t() - mu);

    }

  }

  grad_wrt_vcov_mat = precision_matrix * (outer_product_matrix - 2 * outer_product_matrix_outliers) * precision_matrix - n * precision_matrix;

  grad_wrt_mu = precision_matrix * grad_wrt_mu;

  return Rcpp::List::create(Rcpp::Named("mu") = grad_wrt_mu,
                            Rcpp::Named("vcov_matrix") = grad_wrt_vcov_mat);

}



//Fit the LPTN distribution
// [[Rcpp::export]]
Rcpp::List fit_LPTN(arma::mat & X, arma::vec & mu, arma::mat & vcov_mat,
                        double tau, double tol = 0.01, int max_iter = 100, int ncores = 1){

  //Dimensions of X
  int n = X.n_rows;
  int m = X.n_cols;

  arma::vec mu_before(m, fill::zeros);
  arma::vec mu_after(m, fill::zeros);
  for(int i = 0; i < m; i++){

    mu_before(i) = mu(i);

  }

  arma::rowvec mu_temp(m, fill::zeros);  

  double mu_diff = 0.0;

  arma::mat vcov_before(m, m, fill::zeros);
  arma::mat vcov_after(m, m, fill::zeros);
  for(int j = 0; j < m; j++){

    vcov_before.col(j) = vcov_mat.col(j);

  }


  double vcov_diff = 0.0;  

  arma::vec weights(n, fill::ones);

  arma::mat precision_matrix(m, m, fill::zeros);
  double det_covariance_matrix = 1.0;
  arma::colvec chi_squared(n, fill::zeros);

  double F_chisq = cdf_chisq(std::pow(tau, 2), X.n_cols);
  double tau_sq = std::pow(tau, 2);
  double beta = m * std::exp(-0.5 * tau_sq) * std::pow(tau, m) * std::log(tau) / (std::pow(2, 0.5*m) * std::tgamma(0.5*m + 1) * (1 - F_chisq));
  double x = 0.0;

  double mu_w = m / (m + m * (m + 1) / 2);
  double new_diff = 0.0;
  double old_diff = 0.0;

  for(int k = 0; k < max_iter; k++){

    precision_matrix = vcov_before.i();
    det_covariance_matrix = 1 / arma::det(precision_matrix);
    chi_squared = chi_squared_sums(X, mu_before, precision_matrix, ncores);

    for(int i = 0; i < n; i++){

      if(chi_squared(i) > tau_sq){

        x = std::sqrt(chi_squared(i));
        weights(i) = (1 / x) * ((m / x) + (beta + 1)/(x * std::log(x)));

      }

    }

    mu_temp = weighted_mu(X, weights) / n;
    for(int i = 0; i < m; i++){

      mu_after(i) = mu_temp(i);

    }

    chi_squared = chi_squared_sums(X, mu_after, precision_matrix, ncores);
    weights.fill(1.0);

    for(int i = 0; i < n; i++){

      if(chi_squared(i) > tau_sq){

        x = std::sqrt(chi_squared(i));
        weights(i) = (1 / x) * ((m / x) + (beta + 1)/(x * std::log(x)));

      }

    }

    vcov_after =  weighted_vcov(X, mu_after, weights) / n;

    mu_diff = std::sqrt(arma::accu(arma::pow(mu_after - mu_before, 2)));
    vcov_diff = std::sqrt(arma::accu(arma::pow(vcov_after - vcov_before, 2)));
    weights.fill(1.0);

    new_diff = mu_w * mu_diff + (1 - mu_w) * vcov_diff;
    if(k > 0){

      if(std::abs(new_diff/old_diff - 1) < tol){

        break;

      }

    }

    for(int j = 0; j < m; j++){

      mu_before(j) = mu_after(j);
      vcov_before.col(j) = vcov_after.col(j);

    }

    old_diff = new_diff;

  }


  double LL = Gaussian_LPTN_log_likelihood_2(X, mu_after, vcov_after, tau, ncores);

  return Rcpp::List::create(Rcpp::Named("mu") = mu_after,
                            Rcpp::Named("vcov_matrix") = vcov_after,
                            Rcpp::Named("LL") = LL); 


}




//Fit the LPTN distribution
// [[Rcpp::export]]
double fit_LPTN_2(arma::mat & X, arma::vec & mu, arma::mat & vcov_mat,
                        double tau, double tol = 0.05, int max_iter = 100, int ncores = 1){

  //Dimensions of X
  int n = X.n_rows;
  int m = X.n_cols;

  arma::vec mu_before(m, fill::zeros);
  arma::vec mu_after(m, fill::zeros);
  for(int i = 0; i < m; i++){

    mu_before(i) = mu(i);

  }

  arma::rowvec mu_temp(m, fill::zeros);  

  double mu_diff = 0.0;

  arma::mat vcov_before(m, m, fill::zeros);
  arma::mat vcov_after(m, m, fill::zeros);
  for(int j = 0; j < m; j++){

    vcov_before.col(j) = vcov_mat.col(j);

  }


  double vcov_diff = 0.0;  

  arma::vec weights(n, fill::ones);

  arma::mat precision_matrix(m, m, fill::zeros);
  double det_covariance_matrix = 1.0;
  arma::colvec chi_squared(n, fill::zeros);

  double F_chisq = cdf_chisq(std::pow(tau, 2), X.n_cols);
  double tau_sq = std::pow(tau, 2);
  double beta = m * std::exp(-0.5 * tau_sq) * std::pow(tau, m) * std::log(tau) / (std::pow(2, 0.5*m) * std::tgamma(0.5*m + 1) * (1 - F_chisq));
  double x = 0.0;

  double mu_w = m / (m + m * (m + 1) / 2);
  double new_diff = 0.0;
  double old_diff = 0.0;

  for(int k = 0; k < max_iter; k++){

    precision_matrix = vcov_before.i();
    det_covariance_matrix = 1 / arma::det(precision_matrix);
    chi_squared = chi_squared_sums(X, mu_before, precision_matrix, ncores);

    for(int i = 0; i < n; i++){

      if(chi_squared(i) > tau_sq){

        x = std::sqrt(chi_squared(i));
        weights(i) = (1 / x) * ((m / x) + (beta + 1)/(x * std::log(x)));

      }

    }

    mu_temp = weighted_mu(X, weights) / n;
    for(int i = 0; i < m; i++){

      mu_after(i) = mu_temp(i);

    }

    chi_squared = chi_squared_sums(X, mu_after, precision_matrix, ncores);
    weights.fill(1.0);

    for(int i = 0; i < n; i++){

      if(chi_squared(i) > tau_sq){

        x = std::sqrt(chi_squared(i));
        weights(i) = (1 / x) * ((m / x) + (beta + 1)/(x * std::log(x)));

      }

    }

    vcov_after =  weighted_vcov(X, mu_after, weights) / n;

    mu_diff = std::sqrt(arma::accu(arma::pow(mu_after - mu_before, 2)));
    vcov_diff = std::sqrt(arma::accu(arma::pow(vcov_after - vcov_before, 2)));
    weights.fill(1.0);

    new_diff = mu_w * mu_diff + (1 - mu_w) * vcov_diff;
    if(k > 0){

      if(std::abs(new_diff/old_diff - 1) < tol){

        break;

      }

    }

    for(int j = 0; j < m; j++){

      mu_before(j) = mu_after(j);
      vcov_before.col(j) = vcov_after.col(j);

    }

    old_diff = new_diff;

  }

  return Gaussian_LPTN_log_likelihood_2(X, mu_after, vcov_after, tau, ncores);

}





//Test multiple lambda values
// [[Rcpp::export]]
arma::vec get_LL_values_from_tau_vec(arma::mat & X, arma::vec & mu, arma::mat & vcov_mat,
                        arma::vec tau, double tol = 0.025, int max_iter = 100, int ncores = 1){


  int n_out = tau.n_elem;
  arma::vec out(n_out, fill::zeros);

  omp_set_num_threads(ncores);
  #pragma omp parallel for 
  for(int i = 0; i < n_out; i++){

    out(i) = fit_LPTN_2(X, mu, vcov_mat, tau(i), tol = tol, max_iter = max_iter);

  } 

  return out;


}
