/*
 * mosquito_ode.cpp
 *
 *  Created on: 11 Jun 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include "aquatic_mosquito_eqs.h"

integration_function_t create_eqs(AquaticMosquitoModel& model) {
  return [&model](const state_t& x, state_t& dxdt, double t) {
    auto kt = model.K_history.at(t, false); 
    auto K = carrying_capacity(
      t,
      model.model_seasonality,
      model.g0,
      model.g,
      model.h,
      kt,
      model.R_bar,
      model.rainfall_floor
    );
    
    auto beta = eggs_laid(model.beta, model.mum, model.f);
    auto n_larvae = x[get_idx(AquaticState::E)] + x[get_idx(AquaticState::L)];
    
    if(K == 0){
      // If carrying capacity 0 then all aquatic stages removed
      dxdt[get_idx(AquaticState::E)]  = - x[get_idx(AquaticState::E)];
      dxdt[get_idx(AquaticState::L)] = - x[get_idx(AquaticState::L)];
      dxdt[get_idx(AquaticState::P)] = - x[get_idx(AquaticState::P)];
    } else {
      dxdt[get_idx(AquaticState::E)] = beta * (model.total_M) //new eggs
      - x[get_idx(AquaticState::E)] / model.de //growth to late larval stage
      - x[get_idx(AquaticState::E)] * model.mue * (1 + n_larvae / K); //early larval deaths
      
      dxdt[get_idx(AquaticState::L)] = x[get_idx(AquaticState::E)] / model.de //growth from early larval
      - x[get_idx(AquaticState::L)] / model.dl //growth to pupal
      - x[get_idx(AquaticState::L)] * model.mul * (1 + model.gamma * n_larvae / K); //late larval deaths
      
      dxdt[get_idx(AquaticState::P)] = x[get_idx(AquaticState::L)] / model.dl //growth to pupae
      - x[get_idx(AquaticState::P)] / model.dp //growth to adult
      - x[get_idx(AquaticState::P)] * model.mup; // death of pupae
    }
  };
}

AquaticMosquitoModel::AquaticMosquitoModel(
  double beta,
  double de,
  double mue,
  History K_history,
  double gamma,
  double dl,
  double mul,
  double dp,
  double mup,
  size_t total_M,
  bool model_seasonality,
  double g0,
  std::vector<double> g,
  std::vector<double> h,
  double R_bar,
  double mum,
  double f,
  double rainfall_floor
):
  beta(beta),
  de(de),
  mue(mue),
  K_history(K_history),
  gamma(gamma),
  dl(dl),
  mul(mul),
  dp(dp),
  mup(mup),
  total_M(total_M),
  model_seasonality(model_seasonality),
  g0(g0),
  g(g),
  h(h),
  R_bar(R_bar),
  mum(mum),
  f(f),
  rainfall_floor(rainfall_floor)
{}



//[[Rcpp::export]]
Rcpp::XPtr<AquaticMosquitoModel> create_aquatic_mosquito_model(
    double beta,
    double de,
    double mue,
    double K0,
    bool carrying_capacity,
    std::vector<double> K,
    std::vector<double> K_tt,
    double gamma,
    double dl,
    double mul,
    double dp,
    double mup,
    size_t total_M,
    bool model_seasonality,
    double g0,
    std::vector<double> g,
    std::vector<double> h,
    double R_bar,
    double mum,
    double f,
    double rainfall_floor
) {
  
  // Create the carrying capacity object
  History K_history(K_tt.size() + 1, K0);
  K_history.push(K0, 0.0);
  if(carrying_capacity){
    for(int i = 0; i < K_tt.size(); ++i){
      Rcpp::Rcout << "Adding: " << K[i] << " at " << K_tt[i] << "\n";
      K_history.push(K[i], K_tt[i]);
    }
  }
  Rcpp::Rcout << "T1: " << K_history.at(1, false) << "\n";
  Rcpp::Rcout << "T200: " << K_history.at(200, false) << "\n";
  Rcpp::Rcout << "T364: " << K_history.at(364, false) << "\n";
  Rcpp::Rcout << "T365: " << K_history.at(365, false) << "\n";
  Rcpp::Rcout << "T366: " << K_history.at(366, false) << "\n";
  
  auto model = new AquaticMosquitoModel(
    beta,
    de,
    mue,
    K_history,
    gamma,
    dl,
    mul,
    dp,
    mup,
    total_M,
    model_seasonality,
    g0,
    g,
    h,
    R_bar,
    mum,
    f,
    rainfall_floor
  );
  return Rcpp::XPtr<AquaticMosquitoModel>(model, true);
}

//[[Rcpp::export]]
void aquatic_mosquito_model_update(
    Rcpp::XPtr<AquaticMosquitoModel> model,
    size_t total_M,
    double f,
    double mum
) {
  model->total_M = total_M;
  model->f = f;
  model->mum = mum;
}


//[[Rcpp::export]]
Rcpp::XPtr<Solver> create_aquatic_solver(
    Rcpp::XPtr<AquaticMosquitoModel> model,
    std::vector<double> init,
    double r_tol,
    double a_tol,
    size_t max_steps
) {
  return Rcpp::XPtr<Solver>(
    new Solver(init, create_eqs(*model), r_tol, a_tol, max_steps),
    true
  );
}
