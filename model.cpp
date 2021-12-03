#include <TMB.hpp>
#include "ktools.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);

  // data
  DATA_VECTOR(afs); DATA_UPDATE(afs); // afs is updatable
  DATA_VECTOR(afs_u);
  DATA_VECTOR(afs_l);
  DATA_IVECTOR(age);
  DATA_IVECTOR(event);
  DATA_IVECTOR(yob);

  DATA_VECTOR(svw);

  // meta data
  DATA_INTEGER(to_skew);
  DATA_INTEGER(age_term);
  DATA_INTEGER(yob_term);
  DATA_INTEGER(smooth_age);
  DATA_INTEGER(smooth_yob);
  DATA_INTEGER(interval_censor);
  DATA_IVECTOR(to_fit);

  // priors
  DATA_VECTOR(sd_beta);
  DATA_VECTOR(sd_yob);
  DATA_VECTOR(sd_age);

  DATA_VECTOR(shape_prior);
  DATA_VECTOR(skew_prior);

  // Data model - log-logistic parameters
  PARAMETER(intercept);

  Type prior = 0.0;
  prior -= dnorm(intercept, sd_beta(0), sd_beta(1), true);

  // Shape
  PARAMETER(log_shape); 
  prior -= dnorm(log_shape, shape_prior(0), shape_prior(1), true);
  Type shape = exp(log_shape);

  // Skewness *
  PARAMETER(log_skew);
  Type skew = 1.;
  if (to_skew) {
    prior -= dnorm(log_skew, skew_prior(0), skew_prior(1), true);
    skew = exp(log_skew);
  }

  // yob rw2
  PARAMETER        (beta_yob);
  PARAMETER_VECTOR (yob_rw2);
  PARAMETER_VECTOR (yob_phi);
  DATA_SCALAR(ar_scale);

  density::ARk_t<Type> yob_ar2 = density::ARk(yob_phi);

  if (yob_term)
  {
    if (smooth_yob) {
      prior -= dnorm(log( (1 + yob_phi(0)) / (1 - yob_phi(0))), Type(0), Type(1), true);
      prior -= dnorm(log( (1 + yob_phi(1)) / (1 - yob_phi(1))), Type(0), Type(1), true);
      prior += density::SCALE(yob_ar2, ar_scale)(yob_rw2);
    } else {
      prior -= dnorm(beta_yob, sd_beta(0), sd_beta(1), true);
    }
  }

  // age rw2
  PARAMETER        (beta_age);
  PARAMETER_VECTOR (age_rw2);
  PARAMETER        (age_phi);

  density::AR1_t<density::N01<Type> > age_ar1 = density::AR1(age_phi);

  if (age_term)
  {
    if (smooth_age) {
      prior -= dnorm(log( (1 + age_phi) / (1 - age_phi)), Type(0), Type(1), true);
      prior += age_ar1(age_rw2);
    } else {
      prior -= dnorm(beta_age, sd_beta(0), sd_beta(1), true);
    }
  }

  // Data likelihood
  vector<Type> logliki(afs.size());

  for (int i = 0; i < afs.size(); i++) {

    Type eta = intercept;

    if (yob_term) 
    {
      if (smooth_yob)
        eta += yob_rw2(yob(i));
      else
        eta += beta_yob * yob(i);
    }

    if (age_term)
    {
      if (smooth_age)
        eta += age_rw2(age(i));
      else
        eta += beta_age * age(i);
    }

    Type scale = exp(eta);

    if (event(i)) 
    {
      if (interval_censor)
        logliki(i) = log( svw(i) * (
              ktools::St_llogisI(afs_l(i), shape, scale, skew) - 
              ktools::St_llogisI(afs_u(i), shape, scale, skew) ));
      else
        logliki(i) = log( svw(i) * ktools::ft_llogisI(afs(i), shape, scale, skew) );
    } 
    else {
      logliki(i) = log(svw(i) * ktools::St_llogisI(afs(i), shape, scale, skew));
    }
    if (to_fit(i))
      dll -= logliki(i);
  }
  dll += prior;
  REPORT(logliki);
  return dll;
}
