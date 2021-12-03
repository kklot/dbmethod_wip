#include <TMB.hpp>
#include "ktools.hpp"

template<class Type>
vector<Type> rskewlogis (vector<Type> scale, Type shape, Type skew) {
  vector<Type> u_l(scale.size()), u_u(scale.size());
  u_l.fill(0.);
  u_u.fill(1.);
  vector<Type> u = runif(u_l, u_u);
  vector<Type> t1 = pow(u, -Type(1.) / skew) - Type(1.);
  vector<Type> t2 = pow(t1, -Type(1.) / shape);
  vector<Type> o = Type(1.) / scale *t2;
  return o;
}

template<class Type>
Type rskewlogis (Type scale, Type shape, Type skew) {
  Type u = runif(0., 1.);
  Type t1 = pow(u, -Type(1.) / skew) - Type(1.);
  Type t2 = pow(t1, -Type(1.) / shape);
  Type o = Type(1.) / scale * t2;
  return o;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  parallel_accumulator<Type> dll(this);

  // data
  DATA_VECTOR(afs);
  DATA_IVECTOR(age);
  DATA_IVECTOR(event);
  DATA_IVECTOR(yob);

  // priors
  DATA_VECTOR(sd_beta);
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
  Type skew = exp(log_skew);
  prior -= dnorm(log_skew, skew_prior(0), skew_prior(1), true);

  // yob rw2
  PARAMETER_VECTOR (yob_rw2);
  PARAMETER_VECTOR (yob_phi);

  PARAMETER(log_ar_precision);
  prior -= dgamma(exp(log_ar_precision), Type(1), Type(1./1e-5), true) + log_ar_precision;
  Type ar_scale_sd = pow(1 / exp(log_ar_precision), 0.5);

  density::ARk_t<Type> yob_ar2 = density::ARk(yob_phi);

  prior -= dnorm(log( (1 + yob_phi(0)) / (1 - yob_phi(0))), Type(0), Type(1), true);
  prior -= dnorm(log( (1 + yob_phi(1)) / (1 - yob_phi(1))), Type(0), Type(1), true);
  prior += density::SCALE(yob_ar2, ar_scale_sd)(yob_rw2);
  SIMULATE {
    density::SCALE(density::ARk(yob_phi), ar_scale_sd).simulate(yob_rw2);
  }

  // age rw2
  PARAMETER_VECTOR (age_rw2);
  PARAMETER        (age_phi);

  density::AR1_t<density::N01<Type> > age_ar1 = density::AR1(age_phi);

  prior -= dnorm(log( (1 + age_phi) / (1 - age_phi)), Type(0), Type(1), true);
  prior += age_ar1(age_rw2);
  SIMULATE {
    density::AR1(age_phi).simulate(age_rw2);
  }

  // Data likelihood
  vector<Type> predicted(afs.size()), scale(afs.size());
  for (int i = 0; i < afs.size(); i++)
  {
    Type eta = intercept + yob_rw2(yob(i)) + age_rw2(age(i));
    scale(i) = exp(eta);
    if (event(i)) {
      dll -= log(ktools::ft_llogisI(afs(i), shape, scale(i), skew));
    } else {
      dll -= log(ktools::St_llogisI(afs(i), shape, scale(i), skew));
    }
    SIMULATE {
      predicted(i) = rskewlogis(scale(i), shape, skew);
    }
  }

  dll += prior;

  SIMULATE {
    REPORT(scale);
    REPORT(predicted);
    REPORT(yob_rw2);
    REPORT(age_rw2);
    REPORT(skew);
    REPORT(shape);
  }
  return dll;
}
