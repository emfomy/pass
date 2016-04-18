////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @file   genlog/pass.cpp
/// @brief  The main PaSS algorithm for general logistic regression
///
/// @author Mu Yang <emfomy@gmail.com>
///

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @mainpage Particle Swarm Stepwise (PaSS) Algorithm for General Logistic Regression
///
/// @code{.unparsed}
/// ==========================================================================================================================
///
/// Notation:
/// X     : the regressors
/// Y     : the regressand
/// Beta  : the effects
/// P     : the probability of Y=1
/// Theta : the logit function of P
/// lv    : the likelihood value
/// llv   : the log-likelihood value
///
/// ==========================================================================================================================
///
/// Logistic model:
/// P     := exp(X*Beta) ./ 1+exp(X*Beta)
/// Theta := logit(P) = X*Beta
///
/// Update Beta:
/// Theta := X * Beta
/// Eta   := exp(Theta)
/// P     := Eta ./ (1+Eta)
/// 1-P    = 1 ./ (1+Eta)
/// W     := P .* (1-P)
/// Beta  += inv( X'*diag(W)*X ) * X' * (Y-P)
///
/// Compute log-likelihood:
/// lv    := prod( p^y * (1-p)^(1-y) )
/// llv   := log( lv )
///        = sum( y * log(p) + (1-y) * log(1-p) )
///        = sum( y * log(eta) - y * log(1+eta) - (1-y) * log(1+eta) )
///        = Y' * Theta - sum( log(1+eta) )
///
/// Newton-Raphson method:
/// d(llv)/d(Beta)     =   X' * (Y-P)
/// d^2(llv)/d(Beta)^2 = - X' * diag(W) * X
///
/// ==========================================================================================================================
///
/// Select the index in forward step:
/// idx = argmax_{i not in I} llv_hat
/// Theta_hat := Theta_new - Theta = Beta[i] * X[i col]
/// Eta_hat   := Eta_new  ./ Eta   = exp( Theta_hat )
/// llv_hat   := llv_new - llv
///            = ( Y' * Theta_new - Y' * Theta )
///              - ( sum( log(1+eta_new) ) - sum( log(1+eta) ) )
///            = Y' * Theta_hat - sum( log( (1+eta_new)/(1+eta) ) )
///            = Y' * Theta_hat - sum( log( 1 + (eta_hat-1)*p ) )
///
/// Approximate Beta[i] with Newton-Raphson method:
/// Beta[i]   += ( X[i col]'*(Y-P_new) ) / ( X[i col]'*diag(W_new)*X[i col] )
///
/// ==========================================================================================================================
///
/// Select the index in backward step:
/// idx = argmax_{i in I} llv_hat
/// Theta_hat := Theta_new - Theta = -Beta[i] * X[i col]
/// Eta_hat   := Eta_new  ./ Eta   = exp( Theta_hat )
/// llv_hat   := llv_new - llv
///            = ( Y' * Theta_new - Y' * Theta )
///              - ( sum( log(1+eta_new) ) - sum( log(1+eta) ) )
///            = Y' * Theta_hat - sum( log( (1+eta_new)/(1+eta) ) )
///            = Y' * Theta_hat - sum( log( 1 + (eta_hat-1)*p ) )
///
/// ==========================================================================================================================
///
/// References:
/// Chen, R.-B., Huang, C.-C., & Wang, W. (2015). Particle Swarm Stepwise (PaSS) Algorithm for Variable Selection.
///
/// Liu, Z., & Liu, M. (2011). Logistic Regression Parameter Estimation Based on Parallel Matrix Computation.
///   In Q. Zhou (Ed.), Communications in Computer and Information Science (Vol. 164, pp. 268–275).
///   Berlin, Heidelberg: Springer Berlin Heidelberg. http://doi.org/10.1007/978-3-642-24999-0_38
///
/// Singh, S., Kubica, J., Larsen, S., & Sorokina, D. (2013).
///   Parallel Large Scale Feature Selection for Logistic Regression (pp. 1172–1183).
///   Philadelphia, PA: Society for Industrial and Applied Mathematics. http://doi.org/10.1137/1.9781611972795.100
///
/// Barbu, A., She, Y., Ding, L., & Gramajo, G. (2014). Feature Selection with Annealing for Big Data Learning.
///   http://arxiv.org/pdf/1310.288
///
/// ==========================================================================================================================
/// @endcode
///
/// @author  Mu Yang <emfomy@gmail.com>
///

#include "pass.hpp"
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <mkl.h>
#include <omp.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// The log-binomial function
///
/// @param   n  a number
/// @param   k  a number
///
/// @return     return the log-binomial value of n and k
///
static inline float lbinom( const int n, const int k ) {
  int i;
  return (lgammaf_r(n+1, &i) - lgammaf_r(n-k+1, &i) - lgammaf_r(k+1, &i));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  The namespace of PaSS
//
namespace pass {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  The namespace of general logistic regression
//
namespace genlog {

int n;                // scalar, the number of statistical units
int p;                // scalar, the number of total effects
float *X0;            // matrix, n by p, the regressors
float *Y0;            // vector, n by 1, the regressand
bool *I0;             // vector, 1 by p, the chosen indices
float phi0;           // scalar, the criterion value
Parameter parameter;  // the PaSS parameters

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// The PaSS algorithm for linear regression
///
/// @param[in]   pass::n          scalar, the number of statistical units
/// @param[in]   pass::p          scalar, the number of total effects
/// @param[in]   pass::X0         matrix, n by p, the regressors
/// @param[in]   pass::Y0         vector, n by 1, the regressand
/// @param[in]   pass::parameter  the PaSS parameters
///
/// @param[out]  pass::I0         vector, 1 by p, the chosen indices
/// @param[out]  pass::phi0       scalar, the criterion value
///
/// @note        Please call @c srand() before using this routine.
///
void GenLog() {
  // Check parameters
  auto num_thread = omp_get_max_threads();
  auto num_particle = num_thread * parameter.num_particle_thread;

  // ======== Normalize the original data ================================================================================= //

  if ( !parameter.is_normalized ) {
    // Normalize X0
    for ( auto j = 0; j < p; ++j ) {
      cblas_sscal(n, (1.0f/cblas_snrm2(n, X0+j*n, 1)), X0+j*n, 1);
    }

    // Normalize Y0
    cblas_sscal(n, (1.0f/cblas_snrm2(n, Y0, 1)), Y0, 1);
  }

  // ======== Run PaSS ==================================================================================================== //

  // Allocate particles
  auto particle = new Particle[num_particle];
  phi0 = INFINITY;

  // Use openMP parallel
  #pragma omp parallel
  {
    unsigned int tid = omp_get_thread_num();

    for ( auto j = tid; j < num_particle; j+=num_thread ) {
      // Initialize particles
      particle[j].InitializeModel();
      particle[j].ComputeCriterion();
      particle[j].phi_old = particle[j].phi;

      // Copy best model
      if ( phi0 > particle[j].phi ) {
        phi0 = particle[j].phi;
        memcpy(I0, particle[j].I, sizeof(bool) * p);
      }
    }

    #pragma omp barrier

    // Find best model
    for ( auto i = 1u; i < parameter.num_iteration; ++i ) {
      for ( auto j = tid; j < num_particle; j+=num_thread ) {
        // Update model
        int idx;
        particle[j].SelectIndex(idx);
        particle[j].UpdateModel(idx);
        particle[j].ComputeCriterion();

        // Check singularity
        if ( isnan(particle[j].phi) ) {
          particle[j].InitializeModel();
          particle[j].ComputeCriterion();
        }

        // Change status
        if ( particle[j].phi > particle[j].phi_old ) {
          particle[j].status = !particle[j].status;
        }
        if ( particle[j].k <= 1 ) {
          particle[j].status = true;
        }
        if ( particle[j].k >= n-1 || particle[j].k >= p-4 ) {
          particle[j].status = false;
        }

        particle[j].phi_old = particle[j].phi;

        // Copy best model
        if ( phi0 > particle[j].phi ) {
          phi0 = particle[j].phi;
          memcpy(I0, particle[j].I, sizeof(bool) * p);
        }
      }
    }
  }

  // ====================================================================================================================== //

  // Delete memory
  delete[] particle;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// The constructor
///
Particle::Particle() {
  X        = static_cast<float*>(mkl_malloc(n * n     * sizeof(float), 64));
  Y        = static_cast<float*>(mkl_malloc(n         * sizeof(float), 64));
  Beta     = static_cast<float*>(mkl_malloc(n         * sizeof(float), 64));
  Theta    = static_cast<float*>(mkl_malloc(n         * sizeof(float), 64));
  Eta      = static_cast<float*>(mkl_malloc(n         * sizeof(float), 64));
  P        = static_cast<float*>(mkl_malloc(n         * sizeof(float), 64));
  W        = static_cast<float*>(mkl_malloc(n         * sizeof(float), 64));
  M        = static_cast<float*>(mkl_malloc(n*(n+1)/2 * sizeof(float), 64));
  STemp    = static_cast<float*>(mkl_malloc(n         * sizeof(float), 64));

  Idx_lo   = static_cast<int*  >(mkl_malloc(n         * sizeof(int),   64));
  Idx_ol   = static_cast<int*  >(mkl_malloc(p         * sizeof(int),   64));
  Idx_temp = static_cast<int*  >(mkl_malloc(p         * sizeof(int),   64));
  I        = static_cast<bool* >(mkl_malloc(p         * sizeof(bool),  64));

  iseed    = rand();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// The destructor
///
Particle::~Particle() {
  mkl_free(X);
  mkl_free(Y);
  mkl_free(Beta);
  mkl_free(Theta);
  mkl_free(Eta);
  mkl_free(P);
  mkl_free(W);
  mkl_free(M);
  mkl_free(STemp);

  mkl_free(Idx_lo);
  mkl_free(Idx_ol);
  mkl_free(Idx_temp);
  mkl_free(I);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Initialize the model randomly
///
void Particle::InitializeModel() {
  InitializeModel(rand_r(&iseed) % p);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Initialize the model
///
/// @param[in]  idx  the index of the effect
///
void Particle::InitializeModel( const int idx ) {
  // Initialize size
  k = 0;

  // Initialize index
  memset(I, false, sizeof(bool) * p);

  // X[0 col] := 1.0
  for ( auto i = 0; i < n; ++i ) {
    X[i] = 1.0f;
  }

  // Y := Y0
  cblas_scopy(n, Y0, 1, Y, 1);

  // Set status
  status = true;

  // Insert effect
  UpdateModel(idx);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Update the model
///
/// @param[in]  idx  the index of the effect
///
void Particle::UpdateModel( const int idx ) {
  if ( status ) {  // forward step
    // Update size
    k++;

    // Update index
    I[idx] = true;
    Idx_lo[k] = idx;
    Idx_ol[idx] = k;

    // Set Xnew
    auto Xnew = X+k*n;

    // Insert new row of X
    cblas_scopy(n, X0+idx*n, 1, Xnew, 1);

    // insert Beta by zero
    Beta[k] = 0.0f;
  } else {  // Backward step
    // Update index
    I[idx] = false;

    // Find index
    auto j = Idx_ol[idx];

    // Copy index end to index j
    if ( j != k ) {
      cblas_scopy(n, X+k*n, 1, X+j*n, 1);
      Beta[j] = Beta[k];
      Idx_lo[j] = Idx_lo[k];
      Idx_ol[Idx_lo[j]] = j;
    }

    // Update size
    k--;
  }

  // Compute Beta
  ComputeBeta();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Compute Beta
///
void Particle::ComputeBeta() {
  auto kp = k+1;

  // ======== Find Beta using Newton-Raphson's method ===================================================================== //

  do {
    // Theta := X * Beta
    cblas_sgemv(CblasColMajor, CblasNoTrans, n, kp, 1.0f, X, n, Beta, 1, 0.0f, Theta, 1);

    // Eta := exp(Theta)
    vsExp(n, Theta, Eta);

    // P := Eta ./ (1+Eta)
    vsLinearFrac(n, Eta, Eta, 1.0f, 0.0f, 1.0f, 1.0f, P);

    // W := P .* (1-P)
    vsLinearFrac(n, P, Eta, 1.0f, 0.0f, 1.0f, 1.0f, W);

    // -------- Beta += inv(X'*diag(W)*X) * X' * (Y-P) ------------------------------------------------------------------ //

    // M := X' * diag(W) * X
    for ( auto i = 0; i < kp*(kp+1)/2; ++i ) {
      M[i] = 0.0f;
    }
    for ( auto i = 0; i < n; ++i ) {
      cblas_sspr(CblasColMajor, CblasLower, kp, W[i], X+i, n, M);
    }

    // Compute Cholesky decomposition of M
    LAPACKE_spptrf(LAPACK_COL_MAJOR, 'L', kp, M);

    // W := (Y-P)
    vsSub(n, Y, P, W);

    // STemp := X' * W
    cblas_sgemv(CblasColMajor, CblasTrans, n, kp, 1.0f, X, n, W, 1, 0.0f, STemp, 1);

    // Solve STemp = inv(M) * STemp
    cblas_stpsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, kp, M, STemp, 1);

    // Solve STemp = inv(M') * STemp
    cblas_stpsv(CblasColMajor, CblasLower, CblasTrans, CblasNonUnit, kp, M, STemp, 1);

    // Beta += STemp
    vsAdd(kp, Beta, STemp, Beta);

    // ------------------------------------------------------------------------------------------------------------------ //
  } while ( cblas_snrm2(kp, STemp, 1) > sqrt(kp) * 1e-4f );

  // ====================================================================================================================== //
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Select the index of the effect to update
///
/// @param[out]  idx  the index of the effect
///
void Particle::SelectIndex( int& idx ) {
  auto srand = static_cast<float>(rand_r(&iseed)) / RAND_MAX;

  if ( status ) {  // Forward step
    // Idx_temp[0~itemp] := I0 exclude I
    // Idx_temp[0~(p-k)] := complement of I
    int itemp = 0;
    for ( auto i = 0, j = p-k; i < p; ++i ) {
      if ( !I[i] ) {
        if ( I0[i] ) {
          Idx_temp[itemp] = i;
          itemp++;
        } else {
          j--;
          Idx_temp[j] = i;
        }
      }
    }

    int choose;
    if ( itemp ) {
      choose = (srand < parameter.prob_forward_best) + (srand < parameter.prob_forward_best+parameter.prob_forward_improve);
    } else {
      choose = srand < parameter.prob_forward_improve / (parameter.prob_forward_improve+parameter.prob_forward_random);
    }

    switch( choose ) {
      case 2: {  // Choose randomly from best model
        idx = Idx_temp[rand_r(&iseed) % itemp];
        break;
      }
      case 1: {  // Choose most improvement index
        auto llv_temp = -INFINITY;
        for ( auto i = 0; i < p-k; ++i ) {
          // beta := 0
          auto beta = 0.0f;

          // Xnew := X0[Idx_temp[i] col]
          auto Xnew = X0+Idx_temp[i]*n;

          // stemp := Xnew' * Y
          auto stemp = cblas_sdot(n, Xnew, 1, Y, 1);

          // ======== Solve Y = X * Beta for Beta ========================================================================= //

          auto beta_temp = 0.0f;
          do {
            // STemp(Theta_new) := Theta + beta * Xnew
            cblas_scopy(n, Theta, 1, STemp, 1);
            cblas_saxpy(n, beta, Xnew, 1, STemp, 1);

            // Eta(Eta_new) := exp(Theta_new)
            vsExp(n, STemp, Eta);

            // STemp(P_new) := Eta_new ./ (1+Eta_new)
            vsLinearFrac(n, Eta, Eta, 1.0f, 0.0f, 1.0f, 1.0f, STemp);

            // W(W_new) := P_new .* (1-P_new)
            vsLinearFrac(n, STemp, Eta, 1.0f, 0.0f, 1.0f, 1.0f, W);

            // beta += (Xnew'*(Y-P_new)) / (Xnew'*diag(W_new)*Xnew)
            vsMul(n, W, Xnew, W);
            beta_temp = (stemp - cblas_sdot(n, Xnew, 1, STemp, 1)) / cblas_sdot(n, Xnew, 1, W, 1);
            beta += beta_temp;
          } while ( beta_temp > 1e-4f );

          // ======== stemp := Y' * Theta_hat - sum( log( 1 + (eta_hat-1)*p ) ) =========================================== //

          // STemp(Theta_hat) := beta * Xnew
          cblas_saxpby(n, beta, Xnew, 1, 0.0f, STemp, 1);

          // stemp := Y' * Theta_hat
          stemp = cblas_sdot(n, Y, 1, STemp, 1);

          // STemp := log(1+(Eta_hat-1)*P)
          vsExpm1(n, STemp, Eta);
          vsMul(n, Eta, P, W);
          vsLog1p(n, W, STemp);

          // stemp += sum(STemp)
          for ( auto j = 0; j < n; ++j ) {
            stemp += STemp[j];
          }

          // ============================================================================================================== //

          // Check if this value is maximum
          if ( llv_temp < stemp ) {
            llv_temp = stemp;
            idx = Idx_temp[i];
          }
        }
        break;
      }
      case 0: {  // Choose randomly
        idx = Idx_temp[rand_r(&iseed) % (p-k)];
        break;
      }
    }
  } else {  // Backward step
    if ( srand < parameter.prob_backward_improve ) {  // Choose most improvement index
      auto llv_temp = INFINITY;
      for ( auto i = 1; i <= k; ++i ) {
        // ======== stemp := Y' * Theta_hat - sum( log( 1 + (eta_hat-1)*p ) ) ============================================= //

        // STemp(Theta_hat) := -Beta[i] * X[i]
        cblas_saxpby(n, -Beta[i], X+i*n, 1, 0.0f, STemp, 1);

        // stemp := Y' * Theta_hat
        auto stemp = cblas_sdot(n, Y, 1, STemp, 1);

        // STemp := log(1+(Eta_hat-1)*P)
        vsExpm1(n, STemp, Eta);
        vsMul(n, Eta, P, W);
        vsLog1p(n, W, STemp);

        // stemp += sum(STemp)
        for ( auto j = 0; j < n; ++j ) {
          stemp += STemp[j];
        }

          // ============================================================================================================== //

        // Check if this value is minimal
        if ( llv_temp > stemp ) {
          llv_temp = stemp;
          idx = Idx_lo[i];
        }
      }
    } else {  // Choose randomly
      idx = Idx_lo[rand_r(&iseed) % k + 1];
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Compute the criterion value
///
void Particle::ComputeCriterion() {
  // llv := Y' * Theta - sum(log(1+Eta))
  vsLog1p(n, Eta, STemp);
  llv = cblas_sdot(n, Y, 1, Theta, 1);
  for ( auto i = 0; i < n; ++i ) {
    llv += STemp[i];
  }

  // Compute criterion
  switch(parameter.criterion) {
    case AIC: {    // phi := n*log(e^2/n) + 2k
      phi = -2.0f*llv + 2.0f*k;
      break;
    }
    case BIC: {    // phi := n*log(e^2/n) + k*log(n)
      phi = -2.0f*llv + k*logf(n);
      break;
    }
    case EBIC: {   // phi := n*log(e^2/n) + k*log(n) + 2gamma*log(binom(p, k))
      phi = -2.0f*llv + k*logf(n) + 2.0f*parameter.ebic_gamma*lbinom(p, k);
      break;
    }
    case HDBIC: {  // phi := n*log(e^2/n) + k*log(n)*log(p)
      phi = -2.0f*llv + k*logf(n)*logf(p);
      break;
    }
    case HQC: {    // phi := n*log(e^2/n) + 2k*log(log(n))
      phi = -2.0f*llv + 2.0f*k*logf(logf(n));
      break;
    }
    case HDHQC: {  // phi := n*log(e^2/n) + 2k*log(log(n))*log(p)
      phi = -2.0f*llv + 2.0f*k*logf(logf(n))*logf(p);
      break;
    }
  }
}

}  // namespace genlog

}  // namespace pass
