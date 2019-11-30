////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @file   genlin/pass.cpp
/// @brief  The main PaSS algorithm for general linear regression
///
/// @author Mu Yang <emfomy@gmail.com>
///

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @mainpage Particle Swarm Stepwise (PaSS) Algorithm for General Linear Regression
///
/// @code{.unparsed}
/// ==========================================================================================================================
///
/// Notation:
/// X    : the regressors
/// Y    : the regressand
/// Beta : the effects
/// R    : the residual
///
/// ==========================================================================================================================
///
/// Linear model:
/// Y = X * Beta + error
/// R = Y - X * Beta
///
/// ==========================================================================================================================
///
/// Select the index in forward step:
/// idx = argmax_{i not in I} abs( X[i col]' * R  )
///
/// Select the index in backward step:
/// idx = argmin_{i in I} norm( R_{I exclude i} )
///     = argmin_{i in I} norm( R + Beta[i] * X[i col] )
///
/// ==========================================================================================================================
///
/// References:
/// Chen, R.-B., Huang, C.-C., & Wang, W. (2015). Particle Swarm Stepwise (PaSS) Algorithm for Variable Selection.
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
#include <omp.h>
#include <mkl.h>
#include <magma.h>

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

float *R0, *E0, *dX0, *dR0, *dE0;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  The namespace of PaSS
//
namespace pass {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  The namespace of general linear regression
//
namespace genlin {

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
void GenLin() {
  // Check parameters
  auto num_thread = omp_get_max_threads();
  auto num_particle = num_thread * parameter.num_particle_thread;

  // ======== Centralize and normalize the original data ================================================================== //

  if ( !parameter.is_normalized ) {
    // Centralize and normalize X0
    for ( auto j = 0; j < p; ++j ) {
      float stemp = 1.0f;
      #pragma omp simd
      for ( auto i = 0; i < n; ++i ) {
        stemp += X0[i+j*n];
      }
      stemp /= n;
      vsLinearFrac(n, X0+j*n, X0+j*n, 1.0f, -stemp, 0.0f, 1.0f, X0+j*n);
      cblas_sscal(n, (1.0f/cblas_snrm2(n, X0+j*n, 1)), X0+j*n, 1);
    }

    // Centralize and normalize Y0
    {
      float stemp = 1.0f;
      #pragma omp simd
      for ( auto i = 0; i < n; ++i ) {
        stemp += Y0[i];
      }
      stemp /= n;
      vsLinearFrac(n, Y0, Y0, 1.0f, -stemp, 0.0f, 1.0f, Y0);
      cblas_sscal(n, (1.0f/cblas_snrm2(n, Y0, 1)), Y0, 1);
    }
  }

  // ======== Run PaSS ==================================================================================================== //

  // Allocate particles
  phi0 = INFINITY;
  auto particle = new Particle[num_particle];

  // Allocate GPU memory
  magma_smalloc_cpu(&R0, n * num_particle);
  magma_smalloc_cpu(&E0, p * num_particle);
  magma_smalloc(&dX0, n * p);
  magma_smalloc(&dR0, n * num_particle);
  magma_smalloc(&dE0, p * num_particle);
  magma_ssetmatrix(n, p, X0, n, dX0, n);

  // Initialize particles
  #pragma omp parallel for
  for ( auto j = 0u; j < num_particle; ++j ) {
    particle[j].R = R0 + j*n;
    particle[j].E = E0 + j*p;
    particle[j].InitializeModel();
    particle[j].ComputeCriterion();
    particle[j].phi_old = particle[j].phi;

    // Copy best model
    if ( phi0 > particle[j].phi ) {
      phi0 = particle[j].phi;
      memcpy(I0, particle[j].I, sizeof(bool) * p);
    }
  }

  // Find best model
  for ( auto i = 1u; i < parameter.num_iteration; ++i ) {
    // E0 := X0' * R0
    magma_ssetmatrix(n, num_particle, R0, n, dR0, n);
    magma_sgemm(MagmaTrans, MagmaNoTrans, p, num_particle, n, 1.0, dX0, n, dR0, n, 0.0, dE0, p);
    magma_sgetmatrix(p, num_particle, dE0, p, E0, p);

    #pragma omp parallel for
    for ( auto j = 0u; j < num_particle; ++j ) {
      // Update model
      particle[j].SelectOption();
      particle[j].SelectIndex();
      particle[j].UpdateModel();
      particle[j].ComputeCriterion();

      // Check singularity
      if ( std::isnan(particle[j].phi) ) {
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
      if ( particle[j].k >= n || particle[j].k >= p-4 ) {
        particle[j].status = false;
      }

      particle[j].phi_old = particle[j].phi;
    }

    // Copy best model
    int j_best = -1;
    for ( auto j = 0u; j < num_particle; ++j ) {
      if ( phi0 > particle[j].phi ) {
        phi0 = particle[j].phi;
        j_best = j;
      }
    }
    if ( j_best != -1 ) {
      memcpy(I0, particle[j_best].I, sizeof(bool) * p);
    }
  }

  // ====================================================================================================================== //

  // Delete memory
  delete[] particle;
  magma_free_cpu(R0);
  magma_free_cpu(E0);
  magma_free(dX0);
  magma_free(dR0);
  magma_free(dE0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// The constructor
///
Particle::Particle() {
  X        = static_cast<float*>(mkl_malloc(n * n * sizeof(float), 64));
  Y        = static_cast<float*>(mkl_malloc(n     * sizeof(float), 64));
  Beta     = static_cast<float*>(mkl_malloc(n     * sizeof(float), 64));
  Theta    = static_cast<float*>(mkl_malloc(n     * sizeof(float), 64));
  M        = static_cast<float*>(mkl_malloc(n * n * sizeof(float), 64));
  B        = static_cast<float*>(mkl_malloc(n     * sizeof(float), 64));
  D        = static_cast<float*>(mkl_malloc(n     * sizeof(float), 64));

  Idx_lo   = static_cast<int*  >(mkl_malloc(n     * sizeof(int),   64));
  Idx_ol   = static_cast<int*  >(mkl_malloc(p     * sizeof(int),   64));
  Idx_best = static_cast<int*  >(mkl_malloc(p     * sizeof(int),   64));
  I        = static_cast<bool* >(mkl_malloc(p     * sizeof(bool),  64));

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
  mkl_free(M);
  mkl_free(B);
  mkl_free(D);

  mkl_free(Idx_lo);
  mkl_free(Idx_ol);
  mkl_free(Idx_best);
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
  k = 1;

  // Initialize index
  memset(I, false, sizeof(bool) * p);

  // Update index
  I[idx] = true;
  Idx_lo[0] = idx;
  Idx_ol[idx] = 0;

  // X := X0[idx col]
  cblas_scopy(n, X0+idx*n, 1, X, 1);

  // Y := Y0
  cblas_scopy(n, Y0, 1, Y, 1);

  // M := 1 / (X' * X)
  M[0] = 1.0f;

  // Theta := X' * Y
  Theta[0] = cblas_sdot(n, X, 1, Y, 1);

  // Beta := M * Theta
  Beta[0] = M[0] * Theta[0];

  // R = Y - X * Beta
  cblas_scopy(n, Y, 1, R, 1);
  cblas_saxpy(n, -Beta[0], X, 1, R, 1);

  // Set status
  status = true;
  option = 0;
  r_computed = false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Update the model
///
/// @param[in]  idx  the index of the effect
///
void Particle::UpdateModel( const int idx ) {
  this->idx = idx;
  UpdateModel();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Update the model
///
void Particle::UpdateModel() {
  r_computed = false;

  if ( status ) {  // forward step
    // Update size
    auto km = k;
    k++;

    // Update index
    I[idx] = true;
    Idx_lo[km] = idx;
    Idx_ol[idx] = km;

    // Set Xnew
    auto Xnew = X+km*n;

    // Insert new row of X
    cblas_scopy(n, X0+idx*n, 1, Xnew, 1);

    // ======== Solve Y = X * Beta for Beta =============================================================================== //

    // B := X' * Xnew
    cblas_sgemv(CblasColMajor, CblasTrans, n, km, 1, X, n, Xnew, 1, 0.0f, B, 1);

    // D := M * B
    cblas_ssymv(CblasColMajor, CblasUpper, km, 1.0f, M, n, B, 1, 0.0f, D, 1);

    // a := 1 / (Xnew' * Xnew - B' * D)
    auto a = 1.0f / (1.0f - cblas_sdot(km, B, 1, D, 1));

    // insert D by -1.0
    D[km] = -1.0f;

    // insert M by zeros
    for ( auto i = 0; i < k; i++ ) {
      M[km*n+i] = 0.0f;
    }

    // M += a * D * D'
    cblas_ssyr(CblasColMajor, CblasUpper, k, a, D, 1, M, n);

    // insert Theta by Xnew' * Y
    Theta[km] = cblas_sdot(n, Xnew, 1, Y, 1);

    // insert Beta by zero
    Beta[km] = 0.0f;

    // Beta += a * (D' * Theta) * D
    cblas_saxpy(k, a*cblas_sdot(k, D, 1, Theta, 1), D, 1, Beta, 1);

    // ==================================================================================================================== //
  } else {  // Backward step
    // Update size
    k--;

    // Update index
    I[idx] = false;

    // Find index
    auto j = Idx_ol[idx];

    // Copy index end to index j
    // a := M[end, end], b := Beta[end], D := M[end col]
    float a, b;
    if ( j != k ) {
      a = M[j*n+j];
      b = Beta[j];

      cblas_scopy(n, X+k*n, 1, X+j*n, 1);
      Beta[j] = Beta[k];
      Theta[j] = Theta[k];
      Idx_lo[j] = Idx_lo[k];
      Idx_ol[Idx_lo[j]] = j;

      cblas_scopy(j, M+j*n, 1, D, 1);
      cblas_scopy(k-j-1, M+j*n+n+j, n, D+j+1, 1);
      D[j] = M[k*n+j];

      cblas_scopy(j, M+k*n, 1, M+j*n, 1);
      cblas_scopy(k-j-1, M+k*n+j+1, 1, M+j*n+n+j, n);
      M[j*n+j] = M[k*n+k];
    } else {
      a = M[k*n+k];
      b = Beta[k];
      cblas_scopy(k, M+k*n, 1, D, 1);
    }

    // ======== Solve Y = X * Beta for Beta =============================================================================== //

    // M -= 1/a * D * D'
    cblas_ssyr(CblasColMajor, CblasUpper, k, -1.0f/a, D, 1, M, n);

    // Beta -= b/a * D
    cblas_saxpy(k, -b/a, D, 1, Beta, 1);

    // ==================================================================================================================== //
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Compute the residual
///
void Particle::ComputeResidual() {
  if ( !r_computed ) {
    // R = Y - X * Beta
    cblas_scopy(n, Y, 1, R, 1);
    cblas_sgemv(CblasColMajor, CblasNoTrans, n, k, -1.0f, X, n, Beta, 1, 1.0f, R, 1);
    r_computed = true;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Select the updating option
///
void Particle::SelectOption() {
  auto srand = static_cast<float>(rand_r(&iseed)) / RAND_MAX;
  if ( status ) {  // Forward step
    // Idx_best[0~kbest] := I0 exclude I
    // Idx_best[0~(p-k)] := complement of I
    kbest = 0;
    for ( auto i = 0, j = p-k; i < p; ++i ) {
      if ( !I[i] ) {
        if ( I0[i] ) {
          Idx_best[kbest] = i;
          kbest++;
        } else {
          j--;
          Idx_best[j] = i;
        }
      }
    }

    if ( kbest ) {
      option = (srand < parameter.prob_forward_best) + (srand < parameter.prob_forward_best+parameter.prob_forward_improve);
    } else {
      option = srand < parameter.prob_forward_improve / (parameter.prob_forward_improve+parameter.prob_forward_random);
    }
  } else {  // Backward step
    option = (srand < parameter.prob_backward_improve);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Select the index of the effect to update
///
void Particle::SelectIndex() {
  if ( status ) {  // Forward step
    switch( option ) {
      case Option::BEST: {  // Choose randomly from best model
        idx = Idx_best[rand_r(&iseed) % kbest];
        break;
      }
      case Option::IMPROVE: {  // Choose most improvement index
        // E := abs( E )
        vsAbs(p, E, E);

        // Remove selected entries
        for ( auto i = 0; i < k; ++i ) {
          E[Idx_lo[i]] = 0.0;
        }

        // Find maximum element
        idx = cblas_isamax(p, E, 1);
        break;
      }
      case Option::RANDOM: {  // Choose randomly
        idx = Idx_best[rand_r(&iseed) % (p-k)];
        break;
      }
    }
  } else {  // Backward step
    switch( option ) {
      case Option::IMPROVE: {  // Choose most improvement index
        ComputeResidual();
        auto e_temp = INFINITY;
        for ( auto i = 0; i < k; ++i ) {
          // B := R + Beta[i] * X[i col]
          cblas_scopy(n, R, 1, B, 1);
          cblas_saxpy(n, Beta[i], X+i*n, 1, B, 1);

          // stemp = norm(B)
          auto stemp = cblas_snrm2(n, B, 1);

          // Check if this value is minimal
          if ( e_temp > stemp ) {
            e_temp = stemp;
            idx = Idx_lo[i];
          }
        }
        break;
      }
      case Option::RANDOM: {  // Choose randomly
        idx = Idx_lo[rand_r(&iseed) % k];
        break;
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Compute the criterion value
///
void Particle::ComputeCriterion() {
  ComputeResidual();

  // e := norm(R)
  e = cblas_snrm2(n, R, 1);

  // Compute criterion
  switch(parameter.criterion) {
    case AIC: {    // phi := n*log(e^2/n) + 2k
      phi = n*logf(e*e/n) + 2.0f*k;
      break;
    }
    case BIC: {    // phi := n*log(e^2/n) + k*log(n)
      phi = n*logf(e*e/n) + k*logf(n);
      break;
    }
    case HQC: {    // phi := n*log(e^2/n) + 2k*log(log(n))
      phi = n*logf(e*e/n) + 2.0f*k*logf(logf(n));
      break;
    }
    case EBIC: {   // phi := n*log(e^2/n) + k*log(n) + 2gamma*log(binom(p, k))
      phi = n*logf(e*e/n) + k*logf(n) + 2.0f*parameter.ebic_gamma*lbinom(p, k);
      break;
    }
    case HDAIC: {  // phi := n*log(e^2/n) + 2k*log(p)
      phi = n*logf(e*e/n) + 2.0f*k*logf(p);
      break;
    }
    case HDBIC: {  // phi := n*log(e^2/n) + k*log(n)*log(p)
      phi = n*logf(e*e/n) + k*logf(n)*logf(p);
      break;
    }
    case HDHQC: {  // phi := n*log(e^2/n) + 2k*log(log(n))*log(p)
      phi = n*logf(e*e/n) + 2.0f*k*logf(logf(n))*logf(p);
      break;
    }
  }
}

}  // namespace genlin

}  // namespace pass
