////////////////////////////////////////////////////////////////////////////////
// Particle Swarm Stepwise (PaSS) Algorithm                                   //
//                                                                            //
// pass.cpp                                                                   //
// The PaSS algorithm for general linear regression                           //
//                                                                            //
// Author: Mu Yang <emfomy@gmail.com>                                         //
//                                                                            //
// ========================================================================== //
//                                                                            //
// Notation:                                                                  //
// X    : the regressors                                                      //
// Y    : the regressand                                                      //
// Beta : the effects                                                         //
// R    : the residual                                                        //
//                                                                            //
// ========================================================================== //
//                                                                            //
// Linear model:                                                              //
// Y = X * Beta + error                                                       //
// R = Y - X * Beta                                                           //
//                                                                            //
// ========================================================================== //
//                                                                            //
// Select index in forward step:                                              //
// idx = argmax_{i not in I} abs( X[i col]' * R  )                            //
//                                                                            //
// Select index in backward step:                                             //
// idx = argmin_{i in I} norm( R_{I exclude i} )                              //
//     = argmin_{i in I} norm( R + Beta[i] * X[i col] )                       //
//                                                                            //
// ========================================================================== //
//                                                                            //
// References:                                                                //
// Chen, R.-B., Huang, C.-C., & Wang, W. (2013).                              //
//   Particle Swarm Stepwise (PaSS) Algorithm for Variable Selection.         //
////////////////////////////////////////////////////////////////////////////////

#include "pass.hpp"
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <mkl.h>
#include <omp.h>

// The log-binomial function
static inline float lbinom( const int n, const int k ) {
  int i;
  return (lgammaf_r(n+1, &i) - lgammaf_r(n-k+1, &i) - lgammaf_r(k+1, &i));
}

////////////////////////////////////////////////////////////////////////////////
// The namespace pass                                                         //
////////////////////////////////////////////////////////////////////////////////
namespace pass {
int n;                // scalar, the number of statistical units
int p;                // scalar, the number of total effects
float* X0;            // matrix, n by p, the regressors
float* Y0;            // vector, n by 1, the regressand
bool* I0;             // vector, 1 by p, the chosen indices
float phi0;           // scalar, the criterion value
Parameter parameter;  // the PaSS parameters

////////////////////////////////////////////////////////////////////////////////
// The PaSS algorithm for Linear Regression                                   //
//                                                                            //
// Input Global Parameters:                                                   //
// n:         scalar, the number of statistical units                         //
// p:         scalar, the number of total effects                             //
// X0:        matrix, n by p, the regressors                                  //
// Y0:        vector, n by 1, the regressand                                  //
// parameter: the PaSS parameters                                             //
//                                                                            //
// Output Global Variables:                                                   //
// I0:        vector, 1 by p, the chosen indices                              //
// phi0:      scalar, the criterion value                                     //
//                                                                            //
// Note:                                                                      //
// Please call srand before using this routine.                               //
////////////////////////////////////////////////////////////////////////////////
void GenLin() {
  // Check parameters
  auto num_thread = omp_get_max_threads();
  auto num_particle = num_thread * parameter.num_particle_thread;

  ////////////////////////////////////////////////////////////////////////////
  // Centralize and normalize the original data                             //
  ////////////////////////////////////////////////////////////////////////////

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

  ////////////////////////////////////////////////////////////////////////////
  // Run PaSS                                                               //
  ////////////////////////////////////////////////////////////////////////////

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
        if ( particle[j].k >= n || particle[j].k >= p-4 ) {
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

  ////////////////////////////////////////////////////////////////////////////

  // Delete memory
  delete[] particle;
}

////////////////////////////////////////////////////////////////////////////////
// The constructor of Particle                                                //
////////////////////////////////////////////////////////////////////////////////
Particle::Particle() {
  X        = new float[n*n];
  Y        = new float[n];
  Beta     = new float[n];
  Theta    = new float[n];
  M        = new float[n*n];
  R        = new float[n];
  B        = new float[n];
  D        = new float[n];

  Idx_lo   = new int[n];
  Idx_ol   = new int[p];
  Idx_temp = new int[p];
  I        = new bool[p];

  iseed    = rand();
}

////////////////////////////////////////////////////////////////////////////////
// The destructor of Particle                                                 //
////////////////////////////////////////////////////////////////////////////////
Particle::~Particle() {
  delete[] X;
  delete[] Y;
  delete[] Beta;
  delete[] Theta;
  delete[] M;
  delete[] B;
  delete[] D;
  delete[] R;

  delete[] Idx_lo;
  delete[] Idx_ol;
  delete[] Idx_temp;
  delete[] I;
}

////////////////////////////////////////////////////////////////////////////////
// Initialize the model of Particle randomly                                  //
////////////////////////////////////////////////////////////////////////////////
void Particle::InitializeModel() {
  InitializeModel(rand_r(&iseed) % p);
}

////////////////////////////////////////////////////////////////////////////////
// Initialize the model of Particle                                           //
//                                                                            //
// Parameters:                                                                //
// idx: the index of the effect                                               //
////////////////////////////////////////////////////////////////////////////////
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
}

////////////////////////////////////////////////////////////////////////////////
// Update the model of Particle                                               //
//                                                                            //
// Parameters:                                                                //
// idx: the index of the effect                                               //
////////////////////////////////////////////////////////////////////////////////
void Particle::UpdateModel( const int idx ) {
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

    ////////////////////////////////////////////////////////////////////////
    // Solve Y = X * Beta for Beta                                        //
    ////////////////////////////////////////////////////////////////////////

    // B := X' * Xnew
    cblas_sgemv(CblasColMajor, CblasTrans,
                n, km, 1, X, n, Xnew, 1, 0.0f, B, 1);

    // D := M * B
    cblas_ssymv(CblasColMajor, CblasUpper,
                km, 1.0f, M, n, B, 1, 0.0f, D, 1);

    // a := 1 / (Xnew' * Xnew - B' * D)
    auto a = 1.0f / (1.0f - cblas_sdot(km, B, 1, D, 1));

    // insert D by -1.0
    D[km] = -1.0f;

    // insert M by zeros
    for ( auto i = 0; i < k; i++ ) {
      M[km*n+i] = 0.0f;
    }

    // M += a * D * D'
    cblas_ssyr(CblasColMajor, CblasUpper,
               k, a, D, 1, M, n);

    // insert Theta by Xnew' * Y
    Theta[km] = cblas_sdot(n, Xnew, 1, Y, 1);

    // insert Beta by zero
    Beta[km] = 0.0f;

    // Beta += a * (D' * Theta) * D
    cblas_saxpy(k, a*cblas_sdot(k, D, 1, Theta, 1), D, 1, Beta, 1);

    ////////////////////////////////////////////////////////////////////////
  } else {  // backward step
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

    ////////////////////////////////////////////////////////////////////////
    // Solve Y = X * Beta for Beta                                        //
    ////////////////////////////////////////////////////////////////////////

    // M -= 1/a * D * D'
    cblas_ssyr(CblasColMajor, CblasUpper,
               k, -1.0f/a, D, 1, M, n);

    // Beta -= b/a * D
    cblas_saxpy(k, -b/a, D, 1, Beta, 1);

    ////////////////////////////////////////////////////////////////////////
  }
  
  // R = Y - X * Beta
  cblas_scopy(n, Y, 1, R, 1);
  cblas_sgemv(CblasColMajor, CblasNoTrans,
              n, k, -1.0f, X, n, Beta, 1, 1.0f, R, 1);
}

////////////////////////////////////////////////////////////////////////////////
// Select index of the effect to update                                       //
//                                                                            //
// Output Parameters:                                                         //
// idx: the index of the effect                                               //
////////////////////////////////////////////////////////////////////////////////
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
      choose = (srand < parameter.prob_forward_global)
             + (srand < parameter.prob_forward_global+
                        parameter.prob_forward_local);
    } else {
      choose = srand < parameter.prob_forward_local
                     / (parameter.prob_forward_local+
                        parameter.prob_forward_random);
    }

    switch( choose ) {
      case 2: {  // Global best
        idx = Idx_temp[rand_r(&iseed) % itemp];
        break;
      }
      case 1: {  // Local best
        auto e_temp = -INFINITY;
        for ( auto i = 0; i < p; ++i ) {
          if ( !I[i] ) {
            // stemp := abs( X0[i col]' * R )
            auto stemp = fabs(cblas_sdot(n, X0+i*n, 1, R, 1));

            // Check if this value is maximum
            if ( e_temp < stemp ) {
              e_temp = stemp;
              idx = i;
            }
          }
        }
        break;
      }
      case 0: {  // Random
        idx = Idx_temp[rand_r(&iseed) % (p-k)];
        break;
      }
    }
  } else {  // backward step
    if ( srand < parameter.prob_backward_local ) {  // Local best
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
    } else {  // Random
      idx = Idx_lo[rand_r(&iseed) % k];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Compute the criterion value                                                //
////////////////////////////////////////////////////////////////////////////////
void Particle::ComputeCriterion() {
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
    case EBIC: {   // phi := n*log(e^2/n) + k*log(n) + 2gamma*log(p choose k)
      phi = n*logf(e*e/n) + k*logf(n) + 2.0f*parameter.ebic_gamma*lbinom(p, k);
      break;
    }
    case HDBIC: {  // phi := n*log(e^2/n) + k*log(n)*log(p)
      phi = n*logf(e*e/n) + k*logf(n)*logf(p);
      break;
    }
    case HQC: {    // phi := n*log(e^2/n) + 2k*log(log(n))
      phi = n*logf(e*e/n) + 2.0f*k*logf(logf(n));
      break;
    }
    case HDHQC: {  // phi := n*log(e^2/n) + 2k*log(log(n))*log(p)
      phi = n*logf(e*e/n) + 2.0f*k*logf(logf(n))*logf(p);
      break;
    }
  }
}

}
