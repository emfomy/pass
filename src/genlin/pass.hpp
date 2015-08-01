////////////////////////////////////////////////////////////////////////////////
// Particle Swarm Stepwise (PaSS) Algorithm                                   //
//                                                                            //
// genlin.hpp                                                                 //
// The header of the PaSS algorithm for General Linear Model                  //
//                                                                            //
// Author: emfo<emfomy@gmail.com>                                             //
////////////////////////////////////////////////////////////////////////////////

#ifndef PASS_GENLIN_PASS_HPP_

#define PASS_GENLIN_PASS_HPP_

////////////////////////////////////////////////////////////////////////////////
// The namespace pass                                                         //
////////////////////////////////////////////////////////////////////////////////
namespace pass {

// Global variables
extern int n;                      // scalar, the number of statistical units
extern int p;                      // scalar, the number of total effects
extern float* X0;                  // matrix, n by p, the regressors
extern float* Y0;                  // vector, n by 1, the regressand
extern bool* I0;                   // vector, 1 by p, the chosen indices
extern struct Parameter parameter; // the parameters

////////////////////////////////////////////////////////////////////////////////
// The PaSS algorithm for Linear Regression                                   //
//                                                                            //
// Input Parameters:                                                          //
// n:         scalar, the number of statistical units                         //
// p:         scalar, the number of total effects                             //
// X:         matrix, n by p, the regressors                                  //
// Y:         vector, n by 1, the regressand                                  //
// parameter: the parameters                                                  //
//                                                                            //
// Output Global Variables:                                                   //
// I:         vector, 1 by p, the chosen indices                              //
////////////////////////////////////////////////////////////////////////////////
void GenLin();

////////////////////////////////////////////////////////////////////////////////
// The criterions used in the PaSS algorithm                                  //
////////////////////////////////////////////////////////////////////////////////
enum Criterion {
  AIC,   // Akaike information criterion
  BIC,   // Bayesian information criterion
  EBIC,  // Extended Bayesian information criterion
  HDBIC, // High-dimensional Bayesian information criterion
  HDHQ   // High-dimensional Hannan-Quinn information criterion
};

////////////////////////////////////////////////////////////////////////////////
// Change criterion to string                                                 //
//                                                                            //
// Parameters:                                                                //
// criterion:  the criterion                                                  //
////////////////////////////////////////////////////////////////////////////////
static const char* Criterion2String( const Criterion criterion ) {
  switch(criterion) {
    case AIC: {
      return "AIC";
    }
    case BIC: {
      return "BIC";
    }
    case EBIC: {
      return "EBIC";
    }
    case HDBIC: {
      return "HDBIC";
    }
    case HDHQ: {
      return "HDHQ";
    }
    default: {
      return "";
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// The parameters of the PaSS algorithm                                       //
////////////////////////////////////////////////////////////////////////////////
struct Parameter {
  int num_iteration;           // the number of iterations
  float prob_forward_global;   // the probability of forward step: global
  float prob_forward_local;    // the probability of forward step: local
  float prob_forward_random;   // the probability of forward step: random
  float prob_backward_local;   // the probability of backward step: local
  float prob_backward_random;  // the probability of backward step: random
  Criterion criterion;         // the criterion
  float ebic_gamma;            // the penalty parameter for EBIC
  bool is_normalized;          // the data is normalized of not

  // Constructor
  Parameter() {
    num_iteration = 1024;
    prob_forward_global = 0.1;
    prob_forward_local = 0.8;
    prob_forward_random = 0.1;
    prob_backward_local = 0.9;
    prob_backward_random = 0.1;
    criterion = EBIC;
    ebic_gamma = 1.0;
    is_normalized = false;
  }
};

////////////////////////////////////////////////////////////////////////////////
// The structure of a particle                                                //
////////////////////////////////////////////////////////////////////////////////
struct Particle {
  float *X;        // matrix, n by k, the regressors
  float *Y;        // vector, n by 1, the regressand
  float *Beta;     // vector, k by 1, the effects
  float *Theta;    // vector, k by 1, X'*Y
  float *M;        // matrix, k by k, inv( X'*X ), upper general storage
  float *R;        // vector, n by 1, the residual
  float *B;        // vector, n by 1, temporary vecter
  float *D;        // vector, n by 1, temporary vecter
  float e;         // scalar, the norm of R
  float phi;       // scalar, the value given by criterion
  float phi_old;   // scalar, the value given by criterion, past iteration
  float phi_best;  // scalar, the value given by criterion, global best

  int *Idx_lf;     // vector, 1 by k, map local effects to full effects
  int *Idx_fl;     // vector, 1 by p, map full effects to local effects
  bool *I;         // vector, 1 by p, the chosen indices
  bool *I_best;    // vector, 1 by p, the chosen indices, global best
  int k;           // scalar, the number of chosen effects
  int l;           // scalar, the number of chosen indices

  bool status;     // scalar, the status (forward/backward)

  int *Idx_temp;   // vector, p by 1

  // Constructor
  Particle();

  // Destructor
  ~Particle();

  // Initialize model
  void InitializeModel( const int idx );

  // Update model
  void UpdateModel( const int idx );

  // Select the index to add or remove
  void SelectIndex( int& idx );

  // Compute the value given by criterion
  void ComputeCriterion();
};

}

#endif
