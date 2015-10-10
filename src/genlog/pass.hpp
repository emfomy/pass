////////////////////////////////////////////////////////////////////////////////
// Particle Swarm Stepwise (PaSS) Algorithm                                   //
//                                                                            //
// pass.hpp                                                                   //
// The header of the PaSS algorithm for general logistic regression           //
//                                                                            //
// Author: emfo<emfomy@gmail.com>                                             //
////////////////////////////////////////////////////////////////////////////////

#ifndef PASS_GENLOG_PASS_HPP_

#define PASS_GENLOG_PASS_HPP_

////////////////////////////////////////////////////////////////////////////////
// The namespace pass                                                         //
////////////////////////////////////////////////////////////////////////////////
namespace pass {

// Global variables
extern int n;                       // scalar, the number of statistical units
extern int p;                       // scalar, the number of total effects
extern float* X0;                   // matrix, n by p, the regressors
extern float* Y0;                   // vector, n by 1, the regressand
extern bool* I0;                    // vector, 1 by p, the chosen indices
extern float phi0;                  // scalar, the criterion value
extern struct Parameter parameter;  // the PaSS parameters

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
void GenLog();

////////////////////////////////////////////////////////////////////////////////
// The criterions used in the PaSS algorithm                                  //
////////////////////////////////////////////////////////////////////////////////
enum Criterion {
  AIC,   // Akaike information criterion
  BIC,   // Bayesian information criterion
  EBIC,  // Extended Bayesian information criterion
  HDBIC, // High-dimensional Bayesian information criterion
  HQC,   // Hannan-Quinn information criterion
  HDHQC  // High-dimensional Hannan-Quinn information criterion
};

////////////////////////////////////////////////////////////////////////////////
// Change criterion to string                                                 //
//                                                                            //
// Parameters:                                                                //
// criterion:  the criterion                                                  //
////////////////////////////////////////////////////////////////////////////////
const char* Criterion2String( const Criterion criterion );

////////////////////////////////////////////////////////////////////////////////
// The parameters of the PaSS algorithm                                       //
////////////////////////////////////////////////////////////////////////////////
struct Parameter {
  unsigned int num_iteration;        // the number of iterations
  unsigned int num_particle_thread;  // the number of particles per thread
  float prob_forward_global;         // the probability of forward step: global
  float prob_forward_local;          // the probability of forward step: local
  float prob_forward_random;         // the probability of forward step: random
  float prob_backward_local;         // the probability of backward step: local
  float prob_backward_random;        // the probability of backward step: random
  Criterion criterion;               // the criterion
  float ebic_gamma;                  // the penalty parameter for EBIC
  bool is_normalized;                // the data is normalized of not

  // Constructor
  Parameter() {
    num_iteration = 1024;
    num_particle_thread = 16;
    prob_forward_global = 0.1;
    prob_forward_local = 0.8;
    prob_forward_random = 0.1;
    prob_backward_local = 0.9;
    prob_backward_random = 0.1;
    criterion = HDBIC;
    ebic_gamma = 1.0;
    is_normalized = false;
  }
};

////////////////////////////////////////////////////////////////////////////////
// The structure of a particle                                                //
////////////////////////////////////////////////////////////////////////////////
struct Particle {
  float *X;      // matrix, n by k, the regressors
  float *Y;      // vector, n by 1, the regressand
  float *Beta;   // vector, k by 1, the effects
  float *Theta;  // vector, n by 1, X'*Y, the logit of P
  float *Eta;    // vector, n by 1, exp(Theta)
  float *P;      // vector, n by 1, Eta./(1+Eta), the probability of Y=1
  float *W;      // vector, n by 1, P.*(1-P)
  float *M;      // matrix, k by k, X'*diag(W)*X, lower packed storage
  float *Ones;   // vector, n by 1, ones
  float *STemp;  // vector, n by 1, temporary vector
  float llv;     // scalar, the log-likelihood value
  float phi;     // scalar, the criterion value
  float phi_old; // scalar, the criterion value, past iteration

  int *Idx_lf;   // vector, 1 by k, map local effects to full effects
  int *Idx_fl;   // vector, 1 by p, map full effects to local effects
  int *Idx_temp; // vector, 1 by p, workspace
  bool *I;       // vector, 1 by p, the chosen indices
  int k;         // scalar, the number of chosen effects
  int l;         // scalar, the number of chosen indices

  bool status;   // scalar, the status (forward/backward)

  unsigned int iseed;  // scalar, the random seed;

  // Constructor
  Particle();

  // Destructor
  ~Particle();

  // Initialize model
  void InitializeModel();
  void InitializeModel( const int idx );

  // Update model
  void UpdateModel( const int idx );

  // Select the index to add or remove
  void SelectIndex( int& idx );

  // Compute Beta
  void ComputeBeta();

  // Compute the criterion value
  void ComputeCriterion();
};

}

#endif
