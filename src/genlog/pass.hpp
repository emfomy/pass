////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @file   genlog/pass.hpp
/// @brief  The header of PaSS for general logistic regression
///
/// @author Mu Yang <emfomy@gmail.com>
///

#ifndef PASS_GENLOG_PASS_HPP_

#define PASS_GENLOG_PASS_HPP_

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// The namespace of PaSS
///
namespace pass {

extern int n;                       ///< scalar, the number of statistical units
extern int p;                       ///< scalar, the number of total effects
extern float* X0;                   ///< matrix, n by p, the regressors
extern float* Y0;                   ///< vector, n by 1, the regressand
extern bool* I0;                    ///< vector, 1 by p, the chosen indices
extern float phi0;                  ///< scalar, the criterion value
extern struct Parameter parameter;  ///< the PaSS parameters

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  The PaSS algorithm for logistic regression
//
void GenLog();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// The criterions used in the PaSS algorithm
///
enum Criterion {
  AIC,   ///< Akaike information criterion
  BIC,   ///< Bayesian information criterion
  EBIC,  ///< Extended Bayesian information criterion
  HDBIC, ///< High-dimensional Bayesian information criterion
  HQC,   ///< Hannan-Quinn information criterion
  HDHQC  ///< High-dimensional Hannan-Quinn information criterion
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Change criterion to string
///
/// @param[in]  criterion  the criterion
///
/// @return                the name of criterion
///
static inline const char* Criterion2String( const Criterion criterion ) {
  switch(criterion) {
    case AIC:   return "AIC";
    case BIC:   return "BIC";
    case EBIC:  return "EBIC";
    case HDBIC: return "HDBIC";
    case HQC:   return "HQC";
    case HDHQC: return "HDHQC";
    default:    return "!";
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// The parameters of the PaSS algorithm
///
struct Parameter {
  unsigned int num_iteration;        ///< the number of iterations
  unsigned int num_particle_thread;  ///< the number of particles per thread
  float prob_forward_global;         ///< the probability of forward step: global
  float prob_forward_local;          ///< the probability of forward step: local
  float prob_forward_random;         ///< the probability of forward step: random
  float prob_backward_local;         ///< the probability of backward step: local
  float prob_backward_random;        ///< the probability of backward step: random
  Criterion criterion;               ///< the criterion
  float ebic_gamma;                  ///< the penalty parameter for EBIC
  bool is_normalized;                ///< the data is normalized of not

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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// The structure of particle
///
struct Particle {
  float *X;            ///< matrix, n by (k+1), the regressors
  float *Y;            ///< vector, n by 1, the regressand
  float *Beta;         ///< vector, (k+1) by 1, the effects
  float *Theta;        ///< vector, n by 1, X*Beta, the logit of P
  float *Eta;          ///< vector, n by 1, exp(Theta)
  float *P;            ///< vector, n by 1, Eta./(1+Eta), the probability of Y=1
  float *W;            ///< vector, n by 1, P.*(1-P)
  float *M;            ///< matrix, (k+1) by (k+1), X'*diag(W)*X, lower packed storage
  float *STemp;        ///< vector, n by 1, temporary vector
  float llv;           ///< scalar, the log-likelihood value
  float phi;           ///< scalar, the criterion value
  float phi_old;       ///< scalar, the criterion value, past iteration

  int *Idx_lo;         ///< vector, 1 by (k+1), map local effects to original effects
  int *Idx_ol;         ///< vector, 1 by p, map original effects to local effects
  int *Idx_temp;       ///< vector, 1 by p, workspace
  bool *I;             ///< vector, 1 by p, the chosen indices
  int k;               ///< scalar, the number of chosen effects

  bool status;         ///< scalar, the status (forward/backward)

  unsigned int iseed;  ///< scalar, the random seed

  // Constructor
  Particle();

  // Destructor
  ~Particle();

  // Initialize model
  void InitializeModel();
  void InitializeModel( const int idx );

  // Update model
  void UpdateModel( const int idx );

  // Compute Beta
  void ComputeBeta();

  // Select the index to add or remove
  void SelectIndex( int& idx );

  // Compute the criterion value
  void ComputeCriterion();
};

}

#endif