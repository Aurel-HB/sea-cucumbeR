#include <TMB.hpp>

// Fonction pour détecter les NAs
template<class Type>
  bool isNA(Type x){
    return R_IsNA(asDouble(x));
  }

// Modèle spatial Gamma avec prédiction de biomasse
template<class Type>
  Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  // Données
  DATA_VECTOR( c_i );      // Densités observées (kg/m²)
  DATA_VECTOR( area_i );   // Surface (m²) pour chaque observation
  
  // Matrices SPDE
  DATA_SPARSE_MATRIX(M0);
  DATA_SPARSE_MATRIX(M1);
  DATA_SPARSE_MATRIX(M2);
  DATA_SPARSE_MATRIX(A_is); // INLA 'A' projection matrix for unique stations
  
  // Paramètres
  PARAMETER( beta0 );       // Intercept (log-scale)
  PARAMETER( ln_tau );     // Log de la précision spatiale
  PARAMETER( ln_kappa );   // Log du paramètre de lissage spatial
  PARAMETER( ln_alpha );    // Log du paramètre de forme (alpha) de la Gamma
  PARAMETER_VECTOR( omega_s );  // Effets aléatoires spatiaux
  
  // Quantités dérivées
  Type Range = sqrt(8) / exp(ln_kappa);
  Type SigmaE = 1 / sqrt(4 * M_PI * exp(2*ln_tau) * exp(2*ln_kappa));
  Type alpha = exp(ln_alpha);
  
  // Matrice de précision Q pour les effets aléatoires spatiaux
  Eigen::SparseMatrix<Type> Q = (exp(4*ln_kappa)*M0 + 
    Type(2.0)*exp(2*ln_kappa)*M1 + M2)*exp(2*ln_tau);
  
  // Fonction objectif
  vector<Type> jnll_comp(2);
  jnll_comp.setZero();
  
  // Vraisemblance des effets aléatoires spatiaux
  jnll_comp(1) += GMRF(Q)(omega_s);
  
  // Project using bilinear interpolation 
  vector <Type > omega_i( A_is.rows() ); 
  omega_i = A_is * omega_s;
  
  // Vraisemblance des données (Gamma)
  for( int i=0; i<c_i.size(); i++ ){
    if( !isNA(c_i(i)) ) {
      Type ln_mu_i = beta0 + omega_i(i);
      Type mu_i = exp(ln_mu_i);
      Type beta = mu_i / alpha;
      jnll_comp(0) -= dgamma(c_i(i), alpha, beta, true);
    }
  }
  
  // Rapport des résultats
  Type jnll = jnll_comp.sum();
  REPORT(jnll_comp);
  REPORT(jnll);
  REPORT(Range);
  REPORT(SigmaE);
  REPORT(beta0);
  REPORT(omega_s);
  REPORT(alpha);
  REPORT(omega_i);
  REPORT(ln_kappa);
  REPORT(ln_tau);
  ADREPORT(omega_s)
  
  return jnll;
}
