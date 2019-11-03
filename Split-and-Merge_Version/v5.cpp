// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppNumerical.h>
#include <omp.h>

// [[Rcpp::plugins(openmp)]]
using namespace Numer;
using namespace Rcpp;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Useful functions:

// [[Rcpp::export]]
double dinvgamma_cpp( double x, double a, double b){
  return  R::dgamma(1/x,a,1/b,0)/(x*x);
}

// [[Rcpp::export]]
double log_g_over_g(double a, double b){
  return lgamma(a+b)-lgamma(a);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for numerical integration, used for the computation of the Marginal likelihood, needed for Theta
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Marginal_p0_integrand: public Func  // int_sigma int_mu LIK*N(0,1)*dgamma = int_sigma dgamma * marginal N-N model dsigma
{
private:
  double med; 
  double sq;
  double a;   
  double b; 
  int n;
public:
  Marginal_p0_integrand(double med_,  double sq_, double a_, double b_, int n_) : med(med_), sq(sq_), a(a_), b(b_), n(n_) {}
  double operator()(const double& x) const
  {
    return  pow(x,-(n-1)/2)/pow((2*arma::datum::pi),n/2)  *  pow(x+n,-.5) *
      exp(n * n * med * med * .5 / (x*(x+n)) - sq/(2*x)) * dinvgamma_cpp(x,a,b);
  }
};

// [[Rcpp::export]]
double Marginal_p0( double med, double sq, double a, double b, int n, double low, double upp)
{
  Marginal_p0_integrand f(med,sq,a,b,n);
  double err_est; int err_code;
  const double res = integrate(f, low, upp, err_est, err_code);
  return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double log_Marginal_P0(arma::colvec x, double a, double b, double low=0, double upp=200) {
  int    n   = x.n_rows;
  double sq  = accu(x%x), med = mean(x);
  return log(Marginal_p0(med, sq, a, b, n,low,upp));
}

// [[Rcpp::export]]
arma::colvec log_Marginal_P0_grid(arma::colvec grid, double a, double b, double low=0, double upp=200) {
  int    n_g   = grid.n_rows;
  arma::colvec e_grid(n_g);
  for(int h=0; h<n_g; h++ ){
    double x   = grid[h];
    double sq  = (x*x), med = (x);
    e_grid(h) =  log(Marginal_p0(med, sq, a, b, 1,low,upp));
  }
  return e_grid;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Marginal likelihood for P1 decoupled in two functions:
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]] // P1 is a mixture: here we compute the marginal liklhood wrt 1 component
double log_Marginal_P1_single(arma::colvec x, double m1, double V1, double a1, double b1) {
  
  int n         = x.n_rows;
  double Vn     = 1 / ( n + 1/V1 ),
    ssq    = sum(x%x),
    med    = mean(x),
    an     = a1 + n/2,
    mn     = (m1 / V1 + n * med) * Vn,
    bn     = b1 + .5 * ( m1 *  m1 * 1/V1 + ssq - mn * mn / Vn), val=0; 
  val = .5*(log(Vn)-log(V1)) + a1 * log(2*b1) - an*log(2*bn) + lgamma(an) - lgamma(a1) - n/2*log(arma::datum::pi);
  return (val); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}


// [[Rcpp::export]] 
double log_Marginal_P1(arma::colvec x, double m1, double V1, double a1, double b1) {
  arma::colvec val12(2);
  val12(0) = log_Marginal_P1_single( x,  m1, V1, a1, b1) ;
  val12(1) = log_Marginal_P1_single( x, -m1, V1, a1, b1) ;
  double maxlp = max(val12);
  return maxlp + log( exp(val12(0)-maxlp)+ exp(val12(1)-maxlp) ) +log(.5);  //mancava ultimo termine, pirla! è lo 0.5 della mistura
}

// [[Rcpp::export]] 
double log_Marginal_P1old(arma::colvec x, double m0, double V0, double a0, double b0) {
  
  double val = .5*exp( log_Marginal_P1_single( x, m0, V0, a0, b0) ) +
               .5*exp( log_Marginal_P1_single( x, -m0, V0, a0, b0) );
  
  return log(val); //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
}

// [[Rcpp::export]]
arma::colvec log_Marginal_P1_grid(arma::colvec grid, double m1, double V1, double a1, double b1) {
  int    n_g   = grid.n_rows;
  arma::colvec e_grid(n_g);
  
  for(int h=0; h<n_g; h++ ){
    arma::colvec A(1); A.fill(grid(h));
    e_grid(h) =  log_Marginal_P1( A, m1, V1, a1, b1);
  }
  return e_grid;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Functions for generating random variates from NIG
// [[Rcpp::export]]
arma::colvec NIGmix_generator_sample(arma::colvec y, double m1, double V1, double a1, double b1){
  
  int    n1  = y.n_rows;
  double k1  = 1/V1,
    kn  = k1+n1, 
    an  = a1 + n1/2,
    SQ1 = sum(y%y), 
    M1  = mean(y),
    bn, mn, mNeg = -m1;
  
  ///// correction
  double P1 = exp(log_Marginal_P1_single(y,  m1, V1,a1,b1));
  double P2 = exp(log_Marginal_P1_single(y, -m1, V1,a1,b1));  //IO CREDO NON DEBBA CAMBIAR NULLA
  double P1norm = P1/(P1+P2);
  
  if( R::runif(0,1) < P1norm ){
    mn = (k1*m1+M1*n1)/(kn);
    bn = .5*( k1*m1*m1 + SQ1 - mn*mn*kn  ) + b1;
  }else{
    mn = (k1*mNeg+M1*n1)/(kn);
    bn = .5*( k1*mNeg*mNeg + SQ1 - mn*mn*kn  ) + b1;
  }
  arma::colvec TH(2);
  TH[1] = 1/ R::rgamma( an, 1/(bn) );
  TH[0] = R::rnorm(mn, sqrt( TH[1]/kn  ));
  return TH;
}




// [[Rcpp::export]]
arma::colvec NindepIG_generator_sample(arma::colvec y, double m0, double s0, double previous_variance, double a, double b){
  int n1=y.n_elem;
  double sn = s0 * previous_variance / ( n1*s0 + previous_variance);
  double mn = (s0*accu(y)+m0*previous_variance)/(n1*s0+previous_variance);
  arma::colvec TH(2);
  TH[0] = R::rnorm(mn, sqrt( sn ));
  TH[1] = 1/R::rgamma(.5*n1+a, 1/( b+.5*accu((y-TH[0])%(y-TH[0]))) );
  return TH;
}

// [[Rcpp::export]]
arma::colvec NindepIG_generator_sample_GIBBS(arma::colvec y, double m0, double s0, double a, double b, int SUP=20){
  //microGibbs inside
  double mu=0, sigma=0;
  ////////////////////////////////////////////????
  int n1 = y.n_elem;
  ////////////////////////////////////////////////
  for(int k=0;k<SUP;k++){
    if(k==0){sigma = 1/R::rgamma(a,1/b);}
    double sn = s0 * sigma / ( n1*s0 + sigma);
    double mn = (s0*accu(y) + m0*sigma)/(n1*s0+sigma);
    mu        = R::rnorm(mn, sqrt( sn ));
    sigma     = 1/R::rgamma( n1/2+a , 1/( b+.5*accu( (y-mu)%(y-mu)) ) );
  }
  arma::colvec TH(2);
  TH[0] = mu;
  TH[1] = sigma;
  
  return TH;
}


// this gets unique rows from a matrix 
template <typename T>
inline bool rows_equal(const T& lhs, const T& rhs, double tol = 0.00000001) {
  return arma::approx_equal(lhs, rhs, "absdiff", tol);
}

arma::mat unique_rows(const arma::mat& x) {
  unsigned int count = 1, i = 1, j = 1, nr = x.n_rows, nc = x.n_cols;
  arma::mat result(nr, nc);
  result.row(0) = x.row(0);
  
  for ( ; i < nr; i++) {
    bool matched = false;
    if (rows_equal(x.row(i), result.row(0))) continue;
    
    for (j = i + 1; j < nr; j++) {
      if (rows_equal(x.row(i), x.row(j))) {
        matched = true;
        break;
      }
    }
    if (!matched) result.row(count++) = x.row(i);
  }
  
  return result.rows(0, count - 1);
}

// Table for matrices: needed for updating mu and sigma (even if a simpler unique(mu) should suffice )
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List TABLE_matrix(int i, arma::colvec zi, arma::mat mu_sig){ 
  mu_sig.shed_row(i);
  int thisZZ = zi(i);
  zi.shed_row(i);
  
  arma::colvec mu = mu_sig.col(0), sig = mu_sig.col(1);
  arma::uvec  ind = find(zi == thisZZ);
  if( ind.n_rows==0){
    //Rcout << "solo!";
    int Kzi=0;
    return List::create(_["Kzi"] = Kzi);
  }else{
    arma::mat submat = mu_sig.rows(ind), submat2 = unique_rows(submat);
    arma::colvec submu = submat.col(0), 
      subsig = submat2.col(1), //attento, submat2!!
      u_submu = unique(submu); // sorted, because table returns a sorted vector
    arma::mat DDD=join_rows(sort(submat2.col(0)),
                            subsig(arma::sort_index(submat2.col(0))));
    
    NumericVector submuNV = wrap(submu);
    arma::colvec TI       = as<arma::colvec>(table(submuNV));
    
    
    int Kzi = u_submu.n_rows;
    
    return List::create(
      _["Pn"] = DDD,
      _["TI"] = TI,
      _["Kzi"] = Kzi);
  }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////// 
// FULL CONDITIONALS
///////////////////////////////////////

// [[Rcpp::export]]
List pk_ratio(arma::colvec yyy, 
              arma::colvec TH, 
              arma::colvec zi,  
              double a_rho,               double b_rho, 
              double sigma0,              double sigma1, 
              double theta0,              double theta1, 
              double m0,               double s0,
              double a0,               double b0,
              double a1,               double b1, 
              double m1,               double V1,
              int n, int K, 
              arma::colvec zz,
              arma::colvec uTH // we can recover the clustering from the mean
){
  
  R_CheckUserInterrupt();
  
  int n_star = 0;
  
  for(int k = 0; k<K; k++){
    
    arma::uvec indTH = find( TH == uTH[k] ), // super cautelativo, non può esserci situazione opposta per Lemma 1
      indNOTH        = find( TH != uTH[k] ),
      indZ1          = find( zi == 1  ),
      ind12          = arma::intersect( indNOTH, indZ1 ); // obs in P^*_1, not in k-th group
    
    arma::colvec Yk    = yyy.elem( indTH);//ind11 );
    arma::colvec zz_temp  = zz;
    zz_temp.shed_row(k);
    
    int K_k1 = sum( zz_temp ), // how many clusters in P1, excluding k
      nk = Yk.n_rows,
      n_k1 = ind12.n_elem; // how many obs : zi==1, not in cluster k
    
    double logNUM = log_Marginal_P0(Yk, a0,  b0)    + 
      log_g_over_g(b_rho + n - n_k1 - nk, nk) + 
      log_g_over_g(1 - sigma0, nk-1) + 
      log(theta0 + (K-K_k1-1)*sigma0) + 
      log_g_over_g(theta1 + n_k1, nk),
      
      logDEN = log_Marginal_P1(Yk, m1, V1, a1, b1) + 
        log_g_over_g(a_rho + n_k1 , nk) + 
        log_g_over_g(1 - sigma1, nk-1) + 
        log(theta1 + (K_k1)*sigma1) + 
        log_g_over_g(theta0 + n - n_k1 - nk, nk),
        prob_1 = 1.0/(1.0 + exp(logNUM-logDEN));
    
    if( runif(1)[0] < prob_1 ){ zz[k]=1; }else{ zz[k]=0; }    // here we update zz
    arma::colvec X(nk); X.fill(zz[k]); zi.elem(indTH) = X;    // here we update zi
    
    n_star += nk*zz[k];
  }
  double rho = rbeta(1, a_rho + n_star, b_rho + n - n_star)[0]; 
  return List::create(
    _["zi"]    = zi,
    _["rho"]   = rho );
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// idea qui: fare update di pk con indici di cluster e non coi theta



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat Get_THETA_updated(
    arma::mat data,         // should contain: (yy, mu, sigma, zi) 
    double sig0,            // 2PDP discount parameter
    double theta0,          // 2PDP concentration parameter
    double sig1,            // 2PDP discount parameter
    double theta1,          // 2PDP concentration parameter
    double a0, double b0,
    double m0,  double s0, 
    double m1, double V1, 
    double a1, double b1)
{
  arma::colvec yy = data.col(0),
    zi = data.col(3);
  arma::mat musig = data.cols(1,2);
  int N = yy.n_rows;
  arma::rowvec NEW_THETA(2);
  
  for(int i=0; i<N; i++){
    double this_yy  = yy[i], th_process,  sig_process, log_Marginal_lik_value;  
    int    this_zi  = zi[i];
    ////////////////////////////////////////////////////////////////
    List LL = TABLE_matrix(i, zi, musig);
    int  Kzi     = LL["Kzi"];
    ////////////////////////////////////////////////////////////////
    arma::colvec this_yy_colvec(1); this_yy_colvec.fill(this_yy); // so we can use the marginal likelihood function for just one obs
    
    if(this_zi == 0){
      sig_process = sig0;
      th_process  = theta0;
      log_Marginal_lik_value = log_Marginal_P0(this_yy_colvec, a0, b0);
    }else{ // unique values of param and n_i when i is excluded
      sig_process = sig1;
      th_process  = theta1;
      log_Marginal_lik_value = log_Marginal_P1(this_yy_colvec, m1, V1, a1, b1);
    }
    
    if( Kzi == 0 ){
      if(this_zi==0){
        arma::colvec B = NindepIG_generator_sample_GIBBS(this_yy_colvec, m0, s0,  a0,  b0 );
        musig.row(i) = B.t();
      }else{
        arma::colvec A  = NIGmix_generator_sample(  this_yy_colvec, m1, V1,  a1, b1 );
        musig.row(i) =  A.t();
      }
      
    }else{
      
      arma::mat theta_zi  = LL["Pn"]; // estraggo qui perché non c'è in caso 0
      arma::vec nj_i      = LL["TI"];
      int    m            = theta_zi.n_rows;  
      arma::vec etichette = arma::linspace<arma::vec>(1, (m+1), m+1);
      arma::colvec logProbs(m+1), norm_Probs(m+1);  
      
      for(int k=0; k<(m); k++){
        
        logProbs[k] = log(nj_i(k)-sig_process) + R::dnorm(this_yy, 
                          theta_zi(k,0), 
                          sqrt(theta_zi(k,1)),
                          1);
        
      }
      logProbs[m] = log(th_process + sig_process * Kzi ) + log_Marginal_lik_value;
      double maxlp = max(logProbs);
      norm_Probs = exp( logProbs - maxlp) / sum(exp( logProbs - maxlp));
      norm_Probs = norm_Probs/sum(norm_Probs);
      int VALUE  =  RcppArmadillo::sample( etichette, 1, 1, norm_Probs)[0];
      
      
      if(VALUE == m+1){
        if(this_zi==0){
          arma::colvec B = NindepIG_generator_sample_GIBBS( this_yy_colvec, m0, s0,  a0,  b0 );
          NEW_THETA = B.t();
        }else{
          arma::colvec A = NIGmix_generator_sample(  this_yy_colvec, m1, V1,  a1, b1 );
          NEW_THETA = A.t();
        }
      }else{
        NEW_THETA = theta_zi.row(VALUE-1);
      }
      musig.row(i) = NEW_THETA;
    }
  }
  return(musig);
}


/////////////////////////////////ACCELERATION STEP

// [[Rcpp::export]]
arma::mat Accel_step(arma::mat Data, 
                     double m0, double s0, double a0, double b0,
                     double m1, double V1, double a1, double b1){
  arma::colvec yy  = Data.col(0);
  arma::colvec mu  = Data.col(1);
  arma::colvec sig = Data.col(2);
  arma::colvec zi  = Data.col(3);
  arma::colvec l0  = Data.col(4);
  arma::colvec l1  = Data.col(5);
  int LL0 = max(l0);
  int LL1 = max(l1);
  //Rcout << LL0 << "..." <<LL1 <<  "\n";
  // Accelerate Process 0
  if(LL0>0){
    for(int k0 = 0 ; k0 < LL0; k0++){
      arma::uvec    uk0  = find(l0 == (k0+1));
     // arma::colvec prev_sig    = (sig.elem(uk0)); //prendo minimo indice (il primo) e corrispondente valore di sigma che serve per simulare mu
      arma::colvec suby0 = yy.elem(uk0);
      int n0             = suby0.n_rows;
      // normal normal model
      //Rcout << "n0=" << n0 << " suby0 = "<< suby0 << "\n-----------\n";
      arma::colvec mu_fantoccio0(n0), sig_fantoccio0(n0), TH0(2);
      //TH0 = NindepIG_generator_sample(suby0, m0, s0, prev_sig[0], a0, b0);
      TH0 = NindepIG_generator_sample_GIBBS(suby0, m0, s0, a0, b0);
      sig.elem(uk0) = sig_fantoccio0.fill(TH0[1]) ;
      mu.elem(uk0)  = mu_fantoccio0.fill(TH0[0]) ;
    }}
  // Accelerate Process 1
  if(LL1>0){
    for(int k1 = 0 ; k1 < LL1; k1++){
      arma::uvec    uk1    = find(l1 == (k1+1));
      arma::colvec  suby1  = yy.elem(uk1);
      int n1               = suby1.n_rows;
      arma::colvec mu_fantoccio1(n1), sig_fantoccio1(n1), TH(2);
      TH = NIGmix_generator_sample(suby1, m1, V1, a1, b1);
      
      sig.elem(uk1) = sig_fantoccio1.fill(TH[1]) ;
      mu.elem(uk1)  = mu_fantoccio1.fill(TH[0]) ;
  //  Rcout << mu_fantoccio1<< "\n---\n";
    }}
  arma::mat RES  = arma::join_rows(mu,sig);
  return(RES);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////
// UPD m
/////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double log_NIGMIX(double m, double V,
                  arma::colvec unique_mu1, 
                  arma::colvec unique_sig1){
  
  arma::colvec invsig = 1/unique_sig1;
  double L = accu( log( 
    exp( -(unique_mu1-m)%(unique_mu1-m)/(2*V) % invsig) 
                    +  exp( -(unique_mu1+m)%(unique_mu1+m)/(2*V) % invsig) 
  ));
  return(L);    
}


// [[Rcpp::export]]
double ffm(double m, double k, double s){
  double a = pow(m,2*k) * exp(-(m * m)/(2 * s));
  return(a);
}

// [[Rcpp::export]]
double ABS_VAL(double m, double k, double s){
  double A = ffm(m,k,s) + ffm(-m,k,s);
  return(A);
}


// [[Rcpp::export]]
double unnorm_LOGPOST_MH_rev(double m, double V, double s1,
                             arma::colvec unique_mu1, 
                             arma::colvec unique_sig1, 
                             int kappa){
  return log_NIGMIX( m,  V,
                     unique_mu1, 
                     unique_sig1) +2*kappa*log(m)-m*m/(2*s1);
}

// [[Rcpp::export]]
double unnorm_LOGPOST_MH_rev_prior(double m, double s1, int kappa){
  return 2*kappa*log(m)-m*m/(2*s1);
}



// [[Rcpp::export]]
arma::colvec m_MH_mix(double previous_m, double V, double s1,
                      arma::colvec unique_mu1, 
                      arma::colvec unique_sig1, 
                      int kappa,double sigma_m){
  
  arma::colvec prev_m(2), prop_m(2); 
  
  prev_m[0]= previous_m; prev_m[1]=0.0;
  prop_m[0] = R::rnorm(prev_m[0], sqrt(sigma_m)); prop_m[1] = 1;
  
  double logratio = 
    unnorm_LOGPOST_MH_rev(prop_m[0], V, s1,
                          unique_mu1, unique_sig1, kappa)-  
                            unnorm_LOGPOST_MH_rev(prev_m[0], V, s1,
                                                  unique_mu1, unique_sig1, kappa);
  if( R::runif(0,1) < exp(logratio) ){
    return(prop_m);
  }else{
    return(prev_m);
  }
}




// [[Rcpp::export]]
arma::colvec m_MH_mix_prior(double previous_m, 
                            double s1,
                            int kappa,
                            double sigma_m){
  
  arma::colvec prev_m(2), prop_m(2); 
  prev_m[0]= previous_m; prev_m[1]=0.0;
  prop_m[0] = R::rnorm(prev_m[0], sqrt(sigma_m)); prop_m[1] = 1;
  
  double logratio = 
    unnorm_LOGPOST_MH_rev_prior(prop_m[0], s1, kappa)-
    unnorm_LOGPOST_MH_rev_prior(prev_m[0], s1, kappa);
  if( R::runif(0,1) < exp(logratio) ){
    return(prop_m);
  }else{
    return(prev_m);
  }
}





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////  ♫♫♫♫♫  Well, shake it up, baby, now
//////////////////////  ♫♫♫♫♫♫      Split and Merge             
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double log_ewens_PY(arma::colvec nk,double t0, double s0){
  int K = nk.n_elem;
  int n = accu(nk);
  double q1=0, q2=0, q3 = log_g_over_g(t0+1,n-1);
  for(int k=0; k<K-1; k++){
    q1 += log(t0+s0*(k+1));
    q2 += log_g_over_g( 1-s0, nk[k]-1);
  }
  q2 += log_g_over_g(1-s0, nk[K-1]-1);
  return q1+q2-q3;
  }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double log_Marginal_cancelsout0_SPLIT(arma::colvec ydata0, 
                                      arma::colvec all_lab0_cpp,
                                      arma::colvec labs_duegruppisplitted, 
                                      double a0, double b0){
  double logprob = 0;
  for(int j =0; j< 2; j++){
    arma::uvec   this_group_lab = find(all_lab0_cpp==labs_duegruppisplitted(j));
    arma::colvec currenty       = ydata0.elem(this_group_lab);
    logprob += log_Marginal_P0(currenty, a0,b0);
  } 
  return logprob;
}

// [[Rcpp::export]]
double log_Marginal_cancelsout0_MERGE(arma::colvec ydata0,
                                      arma::colvec all_lab0_cpp,
                                      int lab_gruppo_unito, 
                                      double a0, double b0){
  double logprob = 0;
    arma::uvec   this_group_lab = find(all_lab0_cpp==lab_gruppo_unito);
    arma::colvec currenty       = ydata0.elem(this_group_lab);
    logprob = log_Marginal_P0(currenty, a0,b0);
  return logprob;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double log_Marginal_cancelsout1_SPLIT(arma::colvec ydata1, 
                                      arma::colvec all_lab1_cpp,
                                      arma::colvec labs_duegruppisplitted, 
                                      double m1, double V1,
                                      double a1, double b1){
  
  double logprob = 0;
  for(int j =0; j< 2; j++){
    arma::uvec   this_group_lab = find(all_lab1_cpp==labs_duegruppisplitted(j));
    arma::colvec currenty       = ydata1.elem(this_group_lab);
    logprob += log_Marginal_P1(currenty, m1,V1,a1,b1);
  } 
  return logprob;
}

// [[Rcpp::export]]
double log_Marginal_cancelsout1_MERGE(arma::colvec ydata1,
                                      arma::colvec all_lab1_cpp,
                                      int lab_gruppo_unito, 
                                      double m1, double V1,
                                      double a1, double b1){
  
  double logprob = 0;
  arma::uvec   this_group_lab = find(all_lab1_cpp==lab_gruppo_unito);
  arma::colvec currenty       = ydata1.elem(this_group_lab);
  logprob += log_Marginal_P1(currenty, m1,V1,a1,b1);
  return logprob;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
int log_MH0_SPLITnew(arma::colvec yy0, 
                  arma::colvec labels0,
                  arma::colvec labels0new,
                  arma::colvec labsplit,
                  int labmerge,
                  arma::colvec nk,
                  arma::colvec nknew,
                  double log_prob_split,
                  double a0, double b0,double t0, double s0){
  int acc = 0;
  double NUM = log_ewens_PY(nknew,t0,s0) + log_Marginal_cancelsout0_SPLIT(yy0, labels0new,labsplit, a0, b0) ;
  double DEN = log_ewens_PY(nk   ,t0,s0) + log_Marginal_cancelsout0_MERGE(yy0, labels0,labmerge, a0, b0) + log_prob_split;
  
  //Rcout << NUM-DEN;
  
  if(runif(1)[0]<exp(NUM-DEN)){
    acc = 1 ; 
  }
  
  return acc;
}




// [[Rcpp::export]]
int log_MH1_SPLITnew(arma::colvec yy1, 
                  arma::colvec labels1,
                  arma::colvec labels1new,
                  arma::colvec nk,
                  arma::colvec labsplit,
                  int labmerge,
                  arma::colvec nknew,
                  double log_prob_split,
                    double m1, double V1,
                    double a1, double b1,double t1, double s1){
  int acc = 0;
  double NUM = log_ewens_PY(nknew,t1,s1) + log_Marginal_cancelsout1_SPLIT(yy1, labels1new,labsplit, m1, V1, a1, b1) ;
  double DEN = log_ewens_PY(nk   ,t1,s1) + log_Marginal_cancelsout1_MERGE(yy1, labels1,   labmerge, m1, V1, a1, b1) + log_prob_split;
  
  //Rcout << NUM-DEN;
  
  if(runif(1)[0]<exp(NUM-DEN)){
    acc = 1  ;
  }
  
  return acc;
}


// [[Rcpp::export]]
int log_MH0_MERGEnew(arma::colvec yy0, 
                  arma::colvec labels0,
                  arma::colvec labels0new,
                  arma::colvec nk,
                  arma::colvec nknew,
                  arma::colvec labsplit,
                  int labmerge,
                  double log_prob_split,
                  double a0, double b0,double t0, double s0){
  int acc = 0;
  double NUM = log_ewens_PY(nknew,t0,s0) + log_Marginal_cancelsout0_MERGE(yy0, labels0new,labmerge,a0, b0) + log_prob_split;
  double DEN = log_ewens_PY(nk   ,t0,s0) + log_Marginal_cancelsout0_SPLIT(yy0, labels0   ,labsplit,a0, b0)  ;
  
  //Rcout << NUM-DEN;
  
  if(runif(1)[0]<exp(NUM-DEN)){
    acc = 1  ;
  }
  
  return acc;
}


// [[Rcpp::export]]
int log_MH1_MERGEnew(arma::colvec yy1, 
                  arma::colvec labels1,
                  arma::colvec labels1new,
                  arma::colvec nk,
                  arma::colvec nknew,
                  arma::colvec labsplit,
                  int labmerge,
                  double log_prob_split,
                  double m1, double V1,
                  double a1, double b1,double t1, double s1){
  int acc = 0;
  double NUM = log_ewens_PY(nknew,t1,s1) + log_Marginal_cancelsout1_MERGE(yy1, labels1new,labmerge, m1, V1, a1, b1) + log_prob_split;
  double DEN = log_ewens_PY(nk   ,t1,s1) + log_Marginal_cancelsout1_SPLIT(yy1, labels1,   labsplit, m1, V1, a1, b1)  ;
  
   // Rcout << NUM-DEN;
  
  if(runif(1)[0]<exp(NUM-DEN)){
    acc = 1  ;
  }
  
  return acc;
}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// returns a new colum with updated positions memberships

// [[Rcpp::export]]
List Split0(arma::colvec yy, 
            arma::colvec indici_position_original_noij, //posizioni delle osservazioni da allocare
            int i_ind1, 
            int j_ind2, 
            double sig, 
            double a, 
            double b){
  
  int N  = yy.n_rows;
  int K  = indici_position_original_noij.n_rows;
  //Rcout<< N<<"---"<<K;
  arma::colvec  Si(N); Si.fill(-1);   Si(0) = i_ind1;
  arma::colvec  Sj(N); Sj.fill(-1);   Sj(0) = j_ind2;
  
  arma::colvec YSi(N), tempYSi(N); YSi.fill(-999); YSi(0) = yy(i_ind1-1);
  arma::colvec YSj(N), tempYSj(N); YSj.fill(-999); YSj(0) = yy(j_ind2-1);
  
  
  arma::colvec SHind_noij_cpp = indici_position_original_noij-1;
  SHind_noij_cpp = arma::shuffle(SHind_noij_cpp);

  arma::colvec logprobs(K);

  int n_i=1; int n_j=1;
  double currentlogprob, cum_logprob=0;
  

  for(int k=0; k<K; k++){
    
    int inde = SHind_noij_cpp[k];
    tempYSi = YSi; 
    tempYSj = YSj;
    //Rcout<< tempYSi<<"-----\n" << n_i<<"-----\n";
    tempYSi[n_i] = yy[inde];
    tempYSj[n_j] = yy[inde];
    //Rcout << "ok---";
    double logFactor = log_Marginal_P0(tempYSi.rows(0,n_i-1), a,b  )-log_Marginal_P0(tempYSj.rows(0,n_j-1), a,b  );
    
    currentlogprob = log(n_i-sig) + log_Marginal_P0(tempYSi.rows(0,n_i), a,b  ) - 
                     log( 
                       (n_i-sig) * exp(log_Marginal_P0(tempYSi.rows(0,n_i), a,b)) + 
                       (n_j-sig) * exp(log_Marginal_P0(tempYSj.rows(0,n_j), a,b)) * exp(logFactor)    
                        );

    if( runif(1)[0] < exp(currentlogprob)){
      YSi(n_i) = yy(inde);
      Si(n_i)  = inde+1;  
      n_i += 1;
      cum_logprob += currentlogprob;
    }else{
      YSj(n_j) = yy(inde);
      Sj(n_j)  = inde+1;
      n_j += 1;
      cum_logprob += log(1-exp(currentlogprob));
    }
  }

  
  return List::create(
    _["Si"]    = Si.rows(0,n_i-1), // mi ritorna le positions
    _["Sj"]    = Sj.rows(0,n_j-1),
//    _["YSi"]    = YSi.rows(0,n_i-1),
//    _["YSj"]    = YSj.rows(0,n_j-1),
    _["cumlogprob"]   = cum_logprob
    );
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List Merge0(arma::colvec yy, arma::colvec indici_position_original_noij,
            arma::colvec indici_labels_original, //all
            int i_ind1, int j_ind2, int i_lab1, int j_lab2, 
            double sig,  double a, double b){
  
  int N  = yy.n_rows;
  int K  = indici_position_original_noij.n_rows;
  //Rcout<< N<<"---"<<K;
  arma::colvec  Si(N); Si.fill(-1);   Si(0) = i_ind1;
  arma::colvec  Sj(N); Sj.fill(-1);   Sj(0) = j_ind2;
  
  arma::colvec YSi(N), tempYSi(N); YSi.fill(-999); YSi(0) = yy(i_ind1-1);
  arma::colvec YSj(N), tempYSj(N); YSj.fill(-999); YSj(0) = yy(j_ind2-1);
  
  
  arma::colvec SHind_noij_cpp = indici_position_original_noij-1;
  SHind_noij_cpp = arma::shuffle(SHind_noij_cpp);
  
  arma::colvec logprobs(K);
  
  int n_i=1; int n_j=1;
  double currentlogprob, cum_logprob=0;
  
  //Rcout << "ok2, K=" << K;
  
  for(int k=0; k<K; k++){
    
    int inde = SHind_noij_cpp[k];
    tempYSi = YSi; 
    tempYSj = YSj;
    tempYSi[n_i] = yy[inde];
    tempYSj[n_j] = yy[inde];
    double logFactor = log_Marginal_P0(tempYSi.rows(0,n_i-1), a,b  ) - log_Marginal_P0(tempYSj.rows(0,n_j-1), a,b  );
    currentlogprob = log(n_i-sig) + log_Marginal_P0(tempYSi.rows(0,n_i), a,b  ) - 
      log( 
          (n_i-sig) * exp(log_Marginal_P0(tempYSi.rows(0,n_i), a,b)) + 
          (n_j-sig) * exp(log_Marginal_P0(tempYSj.rows(0,n_j), a,b)) * exp(logFactor)    
      );
    
    
    if( indici_labels_original(inde) == i_lab1){
      YSi(n_i) = yy(inde);
      Si(n_i)  = inde+1;
      n_i += 1;
      cum_logprob += currentlogprob;
    }else{
      YSj(n_j) = yy(inde);
      Sj(n_j)  = inde+1;
      n_j += 1;
      cum_logprob += log(1-exp(currentlogprob));
    }
  }
  return List::create(
    _["Si"]   = Si.rows(0,n_i-1),
    _["Sj"]   = Sj.rows(0,n_i-1),
    _["cumlogprob"]   = cum_logprob
  );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List Split1(arma::colvec yy, arma::colvec indici_position_original_noij,
            int i_ind1, int j_ind2, 
            double sig,  double m1, double V1,double a1, double b1){
  
  int N  = yy.n_rows;
  int K  = indici_position_original_noij.n_rows;
  //Rcout<< N<<"---"<<K;
  arma::colvec  Si(N); Si.fill(-1);   Si(0) = i_ind1;
  arma::colvec  Sj(N); Sj.fill(-1);   Sj(0) = j_ind2;
  
  arma::colvec YSi(N), tempYSi(N); YSi.fill(-999); YSi(0) = yy(i_ind1-1);
  arma::colvec YSj(N), tempYSj(N); YSj.fill(-999); YSj(0) = yy(j_ind2-1);
  
  
  arma::colvec SHind_noij_cpp = indici_position_original_noij-1;
  SHind_noij_cpp = arma::shuffle(SHind_noij_cpp);
  
  arma::colvec logprobs(K);
  
  int n_i=1; int n_j=1;
  double currentlogprob, cum_logprob=0;
  
  //Rcout << "ok2, K=" << K;
  
  for(int k=0; k<K; k++){
    
    int inde = SHind_noij_cpp[k];
    tempYSi = YSi; 
    tempYSj = YSj;
    tempYSi[n_i] = yy[inde];
    tempYSj[n_j] = yy[inde];
    //Rcout << "ok---";
    double logFactor = log_Marginal_P1(tempYSi.rows(0,n_i-1), m1, V1, a1, b1)-log_Marginal_P1(tempYSj.rows(0,n_j-1), m1, V1, a1, b1);
    
    currentlogprob = log(n_i-sig) + log_Marginal_P1(tempYSi.rows(0,n_i), m1, V1, a1, b1)- 
      log( 
          (n_i-sig) * exp(log_Marginal_P1(tempYSi.rows(0,n_i), m1, V1, a1, b1)) + 
          (n_j-sig) * exp(log_Marginal_P1(tempYSj.rows(0,n_j), m1, V1, a1, b1)) * exp(logFactor)    
      );
    if( runif(1)[0] < exp(currentlogprob)){
      YSi(n_i) = yy(inde);
      Si(n_i)  = inde+1;
      n_i += 1;
      cum_logprob += currentlogprob;
    }else{
      YSj(n_j) = yy(inde);
      Sj(n_j)  = inde+1;
      n_j += 1;
      cum_logprob += log(1-exp(currentlogprob));
    }
  }
  
  
  return List::create(
    _["Si"]    = Si.rows(0,n_i-1), // mi ritorna le positions
    _["Sj"]    = Sj.rows(0,n_j-1),
//    _["YSi"]    = YSi.rows(0,n_i-1),
//    _["YSj"]    = YSj.rows(0,n_j-1),
    _["cumlogprob"]   = cum_logprob
  );
  
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List Merge1(arma::colvec yy, arma::colvec indici_position_original_noij,
            arma::colvec indici_labels_original, //all
            int i_ind1, int j_ind2, int i_lab1, int j_lab2, 
            double sig,  double m1, double V1, double a1, double b1){
  
  
  int N  = yy.n_rows;
  int K  = indici_position_original_noij.n_rows;
  //Rcout<< N<<"---"<<K;
  arma::colvec  Si(N); Si.fill(-1);   Si(0) = i_ind1;
  arma::colvec  Sj(N); Sj.fill(-1);   Sj(0) = j_ind2;
  
  arma::colvec YSi(N), tempYSi(N); YSi.fill(-1); YSi(0) = yy(i_ind1-1);
  arma::colvec YSj(N), tempYSj(N); YSj.fill(-1); YSj(0) = yy(j_ind2-1);
  
  
  arma::colvec SHind_noij_cpp = indici_position_original_noij-1;
  SHind_noij_cpp = arma::shuffle(SHind_noij_cpp);
  
  arma::colvec logprobs(K);
  
  int n_i=1; int n_j=1;
  double currentlogprob, cum_logprob=0;
  
  //Rcout << "ok2, K=" << K;
  
  for(int k=0; k<K; k++){
    
    int inde = SHind_noij_cpp[k];
    tempYSi = YSi; 
    tempYSj = YSj;
    tempYSi[n_i] = yy[inde];
    tempYSj[n_j] = yy[inde];
    double logFactor = log_Marginal_P1(tempYSi.rows(0,n_i-1), m1, V1, a1, b1)-log_Marginal_P1(tempYSj.rows(0,n_j-1), m1, V1, a1, b1);
    //Rcout << "ok---";
    currentlogprob = log(n_i-sig) + log_Marginal_P1(tempYSi.rows(0,n_i), m1, V1, a1, b1)- 
      log( 
        (n_i-sig) * exp(log_Marginal_P1(tempYSi.rows(0,n_i), m1, V1, a1, b1)) + 
          (n_j-sig) * exp(log_Marginal_P1(tempYSj.rows(0,n_j), m1, V1, a1, b1)) * exp(logFactor)    
      );
    
    if( indici_labels_original(inde) == i_lab1){
      YSi(n_i) = yy(inde);
      Si(n_i)  = inde+1;
      n_i += 1;
      cum_logprob += currentlogprob;
    }else{
      YSj(n_j) = yy(inde);
      Sj(n_j)  = inde+1;
      n_j += 1;
      cum_logprob += log(1-exp(currentlogprob));
    }
  }
  
  return List::create(
    _["cumlogprob"]   = cum_logprob
  );
  
}























































/*
// [[Rcpp::export]]
double log_Marginal_alldata0(arma::colvec ydata0, 
                             arma::colvec all_lab0, 
                             double a0, double b0){
  int MaxLabel = all_lab0.max();
  arma::colvec all_lab0_cpp = all_lab0 - 1;
  double logprob = 0;
  for(int j =0; j< MaxLabel; j++){
    arma::uvec   this_group_lab = find(all_lab0_cpp==j);
    arma::colvec currenty       = ydata0.elem(this_group_lab);
    // Rcout << log_Marginal_P0(currenty, a0,b0) << ".......\n";
    logprob += log_Marginal_P0(currenty, a0,b0);
  } 
  return logprob;
}
*/


/*
 
 // [[Rcpp::export]]
 double log_Marginal_alldata1(arma::colvec ydata1, 
 arma::colvec all_lab1, 
 double m1, double V1,
 double a1, double b1){
 int MaxLabel = all_lab1.max();
 arma::colvec all_lab1_cpp = all_lab1 - 1;
 double logprob = 0;
 for(int j =0; j< MaxLabel; j++){
 arma::uvec   this_group_lab = find(all_lab1_cpp==j);
 arma::colvec currenty       = ydata1.elem(this_group_lab);
 logprob += log_Marginal_P1(currenty, m1,V1,a1,b1);
 }
 return logprob;
 }
*/


/*
 // [[Rcpp::export]]
 int log_MH0_MERGE(arma::colvec yy0, 
 arma::colvec labels0,
 arma::colvec labels0new,
 arma::colvec nk,
 arma::colvec nknew,
 double log_prob_split,
 double a0, double b0,double t0, double s0){
 int acc = 0;
 
 double NUM = log_ewens_PY(nknew,t0,s0) + log_Marginal_alldata0(yy0, labels0new, a0, b0) + 0;
 double DEN = log_ewens_PY(nk   ,t0,s0) + log_Marginal_alldata0(yy0, labels0   , a0, b0) + log_prob_split;
 Rcout << NUM-DEN;
 if(runif(1)[0]<exp(NUM-DEN)){
 acc = 1  ;
 }
 
 return acc;
 }
 */





/*
 // [[Rcpp::export]]
 int log_MH1_MERGE(arma::colvec yy1, 
 arma::colvec labels1,
 arma::colvec labels1new,
 arma::colvec nk,
 arma::colvec nknew,
 double log_prob_split,
 double m1, double V1,
 double a1, double b1,double t1, double s1){
 int acc = 0;
 double NUM = log_ewens_PY(nknew,t1,s1) + log_Marginal_alldata1(yy1, labels1new, m1, V1, a1, b1) + 0;
 double DEN = log_ewens_PY(nk   ,t1,s1) + log_Marginal_alldata1(yy1, labels1   , m1, V1, a1, b1) + log_prob_split;
 
 Rcout << NUM-DEN;
 
 if(runif(1)[0]<exp(NUM-DEN)){
 acc = 1  ;
 }
 
 return acc;
 }
 */




/*
 // [[Rcpp::export]]
 int log_MH0_SPLIT(arma::colvec yy0, 
 arma::colvec labels0,
 arma::colvec labels0new,
 arma::colvec nk,
 arma::colvec nknew,
 double log_prob_split,
 double a0, double b0,double t0, double s0){
 int acc = 0;
 double NUM = log_ewens_PY(nknew,t0,s0) + log_Marginal_alldata0(yy0, labels0new, a0, b0) + log_prob_split;
 double DEN = log_ewens_PY(nk   ,t0,s0) + log_Marginal_alldata0(yy0, labels0   , a0, b0) + 0;
 
 Rcout << NUM-DEN;
 
 if(runif(1)[0]<exp(NUM-DEN)){
 acc = 1 ; 
 }
 
 return acc;
 }
 // Functions for generating random variates from NIG
 // [[Rcpp::export]]
 arma::colvec NIGmix_generator_prior(double m1, double V1, double a1, double b1){
 double mn;
 if( R::runif(0,1) < .5 ){
 mn = m1;
 }else{
 mn = -m1;
 }
 arma::colvec TH(2);
 TH[1] = 1/ R::rgamma( a1, 1/(b1) );
 TH[0] = R::rnorm(mn, sqrt( TH[1]*V1  ));
 return TH;
 }
 
 
 // [[Rcpp::export]]
 double Rgamma_rcpp(double a, double b){
 return  R::rgamma(a,1/b);
 }
 
 
 
  */

