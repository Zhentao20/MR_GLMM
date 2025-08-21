//[[Rcpp::depends("RcppArmadillo")]]
#include <armadillo>
//#include <boost/math/distributions/normal.hpp>//
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iomanip>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat rw(double d){
  arma::mat z = arma::randn <arma::mat> (d,d,arma::distr_param(0,1));
  return(z);
}

// [[Rcpp::export]]
arma::mat gs(arma::mat u, double d){
  arma::mat z = u + rw(d);
  return(z);
}


// [[Rcpp::export]]
List meanii(List meani, int ti, arma::mat Ut, arma::mat lambda, arma::mat Vt, List xi, List b, double bo){
  arma::mat temp;
  for (int k = 0; k < ti; k++) {
    if (bo > 0){
      temp = Ut * lambda * Ut.t();
    }
    else{
      temp = Ut * Vt.t();
    }    
      List xit = xi[k];
      int p = xit.size();
      
      for (int q = 0; q < p; q++){
        double xitq = xit[q];
        arma::mat biq = as<arma::mat>(b[q]);
        temp += xitq * biq;
      }

    meani[k] = temp;
  }
  return meani;
 }



// [[Rcpp::export]]
List Meani(List mi,arma::mat ui){
  // n =5, indicates time
  int n = mi.size();
  List z(n);
  for(int j=0; j<n; j++){
    z[j] = as<arma::mat>(mi[j]) + ui;
  }
  return(z);
}


// [[Rcpp::export]]
arma::vec matrix2vector(arma::mat m, const bool byrow=false){
  if (byrow) {
    return m.as_row();
  } else {
    return m.as_col();
  }
}

// [[Rcpp::export]]
NumericVector mydnorm(arma::vec x, arma::vec means, arma::vec sds){
  int n = x.size() ;
  NumericVector res(n) ;
  for( int i=0; i<n; i++) res[i] = R::dnorm(x[i], means[i], sds[i], TRUE);
  return res ;
}

// [[Rcpp::export]]
double mydnorm1(double x, double means, double sds){
  double res = R::dnorm(x, means, sds, TRUE);
  return res ;
}

// [[Rcpp::export]]
double mydnorm2(double x, double means, double sds) {
  double variance = pow(sds, 2);
  double exponent = -pow(x - means, 2) / (2 * variance);
  double pdf = (1 / (sds * std::sqrt(2 * M_PI))) * std::exp(exponent);
  double res = std::log(pdf);
  return res;
}

// [[Rcpp::export]]
NumericVector mylogit (arma::vec x, arma::vec means){
  int n = x.size() ;
  NumericVector res(n) ;
  for( int i=0; i<n; i++) res[i] = x[i]*log(exp(means[i])/(1+exp(means[i])))+(1-x[i])*log(1-(exp(means[i])/(1+exp(means[i]))));
  return res ;
}

// [[Rcpp::export]]
double ab(List mi, List zi, arma::mat ui, arma::mat sgammat, double d){
  // n indicates time points = 5
  int n = mi.size();
  NumericVector z(n);
  for(int i = 0; i<n;i++){
    z[i] = sum(mylogit(matrix2vector(zi[i]), matrix2vector(mi[i])));
  }
  double result = sum(z) + sum(mydnorm(matrix2vector(ui), matrix2vector(arma::zeros<arma::mat>(d, d)), matrix2vector(sgammat)));
  return result;
//  Rcpp::Rcout << "result of ab: " << result << std::endl;
}

// // [[Rcpp::export]]
// double ab(List mi, List zi, arma::mat ui, arma::mat sgammat, arma::mat seti){
//   int n = mi.size();
//   NumericVector z(n);
//   for(int i = 0; i<n;i++){
//     z[i] = sum(mydnorm(matrix2vector(zi[i]), matrix2vector(mi[i]), matrix2vector(seti)));
//   }
//   double result = sum(z) + sum(mydnorm(matrix2vector(ui), matrix2vector(set.zeros()), matrix2vector(sgammat)));
//   return result;
// }

// [[Rcpp::export]]
arma::cube sampledi (List zi, List xi, int ti, double M, arma::mat Ut, arma::mat lambda, arma::mat Vt, List b, arma::mat sgammat, double d, double bo){
  // double count = 0;
  // double n= 30;
  arma::cube u = arma::cube(d,d,M);
  /*check*/
  u.slice(0) = arma::randn(d,d);
  List meani(ti);
  for (int i = 0; i < ti; i++) {
    meani[i] = arma::mat(d, d, arma::fill::zeros);

  }  
// used to be i < M - 1
  for(int i=0; i< M-1;i++){
    u.slice(i+1) = gs(u.slice(i),d);
    List z = meani;
    List meani = meanii(z, ti, Ut,lambda, Vt, xi, b, bo);
    List meandi = Meani(meani,u.slice(i));
    /*check*/
    double a = ab(meandi, zi, u.slice(i), sgammat,d);
    List meandi1 = Meani(meani,u.slice(i+1));
    double b = ab(meandi1, zi, u.slice(i+1),sgammat,d);
    double R = exp(b-a);
    NumericVector ratio = NumericVector::create(1, R);
    double r = min(ratio);
    NumericVector keep = Rcpp::rbinom(1,1,r);
    if (max(keep)>=1){
      u.slice(i+1)=u.slice(i+1);
    }else{
      u.slice(i+1)=u.slice(i);
    }
  }
  return u;
}


// [[Rcpp::export]]
List sampleall( List A, List Xt, arma::vec t, double M, arma::mat Ut, arma::mat lambda, arma::mat Vt, List b, arma::mat sgammat, double N, double d, double bo) {
  Rcout << "N in sampleall is: " << N << std::endl;
  List samples(N);
  for (int i = 0; i < N; i++) {
    int ti = t(i);
    List xi(ti);
    List Ai = A[i];
   for(int T = 0; T < ti; T++){
     xi[T] = Xt[i*5+T];
   }
   arma::cube samples_i = sampledi(Ai, xi, ti, M, Ut, lambda, Vt, b, sgammat, d,bo);
//   arma::mat test = samples_i.slice(1);
//   Rcpp::Rcout << "test in sampleall: " << test << std::endl;
   samples[i] = samples_i;
  }
  return samples;
}

// [[Rcpp::export]]
double qq(List Ai, List xi, int ti,double M, arma::mat Ut, arma::mat lambda, arma::mat Vt, List b, arma::mat sgammat, arma::cube gammai, double gam, double d, double bo) {
  int T = xi.size();
//  Rcpp::Rcout << "T in qq: " << T << std::endl;
  
  List meani(T);
  int gslice = gammai.n_slices; //100 data replicates for each time (individual subject)
  int M_dim = T*gslice; // so total is 5 * 100 = 500 replicates
  arma::mat zero = arma::zeros<arma::mat>(d, d);
  arma::mat temp;
  arma::cube listM(d,d,M_dim);
  for (int k = 0; k < ti; k++) {
    if (bo > 0){
      temp = Ut * lambda * Ut.t();
    }
    else{
      temp = Ut * Vt.t();
    }    
    List xit = xi[k];
    int p = xit.size();
//    Rcpp::Rcout << "the p in qq is" << p << std::endl;
    
    for (int q = 0; q < p; q++){
      double xitq = xit[q];
      arma::mat biq = as<arma::mat>(b[q]);
      temp += xitq * biq;
    }
    meani[k] = temp;
  }    
  for(int i = 0; i < gslice; i++){
    for (int j = 0; j < ti; j++) {
      arma::mat meanit = meani[j];
      listM.slice(ti*i+j) = meanit + gammai.slice(i);
    }
  }
  arma::vec store = arma::zeros<arma::vec>(M_dim);

  for (int i = 0; i < gslice; i++){
    for (int j = 0; j < ti; j++){
      store(ti*i+j) =  -sum(mylogit(matrix2vector(Ai[j]), matrix2vector(listM.slice(ti*i+j))));
    }
  }
  
  double sumstore = arma::sum(store);
  arma::vec errorstore = arma::zeros<arma::vec>(M_dim);
  for (int i = 0; i < gslice; i++){
    for (int j = 0; j < ti; j++){
      errorstore(ti*i+j) =  -sum(mydnorm(matrix2vector(gammai.slice(i)), matrix2vector(zero), matrix2vector(sgammat)));
    }
  }
  double sumerror = arma::sum(errorstore);
  double res;
  if (bo > 0){
    res = gam * std::pow(arma::norm(trans(Ut)*Ut-trans(Vt)*Vt,"fro"),2);
  }
  else{
    res = 0;   
  }  
  double qi = (sumstore + sumerror)/M + res;
  return qi;
//  Rcout <<  std::fixed << std::setprecision(10) << "The qi is: "<<  qi << std::endl;
  
}


// Define the Q function
// [[Rcpp::export]]
double Q(List A, List Xt,arma::vec t, arma::mat Ut, arma::mat lambda, arma::mat Vt, List b, arma::mat sgammat, List samples, double M, double gam, double d,double bo) {
  int N = A.length();
  arma::cube Q(1, 1, N, arma::fill::none);
  for (int i = 0; i < N; i++) {
    List Ai = A[i];
    int ti = t(i);
    List xi(ti);
    for (int T = 0; T < ti; T++){
      xi[T] =Xt[5*i+T];
    }
    arma::cube gammai = samples[i];
    Q.slice(i) = qq(Ai, xi, ti, M, Ut,lambda, Vt, b, sgammat, gammai, gam, d, bo);
  }
  double result = accu(Q) / (N * d * d );
  return result;
}


// [[Rcpp::export]]
arma::vec vector_minus_mean(arma::vec v) {
  double mean_v = arma::mean(v); // Calculate mean of the vector
  arma::vec result = v - mean_v; // Subtract mean from each element in the vector
  return result;
}

// [[Rcpp::export]]
arma::vec joinAndVectorize(std::vector<arma::mat> mats) {
  arma::mat result = mats[0];
  for(arma::uword i = 0; i < mats.size(); i++) {
    result = arma::join_rows(result, mats[i]);
  }
  arma::vec vecResult = arma::vectorise(result);
  return vecResult;
}

// [[Rcpp::export]]
arma::mat trun(arma::mat b_ages, double s, double d){
  if(s>0){
  arma::mat b_age_keep = b_ages;
  arma::mat absbages = abs(b_ages);
  arma::uvec sorted_indices_age = arma::sort_index(absbages, "descend"); 
  arma::uvec top_indices_age = sorted_indices_age.head(s*d*d);
  b_ages = arma::zeros<arma::mat>(d, d);
  b_ages.elem(top_indices_age) = b_age_keep.elem(top_indices_age);
  }
  return b_ages;
}

// [[Rcpp::export]]
List output(List Y, List X, List A, List Xt, arma::vec t, double M, List bt, arma::mat sgamma0, arma::mat sgammatrue, double eps, double Eps, double tol, double iter, double maxit, double lhs, arma::mat Abar, arma::mat lambda, List b0, double tol1, double iter1, double maxit1, double maxit2, double s, arma::mat utrue, arma::mat vtrue, List btrue, double gam, int r, double N, double p, double d,  double step1, double step2){
  arma::vec Si(r);
  arma::mat Ui(Abar.n_rows, r), Vi(Abar.n_cols, r);

  arma::sp_mat Ab = arma::sp_mat(Abar);

  arma::svds(Ui, Si, Vi, Ab, r);
  
  arma::vec s0 = Si;
  arma::mat U0 = Ui * arma::diagmat(sqrt(s0));
  arma::mat V0 = Vi * arma::diagmat(sqrt(s0));
  // Rcout << "The initial value of s0 is: " << s0 << std::endl;
  // Rcout << "The initial value of V0 is: " << V0 << std::endl;
  
  arma::mat Ut = U0;
  arma::mat Vt = V0;
  arma::mat sgammat = sgammatrue;
  arma::vec record_q;
  arma::vec record_dif_thetaf;
  
  double btn = bt.size();
  List record_dif_b(btn); // Vector of arma::vec to store the differences
  
  for (int i = 0; i < btn; ++i) {
    //    arma::vec vec(10, arma::fill::zeros); // a vector of 10 elements, filled with zeros
    record_dif_b[i] = arma::vec();
  }  
  
  arma::vec record_step_theta;
  arma::vec record_step_b;
  arma::vec record_eps;
  while (eps > tol && iter < maxit) {
    // double qfunction0 = lhs;
    List samples = sampleall(A, Xt, t, M, Ut, lambda, Vt, bt, sgammat, N,d,0);
    double qfunction = Q(A, Xt, t, Ut, lambda, Vt, bt, sgammat, samples, M, gam, d, 0);
    Rcout <<"The result of Q is: " << qfunction << std::endl;
    double eps1 = arma::norm(Ut*Vt.t()-U0*V0.t(),"fro");
    Rcout << std::fixed<< std::setprecision(10) <<"The eps1 is: " << eps1 << std::endl;
    int bn = bt.size();
    arma::vec epsa(bn);
    for (int i = 0; i < bn; ++i){
      arma::mat bti = bt[i];
      arma::mat b0i = b0[i];
      double epsi = arma::norm(bti-b0i,"fro");
      epsa(i) = epsi;
      Rcout << std::fixed << std::setprecision(5) << "The eps" << (i + 2) << " is: " << epsa(i) << std::endl;
    }
    double max_val = eps1;
    for(arma::uword i = 0; i < epsa.size(); i++) {
      max_val = std::max(max_val, epsa(i));
    }
//0211_change eps tp max_val
    Rcout << std::setprecision(10) <<"The eps at beginning of each EM is: " << max_val << std::endl;
    double approxq = 10000;
    Rcout << std::setprecision(10) <<"The approxq at beginning of each EM is: " << approxq << std::endl;
    
    int iter1 = 0;
    Rcout << std::setprecision(10) <<"The iter1 is: " << iter1 << std::endl;
    double label = std::abs(lhs - approxq)/std::abs(approxq);
    Rcout << std::setprecision(10) <<"The abs(lhs - approxq)/abs(approxq) is: " << label << std::endl;
    
    
    while (std::abs(lhs - approxq)/std::abs(approxq) > tol1 && iter1 < maxit1) {
          arma::mat id = arma::mat(d,d).fill(1);
          approxq = Q(A, Xt, t, Ut, lambda, Vt, bt, sgammat, samples, M, gam, d, 0);
          Rcout  << "The approxq at beginning of inner 1 is: " << approxq << std::endl;
          double denom = N * d * d;
          arma::mat hu = arma::zeros(d, r);
          arma::mat hv = arma::zeros(d, r);
          List hb(bn);
          for (int i = 0; i < bn; ++i){
            arma::mat hbi = arma::zeros(d, d);
            hb[i] = hbi;
          }      
          arma::mat hsg = arma::zeros(d, d);
          List hse(N);
          for (int q = 0; q < bn; ++q) {
            arma::vec xq = X[q];
            arma::mat hbq_adjustment = arma::zeros(d, d); // Assuming xq is a vector
            
            for (int i = 0; i < N; ++i) {
              arma::cube samplesi = samples[i];
              arma::mat samplesi_mean = arma::mean(samplesi, 2);
//             for (int m = 0; m < M; ++m){
//               arma::mat samplesim = samplesi.slice(m);
              for (int T = 0; T < t(i); ++T) {
                arma::cube Yt = Y[T];
                arma::mat Yti = arma::mat(Yt.slice(i));
                arma::mat fm = Ut * trans(Vt) + samplesi_mean;
                double xqit;
                for(int inner_q = 0; inner_q < bn; ++inner_q) {
                  arma::vec inner_xq = X[inner_q];
                  arma::vec xqi = inner_xq.subvec(i*t(i),i*t(i)+t(i)-1);
                  xqit = xqi(T);
                  arma::mat bq = as<arma::mat>(bt[inner_q]);
                  fm += xqit * bq;
                }
                arma::vec outer_xqi = xq.subvec(i*t(i),i*t(i)+t(i)-1);
                double outer_xqit = outer_xqi(T);
                hu -= ((Yt.slice(i) - exp(fm)/(id + exp(fm))) * Vt) / (denom*bn)  - (4*gam*Ut*(trans(Ut)*Ut-trans(Vt)*Vt))/ (denom*bn);
                hv -= (trans(Yt.slice(i) - exp(fm)/(id + exp((fm)))) * Ut) / (denom*bn)  - (4*gam*Vt*(trans(Ut)*Ut-trans(Vt)*Vt))/ (denom*bn);
                hbq_adjustment -= ((Yt.slice(i) - exp(fm)/(id + exp(fm)) ) * outer_xqit) / (denom) ;
          // hu -= ((Yt.slice(i) + exp((Ut * trans(Vt)) - samplesi_mean - X0(i) * b_sext - X1(i) * b_aget)/(id + exp((Ut * trans(Vt)) - samplesi_mean - X0(i) * b_sext - X1(i) * b_aget))) * Vt) / denom  - (4*gam*Ut*(trans(Ut)*Ut-trans(Vt)*Vt))/ (denom * M);
          // hv -= (trans(Yt.slice(i) + exp((Ut * trans(Vt)) - samplesi_mean - X0(i) * b_sext - X1(i) * b_aget)/(id + exp((Ut * trans(Vt)) - samplesi_mean - X0(i) * b_sext - X1(i) * b_aget))) * Ut) / denom  - (4*gam*Vt*(trans(Ut)*Ut-trans(Vt)*Vt))/ (denom * M);
          // hb1 -= ((Yt.slice(i) + exp((Ut * trans(Vt)) - samplesi_mean - X0(i) * b_sext - X1(i) * b_aget)/(id + exp((Ut * trans(Vt)) - samplesi_mean - X0(i) * b_sext - X1(i) * b_aget)) ) * X0(i)) / (denom * M) ;
          // hb2 -= ((Yt.slice(i) + exp((Ut * trans(Vt)) - samplesi_mean - X0(i) * b_sext - X1(i) * b_aget)/(id + exp((Ut * trans(Vt)) - samplesi_mean - X0(i) * b_sext - X1(i) * b_aget)) ) * X1(i)) / (denom * M) ;
           }
        }
          hb[q] = hbq_adjustment;
      }
          // arma::mat Ug = hu.submat(0, 0, r-1, r-1);
          // arma::mat Vg = hv.submat(0, 0, r-1, r-1);
          // arma::mat hb1 = hb[0];
          // arma::mat hb3 = hb[2];
          // arma::mat hb5 = hb[4];
          // 
          // arma::mat b1g = hb1.submat(0, 0, 4, 4);
          // arma::mat b3g = hb3.submat(0, 0, 4, 4);
          // arma::mat b5g = hb5.submat(0, 0, 4, 4);
          // 
          // Rcout << "The initial value of sub-hu of inner 1 is: " << Ug << std::endl;
          // Rcout << "The initial value of sub-hv of inner 1 is: " << Vg << std::endl;
          // Rcout << "The initial value of sub-hb1 of inner 1 is: " << b1g << std::endl;
          // Rcout << "The initial value of sub-hb3 of inner 1 is: " << b3g << std::endl;
          // Rcout << "The initial value of sub-hb5 of inner 1 is: " << b5g << std::endl;

      U0 = Ut;
      V0 = Vt;
      b0 = bt;

      double a = step1;
      double c = step2;
      
      std::vector<arma::mat> dif1 = {-a*hu, -a*hv}; 
      std::vector<arma::mat> gradient1 = {hu, hv}; 
      
      std::vector<arma::mat> dif2;
      List trunb(bn);
      for(int i = 0; i < hb.size(); ++i) {
        arma::mat bti = bt[i];
        arma::mat hbi = hb[i];
//        arma::mat mat = trun(bti-(1/c)*hbi,s,d);
        arma::mat mat = trun(bti-c*hbi,s,d);
        dif2.push_back(mat);
        trunb[i] = mat;
      }
      
      std::vector<arma::mat> gradient2;
      for(int i = 0; i < hb.size(); ++i) {
        arma::mat mat = hb[i];
        gradient2.push_back(mat);
      }
      arma::mat res2 =  trans(joinAndVectorize(dif1))*joinAndVectorize(dif1);
      double Res2 = res2(0);
      
      arma::mat res4 =  trans(joinAndVectorize(dif2))*joinAndVectorize(dif2);
      double Res4 = res4(0);
      
      while (Q(A, Xt, t, U0 - a * hu, lambda, V0 - a * hv, bt, sgammat,  samples, M, gam, d, 0) - approxq - (1/(2.0*a)*Res2) > 0) {
        arma::vec new_rowt(1);
        a *= 0.5;
        new_rowt(0) = a;  
        record_step_theta = join_cols(record_step_theta, new_rowt);
        Rcout << "The record_step_theta of inner 2 is: " << record_step_theta << std::endl;
      }
      
      Ut = U0 - a * hu;
      Vt = V0 - a * hv;
      // arma::mat Uu = Ut.submat(0, 0, r-1, r-1);
      // arma::mat Vu = Vt.submat(0, 0, r-1, r-1);
      // Rcout << "The updated U of inner 1 is: " << Uu << std::endl;
      // Rcout << "The updated V of inner 1 is: " << Vu << std::endl;
      if (s >0){
      while (Q(A, Xt, t, Ut,lambda, Vt, trunb, sgammat,  samples, M, gam, d, 0) - approxq - (1/(2.0*c)*Res4)  > 0){
        arma::vec new_rowb(1);
        c *= 0.8;
        new_rowb(0) = c;  
        record_step_b = join_cols(record_step_b, new_rowb);
//        Rcout << "The record_step_b of inner 2 is: " << record_step_b << std::endl;
      }
      for(int i = 0; i < hb.size(); ++i) {
        arma::mat bti = bt[i];
        arma::mat hbi = hb[i];
        arma::mat btruei = btrue[i];
        double nc = c * 2 * (arma::norm(bti - btruei, "fro"));
        arma::mat mat = trun(bti-nc*hbi,s,d);
        dif2.push_back(mat);
        trunb[i] = mat;
      }
      bt = trunb;
      }
      else{
        bt=b0;
      }
      // arma::mat b1 = bt[0];
      // arma::mat b3 = bt[2];
      // arma::mat b5 = bt[4];
      // arma::mat b1u = b1.submat(0, 0, 4, 4);
      // arma::mat b3u = b3.submat(0, 0, 4, 4);
      // arma::mat b5u = b5.submat(0, 0, 4, 4);
      // 
      // Rcout << "The updated b1 of inner 1 is: " << b1u << std::endl;
      // Rcout << "The updated b3 of inner 1 is: " << b3u << std::endl;
      // Rcout << "The updated b5 of inner 1 is: " << b5u << std::endl;
      lhs = Q(A, Xt, t, Ut, lambda, Vt, bt, sgammat, samples, M, gam, d, 0);
      Rcout << std::setprecision(10) <<"The lhs at INNER (test) is: " << lhs << std::endl;
      iter1++;
      double rdif = std::abs(lhs - approxq)/std::abs(approxq);
      Rcout << std::setprecision(10) <<"The abs(lhs - approxq)/abs(approxq) of inner 1 is: " << rdif  << std::endl;
      
      
      if(rdif < tol1){
        break;
      }  
    }
    iter++;

    approxq = 10000;

    arma::vec new_row1(1);
    new_row1(0) = qfunction;  
    record_q = arma::join_cols(record_q, new_row1);
    Rcout << "The record_q is: " << record_q << std::endl;
    Rcout << std::fixed<< std::setprecision(10) <<"The eps at the end is: " << eps << std::endl;
    Rcout << "The iteration of outer loop is: " << iter << std::endl;
    
    arma::vec new_row2(1);
    new_row2(0) = arma::norm(Ut*Vt.t()-utrue*vtrue.t(), "fro");   
    record_dif_thetaf = arma::join_cols(record_dif_thetaf, new_row2);
    Rcout << std::fixed<< std::setprecision(10) <<"The record_dif_thetaf is: " << record_dif_thetaf << std::endl;
    
    for(int i = 0; i < btn; ++i) {
      arma::mat bi = bt[i];
      arma::mat btruei = btrue[i];
      arma::vec new_row77(1);
      new_row77(0) = arma::norm(bi - btruei, "fro");
      arma::vec rdb = record_dif_b[i];
      record_dif_b[i] = arma::join_cols(rdb, new_row77);
      
      arma::vec temp = record_dif_b[i];
//      Rcout << std::fixed << std::setprecision(5) << "The record_dif_b" << (i + 1) << " is: " << temp << std::endl;
    }

    eps1 = arma::norm(Ut*Vt.t()-U0*V0.t(),"fro");
    Rcout << std::fixed<< std::setprecision(10) <<"The eps1 at the end of each EM is: " << eps1 << std::endl;
    
    for (int i = 0; i < bn; ++i){
      arma::mat bti = bt[i];
      arma::mat b0i = b0[i];
      double epsi = arma::norm(bti-b0i,"fro");
      epsa(i) = epsi;
//      Rcout << std::fixed << std::setprecision(5) << "The eps" << (i + 2) << " at the end of each EM is: " << epsa(i) << std::endl;
    }
    
    double Max_val = eps1;
    for(arma::uword i = 0; i < epsa.size(); i++) {
      Max_val = std::max(Max_val, epsa(i));
    }
    Rcout << std::setprecision(10) <<"The eps at the end of each EM is: " <<  Max_val << std::endl;
    arma::vec new_row7(1);
    Rcout << std::fixed<< std::setprecision(10)<<"The eps is: " << eps << std::endl;
    new_row7(0) = Max_val;   
    record_eps = arma::join_cols(record_eps,new_row7);
    Rcout << std::fixed<<std::setprecision(10) <<"The record_eps is: " << record_eps << std::endl;
  }
  
  //second section
  arma::vec Record_q;
  arma::vec Record_dif_thetaf;
  arma::vec Record_step_theta;
  arma::vec Record_step_b;
  arma::vec Record_eps;  
  arma::mat ini = Ut.t()*Vt;
  Rcout << std::fixed<< std::setprecision(10) <<"The ini is: " << ini << std::endl;  
  
  arma::vec diag_elements = arma::diagvec(ini);    // Extract the diagonal elements as a vector
  arma::vec sign_diag_elements = arma::sign(diag_elements);  
  arma::mat sign_diagonal_matrix = arma::diagmat(sign_diag_elements); 
  lambda = arma::diagmat(sign_diagonal_matrix);
  Rcout << std::fixed<< std::setprecision(10) <<"The lambda is: " << lambda << std::endl;
  U0 = (Ut+Vt*lambda)/2;
  //  U0 = Ut;
  Rcout << std::fixed<< std::setprecision(10) <<"The U0k is: " << U0 << std::endl;
  
  b0 = bt;
  
  iter = 0;
  arma::mat Uk = U0;
  
  List bk = bt;
  double bkn = bk.size();
  List Record_dif_b(bkn); 
  
  for (int i = 0; i < bkn; ++i) {
    Record_dif_b[i] = arma::vec();
  }
  
  while (Eps > tol && iter < maxit2) {
    // double qfunction0 = lhs;
    List samples = sampleall(A, Xt, t, M, Uk, lambda, Vt, bk, sgammat,  N, d, 1);
    double qfunction = Q(A, Xt,t, Uk,lambda, Vt, bk, sgammat, samples, M, gam, d, 1);
    Rcout <<"The result of Q is: " << qfunction << std::endl;
    
    double Eps1 = arma::norm(Uk*lambda*Uk.t()-U0*lambda*U0.t(),"fro");
    Rcout << std::fixed<< std::setprecision(10) <<"The Eps1 is: " << Eps1 << std::endl;
    
    int bkn = bk.size();
    arma::vec Epsa(bkn);
    for (int i = 0; i < bkn; ++i){
      arma::mat bki = bk[i];
      arma::mat b0i = b0[i];
      double Epsi = arma::norm(bki-b0i,"fro");
      Epsa(i) = Epsi;
//      Rcout << std::fixed << std::setprecision(5) << "The eps" << (i + 2) << " is: " << Epsa(i) << std::endl;
    }
    
    double MAX_val = Eps1;
    for(arma::uword i = 0; i < Epsa.size(); i++) {
      MAX_val = std::max(MAX_val, Epsa(i));
    }    
    Rcout << std::setprecision(10) <<"The Eps at beginning of each EM is: " << MAX_val << std::endl;
    double approxq = 10000;
    Rcout  <<"The approxq at beginning of each EM is: " << approxq << std::endl;
    
    // arma::vec X0 = vector_minus_mean(X[0]);
    // arma::vec X1 = vector_minus_mean(X[1]);
    //iter1 and ite2 didn't get resey
    int iter1 = 0;
    Rcout  <<"The iter1 is: " << iter1 << std::endl;
    
    while (std::abs(lhs - approxq)/std::abs(approxq) > tol1 && iter1 < maxit1) {
      approxq = Q(A, Xt, t, Uk, lambda, Vt, bk, sgammat,  samples, M, gam, d,1);
      Rcout  << "The approxq at beginning of inner 2 is: " << approxq << std::endl;
      arma::mat id = arma::mat(d,d).fill(1);
      double denom = N * d * d;
      arma::mat hu = arma::zeros(d, r);
      arma::mat hv = arma::zeros(d, r);
      List hb(bkn);
      for (int i = 0; i < bkn; ++i){
        arma::mat hbi = arma::zeros(d, d);
        hb[i] = hbi;
      }      
      //      arma::mat hb1 = arma::zeros(d, d);
      //      arma::mat hb2 = arma::zeros(d, d);
      arma::mat hsg = arma::zeros(d, d);
      //      double N = 200;
      List hse(N);
      for (int q = 0; q < bkn; ++q) {
        arma::vec xq = X[q];
        arma::mat hbq_adjustment = arma::zeros(d, d); // Assuming xq is a vector
        
        for (int i = 0; i < N; ++i) {
          arma::cube samplesi = samples[i];
          arma::mat samplesi_mean = arma::mean(samplesi, 2);
//         for (int m = 0; m < M; ++m){
//            arma::mat samplesim = samplesi.slice(m);
          for (int T = 0; T < t(i); ++T) {
            arma::cube Yt = Y[T];
            arma::mat Yti = arma::mat(Yt.slice(i));
            arma::mat fm = Ut * trans(Vt) + samplesi_mean;
            double xqit;
            for(int inner_q = 0; inner_q < bkn; ++inner_q) {
              arma::vec inner_xq = X[inner_q];
              arma::vec xqi = inner_xq.subvec(i*t(i),i*t(i)+t(i)-1);
              xqit = xqi(T);
              arma::mat bq = as<arma::mat>(bt[inner_q]);
              fm += xqit * bq;
            }
            arma::vec outer_xqi = xq.subvec(i*t(i),i*t(i)+t(i)-1);
            double outer_xqit = outer_xqi(T);
            hu -= ((Yt.slice(i) - exp(fm)/(id + exp(fm))) * Vt) / (denom*bkn)  ;
            hv -= (trans(Yt.slice(i) - exp(fm)/(id + exp((fm)))) * Ut) / (denom*bkn);
            hbq_adjustment -= ((Yt.slice(i) - exp(fm)/(id + exp(fm)) ) * outer_xqit) / (denom) ;
          }
        }
        hb[q] = hbq_adjustment;
      }      
      // arma::mat Ug = hu.submat(0, 0, r-1, r-1);
      // arma::mat Vg = hv.submat(0, 0, r-1, r-1);
      // arma::mat hb1 = hb[0];
      // arma::mat hb3 = hb[2];
      // arma::mat hb5 = hb[4];
      // 
      // arma::mat b1g = hb1.submat(0, 0, 9, 9);
      // arma::mat b3g = hb3.submat(0, 0, 9, 9);
      // arma::mat b5g = hb5.submat(0, 0, 9, 9);
      // 
      // Rcout << "The initial value of sub-hu of inner 1 is: " << Ug << std::endl;
      // Rcout << "The initial value of sub-hv of inner 1 is: " << Vg << std::endl;
      // Rcout << "The initial value of sub-hb1 of inner 1 is: " << b1g << std::endl;
      // Rcout << "The initial value of sub-hb3 of inner 1 is: " << b3g << std::endl;
      // Rcout << "The initial value of sub-hb5 of inner 1 is: " << b5g << std::endl;
      
      U0 = Uk;
      V0 = Vt;
      b0 = bk;

      double a = step1;
      double c = step2;

      std::vector<arma::mat> dif1 = {-a*hu}; 
      std::vector<arma::mat> gradient1 = {hu}; 
      
      std::vector<arma::mat> dif2;
      List trunb(bkn);
      for(int i = 0; i < hb.size(); ++i) {
        arma::mat bki = bk[i];
        arma::mat hbi = hb[i];
        arma::mat mat = trun(bki-c*hbi,s,d);
        dif2.push_back(mat);
        trunb[i] = mat;
      }
      
      std::vector<arma::mat> gradient2;
      for(int i = 0; i < hb.size(); ++i) {
        arma::mat mat = hb[i];
        gradient2.push_back(mat);
      }
      
      arma::mat res2 =  trans(joinAndVectorize(dif1))*joinAndVectorize(dif1);
      double Res2 = res2(0);
      
      arma::mat res4 =  trans(joinAndVectorize(dif2))*joinAndVectorize(dif2);
      double Res4 = res4(0);
      
      while (Q(A, Xt, t, U0 - a * hu, lambda, V0 - a * hv, bk, sgammat, samples, M, gam, d, 0) - approxq -  (1/(2.0*a)*Res2) > 0) {
        
        arma::vec new_rowt(1);
        a *= 0.5;
        new_rowt(0) = a;  
        Record_step_theta = join_cols(Record_step_theta, new_rowt);
        Rcout << "The Record_step_theta of inner 2 is: " << Record_step_theta << std::endl;
      }
      Uk = U0 - a * hu;
      // arma::mat Uu = Uk.submat(0, 0, r-1, r-1);
      // Rcout << "The updated U of inner 1 is: " << Uu << std::endl;
      if (s>0){
      while (Q(A, Xt, t, Uk,lambda, Vt, trunb, sgammat, samples, M, gam, d, 0) - approxq - (1/(2.0*c)*Res4)  > 0){
        arma::vec new_rowb(1);
        c *= 0.8;
        new_rowb(0) = c;  
        Record_step_b = join_cols(Record_step_b, new_rowb);
//        Rcout << "The Record_step_b of inner 2 is: " << Record_step_b << std::endl;
      }
      for(int i = 0; i < hb.size(); ++i) {
        arma::mat bki = bk[i];
        arma::mat hbi = hb[i];
        arma::mat btruei = btrue[i];
        double nc = c * 2 * (arma::norm(bki - btruei, "fro"));
        arma::mat mat = trun(bki-nc*hbi,s,d);
        dif2.push_back(mat);
        trunb[i] = mat;
      }
      bk = trunb;
      }
      else{
          bk=b0;
        }
      // arma::mat b1 = bk[0];
      // arma::mat b3 = bk[2];
      // arma::mat b5 = bk[4];
      // arma::mat b1u = b1.submat(0, 0, 4, 4);
      // arma::mat b3u = b3.submat(0, 0, 4, 4);
      // arma::mat b5u = b5.submat(0, 0, 4, 4);
      // 
      // Rcout << "The updated b1 of inner 2 is: " << b1u << std::endl;
      // Rcout << "The updated b3 of inner 2 is: " << b3u << std::endl;
      // Rcout << "The updated b5 of inner 2 is: " << b5u << std::endl;

      lhs = Q(A, Xt,t, Uk, lambda, Vt, bk, sgammat,  samples, M, gam, d, 1);
      Rcout << std::setprecision(10) <<"The lhs at INNER 2 (test) is: " << lhs << std::endl;
      iter1++;
      double rdif = std::abs(lhs - approxq)/std::abs(approxq);
      Rcout << std::setprecision(10) <<"The abs(lhs - approxq)/abs(approxq) of inner 2 is: " << rdif  << std::endl;
      
      if(rdif < tol1){
        break;
      }
      
    }
    iter++;
    lhs = Q(A, Xt,t, Uk, lambda, Vt, bk, sgammat,  samples, M, gam, d, 1);
    Rcout << std::setprecision(10) <<"The lhs at end of each EM is: " << lhs << std::endl;
    
    approxq = 10000;
    Rcout  <<"The approxq at end of each EM is: " << approxq << std::endl;
    
    arma::vec new_row1(1);
    new_row1(0) = qfunction;  
    Record_q = arma::join_cols(Record_q, new_row1);
    Rcout << "The Record_q is: " << Record_q << std::endl;
    Rcout << std::fixed<< std::setprecision(10) <<"The Eps at the end is: " << Eps << std::endl;
    Rcout << "The iteration of outer loop is: " << iter << std::endl;
    
    arma::vec new_row2(1);
    new_row2(0) = arma::norm(Uk*lambda*Uk.t()-utrue*vtrue.t(), "fro");   
    Record_dif_thetaf = arma::join_cols(Record_dif_thetaf, new_row2);
    Rcout << std::fixed<< std::setprecision(10) <<"The Record_dif_thetaf is: " << Record_dif_thetaf << std::endl;
    
    
    for(int i = 0; i < bkn; ++i) {
      arma::mat bi = bk[i];
      arma::mat btruei = btrue[i];
      arma::vec new_row88(1);
      new_row88(0) = arma::norm(bi - btruei, "fro");
      // Calculate the Frobenius norm of the difference between bi and btruei
      arma::vec RDB = Record_dif_b[i];
      Record_dif_b[i] = arma::join_cols(RDB, new_row88);
      
      arma::vec temp = Record_dif_b[i];
      // Output the difference for this iteration
//      Rcout << std::fixed << std::setprecision(5) << "The Record_dif_b" << (i + 1) << " is: " << temp << std::endl;
    }
    
    Eps1 = arma::norm(Ut*Vt.t()-U0*V0.t(),"fro");
    Rcout << std::fixed<< std::setprecision(10) <<"The Eps1 at the end of each EM is: " << Eps1 << std::endl;
    
    for (int i = 0; i < bkn; ++i){
      arma::mat bti = bk[i];
      arma::mat b0i = b0[i];
      double epsi = arma::norm(bti-b0i,"fro");
      Epsa(i) = epsi;
      Rcout << std::fixed << std::setprecision(5) << "The Eps" << (i + 2) << " at the end of each EM is: " << Epsa(i) << std::endl;
    }
    
    double Max_VAL = Eps1;
    for(arma::uword i = 0; i < Epsa.size(); i++) {
      Max_VAL = std::max(Max_VAL, Epsa(i));
    }
    Rcout << std::setprecision(10) <<"The eps at the end of each EM is: " <<  Max_VAL << std::endl;
    arma::vec new_row7(1);
    Rcout << std::fixed<< std::setprecision(10)<<"The eps is: " << Max_VAL << std::endl;
    new_row7(0) = Max_VAL;   
    record_eps = arma::join_cols(record_eps,new_row7);
    Rcout << std::fixed<<std::setprecision(10) <<"The Record_eps is: " << record_eps << std::endl;
  }
  double ebic;
  double ebica;
  if (s > 0){
    ebic = 2 * N *d *d * lhs + (2 * d * r + d * d * s * p) * (log(d * d * N) + log(d * d * (p+1)));
    ebica = 2 * N * lhs + (2 * d * r + d * d * s * p) * (log(d * d * N) + log(d * d * (p+1)));
  }
  else{
    ebic = 2 * N *d *d * lhs + (2 * d * r) * (log(d * d * N) + log(d * d));
    ebica = 2 * N * lhs + (2 * d * r ) * (log(d * d * N) + log(d * d));
  }
  // double ebic = 2 * N * lhs + (2 * d * r + s) * (log(d * d * N) + log(d * d * (p+1)));
  List z = Rcpp::List::create(Ut,Vt,bt, record_q, record_dif_thetaf, record_dif_b, record_step_theta,record_step_b,Uk,bk, Record_q, Record_dif_thetaf, Record_dif_b, Record_step_theta,Record_step_b,ebic,lambda,ebica);
  return(z);
}

