#include <Rcpp.h>
using namespace Rcpp;

// This is a c++ version of the the forward part of the explicit duration forward backward algorithm
// for a categorical emission distribution
//' @keywords internal
// [[Rcpp::export(rng = false)]]

 List mult_ed_fb_cpp(int m, int n, NumericVector delta, NumericMatrix allprobs, int Mx, IntegerVector Mx2,NumericMatrix gamma, NumericMatrix d, IntegerVector S, IntegerVector S2) {
     int j, t, i, u, uMax, v, k, Len;

     int zer = 0;
     NumericMatrix d2 = clone(d);

     double x;
     NumericMatrix D2(m,n);
     NumericVector dSum(m);

     NumericVector N(n);
     NumericMatrix Norm(m,n);
     NumericMatrix Forward(m,n);
     NumericMatrix StateIn(m,n);
     double Observ = 0;

     NumericMatrix Backward(m, n);
     NumericMatrix B_star(m, n+2);
     IntegerVector VarL(Mx-1);
     for (i = 1; i < Mx; i++) {
         VarL(i-1) = Mx - i;
     }
     int occNcol = Mx * (n-Mx+1) + sum(VarL);
     NumericMatrix Occupancy(m, occNcol);

     IntegerVector lengthID(n);
     for (i = 0; i < n; i++) {
         if (i < n-Mx+1) {
             lengthID(i) = Mx;
         }
         else {
             lengthID(i) = VarL(i - (n-Mx+1));
         }
     }

     IntegerVector endID(n + 2);
     for (i = 0; i < n+2; i++) {
         if (i == 0) {
             endID(i) = 0;
         }
         if (i > 0) {
             if (i < n+1) {
                 endID(i) =  endID(i-1) + lengthID(i-1);
             }
         }
         if (i == n+1) {
             endID(i) = endID(i-1);
         }
     }

     IntegerVector::const_iterator first = S.begin() + 0;

     // forward recursion
     for (t = 0; t <= n-1; t++) {

         uMax = std::min(t+1, Mx+1);

         IntegerVector::const_iterator last = S.begin() + (t+1);
         IntegerVector SSh(first, last);
         Len = SSh.size();
         IntegerVector SShrev(Len);
         for (i = 0; i < Len; i++) {
             SShrev(i) = SSh(Len-1-i);
         }

         IntegerVector::const_iterator first2 = SShrev.begin() + 0;
         IntegerVector::const_iterator last2 = SShrev.begin() + (uMax);
         IntegerVector SShrev2(first2, last2);

         d2 = clone(d);
         for (i = 1; i < uMax; i++) {
             for (j = 0; j < m; j++) {
                 d2(j,i) *= SShrev2(i-1);
             }
         }

         for (j = 0; j < m; j++) {
             dSum(j) = 0;
             for (i = 0; i <= Mx; i++) {
                 dSum(j) += d2(j,i);
             }
         }

         for (j = 0; j < m; j++) {
             if (dSum(j) != 0) {
                 for (i = 0; i <= Mx; i++) {
                     d2(j,i) /= dSum(j);
                 }
             }
         }

         if (t == n-1) {
             for (j = 0; j < m; j++) {
                 for (u = 1; u <= Mx; u++) {
                     x = 0;
                     for (v = u; v < Mx + 1; v++)
                         x += d2(j,v);
                     D2(j,(u-1)) = x;
                 }
                 for (u = Mx + 1; u <= n; u++) {
                     D2(j,(u-1)) = 0;
                 }
             }
         }

         N(t) = 0;
         for (j = 0; j < m; j++) {
             if (t == 0) {
                 Norm(j,0) = log(delta(j)) + log(allprobs(j,0));
             }
             else
             {
                 Norm(j,t) = log(allprobs(j,t)) + log(std::abs(exp(StateIn(j,t)) - exp(Forward(j, (t-1))) + exp(Norm(j, (t-1)))));
             }
             N(t) += exp(Norm(j,t));
         }
         N(t) = log(N(t));
         for (j = 0; j < m; j++) {
             Norm(j,t) -= N(t);
         }

         for (j = 0; j < m; j++) {
             Forward(j,t) = 0;
             Observ = 0;

             if (t < n-1) {
                 for (u = 1; u <= std::min(t+1, Mx2(j)); u++) {
                     Observ += log(allprobs(j,t-u+1)) - N(t-u+1);
                     if (SShrev2(u-1) == 1) {
                         if (u < t + 1) {
                             Forward(j,t) += exp(Observ + log(d2(j,u)) + StateIn(j, (t-u+1)));
                         }
                         else {
                             Forward(j,t) += exp(Observ + log(d2(j,t+1)) + log(delta(j)));
                         }
                     }
                 }
                 Forward(j,t) = log(Forward(j,t));
             }
             else {
                 for (u = 1; u <= std::min(n, Mx2(j)); u++) {
                     Observ += log(allprobs(j,t-u+1)) - N(t-u+1);
                     if (SShrev2(u-1) == 1) {
                         if (u < n) {
                             Forward(j,n-1) += exp(Observ + log(D2(j,u)) + StateIn(j, n-u));
                         }
                         else {
                             Forward(j,n-1) += exp(Observ + log(D2(j,n)) + log(delta(j)));
                         }
                     }
                 }
                 Forward(j,n-1) = log(Forward(j,n-1));
             }
         }
         if (t < n-1){
             for (j = 0; j < m; j++){
                 StateIn(j, t+1) = 0;
                 for (i = 0; i < m; i++){
                     StateIn(j, t+1) += exp(Forward(i,t) + log(gamma(i,j)));
                 }
                 StateIn(j,t+1) = log(StateIn(j, t+1));
             }
         }

     }

     // Backward recursion

     for (t = n-1; t >= 0; t--) {
         if (S(t) == 1) {

             uMax = std::min(n-t, Mx);
             IntegerVector::const_iterator first3 = S2.begin() + t;
             IntegerVector::const_iterator last3 = S2.begin() + (n);
             IntegerVector SShB(first3, last3);

             IntegerVector::const_iterator first4 = SShB.begin() + 0;
             IntegerVector::const_iterator last4 = SShB.begin() + (uMax);
             IntegerVector SShB2(first4, last4);

             d2 = clone(d);
             for (i = 1; i < uMax+1; i++) {
                 for (j = 0; j < m; j++) {
                     d2(j,i) *= SShB2(i-1);
                 }
             }

             for (j = 0; j < m; j++) {
                 dSum(j) = 0;
                 for (i = 0; i <= Mx; i++) {
                     dSum(j) += d2(j,i);
                 }
             }

             for (j = 0; j < m; j++) {
                 if (dSum(j) != 0) {
                     for (i = 0; i <= Mx; i++) {
                         d2(j,i) /= dSum(j);
                     }
                 }
             }

             for (j = 0; j < m; j++) {
                 for (u = 1; u <= Mx; u++) {
                     x = 0;
                     for (v = u; v < Mx + 1; v++)
                         x += d2(j,v);
                     D2(j,(u-1)) = x;
                 }
                 for (u = Mx + 1; u <= n; u++) {
                     D2(j,(u-1)) = 0;
                 }
             }

             for (j = 0; j < m; j++) {
                 B_star(j,t) = 0;
                 Observ = 0;
                 for (u = 1; u <= std::min(n - t, Mx2(j)); u++) {
                     Observ += log(allprobs(j,t+u-1)) - N(t+u-1);
                     if (SShB2(u-1) == 1) {
                         if (u < n-t) {
                             Occupancy(j, endID(t) + (u-1)) = exp(Backward(j,t+u) + Observ + log(d2(j,u)));
                         }
                         else {
                             Occupancy(j, endID(t) + (u-1)) = exp(Observ + log(D2(j,n-1-t)));
                         }
                         B_star(j,t) += Occupancy(j, endID(t) + (u-1));
                     }
                 }
                 B_star(j,t) = log(B_star(j,t));
             }
             for (j = 0; j < m; j++) {
                 Backward(j,t) = 0;
                 for (k = 0; k < m; k++) {
                     Backward(j,t) += exp(B_star(k,t) + log(gamma(j,k)));
                 }
                 Backward(j,t) = log(Backward(j,t));
             }
         }
     }

     List H(n);
     for (i = 0; i < n; i++) {
         if (S(i) == 0) {
             H(i) = zer;
         }
         else {
             NumericMatrix foo(m,lengthID(i));
             for (k = 0; k < lengthID(i); k++) {
                 foo(_,k) = Occupancy(_, k + endID(i));
             }
             H(i) = foo;
         }
     }

     return List::create(N, B_star, H);
 }

