//include Michael's code in the package
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix prGpostFunc(NumericMatrix phi_index, NumericMatrix lambda_index,
                          NumericMatrix phi, NumericMatrix lambda, NumericMatrix omega,
                          NumericVector pii, int FF, int SS, NumericVector cn_i) {

    int p = phi_index.ncol();
    int n = lambda_index.nrow(), q = lambda_index.ncol();
    int houseindex = 0;
    NumericMatrix prGpost(n, FF);

    for (int i = 0; i < n; i++) {
        NumericVector prGposti(FF);
        double updatesum = 0;
        int houseindexend = houseindex + cn_i[i];

        for (int g = 0; g < FF; g++) {
            double prodhouse = 1;
            for (int j = houseindex; j < houseindexend; j++) {
                double sumomegaprodphi = 0;
                for (int m = 0; m < SS; m++) {
                    double phiprod = 1;
                    for (int k = 0; k < p; k++){
                        int phiindexjk = phi_index(j, k);
                        phiprod *= phi(phiindexjk-1, g+(m*FF));
                    }
                    sumomegaprodphi += phiprod*omega(g, m);
                }
                prodhouse *= sumomegaprodphi;
            }

            double prodlambda = 1;
            for (int kq = 0; kq < q; kq++) {
                int lambdaindexjk = lambda_index(i, kq);
                prodlambda *= lambda(lambdaindexjk-1, g);
            }

            prGposti[g] = pii[g]*prodlambda*prodhouse;
            updatesum += prGposti[g];
        }

        for (int g = 0; g < FF; g++) {
            prGpost(i, g) = prGposti[g]/updatesum;
        }

        houseindex += cn_i[i];
    }
    return prGpost;
}

// [[Rcpp::export]]
NumericMatrix prMpostFunc(NumericMatrix phi_index, NumericMatrix phi,
                          NumericMatrix omega, NumericVector rep_G,
                          int FF, int SS) {

    int p = phi_index.ncol();
    int N =  phi_index.nrow();
    NumericMatrix prMpost(N, SS);

    for (int i = 0; i < N; i++) {
        int Gi = rep_G[i] - 1;
        NumericVector prMposti(SS);
        double updatesum = 0;

        for (int m = 0; m < SS; m++) {
            double phiprod = 1;
            for (int k = 0; k < p; k++) {
                int phiindexjk = phi_index(i, k);
                phiprod *= phi(phiindexjk-1, Gi+(m*FF));
            }
            prMposti[m] = omega(Gi, m)*phiprod;
            updatesum += prMposti[m];
        }

        for (int m = 0; m < SS; m++) {
            prMpost(i, m) = prMposti[m]/updatesum;
        }

    }
    return prMpost;
}

// [[Rcpp::export]]
NumericVector prEachComb02Func(NumericMatrix phi_index_02, NumericMatrix lambda_index_02,
                               NumericMatrix phi, NumericMatrix lambda, NumericMatrix omega,
                               NumericVector pii, int FF, int SS, NumericVector cn_i_02) {

    int p = phi_index_02.ncol();
    int n = phi_index_02.nrow()/2;
    int n_lambda_index = lambda_index_02.nrow();
    int q = lambda_index_02.ncol();
    int houseindex = 0;
    NumericVector prEachComb02(n);

    for (int i = 0; i < n; i++) {
        NumericVector prEachCombi(FF);
        double updatesum = 0;
        int houseindexend = houseindex + cn_i_02[i];

        for (int g = 0; g < FF; g++) {
            double prodhouse = 1;
            for (int j = houseindex; j < houseindexend; j++) {
                double sumomegaprodphi = 0;
                for (int m = 0; m < SS; m++) {
                    double phiprod = 1;
                    for (int k = 0; k < p; k++){
                        int phiindexjk = phi_index_02(j, k);
                        phiprod *= phi(phiindexjk-1, g+(m*FF));
                    }
                    sumomegaprodphi += phiprod*omega(g, m);
                }
                prodhouse *= sumomegaprodphi;
            }

            double sumlambda = 0;
            for (int kqq = 0; kqq < n_lambda_index; kqq++) {
                double prodlambda = 1;
                for (int kq = 0; kq < q; kq++) {
                    int lambdaindexjk = lambda_index_02(kqq, kq);
                    prodlambda *= lambda(lambdaindexjk-1, g);
                }
                sumlambda += prodlambda;
            }

            prEachCombi[g] = pii[g]*sumlambda*prodhouse;
            updatesum += prEachCombi[g];
        }

        prEachComb02[i] = updatesum;
        houseindex += cn_i_02[i];
    }
    return prEachComb02;
}

// [[Rcpp::export]]
NumericMatrix prGpost02Func(NumericMatrix phi_index_02, NumericMatrix lambda_index_02,
                            NumericMatrix phi, NumericMatrix lambda, NumericMatrix omega,
                            NumericVector pii, int FF, int SS, NumericVector cn_i_02) {

    int p = phi_index_02.ncol();
    int n = phi_index_02.nrow()/2;
    int n_lambda_index = lambda_index_02.nrow();
    int q = lambda_index_02.ncol();
    int houseindex = 0;
    NumericMatrix prGpost02(n, FF);

    for (int i = 0; i < n; i++) {
        NumericVector prGposti02(FF);
        double updatesum = 0;
        int houseindexend = houseindex + cn_i_02[i];

        for (int g = 0; g < FF; g++) {
            double prodhouse = 1;
            for (int j = houseindex; j < houseindexend; j++) {
                double sumomegaprodphi = 0;
                for (int m = 0; m < SS; m++) {
                    double phiprod = 1;
                    for (int k = 0; k < p; k++){
                        int phiindexjk = phi_index_02(j, k);
                        phiprod *= phi(phiindexjk-1, g+(m*FF));
                    }
                    sumomegaprodphi += phiprod*omega(g, m);
                }
                prodhouse *= sumomegaprodphi;
            }

            double sumlambda = 0;
            for (int kqq = 0; kqq < n_lambda_index; kqq++) {
                double prodlambda = 1;
                for (int kq = 0; kq < q; kq++) {
                    int lambdaindexjk = lambda_index_02(kqq, kq);
                    prodlambda *= lambda(lambdaindexjk-1, g);
                }
                sumlambda += prodlambda;
            }

            prGposti02[g] = pii[g]*sumlambda*prodhouse;
            updatesum += prGposti02[g];
        }

        for (int g = 0; g < FF; g++) {
            prGpost02(i, g) = prGposti02[g]/updatesum;
        }
        houseindex += cn_i_02[i];
    }
    return prGpost02;
}

// [[Rcpp::export]]
NumericMatrix sampleXMiss02Func(NumericMatrix lambda_index_house_miss_02,
                                NumericMatrix phi_index_struc_zero_pass_02, NumericMatrix lambda, NumericVector pii,
                                NumericMatrix phi, NumericMatrix omega, NumericVector G_X_miss_02,
                                NumericVector M_X_miss_02, int FF) {

    int p = phi_index_struc_zero_pass_02.ncol();
    int q = lambda_index_house_miss_02.ncol();
    int n = G_X_miss_02.size();
    int n_options = phi_index_struc_zero_pass_02.nrow()/2;
    NumericMatrix sampleXMiss02(n, n_options);

    for (int i = 0; i < n; i++) {
        NumericVector sampleXMissi02(n_options);
        double updatesum = 0;
        int g = G_X_miss_02[i];
        int m1 = M_X_miss_02[2*i];
        int m2 = M_X_miss_02[(2*i)+1];

        for(int ii = 0; ii < n_options; ii++){
            int option1 = 2*ii;
            int option2 = (2*ii)+1;

            double phiprod1 = 1;
            for (int k = 0; k < p; k++){
                int phiindexjk1 = phi_index_struc_zero_pass_02(option1, k);
                phiprod1 *= phi(phiindexjk1-1, g+(m1*FF));
            }

            double phiprod2 = 1;
            for (int kk = 0; kk < p; kk++){
                int phiindexjk2 = phi_index_struc_zero_pass_02(option2, kk);
                phiprod2 *= phi(phiindexjk2-1, g+(m2*FF));
            }

            double prodlambda = 1;
            for (int kq = 0; kq < q; kq++) {
                int lambdaindexjk = lambda_index_house_miss_02(i, kq);
                prodlambda *= lambda(lambdaindexjk-1, g);
            }

            sampleXMissi02[ii] = pii[g]*prodlambda*phiprod1*omega(g, m1)*phiprod2*omega(g, m2);
            updatesum += sampleXMissi02[ii];
        }

        for (int ii = 0; ii < n_options; ii++) {
            sampleXMiss02(i,ii) = sampleXMissi02[ii]/updatesum;
        }

    }
    return sampleXMiss02;
}
