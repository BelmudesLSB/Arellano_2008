#ifndef ArellanoEconomy_H
#define ArellanoEconomy_H
#include <cmath>
#include <iostream>
#include <mex.h>
#include "Utils.hpp"
#include "omp.h"


class ArellanoEconomy {

public:
    // parameters for this economy:
    double Beta;    // discount factor.
    int Ny;         // number of income states.
    int Nb;         // number of bond states.
    double Bmin;    // minimum bond holdings.
    double Bmax;    // maximum bond holdings.
    int Bzero_idx;  // index of zero bond holdings.
    double R;       // gross interest rate.
    double Rho;     // persistence of income.
    double Sigma;   // standard deviation of income shocks.
    double Theta;   // probability of re entering financial markets.
    double M;       // number of standard deviations for Tauchen (1986).
    double Y_bar;   // Maximum income at default.
    int Max_iter;   // Maximum number of iterations.
    double Tol;     // Tolerance for convergence.

    // pointers to operate on:
    double* Ygrid;  // pointer for grid of income.
    double* Bgrid;  // pointer for grid of bonds.
    double* P;      // pointer for transition matrix.
    double* V;      // pointer for value function.
    double* V_r;    // pointer for value function at repayment.
    double* V_d;    // pointer for value function at default.
    double* Q;      // pointer for price function.
    int* B_p;    // pointer for bond policy function.
    double* D_p;    // pointer for default policy function.
    
    // Construct Economy:
    ArellanoEconomy(double beta, int ny, int nb, double b_min, double b_max, double r, double rho, double sigma, double theta, double m, double y_bar, double* ygrid, double* bgrid, double* p, double* v, double* v_r, double* v_d, double* q, double* d_p, int* b_p, int max_iter, double tol);

    // Construct Economy  bond, income and probability grids:
    int bond_grid();
    int income_prob_grid();
    int value_default();
    void initialize();
    int solve();
    void export_V();
};

#endif
