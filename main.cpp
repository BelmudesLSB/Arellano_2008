#include <iostream>
#include <mex.h>
#include "ArellanoEconomy.hpp"
#include "omp.h"
#include <cstdint>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // Check for the proper number of arguments
    if (nrhs != 1 || nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "One input required.");
    }

    // Read the input parameters from MATLAB
    const mxArray* paramsStruct = prhs[0];
    int nthreads = static_cast<int>(mxGetScalar(mxGetField(paramsStruct, 0, "nthreads")));
    double beta = mxGetScalar(mxGetField(paramsStruct, 0, "beta"));
    int ny = static_cast<int>(mxGetScalar(mxGetField(paramsStruct, 0, "ny")));
    int nb = static_cast<int>(mxGetScalar(mxGetField(paramsStruct, 0, "nb")));
    double b_min = mxGetScalar(mxGetField(paramsStruct, 0, "b_min"));
    double b_max = mxGetScalar(mxGetField(paramsStruct, 0, "b_max"));
    double r = mxGetScalar(mxGetField(paramsStruct, 0, "r"));
    double rho = mxGetScalar(mxGetField(paramsStruct, 0, "rho"));
    double sigma = mxGetScalar(mxGetField(paramsStruct, 0, "sigma"));
    double theta = mxGetScalar(mxGetField(paramsStruct, 0, "theta"));
    double m = mxGetScalar(mxGetField(paramsStruct, 0, "m"));
    double y_bar = mxGetScalar(mxGetField(paramsStruct, 0, "y_bar"));
    int max_iter = static_cast<int>(mxGetScalar(mxGetField(paramsStruct, 0, "max_iter")));
    double tol = mxGetScalar(mxGetField(paramsStruct, 0, "tol"));

    //New! set the number of threads:
    omp_set_num_threads(nthreads);

    // Create vectors to store the output
    mxArray* ygrid = mxCreateDoubleMatrix(ny, 1, mxREAL);
    mxArray* bgrid = mxCreateDoubleMatrix(nb, 1, mxREAL);
    mxArray* p = mxCreateDoubleMatrix(ny * ny, 1, mxREAL);
    mxArray* v = mxCreateDoubleMatrix(ny * nb, 1, mxREAL);
    mxArray* v_r = mxCreateDoubleMatrix(ny * nb, 1, mxREAL);
    mxArray* v_d = mxCreateDoubleMatrix(ny, 1, mxREAL);
    mxArray* q = mxCreateDoubleMatrix(ny * nb, 1, mxREAL);
    mxArray* b_p = mxCreateNumericMatrix(ny * nb, 1, mxINT32_CLASS, mxREAL);
    mxArray* d_p = mxCreateDoubleMatrix(ny * nb, 1, mxREAL);

    // Create pointers to the output vectors
    double* ygridPtr = mxGetPr(ygrid);
    double* bgridPtr = mxGetPr(bgrid);
    double* pPtr = mxGetPr(p);
    double* vPtr = mxGetPr(v);
    double* v_rPtr = mxGetPr(v_r);
    double* v_dPtr = mxGetPr(v_d);
    double* qPtr = mxGetPr(q);
    int32_t* b_pPtr = reinterpret_cast<int32_t*>(mxGetData(b_p));
    double* d_pPtr =  mxGetPr(d_p);


    // Create an instance of the ArellanoEconomy class
    ArellanoEconomy arellanoeconomy(beta, ny, nb, b_min, b_max, r, rho, sigma, theta, m, y_bar, ygridPtr, bgridPtr, pPtr, vPtr, v_rPtr, v_dPtr, qPtr, d_pPtr, b_pPtr, max_iter, tol);

    // Create grids
    arellanoeconomy.initialize();

    // Solve the model
    arellanoeconomy.solve();    
    

    // Set the output vectors in the MATLAB output
    const char* fieldNames[9] = {"Ygrid", "Bgrid", "P", "V", "V_r", "V_d", "Q", "B_p", "D_p"};
    plhs[0] = mxCreateStructMatrix(1, 1, 9, fieldNames);
    mxSetField(plhs[0], 0, "Ygrid", ygrid);
    mxSetField(plhs[0], 0, "Bgrid", bgrid);
    mxSetField(plhs[0], 0, "P", p);
    mxSetField(plhs[0], 0, "V", v);
    mxSetField(plhs[0], 0, "V_r", v_r);
    mxSetField(plhs[0], 0, "V_d", v_d);
    mxSetField(plhs[0], 0, "Q", q);
    mxSetField(plhs[0], 0, "B_p", b_p);
    mxSetField(plhs[0], 0, "D_p", d_p);
}
