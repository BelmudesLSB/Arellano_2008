#include "ArellanoEconomy.hpp"
#include "Utils.hpp"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <omp.h>

using namespace std;

ArellanoEconomy::ArellanoEconomy(double beta, int ny, int nb, double b_min, double b_max, double r, double rho, double sigma, double theta, double m, double y_bar, double* ygrid, double* bgrid, double* p, double* v, double* v_r, double* v_d, double* q, double* d_p, int* b_p, int max_iter, double tol) {
    // Set the parameters:
    Beta = beta;
    Ny = ny;
    Nb = nb;
    Bmin = b_min;
    Bmax = b_max; 
    R = r;
    Rho = rho;
    Sigma = sigma;
    Theta = theta;
    M = m;
    Y_bar = y_bar;
    Max_iter = max_iter;
    Tol = tol;

    // Set the pointers:
    Ygrid = ygrid;
    Bgrid = bgrid;
    P = p;
    V = v;
    V_r = v_r;
    V_d = v_d;
    Q = q;
    D_p = d_p;
    B_p = b_p;
}

// Create the bond grid:
int ArellanoEconomy::bond_grid(){  
    double bstep = (Bmax - Bmin)/(Nb - 1);
    Bzero_idx = static_cast<int>(std::round(-Bmin / bstep));
    #pragma omp parallel 
    {
        #pragma omp parallel for
        for(int i = 0; i < Nb; i++){
            Bgrid[i] = Bmin + i*bstep;
        }
    }


    return EXIT_SUCCESS;
}

// Create the income and probabilty grid using Tauchen (1986):
 int ArellanoEconomy::income_prob_grid(){

    double sigma_y = sqrt(pow(Sigma,2)/(1-pow(Rho,2)));
    double omega = (2*M*sigma_y)/(Ny-1);

    #pragma omp parallel
    {      
        #pragma omp parallel for
        for (int i=0; i<Ny; i++){ 
            Ygrid[i] = (-M*sigma_y)  + omega * i;
        }

        #pragma omp barrier

        // Create the transition vector:
        //This function will create a vector of size N*N. The first N elements will be the first row of the matrix,
        // the next N elements will be the second row of the matrix, and so on. Element [i * N + j] will represent the
        // probability that from state (i) will move to the state (j).
      
        int j;
        #pragma omp for
        for (int i=0; i<Ny; i++){
            for (int j=0; j<Ny; j++){
                // We need to treat endpoints separately:
                if (j==0 || j==Ny-1){
                    if (j==0){
                        P[i*Ny+j] = normalCDF((Ygrid[0]-Rho*Ygrid[i]+omega/2)/Sigma);
                    }
                    else {
                        P[i*Ny+j] = 1-normalCDF((Ygrid[Ny-1]-Rho*Ygrid[i]-omega/2)/Sigma);
                    }
                // Probability distribution inside the y_grid:
                } else {
                    P[i*Ny+j] = normalCDF((Ygrid[j]-Rho*Ygrid[i]+omega/2)/Sigma)-normalCDF((Ygrid[j]-Rho*Ygrid[i]-omega/2)/Sigma);
                }
            }
        }

        #pragma omp barrier
        
        #pragma omp for
        for (int i=0; i<Ny; i++){
            Ygrid[i] = exp(Ygrid[i]);
        }

    }

    return EXIT_SUCCESS;
 }

// Intitialize the value function at default, using value function iteration:
int ArellanoEconomy::value_default(){

    int iter = 0;
    double* Vd_0 = new double[Ny];
    double dV_d = 0;
    bool Convergence_achieved = false;

    #pragma omp parallel for
    for (int i=0; i<Ny; i++){
            Vd_0[i] = (1/(1-Beta))*utility(h(Ygrid[i],Y_bar));
    }

    while (iter < Max_iter && Convergence_achieved == false){

        #pragma omp parallel 
        {   
            int id_thread = omp_get_thread_num(); 
            double temp2, temp3;
            int j;

            #pragma omp for
            for (int i=0; i<Ny; i++){
                temp2 = 0;
                for (int j=0; j<Ny; j++){
                    temp2 += Beta*P[i*Ny+j]*Vd_0[j];
                }
                V_d[i] = utility(h(Ygrid[i],Y_bar)) + temp2; 
            }

            #pragma omp barrier

            temp3 = 0;
            #pragma omp for reduction(max: dV_d)
            for (int i=0; i<Ny; i++){
                temp3 = abs(V_d[i]-Vd_0[i]);
                if (temp3 > dV_d){
                    dV_d = temp3;
                }
            }
            
            if (dV_d < Tol){
                // Delete dynamic memory:
                if (id_thread == 0){
                    delete[] Vd_0;
                    Convergence_achieved = true;
                    }
            } else {
                #pragma omp for
                for (int i=0; i<Ny; i++){
                    Vd_0[i] = V_d[i];
                }
                if (id_thread == 0){
                    dV_d = 0;
                    iter += 1;
                }
            }

        }

    }

    if (Convergence_achieved){
        return EXIT_SUCCESS;
    } else {
        std::cout<<"Value Function at default failure.";
        return EXIT_FAILURE;
    }
        
}

// Run all previous functions, all functions will check whether the results obtain are correct.
void ArellanoEconomy::initialize(){
    if (income_prob_grid() == 0 && bond_grid() == 0 && value_default()==0){
        std::cout << std::endl;
        std::cout << "..................................................................................." << std::endl;
        std::cout << "Initialization Completed Successfully \n";
        std::cout << std::endl;
    } else {
        std::cout << "ERROR: Initialization!!! \n";
    }
}

// Solve the model:
int ArellanoEconomy::solve(){

    // Step 0. 

    // Allocate memory for old (V_0, V_r0, Q0):
    double* V_0 = new double[Ny*Nb];
    double* V_r0 = new double[Ny*Nb];
    double* Q0 = new double[Ny*Nb];
    double* EV_0 = new double[Ny*Nb];
    // Allocate memory for new (V_1, V_r1, Q1, D_p1, B_p1):
    double* V_1 = new double[Ny*Nb];
    double* V_r1 = new double[Ny*Nb];
    double* Q1 = new double[Ny*Nb];
    double* D_p1 = new double[Ny*Nb];
    int* B_p1 = new int[Ny*Nb];
    

    // Step 1. Initialize (V_0, V_r0, Q0) with a simple guess:
    for (int i=0; i<Ny; i++){
        for (int j=0; j<Nb; j++){
            V_r0[i*Nb+j] = 1/(1-Beta)*utility(Ygrid[i]);
            V_0[i*Nb+j] = V_r0[i*Nb+j];
            Q0[i*Nb+j] = 1/R;
        }
    }

    int iter = 0;
    bool Convergence_achieved = false;
       
    while (iter < Max_iter && !Convergence_achieved){

        double dV = 0;
        double dq = 0;
        
        #pragma omp parallel  
        {
            // Step 1.5. Create a matrix for the expected value of the value function:
            
            int id_thread = omp_get_thread_num();
            int num_threads = omp_get_num_threads();  
            double temp1; // Each thread has their own copy of temp1.
            int j, k;     // Each thread has their own copy of j and k.
            
            #pragma omp for
            for (int i=0; i<Ny; i++){
                for (int j=0; j<Nb; j++){
                    temp1 = 0;
                    for (int k=0; k<Ny; k++){
                        temp1 += P[i*Ny+k]*V_0[k*Nb+j];
                    }
                    EV_0[i*Nb+j] = temp1;
                }
            }
            

            j=0; 
            // Step 2. Update V_r and B_p, given Q0 and EV_0:
            double temp2; // Each thread has their own copy of temp2.
            double temp3; // Each thread has their own copy of temp3.
            double temp4; // Each thread has their own copy of temp4.
            int temp5;    // Each thread has their own copy of temp5.
            #pragma omp for schedule(dynamic)
            for (int i=0; i<Ny; i++){
                for (int j=0; j<Nb; j++){
                    temp4 = -1000000;
                    temp5 = 0;
                    for (int k=0; k<Nb; k++){
                        temp2 = Ygrid[i] + Q0[i*Nb+k]*Bgrid[k] - Bgrid[j];
                        temp3 = Beta*EV_0[i*Nb+k];
                        if (temp2 > 0){
                            if (utility(temp2) + temp3 >= temp4){
                                temp4 = utility(temp2) + temp3;
                                temp5 = k;
                            }
                        }
                    }
                    V_r1[i*Nb+j] = temp4;
                    B_p1[i*Nb+j] = temp5;
                }
            }
  
            // Step 3. Update V and D_p:
            double temp6;
            double temp7;
            int temp11 = Bzero_idx;
            j = 0;
            #pragma omp for
            for (int i=0; i<Ny; i++){
                temp7 = V_d[i];
                for (int j=0; j<Nb; j++){
                    temp6 = V_r1[i*Nb+j];
                    if (temp6 >= temp7){
                        V_1[i*Nb+j] = temp6;
                        D_p1[i*Nb+j] = 0.00;
                    } else {
                        V_1[i*Nb+j] = temp7;
                        D_p1[i*Nb+j] = 1.00;
                        B_p1[i*Nb+j] = temp11;
                    }
                }
            }


            // Step 4. Update Q, as the expected probality of repayment:
            double temp8;
            j = 0;
            k = 0;
            #pragma omp for
            for (int i=0; i<Ny; i++){
                for (int j=0; j<Nb; j++){
                    temp8 = 0;
                    for (int k=0; k<Ny; k++){
                        temp8 += (1/R)*(P[i*Ny+k]*(1-D_p1[k*Nb+j]));
                    }
                    Q1[i*Nb+j] = temp8;
                }
            }
            
            // Step 5. Check for convergence:

            double temp9;   // Each thread has their own copy of temp9.
            double temp10;  // Each thread has their own copy of temp10.
            j=0;

            #pragma omp for reduction(max: dV, dq)
            for (int i=0; i<Ny; i++){
                for (int j=0; j<Nb; j++){
                    temp9 = abs(V_1[i*Nb+j]-V_0[i*Nb+j]);
                    temp10 = abs(Q1[i*Nb+j]-Q0[i*Nb+j]);
                    if (dV < temp9 || dq < temp10){
                        dV = temp9;
                        dq = temp10;
                    }
                }
            }

            j = 0;

            if (dV>Tol || dq>Tol){
                #pragma omp for
                for (int i=0; i<Ny; i++){
                    for (int j=0; j<Nb; j++){
                        V_0[i*Nb+j] = V_1[i*Nb+j];
                        V_r0[i*Nb+j] = V_r1[i*Nb+j];
                        Q0[i*Nb+j] = Q1[i*Nb+j];
                        D_p[i*Nb+j] = D_p1[i*Nb+j];
                        B_p[i*Nb+j] = B_p1[i*Nb+j];
                    }
                }
            } else {
                #pragma omp for
                for (int i=0; i<Ny; i++){
                    for (int j=0; j<Nb; j++){
                        V[i*Nb+j] = V_1[i*Nb+j];
                        V_r[i*Nb+j] = V_r1[i*Nb+j];
                        Q[i*Nb+j] = Q1[i*Nb+j];
                        D_p[i*Nb+j] = D_p1[i*Nb+j];
                        B_p[i*Nb+j] = B_p1[i*Nb+j];
                    }
                }

                if (id_thread == 0) {
                    std::cout << "Solution Found: dV = " << dV << " iter = " << iter << " N. Threads: "<< num_threads <<std::endl;
                    std::cout << std::endl;
                    delete [] V_0;
                    delete [] V_r0;
                    delete [] Q0;
                    delete [] EV_0;
                    delete [] V_1;
                    delete [] V_r1;
                    delete [] Q1;
                    delete [] D_p1;
                    delete [] B_p1;
                    Convergence_achieved = true;
                }
            }

        }// End of parallel region.

    iter += 1;
        
    }

    if (Convergence_achieved){

        return EXIT_SUCCESS;

    } else {

        std::cout << " iter = " << iter << std::endl;
        return EXIT_SUCCESS;

    }

} 

// Export v to a csv file:
void ArellanoEconomy::export_V(){
    // Create a file to store the results:
    FILE *f = fopen("results.csv", "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    // Print the results to the file:
    for (int i=0; i<Ny; i++){
        for (int j=0; j<Nb; j++){
            fprintf(f, "%f \n", V[i*Nb+j]);
        }
    }

    // Close the file:
    fclose(f);
}

