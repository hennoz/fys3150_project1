//
//  main.cpp
//  FYS3150_Project1_Mortensen
//
//  Created by Henrik Mortensen on 27/08/2018.
//  Copyright © 2018 Henrik Mortensen. All rights reserved.
//
//  Project 1) Solve the one-dimensional Poisson equation with Dirichlet boundary conditions by rewriting it as a set of linear equations.

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include <time.h>
#include <ctime>

using namespace std;
using namespace arma;

//  Functions to be used in the program
inline double f(double x){return 100.0*exp(-10.0*x);}
inline double exact(double x){return 1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x);}
inline double erfunc(double x, double v){ return log10(fabs((v - exact(x))/exact(x)));}
void general(int n);
void special(int n);
void ludecomp(int n);
//  Begin main program
int main(int argc, const char * argv[]) {

    //  From commandline x => n = 10^N
    int N = atoi(argv[1]); // N is the exponent

    for (int i = 1; i <= N; i++){
        int n = (int) pow(10, i);

        general(n);
        special(n);
        ludecomp(n);


    }
    return 0;
}
// ---------------------------------------------------------------
// GENERAL
void general(int n){
    //  Step length or spacing defined as h = 1/(n+1)
    double h = 1.0/(n + 1);

    //  Create all arrays
    double *x = new double[n+2];
    double *a = new double[n+1];
    double *b = new double[n+2];
    double *c = new double[n+1];
    double *v = new double[n+2];
    double *b_tilde = new double[n+2];
    double *B = new double[n+2];
    double *B_tilde = new double[n+2];

    //  Grid points x_i = ih from x_0 = 0 to x_{n+1}=1
    for (int i = 0; i < n+2; i++){
        x[i] = i*h ;
    }
    //  Fill vectors a, b, c, and b_tilde with explicit values
    for (int i = 0; i < n+2; i++){
        b_tilde[i] = h*h*f(x[i]);
        b[i] = 2;

    }
    for (int i = 0; i < n+1; i++){
        a[i] = -1;
        c[i] = -1;
    }
    //  Initial conditions
    B_tilde[1] = b_tilde[1];
    B[1] = B[n + 1] = 2;
    v[0] = v[n + 1] = 0;

    double k;
//LOOPS GENERAL
    //  START TIME
    clock_t start,end;
    start = clock();
    //  FORWARD SUBSTITUTION
    //-------------------------------------------------------
    for (int i = 2; i < n + 1; i++){
        k = a[i-1]/B[i-1];
        B[i] = b[i] - k*c[i-1] ;
        B_tilde[i] = b_tilde[i] - k*B_tilde[i-1] ;
    }
    //  BACKWARD SUBSTITUTION
    //------------------------------------------------------
    v[n+1] = B_tilde[n+1] / B[n+1];
    for (int i = n; i > 0; i--){
        v[i] = (B_tilde[i] - c[i]*v[i+1]) / B[i];
        //cout << v[i] << endl;
    }
    //  END TIME
    end = clock();
    double t_g = (double(end - start)/CLOCKS_PER_SEC)*1000;
    cout << t_g << " milliseconds to run general case algorithm when n = " << n << endl;
    // WRITE, write data to file
    //-----------------------------------------------------

    ofstream outfile;
    outfile.open("v_data_file.txt");
    for (int i = 0; i < n+2; i++){
        outfile << setw(8) << x[i];
        outfile << setw(20) << setprecision(8) << exact(x[i]);
        outfile << setw(20) << setprecision(8) << v[i];
        outfile << setw(20) << setprecision(16) << erfunc(x[i], v[i]) << endl;
    }
    outfile.close();
    outfile.open("v_data_file" + to_string(n) + ".txt");
    for (int i = 1; i < n+2; i++){
        outfile << setw(20) << setprecision(16) << erfunc(x[i], v[i]) << endl;
    }
    outfile.close();
    delete[]x;delete[]a;
    delete[]b;delete[]c;
    delete[]v;delete[]b_tilde;
    delete[]B_tilde;

}
// ---------------------------------------------------------------
// SPECIAL
void special(int n){

    double h = 1.0/(n + 1);

    //  Create all arrays
    double *    b = new double[n+2];
    double *b_tilde = new double[n+2];
    double *v = new double[n+2];
    double *x = new double[n+2];

    //  Initial conditions
    b[1] = b[n] = 2;
    v[0] = v[n] = 0;

    for (int i = 1; i < n + 1; i++) {
        b[i] = (i + 1.0)/( (double) i);
        }
    for (int i = 0; i < n + 2; i++){
        x[i]= i*h;
        b_tilde[i] = h*h*f(x[i]);
    }
//LOOPS SPECIAL
    //  START TIME
    clock_t start,end;
    start = clock();
    // Forward substitution
    for (int i = 2; i < n+1; i++){
        b_tilde[i] += b_tilde[i-1]/b[i-1];
    }
    // Backward substitution
    for (int i = n; i > 0; i--){
        v[i] = (b_tilde[i] + v[i+1])/b[i];
        //cout << v[i] << endl;
    }
    //  END TIME
    end = clock();
    double t_s = (double(end - start)/CLOCKS_PER_SEC)*1000;
    cout << t_s << " milliseconds to run special case algorithm when n = " << n << endl;

    // WRITE, write data to file
    //-----------------------------------------------------
    ofstream outfile;
    outfile.open("v_data_file_s.txt");
    for (int i = 0; i < n+2; i++){
        outfile << setw(8) << x[i];
        outfile << setw(20) << setprecision(8) << exact(x[i]);
        outfile << setw(20) << setprecision(8) << v[i];
        outfile << setw(20) << setprecision(8) << erfunc(x[i], v[i]) << endl;
    }

    delete[]b;delete[]b_tilde;delete[]v;delete[]x;
}
// ---------------------------------------------------------------
// LU DECOMPSITION
void ludecomp(int n){
    mat A(n,n);
    for ( int i = 0; i < n; i++){
        for ( int j = 0; j < n; j++){
            if ( i == j){A(i,j) = 2;}
            else if ( i == j - 1){A(i,j) = - 1;}
            else if ( i == j + 1){A(i,j) = - 1;}
        }
    }
    vec b_tilde(n);
    double *x = new double[n+1];
    double h = 1.0/(n + 1);
    for (int i = 1; i < n + 1; i++){
        x[i]= i*h;
        b_tilde[i-1] = h*h*f(x[i]);
    }
    vec w;
    vec y;
    //  START TIME
    clock_t start,end;
    start = clock();
    mat L, U;
    lu(L,U,A);
    w = solve(L,b_tilde);
    y = solve(U,w);
    end = clock();
    double t_lu = (double(end - start)/CLOCKS_PER_SEC)*1000;
    cout << t_lu << " milliseconds to run LU-decomp w. Armadillo when n = " << n << endl;
}


