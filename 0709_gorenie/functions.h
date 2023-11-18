#pragma once
#include <iostream>
#include <fstream>
#include "constants.h"
#include "concentration.h"
#include "thermodynamic.h"
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <cmath>
#include "chemkinReader.h"
#include "boost/regex.hpp"
#include "boost/format.hpp"
using namespace std;

#include <ida/ida.h>   
#include <kinsol/kinsol.h>             /* access to KINSOL func., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver  */
#include <cvode/cvode.h>   

#define FTOL   RCONST(1.e-15) /* function tolerance */
#define STOL   RCONST(1.e-15) /* step tolerance     */

#define ZERO   RCONST(0.0)
#define PT25   RCONST(0.25)
#define PT5    RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)
#define TWO    RCONST(2.0)

#define H2 0
#define H 1
#define O2 2
#define O 3
#define OH 4
#define HO2 5
#define H2O 6
#define H2O2 7
#define N2 8

#define Ith(v,i)    NV_Ith_S(v,i-1)
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1)


#define PI     RCONST(3.1415926)
#define E      RCONST(2.7182818)

/* Problem Constants */
extern double k_mol;
extern double Y_N2;
extern double Y_max;
extern double P;
extern double R;
extern double x_center;
extern double koeff_l;
extern double l;
extern long int myiter;
extern long int nniters;
extern double eps_x ;
extern double eps_func;
extern double eps_fr ;
extern double Tstart ;
extern double Tfinish ;
extern const double kB ;
extern const double Angstroem__ ;
extern const double santimetr ;
extern const vector<double> M ;
extern string name_species[9];
extern SUNContext sunctx;
extern void* cvode_mem;

typedef struct {
    realtype* x; //cells
    realtype* T; // temperature
    realtype* Yiprev;
    realtype* Yi; // mass fr in one cell
    realtype* Yinext; //grad mol fr
    realtype* Xiprev; // mol fr in one cell
    realtype* Xi;
    realtype* Xinext;
    realtype* gradX;
    realtype* X_tmp;
    realtype* Y_tmp;
    realtype* Y_left_bound;
    realtype* wk_add;
    realtype* forward;
    realtype* reverse;
    realtype* equilib;
    SUNContext sunctx;
    int Nx;
    int N_m;
    int NEQ;
    int N_centr;
    realtype Tl;
    realtype M;
    realtype T_center;
    realtype* y;  // molar cons
    realtype* ydot; // d(molar cons) / dt
    void* mykmem;
    IO::ChemkinReader* chemkinReader;
} *UserData;

static int num_gas_species = 9;
static int num_react = 22;

static int check_retval(void* retvalvalue, const char* funcname, int opt);

double F_right(IO::ChemkinReader* chemkinReader, double* Yiprev, double* Yi, double* Yinext,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext, double* Xiprev, double* Xi, double* Xinext, double* gradX, double* Y_tmp,
    double* X_tmp, double M, double* ydot, double* wk_add);

double F_rightY(IO::ChemkinReader* chemkinReader, double* Yiprev, double* Yi, double* Yinext,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext, double* Xiprev, double* Xi, double* Xinext, double* gradX, double* Y_tmp,
    double* X_tmp,
    double M, const int k_spec, double* ydot, double* wk_add);

static int func_Y(N_Vector u, N_Vector f, void* user_data);

int Integrate_Y(IO::ChemkinReader* chemkinReader_temp, int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb);

int InitialData(int& Nx, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect, 
    double& M, double Tstart, double Tfinish, double* Ystart, double* Yend);

void Write_to_file2(string str, ofstream& fout, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double M, int N_x, int number);

int integrate_Y_IDA(IO::ChemkinReader* chemkinReader_temp, int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb);

static int func_Y_IDA(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data);

double gauss_func(double A, double mu, double sigma, double x, int k_spec);

void Add_elem(vector<double>& T, vector<double>& Y, vector<double>& x, int& N_x, int& N_center, double b);

double tanh_spec(double A, double mu, double sigma, double x, int k_spec);
double tanh_T(double T_start, double T_finish, double mu, double sigma, double x);
double tanh_spec_minus(double A, double mu, double sigma, double x, int k_spec);

void Init_Data(UserData data, int N_x, vector<double>& x_vect,
    vector<double>& T_vect, IO::ChemkinReader* chemkinReader_temp, int NEQ,
    int N_center, double* Y_leftb);

void integrate_Y_CVODE(IO::ChemkinReader* chemkinReader_temp, int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb);

static int func_Y_CVODE(realtype t, N_Vector y, N_Vector ydot, void* user_data);

int FindBoundary(IO::ChemkinReader* chemkinReader_temp, double* Y_leftbound_inlet, double* Y_leftbound, double* Yinext, double Tl,
    double xi, double xinext, double M);

static int left_bound_solve(N_Vector u, N_Vector f, void* user_data);