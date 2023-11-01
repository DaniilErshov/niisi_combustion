#pragma once
#include <iostream>
#include <fstream>
#include "constants.h"
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

#define FTOL   RCONST(1.e-8) /* function tolerance */
#define STOL   RCONST(1.e-8) /* step tolerance     */

#define ZERO   RCONST(0.0)
#define PT25   RCONST(0.25)
#define PT5    RCONST(0.5)
#define ONE    RCONST(1.0)
#define ONEPT5 RCONST(1.5)
#define TWO    RCONST(2.0)

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
extern double koeff_l;
extern double l;
extern long int myiter;
extern long int nniters;
extern double eps_x ;
extern double eps_fr ;
extern double Tstart ;
extern double Tfinish ;
extern const double kB ;
extern const double Angstroem__ ;
extern const double santimetr ;
extern const vector<double> M ;

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
    realtype* Y_tmp;
    realtype* Y_left_bound;
    realtype* wk_add;
    int Nx;
    int N_m;
    int NEQ;
    int N_centr;
    realtype Tl;
    realtype M;
    realtype T_center;
    N_Vector y;  // molar cons
    N_Vector ydot; // d(molar cons) / dt
    void* mykmem;
    IO::ChemkinReader* chemkinReader;
} *UserData;

static int num_gas_species = 9;
static int num_react = 22;

static int check_retval(void* retvalvalue, const char* funcname, int opt);

void get_Y(double Y_H2O, double& Y_H2, double& Y_O2, double Y_N2);

double get_W(double* Y);

double get_rho(double* Y, double T);

void make_averageY(double* Yavg, double* Yi, double* Yipred);

void Get_mole_fr(double* X, double* Y);

void Get_molar_cons(double* X, double* Y, double T);

double Cp_all(double T, double* Y);

void get_grad(double* gradX, double* Xi, double* Xinext, double x, double xnext);

double Lambda_All(IO::ChemkinReader* chemkinReader, double* Y, double T);

double Dij_func(IO::ChemkinReader* chemkinReader, int i, int j, double T, double* Y);

double Dk_func(IO::ChemkinReader* chemkinReader, int i, double T, double* Y, double* X, int N);

double rhoYkWk(IO::ChemkinReader* chemkinReader, int k, double T, double* Y, double* gradX, double gradT);

double rhoYkVk(IO::ChemkinReader* chemkinReader, int k, double T, double* Y, double* gradX, double gradT);

//here send molar_cons
void chem_vel(double Tcurr, N_Vector y, N_Vector ydot);

double F_right(IO::ChemkinReader* chemkinReader, double* Yi, double* Yinext,
    double* T, double* Xiprev, double* Xi, double* Xinext, double* gradX, double* Y_tmp,
    double M, realtype* x_cells, const int i, N_Vector ydot, double* wk_add);

double F_rightY(IO::ChemkinReader* chemkinReader, double* Yiprev, double* Yi, double* Yinext,
    double* T, double* Xiprev, double* Xi, double* Xinext, double* gradX, double* Y_tmp,
    double M, realtype* x_cells, const int i, const int k_spec, N_Vector ydot, double* wk_add);

void MakeYvectors(double* Yiprev, double* Yi, double* Yinext, double* Y, double* Y_left_bound, int myNx, int i);

static int func_Y(N_Vector u, N_Vector f, void* user_data);

int Integrate_Y(IO::ChemkinReader* chemkinReader_temp, int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb);

int InitialData(int& Nx, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect, 
    double& M, double Tstart, double Tfinish, double* Ystart, double* Yend);

void Write_to_file2(string str, ofstream& fout, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double M, int N_x, int number);

int Find_final_state_IDA(IO::ChemkinReader* chemkinReader_temp, double& Tinitial, double* Y_vect);

static int func_final_state(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data);

int integrate_Y_IDA(IO::ChemkinReader* chemkinReader_temp, int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb);

static int func_Y_IDA(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data);

double get_GradRho(double* Yi, double* Yinext, double x, double xnext, double Ti, double Tinext);

double myget_Hi(int k, double T);
double myget_Cpi(int k, double T);
double my_mol_weight(int k);

void add_toChemVel(double* wk_add, double M, double* Yi, double* Yinext, double x, double xnext, double Ti, double Tinext);
double myget_enthalpy(int num_gas_speciens, double* Y, double T);
double myget_Cp(int num_gas_speciens, double* Y, double T);
double gauss_func(double A, double mu, double sigma, double x);