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
#include <sundials/sundials_math.h>     /* access to SUNRexp               */
#include <sunmatrix/sunmatrix_band.h>  /* access to band SUNMatrix             */
#include <sunlinsol/sunlinsol_band.h>  /* access to band SUNLinearSolver       */
//#include <cvode/cvode.h>   


#define FTOL   RCONST(1.e-7)/* function tolerance */
#define STOL   RCONST(1.e-20) /* step tolerance     */

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

/* Problem Constants */
extern int num_gas_species;
extern int num_react;
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
extern double nevyaz_Y;
extern double nevyaz_T;
extern const double kB ;
extern const double Angstroem__ ;
extern const double santimetr ;
extern double** Dij_res;

extern double* Yi;
extern double* Yiprev;
extern double* Yinext;

extern double* YkVk_r;
extern double* YkVk_l;

extern double* gradX_r;
extern double* gradX_l;

extern double* X_tmp_r;
extern double* X_tmp_l;

extern double* Y_tmp_r;
extern double* Y_tmp_l;

extern double* Xiprev;
extern double* Xi;
extern double* Xinext;

extern double* gradX;

extern double* Y_tmp;
extern double* X_tmp;
extern double* YkVk;

extern double* Sn;
extern double* Hn;
extern double* Cpn;

extern double* forward_arr;
extern double* reverse_arr;
extern double* equilib_arr;
extern double* Y_left_bound;
extern double* wk_add;
extern double* ydot;
extern double*** diff_polynom;
extern double** lambda_polynom;
extern double* mol_weight;
extern int ida_steps;
extern double eps;
extern vector<string> name_species;
extern std::unordered_map<std::string, int> komponents;
extern std::unordered_map<int, string> komponents_str;
extern std::unordered_map<int, std::unordered_map<string, double>> Dij_saved;
extern std::map<std::string, int> mykomponents;
extern std::map<int, string> mykomponents_str;
extern vector<double> x_vect;
extern vector<double> Y_vect;
extern vector<double> T_vect;

typedef struct {
    double Vc_r, Vc_l, rho_r, rho_l;
    int Nx;
    int N_m;
    int NEQ;
    int NEQ_Y;
    int N_centr;
    realtype Tl;
    realtype M;
    realtype T_center;
    IO::ChemkinReader* chemkinReader;
} *UserData;


static int check_retval(void* retvalvalue, const char* funcname, int opt);

double F_right(UserData data,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext);

double F_rightY(UserData data, int k_spec,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext);


int InitialData(int& Nx, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect, double& M, double Tstart, double Tfinish, double* Ystart, double* Yend);

void Write_to_file(string str, ofstream& fout, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, vector<double>& Yp_vect, double M, int N_x, int number);

int integrate_Y_IDA(int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb, double t_fix);

static int func_Y_IDA(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data);

int integrate_All_IDA(int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb, int iter, double t_fix);

static int func_All_IDA(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data);

double get_M(double* Yiprev, double* Yi, double* Yinext,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext, double* Xiprev, double* Xi, double* Xinext, double* gradX, double* Y_tmp, double* X_tmp,
    double M, double* ydot, double* wk_add);

void Add_elem_simple(vector<double>& T, vector<double>& Y, vector<double>& x, int& N_x, int& N_center, double b, int number, int number_start, double& T_center);

void Init_Data(UserData data, int N_x, vector<double>& x_vect,
    vector<double>& T_vect, int NEQ,
    int N_center, double* Y_leftb);

static int func_Y_IDA(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data);

int integrate_Y_IDA(int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb);

int Integrate_Kinsol(int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb, int iter);

static int func_kinsol(N_Vector u, N_Vector f, void* user_data);

void MakeYvectors(UserData data,
    double* Y, int myNx, int i, double Tl);

int Find_final_state_IDA(double& Tinitial, double& Tend, double* Y_vect, double* Y_end);
static int func_final_state(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data);

void makeYstart(double koeff_topl, double* Ystart);

void find_diff_slag(UserData data, double Tcurr, double Tnext, double* Yi, double* Yinext,
    double* Xi, double* Xinext, double* Ykvk_side, double* Y_tmp_side, double* X_tmp_side, double* gradX_side, double& rho_side, double& Vc_side, int i);

std::vector<std::string> splitString(std::string str, char splitter);

template <typename T>
void findValue(const std::vector<T>& data, bool(*condition)(T));

void set_polynom(double** ploynom, std::string name_file, std::string type_polynom);

void set_polynom_diffusion(double*** polynom, std::string name_file, std::string type_polynom);

void MakeYvectorsY(UserData data,
    double* Y, int myNx, int i, double Tl);

double get_M(double Tprev, double T, double Tnext,
    double xprev, double x, double xnext);


int Integrate_Kinsol_dense(int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb, int iter);

static int func_kinsol_dense(N_Vector u, N_Vector f, void* user_data);


void MakeYvectors_dense(UserData data,
    double* Y, int myNx, int i, double Tl);

void MakeYvectors_kins(UserData data,
    double* Y, int myNx, int i, double Tl);


int integrate_All_IDA_M(int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb, int iter, double t_fix);


static int func_All_IDA_M(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data);

void set_Dij_res(double T);


int integrate_All_IDA_dense(int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb, int iter, double t_fix);


static int func_All_IDA_dense(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data);


int Integrate_Kinsol_withoutM(int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb, int iter);

static int func_kinsol_withoutM(N_Vector u, N_Vector f, void* user_data);