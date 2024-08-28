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
#include <kinsol/kinsol.h>  
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

#define Ith(v,i)    NV_Ith_S(v,i-1)
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1)


#define PI     RCONST(3.1415926)

/* Problem Constants */
extern double p_inter;
extern int num_gas_species;
extern int num_react;
extern double k_mol;
extern double Y_N2;
extern double Y_max;
extern double P;
extern double* Ystart;
extern double R;

extern double* Xi_2;
extern double* Xi_3;
extern double* X_inter;
extern string Fuel;


extern double Tstart ;
extern double Tfinish ;
extern double vel_prev;
extern const double kB ;
extern double t_curr;
extern const double Angstroem__ ;
extern const double santimetr ;
extern double** Dij_res;

extern double* Yi;
extern double* Yiprev;
extern double* Yinext;
extern int Nx;

extern double* YkVk_r;
extern double* YkVk_l;

extern double r_inter;
extern double u_inter;
extern double u_inter_2r;
extern double u_inter_3r;
extern double* Yend;

extern double* Y_inter;
extern double* Y_inter_2r;
extern double* Y_inter_3r;

extern double* gradX_r;
extern double* gradX_l;
extern double corrector;

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

extern int add_variable_count;
extern int count_var_in_cell;
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
extern int preinter;
extern vector<string> name_species;
extern std::unordered_map<std::string, int> komponents;
extern std::unordered_map<int, string> komponents_str;

extern vector<double> x_vect;
extern vector<double> drhodt_vect;
extern vector<double>  dTdt_vect;
extern vector<double>  dWdt_vect;

typedef struct Cell_Properties {
    vector <double> Y;
    double T;
    double u;
    double vel;
    double rho;
};


extern Cell_Properties Cell_Properties_inter;
extern Cell_Properties Cell_prouds_inter;
extern Cell_Properties Cell_rval_inter;

extern vector<Cell_Properties> Cell_Properties_vector;
extern vector<Cell_Properties> Cell_prouds_vector;
extern vector<Cell_Properties> Cell_rval_vector;

typedef struct {
    double Vc_r, Vc_l, rho_r, rho_l;
    int Nx;
    int N_m;
    int NEQ;
    int NEQ_Y;
    int N_centr;
    int my_numjac;
    double my_tcur;
    double t;
    double p;
    void* sun_mem;
    realtype Tl;
    realtype M;
    realtype T_center;
    IO::ChemkinReader* chemkinReader;
} *UserData;


double F_right_u(UserData data,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext,
    double uprev, double u, double unext,
    int number_cell);

void resize_koeff_vectors(int N_x);

void updateKoeffs(double* yval, UserData data);

static int check_retval(void* retvalvalue, const char* funcname, int opt);

double F_right_T_g(UserData data,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext,
    double uprev, double u, double unext, int number_cell);

void Make_Xvector(UserData data, double* Y, double* X_mol, int myNx, int i);

double F_right_T_d(UserData data,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext,
    double uprev, double u, double unext, int number_cell);

double F_right_u_inter_r(UserData data,
    double Tprev, double T, double Tnext,
    double u_inter, double u, double u_inter_2r, double u_inter_3r, double h, double p);

double F_rightY(UserData data, int k_spec,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext, 
    double uprev, double u, double unext, int number_cell);

double F_rightY_interface(UserData data, int k_spec, double Vc,
    double T_inter, double xprev, double x, double u, double us, double p, int number_cell);

double  F_right_u_inter(UserData data, double Vc, int k_spec,
    double T_inter, double xprev, double x, double u, double us, double p, int number_cell);

double F_right_T_interfase_r(double* my_X_inter, double T_inter,
    double T_inter_r, double T_inter_2r, double T_inter_3r, double u, double ri_l, double ri, double ri_r, double h, double p);
double F_right_T_interfase(double* my_X_tmp, double T_inter_3l, double T_inter_2l,
    double T_inter, double T_inter_2r, double T_inter_3r, double us, double ri, double h, double p);

void InitialData(int& Nx, vector<double>& x_vect, vector<Cell_Properties>& Cell_Properties_vector, double Tstart
    , double Tfinish, double* Ystart, double* Yend);

void Write_to_file(string str, string type_str, vector<Cell_Properties>& my_Cell_Properties_vector,
    Cell_Properties& my_Cell_Properties_inter);
void Write_drhodt(string str);

int integrate_All_IDA(int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb, int iter, double t_fix);

static int func_All_IDA(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data);

void Add_elem_simple(vector<double>& T, vector<double>& Y, vector<double>& x, int& N_x, int& N_center, double b, int number, int number_start, double& T_center);

void Init_Data(UserData data, int N_x, int NEQ);


int Find_final_state_IDA(double& Tinitial, double& Tend, double* Y_vect, double* Y_end);
static int func_final_state(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data);


void find_diff_slag(UserData data, double Tcurr, double Tnext, double* Yi, double* Yinext,
    double* Xi, double* Xinext, double* Ykvk_side, double* Y_tmp_side, double* X_tmp_side, double* gradX_side, double& rho_side, double& Vc_side, int i, char side);


double get_M(double Tprev, double T, double Tnext,
    double xprev, double x, double xnext, int number_cell);



int integrate_All_IDA_M(int N_x);


static int func_All_IDA_M(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data);

void set_Dij_res(double T);

void makeYstart(double koeff_topl, string fuel, double O2_in, double N2_in, double* Ystart);
void get_Yi(UserData data, double* yval, double* my_Yi, int i);
void MakeYvectors_interface(double* yval);

void set_fuel_rval(UserData data, double T_prev, double T_curr, double T_next,
    double u_prev, double u_curr, double u_next,
    double vel_prev, double vel_curr, double vel_next,
    double rho_prev, double rho_curr, double rho_next,
    double* rval, double* yval, double* ypval, int& i_temp, int i);


void set_interface_l_rval(UserData data, double T_prev, double T_curr, double T_next,
    double u_prev, double u_curr, double u_next,
    double vel_prev, double vel_curr, double vel_next,
    double rho_prev, double rho_curr, double rho_next,
    double T_inter_3l, double T_inter_2l, double T_inter,
    double* rval, double* yval, double* ypval, int& i_temp, double vel_inter, double rho_inter, int i);

void set_interface_rval(UserData data, double T_prev, double T_curr, double T_next,
    double u_prev, double u_curr, double u_next,
    double vel_prev, double vel_curr, double vel_next,
    double rho_prev, double rho_curr, double rho_next,
    double T_inter_3l, double T_inter_2l, double T_inter, double T_inter_2r, double T_inter_3r,
    double* rval, double* yval, double* ypval, int& i_temp, double vel_inter, double p, double rho_inter, double u_inter, int i);



void set_interface_r_rval(UserData data, double T_prev, double T_curr, double T_next,
    double u_prev, double u_curr, double u_next,
    double vel_prev, double vel_curr, double vel_next,
    double rho_prev, double rho_curr, double rho_next,
    double T_inter_3l, double T_inter_2l, double T_inter, double T_inter_2r, double T_inter_3r,
    double rho_inter, double rho_inter_2r, double rho_inter_3r,
    double* rval, double* yval, double* ypval, int& i_temp, double vel_inter, int i);


void set_rval_gas(UserData data, double T_prev, double T_curr, double T_next,
    double u_prev, double u_curr, double u_next,
    double vel_prev, double vel_curr, double vel_next,
    double rho_prev, double rho_curr, double rho_next,
    double* rval, double* yval, double* ypval, int& i_temp, double vel_inter, int i);

void set_p(double tres, double dt, double step, double V, double h);

double get_vel(double* my_X_tmp, double T_inter_3l, double T_inter_2l,
    double T_inter, double T_inter_2r, double T_inter_3r, double h, double p);
double inital_u(UserData data, double Vc, double T_inter,
    double u, double us, int number_cell);

double F_right_rho(UserData data,
    double rho_l, double rho, double rho_r, double r_l, double r, double r_r,
    double u_l, double u, double u_r,
    int number_cell);

double F_right_rho_inter_r(UserData data,
    double rho_inter, double rho, double rho_inter_2r, double rho_inter_3r,
    double u_inter, double u, double u_inter_2r, double u_inter_3r, double r_l, double r, double r_r, double p);

double F_right_rho_interface(double rhog, double vel, double u_inter, double T_inter);


int KinSetIc(int NEQ);

static int func_kinsol(N_Vector u, N_Vector f, void* user_data);
double get_Qg(double Tval, double Tvalr, double TvalrP, double* Yval, double h, double p);
double get_Qd(double TvallM, double Tvall, double Tval, double* Yval, double h, double p);

void MakeYvector(vector<Cell_Properties>& my_Cell_Properties_vector,
    Cell_Properties& my_Cell_Properties_inter, double* Y, int myNx);

void MakePropertiesvectors(vector<Cell_Properties>& my_Cell_Properties_vector,
    Cell_Properties& my_Cell_Properties_inter,
    double* Y, int myNx);

double get_dotM(double Ys, double Tval, double ri, double rhog, double Dfm);
