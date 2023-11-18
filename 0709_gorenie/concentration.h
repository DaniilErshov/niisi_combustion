#pragma once
#include "functions.h"
#include "chemkinReader.h"
#include "boost/regex.hpp"
#include "boost/format.hpp"

double my_mol_weight(int k);

void make_averageY(double* Yavg, double* Yi, double* Yipred);

void Get_mole_fr(double* X, double* Y);

void Get_molar_cons(double* X, double* Y, double T);

void get_Y(double Y_H2O, double& Y_H2, double& Y_O2, double Y_N2);

double get_W(double* Y);

void MakeYvectors(IO::ChemkinReader* chemkinReader_temp, double* x,
    double* Yiprev, double* Yi, double* Yinext,
    double* Y, double* Y_left_bound, int myNx, int i, double Tl, double M);

double get_rho(double* Y, double T);

double get_GradRho(double* Yi, double* Yinext, double x, double xnext, double Ti, double Tinext);

void get_grad(double* gradX, double* Xi, double* Xinext, double x, double xnext);

void chem_vel(double* forward, double* reverse, double* equilib, double* wk_add, double M, double* Yi, double* Yinext, double x, double xnext, double Tcurr, double Tinext, double* y, double* yprime);

void add_toChemVel(double* wk_add, double M, double* Yi, double* Yinext, double x, double xnext, double Ti, double Tinext);

double YkVk(IO::ChemkinReader* chemkinReader, int k, double T, double* Y, double* gradX, double* Xi);

double Dij_func(IO::ChemkinReader* chemkinReader, int i, int j, double T, double* Y);

double Dk_func(IO::ChemkinReader* chemkinReader, int i, double T, double* Y, double* X, int N);