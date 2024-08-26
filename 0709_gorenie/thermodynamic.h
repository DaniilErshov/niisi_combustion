#pragma once
#include "functions_sundials.h"
#include "chemkinReader.h"
#include "boost/regex.hpp"
#include "boost/format.hpp"

double Lambda_All(double* X, double T, char phase);

double get_dHiRT(double* Cp_coef, double T);

double get_dSiR(double* Cp_coef, double T);
double get_dCpi(double* Cp_coef, double T);

double get_Hi(int component_i, double T);

// Specific heat of ith component
double get_Cpi(int component_i, double T);

double get_Si(int component_i, double T);

double get_Cvi(int component_i, double T);

// Enthalpy of the gas
// Y -- mass fractions
double get_enthalpy(int num_species, double* Y, double T, int number_cell);

// P = rho * R * T
// R -- gas constant
// Y -- mass fractions
double get_gas_constant(int num_gas_species, double* Y);

// Specific heat of the gas
double get_Cp(int num_species, double* Y, double T, char phase);

double get_Cv(int num_species, double* Y, double T, int number_cell);
double get_Lambda(int i, double T, char phase);
double get_Lambda5(int i, double T, int number_cell, char side);
double L_d(double T);
double get_rho_d(double T);
double Pf(double T);