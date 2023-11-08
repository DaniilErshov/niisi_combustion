#pragma once
#include "functions.h"
#include "chemkinReader.h"
#include "boost/regex.hpp"
#include "boost/format.hpp"

double Lambda_H2(double T);
double Lambda_O2(double T);
double Lambda_H2O(double T);
double Lambda_N2(double T);

double myget_Cpi(int k, double T);

double myget_Hi(int k, double T);

double myget_enthalpy(int num_gas_speciens, double* Y, double T);

double myget_Cp(int num_gas_speciens, double* Y, double T);

double Cp_all(double T, double* Y);

double Lambda_All(IO::ChemkinReader* chemkinReader, double* Y, double T);
