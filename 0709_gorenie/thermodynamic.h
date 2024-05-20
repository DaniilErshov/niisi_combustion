#pragma once
#include "functions.h"
#include "chemkinReader.h"
#include "boost/regex.hpp"
#include "boost/format.hpp"


double Lambda_All(double* X, double T, int number_cell, char side);

double get_Cp(const vector<double>& Cp_arr, double T, double* Y);

double Lambda_All(double* X, double T, int number_cell, char side);