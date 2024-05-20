#include "thermodynamic.h"
#include "constants.h"

extern vector<vector<double>> Cp_arr;
extern vector<vector<double>> Lambda_arr;
extern vector<vector<vector<double>>> Dij_arr_r;
extern vector<vector<vector<double>>> Dij_arr_l;


double get_Cp(const vector<double>& Cp_arr, double T, double* Y)
{
    double cp_tmp = 0.;
    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        cp_tmp += Cp_arr[k_spec] * Y[k_spec];
    };
    return cp_tmp;
}


double Lambda_All(double* X, double T, int number_cell, char side)
{
    double res1 = 0;
    double res2 = 0;

    for (int i = 0; i < num_gas_species; i++) {
        double lambda_i = get_Lambda5(i, T, number_cell, side);
        res1 += X[i] * lambda_i;
        res2 += X[i] / lambda_i;
    }


    return (res1 + 1. / res2) / 2. * pow(10, -7);
}
