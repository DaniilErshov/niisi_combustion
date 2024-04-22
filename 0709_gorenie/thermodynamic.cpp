#include "thermodynamic.h"
#include "constants.h"

extern vector<vector<double>> Cp_arr;
extern vector<vector<double>> Lambda_arr;
extern vector<vector<vector<double>>> Dij_arr_r;
extern vector<vector<vector<double>>> Dij_arr_l;


double myget_Cpi(int k, double T) {
    return get_Cpi(k, T);
}
double myget_Si(int k, double T) {
    return get_Si(k, T);
}

double myget_Hi(int k, double T) {
    return get_Hi(k, T);
}

double myget_enthalpy(int num_gas_speciens, double* Y, double T) {
    return get_enthalpy(num_gas_species, Y, T);
}

double myget_Cp(int num_gas_speciens, double* Y, double T) {
    return get_Cp(num_gas_species, Y, T);
}

double get_Cp(const vector<double>& Cp_arr, double T, double* Y)
{
    double cp_tmp = 0.;
    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        cp_tmp += Cp_arr[k_spec] * Y[k_spec];
    };
    return cp_tmp;
}

double Lambda_All(const vector<double>& Lambda_arr, double* X)
{
    double res1 = 0;
    double res2 = 0;

    for (int i = 0; i < num_gas_species; i++) {
        res1 += X[i] * Lambda_arr[i];
        res2 += X[i] / Lambda_arr[i];
    }
    return (res1 + 1. / res2) / 2. * pow(10, -7);
}

double Lambda_All(double* X, double T)
{
    double res1 = 0;
    double res2 = 0;

    for (int i = 0; i < num_gas_species; i++) {
        double lambda_i = get_Lambda5(i, T);
        res1 += X[i] * lambda_i;
        res2 += X[i] / lambda_i;
    }


    return (res1 + 1. / res2) / 2. * pow(10, -7);
}
