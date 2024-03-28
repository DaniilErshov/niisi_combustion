#include "thermodynamic.h"
#include "constants.h"

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

double Cp_all(double T, double* Y)
{
    double cp_tmp = 0.;
    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        cp_tmp += get_Cpi(k_spec, T) * Y[k_spec];
    };
    return cp_tmp;
}

double Lambda_All(double* X, double T)
{
    double res1 = 0;
    double res2 = 0;

    for (int i = 0; i < num_gas_species; i++) {
        double lambda_i = get_Lambda(i, T);
        res1 += X[i] * lambda_i;
        res2 += X[i] / lambda_i;
    }


    return (res1 + 1. / res2) / 2. * pow(10, -7);
}
