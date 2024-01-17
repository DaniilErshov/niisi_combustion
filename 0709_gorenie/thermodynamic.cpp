#include "thermodynamic.h"


double Lambda_H2(double T) {
    double res = (0.7952908032065458 - 1.7543618 * log(T) + 0.29620426 * pow(log(T), 2)
        - 0.01106313 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_H(double T) {
    double res = (-10.794768716134467 + 2.97505023 * log(T) - 0.30072807  * pow(log(T), 2)
        + 0.01292112 * pow(log(T), 3));
    return  exp(res);
} 
double Lambda_O2(double T) {
    double res = (-13.972838830946268 + 3.14253445 * log(T) - 0.31033452 * pow(log(T), 2)
        + 0.01354363 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_O(double T) {
    double res = (-10.510253685474222 + 2.23006088 * log(T) - 0.217143 * pow(log(T), 2)
        + 0.00988099 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_OH(double T) {
    double res = (4.29513792933524 - 3.97497223 * log(T) + 0.64269133 * pow(log(T), 2)
        - 0.02873229 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_HO2(double T) {
    double res = (-13.413518582354056 + 2.671689 * log(T) - 0.20882322 * pow(log(T), 2)
        + 0.00782351 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_H2O(double T) {
    double res = (16.566957613965787 - 11.01054783 * log(T) + 1.82299069 * pow(log(T), 2)
        - 0.0899335 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_H2O2(double T) {
    double res = (-10.106230046345894 + 1.12094883 * log(T) + 0.04337092 * pow(log(T), 2)
        - 0.00542249 * pow(log(T), 3));
    return  exp(res);
}

double Lambda_N2(double T) {
    double res = (0.8952973124090282 - 3.25428557 * log(T) + 0.60089845 * pow(log(T), 2)
        - 0.02961411 * pow(log(T), 3));
    return  exp(res);
}


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
        cp_tmp += myget_Cpi(k_spec, T) * Y[k_spec];
    };
    return cp_tmp;
}

double Lambda_All(double* X, double T)
{
    double res = Lambda_H2(T) * X[0] +
        Lambda_H(T) * X[1] +
        Lambda_O2(T) * X[2] +
        Lambda_O(T) * X[3] +
        Lambda_OH(T) * X[4] +
        Lambda_HO2(T) * X[5] +
        Lambda_H2O(T) * X[6] +
        Lambda_H2O2(T) * X[7] +
        Lambda_N2(T) * X[8];

    res += 1. / (X[0] / Lambda_H2(T) +
        X[1] / Lambda_H(T) +
        X[2] / Lambda_O2(T) +
        X[3] / Lambda_O(T) +
        X[4] / Lambda_OH(T) +
        X[5] / Lambda_HO2(T) +
        X[6] / Lambda_H2O(T) +
        X[7] / Lambda_H2O2(T) +
        X[8] / Lambda_N2(T));
    return res / 2. / 100.;
}