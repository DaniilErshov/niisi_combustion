#include "thermodynamic.h"

double Lambda_H2(double T) {
    double res = (1.6695115531878237 - 1.76734716 * log(T) + 0.26422876 * pow(log(T), 2)
        - 0.00857019 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_H(double T) {
    double res = (-0.2054145803692773 - 0.98198603 * log(T) + 0.19866012 * pow(log(T), 2)
        - 0.00836508 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_O2(double T) {
    double res = (-2.4941624666413817 - 1.09865955 * log(T) + 0.21931218 * pow(log(T), 2)
        - 0.00881651 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_O(double T) {
    double res = (-1.198976718323006 - 1.24807805 * log(T) + 0.22249078 * pow(log(T), 2)
        - 0.00891554 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_OH(double T) {
    double res = (2.8044872137363575 - 2.886673 * log(T) + 0.44357776 * pow(log(T), 2)
        - 0.01793352 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_HO2(double T) {
    double res = (0.31835676143360025 - 2.45120819 * log(T) + 0.43651758 * pow(log(T), 2)
        - 0.01962151 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_H2O(double T) {
    double res = (7.8814997570835645 - 6.24408052 * log(T) + 1.0331172 * pow(log(T), 2)
        - 0.04869067 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_H2O2(double T) {
    double res = (0.5250328719150514 - 2.57886269 * log(T) + 0.47538552 * pow(log(T), 2)
        - 0.02237828 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_N2(double T) {
    double res = (2.2276710981921695 - 3.23746048 * log(T) + 0.53623682 * pow(log(T), 2)
        - 0.02439004 * pow(log(T), 3));
    return  exp(res);
}


double myget_Cpi(int k, double T) {
    return get_Cpi(k, T) * pow(10, 3);
}

double myget_Hi(int k, double T) {
    return get_Hi(k, T) * pow(10, 3);
}

double myget_enthalpy(int num_gas_speciens, double* Y, double T) {
    return get_enthalpy(num_gas_species, Y, T) * pow(10, 3);
}

double myget_Cp(int num_gas_speciens, double* Y, double T) {
    return get_Cp(num_gas_species, Y, T) * pow(10, 3);
}

double Cp_all(double T, double* Y)
{
    double cp_tmp = 0.;
    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        cp_tmp += myget_Cpi(k_spec, T) * Y[k_spec];
    };
    return cp_tmp;
}

double Lambda_All(double* Y, double T)
{
    double* X = new double[num_gas_species];
    Get_mole_fr(X, Y);
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