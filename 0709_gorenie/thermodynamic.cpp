#include "thermodynamic.h"


double Lambda_H2(double T) {

    double res = 12.923981024360526 + -2.1362785572160212 * log(T) + 0.36550872867889844 * pow(log(T), 2) + -0.014973016232083858 * pow(log(T), 3);
    return  exp(res) * pow(10, -7);
}
double Lambda_H(double T) {
    double res = 0.6067604738391079 + 3.0031481452279447 * log(T) + -0.30253077202949874 * pow(log(T), 2) + 0.012922296045869641 * pow(log(T), 3);
    return  exp(res) * pow(10, -7);
} 
double Lambda_O2(double T) {
    double res = -2.1932210255373565 + 3.0298108191775883 * log(T) + -0.2951455273958285 * pow(log(T), 2) + 0.012881156880235112 * pow(log(T), 3);
    return  exp(res) * pow(10, -7);
}
double Lambda_O(double T) {
    double res = 2.7817609871583064 + 1.444073619402008 * log(T) + -0.10287650685055189 * pow(log(T), 2) + 0.004400427609296542 * pow(log(T), 3);
    return  exp(res) * pow(10, -7);
}
double Lambda_OH(double T) {
    double res = 14.826340814534335 + -3.5652377916037046 * log(T) + 0.584490887594824 * pow(log(T), 2) + -0.025931782533391123 * pow(log(T), 3);
    return  exp(res) * pow(10, -7);
}
double Lambda_HO2(double T) {
    double res = 5.238482478867441 + -0.4661770101568432 * log(T) + 0.24626382424375984 * pow(log(T), 2) + -0.013989211813681207 * pow(log(T), 3);
    return  exp(res) * pow(10, -7);
}
double Lambda_H2O(double T) {
    double res = 17.212214150754107 + -6.270146571625887 * log(T) + 1.1378516342304246 * pow(log(T), 2) + -0.05712359759599872 * pow(log(T), 3);
    return  exp(res) * pow(10, -7);
}
double Lambda_H2O2(double T) {
    double res = -0.5050153127654742 + 1.9352583603473448 * log(T) + -0.07235914136238857 * pow(log(T), 2) + 5.5703251616329164e-05 * pow(log(T), 3);
    return  exp(res) * pow(10, -7);
}

double Lambda_N2(double T) {
    double res = 7.995479312358585 + -1.3476582134062585 * log(T) + 0.3264331476770797 * pow(log(T), 2) + -0.016473702040685823 * pow(log(T), 3);
    return  exp(res) * pow(10, -7);
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
        cp_tmp += get_Cpi(k_spec, T) * Y[k_spec];
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
    return res / 2. ;
}