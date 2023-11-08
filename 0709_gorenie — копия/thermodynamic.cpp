#include "thermodynamic.h"

double Lambda_H2(double T) {
    double res = (2.0168072101043486 - 1.93742201 * log(T) + 0.29103516 * pow(log(T), 2)
        - 0.00994456 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_O2(double T) {
    double res = (-3.0551023900370224 - 0.85094309 * log(T) + 0.18225493 * pow(log(T), 2)
        - 0.00696862 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_H2O(double T) {
    double res = (11.007354408138713 - 7.62523857 * log(T) + 1.23399258 * pow(log(T), 2)
        - 0.0583371 * pow(log(T), 3));
    return  exp(res);
}
double Lambda_N2(double T) {
    double res = (-3.0551023900370224 - 0.85094309 * log(T) + 0.18225493 * pow(log(T), 2)
        - 0.00696862 * pow(log(T), 3));
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

double Lambda_All(IO::ChemkinReader* chemkinReader, double* Y, double T)
{
    vector<double> lambda_spec(num_gas_species);
    double* X = new double[num_gas_species];
    double Ti;
    double integr_i;
    double sigma_i;
    double lamb1 = 0.;
    double lamb2 = 0.;
    double kBol = 1.380649 * pow(10, -16); //эрг
    double T_eps_i = 0;
    Get_mole_fr(X, Y);
    double sum = 0;

    for (int i = 0; i < num_gas_species; i++) {

        //std::cout << std::endl << std::endl << std::endl << chemkinReader->species()[i].name() << std::endl;
        sigma_i = chemkinReader->species()[i].transport().getCollisionDiameter() / Angstroem__;
        T_eps_i = chemkinReader->species()[i].transport().getPotentialWellDepth() / kB;

        //cout << "sigma_i = " << sigma_i << endl;
        Ti = T / (T_eps_i);
        //cout << "Ti = " << T_eps_i << endl;

        integr_i = 1.157 * pow(Ti, -0.1472);
        //cout << "integr_i = " << integr_i << endl;
        lambda_spec[i] = 8330. * pow(T / (my_mol_weight(i)), 0.5) / pow(sigma_i, 2.) / integr_i;
        //cout << "lambda_spec[i] = " << lambda_spec[i] << endl;
    }
    for (int i = 0; i < num_gas_species; i++) {
        lamb1 += X[i] * lambda_spec[i];
        lamb2 += X[i] / lambda_spec[i];
    }
    delete[] X;
    return 0.5 * (lamb1 + 1. / lamb2) * pow(10, -7);
}
