#include "thermodynamic.h"
//#include "constants.h"

extern vector<vector<double>> Cp_arr;
extern vector<vector<double>> Lambda_arr;
extern vector<vector<vector<double>>> Dij_arr_r;
extern vector<vector<vector<double>>> Dij_arr_l;




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


double get_dHiRT(double* Cp_coef, double T)
{
    double Hi;
    Hi = Cp_coef[0] + Cp_coef[1] / 2. * T + Cp_coef[2] / 3. * T * T
        + Cp_coef[3] / 4. * T * T * T + Cp_coef[4] / 5. * T * T * T * T + Cp_coef[5] / T;
    return Hi * T;
}

double get_dSiR(double* Cp_coef, double T) {
    double Si;
    Si = Cp_coef[0] * log(T) + Cp_coef[1] * T + Cp_coef[2] / 2. * T * T
        + Cp_coef[3] / 3. * T * T * T + Cp_coef[4] / 4. * T * T * T * T + Cp_coef[6];
    return Si;
}
double get_dCpi(double* Cp_coef, double T)
{
    double Cpi;
    Cpi = Cp_coef[0] + Cp_coef[1] * T + Cp_coef[2] * T * T
        + Cp_coef[3] * T * T * T + Cp_coef[4] * T * T * T * T;

    return Cpi;
}

double get_Hi(int component_i, double T, int number_cell)
{
    double Hi;
    int i = component_i;
    if (!flag_use_save_koeffs) {
        if (T > chec.chemkinReader->species()[component_i].thermo().getTCommon())
            Hi = phyc.Cp_coef_hT[i][0] + phyc.Cp_coef_hT[i][1] / 2. * T + phyc.Cp_coef_hT[i][2] / 3. * T * T
            + phyc.Cp_coef_hT[i][3] / 4. * T * T * T + phyc.Cp_coef_hT[i][4] / 5. * T * T * T * T + phyc.Cp_coef_hT[i][5] / T;
        else
            Hi = phyc.Cp_coef_lT[i][0] + phyc.Cp_coef_lT[i][1] / 2. * T + phyc.Cp_coef_lT[i][2] / 3. * T * T
            + phyc.Cp_coef_lT[i][3] / 4. * T * T * T + phyc.Cp_coef_lT[i][4] / 5. * T * T * T * T + phyc.Cp_coef_lT[i][5] / T;

        return Hi * T;
    }
    else {
        return H_arr[number_cell][component_i];
    }
}

// Specific heat of ith component
double get_Cpi(int component_i, double T, int number_cell)
{
    double Cpi;
    int i = component_i;
    if (!flag_use_save_koeffs) {
        if (T > chec.chemkinReader->species()[component_i].thermo().getTCommon())
            Cpi = phyc.Cp_coef_hT[i][0] + phyc.Cp_coef_hT[i][1] * T + phyc.Cp_coef_hT[i][2] * T * T
            + phyc.Cp_coef_hT[i][3] * T * T * T + phyc.Cp_coef_hT[i][4] * T * T * T * T;
        else
            Cpi = phyc.Cp_coef_lT[i][0] + phyc.Cp_coef_lT[i][1] * T + phyc.Cp_coef_lT[i][2] * T * T
            + phyc.Cp_coef_lT[i][3] * T * T * T + phyc.Cp_coef_lT[i][4] * T * T * T * T;
        return Cpi;
    }
    else {
        return Cp_arr[number_cell][component_i];
    }
}

double get_Si(int component_i, double T) {
    double Si;
    int i = component_i;

    if (T > chec.chemkinReader->species()[component_i].thermo().getTCommon())
        Si = phyc.Cp_coef_hT[i][0] * log(T) + phyc.Cp_coef_hT[i][1] * T + phyc.Cp_coef_hT[i][2] / 2. * T * T
        + phyc.Cp_coef_hT[i][3] / 3. * T * T * T + phyc.Cp_coef_hT[i][4] / 4. * T * T * T * T + phyc.Cp_coef_hT[i][6];
    else
        Si = phyc.Cp_coef_lT[i][0] * log(T) + phyc.Cp_coef_lT[i][1] * T + phyc.Cp_coef_lT[i][2] / 2. * T * T
        + phyc.Cp_coef_lT[i][3] / 3. * T * T * T + phyc.Cp_coef_lT[i][4] / 4. * T * T * T * T + phyc.Cp_coef_lT[i][6];
    return Si;
}

double get_Cvi(int component_i, double T, int number_cell)
{
    double Cpi, Cvi;

    Cpi = get_Cpi(component_i, T, number_cell);
    Cvi = Cpi - phyc.kR / mol_weight[component_i];
    return Cvi;
}

// Enthalpy of the gas
// Y -- mass fractions
double get_enthalpy(int num_species, double* Y, double T, int number_cell)
{
    double H_ = 0;

    for (int i = 0; i < num_species; i++)
        H_ += Y[i] * get_Hi(i, T, number_cell);

    return H_;
}

// P = rho * R * T
// R -- gas constant
// Y -- mass fractions
double get_gas_constant(int num_gas_species, double* Y)
{
    // Gas Constant
    double R = 0;

    for (int i = 0; i < num_gas_species; i++)
        R += Y[i] / mol_weight[i];

    R *= phyc.kR;

    return R;
}

// Specific heat of the gas
double get_Cp(int num_species, double* Y, double T, int number_cell)
{
    double Cp;

    Cp = 0.0;
    for (int i = 0; i < num_species; i++)
        Cp += Y[i] * get_Cpi(i, T, number_cell);

    return Cp;
}

double get_Cv(int num_species, double* Y, double T, int number_cell)
{
    double Cv;

    Cv = 0.0;
    for (int i = 0; i < num_species; i++)
        Cv += Y[i] * get_Cvi(i, T, number_cell);

    return Cv;
}

double get_Lambda(int i, double T, int number_cell, char side)
{
    //double lambda_arg;
    //double logt = log(T);
    //lambda_arg = lambda_polynom[i][0] + lambda_polynom[i][1] * logt
    //    + lambda_polynom[i][2] * logt * logt + lambda_polynom[i][3] * logt * logt * logt + lambda_polynom[i][4] * logt * logt * logt * logt;
    //return pow(T, 0.5) * lambda_arg * pow(10, 5);
    if (!flag_use_save_koeffs) {
        double lambda_arg;
        double logt = log(T);
        lambda_arg = lambda_polynom[i][0] + lambda_polynom[i][1] * logt
            + lambda_polynom[i][2] * logt * logt + lambda_polynom[i][3] * logt * logt * logt;
        //cout << "Lambda " << " " << T << " " << name_species[i] << " = " << exp(lambda_arg) * pow(10, 5) << '\n';
        return exp(lambda_arg) * pow(10, 5);
    }
    else {
        if (side == 'r')
            return Lambda_arr_r[number_cell][i];
        if (side == 'l')
            return Lambda_arr_l[number_cell][i];
        if (side == 'c')
            return Lambda_arr[number_cell][i];
    }
}

double get_Lambda5(int i, double T, int number_cell, char side)
{
    //double lambda_arg;
    //double logt = log(T);
    //lambda_arg = lambda_polynom[i][0] + lambda_polynom[i][1] * logt
    //    + lambda_polynom[i][2] * logt * logt + lambda_polynom[i][3] * logt * logt * logt + lambda_polynom[i][4] * logt * logt * logt * logt;
    //return pow(T, 0.5) * lambda_arg * pow(10, 5);
    if (!flag_use_save_koeffs) {
        double lambda_arg;
        double logt = log(T);
        lambda_arg = lambda_polynom[i][0] + lambda_polynom[i][1] * logt
            + lambda_polynom[i][2] * logt * logt + lambda_polynom[i][3] * logt * logt * logt;
        //cout << "Lambda " << " " << T << " " << name_species[i] << " = " << exp(lambda_arg) * pow(10, 5) << '\n';
        return exp(lambda_arg) * pow(10, 5);
    }
    else {
        if (side == 'r')
            return Lambda_arr_r[number_cell][i];
        if (side == 'l')
            return Lambda_arr_l[number_cell][i];
        if (side == 'c')
            return Lambda_arr[number_cell][i];
    }
}
