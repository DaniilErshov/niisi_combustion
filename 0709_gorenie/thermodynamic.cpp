#include "thermodynamic.h"
//#include "constants.h"

extern vector<vector<double>> Cp_arr;
extern vector<vector<double>> Lambda_arr;
extern vector<vector<vector<double>>> Dij_arr_r;
extern vector<vector<vector<double>>> Dij_arr_l;




double Lambda_All(double* X, double T, char phase)
{
    double res1 = 0;
    double res2 = 0;
    if (phase == 'g') {
        for (int i = 0; i < num_gas_species; i++) {
            double lambda_i = get_Lambda(i, T, phase);
            res1 += X[i] * lambda_i;
            res2 += X[i] / lambda_i;
        }
        return (res1 + 1. / res2) / 2. * pow(10, -7);
    }
    else {
        return (0.124e-2);
    }


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

double get_Hi(int component_i, double T)
{
    double Hi;
    int i = component_i;
    if (T > chec.chemkinReader->species()[component_i].thermo().getTCommon())
        Hi = phyc.Cp_coef_hT[i][0] + phyc.Cp_coef_hT[i][1] / 2. * T + phyc.Cp_coef_hT[i][2] / 3. * T * T
        + phyc.Cp_coef_hT[i][3] / 4. * T * T * T + phyc.Cp_coef_hT[i][4] / 5. * T * T * T * T + phyc.Cp_coef_hT[i][5] / T;
    else
        Hi = phyc.Cp_coef_lT[i][0] + phyc.Cp_coef_lT[i][1] / 2. * T + phyc.Cp_coef_lT[i][2] / 3. * T * T
        + phyc.Cp_coef_lT[i][3] / 4. * T * T * T + phyc.Cp_coef_lT[i][4] / 5. * T * T * T * T + phyc.Cp_coef_lT[i][5] / T;

    return Hi * T;
}

// Specific heat of ith component
double get_Cpi(int component_i, double T)
{
    double Cpi;
    int i = component_i;
    if (T > chec.chemkinReader->species()[component_i].thermo().getTCommon())
        Cpi = phyc.Cp_coef_hT[i][0] + phyc.Cp_coef_hT[i][1] * T + phyc.Cp_coef_hT[i][2] * T * T
        + phyc.Cp_coef_hT[i][3] * T * T * T + phyc.Cp_coef_hT[i][4] * T * T * T * T;
    else
        Cpi = phyc.Cp_coef_lT[i][0] + phyc.Cp_coef_lT[i][1] * T + phyc.Cp_coef_lT[i][2] * T * T
        + phyc.Cp_coef_lT[i][3] * T * T * T + phyc.Cp_coef_lT[i][4] * T * T * T * T;
    return Cpi;
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

double get_Cvi(int component_i, double T)
{
    double Cpi, Cvi;

    Cpi = get_Cpi(component_i, T);
    Cvi = Cpi - phyc.kR / mol_weight[component_i];
    return Cvi;
}

// Enthalpy of the gas
// Y -- mass fractions
double get_enthalpy(int num_species, double* Y, double T)
{
    double H_ = 0;

    for (int i = 0; i < num_species; i++)
        H_ += Y[i] * get_Hi(i, T);

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
double get_Cp(int num_species, double* Y, double T, char phase)
{   
    double Cp;
    if (phase == 'g') {
        Cp = 0.0;
        for (int i = 0; i < num_species; i++)
            Cp += Y[i] * get_Cpi(i, T);
    } 
    else {
        double T_test = T;
        Cp = 1374.5000478170068 + 13.496054354630873 * T_test + -0.11496503788479555 * pow(T_test, 2) + 0.0004478249492254928 * pow(T_test, 3) + -7.506490438296317e-07 * pow(T_test, 4) + 4.767065492350397e-10 * pow(T_test, 5);
        Cp *= 1.e-3;
    }
    
   // cout << "Cp " << phase << " = " << Cp << "\n";
    return Cp;
}

double get_Cv(int num_species, double* Y, double T, int number_cell)
{
    double Cv;

    Cv = 0.0;
    for (int i = 0; i < num_species; i++)
        Cv += Y[i] * get_Cvi(i, T);

    return Cv;
}

double get_Lambda(int i, double T, char phase)
{
    //double lambda_arg;
    //double logt = log(T);
    //lambda_arg = lambda_polynom[i][0] + lambda_polynom[i][1] * logt
    //    + lambda_polynom[i][2] * logt * logt + lambda_polynom[i][3] * logt * logt * logt + lambda_polynom[i][4] * logt * logt * logt * logt;
    //return pow(T, 0.5) * lambda_arg * pow(10, 5);
    if (phase == 'g') {
        double lambda_arg;
        double logt = log(T);
        lambda_arg = lambda_polynom[i][0] + lambda_polynom[i][1] * logt
            + lambda_polynom[i][2] * logt * logt + lambda_polynom[i][3] * logt * logt * logt;
        //cout << "Lambda " << " " << T << " " << name_species[i] << " = " << exp(lambda_arg) * pow(10, 5) << '\n';
        return exp(lambda_arg) * pow(10, 5);
    }
    else{
        return 0.56e-2;
    }
}
double Pf(double T) {
    //double constexpr n1 = 0.11670521452767e4;
    //double constexpr n2 = -0.72421316703206e6;
    //double constexpr n3 = -0.17073846940092e2;
    //double constexpr n4 = 0.12020824702470e5;
    //double constexpr n5 = -0.3232555022333e7;
    //double constexpr n6 = 0.14915108613530e2;
    //double constexpr n7 = -0.48232657361591e4;
    //double constexpr n8 = 0.40511340542057e6;
    //double constexpr n9 = -0.23855557567849;
    //double constexpr n10 = 0.65017534844798e3;
    //double theta = T + n9 / (T - n10);
    //double A = theta * theta + n1 * theta + n2;
    //double B = n3 * theta * theta + n4 * theta + n5;
    //double C = n6 * theta * theta + n7 * theta + n8;
    ////double lgP = 8.07131 - 1730.63 / (Tval - 273 + 233.426);
    ////return 133.33 * pow(10, lgP);
    //////std::cout << 1.0e6 * pow(2 * C / (-B + pow(B * B - 4 * A * C, 0.5)), 4) << std::endl;
    //return pow(2 * C / (-B + pow(B * B - 4 * A * C, 0.5)), 4);
    //double jop = pow(2 * C / (-B + pow(B * B - 4 * A * C, 0.5)), 4);
    double res = 24523.160692131387 + 922.0983808780411 * T + -21.500586403548954 * pow(T, 2) + 0.1552381071666568 * pow(T, 3) + -0.0004781197006891655 * pow(T, 4) + 5.439451804534201e-07 * pow(T, 5);
    res /= 1.e6;
    return res;
}
double L_d(double T) {
    //double Larr[] = { 2481115.000000, 2470000.000000, 2451000.000000, 2431000.000000, 
    //    2411000.000000, 2390000.000000, 2368000.000000, 2346000.000000, 2322000.000000, 2298000.000000, 2273000.000000, 2247000.000000, 2219000.000000 };
    //double Tarr[] = { 270, 280, 290, 300, 310, 320, 330, 340, 350, 360, 370, 380, 390 };
    //double T1, T2;
    //T1 = 0;
    //T2 = 0;
    //// return Ld;
    //for (int i = 0; i < 13; i++) {
    //    if (T > Tarr[11]) {
    //        return Larr[11] * 1.e-3;
    //    }
    //    if ((T - Tarr[i]) == 1.e-3) {
    //        return Larr[i] * 1.e-3;
    //    }
    //    if (T < Tarr[i]) {
    //        return (Larr[i] + (Larr[i - 1] - Larr[i]) * (Tarr[i] - T) / (Tarr[i] - Tarr[i - 1])) * 1.0e-3;
    //    }
    double res;
    res = 122143.71957118751 + 6224.984923585149 * T + -45.41337443600466 * pow(T, 2) + 0.15241169103015095 * pow(T, 3) + -0.0002533266591655675 * pow(T, 4) + 1.6413289216643534e-07 * pow(T, 5);
    res /= 1.e3;
    return res;
}

double get_rho_d(double T) {
    double res;
    res = 382.0229931050386 + 8.233888612039767 * T + -0.06055281806170509 * pow(T, 2) + 0.0002021369250560078 * pow(T, 3) + -3.3385255094635553e-07 * pow(T, 4) + 2.1563943686242482e-10 * pow(T, 5);
    res /= 1.e3;
    return 0.633;
    //return res;
}