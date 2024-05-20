#include <cmath>
#include <iostream>
#include "functions.h"

using namespace std;

struct phy_consts phyc;
struct che_consts chec;

extern vector<vector<double>> Cp_arr;
extern vector<vector<double>> Lambda_arr;
extern vector<vector<vector<double>>> Dij_arr_r;
extern vector<vector<vector<double>>> Dij_arr_l;



extern std::unordered_map<std::string, int> komponents;
extern std::unordered_map<int, string> komponents_str;
extern vector<string> name_species;
extern map<string, double> elem_mol_weight;
extern double*** diff_polynom;
extern double** lambda_polynom;
extern double** Dij_func_saved;
extern double* mol_weight;

extern double* Yi;
extern double* Yiprev;
extern double* Yinext;

extern double* YkVk_r;
extern double* YkVk_l;

extern double* gradX_r;
extern double* gradX_l;

extern double* X_tmp_r;
extern double* X_tmp_l;

extern double* Y_tmp_r;
extern double* Y_tmp_l;

extern double* Xiprev;
extern double* Xi;
extern double* Xinext;

extern double* gradX;

extern double* Y_tmp;
extern double* X_tmp;
extern double* YkVk;

extern double* Sn;
extern double* Hn ;
extern double* Cpn;

extern double* forward_arr;
extern double* reverse_arr;
extern double* equilib_arr;
extern double* Y_left_bound;
extern double* wk_add;
extern double* ydot;

extern double*** diff_polynom;
extern double** lambda_polynom;
extern double* mol_weight;
extern double* Ystart;
extern double* Yend;
extern double* X;
extern double* YkVk_res;


void init_consts(int& num_gas_species, int& num_react) {
    std::string name_file = R"(D:\Storage\Daniil\check_read_file\check_read_file\n-heptane-tran.out)";
    std::string cond_str = "CONDUCTIVITIES";
    std::string visc_str = "VISCOSITIES";
    std::string diff_str = "DIFFUSION";

    const std::string chemfile = R"(D:\Storage\Daniil\n-heptane\heptane.inp)";
    const std::string thermfile = R"(D:\Storage\Daniil\n-heptane\term.dat)";
    const std::string transfile = R"(D:\Storage\Daniil\n-heptane\tran.dat)";
    //const std::string chemfile = R"(D:\Storage\Daniil\gorenie\ChemKin_reader\test\chem2.inp)";
    //const std::string thermfile = R"(D:\Storage\Daniil\gorenie\ChemKin_reader\test\therm-abf.dat)";
    //const std::string transfile = R"(D:\Storage\Daniil\gorenie\ChemKin_reader\test\tran.dat)";

    chec.chemkinReader = new IO::ChemkinReader(chemfile, thermfile, transfile);


    chec.chemkinReader->read();
    //chec.chemkinReader->check();
    //cout << "37\n";
    //cout << chec.chemkinReader->reactions()[37] << "\n";

    cout << "16\n";
    cout << chec.chemkinReader->reactions()[16] << "\n";

    cout << "17\n";
    cout << chec.chemkinReader->reactions()[17] << "\n";

    cout << "193\n";
    cout << chec.chemkinReader->reactions()[193] << "\n";



    std::cout << "NEW BLOCK" << std::endl;
    num_react = chec.chemkinReader->reactions().size();
    num_gas_species = chec.chemkinReader->species().size();
    chec.sum_v = new double[num_react];
    //phyc.kR = 8.314472e+7;  // universal gas constant ( erg/K/mole )
    //phyc.kRc = 1.987207e-3; // universal gas constant ( kcal/K/mole )
    //phyc.kTime = 1.e-3;     // dimensional time ( sec )
    //phyc.kLength = 1.e+2;   // dimensional length ( cm )
    //phyc.kMass = 1.e+3;     // dimensional mass ( g )
    //phyc.kTemp = 1.e0;      // dimensional temperature ( K )
    //phyc.kPres = 1.e+7;     // dimensional pressure ( dyn/cm**2 )
    //phyc.kDens = phyc.kMass / pow(phyc.kLength, 3); // dimensional density ( g/cm**3 )
    //phyc.kVel = phyc.kLength / phyc.kTime;           // dimensional velocity ( cm/sec )

    phyc.kR = 8.314472;  // universal gas constant ( erg/K/mole )
    phyc.kRc = 1.987207e-3; // universal gas constant ( kcal/K/mole )
    phyc.kTime = 1;     // dimensional time ( sec )
    phyc.kLength = 1;   // dimensional length ( cm )
    phyc.kMass = 1;     // dimensional mass ( g )
    phyc.kTemp = 1.e0;      // dimensional temperature ( K )
    phyc.kPres = 1.e+7;     // dimensional pressure ( dyn/cm**2 )
    phyc.kDens = phyc.kMass / pow(phyc.kLength, 3); // dimensional density ( g/cm**3 )
    phyc.kVel = phyc.kLength / phyc.kTime;           // dimensional velocity ( cm/sec )

    phyc.kR /= phyc.kMass * pow(phyc.kVel, 2);

    int i_specie = 0;
    for (const auto& specie_i : chec.chemkinReader->species()) {
        komponents[specie_i.name()] = i_specie;
        komponents_str[i_specie] = specie_i.name();
        //cout << specie_i.name() << " = " << komponents[specie_i.name()] << "\n";
        name_species.push_back(specie_i.name());
        i_specie++;
    }

    
    allocate_memory();

    double d = 0.14;
    int k = 0, l = 0;
    double logF_f, logF_core_f, logF_r, logF_core_r;

    double sum1, sum2;
    double Kci;
    double Kpi;
    double dSiR, dHiRT;
    double sumv = 0;
    double k_0_f, k_inf_f, c, m, Pr_f, Fcent, F_f;
    double sum_ThirdBodies;
    bool M_exist = 1;
    auto& species = chec.chemkinReader->species();


    for (int i = 0; i < num_gas_species; i++) {
        auto koeff_vect = chec.chemkinReader->species()[i].thermo().getLowerTemperatureCoefficients();
        cout << komponents_str[i] << "\n";
        for (int j = 0; j < 9; j++) {
            phyc.Cp_coef_lT[i][j] = 0;
            phyc.Cp_coef_hT[i][j] = 0;
            //cout << "lt = " << phyc.Cp_coef_lT[i][j] << "\n";
        }

        for (int j = 0; j < koeff_vect.size(); j++) {
            phyc.Cp_coef_lT[i][j] = koeff_vect[j];
            //cout << "lt = " << phyc.Cp_coef_lT[i][j] << "\n";
        }

        koeff_vect = chec.chemkinReader->species()[i].thermo().getUpperTemperatureCoefficients();
        for (int j = 0; j < koeff_vect.size(); j++) {
            phyc.Cp_coef_hT[i][j] = koeff_vect[j];
            //cout << "ht = " << phyc.Cp_coef_hT[i][j] << "\n";
        }
    }

    phyc.mol_weight = new double[num_gas_species];
    auto elm = chec.chemkinReader->elements();
    for (int i = 0; i < num_gas_species; i++) {
        phyc.mol_weight[i] = 0;
        auto elms = chec.chemkinReader->species()[i].thermo().getElements();
        for (const auto& elm_i : elms) {
            phyc.mol_weight[i] += elem_mol_weight[elm_i.first] * elm_i.second;
        }
    }


    for (int i = 0; i < num_gas_species; i++)
        phyc.mol_weight[i] /= phyc.kMass;

    // 8th coef. is for enthalpy calculation

    GasTransport poly_obj;
    for (int k = 0; k < num_gas_species; k++) {
        double logt = log(300);
        poly_obj.getConductivityPolynomial(k, lambda_polynom[k]);
        for (int k2 = 0; k2 < num_gas_species; k2++) {
            poly_obj.getBinDiffusivityPolynomial(k, k2, diff_polynom[k][k2]);
           // cout << name_species[k] << " " << name_species[k2] << " " << diff_polynom[k][k2][0] << " " << diff_polynom[k][k2][1] << " " << diff_polynom[k][k2][2] << " " << diff_polynom[k][k2][3] << "\n";
        }
    }
    
   /* set_polynom(lambda_polynom, name_file, cond_str);
    set_polynom_diffusion(diff_polynom, name_file, diff_str);*/
    //for (int k = 0; k < num_gas_species; k++) {
    //    for (int k2 = 0; k2 < num_gas_species; k2++) {
    //        cout << name_species[k] << " " << name_species[k2] << " " << diff_polynom[k][k2][0] << " " << diff_polynom[k][k2][1] << " " << diff_polynom[k][k2][2] << " " << diff_polynom[k][k2][3] << "\n";
    //    }
    //}

    for (int component_i = 0; component_i < num_gas_species; component_i++)
    {
        for (int power_i = 0; power_i <= 8; power_i++)
        {
            // T < 1000 K
            phyc.Cp_coef_lT[component_i][power_i] *= phyc.kR / phyc.mol_weight[component_i];
            // T > 1000 K
            phyc.Cp_coef_hT[component_i][power_i] *= phyc.kR / phyc.mol_weight[component_i];
        }
    }
}
   
double get_dHiRT(double* Cp_coef, double T)
{
    double Hi;
    Hi = Cp_coef[0] + Cp_coef[1] / 2.* T + Cp_coef[2] / 3. * T * T
        + Cp_coef[3] / 4. * T * T * T + Cp_coef[4] / 5. * T * T * T* T + Cp_coef[5] / T;
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


void allocate_memory() {
    chem.T_comon = new double [num_gas_species];
    chem.products = new double* [num_react];
    chem.Arrh_params = new double* [num_react];
    chem.isReversible = new bool[num_react];
    chem.has_Third = new bool[num_react];
    chem.ThirdBodies = new double* [num_react];
    chem.has_low = new bool[num_react];
    chem.Arrh_LP_params = new double* [num_react];
    chem.has_Troe = new bool[num_react];
    chem.Troe_params = new double* [num_react];
    chem.M_exist = new bool[num_react];

    for (int i = 0; i < num_react; i++) {
        chem.products[i] = new double[num_gas_species];
        chem.Arrh_params[i] = new double[3];
        chem.ThirdBodies[i] = new double[num_gas_species];
        chem.Arrh_LP_params[i] = new double[3];
        chem.Troe_params[i] = new double[4];
    }

    Ystart = new double[num_gas_species];
    Yend = new double[num_gas_species];
    X = new double[num_gas_species];
    for (int i = 0; i < num_gas_species; i++) {
        Ystart[i] = 0;
        Yend[i] = 0;
    };

    Yi = new double[num_gas_species];
    Yiprev = new double[num_gas_species];
    Yinext = new double[num_gas_species];
    YkVk_r = new double[num_gas_species];
    YkVk_l = new double[num_gas_species];

    gradX_r = new double[num_gas_species];
    gradX_l = new double[num_gas_species];

    X_tmp_r = new double[num_gas_species];
    X_tmp_l = new double[num_gas_species];

    Y_tmp_r = new double[num_gas_species];
    Y_tmp_l = new double[num_gas_species];

    Xiprev = new double[num_gas_species];
    Xi = new double[num_gas_species];
    Xinext = new double[num_gas_species];
    gradX = new double[num_gas_species];
    Y_tmp = new double[num_gas_species];
    X_tmp = new double[num_gas_species];
    YkVk = new double[num_gas_species];

    Sn = new double[9];
    Hn = new double[9];
    Cpn = new double[9];

    ydot = new double[num_gas_species];
    ydot = new double[num_gas_species];
    forward_arr = new double[num_react];
    reverse_arr = new double[num_react];
    equilib_arr = new double[num_react];
    Y_left_bound = new double[num_gas_species];
    wk_add = new double[num_gas_species];

    Dij_res = new double* [num_gas_species];
    phyc.Cp_coef_hT = new double* [num_gas_species];
    phyc.Cp_coef_lT = new double* [num_gas_species];

    for (int i = 0; i < num_gas_species; i++) {
        Dij_res[i] = new double [num_gas_species];
        phyc.Cp_coef_hT[i] = new double[9];
        phyc.Cp_coef_lT[i] = new double[9];
    }

    diff_polynom = new double** [num_gas_species];
    lambda_polynom = new double* [num_gas_species];
    for (int i = 0; i < num_gas_species; i++) {
        lambda_polynom[i] = new double[5];
        diff_polynom[i] = new double* [num_gas_species];
    }
    for (int i = 0; i < num_gas_species; i++) {
        for (int j = 0; j < num_gas_species; j++) {
            diff_polynom[i][j] = new double[5];
        }
    }
}


void free_memory() {
    delete[] chem.isReversible ;
    delete[] chem.has_Third;
    delete[] chem.has_low ;
    delete[] chem.has_Troe;
    delete[] chem.M_exist;

    for (int i = 0; i < num_react; i++) {
        delete[] chem.products[i];
        delete[] chem.Arrh_params[i];
        delete[] chem.ThirdBodies[i];
        delete[]chem.Arrh_LP_params[i] ;
        delete[] chem.Troe_params[i] ;
    }
    delete[] chem.products;
    delete[] chem.Arrh_params;
    delete[] chem.ThirdBodies ;
    delete[] chem.Arrh_LP_params ;
    delete[] chem.Troe_params;

    delete[] Ystart;
    delete[] Yend;
    delete[] X;

    delete[] Yi;
    delete[] Yiprev;
    delete[] Yinext;
    delete[] YkVk_r;
    delete[] YkVk_l;

    delete[] gradX_r ;
    delete[] gradX_l;

    delete[]  X_tmp_r ;
    delete[] X_tmp_l;

    delete[] Y_tmp_r ;
    delete[] Y_tmp_l ;

    delete[] Xiprev;
    delete[] Xi ;
    delete[] Xinext ;
    delete[] gradX;
    delete[] Y_tmp ;
    delete[] X_tmp ;
    delete[] YkVk ;

    delete[] Sn;
    delete[] Hn ;
    delete[] Cpn;

    delete[] ydot ;
    delete[] forward_arr;
    delete[] reverse_arr;
    delete[] equilib_arr;
    delete[] Y_left_bound;
    delete[] wk_add ;

    for (int i = 0; i < num_gas_species; i++) {
        delete[] phyc.Cp_coef_hT[i];
        delete[] phyc.Cp_coef_lT[i];
    }
    delete[] phyc.Cp_coef_hT;
    delete[] phyc.Cp_coef_lT;


    for (int i = 0; i < num_gas_species; i++) {
        delete[] lambda_polynom[i];
    }
    for (int i = 0; i < num_gas_species; i++) {
        for (int j = 0; j < num_gas_species; j++) {
            delete[] diff_polynom[i][j];
        }
    }
    for (int i = 0; i < num_gas_species; i++) {
        delete[] diff_polynom[i];
    }
    delete[] diff_polynom;
}
