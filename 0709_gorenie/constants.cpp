#include <cmath>
#include <iostream>
#include "constants.h"
#include "functions.h"

using namespace std;

struct phy_consts phyc;
struct che_consts chec;
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

void init_consts(int& num_gas_species, int& num_react)
{
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

    phyc.kR /= phyc.kMass * pow (phyc.kVel, 2);

    int i_specie = 0;
    for (const auto& specie_i : chec.chemkinReader->species()) {
        komponents[specie_i.name()] = i_specie;
        komponents_str[i_specie] = specie_i.name();
        cout << specie_i.name() << " = " << komponents[specie_i.name()] << "\n";
        name_species.push_back(specie_i.name());
        i_specie++;
    }
    
    allocate_memory();

    set_polynom(lambda_polynom, name_file, cond_str);

    set_polynom_diffusion(diff_polynom, name_file, diff_str);


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
            cout << "lt = " << phyc.Cp_coef_lT[i][j] << "\n";
        }

        koeff_vect = chec.chemkinReader->species()[i].thermo().getUpperTemperatureCoefficients();
        for (int j = 0; j < koeff_vect.size(); j++) {
            phyc.Cp_coef_hT[i][j] = koeff_vect[j];
            cout << "ht = " << phyc.Cp_coef_hT[i][j] << "\n";
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
        cout << komponents_str[i] << " = " << phyc.mol_weight[i] << endl;
    }


    for (int i = 0; i < num_gas_species; i++)
        phyc.mol_weight[i] /= phyc.kMass;

    // 8th coef. is for enthalpy calculation

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
    Hi = Cp_coef[0] + Cp_coef[1] / 2. * pow(T, 1) + Cp_coef[2] / 3. * pow(T, 2)
        + Cp_coef[3] / 4. * pow(T, 3) + Cp_coef[4] / 5. * pow(T, 4) + Cp_coef[5] / T;
    return Hi * T;
}

double get_dSiR(double* Cp_coef, double T) {
    double Si;
    Si = Cp_coef[0] * log(T) + Cp_coef[1] * pow(T, 1) + Cp_coef[2] / 2. * pow(T, 2)
        + Cp_coef[3] / 3. * pow(T, 3) + Cp_coef[4] / 4. * pow(T, 4) + Cp_coef[6];
    return Si;
}
double get_dCpi(double* Cp_coef, double T)
{
    double Cpi;
    Cpi = Cp_coef[0] + Cp_coef[1] * pow(T, 1) + Cp_coef[2] * pow(T, 2)
        + Cp_coef[3] * pow(T, 3) + Cp_coef[4] * pow(T, 4);

    return Cpi;
}

double get_Hi(int component_i, double T)
{
    double Hi;
    int i = component_i;

    if (T > chec.chemkinReader->species()[component_i].thermo().getTCommon())
        Hi = phyc.Cp_coef_hT[i][0] + phyc.Cp_coef_hT[i][1] / 2. * pow(T, 1) + phyc.Cp_coef_hT[i][2] / 3. * pow(T, 2)
        + phyc.Cp_coef_hT[i][3] / 4. * pow(T, 3) + phyc.Cp_coef_hT[i][4] / 5. * pow(T, 4) + phyc.Cp_coef_hT[i][5] / T;
    else
        Hi = phyc.Cp_coef_lT[i][0] + phyc.Cp_coef_lT[i][1] / 2. * pow(T, 1) + phyc.Cp_coef_lT[i][2] / 3. * pow(T, 2)
        + phyc.Cp_coef_lT[i][3] / 4. * pow(T, 3) + phyc.Cp_coef_lT[i][4] / 5. * pow(T, 4) + phyc.Cp_coef_lT[i][5] / T;

    return Hi * T;
}

// Specific heat of ith component
double get_Cpi(int component_i, double T)
{
    double Cpi;
    int i = component_i;

    if (T > chec.chemkinReader->species()[component_i].thermo().getTCommon())
        Cpi = phyc.Cp_coef_hT[i][0] + phyc.Cp_coef_hT[i][1] * pow(T, 1) + phyc.Cp_coef_hT[i][2] * pow(T, 2)
        + phyc.Cp_coef_hT[i][3] * pow(T, 3) + phyc.Cp_coef_hT[i][4] * pow(T, 4);
    else
        Cpi = phyc.Cp_coef_lT[i][0] + phyc.Cp_coef_lT[i][1] * pow(T, 1) + phyc.Cp_coef_lT[i][2] * pow(T, 2)
        + phyc.Cp_coef_lT[i][3] * pow(T, 3) + phyc.Cp_coef_lT[i][4] * pow(T, 4);


    return Cpi;
}

double get_Si(int component_i, double T) {
    double Si;
    int i = component_i;

    if (T > chec.chemkinReader->species()[component_i].thermo().getTCommon())
        Si = phyc.Cp_coef_hT[i][0] * log(T) + phyc.Cp_coef_hT[i][1] * pow(T, 1) + phyc.Cp_coef_hT[i][2] / 2. * pow(T, 2)
        + phyc.Cp_coef_hT[i][3] / 3. * pow(T, 3) + phyc.Cp_coef_hT[i][4] / 4. * pow(T, 4) + phyc.Cp_coef_hT[i][6];
    else
        Si = phyc.Cp_coef_lT[i][0] * log(T) + phyc.Cp_coef_lT[i][1] * pow(T, 1) + phyc.Cp_coef_lT[i][2] / 2. * pow(T, 2)
        + phyc.Cp_coef_lT[i][3] / 3. * pow(T, 3) + phyc.Cp_coef_lT[i][4] / 4. * pow(T, 4) + phyc.Cp_coef_lT[i][6];
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
double get_Cp(int num_species, double* Y, double T)
{
    double Cp;

    Cp = 0.0;
    for (int i = 0; i < num_species; i++)
        Cp += Y[i] * get_Cpi(i, T);

    return Cp;
}

double get_Cv(int num_species, double* Y, double T)
{
    double Cv;

    Cv = 0.0;
    for (int i = 0; i < num_species; i++)
        Cv += Y[i] * get_Cvi(i, T);

    return Cv;
}

double get_Lambda(int i, double T)
{
    double lambda_arg;

    lambda_arg = lambda_polynom[i][0] + lambda_polynom[i][1] * log(T)
        + lambda_polynom[i][2] * pow(log(T), 2) + lambda_polynom[i][3] * pow(log(T), 3);

    return exp(lambda_arg);
}

void allocate_memory() {
    Ystart = new double[num_gas_species];
    Yend = new double[num_gas_species];
    X = new double[num_gas_species];
    for (int i = 0; i < num_gas_species; i++) {
        Ystart[i] = 0;
        Yend[i] = 0;
    }

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
        lambda_polynom[i] = new double[4];
        diff_polynom[i] = new double* [num_gas_species];
    }
    for (int i = 0; i < num_gas_species; i++) {
        for (int j = 0; j < num_gas_species; j++) {
            diff_polynom[i][j] = new double[4];
        }
    }
}


void free_memory() {
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
