#include <cmath>
#include <iostream>
#include "functions_sundials.h"

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
extern double* Y_inter;
extern double* Y_inter_2r;
extern double* Y_inter_3r;


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
    std::string cond_str = "CONDUCTIVITIES";
    std::string visc_str = "VISCOSITIES";
    std::string diff_str = "DIFFUSION";

    //const std::string thermfile = R"(C:\nisi_serg\info\term.dat)";
    //const std::string transfile = R"(C:\nisi_serg\info\tran.dat)";
    //const std::string chemfile = R"(C:\nisi_serg\info\heptane.inp)";
    const std::string thermfile = R"(D:\Storage\Daniil\n-heptane\term.dat)";
    const std::string transfile = R"(D:\Storage\Daniil\n-heptane\tran.dat)";
    const std::string chemfile = R"(D:\Storage\Daniil\n-heptane\heptane.inp)";
    chec.chemkinReader = new IO::ChemkinReader(chemfile, thermfile, transfile);
    chec.chemkinReader = new IO::ChemkinReader(chemfile, thermfile, transfile);


    chec.chemkinReader->read();
    //chec.chemkinReader->check();
    //cout << "37\n";
    //cout << chec.chemkinReader->reactions()[37] << "\n";
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
    auto& species = chec.chemkinReader->species();


    for (int i = 0; i < num_gas_species; i++) {
        auto koeff_vect = chec.chemkinReader->species()[i].thermo().getLowerTemperatureCoefficients();
        cout << komponents_str[i] << "\n";
        for (int j = 0; j < 9; j++) {
            phyc.Cp_coef_lT[i][j] = 0;
            phyc.Cp_coef_hT[i][j] = 0;
        }

        for (int j = 0; j < koeff_vect.size(); j++) {
            phyc.Cp_coef_lT[i][j] = koeff_vect[j];
        }

        koeff_vect = chec.chemkinReader->species()[i].thermo().getUpperTemperatureCoefficients();
        for (int j = 0; j < koeff_vect.size(); j++) {
            phyc.Cp_coef_hT[i][j] = koeff_vect[j];
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
void set_polynom(double** ploynom, std::string name_file, std::string type_polynom) {
    std::string line;
    std::string koeff_str = "COEFFICIENTS";
    std::ifstream in(name_file); // îêðûâàåì ôàéë äëÿ ÷òåíèÿ
    if (in.is_open())
    {
        while (getline(in, line))
        {
            vector splitted_string = splitString(line, ' ');
            if (std::find(splitted_string.begin(), splitted_string.end(), koeff_str) != splitted_string.end()
                && std::find(splitted_string.begin(), splitted_string.end(), type_polynom) != splitted_string.end())
            {
                getline(in, line);
                while (getline(in, line))
                {
                    vector splitted_string = splitString(line, ' ');
                    if (std::find(splitted_string.begin(), splitted_string.end(), koeff_str) != splitted_string.end())
                        return;


                    if (splitted_string.size() != 0)
                    {

                        string specie = splitted_string[0];
                        std::transform(specie.begin(), specie.end(), specie.begin(), ::toupper);
                        if (komponents.contains(specie))
                        {
                            int i_komp = komponents[specie];
                            cout << "line = " << line << "\n";
                            cout << "i_komp = " << i_komp << "\n";

                            ploynom[i_komp][0] = stod(splitted_string[1]);
                            ploynom[i_komp][1] = stod(splitted_string[2]);
                            ploynom[i_komp][2] = stod(splitted_string[3]);
                            ploynom[i_komp][3] = stod(splitted_string[4]);
                            cout << specie << " " << i_komp << " " << ploynom[i_komp][0]
                                << " " << ploynom[i_komp][1]
                                << " " << ploynom[i_komp][2]
                                << " " << ploynom[i_komp][3] << "\n";
                        }

                    }

                }
            }


        }
    }
    in.close();
}

void set_polynom_diffusion(double*** polynom, std::string name_file, std::string type_polynom) {
    std::string line;
    std::string koeff_str = "COEFFICIENTS";
    std::ifstream in(name_file); // îêðûâàåì ôàéë äëÿ ÷òåíèÿ
    if (in.is_open())
    {
        while (getline(in, line))
        {
            vector splitted_string = splitString(line, ' ');
            if (std::find(splitted_string.begin(), splitted_string.end(), koeff_str) != splitted_string.end()
                && std::find(splitted_string.begin(), splitted_string.end(), type_polynom) != splitted_string.end())
            {
                getline(in, line);
                while (getline(in, line))
                {
                    vector splitted_string = splitString(line, ' ');
                    if (std::find(splitted_string.begin(), splitted_string.end(), koeff_str) != splitted_string.end())
                        return;

                    if (splitted_string.size() != 0)
                    {
                        string specie1 = splitted_string[0];
                        std::transform(specie1.begin(), specie1.end(), specie1.begin(), ::toupper);

                        string specie2 = splitted_string[1];
                        std::transform(specie2.begin(), specie2.end(), specie2.begin(), ::toupper);

                        if (komponents.contains(specie1) && komponents.contains(specie2))
                        {
                            int i_sp1 = komponents[specie1];
                            int i_sp2 = komponents[specie2];
                            polynom[i_sp1][i_sp2][0] = stod(splitted_string[2]);
                            polynom[i_sp1][i_sp2][1] = stod(splitted_string[3]);
                            polynom[i_sp1][i_sp2][2] = stod(splitted_string[4]);
                            polynom[i_sp1][i_sp2][3] = stod(splitted_string[5]);

                            polynom[i_sp2][i_sp1][0] = stod(splitted_string[2]);
                            polynom[i_sp2][i_sp1][1] = stod(splitted_string[3]);
                            polynom[i_sp2][i_sp1][2] = stod(splitted_string[4]);
                            polynom[i_sp2][i_sp1][3] = stod(splitted_string[5]);

                            /*cout << komponents[specie1] << " " << komponents[specie2]
                                << " " << polynom[i_sp1][i_sp2][0]
                                << " " << polynom[i_sp1][i_sp2][1]
                                << " " << polynom[i_sp1][i_sp2][2]
                                << " " << polynom[i_sp1][i_sp2][3] << "\n";*/
                        }

                    }

                }
            }
        }
    }
    in.close();
}
std::vector<std::string> splitString(std::string str, char splitter) {
    std::vector<std::string> result;
    std::string current = "";
    for (int i = 0; i < str.size(); i++) {
        if (str[i] == splitter) {
            if (current != "") {
                result.push_back(current);
                current = "";
            }
            continue;
        }
        current += str[i];
    }
    if (current.size() != 0)
        result.push_back(current);
    return result;
}

template <typename T>
void findValue(const std::vector<T>& data, bool(*condition)(T))
{
    auto result{ std::find_if(begin(data), end(data), condition) };
    if (result == end(data))
        std::cout << "Value not found" << std::endl;
    else
        std::cout << "Value found at position " << (result - begin(data)) << std::endl;
}




void allocate_memory() {
    Ystart = new double[num_gas_species];
    Yend = new double[num_gas_species];

    Xi_2 = new double[num_gas_species];
    Xi_3 = new double[num_gas_species];
    X_inter = new double[num_gas_species];
    Y_inter_2r = new double[num_gas_species];
    Y_inter_3r = new double[num_gas_species];
    X = new double[num_gas_species];
    for (int i = 0; i < num_gas_species; i++) {
        Ystart[i] = 0;
        Yend[i] = 0;
    };
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

    delete[] Ystart;
    delete[] Yend;
    delete[] X;

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
