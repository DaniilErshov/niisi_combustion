#include "functions_sundials.h"

double k_mol = pow(10, 3);
double Y_N2 = 0.745187;
double Y_max = 1 - Y_N2;
double P = 0.101325;
double R = 8.314;
double koeff_l = 0.4;
double l = 10;
double x_center;
long int myiter = 0;
long int nniters;
double eps_func = pow(10, -8);
double* norm;
double Tstart = 280;
double Tfinish = 1000;
double nevyaz_Y;
double nevyaz_T;
double eps_x = pow(10, -6);
double eps_fr = pow(10, -6);
int add_variable_count = 4;
int count_var_in_cell = 0;
int preinter = 50;
int Nx = 150 + preinter;
const double kB = 1.3806504e-23;
const double Angstroem__ = 1.0e-10;
const double santimetr = 1.0e-8;
double p_inter = 1.0012;
double vel_prev;
string Fuel = "NC7H16";

double* Xi_2;
double* Xi_3;
double* X_inter;

double* Yi;
double* Yiprev;
double* Yinext;
double* Y_inter;
double* Y_inter_2r;
double* Y_inter_3r;

double* YkVk_r;
double* YkVk_l;

double* gradX_r;
double* gradX_l;

double* X_tmp_r;
double* X_tmp_l;

double* Y_tmp_r;
double* Y_tmp_l;

double* Xiprev;
double* Xi;
double* Xinext;

double r_inter;

double* gradX;

double* Y_tmp;
double* X_tmp;
double* YkVk;
double** Dij_res;

double* Sn;
double* Hn;
double* Cpn;
double* ydot;

double* forward_arr;
double* reverse_arr;
double* equilib_arr;
double* Y_left_bound;
double* wk_add;
double*** diff_polynom;
double** lambda_polynom;
double* mol_weight;
double* Ystart;
double* Yend;
double eps = pow(10, -8);
double* X;
int jactimes = -1;
bool flag_use_save_koeffs = 0;
bool save_chem_koeffs = 0;
bool update_koeffs = 1;

double u_inter;
double u_inter_2r;
double u_inter_3r;


double ftol;
double stoler;
int num_gas_species;
int num_react ;
int ida_steps;

vector<double> x_vect;
vector<double> Y_vect;
vector<double> T_vect;
vector<double> u_vect;
vector<double> vel_vect;
vector<double> rho_vect;
vector<double> drhodt_vect;
vector<double> VkH2O_vect;
vector<double> VkN2_vect;
vector<double> dTdt_vect;
vector<double> dWdt_vect;

vector<vector<double>> Cp_arr;
vector<vector<double>> H_arr;
vector<vector<double>> Lambda_arr;
vector<vector<double>> Lambda_arr_r;
vector<vector<double>> Lambda_arr_l;
vector<vector<vector<double>>> Dij_arr;
vector<vector<vector<double>>> Dij_arr_r;
vector<vector<vector<double>>> Dij_arr_l;
vector<vector<double>> forward_arr_save;
vector<vector<double>> reverse_arr_save;
vector<string> name_species;

double corrector = 1;
std::unordered_map<std::string, int> komponents{
};

std::unordered_map<int, string> komponents_str{
};

map<string, double> elem_mol_weight{
    {"H", 1.00797},
    {"O", 15.9994},
    {"N", 14.0067},
    {"C", 12.0096},
    {"AR", 39.94},
    {"HE", 4.002602}
};

Cell_Properties Cell_Properties_inter;
Cell_Properties Cell_prouds_inter;
Cell_Properties Cell_rval_inter;

vector<Cell_Properties> Cell_Properties_vector;
vector<Cell_Properties> Cell_prouds_vector;
vector<Cell_Properties> Cell_rval_vector;

int main()
{
    init_consts(num_gas_species, num_react);
    norm = new double;
    double b = 0.1;
    double M;
    double W, rho, Y_H2, Y_O2;
    int N_center;
    int retval;
    double w_dot;
    double* my_x;
    ofstream fout_v; 
    ofstream fout;
    double X_H2O, X_H2, X_O2, X_N2;

    fout_v.open("detail/v_" + to_string(Tstart) + "_.dat");
    string title2 = R"(VARIABLES= "koeff_fuel", "v", "Tend", "norm")";
    fout_v << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout_v << title2 << endl;
    double koeff_topl = 1;
    double tout1 = pow(10, -7);
    double Tend = Tfinish;
    double T_start = Tstart;
    double T_finish = Tfinish;
    double T_center;
    int N_x = Nx;
    cout << "preinter = " << preinter << "\n";

   /* std::ifstream inFile(R"(C:\nisi_serg\var_files\r.dat)");
    N_x = std::count(std::istreambuf_iterator<char>(inFile),
        std::istreambuf_iterator<char>(), '\n') - 1;*/

    x_vect.resize(N_x);
    Y_vect.resize(N_x * num_gas_species); T_vect.resize(N_x); u_vect.resize(N_x);
    vel_vect.resize(N_x); rho_vect.resize(N_x); drhodt_vect.resize(N_x);
    VkH2O_vect.resize(N_x); VkN2_vect.resize(N_x);
    Cell_Properties_vector.resize(N_x); Cell_prouds_vector.resize(N_x); Cell_rval_vector.resize(N_x);
    resize_koeff_vectors(N_x);
    dTdt_vect.resize(N_x);
    dWdt_vect.resize(N_x);
    for (int i = 0; i < N_x; i++) {
        Cell_Properties_vector[i].Y.resize(num_gas_species);
        Cell_prouds_vector[i].Y.resize(num_gas_species);
        Cell_rval_vector[i].Y.resize(num_gas_species);
    }
    Cell_Properties_inter.Y.resize(num_gas_species);
    Cell_prouds_inter.Y.resize(num_gas_species);
    Cell_rval_inter.Y.resize(num_gas_species);

    Write_to_file("initial", "val", Cell_Properties_vector, Cell_Properties_inter);
    M = 0;
    int j_t = 1;
    double vel_interf = 0;
    N_center = InitialData(N_x, x_vect, Cell_Properties_vector, T_start, Tend, Ystart, Yend);
    //KinSet_T_inter(1);
    //preinter -= 15;
    //vector<string> params_vect = { "r", "T", "rho", "u", "vel" };
    //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
    //    params_vect.push_back(komponents_str[k_spec]);
    //}
    //for (const auto& str_i : params_vect) {
    //    std::ifstream in(R"(C:\nisi_serg\var_files\)" + str_i + R"(.dat)"); // окрываем файл для чтения
    //    string line;
    //    if (in.is_open())
    //    {
    //        int i = 0;
    //        cout << str_i << "\n";
    //        while (std::getline(in, line))
    //        {
    //            if (str_i == "r") {
    //                x_vect[i] = std::stod(line);
    //            }
    //            else if (str_i == "rho") {
    //                Cell_Properties_vector[i].rho = std::stod(line);
    //            }
    //            else if (str_i == "T") {
    //                Cell_Properties_vector[i].T = std::stod(line);
    //            }
    //            else if (str_i == "u") {
    //                Cell_Properties_vector[i].u = std::stod(line);
    //            }
    //            else if (str_i == "vel") {
    //                Cell_Properties_vector[i].vel = std::stod(line);
    //            }
    //            else {
    //                Cell_Properties_vector[i].Y[komponents[str_i]] = std::stod(line);
    //            }

    //            if (i == preinter + 1) {
    //                if (str_i == "r") {
    //                    double h = x_vect[preinter] - x_vect[preinter - 1];
    //                    double ri = std::stod(line);
    //                    p_inter = (ri - x_vect[preinter]) / h + 1.0;
    //                }
    //                else if (str_i == "rho") {
    //                    Cell_Properties_inter.rho = std::stod(line);
    //                }
    //                else if (str_i == "T") {
    //                    Cell_Properties_inter.T = std::stod(line);
    //                }
    //                else if (str_i == "u") {
    //                    Cell_Properties_inter.u = std::stod(line);
    //                }
    //                else if (str_i == "vel") {
    //                    Cell_Properties_inter.vel = std::stod(line);
    //                }
    //                else {
    //                    Cell_Properties_inter.Y[komponents[str_i]] = std::stod(line);
    //                }

    //                std::getline(in, line);

    //                if (str_i == "r") {
    //                    x_vect[i] = std::stod(line);
    //                }
    //                else if (str_i == "rho") {
    //                    Cell_Properties_vector[i].rho = std::stod(line);
    //                }
    //                else if (str_i == "T") {
    //                    Cell_Properties_vector[i].T = std::stod(line);
    //                }
    //                else if (str_i == "u") {
    //                    Cell_Properties_vector[i].u = std::stod(line);
    //                }
    //                else if (str_i == "vel") {
    //                    Cell_Properties_vector[i].vel = std::stod(line);
    //                }
    //                else {
    //                    Cell_Properties_vector[i].Y[komponents[str_i]] = std::stod(line);
    //                }
    //            }
    //            i++;
    //        }
    //    }
    //    in.close(); 
    //}

    Write_to_file("initial", "val", Cell_Properties_vector, Cell_Properties_inter);
    Tfinish = Tend;

    double t_Y = pow(10, -7), t_full = pow(10, -10);
    vel_interf = 0;

    KinSetIc(3 + num_gas_species + (Nx - preinter - 2));
    for (int i = 0; i < Nx; i++)
    {
        if (i <= preinter + 4) {
            for (int k = 0; k < num_gas_species; k++) {
                Cell_Properties_vector[i].Y[k] = Cell_Properties_inter.Y[k];
            }
        }
        else {
            for (int k = 0; k < num_gas_species; k++) {
                Cell_Properties_vector[i].Y[k] = Yend[k];
            }
        }
        Cell_Properties_vector[i].rho = get_rho(Cell_Properties_vector[i].Y.data(), Cell_Properties_vector[i].T, 'g');
    }
    KinSetIc(3 + num_gas_species + (Nx - preinter - 2));
    vel_prev = 10;
    ofstream setka;
    setka.open("setka.dat");
    setka << R"(VARIABLES= "r, cm", "h, cm")";

    setka << "TITLE=\"" << "Graphics" << "\"" << endl;
    for (int i = 1; i < N_x; i++)
    {
        setka << x_vect[i] << " " << x_vect[i] - x_vect[i - 1] << "\n";
    }
    setka.close();
    N_center = 0;
    Write_to_file("initial2", "val", Cell_Properties_vector, Cell_Properties_inter);
    integrate_All_IDA_M(N_x, x_vect,
        T_vect, Y_vect, u_vect, vel_interf, N_center, Ystart, 1, t_full);


    free_memory();
    return 0;
}