#include "functions_sundials.h"

double k_mol = pow(10, 3);
double Y_N2 = 0.745187;
double Y_max = 1 - Y_N2;
double P = 0.101325;
double R = 8.314;
double koeff_l = 0.4;
double l = 2;
double x_center;
long int myiter = 0;
long int nniters;
double eps_func = pow(10, -8);
double* norm;
double Tstart = 400;
double Tfinish = 1500;
double nevyaz_Y;
double nevyaz_T;
double eps_x = pow(10, -6);
double eps_fr = pow(10, -6);
const double kB = 1.3806504e-23;
const double Angstroem__ = 1.0e-10;
const double santimetr = 1.0e-8;

double* Yi;
double* Yiprev;
double*  Yinext;

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

double ftol;
double stoler;
int num_gas_species;
int num_react ;
int ida_steps;
vector<double> x_vect;
vector<double> Y_vect;
vector<double> T_vect;
vector<double> u_vect;

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
    int N_x = 300;

    x_vect.resize(N_x);
    Y_vect.resize(N_x * num_gas_species);
    T_vect.resize(N_x);
    u_vect.resize(N_x);
    resize_koeff_vectors(N_x);

    makeYstart(koeff_topl, "NC7H16", 0.21, 0.79, Ystart);
    M = 60 * get_rho(Ystart, Tstart);
    int j_t = 1;
    N_center = InitialData(N_x, x_vect, T_vect, Y_vect, u_vect, M, T_start, Tend, Ystart, Yend);
    Tfinish = Tend;

    double t_Y = pow(10, -7), t_full = pow(10, -8);
    T_center = T_vect[N_center];
    ida_steps = 2000;
    Write_to_file("detail/initial", fout, x_vect,
        T_vect, Y_vect, u_vect, M, N_x, 1);
    integrate_All_IDA_M(N_x, x_vect,
        T_vect, Y_vect, u_vect, M, N_center, Ystart, 1, t_full);
    Write_to_file("detail/Ida_1", fout, x_vect,
        T_vect, Y_vect, u_vect,  M, N_x, 1);


    free_memory();
    return 0;
}