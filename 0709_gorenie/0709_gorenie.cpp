#include "functions_sundials.h"

double k_mol = pow(10, 3);
double P = 0.101325;
double R = 8.314;

double p_inter = 1.5;
double Tstart = 293.15;
double Tfinish = 1200;
int add_variable_count = 4;
int count_var_in_cell = 0;
int preinter = 50;
int Nx = 150 + preinter;

const double kB = 1.3806504e-23;
const double Angstroem__ = 1.0e-10;
const double santimetr = 1.0e-8;


double vel_prev;
double t_curr = 0;
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

double* YkVk_r; double* YkVk_l;

double* gradX_r; double* gradX_l;

double* X_tmp_r; double* X_tmp_l;

double* Y_tmp_r; double* Y_tmp_l;

double* Xiprev; double* Xi; double* Xinext;

double r_inter;
double* gradX;

double* Y_tmp; double* X_tmp;
double* YkVk;
double** Dij_res;

double* Sn; double* Hn; double* Cpn;
double* ydot;

double* forward_arr; double* reverse_arr; double* equilib_arr;
double* Y_left_bound;
double* wk_add;
double*** diff_polynom; double** lambda_polynom;
double* mol_weight;
double* Ystart; double* Yend;
double eps = pow(10, -8);
double* X;

double u_inter; double u_inter_2r; double u_inter_3r;


int num_gas_species; int num_react;

vector<double> x_vect, Y_vect, T_vect, u_vect, vel_vect, rho_vect;
vector<double> drhodt_vect;
vector<double> VkH2O_vect;
vector<double> VkN2_vect;
vector<double> dTdt_vect;
vector<double> dWdt_vect;


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
   
    double Tbegin = Tstart;
    double Tend = Tfinish;
    double T_center;
    int N_x = Nx;

    x_vect.resize(N_x);
    Y_vect.resize(N_x * num_gas_species); T_vect.resize(N_x); u_vect.resize(N_x);
    vel_vect.resize(N_x); rho_vect.resize(N_x); drhodt_vect.resize(N_x);
    VkH2O_vect.resize(N_x); VkN2_vect.resize(N_x);
    Cell_Properties_vector.resize(N_x); Cell_prouds_vector.resize(N_x); Cell_rval_vector.resize(N_x);
    dTdt_vect.resize(N_x); dWdt_vect.resize(N_x);
    for (int i = 0; i < N_x; i++) {
        Cell_Properties_vector[i].Y.resize(num_gas_species);
        Cell_prouds_vector[i].Y.resize(num_gas_species);
        Cell_rval_vector[i].Y.resize(num_gas_species);
    }
    Cell_Properties_inter.Y.resize(num_gas_species); Cell_prouds_inter.Y.resize(num_gas_species); Cell_rval_inter.Y.resize(num_gas_species);


    InitialData(N_x, x_vect, Cell_Properties_vector, Tbegin, Tend, Ystart, Yend);
    Write_to_file("initial", "val", Cell_Properties_vector, Cell_Properties_inter);

    double t_Y = pow(10, -7), t_full = pow(10, -10);
    KinSetIc(3 + num_gas_species + (Nx - preinter - 2));
    for (int i = 0; i < Nx; i++)
    {
        if (i <= preinter) {
            for (int k = 0; k < num_gas_species; k++) {
                Cell_Properties_vector[i].Y[k] = Cell_Properties_inter.Y[k];
            }
        }
        else {
            for (int k = 0; k < num_gas_species; k++) {
                Cell_Properties_vector[i].Y[k] = Yend[k];
            }
        }
    }

    KinSetIc(3 + num_gas_species + (Nx - preinter - 2));
    vel_prev = 0;
    ofstream setka;
    setka.open("setka.dat");
    setka << R"(VARIABLES= "r, cm", "h, cm")";
    setka << "TITLE=\"" << "Graphics" << "\"" << endl;
    for (int i = 1; i < N_x; i++)
    {
        setka << x_vect[i] << " " << x_vect[i] - x_vect[i - 1] << "\n";
    }
    setka.close();


    Write_to_file("initial2", "val", Cell_Properties_vector, Cell_Properties_inter);
    integrate_All_IDA_M(N_x);


    free_memory();
    return 0;
}