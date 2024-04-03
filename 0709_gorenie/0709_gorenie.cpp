#include "functions.h"

double k_mol = pow(10, 3);
double Y_N2 = 0.745187;
double Y_max = 1 - Y_N2;
double P = 0.101325;
double R = 8.314;
double koeff_l = 0.4;
double l = 0.3;
double x_center;
long int myiter = 0;
long int nniters;
double eps_func = pow(10, -8);

double Tstart = 300;
double Tfinish = 2385.4;
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
double* X;
int num_gas_species;
int num_react ;
int ida_steps;
vector<double> x_vect;
vector<double> Y_vect;
vector<double> T_vect;

vector<string> name_species;

std::map<std::string, int> komponents{
};

std::map<int, string> komponents_str{
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
    double b = 0.1;
    double M;
    double W, rho, Y_H2, Y_O2;
    int N_center;
    int retval;
    double w_dot;
    double* my_x;
    ofstream fout;
    double X_H2O, X_H2, X_O2, X_N2;
    double tout1 = pow(10, -7);
    double Tend = 2385.4;

    double T_start = Tstart;
    double T_finish = Tfinish;
    double T_center;
    std::cout << "Chemkin Reader Test\n";

    for (double koeff_topl = 1.0; koeff_topl <= 1.0; koeff_topl += 0.1){
        int N_x = 19;
        x_vect.resize(N_x);
        Y_vect.resize(N_x * num_gas_species);
        T_vect.resize(N_x);
        makeYstart(koeff_topl, Ystart);
        Find_final_state_IDA(T_start, Tend, Ystart, Yend);
        M = 300 * get_rho(Ystart, Tstart);
        int j_t = 1;

        N_center = InitialData(N_x, x_vect, T_vect, Y_vect, M, T_start, Tend, Ystart, Yend);
        Write_to_file("detail_" + to_string(Tstart) + "/initial", fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);

        double t_Y = pow(10, -7), t_full = pow(10, -6);
        T_center = T_vect[N_center];

        ida_steps = 100;
        integrate_Y_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, t_Y);
        Write_to_file("detail/after_Y", fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);

        ida_steps = 30;
        integrate_All_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 1, t_full);
        Write_to_file("detail/Ida_1", fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);


        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.0001, 1, 0, T_center);
        integrate_Y_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, t_Y);

        ida_steps = 30;
        integrate_All_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 2, t_full);
        Write_to_file("detail/Ida_2", fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);

        ida_steps = 30;
        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.0005, 1, 0, T_center);
        integrate_All_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 2, t_full);
        Write_to_file("detail/Ida_3", fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);

        Integrate_Kinsol(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 6);
       
        Write_to_file("detail/KINSOL1_" + to_string(koeff_topl), fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);


        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.0005, 3, 0, T_center);
        cout << "NX = " << x_vect.size() << "\n";
        Integrate_Kinsol(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 6);
        Write_to_file("detail/KINSOL2_" + to_string(koeff_topl), fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);

        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.007, 3 , 0, T_center);
        cout << "NX = " << x_vect.size() << "\n";
        Integrate_Kinsol(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 6);
        Write_to_file("detail/KINSOL3_" + to_string(koeff_topl), fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);
        fout.open("detail_" + to_string((int)Tstart) + "/M2_" + to_string(koeff_topl) + ".dat");
        fout << "M = " << M << "\n";
        fout << "rho = " << get_rho(Ystart, Tstart) << "\n";
        fout << "v = " << M / get_rho(Ystart, Tstart) << "\n";
        fout << "T = " << T_vect[T_vect.size() - 1] << "\n";
        fout.close();

        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.001, 1, 0, T_center);
        cout << "NX = " << x_vect.size() << "\n";
        Integrate_Kinsol(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 6);
        Write_to_file("detail/KINSOL4_" + to_string(koeff_topl), fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);
        fout.open("detail_" + to_string((int)Tstart) + "/M3_" + to_string(koeff_topl) + ".dat");
        fout << "M = " << M << "\n";
        fout << "rho = " << get_rho(Ystart, Tstart) << "\n";
        fout << "v = " << M / get_rho(Ystart, Tstart) << "\n";
        fout << "T = " << T_vect[T_vect.size() - 1] << "\n";
        fout.close();

        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.00007, 3, 0, T_center);
        cout << "NX = " << x_vect.size() << "\n";
        Integrate_Kinsol(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 6);
        Write_to_file("detail/KINSOL4_" + to_string(koeff_topl), fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);
    } 
    free_memory();
    return 0;
}