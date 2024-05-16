﻿#include "functions.h"

double k_mol = pow(10, 3);
double Y_N2 = 0.745187;
double Y_max = 1 - Y_N2;
double P = 0.101325;
double R = 8.314;
double koeff_l = 0.4;
double l = 0.1;
double x_center;
long int myiter = 0;
long int nniters;
double eps_func = pow(10, -8);
double* norm;
double Tstart = 700;
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


double ftol;
double stoler;
int num_gas_species;
int num_react ;
int ida_steps;
vector<double> x_vect;
vector<double> Y_vect;
vector<double> T_vect;

vector<vector<double>> Cp_arr;
vector<vector<double>> Lambda_arr_r;
vector<vector<double>> Lambda_arr_l;
vector<vector<vector<double>>> Dij_arr_r;
vector<vector<vector<double>>> Dij_arr_l;

struct chem_struct chem;
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
    norm = new double[1];
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
 
    //Get_molar_cons(Xi, Yi, 2282.787);
    //chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr, 2282.787, Xi, ydot);

    //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
    //    fout.open("lambdas/lambda_" + komponents_str[k_spec] + ".dat");
    //    string title2 = R"(VARIABLES= "T", "lambda")";
    //    fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    //    fout << title2 << endl;
    //    for (double T_tmp = 300; T_tmp < 3000; T_tmp += 10) {
    //        fout << T_tmp << " " << get_Lambda5(k_spec, T_tmp) << "\n";
    //    }
    //    fout.close();
    //}

    //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
    //    for (int k_spec2 = 0; k_spec2 < num_gas_species; k_spec2++) {
    //        fout.open("diffs/diff_" + komponents_str[k_spec] + "_" + komponents_str[k_spec2] + ".dat");
    //        string title2 = R"(VARIABLES= "T", "diff_)" + komponents_str[k_spec] + "_" + komponents_str[k_spec2] + R"(")";
    //        fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    //        fout << title2 << endl;
    //        for (double T_tmp = 300; T_tmp < 3000; T_tmp += 10) {
    //            fout << T_tmp << " " << Dij_func5(k_spec, k_spec2, T_tmp) << "\n";
    //        }
    //        fout.close();
    //    }
    //}
    //set_Dij_res(1807.936);
    //Get_mole_fr(Xi, Yi);
    //Get_mole_fr(Xinext, Yinext);
    //get_grad(gradX, Xi, Xinext, 0.165625, 0.1657552);
    //make_averageY(Y_tmp, Yi, Yinext);
    //for (int k = 0; k < num_gas_species; k++) {
    //  YkVk_func(k, 1807.936, Yi, gradX, Xi, Yi);
    //}
    fout_v.open("detail/v_" + to_string(Tstart) + "_.dat");
    string title2 = R"(VARIABLES= "koeff_fuel", "v", "Tend", "norm")";
    fout_v << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout_v << title2 << endl;

    for (double koeff_topl = 1.1; koeff_topl <= 1.1; koeff_topl += 0.1)
    {
        double tout1 = pow(10, -7);
        double Tend = 2385.4;
        double T_start = Tstart;
        double T_finish = Tfinish;
        double T_center;
        int N_x = 10;

        x_vect.resize(N_x);
        Y_vect.resize(N_x * num_gas_species);
        T_vect.resize(N_x);
        Cp_arr.resize(N_x);
        Lambda_arr_r.resize(N_x);
        Lambda_arr_l.resize(N_x);
        Dij_arr_r.resize(N_x);
        Dij_arr_l.resize(N_x);
        for (int i = 0; i < N_x; i++) {
            Cp_arr[i].resize(num_gas_species);
            Lambda_arr_r[i].resize(num_gas_species);
            Lambda_arr_l[i].resize(num_gas_species);
            Dij_arr_r[i].resize(num_gas_species);
            Dij_arr_l[i].resize(num_gas_species);
            for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
                Dij_arr_r[i][k_spec].resize(num_gas_species);
                Dij_arr_l[i][k_spec].resize(num_gas_species);
            }
        }
        makeYstart(koeff_topl, "NC7H16", 0.21, 0.79, Ystart);
        Find_final_state_IDA(T_start, Tend, Ystart, Yend);
        Find_final_state_KINSOL(T_start, Tend, Ystart, Yend);
        M = 100 * get_rho(Ystart, Tstart);
        int j_t = 1;
        N_center = InitialData(N_x, x_vect, T_vect, Y_vect, M, T_start, Tend, Ystart, Yend);
        Tfinish = Tend;
        //Write_to_file("detail/initial", fout, x_vect,
        //    T_vect, Y_vect, Y_vect, M, N_x, 1);

        double t_Y = pow(10, -7), t_full = pow(10, -6);
        T_center = T_vect[N_center];
        ida_steps = 5;
        eps = 0;
        integrate_All_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 1, t_full);
        //Write_to_file("detail/Ida_1", fout, x_vect,
        //    T_vect, Y_vect, Y_vect, M, N_x, 1);

        ida_steps = 60;
        integrate_All_IDA_M(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 1, t_full);
        //Write_to_file("detail/Ida_1", fout, x_vect,
        //    T_vect, Y_vect, Y_vect, M, N_x, 1);

        eps = pow(10, -3);
        Integrate_Kinsol(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 6);
        //Write_to_file("detail/KINSOL1", fout, x_vect,
        //    T_vect, Y_vect, Y_vect, M, N_x, 1);
        cout << "v = " << M / get_rho(Ystart, Tstart) << "\n";
        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.03, 2, 1, T_center);
        ida_steps = 40;
        eps = pow(10, -3);
        integrate_All_IDA_M(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 1, t_full);
        Integrate_Kinsol(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 6);
        cout << "v = " << M / get_rho(Ystart, Tstart) << "\n";
        int vect_size_temp = 0;
        int number_epoch = 0;
        int add_cell = 2;
        int add_cell_start = 2;
        ftol = pow(10, -6);
        stoler = pow(10, -20);

        //ida_steps = 60;
        //eps = pow(10, -3);
        //integrate_All_IDA_M(N_x, x_vect,
        //    T_vect, Y_vect, M, N_center, Ystart, 1, t_full);

        while (vect_size_temp < 400 || fabs(vect_size_temp - x_vect.size()) < 10)
        {
            vect_size_temp = x_vect.size();
            cout << to_string(koeff_topl) << " Nx = " << vect_size_temp << "\n";
            if (number_epoch > 8) add_cell = 0;
            if (number_epoch > 6) add_cell_start = 0;
            Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.025, add_cell, add_cell_start, T_center);

            //ida_steps = 10;
            //eps = pow(10, -3);

            //integrate_All_IDA_M(N_x, x_vect,
            //    T_vect, Y_vect, M, N_center, Ystart, 1, t_full);
            Integrate_Kinsol(N_x, x_vect,
                T_vect, Y_vect, M, N_center, Ystart, 6);
            cout << to_string(koeff_topl) << " v = " << M / get_rho(Ystart, Tstart) << "\n";
            cout << to_string(koeff_topl) << " Tend = " << T_vect[T_vect.size() - 1] << "\n";
            cout << " norm = " << norm[0] << "\n";
            /* Write_to_file("detail/KINSOL_" + to_string(koeff_topl) + "_" + to_string(x_vect.size()), fout, x_vect,
                 T_vect, Y_vect, Y_vect, M, N_x, 1);*/

             //fout.open("detail/M_.dat");
             //fout << "M = " << M << "\n";
             //fout << "rho = " << get_rho(Ystart, Tstart) << "\n";
             //fout << "v = " << M / get_rho(Ystart, Tstart) << "\n";
             //fout << "T = " << T_vect[T_vect.size() - 1] << "\n";
             //fout.close();
            number_epoch++;
            add_cell = 1;
        }

        Write_to_file("detail/KINSOL_" + to_string(koeff_topl), fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);

        fout_v << koeff_topl << " " << M / get_rho(Ystart, Tstart) << " " << T_vect[T_vect.size() - 1] << " "
            << norm[0] << "\n";
        cout << "fout_v did\n";
    }
    fout_v.close();

    free_memory();
    return 0;
}