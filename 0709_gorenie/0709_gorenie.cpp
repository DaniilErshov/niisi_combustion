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

double eps_x = pow(10, -6);
double eps_fr = pow(10, -6);
const double kB = 1.3806504e-23;
const double Angstroem__ = 1.0e-10;
const double santimetr = 1.0e-8;
const vector<double> M = { 2., 1.0, 32.0, 16.0, 17.0, 33.0, 18.0, 34.0, 28.0 };
string name_species[9] = { "H2", "H", "O2",
  "O", "OH", "HO2", "H2O", "H2O2", "N2" };

std::map<std::string, int> komponents{
    {"H2", 0},
    {"H", 1},
    {"O2", 2},
    {"O", 3},
    {"OH", 4},
    {"HO2", 5},
    {"H2O", 6},
    {"H2O2", 7},
    {"N2", 8}
};
std::map<int, string> komponents_str{
    {0, "H2"},
    {1, "H"},
    {2, "O2"},
    {3, "O"},
    {4, "OH"},
    {5, "HO2"},
    {6, "H2O"},
    {7, "H2O2"},
    {8, "N2"}
};
std::map<int, std::map<string, double>> Dij_saved;

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
    //T_find();
    int zh;
    int N_x = 19;

    vector<double> x_vect(N_x);
    vector<double> Y_vect(N_x * num_gas_species);
    vector<double> T_vect(N_x);

    double* Ystart = new double[num_gas_species];
    double* Yend = new double[num_gas_species];
    double* X = new double[num_gas_species];
    for (int i = 0; i < num_gas_species; i++) {
        Ystart[i] = 0;
        Yend[i] = 0;
    }

    double T_start = Tstart;
    double T_finish = Tfinish;
    double T_center;
    std::cout << "Chemkin Reader Test\n";

    double koeff_topl = 1.0;

    makeYstart(koeff_topl, Ystart);

    Find_final_state_IDA(T_start, Tend, Ystart, Yend);

    cout << "Tstart = " << T_start << "\n";
    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        cout << "Ystart = " << Ystart[k_spec] << "\n";
    }

    cout << "Tend = " << Tend << "\n";
    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        cout << "Yend = " << Yend[k_spec] << "\n";
    }

    M = 300 * get_rho(Ystart, Tstart);
    cout << "initial_M = " << M << "\n";


    int j_t = 1;

    N_center = InitialData2(N_x, x_vect, T_vect, Y_vect, M, T_start, Tend, Ystart, Yend);
    Write_to_file2("detail/initial", fout, x_vect,
        T_vect, Y_vect, Y_vect, M, N_x, 1);

    double t_Y = pow(10, -5), t_full = 0.5;
    T_center = T_vect[N_center];
    {

        integrate_Y_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, t_Y);
        Write_to_file2("detail/after_Y", fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);


        integrate_All_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 1, t_full);
        Write_to_file2("detail/Ida_1", fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);


        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.1, 1, 0, T_center);
        integrate_Y_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, t_Y);
        integrate_All_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 2, t_full);
        Write_to_file2("detail/Ida_2", fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);

        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.05, 3, 0, T_center);
        //integrate_Y_IDA(N_x, x_vect,
        //    T_vect, Y_vect, M, N_center, Ystart, t_Y);
        integrate_All_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 3, t_full);
        Write_to_file2("detail/Ida_3", fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);

        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.01, 3, 0, T_center);
        Integrate_Kinsol(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 6);
        Write_to_file2("detail/KINSOL1_" + to_string(Tstart) + "_" + to_string(koeff_topl), fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);

        fout.open("detail/M1.dat");
        fout << "M = " << M << "\n";
        fout << "rho = " << get_rho(Ystart, Tstart) << "\n";
        fout << "v = " << M / get_rho(Ystart, Tstart) << "\n";
        fout.close();


        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.001, 3 , 0, T_center);
        Integrate_Kinsol(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 6);
        Write_to_file2("detail/KINSOL2_" + to_string(Tstart) + "_" + to_string(koeff_topl), fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);
        fout.open("detail/M2.dat");
        fout << "M = " << M << "\n";
        fout << "rho = " << get_rho(Ystart, Tstart) << "\n";
        fout << "v = " << M / get_rho(Ystart, Tstart) << "\n";
        fout.close();

        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.001, 2, 0, T_center);
        Integrate_Kinsol(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 6);
        Write_to_file2("detail/KINSOL3_" + to_string(Tstart) + "_" + to_string(koeff_topl), fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);
        fout.open("detail/M3.dat");
        fout << "M = " << M << "\n";
        fout << "rho = " << get_rho(Ystart, Tstart) << "\n";
        fout << "v = " << M / get_rho(Ystart, Tstart) << "\n";
        fout.close();
    }
    return 0;
}