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

double Tstart = 400;
double Tfinish = 2385.4;
double nevyaz_Y;
double nevyaz_T;
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
std::unordered_map<int, std::unordered_map<string, double>> Dij_saved;

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

    for (double koeff_topl = 1.8; koeff_topl <= 2; koeff_topl += 0.1){
        int N_x = 19;
        vector<double> x_vect(N_x);
        vector<double> Y_vect(N_x * num_gas_species);
        vector<double> T_vect(N_x);
        makeYstart(koeff_topl, Ystart);
        Find_final_state_IDA(T_start, Tend, Ystart, Yend);

        M = 300 * get_rho(Ystart, Tstart);
        int j_t = 1;

        N_center = InitialData2(N_x, x_vect, T_vect, Y_vect, M, T_start, Tend, Ystart, Yend);
        Write_to_file2("detail_" + to_string(Tstart) + "/initial", fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);

        double t_Y = pow(10, -7), t_full = pow(10, -8);
        T_center = T_vect[N_center];
        nevyaz_Y = 1;
        nevyaz_T = 30;
        integrate_Y_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, t_Y);
        //Write_to_file2("detail/after_Y", fout, x_vect,
        //    T_vect, Y_vect, Y_vect, M, N_x, 1);


        integrate_All_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 1, t_full);
        //Write_to_file2("detail/Ida_1", fout, x_vect,
        //    T_vect, Y_vect, Y_vect, M, N_x, 1);


        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.0001, 1, 0, T_center);
        integrate_Y_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, t_Y);
        integrate_All_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 2, t_full);
        //Write_to_file2("detail/Ida_2", fout, x_vect,
        //    T_vect, Y_vect, Y_vect, M, N_x, 1);

        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.0005, 1, 0, T_center);
        //integrate_Y_IDA(N_x, x_vect,
        //    T_vect, Y_vect, M, N_center, Ystart, t_Y);
        Integrate_Kinsol(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 6);
        //Write_to_file2("detail/KINSOL0_" + to_string(Tstart) + "_" + to_string(koeff_topl), fout, x_vect,
        //    T_vect, Y_vect, Y_vect, M, N_x, 1);

        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.0005, 3, 0, T_center);
        cout << "NX = " << x_vect.size() << "\n";
        Integrate_Kinsol(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 6);
        //Write_to_file2("detail/KINSOL1_" + to_string(Tstart) + "_" + to_string(koeff_topl), fout, x_vect,
        //    T_vect, Y_vect, Y_vect, M, N_x, 1);

        //fout.open("detail/M1.dat");
        //fout << "M = " << M << "\n";
        //fout << "rho = " << get_rho(Ystart, Tstart) << "\n";
        //fout << "v = " << M / get_rho(Ystart, Tstart) << "\n";
        //fout.close();


        Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.007, 3 , 0, T_center);
        cout << "NX = " << x_vect.size() << "\n";
        Integrate_Kinsol(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 6);
        Write_to_file2("detail_" + to_string((int)Tstart) + "/KINSOL2_" + to_string(koeff_topl), fout, x_vect,
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
        Write_to_file2("detail_" + to_string((int)Tstart) + "/KINSOL3_" + to_string(koeff_topl), fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);
        fout.open("detail_" + to_string((int)Tstart) + "/M3_" + to_string(koeff_topl) + ".dat");
        fout << "M = " << M << "\n";
        fout << "rho = " << get_rho(Ystart, Tstart) << "\n";
        fout << "v = " << M / get_rho(Ystart, Tstart) << "\n";
        fout << "T = " << T_vect[T_vect.size() - 1] << "\n";
        fout.close();

     /*   Add_elem_simple(T_vect, Y_vect, x_vect, N_x, N_center, 0.006, 2, 0, T_center);
        cout << "NX = " << x_vect.size() << "\n";
        Integrate_Kinsol(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 6);
        Write_to_file2("detail_" + to_string((int)Tstart) + "/KINSOL4_" + to_string(koeff_topl), fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);
        fout.open("detail_" + to_string((int)Tstart) + "/M4_" + to_string(koeff_topl) + ".dat");
        fout << "M = " << M << "\n";
        fout << "rho = " << get_rho(Ystart, Tstart) << "\n";
        fout << "v = " << M / get_rho(Ystart, Tstart) << "\n";
        fout << "T = " << T_vect[T_vect.size() - 1] << "\n";
        fout.close();*/
    }
    return 0;
}