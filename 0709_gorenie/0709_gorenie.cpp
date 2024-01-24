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
    int number = 4000;
    int print_value = 13000;
    int cons_flag = 0;
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

    Ystart[H2] = 0;
    Ystart[H] = 0;
    Ystart[O2] = 0;
    Ystart[O] = 0;
    Ystart[OH] = 0;
    Ystart[HO2] = 0;
    Ystart[H2O] = 0;
    Ystart[H2O2] = 0;
    Ystart[N2] = 0;
    Yend[H2] = 1.4563E-002;
    Yend[H] = 1.7917E-003;
    Yend[O2] = 5.4673E-003;
    Yend[O] = 5.9257E-004;
    Yend[OH] = 7.8587E-003;
    Yend[HO2] = 1.2247E-006;
    Yend[H2O] = 3.2394E-001;
    Yend[H2O2] = 1.3135E-007;
    Yend[N2] = 6.4579E-001;


    double koeff_topl = 0.5;

    makeYstart(koeff_topl, Ystart);
    M = 250 * 8.4955E-004;
    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        cout << "Ystart = " << Ystart[k_spec] << "\n";
    }

    Find_final_state_IDA(T_start, Tend, Ystart, Yend);
    cout << "Tend = " << Tend << "\n";

    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        cout << "Yend = " << Yend[k_spec] << "\n";
    }
    cout << "rho0 = " << get_rho(Ystart, Tstart);

    N_center = InitialData2(N_x, x_vect, T_vect, Y_vect, M, Tstart, Tend, Ystart, Yend);
    cout << "N_center = " << N_center << "\n";
    cout << "T N_center = " << T_vect[N_center] << "\n";

    {
        integrate_Y_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart);

        Write_to_file2("A_initial", fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);

        T_center = T_vect[N_center];

        integrate_All_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 1);

        //Integrate_Kinsol(N_x, x_vect,
        //    T_vect, Y_vect, M, N_center, Ystart, 5);

        Write_to_file2("Ida_1", fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);


        cout << "size1 = " << T_vect.size() << "\n";
        Add_elem(T_vect, Y_vect, x_vect, N_x, N_center, 0.01, 3, 0, T_center);
        Write_to_file2("A_add", fout, x_vect,
            T_vect, Y_vect, Y_vect, M, N_x, 1);

        integrate_Y_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart);
        integrate_All_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 2);



        Add_elem(T_vect, Y_vect, x_vect, N_x, N_center, 0.005, 3, 0, T_center);
        integrate_Y_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart);
        integrate_All_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 3);


        Add_elem(T_vect, Y_vect, x_vect, N_x, N_center, 0.001, 5, 0, T_center);
        integrate_All_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 4);

        Add_elem(T_vect, Y_vect, x_vect, N_x, N_center, 0.0001, 3, 0, T_center);
        integrate_Y_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart);
        integrate_All_IDA(N_x, x_vect,
            T_vect, Y_vect, M, N_center, Ystart, 5);
    }
    return 0;
}