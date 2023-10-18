#include "functions.h"

double k_mol = pow(10, 3);
double Y_N2 = 0.745187;
double Y_max = 1 - Y_N2;
double P = 0.101325;
double R = 8.314;
double koeff_l = 0.35;
double l = 0.8;
long int myiter = 0;
long int nniters;
double Tstart = 300;
double Tfinish = 2200;
double eps_x = pow(10, -6);
double eps_fr = pow(10, -6);
const double kB = 1.3806504e-23;
const double Angstroem__ = 1.0e-10;
const double santimetr = 1.0e-8;
const vector<double> M = { 2., 1.0, 32.0, 16.0, 17.0, 33.0, 18.0, 34.0, 28.0 };

int main()
{
    init_consts(num_gas_species, num_react);
    int N_x = 15;
    double b = 0.01;
    double M;
    double W, rho, Y_H2, Y_O2;
    int N_center;
    int retval;
    double w_dot;
    double* my_x;
    ofstream fout;
    double X_H2O, X_H2, X_O2, X_N2;
    double tout1 = pow(10, -7);
    double Tend = 2200;
    int number = 4000;
    int print_value = 13000;
    int cons_flag = 0;
    //T_find();
    int zh;
    vector<double> x_vect(N_x);
    vector<double> T_vect(N_x);
    vector<double> Y_vect(N_x * 9);
    double* Ystart = new double[num_gas_species];
    double* Yend = new double[num_gas_species];


    const std::string chemfile = R"(C:\Users\Данило\source\repos\check_boost3\ChemKin_reader\test\2_easychem.inp)";
    const std::string thermfile = R"(C:\Users\Данило\source\repos\check_boost3\ChemKin_reader\test\therm-abf.dat)";
    const std::string transfile = R"(C:\Users\Данило\source\repos\check_boost3\ChemKin_reader\test\tran.dat)";
    IO::ChemkinReader chemkinReader(chemfile, thermfile, transfile);
    chemkinReader.read();
    chemkinReader.check();

    std::cout << "NEW BLOCK" << std::endl;
    
    for (int i = 0; i < 9; i++)
    {
        cout << "i = " << i << endl;
        std::cout << chemkinReader.species()[i].name() << std::endl;
        //std::cout << chemkinReader.species()[i].thermo().getPhase() << std::endl;
        //std::cout << chemkinReader.species()[i].transport().getCollisionDiameter()/Angstroem__ << std::endl;
        //std::cout << chemkinReader.species()[i].transport().getPotentialWellDepth() / kB << std::endl;
    }
    
    for (int i = 0; i < num_gas_species; i++) {
        Ystart[i] = 0;
        Yend[i] = 0;
    }
    double T_start = Tstart;
    double T_finish = Tfinish;

    Ystart[8] = Yend[8] = Y_N2;
    Yend[6] = 1 - Y_N2;

    Ystart[0] = (1 - Y_N2) * 1. / 9.;
    Ystart[2] = (1 - Y_N2) * 8. / 9.;
    N_center = InitialData(N_x, x_vect, T_vect, Y_vect, M, T_start, T_finish, Ystart, Yend);
    N_x = x_vect.size();

    for (int i = 0; i < N_x * num_gas_species; i++) {
        cout << "Yi = " << Y_vect[i] << endl;
        if ((i + 1) % num_gas_species == 0) cout << endl;
    }

    cout << "after INITIAL DATA N_x = " << N_x << "\n";
    for (int i = 0; i < N_x - 1; i++) {
        if (x_vect[i] <= koeff_l * l && x_vect[i + 1] > koeff_l * l)
            N_center = i;
    }
    Write_to_file2("detail", fout, x_vect,
        T_vect, Y_vect, M, N_x, 1);
    std::cout << "Chemkin Reader Test\n";

    Find_final_state_IDA(&chemkinReader, T_finish, Yend);

    N_center = InitialData(N_x, x_vect, T_vect, Y_vect, M, T_start, T_finish, Ystart, Yend);

    Write_to_file2("detai2", fout, x_vect,
        T_vect, Y_vect, M, N_x, 2);


    integrate_Y_IDA(&chemkinReader, N_x, x_vect,
       T_vect, Y_vect, M, N_center, Ystart);
    //Integrate_Y(&chemkinReader, N_x, x_vect_new, T_vect_new, Y_vect_new, M, N_center, Y);
    //Write_to_file2("detail", fout, x_vect_new, T_vect_new, Y_vect_new, M, N_x, 1);
    return 0;
    //T_find();
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY KINSOL
 *--------------------------------------------------------------------
 */

 /*
  * System function for predator-prey system
  */
