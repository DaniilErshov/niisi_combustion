#include "functions.h"

double k_mol = pow(10, 3);
double Y_N2 = 0.745187;
double Y_max = 1 - Y_N2;
double P = 0.101325;
double R = 8.314;
double koeff_l = 0.4;
double l = 0.6;
double x_center;
long int myiter = 0;
long int nniters;
double eps_func = pow(10, -8);
double Tstart = 300;
double Tfinish = 2357;
double eps_x = pow(10, -6);
double eps_fr = pow(10, -6);
const double kB = 1.3806504e-23;
const double Angstroem__ = 1.0e-10;
const double santimetr = 1.0e-8;
const vector<double> M = { 2., 1.0, 32.0, 16.0, 17.0, 33.0, 18.0, 34.0, 28.0 };

string name_species[9] = { "H2", "H", "O2",
  "O", "OH", "HO2", "H2O", "H2O2", "N2" };

int main()
{
    init_consts(num_gas_species, num_react);
    int N_x = 25;
    double b = 0.00001;
    double M;
    double W, rho, Y_H2, Y_O2;
    int N_center;
    int retval;
    double w_dot;
    double* my_x;
    ofstream fout;
    double X_H2O, X_H2, X_O2, X_N2;
    double tout1 = pow(10, -7);
    double Tend = 2380;
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
    double* X = new double[num_gas_species];


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

    cout << "N_center = " << N_center << "\n";
    cout << "T N_center = " << T_vect[N_center] << "\n";
    Write_to_file2("detail", fout, x_vect,
        T_vect, Y_vect, M, N_x, 1);
    std::cout << "Chemkin Reader Test\n";

   //Find_final_state_IDA(&chemkinReader, T_finish, Yend);

   cout << "Cp = " << Cp_all(300, Ystart) << "\n";
   //Yend[0] = my_mol_weight(0) * 0.01457938;
   //Yend[1] = my_mol_weight(1) * 0.002126551;
   //Yend[2] = my_mol_weight(2) * 0.007320253
   //    ;
   //Yend[3] = my_mol_weight(3) * 7.984792E-4
   //    ;
   //Yend[4] = my_mol_weight(4) * 0.008885982
   //    ;
   //Yend[5] = my_mol_weight(5) * 1.484319E-6
   //    ;
   //Yend[6] = my_mol_weight(6) * 0.3197389
   //    ;
   //Yend[7] = my_mol_weight(7) * 1.874924E-7
   //    ;
   //Yend[8] = my_mol_weight(8) * 0.6465487
   //    ;
   //for (int i = 0; i < num_gas_species; i++) {
   //    Yend[i] /= 24.30197;
   //}
   Yend[H2] = 0.00119584;
   Yend[H] = 7.25484e-05;
   Yend[O2] = 0.00710046;
   Yend[O] = 0.00038046;
   Yend[OH] = 0.00549423;
   Yend[HO2] = 1.65194e-06;
   Yend[H2O] = 0.240568;
   Yend[H2O2] = 1.81769e-07;
   Yend[N2] = 0.745187;
   //Yend[0] = 0.04252821
   //    ;
   //Yend[1] = 0.03886508
   //    ;
   //Yend[2] = 0.02825736

   //    ;
   //Yend[3] = 0.009236109
   //    ;
   //Yend[4] = 0.01344398

   //    ;
   //Yend[5] = 1.276131E-5

   //    ;
   //Yend[6] = 0.2503474
   //    ;
   //Yend[7] = 5.529905E-6
   //    ;
   //Yend[8] = 0.6173036
   //    ;
   //cout << "lambda = " << Lambda_All(Yend, 1681.869);
  /* cout << "lambda = " << Lambda_All(Yend, 2357.545);
   ofstream f_lambda;
   f_lambda.open("lambda_H2.dat");
   f_lambda << "TITLE=\"" << "Graphics" << "\"" << endl;
   f_lambda << R"(VARIABLES= "T", "lambda_H2", "lambda_OH", "lambda_H2O", "lambda_N2")";
   double T_check = 500;
   double dT_check = 25;
   while (T_check < 3000) {
       f_lambda << T_check << " " << Lambda_H2(T_check) <<  " " << Lambda_OH(T_check)
           << " " << Lambda_H2O(T_check) <<  " " << Lambda_N2(T_check) << "\n";
       T_check += dT_check;
   }
   f_lambda.close();*/
  /* double T_check = 500;
   double dT_check = 25;
   for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
       T_check = 500;
       for (int i = 0; i < num_gas_species; i++) {
           Yend[i] = 0;
       }
       Yend[k_spec] = 1;
       fout.open(name_species[k_spec] + ".dat");
       fout << "TITLE=\"" << "Graphics" << "\"" << endl;
       fout << R"(VARIABLES= "T", "lambda")" << endl;
       while (T_check < 3000) {
           fout << T_check << " " << Lambda_All(&chemkinReader, Yend, T_check) << endl;
           T_check += dT_check;
       }
       fout.close();
    }*/

   //cout << "dij = " << Dij_func(&chemkinReader, 0, 1, 300, Ystart) << "\n";
   //double Dcheck = 0;


   //Get_mole_fr(X, Yend);
   //int my_spec = 0;

   //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
   //    if (my_spec != k_spec) {
   //         Dcheck += X[k_spec] / Dij_func(&chemkinReader, my_spec, k_spec, 2357.545, Yend);
   //    }
   //    
   //}
   //cout << "Dkm = " << (1 - Yend[my_spec]) / Dcheck << "\n";

    N_center = InitialData(N_x, x_vect, T_vect, Y_vect, M, T_start, T_finish, Ystart, Yend);
    //Add_elem(T_vect, Y_vect, x_vect, N_x, N_center, b);
    //Add_elem(T_vect, Y_vect, x_vect, N_x, N_center, b);

    Write_to_file2("A_initial", fout, x_vect,
        T_vect, Y_vect, M, N_x, 2);

    cout << "N_center = " << N_center << "\n";
    cout << "T N_center = " << T_vect[N_center] << "\n";

    //find_M(&chemkinReader, N_x, x_vect, T_vect, Y_vect, M, N_center, Ystart);
    cout << "M = " << M << "\n";
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
