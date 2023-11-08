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
double Tfinish = 2380;
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
    int N_x = 20;
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
 /*  Yend[0] = my_mol_weight(0) * 0.01457938;
   Yend[1] = my_mol_weight(1) * 0.002126551;
   Yend[2] = my_mol_weight(2) * 0.007320253
       ;
   Yend[3] = my_mol_weight(3) * 7.984792E-4
       ;
   Yend[4] = my_mol_weight(4) * 0.008885982
       ;
   Yend[5] = my_mol_weight(5) * 1.484319E-6
       ;
   Yend[6] = my_mol_weight(6) * 0.3197389
       ;
   Yend[7] = my_mol_weight(7) * 1.874924E-7
       ;
   Yend[8] = my_mol_weight(8) * 0.6465487
       ;
   for (int i = 0; i < num_gas_species; i++) {
       Yend[i] /= 24.30197;
   }

   cout << "lambda = " << Lambda_All(&chemkinReader, Yend, 2357) << "\n";*/
   double T_check = 500;
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
    }

   cout << "dij = " << Dij_func(&chemkinReader, 0, 1, 300, Ystart) << "\n";
   double Dcheck = 0;

   //{
   //    double T_check = 500;
   //    double dT_check = 25;
   //    ofstream fout16;
   //    double k_0_f[3], k_inf_f[3], k_0_r[3], k_inf_r[3];
   //    double c[3], m[3], d = 0.14;
   //    double Pr_f[3], Pr_r[3];
   //    int k = 0, l = 0;
   //    double logF_f, logF_core_f, logF_r, logF_core_r;

   //    double sum1, sum2;



   //    double* forward = new double[num_react];
   //    double* reverse = new double[num_react];
   //    double* equilib = new double[num_react];
   //    while (T_check < 5000) {


   //        fout.open("15.dat");
   //        fout << "TITLE=\"" << "Graphics" << "\"" << endl;
   //        fout << R"(VARIABLES= "1 / T", "forward[15]", , "reverse[15]")" << endl;
   //        fout16.open("16.dat");
   //        fout16 << "TITLE=\"" << "Graphics" << "\"" << endl;
   //        fout16 << R"(VARIABLES= "1 / T", "forward[16]", , "reverse[16]")" << endl;
   //        while (T_check < 3000) {
   //            fout << T_check << " " << Lambda_All(&chemkinReader, Yend, T_check) << endl;
   //            T_check += dT_check;
   //        }
   //        fout.close();
   //        for (int i = 0; i < num_react; i++) {
   //            if (i != 8 && i != 15 && i != 16) {
   //                forward[i] = chec.kPrex_f[i] * pow(T_check, chec.kPow_f[i])
   //                    * exp(-chec.kE_f[i] / T_check / phyc.kRc);
   //                reverse[i] = chec.kPrex_r[i] * pow(T_check, chec.kPow_r[i])
   //                    * exp(-chec.kE_r[i] / T_check / phyc.kRc);
   //            }
   //            else {
   //                if (i == 8) k = 0;
   //                if (i == 15) k = 1;
   //                if (i == 16) k = 2;

   //                k_inf_f[k] = chec.kPrex_f[i] * pow(T_check, chec.kPow_f[i])
   //                    * exp(-chec.kE_f[i] / T_check / phyc.kRc);
   //                k_inf_r[k] = chec.kPrex_r[i] * pow(T_check, chec.kPow_r[i])
   //                    * exp(-chec.kE_r[i] / T_check / phyc.kRc);
   //                k_0_f[k] = chec.kPrex_f_lp[i] * pow(T_check, chec.kPow_f_lp[i])
   //                    * exp(-chec.kE_f_lp[i] / T_check / phyc.kRc);
   //                k_0_r[k] = chec.kPrex_r_lp[i] * pow(T_check, chec.kPow_r_lp[i])
   //                    * exp(-chec.kE_r_lp[i] / T_check / phyc.kRc);

   //                c[k] = -0.4 - 0.67 * log10(chec.Fcent[i]);
   //                m[k] = 0.75 - 1.27 * log10(chec.Fcent[i]);
   //            }
   //        }

   //        Pr_f[0] = (k_0_f[0] * (1.3 * y[0] + 10 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8])) / k_inf_f[0];
   //        Pr_r[0] = (k_0_r[0] * (1.3 * y[0] + 10 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8])) / k_inf_r[0];
   //        Pr_f[1] = (k_0_f[1] * (3.7 * y[0] + 1.5 * y[8] + 1.2 * y[2] + 7.7 * y[7] + y[1] + y[3] + y[4] + y[5] + y[6])) / k_inf_f[1];
   //        Pr_r[1] = (k_0_r[1] * (3.7 * y[0] + 1.5 * y[8] + 1.2 * y[2] + 7.7 * y[7] + y[1] + y[3] + y[4] + y[5] + y[6])) / k_inf_r[1];
   //        Pr_f[2] = (k_0_f[2] * y[6]) / k_inf_f[2];
   //        Pr_r[2] = (k_0_r[2] * y[6]) / k_inf_r[2];

   //        for (k = 0; k < 3; k++) {
   //            if (k == 0) l = 8;
   //            if (k == 1) l = 15;
   //            if (k == 2) l = 16;

   //            if (Pr_f[k] == 0) forward[l] = k_inf_f[k];
   //            else {
   //                logF_core_f = pow((log10(Pr_f[k]) + c[k]) / (m[k] - d * (log10(Pr_f[k]) + c[k])), 2);
   //                logF_f = pow(1.0 + logF_core_f, -1) * log10(chec.Fcent[l]);
   //                chec.F_f[l] = pow(10, logF_f);
   //                forward[l] = k_inf_f[k] * (Pr_f[k] / (1 + Pr_f[k])) * chec.F_f[l];
   //            }

   //            if (Pr_r[k] == 0) reverse[l] = k_inf_r[k];
   //            else {
   //                logF_core_r = pow((log10(Pr_r[k]) + c[k]) / (m[k] - d * (log10(Pr_r[k]) + c[k])), 2);
   //                logF_r = pow(1.0 + logF_core_r, -1) * log10(chec.Fcent[l]);
   //                chec.F_r[l] = pow(10, logF_r);
   //                reverse[l] = k_inf_r[k] * (Pr_r[k] / (1 + Pr_r[k])) * chec.F_r[l];
   //            }
   //        }
   //        fout << 1. / T_check << " " << forward[15] << " " << reverse[15] << endl;
   //        fout16 << 1. / T_check << " " << forward[16] << " " << reverse[16] << endl;
   //        T_check += dT_check;
   //    }

   //    delete[] forward;
   //    delete[] reverse;
   //    delete[] equilib;
   //    fout.close();
   //    fout16.close();
   //}

   Get_mole_fr(X, Ystart);
   int my_spec = 6;

   for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
       if (my_spec != k_spec) {
            Dcheck += X[k_spec] / Dij_func(&chemkinReader, my_spec, k_spec, 300, Ystart);
       }
       
   }
   cout << "Dkm = " << (1 - Ystart[my_spec]) / Dcheck << "\n";

    N_center = InitialData(N_x, x_vect, T_vect, Y_vect, M, T_start, T_finish, Ystart, Yend);
    //Add_elem(T_vect, Y_vect, x_vect, N_x, N_center, b);
    //Add_elem(T_vect, Y_vect, x_vect, N_x, N_center, b);

    Write_to_file2("detai2", fout, x_vect,
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
