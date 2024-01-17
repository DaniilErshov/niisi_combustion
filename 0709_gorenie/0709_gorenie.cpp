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
    int N_x = 19;
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
    vector<double> x_vect = { 0.0,
 0.05,
 0.098,
 0.1,
 0.1083333,
 0.1166667,
 0.1208333,
 0.125,
 0.1270833,
 0.1291667,
 0.13125,
 0.1333333,
 0.1416667,
 0.15,
 0.1666667,
 0.1833333,
 0.2,
 0.25,
 0.3 };
    vector<double> T_vect = {
        300.0,
 300.0056,
 300.1539,
 300.244,
 301.6694,
 309.4851,
 326.1897,
 375.8476,
 429.2833,
 508.214,
 603.9516,
 704.3042,
 1074.065,
 1343.69,
 1605.0,
 1746.343,
 1838.286,
 1971.049,
 1971.049
    };
    //vector<double> Y_vect = { 0.29577, 1.2287e-12, 0.14789, 6.0102e-11, 2.7955e-10, 9.7936e-09, 1.4416e-07, 4.9837e-09, 0.55635, 0.2956, 1.4137e-11, 0.14792, 2.2377e-09, 1.0592e-08, 5.5992e-07, 7.6102e-06, 2.8671e-07, 0.55647, 0.29275, 2.102e-08, 0.14831, 8.9039e-08, 7.1014e-07, 3.0303e-05, 0.00037888, 1.535e-05, 0.55852, 0.29158, 2.9945e-06, 0.14829, 5.9167e-07, 7.1588e-06, 8.8276e-05, 0.00078592, 1.9662e-05, 0.55922, 0.2614, 0.0010393, 0.14063, 4.2713e-05, 0.00047397, 0.00025766, 0.024346, 4.3767e-05, 0.57177, 0.17956, 0.015469, 0.10487, 0.0013509, 0.0014159, 0.00012256, 0.10631, 6.1994e-06, 0.5909, 0.079078, 0.0417, 0.049969, 0.0061405, 0.0055246, 3.8457e-05, 0.21235, 3.7622e-06, 0.60519, 0.05032, 0.041567, 0.02843, 0.0089082, 0.013172, 1.2731e-05, 0.24653, 4.931e-06, 0.61105, 0.049027, 0.034385, 0.02266, 0.009695, 0.02092, 6.068e-06, 0.2502, 3.0096e-06, 0.6131, 0.050823, 0.028978, 0.020174, 0.0097374, 0.027773, 4.0324e-06, 0.24852, 1.2855e-06, 0.61399, 0.04858, 0.025412, 0.018648, 0.0084828, 0.026405, 3.7945e-06, 0.25562, 1.1861e-06, 0.61685, 0.04858, 0.025412, 0.018648, 0.0084828, 0.026405, 3.7945e-06, 0.25562, 1.1861e-06, 0.61685 };
    vector<double> Y_vect = { 0.02851092,
 6.699181e-16,
 0.2262768,
 2.396985e-12,
 3.263625e-12,
 5.127868e-11,
 7.106282e-10,
 2.353197e-11,
 0.7452123,
 0.02850836,
 6.443994e-15,
 0.2262776,
 7.407263e-11,
 1.026178e-10,
 2.426414e-09,
 3.075666e-08,
 1.120423e-09,
 0.745214,
 0.028471,
 3.500167e-13,
 0.2262878,
 2.187201e-09,
 5.555058e-09,
 1.005447e-07,
 1.250578e-06,
 5.025773e-08,
 0.7452398,
 0.02845818,
 8.245982e-13,
 0.2262908,
 3.632783e-09,
 1.225712e-08,
 1.909989e-07,
 2.370516e-06,
 9.986949e-08,
 0.7452484,
 0.02832538,
 1.140046e-11,
 0.2263119,
 2.828357e-08,
 1.336385e-07,
 2.349612e-06,
 2.732859e-05,
 1.298455e-06,
 0.7453316,
 0.02788103,
 4.832266e-10,
 0.2262928,
 1.828889e-07,
 1.044985e-06,
 2.089102e-05,
 0.0002240663,
 1.155765e-05,
 0.7455685,
 0.02728935,
 2.478752e-08,
 0.2259902,
 7.10704e-07,
 5.446663e-06,
 8.524923e-05,
 0.0008053601,
 4.2856e-05,
 0.7457808,
 0.02613659,
 8.242283e-07,
 0.2243679,
 3.822436e-06,
 4.181224e-05,
 0.0003255065,
 0.003052685,
 0.0001602669,
 0.7459106,
 0.02523548,
 6.831788e-06,
 0.2219418,
 1.399005e-05,
 0.000177199,
 0.0004993084,
 0.006046387,
 0.0002809691,
 0.7457981,
 0.02409714,
 3.148146e-05,
 0.217497,
 4.070112e-05,
 0.0004923556,
 0.000531871,
 0.01138295,
 0.0004037409,
 0.7455228,
 0.02275174,
 9.636661e-05,
 0.2108316,
 8.845408e-05,
 0.000856491,
 0.0004631158,
 0.01936536,
 0.0004168005,
 0.7451301,
 0.02125152,
 0.0002108554,
 0.2022239,
 0.0001834615,
 0.001064418,
 0.0003736634,
 0.0296957,
 0.0003314137,
 0.7446651,
 0.01420762,
 0.001050501,
 0.1495276,
 0.0018498,
 0.001575279,
 0.0001505922,
 0.08807725,
 3.114399e-05,
 0.7435302,
 0.007907429,
 0.001855132,
 0.08962115,
 0.004805238,
 0.003933212,
 6.531316e-05,
 0.147349,
 7.44082e-06,
 0.7444561,
 0.004158028,
 0.001878083,
 0.04798811,
 0.00665292,
 0.007779996,
 2.493008e-05,
 0.1862175,
 8.431525e-06,
 0.745292,
 0.003541522,
 0.001521342,
 0.03770094,
 0.006242671,
 0.009504854,
 1.605074e-05,
 0.1959969,
 8.012257e-06,
 0.7454677,
 0.003348939,
 0.001246234,
 0.03355879,
 0.005519547,
 0.01020875,
 1.251722e-05,
 0.2006078,
 6.845997e-06,
 0.7454906,
 0.00308068,
 0.0008852767,
 0.02825876,
 0.004261894,
 0.01065343,
 8.669572e-06,
 0.2073325,
 4.545396e-06,
 0.7455143,
 0.00308068,
 0.0008852767,
 0.02825876,
 0.004261894,
 0.01065343,
 8.669572e-06,
 0.2073325,
 4.545396e-06,
 0.7455143 };

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

   //Find_final_state_IDA(&chemkinReader, T_finish, Yend);

   //cout << "Cp = " << Cp_all(300, Ystart) << "\n";
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
   //Yend[H2] = 0.00119584;
   //Yend[H] = 7.25484e-05;
   //Yend[O2] = 0.00710046;
   //Yend[O] = 0.00038046;
   //Yend[OH] = 0.00549423;
   //Yend[HO2] = 1.65194e-06;
   //Yend[H2O] = 0.240568;
   //Yend[H2O2] = 1.81769e-07;
   //Yend[N2] = 0.745187;
    Yend[H2] = 1.4563E-002;
    Yend[H] = 1.7917E-003;
    Yend[O2] = 5.4673E-003;
    Yend[O] = 5.9257E-004;
    Yend[OH] = 7.8587E-003;
    Yend[HO2] = 1.2247E-006;
    Yend[H2O] = 3.2394E-001;
    Yend[H2O2] = 1.3135E-007;
    Yend[N2] = 6.4579E-001;


   Ystart[H2] = 0.02851092;
   Ystart[H] = 6.699181e-16;
   Ystart[O2] = 0.2262768;
   Ystart[O] = 2.396985e-12;
   Ystart[OH] = 3.263625e-12;
   Ystart[HO2] = 5.127868e-11;
   Ystart[H2O] = 7.106282e-10;
   Ystart[H2O2] = 2.353197e-11;
   Ystart[N2] = 0.7452123;


   N_center = 13;
   M = 0.232643;

    Write_to_file2("A_initial", fout, x_vect,
        T_vect, Y_vect, Y_vect, M, N_x, 1);

    cout << "N_center = " << N_center << "\n";
    cout << "T N_center = " << T_vect[N_center] << "\n";
    T_center = T_vect[N_center];
    cout << "M = " << M << "\n";

    integrate_All_IDA(N_x, x_vect,
       T_vect, Y_vect, M, N_center, Ystart, 1);

    for (int i = 0; i < Y_vect.size(); i++) {
        if (Y_vect[i] < 0) Y_vect[i] = 0;
    }
    Write_to_file2("Ida_1", fout, x_vect,
        T_vect, Y_vect, Y_vect, M, N_x, 1);

    //Integrate_Kinsol(N_x, x_vect,
    //    T_vect, Y_vect, M, N_center, Ystart, 1);

    //Write_to_file2("Kinsol_1", fout, x_vect,
    //    T_vect, Y_vect, Y_vect, M, N_x, 1);

    cout << "size1 = " << T_vect.size() << "\n";
    Add_elem(T_vect, Y_vect, x_vect, N_x, N_center, 0.1, 3, 0);
    Write_to_file2("A_add", fout, x_vect,
        T_vect, Y_vect, Y_vect, M, N_x, 1);


    for (int i = 0; i < T_vect.size(); i++) {
        if (T_vect[i] == T_center) {
            N_center = i;
            cout << "new = " << N_center << "\n";
        }
    }
    for (int i = 0; i < Y_vect.size(); i++) {
        if (Y_vect[i] < 0) Y_vect[i] = 0;
    }

   /* integrate_Y_IDA(N_x, x_vect,
        T_vect, Y_vect, M, N_center, Ystart);*/
    cout << "size2 = " << T_vect.size() << "\n";
    cout << "M = " << M << "\n";
    for (int i = 0; i < Y_vect.size(); i++) {
        if (Y_vect[i] < 0) Y_vect[i] = 0;
    }
    integrate_All_IDA(N_x, x_vect,
        T_vect, Y_vect, M, N_center, Ystart, 2);

    Add_elem(T_vect, Y_vect, x_vect, N_x, N_center, 0.05, 3, 0);
    for (int i = 0; i < T_vect.size(); i++) {
        if (T_vect[i] == T_center) {
            N_center = i;
            cout << "new = " << N_center << "\n";
        }
    }

    for (int i = 0; i < Y_vect.size(); i++) {
        if (Y_vect[i] < 0) Y_vect[i] = 0;
    }

    integrate_All_IDA(N_x, x_vect,
        T_vect, Y_vect, M, N_center, Ystart, 3);

    Add_elem(T_vect, Y_vect, x_vect, N_x, N_center, 0.01, 5, 0);
    for (int i = 0; i < T_vect.size(); i++) {
        if (T_vect[i] == T_center) {
            N_center = i;
            cout << "new = " << N_center << "\n";
        }
    }
    for (int i = 0; i < Y_vect.size(); i++) {
        if (Y_vect[i] < 0) Y_vect[i] = 0;
    }

    integrate_All_IDA(N_x, x_vect,
        T_vect, Y_vect, M, N_center, Ystart, 4);
    Add_elem(T_vect, Y_vect, x_vect, N_x, N_center, 0.005, 3, 0);

    for (int i = 0; i < Y_vect.size(); i++) {
        if (Y_vect[i] < 0) Y_vect[i] = 0;
    }

    integrate_All_IDA(N_x, x_vect,
        T_vect, Y_vect, M, N_center, Ystart, 5);
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
