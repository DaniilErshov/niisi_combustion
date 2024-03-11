#include <cmath>
#include <iostream>
#include "constants.h"
using namespace std;

struct phy_consts phyc;
struct che_consts chec;

void init_consts(int& num_gas_species, int& num_react)
{
    const std::string chemfile = R"(C:\Users\Данило\source\repos\check_boost3\ChemKin_reader\test\chem2.inp)";
    const std::string thermfile = R"(C:\Users\Данило\source\repos\check_boost3\ChemKin_reader\test\therm-abf.dat)";
    const std::string transfile = R"(C:\Users\Данило\source\repos\check_boost3\ChemKin_reader\test\tran.dat)";
    chec.chemkinReader = new IO::ChemkinReader(chemfile, thermfile, transfile);
    chec.chemkinReader->read();
    chec.chemkinReader->check();
    std::cout << "NEW BLOCK" << std::endl;
    num_react = chec.chemkinReader->reactions().size();
    num_gas_species = chec.chemkinReader->species().size();
    chec.sum_v = new double[num_react];
    //phyc.kR = 8.314472e+7;  // universal gas constant ( erg/K/mole )
    //phyc.kRc = 1.987207e-3; // universal gas constant ( kcal/K/mole )
    //phyc.kTime = 1.e-3;     // dimensional time ( sec )
    //phyc.kLength = 1.e+2;   // dimensional length ( cm )
    //phyc.kMass = 1.e+3;     // dimensional mass ( g )
    //phyc.kTemp = 1.e0;      // dimensional temperature ( K )
    //phyc.kPres = 1.e+7;     // dimensional pressure ( dyn/cm**2 )
    //phyc.kDens = phyc.kMass / pow(phyc.kLength, 3); // dimensional density ( g/cm**3 )
    //phyc.kVel = phyc.kLength / phyc.kTime;           // dimensional velocity ( cm/sec )

    phyc.kR = 8.314472;  // universal gas constant ( erg/K/mole )
    phyc.kRc = 1.987207e-3; // universal gas constant ( kcal/K/mole )
    phyc.kTime = 1;     // dimensional time ( sec )
    phyc.kLength = 1;   // dimensional length ( cm )
    phyc.kMass = 1;     // dimensional mass ( g )
    phyc.kTemp = 1.e0;      // dimensional temperature ( K )
    phyc.kPres = 1.e+7;     // dimensional pressure ( dyn/cm**2 )
    phyc.kDens = phyc.kMass / pow(phyc.kLength, 3); // dimensional density ( g/cm**3 )
    phyc.kVel = phyc.kLength / phyc.kTime;           // dimensional velocity ( cm/sec )

    phyc.kR /= phyc.kMass * pow (phyc.kVel, 2);

    phyc.Cp_coef_hT = new double*[num_gas_species];
    phyc.Cp_coef_lT = new double*[num_gas_species];

    for(int i = 0; i < num_gas_species; i++) {
        phyc.Cp_coef_hT[i] = new double[9];
        phyc.Cp_coef_lT[i] = new double[9];
    }

    phyc.mol_weight = new double[num_gas_species];

    chec.kPrex_f = new double[num_react];
    chec.kPrex_f_lp = new double[num_react];
    chec.kPrex_r = new double[num_react];
    chec.kPrex_r_lp = new double[num_react];
    chec.kPow_f = new double[num_react];
    chec.kPow_f_lp = new double[num_react];
    chec.kPow_r = new double[num_react];
    chec.kPow_r_lp = new double[num_react];
    chec.kE_f = new double[num_react];
    chec.kE_f_lp = new double[num_react];
    chec.kE_r = new double[num_react];
    chec.kE_r_lp = new double[num_react];
    chec.F_f = new double[num_react];
    chec.F_r = new double[num_react];
    chec.Fcent = new double[num_react];
    
    for(int i = 0; i < num_react; i++)
    {
        chec.kPrex_f_lp[i] = 0;
        chec.kPow_f_lp[i] = 0;
        chec.kE_f_lp[i] = 0;
        chec.kPrex_r_lp[i] = 0;
        chec.kPow_r_lp[i] = 0;
        chec.kE_r_lp[i] = 0;
    }

    

    chec.kPrex_f[0] = 1.04e+14 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[0] = 0.0;
    chec.kE_f[0] = 1.529e+01;


    chec.kPrex_f[1] = 5.08e+04 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[1] = 2.67;
    chec.kE_f[1] = 6.292;


    chec.kPrex_f[2] = 4.38e+13 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[2] = 0.0;
    chec.kE_f[2] = 6.99;


    chec.kPrex_f[3] = 2.97e+06 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[3] = 2.02;
    chec.kE_f[3] = 1.34e+01;


    chec.kPrex_f[4] = 4.577e+19 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[4] = -1.4;
    chec.kE_f[4] = 1.044e+02;


    chec.kPrex_f[5] = 6.165e+15 * phyc.kTime / pow(phyc.kLength, 6);
    chec.kPow_f[5] = -0.5;
    chec.kE_f[5] = 0.0;

    chec.kPrex_f[6] = 4.714e+18 * phyc.kTime / pow(phyc.kLength, 6);
    chec.kPow_f[6] = -1.0;
    chec.kE_f[6] = 0.0;


    chec.kPrex_f[7] = 3.5e+22 * phyc.kTime / pow(phyc.kLength, 6);
    chec.kPow_f[7] = -2.0;
    chec.kE_f[7] = 0.0;


    chec.kPrex_f[8] = 4.650e+12 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[8] = 0.44;
    chec.kE_f[8] = 0.0;
    chec.kPrex_f_lp[8] = 1.737e+19 * phyc.kTime / pow(phyc.kLength, 6);
    chec.kPow_f_lp[8] = -1.23;
    chec.kE_f_lp[8] = 0.0;
    chec.Fcent[8] = 0.67;

    chec.kPrex_f[9] = 5.176e+05 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[9] = 2.43;
    chec.kE_f[9] = 5.35e+01;


    chec.kPrex_f[10] = 7.079e+13 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[10] = 0.0;
    chec.kE_f[10] = 2.95e-01;


    chec.kPrex_f[11] = 3.250e+13 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[11] = 0.0;
    chec.kE_f[11] = 0.0;


    chec.kPrex_f[12] = 2.456e+13 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[12] = 0.0;
    chec.kE_f[12] = -4.97e-01;


    chec.kPrex_f[13] = 1.3e+11 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[13] = 0.0;
    chec.kE_f[13] = -1.63;


    chec.kPrex_f[14] = 3.658e+14 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[14] = 0.0;
    chec.kE_f[14] = 1.2e+01;


    chec.kPrex_f[15] = 2.0e+12 * phyc.kTime;
    chec.kPow_f[15] = 0.9;
    chec.kE_f[15] = 4.875e+01;
    chec.kPrex_f_lp[15] = 2.490e+24 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f_lp[15] = -2.3;
    chec.kE_f_lp[15] = 4.875e+01;

    chec.Fcent[15] = 0.43;

    chec.kPrex_f[16] = 2.0e+12 * phyc.kTime;
    chec.kPow_f[16] = 0.9;
    chec.kE_f[16] = 4.875e+01;
    chec.kPrex_f_lp[16] = 1.865e+25 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f_lp[16] = -2.3;
    chec.kE_f_lp[16] = 4.875e+01;

    chec.Fcent[16] = 0.51;

    chec.kPrex_f[17] = 2.410e+13 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[17] = 0.0;
    chec.kE_f[17] = 3.97;

    chec.kPrex_f[18] = 2.150e+10 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[18] = 1.0;
    chec.kE_f[18] = 6.0;


    chec.kPrex_f[19] = 9.550e+06 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[19] = 2.0;
    chec.kE_f[19] = 3.97;


    chec.kPrex_f[20] = 1.740e+12 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[20] = 0.0;
    chec.kE_f[20] = 3.18e-01;


    chec.kPrex_f[21] = 7.590e+13 * phyc.kTime / pow(phyc.kLength, 3);
    chec.kPow_f[21] = 0.0;
    chec.kE_f[21] = 7.269;


    phyc.mol_weight[0] = 2.01594; // H2
    phyc.mol_weight[1] = 1.00797; // H
    phyc.mol_weight[2] = 31.9988; // O2
    phyc.mol_weight[3] = 15.9994; // O
    phyc.mol_weight[4] = 17.00737; // OH
    phyc.mol_weight[5] = 33.00677; // HO2
    phyc.mol_weight[6] = 18.01534; // H20
    phyc.mol_weight[7] = 34.01474; // H2O2
    phyc.mol_weight[8] = 28.0134; // N2

        // Coefficients for NASA polynoms
        // 8th coef is for enthalpy calculation
        // High temperature (T > 1000 K)

        // H2
        phyc.Cp_coef_hT[0][0] = 5.608128010e+05;
        phyc.Cp_coef_hT[0][1] = -8.371504740e+02;
        phyc.Cp_coef_hT[0][2] = 2.975364532e+00;
        phyc.Cp_coef_hT[0][3] = 1.252249124e-03;
        phyc.Cp_coef_hT[0][4] = -3.740716190e-07;
        phyc.Cp_coef_hT[0][5] = 5.936625200e-11;
        phyc.Cp_coef_hT[0][6] = -3.606994100e-15;
        phyc.Cp_coef_hT[0][7] = 5.339824410e+03;
        phyc.Cp_coef_hT[0][8] = -2.202774769E+00;

        // H
        phyc.Cp_coef_hT[1][0] = 6.078774250e+01;
        phyc.Cp_coef_hT[1][1] = -1.819354417e-01;
        phyc.Cp_coef_hT[1][2] = 2.500211817e+00;
        phyc.Cp_coef_hT[1][3] = -1.226512864e-07;
        phyc.Cp_coef_hT[1][4] = 3.732876330e-11;
        phyc.Cp_coef_hT[1][5] = -5.687744560e-15;
        phyc.Cp_coef_hT[1][6] = 3.410210197e-19;
        phyc.Cp_coef_hT[1][7] = 2.547486398e+04;
        phyc.Cp_coef_hT[1][8] = -4.481917770E-01;
    
        // O2
        phyc.Cp_coef_hT[2][0] = -1.037939022e+06;
        phyc.Cp_coef_hT[2][1] = 2.344830282e+03;
        phyc.Cp_coef_hT[2][2] = 1.819732036e+00;
        phyc.Cp_coef_hT[2][3] = 1.267847582e-03;
        phyc.Cp_coef_hT[2][4] = -2.188067988e-07;
        phyc.Cp_coef_hT[2][5] = 2.053719572e-11;
        phyc.Cp_coef_hT[2][6] = -8.193467050e-16;
        phyc.Cp_coef_hT[2][7] = -1.689010929e+04;
        phyc.Cp_coef_hT[2][8] = 1.738716506E+01;

        // O
        phyc.Cp_coef_hT[3][0] = 2.619020262e+05;
        phyc.Cp_coef_hT[3][1] = -7.298722030e+02;
        phyc.Cp_coef_hT[3][2] = 3.317177270e+00;
        phyc.Cp_coef_hT[3][3] = -4.281334360e-04;
        phyc.Cp_coef_hT[3][4] = 1.036104594e-07;
        phyc.Cp_coef_hT[3][5] = -9.438304330e-12;
        phyc.Cp_coef_hT[3][6] = 2.725038297e-16;
        phyc.Cp_coef_hT[3][7] = 3.392428060e+04;
        phyc.Cp_coef_hT[3][8] = -6.679585350E-01;
    
        // OH
        phyc.Cp_coef_hT[4][0] = 1.017393379e+06;
        phyc.Cp_coef_hT[4][1] = -2.509957276e+03;
        phyc.Cp_coef_hT[4][2] = 5.116547860e+00;
        phyc.Cp_coef_hT[4][3] = 1.305299930e-04;
        phyc.Cp_coef_hT[4][4] = -8.284322260e-08;
        phyc.Cp_coef_hT[4][5] = 2.006475941e-11;
        phyc.Cp_coef_hT[4][6] = -1.556993656e-15;
        phyc.Cp_coef_hT[4][7] = 2.019640206e+04;
        phyc.Cp_coef_hT[4][8] = -1.101282337E+01;
    
        // HO2
        phyc.Cp_coef_hT[5][0] = -1.810669724e+06;
        phyc.Cp_coef_hT[5][1] = 4.963192030e+03;
        phyc.Cp_coef_hT[5][2] = -1.039498992e+00;
        phyc.Cp_coef_hT[5][3] = 4.560148530e-03;
        phyc.Cp_coef_hT[5][4] = -1.061859447e-06;
        phyc.Cp_coef_hT[5][5] = 1.144567878e-10;
        phyc.Cp_coef_hT[5][6] = -4.763064160e-15;
        phyc.Cp_coef_hT[5][7] = -3.200817190e+04;
        phyc.Cp_coef_hT[5][8] = 4.066850920E+01;
        
        // H2O
        phyc.Cp_coef_hT[6][0] = 1.034972096e+06;
        phyc.Cp_coef_hT[6][1] = -2.412698562e+03;
        phyc.Cp_coef_hT[6][2] = 4.646110780e+00;
        phyc.Cp_coef_hT[6][3] = 2.291998307e-03;
        phyc.Cp_coef_hT[6][4] = -6.836830480e-07;
        phyc.Cp_coef_hT[6][5] = 9.426468930e-11;
        phyc.Cp_coef_hT[6][6] = -4.822380530e-15;
        phyc.Cp_coef_hT[6][7] = -1.384286509e+04;
        phyc.Cp_coef_hT[6][8] = -7.978148510E+00;
    
        // H2O2
        phyc.Cp_coef_hT[7][0] = 1.489428027e+06;
        phyc.Cp_coef_hT[7][1] = -5.170821780e+03;
        phyc.Cp_coef_hT[7][2] = 1.128204970e+01;
        phyc.Cp_coef_hT[7][3] = -8.042397790e-05;
        phyc.Cp_coef_hT[7][4] = -1.818383769e-08;
        phyc.Cp_coef_hT[7][5] = 6.947265590e-12;
        phyc.Cp_coef_hT[7][6] = -4.827831900e-16;
        phyc.Cp_coef_hT[7][7] = 1.418251038e+04;
        phyc.Cp_coef_hT[7][8] = -4.650855660E+01;

        // N2
        phyc.Cp_coef_hT[8][0] = 5.877124060e+05;
        phyc.Cp_coef_hT[8][1] = -2.239249073e+03;
        phyc.Cp_coef_hT[8][2] = 6.066949220e+00;
        phyc.Cp_coef_hT[8][3] = -6.139685500e-04;
        phyc.Cp_coef_hT[8][4] = 1.491806679e-07;
        phyc.Cp_coef_hT[8][5] = -1.923105485e-11;
        phyc.Cp_coef_hT[8][6] = 1.061954386e-15;
        phyc.Cp_coef_hT[8][7] = 1.283210415e+04;
        phyc.Cp_coef_hT[8][8] = -1.586640027E+01;

        // H2
        phyc.Cp_coef_lT[0][0] = 4.078323210e+04;
        phyc.Cp_coef_lT[0][1] = -8.009186040e+02;
        phyc.Cp_coef_lT[0][2] = 8.214702010e+00;
        phyc.Cp_coef_lT[0][3] = -1.269714457e-02;
        phyc.Cp_coef_lT[0][4] = 1.753605076e-05;
        phyc.Cp_coef_lT[0][5] = -1.202860270e-08;
        phyc.Cp_coef_lT[0][6] = 3.368093490e-12;
        phyc.Cp_coef_lT[0][7] = 2.682484665e+03;
        phyc.Cp_coef_lT[0][8] = -3.043788844E+01;

        // H
        phyc.Cp_coef_lT[1][0] = 0.000000000e+00;
        phyc.Cp_coef_lT[1][1] = 0.000000000e+00;
        phyc.Cp_coef_lT[1][2] = 2.500000000E+00;
        phyc.Cp_coef_lT[1][3] = 0;
        phyc.Cp_coef_lT[1][4] = 0;
        phyc.Cp_coef_lT[1][5] = 0;
        phyc.Cp_coef_lT[1][6] = 0;
        phyc.Cp_coef_lT[1][7] = 2.547370801E+04;
        phyc.Cp_coef_lT[1][8] = -4.466828530E-01;
        
        // O2
        phyc.Cp_coef_lT[2][0] = -3.425563420e+04;
        phyc.Cp_coef_lT[2][1] = 4.847000970e+02;
        phyc.Cp_coef_lT[2][2] = 1.119010961e+00;
        phyc.Cp_coef_lT[2][3] = 4.293889240e-03;
        phyc.Cp_coef_lT[2][4] = -6.836300520e-07;
        phyc.Cp_coef_lT[2][5] = -2.023372700e-09;
        phyc.Cp_coef_lT[2][6] = 1.039040018e-12;
        phyc.Cp_coef_lT[2][7] = -3.391454870e+03;
        phyc.Cp_coef_lT[2][8] = 1.849699470E+01;
        
        // O
        phyc.Cp_coef_lT[3][0] = -7.953611300e+03;
        phyc.Cp_coef_lT[3][1] = 1.607177787e+02;
        phyc.Cp_coef_lT[3][2] = 1.966226438e+00;
        phyc.Cp_coef_lT[3][3] = 1.013670310e-03;
        phyc.Cp_coef_lT[3][4] = -1.110415423e-06;
        phyc.Cp_coef_lT[3][5] =  6.517507500e-10;
        phyc.Cp_coef_lT[3][6] = -1.584779251e-13;
        phyc.Cp_coef_lT[3][7] = 2.840362437e+04;
        phyc.Cp_coef_lT[3][8] = 8.404241820E+00;
        
        // OH
        phyc.Cp_coef_lT[4][0] = -1.998858990e+03;
        phyc.Cp_coef_lT[4][1] = 9.300136160e+01;
        phyc.Cp_coef_lT[4][2] = 3.050854229e+00;
        phyc.Cp_coef_lT[4][3] = 1.529529288e-03;
        phyc.Cp_coef_lT[4][4] = -3.157890998e-06;
        phyc.Cp_coef_lT[4][5] = 3.315446180e-09;
        phyc.Cp_coef_lT[4][6] = -1.138762683e-12;
        phyc.Cp_coef_lT[4][7] = 2.991214235e+03;
        phyc.Cp_coef_lT[4][8] = 4.674110790E+00;

        // HO2
        phyc.Cp_coef_lT[5][0] = -7.598882540e+04;
        phyc.Cp_coef_lT[5][1] = 1.329383918e+03;
        phyc.Cp_coef_lT[5][2] = -4.677388240e+00;
        phyc.Cp_coef_lT[5][3] = 2.508308202e-02;
        phyc.Cp_coef_lT[5][4] = -3.006551588e-05;
        phyc.Cp_coef_lT[5][5] = 1.895600056e-08;
        phyc.Cp_coef_lT[5][6] = -4.828567390e-12;
        phyc.Cp_coef_lT[5][7] = -5.873350960e+03;
        phyc.Cp_coef_lT[5][8] = 5.193602140E+01;

        // H2O
        phyc.Cp_coef_lT[6][0] = -3.947960830e+04;
        phyc.Cp_coef_lT[6][1] = 5.755731020e+02;
        phyc.Cp_coef_lT[6][2] = 9.317826530e-01;
        phyc.Cp_coef_lT[6][3] = 7.222712860e-03;
        phyc.Cp_coef_lT[6][4] = -7.342557370e-06;
        phyc.Cp_coef_lT[6][5] = 4.955043490e-09;
        phyc.Cp_coef_lT[6][6] = -1.336933246e-12;
        phyc.Cp_coef_lT[6][7] = -3.303974310e+04;
        phyc.Cp_coef_lT[6][8] = 1.724205775E+01;
        
        // H2O2
        phyc.Cp_coef_lT[7][0] = -9.279533580e+04;
        phyc.Cp_coef_lT[7][1] = 1.564748385e+03;
        phyc.Cp_coef_lT[7][2] = -5.976460140e+00;
        phyc.Cp_coef_lT[7][3] = 3.270744520e-02;
        phyc.Cp_coef_lT[7][4] = -3.932193260e-05;
        phyc.Cp_coef_lT[7][5] = 2.509255235e-08;
        phyc.Cp_coef_lT[7][6] = -6.465045290e-12;
        phyc.Cp_coef_lT[7][7] = -2.494004728e+04;
        phyc.Cp_coef_lT[7][8] = 5.877174180E+01;
        
        // N2
        phyc.Cp_coef_lT[8][0] = 2.210371497e+04;
        phyc.Cp_coef_lT[8][1] = -3.818461820e+02;
        phyc.Cp_coef_lT[8][2] = 6.082738360e+00;
        phyc.Cp_coef_lT[8][3] = -8.530914410e-03;
        phyc.Cp_coef_lT[8][4] = 1.384646189e-05;
        phyc.Cp_coef_lT[8][5] = -9.625793620e-09;
        phyc.Cp_coef_lT[8][6] = 2.519705809e-12;
        phyc.Cp_coef_lT[8][7] = 7.108460860e+02;
        phyc.Cp_coef_lT[8][8] = -1.076003744E+01;

    for(int i = 0; i < num_gas_species; i++)
            phyc.mol_weight[i] /= phyc.kMass;

    // 8th coef. is for enthalpy calculation
        
    for (int component_i = 0; component_i < num_gas_species; component_i++)
    {
        for (int power_i = 0; power_i <= 8; power_i++)
        {
            // T < 1000 K
            phyc.Cp_coef_lT[component_i][power_i] *= phyc.kR / phyc.mol_weight[component_i];
            // T > 1000 K
            phyc.Cp_coef_hT[component_i][power_i] *= phyc.kR / phyc.mol_weight[component_i];
        }
    }

}

// Enthalpy of ith component

double get_dHiRT(double* Cp_coef, double T)
{
    double Hi;
        Hi = -1. * Cp_coef[0] * pow(T, -1) + Cp_coef[1] * log(T) + Cp_coef[2] * T + Cp_coef[3] * pow(T, 2) / 2 + Cp_coef[4] * pow(T, 3) / 3
        + Cp_coef[5] * pow(T, 4) / 4 + Cp_coef[6] * pow(T, 5) / 5 + Cp_coef[7];
    return Hi;
}

double get_dSiR(double* Cp_coef, double T) {
    double Si;
        Si = -Cp_coef[0] * pow(T, -2) / 2. - Cp_coef[1] * pow(T, -1) + Cp_coef[2] * log(T) + Cp_coef[3] * pow(T, 1) + Cp_coef[4] * pow(T, 2) / 2.
        + Cp_coef[5] * pow(T, 3) / 3. + Cp_coef[6] * pow(T, 4) / 4. + Cp_coef[8];
    return Si;
}
double get_dCpi(double* Cp_coef, double T)
{
    double Cpi;
        Cpi = Cp_coef[0] * pow(T, -2) + Cp_coef[1] * pow(T, -1) + Cp_coef[2] + Cp_coef[3] * pow(T, 1) + Cp_coef[4] * pow(T, 2)
        + Cp_coef[5] * pow(T, 3) + Cp_coef[6] * pow(T, 4);

    return Cpi;
}

double get_Hi(int component_i, double T)
{
    double Hi;
    int i = component_i;

    if (T > 1000)
        Hi = -1. * phyc.Cp_coef_hT[i][0] * pow(T, -1) + phyc.Cp_coef_hT[i][1] * log(T) + phyc.Cp_coef_hT[i][2] * T + phyc.Cp_coef_hT[i][3] * pow(T, 2) / 2 + phyc.Cp_coef_hT[i][4] * pow(T, 3) / 3
            + phyc.Cp_coef_hT[i][5] * pow(T, 4) / 4 + phyc.Cp_coef_hT[i][6] * pow(T, 5) / 5 + phyc.Cp_coef_hT[i][7];
    else
        Hi = -1. * phyc.Cp_coef_lT[i][0] * pow(T, -1) + phyc.Cp_coef_lT[i][1] * log(T) + phyc.Cp_coef_lT[i][2] * T + phyc.Cp_coef_lT[i][3] * pow(T, 2) / 2 + phyc.Cp_coef_lT[i][4] * pow(T, 3) / 3
            + phyc.Cp_coef_lT[i][5] * pow(T, 4) / 4 + phyc.Cp_coef_lT[i][6] * pow(T, 5) / 5 + phyc.Cp_coef_lT[i][7];
   /* if (T > 1000)
        Hi = -1. * phyc.Cp_coef_hT[i][0] * pow(T, -2) + phyc.Cp_coef_hT[i][1] * log(T) / T + phyc.Cp_coef_hT[i][2] + phyc.Cp_coef_hT[i][3] * T / 2 + phyc.Cp_coef_hT[i][4] * pow(T, 2) / 3
        + phyc.Cp_coef_hT[i][5] * pow(T, 3) / 4 + phyc.Cp_coef_hT[i][6] * pow(T, 4) / 5 + phyc.Cp_coef_hT[i][7] / T;
    else
        Hi = -1. * phyc.Cp_coef_lT[i][0] * pow(T, -2) + phyc.Cp_coef_lT[i][1] * log(T) / T + phyc.Cp_coef_lT[i][2] + phyc.Cp_coef_lT[i][3] * T / 2 + phyc.Cp_coef_lT[i][4] * pow(T, 2) / 3
        + phyc.Cp_coef_lT[i][5] * pow(T, 3) / 4 + phyc.Cp_coef_lT[i][6] * pow(T, 4) / 5 + phyc.Cp_coef_lT[i][7] / T;
    return Hi * T;*/
    return Hi;
}

// Specific heat of ith component
double get_Cpi(int component_i, double T)
{
    double Cpi;
    int i = component_i;

    if (T > 1000)
        Cpi = phyc.Cp_coef_hT[i][0] * pow(T, -2) + phyc.Cp_coef_hT[i][1] * pow(T, -1) + phyc.Cp_coef_hT[i][2] + phyc.Cp_coef_hT[i][3] * pow(T, 1) + phyc.Cp_coef_hT[i][4] * pow(T, 2)
            + phyc.Cp_coef_hT[i][5] * pow(T, 3) + phyc.Cp_coef_hT[i][6] * pow(T, 4);
    else
        Cpi = phyc.Cp_coef_lT[i][0] * pow(T, -2) + phyc.Cp_coef_lT[i][1] * pow(T, -1) + phyc.Cp_coef_lT[i][2] + phyc.Cp_coef_lT[i][3] * pow(T, 1) + phyc.Cp_coef_lT[i][4] * pow(T, 2)
            + phyc.Cp_coef_lT[i][5] * pow(T, 3) + phyc.Cp_coef_lT[i][6] * pow(T, 4);
    
    return Cpi;
}

double get_Si(int component_i, double T) {
    double Si;
    int i = component_i;
    if (T > 1000)
        Si = -phyc.Cp_coef_hT[i][0] * pow(T, -2) / 2. - phyc.Cp_coef_hT[i][1] * pow(T, -1) + phyc.Cp_coef_hT[i][2] * log(T) + phyc.Cp_coef_hT[i][3] * pow(T, 1) + phyc.Cp_coef_hT[i][4] * pow(T, 2) / 2.
        + phyc.Cp_coef_hT[i][5] * pow(T, 3) / 3. + phyc.Cp_coef_hT[i][6] * pow(T, 4) / 4. + phyc.Cp_coef_hT[i][8];
    else
        Si = -phyc.Cp_coef_lT[i][0] * pow(T, -2) / 2. - phyc.Cp_coef_lT[i][1] * pow(T, -1) + phyc.Cp_coef_lT[i][2] * log(T) + phyc.Cp_coef_lT[i][3] * pow(T, 1) + phyc.Cp_coef_lT[i][4] * pow(T, 2) / 2.
        + phyc.Cp_coef_lT[i][5] * pow(T, 3) / 3. + phyc.Cp_coef_lT[i][6] * pow(T, 4) / 4. + phyc.Cp_coef_lT[i][8];
    return Si;
}

double get_Cvi(int component_i, double T)
{
    double Cpi, Cvi;

    Cpi = get_Cpi(component_i, T);
    Cvi = Cpi - phyc.kR / phyc.mol_weight[component_i];
    return Cvi;
}

// Enthalpy of the gas
// Y -- mass fractions
double get_enthalpy(int num_species, double *Y, double T)
{
    double H = 0;
    
    for(int i = 0; i < num_species; i++)
        H += Y[i] * get_Hi(i, T);

    return H;
}

// P = rho * R * T
// R -- gas constant
// Y -- mass fractions
double get_gas_constant(int num_gas_species, double *Y)
{
    // Gas Constant
    double R = 0;

    for(int i = 0; i < num_gas_species; i++)
        R += Y[i] / phyc.mol_weight[i];
    
    R *= phyc.kR;
    
    return R;
}

// Specific heat of the gas
double get_Cp(int num_species, double *Y, double T)
{
    double Cp;
    
    Cp = 0.0;
    for(int i = 0; i < num_species; i++)
        Cp += Y[i] * get_Cpi(i, T);

    return Cp;
}

double get_Cv(int num_species, double *Y, double T)
{
    double Cv;
    
    Cv = 0.0;
    for(int i = 0; i < num_species; i++)
        Cv += Y[i] * get_Cvi(i, T);

    return Cv;
}