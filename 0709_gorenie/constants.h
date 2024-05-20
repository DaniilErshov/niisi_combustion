#pragma once

#include "chemkinReader.h"
#include "GasTransport.h"
void init_consts(int& num_gas_species, int& num_react);
double get_Hi(int component_i, double T, int number_cell);
double get_Cpi(int component_i, double T, int number_cell);
double get_Cvi(int component_i, double T, int number_cell);
double get_Cp(int num_species, double* Y, double T, int number_cell);
double get_enthalpy(int num_species, double *Y, double T, int number_cell);
double get_gas_constant(int num_gas_species, double *Y);
double get_Cv(int num_species, double* Y, double T, int number_cell);
double get_Si(int component_i, double T);
double get_dHiRT(double* Cp_coef, double T);
double get_dSiR(double* Cp_coef, double T);
double get_dCpi(double* Cp_coef, double T);
double get_Lambda(int i, double T, int number_cell, char side);
double get_Lambda5(int i, double T, int number_cell, char side);
void allocate_memory();
void free_memory();

struct phy_consts
{
    double kR;  // universal gas constant ( erg/K/mole )
    double kRc; // universal gas constant ( kcal/K/mole )
    double kTime;     // dimensional time ( sec )
    double kLength;   // dimensional length ( cm )
    double kMass;     // dimensional mass ( g )
    double kTemp;      // dimensional temperature ( K )
    double kPres;     // dimensional pressure ( dyn/cm**2 )
    double kDens; // dimensional density ( g/cm**3 )
    double kVel;           // dimensional velocity ( cm/sec )

    // Coefficients for specific heat Cp and
    // molar enthalpy NASA polynoms
    // 8th coef. is for enthalpy calculation
    double** Cp_coef_lT;
    double** Cp_coef_hT;

    // Molar weights
    double* mol_weight;
};

struct che_consts
{
    double* sum_v;
    const std::string chemfile;
    const std::string thermfile;
    const std::string transfile;
    IO::ChemkinReader* chemkinReader;
};

struct chem_struct
{
    double** products;
    double** Arrh_params;
    bool* isReversible;

    bool* has_Third;
    double ** ThirdBodies;
    bool* has_low;
    double** Arrh_LP_params;

    bool* has_Troe;
    double** Troe_params;
    double* T_comon;
    bool* M_exist;


};

extern phy_consts phyc;
extern che_consts chec;