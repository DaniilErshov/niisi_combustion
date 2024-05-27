#pragma once

#include "chemkinReader.h"
#include "GasTransport.h"
void init_consts(int& num_gas_species, int& num_react);
void allocate_memory();
void free_memory();
void set_polynom(double** ploynom, std::string name_file, std::string type_polynom);

void set_polynom_diffusion(double*** polynom, std::string name_file, std::string type_polynom);

std::vector<std::string> splitString(std::string str, char splitter);

template <typename T>
void findValue(const std::vector<T>& data, bool(*condition)(T));

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

extern phy_consts phyc;
extern che_consts chec;