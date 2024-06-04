#include "concentration.h"
extern vector<vector<double>> forward_arr_save;
extern vector<vector<double>> reverse_arr_save;

double my_mol_weight(int k) {
    return phyc.mol_weight[k];
}

void make_averageY(double* Yavg, double* Yi, double* Yipred) {
    for (int i = 0; i < num_gas_species; i++) {
        Yavg[i] = (Yi[i] + Yipred[i]) / 2.;
    }
}

void Get_mole_fr(double* X, double* Y)
{
    double W = 0;
    W = get_W(Y);
    for (int i = 0; i < num_gas_species; i++) {
        X[i] = Y[i] * W / (my_mol_weight(i));
    }
}

void moleFraction_to_massFraction(double* X, double* Y)
{
    double W = 0;
    for (int i = 0; i < num_gas_species; i++) {
        W += my_mol_weight(i) * X[i];
    }
    cout << "W = " << W << "\n";
    for (int i = 0; i < num_gas_species; i++) {
        Y[i] = X[i] * (my_mol_weight(i)) / W;
    }
}

void Get_molar_cons(double* X, double* Y, double T)
{
    double W = 0;
    double rho = get_rho(Y, T);
    for (int i = 0; i < num_gas_species; i++) {
        X[i] = Y[i] * rho / my_mol_weight(i);
    }
}

void get_Y(double Y_H2O, double& Y_H2, double& Y_O2, double Y_N2)
{
    Y_H2 = (1 - Y_H2O - Y_N2) * 1. / 9.;
    Y_O2 = (1 - Y_H2O - Y_N2) * 8. / 9.;
}

double get_W(double* Y)
{
    double W = 0;
    for (int i = 0; i < num_gas_species; i++) {
        W += Y[i] / (my_mol_weight(i));
    }
    return 1. / W;
}

double get_rho(double* Y, double T) {
    return P * get_W(Y) / R / T;
}

double get_GradRho(double* Yi, double* Yinext, double x, double xnext, double Ti, double Tinext) {
    return (get_rho(Yinext, Tinext) - get_rho(Yi, Ti)) / (xnext - x);
}

void get_grad(double* gradX, double* Xi, double* Xinext, double x, double xnext) {
    for (int i = 0; i < num_gas_species; i++) {
        gradX[i] = (Xinext[i] - Xi[i]) / (xnext - x);
    }
}

void chem_vel(double* Sn, double* Hn, double* forward, double* reverse, double* equilib, double Tcurr, double* y, double* yprime, int num_cell) {

    double d = 0.14;
    int k = 0, l = 0;
    double logF_f, logF_core_f, logF_r, logF_core_r;

    double sum1, sum2;
    double Kci;
    double Kpi;
    double dSiR, dHiRT;
    double sumv = 0;
    double k_0_f, k_inf_f, c, m, Pr_f, Fcent, F_f;
    double sum_ThirdBodies;
    bool M_exist = 1;
    auto& species = chec.chemkinReader->species();
    if (update_koeffs == 0 && save_chem_koeffs == 1) {
        for (int i = 0; i < num_react; i++) {
            forward[i] = forward_arr_save[num_cell][i];
            reverse[i] = reverse_arr_save[num_cell][i];
        }
    }
    else {
        //define TROE pressure dependence
        for (int i = 0; i < num_react; i++) {
            auto& Arrhenius = chec.chemkinReader->reactions()[i].getArrhenius();
            /*cout << "chec.chemkinReader->reactions()[i].hasLOW() = " << chec.chemkinReader->reactions()[i].hasLOW()
                << " " << chem.has_low[i] << "\n";*/
            if (!chec.chemkinReader->reactions()[i].hasLOW())
            {
                forward[i] = Arrhenius.A * pow(Tcurr, Arrhenius.n)
                    * exp(-Arrhenius.E * pow(10, -3) / Tcurr / phyc.kRc);
                int jop = 0;
                //cout << "Arrhenius.A = " << Arrhenius.A << " " << chem.Arrh_params[i][0] << "\n";
                //cout << "Arrhenius.A = " << Arrhenius.n << " " << chem.Arrh_params[i][1] << "\n";
                //cout << "Arrhenius.A = " << Arrhenius.E << " " << chem.Arrh_params[i][2] << "\n";
            }
            else {
                auto& Arrhenius_LP = chec.chemkinReader->reactions()[i].getLOW();
                k_inf_f = Arrhenius.A * pow(Tcurr, Arrhenius.n)
                    * exp(-Arrhenius.E * pow(10, -3) / Tcurr / phyc.kRc);

                k_0_f = Arrhenius_LP[0] * pow(Tcurr, Arrhenius_LP[1])
                    * exp(-Arrhenius_LP[2] * pow(10, -3) / Tcurr / phyc.kRc);

                //cout << "ArrheniusLP.A = " << Arrhenius_LP[0] << " " << chem.Arrh_LP_params[i][0] << "\n";
                //cout << "ArrheniusLP.A = " << Arrhenius_LP[1] << " " << chem.Arrh_LP_params[i][1] << "\n";
                //cout << "ArrheniusLP.A = " << Arrhenius_LP[2] << " " << chem.Arrh_LP_params[i][2] << "\n";

                auto& ThirdBodies = chec.chemkinReader->reactions()[i].getThirdBodies();

                M_exist = 1;
                Pr_f = 0;
                for (const auto& thridSpecie : ThirdBodies) {
                    if (thridSpecie.first == "!M")
                        M_exist = 0;
                }
                if (M_exist) {
                    for (const auto& specie : ThirdBodies) {
                        Pr_f += (specie.second - 1) * y[komponents[specie.first]];
                    }
                    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
                        Pr_f += y[k_spec];
                    }
                }

                if (!M_exist) {
                    //cout << "\n\n\nM_exist = " << M_exist << "\n";
                    for (const auto& specie : ThirdBodies) {
                        if (specie.first != "!M") {
                            //cout << "!M = " << specie.first << "\n";
                            Pr_f = y[komponents[specie.first]];
                        }
                    }
                }
                F_f = 1;
                Pr_f *= k_0_f / k_inf_f;
                if (chec.chemkinReader->reactions()[i].hasTROE() && Pr_f > pow(10, -30))
                {
                    const auto& Troe_koeff = chec.chemkinReader->reactions()[i].getTROE();
                    double alpha = Troe_koeff[0];
                    double T3 = Troe_koeff[1];
                    double T1 = Troe_koeff[2];
                    double T2 = Troe_koeff[3];
                    //cout << "alpha= " << alpha << " " << chem.Troe_params[i][0] << "\n";
                    //cout << "T3 = " << T3 << " " << chem.Troe_params[i][1] << "\n";
                    //cout << "T1 = " << T1 << " " << chem.Troe_params[i][2] << "\n";
                    //cout << "T2 = " << T2 << " " << chem.Troe_params[i][3] << "\n";

                    if (T2 == 0)
                        Fcent = (1 - alpha) * exp(-Tcurr / T3) + alpha * exp(-Tcurr / T1);
                    else
                        Fcent = (1 - alpha) * exp(-Tcurr / T3) + alpha * exp(-Tcurr / T1) + exp(-T2 / Tcurr);
                    c = -0.4 - 0.67 * log10(Fcent);
                    m = 0.75 - 1.27 * log10(Fcent);
                    d = 0.14;
                    double logPr_f = log10(Pr_f);
                    logF_core_f = pow((logPr_f + c) / (m - d * (logPr_f + c)), 2);
                    logF_f = log10(Fcent) / (1.0 + logF_core_f);
                    F_f = pow(10, logF_f);
                }


                forward[i] = k_inf_f * (Pr_f / (1 + Pr_f)) * F_f;
            }
        }

        //find reverce constant
        int komponent;
        double koeff_komp;
        double Sn_tmp;
        for (int i = 0; i < num_react; i++) {
            if (chec.chemkinReader->reactions()[i].isReversible()) {
                auto& prod = chec.chemkinReader->reactions()[i].getProducts();
                auto& react = chec.chemkinReader->reactions()[i].getReactants();
                sumv = 0;
                dSiR = 0;
                dHiRT = 0;
                for (int koeff_i = 0; koeff_i < 9; koeff_i++) {
                    Sn[koeff_i] = 0;
                    Hn[koeff_i] = 0;
                }
                for (const auto& iti : prod) {
                    komponent = komponents[iti.first];
                    koeff_komp = iti.second;
                    if (Tcurr >= chec.chemkinReader->species()[komponent].thermo().getTCommon())
                        for (int koeff_i = 0; koeff_i < 9; koeff_i++) {
                            Sn_tmp = phyc.Cp_coef_hT[komponent][koeff_i] * koeff_komp / phyc.kR * my_mol_weight(komponent);
                            Sn[koeff_i] += Sn_tmp;
                            Hn[koeff_i] += Sn_tmp / Tcurr;
                        }
                    else {
                        for (int koeff_i = 0; koeff_i < 9; koeff_i++) {
                            Sn_tmp = phyc.Cp_coef_lT[komponent][koeff_i] * koeff_komp / phyc.kR * my_mol_weight(komponent);
                            Sn[koeff_i] += Sn_tmp;
                            Hn[koeff_i] += Sn_tmp / Tcurr;
                        }
                    }
                    sumv += koeff_komp;
                }

                for (const auto& iti : react) {
                    komponent = komponents[iti.first];
                    koeff_komp = iti.second;
                    if (Tcurr >= chec.chemkinReader->species()[komponent].thermo().getTCommon())
                        for (int koeff_i = 0; koeff_i < 9; koeff_i++) {
                            Sn_tmp = phyc.Cp_coef_hT[komponent][koeff_i] * koeff_komp / phyc.kR * my_mol_weight(komponent);
                            Sn[koeff_i] -= Sn_tmp;
                            Hn[koeff_i] -= Sn_tmp / Tcurr;
                        }
                    else {
                        for (int koeff_i = 0; koeff_i < 9; koeff_i++) {
                            Sn_tmp = phyc.Cp_coef_lT[komponent][koeff_i] * koeff_komp / phyc.kR * my_mol_weight(komponent);
                            Sn[koeff_i] -= Sn_tmp;
                            Hn[koeff_i] -= Sn_tmp / Tcurr;
                        }
                    }
                    sumv -= koeff_komp;
                }
                dHiRT = get_dHiRT(Hn, Tcurr);
                dSiR = get_dSiR(Sn, Tcurr);
                Kpi = exp(dSiR - dHiRT);
                Kci = Kpi * pow(P / phyc.kR / Tcurr, sumv);
                reverse[i] = forward[i] / Kci;
            }
            else {
                //cout << chec.chemkinReader->reactions()[i] << "\n";
                reverse[i] = 0;
            }
            //cout << "i = " << i << "\n";
            //cout << "forward true = " << forward[i] << "\n";
            //cout << "reverse true =  " << reverse[i] << "\n";
            //cout << "\n";
        }
    }
       
    if (update_koeffs == 1 && save_chem_koeffs == 1) {
        for (int i = 0; i < num_react; i++) {
            forward_arr_save[num_cell][i] = forward[i];
            reverse_arr_save[num_cell][i] = reverse[i];
        }
    }


    //find equilib konstant with THIRDBODY
    for (int i = 0; i < num_react; i++) {
        auto& prod = chec.chemkinReader->reactions()[i].getProducts();
        auto& react = chec.chemkinReader->reactions()[i].getReactants();
        auto& ThirdBodies = chec.chemkinReader->reactions()[i].getThirdBodies();

        double forw_slag = 1, rev_slag = 1;
        double mnozh_M = 1;
        for (const auto& iti : react) {
            int sp_i = komponents[iti.first];
            forw_slag *= pow(y[sp_i], iti.second);
        }

        for (const auto& iti : prod) {
            int sp_i = komponents[iti.first];
            rev_slag *= pow(y[sp_i], iti.second);
        }
        //cout << "forw true = " << forw_slag << "\n";
        //cout << "rev_slag true = " << rev_slag << "\n\n";
        equilib[i] = (forward[i] * forw_slag - reverse[i] * rev_slag);
        if (chec.chemkinReader->reactions()[i].hasThirdBody()) {
            //cout << "i = " << i << "\n";
            mnozh_M = 0;
            for (const auto& specie : ThirdBodies) {
                mnozh_M += (specie.second - 1) * y[komponents[specie.first]];
                //cout << specie.first << " = " << specie.second << "\n";
            }
            for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
                mnozh_M += y[k_spec];
                //cout << specie.name() << " = " << 1 << "\n";
            }

            equilib[i] *= mnozh_M;
        }
        //cout << "forward " << i << " = " << forward[i] * forw_slag * mnozh_M << "\n";
        //cout << "reverse " << i << " = " << reverse[i] * rev_slag * mnozh_M << "\n";
        //cout << " equilib " << i << " = " << equilib[i] << "\n\n\n";
    }

    //find totall yprime
    for (int sp_i = 0; sp_i < num_gas_species; sp_i++) {
        yprime[sp_i] = 0;
    }

    for (int r_i = 0; r_i < num_react; r_i++) {
        auto& prod = chec.chemkinReader->reactions()[r_i].getProducts();
        auto& react = chec.chemkinReader->reactions()[r_i].getReactants();
        for (const auto& iti : prod) {
            yprime[komponents[iti.first]] += equilib[r_i] * iti.second;
        }
        for (const auto& iti : react) {
            yprime[komponents[iti.first]] -= equilib[r_i] * iti.second;

        }
    }
}

double YkVk_func(int k, double T, double* Y, double* gradX, double* Xi, double* Y_average) {
    double sum = 0.;
    double W = get_W(Y_average);
    double rho = get_rho(Y_average, T);
    double YkVk = 0;
    double Dkm = 0;
    double Dij = 0;
    double specie1, specie2;


    for (int j = 0; j < num_gas_species; j++) {
        if (j != k) {
            sum += Xi[j]
                / Dij_res[k][j];
        }
    }

    Dkm = (1. - Y_average[k]) / sum;
    //cout << "Dkm for " << name_species[k] << " = " << Dkm << "\n";
    return  -my_mol_weight(k) / W * Dkm * gradX[k];
}


double Dij_func5(int i, int j, double T, int number_cell)
{
    //double res;
    //double logt = log(T);
    //res = diff_polynom[i][j][0] + diff_polynom[i][j][1] * logt + diff_polynom[i][j][2] * logt * logt + diff_polynom[i][j][3] * logt * logt * logt
    //    + diff_polynom[i][j][4] * logt * logt * logt * logt;
    ////cout << "res = " << res << "\n";
    //return  res * pow(T, 1.5) / pow(10, 1);

    if (!flag_use_save_koeffs) {
        double res;
        double logt = log(T);
        res = diff_polynom[i][j][0] + diff_polynom[i][j][1] * logt + diff_polynom[i][j][2] * logt * logt + diff_polynom[i][j][3] * logt * logt * logt;
        //cout << "res = " << res << "\n";
        return  exp(res) / P / 100.;
    }
    else {
        return Dij_arr[number_cell][i][j];
    }
}