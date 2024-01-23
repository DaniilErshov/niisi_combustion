#include "concentration.h"

double my_mol_weight(int k) {
    //cout << "mol_weight = " << phyc.mol_weight[k] << "\n";
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
    //cout << "W = " << W << "\n";
    for (int i = 0; i < num_gas_species; i++) {
        //cout << "mol_weight =  " << phyc.mol_weight[i] << endl;
        X[i] = Y[i] * W / (my_mol_weight(i));
    }
}

void Get_molar_cons(double* X, double* Y, double T)
{
    double W = 0;
    double rho2 = get_rho(Y, T);
    double rho = P / get_gas_constant(num_gas_species, Y) / T;
    for (int i = 0; i < num_gas_species; i++) {
        //cout << "mol_weight =  " << my_mol_weight(i) << endl;
        //X[i] = Y[i] * rho2 / my_mol_weight(i) * pow(10, 6);
        X[i] = Y[i] * rho2 / my_mol_weight(i);
        //if (Y[i] < pow(10, -20)) X[i] = 0;
        //cout << "X[i] =  " << X[i] << endl;
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

void add_toChemVel(double* wk_add, double M, double* Yi, double* Yinext, double x, double xnext, double Ti, double Tinext) {
    double grad_rho = get_GradRho(Yi, Yinext, x, xnext, Ti, Tinext);
    double rho = get_rho(Yi, Ti);
    for (int i = 0; i < num_gas_species; i++) {
        wk_add[i] = grad_rho * Yi[i] / (my_mol_weight(i)) * M / rho;
        //cout << "wk_add = " << wk_add[i] << "\n";
    }
}

//here send molar_cons
void chem_vel(double* forward, double* reverse, double* equilib, double Tcurr, double* y, double* yprime) {

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

    //define TROE pressure dependence
    for (int i = 0; i < num_react; i++) {
        auto& Arrhenius = chec.chemkinReader->reactions()[i].getArrhenius();
        if (!chec.chemkinReader->reactions()[i].hasLOW())
            forward[i] = Arrhenius.A * pow(Tcurr, Arrhenius.n)
            * exp(-Arrhenius.E * pow(10, -3) / Tcurr / phyc.kRc);
        else {
            auto& Arrhenius_LP = chec.chemkinReader->reactions()[i].getLOW();

            k_inf_f = Arrhenius.A * pow(Tcurr, Arrhenius.n)
                * exp(-Arrhenius.E * pow(10, -3) / Tcurr / phyc.kRc);

            k_0_f = Arrhenius_LP[0] * pow(Tcurr, Arrhenius_LP[1])
                * exp(-Arrhenius_LP[2] * pow(10, -3) / Tcurr / phyc.kRc);

            //cout << "Arrhenius_LP[0] = " << Arrhenius_LP[0] << "\n";
            //cout << "Arrhenius_LP[1] = " << Arrhenius_LP[1] << "\n";
            //cout << "Arrhenius_LP[2] = " << Arrhenius_LP[2] << "\n";

            Fcent = chec.chemkinReader->reactions()[i].getTROE_center()[0];
            c = -0.4 - 0.67 * log10(Fcent);
            m = 0.75 - 1.27 * log10(Fcent);

            auto& ThirdBodies = chec.chemkinReader->reactions()[i].getThirdBodies();

            M_exist = 1;
            Pr_f = 0;
            for (const auto& thridSpecie : ThirdBodies) {
                if (thridSpecie.first == "!M")
                    M_exist = 0;
            }

            if (M_exist) {
                for (const auto& specie : ThirdBodies) {
                    Pr_f += specie.second * y[komponents[specie.first]];
                    /*cout << "in THRID specie.first = " << specie.first << "\n";
                    cout << "in THRID specie.second = " << specie.second << "\n";*/
                }
                for (const auto& specie : species) {

                    if (!ThirdBodies.contains(specie.name())) {
                        //cout << "AFTER THRID specie.first = " << specie.name() << "\n";
                        Pr_f += y[komponents[specie.name()]];
                    }
                }
            }

            if (!M_exist) {
                for (const auto& specie : ThirdBodies) {
                    if (specie.first != "!M") {
                        //cout << "!M = " << specie.first << "\n";
                        Pr_f = y[komponents[specie.first]];
                    }
                }
            }

            Pr_f *= k_0_f / k_inf_f;
            logF_core_f = pow((log10(Pr_f) + c) / (m - d * (log10(Pr_f) + c)), 2);
            logF_f = pow(1.0 + logF_core_f, -1) * log10(Fcent);
            F_f = pow(10, logF_f);
            forward[i] = k_inf_f * (Pr_f / (1 + Pr_f)) * 1.;
        }
    }

    //find reverce constant
    for (int i = 0; i < num_react; i++) {
        auto& prod = chec.chemkinReader->reactions()[i].getProducts();
        auto& react = chec.chemkinReader->reactions()[i].getReactants();
        sumv = 0;
        dSiR = 0;
        dHiRT = 0;
        for (const auto& iti : prod) {
            dSiR += myget_Si(komponents[iti.first], Tcurr) * iti.second / phyc.kR * my_mol_weight(komponents[iti.first]);
            dHiRT += myget_Hi(komponents[iti.first], Tcurr) * iti.second / phyc.kR / Tcurr * my_mol_weight(komponents[iti.first]);
            sumv += iti.second;
        }
        for (const auto& iti : react) {
            dSiR -= myget_Si(komponents[iti.first], Tcurr) * iti.second / phyc.kR * my_mol_weight(komponents[iti.first]);
            dHiRT -= myget_Hi(komponents[iti.first], Tcurr) * iti.second / phyc.kR / Tcurr * my_mol_weight(komponents[iti.first]);
            sumv -= iti.second;
        }
        Kpi = exp(dSiR - dHiRT);
        Kci = Kpi * pow(P / phyc.kR / Tcurr, sumv);
        reverse[i] = forward[i] / Kci;
    }


    //find equilib konstant with THIRDBODY
    for (int i = 0; i < num_react; i++) {
        auto& prod = chec.chemkinReader->reactions()[i].getProducts();
        auto& react = chec.chemkinReader->reactions()[i].getReactants();
        auto& ThirdBodies = chec.chemkinReader->reactions()[i].getThirdBodies();

        double forw_slag = 1, rev_slag = 1;
        double mnozh_M = 0;

        for (const auto& iti : react) {
            forw_slag *= pow(y[komponents[iti.first]], iti.second);
        }

        for (const auto& iti : prod) {
            rev_slag *= pow(y[komponents[iti.first]], iti.second);
        }

        equilib[i] = (forward[i] * forw_slag - reverse[i] * rev_slag);


        if (chec.chemkinReader->reactions()[i].hasThirdBody()) {

            for (const auto& specie : ThirdBodies) {
                mnozh_M += specie.second * y[komponents[specie.first]];
            }

            for (const auto& specie : species) {
                if (!ThirdBodies.contains(specie.name()))
                    mnozh_M += y[komponents[specie.name()]];
            }
            equilib[i] *= mnozh_M;
        }
    }

    //find totall yprime
    for (int sp_i = 0; sp_i < num_gas_species; sp_i++) {
        yprime[sp_i] = 0;
    }
    for (int r_i = 0; r_i < num_react; r_i++) {
        auto& prod = chec.chemkinReader->reactions()[r_i].getProducts();
        auto& react = chec.chemkinReader->reactions()[r_i].getReactants();
        for (const auto& iti : prod) {
            yprime[komponents[iti.first]] += equilib[r_i];
        }
        for (const auto& iti : react) {
            yprime[komponents[iti.first]] -= equilib[r_i];
        }
    }
    return;
}

double YkVk(int k, double T, double* Y, double* gradX, double* Xi) {
    double sum = 0.;
    double W = get_W(Y);
    double rho = get_rho(Y, T);
    double YkVk = 0;
    double Dkm = 0;
    for (int j = 0; j < num_gas_species; j++) {
        if (j != k) {
            sum += Xi[j]
                / Dij_func(k, j, T, Y);
        }
    }
    //cout << "T = " << T << "\n";
    Dkm = (1. - Y[k]) / sum;
    //cout << "Dkm for " << name_species[k] << " = " << Dkm << "\n";
    return  -my_mol_weight(k) / W * Dkm * gradX[k];
}

double Dij_func(int i, int j, double T, double* Y)
{
    double res;
    //double D_res = 3. / 16. * pow(2 * 3.14156 * pow(kB, 3), 0.5) / 3.14156;
    //cout << "konst = " << D_res << "\n";
    //double D_res = 1.858 * pow(10, -3);
    //double mij = my_mol_weight(i) * my_mol_weight(j) / (my_mol_weight(i) + my_mol_weight(j));
    //double rho = get_rho(Y, T);
    //double sigma_ij = 0.5 * (chec.chemkinReader->species()[i].transport().getCollisionDiameter() +
    //    chec.chemkinReader->species()[j].transport().getCollisionDiameter()) / Angstroem__;
    //double eps_ij = pow(chec.chemkinReader->species()[i].transport().getPotentialWellDepth() *
    //    chec.chemkinReader->species()[j].transport().getPotentialWellDepth(), 0.5);
    //double T_ij = kB * T / eps_ij;
    //double integ_ij = 1.074 * pow(T_ij, -0.1604);
    //D_res *= pow(pow(T, 3) / mij, 0.5) / 1. / pow(sigma_ij, 2.) / integ_ij;

    /*cout << "rho = " << rho << "\n";
    cout << "sigmaij = " << sigma_ij << "\n";
    cout << "T_ij = " << T_ij << "\n";*/
    //1.033227
    if (i == H2 && j == H2)
        res = -5.511146184707643 + 0.2539781092813034 * log(T) + 0.18583890855037175 * pow(log(T), 2) + -0.008241706344815997 * pow(log(T), 3);
    if (i == H2 && j == H)
        res = -8.9362260290272 + 1.8410115366188027 * log(T) + -0.03094220816699436 * pow(log(T), 2) + 0.0015796496465732175 * pow(log(T), 3);
    if (i == H2 && j == O2)
        res = -14.310085128291231 + 3.6690223700606657 * log(T) + -0.2839795306356031 * pow(log(T), 2) + 0.013233578082684772 * pow(log(T), 3);
    if (i == H2 && j == O)
        res = -10.147543764389011 + 2.017699491897917 * log(T) + -0.051274947561512144 * pow(log(T), 2) + 0.0023562952003273463 * pow(log(T), 3);
    if (i == H2 && j == OH)
        res = -8.70312790010006 + 1.4615239004667344 * log(T) + 0.01950960336373571 * pow(log(T), 2) + -0.0006262939959915427 * pow(log(T), 3);
    if (i == H2 && j == HO2)
        res = -1.7782964142459061 + -1.640805929718863 * log(T) + 0.46241422422254774 * pow(log(T), 2) + -0.021592849128565184 * pow(log(T), 3);
    if (i == H2 && j == H2O)
        res = -13.016391816092835 + 2.930714012335867 * log(T) + -0.14903413139106828 * pow(log(T), 2) + 0.005707134424530025 * pow(log(T), 3);
    if (i == H2 && j == H2O2)
        res = -6.3083313652087165 + 0.2706940240291075 * log(T) + 0.19445068895086678 * pow(log(T), 2) + -0.009111741098882326 * pow(log(T), 3);
    if (i == H2 && j == N2)
        res = -7.905566719717262 + 0.9462210268795155 * log(T) + 0.09784581392968299 * pow(log(T), 2) + -0.004543834265591781 * pow(log(T), 3);

    if (i == H && j == H2)
        res = -8.9362260290272 + 1.8410115366188027 * log(T) + -0.03094220816699436 * pow(log(T), 2) + 0.0015796496465732175 * pow(log(T), 3);
    if (i == H && j == H)
        res = -4.929754835901546 + 0.1489918551651856 * log(T) + 0.2271857817045916 * pow(log(T), 2) + -0.011264700551599782 * pow(log(T), 3);
    if (i == H && j == O2)
        res = -14.350551114195508 + 3.7100222071029965 * log(T) + -0.2653292086031062 * pow(log(T), 2) + 0.01138788157521553 * pow(log(T), 3);
    if (i == H && j == O)
        res = -9.38520675718858 + 1.8553598026009885 * log(T) + -0.01960828550034363 * pow(log(T), 2) + 0.0005335104136760319 * pow(log(T), 3);
    if (i == H && j == OH)
        res = -8.766209786392151 + 1.5677664329219179 * log(T) + 0.024323748588106164 * pow(log(T), 2) + -0.001675495871676362 * pow(log(T), 3);
    if (i == H && j == HO2)
        res = -12.354062968202447 + 2.906750666170543 * log(T) + -0.157862666217465 * pow(log(T), 2) + 0.006606333855315366 * pow(log(T), 3);
    if (i == H && j == H2O)
        res = -9.511886880913618 + 1.4891063220624936 * log(T) + 0.07215635902790739 * pow(log(T), 2) + -0.005150198332022404 * pow(log(T), 3);
    if (i == H && j == H2O2)
        res = -10.54094712456655 + 2.1320902669330657 * log(T) + -0.04793919160555328 * pow(log(T), 2) + 0.0014266659754885724 * pow(log(T), 3);
    if (i == H && j == N2)
        res = -13.38432100267693 + 3.3244669108028604 * log(T) + -0.21649370218250596 * pow(log(T), 2) + 0.009331507226679664 * pow(log(T), 3);

    if (i == O2 && j == H2)
        res = -14.310085128291231 + 3.6690223700606657 * log(T) + -0.2839795306356031 * pow(log(T), 2) + 0.013233578082684772 * pow(log(T), 3);
    if (i == O2 && j == H)
        res = -14.350551114195508 + 3.7100222071029965 * log(T) + -0.2653292086031062 * pow(log(T), 2) + 0.01138788157521553 * pow(log(T), 3);
    if (i == O2 && j == O2)
        res = -15.374954056917934 + 3.454551196322449 * log(T) + -0.24197532594327462 * pow(log(T), 2) + 0.010813083178420113 * pow(log(T), 3);
    if (i == O2 && j == O)
        res = -5.0548042183230715 + -0.6638248591011925 * log(T) + 0.3295643906182689 * pow(log(T), 2) + -0.015569375382666562 * pow(log(T), 3);
    if (i == O2 && j == OH)
        res = -9.089627987048317 + 1.0198377327256425 * log(T) + 0.09489064642890749 * pow(log(T), 2) + -0.004700528828547379 * pow(log(T), 3);
    if (i == O2 && j == HO2)
        res = -14.447967931537217 + 3.1072835534477634 * log(T) + -0.19926055394469444 * pow(log(T), 2) + 0.009078842248570751 * pow(log(T), 3);
    if (i == O2 && j == H2O)
        res = -20.44353658899271 + 5.366243873443004 * log(T) + -0.4705321692852626 * pow(log(T), 2) + 0.020004375392573104 * pow(log(T), 3);
    if (i == O2 && j == H2O2)
        res = -12.891362690881216 + 2.4386679793380304 * log(T) + -0.10367637372857554 * pow(log(T), 2) + 0.004521099254508942 * pow(log(T), 3);
    if (i == O2 && j == N2)
        res = -11.015455271431685 + 1.5958431311836634 * log(T) + 0.021143926753465363 * pow(log(T), 2) + -0.0015605025179286983 * pow(log(T), 3);

    if (i == O && j == H2)
        res = -10.147543764389011 + 2.017699491897917 * log(T) + -0.051274947561512144 * pow(log(T), 2) + 0.0023562952003273463 * pow(log(T), 3);
    if (i == O && j == H)
        res = -9.38520675718858 + 1.8553598026009885 * log(T) + -0.01960828550034363 * pow(log(T), 2) + 0.0005335104136760319 * pow(log(T), 3);
    if (i == O && j == O2)
        res = -5.0548042183230715 + -0.6638248591011925 * log(T) + 0.3295643906182689 * pow(log(T), 2) + -0.015569375382666562 * pow(log(T), 3);
    if (i == O && j == O)
        res = 1.7049569953091948 + -0.6352900729999933 * log(T) + 0.046768145619869916 * pow(log(T), 2) + 0.006954832942134649 * pow(log(T), 3);
    if (i == O && j == OH)
        res = -6.6698176135505225 + 0.20299018148130044 * log(T) + 0.20570619144344054 * pow(log(T), 2) + -0.009713322976378946 * pow(log(T), 3);
    if (i == O && j == HO2)
        res = -14.484450387286946 + 3.323380067866044 * log(T) + -0.23046418934262425 * pow(log(T), 2) + 0.010555189014602645 * pow(log(T), 3);
    if (i == O && j == H2O)
        res = -11.446189422928509 + 1.8513912862051816 * log(T) + 0.011510555733839906 * pow(log(T), 2) + -0.002087083387809654 * pow(log(T), 3);
    if (i == O && j == H2O2)
        res = -14.747776943808258 + 3.4134174882546118 * log(T) + -0.24097858514496562 * pow(log(T), 2) + 0.010966034441367012 * pow(log(T), 3);
    if (i == O && j == N2)
        res = -9.74988434012846 + 1.3160847258148212 * log(T) + 0.05129180803153138 * pow(log(T), 2) + -0.002587708555491751 * pow(log(T), 3);

    if (i == OH && j == H2)
        res = -8.70312790010006 + 1.4615239004667344 * log(T) + 0.01950960336373571 * pow(log(T), 2) + -0.0006262939959915427 * pow(log(T), 3);
    if (i == OH && j == H)
        res = -8.766209786392151 + 1.5677664329219179 * log(T) + 0.024323748588106164 * pow(log(T), 2) + -0.001675495871676362 * pow(log(T), 3);
    if (i == OH && j == O2)
        res = -9.089627987048317 + 1.0198377327256425 * log(T) + 0.09489064642890749 * pow(log(T), 2) + -0.004700528828547379 * pow(log(T), 3);
    if (i == OH && j == O)
        res = -6.6698176135505225 + 0.20299018148130044 * log(T) + 0.20570619144344054 * pow(log(T), 2) + -0.009713322976378946 * pow(log(T), 3);
    if (i == OH && j == OH)
        res = -11.762500496673484 + 2.362563252860093 * log(T) + -0.09929467159779261 * pow(log(T), 2) + 0.004593329624206208 * pow(log(T), 3);
    if (i == OH && j == HO2)
        res = -15.112559943913636 + 3.564706992237448 * log(T) + -0.2623070237068407 * pow(log(T), 2) + 0.011949237563494641 * pow(log(T), 3);
    if (i == OH && j == H2O)
        res = -15.526477331225582 + 3.5994313791259347 * log(T) + -0.23785803173936182 * pow(log(T), 2) + 0.00971512631753274 * pow(log(T), 3);
    if (i == OH && j == H2O2)
        res = -9.145527816800705 + 1.055554163916615 * log(T) + 0.08814127782415068 * pow(log(T), 2) + -0.00432270944915418 * pow(log(T), 3);
    if (i == OH && j == N2)
        res = -8.513593083331655 + 0.7689792010800627 * log(T) + 0.12994817337380488 * pow(log(T), 2) + -0.006323126550097702 * pow(log(T), 3);

    if (i == HO2 && j == H2)
        res = -1.7782964142459061 + -1.640805929718863 * log(T) + 0.46241422422254774 * pow(log(T), 2) + -0.021592849128565184 * pow(log(T), 3);
    if (i == HO2 && j == H)
        res = -12.354062968202447 + 2.906750666170543 * log(T) + -0.157862666217465 * pow(log(T), 2) + 0.006606333855315366 * pow(log(T), 3);
    if (i == HO2 && j == O2)
        res = -14.447967931537217 + 3.1072835534477634 * log(T) + -0.19926055394469444 * pow(log(T), 2) + 0.009078842248570751 * pow(log(T), 3);
    if (i == HO2 && j == O)
        res = -14.484450387286946 + 3.323380067866044 * log(T) + -0.23046418934262425 * pow(log(T), 2) + 0.010555189014602645 * pow(log(T), 3);
    if (i == HO2 && j == OH)
        res = -15.112559943913636 + 3.564706992237448 * log(T) + -0.2623070237068407 * pow(log(T), 2) + 0.011949237563494641 * pow(log(T), 3);
    if (i == HO2 && j == HO2)
        res = -6.396994466262594 + -0.3137021403151278 * log(T) + 0.2832900222275383 * pow(log(T), 2) + -0.013534288538253558 * pow(log(T), 3);
    if (i == HO2 && j == H2O)
        res = -16.904499288765518 + 3.871094276993996 * log(T) + -0.2601732345291718 * pow(log(T), 2) + 0.010155754345279155 * pow(log(T), 3);
    if (i == HO2 && j == H2O2)
        res = -15.18333918096658 + 3.435052340049357 * log(T) + -0.24731709560900011 * pow(log(T), 2) + 0.011364400992810081 * pow(log(T), 3);
    if (i == HO2 && j == N2)
        res = -5.629772768658993 + -0.5820157824311425 * log(T) + 0.3138295751666163 * pow(log(T), 2) + -0.014656725272155391 * pow(log(T), 3);

    if (i == H2O && j == H2)
        res = -13.016391816092835 + 2.930714012335867 * log(T) + -0.14903413139106828 * pow(log(T), 2) + 0.005707134424530025 * pow(log(T), 3);
    if (i == H2O && j == H)
        res = -9.511886880913618 + 1.4891063220624936 * log(T) + 0.07215635902790739 * pow(log(T), 2) + -0.005150198332022404 * pow(log(T), 3);
    if (i == H2O && j == O2)
        res = -20.44353658899271 + 5.366243873443004 * log(T) + -0.4705321692852626 * pow(log(T), 2) + 0.020004375392573104 * pow(log(T), 3);
    if (i == H2O && j == O)
        res = -11.446189422928509 + 1.8513912862051816 * log(T) + 0.011510555733839906 * pow(log(T), 2) + -0.002087083387809654 * pow(log(T), 3);
    if (i == H2O && j == OH)
        res = -15.526477331225582 + 3.5994313791259347 * log(T) + -0.23785803173936182 * pow(log(T), 2) + 0.00971512631753274 * pow(log(T), 3);
    if (i == H2O && j == HO2)
        res = -16.904499288765518 + 3.871094276993996 * log(T) + -0.2601732345291718 * pow(log(T), 2) + 0.010155754345279155 * pow(log(T), 3);
    if (i == H2O && j == H2O)
        res = -10.790827417374365 + 0.7708550997958658 * log(T) + 0.22306913882903492 * pow(log(T), 2) + -0.013299383714021812 * pow(log(T), 3);
    if (i == H2O && j == H2O2)
        res = -9.464659246290008 + 0.7661161031750491 * log(T) + 0.1704778481359649 * pow(log(T), 2) + -0.009717060407461148 * pow(log(T), 3);
    if (i == H2O && j == N2)
        res = -11.373998705334952 + 1.5419260905886307 * log(T) + 0.06410434615882055 * pow(log(T), 2) + -0.004832886276023816 * pow(log(T), 3);

    if (i == H2O2 && j == H2)
        res = -6.3083313652087165 + 0.2706940240291075 * log(T) + 0.19445068895086678 * pow(log(T), 2) + -0.009111741098882326 * pow(log(T), 3);
    if (i == H2O2 && j == H)
        res = -10.54094712456655 + 2.1320902669330657 * log(T) + -0.04793919160555328 * pow(log(T), 2) + 0.0014266659754885724 * pow(log(T), 3);
    if (i == H2O2 && j == O2)
        res = -12.891362690881216 + 2.4386679793380304 * log(T) + -0.10367637372857554 * pow(log(T), 2) + 0.004521099254508942 * pow(log(T), 3);
    if (i == H2O2 && j == O)
        res = -14.747776943808258 + 3.4134174882546118 * log(T) + -0.24097858514496562 * pow(log(T), 2) + 0.010966034441367012 * pow(log(T), 3);
    if (i == H2O2 && j == OH)
        res = -9.145527816800705 + 1.055554163916615 * log(T) + 0.08814127782415068 * pow(log(T), 2) + -0.00432270944915418 * pow(log(T), 3);
    if (i == H2O2 && j == HO2)
        res = -15.18333918096658 + 3.435052340049357 * log(T) + -0.24731709560900011 * pow(log(T), 2) + 0.011364400992810081 * pow(log(T), 3);
    if (i == H2O2 && j == H2O)
        res = -9.464659246290008 + 0.7661161031750491 * log(T) + 0.1704778481359649 * pow(log(T), 2) + -0.009717060407461148 * pow(log(T), 3);
    if (i == H2O2 && j == H2O2)
        res = -11.443137641163506 + 1.8419122623964752 * log(T) + -0.02240827682289728 * pow(log(T), 2) + 0.0008291037065618657 * pow(log(T), 3);
    if (i == H2O2 && j == N2)
        res = -8.910408124794406 + 0.7929047212992318 * log(T) + 0.12225616772234621 * pow(log(T), 2) + -0.005788832928797144 * pow(log(T), 3);

    if (i == N2 && j == H2)
        res = -7.905566719717262 + 0.9462210268795155 * log(T) + 0.09784581392968299 * pow(log(T), 2) + -0.004543834265591781 * pow(log(T), 3);
    if (i == N2 && j == H)
        res = -13.38432100267693 + 3.3244669108028604 * log(T) + -0.21649370218250596 * pow(log(T), 2) + 0.009331507226679664 * pow(log(T), 3);
    if (i == N2 && j == O2)
        res = -11.015455271431685 + 1.5958431311836634 * log(T) + 0.021143926753465363 * pow(log(T), 2) + -0.0015605025179286983 * pow(log(T), 3);
    if (i == N2 && j == O)
        res = -9.74988434012846 + 1.3160847258148212 * log(T) + 0.05129180803153138 * pow(log(T), 2) + -0.002587708555491751 * pow(log(T), 3);
    if (i == N2 && j == OH)
        res = -8.513593083331655 + 0.7689792010800627 * log(T) + 0.12994817337380488 * pow(log(T), 2) + -0.006323126550097702 * pow(log(T), 3);
    if (i == N2 && j == HO2)
        res = -5.629772768658993 + -0.5820157824311425 * log(T) + 0.3138295751666163 * pow(log(T), 2) + -0.014656725272155391 * pow(log(T), 3);
    if (i == N2 && j == H2O)
        res = -11.373998705334952 + 1.5419260905886307 * log(T) + 0.06410434615882055 * pow(log(T), 2) + -0.004832886276023816 * pow(log(T), 3);
    if (i == N2 && j == H2O2)
        res = -8.910408124794406 + 0.7929047212992318 * log(T) + 0.12225616772234621 * pow(log(T), 2) + -0.005788832928797144 * pow(log(T), 3);
    if (i == N2 && j == N2)
        res = -14.866882292912967 + 3.294249518757591 * log(T) + -0.22632810753009533 * pow(log(T), 2) + 0.010351803087595868 * pow(log(T), 3);


    return  exp(res);
}

double Dk_func(int i, double T, double* Y, double* X, int N) {
    double sum = 0;
    for (int j = 0; j < num_gas_species; j++) {
        if (i != j) sum += X[j] / Dij_func(i, j, T, Y);
        cout << "Dij = " << Dij_func(i, j, T, Y) << "\n";
    }
    return (1 - X[i]) / sum;
}