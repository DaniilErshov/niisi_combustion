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
            if (Pr_f < 0) Pr_f = 1;
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

    //for (int r_i = 0; r_i < num_react; r_i++)
    //{
    //    cout << "forward " << r_i << " = " << forward[r_i] << "\n";
    //    cout << "reverse " << r_i << " = " << reverse[r_i] << "\n";
    //}
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
    int T_int = (int)round(T);;
    double Dij = 0;
    double specie1, specie2;
    bool flag_Add = 1;

   /* for (const auto& item : Dij_saved)
    {
        if (item.first > T_int - 5 && item.first < T_int + 5)
        {
            flag_Add = 0;
            break;
        }
    }
    
    for (int j = 0; j < num_gas_species; j++) {
        specie1 = j, specie2 = k;

        if (k < j) {
            specie1 = k;
            specie2 = j;
        }

        if(flag_Add)
        {
            Dij = Dij_func(k, j, T, Y);
            Dij_saved[T_int].insert(make_pair(komponents_str[specie1] + " " + komponents_str[specie2], Dij));
        }
        else {
                Dij = Dij_func(k, j, T, Y);
                Dij_saved[T_int].insert(make_pair(komponents_str[specie1] + " " + komponents_str[specie2], Dij));
        }

        if (j != k) {
            sum += Xi[j]
                / Dij;
        }*/

    for (int j = 0; j < num_gas_species; j++) {
        if (j != k) {
            sum += Xi[j]
                / Dij_func(k, j, T, Y);
        }
    }

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
        res = -9.595683577746838 + 1.8730595700658685 * log(T) + -0.028671856143716555 * pow(log(T), 2) + 0.0012484323069923454 * pow(log(T), 3);
    if (i == H2 && j == H)
        res = -10.823402227154023 + 2.505862086122786 * log(T) + -0.10899263274795905 * pow(log(T), 2) + 0.004634086713122515 * pow(log(T), 3);
    if (i == H2 && j == O2)
        res = -12.085938093516864 + 2.644538276663122 * log(T) + -0.13075668402762797 * pow(log(T), 2) + 0.0057369590874055965 * pow(log(T), 3);
    if (i == H2 && j == O)
        res = -10.96445026658662 + 2.31655161587527 * log(T) + -0.08841479655712878 * pow(log(T), 2) + 0.003924864554779722 * pow(log(T), 3);
    if (i == H2 && j == OH)
        res = -10.962295316340866 + 2.3138027128833354 * log(T) + -0.08795138462368206 * pow(log(T), 2) + 0.003899697174905001 * pow(log(T), 3);
    if (i == H2 && j == HO2)
        res = -12.072494702068258 + 2.6379862876294524 * log(T) + -0.12976343474521723 * pow(log(T), 2) + 0.005687349182370516 * pow(log(T), 3);
    if (i == H2 && j == H2O)
        res = -16.978056526030944 + 4.51190545379321 * log(T) + -0.36016652392873316 * pow(log(T), 2) + 0.015132253414787578 * pow(log(T), 3);
    if (i == H2 && j == H2O2)
        res = -12.09285843409546 + 2.64652175728679 * log(T) + -0.13099421965976094 * pow(log(T), 2) + 0.005746154756742783 * pow(log(T), 3);
    if (i == H2 && j == N2)
        res = -12.16602344816801 + 2.66762349776338 * log(T) + -0.1346499533182834 * pow(log(T), 2) + 0.0059448158469462 * pow(log(T), 3);

    if (i == H && j == H2)
        res = -10.823402227154023 + 2.505862086122786 * log(T) + -0.10899263274795905 * pow(log(T), 2) + 0.004634086713122515 * pow(log(T), 3);
    if (i == H && j == H)
        res = -14.090977866177216 + 3.8971936827869396 * log(T) + -0.28434119945936825 * pow(log(T), 2) + 0.012013387141699305 * pow(log(T), 3);
    if (i == H && j == O2)
        res = -15.802031116396073 + 4.252016339625318 * log(T) + -0.33308000643226254 * pow(log(T), 2) + 0.014218477696902987 * pow(log(T), 3);
    if (i == H && j == O)
        res = -13.943114318627917 + 3.664583759907452 * log(T) + -0.25937222598545223 * pow(log(T), 2) + 0.011137922638393416 * pow(log(T), 3);
    if (i == H && j == OH)
        res = -13.946943654278872 + 3.6651134051658776 * log(T) + -0.25938519164491014 * pow(log(T), 2) + 0.011135800704483473 * pow(log(T), 3);
    if (i == H && j == HO2)
        res = -15.811695558573987 + 4.255584766504236 * log(T) + -0.33352963129293595 * pow(log(T), 2) + 0.014237153382348433 * pow(log(T), 3);
    if (i == H && j == H2O)
        res = -17.834461751289805 + 4.803494365646052 * log(T) + -0.36773333759636306 * pow(log(T), 2) + 0.014301259769101626 * pow(log(T), 3);
    if (i == H && j == H2O2)
        res = -15.84909197009192 + 4.2711590675119835 * log(T) + -0.3357011355561007 * pow(log(T), 2) + 0.01433770069442436 * pow(log(T), 3);
    if (i == H && j == N2)
        res = -15.448438512831059 + 4.09830128653488 * log(T) + -0.31340352474286276 * pow(log(T), 2) + 0.013379140840107218 * pow(log(T), 3);

    if (i == O2 && j == H2)
        res = -12.085938093516864 + 2.644538276663122 * log(T) + -0.13075668402762797 * pow(log(T), 2) + 0.0057369590874055965 * pow(log(T), 3);
    if (i == O2 && j == H)
        res = -15.802031116396073 + 4.252016339625318 * log(T) + -0.33308000643226254 * pow(log(T), 2) + 0.014218477696902987 * pow(log(T), 3);
    if (i == O2 && j == O2)
        res = -15.064699548887624 + 3.2511958983607867 * log(T) + -0.20498243765965019 * pow(log(T), 2) + 0.008761094852578755 * pow(log(T), 3);
    if (i == O2 && j == O)
        res = -13.95718719250263 + 3.0031415854371524 * log(T) + -0.17388814880063627 * pow(log(T), 2) + 0.007459071273965279 * pow(log(T), 3);
    if (i == O2 && j == OH)
        res = -13.936329787023418 + 2.9858334987758277 * log(T) + -0.17150429684673577 * pow(log(T), 2) + 0.0073503703523701 * pow(log(T), 3);
    if (i == O2 && j == HO2)
        res = -15.06261482724852 + 3.2467361576031717 * log(T) + -0.204312124046108 * pow(log(T), 2) + 0.00872785428285383 * pow(log(T), 3);
    if (i == O2 && j == H2O)
        res = -20.255740108744483 + 5.149760554420596 * log(T) + -0.4234530364335103 * pow(log(T), 2) + 0.017129963709010876 * pow(log(T), 3);
    if (i == O2 && j == H2O2)
        res = -15.085688569886855 + 3.2535474506198687 * log(T) + -0.20528531845858441 * pow(log(T), 2) + 0.00877386292018964 * pow(log(T), 3);
    if (i == O2 && j == N2)
        res = -14.761930889470714 + 3.1336781023010816 * log(T) + -0.1899979031467854 * pow(log(T), 2) + 0.008123797763855078 * pow(log(T), 3);

    if (i == O && j == H2)
        res = -10.96445026658662 + 2.31655161587527 * log(T) + -0.08841479655712878 * pow(log(T), 2) + 0.003924864554779722 * pow(log(T), 3);
    if (i == O && j == H)
        res = -13.943114318627917 + 3.664583759907452 * log(T) + -0.25937222598545223 * pow(log(T), 2) + 0.011137922638393416 * pow(log(T), 3);
    if (i == O && j == O2)
        res = -13.95718719250263 + 3.0031415854371524 * log(T) + -0.17388814880063627 * pow(log(T), 2) + 0.007459071273965279 * pow(log(T), 3);
    if (i == O && j == O)
        res = -12.496609104456216 + 2.586634464142541 * log(T) + -0.11917232376009453 * pow(log(T), 2) + 0.00506381152216262 * pow(log(T), 3);
    if (i == O && j == OH)
        res = -12.495819336613021 + 2.5796546548115216 * log(T) + -0.11815263499394571 * pow(log(T), 2) + 0.005014464713454417 * pow(log(T), 3);
    if (i == O && j == HO2)
        res = -13.979002942570071 + 3.010110358477627 * log(T) + -0.17483718695172576 * pow(log(T), 2) + 0.007502030323316089 * pow(log(T), 3);
    if (i == O && j == H2O)
        res = -18.327443048276123 + 4.665597880321536 * log(T) + -0.3723444263324898 * pow(log(T), 2) + 0.015363471228681924 * pow(log(T), 3);
    if (i == O && j == H2O2)
        res = -13.980006510799385 + 3.008101642728163 * log(T) + -0.1744884106255608 * pow(log(T), 2) + 0.007482803135234167 * pow(log(T), 3);
    if (i == O && j == N2)
        res = -13.65151744315106 + 2.8745311036219334 * log(T) + -0.15705613899615692 * pow(log(T), 2) + 0.006724808757543971 * pow(log(T), 3);

    if (i == OH && j == H2)
        res = -10.962295316340866 + 2.3138027128833354 * log(T) + -0.08795138462368206 * pow(log(T), 2) + 0.003899697174905001 * pow(log(T), 3);
    if (i == OH && j == H)
        res = -13.946943654278872 + 3.6651134051658776 * log(T) + -0.25938519164491014 * pow(log(T), 2) + 0.011135800704483473 * pow(log(T), 3);
    if (i == OH && j == O2)
        res = -13.936329787023418 + 2.9858334987758277 * log(T) + -0.17150429684673577 * pow(log(T), 2) + 0.0073503703523701 * pow(log(T), 3);
    if (i == OH && j == O)
        res = -12.495819336613021 + 2.5796546548115216 * log(T) + -0.11815263499394571 * pow(log(T), 2) + 0.005014464713454417 * pow(log(T), 3);
    if (i == OH && j == OH)
        res = -12.507695121365083 + 2.5780926387275955 * log(T) + -0.11793212300198015 * pow(log(T), 2) + 0.005004151736241927 * pow(log(T), 3);
    if (i == OH && j == HO2)
        res = -13.965421504425606 + 2.9958547081299414 * log(T) + -0.17288753643194013 * pow(log(T), 2) + 0.007413654243839551 * pow(log(T), 3);
    if (i == OH && j == H2O)
        res = -18.362246558359548 + 4.673921950637226 * log(T) + -0.37355746219005515 * pow(log(T), 2) + 0.015421671202477949 * pow(log(T), 3);
    if (i == OH && j == H2O2)
        res = -14.000324577738853 + 3.008696506402565 * log(T) + -0.1747036358379534 * pow(log(T), 2) + 0.007498873729264166 * pow(log(T), 3);
    if (i == OH && j == N2)
        res = -13.644048892967898 + 2.863264495204685 * log(T) + -0.15551145249301063 * pow(log(T), 2) + 0.006654564410578989 * pow(log(T), 3);

    if (i == HO2 && j == H2)
        res = -12.072494702068258 + 2.6379862876294524 * log(T) + -0.12976343474521723 * pow(log(T), 2) + 0.005687349182370516 * pow(log(T), 3);
    if (i == HO2 && j == H)
        res = -15.811695558573987 + 4.255584766504236 * log(T) + -0.33352963129293595 * pow(log(T), 2) + 0.014237153382348433 * pow(log(T), 3);
    if (i == HO2 && j == O2)
        res = -15.06261482724852 + 3.2467361576031717 * log(T) + -0.204312124046108 * pow(log(T), 2) + 0.00872785428285383 * pow(log(T), 3);
    if (i == HO2 && j == O)
        res = -13.979002942570071 + 3.010110358477627 * log(T) + -0.17483718695172576 * pow(log(T), 2) + 0.007502030323316089 * pow(log(T), 3);
    if (i == HO2 && j == OH)
        res = -13.965421504425606 + 2.9958547081299414 * log(T) + -0.17288753643194013 * pow(log(T), 2) + 0.007413654243839551 * pow(log(T), 3);
    if (i == HO2 && j == HO2)
        res = -15.078958090702335 + 3.2507009546651426 * log(T) + -0.2049167645970397 * pow(log(T), 2) + 0.008758179466114152 * pow(log(T), 3);
    if (i == HO2 && j == H2O)
        res = -19.546316861309116 + 4.91126134610257 * log(T) + -0.39704157484667113 * pow(log(T), 2) + 0.016158555403636643 * pow(log(T), 3);
    if (i == HO2 && j == H2O2)
        res = -15.083506979207977 + 3.2495799335938895 * log(T) + -0.2047782749383369 * pow(log(T), 2) + 0.008752709247197427 * pow(log(T), 3);
    if (i == HO2 && j == N2)
        res = -14.770759062659032 + 3.134600441523139 * log(T) + -0.19015870870343127 * pow(log(T), 2) + 0.008132833232204068 * pow(log(T), 3);

    if (i == H2O && j == H2)
        res = -16.978056526030944 + 4.51190545379321 * log(T) + -0.36016652392873316 * pow(log(T), 2) + 0.015132253414787578 * pow(log(T), 3);
    if (i == H2O && j == H)
        res = -17.834461751289805 + 4.803494365646052 * log(T) + -0.36773333759636306 * pow(log(T), 2) + 0.014301259769101626 * pow(log(T), 3);
    if (i == H2O && j == O2)
        res = -20.255740108744483 + 5.149760554420596 * log(T) + -0.4234530364335103 * pow(log(T), 2) + 0.017129963709010876 * pow(log(T), 3);
    if (i == H2O && j == O)
        res = -18.327443048276123 + 4.665597880321536 * log(T) + -0.3723444263324898 * pow(log(T), 2) + 0.015363471228681924 * pow(log(T), 3);
    if (i == H2O && j == OH)
        res = -18.362246558359548 + 4.673921950637226 * log(T) + -0.37355746219005515 * pow(log(T), 2) + 0.015421671202477949 * pow(log(T), 3);
    if (i == H2O && j == HO2)
        res = -19.546316861309116 + 4.91126134610257 * log(T) + -0.39704157484667113 * pow(log(T), 2) + 0.016158555403636643 * pow(log(T), 3);
    if (i == H2O && j == H2O)
        res = -15.281222624159534 + 2.4316573121063683 * log(T) + 0.019636154929359474 * pow(log(T), 2) + -0.0050564602502904815 * pow(log(T), 3);
    if (i == H2O && j == H2O2)
        res = -19.576228499862463 + 4.921474923914928 * log(T) + -0.3984527651065341 * pow(log(T), 2) + 0.01622419122148494 * pow(log(T), 3);
    if (i == H2O && j == N2)
        res = -20.080205379216675 + 5.097944355010822 * log(T) + -0.4195888814625673 * pow(log(T), 2) + 0.017072929263335423 * pow(log(T), 3);

    if (i == H2O2 && j == H2)
        res = -12.09285843409546 + 2.64652175728679 * log(T) + -0.13099421965976094 * pow(log(T), 2) + 0.005746154756742783 * pow(log(T), 3);
    if (i == H2O2 && j == H)
        res = -15.84909197009192 + 4.2711590675119835 * log(T) + -0.3357011355561007 * pow(log(T), 2) + 0.01433770069442436 * pow(log(T), 3);
    if (i == H2O2 && j == O2)
        res = -15.085688569886855 + 3.2535474506198687 * log(T) + -0.20528531845858441 * pow(log(T), 2) + 0.00877386292018964 * pow(log(T), 3);
    if (i == H2O2 && j == O)
        res = -13.980006510799385 + 3.008101642728163 * log(T) + -0.1744884106255608 * pow(log(T), 2) + 0.007482803135234167 * pow(log(T), 3);
    if (i == H2O2 && j == OH)
        res = -14.000324577738853 + 3.008696506402565 * log(T) + -0.1747036358379534 * pow(log(T), 2) + 0.007498873729264166 * pow(log(T), 3);
    if (i == H2O2 && j == HO2)
        res = -15.083506979207977 + 3.2495799335938895 * log(T) + -0.2047782749383369 * pow(log(T), 2) + 0.008752709247197427 * pow(log(T), 3);
    if (i == H2O2 && j == H2O)
        res = -19.576228499862463 + 4.921474923914928 * log(T) + -0.3984527651065341 * pow(log(T), 2) + 0.01622419122148494 * pow(log(T), 3);
    if (i == H2O2 && j == H2O2)
        res = -15.088262647119528 + 3.248091185948195 * log(T) + -0.2045212256705617 * pow(log(T), 2) + 0.008738270450599328 * pow(log(T), 3);
    if (i == H2O2 && j == N2)
        res = -14.78368817690427 + 3.136927916974964 * log(T) + -0.1904462648473799 * pow(log(T), 2) + 0.008144382463101313 * pow(log(T), 3);

    if (i == N2 && j == H2)
        res = -12.16602344816801 + 2.66762349776338 * log(T) + -0.1346499533182834 * pow(log(T), 2) + 0.0059448158469462 * pow(log(T), 3);
    if (i == N2 && j == H)
        res = -15.448438512831059 + 4.09830128653488 * log(T) + -0.31340352474286276 * pow(log(T), 2) + 0.013379140840107218 * pow(log(T), 3);
    if (i == N2 && j == O2)
        res = -14.761930889470714 + 3.1336781023010816 * log(T) + -0.1899979031467854 * pow(log(T), 2) + 0.008123797763855078 * pow(log(T), 3);
    if (i == N2 && j == O)
        res = -13.65151744315106 + 2.8745311036219334 * log(T) + -0.15705613899615692 * pow(log(T), 2) + 0.006724808757543971 * pow(log(T), 3);
    if (i == N2 && j == OH)
        res = -13.644048892967898 + 2.863264495204685 * log(T) + -0.15551145249301063 * pow(log(T), 2) + 0.006654564410578989 * pow(log(T), 3);
    if (i == N2 && j == HO2)
        res = -14.770759062659032 + 3.134600441523139 * log(T) + -0.19015870870343127 * pow(log(T), 2) + 0.008132833232204068 * pow(log(T), 3);
    if (i == N2 && j == H2O)
        res = -20.080205379216675 + 5.097944355010822 * log(T) + -0.4195888814625673 * pow(log(T), 2) + 0.017072929263335423 * pow(log(T), 3);
    if (i == N2 && j == H2O2)
        res = -14.78368817690427 + 3.136927916974964 * log(T) + -0.1904462648473799 * pow(log(T), 2) + 0.008144382463101313 * pow(log(T), 3);
    if (i == N2 && j == N2)
        res = -14.532106736247874 + 3.047027680004394 * log(T) + -0.17938456304379752 * pow(log(T), 2) + 0.007690373469045775 * pow(log(T), 3);

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