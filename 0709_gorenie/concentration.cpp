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
void chem_vel(double* forward, double* reverse, double* equilib, double* wk_add, double M, double* Yi, double* Yinext, double x, double xnext, double Tcurr, double Tinext, double* y, double* yprime) {
    double k_0_f[3], k_inf_f[3], k_0_r[3], k_inf_r[3];
    double c[3], m[3], d = 0.14;
    double Pr_f[3], Pr_r[3];
    int k = 0, l = 0;
    double logF_f, logF_core_f, logF_r, logF_core_r;

    double sum1, sum2;
    double Kci;
    double Kpi;
    double dSiR, dHiRT;
    double sumv = 0;
    for (int i = 0; i < num_react; i++) {
        if (i != 8 && i != 15 && i != 16) {
            forward[i] = chec.kPrex_f[i] * pow(Tcurr, chec.kPow_f[i])
                * exp(-chec.kE_f[i] / Tcurr / phyc.kRc);
            reverse[i] = chec.kPrex_r[i] * pow(Tcurr, chec.kPow_r[i])
                * exp(-chec.kE_r[i] / Tcurr / phyc.kRc);
        }
        else {
            if (i == 8) k = 0;
            if (i == 15) k = 1;
            if (i == 16) k = 2;

            k_inf_f[k] = chec.kPrex_f[i] * pow(Tcurr, chec.kPow_f[i])
                * exp(-chec.kE_f[i] / Tcurr / phyc.kRc);
            k_inf_r[k] = chec.kPrex_r[i] * pow(Tcurr, chec.kPow_r[i])
                * exp(-chec.kE_r[i] / Tcurr / phyc.kRc);
            k_0_f[k] = chec.kPrex_f_lp[i] * pow(Tcurr, chec.kPow_f_lp[i])
                * exp(-chec.kE_f_lp[i] / Tcurr / phyc.kRc);
            k_0_r[k] = chec.kPrex_r_lp[i] * pow(Tcurr, chec.kPow_r_lp[i])
                * exp(-chec.kE_r_lp[i] / Tcurr / phyc.kRc);

            c[k] = -0.4 - 0.67 * log10(chec.Fcent[i]);
            m[k] = 0.75 - 1.27 * log10(chec.Fcent[i]);
        }
    }

    Pr_f[0] = (k_0_f[0] * (1.3 * y[0] + 10 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8])) / k_inf_f[0];
    Pr_r[0] = (k_0_r[0] * (1.3 * y[0] + 10 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8])) / k_inf_r[0];
    Pr_f[1] = (k_0_f[1] * (3.7 * y[0] + 1.5 * y[8] + 1.2 * y[2] + 7.7 * y[7] + y[1] + y[3] + y[4] + y[5] + y[6])) / k_inf_f[1];
    Pr_r[1] = (k_0_r[1] * (3.7 * y[0] + 1.5 * y[8] + 1.2 * y[2] + 7.7 * y[7] + y[1] + y[3] + y[4] + y[5] + y[6])) / k_inf_r[1];
    Pr_f[2] = (k_0_f[2] * y[6]) / k_inf_f[2];
    Pr_r[2] = (k_0_r[2] * y[6]) / k_inf_r[2];

    for (k = 0; k < 3; k++) {
        if (k == 0) l = 8;
        if (k == 1) l = 15;
        if (k == 2) l = 16;

        if (Pr_f[k] == 0) forward[l] = k_inf_f[k];
        else {
            logF_core_f = pow((log10(Pr_f[k]) + c[k]) / (m[k] - d * (log10(Pr_f[k]) + c[k])), 2);
            logF_f = pow(1.0 + logF_core_f, -1) * log10(chec.Fcent[l]);
            chec.F_f[l] = pow(10, logF_f);
            //forward[l] = k_inf_f[k] * (Pr_f[k] / (1 + Pr_f[k])) * chec.F_f[l];
            forward[l] = k_inf_f[k] * (Pr_f[k] / (1 + Pr_f[k])) * 1.;
        }

        if (Pr_r[k] == 0) reverse[l] = k_inf_r[k];
        else {
            logF_core_r = pow((log10(Pr_r[k]) + c[k]) / (m[k] - d * (log10(Pr_r[k]) + c[k])), 2);
            logF_r = pow(1.0 + logF_core_r, -1) * log10(chec.Fcent[l]);
            chec.F_r[l] = pow(10, logF_r);
            //reverse[l] = k_inf_r[k] * (Pr_r[k] / (1 + Pr_r[k])) * chec.F_r[l];
            reverse[l] = k_inf_r[k] * (Pr_r[k] / (1 + Pr_r[k])) * 1;
        }
        //forward[l] = k_inf_f[k];
        //reverse[l] = k_inf_r[k];
    }
    for (int i = 0; i < num_react; i++) {
        auto& prod = chec.chemkinReader->reactions()[i].getProducts();
        auto& react = chec.chemkinReader->reactions()[i].getReactants();
        sumv = 0;
        dSiR = 0;
        dHiRT = 0;
        for (const auto& iti : prod) {
            //cout << "name = " << iti.first << " = " << iti.second << '\n';
            //cout << "S = " << myget_Si(komponents[iti.first], Tcurr) * my_mol_weight(komponents[iti.first]) << "\n";
            //cout << "H = " << myget_Hi(komponents[iti.first], Tcurr) * my_mol_weight(komponents[iti.first]) << "\n";
            dSiR += myget_Si(komponents[iti.first], Tcurr) * iti.second / phyc.kR * my_mol_weight(komponents[iti.first]);
            dHiRT += myget_Hi(komponents[iti.first], Tcurr) * iti.second / phyc.kR / Tcurr * my_mol_weight(komponents[iti.first]);
            sumv += iti.second;
        }
        for (const auto& iti : react) {
            //cout << "name = " << iti.first << " = " << iti.second << '\n';
            //cout << "S = " << myget_Si(komponents[iti.first], Tcurr) * my_mol_weight(komponents[iti.first]) << "\n";
            //cout << "H = " << myget_Hi(komponents[iti.first], Tcurr) * my_mol_weight(komponents[iti.first]) << "\n";
            dSiR -= myget_Si(komponents[iti.first], Tcurr) * iti.second / phyc.kR * my_mol_weight(komponents[iti.first]);
            dHiRT -= myget_Hi(komponents[iti.first], Tcurr) * iti.second / phyc.kR / Tcurr * my_mol_weight(komponents[iti.first]);
            sumv -= iti.second;
        }
        Kpi = exp(dSiR - dHiRT);
        Kci = Kpi * pow(P / phyc.kR / Tcurr, sumv);
        /*cout << "sumv = " << sumv << "\n";
        cout << "Kpi = " << Kpi << "\n";
        cout << "Kci = " << Kci << "\n";
        cout << "forward = " << forward[i] << "\n";*/
        reverse[i] = forward[i] / Kci;
        //cout << "reverse[i] = " << reverse[i] << "\n";
        //cout << "\n\n\n";
        //cout << "rho = " << get_rho(Yi, Tcurr) << "\n";
    }

    equilib[0] = forward[0] * y[1] * y[2] - reverse[0] * y[3] * y[4];
    equilib[1] = forward[1] * y[0] * y[3] - reverse[1] * y[1] * y[4];
    equilib[2] = forward[2] * y[0] * y[4] - reverse[2] * y[1] * y[6];
    equilib[3] = forward[3] * y[6] * y[3] - reverse[3] * y[4] * y[4];
    equilib[4] = (forward[4] * y[0] - reverse[4] * y[1] * y[1]) *
        (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]);
    equilib[5] = (forward[5] * y[3] * y[3] - reverse[5] * y[2]) *
        (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]);
    equilib[6] = (forward[6] * y[3] * y[1] - reverse[6] * y[4]) *
        (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]);
    equilib[7] = (forward[7] * y[1] * y[4] - reverse[7] * y[6]) *
        (0.73 * y[0] + 3.65 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]);
    equilib[8] = forward[8] * y[1] * y[2] - reverse[8] * y[5];
    equilib[9] = forward[9] * y[0] * y[2] - reverse[9] * y[1] * y[5];
    equilib[10] = forward[10] * y[5] * y[1] - reverse[10] * y[4] * y[4];
    equilib[11] = forward[11] * y[5] * y[3] - reverse[11] * y[4] * y[2];
    equilib[12] = forward[12] * y[5] * y[4] - reverse[12] * y[6] * y[2];
    equilib[13] = forward[13] * y[5] * y[5] - reverse[13] * y[7] * y[2];
    equilib[14] = forward[14] * y[5] * y[5] - reverse[14] * y[7] * y[2];
    equilib[15] = forward[15] * y[7] - reverse[15] * y[4] * y[4];
    equilib[16] = forward[16] * y[7] - reverse[16] * y[4] * y[4];
    equilib[17] = forward[17] * y[7] * y[1] - reverse[17] * y[6] * y[4];
    equilib[18] = forward[18] * y[7] * y[1] - reverse[18] * y[5] * y[0];
    equilib[19] = forward[19] * y[7] * y[3] - reverse[19] * y[4] * y[5];
    equilib[20] = forward[20] * y[7] * y[4] - reverse[20] * y[5] * y[6];
    equilib[21] = forward[21] * y[7] * y[4] - reverse[21] * y[5] * y[6];

    yprime[0] = -equilib[1] - equilib[2] - equilib[4] - equilib[9] + equilib[18];
    yprime[1] = -equilib[0] + equilib[4] - equilib[6] - equilib[7] - equilib[8] -
        equilib[10] - equilib[17] - yprime[0];
    yprime[2] = -equilib[0] + equilib[5] - equilib[8] - equilib[9] + equilib[11] +
        equilib[12] + equilib[13] + equilib[14];
    yprime[3] = equilib[0] - equilib[1] - equilib[3] - 2. * equilib[5] - equilib[6] -
        equilib[11] - equilib[19];
    yprime[4] = equilib[0] + equilib[1] - equilib[2] + 2. * equilib[3] + equilib[6] -
        equilib[7] + 2. * equilib[10] + equilib[11] - equilib[12] + 2. * equilib[15] +
        2. * equilib[16] + equilib[17] + equilib[19] - equilib[20] - equilib[21];
    yprime[5] = equilib[8] + equilib[9] - equilib[10] - equilib[11] - equilib[12] -
        2. * equilib[13] - 2. * equilib[14] + equilib[18] + equilib[19] + equilib[20] +
        equilib[21];
    yprime[6] = equilib[2] - equilib[3] + equilib[7] + equilib[12] + equilib[17] +
        equilib[20] + equilib[21];
    yprime[7] = equilib[13] + equilib[14] - equilib[15] - equilib[16] - equilib[17] -
        equilib[18] - equilib[19] - equilib[20] - equilib[21];
    yprime[8] = 0;
    {
       /* cout << "Tcurr = " << Tcurr << "\n";
        cout << "Cp = " << Cp_all(Tcurr, Yi);
        for (int i = 0; i < num_gas_species; i++) {
            cout << "name = " << name_species[i] << " = " << y[i] << "\n";
        }
        cout << "\n\n\n\nT = " << Tcurr << "\n";
        cout << "forward " << 0 << " = " << forward[0] * y[1] * y[2] << "\n";
        cout << "reverrse " << 0 << " = " << reverse[0] * y[3] * y[4] << "\n\n";

        cout << "forward " << 1 << " = " << forward[1] * y[0] * y[3] << "\n";
        cout << "reverrse " << 1 << " = " << reverse[1] * y[1] * y[4] << "\n\n";

        cout << "forward " << 2 << " = " << forward[2] * y[0] * y[4] << "\n";
        cout << "reverrse " << 2 << " = " << reverse[2] * y[1] * y[6] << "\n\n";

        cout << "forward " << 3 << " = " << forward[3] * y[6] * y[3] << "\n";
        cout << "reverrse " << 3 << " = " << reverse[3] * y[4] * y[4] << "\n\n";

        cout << "forward " << 4 << " = " << (forward[4] * y[0]) *
            (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) << "\n";
        cout << "reverrse " << 4 << " = " << (reverse[4] * y[1] * y[1]) *
            (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) << "\n\n";

        cout << "forward " << 5 << " = " << (forward[5] * y[3] * y[3]) *
            (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) << "\n";

        cout << "reverrse " << 5 << " = " << (reverse[5] * y[2]) *
            (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) << "\n\n";

        cout << "forward " << 6 << " = " << (forward[6] * y[3] * y[1]) *
            (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) << "\n";

        cout << "reverrse " << 6 << " = " << (reverse[6] * y[4]) *
            (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) << "\n\n";

        cout << "forward " << 7 << " = " << (forward[7] * y[1] * y[4]) *
            (0.73 * y[0] + 3.65 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) << "\n";

        cout << "reverrse " << 7 << " = " << (reverse[7] * y[6]) *
            (0.73 * y[0] + 3.65 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) << "\n\n";

        cout << "forward " << 8 << " = " << forward[8] * y[1] * y[2] << "\n";
        cout << "reverrse " << 8 << " = " << reverse[8] * y[5] << "\n\n";

        cout << "forward " << 9 << " = " << forward[9] * y[0] * y[2] << "\n";
        cout << "reverrse " << 9 << " = " << reverse[9] * y[1] * y[5] << "\n\n";

        cout << "forward " << 10 << " = " << forward[10] * y[5] * y[1] << "\n";
        cout << "reverrse " << 10 << " = " << reverse[10] * y[4] * y[4] << "\n\n";

        cout << "forward " << 11 << " = " << forward[11] * y[5] * y[3] << "\n";
        cout << "reverrse " << 11 << " = " << reverse[11] * y[4] * y[2] << "\n\n";

        cout << "forward " << 12 << " = " << forward[12] * y[5] * y[4] << "\n";
        cout << "reverrse " << 12 << " = " << reverse[12] * y[6] * y[2] << "\n";

        cout << "forward " << 13 << " = " << forward[13] * y[5] * y[5] << "\n";
        cout << "reverrse " << 13 << " = " << reverse[13] * y[7] * y[2] << "\n\n";

        cout << "forward " << 14 << " = " << forward[14] * y[5] * y[5] << "\n";
        cout << "reverrse " << 14 << " = " << reverse[14] * y[7] * y[2] << "\n\n";

        cout << "forward " << 15 << " = " << forward[15] * y[7] << "\n";
        cout << "reverrse " << 15 << " = " << reverse[15] * y[4] * y[4] << "\n\n";

        cout << "forward " << 16 << " = " << forward[16] * y[7] << "\n";
        cout << "reverrse " << 16 << " = " << reverse[16] * y[4] * y[4] << "\n\n";

        cout << "forward " << 17 << " = " << forward[17] * y[7] * y[1] << "\n";
        cout << "reverrse " << 17 << " = " << reverse[17] * y[6] * y[4] << "\n\n";

        cout << "forward " << 18 << " = " << forward[18] * y[7] * y[1] << "\n";
        cout << "reverrse " << 18 << " = " << reverse[18] * y[5] * y[0] << "\n\n";


        cout << "forward " << 19 << " = " << forward[19] * y[7] * y[3] << "\n";

        cout << "reverrse " << 19 << " = " << reverse[19] * y[4] * y[5] << "\n\n";

        cout << "forward " << 20 << " = " << forward[20] * y[7] * y[4] << "\n";
        cout << "reverrse " << 20 << " = " << reverse[20] * y[5] * y[6] << "\n\n";

        cout << "forward " << 21 << " = " << forward[21] * y[7] * y[4] << "\n";
        cout << "reverrse " << 21 << " = " << reverse[21] * y[5] * y[6] << "\n\n";*/
    }

    //add_toChemVel(wk_add, M, Yi, Yinext, x, xnext, Tcurr, Tinext);
    //for (int i = 0; i < num_gas_species; i++) {
    //    cout << name_species[i] << " = " << Yi[i] << "\n";
    //}
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