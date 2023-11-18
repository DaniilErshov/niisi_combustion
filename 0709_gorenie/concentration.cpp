#include "concentration.h"

double my_mol_weight(int k) {
    //cout << "mol_weight = " << phyc.mol_weight[k] << "\n";
    return phyc.mol_weight[k] * pow(10, 3);
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
    double rho2 = P * get_W(Y) / R / T;
    double rho = P / get_gas_constant(num_gas_species, Y) / T;
    for (int i = 0; i < num_gas_species; i++) {
        //cout << "mol_weight =  " << my_mol_weight(i) << endl;
        X[i] = Y[i] * rho2 / my_mol_weight(i) * pow(10, 6);
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

void MakeYvectors(IO::ChemkinReader* chemkinReader_temp, double* x,
    double* Yiprev, double* Yi, double* Yinext,
    double* Y, double* Y_left_bound, int myNx, int i, double Tl, double M) {
    //cout << "i = " << i << "\n";
    for (int j = 0; j < num_gas_species; j++) {

        Yi[j] = Y[j + (i - 1) * num_gas_species];

        if (i == 1) { 
            FindBoundary(chemkinReader_temp, Y_left_bound, Yiprev, Yi, Tl,
                x[i - 1], x[i], M);
            //Yiprev[j] = Y_left_bound[j];
        }
        else Yiprev[j] = Y[j + (i - 2) * num_gas_species];
        //cout << "k = " << name_species[j - 1] << "\n";
        //cout << " Yiprev[j] = " << Yiprev[j - 1] << "\n";
        //cout << " Yi[j] = " << Yi[j - 1] << "\n";

        if (i == myNx - 2) Yinext[j] = Yi[j];
        else  Yinext[j] = Y[j + (i) * num_gas_species];
        //cout << " Yinext[j] = " << Yinext[j - 1] << "\n\n";
    }
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

        if (Pr_f[k] <= 0) forward[l] = k_inf_f[k];
        else {
            logF_core_f = pow((log10(Pr_f[k]) + c[k]) / (m[k] - d * (log10(Pr_f[k]) + c[k])), 2);
            logF_f = pow(1.0 + logF_core_f, -1) * log10(chec.Fcent[l]);
            chec.F_f[l] = pow(10, logF_f);
            forward[l] = k_inf_f[k] * (Pr_f[k] / (1 + Pr_f[k])) * chec.F_f[l];
        }

        if (Pr_r[k] <= 0) reverse[l] = k_inf_r[k];
        else {
            logF_core_r = pow((log10(Pr_r[k]) + c[k]) / (m[k] - d * (log10(Pr_r[k]) + c[k])), 2);
            logF_r = pow(1.0 + logF_core_r, -1) * log10(chec.Fcent[l]);
            chec.F_r[l] = pow(10, logF_r);
            reverse[l] = k_inf_r[k] * (Pr_r[k] / (1 + Pr_r[k])) * chec.F_r[l];
        }
        //forward[l] = k_inf_f[k];
        //reverse[l] = k_inf_r[k];
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
        //cout << "\n\n\n\nT = " << Tcurr << "\n";
        //cout << "forward " << 0 << " = " << forward[0] * y[1] * y[2] * pow(10, -3) << "\n";
        //cout << "reverrse " << 0 << " = " << reverse[0] * y[3] * y[4] * pow(10, -3) << "\n\n";

        //cout << "forward " << 1 << " = " << forward[1] * y[0] * y[3] * pow(10, -3) << "\n";
        //cout << "reverrse " << 1 << " = " << reverse[1] * y[1] * y[4] * pow(10, -3) << "\n\n";

        //cout << "forward " << 2 << " = " << forward[2] * y[0] * y[4] * pow(10, -3) << "\n";
        //cout << "reverrse " << 2 << " = " << reverse[2] * y[1] * y[6] * pow(10, -3) << "\n\n";

        //cout << "forward " << 3 << " = " << forward[3] * y[6] * y[3] * pow(10, -3) << "\n";
        //cout << "reverrse " << 3 << " = " << reverse[3] * y[4] * y[4] * pow(10, -3) << "\n\n";

        //cout << "forward " << 4 << " = " << (forward[4] * y[0]) *
        //    (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) * pow(10, -3) << "\n";
        //cout << "reverrse " << 4 << " = " << (reverse[4] * y[1] * y[1]) *
        //    (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) * pow(10, -3) << "\n\n";

        //cout << "forward " << 5 << " = " << (forward[5] * y[3] * y[3]) *
        //    (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) * pow(10, -3) << "\n";

        //cout << "reverrse " << 5 << " = " << (reverse[5] * y[2]) *
        //    (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) * pow(10, -3) << "\n\n";

        //cout << "forward " << 6 << " = " << (forward[6] * y[3] * y[1]) *
        //    (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) * pow(10, -3) << "\n";

        //cout << "reverrse " << 6 << " = " << (reverse[6] * y[4]) *
        //    (2.5 * y[0] + 12 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) * pow(10, -3) << "\n\n";

        //cout << "forward " << 7 << " = " << (forward[7] * y[1] * y[4]) *
        //    (0.73 * y[0] + 3.65 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) * pow(10, -3) << "\n";

        //cout << "reverrse " << 7 << " = " << (reverse[7] * y[6]) *
        //    (0.73 * y[0] + 3.65 * y[6] + y[1] + y[2] + y[3] + y[4] + y[5] + y[7] + y[8]) * pow(10, -3) << "\n\n";

        //cout << "forward " << 8 << " = " << forward[8] * y[1] * y[2] * pow(10, -3) << "\n";
        //cout << "reverrse " << 8 << " = " << reverse[8] * y[5] * pow(10, -3) << "\n\n";

        //cout << "forward " << 9 << " = " << forward[9] * y[0] * y[2] * pow(10, -3) << "\n";
        //cout << "reverrse " << 9 << " = " << reverse[9] * y[1] * y[5] * pow(10, -3) << "\n\n";

        //cout << "forward " << 10 << " = " << forward[10] * y[5] * y[1] * pow(10, -3) << "\n";
        //cout << "reverrse " << 10 << " = " << reverse[10] * y[4] * y[4] * pow(10, -3) << "\n\n";

        //cout << "forward " << 11 << " = " << forward[11] * y[5] * y[3] * pow(10, -3) << "\n";
        //cout << "reverrse " << 11 << " = " << reverse[11] * y[4] * y[2] * pow(10, -3) << "\n\n";

        //cout << "forward " << 12 << " = " << forward[12] * y[5] * y[4] * pow(10, -3) << "\n";
        //cout << "reverrse " << 12 << " = " << reverse[12] * y[6] * y[2] * pow(10, -3) << "\n";

        //cout << "forward " << 13 << " = " << forward[13] * y[5] * y[5] * pow(10, -3) << "\n";
        //cout << "reverrse " << 13 << " = " << reverse[13] * y[7] * y[2] * pow(10, -3) << "\n\n";

        //cout << "forward " << 14 << " = " << forward[14] * y[5] * y[5] * pow(10, -3) << "\n";
        //cout << "reverrse " << 14 << " = " << reverse[14] * y[7] * y[2] * pow(10, -3) << "\n\n";

        //cout << "forward " << 15 << " = " << forward[15] * y[7] * pow(10, -3) << "\n";
        //cout << "reverrse " << 15 << " = " << reverse[15] * y[4] * y[4] * pow(10, -3) << "\n\n";

        //cout << "forward " << 16 << " = " << forward[16] * y[7] * pow(10, -3) << "\n";
        //cout << "reverrse " << 16 << " = " << reverse[16] * y[4] * y[4] * pow(10, -3) << "\n\n";

        //cout << "forward " << 17 << " = " << forward[17] * y[7] * y[1] * pow(10, -3) << "\n";
        //cout << "reverrse " << 17 << " = " << reverse[17] * y[6] * y[4] * pow(10, -3) << "\n\n";

        //cout << "forward " << 18 << " = " << forward[18] * y[7] * y[1] * pow(10, -3) << "\n";
        //cout << "reverrse " << 18 << " = " << reverse[18] * y[5] * y[0] * pow(10, -3) << "\n\n";


        //cout << "forward " << 19 << " = " << forward[19] * y[7] * y[3] * pow(10, -3) << "\n";

        //cout << "reverrse " << 19 << " = " << reverse[19] * y[4] * y[5] * pow(10, -3) << "\n\n";

        //cout << "forward " << 20 << " = " << forward[20] * y[7] * y[4] * pow(10, -3) << "\n";
        //cout << "reverrse " << 20 << " = " << reverse[20] * y[5] * y[6] * pow(10, -3) << "\n\n";

        //cout << "forward " << 21 << " = " << forward[21] * y[7] * y[4] * pow(10, -3) << "\n";
        //cout << "reverrse " << 21 << " = " << reverse[21] * y[5] * y[6] * pow(10, -3) << "\n\n";
    }

    //add_toChemVel(wk_add, M, Yi, Yinext, x, xnext, Tcurr, Tinext);

    for (int i = 0; i < num_gas_species; i++) {
        yprime[i] *= pow(10, -3);
        //cout << "velo = " << name_species[i] << " = " << yprime[i] << "\n";
        //yprime[i] += wk_add[i];
    }
}

double YkVk(IO::ChemkinReader* chemkinReader, int k, double T, double* Y, double* gradX, double* Xi) {
    double sum = 0.;
    double W = get_W(Y);
    double rho = get_rho(Y, T);
    double YkVk = 0;
    double Dkm = 0;
    for (int j = 0; j < num_gas_species; j++) {
        if (j != k) {
            sum += Xi[j]
                / Dij_func(chemkinReader, k, j, T, Y);
        }
    }
    Dkm = (1. - Y[k]) / sum;
    //cout << "Dkm for " << name_species[k] << " = " << Dkm << "\n";
    return  -my_mol_weight(k) / W * Dkm * gradX[k];
}

double Dij_func(IO::ChemkinReader* chemkinReader, int i, int j, double T, double* Y)
{
    //double D_res = 3. / 16. * pow(2 * 3.14156 * pow(kB, 3), 0.5) / 3.14156;
    //cout << "konst = " << D_res << "\n";
    double D_res = 1.858 * pow(10, -3);
    double mij = my_mol_weight(i) * my_mol_weight(j) / (my_mol_weight(i) + my_mol_weight(j));
    double rho = get_rho(Y, T);
    double sigma_ij = 0.5 * (chemkinReader->species()[i].transport().getCollisionDiameter() +
        chemkinReader->species()[j].transport().getCollisionDiameter()) / Angstroem__;
    double eps_ij = pow(chemkinReader->species()[i].transport().getPotentialWellDepth() *
        chemkinReader->species()[j].transport().getPotentialWellDepth(), 0.5);
    double T_ij = kB * T / eps_ij;
    double integ_ij = 1.074 * pow(T_ij, -0.1604);
    /*cout << "rho = " << rho << "\n";
    cout << "sigmaij = " << sigma_ij << "\n";
    cout << "T_ij = " << T_ij << "\n";*/
    D_res *= pow(pow(T, 3) / mij, 0.5) / 1. / pow(sigma_ij, 2.) / integ_ij;
    //1.033227
    return D_res;
}

double Dk_func(IO::ChemkinReader* chemkinReader, int i, double T, double* Y, double* X, int N) {
    double sum = 0;
    for (int j = 0; j < num_gas_species; j++) {
        if (i != j) sum += X[j] / Dij_func(chemkinReader, i, j, T, Y);
        cout << "Dij = " << Dij_func(chemkinReader, i, j, T, Y) << "\n";
    }
    return (1 - X[i]) / sum;
}