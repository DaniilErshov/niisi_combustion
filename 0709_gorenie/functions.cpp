#include "functions.h"

static int check_retval(void* retvalvalue, const char* funcname, int opt)
{
    int* errretval;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && retvalvalue == NULL) {
        fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
        return(1);
    }

    /* Check if retval < 0 */
    else if (opt == 1) {
        errretval = (int*)retvalvalue;
        if (*errretval < 0) {
            fprintf(stderr,
                "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
                funcname, *errretval);
            return(1);
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && retvalvalue == NULL) {
        fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
        return(1);
    }

    return(0);
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
        W += Y[i] / (phyc.mol_weight[i] * k_mol);
    }
    return 1. / W;
}

double get_rho(double* Y, double T) {
    return P * get_W(Y) / R / T;
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
        X[i] = Y[i] * W / (phyc.mol_weight[i] * k_mol);
    }
}

void Get_molar_cons(double* X, double* Y, double T)
{
    double W = 0;
    double rho = P * get_W(Y) / R / T;
    //cout << "W = " << W << "\n";
    for (int i = 0; i < num_gas_species; i++) {
        //cout << "mol_weight =  " << phyc.mol_weight[i] << endl;
        X[i] = Y[i] * rho / (phyc.mol_weight[i] * k_mol);
    }
}

double Cp_all(double T, double* Y)
{
    double cp_tmp = 0.;
    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        cp_tmp += get_Cpi(k_spec, T) * Y[k_spec];
    };
    return cp_tmp;
}

void get_grad(double* gradX, double* Xi, double* Xinext, double x, double xnext) {
    for (int i = 0; i < num_gas_species; i++) {
        gradX[i] = (Xinext[i] - Xi[i]) / (xnext - x);
    }
}

double Lambda_All(IO::ChemkinReader* chemkinReader, double* Y, double T)
{
    vector<double> lambda_spec(num_gas_species);
    double* X = new double[num_gas_species];;
    double Ti;
    double integr_i;
    double sigma_i;
    double lamb1 = 0.;
    double lamb2 = 0.;
    Get_mole_fr(X, Y);
    for (int i = 0; i < num_gas_species; i++) {

        //std::cout << std::endl << std::endl << std::endl << chemkinReader->species()[i].name() << std::endl;
        sigma_i = chemkinReader->species()[i].transport().getCollisionDiameter() / Angstroem__;
        //cout << "sigma_i = " << sigma_i << endl;
        Ti = kB * T / (chemkinReader->species()[i].transport().getPotentialWellDepth());
        //cout << "Ti = " << Ti << endl;
        integr_i = 1.157 / pow(Ti, 0.1472);
        //cout << "integr_i = " << integr_i << endl;
        lambda_spec[i] = 8330. * pow(T / M[i], 0.5) / pow(sigma_i, 2.) / integr_i * pow(10, -7);
        //cout << "lambda_spec[i] = " << lambda_spec[i] << endl;
    }
    for (int i = 0; i < num_gas_species; i++) {
        lamb1 += X[i] * lambda_spec[i];
        lamb2 += X[i] / lambda_spec[i];
    }
    delete[] X;
    return 0.5 * (lamb1 + 1. / lamb2);
}

double Dij_func(IO::ChemkinReader* chemkinReader, int i, int j, double T, double* Y)
{
    double D_res = 3. / 16. * pow(2 * 3.14156 * pow(kB, 3), 0.5) / 3.14156;
    //cout << "konst = " << D_res << "\n";
    double mij = M[i] * M[j] / (M[i] + M[j]);
    double rho = P * get_W(Y) / R / T;
    double sigma_ij = 0.5 * (chemkinReader->species()[i].transport().getCollisionDiameter() +
        chemkinReader->species()[j].transport().getCollisionDiameter()) / Angstroem__;
    double eps_ij = pow(chemkinReader->species()[i].transport().getPotentialWellDepth() *
        chemkinReader->species()[j].transport().getPotentialWellDepth(), 0.5);
    double T_ij = kB * T / eps_ij;
    double integ_ij = 1.074 * pow(T_ij, -0.1604);
    //cout << "rho = " << rho << "\n";
    //cout << "sigmaij = " << sigma_ij << "\n";
    //cout << "T_ij = " << T_ij << "\n";
    D_res *= pow(pow(T, 3) / mij, 0.5) / pow(10, 5) / pow(sigma_ij, 2.) / integ_ij;
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

double rhoYkWk(IO::ChemkinReader* chemkinReader, int k, double T, double* Y, double* gradX, double gradT) {

    //std::cout << std::endl << std::endl << std::endl << chemkinReader->species()[i].name() << std::endl;
    double sigma_i = chemkinReader->species()[k].transport().getCollisionDiameter() / Angstroem__;
    //cout << "sigma_i = " << sigma_i << endl;
    double Ti = kB * T / (chemkinReader->species()[k].transport().getPotentialWellDepth());
    //cout << "Ti = " << Ti << endl;
    double integr_i = 1.157 / pow(Ti, 0.1472);
    //cout << "integr_i = " << integr_i << endl;
    double lambda_spec = 8330. * pow(T / M[k], 0.5) / pow(sigma_i, 2.) / integr_i * pow(10, -7);
    double rho = P * get_W(Y) / R / T;
    double Dk = lambda_spec / rho / (get_Cpi(k, T) * pow(10, 3));
    //cout << "lambda_spec = " << lambda_spec << "\n";
    cout << "from Wk Dk = " << Dk << "\n";
    return Dk * (gradT / T);
}

double rhoYkVk(IO::ChemkinReader* chemkinReader, int k, double T, double* Y, double* gradX, double gradT) {
    double sum = 0.;
    double W = get_W(Y);
    double rho = get_rho(Y, T);
    double YkVk = 0;
    for (int j = 0; j < num_gas_species; j++) {
        if (j != k) {
            sum += (phyc.mol_weight[j] * k_mol)
                * Dij_func(chemkinReader, k, j, T, Y)
                * gradX[j];
            //cout << "J = " << j << "\n";
            //cout << "Dij = " << Dij_func(chemkinReader, k, j, T, Y) << "\n";
            //cout << "GradX = " << gradX[j] << "\n";
            //cout << "molw = " << (phyc.mol_weight[j] * k_mol) << "\n\n\n\n";
        }
    }
    cout << "rhoYkV0k = " << sum * (phyc.mol_weight[k] * k_mol) / pow(W, 2) << "\n";
    cout << "rhoYkW0k = " << rhoYkWk(chemkinReader, k, T, Y, gradX, gradT) << "\n";
    return sum * (phyc.mol_weight[k] * k_mol) / pow(W, 2) + rhoYkWk(chemkinReader, k, T, Y, gradX, gradT);
}
//here send molar_cons
void chem_vel(double Tcurr, N_Vector y, N_Vector ydot) {
    double k_0_f[3], k_inf_f[3], k_0_r[3], k_inf_r[3];
    double c[3], m[3], d = 0.14;
    double Pr_f[3], Pr_r[3];
    int k = 0, l = 0;
    double logF_f, logF_core_f, logF_r, logF_core_r;
    double sum1, sum2;
    double* forward = new double[num_react];
    double* reverse = new double[num_react];
    double* equilib = new double[num_react];
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

    Pr_f[0] = (k_0_f[0] * (1.3 * Ith(y, 1) + 10 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9))) / k_inf_f[0];
    Pr_r[0] = (k_0_r[0] * (1.3 * Ith(y, 1) + 10 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9))) / k_inf_r[0];
    Pr_f[1] = (k_0_f[1] * (3.7 * Ith(y, 1) + 1.5 * Ith(y, 9) + 1.2 * Ith(y, 3) + 7.7 * Ith(y, 8) + Ith(y, 2) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 7))) / k_inf_f[1];
    Pr_r[1] = (k_0_r[1] * (3.7 * Ith(y, 1) + 1.5 * Ith(y, 9) + 1.2 * Ith(y, 3) + 7.7 * Ith(y, 8) + Ith(y, 2) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 7))) / k_inf_r[1];
    Pr_f[2] = (k_0_f[2] * Ith(y, 7)) / k_inf_f[2];
    Pr_r[2] = (k_0_r[2] * Ith(y, 7)) / k_inf_r[2];

    for (k = 0; k < 3; k++) {
        if (k == 0) l = 8;
        if (k == 1) l = 15;
        if (k == 2) l = 16;

        if (Pr_f[k] == 0) forward[l] = k_inf_f[k];
        else {
            logF_core_f = pow((log10(Pr_f[k]) + c[k]) / (m[k] - d * (log10(Pr_f[k]) + c[k])), 2);
            logF_f = pow(1.0 + logF_core_f, -1) * log10(chec.Fcent[l]);
            chec.F_f[l] = pow(10, logF_f);
            forward[l] = k_inf_f[k] * (Pr_f[k] / (1 + Pr_f[k])) * chec.F_f[l];
        }

        if (Pr_r[k] == 0) reverse[l] = k_inf_r[k];
        else {
            logF_core_r = pow((log10(Pr_r[k]) + c[k]) / (m[k] - d * (log10(Pr_r[k]) + c[k])), 2);
            logF_r = pow(1.0 + logF_core_r, -1) * log10(chec.Fcent[l]);
            chec.F_r[l] = pow(10, logF_r);
            reverse[l] = k_inf_r[k] * (Pr_r[k] / (1 + Pr_r[k])) * chec.F_r[l];
        }
    }

    equilib[0] = forward[0] * Ith(y, 2) * Ith(y, 3) - reverse[0] * Ith(y, 4) * Ith(y, 5);
    equilib[1] = forward[1] * Ith(y, 1) * Ith(y, 4) - reverse[1] * Ith(y, 2) * Ith(y, 5);
    equilib[2] = forward[2] * Ith(y, 1) * Ith(y, 5) - reverse[2] * Ith(y, 2) * Ith(y, 7);
    equilib[3] = forward[3] * Ith(y, 7) * Ith(y, 4) - reverse[3] * Ith(y, 5) * Ith(y, 5);
    equilib[4] = (forward[4] * Ith(y, 1) - reverse[4] * Ith(y, 2) * Ith(y, 2)) *
        (2.5 * Ith(y, 1) + 12 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9));
    equilib[5] = (forward[5] * Ith(y, 4) * Ith(y, 4) - reverse[5] * Ith(y, 3)) *
        (2.5 * Ith(y, 1) + 12 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9));
    equilib[6] = (forward[6] * Ith(y, 4) * Ith(y, 2) - reverse[6] * Ith(y, 5)) *
        (2.5 * Ith(y, 1) + 12 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9));
    equilib[7] = (forward[7] * Ith(y, 2) * Ith(y, 5) - reverse[7] * Ith(y, 7)) *
        (0.73 * Ith(y, 1) + 3.65 * Ith(y, 7) + Ith(y, 2) + Ith(y, 3) + Ith(y, 4) + Ith(y, 5) + Ith(y, 6) + Ith(y, 8) + Ith(y, 9));
    equilib[8] = forward[8] * Ith(y, 2) * Ith(y, 3) - reverse[8] * Ith(y, 6);
    equilib[9] = forward[9] * Ith(y, 1) * Ith(y, 3) - reverse[9] * Ith(y, 2) * Ith(y, 6);
    equilib[10] = forward[10] * Ith(y, 6) * Ith(y, 2) - reverse[10] * Ith(y, 5) * Ith(y, 5);
    equilib[11] = forward[11] * Ith(y, 6) * Ith(y, 4) - reverse[11] * Ith(y, 5) * Ith(y, 3);
    equilib[12] = forward[12] * Ith(y, 6) * Ith(y, 5) - reverse[12] * Ith(y, 7) * Ith(y, 3);
    equilib[13] = forward[13] * Ith(y, 6) * Ith(y, 6) - reverse[13] * Ith(y, 8) * Ith(y, 3);
    equilib[14] = forward[14] * Ith(y, 6) * Ith(y, 6) - reverse[14] * Ith(y, 8) * Ith(y, 3);
    equilib[15] = forward[15] * Ith(y, 8) - reverse[15] * Ith(y, 5) * Ith(y, 5);
    equilib[16] = forward[16] * Ith(y, 8) - reverse[16] * Ith(y, 5) * Ith(y, 5);
    equilib[17] = forward[17] * Ith(y, 8) * Ith(y, 2) - reverse[17] * Ith(y, 7) * Ith(y, 5);
    equilib[18] = forward[18] * Ith(y, 8) * Ith(y, 2) - reverse[18] * Ith(y, 6) * Ith(y, 1);
    equilib[19] = forward[19] * Ith(y, 8) * Ith(y, 4) - reverse[19] * Ith(y, 5) * Ith(y, 6);
    equilib[20] = forward[20] * Ith(y, 8) * Ith(y, 5) - reverse[20] * Ith(y, 6) * Ith(y, 7);
    equilib[21] = forward[21] * Ith(y, 8) * Ith(y, 5) - reverse[21] * Ith(y, 6) * Ith(y, 7);


    Ith(ydot, 1) = -equilib[1] - equilib[2] - equilib[4] - equilib[9] + equilib[18];
    Ith(ydot, 2) = -equilib[0] + equilib[4] - equilib[6] - equilib[7] - equilib[8] -
        equilib[10] - equilib[17] - Ith(ydot, 1);
    Ith(ydot, 3) = -equilib[0] + equilib[5] - equilib[8] - equilib[9] + equilib[11] +
        equilib[12] + equilib[13] + equilib[14];
    Ith(ydot, 4) = equilib[0] - equilib[1] - equilib[3] - 2. * equilib[5] - equilib[6] -
        equilib[11] - equilib[19];
    Ith(ydot, 5) = equilib[0] + equilib[1] - equilib[2] + 2. * equilib[3] + equilib[6] -
        equilib[7] + 2. * equilib[10] + equilib[11] - equilib[12] + 2. * equilib[15] +
        2. * equilib[16] + equilib[17] + equilib[19] - equilib[20] - equilib[21];
    Ith(ydot, 6) = equilib[8] + equilib[9] - equilib[10] - equilib[11] - equilib[12] -
        2. * equilib[13] - 2. * equilib[14] + equilib[18] + equilib[19] + equilib[20] +
        equilib[21];
    Ith(ydot, 7) = equilib[2] - equilib[3] + equilib[7] + equilib[12] + equilib[17] +
        equilib[20] + equilib[21];
    Ith(ydot, 8) = equilib[13] + equilib[14] - equilib[15] - equilib[16] - equilib[17] -
        equilib[18] - equilib[19] - equilib[20] - equilib[21];
    Ith(ydot, 9) = 0;
}

double F_right(IO::ChemkinReader* chemkinReader, double* Yi, double* Yinext,
    double* T, double* Xiprev, double* Xi, double* Xinext, double* gradX, double* Y_tmp,
    double M, realtype* x_cells, const int i, N_Vector ydot)
{
    double h_left = x_cells[i] - x_cells[i - 1];
    double h = x_cells[i + 1] - x_cells[i];
    double Cp = Cp_all(T[i], Yi) * pow(10, 3);
    double rho = get_rho(Yi, T[i]);
    double W = get_W(Yi);
    Get_mole_fr(Xi, Yi);
    Get_mole_fr(Xinext, Yinext);
    get_grad(gradX, Xi, Xinext, x_cells[i], x_cells[i + 1]);
    double dTdx = (h_left / h / (h + h_left) * T[i + 1] + (h - h_left) / h / h_left * T[i] - h / h_left / (h + h_left) * T[i - 1]);
    double slag_diff = 0;
    double slag_chem = 0.;
    for (int k = 0; k < num_gas_species; k++) {
        cout << " Cpi = " << get_Cpi(k, T[i]) * pow(10, 3) << "\n";
        cout << " rhoYkVk = " << rhoYkVk(chemkinReader, k, T[i], Yi, gradX, dTdx) << "\n";
        cout << "dTdx = " << dTdx << "\n";
        slag_diff += rhoYkVk(chemkinReader, k, T[i], Yi, gradX, dTdx) * get_Cpi(k, T[i]) * pow(10, 3) * dTdx;
        slag_chem += Ith(ydot, k + 1) * get_Hi(k, T[i]) * phyc.mol_weight[k] * pow(10, 6);
        //cout << "Vk_func000 = " << Vk_func(chemkinReader, k, Xi[k], T[i], Yi,
        //    Xi, gradX, dTdx) << "\n";
    }
    cout << " slag_diff = " << slag_diff << "\n";
    //cout << " slag_chem = " << slag_chem << "\n";
    //cout << "dTdx = " << Cp * M * dTdx << "\n";

    return -(2. / (h + h_left)) *
        (Lambda_All(chemkinReader, Yi, (T[i + 1] + T[i]) / 2.) * (T[i + 1] - T[i]) / h
            - Lambda_All(chemkinReader, Yi, (T[i] + T[i - 1]) / 2.) * (T[i] - T[i - 1]) / h_left)
        + Cp * M * dTdx + slag_diff + slag_chem;

}

double F_rightY(IO::ChemkinReader* chemkinReader, double* Yiprev, double* Yi, double* Yinext,
    double* T, double* Xiprev, double* Xi, double* Xinext, double* gradX, double* Y_tmp,
    double M, realtype* x_cells, const int i, const int k_spec, N_Vector ydot)
{

    double h_left = x_cells[i] - x_cells[i - 1];
    double h = x_cells[i + 1] - x_cells[i];
    double x_l = (x_cells[i] + x_cells[i - 1]) / 2.;
    double x_r = (x_cells[i + 1] + x_cells[i]) / 2.;
    double rho = 0.;
    double rhoYkVk_r = 1.;
    double rhoYkVk_l = 1.;
    double slag_diff = 0.;
    double slag_chem = 0.;
    double dTdx = (h_left / h / (h + h_left) * T[i + 1] + (h - h_left) / h / h_left * T[i] - h / h_left / (h + h_left) * T[i - 1]);
    Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);
    make_averageY(Y_tmp, Yi, Yinext);
    get_grad(gradX, Xi, Xinext, x_cells[i], x_cells[i + 1]);
    rhoYkVk_r = rhoYkVk(chemkinReader, k_spec, T[i] / 2. + T[i + 1] / 2., Y_tmp, gradX, dTdx);


    make_averageY(Y_tmp, Yi, Yiprev);
    get_grad(gradX, Xiprev, Xi, x_cells[i - 1], x_cells[i]);
    rhoYkVk_l = rhoYkVk(chemkinReader, k_spec, T[i] / 2. + T[i - 1] / 2., Y_tmp, gradX, dTdx);


    slag_diff = (rhoYkVk_r - rhoYkVk_l) / (x_r - x_l);

    slag_chem = -phyc.mol_weight[k_spec] * k_mol * Ith(ydot, k_spec + 1);
    //cout << "\n\n\nY slag_diff" << slag_diff << "\n";
    //cout << "Y slag_chem" << slag_chem << "\n";
    //cout << "Y slag M = " << M * (h_left / h / (h + h_left) * Yinext[k_spec] + (h - h_left) / h / h_left * Yi[k_spec] - h / h_left / (h + h_left) * Yiprev[k_spec]) << "\n";
    return M * (h_left / h / (h + h_left) * Yinext[k_spec] + (h - h_left) / h / h_left * Yi[k_spec] - h / h_left / (h + h_left) * Yiprev[k_spec])
        + slag_diff + slag_chem;
}

void MakeYvectors(double* Yiprev, double* Yi, double* Yinext, double* Y, double* Y_left_bound, int myNx, int i) {
    //cout << "i = " << i << "\n";
    for (int j = 1; j <= num_gas_species; j++) {

        Yi[j - 1] = Y[j + (i - 1) * num_gas_species];

        if (i == 1) Yiprev[j - 1] = Y_left_bound[j - 1];
        else Yiprev[j - 1] = Y[j + (i - 2) * num_gas_species];
        //cout << " Yiprev[j] = " << Yiprev[j - 1] << "\n";
        //cout << " Yi[j] = " << Yi[j - 1] << "\n";

        if (i == myNx - 2) Yinext[j - 1] = Yi[j - 1];
        else  Yinext[j - 1] = Y[j + (i)*num_gas_species];
        //cout << " Yinext[j] = " << Yinext[j - 1] << "\n\n";
    }
}

static int func_Y(N_Vector u, N_Vector f, void* user_data)
{
    realtype* Y, * fdata;
    realtype x1, l1, L1, x2, l2, L2;
    realtype* x_cells, * T_vect, * Yi, * Yiprev, * Yinext;
    UserData data;
    double h_left, h;
    double tmp;
    double M;
    data = (UserData)user_data;
    x_cells = data->x;
    T_vect = data->T;
    Yi = data->Yi;
    Yiprev = data->Yiprev;
    Yinext = data->Yinext;
    int myNx = data->Nx;
    int myNeq = data->NEQ;
    int myNm = data->N_m;

    Y = N_VGetArrayPointer(u);
    fdata = N_VGetArrayPointer(f);
    double T_cer = data->T_center;

    //cout << "MyNx = " << myNx << "\n";
    //cout << "MyNm = " << myNm << "\n";
    data->M = Y[myNm];
    int j = 1;
    cout << "M = " << data->M << "\n";
    //Get_mole_fr(data->X_mol, Yi);

    int fl = 0;
    int Ncentr = data->N_centr;
    MakeYvectors(Yiprev, Yi, Yinext, Y, data->Y_left_bound, myNx, Ncentr);
    Get_molar_cons(data->Xi, data->Yi, T_vect[Ncentr]);
    for (int k = 0; k < num_gas_species; k++) Ith(data->y, k + 1) = data->Xi[k];
    chem_vel(T_vect[Ncentr], data->y, data->ydot);


    fdata[0] = F_right(data->chemkinReader, Yi, Yinext,
        data->T, data->Xiprev, data->Xi, data->Xinext, data->gradX, data->Y_tmp,
        data->M, x_cells, Ncentr, data->ydot);

    //cout << "fdata[0]" << fdata[0] << "\n";
    double sum = fdata[0];
    for (int i = 1; i < myNx - 1; i++) {

        MakeYvectors(Yiprev, Yi, Yinext, Y, data->Y_left_bound, myNx, i);
        Get_molar_cons(data->Xi, data->Yi, T_vect[i]);
        for (int k = 0; k < num_gas_species; k++) Ith(data->y, k + 1) = data->Xi[k];
        chem_vel(T_vect[i], data->y, data->ydot);


        for (int k_spec = 1; k_spec <= num_gas_species; k_spec++) {
            fdata[k_spec + (i - 1) * num_gas_species] = F_rightY(data->chemkinReader, Yiprev, Yi, Yinext,
                T_vect, data->Xiprev, data->Xi, data->Xinext, data->gradX, data->Y_tmp,
                data->M, x_cells, i, k_spec - 1, data->ydot);
            sum += fdata[k_spec + (i - 1) * num_gas_species];
            //if (fdata[k_spec + (i - 1) * num_gas_species] > pow(10, 10)) cout << "AAAAAAAAAAAAAAAAAAAAAAAaaaaaa\n";
            //cout << "fdata [" << k_spec + (i - 1) * num_gas_species << "] = " << fdata[k_spec + (i - 1) * num_gas_species] << "\n";
        }
        //cout << "j - 1 = " << j - 1 << "\n";
        //cout << "fdata " << i << " = " << fdata[i] << "\n";
    }
    cout << "sum = " << sum << "\n";
    return(0);
}

static int funcFinalState(N_Vector u, N_Vector f, void* user_data)
{
    realtype* Y, * fdata;
    realtype x1, l1, L1, x2, l2, L2;
    realtype* x_cells, * T_vect, * Yi, * Yiprev, * Yinext;
    UserData data;
    double h_left, h;
    double tmp;
    double Temp;
    data = (UserData)user_data;
    x_cells = data->x;
    T_vect = data->T;
    Yi = data->Yi;
    Yiprev = data->Yiprev;
    Yinext = data->Yinext;
    int myNx = data->Nx;
    int myNeq = data->NEQ;
    int myNm = data->N_m;

    Y = N_VGetArrayPointer(u);
    fdata = N_VGetArrayPointer(f);

    for (int j = 0; j < num_gas_species - 1; j++) {
        Yi[j] = Y[j];
    }
    Temp = Y[num_gas_species - 1];

    Get_molar_cons(data->Xi, Yi, Temp);
    for (int k = 0; k < num_gas_species - 1; k++) Ith(data->y, k + 1) = data->Xi[k];
    chem_vel(Temp, data->y, data->ydot);

    for (int i = 0; i < num_gas_species - 1; i++) {
        fdata[i] = Ith(data->ydot, i + 1);
    }
    double sum = 0;
    for (int i = 0; i < num_gas_species - 1; i++) {
        sum += get_Hi(i, Temp) * Ith(data->ydot, i + 1) * phyc.mol_weight[i];
    }

    fdata[num_gas_species - 1] = sum;
    cout << "Temp = " << Temp << "\n";
    return(0);
}

int Integrate_Y(IO::ChemkinReader* chemkinReader_temp, int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb)
{
    SUNContext sunctx;
    UserData data;
    realtype fnormtol, scsteptol;
    N_Vector res_vect, s, c;
    int glstr, mset, retval;
    void* kmem;
    SUNMatrix J;
    SUNLinearSolver LS;
    //int NEQ_T = N_x - 2;
    int NEQ_Y = num_gas_species * (N_x - 2);
    int NEQ = 1 + NEQ_Y;

    res_vect = NULL;
    s = c = NULL;
    kmem = NULL;
    J = NULL;
    LS = NULL;
    data = NULL;

    /* Create the SUNDIALS context that all SUNDIALS objects require */
    retval = SUNContext_Create(NULL, &sunctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

    /* User data */

    data = (UserData)malloc(sizeof * data);

    /* Create serial vectors of length NEQ */

    res_vect = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)res_vect, "N_VNew_Serial", 0)) return(1);

    s = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)s, "N_VNew_Serial", 0)) return(1);

    c = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)c, "N_VNew_Serial", 0)) return(1);

    N_VConst(ONE, s); /* no scaling */

    //data  fill
    data->Nx = N_x;
    data->x = new realtype[N_x];

    data->Yi = new realtype[num_gas_species];
    data->Yiprev = new realtype[num_gas_species];
    data->Yinext = new realtype[num_gas_species];

    data->Xiprev = new realtype[num_gas_species];
    data->Xi = new realtype[num_gas_species];
    data->Xinext = new realtype[num_gas_species];
    data->gradX = new realtype[num_gas_species];
    data->Y_tmp = new realtype[num_gas_species];

    data->Y_left_bound = new realtype[num_gas_species];
    data->T = new realtype[N_x];
    data->NEQ = NEQ;
    data->Tl = T_vect[0];
    data->T_center = T_vect[N_center];
    cout << "T_center = " << data->T_center << "\n";
    data->N_centr = N_center;
    data->chemkinReader = chemkinReader_temp;
    data->y = N_VNew_Serial(num_gas_species, sunctx);
    data->ydot = N_VNew_Serial(num_gas_species, sunctx);
    data->N_m = 1 - 1;

    int j = 0;
    for (int i = 0; i < N_x; i++) {
        data->x[i] = x_vect[i];
        data->T[i] = T_vect[i];
    }
    for (int i = 0; i < num_gas_species; i++) {
        Ith(data->y, i + 1) = 5;
        Ith(data->ydot, i + 1) = 5;
        data->Y_left_bound[i] = Y_leftb[i];
    }

    Ith(c, 1) = ONE;
    Ith(res_vect, 1) = M;
    int k = 2;
    for (int i = 1; i < NEQ; i++) {
        Ith(c, k) = 0.0;   /* no constraint on x1 */
        Ith(res_vect, k) = Y_vect[num_gas_species + i - 1];
        cout << "k = " << k << "   = " << Ith(res_vect, k) << endl;
        if ((i) % num_gas_species == 0) cout << endl;
        k++;
    }

    fnormtol = FTOL; scsteptol = STOL;


    kmem = KINCreate(sunctx);
    if (check_retval((void*)kmem, "KINCreate", 0)) return(1);

    retval = KINSetUserData(kmem, data);
    if (check_retval(&retval, "KINSetUserData", 1)) return(1);
    retval = KINSetConstraints(kmem, c);
    if (check_retval(&retval, "KINSetConstraints", 1)) return(1);
    retval = KINSetFuncNormTol(kmem, fnormtol);
    if (check_retval(&retval, "KINSetFuncNormTol", 1)) return(1);
    retval = KINSetScaledStepTol(kmem, scsteptol);
    if (check_retval(&retval, "KINSetScaledStepTol", 1)) return(1);

    retval = KINInit(kmem, func_Y, res_vect);
    if (check_retval(&retval, "KINInit", 1)) return(1);

    /* Create dense SUNMatrix */
    J = SUNDenseMatrix(NEQ, NEQ, sunctx);
    if (check_retval((void*)J, "SUNDenseMatrix", 0)) return(1);

    /* Create dense SUNLinearSolver object */
    LS = SUNLinSol_Dense(res_vect, J, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return(1);

    /* Attach the matrix and linear solver to KINSOL */
    retval = KINSetLinearSolver(kmem, LS, J);
    if (check_retval(&retval, "KINSetLinearSolver", 1)) return(1);

    glstr = 0;
    mset = 500;

    data->mykmem = kmem;
    retval = KINSol(kmem, res_vect, glstr, s, s);
    if (check_retval(&retval, "KINSol", 1)) return(1);
    //PrintFinalStats(kmem);
    cout << "retval = " << retval << "\n";

    int myNx = data->Nx;
    int myNeq = data->NEQ;
    int myNm = data->N_m;

    T_vect[0] = data->Tl;
    //cout << "MyNx = " << myNx << "\n";
    //cout << "MyNm = " << myNm << "\n";
    M = Ith(res_vect, myNm + 1);
    k = 2;
    for (int i = 1; i < NEQ; i++) {
        Y_vect[num_gas_species + i - 1] = Ith(res_vect, k);
        k++;
    }
    for (int i = NEQ; i < N_x; i++) {
        Y_vect[i] = Y_vect[i - num_gas_species];
    }

    /* Free memory */
    printf("\nFinal statsistics:\n");
    retval = KINPrintAllStats(kmem, stdout, SUN_OUTPUTFORMAT_TABLE);
    N_VDestroy(res_vect);
    N_VDestroy(s);
    N_VDestroy(c);
    KINFree(&kmem);
    SUNLinSolFree(LS);
    SUNMatDestroy(J);
    free(data);
    SUNContext_Free(&sunctx);
    return 0;
}

int FindFinalState(IO::ChemkinReader* chemkinReader_temp, double Tinitial, double* Y_vect)
{
    SUNContext sunctx;
    UserData data;
    realtype fnormtol, scsteptol;
    N_Vector res_vect, s, c;
    int glstr, mset, retval;
    void* kmem;
    SUNMatrix J;
    SUNLinearSolver LS;
    //int NEQ_T = N_x - 2;
    int NEQ_Y = num_gas_species - 1;
    int NEQ = 1 + NEQ_Y;

    data->Xi = new realtype[num_gas_species];
    data->Yi = new realtype[num_gas_species];
    data->y = N_VNew_Serial(num_gas_species, sunctx);
    data->ydot = N_VNew_Serial(num_gas_species, sunctx);

    res_vect = NULL;
    s = c = NULL;
    kmem = NULL;
    J = NULL;
    LS = NULL;
    data = NULL;

    /* Create the SUNDIALS context that all SUNDIALS objects require */
    retval = SUNContext_Create(NULL, &sunctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

    /* User data */

    data = (UserData)malloc(sizeof * data);

    /* Create serial vectors of length NEQ */

    res_vect = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)res_vect, "N_VNew_Serial", 0)) return(1);

    s = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)s, "N_VNew_Serial", 0)) return(1);

    c = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)c, "N_VNew_Serial", 0)) return(1);

    N_VConst(ONE, s); /* no scaling */

    for (int i = 1; i < num_gas_species; i++) {
        Ith(c, i) = 0.0;   /* no constraint on x1 */
        Ith(res_vect, i) = Y_vect[i - 1];
    }
    Ith(c, num_gas_species) = ONE;
    Ith(res_vect, num_gas_species) = Tinitial;
    fnormtol = FTOL; scsteptol = STOL;


    kmem = KINCreate(sunctx);
    if (check_retval((void*)kmem, "KINCreate", 0)) return(1);

    retval = KINSetUserData(kmem, data);
    if (check_retval(&retval, "KINSetUserData", 1)) return(1);
    retval = KINSetConstraints(kmem, c);
    if (check_retval(&retval, "KINSetConstraints", 1)) return(1);
    retval = KINSetFuncNormTol(kmem, fnormtol);
    if (check_retval(&retval, "KINSetFuncNormTol", 1)) return(1);
    retval = KINSetScaledStepTol(kmem, scsteptol);
    if (check_retval(&retval, "KINSetScaledStepTol", 1)) return(1);

    retval = KINInit(kmem, funcFinalState, res_vect);
    if (check_retval(&retval, "KINInit", 1)) return(1);

    /* Create dense SUNMatrix */
    J = SUNDenseMatrix(NEQ, NEQ, sunctx);
    if (check_retval((void*)J, "SUNDenseMatrix", 0)) return(1);

    /* Create dense SUNLinearSolver object */
    LS = SUNLinSol_Dense(res_vect, J, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return(1);

    /* Attach the matrix and linear solver to KINSOL */
    retval = KINSetLinearSolver(kmem, LS, J);
    if (check_retval(&retval, "KINSetLinearSolver", 1)) return(1);

    glstr = 0;
    mset = 500;

    data->mykmem = kmem;
    retval = KINSol(kmem, res_vect, glstr, s, s);
    if (check_retval(&retval, "KINSol", 1)) return(1);
    //PrintFinalStats(kmem);
    cout << "retval = " << retval << "\n";

    /* Free memory */
    printf("\nFinal statsistics:\n");
    retval = KINPrintAllStats(kmem, stdout, SUN_OUTPUTFORMAT_TABLE);
    N_VDestroy(res_vect);
    N_VDestroy(s);
    N_VDestroy(c);
    KINFree(&kmem);
    SUNLinSolFree(LS);
    SUNMatDestroy(J);
    free(data);
    SUNContext_Free(&sunctx);
    return 0;
}

int InitialData(int& Nx, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect, double& M)
{
    double h = l / (Nx - 1);
    double x_start = koeff_l * l - l / 8.;
    double x_finish = koeff_l * l + l / 8.;
    int dN = (x_finish - x_start) / h;
    cout << "dN = " << dN << "\n";
    double j = 0;
    M = 2000 * 0.000871523;
    double Y_H2, Y_O2;
    cout << "M = " << M << "\n";
    double W;

    for (int i = 0; i < Nx; i++) {
        x_vect[i] = h * i;
    }

    T_vect[0] = T_start;
    for (int k = 0; k < 9; k++)
        for (int i = 0; i < Nx; i++)
            Y_vect[k * Nx + i] = 0;
    j = 0;
    for (int i = 0; i < Nx; i++) {
        //Y_vect[8 * Nx + i] = Y_N2;
        if (x_vect[i] <= x_start)
        {
            //Y_vect[6 * Nx + i] = 0;
            T_vect[i] = T_start;
        }
        else if (x_vect[i] >= x_finish)
        {
            //Y_vect[6 * Nx + i] = (1. - Y_N2);
            T_vect[i] = T_finish;
        }
        else {
            //Y_vect[6 * Nx + i] = (1. - Y_N2) / dN * j;
            T_vect[i] = (T_finish - T_start) / dN * j + T_start;
            j++;
        }
        //get_Y(Y_vect[6 * Nx + i], Y_H2, Y_O2, Y_N2);
        //Y_vect[i] = Y_H2;
        //Y_vect[2 * Nx + i] = Y_O2;
    }
    j = 0;
    for (int i = 0; i < Nx; i++)
    {
        Y_vect[8 + i * num_gas_species] = Y_N2;
        if (x_vect[i] <= x_start) Y_vect[6 + i * num_gas_species] = 0.;
        else if (x_vect[i] >= x_finish) Y_vect[6 + i * num_gas_species] = (1. - Y_N2);
        else {
            Y_vect[6 + i * num_gas_species] = (1. - Y_N2) / dN * j;
            j++;
        }

        get_Y(Y_vect[6 + i * num_gas_species], Y_H2, Y_O2, Y_N2);
        Y_vect[i * num_gas_species] = Y_H2;
        Y_vect[2 + i * num_gas_species] = Y_O2;
    }
    /* for (int i = 0; i < Nx; i++) {
        T_vect[i] = 1111.31 * tanh((x_vect[i] - koeff_l * l) / 0.064) + 1111.31 + 293.;
    }*/

    for (int i = 0; i < Nx - 1; i++) {
        if (x_vect[i] <= koeff_l * l && x_vect[i + 1] > koeff_l * l)
            return i;
    }
}

void SetParametres1()
{
    A = 6.85 * pow(10, 12);
    Ea = 46.37 * 293.;
    l = 0.4;
    koeff_l = 0.43;
}

void Write_to_file2(string str, ofstream& fout, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double M, int N_x, int number) {
    double x_start, x_finish, D;
    double rho, W, Y_H2, Y_O2;
    fout.open(str + ".dat");
    string title = (boost::format(R"(VARIABLES= "x", "T%d", "Y_H2%d", "Y_O2%d", "Y_H2O%d", "Y_N2%d")") % number % number % number % number % number).str();
    fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout << title << endl;
    for (int i = 0; i < N_x; i++) {
        fout << x_vect[i] << "  " << T_vect[i] << " " << Y_vect[i * num_gas_species]
            << " " << Y_vect[2 + i * num_gas_species]
            << " " << Y_vect[6 + i * num_gas_species]
            << " " << Y_vect[8 + i * num_gas_species] << endl;
    }
    fout.close();
}