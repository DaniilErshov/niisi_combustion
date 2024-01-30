#include "functions.h"

int flag = 0;
#define RTOL  RCONST(1.0e-5)
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(pow(10, -6))     /* first output time      */
#define TMULT RCONST(1.006)     /* output time factor     */
#define NOUT  1800
#define ZERO  RCONST(0.0)

int InitialData(int& Nx, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect, double& M, double Tstart, double Tfinish, double* Ystart, double* Yend)
{
    double h = l / (Nx - 1);
    x_center = 4 * h;
    double x_start = 0.1;
    double x_finish = 0.2;
    int dN = (x_finish - x_start) / h;
    cout << "dN = " << dN << "\n";
    double j = 0;
    M = 300 * 8.4955E-004;
    double Y_H2, Y_O2;
    cout << "M = " << M << "\n";
    double W;
    double sumY = 0;
    double* Yi = new double[num_gas_species];
    T_vect[0] = Tstart;

    x_vect[0] = 0; x_vect[1] = 0.0500; x_vect[2] = 0.0980; x_vect[3] = 0.1000; x_vect[4] = 0.1167; x_vect[5] = 0.1333; x_vect[6] = 0.1500; x_vect[7] = 0.1667; 
    x_vect[8] = 0.1833; x_vect[9] = 0.2000;
    x_vect[10] = 0.2500; x_vect[11] = 0.3000;

    for (int i = 0; i < Nx; i++) {
        if (x_vect[i] <= x_start)
        {
            T_vect[i] = Tstart;
        }
        else if (x_vect[i] >= x_finish)
        {
            T_vect[i] = Tfinish;
        }
        else {
            T_vect[i] = (Tfinish - Tstart) / (x_finish - x_start) * (x_vect[i] - x_start) + Tstart;
            j++;
        }
    }
    double Wsmes = 0;
    for (int i = 0; i < Nx; i++) {
        sumY = 0;
        Wsmes = 0;

        for (int k = 0; k < num_gas_species; k++) {
            Wsmes += Y_vect[k + i * num_gas_species] * my_mol_weight(k);
        }

        for (int k = 0; k < num_gas_species; k++) {

            Y_vect[k + i * num_gas_species] *= my_mol_weight(k) / Wsmes;
        }
    }
    double T_min = Tfinish;
    double Tcenter = (Tstart + Tfinish) / 2.;
    double i_center = 0;
    for (int i = 0; i < Nx - 1; i++) {
        if (abs(T_vect[i] - Tcenter) < T_min) {
            i_center = i;
            T_min = abs(T_vect[i] - Tcenter);
        }
    }
    delete[] Yi;
    return i_center;
}

int InitialData2(int& Nx, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect, double& M, double Tstart, double Tfinish, double* Ystart, double* Yend)
{
    double h = l / (Nx - 1);
    double x_start = 0.1;
    double x_finish = 0.2;
    int dN = (x_finish - x_start) / h;

    double j = 0;
    cout << "M = " << M << "\n";

    for (int i = 0; i < Nx; i++) {
        x_vect[i] = h * i;
    }
    x_vect = { 0.0,
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

    for (int i = 0; i < Nx; i++) {

        if (x_vect[i] <= x_start)
        {
            T_vect[i] = Tstart;
        }
        else if (x_vect[i] >= x_finish)
        {
            T_vect[i] = Tfinish;
        }
        else {
            T_vect[i] = (Tfinish - Tstart) / (x_finish - x_start) * (x_vect[i] - x_start) + Tstart;
            j++;
        }
    }


    for (int i = 0; i < Nx; i++)
    {
        if (x_vect[i] < x_start) {
            for (int k = 0; k < num_gas_species; k++) {
                Y_vect[k + i * num_gas_species] = Ystart[k];
            }

        }
        else if (x_vect[i] > x_finish) {
            for (int k = 0; k < num_gas_species; k++) {
                Y_vect[k + i * num_gas_species] = Yend[k];
            }

        }
        else {
            for (int k = 0; k < num_gas_species; k++) {
                Y_vect[k + i * num_gas_species] = Ystart[k] + (Yend[k] - Ystart[k]) / (x_finish - x_start) * (x_vect[i] - x_start);
            }
            j++;
        }
    }

    double T_min = Tfinish;
    double Tcenter = (Tstart + Tfinish) / 2.;
    double i_center = 0;
    for (int i = 0; i < Nx - 1; i++) {
        if (abs(T_vect[i] - Tcenter) < T_min) {
            i_center = i;
            T_min = abs(T_vect[i] - Tcenter);
        }
    }
    return i_center;
}

bool Add_elem(vector<double>& T, vector<double>& Y, vector<double>& x, int& N_x, int& N_center, double b, int number, int number_start, double T_center, int& j_t)
{
    double Y_max = 5 * pow(10,-3), Y_min = 0.;
    Y_min = 1; Y_max = 0;
    double T_max = 0, T_min = T[0];
    double Yi[9]{};
    int Nx_add = 0;
    for (int i = 0; i < N_x; i++)
    {
        if (T[i] > T_max) T_max = T[i];
        if (T[i] < T_min) T_min = T[i];
        if (Y[H2O2 + num_gas_species * j_t] > Y_max) Y_max = Y[H2O2 + num_gas_species * j_t];
        if (Y[H2O2 + num_gas_species * j_t] < Y_min) Y_min = Y[H2O2 + num_gas_species * j_t];
    }
    if (j_t == N_x - 2) j_t = 1;
    while (j_t < N_x - 2)
    {
        if (fabs(T[j_t] - T[j_t - 1]) > b * (T_max - T_min) || fabs(Y[H2O2 + num_gas_species * j_t] - Y[H2O2 + num_gas_species * (j_t - 1)]) > b * (Y_max - Y_min))
        {
            T.insert(T.begin() + j_t, (T[j_t] + T[j_t - 1]) / 2.);
            for (int k = 0; k < num_gas_species; k++) {
                Y.insert(Y.begin() + k + num_gas_species * j_t, (Y[k + num_gas_species * j_t + k] + Y[k + num_gas_species * (j_t - 1)]) / 2.);
            }
            x.insert(x.begin() + j_t, (x[j_t] + x[j_t - 1]) / 2.);
            N_x++;
            j_t++;
            Nx_add++;
        }
        j_t++;

        if (Nx_add == 10) {
            N_x = x.size();
            for (int i = 0; i < T.size(); i++) {
                if (T[i] == T_center) {
                    N_center = i;
                }
            }
            for (int i = 0; i < Y.size(); i++) {
                if (Y[i] < 0) Y[i] = 0;
            }

            return true;
        }
        //cout << "j_t = " << j_t << "\n";
    }
    if (Nx_add > 0) {
        N_x = x.size();
        for (int i = 0; i < T.size(); i++) {
            if (T[i] == T_center) {
                N_center = i;
            }
        }
        for (int i = 0; i < Y.size(); i++) {
            if (Y[i] < 0) Y[i] = 0;
        }
        return true;
    }
    if (Nx_add == 0)
    {
        for (int k = 0; k < number; k++) {
            T.push_back(T[N_x - 1]);
            x.push_back(x[N_x - 1] + 1.5 * (x[N_x - 1] - x[N_x - 2]));
            for (int i = 0; i < num_gas_species; i++) {
                Yi[i] = Y[Y.size() - num_gas_species + i];
            }
            for (int i = 0; i < num_gas_species; i++) {
                Y.push_back(Yi[i]);
            }
            N_x++;
        }

        if (number_start == 1) {
            for (int i = 0; i < num_gas_species; i++) {
                Yi[i] = Y[i];
            }
            for (int i = num_gas_species - 1; i >= 0; i--) {
                Y.insert(Y.begin(), Yi[i]);
            }
            T.insert(T.begin(), T[0]);
            x.insert(x.begin(), x[0] - (x[1] - x[0]));

        }
        N_x = x.size();

        for (int i = 0; i < T.size(); i++) {
            if (T[i] == T_center) {
                N_center = i;
            }
        }

        for (int i = 0; i < Y.size(); i++) {
            if (Y[i] < 0) Y[i] = 0;
        }
        return false;
    }


}

void Add_elem_simple(vector<double>& T, vector<double>& Y, vector<double>& x, int& N_x, int& N_center, double b, int number, int number_start, double T_center)
{
    double Y_max = 5 * pow(10, -3), Y_min = 0.;
    Y_min = 1; Y_max = 0;
    double T_max = 0, T_min = T[0];
    double Yi[9]{};
    int Nx_add = 0;
    int j_t = 1;
    for (int i = 0; i < N_x; i++)
    {
        if (T[i] > T_max) T_max = T[i];
        if (T[i] < T_min) T_min = T[i];
        if (Y[H2O2 + num_gas_species * j_t] > Y_max) Y_max = Y[H2O2 + num_gas_species * j_t];
        if (Y[H2O2 + num_gas_species * j_t] < Y_min) Y_min = Y[H2O2 + num_gas_species * j_t];
    }
    while (j_t < N_x - 2)
    {
        if (fabs(T[j_t] - T[j_t - 1]) > b * (T_max - T_min) || fabs(Y[H2O2 + num_gas_species * j_t] - Y[H2O2 + num_gas_species * (j_t - 1)]) > b * (Y_max - Y_min))
        {
            T.insert(T.begin() + j_t, (T[j_t] + T[j_t - 1]) / 2.);
            for (int k = 0; k < num_gas_species; k++) {
                Y.insert(Y.begin() + k + num_gas_species * j_t, (Y[k + num_gas_species * j_t + k] + Y[k + num_gas_species * (j_t - 1)]) / 2.);
            }
            x.insert(x.begin() + j_t, (x[j_t] + x[j_t - 1]) / 2.);
            N_x++;
            j_t++;

        }
        j_t++;

        //cout << "j_t = " << j_t << "\n";
    }
    for (int k = 0; k < number; k++) {
        T.push_back(T[N_x - 1]);
        x.push_back(x[N_x - 1] + 2 * (x[N_x - 1] - x[N_x - 2]));
        for (int i = 0; i < num_gas_species; i++) {
            Yi[i] = Y[Y.size() - num_gas_species + i];
        }
        for (int i = 0; i < num_gas_species; i++) {
            Y.push_back(Yi[i]);
        }
        N_x++;
    }

    if (number_start == 1) {
        for (int i = 0; i < num_gas_species; i++) {
            Yi[i] = Y[i];
        }
        for (int i = num_gas_species - 1; i >= 0; i--) {
            Y.insert(Y.begin(), Yi[i]);
        }
        T.insert(T.begin(), T[0]);
        x.insert(x.begin(), x[0] - (x[1] - x[0]));

    }
    N_x = x.size();

    for (int i = 0; i < T.size(); i++) {
        if (T[i] == T_center) {
            N_center = i;
        }
    }

    for (int i = 0; i < Y.size(); i++) {
        if (Y[i] < 0) Y[i] = 0;
    }
}

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

void MakeYvectors(UserData data,
    double* Y, int myNx, int i, double Tl) {
    //cout << "i = " << i << "\n";
    for (int j = 0; j < num_gas_species; j++) {

        data->Yi[j] = Y[j + (i - 1) * num_gas_species];
        if (i == 1) data->Yiprev[j] = data->Y_left_bound[j];
        else  data->Yiprev[j] = Y[j + (i - 2) * num_gas_species];
         /*cout << "k = " << name_species[j] << "\n";
         cout << " Yiprev[j] = " << data->Yiprev[j] << "\n";
         cout << " Yi[j] = " << data->Yi[j] << "\n";*/

        if (i == myNx - 2)  data->Yinext[j] = data->Yi[j];
        else  data->Yinext[j] = Y[j + (i) * num_gas_species];
        //cout << " Yinext[j] = " << data->Yinext[j] << "\n\n";
    }
}

double get_M(double* Yiprev, double* Yi, double* Yinext,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext, double* Xiprev, double* Xi, double* Xinext, double* gradX, double* Y_tmp, double* X_tmp,
    double M, double* ydot, double* wk_add)
{
    double h_left = x - xprev;
    double h = xnext - x;
    double Cp = Cp_all(T, Yi);
    double rho = get_rho(Yi, T);
    double W = get_W(Yi);
    Get_mole_fr(Xiprev, Yiprev);
    Get_mole_fr(Xi, Yi);
    Get_mole_fr(Xinext, Yinext);
    get_grad(gradX, Xi, Xinext, x, xnext);
    double dTdx = (h_left / h / (h + h_left) * Tnext + (h - h_left) / h / h_left * T - h / h_left / (h + h_left) * Tprev);
    double slag_diff = 0;
    double slag_chem = 0.;

    double Vc = 0;
    double YkVk_slag = 0;
    for (int k = 0; k < num_gas_species; k++) {
        YkVk_slag = YkVk(k, T, Yi, gradX, Xi);
        Vc -= YkVk_slag;
    }

    for (int k = 0; k < num_gas_species; k++) {
        YkVk_slag = YkVk(k, T, Yi, gradX, Xi) + Vc * Yi[k];
        slag_diff += rho * YkVk_slag * myget_Cpi(k, T) * dTdx;
        slag_chem += ydot[k] * myget_Hi(k, T) * my_mol_weight(k);
    }

    //cout << "M  T = " << T << "\n";
    //cout << "M  slag_chem = " << slag_chem << "\n";
    //for (int k = 0; k < num_gas_species; k++) {
    //    cout << "k name = " << name_species[k] << "\n";
    //    cout << "M _rho_ d[X]/dt = " << ydot[k] << "\n";
    //}

    //cout << "\nM  slag_diff = " << slag_diff << "\n\n\n\n";
    make_averageY(X_tmp, Xi, Xinext);
    double lambda_right = Lambda_All(X_tmp, (Tnext + T) / 2.);
    //cout << "lambda = " << lambda_right << "\n";
    //cout << "Cp = " << Cp << "\n";
    make_averageY(X_tmp, Xi, Xiprev);
    double lambda_left = Lambda_All(X_tmp, (T + Tprev) / 2.);
    return ((2. / (h + h_left)) *
        (lambda_right * (Tnext - T) / h
            - lambda_left * (T - Tprev) / h_left)
        - slag_chem - slag_diff) / Cp / dTdx;

}

double F_right(UserData data,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext)
{
    double h_left = x - xprev;
    double h = xnext - x;
    double Cp = Cp_all(T, data->Yi);
    double rho = get_rho(data->Yi, T);
    double W = get_W(data->Yi);
    get_grad(data->gradX, data->Xi, data->Xinext, x, xnext);
    double dTdx = (h_left / h / (h + h_left) * Tnext + (h - h_left) / h / h_left * T - h / h_left / (h + h_left) * Tprev);
    double slag_diff = 0;
    double slag_chem = 0.;

    double Vc = 0;
    double YkVk_slag = 0;
    for (int k = 0; k < num_gas_species; k++) {
        data->YkVk[k] = YkVk(k, T, data->Yi, data->gradX, data->Xi);
        Vc -= data->YkVk[k];
    }

    for (int k = 0; k < num_gas_species; k++) {
        slag_diff += rho * data->YkVk[k] * myget_Cpi(k, T)  * dTdx;
        slag_chem += data->ydot[k] * myget_Hi(k, T) * my_mol_weight(k);
    }

    make_averageY(data->X_tmp, data->Xi, data->Xinext);
    double lambda_right = Lambda_All(data->X_tmp, (Tnext + T) / 2.);
    make_averageY(data->X_tmp, data->Xi, data->Xiprev);
    double lambda_left = Lambda_All(data->X_tmp, (T + Tprev) / 2.);

    return -(2. / (h + h_left)) *
        (lambda_right * (Tnext - T) / h
            - lambda_left * (T - Tprev) / h_left)
        + Cp * data->M * dTdx + slag_chem + slag_diff;

}

double F_rightY(UserData data, int k_spec,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext)
{

    double h_left = x - xprev;
    double h = xnext - x;
    double x_l = (x + xprev) / 2.;
    double x_r = (xnext + x) / 2.;
    double rho = get_rho(data->Yi, T);
    double YkVk_r = 1.;
    double YkVk_l = 1.;
    double rhoYkVk_r = 1.;
    double rhoYkVk_l = 1.;
    double slag_diff = 0.;
    double slag_chem = 0.;
    double dTdx = (h_left / h / (h + h_left) * Tnext + (h - h_left) / h / h_left * T - h / h_left / (h + h_left) * Tprev);

    rhoYkVk_r = data->rho_r * (data->YkVk_r[k_spec] + data->Y_tmp_r[k_spec] * data->Vc_r);
    //cout << "data->YkVk[k_spec] =  " << data->YkVk[k_spec] << "\n";


    rhoYkVk_l = data->rho_l * (data->YkVk_l[k_spec] + data->Y_tmp_l[k_spec] * data->Vc_l);
    //cout << "\n\n";

    slag_diff = (rhoYkVk_r - rhoYkVk_l) / (x_r - x_l);
    slag_chem = -my_mol_weight(k_spec) * data->ydot[k_spec];
    //cout << "slag_chem =  " << slag_chem << "\n";
    //cout << "slag_diff =  " << slag_diff << "\n";

    return data->M * (h_left / h / (h + h_left) * data->Yinext[k_spec] + (h - h_left) / h / h_left * data->Yi[k_spec] - h / h_left / (h + h_left) * data->Yiprev[k_spec])
      + slag_chem + slag_diff;
}

void Write_to_file2(string str, ofstream& fout, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, vector<double>& Yp_vect, double M, int N_x, int number) {
    double x_start, x_finish, D;
    double rho;
    fout.open(str + ".dat");
    string title = (boost::format(R"(VARIABLES= "x", "T%d", "Y_H2_%d", 
"Y_H_%d", "Y_O2_%d", "Y_O_%d", "Y_OH_%d", "Y_HO2_%d", "Y_H2O_%d", "Y_H2O2_%d", "Y_N2_%d")") % number % number % number % number % number
% number % number % number % number % number).str();
    fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout << title << endl;
    double* Yi = new double[num_gas_species];
    for (int i = 0; i < N_x; i++) {
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            Yi[k_spec] = Y_vect[k_spec + i * num_gas_species];
        }
        rho = get_rho(Yi, T_vect[i]);

        if (number == 1)
        fout << x_vect[i] << "  " << T_vect[i] << " " <<Y_vect[i * num_gas_species]
            << " " << Y_vect[1 + i * num_gas_species]
            << " " << Y_vect[2 + i * num_gas_species]
            << " " << Y_vect[3 + i * num_gas_species]
            << " " << Y_vect[4 + i * num_gas_species]
            << " " << Y_vect[5 + i * num_gas_species]
            << " " << Y_vect[6 + i * num_gas_species]
            << " " << Y_vect[7 + i * num_gas_species]
            << " " << Y_vect[8 + i * num_gas_species] << endl; 
        if (number == 2)
            fout << x_vect[i] << "  " << T_vect[i] << " " << rho * Yp_vect[i * num_gas_species]
            << " " << rho * Yp_vect[1 + i * num_gas_species]
            << " " << rho * Yp_vect[2 + i * num_gas_species]
            << " " << rho * Yp_vect[3 + i * num_gas_species]
            << " " << rho * Yp_vect[4 + i * num_gas_species]
            << " " << rho * Yp_vect[5 + i * num_gas_species]
            << " " << rho * Yp_vect[6 + i * num_gas_species]
            << " " << rho * Yp_vect[7 + i * num_gas_species]
            << " " << rho * Yp_vect[8 + i * num_gas_species] << endl;
    }
   
    fout.close();
    delete[] Yi;
}

void Init_Data(UserData data, int N_x, vector<double>& x_vect,
    vector<double>& T_vect, int NEQ, int NEQ_Y,
    int N_center, double* Y_leftb) {
    data->Nx = N_x;
    data->x = new realtype[N_x];

    data->Yi = new realtype[num_gas_species];
    data->Yiprev = new realtype[num_gas_species];
    data->Yinext = new realtype[num_gas_species];
    data->YkVk_r = new realtype[num_gas_species];
    data->YkVk_l = new realtype[num_gas_species];

    data->gradX_r = new realtype[num_gas_species];
    data->gradX_l = new realtype[num_gas_species];

    data->X_tmp_r = new realtype[num_gas_species];
    data->X_tmp_l = new realtype[num_gas_species];

    data->Y_tmp_r = new realtype[num_gas_species];
    data->Y_tmp_l = new realtype[num_gas_species];

    data->Xiprev = new realtype[num_gas_species];
    data->Xi = new realtype[num_gas_species];
    data->Xinext = new realtype[num_gas_species];
    data->gradX = new realtype[num_gas_species];
    data->Y_tmp = new realtype[num_gas_species];
    data->X_tmp = new realtype[num_gas_species];
    data->YkVk = new realtype[num_gas_species];

    data->forward = new double[num_react];
    data->reverse = new double[num_react];
    data->equilib = new double[num_react];

    data->Y_left_bound = new realtype[num_gas_species];
    data->wk_add = new realtype[num_gas_species];
    data->T = new realtype[N_x];
    data->NEQ = NEQ;
    data->NEQ_Y = NEQ_Y;
    data->Tl = T_vect[0];
    data->T_center = T_vect[N_center];
    cout << "T_center = " << data->T_center << "\n";
    data->N_centr = N_center;
    data->y = new realtype[num_gas_species];
    data->ydot = new realtype[num_gas_species];
    data->N_m = 1 - 1;

    for (int i = 0; i < N_x; i++) {
        data->x[i] = x_vect[i];
        data->T[i] = T_vect[i];
    }

    for (int i = 0; i < num_gas_species; i++) {
        data->y[i] = 0;
        data->ydot[i] = 0;
        data->Y_left_bound[i] = Y_leftb[i];
    }
}

int integrate_Y_IDA(int N_x, vector<double>& x_vect,
        vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb, double t_fix){
    void* mem;
    N_Vector yy, yp, avtol, cons;
    realtype rtol, * yval, * ypval, * atval, *consval;
    realtype t0, tout, tret;
    int iout, retval, retvalr;
    SUNMatrix A;
    SUNLinearSolver LS;
    SUNNonlinearSolver NLS;
    SUNContext ctx;
    UserData data;
    data = (UserData)malloc(sizeof * data);
    int NEQ_Y = num_gas_species * (N_x - 2);
    int NEQ = 0 + NEQ_Y;
    vector<double> Yp_vect(N_x * num_gas_species);
    
    int j = 0;
    
    mem = NULL;
    cons = yy = yp = avtol = NULL;
    yval = ypval = atval = NULL;
    A = NULL;
    LS = NULL;
    NLS = NULL;
    /* Create SUNDIALS context */
    retval = SUNContext_Create(NULL, &ctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return(1);
    
    Init_Data(data, N_x, x_vect,
        T_vect, NEQ, NEQ_Y,
        N_center, Y_leftb);


    /* Allocate N-vectors. */
    yy = N_VNew_Serial(NEQ, ctx);
    if (check_retval((void*)yy, "N_VNew_Serial", 0)) return(1);
    
    yp = N_VClone(yy);
    if (check_retval((void*)yp, "N_VNew_Serial", 0)) return(1);
    avtol = N_VClone(yy);
    if (check_retval((void*)avtol, "N_VNew_Serial", 0)) return(1);
    cons = N_VClone(yy);
    if (check_retval((void*)cons, "N_VNew_Serial", 0)) return(1);
        
    
    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rtol = RCONST(1.0e-4);
    atval = N_VGetArrayPointer(avtol);
    consval = N_VGetArrayPointer(cons);
    int number_spec = 0;
    for (int i = 0; i < NEQ; i++) {
        consval[i] = 1.0;   /*constraint*/
        yval[i] = Y_vect[i + num_gas_species];
        ypval[i] = 0;
        atval[i] = RCONST(1.0e-5);
    }
    /* Create and initialize  y, y', and absolute tolerance vectors. */
    t0 = ZERO;
    data->M = M;
    double rho;
    ofstream f_ida;
    double sumY = 0;
        
       
    MakeYvectors(data,
        yval, N_x, 1, T_vect[0]);

    Get_molar_cons(data->Xi, data->Yiprev, T_vect[0]);

    chem_vel(data->forward, data->reverse, data->equilib, T_vect[0], data->Xi, data->ydot);
    Get_mole_fr(data->Xiprev, data->Yiprev); Get_mole_fr(data->Xi, data->Yi); Get_mole_fr(data->Xinext, data->Yinext);

    cout << "M_in_IDA_Y = " << M << "\n";
    cout << "Nx _ Y  =  " << data->Nx << "\n";
    for (int i = 1; i < data->Nx - 1; i++) {
        MakeYvectors(data, 
            yval, N_x, i, T_vect[i]);

        Get_molar_cons(data->Xi, data->Yi, T_vect[i]);

        chem_vel(data->forward, data->reverse, data->equilib, T_vect[i], data->Xi, data->ydot);
        Get_mole_fr(data->Xiprev, data->Yiprev); Get_mole_fr(data->Xi, data->Yi); Get_mole_fr(data->Xinext, data->Yinext);
        rho = get_rho(data->Yi, T_vect[i]);
        sumY = 0;
        
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            ypval[k_spec + (i - 1) * num_gas_species] = -F_rightY(data, k_spec,
                data->T[i - 1], data->T[i], data->T[i + 1],
                data->x[i - 1], data->x[i], data->x[i + 1]) / rho;
        }
    
    }
    for (int i = 0; i < NEQ; i++) {
        Yp_vect[i] = ypval[i];
    }
    
    Write_to_file2("A_proizvod_initial", f_ida, x_vect,
        T_vect, Y_vect, Yp_vect, M, N_x, 2);
    
    /* Call IDACreate and IDAInit to initialize IDA memory */
    mem = IDACreate(ctx);
    if (check_retval((void*)mem, "IDACreate", 0)) return(1);
    
    retval = IDAInit(mem, func_Y_IDA, t0, yy, yp);
    if (check_retval(&retval, "IDAInit", 1)) return(1);
    /* Call IDASVtolerances to set tolerances */
    
    retval = IDASVtolerances(mem, rtol, avtol);
    if (check_retval(&retval, "IDASVtolerances", 1)) return(1);
    
    retval = IDASetConstraints(mem, cons);
    retval = IDASetUserData(mem, data);
    if (check_retval(&retval, "IDASetUserData", 1)) return(1);
    retval = IDASetMaxNumSteps(mem, 1000);
    /* Create dense SUNMatrix for use in linear solves */
    A = SUNDenseMatrix(NEQ, NEQ, ctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) return(1);
    
    /* Create dense SUNLinearSolver object */
    LS = SUNLinSol_Dense(yy, A, ctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return(1);
    
    /* Attach the matrix and linear solver */
    retval = IDASetLinearSolver(mem, LS, A);
    if (check_retval(&retval, "IDASetLinearSolver", 1)) return(1);
    
    NLS = SUNNonlinSol_Newton(yy, ctx);
    if (check_retval((void*)NLS, "SUNNonlinSol_Newton", 0)) return(1);
    
    /* Attach the nonlinear solver */
    retval = IDASetNonlinearSolver(mem, NLS);
    if (check_retval(&retval, "IDASetNonlinearSolver", 1)) return(1);
    //IDASetInitStep(mem, pow(10, -7));
    /* In loop, call IDASolve, print results, and test for error.
        Break out of loop when NOUT preset output times have been reached. */
    
    iout = 0;
    double tout1 = t_fix;
    tout = tout1;
    int iend = 2000000;
    int number = 10;
    ofstream fout;
    ofstream foutw;
    double Y_H2, Y_O2;
    double W, w_dot;
    double sum_Y = 0;
    double* Yvect_rho = new double[num_gas_species];
    
    for (int i = 0; i < num_gas_species * N_x; i++)
        Yp_vect[i] = 0;
    double sum = 100;
    while (sum > 1) {
        retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);

        for (int i = 0; i < NEQ; i++) {
            Y_vect[i + num_gas_species] = yval[i];
            Yp_vect[i + num_gas_species] = ypval[i];
        }
    
        for (int i = NEQ; i < NEQ + num_gas_species - 1; i++) {
            Y_vect[i] = Y_vect[i - num_gas_species];
        }
        for (int i = NEQ; i < N_x * num_gas_species; i++) {
            Yp_vect[i] = 0;        
        }

        sum = 0;
        for (int i = 0; i < N_x - 1; i++) {
            for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
                data->Yi[k_spec] = Y_vect[k_spec + i * num_gas_species];
            }
            rho = get_rho(data->Yi, T_vect[i]);
            sum += rho * (abs(Yp_vect[i * num_gas_species])
                + abs(Yp_vect[1 + i * num_gas_species])
                + abs(Yp_vect[2 + i * num_gas_species])
                + abs(Yp_vect[3 + i * num_gas_species])
                + abs(Yp_vect[4 + i * num_gas_species])
                + abs(Yp_vect[5 + i * num_gas_species])
                + abs(Yp_vect[6 + i * num_gas_species])
                + abs(Yp_vect[7 + i * num_gas_species])
                + abs(Yp_vect[8 + i * num_gas_species]));
        }
    
        iout++;
        tout += tout1;
    }
    cout << "IDA_Y sum_Y = " << sum << "\n";
    for (int i = 0; i < NEQ; i++) {
        Y_vect[i + num_gas_species] = yval[i];
        if (Y_vect[i + num_gas_species] < 0) Y_vect[i + num_gas_species] = 0;
    }

    for (int i = NEQ ; i <  NEQ + num_gas_species ; i++) {
        Y_vect[i + num_gas_species] = Y_vect[i];
    }

    //Write_to_file2("detail/detail" + to_string(tout * pow(10, 6)), f_ida, x_vect,
    //    T_vect, Y_vect, Yp_vect, M, N_x, 1);

    //Write_to_file2("proizvod/proizvod" + to_string(tout * pow(10, 6)), f_ida, x_vect,
    //    T_vect, Y_vect, Yp_vect, M, N_x, 2);
    //if (check_retval(&retval, "IDASolve", 1)) return(1);

    /* Print final statistics to the screen */
    cout << "\nFinal Statistics:\n";
    retval = IDAPrintAllStats(mem, stdout, SUN_OUTPUTFORMAT_TABLE);
    
    /* Free memory */
    delete[](Yvect_rho);
    IDAFree(&mem);
    SUNNonlinSolFree(NLS);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(avtol);
    N_VDestroy(yy);
    N_VDestroy(yp);
    N_VDestroy(cons);
    SUNContext_Free(&ctx);
    
    return(retval);
}
    
static int func_Y_IDA(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
    realtype* yval, * ypval, * rval;
    UserData data;
    realtype* x_cells, * T_vect, * Y_vect, * Tp_vect, * Yp_vect, * Yiprev, * Yinext;
    double* Yi;
    double Temp;
    double M;
    double sumY = 0;
    data = (UserData)user_data;
    T_vect = data->T;
    x_cells = data->x;
    int j;
    Yi = data->Yi;
    x_cells = data->x;
    T_vect = data->T;
    Yi = data->Yi;
    Yiprev = data->Yiprev;
    Yinext = data->Yinext;
    int myNx = data->Nx;
    int myNeq = data->NEQ;
    int myNm = data->N_m;
    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rval = N_VGetArrayPointer(rr);
    int Ncentr = data->N_centr;
    double rho_r, rho_l;
    double rho = get_rho(Yi, T_vect[0]);
    double Vc = 0;
    double Vc_r, Vc_l;
    double rhoYkVk_r, rhoYkVk_l;
    MakeYvectors(data, yval, myNx, 1, T_vect[0]);
    Get_mole_fr(data->Xiprev, data->Yiprev);
    Get_mole_fr(data->Xi, data->Yi);
    Get_mole_fr(data->Xinext, data->Yinext);

    make_averageY(data->Y_tmp, Yiprev, Yi);
    make_averageY(data->X_tmp, data->Xiprev, data->Xi);

    get_grad(data->gradX, data->Xiprev, data->Xi, x_cells[0], x_cells[1]);

    //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
    //    Vc -= YkVk(k_spec, T_vect[0], data->Y_tmp, data->gradX, data->X_tmp);
    //}
    for (int i = 1; i < myNx - 1; i++) {
        //cout << "\n\n\n\nin func i = " << i << "\n";
        MakeYvectors(data, yval, myNx, i, T_vect[i]);
        Get_molar_cons(data->Xi, Yi, T_vect[i]);
        chem_vel(data->forward, data->reverse, data->equilib, T_vect[i], data->Xi, data->ydot);

        
        Get_mole_fr(data->Xiprev, data->Yiprev); Get_mole_fr(data->Xi, data->Yi); Get_mole_fr(data->Xinext, data->Yinext);


        find_diff_slag(data, data->Yi, data->Yinext, data->Xi, data->Xinext, data->YkVk_r, data->Y_tmp_r, data->X_tmp_r, data->gradX_r, data->rho_r, data->Vc_r, i);

        find_diff_slag(data, data->Yiprev, data->Yi, data->Xiprev, data->Xi, data->YkVk_l, data->Y_tmp_l, data->X_tmp_l, data->gradX_l, data->rho_l, data->Vc_l, i - 1);

        rho = get_rho(Yi, T_vect[i]);
        sumY = 0;
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            //cout << "name = " << name_species[k_spec] << "\n";
            //if (k_spec != N2)
            rval[k_spec + (i - 1) * num_gas_species] = rho * ypval[k_spec + (i - 1) * num_gas_species] + F_rightY(data, k_spec,
                data->T[i - 1], data->T[i], data->T[i + 1],
                data->x[i - 1], data->x[i], data->x[i + 1]);


            //cout << "yval = " << yval[k_spec + (i - 1) * num_gas_species] << "\n";
            //cout << "ypval = " << ypval[k_spec + (i - 1) * num_gas_species] << "\n\n\n\n\n";
            //if (k_spec == N2) {
            //    rval[k_spec + (i)*num_gas_species] = sumY - 1.;
            //}
        }
        //cout << "yval = " << yval[num_gas_species] << "\n";
    }
    return(0);
}

int integrate_All_IDA(int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb, int iter, double t_fix) {

    void* mem;
    N_Vector yy, yp, avtol, cons;
    realtype rtol, * yval, * ypval, * atval, * consval;
    realtype t0, tout, tret;
    int iout, retval, retvalr;
    SUNMatrix A;
    SUNLinearSolver LS;
    SUNNonlinearSolver NLS;
    SUNContext ctx;
    UserData data;
    data = (UserData)malloc(sizeof * data);
    int NEQ_Y = num_gas_species * (N_x - 2);
    int NEQ = NEQ_Y + N_x - 2;
    vector<double> Yp_vect(num_gas_species * (N_x));
    vector<double> Tp_vect(N_x);
    int j = 0;

    mem = NULL;
    cons = yy = yp = avtol = NULL;
    yval = ypval = atval = NULL;
    A = NULL;
    LS = NULL;
    NLS = NULL;
    /* Create SUNDIALS context */
    retval = SUNContext_Create(NULL, &ctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

    Init_Data(data, N_x, x_vect,
        T_vect, NEQ, NEQ_Y,
        N_center, Y_leftb);


    /* Allocate N-vectors. */
    yy = N_VNew_Serial(NEQ, ctx);
    if (check_retval((void*)yy, "N_VNew_Serial", 0)) return(1);

    yp = N_VClone(yy);
    if (check_retval((void*)yp, "N_VNew_Serial", 0)) return(1);
    avtol = N_VClone(yy);
    if (check_retval((void*)avtol, "N_VNew_Serial", 0)) return(1);
    cons = N_VClone(yy);
    if (check_retval((void*)cons, "N_VNew_Serial", 0)) return(1);


    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rtol = RCONST(1.0e-3);
    atval = N_VGetArrayPointer(avtol);
    consval = N_VGetArrayPointer(cons);
    int number_spec = 0;

    data->M = M;
    cout << "M_in_IDA = " << M << "\n";
    for (int i = 0; i < NEQ_Y; i++) {
        consval[i] = 1.0;   /*constraint*/
        yval[i] = Y_vect[i + num_gas_species];
        ypval[i] = 0;
        number_spec = i % num_gas_species;
        if (number_spec == H2 || number_spec == O2 || number_spec == H2O || number_spec == N2)
            atval[i] = RCONST(1.0e-3);
        if (number_spec == H || number_spec == O || number_spec == OH)
            atval[i] = RCONST(1.0e-4);
        if (number_spec == HO2 || number_spec == H2O2)
            atval[i] = RCONST(1.0e-5);
        //cout << "i = " << i << "   = " << yval[i] << endl;
        //if ((i + 1) % num_gas_species == 0) cout << endl;
    }
    MakeYvectors(data,
        yval, N_x, N_center, T_vect[N_center]);

    Get_molar_cons(data->Xi, data->Yi, T_vect[N_center]);
    chem_vel(data->forward, data->reverse, data->equilib, T_vect[N_center], data->Xi, data->ydot);

    Get_mole_fr(data->Xiprev, data->Yiprev); Get_mole_fr(data->Xi, data->Yi); Get_mole_fr(data->Xinext, data->Yinext);
    for (int i = NEQ_Y; i < NEQ; i++) {
        consval[i] = 1.0;   /*constraint*/
        atval[i] = RCONST(1.0e0);
        yval[i] = T_vect[i - NEQ_Y + 1];
        if (i - NEQ_Y + 1 == N_center) yval[i] = data->M;
        //cout << "i = " << i << "   = " << yval[i] << endl;
     
    }
    /* Create and initialize  y, y', and absolute tolerance vectors. */
    t0 = ZERO;

    double rho;
    ofstream f_ida;
    double sumY = 0;

    flag = 1;

    cout << "Nx_ALL = " << data->Nx << "\n";
    for (int i = 1; i < data->Nx - 1; i++) {

        MakeYvectors(data,
            yval, N_x, i, T_vect[i]);

        Get_molar_cons(data->Xi, data->Yi, T_vect[i]);
        chem_vel(data->forward, data->reverse, data->equilib, T_vect[i], data->Xi, data->ydot);

        Get_mole_fr(data->Xiprev, data->Yiprev); Get_mole_fr(data->Xi, data->Yi); Get_mole_fr(data->Xinext, data->Yinext);
        rho = get_rho(data->Yi, T_vect[i]);

        double Cp = get_Cp(num_gas_species, data->Yi, T_vect[i]);

        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            ypval[k_spec + (i - 1) * num_gas_species] = -F_rightY(data, k_spec,
                data->T[i - 1], data->T[i], data->T[i + 1],
                data->x[i - 1], data->x[i], data->x[i + 1]
                ) / rho;
        }
        if (i != N_center) {
            ypval[i + NEQ_Y - 1] = -F_right(data,
                data->T[i - 1], data->T[i], data->T[i + 1],
                data->x[i - 1], data->x[i], data->x[i + 1]) / rho / Cp;
        }
        else {
            ypval[i + NEQ_Y - 1] = -F_right(data,
                data->T[i - 1], data->T[i], data->T[i + 1],
                data->x[i - 1], data->x[i], data->x[i + 1]);
        }
             
    }

    //cout << NEQ << " = " << Ith(res_vect, NEQ + 1) << endl;

    /* Call IDACreate and IDAInit to initialize IDA memory */
    mem = IDACreate(ctx);
    if (check_retval((void*)mem, "IDACreate", 0)) return(1);

    retval = IDAInit(mem, func_All_IDA, t0, yy, yp);
    if (check_retval(&retval, "IDAInit", 1)) return(1);
    /* Call IDASVtolerances to set tolerances */

    retval = IDASVtolerances(mem, rtol, avtol);
    if (check_retval(&retval, "IDASVtolerances", 1)) return(1);

    retval = IDASetConstraints(mem, cons);
    retval = IDASetUserData(mem, data);
    if (check_retval(&retval, "IDASetUserData", 1)) return(1);
    retval = IDASetMaxNumSteps(mem, 1000);
    /* Create dense SUNMatrix for use in linear solves */
    A = SUNDenseMatrix(NEQ, NEQ, ctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) return(1);

    /* Create dense SUNLinearSolver object */
    LS = SUNLinSol_Dense(yy, A, ctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return(1);

    /* Attach the matrix and linear solver */
    retval = IDASetLinearSolver(mem, LS, A);
    if (check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

    NLS = SUNNonlinSol_Newton(yy, ctx);
    if (check_retval((void*)NLS, "SUNNonlinSol_Newton", 0)) return(1);

    /* Attach the nonlinear solver */
    retval = IDASetNonlinearSolver(mem, NLS);
    if (check_retval(&retval, "IDASetNonlinearSolver", 1)) return(1);
    //IDASetInitStep(mem, pow(10, -7));

    /* In loop, call IDASolve, print results, and test for error.
        Break out of loop when NOUT preset output times have been reached. */

    iout = 0;
    double tout1 = t_fix;
    tout = tout1;
    int iend = 2000000;
    int number = 10;
    ofstream fout;
    ofstream foutw;
    double Y_H2, Y_O2;
    double W, w_dot;
    double sum_Y = 0;
    double alpha = 1;

    for (int i = 0; i < num_gas_species * N_x; i++)
        Yp_vect[i] = 0;
    double nevyaz_Y = 1000, nevyaz_T = 1000;

    data->N_m = 0;
    retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);
    if (iout % number == 0) {
        for (int i = 0; i < NEQ_Y; i++) {
                Y_vect[i + num_gas_species] = yval[i];
                if (Y_vect[i + num_gas_species] < 0) Y_vect[i + num_gas_species] = 0;
                Yp_vect[i + num_gas_species] = ypval[i];
            
        }
        for (int i = NEQ_Y; i < N_x * num_gas_species; i++) {
            Yp_vect[i] = 0;
        }

        Tp_vect[0] = Tp_vect[N_x - 1] = 0;

        for (int i = NEQ_Y; i < NEQ; i++) {
            Tp_vect[i - NEQ_Y + 1] = ypval[i];
        }

        for (int i = 0; i < num_gas_species; i++)
            Y_vect[(num_gas_species) * (N_x - 1) + i] = Y_vect[(num_gas_species) * (N_x - 2) + i];

        for (int i = NEQ_Y; i < NEQ; i++) {
            T_vect[i - NEQ_Y + 1] = yval[i];
            if (i - NEQ_Y + 1 == N_center) {
                T_vect[i - NEQ_Y + 1] = data->T_center;
                M = yval[i];
            }
        }
        T_vect[N_x - 1] = T_vect[N_x - 2];

        nevyaz_Y = 0, nevyaz_T = 0.;
        for (int i = 0; i < N_x; i++) {
            for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
                data->Yi[k_spec] = Y_vect[k_spec + i * num_gas_species];
            }
            rho = get_rho(data->Yi, T_vect[i]);
            nevyaz_T += rho * (abs(Tp_vect[i]));
            nevyaz_Y += rho * (abs(Yp_vect[i * num_gas_species])
                + abs(Yp_vect[1 + i * num_gas_species])
                + abs(Yp_vect[2 + i * num_gas_species])
                + abs(Yp_vect[3 + i * num_gas_species])
                + abs(Yp_vect[4 + i * num_gas_species])
                + abs(Yp_vect[5 + i * num_gas_species])
                + abs(Yp_vect[6 + i * num_gas_species])
                + abs(Yp_vect[7 + i * num_gas_species])
                + abs(Yp_vect[8 + i * num_gas_species]));
        }
        cout << "tout = " << tout << "\n";
        cout << "nevyaz_Y  = " << nevyaz_Y << "\n";
        cout << "nevyaz_T  = " << nevyaz_T << "\n";
        cout << "M = " << M << "\n";
        if (check_retval(&retval, "IDASolve", 1)) return(1);
        iout++;
        tout1 *= 1.2;
        tout += tout1;
    }

    //cout << "\n\n\n\n";
    /* Print final statistics to the screen */
    cout << "\nFinal Statistics:\n";

    retval = IDAPrintAllStats(mem, stdout, SUN_OUTPUTFORMAT_TABLE);
    /*Write_to_file2("detail/detail" + to_string(iter), f_ida, x_vect,
                            T_vect, Y_vect, Yp_vect, M, N_x, 1);*/
    /* Free memory */
    free(data);
    IDAFree(&mem);
    SUNNonlinSolFree(NLS);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(avtol);
    N_VDestroy(yy);
    N_VDestroy(yp);
    N_VDestroy(cons);
    SUNContext_Free(&ctx);
    return(retval);
}

void find_diff_slag(UserData data, double* Yi, double* Yinext, 
    double* Xi, double* Xinext, double* Ykvk_side, double* Y_tmp_side, double* X_tmp_side, double* gradX_side, double& rho_side, double& Vc_side, int i) {
    make_averageY(Y_tmp_side, Yi, Yinext);
    make_averageY(X_tmp_side, Xi, Xinext);
    get_grad(gradX_side, Xi, Xinext, data->x[i], data->x[i + 1]);
    rho_side = get_rho(Y_tmp_side, (data->T[i] + data->T[i + 1]) / 2.);
    Vc_side = 0;
    for (int k = 0; k < num_gas_species; k++) {
        Ykvk_side[k] = YkVk(k, (data->T[i] + data->T[i + 1]) / 2., Y_tmp_side, gradX_side, X_tmp_side);
        Vc_side -= Ykvk_side[k];
    }
}

static int func_All_IDA(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
    realtype* yval, * ypval, * rval;
    UserData data;
    realtype* x_cells, * T_vect, * Y_vect, * Tp_vect, * Yp_vect, * Yiprev, * Yinext;
    double* Yi;
    double Temp;
    double M;
    double sumY = 0;
    int i_T;
    data = (UserData)user_data;
    double Cp;
    double T_curr, T_prev, T_next;
    T_vect = data->T;
    x_cells = data->x;
    int j;
    Yi = data->Yi;
    x_cells = data->x;
    T_vect = data->T;
    Yi = data->Yi;
    Yiprev = data->Yiprev;
    Yinext = data->Yinext;
    int myNx = data->Nx;
    int myNeq = data->NEQ;
    int myNm = data->N_m;
    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rval = N_VGetArrayPointer(rr);
    int Ncentr = data->N_centr;
    double alpha = 1;
    double rho = get_rho(Yi, data->Tl);
    double Vc = 0;

    data->M = yval[data->N_centr + data->NEQ_Y - 1];
    for (int i = 1; i < myNx - 1; i++) {

        i_T = i + data->NEQ_Y - 1;

        if (i == data->N_centr) {
            T_curr = data->T_center;
            T_prev = yval[i_T - 1];
            T_next  = yval[i_T + 1];
        }
        else if (i == data->N_centr - 1) {
            T_curr = yval[i_T];
            T_prev = yval[i_T - 1];
            T_next = data->T_center;
        }
        else if (i == data->N_centr + 1) {
            T_curr = yval[i_T];
            T_prev = data->T_center;
            T_next = yval[i_T + 1];
        }
        else {
            T_curr = yval[i_T];
            T_prev = yval[i_T - 1];
            T_next = yval[i_T + 1];
        }
        if (i == myNx - 2) T_next = T_curr;
        if (i == 1) T_prev = data->Tl;
        //cout << "in func i = " << i << "\n";
        MakeYvectors(data, yval, myNx, i, data->Tl);

        Get_molar_cons(data->Xi, Yi, T_curr);
        chem_vel(data->forward, data->reverse, data->equilib,
            T_curr,
            data->Xi, data->ydot);


        Get_mole_fr(data->Xiprev, data->Yiprev); Get_mole_fr(data->Xi, data->Yi); Get_mole_fr(data->Xinext, data->Yinext);

        find_diff_slag(data, data->Yi, data->Yinext, data->Xi, data->Xinext, data->YkVk_r, data->Y_tmp_r, data->X_tmp_r, data->gradX_r, data->rho_r, data->Vc_r, i);

        find_diff_slag(data, data->Yiprev, data->Yi, data->Xiprev, data->Xi, data->YkVk_l, data->Y_tmp_l, data->X_tmp_l, data->gradX_l, data->rho_l, data->Vc_l, i - 1);

        rho = get_rho(Yi, T_curr);
        sumY = 0;
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
                rval[k_spec + (i - 1) * num_gas_species] = rho * ypval[k_spec + (i - 1) * num_gas_species] + F_rightY(data, k_spec,
                    T_prev, T_curr, T_next,
                    data->x[i - 1], data->x[i], data->x[i + 1]);
               
        }
        Cp = get_Cp(num_gas_species, data->Yi, T_curr);
        if (i != data->N_centr)
            rval[i_T] = rho * Cp * ypval[i_T] + F_right(data,
                T_prev, T_curr, T_next,
                data->x[i - 1], data->x[i], data->x[i + 1]
                );
        else 
            rval[i_T] = F_right(data,
                T_prev, T_curr, T_next,
                data->x[i - 1], data->x[i], data->x[i + 1]);
    }
    return(0);
}

int Integrate_Kinsol(int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb, int iter)
{
    SUNContext sunctx;
    realtype fnormtol, scsteptol;
    N_Vector res_vect, s, c;
    int glstr, mset, retval;
    void* kmem;
    SUNMatrix J;
    SUNLinearSolver LS;
    realtype* yval, * cval, * sval;
    UserData data;
    data = (UserData)malloc(sizeof * data);

    int NEQ_Y = num_gas_species * (N_x - 2);
    int NEQ = NEQ_Y + N_x - 2;

    res_vect = NULL;
    s = c = NULL;
    kmem = NULL;
    J = NULL;
    LS = NULL;
    data->M = M;
    /* Create the SUNDIALS context that all SUNDIALS objects require */
    retval = SUNContext_Create(NULL, &sunctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return(1);
    res_vect = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)res_vect, "N_VNew_Serial", 0)) return(1);

    s = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)s, "N_VNew_Serial", 0)) return(1);

    c = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)c, "N_VNew_Serial", 0)) return(1);

    N_VConst(ONE, s); /* no scaling */

    /* User data */
    Init_Data(data, N_x, x_vect,
        T_vect, NEQ, NEQ_Y,
        N_center, Y_leftb);

    yval = N_VGetArrayPointer(res_vect);
    cval = N_VGetArrayPointer(c);
    sval = N_VGetArrayPointer(s);
    int number_spec = 0;
    for (int i = 0; i < NEQ_Y; i++) {
        yval[i] = Y_vect[i + num_gas_species];
        cval[i] = 0.0;
        number_spec = i % num_gas_species;
        if (number_spec == H2 || number_spec == O2 || number_spec == H2O || number_spec == N2)
            sval[i] = pow(10, -2);
        if (number_spec == H || number_spec == O || number_spec == OH)
            sval[i] = pow(10, -1);
        if (number_spec == HO2 || number_spec == H2O2)
            sval[i] = 1.;
    }

    for (int i = NEQ_Y; i < NEQ; i++) {
        yval[i] = T_vect[i - NEQ_Y + 1];
        cval[i] = 1.0; 
        sval[i] = pow(10, -5);
        if (i - NEQ_Y + 1 == N_center) yval[i] = data->M;
    }
    fnormtol = FTOL; scsteptol = STOL;

   /* for (int i = 1; i < data->Nx - 1; i++) {

        MakeYvectors(data,
            yval, N_x, i, T_vect[i]);

        Get_molar_cons(data->Xi, data->Yi, T_vect[i]);


        chem_vel(data->forward, data->reverse, data->equilib, T_vect[i], data->Xi, data->ydot);
        Get_mole_fr(data->Xiprev, data->Yiprev); Get_mole_fr(data->Xi, data->Yi); Get_mole_fr(data->Xinext, data->Yinext);


        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            cout << "F_Y = " << F_rightY(data, k_spec,
              T_vect[i - 1], T_vect[i], T_vect[i + 1],
                x_vect[i - 1], x_vect[i], x_vect[i + 1]
            ) << "\n";
        }
        cout << "F_ = " << - F_right(data,
            data->T[i - 1], data->T[i], data->T[i + 1],
            data->x[i - 1], data->x[i], data->x[i + 1]) << "\n";

    }*/

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

    retval = KINInit(kmem, func_kinsol, res_vect);
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
    mset = 1000;

    retval = KINSol(kmem, res_vect, glstr, s, s);
    if (check_retval(&retval, "KINSol", 1)) return(1);
    cout << "retval = " << retval << "\n";
    for (int i = 0; i < NEQ_Y; i++) {
        Y_vect[i + num_gas_species] = yval[i];
    }

    for (int i = 0; i < num_gas_species; i++)
        Y_vect[(num_gas_species) * (N_x - 1) + i] = Y_vect[(num_gas_species) * (N_x - 2) + i];

    for (int i = NEQ_Y; i < NEQ; i++) {
        T_vect[i - NEQ_Y + 1] = yval[i];
        if (i - NEQ_Y + 1 == N_center) {
            T_vect[i - NEQ_Y + 1] = data->T_center;
            M = yval[i];
        }
    }
    T_vect[N_x - 1] = T_vect[N_x - 2];
    cout << "M = " << M << "\n";
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

static int func_kinsol(N_Vector u, N_Vector f, void* user_data)
    {
        realtype* yval, * ypval, * rval;
        UserData data;
        realtype* x_cells, * T_vect, * Y_vect, * Tp_vect, * Yp_vect, * Yiprev, * Yinext;
        double* Yi;
        double Temp;
        double M;
        double sumY = 0;
        int i_T;
        data = (UserData)user_data;
        double Cp;
        double T_curr, T_prev, T_next;
        T_vect = data->T;
        x_cells = data->x;
        int j;
        Yi = data->Yi;
        x_cells = data->x;
        T_vect = data->T;
        Yi = data->Yi;
        Yiprev = data->Yiprev;
        Yinext = data->Yinext;
        int myNx = data->Nx;
        int myNeq = data->NEQ;
        int myNm = data->N_m;
        yval = N_VGetArrayPointer(u);
        rval = N_VGetArrayPointer(f);
        int Ncentr = data->N_centr;
        double rho;
        data->M = yval[data->N_centr + data->NEQ_Y - 1];
        //cout << "M = " << data->M << "\n";

        for (int i = 1; i < myNx - 1; i++) {
            //cout << "\n\n\ni =  " << i << "\n";
            i_T = i + data->NEQ_Y - 1;

            if (i == data->N_centr) {
                T_curr = data->T_center;
                T_prev = yval[i_T - 1];
                T_next = yval[i_T + 1];
            }
            else if (i == data->N_centr - 1) {
                T_curr = yval[i_T];
                T_prev = yval[i_T - 1];
                T_next = data->T_center;
            }
            else if (i == data->N_centr + 1) {
                T_curr = yval[i_T];
                T_prev = data->T_center;
                T_next = yval[i_T + 1];
            }
            else {
                T_curr = yval[i_T];
                T_prev = yval[i_T - 1];
                T_next = yval[i_T + 1];
            }
            if (i == myNx - 2) T_next = T_curr;
            if (i == 1) T_prev = data->Tl;
            //cout << "in func i = " << i << "\n";
            MakeYvectors(data, yval, myNx, i, data->Tl);
            Get_molar_cons(data->Xi, Yi, T_curr);

            chem_vel(data->forward, data->reverse, data->equilib, 
                T_curr, 
                data->Xi, data->ydot);

            Get_mole_fr(data->Xiprev, data->Yiprev); Get_mole_fr(data->Xi, data->Yi); Get_mole_fr(data->Xinext, data->Yinext);

            find_diff_slag(data, data->Yi, data->Yinext, data->Xi, data->Xinext, data->YkVk_r, data->Y_tmp_r, data->X_tmp_r, data->gradX_r, data->rho_r, data->Vc_r, i);

            find_diff_slag(data, data->Yiprev, data->Yi, data->Xiprev, data->Xi, data->YkVk_l, data->Y_tmp_l, data->X_tmp_l, data->gradX_l, data->rho_l, data->Vc_l, i - 1);

            rho = get_rho(Yi, T_curr);
            sumY = 0;
            for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
                //cout << "name = " << name_species[k_spec] << "\n";
                //cout << "name = " << Yi[k_spec] << "\n";
                rval[k_spec + (i - 1) * num_gas_species] = F_rightY(data, k_spec,
                    T_prev, T_curr, T_next,
                    data->x[i - 1], data->x[i], data->x[i + 1]) /rho ;
                //cout << "rval = " << rval[k_spec + (i - 1) * num_gas_species] << "\n";
            }
            //cout << "\n\n\n";
            Cp = get_Cp(num_gas_species, data->Yi, T_curr);

            if (i != data->N_centr)
                rval[i_T] = F_right(data,
                    T_prev, T_curr, T_next,
                    data->x[i - 1], data->x[i], data->x[i + 1]) / rho / Cp;
            else
                rval[i_T] = F_right(data,
                    T_prev, T_curr, T_next,
                    data->x[i - 1], data->x[i], data->x[i + 1]);
            //cout << "rval_T = " << rval[i_T] << "\n\n";
        }
        //cout << "\n\n\n\n\n";
        int jop;
        return(0);
    }

int Find_final_state_IDA(double& Tinitial, double& Tend, double* Y_vect, double* Y_end)
{
    void* mem;
    N_Vector yy, yp, avtol, cons;
    realtype rtol, * yval, * ypval, * atval;
    realtype t0, tout, tret;
    int iout, retval, retvalr;
    SUNMatrix A;
    SUNLinearSolver LS;
    SUNNonlinearSolver NLS;
    SUNContext ctx;
    UserData data;
    data = (UserData)malloc(sizeof * data);
    int NEQ_Y = num_gas_species;
    int NEQ = NEQ_Y + 1;

    mem = NULL;
    cons = yy = yp = avtol = NULL;
    yval = ypval = atval = NULL;
    A = NULL;
    LS = NULL;
    NLS = NULL;
    /* Create SUNDIALS context */
    retval = SUNContext_Create(NULL, &ctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

    data->y = new realtype[num_gas_species];
    data->ydot = new realtype[num_gas_species];
    data->Xi = new realtype[num_gas_species];
    data->Yi = new realtype[num_gas_species];
    data->Yiprev = new realtype[num_gas_species];
    data->forward = new double[num_react];
    data->reverse = new double[num_react];
    data->equilib = new double[num_react];

    data->NEQ = NEQ;
    data->Tl = Tinitial;

    for (int i = 0; i < num_gas_species; i++) {
        data->Yi[i] = Y_vect[i];
        data->Yiprev[i] = Y_vect[i];
        //cout << i << " = " << data->Yi[i] << endl;
    }

    /* Allocate N-vectors. */
    yy = N_VNew_Serial(NEQ, ctx);
    if (check_retval((void*)yy, "N_VNew_Serial", 0)) return(1);
    yp = N_VClone(yy);
    if (check_retval((void*)yp, "N_VNew_Serial", 0)) return(1);
    avtol = N_VClone(yy);
    if (check_retval((void*)avtol, "N_VNew_Serial", 0)) return(1);
    cons = N_VClone(yy);
    if (check_retval((void*)cons, "N_VNew_Serial", 0)) return(1);

    /* Create and initialize  y, y', and absolute tolerance vectors. */
    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rtol = RCONST(1.0e-12);
    atval = N_VGetArrayPointer(avtol);
    /* Integration limits */
    t0 = ZERO;

    //cout << NEQ << " = " << Ith(res_vect, NEQ + 1) << endl;
    for (int i = 1; i < NEQ; i++) {
        Ith(avtol, i) = RCONST(1.0e-10);
        Ith(yy, i) = Y_vect[i - 1];
    }

    Ith(avtol, NEQ) = RCONST(1.0e-10);
    Ith(yy, NEQ) = Tinitial * 5;

    double sumY = 0;
    Get_molar_cons(data->Xi, Y_vect, Tinitial);

    chem_vel(data->forward, data->reverse, data->equilib, Tinitial, data->Xi, data->ydot);

    double rho = get_rho(Y_vect, Tinitial);

    for (int i = 1; i < NEQ; i++) {
        Ith(yp, i) = data->ydot[i - 1];
    }

    double sum = 0;
    double Cp = get_Cp(num_gas_species, Y_vect, Tinitial);
    sum = get_enthalpy(num_gas_species, Y_vect, Tinitial);
    double sum2 = get_enthalpy(num_gas_species, Y_vect, Ith(yy, NEQ));

    Ith(yp, NEQ) = (sum - sum2) / rho / Cp;

    cout << "Ith(yp, NEQ) = " << Ith(yp, NEQ) << "\n";
    /* Call IDACreate and IDAInit to initialize IDA memory */
    mem = IDACreate(ctx);
    if (check_retval((void*)mem, "IDACreate", 0)) return(1);

    retval = IDAInit(mem, func_final_state, t0, yy, yp);
    if (check_retval(&retval, "IDAInit", 1)) return(1);
    /* Call IDASVtolerances to set tolerances */

    retval = IDASVtolerances(mem, rtol, avtol);
    if (check_retval(&retval, "IDASVtolerances", 1)) return(1);
    for (int i = 1; i <= NEQ; i++) {
        Ith(cons, i) = 1.0;
    }

    retval = IDASetConstraints(mem, cons);
    retval = IDASetUserData(mem, data);
    if (check_retval(&retval, "IDASetUserData", 1)) return(1);
    retval = IDASetMaxNumSteps(mem, 40000);
    /* Create dense SUNMatrix for use in linear solves */
    A = SUNDenseMatrix(NEQ, NEQ, ctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) return(1);

    /* Create dense SUNLinearSolver object */
    LS = SUNLinSol_Dense(yy, A, ctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return(1);

    /* Attach the matrix and linear solver */
    retval = IDASetLinearSolver(mem, LS, A);
    if (check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

    NLS = SUNNonlinSol_Newton(yy, ctx);
    if (check_retval((void*)NLS, "SUNNonlinSol_Newton", 0)) return(1);

    /* Attach the nonlinear solver */
    retval = IDASetNonlinearSolver(mem, NLS);
    if (check_retval(&retval, "IDASetNonlinearSolver", 1)) return(1);

    /* In loop, call IDASolve, print results, and test for error.
       Break out of loop when NOUT preset output times have been reached. */

    iout = 0;
    double tout1 = pow(10, 1);
    tout = tout1;
    int iend = 500;
    int number = 100;
    ofstream fout;
    double Y_H2, Y_O2;
    double W, w_dot;
    double sum_Y = 0;

    while (iout < iend) {
        retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);
        if (check_retval(&retval, "IDASolve", 1)) return(1);

        iout++;
        tout += tout1;
    }
    for (int i = 0; i < num_gas_species; i++) {
        Y_end[i] = Ith(yy, i + 1);
    }

    Tend = Ith(yy, num_gas_species + 1);
    /* Print final statistics to the screen */
    cout << "\nFinal Statistics:\n";
    retval = IDAPrintAllStats(mem, stdout, SUN_OUTPUTFORMAT_TABLE);

    /* Print final statistics to a file in CSV format */
    //FID = fopen("idaRoberts_dns_stats.csv", "w");
    //retval = IDAPrintAllStats(mem, FID, SUN_OUTPUTFORMAT_CSV);
    //fclose(FID);

    /* check the solution error */
    //retval = check_ans(yy, tret, rtol, avtol);

    /* Free memory */
    IDAFree(&mem);
    SUNNonlinSolFree(NLS);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(avtol);
    N_VDestroy(yy);
    N_VDestroy(yp);
    N_VDestroy(cons);
    SUNContext_Free(&ctx);

    return(retval);
}

static int func_final_state(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
    realtype* yval, * ypval, * rval;
    UserData data;
    realtype* x_cells, * T_vect, * Y_vect, * Tp_vect, * Yp_vect;
    double* Yi;
    double Temp;
    double M;
    data = (UserData)user_data;
    T_vect = data->T;
    x_cells = data->x;
    int j;
    Yi = data->Yi;
    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rval = N_VGetArrayPointer(rr);

    for (int j = 0; j < num_gas_species; j++) {
        Yi[j] = yval[j];
    }
    Temp = yval[num_gas_species];
    Get_molar_cons(data->Xi, Yi, Temp);
    chem_vel(data->forward, data->reverse, data->equilib, Temp, data->Xi, data->ydot);

    double rho = get_rho(Yi, Temp);
    for (j = 0; j < num_gas_species; j++) {
        rval[j] = rho * ypval[j] - data->ydot[j] * phyc.mol_weight[j];
    }
    double sum = 0;

    double Cp = get_Cp(num_gas_species, Yi, Temp);
    double sum2;

    sum2 = get_enthalpy(num_gas_species, data->Yiprev, data->Tl);

    sum = get_enthalpy(num_gas_species, Yi, Temp);

    rval[num_gas_species] = rho * Cp * ypval[num_gas_species] + (sum - sum2);
    return(0);
}

void makeYstart(double koeff_topl, double* Ystart) {

    Ystart[0] = koeff_topl * 2.0 / (koeff_topl * 2. + 1. + 3.76);
    Ystart[2] = 1.0 / (koeff_topl * 2. + 1. + 3.76);
    Ystart[8] = 3.76 / (koeff_topl * 2. + 1. + 3.76);
    double Wsmes = 0;

    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        Wsmes += Ystart[k_spec] * phyc.mol_weight[k_spec];
    }

    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        Ystart[k_spec] *= phyc.mol_weight[k_spec] / Wsmes;
        cout << "Yi = " << Ystart[k_spec] << "\n";
    }
}