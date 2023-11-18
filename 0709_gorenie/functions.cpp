#include "functions.h"

int flag = 0;
#define RTOL  RCONST(1.0e-3)
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
    //M = 0.305867299598755;
    double Y_H2, Y_O2;
    cout << "M = " << M << "\n";
    double W;
    double sumY = 0;
    double A_OH = 0.008;
    double A_H2O = 0.22;
    double A_H2O2 = 8 * pow(10, -4);
    double* Yi = new double[num_gas_species];
    for (int i = 0; i < Nx; i++) {
        x_vect[i] = h * i;
    }

    T_vect[0] = Tstart;
    for (int k = 0; k < 9; k++)
        for (int i = 0; i < Nx; i++)
            Y_vect[k * Nx + i] = 0;
    j = 0;
    //for (int i = 0; i < Nx; i++) {
    //    //T_vect[i] = tanh_T(Tstart, Tfinish, x_center + 0.1 * h, 0.025, x_vect[i]);
    //    if (x_vect[i] <= x_start)
    //    {
    //        T_vect[i] = Tstart;
    //    }
    //    else if (x_vect[i] >= x_finish)
    //    {
    //        T_vect[i] = Tfinish;
    //    }
    //    else {
    //        T_vect[i] = (Tfinish - Tstart) / (x_finish - x_start) * (x_vect[i] - x_start) + Tstart;
    //        j++;
    //    }
    //}
    //j = 0;
    //x_start += 0;
    //x_finish += 0;
    //for (int i = 0; i < Nx; i++)
    //{
    //    sumY = 0;
    //    if (x_vect[i] < x_start) {
    //        for (int k = 0; k < num_gas_species; k++) {
    //            Y_vect[k + i * num_gas_species] = Ystart[k];
    //        }

    //    }
    //    else if (x_vect[i] > x_finish) {
    //        for (int k = 0; k < num_gas_species; k++) {
    //            Y_vect[k + i * num_gas_species] = Yend[k];
    //        }

    //    }
    //    else {
    //        for (int k = 0; k < num_gas_species; k++) {
    //            Y_vect[k + i * num_gas_species] = Ystart[k] + (Yend[k] - Ystart[k]) / (x_finish - x_start) * (x_vect[i] - x_start);
    //        }
    //        j++;
    //    }
    //    for (int k = 0; k < num_gas_species; k++) {

    //    }
    //}
    //double Y_OH, Y_H2O2, Y_HO2;
    //for (int i = 0; i < Nx; i++)
    //{
    //    Y_OH = tanh_spec(Yend[OH], x_center + 0.25 * h, 0.02, x_vect[i], 4);
    //    Y_H2O2 = gauss_func(A_H2O2, x_center + 2.5 * h, h / 2., x_vect[i], 7);
    //    Y_HO2 = gauss_func(A_H2O2, x_center + 2.5 * h, h / 2., x_vect[i], 5);
    //    Y_vect[H2O + i * num_gas_species] = tanh_spec(Yend[H2O], x_center + 0.18 * h, 0.05, x_vect[i], 6);
    //    Y_vect[H2 + i * num_gas_species] = tanh_spec_minus(Ystart[H2], x_center + 0.18 * h, 0.06, x_vect[i], H2);
    //    Y_vect[O2 + i * num_gas_species] = tanh_spec_minus(Ystart[O2], x_center + 0.19 * h, 0.08, x_vect[i], H2);
    //    Y_vect[OH + i * num_gas_species] = Y_OH;
    //    Y_vect[H2O2 + i * num_gas_species] = Y_H2O2;
    //    Y_vect[HO2 + i * num_gas_species] = Y_HO2;
    //    Y_vect[N2 + i * num_gas_species] = Ystart[N2];
    //    sumY = 0;
    //    
    //    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
    //            sumY += Y_vect[k_spec + i * num_gas_species];
    //    }
    //    if (sumY > 1.) {
    //        if (Y_vect[O2 + i * num_gas_species] + 1 - sumY > 0)
    //            Y_vect[O2 + i * num_gas_species] += 1 - sumY;
    //    }
    //    else {
    //        Y_vect[H2O + i * num_gas_species] += 1 - sumY;
    //    }
    //    /*Y_vect[H2 + i * num_gas_species] = (1. - sumY) * 1. / 9.;
    //    Y_vect[O2 + i * num_gas_species] = (1. - sumY) * 8. / 9.;*/
    //    
    //}

    /*for (int i = 0; i < Nx - 1; i++) {
        if (x_vect[i] <= x_center && x_vect[i + 1] > x_center)
            return i;
    }*/
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
        }


        for (int k = 0; k < num_gas_species; k++) {
            Wsmes += Y_vect[k + i * num_gas_species] * my_mol_weight(k);
        }

        for (int k = 0; k < num_gas_species; k++) {

            Y_vect[k + i * num_gas_species] *= my_mol_weight(k) / Wsmes;
        }
    }
    Wsmes = 0;
    for (int k = 0; k < num_gas_species; k++) {
        Wsmes += Ystart[k] * my_mol_weight(k);
    }
    for (int k = 0; k < num_gas_species; k++) {

        Ystart[k] *= my_mol_weight(k) / Wsmes;
    }
    Wsmes = 0;
    for (int k = 0; k < num_gas_species; k++) {
        Wsmes += Yend[k] * my_mol_weight(k);
    }
    for (int k = 0; k < num_gas_species; k++) {

        Yend[k] *= my_mol_weight(k) / Wsmes;
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

void Add_elem(vector<double>& T, vector<double>& Y, vector<double>& x, int& N_x, int& N_center, double b)
{
    int j_t = 1;
    double Y_max = 0, Y_min = Y[N_x - 1];
    j_t = 1;
    Y_min = 0; Y_max = Y_N2;
    double T_max = 0, T_min = T[0];

    for (int i = 0; i < N_x; i++)
    {
        if (T[i] > T_max) T_max = T[i];
        if (T[i] < T_min) T_min = T[i];
    }

    while (j_t < N_x - 2)
    {
        if (fabs(Y[7 + num_gas_species * j_t] - Y[7 + num_gas_species * (j_t - 1)]) > b * (Y_max - Y_min) || fabs(T[j_t] - T[j_t - 1]) > b * (T_max - T_min))
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
    for (int i = 0; i < N_x - 1; i++) {
        if (x[i] <= x_center && x[i + 1] > x_center)
            N_center = i;
    }
    cout << "Add N_center = " << N_center << "\n";
    cout << "Add T N_center = " << T[N_center] << "\n";
    cout << "N_x = " << T.size() << "\n";
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

double gauss_func(double A, double mu, double sigma, double x, int k_spec) {

        return A * exp(-pow(x - mu, 2) / 2 / pow(sigma, 2));
}

double tanh_spec(double A, double mu, double sigma, double x, int k_spec) {

    return A / 2. * tanh((x - mu) / sigma) + A / 2.;
}

double tanh_spec_minus(double A, double mu, double sigma, double x, int k_spec) {

    return A / 2. * tanh(-(x - mu) / sigma) + A / 2.;
}

double tanh_T(double T_start, double T_finish, double mu, double sigma, double x) {

    return (T_finish - T_start) / 2. * tanh((x - mu) / sigma) + (T_finish - T_start) / 2. + T_start;
}

double F_right(IO::ChemkinReader* chemkinReader, double* Yiprev, double* Yi, double* Yinext,
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
        YkVk_slag = YkVk(chemkinReader, k, T, Yi, gradX, Xi);
        Vc -= YkVk_slag;
    }

    for (int k = 0; k < num_gas_species; k++) {
        YkVk_slag = YkVk(chemkinReader, k, T, Yi, gradX, Xi) + Vc * Yi[k];
        slag_diff += rho * YkVk_slag * myget_Cpi(k, T)  * dTdx;
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

    return -(2. / (h + h_left)) *
        (lambda_right * (Tnext - T) / h
            - lambda_left * (T - Tprev) / h_left)
        + Cp * M * dTdx + slag_chem + slag_diff;

}

double F_rightY(IO::ChemkinReader* chemkinReader, double* Yiprev, double* Yi, double* Yinext,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext, double* Xiprev, double* Xi, double* Xinext, double* gradX, double* Y_tmp,
    double* X_tmp,
    double M, const int k_spec, double* ydot, double* wk_add)
{

    double h_left = x - xprev;
    double h = xnext - x;
    double x_l = (x + xprev) / 2.;
    double x_r = (xnext + x) / 2.;
    double rho = get_rho(Yi, T);
    double YkVk_r = 1.;
    double YkVk_l = 1.;
    double rhoYkVk_r = 1.;
    double rhoYkVk_l = 1.;
    double slag_diff = 0.;
    double slag_chem = 0.;
    double dTdx = (h_left / h / (h + h_left) * Tnext + (h - h_left) / h / h_left * T - h / h_left / (h + h_left) * Tprev);
    Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);
    make_averageY(Y_tmp, Yi, Yinext);
    make_averageY(X_tmp, Xi, Xinext);
    get_grad(gradX, Xi, Xinext, x, xnext);

    rho = get_rho(Y_tmp, T / 2. + Tnext / 2.);
    double YkVk_slag = 0;
    double Vc = 0;
    for (int k = 0; k < num_gas_species; k++) {
        YkVk_slag = YkVk(chemkinReader, k, T / 2. + Tnext / 2., Y_tmp, gradX, X_tmp);
        Vc -= YkVk_slag;
    }
    rhoYkVk_r = rho * (YkVk(chemkinReader, k_spec, T / 2. + Tnext/ 2., Y_tmp, gradX, X_tmp) + Y_tmp[k_spec] * Vc);


    make_averageY(Y_tmp, Yi, Yiprev);
    make_averageY(X_tmp, Xi, Xiprev);
    get_grad(gradX, Xiprev, Xi, xprev, x);
    rho = get_rho(Y_tmp, T / 2. + Tprev / 2.);
    Vc = 0;
    for (int k = 0; k < num_gas_species; k++) {
        YkVk_slag = YkVk(chemkinReader, k, T / 2. + Tprev / 2., Y_tmp, gradX, X_tmp);
        Vc -= YkVk_slag;
    }

    rhoYkVk_l = rho * (YkVk(chemkinReader, k_spec, T / 2. + Tprev / 2., Y_tmp, gradX, X_tmp) + Y_tmp[k_spec] * Vc);
    //cout << "\n\n";

    slag_diff = (rhoYkVk_r - rhoYkVk_l) / (x_r - x_l);
    slag_chem = -my_mol_weight(k_spec) * ydot[k_spec];

    //cout << "k = " << name_species[k_spec] << "\n";
    //cout << "d[Xk]/dt = " << ydot[k_spec] << "\n";
    ////cout << "MdY/dx = " << M * (h_left / h / (h + h_left) * Yinext[k_spec] + (h - h_left) / h / h_left * Yi[k_spec] - h / h_left / (h + h_left) * Yiprev[k_spec]) << "\n";
    //cout << "Y slag_diff = " << slag_diff << "\n";
    //cout << "Y slag_chem = " << slag_chem << "\n";
    //cout << "\n\n";
    //if (flag == 1) {
    //   // cout << "Y slag_diff = " << slag_diff << "\n";
    //    //cout << "my_mol_weight  = " << my_mol_weight(k_spec) << "\n";
    //    cout << "rho = const d[Xk]/dt = " << Ith(ydot, k_spec + 1) << "\n";
    //    cout << "drho/dt * Yk / Wk  = " << wk_add[k_spec] << "\n";
    //    cout << "MdY/dx = " << M * (h_left / h / (h + h_left) * Yinext[k_spec] + (h - h_left) / h / h_left * Yi[k_spec] - h / h_left / (h + h_left) * Yiprev[k_spec]) << "\n";
    //    cout << "\n\n";
    //}
    return M * (h_left / h / (h + h_left) * Yinext[k_spec] + (h - h_left) / h / h_left * Yi[k_spec] - h / h_left / (h + h_left) * Yiprev[k_spec])
      + slag_chem + slag_diff;
}

void Write_to_file2(string str, ofstream& fout, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double M, int N_x, int number) {
    double x_start, x_finish, D;
    double rho, W, Y_H2, Y_O2;
    fout.open(str + ".dat");
    string title = (boost::format(R"(VARIABLES= "x", "T%d", "Y_H2_%d", 
"Y_H_%d", "Y_O2_%d", "Y_O_%d", "Y_OH_%d", "Y_HO2_%d", "Y_H2O_%d", "Y_H2O2_%d", "Y_N2_%d")") % number % number % number % number % number
% number % number % number % number % number).str();
    fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout << title << endl;
    for (int i = 0; i < N_x; i++) {
        fout << x_vect[i] << "  " << T_vect[i] << " " << Y_vect[i * num_gas_species]
            << " " << Y_vect[1 + i * num_gas_species]
            << " " << Y_vect[2 + i * num_gas_species]
            << " " << Y_vect[3 + i * num_gas_species]
            << " " << Y_vect[4 + i * num_gas_species]
            << " " << Y_vect[5 + i * num_gas_species]
            << " " << Y_vect[6 + i * num_gas_species]
            << " " << Y_vect[7 + i * num_gas_species]
            << " " << Y_vect[8 + i * num_gas_species] << endl;
    }
    fout.close();
}

void Init_Data(UserData data, int N_x, vector<double>& x_vect,
    vector<double>& T_vect, IO::ChemkinReader* chemkinReader_temp, int NEQ,
    int N_center, double* Y_leftb) {
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
    data->X_tmp = new realtype[num_gas_species];

    data->forward = new double[num_react];
    data->reverse = new double[num_react];
    data->equilib = new double[num_react];

    data->Y_left_bound = new realtype[num_gas_species];
    data->wk_add = new realtype[num_gas_species];
    data->T = new realtype[N_x];
    data->NEQ = NEQ;
    data->Tl = T_vect[0];
    data->T_center = T_vect[N_center];
    cout << "T_center = " << data->T_center << "\n";
    data->N_centr = N_center;
    data->chemkinReader = chemkinReader_temp;
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

void integrate_Y_CVODE(IO::ChemkinReader* chemkinReader_temp, int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb)
{
    SUNContext sunctx;
    void* cvode_mem;
    N_Vector yy, yp, avtol, cons;
    realtype rtol, * yval, * ypval, * atval, * consval;
    realtype t0, tret;
    int iout, retval, retvalr;

    SUNNonlinearSolver NLS;
    double t, tout;
    N_Vector y;
    N_Vector abstol;
    SUNMatrix A;
    SUNLinearSolver LS;
    UserData data;
    data = (UserData)malloc(sizeof * data);
    int NEQ_Y = num_gas_species * (N_x - 2);
    int NEQ = 0 + NEQ_Y;

    int j = 0;

    cvode_mem = NULL;
    cons = yy = yp = avtol = NULL;
    yval = ypval = atval = NULL;
    A = NULL;
    LS = NULL;
    NLS = NULL;

    /* Create SUNDIALS context */
    retval = SUNContext_Create(NULL, &sunctx);

    //data->sunctx = sunctx;
    if (check_retval(&retval, "SUNContext_Create", 1)) return;
    /* Initial conditions */
    y = N_VNew_Serial(NEQ, sunctx);
    cons = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)y, "N_VNew_Serial", 0)) return;
    abstol = N_VNew_Serial(NEQ, sunctx);

    for (int i = 0; i < NEQ; i++) {
        Ith(abstol, i + 1) =  pow(10, -6);
        Ith(cons, i + 1) = 1.0;
    }

    Init_Data(data, N_x, x_vect,
        T_vect, chemkinReader_temp, NEQ,
        N_center, Y_leftb);

    for (int i = 0; i < N_x; i++) {
        cout << "T = " << data->T[i] << "\n";
    }

    yval = N_VGetArrayPointer(y);
    for (int i = 0; i < NEQ; i++) {
        if (i % num_gas_species == 0)
            cout << "\nnn i = " << i / num_gas_species + 1 << "\n";
        yval[i] = Y_vect[i + num_gas_species];
        cout << "yval = " << name_species[i % num_gas_species] << " = " << yval[i] << "\n";
    }

    cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (check_retval((void*)cvode_mem, "CVodeCreate", 0)) return;

    retval = CVodeInit(cvode_mem, func_Y_CVODE, T0, y);
    if (check_retval(&retval, "CVodeInit", 1)) return;

    retval = CVodeSVtolerances(cvode_mem, RTOL, abstol);
    if (check_retval(&retval, "CVodeSVtolerances", 1)) return;
    retval = CVodeSetConstraints(cvode_mem, cons);
    
    A = SUNDenseMatrix(NEQ, NEQ, sunctx);
    if (check_retval((void*)A, "SUNDenseMatrix", 0)) return;

    /* Create dense SUNLinearSolver object for use by CVode */
    LS = SUNLinSol_Dense(y, A, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_Dense", 0)) return;

    retval = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return;

    retval = CVodeSetMaxNumSteps(cvode_mem, 10000);
    retval = CVodeSetMaxHnilWarns(cvode_mem, 50);
    retval = CVodeSetUserData(cvode_mem, data);
    if (check_retval(&retval, "CVODESetUserData", 1)) return;
    /* Create and initialize  y, y', and absolute tolerance vectors. */
    t0 = ZERO;
    data->M = M;
    double rho;
    flag = 1;
    ofstream f_ida;
    double sumY = 0;

    iout = 0;
    double dt =  pow(10, -6);
    tout = dt;
    int iend = 300;
    int number = 100;
    ofstream fout;
    double sum_Y = 0;
    while (iout < iend) {
        retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        cout << "tout = " << tout << "\n";
        for (int i = 0; i < NEQ; i++) {
            Y_vect[num_gas_species + i] = yval[i];
        }

        for (int i = NEQ; i < NEQ + num_gas_species; i++) {
            Y_vect[num_gas_species + i - 1] = Y_vect[i - 1];
        }

        Write_to_file2("detail" + to_string(tout * pow(10, 6)), f_ida, x_vect,
            T_vect, Y_vect, M, N_x, 2);

        if (check_retval(&retval, "IDASolve", 1)) return;

        iout++;
        tout += dt;
    }

    /* Print final statistics to the screen */
    cout << "\nFinal Statistics:\n";
    retval = CVodePrintAllStats(cvode_mem, stdout, SUN_OUTPUTFORMAT_TABLE);

    /* Print final statistics to a file in CSV format */
    //FID = fopen("idaRoberts_dns_stats.csv", "w");
    //retval = IDAPrintAllStats(mem, FID, SUN_OUTPUTFORMAT_CSV);
    //fclose(FID);

    /* check the solution error */
    //retval = check_ans(yy, tret, rtol, avtol);

    /* Free memory */
    free(data);
    N_VDestroy(y);                            /* Free y vector */
    N_VDestroy(abstol);                       /* Free abstol vector */
    CVodeFree(&cvode_mem);                    /* Free CVODE memory */
    SUNLinSolFree(LS);                        /* Free the linear solver memory */
    SUNMatDestroy(A);                         /* Free the matrix memory */
    SUNContext_Free(&sunctx);                 /* Free the SUNDIALS context */

    return;
}

static int func_Y_CVODE(realtype t, N_Vector y, N_Vector ydot, void* user_data)
{
    realtype* yval, * ypval, * rval;
    UserData data;
    realtype* x_cells, * T_vect, * Y_vect, * Tp_vect, * Yp_vect, * Yiprev, * Yinext;
    double* Yi;
    double Temp;
    double M;
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
    int Ncentr = data->N_centr;

    yval = N_VGetArrayPointer(y);
    ypval = N_VGetArrayPointer(ydot);
    double rho;
    for (int i = 1; i < myNx - 1; i++) {
        //cout << "\n\n\n\nin func i = " << i << "\n";
        MakeYvectors(data->chemkinReader, x_cells, Yiprev, Yi, Yinext, yval, data->Y_left_bound, myNx, i, T_vect[0], data->M);
        Get_molar_cons(data->Xi, Yi, T_vect[i]);

        chem_vel(data->forward, data->reverse, data->equilib, data->wk_add, data->M,
            Yi, Yinext,
            x_cells[i], x_cells[i + 1],
            T_vect[i], T_vect[i + 1],
            data->Xi, data->ydot);

        rho = get_rho(Yi, T_vect[i]);
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
           //cout << "yval = " << yval[k_spec] << "\n";
           ypval[k_spec + (i - 1) * num_gas_species] = -F_rightY(data->chemkinReader, 
                data->Yiprev, data->Yi, data->Yinext,
                data->T[i - 1], data->T[i], data->T[i + 1],
                data->x[i - 1], data->x[i], data->x[i + 1],
                data->Xiprev, data->Xi, data->Xinext, data->gradX, data->Y_tmp, data->X_tmp,
                data->M, k_spec, data->ydot, data->wk_add) / rho;
           //cout << "ypval = " << ypval[k_spec] << "\n";
        }
        int jop = 0;
    }
    return(0);
}

int FindBoundary(IO::ChemkinReader* chemkinReader_temp, double* Y_leftbound_inlet, double* Y_leftbound, double* Yinext, double Tl,
    double xi, double xinext, double M)
{
    SUNContext sunctx;
    UserData data;
    realtype fnormtol, scsteptol;
    N_Vector res_vect, s, c;
    int glstr, mset, retval;
    void* kmem;
    SUNMatrix J;
    SUNLinearSolver LS;
    int NEQ = num_gas_species;
    res_vect = NULL;
    s = c = NULL;
    kmem = NULL;
    J = NULL;
    LS = NULL;
    data = NULL;
    //data->sunctx = sunctx;
    /* Create the SUNDIALS context that all SUNDIALS objects require */
    retval = SUNContext_Create(NULL, &sunctx);
    if (check_retval(&retval, "SUNContext_Create", 1)) return(1);

    /* Create serial vectors of length NEQ */

    res_vect = N_VNew_Serial(NEQ, sunctx);

    if (check_retval((void*)res_vect, "N_VNew_Serial", 0)) return(1);

    s = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)s, "N_VNew_Serial", 0)) return(1);

    c = N_VNew_Serial(NEQ, sunctx);
    if (check_retval((void*)c, "N_VNew_Serial", 0)) return(1);

    //SetInitialGuess1(res_vect, data, NEQ);

    N_VConst(ONE, s); /* no scaling */
    data = (UserData)malloc(sizeof * data);

    data->chemkinReader = chemkinReader_temp;
    data->Tl = Tl;
    data->Xi = new realtype[num_gas_species];
    data->Xinext = new realtype[num_gas_species];
    data->gradX = new realtype[num_gas_species];
    data->Y_left_bound = new realtype[num_gas_species];
    data->x = new realtype[2];
    data->Y_tmp = new realtype[num_gas_species];
    data->X_tmp = new realtype[num_gas_species];
    data->Yi = new realtype[num_gas_species];
    data->Yinext = new realtype[num_gas_species];
    data->x[0] = xi;
    data->x[1] = xinext;
    data->M = M;

    for (int k = 0; k < NEQ; k++) {
        Ith(res_vect, k + 1) = Y_leftbound_inlet[k];
        Ith(c, k + 1) = ONE;
        data->Y_left_bound[k] = Y_leftbound_inlet[k];
        data->Yi[k] = Y_leftbound_inlet[k];
        data->Yinext[k] = Yinext[k];
        cout << "Yleft = " << Ith(res_vect, k + 1) << "\n";
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

    retval = KINInit(kmem, left_bound_solve, res_vect);
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

    /* Print out the problem size, solution parameters, initial guess. */
    //PrintHeader(fnormtol, scsteptol);

    /* --------------------------- */

    glstr = 0;
    retval = KINSol(kmem, res_vect, glstr, s, s);

    for (int i = 0; i < num_gas_species; i++) {
        Y_leftbound[i] = Ith(res_vect, i + 1);
    }
    if (check_retval(&retval, "KINSol", 1)) return(1);
    cout << "T res = " << Ith(res_vect, 1) << endl;
    //PrintFinalStats(kmem);
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
    //SUNContext_Free(&sunctx);
    return 0;
}

static int left_bound_solve(N_Vector u, N_Vector f, void* user_data)
{
    int k = 0;
    realtype* Y, * fdata; 
    UserData data;
    data = (UserData)user_data;
    
    Y = N_VGetArrayPointer(u);
    fdata = N_VGetArrayPointer(f);

    double sum = 0.;
    double W = get_W(Y);
    double YkVk_slag = 0;
    double Vc = 0;
    double Dkm = 0;
    Get_mole_fr(data->Xi, Y);
    Get_mole_fr(data->Xinext, data->Yinext);

    get_grad(data->gradX, data->Xi, data->Xinext, data->x[0], data->x[1]);

    make_averageY(data->Y_tmp, Y, data->Yinext);
    make_averageY(data->X_tmp, data->Xi, data->Xinext);

    double rho = get_rho(data->Y_tmp, data->Tl);

    for (int k = 0; k < num_gas_species; k++) {
        sum = 0;
        for (int j = 0; j < num_gas_species; j++) {
            if (j != k) {
                sum += data->X_tmp[j]
                    / Dij_func(data->chemkinReader, k, j, data->Tl, data->Y_tmp);
            }
            YkVk_slag = YkVk(data->chemkinReader, k, data->Tl, data->Y_tmp, data->gradX, data->X_tmp);
            Vc -= YkVk_slag;
        }

        Dkm = (1. - Y[k]) / sum;
        YkVk_slag = -my_mol_weight(k) / W * Dkm * data->gradX[k]
            + Y[k] * Vc;

        fdata[k] = data->Y_left_bound[k] - Y[k] - rho * YkVk_slag / data->M;
    }

    return(0);
}