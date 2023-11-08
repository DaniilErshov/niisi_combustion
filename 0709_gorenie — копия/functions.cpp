#include "functions.h"

int flag = 0;

int InitialData(int& Nx, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect, double& M, double Tstart, double Tfinish, double* Ystart, double* Yend)
{
    double h = l / (Nx - 1);
    x_center = 8 * h;
    double x_start = x_center -  2 * h;
    double x_finish = x_center +  h;
    int dN = (x_finish - x_start) / h;
    cout << "dN = " << dN << "\n";
    double j = 0;
   //M = 5522.86532 *  get_rho(Ystart, Tstart);
    M = 0.454534374246507;
    double Y_H2, Y_O2;
    cout << "M = " << M << "\n";
    double W;
    double sumY = 0;
    double A_OH = 0.0055;
    double A_H2O2 = 6.5 * pow(10, -4);
    for (int i = 0; i < Nx; i++) {
        x_vect[i] = h * i;
    }

    T_vect[0] = Tstart;
    for (int k = 0; k < 9; k++)
        for (int i = 0; i < Nx; i++)
            Y_vect[k * Nx + i] = 0;
    j = 0;
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
    j = 0;
    x_start += h;
    x_finish += h;
    for (int i = 0; i < Nx; i++)
    {
        sumY = 0;
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
    double Y_OH, Y_H2O2, Y_HO2;
    for (int i = 0; i < Nx; i++)
    {
        Y_OH = tanh_OH(Yend[4], x_center, 0.04, x_vect[i], 4);
        Y_H2O2 = gauss_func(A_H2O2, (x_finish - x_start) / 2. + x_start + -1. * h, (x_finish - x_start) / 6., x_vect[i], 7);
        Y_HO2 = gauss_func(A_H2O2, (x_finish - x_start) / 2. + x_start + -1. * h, (x_finish - x_start) / 6., x_vect[i], 5);

        if (Y_vect[2 + i * num_gas_species] - Y_OH > 0.) {
            Y_vect[2 + i * num_gas_species] -= Y_OH;
            Y_vect[4 + i * num_gas_species] = Y_OH;
        }
        else if (Y_vect[6 + i * num_gas_species] - Y_OH > 0.) {
            Y_vect[6 + i * num_gas_species] -= Y_OH;
            Y_vect[4 + i * num_gas_species] = Y_OH;
        }

        if (Y_vect[2 + i * num_gas_species] - Y_H2O2 > 0.) {
            Y_vect[2 + i * num_gas_species] -= Y_H2O2;
            Y_vect[7 + i * num_gas_species] = Y_H2O2;
        }
        else if (Y_vect[6 + i * num_gas_species] - Y_H2O2 > 0.) {
            Y_vect[6 + i * num_gas_species] -= Y_H2O2;
            Y_vect[7 + i * num_gas_species] = Y_H2O2;
        }

        if (Y_vect[2 + i * num_gas_species] - Y_HO2 > 0.) {
            Y_vect[2 + i * num_gas_species] -= Y_HO2;
            Y_vect[5 + i * num_gas_species] = Y_HO2;
        }
        else if (Y_vect[6 + i * num_gas_species] - Y_HO2 > 0.) {
            Y_vect[6 + i * num_gas_species] -= Y_HO2;
            Y_vect[5 + i * num_gas_species] = Y_HO2;
        }

    }

    for (int i = 0; i < Nx - 1; i++) {
        if (x_vect[i] <= x_center && x_vect[i + 1] > x_center)
            return i;
    }

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

double tanh_OH(double A, double mu, double sigma, double x, int k_spec) {

    return A / 2. * tanh((x - mu) / sigma) + A / 2.;
}

double F_right(IO::ChemkinReader* chemkinReader, double* Yi, double* Yinext,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext, double* Xiprev, double* Xi, double* Xinext, double* gradX, double* Y_tmp, 
    double M, double* ydot, double* wk_add)
{
    double h_left = x - xprev;
    double h = xnext - x;
    double Cp = Cp_all(T, Yi);
    double rho = get_rho(Yi, T);
    double W = get_W(Yi);
    Get_mole_fr(Xi, Yi);
    Get_mole_fr(Xinext, Yinext);
    get_grad(gradX, Xi, Xinext, x, xnext);
    double dTdx = (h_left / h / (h + h_left) * Tnext + (h - h_left) / h / h_left * T - h / h_left / (h + h_left) * Tprev);
    double slag_diff = 0;
    double slag_chem = 0.;
    /*cout << "Tcurr = " << T << "\n";
    cout << "Cp = " << Cp << "\n";
    cout << "Lambda =" << Lambda_All(chemkinReader, Yi, (Tnext + T) / 2.) << "\n";*/
    for (int k = 0; k < num_gas_species; k++) {
        slag_diff += rhoYkVk(chemkinReader, k, T, Yi, gradX, dTdx) * myget_Cpi(k, T)  * dTdx;
        slag_chem += ydot[k] * myget_Hi(k, T) * my_mol_weight(k);
        
        /*cout << "k = " << k << "\n";
        cout << "wdot = " << Ith(ydot, k + 1) << "\n";
        cout << "y = " << Yi[k] << "\n\n\n\n";*/
    }
    /*cout << "M  T = " << T << "\n";
    cout << "M  slag_chem = " << slag_chem << "\n";
    cout << "M  lambda = " << Lambda_All(chemkinReader, Yi, (Tnext + T) / 2.) << "\n";
    cout << "M Cp = " << Cp << "\n";
    for (int k = 0; k < num_gas_species; k++) {
        cout << "k name = " << name_species[k] << "\n";
        cout << "M d[X]/dt = " << Ith(ydot, k + 1) + wk_add[k] << "\n";
        cout << "h  = " << myget_Hi(k, T) << "\n";
    }
    cout << "M  slag_diff = " << slag_diff << "\n\n\n\n";*/

    /*cout << "\n\n\nM  slag_diff = " << slag_diff << "\n";*/
    if (flag == 1) {
        cout << "M  T = " << T << "\n";
        cout << "M  slag_chem = " << slag_chem << "\n";
        cout << "M  lambda = " << -(2. / (h + h_left)) *
            (Lambda_All(chemkinReader, Yi, (Tnext + T) / 2.) * (Tnext- T) / h
                - Lambda_All(chemkinReader, Yi, (T + Tprev) / 2.) * (T - Tprev) / h_left) << "\n";
        cout << "M  Cp = " << Cp * M * dTdx << "\n";
    }
    return -(2. / (h + h_left)) *
        (Lambda_All(chemkinReader, Yi, (Tnext + T) / 2.) * (Tnext - T) / h
            - Lambda_All(chemkinReader, Yi, (T + Tprev) / 2.) * (T - Tprev) / h_left)
        + Cp * M * dTdx + slag_chem + slag_diff;

}

double F_rightY(IO::ChemkinReader* chemkinReader, double* Yiprev, double* Yi, double* Yinext,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext, double* Xiprev, double* Xi, double* Xinext, double* gradX, double* Y_tmp,
    double M, const int k_spec, double* ydot, double* wk_add)
{

    double h_left = x - xprev;
    double h = xnext - x;
    double x_l = (x + xprev) / 2.;
    double x_r = (xnext + x) / 2.;
    double rho = get_rho(Yi, T);
    double rhoYkVk_r = 1.;
    double rhoYkVk_l = 1.;
    double slag_diff = 0.;
    double slag_chem = 0.;
    double dTdx = (h_left / h / (h + h_left) * Tnext + (h - h_left) / h / h_left * T - h / h_left / (h + h_left) * Tprev);
    Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);
    make_averageY(Y_tmp, Yi, Yinext);
    get_grad(gradX, Xi, Xinext, x, xnext);
    rhoYkVk_r = rhoYkVk(chemkinReader, k_spec, T / 2. + Tnext/ 2., Y_tmp, gradX, dTdx);


    make_averageY(Y_tmp, Yi, Yiprev);
    get_grad(gradX, Xiprev, Xi, xprev, x);
    rhoYkVk_l = rhoYkVk(chemkinReader, k_spec, T / 2. + Tprev / 2., Y_tmp, gradX, dTdx);
    //cout << "\n\n";

    slag_diff = (rhoYkVk_r - rhoYkVk_l) / (x_r - x_l);
    slag_chem = -my_mol_weight(k_spec) * ydot[k_spec];

    /*cout << "k = " << name_species[k_spec] << "\n";
    cout << "rho = const d[Xk]/dt = " << Ith(ydot, k_spec + 1) << "\n";
    cout << "drho/dt * Yk / Wk  = " << wk_add[k_spec] << "\n";
    cout << "MdY/dx = " << M * (h_left / h / (h + h_left) * Yinext[k_spec] + (h - h_left) / h / h_left * Yi[k_spec] - h / h_left / (h + h_left) * Yiprev[k_spec]) << "\n";
    cout << "Y slag_diff = " << slag_diff << "\n";
    cout << "Y slag_chem = " << slag_chem << "\n";
    cout << "Y rho = " << rho << "\n";
    cout << "\n\n";*/
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

int integrate_Y_IDA(IO::ChemkinReader* chemkinReader_temp, int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, double& M, int N_center, double* Y_leftb)
{
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
    int NEQ = 1 + NEQ_Y;
    vector<double> Yp_vect(N_x * 9);

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

    double* Yend = new double[num_gas_species];
    Yend[0] = my_mol_weight(0) * 0.124713
        ;
    Yend[1] = my_mol_weight(1) * 0.0314856
        ;
    Yend[2] = my_mol_weight(2) * 0.09486124
        ;
    Yend[3] = my_mol_weight(3) * 0.004672889
        ;
    Yend[4] = my_mol_weight(4) * 0.003471107
        ;
    Yend[5] = my_mol_weight(5) * 7.41308E-5
        ;
    Yend[6] = my_mol_weight(6) * 0.1357514
        ;
    Yend[7] = my_mol_weight(7) * 3.619971E-6
        ;
    Yend[8] = my_mol_weight(8) * 0.604967
        ;

    for (int i = 0; i < num_gas_species; i++) {
        Yend[i] /= 22.84775;
    }
    Get_molar_cons(data->Xi, Yend, 1206.609);
    for (int i = 0; i < num_gas_species; i++) {
        cout << data->Xi[i] << "\n";
    }
    //chem_vel(data->wk_add, M, data->Yi, data->Yinext, data->x[i], data->x[i + 1], T_vect[i], T_vect[i + 1], data->Xi, data->ydot);
    //chem_vel(1206.609, data->Xi, data->ydot);

    /* Allocate N-vectors. */
    yy = N_VNew_Serial(NEQ, ctx);
    if (check_retval((void*)yy, "N_VNew_Serial", 0)) return(1);

    yp = N_VClone(yy);
    if (check_retval((void*)yp, "N_VNew_Serial", 0)) return(1);
    avtol = N_VClone(yy);
    if (check_retval((void*)avtol, "N_VNew_Serial", 0)) return(1);
    cons = N_VClone(yy);
    if (check_retval((void*)cons, "N_VNew_Serial", 0)) return(1);
    
    Ith(cons, 1) = ONE;
    Ith(avtol, 1) = RCONST(1.0e-6);
    Ith(yy, 1) = M;
    Ith(yp, 1) = 0;

    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rtol = RCONST(1.0e-8);
    atval = N_VGetArrayPointer(avtol);
    consval = N_VGetArrayPointer(cons);

    for (int i = 1; i < NEQ; i++) {
        consval[i] = 1.0;   /*constraint*/
        yval[i] = Y_vect[num_gas_species + i - 1];
        atval[i] = RCONST(1.0e-8);
        cout << "i = " << i << "   = " << yval[i] << endl;
        if ((i) % num_gas_species == 0) cout << endl;
    }
    /* Create and initialize  y, y', and absolute tolerance vectors. */
    t0 = ZERO;
    data->M = Ith(yy, 1);
    double rho;
    flag = 1;

    for (int i = 1; i < data->Nx - 1; i++) {
        MakeYvectors(data->Yiprev, data->Yi, data->Yinext, yval, data->Y_left_bound, N_x, i);
        Get_molar_cons(data->Xi, data->Yi, T_vect[i]);
        chem_vel(data->wk_add, M, data->Yi, data->Yinext, data->x[i], data->x[i + 1], T_vect[i], T_vect[i + 1], data->Xi, data->ydot);

        cout << "\n\n\n\ni = " << i << "\n";
        cout << "M = " << data->M << "\n";
        cout << "T = " << T_vect[i] << "\n";
        rho = get_rho(data->Yi, T_vect[i]);
        for (int k_spec = 1; k_spec <= num_gas_species; k_spec++) {
           ypval[k_spec + (i - 1) * num_gas_species] = -F_rightY(data->chemkinReader, data->Yiprev, data->Yi, data->Yinext,
               data->T[i - 1], data->T[i], data->T[i + 1],
               data->x[i - 1], data->x[i], data->x[i + 1],
               data->Xiprev, data->Xi, data->Xinext, data->gradX, data->Y_tmp,
               data->M, k_spec - 1, data->ydot, data->wk_add) / rho;
           cout << "yval = " << yval[k_spec + (i - 1) * num_gas_species] << "\n";
           cout << "ypval = " << ypval[k_spec + (i - 1) * num_gas_species] << "\n\n\n\n\n";
        }

    }
    //cout << NEQ << " = " << Ith(res_vect, NEQ + 1) << endl;

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
    retval = IDASetMaxNumSteps(mem, 2000);
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
    double tout1 = pow(10, -7);
    tout = tout1;
    int iend = 500;
    int number = 100;
    ofstream fout;
    double Y_H2, Y_O2;
    double W, w_dot;
    double sum_Y = 0;
    ofstream f_ida;
    for (int i = 0; i < num_gas_species * N_x; i++)
        Yp_vect[i] = 0;

    while (iout < iend) {
        retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);
        cout << "tout = " << tout << "\n";
        for (int i = 1; i < NEQ; i++) {
            Y_vect[num_gas_species + i - 1] = yval[i];
            Yp_vect[num_gas_species + i - 1] = ypval[i];
        }

        for (int i = NEQ; i < NEQ + num_gas_species; i++) {
            Y_vect[num_gas_species + i - 1] = Y_vect[i - 1];
        }
        cout << "M = " << yval[0] << "\n";
        
        Write_to_file2("detai" + to_string(tout * pow(10, 6)), f_ida, x_vect,
            T_vect, Y_vect, M, N_x, 2);
        Write_to_file2("proizvod" + to_string(tout * pow(10, 6)), f_ida, x_vect,
            T_vect, Yp_vect, M, N_x, 2);
        if (check_retval(&retval, "IDASolve", 1)) return(1);

        iout++;
        tout += tout1;
    }

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

static int func_Y_IDA(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
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
    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rval = N_VGetArrayPointer(rr);
    int Ncentr = data->N_centr;

    MakeYvectors(Yiprev, Yi, Yinext, yval, data->Y_left_bound, myNx, Ncentr);
    Get_molar_cons(data->Xi, data->Yi, T_vect[Ncentr]);
    Get_molar_cons(data->Xiprev, data->Yiprev, T_vect[Ncentr - 1]);
    Get_molar_cons(data->Xinext, data->Yinext, T_vect[Ncentr + 1]);

    //for (int j = 1; j <= num_gas_species; j++) {
    //    cout << "in IDA Yiprev[j] = " << Yiprev[j - 1] << "\n";
    //    cout << "in IDA Yi[j] = " << Yi[j - 1] << "\n";
    //    cout << "in IDA Yinext[j] = " << Yinext[j - 1] << "\n\n";
    //}

    chem_vel(T_vect[Ncentr], data->Xi, data->ydot);
    double rho = get_rho(Yi, T_vect[Ncentr]);
    data->M = yval[0];
    rval[0] = F_right(data->chemkinReader, Yi, Yinext,
        data->T[Ncentr - 1], data->T[Ncentr], data->T[Ncentr + 1],
        x_cells[Ncentr - 1], x_cells[Ncentr], x_cells[Ncentr + 1],
        data->Xiprev, data->Xi, data->Xinext, data->gradX, data->Y_tmp,
        data->M, data->ydot, data->wk_add);
    cout << "fdata[0]" << rval[0] << "\n";
    cout <<  "rval low M = ";
    cout << fixed << setprecision(15) << yval[0] << "\n";

   for (int i = 1; i < myNx - 1; i++) {
        //cout << "\n\n\n\nin func i = " << i << "\n";
        MakeYvectors(Yiprev, Yi, Yinext, yval, data->Y_left_bound, myNx, i);
        Get_molar_cons(data->Xi, Yi, T_vect[i]);

        chem_vel(data->wk_add, M, Yi, Yinext, x_cells[i], x_cells[i + 1], T_vect[i], T_vect[i + 1], data->Xi, data->ydot);
        flag = 0;
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            if (data->ydot[k_spec] != data->ydot[k_spec]){
                flag = 1;
            }
        }
        rho = get_rho(Yi, T_vect[i]);
        /*if (flag == 1) {
            cout << "fdata[0]" << rval[0] << "\n";
            cout << "\n\n\n\nin func i = " << i << "\n";
            cout << "M = " << data->M << "\n";
        }*/
        
        for (int k_spec = 1; k_spec <= num_gas_species; k_spec++) {
           /* if (flag == 1) {
                cout << "T = " << T_vect[i] << "\n";
                cout << "k_spec = " << name_species[k_spec - 1] << "\n";
                cout << "Y = " << yval[k_spec + (i - 1) * num_gas_species] << "\n";
                cout << "dY/dt = " << ypval[k_spec + (i - 1) * num_gas_species] << "\n";
            }*/
            rval[k_spec + (i - 1) * num_gas_species] = ypval[k_spec + (i - 1) * num_gas_species] + F_rightY(data->chemkinReader, data->Yiprev, data->Yi, data->Yinext,
                data->T[i - 1], data->T[i], data->T[i + 1],
                data->x[i - 1], data->x[i], data->x[i + 1],
                data->Xiprev, data->Xi, data->Xinext, data->gradX, data->Y_tmp,
                data->M, k_spec - 1, data->ydot, data->wk_add) / rho;
            //cout << "yval = " << yval[k_spec + (i - 1) * num_gas_species] << "\n";
            //cout << "ypval = " << ypval[k_spec + (i - 1) * num_gas_species] << "\n\n\n\n\n";
        }
        if (flag == 1) {
            flag = 0;
        }
        //cout << "yval = " << yval[num_gas_species] << "\n";
    }
    return(0);
}

