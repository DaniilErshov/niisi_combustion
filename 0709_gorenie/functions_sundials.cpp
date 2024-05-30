#include "functions_sundials.h"

int flag = 0;
#define RTOL  RCONST(1.0e-5)
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(pow(10, -6))     /* first output time      */
#define TMULT RCONST(1.006)     /* output time factor     */
#define NOUT  1800
#define ZERO  RCONST(0.0)
int InitialData(int& Nx, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect, 
    vector<double>& u_vect, double& M, double Tstart, double Tfinish, double* Ystart, double* Yend)
{
    double h = l / (Nx - 1);
    double x_start = 0.5;
    //double x_start = 0.1;
    //double x_finish = 0.2;

    double j = 0;
    //std::cout << "M = " << M << "\n";

    for (int i = 0; i < Nx; i++) {
        x_vect[i] = h * i;
    }

    for (int i = 0; i < Nx; i++) {

        if (x_vect[i] <= x_start)
        {
            T_vect[i] = Tstart;
        }
        else {
            //T_vect[i] = (Tfinish - Tstart) / (x_finish - x_start) * (x_vect[i] - x_start) + Tstart;
            T_vect[i] = Tfinish;
            j++;
        }
        u_vect[i] = 0;
    }

    double O2_in = 0.21;
    double N2_in = 0.79;

    double koeff_O2 = 1.;

    double koeff_N2 = N2_in / O2_in * koeff_O2;

    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        Xi[k_spec] = 0;
    }

    Xi[komponents["O2"]] = koeff_O2 / (koeff_O2 + koeff_N2);
    Xi[komponents["N2"]] = koeff_N2 / (koeff_O2 + koeff_N2);
    moleFraction_to_massFraction(Xi, Yi);

    for (int i = 0; i < Nx; i++)
    {
        if (x_vect[i] < x_start) {
            for (int k = 0; k < num_gas_species; k++) {
                Y_vect[k + i * num_gas_species] = 0;
                //Y_vect[k + i * num_gas_species] = Ystart[k] + (Yend[k] - Ystart[k]) / (x_finish - x_start) * (x_vect[i] - x_start);
            }
            Y_vect[komponents["NC7H16"] + i * num_gas_species] = 1;
        }
        else {
            for (int k = 0; k < num_gas_species; k++) {
                Y_vect[k + i * num_gas_species] = Yi[k];
                //Y_vect[k + i * num_gas_species] = Ystart[k] + (Yend[k] - Ystart[k]) / (x_finish - x_start) * (x_vect[i] - x_start);
            }
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

void Add_elem_simple(vector<double>& T, vector<double>& Y, vector<double>& x, int& N_x, int& N_center, double b, int number, int number_start, double& T_center)
{
    double T_max = 0, T_min = T[0];
    int Nx_add = 0;
    int j_t = 1;
    vector<double> Ymax;
    vector<double> Ymin;
    Ymax.resize(num_gas_species);
    Ymin.resize(num_gas_species);
    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        Ymax[k_spec] = 0;
        Ymin[k_spec] = 1.0;
    }
    for (int i = 0; i < N_x; i++)
    {
        if (T[i] > T_max) T_max = T[i];
        if (T[i] < T_min) T_min = T[i];
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            if (Y[k_spec + num_gas_species * i] > Ymax[k_spec]) Ymax[k_spec] = Y[k_spec + num_gas_species * i];
            if (Y[k_spec + num_gas_species * i] < Ymin[k_spec]) Ymin[k_spec] = Y[k_spec + num_gas_species * i];
        }
    }
    bool need_add_cell = 0;
    while (j_t < N_x - 2)
    {
        need_add_cell = 0;
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            if (fabs(Y[k_spec + num_gas_species * j_t] - Y[k_spec + num_gas_species * (j_t - 1)]) > b * (Ymax[k_spec] - Ymin[k_spec]))
            {
                need_add_cell = 1;
            }
        }
        if (fabs(T[j_t] - T[j_t - 1]) > b * (T_max - T_min))
        {
            need_add_cell = 1;
        }
        if (need_add_cell) {
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
        x.push_back(x[N_x - 1] + 1.5 * (x[N_x - 1] - x[N_x - 2]));
        for (int i = 0; i < num_gas_species; i++) {
            Yi[i] = Y[Y.size() - num_gas_species + i];
        }
        for (int i = 0; i < num_gas_species; i++) {
            Y.push_back(Yi[i]);
        }
        N_x++;
    }

    for (int k = 0; k < number_start; k++) {
        for (int i = 0; i < num_gas_species; i++) {
            Yi[i] = Y[i];
        }
        for (int i = num_gas_species - 1; i >= 0; i--) {
            Y.insert(Y.begin(), Yi[i]);
        }
        T.insert(T.begin(), T[0]);
        x.insert(x.begin(), x[0] - 1.5 * (x[1] - x[0]));

    }
    N_x = x.size();
    resize_koeff_vectors(N_x);
    //T_center = (Tstart + Tfinish) / 2.;
    T_min = 1000;
    for (int i = 0; i < T.size(); i++) {
        if (abs(T[i] - T_center) < T_min) {
            N_center = i;
            T_min = abs(T_vect[i] - T_center);
        }
    }
    T_center = T[N_center];
    cout << "Tcenter = " << T[N_center] << "\n";
    cout << "Ncenter = " << N_center << "\n";
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

double F_right(UserData data,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext, 
    double uprev, double u, double unext, int number_cell)
{
    double h_left = x - xprev;
    double h = xnext - x;
    double Cp = get_Cp(num_gas_species, Yi, T, number_cell);
    double rho = get_rho(Yi, T);
    double W = get_W(Yi);
    get_grad(gradX, Xi, Xinext, x, xnext);
    double dTdx = (h_left / h / (h + h_left) * Tnext + (h - h_left) / h / h_left * T - h / h_left / (h + h_left) * Tprev);
    double slag_diff = 0;
    double slag_chem = 0.;

    double Vc = 0;
    double YkVk_slag = 0;
    set_Dij_res(T, number_cell, 'c');
    for (int k = 0; k < num_gas_species; k++) {
        YkVk[k] = YkVk_func(k, T, Yi, gradX, Xi, Yi);
        Vc -= YkVk[k];
    }
    for (int i = 0; i < 9; i++) {
        Cpn[i] = 0;
        Hn[i] = 0;
    }

    for (int k = 0; k < num_gas_species; k++) {
        if (T > chec.chemkinReader->species()[k].thermo().getTCommon())
            for (int i = 0; i < 9; i++) {
                Hn[i] += ydot[k] * phyc.Cp_coef_hT[k][i] * my_mol_weight(k);
                Cpn[i] += rho * (YkVk[k] + Yi[k] * Vc) * phyc.Cp_coef_hT[k][i] * dTdx;
            }
        else {
            for (int i = 0; i < 9; i++) {
                Hn[i] += ydot[k] * phyc.Cp_coef_lT[k][i] * my_mol_weight(k);
                Cpn[i] += rho * (YkVk[k] + Yi[k] * Vc) * phyc.Cp_coef_lT[k][i] * dTdx;
            }
        }
        //slag_diff += rho * data->YkVk[k] * myget_Cpi(k, T)  * dTdx;
        //slag_chem += data->ydot[k] * myget_Hi(k, T) * my_mol_weight(k);
    }

    //slag_diff = get_dCpi(Cpn, T);
    slag_diff = 0;
    slag_chem = get_dHiRT(Hn, T);

    make_averageY(X_tmp, Xi, Xinext);
    double lambda_right = Lambda_All(X_tmp, (Tnext + T) / 2., number_cell, 'r');
    make_averageY(X_tmp, Xi, Xiprev);
    double lambda_left = Lambda_All(X_tmp, (T + Tprev) / 2., number_cell, 'l');
    double M = get_rho(Yi, T) * u;
    return (2. / (h + h_left)) / x / x *
        (lambda_right * (Tnext - T) / h * pow((xnext + x) / 2., 2) 
            - lambda_left * (T - Tprev) / h_left * pow((x + xprev) / 2., 2))
        - Cp * M * dTdx + slag_chem + slag_diff;

}

double F_rightY(UserData data, int k_spec,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext, 
    double uprev, double u, double unext, int i)
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

    rhoYkVk_r = data->rho_r * (YkVk_r[k_spec] + Y_tmp_r[k_spec] * data->Vc_r);
    //cout << "data->YkVk[k_spec] =  " << data->YkVk[k_spec] << "\n";


    rhoYkVk_l = data->rho_l * (YkVk_l[k_spec] + Y_tmp_l[k_spec] * data->Vc_l);
    //cout << "\n\n";

    slag_diff = (rhoYkVk_r * pow(x_r, 2) - rhoYkVk_l * pow(x_l, 2)) / (x_r - x_l) / x / x;
    slag_chem = my_mol_weight(k_spec) * ydot[k_spec];
    //cout << "slag_chem =  " << slag_chem << "\n";
    //cout << "slag_diff =  " << slag_diff << "\n";
    double M = get_rho(Yi, T) * u;
    return -M * (h_left / h / (h + h_left) * Yinext[k_spec] + (h - h_left) / h / h_left * Yi[k_spec] - h / h_left / (h + h_left) * Yiprev[k_spec])
        + slag_chem + slag_diff;
}

void Write_to_file(string str, ofstream& fout, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, vector<double>& u_vect, double M, int N_x, int number) {
    double x_start, x_finish, D;
    double rho;
    fout.open(str + ".dat");
    string title2 = R"(VARIABLES= "x", "T", "u",)";
    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        title2.append(R"(, "Y_)" + komponents_str[k_spec] + R"(")");
    }
    fout << "TITLE=\"" << "Graphics" << "\"" << endl;
    fout << title2 << endl;
    for (int i = 0; i < N_x; i++) {
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            Yi[k_spec] = Y_vect[k_spec + i * num_gas_species];
        }
        rho = get_rho(Yi, T_vect[i]);

        if (number == 1) {
            fout << x_vect[i] << "  " << T_vect[i] << " " << u_vect[i];
            for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
                fout << " " << Y_vect[k_spec + i * num_gas_species];
            }
            fout << endl;
        }
    }
    fout.close();
}

void Init_Data(UserData data, int N_x, vector<double>& x_vect,
    vector<double>& T_vect, int NEQ, int NEQ_Y,
    int N_center, double* Y_leftb) {
    data->Nx = N_x;

    data->NEQ = NEQ;
    data->NEQ_Y = NEQ_Y;
    data->Tl = T_vect[0];
    data->T_center = T_vect[N_center];
    std::cout << "T_center = " << data->T_center << "\n";
    data->N_centr = N_center;
    data->N_m = 1 - 1;

    for (int i = 0; i < num_gas_species; i++) {
        ydot[i] = 0;
        Y_left_bound[i] = Y_leftb[i];
    }
}

void find_diff_slag(UserData data, double Tcurr, double Tnext, double* Yi, double* Yinext,
    double* Xi, double* Xinext, double* Ykvk_side, double* Y_tmp_side, double* X_tmp_side, double* gradX_side, double& rho_side, double& Vc_side, int i, char side) {
    make_averageY(Y_tmp_side, Yi, Yinext);
    make_averageY(X_tmp_side, Xi, Xinext);
    get_grad(gradX_side, Xi, Xinext, x_vect[i], x_vect[i + 1]);
    rho_side = get_rho(Y_tmp_side, (Tcurr + Tnext) / 2.);
    Vc_side = 0;
    set_Dij_res((Tcurr + Tnext) / 2., i, side);
    for (int k = 0; k < num_gas_species; k++) {
        //cout << "k = " << komponents_str[k] << "\n";
        int jop;
        Ykvk_side[k] = YkVk_func(k, (Tcurr + Tnext) / 2., Yi, gradX_side, X_tmp_side, Y_tmp_side);
        Vc_side -= Ykvk_side[k];
        //cout << "k = " << komponents_str[k] << "  Ykvk = " << Ykvk_side[k] << "\n";
    }

}

void MakeYvectors_kins(UserData data,
    double* Y, int myNx, int i, double Tl) {
    //cout << "i = " << i << "\n";
    for (int j = 0; j < num_gas_species; j++) {

        Yi[j] = Y[j + (i - 1) * (num_gas_species + 1) + i - 1];
        if (i == 1) {
            Yiprev[j] = 0;
            Yiprev[komponents["NC7H16"]] = 1;
        }

        else Yiprev[j] = Y[j + (i - 2) * (num_gas_species + 1) + i - 2];

        if (i == myNx - 2)  Yinext[j] = Yi[j];
        else  Yinext[j] = Y[j + (i) * (num_gas_species + 1) + i];
    }
}

void makeYstart(double koeff_topl, string fuel, double O2_in, double N2_in, double* Ystart) {

    /* Ystart[0] = koeff_topl * 2.0 / (koeff_topl * 2. + 1. + 3.76);
     Ystart[2] = 1.0 / (koeff_topl * 2. + 1. + 3.76);
     Ystart[8] = 3.76 / (koeff_topl * 2. + 1. + 3.76);
     double Wsmes = 0;

     for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
         Wsmes += Ystart[k_spec] * phyc.mol_weight[k_spec];
     }

     for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
         Ystart[k_spec] *= phyc.mol_weight[k_spec] / Wsmes;
         cout << "Yi = " << Ystart[k_spec] << "\n";
     }*/
    double n = 7;
    double m = 16;
    double koeff_O2 = (n + 0.25 * m);

    double koeff_N2 = N2_in / O2_in * koeff_O2;
    double koeff_fuel = koeff_topl;
    cout << "koeff_fuel = " << koeff_fuel << "\n";
    for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        Xi[k_spec] = 0;
    }

    Xi[komponents["O2"]] = koeff_O2 / (koeff_O2 + koeff_N2 + koeff_fuel);
    Xi[komponents["N2"]] = koeff_N2 / (koeff_O2 + koeff_N2 + koeff_fuel);
    Xi[komponents[fuel]] = koeff_fuel / (koeff_O2 + koeff_N2 + koeff_fuel);

    cout << "mol O2 = " << Xi[komponents["O2"]] << "\n";
    cout << "mol N2 = " << Xi[komponents["N2"]] << "\n";
    cout << "mol " << fuel << " = " << Xi[komponents[fuel]] << "\n";

    moleFraction_to_massFraction(Xi, Ystart);


    cout << "O2 = " << Ystart[komponents["O2"]] << "\n";
    cout << "N2 = " << Ystart[komponents["N2"]] << "\n";
    cout << fuel << " = " << Ystart[komponents[fuel]] << "\n";
}

int integrate_All_IDA_M(int N_x, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, 
    vector<double>& u_vect, double& M, int N_center, double* Y_leftb, int iter, double t_fix) {

    void* mem;
    N_Vector yy, yp, avtol, cons, id;
    realtype rtol, * yval, * ypval, * atval, * consval, * id_val;
    realtype t0, tout, tret;
    int iout, retval, retvalr;
    SUNMatrix A;
    SUNLinearSolver LS;
    SUNNonlinearSolver NLS;
    SUNContext ctx;
    UserData data;
    data = (UserData)malloc(sizeof * data);
    int NEQ_Y = num_gas_species * (N_x - 2);
    int NEQ = (N_x - 2) * ((num_gas_species + 2));
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


    data->my_tcur = 0;
    data->my_numjac = -1;
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
    id = N_VClone(yy);
 

    yval = N_VGetArrayPointer(yy);
    id_val = N_VGetArrayPointer(id);
    ypval = N_VGetArrayPointer(yp);
    rtol = RCONST(1.0e-3);
    atval = N_VGetArrayPointer(avtol);
    consval = N_VGetArrayPointer(cons);
    int number_spec = 0;

    data->M = M;
    std::cout << "M_in_IDA = " << M << "\n";

    int i_temp = 0;
    for (int i = 1; i < N_x - 1; i++) {
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            consval[i_temp] = 0.0;   /*constraint*/
            yval[i_temp] = Y_vect[i * num_gas_species + k_spec];
            atval[i_temp] = pow(10, -5);
            id_val[i_temp] = 1.;
            //cout << "Yi = " << i_temp << "   = " << yval[i_temp] << endl;
            i_temp++;
        }
        yval[i_temp] = T_vect[i];
        consval[i_temp] = 1.0;
        atval[i_temp] = pow(10, 0);
        id_val[i_temp] = 1.;
        //cout << "Ti = " << i_temp << "   = " << yval[i_temp] << endl << endl;
        i_temp++;
        yval[i_temp] = u_vect[i];
        atval[i_temp] = pow(10, -3);
        consval[i_temp] = 0.0;
        id_val[i_temp] = 0;
        //cout << "Mi = " << i_temp << "   = " << yval[i_temp] << endl << endl;
        i_temp++;

        //if ((i + 1) % num_gas_species == 0) cout << endl;
    }

    /* Create and initialize  y, y', and absolute tolerance vectors. */
    t0 = ZERO;

    double rho;
    ofstream f_ida;
    double sumY = 0;

    flag = 1;

    cout << "Nx_ALL = " << data->Nx << "\n";

    i_temp = 0;
    int i_T, i_u;
    double T_prev, T_curr, T_next;
    double u_prev, u_curr, u_next;
    for (int i = 1; i < data->Nx - 1; i++) {

        i_T = i * (num_gas_species + 2) - 2;
        i_u = i * (num_gas_species + 2) - 1;
        //cout << "in func i = " << i << "\n";
        //cout << "M = " << data->M << "\n\n";
        if (i == N_x - 2) {
            T_curr = yval[i_T];
            T_prev = yval[i_T - (num_gas_species + 2)];
            T_next = T_curr;

            u_curr = yval[i_u];
            u_prev = yval[i_u - (num_gas_species + 2)];
            u_next = u_curr;

        }
        else if (i == 1) {
            T_prev = data->Tl;
            T_curr = yval[i_T];
            T_next = yval[i_T + (num_gas_species + 2)];

            u_curr = yval[i_u];
            u_prev = u_curr;
            u_next = yval[i_u + (num_gas_species + 2)];
        }
        else {
            T_curr = yval[i_T];
            T_prev = yval[i_T - (num_gas_species + 2)];
            T_next = yval[i_T + (num_gas_species + 2)];

            u_curr = yval[i_u];
            u_prev = yval[i_u - (num_gas_species + 2)];
            u_next = yval[i_u + (num_gas_species + 2)];
        };

        //cout << "T_prev" << T_prev << "\n";
        //cout << "T_curr" << T_curr << "\n";
        //cout << "T_next" << T_next << "\n\n";
        MakeYvectors_kins(data, yval, N_x, i, data->Tl);
        //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        //    cout << "Yi = " << Yi[k_spec] << "\n";
        //}
        Get_molar_cons(Xi, Yi, T_curr);

        chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr,
            T_curr, Xi, ydot, i);


        Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);

        find_diff_slag(data, T_curr, T_next, Yi, Yinext, Xi, Xinext, YkVk_r, Y_tmp_r, X_tmp_r, gradX_r, data->rho_r, data->Vc_r, i, 'r');

        find_diff_slag(data, T_prev, T_curr, Yiprev, Yi, Xiprev, Xi, YkVk_l, Y_tmp_l, X_tmp_l, gradX_l, data->rho_l, data->Vc_l, i - 1, 'l');

        rho = get_rho(Yi, T_curr);
        sumY = 0;
        double dWdt = 0;
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            int jop;
            ypval[i_temp] = F_rightY(data, k_spec,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1], 
                u_prev, u_curr, u_next, i) / rho;
            //cout << "i_temp initial  = " << i_temp << " Y  = " << ypval[i_temp] << "\n";
            dWdt += ypval[i_temp] / phyc.mol_weight[k_spec];
            i_temp++;
        }
        dWdt *= -pow(get_W(Yi), 2);
        double Cp = get_Cp(num_gas_species, Yi, T_curr, i);

        //cout << "i_temp = " << i_temp << " T  = " << yval[i_temp] << "\n\n";
        ypval[i_temp] = F_right(data,
            T_prev, T_curr, T_next,
            x_vect[i - 1], x_vect[i], x_vect[i + 1], 
            u_prev, u_curr, u_next, i) / rho / Cp;
        double dTdt = ypval[i_temp];
        //cout << "i_temp initial = " << i_temp << " = " << ypval[i_temp] << "\n";
        i_temp++;
        double drhodt = dWdt * P / R / T_curr - P * get_W(Yi) / R / pow(T_curr, 2) * dTdt;

        ypval[i_temp] = 0;
        cout << "ind = " << i_temp << " last eq = " << drhodt - F_right_rho(data,
            T_prev, T_curr, T_next,
            x_vect[i - 1], x_vect[i], x_vect[i + 1],
            u_prev, u_curr, u_next, i) << "\n\n";

        i_temp++;
    }

    /* Call IDACreate and IDAInit to initialize IDA memory */
    mem = IDACreate(ctx);
    if (check_retval((void*)mem, "IDACreate", 0)) return(1);

    retval = IDAInit(mem, func_All_IDA_M, t0, yy, yp);
    if (check_retval(&retval, "IDAInit", 1)) return(1);
    /* Call IDASVtolerances to set tolerances */
    data->sun_mem = mem;
    retval = IDASVtolerances(mem, rtol, avtol);
    if (check_retval(&retval, "IDASVtolerances", 1)) return(1);
    retval = IDASetConstraints(mem, cons);
    retval = IDASetUserData(mem, data);
    if (check_retval(&retval, "IDASetUserData", 1)) return(1);
    //retval = IDASetMaxNumSteps(mem, 1000);
    /* Create dense SUNMatrix for use in linear solves */
    int mu = 2 * (num_gas_species + 2);
    int ml = 2 * (num_gas_species + 2);

    A = SUNBandMatrix(NEQ, mu, ml, ctx);
    if (check_retval((void*)A, "SUNBandMatrix", 0)) return(1);

    /* Create banded SUNLinearSolver object */
    LS = SUNLinSol_Band(yy, A, ctx);
    if (check_retval((void*)LS, "SUNLinSol_Band", 0)) return(1);

    /* Attach the matrix and linear solver */
    retval = IDASetLinearSolver(mem, LS, A);
    if (check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

    NLS = SUNNonlinSol_Newton(yy, ctx);
    if (check_retval((void*)NLS, "SUNNonlinSol_Newton", 0)) return(1);

    /* Attach the nonlinear solver */
    retval = IDASetNonlinearSolver(mem, NLS);
    if (check_retval(&retval, "IDASetNonlinearSolver", 1)) return(1);
    //IDASetInitStep(mem, pow(10, -7));

    //IDASetMaxOrd(mem, 1);
    //IDASetInitStep(mem, pow(10, -7));
    tout = RCONST(0.00001);
    retval = IDASetId(mem, id);
    if (check_retval(&retval, "IDASetId", 1)) return(1);

    retval = IDACalcIC(mem, 1, tout);
    if (check_retval(&retval, "IDACalcIC", 1)) return(1);
    /* In loop, call IDASolve, print results, and test for error.
        Break out of loop when NOUT preset output times have been reached. */

    iout = 0;
    double tout1 = t_fix;
    tout = tout1;
    int iend = 2000000;
    int number = 1;
    ofstream fout;
    ofstream foutw;
    double Y_H2, Y_O2;
    double W, w_dot;
    double sum_Y = 0;
    double alpha = 1;


    data->N_m = 0;
    while (iout < ida_steps) {
        retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);

        cout << "t = " << tout << "\n";
        if (check_retval(&retval, "IDASolve", 1)) return(1);
        iout++;
        tout1 *= 1.0;
        tout += tout1;
    }
    i_temp = 0;
    for (int i = 1; i < N_x - 1; i++) {
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            Y_vect[i * num_gas_species + k_spec] = yval[i_temp] - eps;
            i_temp++;
        }
        T_vect[i] = yval[i_temp];
        i_temp++;
        u_vect[i] = yval[i_temp];
        //cout << "i = " << i << "M = " << yval[i_temp] << "\n";
        i_temp++;
    }
    for (int i = 0; i < num_gas_species; i++) {
        Y_vect[(num_gas_species) * (N_x - 1) + i] = Y_vect[(num_gas_species) * (N_x - 2) + i];
    }

    T_vect[N_x - 1] = T_vect[N_x - 2];

    std::cout << "tout = " << tout << "\n";
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

static int func_All_IDA_M(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
    realtype* yval, * ypval, * rval;
    UserData data;
    double Temp;
    double M;
    double sumY = 0;
    int i_T, i_u;
    data = (UserData)user_data;
    double Cp;
    double T_curr, T_prev, T_next;
    int j;
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
    int i_temp = 0;
    double u_curr, u_next, u_prev;

    for (int i = 1; i < myNx - 1; i++) {
        //data->M = yval[i * (num_gas_species + 2) - 1];
        i_T = i * (num_gas_species + 2) - 2;
        i_u = i * (num_gas_species + 2) - 1;
        //cout << "in func i = " << i << "\n";
        //cout << " T_curr = " << yval[i_T] << "\n\n";
        bool nan_flag = 0;
        if (i == myNx - 2) {
            T_curr = yval[i_T];
            T_prev = yval[i_T - (num_gas_species + 2)];
            T_next = T_curr;

            u_curr = yval[i_u];
            u_prev = yval[i_u - (num_gas_species + 2)];
            u_next = u_curr;

        }
        else if (i == 1) {
            T_curr = yval[i_T];
            T_prev = T_curr;
            T_next = yval[i_T + (num_gas_species + 2)];

            u_curr = yval[i_u];
            u_prev = 0;
            u_next = yval[i_u + (num_gas_species + 2)];
        }
        else {
            T_curr = yval[i_T];
            T_prev = yval[i_T - (num_gas_species + 2)];
            T_next = yval[i_T + (num_gas_species + 2)];

            u_curr = yval[i_u];
            u_prev = yval[i_u - (num_gas_species + 2)];
            u_next = yval[i_u + (num_gas_species + 2)];
        };
        //cout << "T_prev" << T_prev << "\n";
        //cout << "T_curr" << T_curr << "\n";
        //cout << "T_next" << T_next << "\n\n";
        MakeYvectors_kins(data, yval, myNx, i, data->Tl);
        Get_molar_cons(Xi, Yi, T_curr);

        chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr,
            T_curr, Xi, ydot, i);


        Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);

        find_diff_slag(data, T_curr, T_next, Yi, Yinext, Xi, Xinext, YkVk_r, Y_tmp_r, X_tmp_r, gradX_r, data->rho_r, data->Vc_r, i, 'r');

        find_diff_slag(data, T_prev, T_curr, Yiprev, Yi, Xiprev, Xi, YkVk_l, Y_tmp_l, X_tmp_l, gradX_l, data->rho_l, data->Vc_l, i - 1, 'l');

        rho = get_rho(Yi, T_curr);
        sumY = 0;
        double dWdt = 0;
        int i_temp_start = i_temp;
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            rval[i_temp] = rho * ypval[i_temp] - F_rightY(data, k_spec,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1], 
                u_prev, u_curr, u_next, i);

           /* cout << "rval Y i_temp = " << i_temp << " = " << rval[i_temp] << "\n";
            cout << "yval Y i_temp = " << i_temp << " = " << yval[i_temp] << "\n";
            cout << "ypval Y i_temp = " << i_temp << " = " << ypval[i_temp] << "\n\n\n\n\n";*/

            if (rval[i_temp] != rval[i_temp]) nan_flag = 1;
            dWdt += ypval[i_temp] / phyc.mol_weight[k_spec];
            i_temp++;
        }
        dWdt *= -pow(get_W(Yi), 2);

        Cp = get_Cp(num_gas_species, Yi, T_curr, i);

        double dTdt = ypval[i_temp];
        rval[i_temp] = rho * Cp * ypval[i_temp] - F_right(data,
            T_prev, T_curr, T_next,
            x_vect[i - 1], x_vect[i], x_vect[i + 1], 
            u_prev, u_curr, u_next, i);

       /* cout << "rval T i_temp = " << i_temp << " = " << rval[i_temp] << "\n";
        cout << "yval T i_temp = " << i_temp << " = " << yval[i_temp] << "\n";
        cout << "ypval T i_temp = " << i_temp << " = " << ypval[i_temp] << "\n\n\n\n\n";*/

        if (rval[i_temp] != rval[i_temp]) nan_flag = 1;
        i_temp++;
        double drhodt = dWdt * P / R / T_curr - P * get_W(Yi) / R / pow(T_curr, 2) * dTdt;

        rval[i_temp] = drhodt - F_right_rho(data,
            T_prev, T_curr, T_next,
            x_vect[i - 1], x_vect[i], x_vect[i + 1], 
            u_prev, u_curr, u_next, i);

       /* cout << "rval u i_temp = " << i_temp << " = " << rval[i_temp] << "\n";
        cout << "yval u i_temp = " << i_temp << " = " << yval[i_temp] << "\n";
        cout << "ypval u i_temp = " << i_temp << " = " << ypval[i_temp] << "\n\n\n\n\n";*/

        int i_temp_finish = i_temp;
        if (rval[i_temp] != rval[i_temp]) nan_flag = 1;
        if (nan_flag) {
            cout << "i = " << i << "\n";
            for (int i_temp2 = i_temp_start; i_temp2 <= i_temp_finish; i_temp2++) {
                cout << "rval i_temp = " << i_temp << " = " << rval[i_temp] << "\n";
                cout << "yval i_temp = " << i_temp << " = " << yval[i_temp] << "\n";
                cout << "ypval i_temp = " << i_temp << " = " << ypval[i_temp] << "\n\n\n\n\n";
            }
        }

    }
    return 0;
}

double F_right_rho(UserData data,
    double Tprev, double T, double Tnext, double xprev, double x, double xnext, 
    double uprev, double u, double unext, 
    int number_cell){
    double h_left = x - xprev;
    double h = xnext - x;

    double rho_r = get_rho(Yinext, Tnext);
    double rho = get_rho(Yi, T);
    double rho_l = get_rho(Yiprev, Tprev);

    double d_r2_rho_u_dr = (h_left / h / (h + h_left) * xnext * xnext * rho_r * unext + 
        (h - h_left) / h / h_left * x * x * rho * u - 
        h / h_left / (h + h_left) * xprev * xprev * rho_l * uprev);

    return -1. / x / x * d_r2_rho_u_dr;

}

void set_Dij_res(double T, int number_cell, char side) {
    if (!flag_use_save_koeffs) {
        for (unsigned short int i = 0; i < num_gas_species; ++i) {
            for (unsigned short int j = i; j < num_gas_species; ++j) {
                Dij_res[i][j] = Dij_func5(i, j, T, number_cell);
                Dij_res[j][i] = Dij_res[i][j];
            }
        }
    }
    else {
        if (side == 'r') {
            for (unsigned short int i = 0; i < num_gas_species; ++i) {
                for (unsigned short int j = i; j < num_gas_species; ++j) {
                    Dij_res[i][j] = Dij_arr_r[number_cell][i][j];
                    Dij_res[j][i] = Dij_res[i][j];
                }
            }
        }

        if (side == 'c') {
            for (unsigned short int i = 0; i < num_gas_species; ++i) {
                for (unsigned short int j = i; j < num_gas_species; ++j) {
                    Dij_res[i][j] = Dij_arr[number_cell][i][j];
                    Dij_res[j][i] = Dij_res[i][j];
                }
            }
        }

        if (side == 'l') {
            for (unsigned short int i = 0; i < num_gas_species; ++i) {
                for (unsigned short int j = i; j < num_gas_species; ++j) {
                    Dij_res[i][j] = Dij_arr_l[number_cell + 1][i][j];
                    Dij_res[j][i] = Dij_res[i][j];
                }
            }
        }
    }
}

void resize_koeff_vectors(int N_x) {
    update_koeffs = 1;
    Cp_arr.resize(N_x);
    forward_arr_save.resize(N_x);
    reverse_arr_save.resize(N_x);
    H_arr.resize(N_x);
    Lambda_arr.resize(N_x);
    Lambda_arr_r.resize(N_x);
    Lambda_arr_l.resize(N_x);
    Dij_arr.resize(N_x);
    Dij_arr_r.resize(N_x);
    Dij_arr_l.resize(N_x);

    for (int i = 0; i < N_x; i++) {
        forward_arr_save[i].resize(num_react);
        reverse_arr_save[i].resize(num_react);
        reverse_arr_save.resize(N_x);
        Cp_arr[i].resize(num_gas_species);
        H_arr[i].resize(num_gas_species);
        Lambda_arr[i].resize(num_gas_species);
        Lambda_arr_r[i].resize(num_gas_species);
        Lambda_arr_l[i].resize(num_gas_species);

        Dij_arr[i].resize(num_gas_species);
        Dij_arr_r[i].resize(num_gas_species);
        Dij_arr_l[i].resize(num_gas_species);
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            Dij_arr[i][k_spec].resize(num_gas_species);
            Dij_arr_r[i][k_spec].resize(num_gas_species);
            Dij_arr_l[i][k_spec].resize(num_gas_species);
        }
    }
}

void updateKoeffs(double* yval, UserData data) {
    int i_T;
    double T_curr, T_prev, T_next;
    flag_use_save_koeffs = 0;
    for (int i = 1; i < data->Nx - 1; i++) {

        i_T = i * (num_gas_species + 2) - 2;
        if (i == data->Nx - 2) {
            T_curr = yval[i_T];
            T_prev = yval[i_T - (num_gas_species + 2)];
            T_next = T_curr;

        }
        else if (i == 1) {
            T_prev = data->Tl;
            T_curr = yval[i_T];
            T_next = yval[i_T + (num_gas_species + 2)];
        }
        else {
            T_curr = yval[i_T];
            T_prev = yval[i_T - (num_gas_species + 2)];
            T_next = yval[i_T + (num_gas_species + 2)];
        }
        if (save_chem_koeffs == 1) {
            MakeYvectors_kins(data, yval, data->Nx, i, data->Tl);
                //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
                //    cout << "Yi = " << Yi[k_spec] << "\n";
                //}
            for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
                Yi[k_spec] -= eps;
                if (i > 1) Yiprev[k_spec] -= eps;
                Yinext[k_spec] -= eps;
                //cout << "Yi = " << Yi[k_spec] << "\n";
            }
            Get_molar_cons(Xi, Yi, T_curr);
            chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr,
                T_curr, Xi, ydot, i);
        }
       

        for (int k_spec = 0; k_spec < num_gas_species; ++k_spec) {
            Cp_arr[i][k_spec] = get_Cpi(k_spec, T_curr, i);
            H_arr[i][k_spec] = get_Hi(k_spec, T_curr, i);

            Lambda_arr[i][k_spec] = get_Lambda5(k_spec, T_curr, i, 'c');
            Lambda_arr_r[i][k_spec] = get_Lambda5(k_spec, (T_curr + T_next) / 2., i, 'r');
            Lambda_arr_l[i][k_spec] = get_Lambda5(k_spec, (T_curr + T_prev) / 2., i, 'l');

            for (int k_spec2 = k_spec; k_spec2 < num_gas_species; ++k_spec2) {
                Dij_arr[i][k_spec][k_spec2] = Dij_func5(k_spec, k_spec2, T_curr, i);
                Dij_arr[i][k_spec2][k_spec] = Dij_arr[i][k_spec][k_spec2];

                Dij_arr_r[i][k_spec][k_spec2] = Dij_func5(k_spec, k_spec2, (T_curr + T_next) / 2., i);
                Dij_arr_r[i][k_spec2][k_spec] = Dij_arr_r[i][k_spec][k_spec2];

                Dij_arr_l[i][k_spec][k_spec2] = Dij_func5(k_spec, k_spec2, (T_curr + T_prev) / 2., i);
                Dij_arr_l[i][k_spec2][k_spec] = Dij_arr_l[i][k_spec][k_spec2];
            }
        }
    }
    flag_use_save_koeffs = 1;
    update_koeffs = 0;
}
