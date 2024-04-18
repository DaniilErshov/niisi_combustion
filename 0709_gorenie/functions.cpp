#include "functions.h"

int flag = 0;
#define RTOL  RCONST(1.0e-5)
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(pow(10, -6))     /* first output time      */
#define TMULT RCONST(1.006)     /* output time factor     */
#define NOUT  1800
#define ZERO  RCONST(0.0)
double get_M(double Tprev, double T, double Tnext, 
    double xprev, double x, double xnext)
{
        double h_left = x - xprev;
        double h = xnext - x;
        double Cp = Cp_all(T, Yi);
        double rho = get_rho(Yi, T);
        double W = get_W(Yi);
        get_grad(gradX, Xi, Xinext, x, xnext);
        double dTdx = (h_left / h / (h + h_left) * Tnext + (h - h_left) / h / h_left * T - h / h_left / (h + h_left) * Tprev);
        double slag_diff = 0;
        double slag_chem = 0.;

        double Vc = 0;
        double YkVk_slag = 0;
        set_Dij_res(T);
        for (int k = 0; k < num_gas_species; k++) {
            YkVk[k] = YkVk_func(k, T, Yi, gradX, Xi);
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
                    Cpn[i] += rho * YkVk[k] * phyc.Cp_coef_hT[k][i] * dTdx;
                }
                else {
                for (int i = 0; i < 9; i++) {
                    Hn[i] += ydot[k] * phyc.Cp_coef_lT[k][i] * my_mol_weight(k);
                    Cpn[i] += rho * YkVk[k] * phyc.Cp_coef_lT[k][i] * dTdx;
                }
            }
            //slag_diff += rho * data->YkVk[k] * myget_Cpi(k, T)  * dTdx;
            //slag_chem += data->ydot[k] * myget_Hi(k, T) * my_mol_weight(k);
        }

        slag_diff = get_dCpi(Cpn, T);
        slag_chem = get_dHiRT(Hn, T);

        make_averageY(X_tmp, Xi, Xinext);
        double lambda_right = Lambda_All(X_tmp, (Tnext + T) / 2.);
        make_averageY(X_tmp, Xi, Xiprev);
        double lambda_left = Lambda_All(X_tmp, (T + Tprev) / 2.);

        return ((2. / (h + h_left)) *
            (lambda_right * (Tnext - T) / h
                - lambda_left * (T - Tprev) / h_left) - slag_chem - slag_diff) / Cp / dTdx;
}

int InitialData(int& Nx, vector<double>& x_vect, vector<double>& T_vect, vector<double>& Y_vect, double& M, double Tstart, double Tfinish, double* Ystart, double* Yend)
{
    double h = l / (Nx - 1);
    //double x_start = l / 3.;
    //double x_finish = 2. * l / 3.;
    double x_start = 0.1;
    double x_finish = 0.2;
    int dN = (x_finish - x_start) / h;

    double j = 0;
    std::cout << "M = " << M << "\n";

    //for (int i = 0; i < Nx; i++) {
    //    x_vect[i] = h * i;
    //}
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
    T_center = (Tstart + Tfinish) / 2.;
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

void MakeYvectors(UserData data,
    double* Y, int myNx, int i, double Tl) {
    //cout << "i = " << i << "\n";
    for (int j = 0; j < num_gas_species; j++) {

        Yi[j] = Y[j + (i - 1) * (num_gas_species) + i - 1];
        if (i == 1) Yiprev[j] = Y_left_bound[j];
        else Yiprev[j] = Y[j + (i - 2) * (num_gas_species) + i - 2];

        if (i == myNx - 2)  Yinext[j] = Yi[j];
        else  Yinext[j] = Y[j + (i) * (num_gas_species) + i];
    }
}

void MakeYvectorsY(UserData data,
    double* Y, int myNx, int i, double Tl) {
    //cout << "i = " << i << "\n";
    for (int j = 0; j < num_gas_species; j++) {

        Yi[j] = Y[j + (i - 1) * num_gas_species ];
        if (i == 1) Yiprev[j] = Y_left_bound[j];
        else Yiprev[j] = Y[j + (i - 2) * num_gas_species];

        if (i == myNx - 2)  Yinext[j] = Yi[j];
        else  Yinext[j] = Y[j + (i)*num_gas_species];
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
        YkVk_slag = YkVk_func(k, T, Yi, gradX, Xi);
        Vc -= YkVk_slag;
    }

    for (int k = 0; k < num_gas_species; k++) {
        YkVk_slag = YkVk_func(k, T, Yi, gradX, Xi) + Vc * Yi[k];
        slag_diff += rho * YkVk_slag * myget_Cpi(k, T) * dTdx;
        slag_chem += ydot[k] * myget_Hi(k, T) * my_mol_weight(k);
    }

    make_averageY(X_tmp, Xi, Xinext);
    double lambda_right = Lambda_All(X_tmp, (Tnext + T) / 2.);

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
    double Cp = Cp_all(T, Yi);
    double rho = get_rho(Yi, T);
    double W = get_W(Yi);
    get_grad(gradX, Xi, Xinext, x, xnext);
    double dTdx = (h_left / h / (h + h_left) * Tnext + (h - h_left) / h / h_left * T - h / h_left / (h + h_left) * Tprev);
    double slag_diff = 0;
    double slag_chem = 0.;

    double Vc = 0;
    double YkVk_slag = 0;
    set_Dij_res(T);
    for (int k = 0; k < num_gas_species; k++) {
        YkVk[k] = YkVk_func(k, T, Yi, gradX, Xi);
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
                Cpn[i] += rho * YkVk[k] * phyc.Cp_coef_hT[k][i] * dTdx;
            }
            else {
            for (int i = 0; i < 9; i++) {
                Hn[i] += ydot[k] * phyc.Cp_coef_lT[k][i] * my_mol_weight(k);
                Cpn[i] += rho * YkVk[k] * phyc.Cp_coef_lT[k][i] * dTdx;
            }
        }
        //slag_diff += rho * data->YkVk[k] * myget_Cpi(k, T)  * dTdx;
        //slag_chem += data->ydot[k] * myget_Hi(k, T) * my_mol_weight(k);
    }

    slag_diff = get_dCpi(Cpn, T);
    slag_chem = get_dHiRT(Hn, T);

    make_averageY(X_tmp, Xi, Xinext);
    double lambda_right = Lambda_All(X_tmp, (Tnext + T) / 2.);
    make_averageY(X_tmp, Xi, Xiprev);
    double lambda_left = Lambda_All(X_tmp, (T + Tprev) / 2.);

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

    slag_diff = (rhoYkVk_r - rhoYkVk_l) / (x_r - x_l);
    slag_chem = -my_mol_weight(k_spec) * ydot[k_spec];
    //cout << "slag_chem =  " << slag_chem << "\n";
    //cout << "slag_diff =  " << slag_diff << "\n";

    return data->M * (h_left / h / (h + h_left) * Yinext[k_spec] + (h - h_left) / h / h_left * Yi[k_spec] - h / h_left / (h + h_left) * Yiprev[k_spec])
      + slag_chem + slag_diff;
}

void Write_to_file(string str, ofstream& fout, vector<double>& x_vect,
    vector<double>& T_vect, vector<double>& Y_vect, vector<double>& Yp_vect, double M, int N_x, int number) {
    double x_start, x_finish, D;
    double rho;
    fout.open(str + ".dat");
    string title2 = R"(VARIABLES= "x", "T")";
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
            fout << x_vect[i] << "  " << T_vect[i];
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
    rtol = RCONST(1.0e-3);
    atval = N_VGetArrayPointer(avtol);
    consval = N_VGetArrayPointer(cons);
    int number_spec = 0;
    for (int i = 0; i < NEQ; i++) {
        consval[i] = 0.0;   /*constraint*/
        yval[i] = Y_vect[i + num_gas_species];
        ypval[i] = 0;
        atval[i] = RCONST(1.0e-5);
    }
    consval[1] = 1.0;
    /* Create and initialize  y, y', and absolute tolerance vectors. */
    t0 = ZERO;
    data->M = M;
    double rho;
    ofstream f_ida;
    double sumY = 0;
        
       
    MakeYvectorsY(data,
        yval, N_x, 1, T_vect[0]);

    Get_molar_cons(Xi, Yiprev, T_vect[0]);

    chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr, T_vect[0], Xi, ydot);
    Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);

    std::cout << "M_in_IDA_Y = " << M << "\n";
    std::cout << "Nx _ Y  =  " << data->Nx << "\n";
    for (int i = 1; i < data->Nx - 1; i++) {
        MakeYvectorsY(data, 
            yval, N_x, i, T_vect[i]);

        Get_molar_cons(Xi, Yi, T_vect[i]);

        chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr, T_vect[i], Xi, ydot);
        Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);
        rho = get_rho(Yi, T_vect[i]);
        sumY = 0;
        find_diff_slag(data, T_vect[i], T_vect[i + 1], Yi, Yinext, Xi, Xinext, YkVk_r, Y_tmp_r, X_tmp_r, gradX_r, data->rho_r, data->Vc_r, i);

        find_diff_slag(data, T_vect[i - 1], T_vect[i], Yiprev, Yi, Xiprev, Xi, YkVk_l, Y_tmp_l, X_tmp_l, gradX_l, data->rho_l, data->Vc_l, i - 1);
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            ypval[k_spec + (i - 1) * num_gas_species] = -F_rightY(data, k_spec,
                T_vect[i - 1], T_vect[i], T_vect[i + 1],
                x_vect[i - 1], x_vect[i], x_vect[i + 1]) / rho;


        }
    
    }
   
    
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


    int mu = num_gas_species;
    int ml = num_gas_species;

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

    double sum = 100;
    while (iout < 80) {
        retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);

        for (int i = 0; i < NEQ; i++) {
            Y_vect[i + num_gas_species] = yval[i];
        }
    
        for (int i = NEQ; i < NEQ + num_gas_species - 1; i++) {
            Y_vect[i] = Y_vect[i - num_gas_species];
        }
        iout++;
        tout += tout1;
    }
    std::cout << "t = " << tout << "\n";

    for (int i = 0; i < NEQ; i++) {
        Y_vect[i + num_gas_species] = yval[i];
        if (Y_vect[i + num_gas_species] < 0) Y_vect[i + num_gas_species] = 0;
    }

    for (int i = NEQ ; i <  NEQ + num_gas_species ; i++) {
        Y_vect[i + num_gas_species] = Y_vect[i];
    }


    /* Print final statistics to the screen */
    cout << "\nFinal Statistics:\n";
    retval = IDAPrintAllStats(mem, stdout, SUN_OUTPUTFORMAT_TABLE);
    
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
    double Temp;
    double M;
    double sumY = 0;
    data = (UserData)user_data;
    int j;
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
    MakeYvectorsY(data, yval, myNx, 1, T_vect[0]);
    Get_mole_fr(Xiprev, Yiprev);
    Get_mole_fr(Xi, Yi);
    Get_mole_fr(Xinext, Yinext);

    make_averageY(Y_tmp, Yiprev, Yi);
    make_averageY(X_tmp, Xiprev, Xi);

    get_grad(gradX, Xiprev, Xi, x_vect[0], x_vect[1]);

    //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
    //    Vc -= YkVk(k_spec, T_vect[0], data->Y_tmp, data->gradX, data->X_tmp);
    //}
    for (int i = 1; i < myNx - 1; i++) {
        //cout << "\n\n\n\nin func i = " << i << "\n";
        MakeYvectorsY(data, yval, myNx, i, T_vect[i]);
        Get_molar_cons(Xi, Yi, T_vect[i]);
        chem_vel(Sn, Hn,forward_arr, reverse_arr, equilib_arr, T_vect[i], Xi, ydot);

        
        Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);


        find_diff_slag(data, T_vect[i], T_vect[i + 1], Yi, Yinext, Xi, Xinext, YkVk_r, Y_tmp_r, X_tmp_r, gradX_r, rho_r, Vc_r, i);

        find_diff_slag(data, T_vect[i - 1], T_vect[i], Yiprev, Yi, Xiprev, Xi, YkVk_l, Y_tmp_l, X_tmp_l, gradX_l, rho_l, Vc_l, i - 1);

        rho = get_rho(Yi, T_vect[i]);
        sumY = 0;
        //cout << "i = " << i << "n";
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            //cout << "name = " << name_species[k_spec] << "\n";
            //cout << "yval = " << yval[k_spec] << "\n";
            //if (k_spec != N2)
            rval[k_spec + (i - 1) * num_gas_species] = rho * ypval[k_spec + (i - 1) * num_gas_species] + F_rightY(data, k_spec,
                T_vect[i - 1], T_vect[i], T_vect[i + 1],
                x_vect[i - 1], x_vect[i], x_vect[i + 1]);
        }
        //cout << "\n";
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
    std::cout << "M_in_IDA = " << M << "\n";
    int i_temp = 0;
    for (int i = 1; i < N_x - 1; i++) {
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            consval[i_temp] = 0.0;   /*constraint*/
            yval[i_temp] = Y_vect[i * num_gas_species + k_spec] ;
            ypval[i_temp] = 0;
            atval[i_temp] = RCONST(1.0e-6);
            //cout << "Yi = " << i_temp << "   = " << yval[i_temp] << endl;
            i_temp++;
        }
        consval[i_temp] = 1.0;   /*constraint*/
        atval[i_temp] = pow(10, 0);
        yval[i_temp] = T_vect[i];
        ypval[i_temp] = 0;
        //cout << "Ti = " << i_temp << "   = " << yval[i_temp] << endl << endl;
        i_temp++;
        //if ((i + 1) % num_gas_species == 0) cout << endl;
    }

    t0 = ZERO;

    double rho;
    ofstream f_ida;
    double sumY = 0;

    flag = 1;

    cout << "Nx_ALL = " << data->Nx << "\n";
    i_temp = 0;

    for (int i = 1; i < data->Nx - 1; i++) {

        MakeYvectors(data,
            yval, N_x, i, T_vect[i]);
        //cout << "i = " << i << "\n";
        //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        //    cout << "Yiprev = " << Yiprev[k_spec] << "\n";
        //    cout << "Ycurr = " << Yi[k_spec] << "\n";
        //    cout << "Yinext = " << Yinext[k_spec] << "\n\n";
        //}
        //cout << "\n\n";
        //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        //    Yi[k_spec] -= eps;
        //    if (i > 1) Yiprev[k_spec] -= eps;
        //    Yinext[k_spec] -= eps;
        //    //cout << "Yi = " << Yi[k_spec] << "\n";
        //}
        Get_molar_cons(Xi, Yi, T_vect[i]);
        chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr, T_vect[i], Xi, ydot);

        Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);
        rho = get_rho(Yi, T_vect[i]);

        double Cp = get_Cp(num_gas_species, Yi, T_vect[i]);
        find_diff_slag(data, T_vect[i], T_vect[i + 1], Yi, Yinext, Xi, Xinext, YkVk_r, Y_tmp_r, X_tmp_r, gradX_r, data->rho_r, data->Vc_r, i);

        find_diff_slag(data, T_vect[i - 1], T_vect[i], Yiprev, Yi, Xiprev, Xi, YkVk_l, Y_tmp_l, X_tmp_l, gradX_l, data->rho_l, data->Vc_l, i - 1);

        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            ypval[i_temp] = -F_rightY(data, k_spec,
                T_vect[i - 1], T_vect[i], T_vect[i + 1],
                x_vect[i - 1], x_vect[i], x_vect[i + 1]
            ) / rho;
            //cout << "i_temp = " << i_temp << " = " << ypval[i_temp] << "\n";
            i_temp++;

        }
        ypval[i_temp] = -F_right(data,
            T_vect[i - 1], T_vect[i], T_vect[i + 1],
            x_vect[i - 1], x_vect[i], x_vect[i + 1]) / rho / Cp;

        if (i == N_center) ypval[i_temp] = 0;
        //cout << "i_temp T = " << i_temp << " = " << ypval[i_temp] << "\n";
        i_temp++;
        //if (i != N_center) {
        //    ypval[i + NEQ_Y - 1] = -F_right(data,
        //        T_vect[i - 1], T_vect[i], T_vect[i + 1],
        //        x_vect[i - 1], x_vect[i], x_vect[i + 1]) / rho / Cp;
        //}
        //else {
        //    ypval[i + NEQ_Y - 1] = -F_right(data,
        //        T_vect[i - 1], T_vect[i], T_vect[i + 1],
        //        x_vect[i - 1], x_vect[i], x_vect[i + 1]);
        //}

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
    int mu = 2 * (num_gas_species + 1);
    int ml = 2 * (num_gas_species + 1);
    A = SUNBandMatrix(NEQ, mu, ml, ctx);
    if (check_retval((void*)A, "SUNBandMatrix", 0)) return(1);

    /* Create banded SUNLinearSolver object */
    LS = SUNLinSol_Band(yy, A, ctx);
    if (check_retval((void*)LS, "SUNLinSol_Band", 0)) return(1);

    /* Attach the matrix and linear solver */
    retval = IDASetLinearSolver(mem, LS, A);
    if (check_retval(&retval, "IDASetLinearSolver", 1)) return(1);

    IDASetMaxOrd(mem, 1);
    //NLS = SUNNonlinSol_Newton(yy, ctx);
    //if (check_retval((void*)NLS, "SUNNonlinSol_Newton", 0)) return(1);

    ///* Attach the nonlinear solver */
    //retval = IDASetNonlinearSolver(mem, NLS);
    //if (check_retval(&retval, "IDASetNonlinearSolver", 1)) return(1);
    //IDASetInitStep(mem, pow(10, -7));

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

    for (int i = 0; i < num_gas_species * N_x; i++)
        Yp_vect[i] = 0;

    double nevyaz_Y_IDA = 1000, nevyaz_T_IDA = 1000;

    data->N_m = 0;
    while (iout < ida_steps) {
        retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);
        int i_center = N_center * (num_gas_species)+N_center - 1;
        int i_left = i_center - (num_gas_species + 1);
        int i_right = i_center + (num_gas_species + 1);

        MakeYvectors(data, yval, N_x, N_center, data->Tl);
        //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        //    Yi[k_spec] -= eps;
        //    Yiprev[k_spec] -= eps;
        //    Yinext[k_spec] -= eps;
        //    //cout << "Yi = " << Yi[k_spec] << "\n";
        //}
        Get_molar_cons(Xi, Yi, yval[i_center]);

        chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr,
            yval[i_center], Xi, ydot);


        Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);

        M = get_M(yval[i_left], yval[i_center], yval[i_right],
            x_vect[N_center - 1], x_vect[N_center], x_vect[N_center + 1]);

        data->M = M;
        cout << "t = " << tout << "\n";
        cout << "M = " << M << "\n\n";
        if (check_retval(&retval, "IDASolve", 1)) return(1);
        iout++;
        tout1 *= 1;
        tout += tout1;
    }
    i_temp = 0;
    for (int i = 1; i < N_x - 1; i++) {
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            Y_vect[i * num_gas_species + k_spec] = yval[i_temp] - eps;
            if (Y_vect[i * num_gas_species + k_spec] < 0) Y_vect[i * num_gas_species + k_spec] = pow(10, -15);
            i_temp++;
        }
        T_vect[i] = yval[i_temp];
        i_temp++;
    }

    for (int i = 0; i < num_gas_species; i++) {
        Y_vect[(num_gas_species) * (N_x - 1) + i] = Y_vect[(num_gas_species) * (N_x - 2) + i];
    }

    T_vect[N_x - 1] = T_vect[N_x - 2];

    std::cout << "tout = " << tout << "\n";
    std::cout << "end IDA M = " << M << "\n";
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

void find_diff_slag(UserData data, double Tcurr, double Tnext, double* Yi, double* Yinext,
    double* Xi, double* Xinext, double* Ykvk_side, double* Y_tmp_side, double* X_tmp_side, double* gradX_side, double& rho_side, double& Vc_side, int i) {
    make_averageY(Y_tmp_side, Yi, Yinext);
    make_averageY(X_tmp_side, Xi, Xinext);
    get_grad(gradX_side, Xi, Xinext, x_vect[i], x_vect[i + 1]);
    rho_side = get_rho(Y_tmp_side, (Tcurr + Tnext) / 2.);
    Vc_side = 0;
    set_Dij_res((Tcurr + Tnext) / 2.);
    for (int k = 0; k < num_gas_species; k++) {
        Ykvk_side[k] = YkVk_func(k, (Tcurr + Tnext) / 2., Y_tmp_side, gradX_side, X_tmp_side);
        Vc_side -= Ykvk_side[k];
        //cout << "Ykvk = " << Ykvk_side[k] << "\n";
    }

}

static int func_All_IDA(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
    realtype* yval, * ypval, * rval;
    UserData data;
    double Temp;
    double M;
    double sumY = 0;
    int i_T;
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

    //data->M = yval[data->N_centr + data->NEQ_Y - 1];
    int i_temp = 0;
    for (int i = 1; i < myNx - 1; i++) {

        i_T = i * (num_gas_species)+i - 1;

        if (i == myNx - 2) {
            T_curr = yval[i_T];
            T_prev = yval[i_T - (num_gas_species + 1)];
            T_next = T_curr;

        }
        else if (i == 1) {
            T_prev = data->Tl;
            T_curr = yval[i_T];
            T_next = yval[i_T + (num_gas_species + 1)];
        }
        else {
            T_curr = yval[i_T];
            T_prev = yval[i_T - (num_gas_species + 1)];
            T_next = yval[i_T + (num_gas_species + 1)];
        };

        //cout << "in func i = " << i << "\n";
        //cout << "T_prev" << T_prev << "\n";
        //cout << "T_curr" << T_curr << "\n";
        //cout << "T_next" << T_next << "\n\n";
        MakeYvectors(data, yval, myNx, i, data->Tl);
        //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        //    Yi[k_spec] -= eps;
        //    if (i > 1) Yiprev[k_spec] -= eps;
        //    Yinext[k_spec] -= eps;
        //    //cout << "Yi = " << Yi[k_spec] << "\n";
        //}

        Get_molar_cons(Xi, Yi, T_curr);
        
        chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr, T_curr, Xi, ydot);


        Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);

        find_diff_slag(data, T_curr, T_next, Yi, Yinext, Xi, Xinext, YkVk_r, Y_tmp_r, X_tmp_r, gradX_r, data->rho_r, data->Vc_r, i);

        find_diff_slag(data, T_prev, T_curr, Yiprev, Yi, Xiprev, Xi, YkVk_l, Y_tmp_l, X_tmp_l, gradX_l, data->rho_l, data->Vc_l, i - 1);

        rho = get_rho(Yi, T_curr);
        sumY = 0;
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            rval[i_temp] = rho * ypval[i_temp] + F_rightY(data, k_spec,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]
            );
            //cout << "i_temp = " << i_temp << " = " << ypval[i_temp] << "\n";
            i_temp++;
        }
        Cp = get_Cp(num_gas_species, Yi, T_curr);

        rval[i_temp] = rho * Cp * ypval[i_temp] + F_right(data,
            T_prev, T_curr, T_next,
            x_vect[i - 1], x_vect[i], x_vect[i + 1]);
        //cout << "i_temp = " << i_temp << " = " << ypval[i_temp] << "\n";
        if (i == Ncentr)
        {
            rval[i_temp] = yval[i_temp] - data->T_center;
        }



        i_temp++;
        /* if (i != data->N_centr)
             rval[i_T] = rho * Cp * ypval[i_T] + F_right(data,
                 T_prev, T_curr, T_next,
                 x_vect[i - 1], x_vect[i], x_vect[i + 1]
                 );
         else
             rval[i_T] = F_right(data,
                 T_prev, T_curr, T_next,
                 x_vect[i - 1], x_vect[i], x_vect[i + 1]);*/
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


    int NEQ = (N_x - 2) * ((num_gas_species + 2));
    int NEQ_Y = 0;
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
    data->M = M;
    std::cout << "INTEGRATE KINSOL = " << M << "\n";
    int i_temp = 0;
    for (int i = 1; i < N_x - 1; i++) {
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            cval[i_temp] = 1.0;   /*constraint*/
            yval[i_temp] = Y_vect[i * num_gas_species + k_spec] + eps;
            sval[i_temp] = pow(10, 0);
            //cout << "Yi = " << i_temp << "   = " << Y_vect[i * num_gas_species + k_spec] << endl;
            i_temp++;
        }
        yval[i_temp] = T_vect[i];
        cval[i_temp] = 1.0;
        sval[i_temp] = pow(10, -3);
        //cout << "Ti = " << i_temp << "   = " << yval[i_temp] << endl << endl;
        i_temp++;
        yval[i_temp] = M;
        cval[i_temp] = 1.0;
        sval[i_temp] = pow(10, -3);
        //cout << "Mi = " << i_temp << "   = " << yval[i_temp] << endl << endl;
        i_temp++;

        //if ((i + 1) % num_gas_species == 0) cout << endl;
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
    KINSetNumMaxIters(kmem, 1000);
    retval = KINInit(kmem, func_kinsol, res_vect);
    if (check_retval(&retval, "KINInit", 1)) return(1);

    /* Create dense SUNMatrix */
    int mu =  2 * (num_gas_species + 2);
    int ml = 2 * (num_gas_species + 2);

    J = SUNBandMatrix(NEQ, mu, ml, sunctx);
    if (check_retval((void*)J, "SUNBandMatrix", 0)) return(1);

    /* Create banded SUNLinearSolver object */
    LS = SUNLinSol_Band(res_vect, J, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_Band", 0)) return(1);

    /* Attach the matrix and linear solver to KINSOL */
    retval = KINSetLinearSolver(kmem, LS, J);
    if (check_retval(&retval, "KINSetLinearSolver", 1)) return(1);

    glstr = 0;
    mset = 1000;

    retval = KINSol(kmem, res_vect, glstr, s, s);
    if (check_retval(&retval, "KINSol", 1)) return(1);
    cout << "retval = " << retval << "\n";
    i_temp = 0;
    for (int i = 1; i < N_x - 1; i++) {
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            Y_vect[i * num_gas_species + k_spec] = yval[i_temp] - eps;
            i_temp++;
        }
        T_vect[i] = yval[i_temp];
        i_temp++;
        M = yval[i_temp];
        //cout << "i = " << i << "M = " << yval[i_temp] << "\n";
        i_temp++;
    }
    for (int i = 0; i < num_gas_species; i++) {
        Y_vect[(num_gas_species) * (N_x - 1) + i] = Y_vect[(num_gas_species) * (N_x - 2) + i];
    }

    T_vect[N_x - 1] = T_vect[N_x - 2];
    cout << "M = " << M << "\n";
    cout << "N_x = " << N_x << "\n";
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
    double Temp;
    double M;
    double sumY = 0;
    int i_T;
    data = (UserData)user_data;
    double Cp;
    double T_curr, T_prev, T_next;
    int j;
    int myNx = data->Nx;
    int myNeq = data->NEQ;
    int myNm = data->N_m;
    yval = N_VGetArrayPointer(u);
    rval = N_VGetArrayPointer(f);
    int Ncentr = data->N_centr;
    double rho;
    int i_temp = 0;

    for (int i = 1; i < myNx - 1; i++) {

        i_T = i * (num_gas_species + 2) - 2;
        data->M = yval[i * (num_gas_species + 2) - 1];
        //cout << "in func i = " << i << "\n";
        //cout << "M = " << data->M << "\n\n";
        if (i == myNx - 2) {
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
        };

        //cout << "T_prev" << T_prev << "\n";
        //cout << "T_curr" << T_curr << "\n";
        //cout << "T_next" << T_next << "\n\n";
        MakeYvectors_kins(data, yval, myNx, i, data->Tl);

       for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
          Yi[k_spec] -= eps;
          if (i > 1) Yiprev[k_spec] -= eps;
          Yinext[k_spec] -= eps;
          //cout << "Yi = " << Yi[k_spec] << "\n";
       }
        Get_molar_cons(Xi, Yi, T_curr);

        chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr,
            T_curr, Xi, ydot);

        Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);

        find_diff_slag(data, T_curr, T_next, Yi, Yinext, Xi, Xinext, YkVk_r, Y_tmp_r, X_tmp_r, gradX_r, data->rho_r, data->Vc_r, i);

        find_diff_slag(data, T_prev, T_curr, Yiprev, Yi, Xiprev, Xi, YkVk_l, Y_tmp_l, X_tmp_l, gradX_l, data->rho_l, data->Vc_l, i - 1);

        rho = get_rho(Yi, T_curr);
        sumY = 0;
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            rval[i_temp] = F_rightY(data, k_spec,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]
            );
            //cout << "i_temp = " << i_temp << " = " << yval[i_temp] << "\n";
            i_temp++;
        }
        Cp = get_Cp(num_gas_species, Yi, T_curr);

        //cout << "i_temp = " << i_temp << " = " << yval[i_temp] << "\n\n";
        if (i == Ncentr)
        {
            rval[i_temp] = yval[i_temp] - data->T_center;
        }
        else {
            rval[i_temp] = F_right(data,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]);
        }
        i_temp++;

        if (i == Ncentr) {
            rval[i_temp] = F_right(data,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]);
        }
        if (i > Ncentr){
            rval[i_temp] = yval[i_temp] - yval[i_temp - (num_gas_species + 2)];
        }
        if (i < Ncentr) {
            rval[i_temp] = yval[i_temp] - yval[i_temp + (num_gas_species + 2)];
        }
        i_temp++;
    }
    return 0;
}

void MakeYvectors_kins(UserData data,
    double* Y, int myNx, int i, double Tl) {
    //cout << "i = " << i << "\n";
    for (int j = 0; j < num_gas_species; j++) {

        Yi[j] = Y[j + (i - 1) * (num_gas_species + 1) + i - 1];
        if (i == 1) Yiprev[j] = Y_left_bound[j];
        else Yiprev[j] = Y[j + (i - 2) * (num_gas_species + 1) + i - 2];

        if (i == myNx - 2)  Yinext[j] = Yi[j];
        else  Yinext[j] = Y[j + (i) * (num_gas_species + 1)+i];
    }
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

    data->NEQ = NEQ;
    data->Tl = Tinitial;

    for (int i = 0; i < num_gas_species; i++) {
        Yi[i] = Y_vect[i];
        Yiprev[i] = Y_vect[i];
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
    Ith(yy, NEQ) = Tinitial * 8;

    double sumY = 0;
    Get_molar_cons(Xi, Y_vect, Tinitial);

    chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr, Tinitial, Xi, ydot);

    double rho = get_rho(Y_vect, Tinitial);

    for (int i = 1; i < NEQ; i++) {
        Ith(yp, i) = ydot[i - 1];
        cout << Ith(yp, i) << "\n";
    }

    double sum = 0;
    double Cp = get_Cp(num_gas_species, Y_vect, Tinitial);
    cout << "Cp = " << Cp << "\n";
    sum = get_enthalpy(num_gas_species, Y_vect, Tinitial);
    double sum2 = get_enthalpy(num_gas_species, Y_vect, Ith(yy, NEQ));
    cout << "enthalpy = " << sum << "\n";
    cout << "rho = " << get_rho(Y_vect, Tinitial) << "\n";
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
    //for (int i = 1; i <= NEQ; i++) {
    //    Ith(cons, i) = 1.0;
    //}

    //retval = IDASetConstraints(mem, cons);
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
    double tout1 = 2 * pow(10, 1);
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
        cout << komponents_str[i] << " = " << Y_end[i] << "\n";
    }

    Tend = Ith(yy, num_gas_species + 1);
    cout << "Tend = " << Tend << "\n";
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

static int func_final_state(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
    realtype* yval, * ypval, * rval;
    UserData data;
    double Temp;
    double M;
    data = (UserData)user_data;
    int j;
    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rval = N_VGetArrayPointer(rr);

    for (int j = 0; j < num_gas_species; j++) {
        Yi[j] = yval[j];
    }
    Temp = yval[num_gas_species];
    Get_molar_cons(Xi, Yi, Temp);
    chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr, Temp, Xi, ydot);

    double rho = get_rho(Yi, Temp);
    for (j = 0; j < num_gas_species; j++) {
        rval[j] = rho * ypval[j] - ydot[j] * phyc.mol_weight[j];
    }
    double sum = 0;

    double Cp = get_Cp(num_gas_species, Yi, Temp);
    double sum2;

    sum2 = get_enthalpy(num_gas_species, Yiprev, data->Tl);

    sum = get_enthalpy(num_gas_species, Yi, Temp);

    rval[num_gas_species] = rho * Cp * ypval[num_gas_species] + (sum - sum2);
    return(0);
}

void makeYstart(double koeff_topl, double* Ystart) {

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
    Ystart[komponents["H2"]] = 1.988419E-9;
    Ystart[komponents["O2"]] = 0.2184336;
    Ystart[komponents["N2"]] = 0.7193815;
    Ystart[komponents["NC7H16"]] = 0.06218489;
}

void set_polynom(double** ploynom, std::string name_file, std::string type_polynom) {
    std::string line;
    std::string koeff_str = "COEFFICIENTS";
    std::ifstream in(name_file); // окрываем файл для чтения
    if (in.is_open())
    {
        while (getline(in, line))
        {
            vector splitted_string = splitString(line, ' ');
            if (std::find(splitted_string.begin(), splitted_string.end(), koeff_str) != splitted_string.end()
                && std::find(splitted_string.begin(), splitted_string.end(), type_polynom) != splitted_string.end())
            {
                getline(in, line);
                while (getline(in, line))
                {
                    vector splitted_string = splitString(line, ' ');
                    if (std::find(splitted_string.begin(), splitted_string.end(), koeff_str) != splitted_string.end())
                        return;


                    if (splitted_string.size() != 0)
                    {
                        
                        string specie = splitted_string[0];
                        std::transform(specie.begin(), specie.end(), specie.begin(), ::toupper);
                        if (komponents.contains(specie))
                        {                            
                            int i_komp = komponents[specie];
                            cout << "line = " << line << "\n";
                            cout << "i_komp = " << i_komp << "\n";

                            ploynom[i_komp][0] = stod(splitted_string[1]);
                            ploynom[i_komp][1] = stod(splitted_string[2]);
                            ploynom[i_komp][2] = stod(splitted_string[3]);
                            ploynom[i_komp][3] = stod(splitted_string[4]);
                            cout << specie << " " << i_komp << " " << ploynom[i_komp][0]
                                << " " << ploynom[i_komp][1]
                                << " " << ploynom[i_komp][2]
                                << " " << ploynom[i_komp][3] << "\n";
                        }

                    }

                }
            }


        }
    }
    in.close();
}

void set_polynom_diffusion(double*** polynom, std::string name_file, std::string type_polynom) {
    std::string line;
    std::string koeff_str = "COEFFICIENTS";
    std::ifstream in(name_file); // окрываем файл для чтения
    if (in.is_open())
    {
        while (getline(in, line))
        {
            vector splitted_string = splitString(line, ' ');
            if (std::find(splitted_string.begin(), splitted_string.end(), koeff_str) != splitted_string.end()
                && std::find(splitted_string.begin(), splitted_string.end(), type_polynom) != splitted_string.end())
            {
                getline(in, line);
                while (getline(in, line))
                {
                    vector splitted_string = splitString(line, ' ');
                    if (std::find(splitted_string.begin(), splitted_string.end(), koeff_str) != splitted_string.end())
                        return;

                    if (splitted_string.size() != 0)
                    {
                        string specie1 = splitted_string[0];
                        std::transform(specie1.begin(), specie1.end(), specie1.begin(), ::toupper);

                        string specie2 = splitted_string[1];
                        std::transform(specie2.begin(), specie2.end(), specie2.begin(), ::toupper);

                        if (komponents.contains(specie1) && komponents.contains(specie2))
                        {
                            int i_sp1 = komponents[specie1];
                            int i_sp2 = komponents[specie2];
                            polynom[i_sp1][i_sp2][0] = stod(splitted_string[2]);
                            polynom[i_sp1][i_sp2][1] = stod(splitted_string[3]);
                            polynom[i_sp1][i_sp2][2] = stod(splitted_string[4]);
                            polynom[i_sp1][i_sp2][3] = stod(splitted_string[5]);

                            polynom[i_sp2][i_sp1][0] = stod(splitted_string[2]);
                            polynom[i_sp2][i_sp1][1] = stod(splitted_string[3]);
                            polynom[i_sp2][i_sp1][2] = stod(splitted_string[4]);
                            polynom[i_sp2][i_sp1][3] = stod(splitted_string[5]);

                            /*cout << komponents[specie1] << " " << komponents[specie2]
                                << " " << polynom[i_sp1][i_sp2][0]
                                << " " << polynom[i_sp1][i_sp2][1]
                                << " " << polynom[i_sp1][i_sp2][2]
                                << " " << polynom[i_sp1][i_sp2][3] << "\n";*/
                        }

                    }

                }
            }
        }
    }
    in.close();
}
std::vector<std::string> splitString(std::string str, char splitter) {
    std::vector<std::string> result;
    std::string current = "";
    for (int i = 0; i < str.size(); i++) {
        if (str[i] == splitter) {
            if (current != "") {
                result.push_back(current);
                current = "";
            }
            continue;
        }
        current += str[i];
    }
    if (current.size() != 0)
        result.push_back(current);
    return result;
}

template <typename T>
void findValue(const std::vector<T>& data, bool(*condition)(T))
{
    auto result{ std::find_if(begin(data), end(data), condition) };
    if (result == end(data))
        std::cout << "Value not found" << std::endl;
    else
        std::cout << "Value found at position " << (result - begin(data)) << std::endl;
}



int Integrate_Kinsol_dense(int N_x, vector<double>& x_vect,
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
    data->M = M;
    std::cout << "INTEGRATE KINSOL = " << M << "\n";
    int i_temp = 0;
    for (int i = 1; i < N_x - 1; i++) {
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            cval[i_temp] = 1.0;   /*constraint*/
            yval[i_temp] = Y_vect[i * num_gas_species + k_spec] + eps;
            sval[i_temp] = pow(10, 0);
            cout << "Yi = " << i_temp << "   = " << yval[i_temp] << endl;
            i_temp++;
        }
        yval[i_temp] = T_vect[i];
        if (i == N_center)  yval[i_temp] = M;
        cval[i_temp] = 1.0;
        sval[i_temp] = pow(10, -5);
        cout << "Ti = " << i_temp << "   = " << yval[i_temp] << endl << endl;
        i_temp++;
        //if ((i + 1) % num_gas_species == 0) cout << endl;
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

    retval = KINInit(kmem, func_kinsol_dense, res_vect);
    if (check_retval(&retval, "KINInit", 1)) return(1);

    /* Create dense SUNMatrix */
    int mu = 2 * (num_gas_species);
    int ml = 2 * (num_gas_species);

    J = SUNDenseMatrix(NEQ, NEQ, sunctx);
    if (check_retval((void*)J, "SUNBandMatrix", 0)) return(1);

    /* Create banded SUNLinearSolver object */
    LS = SUNLinSol_Dense(res_vect, J, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_Band", 0)) return(1);

    /* Attach the matrix and linear solver to KINSOL */
    retval = KINSetLinearSolver(kmem, LS, J);
    if (check_retval(&retval, "KINSetLinearSolver", 1)) return(1);

    glstr = 0;
    mset = 1000;

    retval = KINSol(kmem, res_vect, glstr, s, s);
    if (check_retval(&retval, "KINSol", 1)) return(1);
    cout << "retval = " << retval << "\n";
    i_temp = 0;
    for (int i = 1; i < N_x - 1; i++) {
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            Y_vect[i * num_gas_species + k_spec] = yval[i_temp] - eps;
            i_temp++;
        }
        if (i != N_center)
            T_vect[i] = yval[i_temp];
        else 
            M = yval[i_temp];
        i_temp++;
    }

    for (int i = 0; i < num_gas_species; i++) {
        Y_vect[(num_gas_species) * (N_x - 1) + i] = Y_vect[(num_gas_species) * (N_x - 2) + i];
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

static int func_kinsol_dense(N_Vector u, N_Vector f, void* user_data)
{
    realtype* yval, * ypval, * rval;
    UserData data;
    double Temp;
    double M;
    double sumY = 0;
    int i_T;
    data = (UserData)user_data;
    double Cp;
    double T_curr, T_prev, T_next;
    int j;
    int myNx = data->Nx;
    int myNeq = data->NEQ;
    int myNm = data->N_m;
    yval = N_VGetArrayPointer(u);
    rval = N_VGetArrayPointer(f);
    int Ncentr = data->N_centr;
    double rho;
    data->M = yval[(Ncentr) * (num_gas_species + 1) - 1];
    //cout << "M = " << data->M << "\n";
    //cout << "number = " << (Ncentr - 1) * (num_gas_species + 1) << "\n";
    int i_temp = 0;
    for (int i = 1; i < myNx - 1; i++) {

        i_T = i * (num_gas_species + 1) - 1;

        if (i == myNx - 2) {
            T_curr = yval[i_T];
            T_prev = yval[i_T - (num_gas_species + 1)];
            T_next = T_curr;

        }
        else if (i == 1) {
            T_prev = data->Tl;
            T_curr = yval[i_T];
            T_next = yval[i_T + (num_gas_species + 1)];
        }
        else {
            T_curr = yval[i_T];
            T_prev = yval[i_T - (num_gas_species + 1)];
            T_next = yval[i_T + (num_gas_species + 1)];
        };
        if (i == Ncentr) T_curr = data->T_center;
        if (i == Ncentr - 1) T_next = data->T_center;
        if (i == Ncentr + 1) T_prev = data->T_center;

        //cout << "in func i = " << i << "\n";
        //cout << "T_prev" << T_prev << "\n";
        //cout << "T_curr" << T_curr << "\n";
        //cout << "T_next" << T_next << "\n\n";
        MakeYvectors(data, yval, myNx, i, data->Tl);
        //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        //    cout << "Yprev = " << Yiprev[k_spec] << "\n";
        //}
        //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        //    cout << "Yi = " << Yi[k_spec] << "\n";
        //}
        //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        //    cout << "Ynext = " << Yinext[k_spec] << "\n";
        //}
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            Yi[k_spec] -= eps;
            if (i > 1) Yiprev[k_spec] -= eps;
            Yinext[k_spec] -= eps;
            //cout << "Yi = " << Yi[k_spec] << "\n";
        }
        Get_molar_cons(Xi, Yi, T_curr);

        chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr,
            T_curr, Xi, ydot);


        Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);

        find_diff_slag(data, T_curr, T_next, Yi, Yinext, Xi, Xinext, YkVk_r, Y_tmp_r, X_tmp_r, gradX_r, data->rho_r, data->Vc_r, i);

        find_diff_slag(data, T_prev, T_curr, Yiprev, Yi, Xiprev, Xi, YkVk_l, Y_tmp_l, X_tmp_l, gradX_l, data->rho_l, data->Vc_l, i - 1);

        rho = get_rho(Yi, T_curr);
        sumY = 0;
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            rval[i_temp] = F_rightY(data, k_spec,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]
            );
            //cout << "i_temp = " << i_temp << " = " << yval[i_temp] << "\n";
            i_temp++;
        }
        Cp = get_Cp(num_gas_species, Yi, T_curr);

        //cout << "i_temp = " << i_temp << " = " << yval[i_temp] << "\n\n";
            rval[i_temp] = F_right(data,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]);
        i_temp++;

    }
    return 0;
}


void MakeYvectors_dense(UserData data,
    double* Y, int myNx, int i, double Tl) {
    //cout << "i = " << i << "\n";
    for (int j = 0; j < num_gas_species; j++) {

        Yi[j] = Y[j + (i - 1) * num_gas_species];
        if (i == 1) Yiprev[j] = Y_left_bound[j];
        else Yiprev[j] = Y[j + (i - 2) * num_gas_species];

        if (i == myNx - 2)  Yinext[j] = Yi[j];
        else  Yinext[j] = Y[j + (i)*num_gas_species];
    }
}



int integrate_All_IDA_M(int N_x, vector<double>& x_vect,
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
    std::cout << "M_in_IDA = " << M << "\n";

    int i_temp = 0;
    for (int i = 1; i < N_x - 1; i++) {
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            consval[i_temp] = 0.0;   /*constraint*/
            yval[i_temp] = Y_vect[i * num_gas_species + k_spec] + eps;
            atval[i_temp] = pow(10, -5);
            //cout << "Yi = " << i_temp << "   = " << yval[i_temp] << endl;
            i_temp++;
        }
        yval[i_temp] = T_vect[i];
        consval[i_temp] = 1.0;
        atval[i_temp] = pow(10, 0);
        //cout << "Ti = " << i_temp << "   = " << yval[i_temp] << endl << endl;
        i_temp++;
        yval[i_temp] = M;
        atval[i_temp] = pow(10, -3);
        consval[i_temp] = 1.0;
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
    int i_T;
    double T_prev, T_curr, T_next;
    for (int i = 1; i < data->Nx - 1; i++) {

        i_T = i * (num_gas_species + 2) - 2;
        data->M = yval[i * (num_gas_species + 2) - 1];
        //cout << "in func i = " << i << "\n";
        //cout << "M = " << data->M << "\n\n";
        if (i == N_x - 2) {
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
        };

        //cout << "T_prev" << T_prev << "\n";
        //cout << "T_curr" << T_curr << "\n";
        //cout << "T_next" << T_next << "\n\n";
        MakeYvectors_kins(data, yval, N_x, i, data->Tl);
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            Yi[k_spec] -= eps;
            Yiprev[k_spec] -= eps;
            Yinext[k_spec] -= eps;
            //cout << "Yi = " << Yi[k_spec] << "\n";
        }
        //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        //    cout << "Yi = " << Yi[k_spec] << "\n";
        //}
        Get_molar_cons(Xi, Yi, T_curr);

        chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr,
            T_curr, Xi, ydot);


        Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);

        find_diff_slag(data, T_curr, T_next, Yi, Yinext, Xi, Xinext, YkVk_r, Y_tmp_r, X_tmp_r, gradX_r, data->rho_r, data->Vc_r, i);

        find_diff_slag(data, T_prev, T_curr, Yiprev, Yi, Xiprev, Xi, YkVk_l, Y_tmp_l, X_tmp_l, gradX_l, data->rho_l, data->Vc_l, i - 1);

        rho = get_rho(Yi, T_curr);
        sumY = 0;
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            ypval[i_temp] = -F_rightY(data, k_spec,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]
            ) / rho;
            //cout << "i_temp = " << i_temp << " = " << ypval[i_temp] << "\n";
            i_temp++;
        }
        double Cp = get_Cp(num_gas_species, Yi, T_curr);

        //cout << "i_temp = " << i_temp << " = " << yval[i_temp] << "\n\n";
        if (i == N_center)
        {
            ypval[i_temp] = 0;
        }
        else {
            ypval[i_temp] = -F_right(data,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]) / rho / Cp;
        }
        //cout << "i_temp = " << i_temp << " = " << ypval[i_temp] << "\n";
        i_temp++;
        if (i == N_center) {
            ypval[i_temp] = 0;
        }
        if (i > N_center) {
            ypval[i_temp] = 0;
        }
        if (i < N_center) {
            ypval[i_temp] = 0;
        }
        //cout << "i_temp = " << i_temp << " = " << ypval[i_temp] << "\n";
        i_temp++;
    }
    //cout << NEQ << " = " << Ith(res_vect, NEQ + 1) << endl;

    /* Call IDACreate and IDAInit to initialize IDA memory */    
    mem = IDACreate(ctx);
    if (check_retval((void*)mem, "IDACreate", 0)) return(1);

    retval = IDAInit(mem, func_All_IDA_M, t0, yy, yp);
    if (check_retval(&retval, "IDAInit", 1)) return(1);
    /* Call IDASVtolerances to set tolerances */

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

    IDASetMaxOrd(mem, 1);
    IDASetInitStep(mem, pow(10, -7));
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
        M = yval[i_temp];
        //cout << "i = " << i << "M = " << yval[i_temp] << "\n";
        i_temp++;
    }
    for (int i = 0; i < num_gas_species; i++) {
        Y_vect[(num_gas_species) * (N_x - 1) + i] = Y_vect[(num_gas_species) * (N_x - 2) + i];
    }

    T_vect[N_x - 1] = T_vect[N_x - 2];

    std::cout << "tout = " << tout << "\n";
    std::cout << "M = " << M << "\n";
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
    int i_T;
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
    for (int i = 1; i < myNx - 1; i++) {
        data->M = yval[i * (num_gas_species + 2) - 1];
        i_T = i * (num_gas_species + 2) - 2;
        //cout << "in func i = " << i << "\n";
        //cout << "M = " << data->M << "\n\n";
        if (i == myNx - 2) {
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
        };

        //cout << "T_prev" << T_prev << "\n";
        //cout << "T_curr" << T_curr << "\n";
        //cout << "T_next" << T_next << "\n\n";
        MakeYvectors_kins(data, yval, myNx, i, data->Tl);
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
            T_curr, Xi, ydot);


        Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);

        find_diff_slag(data, T_curr, T_next, Yi, Yinext, Xi, Xinext, YkVk_r, Y_tmp_r, X_tmp_r, gradX_r, data->rho_r, data->Vc_r, i);

        find_diff_slag(data, T_prev, T_curr, Yiprev, Yi, Xiprev, Xi, YkVk_l, Y_tmp_l, X_tmp_l, gradX_l, data->rho_l, data->Vc_l, i - 1);

        rho = get_rho(Yi, T_curr);
        sumY = 0;
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            rval[i_temp] = rho * ypval[i_temp] + F_rightY(data, k_spec,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]
            );
            //cout << "i_temp = " << i_temp << " = " << rval[i_temp] << "\n";
            i_temp++;
        }
        Cp = get_Cp(num_gas_species, Yi, T_curr);

        if (i != Ncentr) {
            rval[i_temp] = rho * Cp * ypval[i_temp] + F_right(data,
                    T_prev, T_curr, T_next,
                    x_vect[i - 1], x_vect[i], x_vect[i + 1]);
        }
        else {
            rval[i_temp] = yval[i_temp] - data->T_center;
        }
       
        //cout << "i_temp = " << i_temp << " = " << rval[i_temp] << "\n\n";
        i_temp++;

        if (i == Ncentr) {
            rval[i_temp] = F_right(data,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]);
        }
        if (i > Ncentr) {
            rval[i_temp] = -(yval[i_temp] - yval[i_temp - (num_gas_species + 2)]) / (x_vect[i] - x_vect[i - 1]);
        }
        if (i < Ncentr) {
            rval[i_temp] = -(yval[i_temp + (num_gas_species + 2)] - yval[i_temp]) / (x_vect[i + 1] - x_vect[i]);
        }
        //cout << "i_temp = " << i_temp << " = " << rval[i_temp] << "\n\n";
        i_temp++;
    }
    return 0;
}

void set_Dij_res(double T) {
    for (unsigned short int i = 0; i < num_gas_species; ++i) {
        for (unsigned short int j = i; j < num_gas_species; ++j) {
            Dij_res[i][j] = Dij_func5(i, j, T);
            Dij_res[j][i] = Dij_res[i][j];
        }
    }
}



int integrate_All_IDA_dense(int N_x, vector<double>& x_vect,
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
        consval[i] = 0.0;   /*constraint*/
        yval[i] = Y_vect[i + num_gas_species];
        ypval[i] = 0;
        number_spec = i % num_gas_species;
        atval[i] = RCONST(1.0e-5);
        //cout << "i = " << i << "   = " << yval[i] << endl;
        //if ((i + 1) % num_gas_species == 0) cout << endl;
    }
    consval[1] = 1.;

    MakeYvectors(data,
        yval, N_x, N_center, T_vect[N_center]);

    Get_molar_cons(Xi, Yi, T_vect[N_center]);
    chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr, T_vect[N_center], Xi, ydot);

    Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);
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

        Get_molar_cons(Xi, Yi, T_vect[i]);
        chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr, T_vect[i], Xi, ydot);

        Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);
        rho = get_rho(Yi, T_vect[i]);

        double Cp = get_Cp(num_gas_species, Yi, T_vect[i]);
        find_diff_slag(data, T_vect[i], T_vect[i + 1], Yi, Yinext, Xi, Xinext, YkVk_r, Y_tmp_r, X_tmp_r, gradX_r, data->rho_r, data->Vc_r, i);

        find_diff_slag(data, T_vect[i - 1], T_vect[i], Yiprev, Yi, Xiprev, Xi, YkVk_l, Y_tmp_l, X_tmp_l, gradX_l, data->rho_l, data->Vc_l, i - 1);
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            ypval[k_spec + (i - 1) * num_gas_species] = -F_rightY(data, k_spec,
                T_vect[i - 1], T_vect[i], T_vect[i + 1],
                x_vect[i - 1], x_vect[i], x_vect[i + 1]
            ) / rho;
        }
        if (i != N_center) {
            ypval[i + NEQ_Y - 1] = -F_right(data,
                T_vect[i - 1], T_vect[i], T_vect[i + 1],
                x_vect[i - 1], x_vect[i], x_vect[i + 1]) / rho / Cp;
        }
        else {
            ypval[i + NEQ_Y - 1] = -F_right(data,
                T_vect[i - 1], T_vect[i], T_vect[i + 1],
                x_vect[i - 1], x_vect[i], x_vect[i + 1]);
        }

    }
    //cout << NEQ << " = " << Ith(res_vect, NEQ + 1) << endl;

    /* Call IDACreate and IDAInit to initialize IDA memory */
    mem = IDACreate(ctx);
    if (check_retval((void*)mem, "IDACreate", 0)) return(1);

    retval = IDAInit(mem, func_All_IDA_dense, t0, yy, yp);
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
    int number = 1;
    ofstream fout;
    ofstream foutw;
    double Y_H2, Y_O2;
    double W, w_dot;
    double sum_Y = 0;
    double alpha = 1;

    for (int i = 0; i < num_gas_species * N_x; i++)
        Yp_vect[i] = 0;

    data->N_m = 0;
    while (iout < ida_steps) {
        retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);

        if (check_retval(&retval, "IDASolve", 1)) return(1);
        iout++;
        tout1 *= 1.0;
        cout << "tout = " << tout << "\n";
        cout << "M = " << M << "\n";
        tout += tout1;
    }

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

    for (int i = 0; i < N_x; i++) {
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            Yi[k_spec] = Y_vect[k_spec + i * num_gas_species];
        }
       
    }
    cout << "tout = " << tout << "\n";
    cout << "M = " << M << "\n";
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


static int func_All_IDA_dense(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void* user_data)
{
    realtype* yval, * ypval, * rval;
    UserData data;
    double Temp;
    double M;
    double sumY = 0;
    int i_T;
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

    data->M = yval[data->N_centr + data->NEQ_Y - 1];
    for (int i = 1; i < myNx - 1; i++) {

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
        MakeYvectors_dense(data, yval, myNx, i, data->Tl);

        Get_molar_cons(Xi, Yi, T_curr);
        chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr,
            T_curr,
            Xi, ydot);


        Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);

        find_diff_slag(data, T_curr, T_next, Yi, Yinext, Xi, Xinext, YkVk_r, Y_tmp_r, X_tmp_r, gradX_r, data->rho_r, data->Vc_r, i);

        find_diff_slag(data, T_prev, T_curr, Yiprev, Yi, Xiprev, Xi, YkVk_l, Y_tmp_l, X_tmp_l, gradX_l, data->rho_l, data->Vc_l, i - 1);

        rho = get_rho(Yi, T_curr);
        sumY = 0;
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            rval[k_spec + (i - 1) * num_gas_species] = rho * ypval[k_spec + (i - 1) * num_gas_species] + F_rightY(data, k_spec,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]);

        }
        Cp = get_Cp(num_gas_species, Yi, T_curr);
        if (i != data->N_centr)
            rval[i_T] = rho * Cp * ypval[i_T] + F_right(data,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]
            );
        else
            rval[i_T] = F_right(data,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]);
    }
    return(0);
}




int Integrate_Kinsol_withoutM(int N_x, vector<double>& x_vect,
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
    data->M = M;
    std::cout << "INTEGRATE KINSOL = " << M << "\n";
    int i_temp = 0;
    for (int i = 1; i < N_x - 1; i++) {
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            cval[i_temp] = 1.0;   /*constraint*/
            yval[i_temp] = Y_vect[i * num_gas_species + k_spec] + eps;
            sval[i_temp] = pow(10, 0);
            //cout << "Yi = " << i_temp << "   = " << yval[i_temp] << endl;
            i_temp++;
        }
        yval[i_temp] = T_vect[i];
        cval[i_temp] = 1.0;
        sval[i_temp] = pow(10, -5);
        //cout << "Ti = " << i_temp << "   = " << yval[i_temp] << endl << endl;
        i_temp++;
        //if ((i + 1) % num_gas_species == 0) cout << endl;
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

    retval = KINInit(kmem, func_kinsol_withoutM, res_vect);
    if (check_retval(&retval, "KINInit", 1)) return(1);

    /* Create dense SUNMatrix */
    int mu = 2 * (num_gas_species + 1);
    int ml = 2 * (num_gas_species + 1);

    J = SUNBandMatrix(NEQ, mu, ml, sunctx);
    if (check_retval((void*)J, "SUNBandMatrix", 0)) return(1);

    /* Create banded SUNLinearSolver object */
    LS = SUNLinSol_Band(res_vect, J, sunctx);
    if (check_retval((void*)LS, "SUNLinSol_Band", 0)) return(1);

    /* Attach the matrix and linear solver to KINSOL */
    retval = KINSetLinearSolver(kmem, LS, J);
    if (check_retval(&retval, "KINSetLinearSolver", 1)) return(1);

    glstr = 0;
    mset = 1000;

    retval = KINSol(kmem, res_vect, glstr, s, s);
    if (check_retval(&retval, "KINSol", 1)) return(1);
    cout << "retval = " << retval << "\n";
    i_temp = 0;
    for (int i = 1; i < N_x - 1; i++) {
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            Y_vect[i * num_gas_species + k_spec] = yval[i_temp] - eps;
            i_temp++;
        }
        T_vect[i] = yval[i_temp];
        i_temp++;
    }

    for (int i = 0; i < num_gas_species; i++) {
        Y_vect[(num_gas_species) * (N_x - 1) + i] = Y_vect[(num_gas_species) * (N_x - 2) + i];
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

static int func_kinsol_withoutM(N_Vector u, N_Vector f, void* user_data)
{
    realtype* yval, * ypval, * rval;
    UserData data;
    double Temp;
    double M;
    double sumY = 0;
    int i_T;
    data = (UserData)user_data;
    double Cp;
    double T_curr, T_prev, T_next;
    int j;
    int myNx = data->Nx;
    int myNeq = data->NEQ;
    int myNm = data->N_m;
    yval = N_VGetArrayPointer(u);
    rval = N_VGetArrayPointer(f);
    int Ncentr = data->N_centr;
    double rho;
    int i_temp = 0;
    for (int i = 1; i < myNx - 1; i++) {

        i_T = i * (num_gas_species + 1) - 1;

        if (i == myNx - 2) {
            T_curr = yval[i_T];
            T_prev = yval[i_T - (num_gas_species + 1)];
            T_next = T_curr;

        }
        else if (i == 1) {
            T_prev = data->Tl;
            T_curr = yval[i_T];
            T_next = yval[i_T + (num_gas_species + 1)];
        }
        else {
            T_curr = yval[i_T];
            T_prev = yval[i_T - (num_gas_species + 1)];
            T_next = yval[i_T + (num_gas_species + 1)];
        };

        //cout << "in func i = " << i << "\n";
        //cout << "T_prev" << T_prev << "\n";
        //cout << "T_curr" << T_curr << "\n";
        //cout << "T_next" << T_next << "\n\n";
        MakeYvectors(data, yval, myNx, i, data->Tl);
        // for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        //    cout << "Yprev = " << Yiprev[k_spec] << "\n";
        //}
        //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        //    cout << "Yi = " << Yi[k_spec] << "\n";
        //}
        //for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
        //    cout << "Ynext = " << Yinext[k_spec] << "\n";
        //}
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            Yi[k_spec] -= eps;
            if (i > 1) Yiprev[k_spec] -= eps;
            Yinext[k_spec] -= eps;
            //cout << "Yi = " << Yi[k_spec] << "\n";
        }

        Get_molar_cons(Xi, Yi, T_curr);

        chem_vel(Sn, Hn, forward_arr, reverse_arr, equilib_arr,
            T_curr, Xi, ydot);


        Get_mole_fr(Xiprev, Yiprev); Get_mole_fr(Xi, Yi); Get_mole_fr(Xinext, Yinext);

        find_diff_slag(data, T_curr, T_next, Yi, Yinext, Xi, Xinext, YkVk_r, Y_tmp_r, X_tmp_r, gradX_r, data->rho_r, data->Vc_r, i);

        find_diff_slag(data, T_prev, T_curr, Yiprev, Yi, Xiprev, Xi, YkVk_l, Y_tmp_l, X_tmp_l, gradX_l, data->rho_l, data->Vc_l, i - 1);

        rho = get_rho(Yi, T_curr);
        sumY = 0;
        for (int k_spec = 0; k_spec < num_gas_species; k_spec++) {
            rval[i_temp] = F_rightY(data, k_spec,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]
            );
            //cout << "i_temp = " << i_temp << " = " << yval[i_temp] << "\n";
            i_temp++;
        }
        Cp = get_Cp(num_gas_species, Yi, T_curr);

        //cout << "i_temp = " << i_temp << " = " << yval[i_temp] << "\n\n";
        if (i == Ncentr)
        {
            rval[i_temp] = yval[i_temp] - data->T_center;
        }
        else {
            rval[i_temp] = F_right(data,
                T_prev, T_curr, T_next,
                x_vect[i - 1], x_vect[i], x_vect[i + 1]);
        }

        i_temp++;

    }
    return 0;
}