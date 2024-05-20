//! @file GasTransport.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "GasTransport.h"
#include "MMCollisionInt.h"
#include "polyfit.h"


//! polynomial degree used for fitting collision integrals
//! except in CK mode, where the degree is 6.
#define COLL_INT_POLY_DEGREE 8
#define minTemp 300
#define maxTemp 3000
#define Press_diff 101325


void resize2D(vector<vector<double>>& vect,int dim1, int dim2) {
    vect.resize(dim1);
    for (int i = 0; i < dim1; i++) {
        vect[i].resize(dim2);
        for (int j = 0; j < dim2; j++) {
            vect[i][j] = 0.;
        }
    }
}


void GasTransport::init(int mode, int m_nsp)
{
    m_mode = mode;

    // set up Monchick and Mason collision integrals
    setupCollisionParameters();
    setupCollisionIntegral();

    m_molefracs.resize(num_gas_species);
    m_spwork.resize(num_gas_species);
    m_visc.resize(num_gas_species);
    m_sqvisc.resize(num_gas_species);
    resize2D(m_phi, num_gas_species, num_gas_species);
    resize2D(m_bdiff, num_gas_species, num_gas_species);


    resize2D(m_wratjk, num_gas_species, num_gas_species);
    resize2D(m_wratkj1, num_gas_species, num_gas_species);
    for (size_t j = 0; j < num_gas_species; j++) {
        for (size_t k = j; k < num_gas_species; k++) {
            m_wratjk[j][k] = sqrt(phyc.mol_weight[j] / phyc.mol_weight[k]);
            m_wratjk[k][j] = sqrt(m_wratjk[j][k]);
            m_wratkj1[j][k] = sqrt(1.0 + phyc.mol_weight[k] / phyc.mol_weight[j]);
        }
    }
}

void GasTransport::setupCollisionParameters()
{


    resize2D(m_epsilon, num_gas_species, num_gas_species);
    resize2D(m_delta, num_gas_species, num_gas_species);
    resize2D(m_reducedMass, num_gas_species, num_gas_species);
    resize2D(m_dipole, num_gas_species, num_gas_species);
    resize2D(m_diam, num_gas_species, num_gas_species);

    m_crot.resize(num_gas_species);
    m_zrot.resize(num_gas_species);
    m_polar.resize(num_gas_species, false);
    m_alpha.resize(num_gas_species, 0.0);
    m_poly.resize(num_gas_species);
    m_star_poly_uses_actualT.resize(num_gas_species);
    m_sigma.resize(num_gas_species);
    m_eps.resize(num_gas_species);
    m_w_ac.resize(num_gas_species);
    m_disp.resize(num_gas_species, 0.0);
    m_quad_polar.resize(num_gas_species, 0.0);
    getTransportData();

    for (size_t i = 0; i < num_gas_species; i++) {
        m_poly[i].resize(num_gas_species);
        m_star_poly_uses_actualT[i].resize(num_gas_species);
    }

    double f_eps, f_sigma;

    for (size_t i = 0; i < num_gas_species; i++) {
        for (size_t j = i; j < num_gas_species; j++) {
            // the reduced mass
            m_reducedMass[i][j] = phyc.mol_weight[i] * phyc.mol_weight[j] / (Avogadro * (phyc.mol_weight[i] + phyc.mol_weight[j]));

            // hard-sphere diameter for (i,j) collisions
            m_diam[i][j] = 0.5 * (m_sigma[i] + m_sigma[j]);

            // the effective well depth for (i,j) collisions
            m_epsilon[i][j] = sqrt(m_eps[i] * m_eps[j]);

            // the effective dipole moment for (i,j) collisions
            m_dipole[i][j] = sqrt(m_dipole[i][i] * m_dipole[j][j]);

            // reduced dipole moment delta* (nondimensional)
            double d = m_diam[i][j];
            m_delta[i][j] = 0.5 * m_dipole[i][j] * m_dipole[i][j]
                / (4 * Pi * epsilon_0 * m_epsilon[i][j] * d * d * d);
            makePolarCorrections(i, j, f_eps, f_sigma);
            m_diam[i][j] *= f_sigma;
            m_epsilon[i][j] *= f_eps;

            // properties are symmetric
            m_reducedMass[j][i] = m_reducedMass[i][j];
            m_diam[j][i] = m_diam[i][j];
            m_epsilon[j][i] = m_epsilon[i][j];
            m_dipole[j][i] = m_dipole[i][j];
            m_delta[j][i] = m_delta[i][j];
        }
    }
}

void GasTransport::setupCollisionIntegral()
{
    double tstar_min = 1.e8, tstar_max = 0.0;
    for (size_t i = 0; i < num_gas_species; i++) {
        for (size_t j = i; j < num_gas_species; j++) {
            // The polynomial fits of collision integrals vs. T*
            // will be done for the T* from tstar_min to tstar_max
            tstar_min = std::min(tstar_min, Boltzmann * minTemp/ m_epsilon[i][j]);
            tstar_max = std::max(tstar_max, Boltzmann * maxTemp / m_epsilon[i][j]);
        }
    }
    // Chemkin fits the entire T* range in the Monchick and Mason tables,
    // so modify tstar_min and tstar_max if in Chemkin compatibility mode
    if (m_mode == CK_Mode) {
        tstar_min = 0.101;
        tstar_max = 99.9;
    }

    // initialize the collision integral calculator for the desired T* range
    MMCollisionInt integrals;
    integrals.init(tstar_min, tstar_max, m_log_level);
    fitCollisionIntegrals(integrals);
    // make polynomial fits
    fitProperties(integrals);
}

void GasTransport::getTransportData()
{
    auto& Species = chec.chemkinReader->species();
    for (size_t k = 0; k < num_gas_species; k++) {

        if (Species[k].transport().getMoleculeIndex() == 0) {
            m_crot[k] = 0.0;
        }
        else if (Species[k].transport().getMoleculeIndex() == 1) {
            m_crot[k] = 1.0;
        }
        else if (Species[k].transport().getMoleculeIndex() == 2) {
            m_crot[k] = 1.5;
        }

        m_sigma[k] = Species[k].transport().getCollisionDiameter();
        m_eps[k] = Species[k].transport().getPotentialWellDepth();
        m_dipole[k][k] = Species[k].transport().getDipoleMoment();
        m_polar[k] = (Species[k].transport().getDipoleMoment() > 0);
        m_alpha[k] = Species[k].transport().getPolarizability();
        m_zrot[k] = Species[k].transport().getRotRelaxationNumber();

        //cout << "\nmol_weight = " << phyc.mol_weight[k] << "\n";
        //cout << "m_sigma[k] = " << m_sigma[k] << "\n";
        //cout << "m_eps[k] = " << m_eps[k] << "\n";
        //cout << "m_dipole[k][k] = " << m_dipole[k][k] << "\n";
        //cout << "m_polar[k] = " << m_polar[k] << "\n";
        //cout << "m_alpha[k] = " << m_alpha[k] << "\n";
        //cout << "m_zrot[k] = " << m_zrot[k] << "\n\n";
    }
}

void GasTransport::makePolarCorrections(size_t i, size_t j,
    double& f_eps, double& f_sigma)
{
    // no correction if both are nonpolar, or both are polar
    if (m_polar[i] == m_polar[j]) {
        f_eps = 1.0;
        f_sigma = 1.0;
        return;
    }

    // corrections to the effective diameter and well depth
    // if one is polar and one is non-polar
    size_t kp = (m_polar[i] ? i : j); // the polar one
    size_t knp = (i == kp ? j : i); // the nonpolar one
    double d3np, d3p, alpha_star, mu_p_star, xi;
    d3np = pow(m_sigma[knp], 3);
    d3p = pow(m_sigma[kp], 3);
    alpha_star = m_alpha[knp] / d3np;
    mu_p_star = m_dipole[kp][kp] / sqrt(4 * Pi * epsilon_0 * d3p * m_eps[kp]);
    xi = 1.0 + 0.25 * alpha_star * mu_p_star * mu_p_star *
        sqrt(m_eps[kp] / m_eps[knp]);
    f_sigma = pow(xi, -1.0 / 6.0);
    f_eps = xi * xi;
}

void GasTransport::fitCollisionIntegrals(MMCollisionInt& integrals)
{
    // Chemkin fits to sixth order polynomials
    int degree = (m_mode == CK_Mode ? 6 : COLL_INT_POLY_DEGREE);
    vector<double> fitlist;
    m_omega22_poly.clear();
    m_astar_poly.clear();
    m_bstar_poly.clear();
    m_cstar_poly.clear();
    for (size_t i = 0; i < num_gas_species; i++) {
        for (size_t j = i; j < num_gas_species; j++) {
            // Chemkin fits only delta* = 0
            double dstar = (m_mode != CK_Mode) ? m_delta[i][j] : 0.0;

            // if a fit has already been generated for delta* = m_delta(i,j),
            // then use it. Otherwise, make a new fit, and add m_delta(i,j) to
            // the list of delta* values for which fits have been done.

            // 'find' returns a pointer to end() if not found
            auto dptr = find(fitlist.begin(), fitlist.end(), dstar);
            if (dptr == fitlist.end()) {
                vector<double> ca(degree + 1), cb(degree + 1), cc(degree + 1);
                vector<double> co22(degree + 1);
                integrals.fit(degree, dstar, ca.data(), cb.data(), cc.data());
                integrals.fit_omega22(degree, dstar, co22.data());
                m_omega22_poly.push_back(co22);
                m_astar_poly.push_back(ca);
                m_bstar_poly.push_back(cb);
                m_cstar_poly.push_back(cc);
                m_poly[i][j] = static_cast<int>(m_astar_poly.size()) - 1;
                m_star_poly_uses_actualT[i][j] = 0;
                fitlist.push_back(dstar);
            }
            else {
                // delta* found in fitlist, so just point to this polynomial
                m_poly[i][j] = static_cast<int>((dptr - fitlist.begin()));
                m_star_poly_uses_actualT[i][j] = 0;
            }
            m_poly[j][i] = m_poly[i][j];
            m_star_poly_uses_actualT[j][i] = m_star_poly_uses_actualT[i][j];
        }
    }
}

void GasTransport::fitProperties(MMCollisionInt& integrals)
{
    // number of points to use in generating fit data
    const size_t np = 400;
    //ATTENTION
    m_mode = CK_Mode;
    //ATTENTION
    int degree = (m_mode == CK_Mode ? 3 : 4);
    double dt = (maxTemp - minTemp) / (np - 1);
    vector<double> tlog(np), spvisc(np), spcond(np);
    vector<double> w(np), w2(np);

    m_visccoeffs.clear();
    m_condcoeffs.clear();

    // generate array of log(t) values
    for (size_t n = 0; n < np; n++) {
        double t = minTemp + dt * n;
        tlog[n] = log(t);
    }

    // vector of polynomial coefficients
    vector<double> c(degree + 1), c2(degree + 1);

    // fit the pure-species viscosity and thermal conductivity for each species
    double visc, err, relerr,
        mxerr = 0.0, mxrelerr = 0.0, mxerr_cond = 0.0, mxrelerr_cond = 0.0;


    for (size_t k = 0; k < num_gas_species; k++) {
        double tstar = Boltzmann * 298.0 / m_eps[k];
        // Scaling factor for temperature dependence of z_rot. [Kee2003] Eq.
        // 12.112 or [Kee2017] Eq. 11.115
        double fz_298 = 1.0 + pow(Pi, 1.5) / sqrt(tstar) * (0.5 + 1.0 / tstar) +
            (0.25 * Pi * Pi + 2) / tstar;

        for (size_t n = 0; n < np; n++) {
            double t = minTemp + dt * n;
            //vector<double> cp_R_all(m_thermo->nSpecies());
            //m_thermo->getCp_R_ref(&cp_R_all[0]);
            //double cp_R = cp_R_all[k];
            double cp_R = get_Cpi(k, t, 0);
            //std::cout << "Cp_R = " << cp_R << "\n";
            tstar = Boltzmann * t / m_eps[k];
            double sqrt_T = sqrt(t);
            double om22 = integrals.omega22(tstar, m_delta[k][k]);
            double om11 = integrals.omega11(tstar, m_delta[k][k]);

            // self-diffusion coefficient, without polar corrections
            double diffcoeff = 3.0 / 16.0 * sqrt(2.0 * Pi / m_reducedMass[k][k]) *
                pow((Boltzmann * t), 1.5) /
                (Pi * m_sigma[k] * m_sigma[k] * om11);

            // viscosity
            visc = 5.0 / 16.0 * sqrt(Pi * phyc.mol_weight[k] * Boltzmann * t / Avogadro) /
                (om22 * Pi * m_sigma[k] * m_sigma[k]);

            // thermal conductivity
            double f_int = phyc.mol_weight[k] / (GasConstant * t) * diffcoeff / visc;
            double cv_rot = m_crot[k];
            double A_factor = 2.5 - f_int;
            double fz_tstar = 1.0 + pow(Pi, 1.5) / sqrt(tstar) * (0.5 + 1.0 / tstar) +
                (0.25 * Pi * Pi + 2) / tstar;
            double B_factor = m_zrot[k] * fz_298 / fz_tstar + 2.0 / Pi * (5.0 / 3.0 * cv_rot + f_int);
            double c1 = 2.0 / Pi * A_factor / B_factor;
            double cv_int = cp_R - 2.5 - cv_rot;
            double f_rot = f_int * (1.0 + c1);
            double f_trans = 2.5 * (1.0 - c1 * cv_rot / 1.5);
            double cond = (visc / phyc.mol_weight[k]) * GasConstant * (f_trans * 1.5
                + f_rot * cv_rot + f_int * cv_int);

            if (m_mode == CK_Mode) {
                spvisc[n] = log(visc);
                spcond[n] = log(cond);
                w[n] = -1.0;
                w2[n] = -1.0;
            }
            else {
                // the viscosity should be proportional approximately to
                // sqrt(T); therefore, visc/sqrt(T) should have only a weak
                // temperature dependence. And since the mixture rule requires
                // the square root of the pure-species viscosity, fit the square
                // root of (visc/sqrt(T)) to avoid having to compute square
                // roots in the mixture rule.
                spvisc[n] = sqrt(visc / sqrt_T);

                // the pure-species conductivity scales approximately with
                // sqrt(T). Unlike the viscosity, there is no reason here to fit
                // the square root, since a different mixture rule is used.
                spcond[n] = cond / sqrt_T;
                w[n] = 1.0 / (spvisc[n] * spvisc[n]);
                w2[n] = 1.0 / (spcond[n] * spcond[n]);
            }
        }
        polyfit(np, degree, tlog.data(), spvisc.data(), w.data(), c.data());
        polyfit(np, degree, tlog.data(), spcond.data(), w2.data(), c2.data());

        // evaluate max fit errors for viscosity
        for (size_t n = 0; n < np; n++) {
            double val, fit;
            if (m_mode == CK_Mode) {
                val = exp(spvisc[n]);
                fit = exp(poly3(tlog[n], c.data()));
            }
            else {
                double sqrt_T = exp(0.5 * tlog[n]);
                val = sqrt_T * pow(spvisc[n], 2);
                fit = sqrt_T * pow(poly4(tlog[n], c.data()), 2);
            }
            err = fit - val;
            relerr = err / val;
            mxerr = std::max(mxerr, fabs(err));
            mxrelerr = std::max(mxrelerr, fabs(relerr));
        }

        // evaluate max fit errors for conductivity
        for (size_t n = 0; n < np; n++) {
            double val, fit;
            if (m_mode == CK_Mode) {
                val = exp(spcond[n]);
                fit = exp(poly3(tlog[n], c2.data()));
            }
            else {
                double sqrt_T = exp(0.5 * tlog[n]);
                val = sqrt_T * spcond[n];
                fit = sqrt_T * poly4(tlog[n], c2.data());
            }
            err = fit - val;
            relerr = err / val;
            mxerr_cond = std::max(mxerr_cond, fabs(err));
            mxrelerr_cond = std::max(mxrelerr_cond, fabs(relerr));
        }
        m_visccoeffs.push_back(c);
        m_condcoeffs.push_back(c2);
    }

    fitDiffCoeffs(integrals);
}

void GasTransport::fitDiffCoeffs(MMCollisionInt& integrals)
{
    // number of points to use in generating fit data
    const size_t np = 400;
    int degree = (m_mode == CK_Mode ? 3 : 4);
    double dt = (maxTemp - minTemp) / (np - 1);
    vector<double> tlog(np);
    vector<double> w(np), w2(np);

    // generate array of log(t) values
    for (size_t n = 0; n < np; n++) {
        double t = minTemp + dt * n;
        tlog[n] = log(t);
    }

    // vector of polynomial coefficients
    vector<double> c(degree + 1), c2(degree + 1);
    double err, relerr,
        mxerr = 0.0, mxrelerr = 0.0;

    vector<double> diff(np + 1);
    m_diffcoeffs.clear();
    for (size_t k = 0; k < num_gas_species; k++) {
        for (size_t j = k; j < num_gas_species; j++) {
            for (size_t n = 0; n < np; n++) {
                double t = minTemp + dt * n;
                double eps = m_epsilon[j][k];
                double tstar = Boltzmann * t / eps;
                double sigma = m_diam[j][k];
                double om11 = integrals.omega11(tstar, m_delta[j][k]);
                double diffcoeff = 3.0 / 16.0 * sqrt(2.0 * Pi / m_reducedMass[k][j])
                    * pow(Boltzmann * t, 1.5) / (Pi * sigma * sigma * om11);
                // 2nd order correction
                // NOTE: THIS CORRECTION IS NOT APPLIED
                //double fkj, fjk;
                //getBinDiffCorrection(t, integrals, k, j, 1.0, 1.0, fkj, fjk);

                if (m_mode == CK_Mode) {
                    diff[n] = log(diffcoeff);
                    w[n] = -1.0;
                }
                else {
                    diff[n] = diffcoeff / pow(t, 1.5);
                    w[n] = 1.0 / (diff[n] * diff[n]);
                }
            }
            polyfit(np, degree, tlog.data(), diff.data(), w.data(), c.data());

            for (size_t n = 0; n < np; n++) {
                double val, fit;
                if (m_mode == CK_Mode) {
                    val = exp(diff[n]);
                    fit = exp(poly3(tlog[n], c.data()));
                }
                else {
                    double t = exp(tlog[n]);
                    double pre = pow(t, 1.5);
                    val = pre * diff[n];
                    fit = pre * poly4(tlog[n], c.data());
                }
                err = fit - val;
                relerr = err / val;
                mxerr = std::max(mxerr, fabs(err));
                mxrelerr = std::max(mxrelerr, fabs(relerr));
            }
            m_diffcoeffs.push_back(c);
        }
    }
}

void GasTransport::getBinDiffCorrection(double t, MMCollisionInt& integrals,
    size_t k, size_t j, double xk, double xj, double& fkj, double& fjk)
{
    double w1 = phyc.mol_weight[k];
    double w2 = phyc.mol_weight[j];
    double wsum = w1 + w2;
    double wmwp = (w1 - w2) / wsum;
    double sqw12 = sqrt(w1 * w2);
    double sig1 = m_sigma[k];
    double sig2 = m_sigma[j];
    double sig12 = 0.5 * (m_sigma[k] + m_sigma[j]);
    double sigratio = sig1 * sig1 / (sig2 * sig2);
    double sigratio2 = sig1 * sig1 / (sig12 * sig12);
    double sigratio3 = sig2 * sig2 / (sig12 * sig12);
    double tstar1 = Boltzmann * t / m_eps[k];
    double tstar2 = Boltzmann * t / m_eps[j];
    double tstar12 = Boltzmann * t / sqrt(m_eps[k] * m_eps[j]);
    double om22_1 = integrals.omega22(tstar1, m_delta[k][k]);
    double om22_2 = integrals.omega22(tstar2, m_delta[j][j]);
    double om11_12 = integrals.omega11(tstar12, m_delta[k][j]);
    double astar_12 = integrals.astar(tstar12, m_delta[k][j]);
    double bstar_12 = integrals.bstar(tstar12, m_delta[k][j]);
    double cstar_12 = integrals.cstar(tstar12, m_delta[k][j]);

    double cnst = sigratio * sqrt(2.0 * w2 / wsum) * 2.0 * w1 * w1 / (wsum * w2);
    double p1 = cnst * om22_1 / om11_12;

    cnst = (1.0 / sigratio) * sqrt(2.0 * w1 / wsum) * 2.0 * w2 * w2 / (wsum * w1);
    double p2 = cnst * om22_2 / om11_12;
    double p12 = 15.0 * wmwp * wmwp + 8.0 * w1 * w2 * astar_12 / (wsum * wsum);

    cnst = (2.0 / (w2 * wsum)) * sqrt(2.0 * w2 / wsum) * sigratio2;
    double q1 = cnst * ((2.5 - 1.2 * bstar_12) * w1 * w1 + 3.0 * w2 * w2
        + 1.6 * w1 * w2 * astar_12);

    cnst = (2.0 / (w1 * wsum)) * sqrt(2.0 * w1 / wsum) * sigratio3;
    double q2 = cnst * ((2.5 - 1.2 * bstar_12) * w2 * w2 + 3.0 * w1 * w1
        + 1.6 * w1 * w2 * astar_12);
    double q12 = wmwp * wmwp * 15.0 * (2.5 - 1.2 * bstar_12)
        + 4.0 * w1 * w2 * astar_12 * (11.0 - 2.4 * bstar_12) / (wsum * wsum)
        + 1.6 * wsum * om22_1 * om22_2 / (om11_12 * om11_12 * sqw12)
        * sigratio2 * sigratio3;

    cnst = 6.0 * cstar_12 - 5.0;
    fkj = 1.0 + 0.1 * cnst * cnst *
        (p1 * xk * xk + p2 * xj * xj + p12 * xk * xj) /
        (q1 * xk * xk + q2 * xj * xj + q12 * xk * xj);
    fjk = 1.0 + 0.1 * cnst * cnst *
        (p2 * xk * xk + p1 * xj * xj + p12 * xk * xj) /
        (q2 * xk * xk + q1 * xj * xj + q12 * xk * xj);
}

void GasTransport::getViscosityPolynomial(size_t i, double* coeffs) const
{
    for (int k = 0; k < (m_mode == CK_Mode ? 4 : 5); k++) {
        coeffs[k] = m_visccoeffs[i][k];
    }
}

void GasTransport::getConductivityPolynomial(size_t i, double* coeffs) const
{
    for (int k = 0; k < (m_mode == CK_Mode ? 4 : 5); k++) {
        coeffs[k] = m_condcoeffs[i][k];
    }
}

void GasTransport::getBinDiffusivityPolynomial(size_t i, size_t j, double* coeffs) const
{
    size_t mi = (j >= i ? i : j);
    size_t mj = (j >= i ? j : i);
    size_t ic = 0;
    for (size_t ii = 0; ii < mi; ii++) {
        ic += num_gas_species - ii;
    }
    ic += mj - mi;

    for (int k = 0; k < (m_mode == CK_Mode ? 4 : 5); k++) {
        coeffs[k] = m_diffcoeffs[ic][k];
    }
}

void GasTransport::getCollisionIntegralPolynomial(size_t i, size_t j,
    double* astar_coeffs,
    double* bstar_coeffs,
    double* cstar_coeffs) const
{
    for (int k = 0; k < (m_mode == CK_Mode ? 6 : COLL_INT_POLY_DEGREE) + 1; k++) {
        astar_coeffs[k] = m_astar_poly[m_poly[i][j]][k];
        bstar_coeffs[k] = m_bstar_poly[m_poly[i][j]][k];
        cstar_coeffs[k] = m_cstar_poly[m_poly[i][j]][k];
    }
}

void GasTransport::setViscosityPolynomial(size_t i, double* coeffs)
{
    for (int k = 0; k < (m_mode == CK_Mode ? 4 : 5); k++) {
        m_visccoeffs[i][k] = coeffs[k];
    }

    m_visc_ok = false;
    m_spvisc_ok = false;
    m_viscwt_ok = false;
    m_bindiff_ok = false;
    m_temp = -1;
}

void GasTransport::setConductivityPolynomial(size_t i, double* coeffs)
{
    for (int k = 0; k < (m_mode == CK_Mode ? 4 : 5); k++) {
        m_condcoeffs[i][k] = coeffs[k];
    }

    m_visc_ok = false;
    m_spvisc_ok = false;
    m_viscwt_ok = false;
    m_bindiff_ok = false;
    m_temp = -1;
}

void GasTransport::setBinDiffusivityPolynomial(size_t i, size_t j, double* coeffs)
{
    size_t mi = (j >= i ? i : j);
    size_t mj = (j >= i ? j : i);
    size_t ic = 0;
    for (size_t ii = 0; ii < mi; ii++) {
        ic += num_gas_species - ii;
    }
    ic += mj - mi;

    for (int k = 0; k < (m_mode == CK_Mode ? 4 : 5); k++) {
        m_diffcoeffs[ic][k] = coeffs[k];
    }

    m_visc_ok = false;
    m_spvisc_ok = false;
    m_viscwt_ok = false;
    m_bindiff_ok = false;
    m_temp = -1;
}

void GasTransport::setCollisionIntegralPolynomial(size_t i, size_t j,
    double* astar_coeffs,
    double* bstar_coeffs,
    double* cstar_coeffs, bool actualT)
{
    size_t degree = (m_mode == CK_Mode ? 6 : COLL_INT_POLY_DEGREE);
    vector<double> ca(degree + 1), cb(degree + 1), cc(degree + 1);

    for (size_t k = 0; k < degree + 1; k++) {
        ca[k] = astar_coeffs[k];
        cb[k] = bstar_coeffs[k];
        cc[k] = cstar_coeffs[k];
    }

    m_astar_poly.push_back(ca);
    m_bstar_poly.push_back(cb);
    m_cstar_poly.push_back(cc);
    m_poly[i][j] = static_cast<int>(m_astar_poly.size()) - 1;
    m_poly[j][i] = m_poly[i][j];
    if (actualT) {
        m_star_poly_uses_actualT[i][j] = 1;
        m_star_poly_uses_actualT[j][i] = m_star_poly_uses_actualT[i][j];
    }

    m_visc_ok = false;
    m_spvisc_ok = false;
    m_viscwt_ok = false;
    m_bindiff_ok = false;
    m_temp = -1;
}