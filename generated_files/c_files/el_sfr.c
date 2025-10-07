/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/el_sfr/el_sfr.c
 * \date        08/2025
 * \author      Ezequiel Lozano
 * \brief       Compute the star formation rate for a given gas cell.
 * \details     This file contains the routines to compute the star formation rate, according to our
 *              star formation model. The model evolves a set of ODEs which describe the mass
 *              exchange between the different phases of hydrogen (ionized, atomic, and molecular),
 *              and stars, metals, and dust.
 *              contains functions:
 *                double *normalize_vector(double *vec, int size)
 *                int has_negative(double *vec, int size)
 *                void *read_ftable(const char *file_path, const int n_rows, const int n_cols)
 *                static double interpolate1D(double x, const void *dtable)
 *                static double interpolate2D(double x, double y, const void *dtable)
 *                static double J21(double z)
 *                static int sf_ode(double t, const double y[], double f[], void *parameters)
 *                static int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *parameters)
 *                static void integrate_ode(double *ic, double *parameters, double it)
 *                static double compute_gas_sfr(const int index)
 *                double rate_of_star_formation(const int index)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <libgen.h>

#ifdef TESTING
#include "./el_sfr.h"
#else /* #ifdef TESTING */
#include "../allvars.h"
#include "../proto.h"
#endif /* #ifdef TESTING */

#ifdef EL_SFR

/***************************************************************************************************
 * Auxilary Functions
 ***************************************************************************************************/

/*! \brief Renormalize a vector of fractions.
 *
 *  Sets every element that is negative to 0, then renormalizes the vector such
 *  that the sum equals 1. If all elements are zero or negative, the resulting
 *  vector will contain NaN values due to division by zero.
 *
 *  This function mutates the vector in-place.
 *
 *  \param[in,out] vec Pointer to the array of doubles to be normalized.
 *  \param[in] size Size of the array. Must be > 0.
 *
 *  \return Pointer to the normalized vector (same as input vec).
 */
double *normalize_vector(double *vec, int size)
{
    double sum = 0.0;

    for (int i = 0; i < size; i++)
    {
        vec[i] = fmax(vec[i], 0.0);
        sum += vec[i];
    }

    double factor = 1.0 / sum;

    for (int i = 0; i < size; i++)
    {
        vec[i] *= factor;
    }

    return vec;
}

/*! \brief Check for negative elements.
 *
 *  Returns 1 if any element of the array is negative, 0 otherwise.
 *
 *  \param[in] vec Pointer to the array of doubles.
 *  \param[in] size Size of the array.
 *
 *  \return 1 if the array has any negative elements, 0 otherwise.
 */
int has_negative(double *vec, int size)
{
    for (int i = 0; i < size; i++)
    {
        if (vec[i] < 0.0)
        {
            return 1;
        }
    }

    return 0;
}

/*! \brief Load a text file into memory as a C array.
 *
 *  Read a text file with space separated values and store the numbers as doubles in a C array.
 *
 *  The value on row i and column j will be stored as element i * n_cols + j in the array.
 *
 *  \param[in] file_path Path to the file.
 *  \param[in] n_rows Number of rows in the table.
 *  \param[in] n_cols Number of columns in the table.
 *
 *  \return The values of the text file in a `data_table` struct.
 */
void *read_ftable(const char *file_path, const int n_rows, const int n_cols)
{
    FILE *file_ptr = fopen(file_path, "r");
    if (!file_ptr)
    {
        fprintf(stderr, "Error: could not open file %s\n", file_path);
        return NULL;
    }

    int count = n_rows * n_cols;
    double *data = (double *)malloc(count * sizeof(double));
    for (int i = 0; i < count; i++)
    {
        if (fscanf(file_ptr, "%lf", &data[i]) != 1)
        {
            #ifndef TESTING
            terminate("Error in read_ftable(): failed to read data at index %d\n", i);
            #else /* #ifndef TESTING */
            fprintf(stderr, "Error in read_ftable(): failed to read data at index %d\n", i);
            #endif /* #ifndef TESTING */

            free(data);
            fclose(file_ptr);

            return NULL;
        }
    }

    fclose(file_ptr);

    data_table *dtable = malloc(sizeof(data_table));
    dtable->data = data;
    dtable->n_rows = n_rows;
    dtable->n_cols = n_cols;

    return (void *)dtable;
}

/*! \brief Use linear interpolation to approximate f(x).
 *
 *  If the value is out of range, the boundaries of the function are used (flat extrapolation).
 *
 *  \param[in] x Value at which the interpolation function will be evaluated.
 *  \param[in] dtable Path to the table containing the values of f(xi), for several xi.
 *
 *  \return The interpolation function evaluated at `x`.
 */
#ifdef TESTING
double interpolate1D(double x, const char *table_path, const int NROWS, const int NCOLS)
{
    void *dtable = read_ftable(table_path, NROWS, NCOLS);
#else  /* #ifdef TESTING */
static double interpolate1D(double x, const void *dtable)
{
#endif /* #ifdef TESTING */

    data_table *interp_table = (data_table *)dtable;
    double *data = interp_table->data;
    int n_rows = interp_table->n_rows;
    int n_cols = interp_table->n_cols;

    double *xa = malloc(n_rows * sizeof(double));
    double *ya = malloc(n_rows * sizeof(double));

    for (size_t i = 0; i < n_rows; ++i)
    {
        xa[i] = data[i * n_cols];
    }

    for (size_t i = 0; i < n_rows; ++i)
    {
        ya[i] = data[i * n_cols + 1];
    }

    double x1, x2;
    int idx_x1, idx_x2;

    if (x < xa[0])
    {
        x = xa[0];
        x1 = xa[0];
        x2 = xa[1];
        idx_x1 = 0;
        idx_x2 = 1;
    }
    else if (x > xa[n_rows - 1])
    {
        x = xa[n_rows - 1];
        x1 = xa[n_rows - 2];
        x2 = xa[n_rows - 1];
        idx_x1 = n_rows - 2;
        idx_x2 = n_rows - 1;
    }
    else
    {
        for (size_t i = 1; i < n_rows; ++i)
        {
            if (xa[i - 1] <= x && x <= xa[i])
            {
                x1 = xa[i - 1];
                x2 = xa[i];
                idx_x1 = i - 1;
                idx_x2 = i;
                break;
            }
        }
    }

    double fQ1 = ya[idx_x1];
    double fQ2 = ya[idx_x2];

    double coeff = 1 / (x2 - x1);
    double term_01 = fQ1 * (x2 - x);
    double term_02 = fQ2 * (x - x1);

    free(xa);
    free(ya);

    return coeff * (term_01 + term_02);
}

/*! \brief Use bilinear interpolation to approximate f(x, y).
 *
 *  If the values are out of range, the boundaries of the function are used (flat extrapolation).
 *
 *  \param[in] x X coordinate at which the interpolation function will be evaluated.
 *  \param[in] y Y coordinate at which the interpolation function will be evaluated.
 *  \param[in] dtable Path to the table containing the values of f(xi, yi), for many xi and yi.
 *
 *  \return The interpolation function evaluated at `x` and `y`.
 */
#ifdef TESTING
double interpolate2D(double x, double y, const char *table_path)
{
    void *dtable = read_ftable(table_path, ETA_NROWS, ETA_NCOLS);
#else  /* #ifdef TESTING */
static double interpolate2D(double x, double y, const void *dtable)
{
#endif /* #ifdef TESTING */

    data_table *interp_table = (data_table *)dtable;
    double *data = interp_table->data;
    int n_rows = interp_table->n_rows;
    int n_cols = interp_table->n_cols;

    int nx = n_rows - 1;
    int ny = n_cols - 1;

    double *xa = malloc(nx * sizeof(double));
    double *ya = malloc(ny * sizeof(double));
    double *za = malloc(nx * ny * sizeof(double));

    for (size_t i = 0; i < nx; ++i)
    {
        xa[i] = data[(i + 1) * n_cols];
    }

    for (size_t i = 0; i < ny; ++i)
    {
        ya[i] = data[i + 1];
    }

    for (size_t i = 0; i < nx; ++i)
    {
        for (size_t j = 0; j < ny; ++j)
        {
            za[i * ny + j] = data[(i + 1) * n_cols + (j + 1)];
        }
    }

    double x1, x2, y1, y2;
    int idx_x1, idx_x2, idx_y1, idx_y2;

    if (x < xa[0])
    {
        x = xa[0];
        x1 = xa[0];
        x2 = xa[1];
        idx_x1 = 0;
        idx_x2 = 1;
    }
    else if (x > xa[nx - 1])
    {
        x = xa[nx - 1];
        x1 = xa[nx - 2];
        x2 = xa[nx - 1];
        idx_x1 = nx - 2;
        idx_x2 = nx - 1;
    }
    else
    {
        for (size_t i = 1; i < nx; ++i)
        {
            if (xa[i - 1] <= x && x <= xa[i])
            {
                x1 = xa[i - 1];
                x2 = xa[i];
                idx_x1 = i - 1;
                idx_x2 = i;
                break;
            }
        }
    }

    if (y < ya[0])
    {
        y = ya[0];
        y1 = ya[0];
        y2 = ya[1];
        idx_y1 = 0;
        idx_y2 = 1;
    }
    else if (y > ya[ny - 1])
    {
        y = ya[ny - 1];
        y1 = ya[ny - 2];
        y2 = ya[ny - 1];
        idx_y1 = ny - 2;
        idx_y2 = ny - 1;
    }
    else
    {
        for (size_t i = 1; i < ny; ++i)
        {
            if (ya[i - 1] <= y && y <= ya[i])
            {
                y1 = ya[i - 1];
                y2 = ya[i];
                idx_y1 = i - 1;
                idx_y2 = i;
                break;
            }
        }
    }

    double fQ11 = za[idx_x1 * ny + idx_y1];
    double fQ21 = za[idx_x2 * ny + idx_y1];
    double fQ12 = za[idx_x1 * ny + idx_y2];
    double fQ22 = za[idx_x2 * ny + idx_y2];

    double coeff = 1 / ((x2 - x1) * (y2 - y1));
    double term_01 = fQ11 * (x2 - x) * (y2 - y);
    double term_02 = fQ21 * (x - x1) * (y2 - y);
    double term_03 = fQ12 * (x2 - x) * (y - y1);
    double term_04 = fQ22 * (x - x1) * (y - y1);

    free(xa);
    free(ya);
    free(za);

    return coeff * (term_01 + term_02 + term_03 + term_04);
}

/*! \brief Compute the LW background radiation intensity.
 *
 *  \param[in] z Redshift.
 *
 *  \return The LW background radiation intensity in [10^(-21) erg s^(-1) cm^(-2) Hz^(-1) sr^(-1)].
 */
#ifdef TESTING
double J21(double z)
#else  /* #ifdef TESTING */
static double J21(double z)
#endif /* #ifdef TESTING */
{
    double LWB_A = 2.119;
    double LWB_B = -0.1117;
    double LWB_C = -0.002782;

    double z1 = 1.0 + z;
    double logJ21 = LWB_A + LWB_B * z1 + LWB_C * z1 * z1;

    return pow(10, logJ21);
}

/***************************************************************************************************
 * ODE functions
 ***************************************************************************************************/

/*! \brief Evaluate the systems of equations.
 *
 *  Evaluate the six ODEs of the model, using the following variables:
 *
 *  Ionized gas fraction:    fi(t) = Mi(t) / MC --> y[0]
 *  Atomic gas fraction:     fa(t) = Ma(t) / MC --> y[1]
 *  Molecular gas fraction:  fm(t) = Mm(t) / MC --> y[2]
 *  Stellar fraction:        fs(t) = Ms(t) / MC --> y[3]
 *  Metal fraction:          fZ(t) = MZ(t) / MC --> y[4]
 *  Dust fraction:           fd(t) = Md(t) / MC --> y[5]
 *
 *  where MC = Mi(t) + Ma(t) + Mm(t) + Ms(t) + MZ(t) + Md(t) is the total density of the gas cell,
 *  and each equation has units of Myr^(-1).
 *
 *  \param[in] t Unused variable to comply with the `gsl_odeiv2_driver_alloc_y_new()` API.
 *  \param[in] y Values of the variables at which the ODEs will be evaluated.
 *  \param[out] f Where the results of evaluating the ODEs will be stored.
 *  \param[in] parameters Parameters for the ODEs.
 *
 *  \return Constant `GSL_SUCCESS`, to confirm that the computation was successful.
 */
static int sf_ode(double t, const double y[], double f[], void *parameters)
{
    (void)(t);

    /**********************************************************************************************
     * Initial conditions
     **********************************************************************************************/

    double fi = y[0];
    double fa = y[1];
    double fm = y[2];
    double fs = y[3];
    double fZ = y[4];
    double fd = y[5];

    /**********************************************************************************************
     * Parameters
     **********************************************************************************************
     *
     * rho_c: Total cell density           [mp * cm^(-3)]
     * UVB:   UVB photoionization rate     [Myr^(-1)]
     * LWB:   LWB photodissociation rate   [Myr^(-1)]
     * eta_d: Photodissociation efficiency [dimensionless]
     * eta_i: Photoionization efficiency   [dimensionless]
     * R:     Mass recycling fraction      [dimensionless]
     * Zsn:   Metals recycling fraction    [dimensionless]
     * h:     Column height                [cm]
     **********************************************************************************************/

    double *p = (double *)parameters;
    double rho_c = p[0];
    double UVB = p[1];
    double LWB = p[2];
    double eta_d = p[3];
    double eta_i = p[4];
    double R = p[5];
    double Zsn = p[6];
    double h = p[7];

    /**********************************************************************************************
     * Auxiliary equations
     **********************************************************************************************/

    /* Star formation rate [Myr^(-1)] */
    double psi = ODE_CS * fm * sqrt(rho_c);

    /*******************
     * Partial fraction
     *******************/

    /* Neutral fraction */
    double fn = fa + fm;

    /* Gas fraction */
    double fg = fn + fi;

    /* Metal fraction */
    double Zt = fZ + fd;

    /************
	 * Shielding
     ************/

    /* Dust optical depth */
    double tau_d = ODE_CSD * h * rho_c * fn * Zt;

	/* Dust shielding */
    double sd = exp(-tau_d);

    /* Molecular self-shielding optical depth */
    double xp1 = ODE_CSH2 * h * fm * rho_c + 1.0;
    double xp1_2 = xp1 * xp1;
    double sq_xp1 = sqrt(xp1);
    double tau_sh = 8.5e-4 * sq_xp1;

    /* Molecular self-shielding */
    double exp_xp1 = exp(-tau_sh);
    double sh2 = ((1 - WH2) / xp1_2) + (WH2 * exp_xp1 / sq_xp1);

    /* Combined shielding  */
	double e_diss = sd * sh2;

    /*****************
     * Net ionization
     *****************/

    /* Recombination [Myr^(-1)] */
    double recomb = ODE_CREC * fi * fi * rho_c;

    /* Ionization optical depth */
	double tau_ion = ODE_CTION * fa * rho_c * h;

    /* Stellar ionization [Myr^(-1)] */
	double ion_prob = -expm1(-tau_ion);
    double s_ion = sd * eta_i * ion_prob * psi;

    /* UVB photoionization [Myr^(-1)] */
    double uvb = sd * UVB * fa;

    double net_ionization = uvb + s_ion - recomb;

    /*************************
     * Stellar gas production
     *************************/

    /* Gas and metals produced at stellar death [Myr^(-1)] */
    double R_psi = R * psi;

    double s_gas_production = fma(R_psi, -Zsn, R_psi);

    /*******************
     * Net dissociation
     *******************/

    /* Condensation [Myr^(-1)] */
    double Zt_eff = Zt + ZEFF;
    double cond = ODE_CCOND * fa * rho_c * fg * Zt_eff;

    /* Dissociation optical depth */
    double tau_diss = ODE_CTDISS * fm * rho_c * h;

    /* Stellar dissociation [Myr^(-1)] */
	double diss_prob = -expm1(-tau_diss);
    double s_diss = e_diss * eta_d * diss_prob * psi;

    /* LW background dissociation [Myr^(-1)] */
	double lwb = e_diss * LWB * fm;

    double net_dissociation = lwb + s_diss - cond;

    /******************
     * Net dust growth
     ******************/

    /* Dust growth  [Myr^(-1)] */
    double dg = ODE_CDG * fZ * fd * fn * fn * rho_c;

    /* Dust destruction [Myr^(-1)] */
    double dd = fd * INV_T_DD;

    double net_dust_growth = dg - dd;

    /**********************************************************************************************
     * Evaluate the ODE system
     **********************************************************************************************/

    f[0] = net_ionization + s_gas_production;
    f[1] = net_dissociation - net_ionization;
    f[2] = -net_dissociation - psi;
    f[3] = fma(R, -psi, psi);
    f[4] = fma(Zsn, R_psi, -net_dust_growth);
    f[5] = net_dust_growth;

    return GSL_SUCCESS;
}

// JACOBIAN_START
// JACOBIAN_END

/*! \brief Solve the system of ODEs, using numerical integration.
 *
 *  ICs:
 *
 *  Ionized gas fraction:    fi = Mi / MC [dimensionless]
 *  Atomic gas fraction:     fa = Ma / MC [dimensionless]
 *  Molecular gas fraction:  fm = Mm / MC [dimensionless]
 *  Stellar fraction:        fs = Ms / MC [dimensionless]
 *  Metal fraction:          fZ = MZ / MC [dimensionless]
 *  Dust fraction:           fd = Md / MC [dimensionless]
 *
 *  Parameters
 *
 *  rho_C: Total cell density           [mp * cm^(-3)]
 *  UVB:   UVB photoionization rate     [Myr^(-1)]
 *  LWB:   LWB photodissociation rate   [Myr^(-1)]
 *  eta_d: Photodissociation efficiency [dimensionless]
 *  eta_i: Photoionization efficiency   [dimensionless]
 *  R:     Mass recycling fraction      [dimensionless]
 *  Zsn:   Metals recycling fraction    [dimensionless]
 *  h:     Column height                [cm]
 *
 *  \param[in] ic Initial conditions.
 *  \param[in] parameters Parameters for the ODEs.
 *  \param[in] it Integration time in Myr.
 *
 *  \return void.
 */
#ifdef TESTING
void integrate_ode(double *ic, double *parameters, double it)
#else  /* #ifdef TESTING */
static void integrate_ode(double *ic, double *parameters, double it)
#endif /* #ifdef TESTING */
{
    /* Setup and run ODE integration */
    double t0 = 0.0;
    gsl_odeiv2_system sys = {sf_ode, NULL, N_EQU, parameters};
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msadams, it * 1e-4, 1e-10, 0.0);

    int status = gsl_odeiv2_driver_apply(driver, &t0, it, ic);

#ifndef TESTING
    if (status != GSL_SUCCESS)
    {
        terminate("GSL ERROR: Could not integrate, status = %d\n", status);
    }
#endif /* #ifndef TESTING */

    gsl_odeiv2_driver_free(driver);

#ifndef TESTING
    /* Check that none of the outputs are negative */
    if (has_negative(ic, N_EQU))
    {
        warn("WARNING: Negative equation outputs: \nfi = %.9lf, fa = %.9lf, fm = %.9lf, fs = %.9lf, fZ = %.9lf, fd = %.9lf ", ic[0], ic[1], ic[2], ic[3], ic[4], ic[5]);

        /* Set to 0 the fractions that are negative and renormalize */
        normalize_vector(ic, N_EQU);
    }
#endif /* #ifndef TESTING */
}

#ifndef TESTING

/*! \brief Compute the star formation rate.
 *
 *  Compute the star formation rate using the stellar mass fraction of the cell.
 *
 *  \param[in] index Index of the gas cell in question.
 *
 *  \return The star formation rate in [Mₒ yr^(-1)].
 */
static double compute_gas_sfr(const int index)
{
    double integration_time = SphP[index].integration_time * 1000000; // [yr]
    double cell_mass = P[index].Mass * M_COSMO;                       // [Mₒ]

    return SphP[index].ODE_fractions[3] * cell_mass / integration_time;
}

/*! \brief Compute the star formation rate.
 *
 *  Compute the star formation rate solving a set of ODEs.
 *
 *  \param[in] index Index of the gas cell in question.
 *
 *  \return The star formation rate in [Mₒ yr^(-1)].
 */
double rate_of_star_formation(const int index)
{
    /* Check for cases when the SFR should be 0, to shortcut the computation */
    if (!P[index].TimeBinHydro)
    {
        return 0.0;
    }
    if (All.ComovingIntegrationOn && SphP[index].Density < All.OverDensThresh)
    {
        return 0.0;
    }
    if (!All.ComovingIntegrationOn && SphP[index].Density * All.cf_a3inv < All.PhysDensThresh)
    {
        return 0.0;
    }

    /* Check if we are still in the same global timestep, to shortcut the computation */
    if (!isnan(SphP[index].parameter_a))
    {
        if (All.Time == SphP[index].parameter_a)
        {
            return compute_gas_sfr(index);
        }
    }

    /**********************************************************************************************
     * Compute the Integration time
     **********************************************************************************************/

    /* Integration time [Myr] */
    double integration_time = (((integertime)1) << P[index].TimeBinHydro) * All.Timebase_interval;
    if (integration_time <= 0.0)
    {
        return 0.0;
    }
    integration_time *= T_MYR * All.cf_atime / All.cf_time_hubble_a;

    /* Integration time as log10(it [yr]) for the interpolation tables */
    double log_it = log10(integration_time * 1.0e6);

    /* Store the integration time */
    SphP[index].integration_time = integration_time;

    /**********************************************************************************************
     * Redshift
     **********************************************************************************************/

    double redshift = All.cf_redshift;

    if (!All.ComovingIntegrationOn)
    {
        double t = All.Time; // Gyr
        double t2 = t * t;
        double t3 = t2 * t;
        double at_01 = 0.0185084;
        double at_02 = 0.00151085;
        double at_03 = -0.0207921;
        double at_04 = 0.149462;

        double a = at_01 + at_02 * t3 + at_03 * t2 + at_04 * t;

        /* Scale factor to redshift */
        redshift = (1.0 / a) - 1.0;
    }
    
    /**********************************************************************************************
     * Compute the ODE parameters
     *
     * rho_C: Total cell density           [mp * cm^(-3)]
     * UVB:   UVB photoionization rate     [Myr^(-1)]
     * LWB:   LWB photodissociation rate   [Myr^(-1)]
     * eta_d: Photodissociation efficiency [dimensionless]
     * eta_i: Photoionization efficiency   [dimensionless]
     * R:     Mass recycling fraction      [dimensionless]
     * Zsn:   Metals recycling fraction    [dimensionless]
     * h:     Column height                [cm]
     **********************************************************************************************/

    /* Cell density [mp * cm^(-3)] */
    double rhoC = SphP[index].Density * RHO_COSMO;

    /* UVB photoionization rate [Myr^(-1)] */
    double UVB = interpolate1D(redshift, All.UVB_TABLE_DATA);

    /* LWB photodissociation rate [Myr^(-1)] */
    double LWB = ABEL97 * J21(redshift);

    /* Metallicity [dimensionless] */
    double Z = fmax(0.0, SphP[index].Metallicity);

    /* Photodissociation efficiency [dimensionless] */
    double eta_d = interpolate2D(log_it, Z, All.ETA_D_TABLE_DATA);

    /* Photoionization efficiency [dimensionless] */
    double eta_i = interpolate2D(log_it, Z, All.ETA_I_TABLE_DATA);

    /* Mass recycling fraction [dimensionless] */
    double R = interpolate1D(Z, All.R_TABLE_DATA);

    /* Metals recycling fraction [dimensionless] */
    double Zsn = interpolate1D(Z, All.ZSN_TABLE_DATA);

    /* Column height [cm] */
    double M = P[index].Mass * M_CGS;
    
    /* Cell volume [cm^3] */
    double V = M / rhoC;

    /* Column height [cm] */
    double h = cbrt(V / (4 * M_PI / 3));

    double parameters[] = {rhoC, UVB, LWB, eta_d, eta_i, R, Zsn, h};

    /* Store the ODE parameters */
    SphP[index].tau_S = ODE_CS / sqrt(rhoC);
    SphP[index].parameter_a = All.Time;
    SphP[index].parameter_z = redshift;
    SphP[index].parameter_rhoC = rhoC;
    SphP[index].parameter_UVB = UVB;
    SphP[index].parameter_LWB = LWB;
    SphP[index].parameter_Z = Z;
    SphP[index].parameter_eta_d = eta_d;
    SphP[index].parameter_eta_i = eta_i;
    SphP[index].parameter_R = R;
    SphP[index].parameter_Zsn = Zsn;
    SphP[index].parameter_h = h;

    /**********************************************************************************************
     * Compute the initial conditions
     **********************************************************************************************/

    double fi, fa, fm, fs, fZ, fd;

    double nhp, nh, df;

    get_arepo_fraction(index, &nhp, &nh);

    /* Ionized gas mass fraction [dimensionless] */
    fi = (1 - Z) * nhp / (nhp + nh);

    /* Atomic gas mass fraction [dimensionless] */
    fa = (1 - Z) * nh / (nhp + nh);

    /* Molecular gas mass fraction [dimensionless] */
    fm = 0.0;

    /* Stellar mass fraction [dimensionless] */
    fs = 0.0;

    /* Dust fraction in the metals [dimensionless] */
    df = ODE_CXD * fa;

    /* Metal mass fraction [dimensionless] */
    fZ = Z * (1 - df);

    /* Dust mass fraction [dimensionless] */
    fd = Z * df;

    double ic[] = {fi, fa, fm, fs, fZ, fd};

    /**********************************************************************************************
     * Integrate the ODEs
     **********************************************************************************************/

    integrate_ode(ic, parameters, integration_time);

    /* Store the results after solving the ODEs */
    for (size_t i = 0; i < N_EQU; ++i)
    {
        SphP[index].ODE_fractions[i] = ic[i];
    }

    /* Return the star formation rate in [Mₒ yr^(-1)] */
    return compute_gas_sfr(index);
}

#endif /* #ifndef TESTING */

#endif /* #ifdef EL_SFR */
