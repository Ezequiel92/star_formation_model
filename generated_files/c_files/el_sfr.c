/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/el_sfr/el_sfr.c
 * \date        11/2026
 * \author      Ezequiel Lozano
 * \brief       Compute the star formation rate for a given gas cell.
 * \details     This file contains the routines to compute the star formation rate, according to our
 *              star formation model. The model evolves a set of ODEs which describe the mass
 *              exchange between the different phases of hydrogen (ionized, atomic, and molecular),
 *              and stars, metals, and dust.
 *              contains functions:
 *                double *normalize_array(double *arr, int size)
 *                int has_negative(double *arr, int size)
 *                void *read_ftable(const char *file_path, const int n_rows, const int n_cols)
 *                static int find_interval(double *x, const double *axis, const int n, const int stride)
 *                static double linear_interpolation(double x, const data_table *table)
 *                static double bilinear_interpolation(double x, double y, const data_table *table)
 *                static const data_table *cached_table(const char *table_path, const int n_rows, const int n_cols)
 *                static double interpolate1D(double x, const void *dtable)
 *                static double interpolate2D(double x, double y, const void *dtable)
 *                static double J21(double z)
 *                static int sf_ode(double t, const double y[], double f[], void *parameters)
 *                static int integrate_ode(double *ic, double *parameters, double it)
 *                static void two_phase_fractions(const int index, double *fi_gas, double *fn_gas)
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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#ifdef TESTING

#include <string.h>
#include "./el_sfr.h"

/* Maximum number of tables kept in memory by `cached_table()` */
#define N_CACHED_TABLES 8
/* Maximum length of the paths of the tables kept in memory by `cached_table()` */
#define CACHED_PATH_LENGTH 4096

#else /* #ifdef TESTING */

#include "../allvars.h"
#include "../proto.h"

#endif /* #ifdef TESTING */

#ifdef EL_SFR

/***************************************************************************************************
 * Auxiliary Functions
 ***************************************************************************************************/

/*! \brief Renormalize an array of fractions.
 *
 *  Sets every element that is negative to 0, then renormalizes the array such
 *  that the sum of its elements equals 1. If all elements are zero or negative,
 *  the resulting array will contain NaN values due to division by zero.
 *
 *  This function mutates the array in-place.
 *
 *  \param[in,out] arr Pointer to the array of doubles.
 *  \param[in] size Size of the array. Must be > 0.
 *
 *  \return Pointer to the normalized array (same as input arr).
 */
double *normalize_array(double *arr, int size)
{
    double sum = 0.0;

    for (int i = 0; i < size; i++)
    {
        arr[i] = fmax(arr[i], 0.0);
        sum += arr[i];
    }

    double factor = 1.0 / sum;

    for (int i = 0; i < size; i++)
    {
        arr[i] *= factor;
    }

    return arr;
}

/*! \brief Check for negative elements.
 *
 *  Returns 1 if any element of the array is negative, 0 otherwise.
 *
 *  \param[in] arr Pointer to the array of doubles.
 *  \param[in] size Size of the array. Must be > 0.
 *
 *  \return 1 if the array has any negative elements, 0 otherwise.
 */
int has_negative(double *arr, int size)
{
    for (int i = 0; i < size; i++)
    {
        if (arr[i] < 0.0)
        {
            return 1;
        }
    }

    return 0;
}

/*! \brief Load a text file into memory as a C array.
 *
 *  Read a text file with space-separated values and store the numbers as doubles in a C array.
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

    int count    = n_rows * n_cols;
    double *data = (double *)malloc(count * sizeof(double));
    if (!data)
    {
        fprintf(stderr, "Error in read_ftable(): could not allocate memory for %s\n", file_path);
        fclose(file_ptr);
        return NULL;
    }

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
    if (!dtable)
    {
        fprintf(stderr, "Error in read_ftable(): could not allocate memory for %s\n", file_path);
        free(data);
        return NULL;
    }

    dtable->data   = data;
    dtable->n_rows = n_rows;
    dtable->n_cols = n_cols;

    return (void *)dtable;
}

/*! \brief Find the interval of a sorted axis that contains a value.
 *
 *  The axis values are `axis[i * stride]`, for i = 0, ..., n - 1, sorted in increasing order. The
 *  stride allows the use of a row or a column of a table as the axis, without copying it. Values
 *  out of range are moved to the closest boundary (flat extrapolation).
 *
 *  \param[in,out] x Value to locate. If it is out of range, it is replaced by the closest boundary.
 *  \param[in] axis Pointer to the first value of the axis.
 *  \param[in] n Number of values in the axis. Must be >= 2.
 *  \param[in] stride Distance between consecutive values of the axis in memory.
 *
 *  \return Index i of the interval, such that axis[i * stride] <= x <= axis[(i + 1) * stride],
 *          with 0 <= i <= n - 2.
 */
static int find_interval(double *x, const double *axis, const int n, const int stride)
{
    if (*x <= axis[0])
    {
        *x = axis[0];
        return 0;
    }

    if (*x >= axis[(n - 1) * stride])
    {
        *x = axis[(n - 1) * stride];
        return n - 2;
    }

    int low  = 0;
    int high = n - 1;

    while (high - low > 1)
    {
        int mid = (low + high) / 2;

        if (axis[mid * stride] <= *x)
        {
            low = mid;
        }
        else
        {
            high = mid;
        }
    }

    return low;
}

/*! \brief Use linear interpolation to approximate f(x), from a table in memory.
 *
 *  The table has two columns, xi in the first one (sorted in increasing order) and f(xi) in the
 *  second one. If `x` is out of range, the boundaries of the function are used (flat extrapolation).
 *
 *  \param[in] x Value at which the interpolation function will be evaluated.
 *  \param[in] table Table with the values of f(xi), for several xi.
 *
 *  \return The interpolation function evaluated at `x`.
 */
static double linear_interpolation(double x, const data_table *table)
{
    const double *data = table->data;
    const int n_rows   = table->n_rows;
    const int n_cols   = table->n_cols;

    /* xi are the values of the first column */
    int i = find_interval(&x, data, n_rows, n_cols);

    double x1 = data[i * n_cols];
    double x2 = data[(i + 1) * n_cols];

    double fQ1 = data[i * n_cols + 1];
    double fQ2 = data[(i + 1) * n_cols + 1];

    double coeff   = 1.0 / (x2 - x1);
    double term_01 = fQ1 * (x2 - x);
    double term_02 = fQ2 * (x - x1);

    return coeff * (term_01 + term_02);
}

/*! \brief Use bilinear interpolation to approximate f(x, y), from a table in memory.
 *
 *  The xi values are in the first column (from the second row onwards), the yi values are in the
 *  first row (from the second column onwards), and f(xi, yi) fill the rest of the table. Both xi
 *  and yi are sorted in increasing order. If the values are out of range, the boundaries of the
 *  function are used (flat extrapolation).
 *
 *  \param[in] x X coordinate at which the interpolation function will be evaluated.
 *  \param[in] y Y coordinate at which the interpolation function will be evaluated.
 *  \param[in] table Table with the values of f(xi, yi), for many xi and yi.
 *
 *  \return The interpolation function evaluated at `x` and `y`.
 */
static double bilinear_interpolation(double x, double y, const data_table *table)
{
    const double *data = table->data;
    const int n_cols   = table->n_cols;

    /* Number of values of xi and yi */
    const int nx = table->n_rows - 1;
    const int ny = n_cols - 1;

    /* xi are the values of the first column, from the second row onwards */
    int i = find_interval(&x, data + n_cols, nx, n_cols);

    /* yi are the values of the first row, from the second column onwards */
    int j = find_interval(&y, data + 1, ny, 1);

    double x1 = data[(i + 1) * n_cols];
    double x2 = data[(i + 2) * n_cols];
    double y1 = data[j + 1];
    double y2 = data[j + 2];

    /* f(xi, yj) is the element (i + 1, j + 1) of the table */
    double fQ11 = data[(i + 1) * n_cols + (j + 1)];
    double fQ21 = data[(i + 2) * n_cols + (j + 1)];
    double fQ12 = data[(i + 1) * n_cols + (j + 2)];
    double fQ22 = data[(i + 2) * n_cols + (j + 2)];

    double coeff   = 1.0 / ((x2 - x1) * (y2 - y1));
    double term_01 = fQ11 * (x2 - x) * (y2 - y);
    double term_02 = fQ21 * (x - x1) * (y2 - y);
    double term_03 = fQ12 * (x2 - x) * (y - y1);
    double term_04 = fQ22 * (x - x1) * (y - y1);

    return coeff * (term_01 + term_02 + term_03 + term_04);
}

#ifdef TESTING

/*! \brief Get a table from a text file, reading the file only the first time.
 *
 *  The tables are kept in memory until the program ends, so each file is read only once, no matter
 *  how many times the interpolation functions are called. The cached tables are not updated if
 *  the file changes, so the library has to be reloaded if the tables change.
 *
 *  \param[in] table_path Path to the file.
 *  \param[in] n_rows Number of rows in the table.
 *  \param[in] n_cols Number of columns in the table.
 *
 *  \return The table, or NULL if the file could not be read or the cache is full.
 */
static const data_table *cached_table(const char *table_path, const int n_rows, const int n_cols)
{
    /* Tables already read, and their paths */
    static const data_table *tables[N_CACHED_TABLES];
    static char paths[N_CACHED_TABLES][CACHED_PATH_LENGTH];
    static int n_tables = 0;

    for (int i = 0; i < n_tables; i++)
    {
        if (strcmp(paths[i], table_path) == 0 && tables[i]->n_rows == n_rows && tables[i]->n_cols == n_cols)
        {
            return tables[i];
        }
    }

    if (n_tables == N_CACHED_TABLES || strlen(table_path) >= CACHED_PATH_LENGTH)
    {
        fprintf(stderr, "Error in cached_table(): could not cache the table %s\n", table_path);
        return NULL;
    }

    const data_table *table = (const data_table *)read_ftable(table_path, n_rows, n_cols);
    if (!table)
    {
        return NULL;
    }

    strcpy(paths[n_tables], table_path);
    tables[n_tables] = table;
    n_tables++;

    return table;
}

#endif /* #ifdef TESTING */

/*! \brief Use linear interpolation to approximate f(x).
 *
 *  If the value is out of range, the boundaries of the function are used (flat extrapolation).
 *
 *  \param[in] x Value at which the interpolation function will be evaluated.
 *  \param[in] dtable Table containing the values of f(xi), for several xi.
 *
 *  \return The interpolation function evaluated at `x`.
 */
#ifdef TESTING
double interpolate1D(double x, const char *table_path, const int NROWS, const int NCOLS)
{
    const data_table *table = cached_table(table_path, NROWS, NCOLS);
    if (!table)
    {
        return NAN;
    }

    return linear_interpolation(x, table);
}
#else  /* #ifdef TESTING */
static double interpolate1D(double x, const void *dtable)
{
    return linear_interpolation(x, (const data_table *)dtable);
}
#endif /* #ifdef TESTING */

/*! \brief Use bilinear interpolation to approximate f(x, y).
 *
 *  If the values are out of range, the boundaries of the function are used (flat extrapolation).
 *
 *  \param[in] x X coordinate at which the interpolation function will be evaluated.
 *  \param[in] y Y coordinate at which the interpolation function will be evaluated.
 *  \param[in] dtable Table containing the values of f(xi, yi), for many xi and yi.
 *
 *  \return The interpolation function evaluated at `x` and `y`.
 */
#ifdef TESTING
double interpolate2D(double x, double y, const char *table_path)
{
    const data_table *table = cached_table(table_path, ETA_NROWS, ETA_NCOLS);
    if (!table)
    {
        return NAN;
    }

    return bilinear_interpolation(x, y, table);
}
#else  /* #ifdef TESTING */
static double interpolate2D(double x, double y, const void *dtable)
{
    return bilinear_interpolation(x, y, (const data_table *)dtable);
}
#endif /* #ifdef TESTING */

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
    /* Redshift restricted to the validity range of the fit */
    double zeff = fmax(z, LWB_Z_MIN);

    double z1     = 1.0 + zeff;
    double logJ21 = LWB_A + LWB_B * z1 + LWB_C * z1 * z1;

    return pow(10, logJ21);
}

/***************************************************************************************************
 * ODE functions
 ***************************************************************************************************/

/*! \brief Evaluate the system of equations of our model.
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
 *  where MC = Mi(t) + Ma(t) + Mm(t) + Ms(t) + MZ(t) + Md(t) is the total mass of the gas cell,
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
    double fd = y[5];

    /**********************************************************************************************
     * Parameters
     **********************************************************************************************
     *
     * rho_C:  Total cell density           [mp * cm^(-3)]
     * UVB:    UVB photoionization rate     [Myr^(-1)]
     * LWB:    LWB photodissociation rate   [Myr^(-1)]
     * eta_d:  Photodissociation efficiency [dimensionless]
     * eta_i:  Photoionization efficiency   [dimensionless]
     * R:      Mass recycling fraction      [dimensionless]
     * Zsn:    Metallicity of SN ejecta     [dimensionless]
     * h:      Cell size                    [cm]
     * fd_max: Maximum dust fraction        [dimensionless]
     * c_dg:   Dust growth constant         [Myr^(-1) * mp^(-1) * cm^3]
     **********************************************************************************************/

    double *p     = (double *)parameters;
    double rho_c  = p[0];
    double UVB    = p[1];
    double LWB    = p[2];
    double eta_d  = p[3];
    double eta_i  = p[4];
    double R      = p[5];
    double Zsn    = p[6];
    double h      = p[7];
    double fd_max = p[8];
    double c_dg   = p[9];

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

    /************
	 * Shielding
     ************/

    /* Shielding length [cm] */
    double Lsh = ODE_CLSH / sqrt(rho_c);

    /* Dust optical depth */
    double tau_d = ODE_CSD * fd * rho_c * Lsh;

	/* Dust shielding */
    double sd = exp(-tau_d);

    /* Molecular self-shielding optical depth */
    double xp1    = ODE_CSH2 * fmax(fm, 0.0) * rho_c * Lsh + 1.0;
    double xp1_2  = xp1 * xp1;
    double sq_xp1 = sqrt(xp1);
    double tau_sh = ODE_EXPFAC * sq_xp1;

    /* Molecular self-shielding */
    double exp_xp1 = exp(-tau_sh);
    double sh2     = ((1 - WH2) / xp1_2) + (WH2 * exp_xp1 / sq_xp1);

    /* Combined shielding  */
	double e_diss = sd * sh2;

    /*****************
     * Net ionization
     *****************/

    /* Recombination [Myr^(-1)] */
    double recomb = ODE_CREC * fi * fmax(fi, 0.0) * rho_c;

    /* Ionization optical depth */
	double tau_ion = ODE_CTION * fa * rho_c * h;

    /* Stellar ionization [Myr^(-1)] */
	double ion_prob = -expm1(-tau_ion);
    double s_ion    = sd * eta_i * ion_prob * psi;

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
    double cond = ODE_CCOND * fa * rho_c * (fd + FD_EFF);

    /* Dissociation optical depth */
    double tau_diss = ODE_CTDISS * fm * rho_c * h;

    /* Stellar dissociation [Myr^(-1)] */
	double diss_prob = -expm1(-tau_diss);
    double s_diss    = e_diss * eta_d * diss_prob * psi;

    /* LW background dissociation [Myr^(-1)] */
	double lwb = e_diss * LWB * fm;

    double net_dissociation = lwb + s_diss - cond;

    /******************
     * Net dust growth
     ******************/

    /* Dust yield by SNe [Myr^(-1)] */
    double dg_SNe = DELTA_DUST * Zsn * R_psi;

    /* Gas-phase fraction of the key elements [dimensionless] */
    double xi = fmax(0.0, 1.0 - fd / fd_max);

    /* Dust accretion [Myr^(-1)] */
    double dg_acc = c_dg * xi * fd * fn * rho_c;

    /* Dust destruction [Myr^(-1)] */
    double dd = ODE_CDD * fd * psi;

    double net_dust_growth = dg_SNe + dg_acc - dd;

    /**********************************************************************************************
     * ODE system
     **********************************************************************************************/

    f[0] = net_ionization + s_gas_production;
    f[1] = net_dissociation - net_ionization;
    f[2] = -net_dissociation - psi;
    f[3] = fma(R, -psi, psi);
    f[4] = fma(Zsn, R_psi, -net_dust_growth);
    f[5] = net_dust_growth;

    return GSL_SUCCESS;
}

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
 *  rho_C:  Total cell density           [mp * cm^(-3)]
 *  UVB:    UVB photoionization rate     [Myr^(-1)]
 *  LWB:    LWB photodissociation rate   [Myr^(-1)]
 *  eta_d:  Photodissociation efficiency [dimensionless]
 *  eta_i:  Photoionization efficiency   [dimensionless]
 *  R:      Mass recycling fraction      [dimensionless]
 *  Zsn:    Metallicity of SN ejecta     [dimensionless]
 *  h:      Cell size                    [cm]
 *  fd_max: Maximum dust fraction        [dimensionless]
 *  c_dg:   Dust growth constant         [Myr^(-1) * mp^(-1) * cm^3]
 *
 *  \param[in] ic Initial conditions.
 *  \param[in] parameters Parameters for the ODEs.
 *  \param[in] it Integration time in Myr.
 *
 *  \return The status of the GSL integration (GSL_SUCCESS if it was successful). Outside of
 *          TESTING, a failed integration terminates the run, printing its inputs.
 */
#ifdef TESTING
int integrate_ode(double *ic, double *parameters, double it)
#else  /* #ifdef TESTING */
static int integrate_ode(double *ic, double *parameters, double it)
#endif /* #ifdef TESTING */
{
#ifndef TESTING
    /* Copy of the initial conditions, to report them if the integration fails */
    double ic0[N_EQU];
    for (size_t i = 0; i < N_EQU; ++i)
    {
        ic0[i] = ic[i];
    }
#endif /* #ifndef TESTING */

    /* Set up and run the ODE integration */
    double t0 = 0.0;
    gsl_odeiv2_system sys = {sf_ode, NULL, N_EQU, parameters};
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
        &sys,
        gsl_odeiv2_step_msadams,
        it * 1e-4,
        1e-10,
        0.0
    );

    int status = gsl_odeiv2_driver_apply(driver, &t0, it, ic);

#ifndef TESTING
    if (status != GSL_SUCCESS)
    {
        double *p = (double *)parameters;
        terminate(
            "GSL ERROR: Could not integrate, status = %d, reached t = %.17g of it = %.17g [Myr]\n"
            "ic = {%.17g, %.17g, %.17g, %.17g, %.17g, %.17g}\n"
            "parameters = {%.17g, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g, %.17g}\n"
            "state at failure = {%.17g, %.17g, %.17g, %.17g, %.17g, %.17g}\n",
            status, t0, it,
            ic0[0], ic0[1], ic0[2], ic0[3], ic0[4], ic0[5],
            p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9],
            ic[0], ic[1], ic[2], ic[3], ic[4], ic[5]
        );
    }
#endif /* #ifndef TESTING */

    gsl_odeiv2_driver_free(driver);

#ifndef TESTING

    /* Check that none of the outputs are negative */
    if (has_negative(ic, N_EQU))
    {
        warn(
            "WARNING: Negative equation outputs: \nfi = %.9lf, fa = %.9lf, fm = %.9lf, fs = %.9lf, fZ = %.9lf, fd = %.9lf ",
            ic[0], ic[1], ic[2], ic[3], ic[4], ic[5]
        );

        /* Set to 0 the fractions that are negative and renormalize */
        normalize_array(ic, N_EQU);
    }

#endif /* #ifndef TESTING */

    return status;
}

#ifndef TESTING

/*! \brief Compute the ionized and neutral fractions of the gas in a star-forming cell.
 *
 *  In AREPO a star-forming cell follows the two-phase ISM model of Springel & Hernquist
 *  (2003, https://doi.org/10.1046/j.1365-8711.2003.06206.x, SH03): cold clouds (T_c = TempClouds)
 *  with mass fraction x, embedded in a hot phase with mass fraction 1 - x.
 *
 *  \param[in] index Index of the gas cell in question.
 *  \param[out] fi_gas Ionized fraction of the gas [dimensionless].
 *  \param[out] fn_gas Neutral fraction of the gas [dimensionless].
 *
 *  \return void
 */
static void two_phase_fractions(const int index, double *fi_gas, double *fn_gas)
{
    /* Physical density [internal units] */
    double rho = SphP[index].Density * All.cf_a3inv;

    /* Electron abundance guess for the iteration in `GetCoolingTime()` */
    double ne = SphP[index].Ne;

    /**********************************************************************************************
    * Cold mass fraction x of SH03
    ***********************************************************************************************/

    double x, tsfr, factorEVP;
    calc_egyeff(index, rho, &ne, &x, &tsfr, &factorEVP);
    x = fmin(fmax(x, 0.0), 1.0);

    /**********************************************************************************************
    * Ionized hydrogen abundance of the hot phase
    ***********************************************************************************************/

    double nhp_hot;
    getHp(&nhp_hot);
    nhp_hot = fmin(fmax(nhp_hot, 0.0), 1.0);

    /**********************************************************************************************
    * Neutral and ionized fractions of the star-forming gas
    ***********************************************************************************************/

    /* Ionized fraction of the gas */
    /* Assuming the cold phase is fully neutral (T_c ~ 1e3 K) */
    *fi_gas = (1.0 - x) * nhp_hot;

    /* Neutral fraction of the gas */
    /* Assuming the cold phase is fully neutral (T_c ~ 1e3 K) */
    *fn_gas = x + (1.0 - x) * (1.0 - nhp_hot);
}

/*! \brief Compute the star formation rate.
 *
 *  Compute the gross star formation rate, i.e. the mass formed by psi during the integration,
 *
 *  \param[in] index Index of the gas cell in question.
 *
 *  \return The star formation rate in [Msun yr^(-1)].
 */
static double compute_gas_sfr(const int index)
{
    double integration_time = SphP[index].integration_time * 1000000; // [yr]
    double cell_mass = P[index].Mass * M_COSMO;                       // [Msun]

    /* Gross stellar fraction: fs(it) = (1 - R) * integral of psi, as R is constant during the integration */
    double fs_gross = SphP[index].ODE_fractions[3] / (1.0 - SphP[index].parameter_R);

    return fs_gross * cell_mass / integration_time;
}

/*! \brief Compute the star formation rate.
 *
 *  Compute the star formation rate solving the set of ODEs of our model.
 *
 *  \param[in] index Index of the gas cell in question.
 *
 *  \return The star formation rate in [Msun yr^(-1)].
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
     * Compute the integration time
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

    /**********************************************************************************************
     * Compute the ODE parameters
     *
     * rho_C:  Total cell density           [mp * cm^(-3)]
     * UVB:    UVB photoionization rate     [Myr^(-1)]
     * LWB:    LWB photodissociation rate   [Myr^(-1)]
     * eta_d:  Photodissociation efficiency [dimensionless]
     * eta_i:  Photoionization efficiency   [dimensionless]
     * R:      Mass recycling fraction      [dimensionless]
     * Zsn:    Metallicity of SN ejecta     [dimensionless]
     * h:      Cell size                    [cm]
     * fd_max: Maximum dust fraction        [dimensionless]
     * c_dg:   Dust growth constant         [Myr^(-1) * mp^(-1) * cm^3]
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

    /* Below the first age of the η tables, the photon rates Q are constant,
       so η = c * ∫ Q dt is proportional to the integration time */
    if (log_it < ETA_LOG_AGE_MIN)
    {
        double eta_scale = pow(10.0, log_it - ETA_LOG_AGE_MIN);
        eta_d *= eta_scale;
        eta_i *= eta_scale;
    }

    /* Mass recycling fraction [dimensionless] */
    double R = interpolate1D(Z, All.R_TABLE_DATA);

    /* Metallicity of the supernova ejecta [dimensionless] */
    double Zsn = interpolate1D(Z, All.ZSN_TABLE_DATA);

    /* Cell mass [mp] */
    double M = P[index].Mass * M_CGS;

    /* Cell volume [cm^3] */
    double V = M / rhoC;

    /* Cell size [cm] (diameter of the sphere with the same volume as the cell) */
    double h = 2.0 * cbrt(V / (4 * M_PI / 3));

    /* Indices of Si and C in MetalsFraction, looked up only once (they never change) */
    static int idx_Si = -1;
    static int idx_C  = -1;
    if (idx_Si < 0 || idx_C < 0)
    {
        idx_Si = element_index("Silicon");
        idx_C  = element_index("Carbon");
    }

    /* Silicon mass fraction [dimensionless] */
    double X_Si = fmax(0.0, SphP[index].MetalsFraction[idx_Si]);

    /* Carbon mass fraction [dimensionless] */
    double X_C = fmax(0.0, SphP[index].MetalsFraction[idx_C]);

    /* Clean the metal fractions from numerical noise */
    if (Z < Z_FLOOR)
    {
        /* Pristine gas: below the floor, the values are numerical noise */
        Z = 0.0;
        X_Si = 0.0;
        X_C = 0.0;
    }
    else if (X_Si + X_C > Z)
    {
        /* Si and C cannot be more than all the metals */
        double X_scale = Z / (X_Si + X_C);
        X_Si *= X_scale;
        X_C *= X_scale;
    }

    /* Maximum silicate and graphite mass fractions [dimensionless] */
    double y_Si = ODE_CMAXSI * X_Si;
    double y_C  = ODE_CMAXC * X_C;

    /* Maximum dust mass fraction [dimensionless] (capped at Z) */
    double fd_max = fmin(y_Si + y_C, Z);

    /* Dust growth constant [Myr^(-1) * mp^(-1) * cm^3] */
    double c_dg = 0.0;
    if (fd_max > 0.0)
    {
        c_dg = (ODE_CACCSI * X_Si * y_Si + ODE_CACCC * X_C * y_C) / (y_Si + y_C);
    }

    /* Maximum dust mass fraction for the ODEs [dimensionless] */
    double fd_max_ode = (fd_max > 0.0) ? fd_max : 1.0;

    /* Ionized and neutral fractions of the gas, from the two-phase ISM model of AREPO (SH03) */
    double fi_gas, fn_gas;
    two_phase_fractions(index, &fi_gas, &fn_gas);

    double parameters[] = {rhoC, UVB, LWB, eta_d, eta_i, R, Zsn, h, fd_max_ode, c_dg};

    /* Store the ODE parameters */
    SphP[index].tau_S            = 1.0 / (ODE_CS * sqrt(rhoC));
    SphP[index].parameter_a      = All.Time;
    SphP[index].parameter_z      = redshift;
    SphP[index].parameter_rhoC   = rhoC;
    SphP[index].parameter_UVB    = UVB;
    SphP[index].parameter_LWB    = LWB;
    SphP[index].parameter_Z      = Z;
    SphP[index].parameter_eta_d  = eta_d;
    SphP[index].parameter_eta_i  = eta_i;
    SphP[index].parameter_R      = R;
    SphP[index].parameter_Zsn    = Zsn;
    SphP[index].parameter_h      = h;
    SphP[index].parameter_fd_max = fd_max;
    SphP[index].parameter_c_dg   = c_dg;
    SphP[index].parameter_fi_gas = fi_gas;
    SphP[index].parameter_fn_gas = fn_gas;

    /**********************************************************************************************
     * Compute the initial conditions
     **********************************************************************************************/

    double fi, fa, fm, fs, fZ, fd;

    /* Ionized gas mass fraction [dimensionless] */
    fi = (1 - Z) * fi_gas;

    /* Atomic gas mass fraction [dimensionless] (all the neutral gas starts as atomic) */
    fa = (1 - Z) * fn_gas;

    /* Molecular gas mass fraction [dimensionless] */
    fm = 0.0;

    /* Stellar mass fraction [dimensionless] */
    fs = 0.0;

    /* Dust mass fraction [dimensionless] */
    double fZ_prev = SphP[index].ODE_fractions[4];
    double fd_prev = SphP[index].ODE_fractions[5];

    if (isnan(fZ_prev) || isnan(fd_prev) || fZ_prev + fd_prev <= 0.0)
    {
        /* First entry: Hirashita 2012, eq. 18, silicate + graphite, capped at Z */
    	fd = fmin(Z, ODE_CXSI * X_Si + ODE_CXC * X_C);
    }
    else
    {
        /* Later entries: dust-to-metal ratio of the previous integration, capped at fd_max */
        fd = fmin(fd_prev / (fZ_prev + fd_prev) * Z, fd_max);
    }

    /* Metal mass fraction [dimensionless] */
    fZ = Z - fd;

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

    /* Return the star formation rate in [Msun yr^(-1)] */
    return compute_gas_sfr(index);
}

#endif /* #ifndef TESTING */

#endif /* #ifdef EL_SFR */
