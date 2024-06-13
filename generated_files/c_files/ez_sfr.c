/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/ez_sfr/ez_sfr.c
 * \date        06/2024
 * \author      Ezequiel Lozano
 * \brief       Compute the star formation rate for a given gas cell.
 * \details     This file contains the routines to compute the star formation rate, according to our
 *              star formation model. The model evolves a set of ODEs which describe the mass
 *              exchange between the different phases of Hydrogen (ionized, atomic and molecular),
 *              and stars.
 *              contains functions:
 *                void *read_ftable(const char *file_path, const int n_rows, const int n_cols)
 *                static double interpolate1D(double x, const void *table)
 *                static double interpolate2D(double x, double y, const void *table)
 *                static int sf_ode(double t, const double y[], double f[], void *parameters)
 *                static int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *parameters)
 *                static void integrate_ode(const double *ic, double *parameters, double it, double *fractions)
 *                static double compute_sfr(const int index)
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
#include "./ez_sfr.h"
#else /* #ifdef TESTING */
#include "../allvars.h"
#include "../proto.h"
#endif /* #ifdef TESTING */

#ifdef EZ_SFR

/*! \brief Loads a text file into memory, as a C array.
 *
 *  Read a text file with space separated values, and store the numbers as doubles in a C array.
 *  Each line in the file must be at most 1024 characters long.
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
  if(!file_ptr)
    {
      fprintf(stderr, "Error: could not open file %s\n", file_path);
      return NULL;
    }

  int count    = n_rows * n_cols;
  double *data = (double *)malloc(count * sizeof(double));
  for(int i = 0; i < count; i++)
    {
      fscanf(file_ptr, "%lf", &data[i]);
    }

  fclose(file_ptr);

  data_table *table = malloc(sizeof(data_table));
  table->data       = data;
  table->n_rows     = n_rows;
  table->n_cols     = n_cols;

  return (void *)table;
}

/*! \brief Use linear interpolation to approximate f(x).
 *
 *  If the value is out of range, the boundaries of the function are used (constant extrapolation).
 *
 *  \param[in] x Value at which the interpolation function will be evaluated.
 *  \param[in] table_path Path to the table containing the values of f(xi), for many xi.
 *
 *  \return The interpolation function evaluated at `x`.
 */
#ifdef TESTING
double interpolate1D(double x, const char *table_path)
{
  void *table = read_ftable(table_path, R_NROWS, R_NCOLS);
#else  /* #ifdef TESTING */
static double interpolate1D(double x, const void *table)
{
#endif /* #ifdef TESTING */

  data_table *interp_table = (data_table *)table;
  double *data             = interp_table->data;
  int n_rows               = interp_table->n_rows;
  int n_cols               = interp_table->n_cols;

  double *xa = malloc(n_rows * sizeof(double));
  double *ya = malloc(n_rows * sizeof(double));

  for(size_t i = 0; i < n_rows; ++i)
    {
      xa[i] = data[i * n_cols];
    }

  for(size_t i = 0; i < n_rows; ++i)
    {
      ya[i] = data[i * n_cols + 1];
    }

  double x1, x2;
  int idx_x1, idx_x2;

  if(x < xa[0])
    {
      x      = xa[0];
      x1     = xa[0];
      x2     = xa[1];
      idx_x1 = 0;
      idx_x2 = 1;
    }
  else if(x > xa[n_rows - 1])
    {
      x      = xa[n_rows - 1];
      x1     = xa[n_rows - 2];
      x2     = xa[n_rows - 1];
      idx_x1 = n_rows - 2;
      idx_x2 = n_rows - 1;
    }
  else
    {
      for(size_t i = 1; i < n_rows; ++i)
        {
          if(xa[i - 1] <= x && x <= xa[i])
            {
              x1     = xa[i - 1];
              x2     = xa[i];
              idx_x1 = i - 1;
              idx_x2 = i;
              break;
            }
        }
    }

  double fQ1 = ya[idx_x1];
  double fQ2 = ya[idx_x2];

  double coeff   = 1 / (x2 - x1);
  double term_01 = fQ1 * (x2 - x);
  double term_02 = fQ2 * (x - x1);

  free(xa);
  free(ya);

  return coeff * (term_01 + term_02);
}

/*! \brief Use bilinear interpolation to approximate f(x,y).
 *
 *  If the values are out of range, the boundaries of the function are used (constant extrapolation).
 *
 *  \param[in] x X coordinate at which the interpolation function will be evaluated.
 *  \param[in] y Y coordinate at which the interpolation function will be evaluated.
 *  \param[in] table_path Path to the table containing the values of f(xi,yi), for many xi and yi.
 *
 *  \return The interpolation function evaluated at `x` and `y`.
 */
#ifdef TESTING
double interpolate2D(double x, double y, const char *table_path)
{
  void *table = read_ftable(table_path, ETA_NROWS, ETA_NCOLS);
#else  /* #ifdef TESTING */
static double interpolate2D(double x, double y, const void *table)
{
#endif /* #ifdef TESTING */

  data_table *interp_table = (data_table *)table;
  double *data             = interp_table->data;
  int n_rows               = interp_table->n_rows;
  int n_cols               = interp_table->n_cols;

  int nx = n_rows - 1;
  int ny = n_cols - 1;

  double *xa = malloc(nx * sizeof(double));
  double *ya = malloc(ny * sizeof(double));
  double *za = malloc(nx * ny * sizeof(double));

  for(size_t i = 0; i < nx; ++i)
    {
      xa[i] = data[(i + 1) * n_cols];
    }

  for(size_t i = 0; i < ny; ++i)
    {
      ya[i] = data[i + 1];
    }

  for(size_t i = 0; i < nx; ++i)
    {
      for(size_t j = 0; j < ny; ++j)
        {
          za[i * ny + j] = data[(i + 1) * n_cols + (j + 1)];
        }
    }

  double x1, x2, y1, y2;
  int idx_x1, idx_x2, idx_y1, idx_y2;

  if(x < xa[0])
    {
      x      = xa[0];
      x1     = xa[0];
      x2     = xa[1];
      idx_x1 = 0;
      idx_x2 = 1;
    }
  else if(x > xa[nx - 1])
    {
      x      = xa[nx - 1];
      x1     = xa[nx - 2];
      x2     = xa[nx - 1];
      idx_x1 = nx - 2;
      idx_x2 = nx - 1;
    }
  else
    {
      for(size_t i = 1; i < nx; ++i)
        {
          if(xa[i - 1] <= x && x <= xa[i])
            {
              x1     = xa[i - 1];
              x2     = xa[i];
              idx_x1 = i - 1;
              idx_x2 = i;
              break;
            }
        }
    }

  if(y < ya[0])
    {
      y      = ya[0];
      y1     = ya[0];
      y2     = ya[1];
      idx_y1 = 0;
      idx_y2 = 1;
    }
  else if(y > ya[ny - 1])
    {
      y      = ya[ny - 1];
      y1     = ya[ny - 2];
      y2     = ya[ny - 1];
      idx_y1 = ny - 2;
      idx_y2 = ny - 1;
    }
  else
    {
      for(size_t i = 1; i < ny; ++i)
        {
          if(ya[i - 1] <= y && y <= ya[i])
            {
              y1     = ya[i - 1];
              y2     = ya[i];
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

  double coeff   = 1 / ((x2 - x1) * (y2 - y1));
  double term_01 = fQ11 * (x2 - x) * (y2 - y);
  double term_02 = fQ21 * (x - x1) * (y2 - y);
  double term_03 = fQ12 * (x2 - x) * (y - y1);
  double term_04 = fQ22 * (x - x1) * (y - y1);

  free(xa);
  free(ya);
  free(za);

  return coeff * (term_01 + term_02 + term_03 + term_04);
}

/*! \brief Evaluate the systems of equations.
 *
 *  Evaluate the four ODEs of the model, using the following variables:
 *
 *  Ionized gas fraction:    fi(t) = Mi(t) / MC --> y[0]
 *  Atomic gas fraction:     fa(t) = Ma(t) / MC --> y[1]
 *  Molecular gas fraction:  fm(t) = Mm(t) / MC --> y[2]
 *  Stellar fraction:        fs(t) = Ms(t) / MC --> y[3]
 *
 *  where MC = Mi(t) + Ma(t) + Mm(t) + Ms(t) is the total density of the gas cell,
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

  /*
   * Destructure the parameters
   *
   * rho_C: Total cell density                                 [mp * cm⁻³]
   * Z:     Metallicity                                        [dimensionless]
   * eta_d: Photodissociation efficiency of Hydrogen molecules [dimensionless]
   * eta_i: Photoionization efficiency of Hydrogen atoms       [dimensionless]
   * R:     Mass recycling fraction                            [dimensionless]
   */
  double *p    = (double *)parameters;
  double rho_C = p[0];
  double Z     = p[1];
  double eta_d = p[2];
  double eta_i = p[3];
  double R     = p[4];

  /* Compute auxiliary equations */
  double tau_R           = ODE_CR / (y[0] * rho_C);
  double tau_C           = ODE_CC / (rho_C * (Z + ZEFF) * (1.0 - y[3]));
  double tau_S           = ODE_CS / sqrt(rho_C);
  double recombination   = y[0] / tau_R;
  double cloud_formation = y[1] / tau_C;
  double sfr             = y[2] / tau_S;

  /* Evaluate the ODE system */
  f[0] = -recombination + (eta_i + R) * sfr;
  f[1] = -cloud_formation + recombination + (eta_d - eta_i) * sfr;
  f[2] = cloud_formation - (1 + eta_d) * sfr;
  f[3] = (1 - R) * sfr;

  return GSL_SUCCESS;
}

/*! \brief Evaluate the Jacobian of the systems of equations.
*
*  Evaluate the Jacobian matrix of the model, using the following variables:
*
*  Ionized gas fraction:    fi(t) = Mi(t) / MC --> y[0]
*  Atomic gas fraction:     fa(t) = Ma(t) / MC --> y[1]
*  Molecular gas fraction:  fm(t) = Mm(t) / MC --> y[2]
*  Stellar fraction:        fs(t) = Ms(t) / MC --> y[3]
*
*  where MC = Mi(t) + Ma(t) + Mm(t) + Ms(t) is the total density of the gas cell,
*  and each equation has units of Myr^(-1).
*
*  \param[in] t Unused variable to comply with the `gsl_odeiv2_driver_alloc_y_new()` API.
*  \param[in] y Values of the variables at which the Jacobian will be evaluated.
*  \param[out] dfdy Where the results of evaluating the Jacobian will be stored.
*  \param[out] dfdt Where the results of evaluating the time derivatives will be stored.
*  \param[in] parameters Parameters for the Jacobian.
*
*  \return Constant `GSL_SUCCESS`, to confirm that the computation was successful.
*/
static int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *parameters)
{
	(void)(t);

    /*
	* Destructure the parameters
	*
	* rho_C: Total cell density [mp * cm⁻³]
	* Z:     Metallicity [dimensionless]
	* eta_d: Photodissociation efficiency of Hydrogen molecules [dimensionless]
	* eta_i: Photoionization efficiency of Hydrogen atoms [dimensionless]
	* R:     Mass recycling fraction [dimensionless]
	*/
	double *p    = (double *)parameters;
	double rho_C = p[0];
	double Z     = p[1];
	double eta_d = p[2];
	double eta_i = p[3];
	double R     = p[4];

	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 4, 4);
	gsl_matrix *m = &dfdy_mat.matrix;

    double aux_var = sqrt(rho_C);

	gsl_matrix_set(m, 0, 0, -16.409952 * y[0] * rho_C);
	gsl_matrix_set(m, 0, 1, 0);
	gsl_matrix_set(m, 0, 2, 0.019428762831580126 * (eta_i + R) * aux_var);
	gsl_matrix_set(m, 0, 3, 0);

	gsl_matrix_set(m, 1, 0, 16.409952 * y[0] * rho_C);
	gsl_matrix_set(m, 1, 1, -17.393952755905513 * (1.0 - y[3]) * rho_C * (1.27e-5 + Z));
	gsl_matrix_set(m, 1, 2, 0.019428762831580126 * (eta_d - eta_i) * aux_var);
	gsl_matrix_set(m, 1, 3, 17.393952755905513 * y[1] * rho_C * (1.27e-5 + Z));

	gsl_matrix_set(m, 2, 0, 0);
	gsl_matrix_set(m, 2, 1, 17.393952755905513 * (1.0 - y[3]) * rho_C * (1.27e-5 + Z));
	gsl_matrix_set(m, 2, 2, -0.019428762831580126 * (1.0 + eta_d) * aux_var);
	gsl_matrix_set(m, 2, 3, -17.393952755905513 * y[1] * rho_C * (1.27e-5 + Z));

	gsl_matrix_set(m, 3, 0, 0);
	gsl_matrix_set(m, 3, 1, 0);
	gsl_matrix_set(m, 3, 2, 0.019428762831580126 * (1.0 - R) * aux_var);
	gsl_matrix_set(m, 3, 3, 0);

	dfdt[0] = 0;
	dfdt[1] = 0;
	dfdt[2] = 0;
	dfdt[3] = 0;

	return GSL_SUCCESS;
}

/*! \brief Solve the system of ODEs, using numerical integration.
 *
 * ICs:
 *
 * fi: Ionized gas fraction   [dimensionless]
 * fa: Atomic gas fraction    [dimensionless]
 * fm: Molecular gas fraction [dimensionless]
 * fs: Stellar fraction       [dimensionless]
 *
 * Parameters
 *
 * rho_C: Total cell density                                 [mp * cm⁻³]
 * Z:     Metallicity                                        [dimensionless]
 * eta_d: Photodissociation efficiency of Hydrogen molecules [dimensionless]
 * eta_i: Photoionization efficiency of Hydrogen atoms       [dimensionless]
 * R:     Mass recycling fraction                            [dimensionless]
 *
 *  \param[in] ic Initial conditions.
 *  \param[in] parameters Parameters for the ODEs.
 *  \param[in] it Integration time in Myr.
 *  \param[out] fractions Array to store the resulting fractions after integration.
 *
 *  \return void.
 */
#ifdef TESTING
void integrate_ode(const double *ic, double *parameters, double it, double *fractions)
#else  /* #ifdef TESTING */
static void integrate_ode(const double *ic, double *parameters, double it, double *fractions)
#endif /* #ifdef TESTING */
{
  /* Initial conditions */
  double fi = ic[0];
  double fa = ic[1];
  double fm = ic[2];
  double fs = ic[3];

#ifndef TESTING
  /* Check that none of the inputs are too negative */
  if(fi < -1e-8 || fa < -1e-8 || fm < -1e-8 || fs < -1e-8)
    {
      terminate("ERROR: Negative equation inputs: \nfi = %.9lf, fa = %.9lf, fm = %.9lf, fs = %.9lf", fi, fa, fm, fs);
    }
#endif /* #ifndef TESTING */

  /* Set to 0 the negative fractions */
  fi           = fmax(0.0, fi);
  fa           = fmax(0.0, fa);
  fm           = fmax(0.0, fm);
  fs           = fmax(0.0, fs);
  double total = fi + fa + fm + fs;

  /* Initialize integration variables */
  double t0  = 0.0;
  double i_f = 0.0;
  double a_f = 0.0;
  double m_f = 0.0;
  double s_f = 0.0;

  /* Save mean cell density */
  double rhoC = parameters[0];

  gsl_odeiv2_system sys     = {sf_ode, jacobian, N_EQU, parameters};
  gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_bsimp, it * 1e-3, 1e-8, 0.0);

  for(size_t i = 0; i < DIVISIONS; ++i)
    {
      fractions[0] = fi / total;
      fractions[1] = fa / total;
      fractions[2] = fm / total;
      fractions[3] = fs / total;

      parameters[0] = rhoC * F_RHO[i];

      gsl_odeiv2_driver_reset(driver);
      int status = gsl_odeiv2_driver_apply(driver, &t0, it, fractions);

      i_f += fractions[0] * PDF[i];
      a_f += fractions[1] * PDF[i];
      m_f += fractions[2] * PDF[i];
      s_f += fractions[3] * PDF[i];

      t0 = 0.0;

#ifndef TESTING
      if(status != GSL_SUCCESS)
        {
          terminate("GSL ERROR: Could not integrate, status = %d\n", status);
        }

#endif /* #ifndef TESTING */
    }

  gsl_odeiv2_driver_free(driver);

#ifndef TESTING
  /* Check that none of the outputs are too negative */
  if(i_f < -1e-8 || a_f < -1e-8 || m_f < -1e-8 || s_f < -1e-8)
    {
      warn("WARNING: Negative equation outputs: \nfi = %.9lf, fa = %.9lf, fm = %.9lf, fs = %.9lf", i_f, a_f, m_f, s_f);
    }
#endif /* #ifndef TESTING */

  /* Set to 0 the negative fractions */
  fractions[0] = fmax(0.0, i_f);
  fractions[1] = fmax(0.0, a_f);
  fractions[2] = fmax(0.0, m_f);
  fractions[3] = fmax(0.0, s_f);

  /* Renormalize IC */
  total = fractions[0] + fractions[1] + fractions[2] + fractions[3];
  fractions[0] /= total;
  fractions[1] /= total;
  fractions[2] /= total;
  fractions[3] /= total;
}

#ifndef TESTING

/*! \brief Compute the star formation rate.
 *
 *  Compute the star formation rate using the star mass fraction of the cell.
 *
 *  \param[in] index Index of the gas cell in question.
 *
 *  \return The star formation rate in [Mₒ yr^(-1)].
 */
static double compute_sfr(const int index)
{
  double integration_time = SphP[index].accu_integration_time * 1000000;  // [yr]
  double cell_mass        = P[index].Mass * M_COSMO;                      // [Mₒ]

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
  if(!P[index].TimeBinHydro)
    {
      return 0.0;
    }
  if(All.ComovingIntegrationOn && SphP[index].Density < All.OverDensThresh)
    {
      return 0.0;
    }
  if(!All.ComovingIntegrationOn && SphP[index].Density * All.cf_a3inv < All.PhysDensThresh)
    {
      return 0.0;
    }

  /*************************************************************************************************
   * Compute the time parameters
   *************************************************************************************************/

  /* Current time [Myr] */
  double current_time = All.Time * T_MYR;

  /* Delta time [Myr] */
  double delta_time = 0.0;
  if(!isnan(SphP[index].current_time))
    {
      delta_time = current_time - SphP[index].current_time;
      if(delta_time == 0.0)
        {
          return compute_sfr(index);
        }
    }

  /* Integration time [Myr] */
  double integration_time = (((integertime)1) << P[index].TimeBinHydro) * All.Timebase_interval;
  if(integration_time <= 0.0)
    {
      return 0.0;
    }
  integration_time *= T_MYR * All.cf_atime / All.cf_time_hubble_a;

  /* Accumulated integration time [Myr] */
  double accu_integration_time = integration_time;
  if(!isnan(SphP[index].accu_integration_time) && delta_time < SphP[index].tau_S)
    {
      accu_integration_time += SphP[index].accu_integration_time;
    }

  /* Store the integration time parameters */
  SphP[index].current_time          = current_time;
  SphP[index].delta_time            = delta_time;
  SphP[index].integration_time      = integration_time;
  SphP[index].accu_integration_time = accu_integration_time;

  /*************************************************************************************************
   * Compute the ODE parameters
   *************************************************************************************************/

  /* Cell density [cm⁻³] (only hydrogen) */
  double rhoC = SphP[index].Density * RHO_COSMO * HYDROGEN_MASSFRAC;

  /* Metallicity [dimensionless] */
  double Z = (SphP[index].Metallicity > 0.0) ? SphP[index].Metallicity : 0.0;

  /* Photodissociation efficiency [dimensionless] */
  double eta_d = interpolate2D(accu_integration_time, Z, All.ETA_D_TABLE_DATA);

  /* Photoionization efficiency [dimensionless] */
  double eta_i = interpolate2D(accu_integration_time, Z, All.ETA_I_TABLE_DATA);

  /* Mass recycling fraction [dimensionless] */
  double R = interpolate1D(Z, All.R_TABLE_DATA);

  double parameters[] = {rhoC, Z, eta_d, eta_i, R};

  /* Store the ODE parameters */
  SphP[index].parameter_rhoC  = rhoC;
  SphP[index].parameter_Z     = Z;
  SphP[index].parameter_eta_d = eta_d;
  SphP[index].parameter_eta_i = eta_i;
  SphP[index].parameter_R     = R;

  /*************************************************************************************************
   * Compute the initial conditions
   *************************************************************************************************/

  double fi, fa, fm, fs;

  if(!isnan(SphP[index].ODE_fractions[0]) && delta_time < SphP[index].tau_S)
    {
      /* Ionized gas mass fraction [dimensionless] */
      fi = SphP[index].ODE_fractions[0];
      /* Atomic gas mass fraction [dimensionless] */
      fa = SphP[index].ODE_fractions[1];
      /* Molecular gas mass fraction [dimensionless] */
      fm = SphP[index].ODE_fractions[2];
      /* Stellar mass fraction [dimensionless] */
      fs = SphP[index].ODE_fractions[3];
    }
  else
    {
      double nhp, nh;
      get_arepo_fraction(index, &nhp, &nh);
      /* Ionized gas mass fraction [dimensionless] */
      fi = nhp / (nhp + nh);
      /* Atomic gas mass fraction [dimensionless] */
      fa = 1.0 - fi;
      /* Molecular gas mass fraction [dimensionless] */
      fm = 0.0;
      /* Stellar mass fraction [dimensionless] */
      fs = 0.0;
    }

  const double ic[] = {fi, fa, fm, fs};

  /*************************************************************************************************
   * Integrate the ODEs
   *************************************************************************************************/

  double fractions[4] = {0.0};
  integrate_ode(ic, parameters, integration_time, fractions);

  /* Store the star formation time scale (τ_star) */
  SphP[index].tau_S = ODE_CS / sqrt((1 - fractions[3]) * rhoC);

  /* Store the results after solving the ODEs */
  for(size_t i = 0; i < N_EQU; ++i)
    {
      SphP[index].ODE_fractions[i] = fractions[i];
    }

  /* Return the star formation rate in [Mₒ yr^(-1)] */
  return compute_sfr(index);
}

#endif /* #ifndef TESTING */

#endif /* #ifdef EZ_SFR */
