/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *
 * \file        src/ez_sfr/ez_sfr.c
 * \date        12/2022
 * \author      Ezequiel Lozano
 * \brief       Compute the star formation rate for a given gas cell.
 * \details     This file contains the routines to compute the star formation rate, according to our
 *              star formation model. The model evolves a set of ODEs which describe the mass
 *              exchange between the different phases of neutral Hydrogen (atomic and molecular),
 *              and stars.
 *              contains functions:
 *                static double *read_ftable(const char *filepath, const int n_rows, const int n_cols)
 *                double interpolate(double x, double y, const void *interpolation_function)
 *                void *get_eta_interp_function(const char *eta_table)
 *                static int sf_ode(double t, const double y[], double f[], void *ode_params)
 *                static int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *ode_params)
 *                void integrate_ode(const double *ic, double *parameters, double it, double *fractions)
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

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline2d.h>
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
 *  Read a text file with space separated values, and stores the numerical values as
 *  doubles in a C array. Each line in the file must be at most 1024 characters long.
 *
 *  The value on row i and column j will be stored as element i * n_cols + j in the array.
 *
 *  \param[in] filepath Path to the file.
 *  \param[in] n_rows Number of rows in the table.
 *  \param[in] n_cols Number of columns in the table.
 *
 *  \return Values in the table as an array.
 */
#ifdef TESTING
double *read_ftable(const char *filepath, const int n_rows, const int n_cols)
#else  /* #ifdef TESTING */
static double *read_ftable(const char *filepath, const int n_rows, const int n_cols)
#endif /* #ifdef TESTING */
{
  FILE *file_ptr = fopen(filepath, "r");
  double *data   = malloc(n_rows * n_cols * sizeof(double));
  char line[1024];
  char *scan;
  int offset;

  for(size_t i = 0; i < n_rows; ++i)
    {
      fgets(line, sizeof line, file_ptr);
      scan   = line;
      offset = 0;

      for(size_t j = 0; j < n_cols; ++j)
        {
          sscanf(scan, "%lf%n", &(data[i * n_cols + j]), &offset);
          scan += offset;
        }
    }

  fclose(file_ptr);

  return data;
}

/*! \brief Use linear interpolation to compute f(x).
 *
 *  If the value is out of range, the function at its boundaries is used.
 *
 *  \param[in] x Value at which the interpolation function will be evaluated.
 *  \param[in] R_table Path to the table containing the values of x for every y.
 *
 *  \return The interpolation function.
 */
#ifdef TESTING
double interpolate1D(const double x, const char *R_table)
{
  double *R_data = read_ftable(table, R_NROWS, R_NCOLS);
#else  /* #ifdef TESTING */
static double interpolate1D(const double x, const double *R_data)
{
#endif /* #ifdef TESTING */

  double *xa = malloc(R_NROWS * sizeof(double));
  double *ya = malloc(R_NROWS * sizeof(double));

  for(size_t i = 0; i < R_NROWS; ++i)
    {
      xa[i] = R_data[i * R_NCOLS];
    }

  for(size_t i = 0; i < R_NROWS; ++i)
    {
      ya[i] = R_data[i * R_NCOLS + 1];
    }

  double x1, x2;
  int idx_x1, idx_x2;

  if(x < xa[0])
    {
      x1     = xa[0];
      x2     = xa[1];
      idx_x1 = 0;
      idx_x2 = 1;
    }
  else if(x > xa[R_NROWS - 1])
    {
      x1     = xa[R_NROWS - 2];
      x2     = xa[R_NROWS - 1];
      idx_x1 = R_NROWS - 2;
      idx_x2 = R_NROWS - 1;
    }
  else
    {
      for(size_t i = 1; i < R_NROWS; ++i)
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

/*! \brief Use bilinear interpolation to compute f(x,y).
 *
 *  If the values are out of range, the function at its boundaries is used.
 *
 *  \param[in] x X coordinate at which the interpolation function will be evaluated.
 *  \param[in] y Y coordinate at which the interpolation function will be evaluated.
 *  \param[in] eta_table Path to the table containing the values of z for every x and y.
 *
 *  \return The interpolation function.
 */
#ifdef TESTING
double interpolate2D(const double x, const double y, const char *eta_table)
{
  double *eta_data = read_ftable(eta_table, ETA_NROWS, ETA_NCOLS);
#else  /* #ifdef TESTING */
static double interpolate2D(const double x, const double y, const double *eta_data)
{
#endif /* #ifdef TESTING */

  int nx = ETA_NROWS - 1;
  int ny = ETA_NCOLS - 1;

  double *xa = malloc(nx * sizeof(double));
  double *ya = malloc(ny * sizeof(double));
  double *za = malloc(nx * ny * sizeof(double));

  for(size_t i = 0; i < nx; ++i)
    {
      xa[i] = eta_data[(i + 1) * ETA_NCOLS];
    }

  for(size_t i = 0; i < ny; ++i)
    {
      ya[i] = eta_data[i + 1];
    }

  for(size_t i = 0; i < nx; ++i)
    {
      for(size_t j = 0; j < ny; ++j)
        {
          za[i * ny + j] = eta_data[(i + 1) * ETA_NCOLS + (j + 1)];
        }
    }

  double x1, x2, y1, y2;
  int idx_x1, idx_x2, idx_y1, idx_y2;

  if(x < xa[0])
    {
      x1     = xa[0];
      x2     = xa[1];
      idx_x1 = 0;
      idx_x2 = 1;
    }
  else if(x > xa[nx - 1])
    {
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
      y1     = ya[0];
      y2     = ya[1];
      idx_y1 = 0;
      idx_y2 = 1;
    }
  else if(y > ya[ny - 1])
    {
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

/*! \brief Evaluate the systems of equations for the star formation model.
 *
 *  Evaluate the four ODEs of the model, using the following variables:
 *
 *  Ionized gas fraction:    if(t) = ρi(t) / ρC --> y[0]
 *  Atomic gas fraction:     af(t) = ρa(t) / ρC --> y[1]
 *  Molecular gas fraction:  mf(t) = ρm(t) / ρC --> y[2]
 *  Stellar fraction:        sf(t) = ρs(t) / ρC --> y[3]
 *
 *  where ρC = ρi(t) + ρa(t) + ρm(t) + ρs(t), and each equation has units of Myr^(-1).
 *
 *  \param[in] t Unused variable to comply with the gsl_odeiv2_driver_alloc_y_new() API.
 *  \param[in] y Values of the variables at which the ODEs will be evaluated.
 *  \param[out] f Where the results of evaluating the ODEs will be stored.
 *  \param[in] ode_params Parameters for the ODEs.
 *
 *  \return Constant GSL_SUCCESS, to confirm that the computation was successful.
 */
static int sf_ode(double t, const double y[], double f[], void *ode_params)
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
  double *parameters = (double *)ode_params;
  double rho_C       = parameters[0];
  double Z           = parameters[1];
  double eta_d       = parameters[2];
  double eta_i       = parameters[3];
  double R           = parameters[4];

  /* Compute auxiliary equations */
  double tau_R           = ODE_CR / (y[0] * rho_C);
  double tau_C           = ODE_CC / ((1 - y[3]) * rho_C * (Z + ZEFF));
  double tau_S           = ODE_CS / sqrt((1 - y[3]) * rho_C);
  double recombination   = y[0] / tau_R;
  double cloud_formation = y[1] / tau_C;
  double sfr             = (AW * y[1] + MW * y[2]) / tau_S;

  /* Evaluate ODE system */
  f[0] = -recombination + (eta_i + R) * sfr;
  f[1] = -cloud_formation + recombination + (eta_d - eta_i) * sfr;
  f[2] = cloud_formation - (1 + eta_d) * sfr;
  f[3] = (1 - R) * sfr;

  return GSL_SUCCESS;
}

/*! \brief Evaluate the Jacobian matrix of the star formation model.
 *
 *  Evaluate the Jacobian matrix of the model, using the following variables:
 *
 *  Ionized gas fraction:    if(t) = ρi(t) / ρC --> y[0]
 *  Atomic gas fraction:     af(t) = ρa(t) / ρC --> y[1]
 *  Molecular gas fraction:  mf(t) = ρm(t) / ρC --> y[2]
 *  Stellar fraction:        sf(t) = ρs(t) / ρC --> y[3]
 *
 *  where ρC = ρi(t) + ρa(t) + ρm(t) + ρs(t), and each equation has units of Myr^(-1).
 *
 *  \param[in] t Unused variable to conform to the gsl_odeiv2_driver_alloc_y_new() API.
 *  \param[in] y Values of the variables at which the Jacobian will be evaluated.
 *  \param[out] dfdy Where the results of evaluating the Jacobian will be stored.
 *  \param[out] dfdt Where the result of evaluating the time derivatives will be stored.
 *  \param[in] ode_params Parameters for the Jacobian.
 *
 *  \return Constant GSL_SUCCESS, to confirm that the computation was successful.
 */
static int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *ode_params)
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
  double *parameters = (double *)ode_params;
  double rho_C       = parameters[0];
  double Z           = parameters[1];
  double eta_d       = parameters[2];
  double eta_i       = parameters[3];
  double R           = parameters[4];

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 4, 4);
  gsl_matrix *m            = &dfdy_mat.matrix;

  /* Compute once operations that repeat in the Jacobian*/
  double aux_var = sqrt((1.0 - y[3]) * rho_C);

  gsl_matrix_set(m, 0, 0, -16.346836800000002 * y[0] * rho_C);
  gsl_matrix_set(m, 0, 1, 0);
  gsl_matrix_set(m, 0, 2, 0.0003885752566316025 * (eta_i + R) * aux_var);
  gsl_matrix_set(m, 0, 3, -0.00019428762831580126 * (1.0 / aux_var) * (eta_i + R) * y[2] * rho_C);

  gsl_matrix_set(m, 1, 0, 16.346836800000002 * y[0] * rho_C);
  gsl_matrix_set(m, 1, 1, 0.2826053731343284 * (1 - y[3]) * (-1.34e-5 - Z) * rho_C);
  gsl_matrix_set(m, 1, 2, 0.0003885752566316025 * (-1 * eta_i + eta_d) * aux_var);
  gsl_matrix_set(m, 1, 3,
                 -0.2826053731343284 * (-1.34e-5 - Z) * y[1] * rho_C -
                     0.00019428762831580126 * (-1 * eta_i + eta_d) * (1.0 / aux_var) * y[2] * rho_C);

  gsl_matrix_set(m, 2, 0, 0);
  gsl_matrix_set(m, 2, 1, 0.2826053731343284 * (1.34e-5 + Z) * (1 - y[3]) * rho_C);
  gsl_matrix_set(m, 2, 2, 0.0003885752566316025 * (-1 - eta_d) * aux_var);
  gsl_matrix_set(
      m, 2, 3,
      -0.2826053731343284 * (1.34e-5 + Z) * y[1] * rho_C - 0.00019428762831580126 * (-1 - eta_d) * (1.0 / aux_var) * y[2] * rho_C);

  gsl_matrix_set(m, 3, 0, 0);
  gsl_matrix_set(m, 3, 1, 0);
  gsl_matrix_set(m, 3, 2, 0.0003885752566316025 * (1 - R) * aux_var);
  gsl_matrix_set(m, 3, 3, -0.00019428762831580126 * (1 - R) * (1.0 / aux_var) * y[2] * rho_C);

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
 * i0: Ionized gas fraction (of the total cell density) [dimensionless]
 * a0: Atomic gas fraction (of the total cell density) [dimensionless]
 * m0: Molecular gas fraction (of the total cell density) [dimensionless]
 * s0: Stellar fraction (of the total cell density) [dimensionless]
 *
 * Parameters
 *
 * rho_C: Total cell density [mp * cm⁻³]
 * Z:     Metallicity [dimensionless]
 * eta_d: Photodissociation efficiency of Hydrogen molecules [dimensionless]
 * eta_i: Photoionization efficiency of Hydrogen atoms [dimensionless]
 * R:     Mass recycling fraction [dimensionless]
 *
 *  \param[in] ic Initial conditions.
 *  \param[in] parameters Parameters for the ODEs.
 *  \param[in] it Integration time in Myr.
 *  \param[out] fractions Array to store the resulting fraction after integration.
 *
 *  \return void.
 */
#ifdef TESTING
int integrate_ode(const double *ic, double *parameters, double it, double *fractions)
#else  /* #ifdef TESTING */
static int integrate_ode(const double *ic, double *parameters, double it, double *fractions)
#endif /* #ifdef TESTING */
{
  /* Initial conditions */
  double i0 = ic[0];
  double a0 = ic[1];
  double m0 = ic[2];
  double s0 = ic[3];

  double rho = parameters[0];

  /* Initialize integration variables */
  double t0  = 0.0;
  double i_f = 0.0;
  double a_f = 0.0;
  double m_f = 0.0;
  double s_f = 0.0;

  gsl_odeiv2_system sys     = {sf_ode, jacobian, 4, parameters};
  gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_msbdf, it * 1e-3, ABS_TOL, REL_TOL);

  /* Use the density PDF to increase the effective resolution of the computation */
  for(size_t i = 0; i < DIVISIONS; ++i)
    {
      double y[]    = {i0, a0, m0, s0};
      parameters[0] = rho * F_RHO[i];

      gsl_odeiv2_driver_reset(driver);
      gsl_odeiv2_driver_apply(driver, &t0, it, y);

      i_f += y[0] * PDF[i];
      a_f += y[1] * PDF[i];
      m_f += y[2] * PDF[i];
      s_f += y[3] * PDF[i];

      t0 = 0.0;
    }

  gsl_odeiv2_driver_free(driver);

  fractions[0] = i_f;
  fractions[1] = a_f;
  fractions[2] = m_f;
  fractions[3] = s_f;

  return 0;
}

#ifndef TESTING

/*! \brief Compute the star formation rate.
 *
 *  The SFR will be used to calculate the probability of forming a star population.
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
   * Compute the integration time
   *************************************************************************************************/

  double current_time = All.Time * T_MYR;
  double delta_time   = current_time - SphP[index].C_TIME;

  /* Integration time [im Myr)] */
  double it = (((integertime)1) << P[index].TimeBinHydro) * All.Timebase_interval * T_MYR * All.cf_atime / All.cf_time_hubble_a;
  if(it <= 0.0)
    {
      return 0.0;
    }
  double accu_it = it;
  if(!isnan(SphP[index].ACCU_IT) && delta_time < SphP[index].TAU_S)
    { /* If this is not the first time for this particle */
      accu_it = SphP[index].ACCU_IT + it;
    }

  /* Store the time parameters */
  SphP[index].IT      = it;
  SphP[index].ACCU_IT = accu_it;
  SphP[index].C_TIME  = current_time;
  SphP[index].D_TIME  = delta_time;

  /*************************************************************************************************
   * Compute the parameters
   *************************************************************************************************/

  /* Density in cm⁻³ */
  double density = SphP[index].Density * RHO_COSMO;

  /* Metallicity [dimensionless] */
  double metallicity = 0.0;
  if(SphP[index].Metallicity > 0.0)
    {
      metallicity = SphP[index].Metallicity;
    }

  /* Compute the photodissociation efficiency */
  double eta_d = interpolate2D(accu_it, metallicity, ETA_D_DATA);
  double eta_i = interpolate2D(accu_it, metallicity, ETA_I_DATA);

  /* Mass recycling fraction */
  double R = interpolate1D(metallicity, R_DATA);

  double parameters[] = {density, metallicity, eta_d, eta_i, R};

  /* Store the parameters */
  SphP[index].PAR_rho_C = density;
  SphP[index].PAR_Z     = metallicity;
  SphP[index].PAR_eta_d = eta_d;
  SphP[index].PAR_eta_i = eta_i;
  SphP[index].PAR_R     = R;

  /*************************************************************************************************
   * Compute the ICs
   *************************************************************************************************/

  /* Atomic gas fraction [dimensionless] */
  double a_f = 0.0;
  get_neutral_fraction(index, &af);
  if(!isnan(SphP[index].atomicFraction) && delta_time < SphP[index].TAU_S)
    { /* If this is not the first time for this particle */
      a_f = SphP[index].atomicFraction;
    }

  /* Ionized gas fraction */
  double i_f = 1.0 - a_f;
  get_neutral_fraction(index, &af);
  if(!isnan(SphP[index].ionizedFraction) && delta_time < SphP[index].TAU_S)
    { /* If this is not the first time for this particle */
      i_f = SphP[index].ionizedFraction;
    }

  /* Fraction of the neutral density that is molecular Hydrogen [dimensionless] */
  double m_f = 0.0;
  if(!isnan(SphP[index].molecularFraction) && delta_time < SphP[index].TAU_S)
    { /* If this is not the first time for this particle */
      m_f = SphP[index].molecularFraction;
    }

  const double ic[] = {i_f, a_f, m_f, 1.0 - i_f - a_f - m_f};

  /* Store stellar time parameter */
  SphP[index].TAU_S = ODE_CS / sqrt((1 - fractions[3]) * density);

  double fractions[4] = {0.0};
  integrate_ode(ic, parameters, it, fractions);

  /* Store the results after solving the ODEs */
  SphP[index].ionizedFraction   = fractions[0];
  SphP[index].atomicFraction    = fractions[1];
  SphP[index].molecularFraction = fractions[2];
  SphP[index].stellarFraction   = fractions[3];

  /* Return the star formation rate in [Mₒ yr^(-1)] */
  return fractions[3] * (P[index].Mass * M_COSMO) / accu_it;
}

#endif /* #ifndef TESTING */

#endif /* #ifdef EZ_SFR */