typedef struct InterpFunc2D
{
  gsl_spline2d *spline;
  gsl_interp_accel *xacc;
  gsl_interp_accel *yacc;
  double x_min;
  double x_max;
  double y_min;
  double y_max;
} InterpFunc2D;

typedef struct InterpFunc
{
  gsl_spline2d *spline;
  gsl_interp_accel *acc;
  double min;
  double max;
} InterpFunc2D;

typedef struct DataTable
{
  double *data;  // Values of the table
  int n_rows;    // Number of rows in the table
  int n_cols;    // Number of columns in the table
} DataTable;

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
double *read_ftable(const char *filepath, const int n_rows, const int n_cols)
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

// double *read_numbers_from_file(const char *file_path, int rows, int cols) {
//     FILE *fp = fopen(file_path, "r");
//     if (!fp) {
//         printf("Error: could not open file %s\n", file_path);
//         return NULL;
//     }

//     int count = rows * cols;
//     double *numbers = (double *) malloc(count * sizeof(double));
//     for (int i = 0; i < count; i++) {
//         fscanf(fp, "%lf", &numbers[i]);
//     }

//     fclose(fp);
//     return numbers;
// }


/*! \brief Evaluate a GSL interpolation function, using flat extrapolation.
 *
 *  Evaluate a bilinear interpolation function from the library GSL at a given value.
 *
 *  If the values are out of range, the function at its boundaries is used.
 *
 *  \param[in] x X coordinate at which the interpolation function will be evaluated.
 *  \param[in] y Y coordinate at which the interpolation function will be evaluated.
 *  \param[in] interpolation_function Interpolation function.
 *
 *  \return The result of evaluating the interpolation function.
 */
double eval_interpolation2D(double x, double y, const void *interpolation_function)
{
  double z;
  InterpFunc2D *interp = (InterpFunc2D *)interpolation_function;

  if(x < interp->x_min)
    {
      if(y < interp->y_min)
        {
          z = gsl_spline2d_eval(interp->spline, interp->x_min, interp->y_min, interp->xacc, interp->yacc);
        }
      else if(y > interp->y_max)
        {
          z = gsl_spline2d_eval(interp->spline, interp->x_min, interp->y_max, interp->xacc, interp->yacc);
        }
      else
        {
          z = gsl_spline2d_eval(interp->spline, interp->x_min, y, interp->xacc, interp->yacc);
        }
    }
  else if(x > interp->x_max)
    {
      if(y < interp->y_min)
        {
          z = gsl_spline2d_eval(interp->spline, interp->x_max, interp->y_min, interp->xacc, interp->yacc);
        }
      else if(y > interp->y_max)
        {
          z = gsl_spline2d_eval(interp->spline, interp->x_max, interp->y_max, interp->xacc, interp->yacc);
        }
      else
        {
          z = gsl_spline2d_eval(interp->spline, interp->x_max, y, interp->xacc, interp->yacc);
        }
    }
  else
    {
      if(y < interp->y_min)
        {
          z = gsl_spline2d_eval(interp->spline, x, interp->y_min, interp->xacc, interp->yacc);
        }
      else if(y > interp->y_max)
        {
          z = gsl_spline2d_eval(interp->spline, x, interp->y_max, interp->xacc, interp->yacc);
        }
      else
        {
          z = gsl_spline2d_eval(interp->spline, x, y, interp->xacc, interp->yacc);
        }
    }

  return z;
}

/*! \brief Construct a bilinear interpolation functions.
 *
 *  \param[in] table Path to the table containing the values of z for every x and y.
 *
 *  \return The interpolation function.
 */
void *get_interpolation_function(const char *table)
{
  double *data = read_ftable(table, NROWS, NCOLS);

  int nx = NROWS - 1;
  int ny = NCOLS - 1;

  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  gsl_spline2d *spline       = gsl_spline2d_alloc(T, nx, ny);
  gsl_interp_accel *xacc     = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc     = gsl_interp_accel_alloc();

  double *xa = malloc(nx * sizeof(double));
  double *ya = malloc(ny * sizeof(double));
  double *za = malloc(nx * ny * sizeof(double));

  for(size_t i = 0; i < nx; ++i)
    {
      xa[i] = data[(i + 1) * NCOLS];
    }

  for(size_t i = 0; i < ny; ++i)
    {
      ya[i] = data[i + 1];
    }

  for(size_t i = 0; i < nx; ++i)
    {
      for(size_t j = 0; j < ny; ++j)
        {
          gsl_spline2d_set(spline, za, i, j, data[(i + 1) * NCOLS + (j + 1)]);
        }
    }

  gsl_spline2d_init(spline, xa, ya, za, nx, ny);

  InterpFunc2D *interpolation_function = (InterpFunc2D *)malloc(sizeof(InterpFunc2D));

  interpolation_function->spline = spline;
  interpolation_function->xacc   = xacc;
  interpolation_function->yacc   = yacc;
  interpolation_function->x_min  = xa[0];
  interpolation_function->x_max  = xa[nx - 1];
  interpolation_function->y_min  = ya[0];
  interpolation_function->y_max  = ya[ny - 1];

  free(xa);
  free(ya);
  free(za);

  return (void *)interpolation_function;
}

/*! \brief Use linear interpolation to compute f(x).
 *
 *  If the value is out of range, the function at its boundaries is used.
 *
 *  \param[in] x Value at which the interpolation function will be evaluated.
 *  \param[in] table Path to the table containing the values of x for every y.
 *
 *  \return The interpolation function.
 */
double interpolate1D(double x, const char *table)
{
  double *data = read_ftable(table, NROWS, 2);
  double *xa   = malloc(NROWS * sizeof(double));
  double *ya   = malloc(NROWS * sizeof(double));

  for(size_t i = 0; i < NROWS; ++i)
    {
      xa[i] = eta_data[i * NCOLS];
    }

  for(size_t i = 0; i < NROWS; ++i)
    {
      ya[i] = eta_data[i * NCOLS + 1];
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
  else if(x > xa[NROWS - 1])
    {
      x      = xa[NROWS - 1];
      x1     = xa[NROWS - 2];
      x2     = xa[NROWS - 1];
      idx_x1 = NROWS - 2;
      idx_x2 = NROWS - 1;
    }
  else
    {
      for(size_t i = 1; i < NROWS; ++i)
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
 *  \param[in] table Path to the table containing the values of z for every x and y.
 *
 *  \return The interpolation function.
 */
double interpolate2D(double x, double y, const char *table)
{
  double *eta_data = read_ftable(table, NROWS, NCOLS);

  int nx = NROWS - 1;
  int ny = NCOLS - 1;

  double *xa = malloc(nx * sizeof(double));
  double *ya = malloc(ny * sizeof(double));
  double *za = malloc(nx * ny * sizeof(double));

  for(size_t i = 0; i < nx; ++i)
    {
      xa[i] = data[(i + 1) * NCOLS];
    }

  for(size_t i = 0; i < ny; ++i)
    {
      ya[i] = data[i + 1];
    }

  for(size_t i = 0; i < nx; ++i)
    {
      for(size_t j = 0; j < ny; ++j)
        {
          za[i * ny + j] = data[(i + 1) * NCOLS + (j + 1)];
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

/*! \brief Evaluates a n-dimensional function using linear interpolation.
 *
 *  Return the value of a function F at a given point x, using several known F(x_i).
 *
 *  \param[in] x Point at which the function is not known.
 *  \param[in] ode_table Array with the known values, plus metadata.
 *  \param[in] out_idx Index of the output column which will be interpolated, goes from 0 to number of outputs - 1.
 *  \param[in] n_dims Number of dimensions of the problem.
 *
 *  \return An approximation of F(x).
 */
double interpolate(const double *x, void *ode_table, const int out_idx, const int n_dims)
{
  DataTable *table = (DataTable *)ode_table;

  const double *data = table->data;    // Values in the interpolation table
  const int n_rows   = table->n_rows;  // Number of rows in the table
  const int n_cols   = table->n_cols;  // Number of columns in the table
  const int n_vert   = 1UL << n_dims;  // Number of vertices of a n-rectangle (= 2^n_dims)
  const int n_inp1   = n_dims + 1;     // Number of input parameters plus 1 (= number of columns in the reduced table)

  /***********************************************************************************************
   * Find the vertices of the n-rectangle surrounding the target point
   ***********************************************************************************************/

  double idx_low[n_dims];   // Grid indexes of the low bounding hyperplanes on each dimension
  double idx_high[n_dims];  // Grid indexes of the upper bounding hyperplanes in each dimension
  double val_low[n_dims];   // Input value corresponding to each low bounding hyperplanes in each dimension
  double val_high[n_dims];  // Input value corresponding to each upper bounding hyperplane in each dimension

  for(size_t i = 0; i < n_dims; ++i)
    {
      if(x[i] < data[i])
        {
          idx_low[i]  = 0;
          idx_high[i] = 0;
        }
      else if(x[i] > data[(n_rows - 1) * n_cols + i])
        {
          idx_low[i]  = N_GRID[i] - 1;
          idx_high[i] = N_GRID[i] - 1;
        }
      else
        {
          idx_low[i]  = floor((*INV_FUN[i])(x));
          idx_high[i] = ceil((*INV_FUN[i])(x));
        }

      val_low[i]  = (*FUN[i])(idx_low);
      val_high[i] = (*FUN[i])(idx_high);

      /* If the point coincides with a vertex, collapse that dimension */
      if(DEQU(x[i], val_low[i]) && DNEQ(x[i], val_high[i]))
        {
          val_high[i] = val_low[i];
          idx_high[i] = idx_low[i];
        }
      else if(DEQU(x[i], val_high[i]) && DNEQ(x[i], val_low[i]))
        {
          val_low[i] = val_high[i];
          idx_low[i] = idx_high[i];
        }
    }

  /*************************************************************************************************
   * Find the rows corresponding to the vertices
   ************************************************************************************************/

  double vertices[n_vert * n_inp1];  // Reduce table, with only the vertices and a single output parameter
  size_t vertex_idx[n_dims];         // Grid index in each dimension for the vertex
  size_t vertex_row;                 // Table row of the vertex

  for(size_t i = 0; i < n_vert; ++i)
    {
      size_t bit_i     = i;
      size_t dimension = 0;

      /* This runs up to the most significant bit set to 1 */
      while(bit_i > 0)
        {
          if(bit_i & 1)
            {
              /* If current bit is 1 */
              vertex_idx[dimension] = (size_t)idx_high[dimension];
            }
          else
            {
              /* If current bit is 0 */
              vertex_idx[dimension] = (size_t)idx_low[dimension];
            }

          bit_i = bit_i >> 1;
          dimension++;
        }

      /* The rest of the bits are 0, so use idx_low */
      for(; dimension < n_dims; ++dimension)
        {
          vertex_idx[dimension] = (size_t)idx_low[dimension];
        }

      /* Compute the row of the vertex */
      vertex_row = vertex_idx[n_dims - 1];
      for(size_t j = 0; j < n_dims - 1; ++j)
        {
          size_t mul_ngrid = 1;
          for(size_t k = j + 1; k < n_dims; ++k)
            {
              mul_ngrid *= N_GRID[k];
            }
          vertex_row += vertex_idx[j] * mul_ngrid;
        }

      /* Reduce the table to only have the vertices and a single output parameter */
      for(size_t j = 0; j < n_dims; ++j)
        {
          vertices[i * n_inp1 + j] = data[vertex_row * n_cols + j];
        }
      vertices[i * n_inp1 + n_dims] = data[vertex_row * n_cols + n_dims + out_idx];
    }

  /*************************************************************************************************
   * Compute the value for the target point using n-linear interpolation
   * See Zhang et al. (2021), https://doi.org/10.1145/3423184 and references therein
   ************************************************************************************************/

  double interp_val   = 0.0;  // Value for the target point
  double weight       = 1.0;  // Interpolation weight (fractional volume of the sub n-rectangles)
  size_t degen_factor = 0.0;  // Degeneracy factor (for when the target point lies in an edge or face of the n-rectangle)

  for(size_t i = 0; i < n_vert; ++i)
    {
      for(size_t j = 0; j < n_dims; ++j)
        {
          if(DNEQ(val_low[j], val_high[j]))
            {
              weight *= DEQU(vertices[i * n_inp1 + j], val_low[j]) ? (val_high[j] - x[j]) : (x[j] - val_low[j]);
              weight /= (val_high[j] - val_low[j]);
            }
          else
            {
              degen_factor++;
            }
        }

      interp_val += vertices[i * n_inp1 + n_dims] * weight / (0 | 1UL << degen_factor);
      weight       = 1.0;
      degen_factor = 0;
    }

  return interp_val;
}

/*! \brief Evaluate a GSL interpolation function, using flat extrapolation.
 *
 *  Evaluate a interpolation function from the library GSL at a given value.
 *
 *  If the value is out of range, the function at its boundaries is used.
 *
 *  \param[in] interp Interpolation function.
 *  \param[in] x Value at which the interpolation function will be evaluated.
 *
 *  \return The result of evaluating the interpolation function.
 */
static double eval_interp(const InterpFunc *interp, const double x)
{
  double y;

  if(x < interp->min)
    {
      y = gsl_spline_eval(interp->spline, interp->min, interp->acc);
    }
  else if(x > interp->max)
    {
      y = gsl_spline_eval(interp->spline, interp->max, interp->acc);
    }
  else
    {
      y = gsl_spline_eval(interp->spline, x, interp->acc);
    }

  return y;
}
