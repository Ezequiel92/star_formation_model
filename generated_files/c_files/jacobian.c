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
	
    double aux_var = sqrt((1e2 * y[1] + 1e5 * y[2]) * rho_C);

	gsl_matrix_set(m, 0, 0, -0.16409952000000003 * y[0] * rho_C);
	gsl_matrix_set(m, 0, 1, (0.03885752566316025 * (eta_i + R) * y[2] * rho_C) / (2 * aux_var));
	gsl_matrix_set(m, 0, 2, (38.857525663160246 * (eta_i + R) * y[2] * rho_C) / (2 * aux_var) + 0.0003885752566316025 * (eta_i + R) * aux_var);
	gsl_matrix_set(m, 0, 3, 0);

	gsl_matrix_set(m, 1, 0, 0.16409952000000003 * y[0] * rho_C);
	gsl_matrix_set(m, 1, 1, (0.03885752566316025 * (- eta_i + eta_d) * y[2] * rho_C) / (2 * aux_var) - 0.017393952755905513 * (1.27e-5 + Z) * y[1] * rho_C + 0.00017393952755905513 * (-100.0 * y[1] - 99999.99999999999 * y[2]) * (1.27e-5 + Z) * rho_C);
	gsl_matrix_set(m, 1, 2, (38.857525663160246 * (- eta_i + eta_d) * y[2] * rho_C) / (2 * aux_var) + 0.0003885752566316025 * (- eta_i + eta_d) * aux_var - 17.39395275590551 * (1.27e-5 + Z) * y[1] * rho_C);
	gsl_matrix_set(m, 1, 3, 0);

	gsl_matrix_set(m, 2, 0, 0);
	gsl_matrix_set(m, 2, 1, (0.03885752566316025 * (-1 - eta_d) * y[2] * rho_C) / (2 * aux_var) + 0.00017393952755905513 * (100.0 * y[1] + 99999.99999999999 * y[2]) * (1.27e-5 + Z) * rho_C + 0.017393952755905513 * (1.27e-5 + Z) * y[1] * rho_C);
	gsl_matrix_set(m, 2, 2, (38.857525663160246 * (-1 - eta_d) * y[2] * rho_C) / (2 * aux_var) + 0.0003885752566316025 * (-1 - eta_d) * aux_var + 17.39395275590551 * (1.27e-5 + Z) * y[1] * rho_C);
	gsl_matrix_set(m, 2, 3, 0);

	gsl_matrix_set(m, 3, 0, 0);
	gsl_matrix_set(m, 3, 1, (0.03885752566316025 * (1 - R) * y[2] * rho_C) / (2 * aux_var));
	gsl_matrix_set(m, 3, 2, (38.857525663160246 * (1 - R) * y[2] * rho_C) / (2 * aux_var) + 0.0003885752566316025 * (1 - R) * aux_var);
	gsl_matrix_set(m, 3, 3, 0);

	dfdt[0] = 0;
	dfdt[1] = 0;
	dfdt[2] = 0;
	dfdt[3] = 0;

	return GSL_SUCCESS;
}
