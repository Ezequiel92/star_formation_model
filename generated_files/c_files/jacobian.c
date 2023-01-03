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
	gsl_matrix *m = &dfdy_mat.matrix;
	
    /* Compute once operations that repeat in the Jacobian*/
    double aux_var = sqrt((1.0 - y[3]) * rho_C);

	gsl_matrix_set(m, 0, 0, -16.346836800000002 * y[0] * rho_C);
	gsl_matrix_set(m, 0, 1, 0);
	gsl_matrix_set(m, 0, 2, 0.0003885752566316025 * (eta_i + R) * aux_var);
	gsl_matrix_set(m, 0, 3, -0.00019428762831580126 * (1.0 / aux_var) * (eta_i + R) * y[2] * rho_C);

	gsl_matrix_set(m, 1, 0, 16.346836800000002 * y[0] * rho_C);
	gsl_matrix_set(m, 1, 1, 0.2826053731343284 * (1 - y[3]) * (-1.34e-5 - Z) * rho_C);
	gsl_matrix_set(m, 1, 2, 0.0003885752566316025 * (-1 * eta_i + eta_d) * aux_var);
	gsl_matrix_set(m, 1, 3, -0.2826053731343284 * (-1.34e-5 - Z) * y[1] * rho_C - 0.00019428762831580126 * (-1 * eta_i + eta_d) * (1.0 / aux_var) * y[2] * rho_C);

	gsl_matrix_set(m, 2, 0, 0);
	gsl_matrix_set(m, 2, 1, 0.2826053731343284 * (1.34e-5 + Z) * (1 - y[3]) * rho_C);
	gsl_matrix_set(m, 2, 2, 0.0003885752566316025 * (-1 - eta_d) * aux_var);
	gsl_matrix_set(m, 2, 3, -0.2826053731343284 * (1.34e-5 + Z) * y[1] * rho_C - 0.00019428762831580126 * (-1 - eta_d) * (1.0 / aux_var) * y[2] * rho_C);

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
