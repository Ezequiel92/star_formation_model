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

    double i_f = fmax(0.0, y[0]);
    double a_f = fmax(0.0, y[1]);
    double m_f = fmax(0.0, y[2]);
    double s_f = fmax(0.0, y[3]);
	
    double aux_var = sqrt((1.0 - s_f) * rho_C);

	gsl_matrix_set(m, 0, 0, -16.409952 * i_f * rho_C);
	gsl_matrix_set(m, 0, 1, 0);
	gsl_matrix_set(m, 0, 2, 0.0003885752566316025 * (eta_i + R) * aux_var);
	gsl_matrix_set(m, 0, 3, (-0.0003885752566316025 * (eta_i + R) * m_f * rho_C) / (2 * aux_var));

	gsl_matrix_set(m, 1, 0, 16.409952 * i_f * rho_C);
	gsl_matrix_set(m, 1, 1, 0.17393952755905515 * (- a_f - m_f) * (1.27e-5 + Z) * rho_C - 0.17393952755905515 * (1.27e-5 + Z) * a_f * rho_C);
	gsl_matrix_set(m, 1, 2, 0.0003885752566316025 * (- eta_i + eta_d) * aux_var - 0.17393952755905515 * (1.27e-5 + Z) * a_f * rho_C);
	gsl_matrix_set(m, 1, 3, (-0.0003885752566316025 * (- eta_i + eta_d) * m_f * rho_C) / (2 * aux_var));

	gsl_matrix_set(m, 2, 0, 0);
	gsl_matrix_set(m, 2, 1, 0.17393952755905515 * (a_f + m_f) * (1.27e-5 + Z) * rho_C + 0.17393952755905515 * (1.27e-5 + Z) * a_f * rho_C);
	gsl_matrix_set(m, 2, 2, 0.0003885752566316025 * (-1 - eta_d) * aux_var + 0.17393952755905515 * (1.27e-5 + Z) * a_f * rho_C);
	gsl_matrix_set(m, 2, 3, (-0.0003885752566316025 * (-1 - eta_d) * m_f * rho_C) / (2 * aux_var));

	gsl_matrix_set(m, 3, 0, 0);
	gsl_matrix_set(m, 3, 1, 0);
	gsl_matrix_set(m, 3, 2, 0.0003885752566316025 * (1 - R) * aux_var);
	gsl_matrix_set(m, 3, 3, (-0.0003885752566316025 * (1 - R) * m_f * rho_C) / (2 * aux_var));

	dfdt[0] = 0;
	dfdt[1] = 0;
	dfdt[2] = 0;
	dfdt[3] = 0;

	return GSL_SUCCESS;
}
