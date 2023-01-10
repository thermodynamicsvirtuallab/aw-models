#include "../definitions_and_headers.h"

double x_to_m ( double x, double xw );
double polynomial ( double x, double *coefs, int degree );
void fit_polynomial ( double *x_data, double *y_data,
		double *coefs, int size, int degree );
int aw_is_sane ( double aw, System *data);
void free_ptrs ( double *p0, double *p1, double *p2, double *p3, double *p4,
		double *p5, double *p6, double *p7, double *p8, double *p9);



/*
 * This function compares the water activities with the values
 * predicted using zdanovskii's relation
 */
void check_zdanovskii ( System *data, info *user_data, double *errors ) {

	int i, j, n, p, iter, degree_01, degree_02;
	double xw, aw, m_1, m_2, m_01, m_02, phi_calc, phi_real, S, m_std;
	double *m_1st_vec, *m_2nd_vec, *aw_vec;
	double *m_01_std, *m_02_std, *aw_01_std, *aw_02_std;
	double *K_aw_to_m_1st;
	double *K_aw_to_m_2nd;
	double *K_m_1st_to_aw;
	polynomial_function aw_to_m_1st, aw_to_m_2nd, m_1st_to_aw;

	aw_to_m_1st = &polynomial;
	aw_to_m_2nd = &polynomial;
	m_1st_to_aw = &polynomial;

	n = data->description.dataset_size;
	p = data->description.n_of_comps;

	/* allocating memory for molalities and water activities... */
		/* ...of data to be analyzed */
	m_1st_vec = malloc ( n * sizeof (double) );
	m_2nd_vec = malloc ( n * sizeof (double) );
	aw_vec = malloc ( n * sizeof (double) );
	data->x_and_aw.aw_calc = malloc ( n * sizeof (double) );

		/* ...of reference data... */
			/* ...of the first component */
	m_01_std = malloc ( data->x_and_aw.n_zdan[0] * sizeof (double) );
	aw_01_std = malloc ( data->x_and_aw.n_zdan[0] * sizeof (double) );

	for ( i = 0; i < data->x_and_aw.n_zdan[0]; i++ ) {
		m_01_std[i] = x_to_m ( data->x_and_aw.x_zdan[0][i],
				1 - data->x_and_aw.x_zdan[0][i] );
		aw_01_std[i] = data->x_and_aw.aw_zdan[0][i];
	}

			/* ...of the second component */
	m_02_std = malloc ( data->x_and_aw.n_zdan[1] * sizeof (double) );
	aw_02_std = malloc ( data->x_and_aw.n_zdan[1] * sizeof (double) );

	for ( i = 0; i < data->x_and_aw.n_zdan[1]; i++ ) {
		m_02_std[i] = x_to_m ( data->x_and_aw.x_zdan[1][i],
				1 - data->x_and_aw.x_zdan[1][i] );
		aw_02_std[i] = data->x_and_aw.aw_zdan[1][i];
	}

	/* converting data from mol fraction to molality */
	for ( i = 0; i < n; i++ ) {
		xw = 1;
		for ( j = 0; j < p; j++ ) {
			xw -= data->x_and_aw.x[i][j];
		}
		m_1st_vec[i] = x_to_m ( data->x_and_aw.x[i][0], xw );
		m_2nd_vec[i] = x_to_m ( data->x_and_aw.x[i][1], xw );
		aw_vec[i] = data->x_and_aw.aw[i];
	}

	/* get relations between molality and water activity */

	if ( data->x_and_aw.n_zdan[0] < 10 ) {
		degree_01 = (int) sqrt (data->x_and_aw.n_zdan[0] );
	} else {
		degree_01 = DEG_POLY_ZDAN;
	}

	if ( data->x_and_aw.n_zdan[1] < 10 ) {
		degree_02 = (int) sqrt (data->x_and_aw.n_zdan[0] );
	} else {
		degree_02 = DEG_POLY_ZDAN;
	}

	K_m_1st_to_aw = malloc ( degree_01 * sizeof (double) );
	K_aw_to_m_1st = malloc ( degree_01 * sizeof (double) );
	K_aw_to_m_2nd = malloc ( degree_02 * sizeof (double) );

	fit_polynomial ( m_01_std, aw_01_std, K_m_1st_to_aw,
			data->x_and_aw.n_zdan[0], degree_01 );
	fit_polynomial ( aw_01_std, m_01_std, K_aw_to_m_1st,
			data->x_and_aw.n_zdan[0], degree_01 );
	fit_polynomial ( aw_02_std, m_02_std, K_aw_to_m_2nd,
			data->x_and_aw.n_zdan[1], degree_02 );

	/* apply zdanovskii model */
	for ( i = 0; i < n; i++ ) {
		xw = 1;
		for ( j = 0; j < p; j++ ) {
			xw -= data->x_and_aw.x[i][j];
		}
		aw = xw; /* Using Raoult as initial guess */

		if ( aw_is_sane ( data->x_and_aw.aw[i], data ) != TRUE ) {
			fprintf ( stderr, "Water activity in mixture not in " );
			fprintf ( stderr, "acceptable range. Aborting...\n" );
			free_ptrs ( m_1st_vec, m_2nd_vec, aw_vec,
					m_01_std, aw_01_std, m_02_std,
					aw_02_std, K_m_1st_to_aw, K_aw_to_m_1st,
					K_aw_to_m_2nd );
			free (errors);
			finalize ( &data->description, &data->x_and_aw, user_data );
			exit (35);
		}

		m_01 = aw_to_m_1st ( aw, K_aw_to_m_1st, degree_01 );
		m_02 = aw_to_m_2nd ( aw, K_aw_to_m_2nd, degree_02 );
		m_1 = x_to_m ( data->x_and_aw.x[i][0], xw );
		m_2 = x_to_m ( data->x_and_aw.x[i][1], xw );
		S = ( m_1 / m_01 ) + ( m_2 / m_02 );

		iter = 0;
		while ( fabs ( S - 1 ) < TOL_ZDAN ) {
			m_01 = S * m_01;
			aw = m_1st_to_aw ( m_01, K_m_1st_to_aw, degree_01 );
			m_02 = aw_to_m_2nd ( aw, K_aw_to_m_2nd, degree_02 );
			S = ( m_1 / m_01 ) + ( m_2 / m_02 );
			iter ++;
			if ( iter >= MAX_ITER_ZDAN ) {
				fprintf ( stderr,
					"Maximum number of iterations (%d) ",
					MAX_ITER_ZDAN );
				fprintf ( stderr,
					"now surpassed. Aborting...\n" );
				free_ptrs ( m_1st_vec, m_2nd_vec, aw_vec,
					m_01_std, aw_01_std, m_02_std,
					aw_02_std, K_m_1st_to_aw, K_aw_to_m_1st,
					K_aw_to_m_2nd );
				free (errors);
				finalize ( &data->description,
						&data->x_and_aw, user_data );
				exit (39);
			}
		}

		aw = m_1st_to_aw ( m_01, K_m_1st_to_aw, degree_01 );
		data->x_and_aw.aw_calc[i] = aw;
		phi_calc = log (aw) / log (xw);
		if ( data->description.has_aw_data == TRUE ) {
			phi_real = log (data->x_and_aw.aw[i]) / log (xw);
		} else {
			m_std = ( data->x_and_aw.aw[i] ) /
				( ( 1 - data->x_and_aw.aw[i] ) * KGS_IN_MOL_WATER );
			aw = m_1st_to_aw ( m_std, K_m_1st_to_aw, degree_01 );
			phi_real = log (aw) / log (xw);
		}
		errors[i] = ( phi_calc - phi_real );
	}

	free_ptrs ( m_1st_vec, m_2nd_vec, aw_vec, m_01_std, aw_01_std, m_02_std,
			aw_02_std, K_m_1st_to_aw, K_aw_to_m_1st, K_aw_to_m_2nd );
}

void print_zdanovskii ( System *data, info *user_data, double *errors ) {

	int i, n;
	double R_squared, R_squared_aw;

	n = data->description.dataset_size;

	user_data->cost = 0;
	for ( i = 0; i < n; i++ ) {
		user_data->cost += pow ( errors[i], 2 );
	}

	user_data->cost = sqrt ( user_data->cost / n );

	R_squared = get_R_squared_check ( errors, data );
	R_squared_aw = get_R_squared_aw_check ( errors, data );

	if ( user_data->is_all == FALSE ) {
		fprintf ( stdout, "coeff. of determination      = %f\n",
				R_squared );
		fprintf ( stdout, "coeff. of determination (aw) = %f\n",
				R_squared_aw );
		fprintf ( stdout, "final cost:           |f(x)| = %f\n",
				user_data->cost );
	}

}

void save_zdanovskii ( System *data, info *user_data ) {

	int i, j, lines, comps;
	double xw, phi_calc, phi_exp;
	char *filename;
	FILE *results_file;

	filename = user_data->filename_new_results;
	results_file = fopen ( filename, "w" );

	lines = data->description.dataset_size;
	comps = data->description.n_of_comps;

	fprintf ( results_file, "phi_calc,phi_exp," );
	if ( user_data->aw_in_results == TRUE ) {
		fprintf ( results_file, "aw_calc,aw_exp," );
	}

	for ( i = 0; i < comps - 1; i++ ) {
		fprintf ( results_file, "%s,", data->description.components[i] );
	}
	fprintf ( results_file, "%s\n", data->description.components[comps-1] );

	for ( i = 0; i < lines; i++ ) {
		xw = 1;
		for ( j = 0; j < comps; j++ ) {
			xw -= data->x_and_aw.x[i][j];
		}
		phi_exp = log (data->x_and_aw.aw[i]) / log (xw);
		phi_calc = log (data->x_and_aw.aw_calc[i]) / log(xw);
		fprintf ( results_file, "%f,%f,", phi_calc, phi_exp );
		if ( user_data->aw_in_results == TRUE ) {
			fprintf ( results_file, "%f,%f,",
				data->x_and_aw.aw_calc[i], data->x_and_aw.aw[i] );
		}
		for ( j = 0; j < comps - 1; j++ ) {
			fprintf ( results_file, "%f,", data->x_and_aw.x[i][j] );
		}
		fprintf ( results_file, "%f\n", data->x_and_aw.x[i][comps-1] );
	}
	fclose (results_file);
}

/* this is the polynomial as needed for GSL non-linear fitting */
int polynomial_as_gsl ( const gsl_vector *K, void *params, gsl_vector *f ) {

	int size, degree, i, j;
	double real, calc;
	struct_for_zdan *data;

	data = (struct_for_zdan *) params;
	size = data->size;
	degree = data->degree;

	for ( i = 0; i < size; i++ ) {
		calc = 0;
		for ( j = 0; j < degree; j++ ) {
			calc += gsl_vector_get ( K, j ) *
				pow ( data->x[i], j );
		}
		real = data->y[i];
		gsl_vector_set ( f, i, real - calc );
	}

	return GSL_SUCCESS;
}

/* get polynomial coeficients from data  */
void fit_polynomial ( double *x_data, double *y_data,
		double *coefs, int size, int degree ) {

	int i, n, info;
	const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_workspace *w;
	gsl_multifit_nlinear_fdf fdf;
	gsl_multifit_nlinear_parameters fdf_params =
		gsl_multifit_nlinear_default_parameters ();
	gsl_vector_view x;
	callback_function callback;
	double *x_init;
	const double xtol = 1e-9;
	const double gtol = 1e-9;
	const double ftol = 0;
	struct_for_zdan data_for_func;

	data_for_func.size = size;
	data_for_func.degree = degree;
	data_for_func.x = x_data;
	data_for_func.y = y_data;

	x_init = malloc ( data_for_func.degree * sizeof (double) );

	n = size;
	fdf.f = polynomial_as_gsl;
	fdf.df = NULL;
	fdf.fvv = NULL;
	callback = NULL;
	fdf.n = n;
	fdf.p = data_for_func.degree;
	fdf.params = &data_for_func;

	x_init[0] = 1;
	for ( i = 1; i < data_for_func.degree; i++ ) {
		x_init[i] = 0;
	}
	x = gsl_vector_view_array ( x_init, data_for_func.degree );
	w = gsl_multifit_nlinear_alloc ( T, &fdf_params, n, data_for_func.degree );
	gsl_multifit_nlinear_winit ( &x.vector, NULL, &fdf, w );

	gsl_multifit_nlinear_driver ( 100, xtol, gtol, ftol,
			callback, NULL, &info, w);

	for ( i = 0; i < data_for_func.degree; i++ ) {
		coefs[i] = gsl_vector_get ( w->x, i );
	}

	gsl_multifit_nlinear_free (w);
	free (x_init);

}


/* auxiliary functions for property conversion */
double x_to_m ( double x, double xw ) {

	return x / ( xw * KGS_IN_MOL_WATER );

}

double polynomial ( double x, double *coefs, int degree ) {

	int i;
	double y = 0;

	for ( i = 0; i < degree; i++ ) {
		y += coefs[i] * pow ( x, i );
	}

	return y;
}

/*
 * This funtion checks if the binary mixture data are sufficient to
 * calculate the water activity of the ternary mixture.
 */
int aw_is_sane ( double aw, System *data ) {

	double aw_min;
	int n_1st, n_2nd, i, is_sane;

	n_1st = data->x_and_aw.n_zdan[0];
	n_2nd = data->x_and_aw.n_zdan[1];

	aw_min = 1.0;
	for ( i = 0; i < n_1st; i++ ) {
		if ( aw_min > data->x_and_aw.aw_zdan[0][i] ) {
			aw_min = data->x_and_aw.aw_zdan[0][i];
		}
	}
	for ( i = 0; i < n_2nd; i++ ) {
		if ( aw_min > data->x_and_aw.aw_zdan[1][i] ) {
			aw_min = data->x_and_aw.aw_zdan[1][i];
		}
	}

	if ( aw < aw_min ) {
		is_sane = FALSE;
	} else {
		is_sane = TRUE;
	}

	return is_sane;
}

/*
 * This function frees all pointers used in the zdanovskii's model.
 */

void free_ptrs ( double *p0, double *p1, double *p2, double *p3, double *p4,
		double *p5, double *p6, double *p7, double *p8, double *p9) {

	free (p0);
	free (p1);
	free (p2);
	free (p3);
	free (p4);
	free (p5);
	free (p6);
	free (p7);
	free (p8);
	free (p9);

}
