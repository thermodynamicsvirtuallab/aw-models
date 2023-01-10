#include "../definitions_and_headers.h"

/* This function is an implementation of Norrish's model for predicting water
 * activity in solutions, to be used for fitting.
 */

int phi_norrish ( const gsl_vector *K, void *params, gsl_vector * f ) {

	System *data;
	double phi_i_calc, phi_i_real;
	double xw, sumxiki;
	int n, c;
	int i, j;

	data = (System *) params;

	n = data->description.dataset_size;
	c = data->description.n_of_comps;

	for ( i = 0; i < n; i++ ) {
		xw = 1;
		sumxiki = 0;
		for ( j = 0; j < c; j++ ) {
			xw -= data->x_and_aw.x[i][j];
			sumxiki += sqrt ( fabs ( gsl_vector_get ( K, j ) ) )
				* data->x_and_aw.x[i][j];
		}
		sumxiki = sumxiki * sumxiki;
		xw = log (xw);
		phi_i_calc = ( xw + sumxiki ) / xw;
		if ( data->description.has_aw_data == TRUE ) {
			phi_i_real = log ( data->x_and_aw.aw[i] ) / xw;
		} else {
			phi_i_real = gsl_vector_get ( K, 0 ) * data->x_and_aw.aw[i];
			phi_i_real = phi_i_real * phi_i_real;
			phi_i_real += 1 - data->x_and_aw.aw[i];
			phi_i_real = phi_i_real / xw;
		}
		gsl_vector_set ( f, i, phi_i_calc - phi_i_real );
	}

	return GSL_SUCCESS;
}


/* function to call every iteration */
void callback_norrish ( const size_t iter, void *params,
		const gsl_multifit_nlinear_workspace *w ) {

	gsl_vector *f = gsl_multifit_nlinear_residual (w);
	gsl_vector *x = gsl_multifit_nlinear_position (w);
	size_t size, i;

	fprintf ( stderr, "Iteration number %2zu:\n", iter );
	fprintf ( stderr, "\t|f(x)| = %.4f\n", gsl_blas_dnrm2 (f) / f->size );
	size = x->size;
	for ( i = 0; i < size; i++ ) {
		fprintf ( stderr, "\tK_%d = %.4f\n", (int) i,
				gsl_vector_get ( x, i ) );
	}
	fprintf ( stderr, "\n" );

}

void print_norrish ( gsl_matrix *covar, gsl_multifit_nlinear_workspace *w,
		int status, double chisq0, double chisq,
		System *data, info *user_data ) {

	int p, n;
	size_t i;
	double correction, R_squared, R_squared_aw;

	p = data->description.n_of_comps;
	n = data->description.dataset_size;

	correction = GSL_MAX_DBL ( 1, sqrt ( chisq / (n - p) ) );

	if ( user_data->is_all == FALSE ) {
		fprintf ( stdout, "Results Obtained:\n" );
		for ( i = 0; i < p; i++ ) {
			fprintf ( stdout, "\tK_%d = %.5e\t+/-\t%.5e\t", (int) i,
				gsl_vector_get ( w->x, i), correction *
				sqrt ( gsl_matrix_get ( covar, i, i ) ) );
			fprintf ( stdout, "(%s)\n",
					data->description.components[i] );
			/* add name of the solute */
		}

		R_squared = get_R_squared ( w, data );
		R_squared_aw = get_R_squared_aw ( w, data );

		fprintf ( stdout, "initial cost:         |f(x)| = %f\n",
				sqrt ( chisq0 / n ) );
		fprintf ( stdout, "final cost:           |f(x)| = %f\n",
				sqrt ( chisq / n ) );
		fprintf ( stdout, "adj. coeff. of determination = %f\n",
				1 - ( 1 - R_squared ) * ( n - 1 ) / ( n - p - 1 ) );
		fprintf ( stdout, "coeff. of determination (aw) = %f\n",
				1 - ( 1 - R_squared_aw ) * ( n - 1 ) / ( n - p - 1 ) );
		fprintf ( stdout, "Exit status is \"%s\".\n\n",
				gsl_strerror (status) );
	}

	user_data->cost = sqrt ( chisq / n );

}

void save_norrish ( System *data, info *user_data,
		gsl_multifit_nlinear_workspace *w ) {

	int i, j, lines, comps;
	double xw, sumxiki, phi_calc, phi_exp;
	char *filename;
	FILE *results_file;
	gsl_vector *x = gsl_multifit_nlinear_position (w);

	lines = data->description.dataset_size;
	comps = data->description.n_of_comps;

	filename = user_data->filename_new_results;
	results_file = fopen ( filename, "w" );

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
		sumxiki = 0;
		for ( j = 0; j < comps; j++ ) {
			xw -= data->x_and_aw.x[i][j];
			sumxiki += data->x_and_aw.x[i][j] *
				sqrt ( fabs ( gsl_vector_get ( x, j ) ) );
		}
		sumxiki = sumxiki * sumxiki;
		xw = log(xw);
		phi_calc = ( xw + sumxiki) / xw;
		phi_exp = log ( data->x_and_aw.aw[i] ) / xw;
		fprintf ( results_file, "%f,%f,", phi_calc, phi_exp );
		if ( user_data->aw_in_results == TRUE ) {
			fprintf ( results_file, "%f,%f,",
				exp ( xw + sumxiki ), data->x_and_aw.aw[i] );
		}
		for ( j = 0; j < comps - 1; j++ ) {
			fprintf ( results_file, "%f,", data->x_and_aw.x[i][j] );
		}
		fprintf ( results_file, "%f\n", data->x_and_aw.x[i][comps-1] );
	}

	fclose (results_file);

}
