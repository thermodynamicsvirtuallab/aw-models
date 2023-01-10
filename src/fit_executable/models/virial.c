#include "../definitions_and_headers.h"

/* This function is an implementation of a virial model for predicting water
 * activity in solutions, to be used for fitting.
 */

int phi_virial ( const gsl_vector *K, void *params, gsl_vector * f ) {

	System *data;
	double phi_i_calc, phi_i_real;
	double xw, sumxiki, ln_aw_mono, m_mono;
	int n, c_mono;
	int i, j, k, counter;

	data = (System *) params;

	n = data->description.dataset_size;
	c_mono = data->description.n_of_comps;

	for ( i = 0; i < n; i++ ) {
		xw = 1;
		sumxiki = 0;
		counter = 0;
		for ( j = 0; j < c_mono; j++ ) {
			xw -= data->x_and_aw.x[i][j];
		}
		for ( j = 0; j < c_mono; j++ ) {
			sumxiki += ( data->x_and_aw.x[i][j]
					/ xw ) *
				( 1 + 2 * gsl_vector_get ( K, counter ) );
			counter++;
			for ( k = 0; k < j; k++ ) {
				sumxiki += 3 * ( data->x_and_aw.x[i][j] *
					data->x_and_aw.x[i][k] /
					pow ( xw, 2 ) ) *
					gsl_vector_get ( K, counter );
				counter++;
			}
		}
		xw = log(xw);
		phi_i_calc = sumxiki / xw;
		if ( data->description.has_aw_data == TRUE ) {
			phi_i_real = log ( data->x_and_aw.aw[i] ) / xw;
		} else {
			m_mono = data->x_and_aw.aw[i] /
				( ( 1 - data->x_and_aw.aw[i] ) * KGS_IN_MOL_WATER );
			ln_aw_mono = m_mono * ( 1 + 2 *
					gsl_vector_get ( K, 0 ) );
			phi_i_real = ln_aw_mono / xw;
		}      /*
			* This is weird, but instead of graciously
			* informing the water activity of a n-ary
			* solution, or some equivalent property as
			* BPE or phi or Pvap, as one would expect,
			* a few studies inform the molality of a
			* binary solution with equal activity. Thus
			* we can't actually obtain a phi_real, but
			* being the composition of the binary mixture
			* used as reference a subset of the n-ary
			* mixture of interest, the model can be, very
			* easily, adapted. Crisis, therefore, averted.
			*/
		gsl_vector_set ( f, i, phi_i_calc - phi_i_real );
	}

	return GSL_SUCCESS;
}

/* function to call every iteration */
void callback_virial ( const size_t iter, void *params,
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

void print_virial ( gsl_matrix *covar, gsl_multifit_nlinear_workspace *w,
		int status, double chisq0, double chisq,
		System *data, info *user_data ) {

	int p, n;
	size_t i, j, counter;
	double correction, R_squared, R_squared_aw;

	p = data->description.n_of_comps;
	n = data->description.dataset_size;

	correction = GSL_MAX_DBL ( 1, sqrt ( chisq / (n - p) ) );

	if ( user_data->is_all == FALSE ) {
		fprintf ( stdout, "Results Obtained:\n" );

		counter = 0;

		for ( i = 0; i < p; i++ ) {
			fprintf ( stdout, "\tb_%d = %.5e\t+/-\t%.5e\t", (int) counter,
				gsl_vector_get ( w->x, counter),
				correction * sqrt ( gsl_matrix_get
					( covar, counter, counter ) ) );
			fprintf ( stdout, "(%s)\n", data->description.components[i] );
			/* add solute name */
			counter++;
			for ( j = 0; j < i; j++ ) {
				fprintf ( stdout, "\tc_%d = %.5e\t+/-\t%.5e\t",
					(int) counter, gsl_vector_get ( w->x, counter),
					correction * sqrt ( gsl_matrix_get
						( covar, counter, counter ) ) );
				fprintf ( stdout, "(%s/%s)\n",
					data->description.components[i],
					data->description.components[j] );
				/* add solute pair name */
				counter++;
			}
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

void save_virial ( System *data, info *user_data,
		gsl_multifit_nlinear_workspace *w ) {

	int i, j, k, lines, comps, counter;
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
		counter = 0;
		for ( j = 0; j < comps; j++ ) {
			xw -= data->x_and_aw.x[i][j];
			sumxiki += ( data->x_and_aw.x[i][j]
					/ xw ) *
				( 1 + 2 * gsl_vector_get ( x, counter ) );
			counter++;
			for ( k = 0; k < j; k++ ) {
				sumxiki += 3 * ( data->x_and_aw.x[i][j] *
					data->x_and_aw.x[i][k] /
					pow ( xw, 2 ) ) *
					gsl_vector_get ( x, counter );
				counter++;
			}
		}
		xw = log(xw);
		phi_calc = sumxiki / xw;
		phi_exp = log ( data->x_and_aw.aw[i] ) / xw;
		fprintf ( results_file, "%f,%f,", phi_calc, phi_exp );
		if ( user_data->aw_in_results == TRUE ) {
			fprintf ( results_file, "%f,%f,",
				exp (sumxiki), data->x_and_aw.aw[i] );
		}
		for ( j = 0; j < comps - 1; j++ ) {
			fprintf ( results_file, "%f,", data->x_and_aw.x[i][j] );
		}
		fprintf ( results_file, "%f\n", data->x_and_aw.x[i][comps-1] );
	}

	fclose (results_file);
}

