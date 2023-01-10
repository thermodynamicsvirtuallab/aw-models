#include "../definitions_and_headers.h"

/*
 * This function implements a simplified version of the UNIQUAC model,
 * assuming that for all component i the respective parameters r_i and
 * q_i are equal, and that for all pairs of components i and j, the
 * parameter u_ij simplifies to the geometric mean between u_ii and u_jj.
 *
 * You may ask why the complete model isn't used: it's because the data
 * sets are, usually, small, and it isn't reasonable to fit a model with
 * 18 or so parameters to a data set with 10 or 12 points.
 *
 * You may also ask why these simplifications, specifically, were chosen.
 * r_i = q_i appears, according to the original Abrams & Prausnitz 1975
 * paper, to be one of the most commonly used simplifications. In fact,
 * r_i and q_i are strongly correlated (both are a measure of molecular
 * "size"), therefore this simplication is sane. u_ij = sqrt( u_ii * u_jj )
 * is another simplification mentioned in the article, specially atractive
 * because it reduces the quadratic dependence of the number of parameters
 * to fit with the number of components to a linear one, better suited for
 * the small data sets used.
 *
 * One could also question our choice of simplification for u_ij, pointing
 * that, in the original paper, is also mentioned u_ij = u_ii. We opted not
 * to assert this, to maintain the dependence of u_ij in the two components.
 */

int phi_uniquac ( const gsl_vector *K, void *params, gsl_vector * f ) {

	System *data;
	long double phi_i_calc, phi_i_real;
	long double ln_gamma_w;
	long double x_w, x_j, x_k;
	long double sumxjlj, sumqjxj, sumthetajtaujw, sumthetaktaukj, sumsum;
	long double l_w, l_j, r_w, r_j, q_w, q_j, q_k;
	long double theta_w, theta_j, theta_k, Phi_w;
	long double u_jj, u_jw, u_wj, u_kj, u_kk, u_ww, tau_jw, tau_kj, tau_wj;
	int n, p;
	int i, j, k;

	data = (System *) params;

	n = data->description.dataset_size;
	p = data->description.n_of_comps;

	for ( i = 0; i < n; i++ ) {
		x_w = 1;
		for ( j = 0; j < p; j++ ) {
			x_w -= data->x_and_aw.x[i][j];
		}      /* In this model the water is seen as another component;
			* therefore, we also need its molar fraction.
			*/

		r_w = fabs ( gsl_vector_get ( K, 0 ) );
		l_w = 1 - r_w;
			/*
			* l = (z/2) * (r-q) - (r-1)
			* :. l = -(r-1) = 1 - r
			*/
		sumxjlj = x_w * l_w;

		q_w = r_w;
		sumqjxj = q_w * x_w;

		for ( j = 0; j < p; j++ ) {

			x_j = data->x_and_aw.x[i][j];
			r_j = fabs ( gsl_vector_get ( K, j + 1 ) );
			l_j = 1 - r_j;
			sumxjlj += x_j * l_j;
			/* here we dealt with the l's and their sum */

			q_j = r_j;
			sumqjxj += q_j * x_j;
			/* needed for theta_j */
		}

		theta_w = ( q_w * x_w ) / sumqjxj;
		Phi_w = theta_w; /* because q_i = r_i */

		sumthetajtaujw = theta_w; /* = theta_w * tau_ww = theta_w * 1 */
		for ( j = 0; j < p; j++ ) {
			q_j = fabs ( gsl_vector_get ( K, j + 1 ) );
			x_j = data->x_and_aw.x[i][j];
			theta_j = ( q_j * x_j ) / sumqjxj;
			u_ww = fabs ( gsl_vector_get ( K, p + 1 ) );
			u_jj = fabs ( gsl_vector_get ( K, p + j + 2 ) );
			u_jw = sqrt ( u_jj * u_ww );
			tau_jw = exp ( - ( u_jw - u_ww ) / ( R * TEMP ) );
			sumthetajtaujw += theta_j * tau_jw;
		}

		sumsum = 0;
		for ( j = -1; j < p; j++ ) {

			sumthetaktaukj = theta_w;
			for ( k = -1; k < p; k++ ) {
				if ( k == -1 ) {
					x_k = x_w;
				} else {
					x_k = data->x_and_aw.x[i][k];
				}
				q_k = fabs ( gsl_vector_get ( K, k + 1 ) );
				theta_k = ( q_k * x_k ) / sumqjxj;
				u_kk = fabs ( gsl_vector_get ( K, p + k + 2 ) );
				u_jj = fabs ( gsl_vector_get ( K, p + j + 2 ) );
				u_kj = sqrt ( u_jj * u_kk );
				tau_kj = exp ( - ( u_kj - u_jj ) / ( R * TEMP ) );
				sumthetaktaukj += theta_k * tau_kj;
			}

			q_j = fabs ( gsl_vector_get ( K, j + 1 ) );
			if ( j == -1 ) {
				x_j = x_w;
			} else {
				x_j = data->x_and_aw.x[i][j];
			}
			theta_j = ( q_j * x_j ) / sumqjxj;
			u_ww = fabs ( gsl_vector_get ( K, p ) );
			u_jj = fabs ( gsl_vector_get ( K, p + j + 2 ) );
			u_wj = sqrt ( u_ww * u_jj );
			tau_wj = exp ( - ( u_wj - u_jj ) / ( R * TEMP ) );

			sumsum += ( theta_j * tau_wj ) / sumthetaktaukj;

		}

		ln_gamma_w = log ( Phi_w / x_w );
			/* "ln_gamma_w +=
			* ( Z_COORD / 2 ) * q_w * log ( theta_w / Phi_w );"
			* is a expression that would be here if we did not
			* assume r_j = q_j; however, because we assumed this,
			* the expression between parenthesis inside the log
			* equals  one, and this forces the whole right-hand
			* side to 0.
			*/
		ln_gamma_w += l_w;
		ln_gamma_w -= ( Phi_w / x_w ) * sumxjlj;
		ln_gamma_w -= q_w * ( 1 - log ( sumthetajtaujw ) - sumsum );
		phi_i_calc = ( ln_gamma_w + log (x_w) ) / log (x_w);

		if ( data->description.has_aw_data == TRUE ) {
			phi_i_real = log (data->x_and_aw.aw[i]) / log (x_w);
		} else {
			/*
			* This is exactly the algorithm programmed above,
			* for a binary mixture. For the reasons why we need
			* implement this ugly thing, see the file ./virial.c,
			* specifically the function "phi_virial".
			*/
			x_j = data->x_and_aw.aw[i];
			x_w = 1 - x_j;

			r_w = fabs ( gsl_vector_get ( K, 0 ) );
			q_w = r_w;
			l_w = 1 - r_w;

			r_j = fabs ( gsl_vector_get ( K, 1 ) );
			q_j = r_j;
			l_j = 1 - r_j;

			sumxjlj = ( x_w * l_w ) + ( x_j * l_j );
			sumqjxj = ( q_w * x_w ) + ( q_j * x_j );
			theta_w = ( q_w * x_w ) / sumqjxj;
			Phi_w = theta_w;

			sumthetajtaujw = theta_w;
			theta_j = ( q_j * x_j ) / sumqjxj;
			u_ww = fabs ( gsl_vector_get ( K, p + 1 ) );
			u_jj = fabs ( gsl_vector_get ( K, p + 2 ) );
			u_jw = sqrt ( u_jj * u_ww );
			tau_jw = exp ( - ( u_jw - u_ww ) / ( R * TEMP ) );
			sumthetajtaujw += theta_j * tau_jw;

			sumsum = 0;
			for ( j = -1; j < 1; j++ ) {

				sumthetaktaukj = theta_w;
				for ( k = -1; k < 1; k++ ) {
					if ( k == -1 ) {
						x_k = x_w;
					} else {
						x_k = data->x_and_aw.aw[i];
					}
					q_k = fabs ( gsl_vector_get ( K, k + 1 ) );
					theta_k = ( q_k * x_k ) / sumqjxj;
					u_kk = fabs
						( gsl_vector_get ( K, p + k + 2 ) );
					u_jj = fabs
						( gsl_vector_get ( K, p + j + 2 ) );
					u_kj = sqrt ( u_jj * u_kk );
					tau_kj = exp
						( - ( u_kj - u_jj ) / ( R * TEMP ) );
					sumthetaktaukj += theta_k * tau_kj;
				}

				q_j = fabs ( gsl_vector_get ( K, j + 1 ) );
				if ( j == -1 ) {
					x_j = x_w;
				} else {
					x_j = data->x_and_aw.x[i][j];
				}
				theta_j = ( q_j * x_j ) / sumqjxj;
				u_ww = fabs ( gsl_vector_get ( K, p ) );
				u_jj = fabs ( gsl_vector_get ( K, p + j + 2 ) );
				u_wj = sqrt ( u_ww * u_jj );
				tau_wj = exp ( - ( u_wj - u_jj ) / ( R * TEMP ) );

				sumsum += ( theta_j * tau_wj ) / sumthetaktaukj;

			}


			ln_gamma_w = log ( Phi_w / x_w );
			ln_gamma_w += l_w;
			ln_gamma_w -= ( Phi_w / x_w ) * sumxjlj;
			ln_gamma_w -= q_w * ( 1 - log ( sumthetajtaujw ) - sumsum );

			phi_i_real = ( ln_gamma_w +
				log ( 1 - data->x_and_aw.aw[i] ) ) / log (x_w);
		}
		gsl_vector_set ( f, i, phi_i_calc - phi_i_real );
	}

	return GSL_SUCCESS;
}

/* function to call every iteration */
void callback_uniquac ( const size_t iter, void *params,
		const gsl_multifit_nlinear_workspace *w ) {

	gsl_vector *f = gsl_multifit_nlinear_residual (w);
	gsl_vector *x = gsl_multifit_nlinear_position (w);
	size_t size, i;

	fprintf ( stderr, "Iteration number %2zu:\n", iter );
	fprintf ( stderr, "\t|f(x)| = %.4f\n", gsl_blas_dnrm2 (f) / f->size );
	size = x->size;
	for ( i = 0; i < size; i++ ) {
		fprintf ( stderr, "\tK_%d = %.7e\n", (int) i,
				fabs ( gsl_vector_get ( x, i ) ) );
	}
	fprintf ( stderr, "\n" );

}

void print_uniquac ( gsl_matrix *covar, gsl_multifit_nlinear_workspace *w,
		int status, double chisq0, double chisq,
		System *data, info *user_data ) {

	int p, n;
	int i;
	double correction, R_squared, R_squared_aw;

	p = data->description.n_of_comps;
	n = data->description.dataset_size;

	correction = GSL_MAX_DBL ( 1, sqrt ( chisq / (n - p) ) );

	if ( user_data->is_all == FALSE ) {
		fprintf ( stdout, "Results Obtained:\n" );

		for ( i = -1; i < p; i++ ) {
			fprintf ( stdout, "\tq_%d  = %.5e\t+/-\t%.5e\t",
				i + 1, fabs ( gsl_vector_get ( w->x, i + 1) ),
				correction * sqrt ( gsl_matrix_get
					( covar, i + 1, i + 1 ) ) );
			if ( i != -1 ) {
				fprintf ( stdout, "(%s)\n",
					data->description.components[i] );
			} else {
				fprintf ( stdout, "(water)\n" );
			}
			/* add solute name */
		}
		for ( i = -1; i < p; i++ ) {
			fprintf ( stdout, "\tu_%d%d = %.5e\t+/-\t%.5e\t",
				i + 1, i + 1,
				fabs ( gsl_vector_get ( w->x, p + i + 2 ) ),
				correction * sqrt ( gsl_matrix_get
					( covar, p + i + 2, p + i + 2 ) ) );
			if ( i != -1 ) {
				fprintf ( stdout, "(%s)\n",
					data->description.components[i] );
			} else {
				fprintf ( stdout, "(water)\n" );
			}
			/* add solute name */
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

void save_uniquac ( System *data, info *user_data,
		gsl_multifit_nlinear_workspace *w ) {

	int i, j, k, lines, comps;
	long double phi_calc, phi_real;
	long double ln_gamma_w;
	long double x_w, x_j, x_k;
	long double sumxjlj, sumqjxj, sumthetajtaujw, sumthetaktaukj, sumsum;
	long double l_w, l_j, r_w, r_j, q_w, q_j, q_k;
	long double theta_w, theta_j, theta_k, Phi_w;
	long double u_jj, u_jw, u_wj, u_kj, u_kk, u_ww, tau_jw, tau_kj, tau_wj;
	char *filename;
	FILE *results_file;
	gsl_vector *x = gsl_multifit_nlinear_position (w);

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
	fprintf ( results_file, "%s\n", data->description.components[i] );

	/*
	* Shamelessly copied from the function at the top of this file.
	* I'm sorry.
	*/
	for ( i = 0; i < lines ; i++ ) {
		x_w = 1;
		for ( j = 0; j < comps; j++ ) {
			x_w -= data->x_and_aw.x[i][j];
		}      /* In this model the water is seen as another component;
			* therefore, we also need its molar fraction.
			*/

		r_w = fabs ( gsl_vector_get ( x, 0 ) );
		l_w = 1 - r_w;
			/*
			* l = (z/2) * (r-q) - (r-1)
			* :. l = -(r-1) = 1 - r
			*/
		sumxjlj = x_w * l_w;

		q_w = r_w;
		sumqjxj = q_w * x_w;

		for ( j = 0; j < comps; j++ ) {

			x_j = data->x_and_aw.x[i][j];
			r_j = fabs ( gsl_vector_get ( x, j + 1 ) );
			l_j = 1 - r_j;
			sumxjlj += x_j * l_j;
			/* here we dealt with the l's and their sum */

			q_j = r_j;
			sumqjxj += q_j * x_j;
			/* needed for theta_j */
		}

		theta_w = ( q_w * x_w ) / sumqjxj;
		Phi_w = theta_w; /* because q_i = r_i */

		sumthetajtaujw = theta_w; /* = theta_w * tau_ww = theta_w * 1 */
		for ( j = 0; j < comps; j++ ) {
			q_j = fabs ( gsl_vector_get ( x, j + 1 ) );
			x_j = data->x_and_aw.x[i][j];
			theta_j = ( q_j * x_j ) / sumqjxj;
			u_ww = fabs ( gsl_vector_get ( x, comps + 1 ) );
			u_jj = fabs ( gsl_vector_get ( x, comps + j + 2 ) );
			u_jw = sqrt ( u_jj * u_ww );
			tau_jw = exp ( - ( u_jw - u_ww ) / ( R * TEMP ) );
			/*fprintf ( stderr, "tau_jw = %Lf\n", tau_jw );*/
			sumthetajtaujw += theta_j * tau_jw;
		}

		sumsum = 0;
		for ( j = -1; j < comps; j++ ) {

			sumthetaktaukj = theta_w;
			for ( k = -1; k < comps; k++ ) {
				if ( k == -1 ) {
					x_k = x_w;
				} else {
					x_k = data->x_and_aw.x[i][k];
				}
				q_k = fabs ( gsl_vector_get ( x, k + 1 ) );
				theta_k = ( q_k * x_k ) / sumqjxj;
				u_kk = fabs ( gsl_vector_get ( x, comps + k + 2 ) );
				u_jj = fabs ( gsl_vector_get ( x, comps + j + 2 ) );
				u_kj = sqrt ( u_jj * u_kk );
				tau_kj = exp ( - ( u_kj - u_jj ) / ( R * TEMP ) );
				sumthetaktaukj += theta_k * tau_kj;
			}

			q_j = fabs ( gsl_vector_get ( x, j + 1 ) );
			if ( j == -1 ) {
				x_j = x_w;
			} else {
				x_j = data->x_and_aw.x[i][j];
			}
			theta_j = ( q_j * x_j ) / sumqjxj;
			u_ww = fabs ( gsl_vector_get ( x, comps ) );
			u_jj = fabs ( gsl_vector_get ( x, comps + j + 2 ) );
			u_wj = sqrt ( u_ww * u_jj );
			tau_wj = exp ( - ( u_wj - u_jj ) / ( R * TEMP ) );

			sumsum += ( theta_j * tau_wj ) / sumthetaktaukj;

		}

		ln_gamma_w = log ( Phi_w / x_w );
			/* "ln_gamma_w +=
			* ( Z_COORD / 2 ) * q_w * log ( theta_w / Phi_w );"
			* is a expression that would be here if we did not
			* assume r_j = q_j; however, because we assumed this,
			* the expression between parenthesis inside the log
			* equals  one, and this forces the whole right-hand
			* side to 0.
			*/
		ln_gamma_w += l_w;
		ln_gamma_w -= ( Phi_w / x_w ) * sumxjlj;
		ln_gamma_w -= q_w * ( 1 - log ( sumthetajtaujw ) - sumsum );

		phi_calc = ( ln_gamma_w + log (x_w) ) / log (x_w);
		phi_real = log (data->x_and_aw.aw[i]) / log (x_w);
		fprintf ( results_file, "%Lf,%Lf,", phi_calc, phi_real );
		if ( user_data->aw_in_results == TRUE ) {
			fprintf ( results_file, "%Lf,%f,",
				exp (ln_gamma_w) * x_w, data->x_and_aw.aw[i] );
		}
		for ( j = 0; j < comps - 1; j++ ) {
			fprintf ( results_file, "%f,", data->x_and_aw.x[i][j] );
		}
		fprintf ( results_file, "%f\n", data->x_and_aw.x[i][comps-1] );
	}
	fclose (results_file);
}
