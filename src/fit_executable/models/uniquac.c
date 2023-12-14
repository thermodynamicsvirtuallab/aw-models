#include "../definitions_and_headers.h"
#include "../unifac_header.h"

double uni_ind ( int i, int j, int n, const gsl_vector *K );

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
	long double sumrjxj;
	long double l_w, l_j, r_w, r_j, q_w, q_j, q_k;
	long double theta_w, theta_j, theta_k, Phi_w;
	long double tau_jw, tau_kj, tau_wj;
	long double A_wj, A_jw, A_kj;
	int n, p;
	int i, j, k;

	data = (System *) params;

	n = data->description.dataset_size;
	p = data->description.n_of_comps;

	ln_gamma_w = 0.0;
	x_w = 0.0;
	for ( i = 0; i < n; i++ ) {
		x_w = 1;
		for ( j = 0; j < p; j++ ) {
			x_w -= data->x_and_aw.x[i][j];
		}      /* In this model the water is seen as another component;
			* therefore, we also need its molar fraction.
			*/

		r_w = R_WATER;
		q_w = Q_WATER;
		l_w = (Z_COORD / 2) * (r_w - q_w) - (r_w - 1);
		sumxjlj = x_w * l_w;

		sumqjxj = q_w * x_w;
		sumrjxj = r_w * x_w;

		for ( j = 0; j < p; j++ ) {

			x_j = data->x_and_aw.x[i][j];
			r_j = data->description.r_vals[j];
			q_j = data->description.q_vals[j];
			l_j = (Z_COORD / 2) * (r_j - q_j) - (r_j - 1);
			sumxjlj += x_j * l_j;
			/* here we dealt with the l's and their sum */

			sumqjxj += q_j * x_j;
			/* needed for theta_j */
			sumrjxj += r_j * x_j;
			/* needed for Phi_j */
		}

		theta_w = ( q_w * x_w ) / sumqjxj;
		Phi_w = (r_w * x_w) / sumrjxj;

		sumthetajtaujw = theta_w; /* = theta_w * tau_ww = theta_w * 1 */
		for ( j = 0; j < p; j++ ) {
			r_j = data->description.r_vals[j];
			q_j = data->description.q_vals[j];
			x_j = data->x_and_aw.x[i][j];
			theta_j = ( q_j * x_j ) / sumqjxj;
			A_jw = uni_ind( j+1, 0, p+1, K );
			tau_jw = exp ( - A_jw / ( R * data->description.temp ) );
			sumthetajtaujw += theta_j * tau_jw;
		}

		sumsum = 0;
		for ( j = -1; j < p; j++ ) {

			sumthetaktaukj = 0;
			for ( k = -1; k < p; k++ ) {
				if ( k == -1 ) {
					x_k = x_w;
					q_k = Q_WATER;
				} else {
					x_k = data->x_and_aw.x[i][k];
					q_k = data->description.q_vals[k];
				}
				theta_k = ( q_k * x_k ) / sumqjxj;
				A_kj = uni_ind ( k + 1, j + 1, p + 1, K );
				tau_kj = exp ( - A_kj /
						( R * data->description.temp ) );
				sumthetaktaukj += theta_k * tau_kj;
			}

			if ( j == -1 ) {
				x_j = x_w;
				r_j = R_WATER;
				q_j = Q_WATER;
			} else {
				x_j = data->x_and_aw.x[i][j];
				r_j = data->description.r_vals[j];
				q_j = data->description.q_vals[j];
			}
			theta_j = ( q_j * x_j ) / sumqjxj;
			A_wj = uni_ind ( 0, j + 1, p + 1, K );
			tau_wj = exp ( - A_wj / ( R * data->description.temp ) );

			sumsum += ( theta_j * tau_wj ) / sumthetaktaukj;

		}

		ln_gamma_w = log ( Phi_w / x_w );
		ln_gamma_w += ( Z_COORD / 2 ) * q_w * log ( theta_w / Phi_w );
		ln_gamma_w += l_w;
		ln_gamma_w -= ( Phi_w / x_w ) * sumxjlj;
		ln_gamma_w += q_w * ( 1 - log ( sumthetajtaujw ) - sumsum );
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

			r_w = R_WATER;
			q_w = Q_WATER;
			l_w = (Z_COORD / 2) * (r_w - q_w) - (r_w - 1);

			r_j = data->description.r_vals[0];
			q_j = data->description.q_vals[0];
			l_j = (Z_COORD / 2) * (r_j - q_j) - (r_j - 1);

			sumxjlj = ( x_w * l_w ) + ( x_j * l_j );
			sumqjxj = ( q_w * x_w ) + ( q_j * x_j );
			theta_w = ( q_w * x_w ) / sumqjxj;
			Phi_w = theta_w;

			sumthetajtaujw = theta_w;
			theta_j = ( q_j * x_j ) / sumqjxj;
			A_jw = uni_ind ( 1, 0, p + 1, K );
			tau_jw = exp ( - A_jw / ( R * data->description.temp ) );
			sumthetajtaujw += theta_j * tau_jw;

			sumsum = 0;
			for ( j = -1; j < 1; j++ ) {

				sumthetaktaukj = 0;
				for ( k = -1; k < 1; k++ ) {
					if ( k == -1 ) {
						x_k = x_w;
						q_k = Q_WATER;
					} else {
						x_k = data->x_and_aw.aw[i];
						q_k = data->description.q_vals[0];
					}
					theta_k = ( q_k * x_k ) / sumqjxj;
					A_kj = uni_ind ( k+1, j+1, p+1, K );
					tau_kj = exp ( - A_kj /
						( R * data->description.temp ) );
					sumthetaktaukj += theta_k * tau_kj;
				}

				if ( j == -1 ) {
					x_j = x_w;
					r_j = R_WATER;
					q_j = Q_WATER;
				} else {
					x_j = data->x_and_aw.x[i][j];
					r_j = data->description.r_vals[0];
					q_j = data->description.q_vals[0];
				}
				theta_j = ( q_j * x_j ) / sumqjxj;
				A_wj = uni_ind ( 0, j+1, p+1, K );
				tau_wj = exp ( - A_wj /
						( R * data->description.temp ) );

				sumsum += ( theta_j * tau_wj ) / sumthetaktaukj;

			}


			ln_gamma_w = log ( Phi_w / x_w );
			ln_gamma_w += ( Z_COORD / 2 ) * q_w *
						log ( theta_w / Phi_w );
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
				gsl_vector_get ( x, i ) );
	}
	fprintf ( stderr, "\n" );

}

void print_uniquac ( gsl_matrix *covar, gsl_multifit_nlinear_workspace *w,
		int status, double chisq0, double chisq,
		System *data, info *user_data ) {

	int p, n;
	int i, j, counter;
	double correction, R_squared, R_squared_aw;

	p = data->description.n_of_comps;
	n = data->description.dataset_size;

	correction = GSL_MAX_DBL ( 1, sqrt ( chisq / (n - p) ) );

	if ( user_data->is_all == FALSE ) {
		fprintf ( stdout, "Results Obtained:\n" );

		for ( i = -1; i < p; i++ ) {
			if ( i != -1 ) {
				fprintf ( stdout, "\tr_%d  = %.5e\t",
					i + 1, data->description.r_vals[i]);
				fprintf ( stdout, "(%s)\n",
					data->description.components[i] );
				fprintf ( stdout, "\tq_%d  = %.5e\t",
					i + 1, data->description.q_vals[i]);
				fprintf ( stdout, "(%s)\n",
					data->description.components[i] );
			} else {
				fprintf ( stdout, "\tr_%d  = %.5e\t",
					i + 1, R_WATER);
				fprintf ( stdout, "(water)\n" );
				fprintf ( stdout, "\tq_%d  = %.5e\t",
					i + 1, Q_WATER);
				fprintf ( stdout, "(water)\n" );
			}
			/* add solute name */
		}
		counter = 0;
		for ( i = -1; i < p; i++ ) {
			for ( j = -1; j < p; j++ ) {
				if ( i != j ) {
					fprintf ( stdout,
						"\tA_%d%d = %.5e\t+/-\t%.5e\t",
						i + 1, j + 1,
						uni_ind ( i + 1, j + 1,
							p + 1, w->x ),
					correction * sqrt ( gsl_matrix_get
						( covar, counter, counter ) ) );
					counter++;
				} else {
					fprintf ( stdout,
						"\tA_%d%d = 0.0\t\t+/-",
						i + 1, j + 1);
					fprintf ( stdout, "\t\t0.0\t");
				}
				if ( i != -1 ) {
					fprintf ( stdout, "(%s/",
						data->description.components[i]);
				} else {
					fprintf ( stdout, "(water/" );
				}
				if ( j != -1 ) {
					fprintf ( stdout, "%s)\n",
						data->description.components[j]);
				} else {
					fprintf ( stdout, "water)\n" );
				}
			/* add solute name */
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

void save_uniquac ( System *data, info *user_data,
		gsl_multifit_nlinear_workspace *w ) {

	int i, j, k, n, p;
	long double phi_i_calc, phi_i_real;
	long double ln_gamma_w;
	long double x_w, x_j, x_k;
	long double sumxjlj, sumrjxj, sumqjxj, sumthetajtaujw;
	long double sumthetaktaukj, sumsum;
	long double l_w, l_j, r_w, r_j, q_w, q_j, q_k;
	long double theta_w, theta_j, theta_k, Phi_w;
	long double tau_jw, tau_kj, tau_wj;
	long double A_jw, A_wj, A_kj;
	char *filename;
	FILE *results_file;

	filename = user_data->filename_new_results;
	results_file = fopen ( filename, "w" );

	n = data->description.dataset_size;
	p = data->description.n_of_comps;

	fprintf ( results_file, "phi_calc,phi_exp," );
	if ( user_data->aw_in_results == TRUE ) {
		fprintf ( results_file, "aw_calc,aw_exp,xw," );
	}

	for ( i = 0; i < p - 1; i++ ) {
		fprintf ( results_file, "%s,", data->description.components[i] );
	}
	fprintf ( results_file, "%s\n", data->description.components[i] );

	/*
	* Shamelessly copied from the function at the top of this file.
	* I'm sorry.
	*/
	for ( i = 0; i < n; i++ ) {
		x_w = 1;
		for ( j = 0; j < p; j++ ) {
			x_w -= data->x_and_aw.x[i][j];
		}      /* In this model the water is seen as another component;
			* therefore, we also need its molar fraction.
			*/

		r_w = R_WATER;
		q_w = Q_WATER;
		l_w = (Z_COORD / 2) * (r_w - q_w) - (r_w - 1);
		sumxjlj = x_w * l_w;

		sumqjxj = q_w * x_w;
		sumrjxj = r_w * x_w;

		for ( j = 0; j < p; j++ ) {

			x_j = data->x_and_aw.x[i][j];
			r_j = data->description.r_vals[j];
			q_j = data->description.q_vals[j];
			l_j = (Z_COORD / 2) * (r_j - q_j) - (r_j - 1);
			sumxjlj += x_j * l_j;
			/* here we dealt with the l's and their sum */

			sumqjxj += q_j * x_j;
			/* needed for theta_j */
			sumrjxj += r_j * x_j;
			/* needed for Phi_j */
		}

		theta_w = ( q_w * x_w ) / sumqjxj;
		Phi_w = (r_w * x_w) / sumrjxj;

		sumthetajtaujw = theta_w; /* = theta_w * tau_ww = theta_w * 1 */
		for ( j = 0; j < p; j++ ) {
			r_j = data->description.r_vals[j];
			q_j = data->description.q_vals[j];
			x_j = data->x_and_aw.x[i][j];
			theta_j = ( q_j * x_j ) / sumqjxj;
			A_jw = uni_ind( j + 1, 0, p + 1, w->x );
			tau_jw = exp ( - A_jw / ( R * data->description.temp ) );
			sumthetajtaujw += theta_j * tau_jw;
		}

		sumsum = 0;
		for ( j = -1; j < p; j++ ) {

			sumthetaktaukj = 0;
			for ( k = -1; k < p; k++ ) {
				if ( k == -1 ) {
					x_k = x_w;
					q_k = Q_WATER;
				} else {
					x_k = data->x_and_aw.x[i][k];
					q_k = data->description.q_vals[k];
				}
				theta_k = ( q_k * x_k ) / sumqjxj;
				A_kj = uni_ind ( k + 1, j + 1, p + 1, w->x );
				tau_kj = exp ( - A_kj /
						( R * data->description.temp ) );
				sumthetaktaukj += theta_k * tau_kj;
			}

			if ( j == -1 ) {
				x_j = x_w;
				r_j = R_WATER;
				q_j = Q_WATER;
			} else {
				x_j = data->x_and_aw.x[i][j];
				r_j = data->description.r_vals[j];
				q_j = data->description.q_vals[j];
			}
			theta_j = ( q_j * x_j ) / sumqjxj;
			A_wj = uni_ind ( 0, j + 1, p + 1, w->x );
			tau_wj = exp ( - A_wj /
					( R * data->description.temp ) );

			sumsum += ( theta_j * tau_wj ) / sumthetaktaukj;

		}

		ln_gamma_w = log ( Phi_w / x_w );
		ln_gamma_w += ( Z_COORD / 2 ) * q_w * log ( theta_w / Phi_w );
		ln_gamma_w += l_w;
		ln_gamma_w -= ( Phi_w / x_w ) * sumxjlj;
		ln_gamma_w += q_w * ( 1 - log ( sumthetajtaujw ) - sumsum );
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

			r_w = R_WATER;
			q_w = Q_WATER;
			l_w = (Z_COORD / 2) * (r_w - q_w) - (r_w - 1);

			r_j = data->description.r_vals[0];
			q_j = data->description.q_vals[0];
			l_j = (Z_COORD / 2) * (r_j - q_j) - (r_j - 1);

			sumxjlj = ( x_w * l_w ) + ( x_j * l_j );
			sumqjxj = ( q_w * x_w ) + ( q_j * x_j );
			theta_w = ( q_w * x_w ) / sumqjxj;
			Phi_w = theta_w;

			sumthetajtaujw = theta_w;
			theta_j = ( q_j * x_j ) / sumqjxj;
			A_jw = uni_ind ( 1, 0, p + 1, w->x );
			tau_jw = exp ( - A_jw /
					( R * data->description.temp ) );
			sumthetajtaujw += theta_j * tau_jw;

			sumsum = 0;
			for ( j = -1; j < 1; j++ ) {

				sumthetaktaukj = theta_w;
				for ( k = -1; k < 1; k++ ) {
					if ( k == -1 ) {
						x_k = x_w;
						q_k = Q_WATER;
					} else {
						x_k = data->x_and_aw.aw[i];
						q_k = data->description.q_vals[0];
					}
					theta_k = ( q_k * x_k ) / sumqjxj;
					A_kj = uni_ind ( k+1, j+1, p+1, w->x );
					tau_kj = exp ( - A_kj /
						( R * data->description.temp ) );
					sumthetaktaukj += theta_k * tau_kj;
				}

				if ( j == -1 ) {
					x_j = x_w;
					r_j = R_WATER;
					q_j = Q_WATER;
				} else {
					x_j = data->x_and_aw.x[i][j];
					r_j = data->description.r_vals[0];
					q_j = data->description.q_vals[0];
				}
				theta_j = ( q_j * x_j ) / sumqjxj;
				A_wj = uni_ind ( 0, j+1, p+1, w->x );
				tau_wj = exp ( - A_wj /
						( R * data->description.temp ) );

				sumsum += ( theta_j * tau_wj ) / sumthetaktaukj;

			}


			ln_gamma_w = log ( Phi_w / x_w );
			ln_gamma_w += ( Z_COORD / 2 ) * q_w *
						log ( theta_w / Phi_w );
			ln_gamma_w += l_w;
			ln_gamma_w -= ( Phi_w / x_w ) * sumxjlj;
			ln_gamma_w -= q_w * ( 1 - log ( sumthetajtaujw ) - sumsum );

			phi_i_real = ( ln_gamma_w +
				log ( 1 - data->x_and_aw.aw[i] ) ) / log (x_w);
		}

		fprintf ( results_file, "%Lf,%Lf,", phi_i_calc, phi_i_real );
		if ( user_data->aw_in_results == TRUE ) {
			fprintf ( results_file, "%Lf,%f,%Lf,",
				exp (ln_gamma_w) * x_w,
				data->x_and_aw.aw[i], x_w);
		}
		for ( j = 0; j < p - 1; j++ ) {
			fprintf ( results_file, "%f,", data->x_and_aw.x[i][j] );
		}
		fprintf ( results_file, "%f\n", data->x_and_aw.x[i][p-1] );
	}
	fclose (results_file);
}

/*
 * This function returns the value of A_ij in the gsl_vector K, when i != j;
 * otherwise, it returns the defined value of A_ij, zero.
 */
double uni_ind ( int i, int j, int n, const gsl_vector *K ) {

	double retval;

	retval = 0;
	if ( i > j ) {
		retval = gsl_vector_get ( K, i - 1 + ( n - 1 ) * j );
	} else if ( i < j ) {
		retval = gsl_vector_get ( K, i+ ( n - 1 ) * j );
	}

	return retval;
}
