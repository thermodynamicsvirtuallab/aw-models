#include "../definitions_and_headers.h"

/*
 * This function compares the values of phi predicted with Raoult's model
 * with the values of phi experimentally obtained.
 */
void check_caurie ( System *data, info *user_data, double *errors ) {

	int i, j, k, l, o, n, p;
	double xw, prod_aw, *m, phi_real, phi_calc;

	n = data->description.dataset_size;
	p = data->description.n_of_comps;

	data->x_and_aw.aw_calc = malloc ( n * sizeof(double) );

	m = malloc ( p * sizeof (double) );
	for ( i = 0; i < p; i++ ) {
		m[i] = 0;
	}

	for ( i = 0; i < n; i++ ) {
		xw = 1;
		for ( j = 0; j < p; j++ ) {
			xw -= data->x_and_aw.x[i][j];
		}
		prod_aw = 1;
		for ( j = 0; j < p; j++ ) {
			prod_aw *= xw / ( xw + data->x_and_aw.x[i][j] );
		}
		/*
		* obtaining the water activity of a binary solution
		* in which the solute maintains constant its molality
		*/
		if ( p > 1 ) {
			for ( j = 0; j < p; j++ ) {
				m[j] = data->x_and_aw.x[i][j] /
					( xw * KGS_IN_MOL_WATER );
			}
			for ( k = 0; k < p; k++ ) {
				for ( l = 0; l < k; l++ ) {
					prod_aw -= ( p / FAC_B_CAURIE ) *
						m[k] * m[l];
				}
			}
			/*
			* Caurie's correction for ternary solutions
			*/
			if ( p > 2 ) {
				for ( k = 0; k < p; k++ ) {
					for ( l = 0; l < k; l++ ) {
						for ( o = 0; o < l; o++ ) {
							prod_aw -=
								( ( p + 1 ) /
								FAC_C_CAURIE ) *
								m[k] * m[l] * m[o];
						}
					}
				}
				/*
				* Caurie's correction for n-ary solutions
				*/
			}
		}
		data->x_and_aw.aw_calc[i] = prod_aw;
		phi_calc = log (prod_aw) / log (xw);
		if ( data->description.has_aw_data == TRUE ) {
			phi_real = log (data->x_and_aw.aw[i]) / log (xw);
		} else {
			/* This is for non-aw (isopiestic) data.
			* Check ./virial.c for detailed explanations.
			*/
			phi_real = ( 1 - log (data->x_and_aw.aw[i]) ) /
				log (xw);
		}
		errors[i] = phi_calc - phi_real;
	}

	free (m);
}

void print_caurie ( System *data, info *user_data, double *errors ) {

	int i, n;
	double cost, R_squared, R_squared_aw;

	n = data->description.dataset_size;

	cost = 0;
	for ( i = 0; i < n; i++ ) {
		cost += pow ( errors[i], 2 );
	}

	user_data->cost = sqrt ( cost / n );

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

void save_caurie ( System *data, info *user_data ) {

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

