#include "../definitions_and_headers.h"

/*
 * This function compares the values of phi predicted with Raoult's model
 * with the values of phi experimentally obtained.
 *     phi, here, is constant and equals 1; phi = ln(aw_raoult)/ln(xw) =
 * ln(xw)/ln(xw) = 1.
 */
void check_raoult ( System *data, info *user_data, double *errors ) {

	int i, j, n, p;
	double xw;

	n = data->description.dataset_size;
	p = data->description.n_of_comps;

	for ( i = 0; i < n; i++ ) {
		xw = 1;
		for ( j = 0; j < p; j++ ) {
			xw -= data->x_and_aw.x[i][j];
		}
		if ( data->description.has_aw_data == TRUE ) {
			errors[i] = ( 1 - ( log(data->x_and_aw.aw[i]) / log(xw) ) );
		} else {
			/* For when isopiestic data (not aw, directly)
			* is informed (check ./virial.c).
			*/
			errors[i] = ( 1 - log ( 1 - data->x_and_aw.aw[i] ) /
					log (xw) );
		}
	}
}

void print_raoult ( System *data, info *user_data, double *errors ) {

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

void save_raoult ( System *data, info *user_data ) {

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
		phi_calc = 1.0;
		fprintf ( results_file, "%f,%f,", phi_calc, phi_exp );
		if ( user_data->aw_in_results == TRUE ) {
			fprintf ( results_file, "%f,%f,",
				xw, data->x_and_aw.aw[i] );
		}
		for ( j = 0; j < comps - 1; j++ ) {
			fprintf ( results_file, "%f,", data->x_and_aw.x[i][j] );
		}
		fprintf ( results_file, "%f\n", data->x_and_aw.x[i][comps-1] );
	}

	fclose (results_file);

}
