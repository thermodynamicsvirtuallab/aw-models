#include "definitions_and_headers.h"

/*
 * This function performs an analysis of non-parametric models, as
 * Raoult's, Caurie's or Zdanovskii's.
 */

int check_model ( System *data, info *user_data ) {

	check_function check;
	print_function_for_check print;
	save_function_for_check save;
	int n;
	double *errors;

	n = data->description.dataset_size;
	errors = malloc ( n * sizeof (double) );
	/* array in which the absolute values
	* of the errors shall be stored
	*/

	if ( strcmp ( user_data->model, "raoult" ) == TRUE ) {
		check = &check_raoult;
		print = &print_raoult;
		save = &save_raoult;
	} else if ( strcmp ( user_data->model, "caurie" ) == TRUE ) {
		check = &check_caurie;
		print = &print_caurie;
		save = &save_caurie;
	} else if ( strcmp ( user_data->model, "zdanovskii" ) == TRUE ) {
		check = &check_zdanovskii;
		print = &print_zdanovskii;
		save = &save_zdanovskii;
	} else {
		fprintf ( stderr, "Model unknown. Assuming Raoult's Law.\n" );
		check = &check_raoult;
		print = &print_raoult;
		save = &save_raoult;
	}

	check ( data, user_data, errors );
	print ( data, user_data, errors );

	if ( user_data->save_new_results == TRUE ) {
		save ( data, user_data );
	}

	free (errors);

	return 0;
}

