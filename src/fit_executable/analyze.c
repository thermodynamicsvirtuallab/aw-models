#include "definitions_and_headers.h"

void analyze_all_models ( System *data, info *user_data ) {

	free (user_data->model);
	user_data->model = malloc ( strlen ("norrish") + 1 );
	strcpy ( user_data->model, "norrish" );
	fit_to_model ( data, user_data );
	if ( fpclassify ( user_data->cost ) == FP_ZERO ||
			fpclassify ( user_data->cost ) == FP_NORMAL ) {
		fprintf ( stdout, "Final cost for Norrish's model:   \t\t%f\n",
			user_data->cost );
	} else {
		fprintf ( stderr,
			"\t\tFinal cost for Norrish's model not available.\n" );
	}

	free (user_data->model);
	user_data->model = malloc ( strlen ("virial") + 1 );
	strcpy ( user_data->model, "virial" );
	fit_to_model ( data, user_data );
	if ( fpclassify ( user_data->cost ) == FP_ZERO ||
			fpclassify ( user_data->cost ) == FP_NORMAL ) {
		fprintf ( stdout, "Final cost for Virial model:      \t\t%f\n",
			user_data->cost );
	} else {
		fprintf ( stderr,
			"\t\tFinal cost for Virial model not available.\n" );
	}

	free (user_data->model);
	user_data->model = malloc ( strlen ("uniquac") + 1 );
	strcpy ( user_data->model, "uniquac" );
	fit_to_model ( data, user_data );
	if ( fpclassify ( user_data->cost ) == FP_ZERO ||
			fpclassify ( user_data->cost ) == FP_NORMAL ) {
		fprintf ( stdout, "Final cost for UNIQUAC model:      \t\t%f\n",
			user_data->cost );
	} else {
		fprintf ( stderr,
			"\t\tFinal cost for UNIQUAC model not available.\n" );
	}

	free (user_data->model);
	user_data->model = malloc ( strlen ("raoult") + 1 );
	strcpy ( user_data->model, "raoult" );
	check_model ( data, user_data );
	if ( fpclassify ( user_data->cost ) == FP_ZERO ||
			fpclassify ( user_data->cost ) == FP_NORMAL ) {
		fprintf ( stdout, "Final cost for Raoult's model:    \t\t%f\n",
			user_data->cost );
	} else {
		fprintf ( stderr,
			"\t\tFinal cost for Raoult's model not available.\n" );
	}

	free (user_data->model);
	user_data->model = malloc ( strlen ("caurie") + 1 );
	strcpy ( user_data->model, "caurie" );
	check_model ( data, user_data );
	if ( fpclassify ( user_data->cost ) == FP_ZERO ||
			fpclassify ( user_data->cost ) == FP_NORMAL ) {
		fprintf ( stdout, "Final cost for Caurie's model:    \t\t%f\n",
			user_data->cost );
	} else {
		fprintf ( stderr,
			"\t\tFinal cost for Caurie's model not available.\n" );
	}

	if ( user_data->gave_filenames == TRUE &&
			user_data->not_zdan != TRUE ) {
		free (user_data->model);
		user_data->model = malloc ( strlen ("zdanovskii") + 1 );
		strcpy ( user_data->model, "zdanovskii" );
		check_model ( data, user_data );
		if ( fpclassify ( user_data->cost ) == FP_ZERO ||
				fpclassify ( user_data->cost ) == FP_NORMAL ) {
			fprintf ( stdout, "Final cost for Zdanovskii's model:\t\t%f\n",
				user_data->cost );
		} else {
			fprintf ( stderr,
			"\t\tFinal cost for Zdanovskii's model not available.\n" );
		}
	} else if ( user_data->not_zdan != TRUE ) {
		fprintf ( stderr, "\tWARNING: Not possible to compute" );
		fprintf ( stderr, " Zdanovskii's results \n" );
	}

}
