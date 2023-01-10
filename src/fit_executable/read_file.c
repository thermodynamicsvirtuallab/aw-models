# include "definitions_and_headers.h"

/*
 * This function reads data from a file and stores it in a struct of type Data.
 * It also stores the number of data points, the number of components and their
 * names in another struct, of type Metadata. Its arguments are the name of the
 * .csv file in which the data are stored, and pointers to the structs.
 */

void initialize ( char *filename, Metadata *system_description,
		Data *system, info *user_data ) {

	int lines, columns, line_is_empty, i, j;
	int lines_zdan[2];
	char tmpchar, str[256];
	char **composition;
	double **x_values;
	double *y_values;
	double avg_aw, avg_xw;
	FILE *file;

	file = fopen ( filename, "r" );
	lines = -1;
	columns = 0;
	line_is_empty = TRUE;
	do {
		tmpchar = fgetc ( file );
		if ( (tmpchar == '\n') && (line_is_empty == FALSE) ) {
			lines ++;
			line_is_empty = TRUE;
		} else if ( tmpchar == ',' ) {
			columns ++;
			line_is_empty = FALSE;
		}
	} while ( tmpchar != EOF );
	fclose ( file );
	columns = columns / lines;
		/*
		* Here we get the number of components and
		* data points in the file (lines and columns)
		*/

	composition = malloc ( columns * sizeof(char *) );

	x_values = malloc ( lines * sizeof(double *) );
	y_values = malloc ( lines * sizeof(double) );
	for ( i = 0; i < lines; i++ ) {
		x_values[i] = malloc ( columns * sizeof (double) );
	}

	if ( strcmp ( user_data->model, "zdanovskii" ) == TRUE ||
			( strcmp ( user_data->model, "all" ) == TRUE &&
			user_data->gave_filenames == TRUE ) ) {
		system->x_zdan = malloc ( 2 * sizeof (double *) );
		system->aw_zdan = malloc ( 2 * sizeof (double *) );
		system->n_zdan = malloc ( 2 * sizeof (int) );
		for ( i = 0; i < 2; i++ ) {

			/* read number of lines in file to allocate memory */
			file = fopen ( user_data->files_zdan[i], "r" );
			lines_zdan[i] = -1;
			line_is_empty = TRUE;
			do {
				tmpchar = fgetc (file);
				if ( tmpchar == '\n' && line_is_empty == FALSE ) {
					lines_zdan[i]++;
					line_is_empty = TRUE;
				} else if ( tmpchar == ',' ) {
					line_is_empty = FALSE;
				}
			} while ( tmpchar != EOF );
			fclose (file);
			system->x_zdan[i] = malloc ( lines_zdan[i] *
					sizeof (double) );
			system->aw_zdan[i] = malloc ( lines_zdan[i] *
					sizeof (double) );
			system->n_zdan[i] = lines_zdan[i];


			/* read file data to data structures */
			file = fopen ( user_data->files_zdan[i], "r" );
			fscanf ( file, "%127[^,\n]", str );
			fscanf ( file, "%*c" );
			fscanf ( file, "%127[^,\n]", str);
			fscanf ( file, "%*c"); /* get rid of names */
			for ( j = 0; j < lines_zdan[i]; j++ ) {
				fscanf ( file, "%127[^,\n]", str );
				fscanf ( file, "%*c" );
				system->aw_zdan[i][j] = strtod ( str, NULL );
				fscanf ( file, "%127[^,\n]", str );
				fscanf ( file, "%*c" );
				system->x_zdan[i][j] = strtod ( str, NULL );
			}
			fclose (file);
		}

	} else {
		system->n_zdan = NULL;
		system->x_zdan = NULL;
		system->aw_zdan = NULL;
		user_data->files_zdan = NULL;
	}

	system->aw_calc = NULL;


		/*
		* And now our structures are initialized
		*/

	file = fopen ( filename, "r" );

	fscanf(file, "%127[^,\n]", str);
	fscanf(file, "%*c");

	for ( i = 0; i < columns; i ++ ) {
		fscanf(file, "%127[^,\n]", str);
		fscanf(file, "%*c");
		composition[i] = malloc ( (strlen(str) + 1) * sizeof (char) );
		strcpy ( composition[i], str );
	}

	for ( i = 0; i < lines; i++ ) {
		fscanf(file, "%127[^,\n]", str);
		fscanf(file, "%*c");
		y_values[i] = strtod ( str, NULL );
		for ( j = 0; j < columns; j++ ) {
			fscanf(file, "%127[^,\n]", str);
			fscanf(file, "%*c");
			x_values[i][j] = strtod ( str, NULL );
		}
	}

	fclose ( file );

	avg_aw = 0;
	avg_xw = 0;
	for ( i = 0; i < lines; i++ ) {
		avg_aw += y_values[i] / lines;
		avg_xw += 1.0 / lines;
		for ( j = 0; j < columns; j++ ) {
			avg_xw -= x_values[i][j] / lines;
		}
	}
		/*
		* And now we have read the data from the file
		*/

	system_description->dataset_size = lines;
	system_description->n_of_comps = columns;
	system_description->has_aw_data = user_data->has_aw_data;
	system_description->components = composition;

	system->x = x_values;
	system->aw = y_values;

	fprintf ( stdout, "%d lines read from file \"%s\" read.\n",
			lines, filename );
	fprintf ( stdout, "Average activity: avg(aw) = %f\n", avg_aw );
	fprintf ( stdout, "Average dilution: avg(xw) = %f\n", avg_xw );

		/*
		* And now our data and metadata are stored in our structs
		*/

}

/*
 * Initialize solver for different models.
 */
void init_data ( char *model, double *x_init, int p, info *user_data ) {

	int i, n;

	if ( user_data->K_number != p ) {
		if ( user_data->K_number != 0 ) {
			fprintf ( stderr,
				"Testing impossible; not enough parameters.\n" );
			fprintf ( stderr,
				"Fitting as normally\n" );
		}
		if ( strcmp ( model, "uniquac" ) == TRUE ) {
			n = ( p / 2 );
			for ( i = 0; i < n; i++ ) {
				x_init[i] = 10;
				x_init[n+i] = 1000;
			}
		} else {
			for ( i = 0; i < p; i++ ) {
				x_init[i] = 1.0;
			}
		}
	} else {
		for ( i = 0; i < p; i++ ) {
			x_init[i] = user_data->K[i];
		}
	}
}

/*
 * This function cleans up the memory utilized,
 * freeing pointers to data and metadata.
 */

void finalize ( Metadata *system_description, Data *system,
		info *user_data ) {

	int i, lines = system_description->dataset_size;
	int components = system_description->n_of_comps;

	for ( i = 0; i < lines; i++ ) {
		free (system->x[i]);
	}

	if ( strcmp ( user_data->model, "zdanovskii" ) == TRUE ) {
		free (system->x_zdan[0]);
		free (system->x_zdan[1]);
		free (system->x_zdan);
		free (system->aw_zdan[0]);
		free (system->aw_zdan[1]);
		free (system->aw_zdan);
		free (system->n_zdan);
	}

	if ( system->aw_calc != NULL ) {
		free (system->aw_calc);
	}

	if ( user_data->files_zdan != NULL ) {
		for ( i = 0; i < 2; i++ ) {
			if ( user_data->files_zdan[i] != NULL ) {
				free (user_data->files_zdan[i]);
			}
		}
		free (user_data->files_zdan);
	}

	if ( user_data->save_new_results == TRUE ) {
		free (user_data->filename_new_results);
	}

	free (system->x);
	free (system->aw);

	for ( i = 0; i < components; i++ ) {
		free (system_description->components[i]);
	}

	if ( user_data->K != NULL ) {
		free (user_data->K);
	}

	free (system_description->components);
	free (user_data->model);
	free (user_data->filename);

}
