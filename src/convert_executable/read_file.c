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
	char tmpchar, str[256];
	char **composition;
	double **x_values;
	double *y_values;
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
		y_values[i] = strtod ( str, NULL ) + user_data->temperature_to_add;
		for ( j = 0; j < columns; j++ ) {
			fscanf(file, "%127[^,\n]", str);
			fscanf(file, "%*c");
			x_values[i][j] = strtod ( str, NULL );
		}
	}

	fclose ( file );
		/*
		* And now we have read the data from the file
		*/

	system_description->dataset_size = lines;
	system_description->n_of_comps = columns;
	system_description->components = composition;

	system->x = x_values;
	system->aw = y_values;
		/*
		* And now our data and metadata are stored in our structs
		*/

}

/*
* This function cleans up the memory utilized,
* freeing pointers to data and metadata.
*/

void finalize ( Metadata *system_description, Data *new_system,
		Data *system, info *user_data ) {

	int i, lines = system_description->dataset_size;
	int components = system_description->n_of_comps;

	for ( i = 0; i < lines; i++ ) {
		free (system->x[i]);
		free (new_system->x[i]);
	}

	free (system->x);
	free (system->aw);
	free (new_system->x);
	free (new_system->aw);

	for ( i = 0; i < components; i++ ) {
		free (system_description->components[i]);
	}
	free (system_description->components);

	if ( user_data->filename != NULL ) {
		free (user_data->filename);
	}

	if ( user_data->new_filename != NULL ) {
		free (user_data->new_filename);
	}

	if ( user_data->y_property != NULL ) {
		free (user_data->y_property);
	}

	if ( user_data->x_property != NULL ) {
		free (user_data->x_property);
	}

	if ( user_data->molar_masses != NULL ) {
		free (user_data->molar_masses);
	}

	if ( user_data->isopiestic_property != NULL ) {
		free (user_data->isopiestic_property);
	}


}
