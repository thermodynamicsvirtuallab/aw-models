#include "definitions_and_headers.h"

void write_to_file ( Metadata *system_description, Data *new_system,
		info *user_data ) {

	int i, j, lines, comps;
	char *filename;
	FILE *new_file;

	filename = user_data->new_filename;

	new_file = fopen ( filename, "w" );

	lines = system_description->dataset_size;
	comps = system_description->n_of_comps;

	fprintf ( new_file, "aw," );

	for ( i = 0; i < comps - 1; i++ ) {
		fprintf ( new_file, "%s,", system_description->components[i] );
	}
	/* write component names to the first line */

	fprintf ( new_file, "%s\n", system_description->components[comps-1] );

	for ( i = 0; i < lines; i++ ) {
		fprintf ( new_file, "%.9f,", new_system->aw[i] );
		for ( j = 0; j < comps - 1; j++ ) {
			fprintf ( new_file, "%.9f,", new_system->x[i][j] );
		}
		fprintf ( new_file, "%.9f\n", new_system->x[i][j] );
		/* write obtained water activities and molar fractions */
	}

	fclose (new_file);

}

