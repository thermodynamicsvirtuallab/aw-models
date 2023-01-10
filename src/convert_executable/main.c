#include "definitions_and_headers.h"

int main ( int argc, char **argv ) {

	Metadata system_description;
	Data system, new_system;
	info user_data;

	getargs ( argc, argv, &user_data );
	/* read command line arguments */

	initialize ( user_data.filename, &system_description,
			&system, &user_data );
			/* read file to data structure */

	new_system = convert ( &system_description, &system, &user_data );
	/* convert old properties to the desired ones,
	* storing them in a new data structure
	*/

	write_to_file ( &system_description, &new_system, &user_data );
	/* write new data structure to file */

	finalize ( &system_description, &new_system, &system, &user_data );
	/* free remaining malloc'd blocks */

	return 0;

}
