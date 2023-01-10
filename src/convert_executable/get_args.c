#include "definitions_and_headers.h"

void print_usage (void) {

	fprintf ( stderr, "Usage: ConvertWaterActivity [OPTIONS]\n" );
	fprintf ( stderr, "  Options:\n" );
	fprintf ( stderr, "    -h\n" );
	fprintf ( stderr, "      Show this help.\n" );
	fprintf ( stderr, "    -q\n" );
	fprintf ( stderr, "      Quiet; hide verbose output.\n" );
	fprintf ( stderr, "    -f <filename>\n" );
	fprintf ( stderr, "      File in which the old data are stored.\n" );
	fprintf ( stderr, "    -n <new filename>\n" );
	fprintf ( stderr, "      File in which the new data shall be stored.\n" );
	fprintf ( stderr, "    -x <property>\n" );
	fprintf ( stderr, "      Property to convert to molar fraction.\n" );
	fprintf ( stderr, "      Can be one of \"molality\", \"mass_fraction\",\n" );
	fprintf ( stderr, "      \"mass_concentration\"" );
	fprintf ( stderr,      " or \"molar_concentration\".\n" );
	fprintf ( stderr, "      If set, the following options are valid:\n" );
	fprintf ( stderr, "        -m <molar mass 1> ... <molar mass n>\n" );
	fprintf ( stderr, "          If the property replacing molar fraction\n" );
	fprintf ( stderr, "          is mass fraction, the the molar masses of\n" );
	fprintf ( stderr, "          each component shall be informed through\n" );
	fprintf ( stderr, "          this option, in g/mol.\n" );
	fprintf ( stderr, "    -y <property>\n" );
	fprintf ( stderr, "      Property to convert to water activity.\n" );
	fprintf ( stderr, "      Can be one of \"vapour_pressure\", \"BPE\",\n" );
	fprintf ( stderr, "      \"FPD\", \"phi\" or \"osmotic_coefficient\",\n" );
	fprintf ( stderr, "      \"isopiestic_nacl\" (or \"nacl\"), \n" );
	fprintf ( stderr, "      \"isopiestic_sucrose\" (or \"sucrose\") or \n" );
	fprintf ( stderr, "      \"molality\".\n" );
	fprintf ( stderr, "      If set, the following options are valid:\n" );
	fprintf ( stderr, "        -k\n" );
	fprintf ( stderr, "          Indicates pressure in kPa, not Pa.\n" );
	fprintf ( stderr, "        -M\n" );
	fprintf ( stderr, "          Indicates pressure in mmHg, not Pa.\n" );
	fprintf ( stderr, "        -t <temperature>\n" );
	fprintf ( stderr, "          If the property replacing water activity\n" );
	fprintf ( stderr, "          is the vapour pressure (-y pressure) of\n" );
	fprintf ( stderr, "          the solution, then its temperature shall\n" );
	fprintf ( stderr, "          be informed (see also: -T). Otherwise the\n" );
	fprintf ( stderr, "          program assumes T = 298.15 K.\n" );
	fprintf ( stderr, "        -P <standard pressure in Pa>\n" );
	fprintf ( stderr, "          If informed the vapour pressure won't\n" );
	fprintf ( stderr, "          be obtained through Arden-Buck, it will,\n" );
	fprintf ( stderr, "          instead, be defined by the user.\n" );
	fprintf ( stderr, "        -C\n" );
	fprintf ( stderr, "          If passed, the program will assume\n" );
	fprintf ( stderr, "          temperatures in oC; otherwise, in K.\n" );
	fprintf ( stderr, "        -p <pressure>\n" );
	fprintf ( stderr, "          If the property replacing water activity is\n" );
	fprintf ( stderr, "          freezing point depression or boiling point\n" );
	fprintf ( stderr, "          elevation, the pressure shall be informed\n" );
	fprintf ( stderr, "          through this option. Otherwise assumed to be\n" );
	fprintf ( stderr, "          P = 101.325 kPa.\n" );
	fprintf ( stderr, "        -I <property>\n" );
	fprintf ( stderr, "          Property in which isopiestic NaCl data\n" );
	fprintf ( stderr, "          are stored. The options are the same for\n" );
	fprintf ( stderr, "          \"-x\", above. If not informed, will assume\n" );
	fprintf ( stderr, "          molar fraction.\n" );

}

void getargs ( int argc, char **argv, info *user_data ) {

	int opt, gave_file, gave_new_file, gave_C_yet, index, comps;
	double next;
	double *tmp_mol;

	gave_file = FALSE;
	gave_new_file = FALSE;
	gave_C_yet = FALSE;

	user_data->filename = NULL;
	user_data->new_filename = NULL;
	user_data->y_property = NULL;
	user_data->temperature = 298.15;
	user_data->temperature_to_add = 0;
	user_data->pressure = 101325;
	user_data->pressure_factor = 1;
	user_data->x_property = NULL;
	user_data->molar_masses = NULL;
	user_data->quiet = FALSE;
	user_data->isopiestic_property = NULL;
	user_data->gave_std_pressure = FALSE;

	if ( argc < 2 ) {
		fprintf ( stderr, "No arguments provided; for help, use -h\n" );
		fprintf ( stderr, "Aborting...\n" );
		exit (24);
	}

	user_data->quiet = FALSE;

	while ( ( opt = getopt ( argc, argv, "hqCkMf:n:x:I:m:y:t:p:P:" ) ) != -1 ) {
		switch (opt) {
			case 'h':
				print_usage ();
				exit (23);
				break;
			case 'q':
				user_data->quiet = TRUE;
				break;
			case 'f':
				/* read filename from user (original data) */
				gave_file = TRUE;
				if ( access( optarg, F_OK|R_OK ) == TRUE ) {
					user_data->filename = malloc (
						( strlen (optarg) + 1 ) *
						sizeof (char)
					);
					strcpy ( user_data->filename, optarg );
				} else {
					fprintf ( stderr,
					"File not found or wrong permissions. " );
					fprintf ( stderr, "Aborting...\n" );
					free(user_data->filename);
					exit (34);
				}
				break;
			case 'k':
				user_data->pressure_factor = 1000;
				break;
			case 'M':
				user_data->pressure_factor = 133.322387415;
				break;
			case 'n':
				/* read new filename from user (new data) */
				gave_new_file = TRUE;
				user_data->new_filename = malloc (
					( strlen (optarg) + 1 ) *
					sizeof (char) );
				strcpy ( user_data->new_filename, optarg );
				break;
			case 'x':
				/* read property to convert to molar fraction */
				user_data->x_property = malloc (
					( strlen (optarg) + 1 ) *
					sizeof (char) );
				strcpy ( user_data->x_property, optarg );
				break;
			case 'I':
				/* read property to convert to molar fraction */
				user_data->isopiestic_property = malloc (
					( strlen (optarg) + 1 ) *
					sizeof (char) );
				strcpy ( user_data->isopiestic_property, optarg );
				break;
			case 'm':
				/* read molar mass of each component */
				index = optind - 1;
				comps = 0;
				while ( index < argc ) {
					if ( argv[index][0] != '-' ) {
						next = strtod ( argv[index], NULL );
						comps++;
						tmp_mol =
							realloc (
							user_data->molar_masses,
							comps * sizeof (double) );
						if ( tmp_mol == NULL ) {
							fprintf ( stderr,
								"Memory error\n" );
							exit (56);
						}
						user_data->molar_masses = tmp_mol;
						user_data->molar_masses[comps-1]
							= next;
					} else break;
					index++;
				}
				break;
			case 'y':
				/* read property to convert to water activity */
				user_data->y_property = malloc (
					( strlen (optarg) + 1 ) *
					sizeof (char) );
				strcpy ( user_data->y_property, optarg );
				break;
			case 't':
				user_data->temperature = strtod ( optarg, NULL );
				break;
			case 'C':
				/* read if temperature in C or K */
				gave_C_yet = TRUE;
				break;
			case 'p':
				user_data->pressure = strtod ( optarg, NULL );
				break;
			case 'P':
				user_data->std_pressure = strtod ( optarg, NULL );
				user_data->gave_std_pressure = FALSE;
				break;
		}

	}

	if ( gave_file == FALSE || gave_new_file == FALSE ) {
		fprintf ( stderr, "No filename (or only one) given, aborting...\n" );
		free (user_data->filename);
		free (user_data->new_filename);
		free (user_data->y_property);
		free (user_data->x_property);
		free (user_data->molar_masses);
		exit (24);
	}
	if ( gave_C_yet == TRUE ) {
		/* convert Kelvin (default) to Celsius */
		user_data->temperature_to_add = KELVIN_TO_CELSIUS;
	}
}

