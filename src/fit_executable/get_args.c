#include "definitions_and_headers.h"

/* print accepted arguments and options */
void print_usage (void) {

	fprintf ( stderr, "Usage: FitWaterActivity [OPTIONS]\n" );
	fprintf ( stderr, "    Options:\n" );
	fprintf ( stderr, "        -f <filename>\n" );
	fprintf ( stderr, "          File in which the activity data are stored\n" );
	fprintf ( stderr, "        -F <filename to store predicted values>\n" );
	fprintf ( stderr, "          File to store the calculated values.\n" );
	fprintf ( stderr, "        -K <param 1> <param 2> <param 3> ...\n" );
	fprintf ( stderr, "          \"param n\" is the n-th previously\n" );
	fprintf ( stderr, "          obtained parameter (K_n for norrish, b_n or\n" );
	fprintf ( stderr, "          c_nm for virial and q_n/u_nn for UNIQUAC).\n" );
	fprintf ( stderr, "          passed, this option will skip the fitting\n" );
	fprintf ( stderr, "          step, and replace the parameters that would\n" );
	fprintf ( stderr, "          be obtained through regression with\n" );
	fprintf ( stderr, "          user-informed ones. This is useful to test\n" );
	fprintf ( stderr, "          predictive abilities of the regression\n" );
	fprintf ( stderr, "          (training and test data) and also to set\n" );
	fprintf ( stderr, "          the initial values of the iterative process.\n" );
	fprintf ( stderr, "        -m <model>\n" );
	fprintf ( stderr, "          Model: can be one of norrish, virial, \n" );
	fprintf ( stderr, "          UNIQUAC, caurie, raoult, zdanovskii \n" );
	fprintf ( stderr, "          or all, which compares all other models. \n" );
	fprintf ( stderr, "        -M <maximum number of iterations>\n" );
	fprintf ( stderr, "          Maximum number of iterations; defaults to \n" );
	fprintf ( stderr, "          MAX_ITER, defined as 100, if not informed. \n" );
	fprintf ( stderr, "          Useful mainly for the UNIQUAC model, which\n" );
	fprintf ( stderr, "          demands more iterations for convergence.\n" );
	fprintf ( stderr, "        -Z <1st filename> <2nd filename>\n" );
	fprintf ( stderr, "          If the model chosen is \"zdanovskii\",\n" );
	fprintf ( stderr, "          these two files must store the isotherms\n" );
	fprintf ( stderr, "          of the two components of the ternary mixture\n" );
	fprintf ( stderr, "          in a binary solution.\n" );
	fprintf ( stderr, "        -A\n" );
	fprintf ( stderr, "          Store also the predicted and real values of\n" );
	fprintf ( stderr, "          water activity in the file passed to option\n" );
	fprintf ( stderr, "          -F.\n" );
	fprintf ( stderr, "        -E\n" );
	fprintf ( stderr, "          If passed, we assume that the data marked\n" );
	fprintf ( stderr, "          \"aw\" are actually molar fraction data of\n" );
	fprintf ( stderr, "          binary solutions in osmotic equilibrium.\n" );
	fprintf ( stderr, "          This, obviously, makes sense only for n-ary\n" );
	fprintf ( stderr, "          solutions. Using with option \"-F\" will\n" );
	fprintf ( stderr, "          lead to unpredictable gibberish, for now.\n" );
	fprintf ( stderr, "        -h\n" );
	fprintf ( stderr, "          Show this help\n" );
	fprintf ( stderr, "        -O\n" );
	fprintf ( stderr, "          Fit models that need only one one file,\n" );
	fprintf ( stderr, "          currently all but \"zdanovskii\". Useful\n" );
	fprintf ( stderr, "          for combination with option \"-m all\".\n" );
	fprintf ( stderr, "        -q\n" );
	fprintf ( stderr, "          Quiet; hide verbose output\n" );

}

/* read arguments and options from user */
void getargs ( int argc, char **argv, info *user_data ) {

	int opt, gave_file, index, count, K;
	double next, *tmp_K;
	char *endptr_K;

	opterr = 0; /* This is a very ugly hack; if the user-informed
			* fitting parameters passed to the option '-K'
			* are negative, getopt will interpret any value
			* after the first as an option; this is a consequence
			* of the (other) hack we used to pass several arguments
			* to an option, that has also been used in the convert
			* source code (../convert_executable/get_args.c) to
			* supply the molar masses of the components of a solution;
			* unfortunately, here that solution won't work optimally,
			* due to the possibility of negative arguments. Luckily,
			* however, the absence of an "invalid option" warning is
			* not a problem in this situation.
			*/

	gave_file = FALSE;
	user_data->gave_filenames = FALSE;

	if ( argc < 2 ) {
		fprintf ( stderr, "No arguments provided; for help, use -h\n" );
		fprintf ( stderr, "Aborting...\n" );
		exit (24);
	}

	user_data->model = malloc ( 8 * sizeof (char) );
	strcpy ( user_data->model, "raoult" );
	user_data->quiet = FALSE;
	user_data->cost = 0;
	user_data->K = NULL;
	user_data->K_number = 0;
	user_data->is_all = FALSE;
	user_data->not_zdan = FALSE;
	user_data->files_zdan = NULL;
	user_data->save_new_results = FALSE;
	user_data->filename_new_results = NULL;
	user_data->has_aw_data = TRUE;
	user_data->max_iter = MAX_ITER;
	user_data->aw_in_results = FALSE;

	while ( ( opt = getopt ( argc, argv, "hqf:F:m:Z:M:K:OEA" ) ) != -1 ) {
		switch (opt) {
			case 'h':
				free (user_data->model);
				print_usage ();
				exit (23);
				break;
			case 'q':
				user_data->quiet = TRUE;
				break;
			case 'O':
				user_data->not_zdan = TRUE;
				break;
			case 'E':
				user_data->has_aw_data = FALSE;
				break;
			case 'A':
				user_data->aw_in_results = TRUE;
				break;
			case 'Z':
				user_data->files_zdan =
					malloc ( 2 * sizeof (char *) );
				index = optind - 1;
				count = 0;
				while ( index < argc ) {
					if ( argv[index][0] != '-' ) {
						user_data->files_zdan[count] =
							malloc ( strlen (argv[index])
									+ 1 );
						strcpy ( user_data->files_zdan[count],
								argv[index] );
						count++;
					} else break;
					index ++;
				}
				if ( count >= 2 ) {
					user_data->gave_filenames = TRUE;
				}
				break;
			case 'f':
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
					free(user_data->model);
					exit (34);
				}
				break;
			case 'F':
				user_data->save_new_results = TRUE;
				user_data->filename_new_results =
					malloc ( ( strlen (optarg) + 1 ) *
					sizeof (char) );
				strcpy ( user_data->filename_new_results,
						optarg );
				break;
			case 'm':
				free(user_data->model);
				user_data->model = malloc (
					( strlen (optarg) + 1 ) * sizeof (char)
				);
				strcpy ( user_data->model, optarg );
				if ( strcmp ( user_data->model, "all" ) == TRUE ) {
					user_data->is_all = TRUE;
					user_data->quiet = TRUE;
				}
				break;
			case 'M':
				user_data->max_iter = atoi (optarg);
				break;
			case 'K':
				/*
				* This option reads the values of the fitting
				* parameters informed by the user; this is required
				* by the analysis of test and training data
				* differences.
				*/
				index = optind - 1;
				K = 0;
				while ( index < argc ) {
					errno = 0;
					next = strtod ( argv[index], &endptr_K );
					/* Unlike the other use of this solution
					* to pass multiple arguments to a single
					* option, here we can't assume positive
					* arguments; therefore, to check if the
					* next string in *argv is an option or an
					* argument, we simply check if it's possible
					* to convert it to a double.*/
					if ( errno == 0 && *endptr_K == '\0' ) {
						K++;
						tmp_K = realloc ( user_data->K,
							K * sizeof (double) );
						if ( tmp_K == NULL ) {
							fprintf ( stderr,
								"Memory error\n" );
							exit (56);
						}
						user_data->K = tmp_K;
						user_data->K[K-1] = next;
						user_data->K_number = K;
					} else break;
					index++;
				}
				break;
		}

	}

	if ( gave_file == FALSE ) {
		fprintf ( stderr, "No filename given, aborting...\n" );
		free (user_data->model);
		exit (24);
	}

	if ( user_data->gave_filenames == FALSE &&
			strcmp ( user_data->model, "zdanovskii" ) == TRUE ) {
		fprintf ( stderr, "No isotherms given for Zdanovskii Relation." );
		fprintf ( stderr, " Aborting.\n" );
		free (user_data->model);
		if ( user_data->filename != NULL ) {
			free (user_data->filename);
		}
		if ( user_data->files_zdan != NULL ) {
			if ( user_data->files_zdan[0] != NULL ) {
				free (user_data->files_zdan[0]);
			}
			free (user_data->files_zdan);
		}
		exit (24);
	}
}

