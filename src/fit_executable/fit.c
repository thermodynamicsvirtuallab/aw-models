#include "definitions_and_headers.h"

/*
 * NOTICE:
 *
 *
 * A significant portion of this file was adapted from the examples provided in
 * the GNU Scientific Library Reference Manual.
 */

int fit_to_model ( System *data, info *user_data ) {

	size_t n, p;
	const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_workspace *w;
	gsl_multifit_nlinear_fdf fdf;
	gsl_multifit_nlinear_parameters fdf_params =
		gsl_multifit_nlinear_default_parameters ();
	gsl_vector *f;
	gsl_matrix *J, *covar;
	gsl_vector_view x;
	callback_function callback;
	print_function print;
	save_function save;
	const double xtol = 1e-9;
	const double gtol = 1e-9;
	const double ftol = 0.0;
	double *x_init;
	double chisq, chisq0;
	int status, info;

	fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;

	n = data->description.dataset_size;
	p = data->description.n_of_comps;

	save = NULL;
	if ( strcmp ( user_data->model, "norrish" ) == TRUE ) {
		fdf.f = phi_norrish;
		fdf.p = p;
		callback = &callback_norrish;
		print = &print_norrish;
		save = &save_norrish;
	} else if ( strcmp ( user_data->model, "virial" ) == TRUE ) {
		fdf.f = phi_virial;
		fdf.p = p + ( p * ( p - 1 ) ) / 2;
		callback = &callback_virial;
		print = &print_virial;
		save = &save_virial;
	} else if ( strcmp ( user_data->model, "uniquac" ) == TRUE ) {
		fdf.f = phi_uniquac;
		fdf.p = 2 * ( p + 1 );
		callback = &callback_uniquac;
		print = &print_uniquac;
		save = &save_uniquac;
	} else {
		fprintf (stderr, "Model unknown. Aborting...\n" );
		exit (45);
	}

	fdf.df = NULL;
	fdf.fvv = NULL;
		/*
		* We won't specify the values of the jacobians and hessians,
		* because the data are small enough to make the performance
		* difference negligible
		*/
	fdf.n = n;
	fdf.params = data;

	p = fdf.p;

	if ( n < p ) {
		fprintf ( stderr,
			"Excessive number of parameters/ dataset too small.\n" );
		fprintf ( stderr,
			"Aborting...\n" );
		return 0;
	}

	covar = gsl_matrix_alloc ( p, p );
	x_init = malloc ( p * sizeof (double) );
	init_data ( user_data->model, x_init, p, user_data );

	x = gsl_vector_view_array ( x_init, p );

	w = gsl_multifit_nlinear_alloc ( T, &fdf_params, n, p );
	gsl_multifit_nlinear_winit ( &x.vector, NULL, &fdf, w );

	f = gsl_multifit_nlinear_residual (w);
	gsl_blas_ddot ( f, f, &chisq0 );

	if ( user_data->quiet == TRUE ) {
		callback = NULL;
		/* limit ouput if required */
	}

	status = gsl_multifit_nlinear_driver ( user_data->max_iter, xtol,
			gtol, ftol, callback, NULL, &info, w);

	J = gsl_multifit_nlinear_jac (w);
	gsl_multifit_nlinear_covar ( J, 0.0, covar );

	gsl_blas_ddot ( f, f, &chisq );

	print ( covar, w, status, chisq0, chisq, data, user_data );

	if ( user_data->save_new_results == TRUE && save != NULL ) {
		save ( data, user_data, w );
	}

	gsl_multifit_nlinear_free (w) ;
	gsl_matrix_free (covar);
	free (x_init);

	return 0;
}

