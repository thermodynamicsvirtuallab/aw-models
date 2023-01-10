#include "definitions_and_headers.h"


/*
 * Function declarations:
 */
double arden_buck ( double T );
double ge_and_wang_bpe ( double T, info *user_data );
double ge_and_wang_fpd ( double T, info *user_data );
double boiling_point ( double P );
double delta_H_vap ( double T );
double delta_C_vap ( double T );
double freezing_point ( double P );
double delta_H_fus ( double T );
double delta_C_fus ( double T );
double scatchard_isopiestic ( double x_nacl );
double stokes_isopiestic ( double x_sucrose );


/*
 * Central function for conversion.
 */
Data convert ( Metadata *system_description, Data *system,
		info *user_data ) {

	Data new_system;
	double **x, **x_var;
	double *aw, *y_var;
	int lines, comps, i, j;
	double tmp_xtot_for_molality, tmp_xw_for_phi;
	double tmp_mw_for_mass_fraction, tmp_xw_for_mass_fraction,
		tmp_xtot_for_mass_fraction;
	double pure_pressure;

	lines = system_description->dataset_size;
	comps = system_description->n_of_comps;

	x_var = system->x;
	y_var = system->aw;

	x = malloc ( lines * sizeof (double *) );
	aw = malloc ( lines * sizeof (double) );
	for ( i = 0; i < lines; i ++ ) {
		x[i] = malloc ( comps * sizeof (double) );
	}

	/*
	* * In this part of the file we deal with the
	* * variables associated to concentration, with
	* * the conversion of the properties in molar
	* * molar fractions.
	*/

	if ( user_data->x_property != NULL ) {
		if ( strcmp ( user_data->x_property, "molality" ) == TRUE ) {
			for ( i = 0; i < lines; i++ ) {
				tmp_xtot_for_molality = MOLALITY_TO_FRACTION;
				for ( j = 0; j < comps; j++ ) {
					tmp_xtot_for_molality += x_var[i][j];
				}
				for ( j = 0; j < comps; j++ ) {
					x[i][j] = x_var[i][j] / tmp_xtot_for_molality;
				}
			}
		}
		/*
		* Here we convert molality data in molar fraction data;
		*/

		else if ( strcmp ( user_data->x_property, "mass_fraction" ) == TRUE ) {
			for ( i = 0; i < lines; i++ ) {
				tmp_mw_for_mass_fraction = 1;
				tmp_xtot_for_mass_fraction = 0;
				for ( j = 0; j < comps; j++ ) {
					tmp_mw_for_mass_fraction -= x_var[i][j];
					x[i][j] = x_var[i][j] /
						user_data->molar_masses[j];
					tmp_xtot_for_mass_fraction += x[i][j];
				}
				tmp_xw_for_mass_fraction =
					tmp_mw_for_mass_fraction / MASS_TO_FRACTION;
				tmp_xtot_for_mass_fraction += tmp_xw_for_mass_fraction;
				for ( j = 0; j < comps; j++ ) {
					x[i][j] = x[i][j] / tmp_xtot_for_mass_fraction;
				}
			}
		}
		/*
		* Here we convert mass fraction data in molar fraction data.
		*/

		else if ( strcmp ( user_data->x_property,
				"mass_concentration" ) == TRUE ) {
			for ( i = 0; i < lines; i++ ) {
				tmp_mw_for_mass_fraction = 1;
				tmp_xtot_for_mass_fraction = 0;
				for ( j = 0; j < comps; j++ ) {
					x[i][j] = x_var[i][j] /
						user_data->molar_masses[j];
					tmp_xtot_for_mass_fraction += x[i][j];
				}
				tmp_xw_for_mass_fraction =
					tmp_mw_for_mass_fraction / MASS_TO_FRACTION;
				tmp_xtot_for_mass_fraction += tmp_xw_for_mass_fraction;
				for ( j = 0; j < comps; j++ ) {
					x[i][j] = x[i][j] / tmp_xtot_for_mass_fraction;
				}
			}
		}
		/*
		* Here we convert mass concentration data in molar fraction data.
		*/

		else if ( strcmp ( user_data->x_property,
				"molar_concentration" ) == TRUE ) {
			for ( i = 0; i < lines; i++ ) {
				tmp_mw_for_mass_fraction = 1;
				tmp_xtot_for_mass_fraction = MOLALITY_TO_FRACTION;
				for ( j = 0; j < comps; j++ ) {
					tmp_xtot_for_mass_fraction += x_var[i][j];
				}
				tmp_xw_for_mass_fraction =
					tmp_mw_for_mass_fraction / MASS_TO_FRACTION;
				tmp_xtot_for_mass_fraction += tmp_xw_for_mass_fraction;
				for ( j = 0; j < comps; j++ ) {
					x[i][j] = x_var[i][j] /
						tmp_xtot_for_mass_fraction;
				}
			}
		}
		/*
		* Here we convert molar concentration data in molar fraction data.
		*/
	} else {
		for ( i = 0; i < lines; i++ ) {
			for ( j = 0; j < comps; j++ ) {
				x[i][j] = x_var[i][j];
			}
		}
	}

	/*
	* * In this other half, we proceed to convert data from boiling
	* * point elevation, freeze point depression, vapour pressure and
	* * osmotic coefficient in water activity data.
	*/

	if ( user_data->y_property != NULL ) {
		if ( strcmp ( user_data->y_property, "phi" ) == TRUE ||
			strcmp ( user_data->y_property,
			"osmotic_coefficient" ) == TRUE ) {
			for ( i = 0; i < lines; i++ ) {
				tmp_xw_for_phi = 1;
				for ( j = 0; j < comps; j++ ) {
					tmp_xw_for_phi -= x[i][j];
				}
				aw[i] = pow ( tmp_xw_for_phi, y_var[i] );
			}
		}
		/*
		* Here we convert the osmotic coefficient phi in water activity data;
		* N.b.: the x must be already molar fraction data. So this comes last.
		*/

		else if ( strcmp ( user_data->y_property,
				"vapour_pressure" ) == TRUE ||
				strcmp ( user_data->y_property,
				"vapor_pressure" ) == TRUE ) {
			if ( user_data->gave_std_pressure == TRUE ) {
				pure_pressure = user_data->std_pressure;
			} else {
				pure_pressure =
					arden_buck ( user_data->temperature );
			}
			for ( i = 0; i < lines; i++ ) {
				aw[i] =  user_data->pressure_factor * y_var[i];
				aw[i] /= pure_pressure;
			}
		}      /*
		* Here we convert vapour pressure data in water activity (via
		* ERH); the vapour pressure of pure water is approximated with
		* empirical relations, especifically the Arden Buck equation.
		*/

		else if ( strcmp ( user_data->y_property, "BPE" ) == TRUE ) {
			for ( i = 0; i < lines; i++ ) {
				aw[i] = ge_and_wang_bpe ( y_var[i], user_data );
			}
		}
		/*
		* Here we convert boiling point elevation data in water activity,
		* using the equation of Ge and Wang.
		*/

		else if ( strcmp ( user_data->y_property, "FPD" ) == TRUE ) {
			for ( i = 0; i < lines; i++ ) {
				aw[i] = ge_and_wang_fpd ( y_var[i], user_data );
			}
		}
		/*
		* Here we convert freezing point depression data in water activity,
		* using the equation of Ge and Wang.
		*/

		else if ( strcmp ( user_data->y_property,
					"isopiestic_nacl" ) == TRUE ||
				strcmp ( user_data->y_property,
					"nacl" ) == TRUE ||
				strcmp ( user_data->y_property,
					"isopiestic_sucrose" ) == TRUE ||
				strcmp ( user_data->y_property,
					"sucrose" ) == TRUE ) {
			if ( user_data->isopiestic_property == NULL ) {
			} else if ( strcmp ( user_data->isopiestic_property,
					"molality" ) == TRUE ) {
				for ( i = 0; i < lines; i++ ) {
					y_var[i] = y_var[i] /
						( y_var[i] + MOLALITY_TO_FRACTION );
				}
			} else if ( strcmp ( user_data->isopiestic_property,
					"mass_fraction" ) == TRUE ) {
				for ( i = 0; i < lines; i++ ) {
					y_var[i] = ( y_var[i] / NACL_MASS ) /
						( ( ( 1 - y_var[i] ) /
						MASS_TO_FRACTION ) +
						( y_var[i] / NACL_MASS ) );
				}
			} else if ( strcmp ( user_data->isopiestic_property,
					"mass_concentration" ) == TRUE ) {
				for ( i = 0; i < lines; i++ ) {
					y_var[i] = ( y_var[i] / NACL_MASS ) /
						( ( y_var[i] / NACL_MASS ) +
						( 1 / MASS_TO_FRACTION ) );
				}
			}
			if ( strcmp ( user_data->y_property,
					"sucrose" ) == TRUE ||
				strcmp ( user_data->y_property,
					"isopiestic_sucrose" ) == TRUE ) {
				for ( i = 0; i < lines; i++ ) {
					aw[i] = stokes_isopiestic ( y_var[i] );
				}
			} else {
				for ( i = 0; i < lines; i++ ) {
					aw[i] = scatchard_isopiestic ( y_var[i] );
				}
			}
		}
		/*
		* We use data from Scatchard et. al (1938) to obtain water
		* activity values from isopiestic measurements, i.e. from the
		* NaCl solution in osmotic equilibrium with the solution to be
		* analysed. It's also possible to use data from Stokes &
		* Robinson to obtain water activity from isopiestic sucrose
		* data.
		*/

		else if ( strcmp ( user_data->y_property,
					"molality" ) == TRUE ) {
			for ( i = 0; i < lines; i++ ) {
				tmp_xtot_for_molality = MOLALITY_TO_FRACTION;
				tmp_xtot_for_molality += y_var[i];
				aw[i] = y_var[i] / tmp_xtot_for_molality;
			}

		}
		/*
		* The reason for this seemingly absurd is defined in the file
		* ../fit_executable/models/virial.c, at the function "phi_virial".
		*/

	} else {
		for ( i = 0; i < lines; i++ ) {
			aw[i] = y_var[i];
		}
	}

	/*
	* One last sanity check: we need this because sometimes floating
	* point and empirical relations precision makes a result a little
	* bit wrong, and this should be corrected.
	*/

	for ( i = 0; i < lines; i++ ) {
		if ( user_data->y_property != NULL ) {
			if ( aw[i] > 1.0 ) {
				aw[i] = 1.0;
			} else if ( aw[i] < 0.0 ) {
				aw[i] = 0.0;
			}
		}
		if ( user_data->x_property != NULL ) {
			for ( j = 0; j < comps; j++ ) {
				if ( x[i][j] > 1.0 ) {
					x[i][j] = 1.0;
				} else if ( x[i][j] < 0.0 ) {
					x[i][j] = 0.0;
				}
			}
		}
	}

	new_system.x = x;
	new_system.aw = aw;

	return new_system;
}

/*
 * *
 * * Auxiliary functions
 * *
 */

/* For vapour pressure calculation, more precise than Antoine. */
double arden_buck ( double T ) {
	if ( T > 150.0 ) {
		T = T - 273.15;
	}
	return ( A_A * exp ( ( B_A - T / C_A ) * ( T / ( D_A + T ) ) ) );
}

/* Ge and Wang's relation between boiling point elevation and water activity. */
double ge_and_wang_bpe ( double T, info *user_data ) {

	double aw, aw_log;
	double T_B;

	T_B = boiling_point (user_data->pressure);

	aw_log = - delta_H_vap ( T_B ) * ( ( 1 / T_B ) - ( 1 / T ) );
	aw_log += delta_C_vap ( T ) * ( log ( T / T_B ) - ( ( T - T_B ) / T ) );
	aw_log = aw_log / R_CONST;

	aw = exp ( aw_log );

	return aw;
}

/* Ge and Wang's relation between freezing point depression and water activity. */
double ge_and_wang_fpd ( double T, info *user_data ) {

	double aw, aw_log;
	double T_F;

	T_F = freezing_point (user_data->pressure);

	aw_log = delta_H_fus ( T_F ) * ( ( 1 / T_F ) - ( 1 / T ) );
	aw_log += delta_C_fus ( T ) * ( log ( T / T_F ) - ( ( T_F - T ) / T ) );
	aw_log = aw_log / R_CONST;

	aw = exp ( aw_log );

	return aw;
}

/*
 * Basically the inverse of Antoine.
 */
double boiling_point ( double P ) {

	double T_B;

	T_B = - C_B + B_B / ( A_B - log (P) );

	return T_B;
}

double freezing_point ( double P ) {

	return FREEZING_POINT_WATER;
}

/*
 * The next four functions are polynomial regressions to deal with the relation
 * of some property (such as C_p or Delta_H) with the temperature.
 */

double delta_H_vap ( double T ) {

	double delta_H;

	delta_H = A_H + B_H * T + C_H * pow ( T, 2 ) + D_H * pow ( T, 3 );

	return delta_H;
}

double delta_C_vap ( double T ) {

	double delta_C, C_vap, C_liq;

	C_vap = A_G + B_G * T + C_G * pow ( T, 2 ) + D_G * pow ( T, 3 );
	C_liq = A_L + B_L * T + C_L * pow ( T, 2 ) +
		D_L * pow ( T, 3 ) + E_L * pow ( T, 4 );
	delta_C = C_vap - C_liq;

	return delta_C;
}

double delta_H_fus ( double T ) {

	return DELTA_H_FUS;
}

double delta_C_fus ( double T ) {

	double delta_C, C_liq, C_sol;

	C_liq = A_L + B_L * T + C_L * pow ( T, 2 ) +
		D_L * pow ( T, 3 ) + E_L * pow ( T, 4 );
	C_sol = A_S + B_S * T + C_S * pow ( T, 2 ) + D_S * pow ( T, 3 );
	delta_C = C_liq - C_sol;

	return delta_C;
}

/*
 * This function converts data of molar fraction of NaCl in binary solution
 * to water activity data.
 */

double scatchard_isopiestic ( double x_nacl ) {

	double aw_equiv;

	aw_equiv = 0;
	aw_equiv += A_I * pow ( x_nacl, 4 );
	aw_equiv += B_I * pow ( x_nacl, 3 );
	aw_equiv += C_I * pow ( x_nacl, 2 );
	aw_equiv += D_I * x_nacl;
	aw_equiv += E_I;

	return aw_equiv;

}

/*
 * This function converts data of molar fraction of sucrose in binary solution
 * to water activity data.
 */

double stokes_isopiestic ( double x_sucrose ) {

	double aw_equiv;

	aw_equiv = 0;
	aw_equiv += A_SR * pow ( x_sucrose, 8 );
	aw_equiv += B_SR * pow ( x_sucrose, 7 );
	aw_equiv += C_SR * pow ( x_sucrose, 6 );
	aw_equiv += D_SR * pow ( x_sucrose, 5 );
	aw_equiv += E_SR * pow ( x_sucrose, 4 );
	aw_equiv += F_SR * pow ( x_sucrose, 3 );
	aw_equiv += G_SR * pow ( x_sucrose, 2 );
	aw_equiv += H_SR * x_sucrose;
	aw_equiv += I_SR;

	return aw_equiv;

}
