# include "definitions_and_headers.h"

/*
 * This function reads a dictionary-like data structure, present in "System",
 * and updates the values of the q and r pointers to their previously calculated
 * ones.
 */

void get_q_and_r ( System *system, char *comp_name, long double *q, long double *r ) {


	int i;
	double retval;

	for ( i = 0; i < KNOWN_COMPS; i++ ) {
		if ( strcmp ( system->description.uniquac_keys[i],
					comp_name) == TRUE ) {
			*r = system->description.r_vals[i];
			*q = system->description.q_vals[i];
		}
	}
}

/*
 * This function fills a dictionary-like structure with UNIFAC's R and Q values
 * from the Dortmund Data Bank (available at
 * https://www.ddbst.com/published-parameters-unifac.html, specifically the
 * section "List Of Sub Groups And Their Group Surfaces And Volumes").
 */
void fill_uniquac ( Metadata *system_description ) {

	int i;

	system_description->uniquac_keys = malloc ( KNOWN_COMPS * sizeof(char *) );
	system_description->q_vals = malloc ( KNOWN_COMPS * sizeof(double) );
	system_description->r_vals = malloc ( KNOWN_COMPS * sizeof(double) );

	/* water */
	system_description->uniquac_keys[0] = malloc ( strlen ("water") + 1);
	strcpy (system_description->uniquac_keys[0], "water");
	system_description->q_vals[i] = 1.4000;
	system_description->r_vals[i] = 0.9200;
}

