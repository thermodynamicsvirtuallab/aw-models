# include "definitions_and_headers.h"
# include "unifac_header.h"

/*
 * This function fills a dictionary-like structure with UNIFAC's R and Q values
 * from the Dortmund Data Bank (available at
 * https://www.ddbst.com/published-parameters-unifac.html, specifically the
 * section "List Of Sub Groups And Their Group Surfaces And Volumes").
 */
void fill_uniquac ( Metadata *system_description ) {

	int i;

	system_description->q_vals = malloc ( system_description->n_of_comps *
						sizeof(double) );
	system_description->r_vals = malloc ( system_description->n_of_comps *
						sizeof(double) );

	for ( i = 0; i < system_description->n_of_comps; i++ ) {
		if ( strcmp(system_description->components[i], "water") ){
			system_description->r_vals[i] = R_WATER;
			system_description->q_vals[i] = Q_WATER;
		}
		if ( strcmp(system_description->components[i],
			"1_2_3_4_butanetetrol") ||
				strcmp(system_description->components[i],
					"erythritol") ){
			system_description->r_vals[i] = R_ERYTHRITOL;
			system_description->q_vals[i] = Q_ERYTHRITOL;
		}
		if ( strcmp(system_description->components[i],
			"1_2_3_propanetriol") ||
				strcmp(system_description->components[i],
					"glycerol") ){
			system_description->r_vals[i] = R_GLYCEROL;
			system_description->q_vals[i] = Q_GLYCEROL;
		}
		if ( strcmp(system_description->components[i],
					"1_2_ethanediol") ){
			system_description->r_vals[i] = R_1_2_ETHANEDIOL;
			system_description->q_vals[i] = Q_1_2_ETHANEDIOL;
		}
		if ( strcmp(system_description->components[i],
			"alanine") ||
				strcmp(system_description->components[i],
					"dl_alanine") ){
			system_description->r_vals[i] = R_ALANINE;
			system_description->q_vals[i] = Q_ALANINE;
		}
		if ( strcmp(system_description->components[i],
			"alpha_amino_n_butyric_acid") ||
				strcmp(system_description->components[i],
					"dl_alpha_aminobutyric_acid") ){
			system_description->r_vals[i] =
				R_ALPHA_AMINO_N_BUTYRIC_ACID;
			system_description->q_vals[i] =
				Q_ALPHA_AMINO_N_BUTYRIC_ACID;
		}
		if ( strcmp(system_description->components[i],
			"arabinose") ||
				strcmp(system_description->components[i],
					"ribose") ){
			system_description->r_vals[i] = R_ARABINOSE;
			system_description->q_vals[i] = Q_ARABINOSE;
		}
		if ( strcmp(system_description->components[i],
			"arginine") ||
				strcmp(system_description->components[i],
					"l_arginine") ){
			system_description->r_vals[i] = R_ARGININE;
			system_description->q_vals[i] = Q_ARGININE;
		}
		if ( strcmp(system_description->components[i],
			"dextrose") ||
				strcmp(system_description->components[i],
					"d_glucose") ||
				strcmp(system_description->components[i],
					"glucose") ||
				strcmp(system_description->components[i],
					"mannose") ||
				strcmp(system_description->components[i],
					"galactose") ){
			system_description->r_vals[i] = R_GLUCOSE;
			system_description->q_vals[i] = Q_GLUCOSE;
		}
		if ( strcmp(system_description->components[i],
			"d_fructose") ||
				strcmp(system_description->components[i],
					"fructose") ){
			system_description->r_vals[i] = R_FRUCTOSE;
			system_description->q_vals[i] = Q_FRUCTOSE;
		}
		if ( strcmp(system_description->components[i],
					"dl_methionine") ){
			system_description->r_vals[i] = R_DL_METHIONINE;
			system_description->q_vals[i] = Q_DL_METHIONINE;
		}
		if ( strcmp(system_description->components[i],
			"valine") ||
				strcmp(system_description->components[i],
					"dl_valine") ){
			system_description->r_vals[i] = R_VALINE;
			system_description->q_vals[i] = Q_VALINE;
		}
		if ( strcmp(system_description->components[i],
			"d_mannitol") ||
				strcmp(system_description->components[i],
					"mannitol") ||
				strcmp(system_description->components[i],
					"sorbitol") ){
			system_description->r_vals[i] = R_MANNITOL;
			system_description->q_vals[i] = Q_MANNITOL;
		}
		if ( strcmp(system_description->components[i],
					"glutamic_acid") ){
			system_description->r_vals[i] = R_GLUTAMIC_ACID;
			system_description->q_vals[i] = Q_GLUTAMIC_ACID;
		}
		if ( strcmp(system_description->components[i],
					"glycine") ){
			system_description->r_vals[i] = R_GLYCINE;
			system_description->q_vals[i] = Q_GLYCINE;
		}
		if ( strcmp(system_description->components[i],
					"glycylglycine") ){
			system_description->r_vals[i] = R_GLYCYLGLYCINE;
			system_description->q_vals[i] = Q_GLYCYLGLYCINE;
		}
		if ( strcmp(system_description->components[i],
					"histidine") ){
			system_description->r_vals[i] = R_HISTIDINE;
			system_description->q_vals[i] = Q_HISTIDINE;
		}
		if ( strcmp(system_description->components[i],
					"hydroxy_l_proline") ){
			system_description->r_vals[i] = R_HYDROXY_L_PROLINE;
			system_description->q_vals[i] = Q_HYDROXY_L_PROLINE;
		}
		if ( strcmp(system_description->components[i],
					"lactamide") ){
			system_description->r_vals[i] = R_LACTAMIDE;
			system_description->q_vals[i] = Q_LACTAMIDE;
		}
		if ( strcmp(system_description->components[i],
					"lactose") ){
			system_description->r_vals[i] = R_LACTOSE;
			system_description->q_vals[i] = Q_LACTOSE;
		}
		if ( strcmp(system_description->components[i],
					"l_isoleucine") ){
			system_description->r_vals[i] = R_L_ISOLEUCINE;
			system_description->q_vals[i] = Q_L_ISOLEUCINE;
		}
		if ( strcmp(system_description->components[i],
					"l_proline") ){
			system_description->r_vals[i] = R_L_PROLINE;
			system_description->q_vals[i] = Q_L_PROLINE;
		}
		if ( strcmp(system_description->components[i],
			"serine") ||
				strcmp(system_description->components[i],
					"l_serine") ){
			system_description->r_vals[i] = R_SERINE;
			system_description->q_vals[i] = Q_SERINE;
		}
		if ( strcmp(system_description->components[i],
					"l_threonine") ){
			system_description->r_vals[i] = R_L_THREONINE;
			system_description->q_vals[i] = Q_L_THREONINE;
		}
		if ( strcmp(system_description->components[i],
					"lysine") ){
			system_description->r_vals[i] = R_LYSINE;
			system_description->q_vals[i] = Q_LYSINE;
		}
		if ( strcmp(system_description->components[i],
					"maltitol") ){
			system_description->r_vals[i] = R_MALTITOL;
			system_description->q_vals[i] = Q_MALTITOL;
		}
		if ( strcmp(system_description->components[i],
			"maltose") ||
				strcmp(system_description->components[i],
					"maltose_monohydrate") ){
			system_description->r_vals[i] = R_MALTOSE;
			system_description->q_vals[i] = Q_MALTOSE;
		}
		if ( strcmp(system_description->components[i],
					"raffinose") ){
			system_description->r_vals[i] = R_RAFFINOSE;
			system_description->q_vals[i] = Q_RAFFINOSE;
		}
		if ( strcmp(system_description->components[i],
					"sucrose") ){
			system_description->r_vals[i] = R_SUCROSE;
			system_description->q_vals[i] = Q_SUCROSE;
		}
		if ( strcmp(system_description->components[i],
					"tris_hydroxymethyl_aminomethane") ){
			system_description->r_vals[i] =
				R_TRIS_HYDROXYMETHYL_AMINOMETHANE;
			system_description->q_vals[i] =
				Q_TRIS_HYDROXYMETHYL_AMINOMETHANE;
		}
		if ( strcmp(system_description->components[i],
					"urea") ){
			system_description->r_vals[i] = R_UREA;
			system_description->q_vals[i] = Q_UREA;
		}
		if ( strcmp(system_description->components[i],
					"xylitol") ){
			system_description->r_vals[i] = R_XYLITOL;
			system_description->q_vals[i] = Q_XYLITOL;
		}
		if ( strcmp(system_description->components[i],
					"xylose") ){
			system_description->r_vals[i] = R_XYLOSE;
			system_description->q_vals[i] = Q_XYLOSE;
		}
	}
}

