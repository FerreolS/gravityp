/* $Id: gravi_calib.c,v 1.10 2012/03/23 15:10:40 nazouaoui Exp $
 *
 * This file is part of the GRAVI Pipeline
 * Copyright (C) 2002,2003 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */


/* 
 * This file contains outdated function, not compiled
 * nor included in H files.
 */


int gravi_get_base_from_tel (int t1, int t2)
{
  for ( int i=0 ; i<6; i++ ) {
	if ( GRAVI_BASE_TEL[i][0] == t1 && GRAVI_BASE_TEL[i][1] == t2 ) return i;
	if ( GRAVI_BASE_TEL[i][0] == t2 && GRAVI_BASE_TEL[i][1] == t1 ) return i;
  }
  return 99;
}

int gravi_get_basesign_from_tel (int t1, int t2)
{
  for ( int i=0 ; i<6; i++ ) {
	if ( GRAVI_BASE_TEL[i][0] == t1 && GRAVI_BASE_TEL[i][1] == t2 ) return   1;
	if ( GRAVI_BASE_TEL[i][0] == t2 && GRAVI_BASE_TEL[i][1] == t1 ) return  -1;
  }
  return 99;
}

cpl_error_code gravi_lazer_get_wavelength (gravi_data * lazer_data){
	gravi_msg_function_start(1);

	cpl_table * spectrum_table, * detector_table, * wavelength_table;
	char  * polarisation, * data_x;
	const char * regname;
	cpl_propertylist * primary_hdr;
	int type_data, pol, npol, nregion, reg, nwave, size, wave, nv;

	primary_hdr = gravi_data_get_header (lazer_data);

	for (type_data = 0; type_data < 2; type_data ++) {

		/* Extract the IMAGING_DETECTOR table, OI_WAVELENGTH property list and the
		 * primary header */

		if (type_data == 0){
			detector_table = gravi_data_get_table (lazer_data,
													  GRAVI_IMAGING_DETECTOR_SC_EXT);
			/* Find the spectrum data field */
			spectrum_table = gravi_data_get_table (lazer_data,
														 GRAVI_SPECTRUM_DATA_SC_EXT);
		}
		else {
			detector_table = gravi_data_get_table (lazer_data,
													  GRAVI_IMAGING_DETECTOR_FT_EXT);

			/* Find the spectrum data field */
			spectrum_table = gravi_data_get_table (lazer_data,
														 GRAVI_SPECTRUM_DATA_FT_EXT);
		}

		if ((spectrum_table == NULL) || (primary_hdr == NULL) ||
				(wavelength_table == NULL) || (detector_table == NULL)){

			return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
									  "The preproc data must contrain "
												"SPECTRUM_DATA field");
		}


		/* Get the number of wavelength and the region */
		nwave = cpl_table_get_column_dimension(spectrum_table, "DATA1", 1);
		size = cpl_table_get_column_dimension(spectrum_table, "DATA1", 0);

		nregion = cpl_table_get_nrow(detector_table);
		if (nregion > 24)
			npol = 2;
		else
			npol = 1;


		wavelength_table = gravi_data_get_oi_table(lazer_data,
						GRAVI_OI_WAVE_EXT(type_data),
						GRAVI_INSNAME(type_data,0,npol));

		for (pol = 0; pol < npol; pol++){

			if (pol == 0)
				if (npol == 2)
					polarisation = POLAR_2;
				else {
					polarisation = POLAR_3;
					pol = 1;
				}
			else
				polarisation = POLAR_1;

			cpl_image * mean_pol = cpl_image_new(size, nwave, CPL_TYPE_DOUBLE);
			cpl_image_fill_window (mean_pol, 1, 1, size, nwave, 0.0);

			for (reg = 0; reg < nregion; reg ++){

				regname = cpl_table_get_string (detector_table, "REGANAME", reg);

				if (!strstr(regname, polarisation))
					continue;


				data_x = cpl_sprintf("DATA%d", reg+1);

				cpl_imagelist * spec_data = gravi_table_data_to_imagelist(
						spectrum_table, reg+1);

				cpl_image * mean_spec = cpl_imagelist_collapse_create(spec_data);

				cpl_image_add (mean_pol, mean_spec);
				cpl_imagelist_delete (spec_data);
				cpl_image_delete (mean_spec);

			}

			int ind = 0;
			double val_spec = cpl_image_get (mean_pol, 1, 1, &nv);
			for (wave = 1; wave < nwave; wave ++){
				if (val_spec < cpl_image_get (mean_pol, 1, wave, &nv)){
					val_spec = cpl_image_get (mean_pol, 1, wave, &nv);
					ind = wave;
				}

			}

			float wave_lazer = cpl_table_get_float (wavelength_table,
					"EFF_WAVE", ind, &nv);
			char * key = cpl_sprintf("LASER %s %s", GRAVI_TYPE(type_data), polarisation);
			cpl_propertylist_append_float (primary_hdr, key, wave_lazer);
			cpl_free (key);

			cpl_image_delete (mean_pol);
		}
	}
	
	/* Verbose */
	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}


gravi_data * gravi_compute_disp_old  (gravi_data * vis_data){

	cpl_table * vistable_sc, * wavelenght_sc, * fddl_met;
	cpl_propertylist * plist, * vis_FT_p, * vis_SC_p, * wavesc_plist;
	gravi_data * disp_data;
	int wave, nwave_sc, base, nbase = 6, row , nfile, nv;
	cpl_array * phase_sc_base;
	const cpl_array * phase_sc,
				* phase_met1, * phase_met2, * fddl_sc, * fddl_ft;
	cpl_vector * wave_sc;
	cpl_bivector * fout, * fref;
	int tel_1[6] = {0,0,0,1,1,2}, pol = 0;
	int tel_2[6] = {1,2,3,2,3,3}, size_matrix, size_opd, npol;
	double pos_sc, pos_ft, pos_mean, ph_sc, opd_sc, opd_met,
	ph_met1, ph_met2;
	double time_sc, time_ft, exptime_sc;
	gsl_matrix * phase_disp;
	gsl_vector ** opd_disp;
	char * pol_FT_check[2], * pol_SC_check[2];

	/* Verbose */
	gravi_msg_function_start(1);

	vis_FT_p = gravi_data_get_plist (vis_data, GRAVI_OI_VIS_FT_EXT);
	vis_SC_p = gravi_data_get_plist (vis_data, GRAVI_OI_VIS_SC_EXT);

	/* Check if there is 2 polarization */
	if (!(strcmp(gravi_pfits_get_insname (vis_FT_p), INSNAME_FT_P1)) &&
			(strcmp(gravi_pfits_get_insname (vis_FT_p), INSNAME_FT_P2))){
		pol_FT_check[0] = INSNAME_FT_P1;
		pol_FT_check[1] = INSNAME_FT_P2;

		npol = 2;
//		polFT_test = 1;
	}
	else if (!strcmp(gravi_pfits_get_insname (vis_FT_p), INSNAME_FT)){

		pol_FT_check[0] = INSNAME_FT;
		pol_FT_check[1] = INSNAME_FT;
		npol = 1;
	}

	if (!(strcmp(gravi_pfits_get_insname (vis_SC_p), INSNAME_SC_P1)) &&
			(strcmp(gravi_pfits_get_insname (vis_SC_p), INSNAME_SC_P2))){
		pol_SC_check[0] = INSNAME_SC_P1;
		pol_SC_check[1] = INSNAME_SC_P2;
		npol = 2;
//		polSC_test = 1;
	}
	else if (!strcmp(gravi_pfits_get_insname (vis_SC_p), INSNAME_SC)){
		pol_SC_check[0] = INSNAME_SC;
		pol_SC_check[1] = INSNAME_SC;
		npol = 1;
	}

	/* Get the FDDL_MET table and load the matrix fddl matrix*/
	fddl_met = gravi_data_get_table(vis_data, "FDDL_MET_MEAN");
	nfile=cpl_table_get_nrow(fddl_met);

	disp_data = gravi_data_new(0);

	gravi_data_append_header (disp_data, gravi_data_get_header (vis_data));


//	cpl_table_delete (p2vm_met);
	for (pol = 0; pol < npol; pol ++){
		vistable_sc = gravi_data_get_oi_table (vis_data,
				GRAVI_OI_VIS_SC_EXT, pol_SC_check[pol]);
		wavelenght_sc = gravi_data_get_oi_table(vis_data,
				GRAVI_OI_WAVELENGTH_SC_EXT, pol_SC_check[pol]);
		wavesc_plist = gravi_data_get_oi_plist(vis_data,
				GRAVI_OI_WAVELENGTH_SC_EXT, pol_SC_check[pol]);
		nwave_sc = gravi_pfits_get_nwave (wavesc_plist);
		wave_sc = cpl_vector_new(nwave_sc);

		plist = gravi_data_get_header (vis_data);

	// not used ?	exptime_sc = gravi_pfits_get_dit (plist, GRAVI_SC)*pow(10, 6);
	//	cpl_matrix * matrix = cpl_matrix_new (nbase, ncol);
	//	cpl_matrix_fill (matrix, 0.0);
	//	for (i = 0; i < nbase; i ++){
	//		cpl_matrix_set (matrix, i, tel_1[i], 1);
	//		cpl_matrix_set (matrix, i, tel_2[i], -1);
	//		cpl_matrix_set (matrix, i, 4 + tel_1[i], 1);
	//		cpl_matrix_set (matrix, i, 4 + tel_2[i], -1);
	//	}

		nfile = cpl_table_get_nrow (vistable_sc) / 6;

	//	gsl_vector * bis;
	//	gsl_matrix * A;

		phase_disp = gsl_matrix_calloc (nfile * 6, 8);

		opd_disp = cpl_malloc(nwave_sc * sizeof(gsl_vector *));
		for (wave = 0; wave < nwave_sc; wave ++){

			cpl_vector_set (wave_sc, wave,
					cpl_table_get_float (wavelenght_sc, "EFF_WAVE", wave, &nv));
			opd_disp[wave] = gsl_vector_alloc (nfile * 6);
		}

		/* Get the mean of the metrology phase */
		cpl_array * met_phase = cpl_array_new(4, CPL_TYPE_DOUBLE);
		cpl_array_fill_window (met_phase, 0, 4, 0.0);
		for (row = 0; row < cpl_table_get_nrow(fddl_met); row ++){
			cpl_array_add (met_phase, cpl_table_get_array(fddl_met, "PHASE_FC", row));
		}

		cpl_array_divide_scalar(met_phase, cpl_table_get_nrow(fddl_met));

		for (base = 0; base < 6; base++){

			ph_met1 = cpl_array_get_double (met_phase, tel_1[base], &nv);/////////////////////?????
			ph_met2 = cpl_array_get_double (met_phase, tel_2[base], &nv);/////////////////////?????
			opd_met = (ph_met1 - ph_met2) * 2 * M_PI * LAMBDA_MET;
            
			/* Construction of each mean complex visibilities and the mean
			 * squared visibilities array */


			/* Compute the average phase of FT to evaluate the phase between
			 * the interval [n DIT, (n+1)DIT] */
			for (row = 0; row < nfile; row ++){

	//				printf("time_SC = %d  time_FT = %d\n", time_SC, time_FT);

				phase_sc = cpl_table_get_array (vistable_sc, "VISPHI", row * 6 + base);

				/* Get the visibilities product by the reduced p2vm */

				/* Compute the opd dispertion */

				/* Interpolate the phase ft  */
				phase_sc_base = cpl_array_duplicate (phase_sc);

				if (cpl_error_get_code()){
					cpl_bivector_unwrap_vectors(fout);
					cpl_bivector_unwrap_vectors(fref);
					 cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT,
						"Error during the interpolation of phase ft");
					 return NULL;
				}

				/*Construction of the phase matrix */
				fddl_sc = cpl_table_get_array (fddl_met, "FT_POS", row);
				fddl_ft = cpl_table_get_array (fddl_met, "SC_POS", row);
				pos_sc = cpl_array_get_float (fddl_sc, tel_1[base], &nv);
				pos_ft = cpl_array_get_float (fddl_ft, tel_1[base], &nv);

				pos_mean = (pos_sc + pos_ft) / 2;
				gsl_matrix_set (phase_disp, row * nbase + base, tel_1[base], ph_met1);
				gsl_matrix_set (phase_disp, row * nbase + base, tel_2[base], ph_met2);
				gsl_matrix_set (phase_disp, row * nbase + base, 4 + tel_1[base], pos_mean);

				pos_sc = cpl_array_get_float (fddl_sc, tel_2[base], &nv);
				pos_ft = cpl_array_get_float (fddl_ft, tel_2[base], &nv);

				pos_mean = (pos_sc + pos_ft) / 2;

				gsl_matrix_set (phase_disp, row * nbase + base, 4 + tel_2[base], pos_mean);

				/* Construction  of the opd vector  */
				for (wave = 0; wave < nwave_sc; wave++){
					ph_sc = cpl_array_get_double (phase_sc_base, wave, &nv);
					opd_sc = ph_sc * cpl_vector_get (wave_sc, wave)/ (2 * M_PI);
	//				printf("cpl_vector_get (wave_sc, wave) = %e\n", cpl_vector_get (wave_sc, wave));
	//				printf("wave = %d\n", wave);
	//				printf("ph_sc = %e opd_sc \n", ph_sc,opd_sc);
					gsl_vector_set (opd_disp[wave], base + nbase * row,
													opd_sc);



				}
				if (cpl_error_get_code()){
					cpl_bivector_unwrap_vectors(fout);
					cpl_bivector_unwrap_vectors(fref);
					 cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT,
						"Error during the construction of the opd vector");
					 return NULL;
				}

				cpl_array_delete (phase_sc_base);
			}


		}
		cpl_vector_unwrap (wave_sc);
		cpl_array_delete (met_phase);
		cpl_table * disp_table = cpl_table_new (nwave_sc);
		cpl_table_new_column_array (disp_table, "COEFF",
				CPL_TYPE_DOUBLE, 8);
		cpl_array * coeff;

		gsl_vector * S;
		gsl_vector * work;
		gsl_vector * x;
		gsl_matrix * u = gsl_matrix_calloc (nfile * 6, 8), * V;
		int i;

		for (wave = 0; wave < nwave_sc; wave ++){

			gsl_matrix_memcpy(u, phase_disp);
			x = gsl_vector_alloc (8);
			work = gsl_vector_alloc (8);
			S = gsl_vector_alloc (8);
			coeff = cpl_array_new (8, CPL_TYPE_DOUBLE);
			V = gsl_matrix_alloc (8, 8);

			gsl_linalg_SV_decomp (u, V, S, work);
			gsl_linalg_SV_solve (u, V, S, opd_disp[wave], x);

			for (i = 0; i < 8; i ++) {
				cpl_array_set_double (coeff, i, gsl_vector_get (x, i));
			}

			cpl_table_set_array (disp_table, "COEFF", wave, coeff);

			cpl_array_delete (coeff);
			gsl_vector_free (opd_disp[wave]);
			gsl_matrix_free (V);
			gsl_vector_free (x);
			gsl_vector_free (work);
			gsl_vector_free (S);
		}
		gsl_matrix_free (u);
		gsl_matrix_free (phase_disp);
		cpl_free (opd_disp);


		cpl_propertylist * met_plist = cpl_propertylist_duplicate (wavesc_plist);
		char * name_tb = cpl_sprintf("DISP_%s", pol_SC_check[pol]);
		cpl_propertylist_append_string (met_plist, "EXTNAME",
				name_tb);
		cpl_free (name_tb);
		gravi_data_add (disp_data, met_plist, disp_table);
		cpl_propertylist_delete (met_plist);
	} // End loop on polarisation


	gravi_msg_function_exit(1);
	return disp_data;
}



/*----------------------------------------------------------------------------*/
/**
  @brief   This function calibrate the wavelength using the argon

  @param argon_data  	The input argon data
  @param wave_data		The calibrated wavelength using the gravi_compute_wave function
  @param dark_map		the dark map computed using the function gravi_compute_dark
  @param bad_map		the bad pixwel map computed using the function
  	  	  	  	  	  	gravi_compute_badpix

  @return   the output gravi data containing the wavelength calibration

  The wavelength is compute by fitting the wavelength argon on the wavelength
  already calibrated using the comparison of the measured phase of each spectral
  element with the realized OPD.
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_fit_argon (gravi_data * argon_data, gravi_data * wave_data,
							  gravi_data * profile_map, gravi_data * dark_map,
							  gravi_data * bad_map) {
    gravi_msg_function_start(1);

	cpl_propertylist * primary_hdr, * plist;
	cpl_table * spectrum_table, * imaging_detector, * wave_data_sc,
			  * img_output;
	cpl_matrix * all_coord;
	cpl_array * wavelength;
	cpl_polynomial  * fit2d;
	cpl_vector  * residuals;
	double rechisq;
	// Degrees of wavelength fit polynomial {2 in vertical direction, 4 in disp direction}
	const cpl_size  deg2d[2] = {2, 2};
	const cpl_size  deg1d[1] = {2};
	cpl_vector * coord_X, * coord_Y, * all_wavelength;
	double result;
	int wave;
	double slope;
	int y_corner;
	int nwave, n_region, region;
	char * data_x;
	int nv;
	cpl_bivector * plot;
	cpl_vector *pos;
	double minwave = 0, maxwave = 1e10;
	cpl_image * img_profile, * image_wave, * profile_image;
	cpl_imagelist * imglist_wave;
	cpl_table * profile_table;
	cpl_table * wave_fibre;
	int sizey, sizex, ind;

	/* Check the inputs data */
	if (argon_data == NULL){
        cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT, "The data is NULL");
        return NULL;
    }

    clock_t start, end;
    start = clock();
    cpl_imagelist * imglist = gravi_data_get_cube (argon_data, GRAVI_IMAGING_DATA_SC_EXT);
    cpl_image * img_median = cpl_imagelist_collapse_median_create (imglist);

    cpl_imagelist * imglist_med = cpl_imagelist_new ();
    cpl_imagelist_set (imglist_med, img_median, 0);

    /* Put the data in the output table : dark_map */
    gravi_data_set_cube (argon_data, GRAVI_IMAGING_DATA_SC_EXT, imglist_med);

    gravi_data * output_data = gravi_extract_spectrum(argon_data, profile_map,
    		dark_map, bad_map, NULL);

    end = clock();
    cpl_msg_info (cpl_func, "Execution time gravi_extract_spectrum : %f", (end - start) / (double)CLOCKS_PER_SEC);

	CPLCHECK_NUL("Connot extract the spectrum");

	/* Wavelengths of the argon emission lines (microns) */
	int nlines = 10;
	double  line_wave[] = {/*1.982291,*/
                   1.997118e-6,
                   2.032256e-6,
                   2.062186e-6,
                   /*2.065277e-6,*/
                   /*2.073922e-6,*/
                   /*2.081672e-6,*/
                   2.099184e-6,
                   2.133871e-6,
                   2.154009e-6,
                   2.208321e-6,
                   2.313952e-6,
                   2.385154e-6,
                   2.397306e-6};

	slope = 0.02;
	/* Correction of the wavelength scale in glass to go to vaccuum scale */
	imaging_detector = gravi_data_get_table(output_data, GRAVI_IMAGING_DETECTOR_SC_EXT);
	spectrum_table = gravi_data_get_table (output_data, GRAVI_SPECTRUM_DATA_SC_EXT);

	if ((spectrum_table == NULL) || (imaging_detector == NULL)) {
		gravi_data_delete(output_data);
		cpl_error_set_message(cpl_func,
							  CPL_ERROR_ILLEGAL_OUTPUT, "Data must contain SPECTRUM_DATA");
		return NULL;
	}

	wave_data_sc =  cpl_table_duplicate (gravi_data_get_table (wave_data, GRAVI_WAVE_DATA_SC_EXT));
	double wavestep = cpl_array_get_double (cpl_table_get_array (wave_data_sc, "DATA1", 0), 1, &nv) -
	  cpl_array_get_double (cpl_table_get_array (wave_data_sc, "DATA1", 0), 0, &nv);

	cpl_array ** wave_array;
	double * wave_ ;
//	for (region = 0; region < n_region; region ++){
//
//		data_x = cpl_sprintf("DATA%d", region + 1);
//		wave_array = cpl_table_get_data_array (wave_data_sc, data_x);
//		wave_ = cpl_array_get_data_double (wave_array[0]);
//		for (wave = 0; wave < nwave; wave ++){
//			wave_ [wave] -= wave * wavestep * slope;
//		}
//
//		cpl_free (data_x);
//	}


	/* Get the number of wavelength and the region */
	nwave = cpl_table_get_column_dimension(spectrum_table, "DATA1", 1);
	n_region = cpl_table_get_nrow(imaging_detector);

	/* Loop region */
	coord_X = cpl_vector_new(n_region * nlines);
	coord_Y = cpl_vector_new(n_region * nlines);
	cpl_vector * coord_X_fit = cpl_vector_new(n_region * nlines);
	all_wavelength = cpl_vector_new(n_region * nlines);
	cpl_vector * fitsigm = cpl_vector_new(n_region * nlines);

	cpl_vector_fill (coord_X, -1);
	cpl_vector_fill (coord_Y, -1);
	cpl_vector_fill (all_wavelength, -1);

	int fitwidth, fit_in, comp = 0, list;

	primary_hdr = gravi_data_get_header (argon_data);
	const char * resolution = gravi_pfits_get_resolution (primary_hdr);

	if (! (strcmp(resolution, "LOW") && strcmp(resolution, "MED")) )
		fitwidth = 3;
	else
		fitwidth = 10;

	/* Fit of each emission line position using a quadratic model of their
	 * position as a function of channel number */
	const cpl_array * argon;

	for (region = 0; region < n_region; region ++) {

		data_x = cpl_sprintf("DATA%d", region + 1);
		wave_array = cpl_table_get_data_array (wave_data_sc, data_x);
		wave_ = cpl_array_get_data_double (wave_array[0]);

		for (wave = 0; wave < nwave; wave ++){
			wave_ [wave] -= wave * wavestep * slope;
		}

		argon = cpl_table_get_array (spectrum_table, data_x, 0);

		cpl_free (data_x);

		wave = 0;
		for (list = 0; list < nlines; list ++) {

			while (wave_[wave] < line_wave[list]){
				wave ++;
			}

//			wave += 20;

			if (wave >= nwave){
				cpl_error_set_message(cpl_func,
						CPL_ERROR_ILLEGAL_OUTPUT, "The argon wavelength does "
								"not much with the calibration wavelength ");
				return NULL;
			}

			double flux = 0;
			double coord_flux = 0, coord_err = 0;
			cpl_vector * vector_x=cpl_vector_new(fitwidth*2);
			cpl_vector * vector_y=cpl_vector_new(fitwidth*2);
			int i=0;
			for (fit_in = wave - fitwidth; fit_in < wave + fitwidth; fit_in ++){
				cpl_vector_set(vector_x, i, fit_in);
				cpl_vector_set(vector_y, i, cpl_array_get_double (argon, fit_in, &nv));
				i++;
				flux += cpl_array_get_double (argon, fit_in, &nv);
				coord_flux += fit_in * cpl_array_get_double (argon, fit_in, &nv);
//				printf("region %d; line %d; %d, flux= %g \n", region, list, fit_in, cpl_array_get_double (argon, fit_in, &nv));
			}

			cpl_errorstate prestate = cpl_errorstate_get();
			double x0, sigma, area, offset, mse;
			cpl_vector_fit_gaussian (vector_x, NULL, vector_y, NULL,
								CPL_FIT_ALL, &x0, &sigma, &area,
									 &offset, &mse, NULL, NULL);

			coord_flux =  (coord_flux / flux);
			coord_err = fabs (flux);
//			printf("Barycentre = %g Gaussian FIT :%g (+-%g) \n", coord_flux, x0, sigma);

			if (cpl_error_get_code() == CPL_ERROR_CONTINUE){
		    	cpl_errorstate_set (prestate);
		    	cpl_msg_warning(cpl_func, "The gaussian fit did not converge while fitting argon line");
		    	x0=coord_flux;
		    	sigma=100;

		    }

			CPLCHECK_NUL("Error during the computation");

//			cpl_vector_set (coord_X, list * n_region + region, coord_flux);
			cpl_vector_set (coord_X, list * n_region + region, x0);
			cpl_vector_set (coord_Y, list * n_region + region, region + 1);
//			cpl_vector_set (fitsigm, list * n_region + region, coord_err);
			cpl_vector_set (fitsigm, list * n_region + region, sigma);
			cpl_vector_set (all_wavelength, list * n_region + region, line_wave[list]);

			comp++;
		} /* End loop on list of lines */

	} /* End loop on regions */

	cpl_table_delete (wave_data_sc);
	int i;

	/* Fit of the barycenters of each argon wavelength of all the regions */
	cpl_vector * coord_x_list, * fitsigm_list,
	* coord_y_list, * coord_x_fit = cpl_vector_new (n_region);

	cpl_matrix * matrix3, * matrix2;

	gsl_matrix * X = gsl_matrix_alloc (n_region, 3);
	gsl_matrix * X_bis = gsl_matrix_alloc (n_region, 3);
	gsl_vector * y = gsl_vector_alloc (n_region);
	gsl_vector * w = gsl_vector_alloc (n_region);

	gsl_vector * c = gsl_vector_alloc (3);
	gsl_matrix * cov = gsl_matrix_alloc (3, 3);
	double chisq, result_gsl;

	/* Loop on the list of lines */
	for (list = 0; list < nlines; list++) {

		coord_x_list = cpl_vector_extract (coord_X, list * n_region, (list + 1) * n_region - 1, 1);
		fitsigm_list = cpl_vector_extract (fitsigm, list * n_region, (list + 1) * n_region - 1, 1);
		coord_y_list = cpl_vector_extract (coord_Y, list * n_region, (list + 1) * n_region - 1, 1);

		for (i = 0; i < n_region; i++){
			gsl_matrix_set (X, i, 0, 1.0);
			gsl_matrix_set (X, i, 1, cpl_vector_get (coord_y_list, i));
			gsl_matrix_set (X, i, 2, cpl_vector_get (coord_y_list, i) *
							cpl_vector_get (coord_y_list, i));

			gsl_vector_set (y, i, cpl_vector_get (coord_x_list, i));
			gsl_vector_set (w, i, cpl_vector_get (fitsigm_list, i));
		}

		gsl_matrix_memcpy (X_bis, X);
	    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n_region, 3);
	    gsl_multifit_wlinear (X_bis, w, y, c, cov,
	                          &chisq, work);
	    gsl_multifit_linear_free (work);

		CPLCHECK_NUL("Cannot compute the barycneter of all regions");

		for (region = 0; region < n_region; region ++){

			result_gsl = gsl_vector_get (c, 0) + gsl_vector_get (c, 1) *
					gsl_matrix_get (X, region, 1) + gsl_vector_get (c, 2) *
						gsl_matrix_get (X, region, 2);
			cpl_vector_set (coord_X_fit, list * n_region + region, result_gsl);
			cpl_vector_set (coord_x_fit, region, result_gsl);
		}

		cpl_vector_delete (coord_x_list);
		cpl_vector_delete (fitsigm_list);
		cpl_vector_delete (coord_y_list);

	} /* End loop on the list of lines */

	gsl_matrix_free (cov);
	gsl_matrix_free (X);
	gsl_matrix_free (X_bis);
	gsl_vector_free (y);
	gsl_vector_free (w);
	gsl_vector_free (c);

	cpl_vector_delete (coord_x_fit);

	/* Compute the polynomial model between the position and the
	 * wavelength = F(x, y) where F is polynomial function */

	fit2d = cpl_polynomial_new(2);

	matrix3 = cpl_matrix_wrap(1, nlines * n_region, cpl_vector_get_data(coord_X));
	matrix2 = cpl_matrix_wrap(1, nlines * n_region, cpl_vector_get_data(coord_Y));

	all_coord = cpl_matrix_duplicate (matrix3) ;
	cpl_matrix_append(all_coord, matrix2, 1);

	/* fit CPL */
	cpl_polynomial_fit(fit2d, all_coord, NULL, all_wavelength, NULL,
											 CPL_TRUE, NULL, deg2d);
	residuals = cpl_vector_new(nlines*n_region);
	cpl_vector_fill_polynomial_fit_residual	(residuals, all_wavelength, NULL,
			fit2d, all_coord, &rechisq );

	/* fits gsl */
	X = gsl_matrix_alloc (n_region*nlines, 6);
	X_bis = gsl_matrix_alloc (n_region*nlines, 6);
	y = gsl_vector_alloc (n_region*nlines);
	w = gsl_vector_alloc (n_region*nlines);

	c = gsl_vector_alloc (6);
	cov = gsl_matrix_alloc (6, 6);

	for (list = 0; list < nlines; list++) {

		coord_x_list = cpl_vector_extract (coord_X, list * n_region, (list + 1) * n_region - 1, 1);
		fitsigm_list = cpl_vector_extract (fitsigm, list * n_region, (list + 1) * n_region - 1, 1);
		coord_y_list = cpl_vector_extract (coord_Y, list * n_region, (list + 1) * n_region - 1, 1);

		for (i = 0; i < n_region; i++){
			gsl_matrix_set (X, i+list*n_region, 0, 1.0);
			gsl_matrix_set (X, i+list*n_region, 1, cpl_vector_get (coord_x_list, i));
			gsl_matrix_set (X, i+list*n_region, 2, cpl_vector_get (coord_x_list, i) *
							cpl_vector_get (coord_x_list, i));
			gsl_matrix_set (X, i+list*n_region, 3, cpl_vector_get (coord_y_list, i));
			gsl_matrix_set (X, i+list*n_region, 4, cpl_vector_get (coord_y_list, i) *
							cpl_vector_get (coord_y_list, i));
			gsl_matrix_set (X, i+list*n_region, 5, cpl_vector_get (coord_x_list, i) *
							cpl_vector_get (coord_y_list, i));

			gsl_vector_set (y, i+list*n_region, line_wave[list]);
			gsl_vector_set (w, i+list*n_region, 1./pow(cpl_vector_get (fitsigm_list, i),2));
		}
	}

	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n_region*nlines, 6);
	gsl_matrix_memcpy (X_bis, X);
	gsl_multifit_wlinear (X_bis, w, y, c, cov,
                          &chisq, work);
	gsl_vector * r=gsl_vector_alloc (n_region*nlines);
	gsl_multifit_linear_residuals (X_bis, y, c, r);
    gsl_multifit_linear_free (work);

	for (list = 0; list < nlines; list++) {
		for (i = 0; i < n_region; i++){
//			printf("line %d region %d residual %g\n", list, i, gsl_vector_get(r, i+list*n_region));
		}
	}



	CPLCHECK_NUL("Error in the fit of the wavelength and position");


//	if (PLOT_WAVELENGTH)
//	{
//		vectors=malloc(3 * sizeof(cpl_vector*));
//		vectors[0]=NULL;
//		vectors[1]=all_wavelength;
//		cpl_vector *XYpos=cpl_vector_duplicate(coord_X);
//		cpl_vector_multiply_scalar(XYpos, nlines/5);
//		cpl_vector_add(XYpos, coord_Y);
//		cpl_vector_divide_scalar(XYpos, nlines/5);
//		vectors[0]=XYpos;
//		cpl_bivector *error=cpl_bivector_wrap_vectors( XYpos, cpl_vector_duplicate(residuals));
//		cpl_vector_subtract(residuals, all_wavelength);
//		cpl_vector_multiply_scalar(residuals, -1000000);
//		cpl_vector_multiply_scalar(all_wavelength, 1000000);
//		vectors[2]=residuals;
////					vectors[2]=all_wavelength;
//		if (POSTSCRIPT_PLOT) ps_string=cpl_sprintf("set term png; set output 'plot_wavelength_%s.png';", "SC_Argon");
//		else  ps_string=cpl_sprintf(" ");
//		cpl_plot_vectors(cpl_sprintf("%s set xlabel \"Index\";set ylabel \"Wavelength\";", ps_string)," title 'mesured';", NULL, vectors, 3);
//		cpl_plot_bivector(cpl_sprintf("%s set xlabel \"Index\";set ylabel \"Wavelength\";", ps_string)," title 'error';", NULL, error);
//		free(vectors);
//		cpl_bivector_delete(error);
//		cpl_free(ps_string);
//	}

	cpl_vector_delete(coord_X_fit);
	cpl_vector_delete(all_wavelength);
	cpl_vector_delete(coord_X);
	cpl_vector_delete(coord_Y);
	cpl_matrix_delete (all_coord);
	cpl_matrix_unwrap(matrix2);
	cpl_matrix_unwrap(matrix3);
	cpl_vector_delete(residuals);
	cpl_vector_delete(fitsigm);

	/* Calcul of the new wavelength vector after calibration */

	/* Save the result wavelength on the associated bases depending
	 * of the polarization */
	cpl_array * dimension = cpl_array_new(2, CPL_TYPE_INT);
	cpl_array_set(dimension, 0, 1);
	cpl_array_set(dimension, 1, nwave);

	profile_table = gravi_data_get_table (profile_map, GRAVI_PROFILE_DATA_EXT);
	sizex = cpl_table_get_column_dimension (profile_table, "DATA1", 0);
	sizey = cpl_table_get_column_dimension (profile_table, "DATA1", 1);
	profile_image = cpl_image_new (sizex, sizey, CPL_TYPE_DOUBLE);
	cpl_image_fill_window (profile_image, 1, 1, sizex, sizey, 0.0);
	image_wave = cpl_image_new (sizex, sizey, CPL_TYPE_DOUBLE);
	cpl_image_fill_window (image_wave, 1, 1, sizex, sizey, 0.0);

	img_output = cpl_table_new (1);

	for (region = 0 ; region < n_region; region ++) {

		data_x = cpl_sprintf("DATA%d", region + 1);

		imglist_wave = gravi_table_data_to_imagelist(profile_table, region + 1);
		img_profile = cpl_imagelist_get(imglist_wave, 0);

		cpl_table_new_column_array (img_output, data_x, CPL_TYPE_DOUBLE, nwave);
		cpl_table_set_column_dimensions(img_output, data_x, dimension);

		wavelength = cpl_array_new(nwave, CPL_TYPE_DOUBLE);

		/* Calculate of the new wavelength vector after calibration */
		y_corner = 1;

		cpl_image_add (profile_image, img_profile);
		pos = cpl_vector_new(2);
		for (wave = 0; wave < nwave; wave ++) {

			cpl_vector_set(pos, 0, wave);
			cpl_vector_set(pos, 1, region + 1);
			result=gsl_vector_get (c, 0) +
					gsl_vector_get (c, 1) *	wave +
					gsl_vector_get (c, 2) *	wave*wave +
					gsl_vector_get (c, 3) *	(region+1) +
					gsl_vector_get (c, 4) *	(region+1)*(region+1) +
					gsl_vector_get (c, 5) *	(region+1)*wave;

			//result = cpl_polynomial_eval(fit2d, pos);

			cpl_array_set(wavelength, wave, result);
			y_corner ++;

			for (ind = 0; ind < sizey; ind ++){

			  if (cpl_image_get (img_profile, wave+1, ind+1, &nv) > 0.01)
				cpl_image_set (image_wave, wave+1, ind+1, result);

			  CPLCHECK_NUL("The corner and image_wave");
			}
		}
		cpl_vector_delete(pos);

		/* Get the miniumum and maximum wavelength */
		minwave = CPL_MAX (minwave, cpl_array_get_min(wavelength));
		maxwave = CPL_MIN (maxwave, cpl_array_get_max(wavelength));

		cpl_table_set_array(img_output, data_x, 0, wavelength);

		cpl_imagelist_delete (imglist_wave);
		cpl_array_delete(wavelength);
		cpl_free(data_x);

	} /* End loop on regions */


	gravi_data * wave_map;

	wave_map = gravi_data_new (0);
	plist = gravi_data_get_header (wave_data);
	gravi_data_append_header (wave_map, plist);
	int nb_ext =  gravi_data_get_size (wave_data), ext, type_ext;
	cpl_propertylist * img_plist;

	for (ext = 0; ext < nb_ext; ext++ ){

    	/*
    	 * Load the FT or SC
    	 */
    	img_plist = gravi_data_get_plist_x (wave_data, ext);

		const char * plist_name = gravi_pfits_get_extname (img_plist);
		/* Check if the needed extentions are there */
    	type_ext = gravi_pfits_get_extension_type (img_plist);
		if (!(strcmp (plist_name, "WAVE_FIBRE_SC") &&
				    strcmp (plist_name, GRAVI_IMAGING_DETECTOR_SC_EXT) &&
				      strcmp (plist_name, GRAVI_WAVE_DATA_SC_EXT) &&
				        strcmp (plist_name, "TEST_WAVE"))){
			/*	Load the FT table	*/
			if (type_ext == 2){
				gravi_data_add (wave_map, img_plist,
							cpl_table_duplicate (gravi_data_get_table (wave_data, plist_name)));
				if(! strcmp (plist_name, GRAVI_WAVE_DATA_SC_EXT))
					cpl_propertylist_set_string (gravi_data_get_plist (wave_map,
						GRAVI_WAVE_DATA_SC_EXT),  "EXTNAME", GRAVI_WAVE_ARGON_EXT);
			}

			/*	Load the SC image_list	*/
			else if (type_ext == 3)
				gravi_data_add_cube (wave_map, img_plist,
							cpl_imagelist_duplicate (gravi_data_get_cube (wave_data, plist_name)));
		}
    }

	/* Set the QC parameters */
	primary_hdr = gravi_data_get_header (wave_map);
	cpl_propertylist_append_double (primary_hdr, QC_MINWAVE_SC,	minwave);
	cpl_propertylist_append_double (primary_hdr, QC_MAXWAVE_SC, maxwave);
	cpl_msg_info (cpl_func, "QC_MINWAVE_SC = %e QC_MAXWAVE_SC = %e", minwave, maxwave);

	/* Add the img_output and wave_fibre */
	gravi_data_set_table (wave_map, GRAVI_WAVE_ARGON_EXT, img_output);

	/* Add the image_wave and profile_image */
	imglist_wave = cpl_imagelist_new ();
	cpl_imagelist_set(imglist_wave, image_wave, 0);
	cpl_imagelist_set(imglist_wave, profile_image, 1);

	gravi_data_set_cube (wave_map, "TEST_WAVE", imglist_wave);

	cpl_polynomial_delete (fit2d);
	cpl_array_delete (dimension);
	gravi_data_delete (output_data);

    gravi_msg_function_exit(1);
	return wave_map;
}

gravi_data * gravi_compute_wave_offset (gravi_data * argon_wave, gravi_data * wave_data) {

	cpl_propertylist * primary_hdr, * plist, * wavePlist, * detectorPlist;
	cpl_table * detector, * waveArgon_table, * waveData_table,
			  * img_output;
	cpl_matrix * all_coord;
	cpl_array * wavelength;
	cpl_polynomial  * fit_slope;
	cpl_polynomial  * fit2d;
	cpl_vector  * residuals;
	double rechisq;
	// Degrees of wavelength fit polynomial {2 in vertical direction, 4 in disp direction}
	const cpl_size  deg2d[2] = {4, 4};
	cpl_vector * coord_X, * coord_Y, * all_wavelength;
	double result;
	int wave;
	double slope;
	int y_corner;
	int nwave, n_region, region;
	char * data_x;
	int nv;
//	const cpl_vector ** vectors;
	cpl_bivector * plot;
	cpl_vector *pos;
	double minwave, maxwave;
	cpl_image * image_wave;
	cpl_imagelist * imglist_wave;
	int sizey, sizex, ind, npol;
	/* Tel for base */
	//int tel1[6] = {0,0,0,1,1,2};
	//int tel2[6] = {1,2,3,2,3,3};
	char * baseString [6] = {"12", "13", "14", "23", "24", "34"};
	char * baseString_bis [6] = {"21", "31", "41", "32", "42", "43"};
	const cpl_size   deg = 0, degM = 1;
	cpl_vector * position;
	cpl_matrix * matFit2;


	if ((argon_wave == NULL) || (wave_data == NULL)){
        cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT,
        		                                         "The data is NULL");
        return NULL;
    }

	int nbase = 6, base, type_data;
	primary_hdr = gravi_data_get_header (argon_wave);
	const char * resolArg = gravi_pfits_get_resolution (primary_hdr);
	const char * resolWave = gravi_pfits_get_resolution (gravi_data_get_header (wave_data));

	if (strcmp(resolArg, resolWave)){
		cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
								  "The argon file doesn t have the same"
		"resultion as the input wave");
		return NULL;
	}

	/* Compare the wavelength for each base */

	type_data = 0;

	detector = gravi_data_get_table (argon_wave, GRAVI_IMAGING_DETECTOR_SC_EXT);
	detectorPlist = gravi_data_get_plist (wave_data, GRAVI_IMAGING_DETECTOR_SC_EXT);
	wavePlist = gravi_data_get_plist (wave_data, GRAVI_WAVE_DATA_SC_EXT);

	n_region  = cpl_table_get_nrow (detector) ;

	npol = (n_region > 24)?2:1;

	char * pola = (npol == 2)?POLAR_1:POLAR_3;

	waveArgon_table = gravi_data_get_table (argon_wave, GRAVI_WAVE_ARGON_EXT);

	waveData_table = gravi_data_get_table (wave_data, GRAVI_WAVE_FIBRE_SC_EXT);

	nwave = cpl_table_get_column_depth(waveArgon_table, "DATA1");
	char * base_fiber = cpl_sprintf("BASE_12_%s", pola);
	if (nwave!= cpl_table_get_column_depth(waveData_table, base_fiber)){
		cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
								  "The argon and the wave data do not have"
								  " the same number of wavelength");
		return NULL;
	}
	cpl_array * wave_;
	int comp;

	cpl_imagelist * imglist = gravi_data_get_cube (wave_data, "TEST_WAVE");
	cpl_image * img_profile = cpl_imagelist_get (imglist, 1);
	sizex = cpl_image_get_size_x (img_profile);
	sizey = cpl_image_get_size_y (img_profile);
	image_wave = cpl_image_new (sizex, sizey, CPL_TYPE_DOUBLE);
	cpl_image * image_real = cpl_image_new (sizex, sizey, CPL_TYPE_DOUBLE);
	cpl_image_fill_window (image_real, 1, 1, sizex, sizey, 0.0);
	matFit2 = cpl_matrix_new (2,n_region*nwave);
	all_wavelength = cpl_vector_new (n_region*nwave);
	int regionA, regionD, regionC, regionB;
	cpl_array * waveToCal;

	cpl_table * cooef_table = cpl_table_new (n_region);
	cpl_table_new_column(cooef_table, "SLOPE", CPL_TYPE_DOUBLE);
	cpl_table_new_column(cooef_table, "OFFSET", CPL_TYPE_DOUBLE);

	int startx = gravi_pfits_get_startx(wavePlist) +
			gravi_pfits_get_window_start(primary_hdr) -2,  i;
	if (cpl_error_get_code ()){
		cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT, "Problem to get "
				"the start window detector");
		return NULL;
	}
	int np = 0, nq = 1;
	cpl_vector * vect_index = cpl_vector_new (nwave), * vect;

	for (i = startx; i < startx + nwave; i++){
		if ((i/64)%2 == 0){
			np ++;
			cpl_vector_set (vect_index, i - startx, 0);
		}
		else {
			nq ++;
			cpl_vector_set (vect_index, i - startx, 1);
		}
	}

	cpl_vector * residual_p = cpl_vector_new(np * n_region),
			* residual_q = cpl_vector_new(nq * n_region);

	for (region = 0; region < n_region; region ++) {
		/* Look for the region with the same base none looking for
		 * the polarization */
//		comp = 0;
//
//		wavelength = cpl_array_new (nwave, CPL_TYPE_DOUBLE);
//		waveToCal = cpl_array_new (nwave, CPL_TYPE_DOUBLE);
//		cpl_array_fill_window (wavelength, 0, nwave, 0.0);
//		cpl_array_fill_window (waveToCal, 0, nwave, 0.0);
//
//		regionA = gravi_table_get_regindex(detector, tel1[base] + 1,
//				tel2[base] + 1, "A", pola);
//
//		if (regionA != 0) {
//			data_x = cpl_sprintf ("DATA%d", regionA);
//
//
//			wave_ = cpl_array_duplicate (cpl_table_get_array (waveData_table,
//					data_x, 0));
//			cpl_array_add (waveToCal, wave_);
//			cpl_array_subtract (wave_, cpl_table_get_array (waveArgon_table,
//					data_x, 0));
//			cpl_free (data_x);
//
//			cpl_array_add (wavelength, wave_);
//			cpl_array_delete (wave_);
//			comp ++;
//		}
//
//		regionB = gravi_table_get_regindex(detector, tel1[base] + 1,
//							tel2[base] + 1, "B", pola);
//
//		if (regionB != 0) {
//			data_x = cpl_sprintf ("DATA%d", regionB);
//
//
//			wave_ = cpl_array_duplicate (cpl_table_get_array (waveData_table,
//					data_x, 0));
//			cpl_array_add (waveToCal, wave_);
//			cpl_array_subtract (wave_, cpl_table_get_array (waveArgon_table,
//					data_x, 0));
//			cpl_free (data_x);
//
//			cpl_array_add(wavelength, wave_);
//			cpl_array_delete (wave_);
//			comp ++;
//		}
//
//		regionC = gravi_table_get_regindex(detector, tel1[base] + 1,
//							tel2[base] + 1, "C", pola);
//		if (regionC != 0) {
//			data_x = cpl_sprintf ("DATA%d", regionC);
//
//
//			wave_ = cpl_array_duplicate (cpl_table_get_array (waveData_table,
//					data_x, 0));
//			cpl_array_add (waveToCal, wave_);
//			cpl_array_subtract (wave_, cpl_table_get_array (waveArgon_table,
//					data_x, 0));
//			cpl_free (data_x);
//
//			cpl_array_add(wavelength, wave_);
//			cpl_array_delete (wave_);
//			comp ++;
//		}
//
//		regionD = gravi_table_get_regindex(detector, tel1[base] + 1,
//							tel2[base] + 1, "D", pola);
//		if (regionD != 0) {
//			data_x = cpl_sprintf ("DATA%d", regionD);
//
//
//			wave_ = cpl_array_duplicate (cpl_table_get_array (waveData_table,
//					data_x, 0));
//			cpl_array_add (waveToCal, wave_);
//			cpl_array_subtract (wave_, cpl_table_get_array (waveArgon_table,
//					data_x, 0));
//			cpl_free (data_x);
//
//			cpl_array_add(wavelength, wave_);
//			cpl_array_delete (wave_);
//			comp ++;
//		}
//
//		if (comp != 0){
//			cpl_array_divide_scalar (wavelength, comp);
//			cpl_array_divide_scalar (waveToCal, comp);
//		}
//
//		else {
//			cpl_msg_info (cpl_func, "The base %d %d does "
//					"not exist", tel1[base]+1, tel2[base]+1);
//			cpl_array_delete(wavelength);
//			cpl_array_delete(waveToCal);
//			continue;
//		}

		const char * regname = cpl_table_get_string (detector, "REGNAME", region);
		for (base = 0; base < nbase; base++){
//			printf("baseString[%d] = %s  baseString_bis[%d] = %s\n",base,
//					baseString[base],base,  baseString_bis[base]);
			if (strstr(regname, baseString[base]) ||
					strstr(regname, baseString_bis[base]))
				break;
		}

		if (npol == 2) {
		  if (strstr (regname, POLAR_1))
			pola = POLAR_1;
		  else
			pola = POLAR_2;
		}

		base_fiber = cpl_sprintf ("BASE_%s_%s", baseString[base], pola);
		data_x = cpl_sprintf ("DATA%d", region + 1);
		wavelength = cpl_array_duplicate (cpl_table_get_array (waveData_table,
				base_fiber, 0));
		waveToCal = cpl_array_duplicate (wavelength);
		cpl_array_subtract (wavelength, cpl_table_get_array (waveArgon_table,
				data_x, 0));

		/* Fit to get the slope of the polynomial wavelength */

		cpl_matrix * matrix = cpl_matrix_new (1, nwave);
		position = cpl_vector_wrap (nwave, cpl_array_get_data_double(wavelength));

		for (ind = 0; ind < nwave; ind++) {
			cpl_matrix_set (matFit2,0, region*nwave + ind,
					region);
			cpl_matrix_set (matFit2,1, region*nwave + ind,
					ind);
			cpl_matrix_set(matrix, 0, ind, ind);

		}

//		// fit with cpl_fit_lvmq
//		int
//		cpl_fit_lvmq(matrix, NULL, position,
//		NULL, init_val, val_to_fit, &polystep,
//		NULL, CPL_FIT_LVMQ_TOLERANCE, CPL_FIT_LVMQ_COUNT,
//		CPL_FIT_LVMQ_MAXITER, &mse, &red_chisq, NULL);


		fit_slope =	cpl_polynomial_new(1);
		cpl_polynomial_fit(fit_slope, matrix, NULL, position, NULL,
										CPL_FALSE, &deg, &degM);

		vect = cpl_vector_new (nwave);
		cpl_vector_fill_polynomial_fit_residual	(vect, position, NULL,
				fit_slope, matrix, &rechisq );
		if (cpl_error_get_code ()){
			cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT, "fit the residual");
			return NULL;
		}

		int ip = 0, iq = 0;

		for (wave = 0; wave < nwave; wave ++)
			if (cpl_vector_get (vect_index, wave) == 0) {
				cpl_vector_set (residual_p, ip + np * region ,
						cpl_vector_get (vect, wave));
				ip ++;
			}
			else {
				cpl_vector_set (residual_q, iq + nq * region,
						cpl_vector_get (vect, wave));
				iq ++;
			}



		const cpl_size  pow_slope = 1;
		const cpl_size  pow_slope2 = 0;
		cpl_msg_info (cpl_func, "region %s : y = %e * x + %e", regname,
				cpl_polynomial_get_coeff(fit_slope,&pow_slope ),
				cpl_polynomial_get_coeff(fit_slope,
						&pow_slope2 ));

		cpl_table_set_double (cooef_table, "SLOPE", region,
				cpl_polynomial_get_coeff(fit_slope,&pow_slope ));
		cpl_table_set_double (cooef_table, "OFFSET", region,
				cpl_polynomial_get_coeff(fit_slope,&pow_slope2 ));
		if (cpl_error_get_code()){

			cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT, "fit slope");
			return NULL;
		}

		pos = cpl_vector_new (1);
		for (wave = 0; wave < nwave; wave++){
			cpl_vector_set (pos, 0, wave);
			cpl_vector_set(all_wavelength, region*nwave + wave,
					cpl_array_get_double (waveToCal, wave, &nv) -
					cpl_polynomial_eval(fit_slope, pos));
		}

		cpl_vector_delete (pos);
		cpl_polynomial_delete (fit_slope);
		cpl_matrix_delete (matrix);
		cpl_vector_unwrap (position);
		cpl_array_delete(wavelength);
		cpl_array_delete(waveToCal);
		cpl_free (data_x);
		cpl_free (base_fiber);
	}

	/* Computation of the step */
	double step = cpl_vector_get_median (residual_p) - cpl_vector_get_median (residual_q);


	cpl_msg_info(cpl_func, "step = %g", step);

	/* Fit  */
	fit2d = cpl_polynomial_new(2);
	cpl_polynomial_fit(fit2d, matFit2, NULL, all_wavelength, NULL,
											 CPL_TRUE, NULL, deg2d);
	cpl_matrix_delete (matFit2);

	if (cpl_error_get_code()){

		cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT, "fit 2D");
		return NULL;
	}

	/* Calcul of the new wavelength vector after calibration */

	/* Save the result wavelength on the associated bases depending
	 * of the polarization */
	cpl_array * dimension = cpl_array_new(2, CPL_TYPE_INT);
	cpl_array_set(dimension, 0, 1);
	cpl_array_set(dimension, 1, nwave);

	cpl_image_fill_window (image_wave, 1, 1, sizex, sizey, 0.0);

	img_output = NULL;
	img_output = cpl_table_new (1);

	for (region = 0 ; region < n_region; region ++){

		data_x = cpl_sprintf("DATA%d", region + 1);

		cpl_table_new_column_array (img_output, data_x, CPL_TYPE_DOUBLE, nwave);
		cpl_table_set_column_dimensions(img_output, data_x,
													dimension);
		const char * regname = cpl_table_get_string (detector, "REGNAME", region);

//		printf("regname = %s\n", regname);
		for (base = 0; base < nbase; base++){
//			printf("baseString[%d] = %s  baseString_bis[%d] = %s\n",base,
//					baseString[base],base,  baseString_bis[base]);
			if (strstr(regname, baseString[base]) ||
					strstr(regname, baseString_bis[base]))
				break;
		}
//		printf("base = %d\n", base);
		if (base == nbase){
			cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT,
					"Check the base");
			return NULL;
		}
//		printf("base = %d\n", base);
		wavelength = cpl_array_new(nwave, CPL_TYPE_DOUBLE);
		cpl_vector *wavelength_vect = cpl_vector_new(nwave);
//		vectors[region + 1] = cpl_vector_new (nwave);
		/* Calculate of the new wavelength vector after calibration */
		y_corner = 1;

		pos = cpl_vector_new(2);
		for (wave = 0; wave < nwave; wave ++){


			cpl_vector_set(pos, 0, region + 1);
			cpl_vector_set(pos, 1, wave);

			result = cpl_polynomial_eval(fit2d, pos);
//			result = cpl_vector_get (all_wavelength, base*nwave + wave) +
			result = cpl_vector_get (all_wavelength, region*nwave + wave) -
					cpl_vector_get (vect_index, wave) * step;
			cpl_array_set(wavelength, wave, result);
			cpl_vector_set (wavelength_vect, wave, result );
			y_corner ++;

			for (ind = 0; ind < sizey; ind ++){

				if (cpl_image_get (img_profile, wave+1, ind+1, &nv) > 0.01){

					cpl_image_set (image_wave, wave+1, ind+1, result);

					cpl_image_set (image_real, wave+1, ind+1,
//							cpl_vector_get (all_wavelength, base*nwave + wave));
							cpl_vector_get (all_wavelength, region*nwave + wave));
				}

				if (cpl_error_get_code()){

					cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT, "The corner and "
							"image_wave");
					return NULL;
				}

			}


		}
		cpl_vector_delete(pos);


		/*
		 * Fit the wave_fibre, and add the OFFSET, SLOPE and STEP
		 */
		cpl_matrix * matrix = cpl_matrix_new (1, nwave);
		const cpl_size deg_2=2;

		for (ind = 0; ind < nwave; ind++) {
			cpl_matrix_set(matrix, 0, ind, ind);
		}

		fit_slope =	cpl_polynomial_new(1);
		cpl_polynomial_fit(fit_slope, matrix, NULL, wavelength_vect, NULL,
										CPL_FALSE, &deg, &deg_2);

		pos=cpl_vector_new(1);
		for (wave = 0; wave < nwave; wave ++){
			cpl_vector_set(pos, 0, wave);
			result = cpl_polynomial_eval(fit_slope, pos)-
					cpl_vector_get (vect_index, wave) * step;
			cpl_array_set(wavelength, wave, result);
		}
		cpl_vector_delete(pos);
		cpl_polynomial_delete(fit_slope);

		/* Get the miniumum and maximum wavelength */
		if (region == 0){
			minwave = cpl_array_get_min(wavelength);
			maxwave = cpl_array_get_max(wavelength);
		}
		else {

			if (minwave < cpl_array_get_min(wavelength))
				minwave = cpl_array_get_min(wavelength);


			if (maxwave > cpl_array_get_max(wavelength))
				maxwave = cpl_array_get_max(wavelength);
		}


		cpl_table_set_array(img_output, data_x, 0, wavelength);

		cpl_array_delete(wavelength);
		cpl_free(data_x);
	}

	cpl_vector_delete (all_wavelength);
	gravi_data * wave_map;
	wave_map = gravi_data_new (0);
	plist = gravi_data_get_header (wave_data);
	gravi_data_append_header (wave_map, plist);
	int nb_ext =  gravi_data_get_size (wave_data), ext, type_ext;
	cpl_propertylist * img_plist;

	for (ext = 0; ext < nb_ext; ext++ ){

    	/*
    	 * Load the FT or SC
    	 */
    	img_plist = gravi_data_get_plist_x (wave_data, ext);

		const char * plist_name = gravi_pfits_get_extname (img_plist);
		/* Check if the needed extentions are there */
    	type_ext = gravi_pfits_get_extension_type (img_plist);
		if (!(strcmp (plist_name, GRAVI_WAVE_FIBRE_SC_EXT) &&
				    strcmp (plist_name, GRAVI_IMAGING_DETECTOR_SC_EXT) &&
				      strcmp (plist_name, GRAVI_WAVE_DATA_SC_EXT) &&
				        strcmp (plist_name, "TEST_WAVE") &&
				        strcmp (plist_name, GRAVI_WAVE_DATA_FT_EXT) &&
					    strcmp (plist_name, GRAVI_IMAGING_DETECTOR_FT_EXT))){
			/*	Load the FT table	*/
			if (type_ext == 2)
				gravi_data_add (wave_map, img_plist,
							cpl_table_duplicate (gravi_data_get_table (wave_data, plist_name)));

			/*	Load the SC image_list	*/
			else if (type_ext == 3)
				gravi_data_add_cube (wave_map, img_plist,
							cpl_imagelist_duplicate (gravi_data_get_cube (wave_data, plist_name)));
		}
    }

//	wave_map = gravi_data_duplicate (wave_data);

	/* Add the offset table */
	plist = cpl_propertylist_duplicate (gravi_data_get_plist (wave_map,
			GRAVI_WAVE_DATA_SC_EXT));
	cpl_propertylist_append_string (plist,  "EXTNAME",
										"WAVE_OFFSET");
	gravi_data_add (wave_map, plist, cooef_table);
	cpl_propertylist_delete (plist);

	primary_hdr = gravi_data_get_header (wave_map);

	cpl_polynomial_delete(fit2d);
	cpl_propertylist_append_double (primary_hdr, QC_MINWAVE_SC,	minwave);
	cpl_propertylist_append_double (primary_hdr, QC_MAXWAVE_SC, maxwave);

	cpl_msg_info (cpl_func, "wave corrected : QC_MINWAVE_SC = %e QC_MAXWAVE_SC = %e",
			minwave, maxwave);

	gravi_data_set_table (wave_map, GRAVI_WAVE_DATA_SC_EXT, img_output);

	imglist_wave = cpl_imagelist_new ();
	cpl_imagelist_set(imglist_wave,image_wave , 0);
	cpl_imagelist_set(imglist_wave, cpl_image_duplicate (img_profile), 1);
	cpl_imagelist_set(imglist_wave, image_real, 2);
	gravi_data_set_cube (wave_map, "TEST_WAVE", imglist_wave);

	cpl_array_delete (dimension);

	return wave_map;
}

cpl_error_code gravi_data_save( gravi_data 		  * self,
								cpl_frameset 	  * allframes,
								const char 		  * filename,
								const cpl_parameterlist * parlist,
								cpl_frameset	  * usedframes,
								cpl_frame * frame,
								const char 		  * recipe,
								cpl_propertylist  * applist){
	cpl_frameset * frameset;
	int ext;
	/* Check the inputs */
	if ((filename == NULL) || (self == NULL)){
		cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT,
				                         "one of the inputs at least is NULL");
		return CPL_ERROR_NULL_INPUT;
	}

	/* Create the file and save the primary header. */

	frameset = cpl_frameset_duplicate(usedframes);
//	cpl_frame * frame = cpl_frameset_get_first(frameset);

	if (cpl_dfs_save_propertylist (allframes, self->primary_hdr, parlist,
			frameset, frame, recipe, applist,
			                 NULL, PACKAGE_STRING, filename ) != CPL_ERROR_NONE){
		cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_OUTPUT,
				          "Cannot save the first extension primary header");
		return CPL_ERROR_NULL_INPUT;
	}

	/* Save the extensions */
	for (ext = 0; ext < self->nb_ext; ext ++){
		if (self->exts_tbs[ext] != NULL)
			cpl_table_save(self->exts_tbs[ext], NULL, self->exts_hdrs[ext],
					                                filename, CPL_IO_EXTEND);
		else if (self->exts_imgl[ext] != NULL)
			cpl_imagelist_save (self->exts_imgl[ext], filename,
				   cpl_image_get_type (cpl_imagelist_get (self->exts_imgl[ext], 0)),
				   self->exts_hdrs[ext], CPL_IO_EXTEND);
	}
	cpl_frameset_delete (frameset);

	return CPL_ERROR_NONE;
}

/* 
 * Normalize the OI_FLUX (FLUX and FLUXERR columns) by the mean of FLUX(lbd)
 * Independly for each row.
 */

cpl_error_code gravi_normalize_flux (gravi_data * data)
{
	gravi_msg_function_start(0);

	/* Get the header */
	cpl_propertylist * hdr_data = gravi_data_get_header (data);
	
	/* Loop on FT and SC */
	for (int type_data = 0; type_data < 2; type_data ++) {
	  
	  /* Loop on polarisation */
	  int npol = gravi_pfits_get_pola_num (hdr_data, type_data);
	  for ( int pol= 0 ; pol < npol ; pol++ ) {

		cpl_msg_info (cpl_func, "Normalize the flux of %s (pol %i over %i)",GRAVI_TYPE(type_data),pol+1,npol);

		cpl_table * oi_flux = gravi_data_get_oi_flux (data, type_data, pol, npol);
		cpl_array ** flux    = cpl_table_get_data_array (oi_flux, "FLUX");
		cpl_array ** fluxErr = cpl_table_get_data_array (oi_flux, "FLUXERR");
		cpl_size nrow = cpl_table_get_nrow (oi_flux);
		CPLCHECK_MSG("Cannot get data");

		/* Loop on rows to normalize */
		for (cpl_size row = 0; row < nrow; row++) {
		  double flux_mean = cpl_array_get_mean (flux[row]);
		  cpl_array_divide_scalar (flux[row], flux_mean);
		  cpl_array_divide_scalar (fluxErr[row], flux_mean);
		}
	  
		CPLCHECK_MSG("Cannot divide by the mean flux");
	  } /* End loop on pol */
	} /* End loop on SC/FT */

	gravi_msg_function_exit(0);
	return CPL_ERROR_NONE;
}

/*---------------------------------------------------------------------------*/

cpl_vector * gravi_image_to_vector(cpl_image * img)
{
    cpl_ensure (img, CPL_ERROR_NULL_INPUT, NULL);
	
	cpl_vector * vector;
	cpl_type type_img;
	cpl_image * _image;
	int x, y;
	double * data;

	type_img = cpl_image_get_type (img);

	_image = cpl_image_cast (img, CPL_TYPE_DOUBLE);

	x = cpl_image_get_size_x (img);
	y = cpl_image_get_size_y (img);
	data = cpl_image_get_data_double(_image);
	vector = cpl_vector_new(x*y);
	for (int i = 0 ; i < x*y; i++)
		cpl_vector_set(vector, i, data[i]);

	cpl_image_delete(_image);
	return vector;
}

/*  Transform tangent plane coordinates into spherical
 *  (double precision). All in [rad]. Return RA in
 *
 *  Given:
 *     XI,ETA      dp   tangent plane rectangular coordinates
 *     RAZ,DECZ    dp   spherical coordinates of tangent point
 *
 *  Returned:
 *     RA,DEC      dp   spherical coordinates (0-2pi,+/-pi/2)
 *
 *  From:        sla_DRANRM
 *  P.T.Wallace   Starlink   24 July 1995
 *  Copyright (C) 1995 Rutherford Appleton Laboratory
 */
cpl_error_code gravi_dtps(double xi, double eta, double raz, double decz, double * ra, double * dec)
{
  double sdecz = sin(decz);
  double cdecz = cos(decz);
  double denom = cdecz - eta * sdecz;
  /* Compute new coordinates */
  *dec = atan2 (sdecz+eta*cdecz, sqrt(xi*xi + denom*denom));
  *ra  = atan2 (xi,denom) + raz;

  /* Make ra within 0-2pi */
  *ra = fmod (*ra, CPL_MATH_2PI);
  if ( *ra < 0.0 ) *ra += CPL_MATH_2PI;

  return CPL_ERROR_NONE;
}


int my_gsl_matrix_complex_fprintf(FILE *stream,gsl_matrix_complex *m,char *fmt)
{
        size_t rows=m->size1;
        size_t cols=m->size2;
        size_t row,col,ml;
        int fill;
        char buf[100];
        gsl_vector *maxlen;
        gsl_complex buff;

        maxlen=gsl_vector_alloc(cols);
        for (col=0;col<cols;++col) {
                ml=0;
                for (row=0;row<rows;++row) {
                		buff=gsl_matrix_complex_get(m,row,col);
                        sprintf(buf,"%g + i%g", GSL_REAL (buff), GSL_IMAG(buff));
                        if (strlen(buf)>ml)
                                ml=strlen(buf);
                }
                gsl_vector_set(maxlen,col,ml);
        }

        for (row=0;row<rows;++row) {
                for (col=0;col<cols;++col) {
						buff=gsl_matrix_complex_get(m,row,col);
						sprintf(buf,"%g + i%g", GSL_REAL (buff), GSL_IMAG(buff));
                       fprintf(stream,"%s",buf);
                        fill=gsl_vector_get(maxlen,col)+1-strlen(buf);
                        while (--fill>=0)
                                fprintf(stream," ");
                }
                fprintf(stream,"\n");
        }
        gsl_vector_free(maxlen);
        return 0;
}


/*----------------------------------------------------------------------------*/
/**
  @brief   The aim of this function is to identify the bad pixels and to update
  	  	  	  profile accordingly


  @param profile_map 	  	The FLAT calibration

  Pixel with flat value lower than a threshold are declared as
  bad pixels. the profile is set to 0 for this bad pixel

 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_compute_flat_badpix(gravi_data * flat, gravi_data * dark)
{
	double threshold;
	cpl_ensure_code (flat, CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (dark, CPL_ERROR_NULL_INPUT);
	
	cpl_msg_info(cpl_func, "Compute SC flat bad pixels");
	
	/* get the median flat */
	cpl_image * median_flat = cpl_imagelist_get(
					gravi_data_get_cube(flat, GRAVI_IMAGING_DATA_SC_EXT), 0);

	/* Compute bad pixels from flat*/
	cpl_image * grad = cpl_image_duplicate(median_flat);
	int x, size_x = cpl_image_get_size_x (median_flat);
	int y, size_y = cpl_image_get_size_y (median_flat);

	/* Mean profile */
	cpl_vector * vect_mean = cpl_vector_new_from_image_row (median_flat, 1);
	cpl_vector * vector;
	for (y = 2; y <= size_y; y++) {
		cpl_vector_add (vect_mean, cpl_vector_new_from_image_row (median_flat, y));
	}
	cpl_vector_divide_scalar (vect_mean, size_y);
	if (size_x < 60) { // case for LOW res
		vector = cpl_vector_filter_median_create(vect_mean, 2);
	}
	else{ // case for MED and HIGH
		vector = cpl_vector_filter_median_create (vect_mean, 5);
	}
	threshold = cpl_vector_get_max(vector)/100.;

	/* Threshold from DARK QC */
	cpl_propertylist * p_dark = gravi_data_get_header (dark);

	threshold = 20 * cpl_propertylist_get_double (p_dark, QC_DARKRMS_SC);
//	/* mask gradian */
//	cpl_image_shift(grad, 1, 0);
//	cpl_image_subtract(grad, median_flat);
//	cpl_image_abs(grad);
//	threshold = cpl_image_get_stdev(grad)*10;

	cpl_image_threshold(grad, threshold, threshold, 1., 0.);
	cpl_image * mask = cpl_image_duplicate(grad);

	/* add mask to FLAT data */
	cpl_imagelist * list_mask = cpl_imagelist_new();
	cpl_propertylist * plist = cpl_propertylist_duplicate(
			gravi_data_get_plist(flat, GRAVI_IMAGING_DATA_SC_EXT));
	cpl_propertylist_update_string(plist, "EXTNAME", "FLAT_MASK");
	cpl_imagelist_set(list_mask, grad, 0);
	gravi_data_add_cube(flat, plist, list_mask);
	CPLCHECK_INT("Error mask");

	/* Set badpix in profile */
	cpl_table * profiles = gravi_data_get_table(flat, "PROFILE_DATA");
	int region, nregion = cpl_table_get_ncol(profiles);
	char * regname;
	cpl_array ** array_data;
	double sum;
	cpl_size count;
	CPLCHECK_INT("Error mask");

	printf("%g \n", cpl_image_get(mask, 1, 1, NULL));
	CPLCHECK_INT("Error mask");

	for (region = 1; region <= nregion; region++) {
		regname = cpl_sprintf("DATA%d", region);
		printf("region %d \n", region);
		array_data=cpl_table_get_data_array(profiles, regname);
		count=0;
		for (x = 0; x < size_x; x++) {
			sum=0;
			printf("%d ", x);
			for (y = 0; y < size_y; y++) {
				printf("%d %g ", y, cpl_image_get(grad, x+1, y+1, NULL));
				if (cpl_image_get(grad, x+1, y+1, NULL) > 0.1){
					cpl_array_set(array_data[0], x+size_x*y, 0);
					count++;
					printf("1");
				}
				sum+=cpl_array_get(array_data[0], x+size_x*y, NULL);
			}
			if ( sum != 0 ){
				for (y = 0; y < size_y; y++) {
					cpl_array_set(array_data[0], x+size_x*y,
							cpl_array_get(array_data[0], x+size_x*y, NULL)/sum);
				}
			}
		}
	}
	cpl_msg_info(cpl_func, "Number of flat bad pixels : %lld", count);

	cpl_vector_delete(vect_mean);
	return(CPL_ERROR_NONE);
}



/*----------------------------------------------------------------------------*/
/**
  @brief    Create a yorick batch to execute mira on a OIFITS  file
  @param filename		batch filename
  @param input_file		input OIFITS file
  @param output_file	output FITS file

  @return   0 on success

  The function generate a string of FITS files automatically and randomly
 */
/*----------------------------------------------------------------------------*/

int gravi_write_yorick_batch(const char* filename, const char* input_file, const char* output_file)
{
	gravi_msg_function_start(0);
	FILE *fp;
	//char buffer[MAX_LEN_BUFFER];

	fp = fopen(filename,"w");

	fprintf(fp,
"/* \n"
"* mira-batch.i -\n"
"*\n"
"*	Mira batch generated by gravi_mira recipe \n"
"*/\n"
"include, \"mira.i\";\n"
"dataset = mira_new(\"%s\");\n"
"mira_config, dataset, dim=100, pixelsize=0.4*MIRA_MILLIARCSECOND, xform=\"exact\";\n"
"rgl = rgl_new(\"smoothness\");\n"
"dim = mira_get_dim(dataset);\n"
"img0 = array(double, dim, dim);\n"
"img0(dim/2, dim/2) = 1.0;\n"
"img1 = mira_solve(dataset, img0, maxeval=500, verb=0, xmin=0.0, normalization=1,\n"
"                  regul=rgl, mu=1e6);\n"
"fh = fits_create(\"%s\", overwrite=1, bitpix=-64, dimlist=dimsof(img0)); \n"
"fits_write_header, fh; \n"
"fits_write_array, fh, img1;\n"
"fits_close, fh \n"


"quit" , input_file, output_file);

	//fputs(buffer, fp);
	fclose(fp);

	gravi_msg_function_exit(0);
	return 0;
}


static void finals2000A_read_line(FILE *pFile, char *flag, double *mjd, double *pmx, double *pmy, double *dut)
{
  char line[LINE_SIZE];
  fread(line, LINE_SIZE, 1, pFile);
  *mjd = atof(line+7);
  *pmx = atof(line+18);
  *pmy = atof(line+37);
  *dut = atof(line+58);
  *flag = line[16];
}

static int finals2000A_mjd_first (FILE *pFile)
{
  char flag;
  double mjd, pmx, pmy, dut;
  fseek(pFile, 0, SEEK_SET);
  finals2000A_read_line(pFile, &flag, &mjd, &pmx, &pmy, &dut);
  return mjd;
}

static int finals2000A_mjd_last_type(FILE *pFile, char type)
{
  char flag;
  double mjd, pmx, pmy, dut;
  fseek(pFile, -LINE_SIZE, SEEK_END);
  finals2000A_read_line(pFile, &flag, &mjd, &pmx, &pmy, &dut);
  while( flag != type)
  {
    fseek(pFile, -2*LINE_SIZE, SEEK_CUR);
    finals2000A_read_line(pFile, &flag, &mjd, &pmx, &pmy, &dut);
  }
  fseek(pFile, 0, SEEK_SET);
  return mjd;
}

gravi_data * gravi_eop_load_finals2000A(const char *eop_file)
{
  gravi_msg_function_start(1);

  FILE *pFile;
  cpl_size last_mjd, n_entries;
  double *mjd, *pmx, *pmy, *dut;
  char **flag;
  gravi_data * eop_data = NULL;

  cpl_ensure (eop_file, CPL_ERROR_NULL_INPUT, NULL);
  
  /* Open finals2000A.data file */
  pFile = fopen ((char *)eop_file, "r");
  
  if(pFile != NULL)
  {
	cpl_msg_info (cpl_func, "Load the file: %s", eop_file);
	
    /* Find number of entries to last predicted entry */
    double mjd_P = finals2000A_mjd_last_type (pFile, 'P');
    double mjd_I = finals2000A_mjd_last_type (pFile, 'I');
	double mjd_S = finals2000A_mjd_first (pFile);
    n_entries = mjd_P - mjd_S;
    cpl_msg_info(cpl_func, "Reading %lli earth orientation parameters.", n_entries);
    cpl_msg_info(cpl_func, "  First entry:          MJD=%.1f", mjd_S);
    cpl_msg_info(cpl_func, "  Last IERS entry:      MJD=%.1f", mjd_I);
    cpl_msg_info(cpl_func, "  Last predicted entry: MJD=%.1f", mjd_P);

    /* Create tables */
    cpl_table * eop_table = cpl_table_new (n_entries);

    /* Create columns (filled with zero, thus valid) */
    gravi_table_new_column (eop_table, "MJD", "d", CPL_TYPE_DOUBLE);
    gravi_table_new_column (eop_table, "PMX", "arcsec", CPL_TYPE_DOUBLE);
    gravi_table_new_column (eop_table, "PMY", "arcsec", CPL_TYPE_DOUBLE);
    gravi_table_new_column (eop_table, "DUT", "s", CPL_TYPE_DOUBLE);
	cpl_table_new_column (eop_table, "FLAG", CPL_TYPE_STRING);
	cpl_table_fill_column_window_string (eop_table, "FLAG", 0, n_entries, " ");
	
    mjd = cpl_table_get_data_double (eop_table, "MJD");
    pmx = cpl_table_get_data_double (eop_table, "PMX");
    pmy = cpl_table_get_data_double (eop_table, "PMY");
    dut = cpl_table_get_data_double (eop_table, "DUT");
	flag = cpl_table_get_data_string (eop_table, "FLAG");

    /* Read finals2000A */
    for(int i=0; i<n_entries; i++)
    {
      finals2000A_read_line (pFile, flag[i]+0, mjd+i, pmx+i, pmy+i, dut+i);
    }

    /* close file */
    fclose (pFile);

    /* Build the gravi_data */
    eop_data = gravi_data_new (0);
    gravi_data_add (eop_data, NULL, eop_table);

	/* Create main header */
    cpl_propertylist * header = cpl_propertylist_new();
    cpl_propertylist_append_double (header, "ESO QC EOP_PARAM MJD_S", mjd_S);
    cpl_propertylist_append_double (header, "ESO QC EOP_PARAM MJD_I", mjd_I);
    cpl_propertylist_append_double (header, "ESO QC EOP_PARAM MJD_P", mjd_P);
    cpl_propertylist_append_double (header, "MJD-OBS", mjd_I);
    cpl_propertylist_append_string (header, "ESO PRO CATG", "EOP_PARAM");
    cpl_propertylist_append_string (header, "ESO PRO TECH", "CATALOG");
    cpl_propertylist_append_string (header, "ESO PRO TYPE", "IERS");
    gravi_data_append_header (eop_data, header);
  }
  else
  {
    cpl_msg_warning (cpl_func, "Cannot load the file: %s", eop_file);
    eop_data = NULL;
  }

  gravi_msg_function_exit(1);
  return eop_data;
}

// FLAT high frequency

/* Compute the flat, mean image of the flat multiplied by the profile of this region */
cpl_image * flat_profiled = cpl_image_extract (profile_mean, xmin, ymin, xmax, ymax);
cpl_image_multiply (flat_profiled, profile_crop);
cpl_image * specMean = cpl_image_collapse_window_create (flat_profiled, 1,1,nxc,nyc,0);
cpl_image_delete (flat_profiled);

CPLCHECK_NUL("Compute the flat");

/* Keep only the high frequencies in the FLAT, that is the detector
 * pixel to pixel efficiency difference. So that we remove the detector
 * imprint in the data, while keeping the overall instrumental transmission. */
if ( nxc > 5000 )
{
  cpl_msg_debug (cpl_func, "Spectra has >50 pixels -> flat high frequencies");
  cpl_image * specFlat = cpl_image_duplicate (specMean);
  cpl_mask * kernel = cpl_mask_new (11, 1);
  /* Faulty line, replaced by cpl_mask_not(kernel) */
  // for (int imask=1;imask<12;imask++) cpl_mask_set (kernel, i, 1, CPL_BINARY_1);
  cpl_mask_not (kernel);
  cpl_image_filter_mask (specFlat, specMean, kernel, CPL_FILTER_MEDIAN,CPL_BORDER_FILTER);
  cpl_image_divide (specMean, specFlat);
  cpl_image_delete (specFlat);
  cpl_mask_delete (kernel);
}
else
{
  cpl_image_fill_window (specMean, 1, 1, nxc, 1, 1.0);
  cpl_msg_debug (cpl_func, "Spectra has <50 pixels -> don't flat");
}

CPLCHECK_NUL("Creating spectrum flat");




	/* Fit the dispersion to ARGON */
	gravi_data * wave_data = gravi_compute_argon_wave (argon_data, profile_map, dark_map, badpix_map, parlist);
	CPLCHECK_CLEAN ("Cannot fit argon");

	cpl_table * wave_data_sc = gravi_data_get_table (wave_data, "WAVE_DATA_SC");
	cpl_table * wave_map_sc = gravi_data_get_table (wave_map, "WAVE_DATA_SC");

	cpl_propertylist * plist = gravi_data_get_plist (wave_data, "WAVE_DATA_SC");
	gravi_data_set_propertylist (wave_map, "WAVE_DATA_SC", plist);
	gravi_data_set_table (wave_map, "WAVE_DATA_SC", cpl_table_duplicate (wave_data_sc));

	cpl_propertylist * header = gravi_data_get_header (wave_map);
	plist = gravi_data_get_header (wave_data);
	cpl_propertylist_update_double (header, "ESO QC MINWAVE SC", cpl_propertylist_get_double (plist, "ESO QC MINWAVE SC"));
	cpl_propertylist_update_double (header, "ESO QC MAXWAVE SC", cpl_propertylist_get_double (plist, "ESO QC MAXWAVE SC"));

	gravi_data_save_data (wave_data, "wave_map_argon.fits", CPL_IO_CREATE);
	FREE (gravi_data_delete, wave_data);





gravi_data * gravi_visdata_fromellipse (gravi_data * spectrum_data)
{
  gravi_msg_function_start(1);
  cpl_ensure (spectrum_data, CPL_ERROR_NULL_INPUT, NULL);

  cpl_propertylist * plist = cpl_propertylist_new ();

  /* Create the output */
  gravi_data * output = gravi_data_new (0);
  cpl_propertylist * header = gravi_data_get_header (spectrum_data);
  gravi_data_append_header (output, cpl_propertylist_duplicate (header));

  /* Loop on SC/FT and polarisation*/
  for (int type_data = 0; type_data < 2; type_data ++ ) {
	int npol = gravi_pfits_get_pola_num (header, type_data);
	for (int pol = 0; pol < npol; pol ++) {

	  /* Get the data */
	  cpl_table * spectrum_table = gravi_data_get_spectrum_data (spectrum_data, type_data);
	  cpl_table * img_det = gravi_data_get_imaging_detector (spectrum_data, type_data);
	  
	  cpl_size nwave = cpl_table_get_column_depth (spectrum_table,"DATA1");
	  cpl_size nrow  = cpl_table_get_nrow (spectrum_table);
	  CPLCHECK_NUL ("Cannot get tables");

	  /* Create output table */
	  cpl_size nbase = 6;
	  cpl_table * vis_table = cpl_table_new (nbase*nrow);
	  gravi_table_new_column (vis_table, "TIME", "us", CPL_TYPE_INT);
	  gravi_table_new_column_array (vis_table, "VISDATA", "e", CPL_TYPE_FLOAT_COMPLEX, nwave);

	  /* Create memory for output */
	  float complex** mVis = cpl_malloc (nbase*nrow * sizeof(float complex*));
	  for (cpl_size row = 0; row < nbase*nrow; row ++) 
		mVis[row] = cpl_malloc (nwave * sizeof(float complex));

	  /* Create vectors over row direction, for the ellipse */
	  cpl_vector * vectX = cpl_vector_new (nrow);
	  cpl_vector * vectY = cpl_vector_new (nrow);

	  /* Get pointer to speed up */
	  cpl_array ** tVis = cpl_table_get_data_array (vis_table, "VISDATA");
	  double * datX = cpl_vector_get_data (vectX);
	  double * datY = cpl_vector_get_data (vectY);

	  /* Loop on base */
	  for (int base = 0; base < nbase; base++) {
		cpl_msg_info_overwritable (cpl_func,"Ellipse to type=%s, pol=%i over %i, base=%i",
								   GRAVI_TYPE(type_data), pol+1, npol, base);

		/* Set the time */
		for (cpl_size row = 0; row < nrow; row ++)  {
		  cpl_size irow = row*nbase + base;
		  cpl_table_set (vis_table, "TIME", irow,
						 cpl_table_get (spectrum_table, "TIME", row, NULL));
		}
		
		/* Get the sign of the phase for this baseline */
		int phi_sign = gravi_table_get_phase_sign (img_det, GRAVI_BASE_TEL[base][0]+1, GRAVI_BASE_TEL[base][1]+1);

		/* Get regions */
		int regA = gravi_get_region (img_det, base, 'A', pol);
		int regB = gravi_get_region (img_det, base, 'B', pol);
		int regC = gravi_get_region (img_det, base, 'C', pol);
		int regD = gravi_get_region (img_det, base, 'D', pol);
		cpl_ensure (regA>=0 && regA>=0 && regA>=0 && regA>=0,
					CPL_ERROR_ILLEGAL_INPUT, NULL);

		/* Get pointer to data to speed up */
		double ** datA = gravi_table_get_data_array_double (spectrum_table, GRAVI_DATA[regA]);
		double ** datB = gravi_table_get_data_array_double (spectrum_table, GRAVI_DATA[regB]);
		double ** datC = gravi_table_get_data_array_double (spectrum_table, GRAVI_DATA[regC]);
		double ** datD = gravi_table_get_data_array_double (spectrum_table, GRAVI_DATA[regD]);
		CPLCHECK_NUL ("Cannot get data");

		/* Loop on wavelength because we need the
		 * ellipse vectors in the row direction */
		for (cpl_size wave = 0 ; wave < nwave ; wave++) {

		  /* Create the ellipse vectors as:
		   *    X = C-A   and   Y = D-B   */
		  for (cpl_size row = 0; row < nrow; row ++)  {
			datX[row] = datC[row][wave] - datA[row][wave];
			datY[row] = datD[row][wave] - datB[row][wave];
		  }
		  
		  /* Re-center the ellipse to ease fit */
		  cpl_vector_subtract_scalar (vectY, cpl_vector_get_mean (vectY));
		  cpl_vector_subtract_scalar (vectX, cpl_vector_get_mean (vectX));

		  /* Get the phase by fitting the ellipse */
		  cpl_vector * phase = gravi_vectors_phase_create (vectX, vectY);
		  cpl_vector_multiply_scalar (phase, phi_sign);

		  /* Set the VISDATA = sqrt (X^2 + Y^2) * exp (i.phase) */
		  for (cpl_size row = 0; row < nrow; row ++)  {
                      cpl_size irow = row*nbase + base;
		      mVis[irow][wave] = (float complex)
			  (sqrt (datX[row]*datX[row] + datY[row]*datY[row]) *
			   cexp (1.*I * cpl_vector_get (phase, row)));
		  }

		  FREE (cpl_vector_delete, phase);
		} /* End loop on wave */

		/* Wrap the array into the table */
		for (cpl_size row = 0; row < nrow*nbase; row ++)
		  tVis[row] = cpl_array_wrap_float_complex (mVis[row], nwave);

		/* Free pointer to data array */
		FREE (cpl_free, datA);
		FREE (cpl_free, datB);
		FREE (cpl_free, datD);
		FREE (cpl_free, datD);
		
	  } /* End loop on base */

	  /* Desalocate ellipses vectors */
	  FREE (cpl_vector_delete, vectX);
	  FREE (cpl_vector_delete, vectY);
	  
	  /* Add table to data with EXTNAME */
	  cpl_propertylist_append_string (plist, "EXTNAME", "ELLIPSE_VIS");
	  cpl_propertylist_append_string (plist, "INSNAME", GRAVI_INSNAME(type_data,pol,npol));
	  cpl_propertylist_append_int (plist, "EXTVER", GRAVI_EXTVER(type_data,pol,npol));
	  
	  gravi_data_add (output, plist, vis_table);

	  
	}	/* End loop on pol */
  }	/* End loop on SC/FT */

  gravi_msg_function_exit(1);
  return output;
}

double gravi_pfits_get_lambdamet(const cpl_propertylist * plist)
{
    cpl_errorstate prestate = cpl_errorstate_get();
    double value = cpl_propertylist_get_double(plist, "ESO INS MLC WAVELENG");
    cpl_ensure (cpl_errorstate_is_equal(prestate), cpl_error_get_code(), 0.0);
    return value;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    This function smooth the SC preproc data in the spectral direction
            in place, only for preproc (per region). Note that the data are
			copied during the execution.

  @param data 	  	 The gravi_data to smooth
  @param parlist     The options of the recipe

  @return   cpl_error_code
*/
cpl_error_code gravi_smooth_preproc (gravi_data * data,
									 const cpl_parameterlist * parlist)
{
  cpl_table * detector_table, * spectrum_table;
  cpl_size row, n_row, reg, n_region;
  cpl_array * output = NULL;
  
  /* Verbose */
  gravi_msg_function_start(1);
  cpl_ensure_code (data,    CPL_ERROR_NULL_INPUT);
  cpl_ensure_code (parlist, CPL_ERROR_NULL_INPUT);

  /* Get the nsmooth from parameter list */
  int nsmooth = cpl_parameter_get_int (cpl_parameterlist_find_const (parlist, "gravi.nsmooth_sc"));

  /* Default from spectral resolution */
  if (nsmooth < 0) {
	const char * resolution  = gravi_pfits_get_spec_res (gravi_data_get_header (data));
	
	if ( !(strcmp(resolution, "HIGH")) ) {
	  cpl_msg_info(cpl_func,"Default smoothing is 15 (HIGH)");
	  nsmooth = 15;
	}
	else if ( !(strcmp(resolution, "MEDIUM")) ) {
	  cpl_msg_info(cpl_func,"Default smoothing is 3 (MEDIUM)");
	  nsmooth = 3;
	}
	else if ( !(strcmp(resolution, "LOW")) ) {
	  cpl_msg_info(cpl_func,"Default is no smoothing of data (LOW)");
	  nsmooth = 0;
	} else {
	  cpl_msg_warning(cpl_func,"Unknown spectral resolution thus no smoothing of data");
	  nsmooth = 0;
	}
  }
  
  /* Case no smoothing */
  if ( nsmooth <= 0 ) {
	cpl_msg_info (cpl_func, "End function (no smoothing)");
	return CPL_ERROR_NONE;
  }
  
  /* Search the number of regions for SC */
  detector_table = gravi_data_get_table (data, GRAVI_IMAGING_DETECTOR_SC_EXT);
  n_region = cpl_table_get_nrow (detector_table);

  /* Get the SC table data */
  spectrum_table = gravi_data_get_table (data, GRAVI_SPECTRUM_DATA_SC_EXT);
  n_row = cpl_table_get_nrow (spectrum_table);

  CPLCHECK_MSG ("Cannot get data");

  /* Loop on SC regions */
  for (reg = 0; reg < n_region; reg++) {
	
	/* Load pointer to all rows */
	cpl_msg_info_overwritable(cpl_func,"Smooth region %lld over %lld of SC", reg+1, n_region);
	cpl_array** arrays = cpl_table_get_data_array (spectrum_table, GRAVI_DATA[reg]);

	/* Loop on rows */
	for (row = 0; row < n_row; row++) {
	  /* Compute a smoothed version of the spectra of this row */
	  output = gravi_array_smooth (arrays[row], nsmooth);
	  /* Put it back inplace in the table */
	  cpl_array_delete (arrays[row]);
	  arrays[row] = output;
	}
	/* End loop on rows */
  }
  /* End loop on region */

  /* Add the NSMOOTH parameter */
  cpl_propertylist * hdr_data = gravi_data_get_header (data);
  cpl_propertylist_append_int (hdr_data, "ESO QC P2VM NSMOOTH SC", nsmooth);
  cpl_propertylist_set_comment (hdr_data, "ESO QC P2VM NSMOOTH SC", "nb of smoothed channels in P2VM computation");
    
  /* Verbose */
  gravi_msg_function_exit(1);
  return CPL_ERROR_NONE;
}


cpl_vector * gravi_construction_opd_phase (cpl_table * opl_table,
										   cpl_table * phase_sc,
										   cpl_table * phase_ft,
										   double dit_sc)
{
	int ind_sc, nrow, row, nbase = 6, base;
	double exptime, exptime_sc, opl, time_sc, time_ft, time_metrology;
	int nv, comp, nrow_met, nrow_sc;
	int tel_1[6] = {0,0,0,1,1,2};
	int tel_2[6] = {1,2,3,2,3,3};
	const cpl_array * time_array_ft, * time_array_sc;

    gravi_msg_function_start(0);
	cpl_ensure (opl_table, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (phase_sc,  CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (phase_ft,  CPL_ERROR_NULL_INPUT, NULL);
	
	/* Get the number of acquisitions */
	time_array_ft = cpl_table_get_array (phase_ft, "TIME", 0);
	time_array_sc = cpl_table_get_array (phase_sc, "TIME", 0);
	nrow =  cpl_array_get_size (time_array_ft);
	nrow_sc = cpl_array_get_size (time_array_sc);
	nrow_met = cpl_table_get_nrow (opl_table);

	CPLCHECK_NUL("Cannot get data");
	
	/* allocate the matrices */
	gsl_matrix * A_data = gsl_matrix_calloc (nbase * nrow, 8), * A;
	gsl_matrix * U;
	gsl_vector * bis_data = gsl_vector_alloc (nbase * nrow), * bis;
	gsl_matrix * X = gsl_matrix_alloc (8, 8),
					* V = gsl_matrix_alloc (8, 8);
	gsl_vector * S = gsl_vector_alloc (8);
	gsl_vector * work = gsl_vector_alloc (8), * x = gsl_vector_alloc (8);
	CPLCHECK_NUL("Allocate matrix");

	/* compute the periods of the signals */
	exptime = cpl_array_get_int (time_array_ft, 1, &nv) -
						cpl_array_get_int (time_array_ft, 0, &nv);
	exptime_sc = cpl_array_get_int (time_array_sc, 1, &nv) -
						cpl_array_get_int (time_array_sc, 0, &nv);
	//double exptime_met = cpl_table_get_int (opl_table, "TIME", 1, &nv) -
	//		cpl_table_get_int (opl_table, "TIME", 0, &nv);

	/*
	 *  Extract the OPl metrology and compute the mean of each exposure time
	 */

	for (base = 0; base < nbase; base ++){
		ind_sc = 0;
		int im = 0;
		const cpl_array * ft_array = cpl_table_get_array (phase_ft, "PHASE", base);
		const cpl_array * sc_array = cpl_table_get_array (phase_sc, "PHASE", base);
		char * opl1 = cpl_sprintf("OPL%d", tel_1[base]+1);
		char * opl2 = cpl_sprintf("OPL%d", tel_2[base]+1);


		cpl_array * temp = cpl_array_wrap_double (
				cpl_table_get_data_double (opl_table, opl1), nrow_met);
		cpl_array * opl_met1 = cpl_array_duplicate (temp);
		cpl_array_unwrap (temp);
		cpl_array * opl_met2 = cpl_array_wrap_double (
				cpl_table_get_data_double (opl_table, opl2), nrow_met);
		cpl_array_subtract (opl_met1, opl_met2);
		CPLCHECK_NUL ("Opd diff");

		/*
		 * resample the SC and Metrology at FT period
		 */
		for (row = 0; row < nrow; row ++){ // loop on row FT

			gsl_matrix_set (A_data, base * nrow + row, 6,
					- cpl_array_get_double (ft_array, row, &nv));
			CPLCHECK_NUL ("Set FT phase in matrix");

			time_ft = cpl_array_get_int (time_array_ft, row, &nv);
			if (ind_sc < cpl_array_get_size (sc_array)){
				time_sc = cpl_array_get_int (time_array_sc, ind_sc, &nv);
			}
			CPLCHECK_NUL ("Get time sc");

			/* increase ind ft if time_ft in the range of the nex DIT */
			while(time_ft  > (time_sc+exptime_sc/2.)){
				ind_sc++;
				if (ind_sc < cpl_array_get_size (sc_array))
					time_sc = cpl_array_get_int (time_array_sc, ind_sc, &nv);
				else
					break;
			}
			CPLCHECK_NUL ("Get next sc dit");

			/* get phi sc */
			double phi_sc=0;
			if (ind_sc > 0){
				if(ind_sc < nrow_sc) {
					phi_sc=cpl_array_get_double (sc_array, ind_sc, &nv);
				}
				else {
					phi_sc=cpl_array_get_double (sc_array, nrow_sc - 1, &nv);
				}
			}
			else {
				phi_sc = cpl_array_get_double (sc_array, 0, &nv);
			}
			CPLCHECK_NUL ("Get phi_sc");

			gsl_matrix_set (A_data, base * nrow + row, 7, phi_sc);
			gsl_matrix_set (A_data, base * nrow + row, base, 1.0);


			/* Metrology case */
			opl = 0;
			comp = 0;

			if (im < nrow_met)
			time_metrology = cpl_table_get_int (opl_table, "TIME", im, &nv);

			while ((time_metrology  < (time_ft + exptime))){ //((time_metrology + exptime_met)< time_ft){

				if (im < nrow_met) {
					opl += cpl_array_get_double (opl_met1, im, &nv); //cpl_vector_get (vector_opd, im);
					comp ++;
				}

				im++;
				if (im < nrow_met) {
					time_metrology = cpl_table_get_int (opl_table, "TIME", im, &nv);
				}
				else {
					break;
				}

				CPLCHECK_NUL("Get time metrology");
			}

			/* average the metrology over the FT DIT */
			if (comp != 0)
				gsl_vector_set (bis_data, base * nrow + row, opl/comp);

			/* if no metrology signal within the FT DIT interpolate the metrology */
			else {
				if (im > 0){
					if(im < nrow_met) {
						opl = (cpl_array_get_double (opl_met1, im, &nv) -
							cpl_array_get_double (opl_met1, im - 1, &nv)) * (time_ft -
							cpl_table_get_int (opl_table, "TIME", im - 1, &nv))
							/(time_metrology -
							cpl_table_get_int (opl_table, "TIME", im - 1, &nv)) +
							cpl_array_get_double (opl_met1, im - 1, &nv);

						gsl_vector_set (bis_data, base * nrow + row, opl);
					}
					else {
						gsl_vector_set (bis_data, base * nrow + row,
							cpl_array_get_double (opl_met1, nrow_met - 1, &nv));
					}
				}
				else {
					gsl_vector_set (bis_data, base * nrow + row,
							cpl_array_get_double (opl_met1, 0, &nv));
				}
				CPLCHECK_NUL("Interpolate the metrology");
			}
		} /* end loop on row FT */

		cpl_array_delete (opl_met1);
		cpl_array_unwrap (opl_met2);
		cpl_free(opl1);
		cpl_free(opl2);
	}

	/*
	 *  filter the data out of the SC integration time
	 */

	long n_row_A=nrow_sc*(int)(dit_sc/exptime+1);
	cpl_vector *i_A_vector = cpl_vector_new(n_row_A);
	int i_A=0;
	int time_mod;
	int t0_sc=cpl_array_get_int (time_array_sc, 0, &nv);
	int tend_sc=cpl_array_get_int (time_array_sc, nrow_sc-1, &nv);

	/* Find the index of the frames within the SC integration */
	for (row=0; row<nrow; row++){
		time_ft = cpl_array_get_int (time_array_ft, row, &nv);
		if ((time_ft >= t0_sc-dit_sc/2) && (time_ft < tend_sc+dit_sc/2)){
			time_mod=((int)(time_ft-(t0_sc-dit_sc/2)) % (int)exptime_sc);
			if ( time_mod >= 0 && time_mod < dit_sc ) {
				cpl_vector_set (i_A_vector, i_A, row);
				i_A++;
			}
		}
	}
	n_row_A=i_A;

	/* copy the frames within SC intergration in the matrix A */
	A = gsl_matrix_alloc (nbase * n_row_A, 8);
	bis = gsl_vector_alloc (nbase * n_row_A);
	cpl_vector *time_A = cpl_vector_new(n_row_A);
	int col;

	for (base = 0; base < nbase; base ++)
		for (i_A = 0; i_A < n_row_A; i_A++){
			row = cpl_vector_get(i_A_vector, i_A);
			for (col = 0; col < 8; col ++)
				gsl_matrix_set (A, base * n_row_A + i_A, col, gsl_matrix_get (A_data, base * nrow + row, col));

			gsl_vector_set (bis, base * n_row_A + i_A, gsl_vector_get (bis_data, base * nrow + row));
			if (base == 0)
				cpl_vector_set (time_A, i_A, cpl_array_get_int (time_array_ft, row, &nv));
		}

	cpl_vector_delete(i_A_vector);

	/*
	 * Solve the linear equation Ax = b to get the coefficients to deduce the OPDs
	 */
	nrow=n_row_A;
	U = gsl_matrix_alloc (nbase * nrow, 8);
	gsl_matrix_memcpy (U, A);


	gsl_linalg_SV_decomp (U, V, S, work);
	gsl_linalg_SV_solve (U, V, S, bis, x);
	cpl_vector * opd_coeff = cpl_vector_new (2);
	cpl_vector_set (opd_coeff, 0, gsl_vector_get (x, 7));
	cpl_vector_set (opd_coeff, 1, gsl_vector_get (x, 6));

	cpl_msg_debug(cpl_func, "Wavelength base %d => SC = %g, FT = %g\n",
			base, -cpl_vector_get(opd_coeff,0)*2*M_PI, cpl_vector_get(opd_coeff,1)*2*M_PI);

	// TODO
	if (PLOT_MET_PHASE_FIT)

	{
		cpl_msg_info(cpl_func, "Plot fit residuals");
		cpl_errorstate prestate = cpl_errorstate_get();
		gsl_matrix_memcpy (U, A);

		const cpl_vector ** vectors=malloc(3 * sizeof(cpl_vector*));
		vectors[0]=NULL;
		cpl_vector *vect_phase_sc_ft=cpl_vector_new(nrow*nbase);
		cpl_vector *vect_dopd_met=cpl_vector_new(nrow*nbase);
		cpl_vector *vect_diff=cpl_vector_new(nrow*nbase);
		for (row = 0; row < nrow*nbase; row ++){
			cpl_vector_set(vect_dopd_met, row, gsl_vector_get (bis, row));
			cpl_vector_set(vect_phase_sc_ft, row, gsl_matrix_get (A, row, 0)*gsl_vector_get (x, 0)+
					gsl_matrix_get (A, row, 1)*gsl_vector_get (x, 1)+
					gsl_matrix_get (A, row, 2)*gsl_vector_get (x, 2)+
					gsl_matrix_get (A, row, 3)*gsl_vector_get (x, 3)+
					gsl_matrix_get (A, row, 4)*gsl_vector_get (x, 4)+
					gsl_matrix_get (A, row, 5)*gsl_vector_get (x, 5)+
					gsl_matrix_get (A, row, 6)*gsl_vector_get (x, 6)+
					gsl_matrix_get (A, row, 7)*gsl_vector_get (x, 7));
			cpl_vector_set(vect_diff, row, cpl_vector_get(vect_dopd_met, row)-cpl_vector_get(vect_phase_sc_ft, row));
		}

		vectors[1]=vect_phase_sc_ft;
		vectors[2]=vect_dopd_met;

		cpl_plot_vectors (cpl_sprintf("set title 'Met fit  case %d a=%g  b=%g'; set xlabel 'time[10-6s]'; set ylabel 'Phase_ft+Phase_sc, dopd_met';", base+1, gsl_vector_get (x, 0), gsl_vector_get (x, 1)),
				 "", "", vectors, 3);
		cpl_plot_vector(cpl_sprintf("set title 'Met fit  case %d a=%g  b=%g'; set xlabel 'time[10-6s]'; set ylabel 'Phase_ft+Phase_sc-dopd_met';", base+1, gsl_vector_get (x, 0), gsl_vector_get (x, 1)),
						"", "", vect_diff);
		cpl_vector_delete(vect_phase_sc_ft);
		cpl_vector_delete(vect_dopd_met);
		cpl_vector_delete(vect_diff);
		cpl_free(vectors);
		cpl_errorstate_set (prestate);
	}

	gsl_matrix_free (A);
	gsl_matrix_free (A_data);
	gsl_matrix_free (U);
	gsl_matrix_free (V);
	gsl_vector_free (x);
	gsl_vector_free (bis);
	gsl_vector_free (bis_data);
	cpl_vector_delete(time_A);
	gsl_matrix_free (X);
	gsl_vector_free (S);
	gsl_vector_free (work);
	
    gravi_msg_function_exit(0);
	return opd_coeff;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    The given output FITS file contain a p2vm table of the metrology
  	        with the values of the transmission, phase and coherence extract
  	        using the p2vm matrix

  @param metrology_table 	The metrology table from a fits file wish all
  	  	  	  	  	  	  	the shutters are open
  @param opl_table			This table will contain the opl of each diode
                            and telescope and must allocated before using
                            this function
  @param mjd_obs			The Modified Julian Date when metrology stated

  @return   The p2vm table of the metrology with the values of the
  transmission, phase and coherence

  The function returns a pointer of cpl_table compute the p2vm table
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_metrology_calibration (cpl_table * metrology_table,
										 cpl_table * opl_table,
										 double mjd_obs)
{
	gravi_msg_function_start(1);
	cpl_ensure (metrology_table, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (opl_table,       CPL_ERROR_NULL_INPUT, NULL);
	
	cpl_table * p2vm_met;

	int row, i, nb_row;
	cpl_vector * vectA, * vectB,
			   * phase, * y_sigma, * init_val;
	cpl_array * volt;
	cpl_array * temp, * coherence_fit, * phase_fit;
	const cpl_array * volt0;
	cpl_matrix * opd_matrix;
	int tel, infos = 0, n_tel;
	cpl_vector ** vect;
	int val_to_fit2[] = {1,1,1,0}, nv;
	double mse, red_chisq;
	cpl_bivector *plot;
	char* ps_string;
	char * opl;
	cpl_array * trans;
	
	/* Extract inputs */

	nb_row = cpl_table_get_nrow (opl_table);
	for (tel = 0; tel < 4; tel ++){
		opl = cpl_sprintf("OPL%d", tel + 1);
		cpl_table_new_column (opl_table, opl, CPL_TYPE_DOUBLE);
		cpl_free (opl);
	}

	/* Define the first ABCD k band and in witch case it starts */
	volt0 = cpl_table_get_array (metrology_table, "VOLT", 0);


	n_tel = cpl_array_get_size (volt0);

	/* Creating the 4 columns (REGNAME, TRANSMISSION, COHERENCE and PHASE) */
	/** Regname **/

	int start_tel = n_tel/2 - 8;// 0
	p2vm_met = cpl_table_new (n_tel / 2 - start_tel);
	cpl_table_new_column_array (p2vm_met,"TRANSMISSION", CPL_TYPE_DOUBLE, 2);
	cpl_table_new_column_array (p2vm_met,"COHERENCE", CPL_TYPE_DOUBLE, 2);
	cpl_table_new_column_array (p2vm_met,"PHASE", CPL_TYPE_DOUBLE, 2);

	/* Create a new table to stoke the OPL of each telescope */
	cpl_table_new_column (opl_table, "TIME", CPL_TYPE_INT);

	if (cpl_error_get_code()){

		printf("error %f  %s and %s  \n", 1.0, cpl_error_get_message(), cpl_error_get_where());
		return NULL;
	}


	for (tel = n_tel/2 - 8; tel < n_tel/2; tel++){

		vectA = cpl_vector_new (nb_row);
		vectB = cpl_vector_new (nb_row);

		/* Get the index of the k band AB */

		/* Define the first AB k band and in witch case it starts */
		for (row = 0; row < nb_row; row ++){
			temp = cpl_array_duplicate (cpl_table_get_array (metrology_table,
														"VOLT", row));
			volt = cpl_array_cast (temp, CPL_TYPE_DOUBLE);
			cpl_array_delete (temp);
			cpl_vector_set (vectA, row, cpl_array_get_double (volt, 2*tel, &nv));
			cpl_vector_set (vectB, row, cpl_array_get_double (volt, 2*tel + 1, &nv));

			cpl_array_delete (volt);
		}

		/* Fit OPD */
		/* Compute the OPD */

		phase = gravi_vectors_phase_create (vectA, vectB);

		if (PLOT_MET_ELLIPSE)
		{
			if (POSTSCRIPT_PLOT) ps_string=cpl_sprintf("set term png; set output 'plot_met_ellipse_T%d.png';", tel);
			else  ps_string=cpl_sprintf(" ");
			plot = cpl_bivector_wrap_vectors (vectA, vectB);
			cpl_plot_bivector (cpl_sprintf("set title 'Met ellipse Tel %d'; set xlabel 'C-A'; set ylabel 'D-B';", tel+1),
					NULL, NULL, plot);

			if (POSTSCRIPT_PLOT) ps_string=cpl_sprintf("set term png; set output 'plot_met_phase_T%d.png';", tel);
			else  ps_string=cpl_sprintf(" ");
			cpl_plot_vector (cpl_sprintf("set title 'Met phase Tel %d'; set xlabel 'index'; set ylabel 'Phase';", tel+1),
					NULL, NULL, phase);

			cpl_free(ps_string);
		}

		
		if (cpl_error_get_code()){
			printf("error %f  %s and %s  \n", 4.0, cpl_error_get_message(), cpl_error_get_where());

			return NULL;
		}
		/* Get the OPD vector from the phase */

		opd_matrix = cpl_matrix_new(nb_row, 1);

		for (i = 0; i < nb_row; i++){

			cpl_matrix_set(opd_matrix, i, 0,
						  cpl_vector_get(phase, i) * LAMBDA_MET/(2*M_PI));
			if (tel >= n_tel/2 - 8){

				if (tel < n_tel/2 - 4){
					opl = cpl_sprintf("OPL%d", tel - (n_tel/2 - 8) + 1);
					cpl_table_set_double (opl_table, opl, i,
						cpl_vector_get (phase, i) * LAMBDA_MET/(2*M_PI));
					cpl_free (opl);
				}
				else {
					opl = cpl_sprintf("OPL%d", tel - (n_tel/2 - 4) + 1);
					double opl_ft = cpl_table_get_double (opl_table, opl, i, &nv);
					cpl_table_set_double (opl_table, opl, i, (cpl_vector_get (phase, i) * LAMBDA_MET/(2*M_PI))- opl_ft);
					cpl_free (opl);
				}

				if (tel == n_tel/2 - 4)
					cpl_table_set_int (opl_table, "TIME", i,
						cpl_table_get_int (metrology_table, "TIME",
													i, &nv) + mjd_obs);
			}
			if (cpl_error_get_code()){
				printf("error %f  %s and %s  \n", 5.0, cpl_error_get_message(), cpl_error_get_where());

				return NULL;
			}
		}

		cpl_vector_delete(phase);

		/* End fit opd */

		/* fit on a central window of 5*lambda width */

		vect = cpl_malloc(2*sizeof(cpl_vector*));
		vect[0] = vectA;
		vect[1] = vectB;
		phase_fit = cpl_array_new(2, CPL_TYPE_DOUBLE);
		coherence_fit = cpl_array_new(2, CPL_TYPE_DOUBLE);
		trans = cpl_array_new(2, CPL_TYPE_DOUBLE);

		for (i = 0; i < 2; i++){

			/* Fit to get the phase and the coherence */


			/* Get the spectrum of the vector region */
			y_sigma = cpl_vector_new(nb_row);
			cpl_vector_fill(y_sigma, 1);

			/* Define and initialize all variables to make a FIT */
			mse = 0;
			red_chisq = 0;
			init_val = cpl_vector_new(4);
			cpl_vector_set(init_val, 0, 1);
			cpl_vector_set(init_val, 1, 1);
			cpl_vector_set(init_val, 2, 1);
			cpl_vector_set(init_val, 3, LAMBDA_MET);

			cpl_errorstate prestate = cpl_errorstate_get();
			cpl_fit_lvmq(opd_matrix, NULL, vect[i],
			y_sigma, init_val, val_to_fit2, &sin_lambda,
			&dfda_sin, CPL_FIT_LVMQ_TOLERANCE, CPL_FIT_LVMQ_COUNT,
			CPL_FIT_LVMQ_MAXITER, &mse, &red_chisq, NULL);
			
			if (cpl_error_get_code()){
				printf("error %f  %s and %s  \n", 6.0, cpl_error_get_message(), cpl_error_get_where());
				return NULL;
			}

			if (!strcmp("The iterative process did not converge",
				cpl_error_get_message())){
				if (infos)
					cpl_msg_info(cpl_func, "The iterative process "
											"did not converge");
				cpl_errorstate_set (prestate);
			}
			if (infos)
				cpl_msg_info(cpl_func, "tel = %d : mse "
						"%g chi2 %g", tel, mse, red_chisq);



			/* Compute the P2VM value */
			cpl_array_set_double (coherence_fit, i,
					sqrt( pow( cpl_vector_get(init_val, 2), 2) +
						  pow( cpl_vector_get(init_val, 1), 2)));

			cpl_array_set_double (phase_fit, i, atan2( cpl_vector_get(init_val, 2),
								cpl_vector_get(init_val, 1)));

			cpl_array_set_double (trans, i, cpl_vector_get(init_val, 0));

			if (cpl_error_get_code()){
				printf("error %f  %s and %s  \n", 7.0, cpl_error_get_message(), cpl_error_get_where());
				return NULL;
			}

			cpl_vector_delete(init_val);
			cpl_vector_delete(y_sigma);
			cpl_vector_delete (vect[i]);
		}

		cpl_free (vect);
		cpl_matrix_delete(opd_matrix);

		cpl_table_set_array (p2vm_met, "COHERENCE", tel - start_tel, coherence_fit);
		cpl_array_subtract_scalar (phase_fit,
				cpl_array_get_double (phase_fit, 0, &nv));
		cpl_table_set_array (p2vm_met, "PHASE", tel - start_tel, phase_fit);
		cpl_table_set_array (p2vm_met, "TRANSMISSION", tel - start_tel, trans);

		cpl_array_delete (trans);
		cpl_array_delete (coherence_fit);
		cpl_array_delete (phase_fit);
	}

	/* Verbose */
	gravi_msg_function_exit(1);
	return p2vm_met;

}



/*----------------------------------------------------------------------------*/
/**
 * @brief Compute mean position of MET and FDDL in each DITs
 * 
 * @param oi_vis:       input/output OIFITS
 * @param p2vmred_data:  
 * @param preproc_data:   
 * 
 * We shall deprecate this function and use instead the already
 * computed FDDL and MET values from the gravi_compute_signals.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_data_mean_metFddl (gravi_data * oi_vis, gravi_data * p2vmred_data,
									    gravi_data *  preproc_data,
									    int thread){

	cpl_table * visMet_data, * fddl_data, * spectrum_sc;
	cpl_propertylist * primary_hdr;
	int nbrow_met, row, nbrow_sc;
	int nv;
	cpl_table * fddl_met_mean;

	gravi_msg_function_start(1);
	cpl_ensure_code (oi_vis,       CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (p2vmred_data,  CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (preproc_data, CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (thread>=0,    CPL_ERROR_ILLEGAL_INPUT);

	/*
	 * Construction of the table containing columns of the mean FT_POS,
	 * SC_POS and the phase FC of each dispersion file
	 */

	fddl_met_mean = cpl_table_new (1);
	cpl_table_new_column_array (fddl_met_mean, "FT_POS", CPL_TYPE_DOUBLE, 4);
	cpl_table_new_column_array (fddl_met_mean, "SC_POS", CPL_TYPE_DOUBLE, 4);
	cpl_table_new_column_array (fddl_met_mean, "PHASE_FC", CPL_TYPE_DOUBLE, 4);
	cpl_table_new_column (fddl_met_mean, "TIME", CPL_TYPE_INT);

	/* Get the data */
	visMet_data = gravi_data_get_table (p2vmred_data, GRAVI_OI_VIS_MET_EXT);

	primary_hdr = gravi_data_get_plist (p2vmred_data, GRAVI_PRIMARY_HDR_EXT);
	spectrum_sc = gravi_data_get_table (preproc_data, GRAVI_SPECTRUM_DATA_SC_EXT);
	fddl_data =  gravi_data_get_table (p2vmred_data, GRAVI_FDDL_EXT);

	CPLCHECK_MSG ("ERROR1");

	nbrow_met = cpl_table_get_nrow (visMet_data);
	nbrow_sc = cpl_table_get_nrow (spectrum_sc);

	CPLCHECK_MSG ("ERROR2");

	int * time = cpl_table_get_data_int (visMet_data, "TIME");
	int exptime_sc = cpl_table_get_int (spectrum_sc, "TIME", 1, &nv) -
			cpl_table_get_int (spectrum_sc, "TIME", 0, &nv);
	cpl_array ** met_arrays = cpl_table_get_data_array (visMet_data, "PHASE_FC");
	cpl_array * met_mean = cpl_array_new (4, CPL_TYPE_DOUBLE);
	cpl_array_fill_window (met_mean, 0, 4, 0.0);
	CPLCHECK_MSG ("ERROR3");

	/*
	 * Compute the mean of the metrology included between the SC DIT
	 */

	int comp = 0;
	for (row = 0; row < nbrow_met; row ++){

		if ((time[row] > cpl_table_get_int (spectrum_sc, "TIME", 0, &nv) - exptime_sc/2) &&
				(time[row] < cpl_table_get_int (spectrum_sc, "TIME", nbrow_sc-1, &nv) + exptime_sc/2)){
			cpl_array_add (met_mean, met_arrays[row]);
			comp ++;
		}
		CPLCHECK_MSG ("ERROR4");
	}

	if (comp != 0)
		cpl_array_divide_scalar (met_mean, comp);

	/*
	 * Save the mean metrology
	 */

	cpl_table_set_array (fddl_met_mean, "PHASE_FC", 0, met_mean);
	cpl_table_set_int (fddl_met_mean, "TIME", 0,
			cpl_table_get_int (spectrum_sc, "TIME", 0, &nv) );
	cpl_array_delete (met_mean);

	int nbrow_fddl = cpl_table_get_nrow (fddl_data);

	/*
	 * Compute the mean FT po and SC pos
	 */

	cpl_array ** fddl_ft = cpl_table_get_data_array (fddl_data, "FT_POS");
	cpl_array ** fddl_sc = cpl_table_get_data_array (fddl_data, "SC_POS");
	time = cpl_table_get_data_int (fddl_data, "TIME");
	cpl_array * ft_pos = cpl_array_new (4, CPL_TYPE_DOUBLE);
	cpl_array * sc_pos = cpl_array_new (4, CPL_TYPE_DOUBLE);
	cpl_array_fill_window (ft_pos, 0, 4, 0.0);
	cpl_array_fill_window (sc_pos, 0, 4, 0.0);
	comp = 0;

	for (row = 0; row < nbrow_fddl; row ++){

		if ((time[row] > cpl_table_get_int (spectrum_sc, "TIME", 0, &nv) - exptime_sc/2) &&
				(time[row] < cpl_table_get_int (spectrum_sc, "TIME", nbrow_sc-1, &nv) + exptime_sc/2)){
			cpl_array_add (ft_pos, fddl_ft[row]);
			cpl_array_add (sc_pos, fddl_sc[row]);
			comp ++;
		}

		CPLCHECK_MSG ("Problem during the compute of the mean of the fddl");
	}

	if (comp != 0) {
		cpl_array_divide_scalar (ft_pos, comp);
		cpl_array_divide_scalar (sc_pos, comp);
	}

	/*
	 * Save the mean FT_POS and SC_POS
	 */
	cpl_table_set_array (fddl_met_mean, "FT_POS", 0, ft_pos);
	cpl_table_set_array (fddl_met_mean, "SC_POS", 0, sc_pos);
	cpl_array_delete (ft_pos);
	cpl_array_delete (sc_pos);

	if (thread == 0){
		cpl_propertylist * plist_img =  cpl_propertylist_duplicate (
				gravi_data_get_plist (p2vmred_data, GRAVI_FDDL_EXT));
    	cpl_propertylist_set_string (plist_img, "EXTNAME", "FDDL_MET_MEAN");
    	gravi_data_add (oi_vis, plist_img, fddl_met_mean);
    	cpl_propertylist_delete (plist_img);
	}
	else {

		cpl_table_insert (gravi_data_get_table (oi_vis, "FDDL_MET_MEAN"), fddl_met_mean, thread);
		cpl_table_delete (fddl_met_mean);
	}
	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}


        /* TEST: Make this cropped profile flux conservative if complex */
        //cpl_msg_info ("TEST:","profile is forced flux conservative don't work for boxcard!!");
        //cpl_msg_info ("TEST:","nxc = %lld nyc = %lld", nxc, nyc);
        //int null_val;
        //for (cpl_size x = 0; x < nxc; x++ ) {
        //    double sum_y  = 0.0;
        //    double sum_yy = 0.0;
        //    for (cpl_size y = 0; y < nyc; y++ ) {
        //        sum_y  += cpl_image_get (profile_crop, x+1, y+1, &null_val);
        //        sum_yy += pow (cpl_image_get (profile_crop, x+1, y+1, &null_val), 2.0);
        //    }
        //    for (cpl_size y = 0; y < nyc; y++ ) {
        //        double current = cpl_image_get (profile_crop, x+1, y+1, &null_val);
        //        cpl_image_set (profile_crop, x+1, y+1, current * sum_y / sum_yy);
        //    }
        //}    

cpl_error_code gravi_p2vm_mean_spectrum (gravi_data * p2vm_map, gravi_data * preproc_data)
{
	gravi_msg_function_start(1);
    cpl_ensure_code (p2vm_map,     CPL_ERROR_NULL_INPUT);
    cpl_ensure_code (preproc_data, CPL_ERROR_NULL_INPUT);

    /* Extract the spectrum of all region for SC only */
	cpl_table * spectrum_table = gravi_data_get_spectrum_data (preproc_data, GRAVI_SC);
    cpl_ensure_code (spectrum_table, CPL_ERROR_ILLEGAL_INPUT);
    
    /* Get sizes */
    cpl_size nrow = cpl_table_get_nrow (spectrum_table);
    cpl_size nreg = gravi_spectrum_get_nregion (spectrum_table);

    /* Build output table */
    cpl_table * output_table = cpl_table_extract (spectrum_table, 0, 1);

    /* Loop on region and rows to integrate the spectrums */
    for (cpl_size reg = 0; reg < nreg ; reg++) {
        cpl_array ** arrays = cpl_table_get_data_array (spectrum_table, GRAVI_DATA[reg]);
        cpl_array *  array  = cpl_table_get_data_array (output_table, GRAVI_DATA[reg])[0];
        cpl_ensure_code (arrays, CPL_ERROR_ILLEGAL_INPUT);
        for (cpl_size row = 1; row < nrow ; row++) cpl_array_add (array, arrays[row]);
        cpl_table_erase_column (output_table, GRAVI_DATAERR[reg]);
    }

    /* Add this table to the P2VM */
	gravi_data_add_table (p2vm_map, NULL, GRAVI_SPECTRUM_DATA_SC_EXT, output_table);

	gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief   Compute the dispersion calibration map (DISP)
 * FIXME: to be done
 */
/*----------------------------------------------------------------------------*/

gravi_data * gravi_compute_disp  (gravi_data * vis_data)
{
	gravi_data * disp_data;
	int ntel = 4, nbase = 6, npol = 2;
	FILE * file;
	gsl_matrix * U, * V;
	gsl_vector * S, * work;

    gravi_msg_function_start(1);
	cpl_ensure (vis_data, CPL_ERROR_NULL_INPUT, NULL);

	/*
	 * Init useful matrix
	 */
	
	/* M_matrix is to go from OPL to OPD */
	double M_tab[24]={	1., -1., 0.0, 0.0,
						1., 0.0, -1., 0.0,
						1., 0.0, 0.0, -1.,
						0.0, 1., -1., 0.0,
						0.0, 1., 0.0, -1.,
						0.0, 0.0, 1., -1.};
	gsl_matrix * M_matrix=gsl_matrix_alloc(6,4);
	memcpy(gsl_matrix_ptr(M_matrix, 0, 0), M_tab, 24*sizeof(double));

	/* M_matrix2 is to go from OPL of FDDL FT and SC to OPD */
	gsl_matrix * M_matrix2=gsl_matrix_alloc(6,8);
	gsl_vector * M_vector=gsl_vector_alloc(6);
	for (int tel=0; tel<ntel; tel++){
		gsl_matrix_get_col(M_vector, M_matrix, tel);
		gsl_matrix_set_col(M_matrix2, tel, M_vector);
		gsl_vector_scale(M_vector, -1.);
		gsl_matrix_set_col(M_matrix2, tel+ntel, M_vector);
	}
	gsl_vector_free(M_vector);

	if (INFO_DEBUG) {
		file = fopen("M_matrix2.txt","w");
		my_gsl_matrix_fprintf (file, M_matrix2, "%g");
		fclose(file);
	}

	/* Compute inverse of M_matrix2	 */

	U = gsl_matrix_alloc (8, 6);
	V = gsl_matrix_alloc (6, 6);
	S = gsl_vector_alloc (6);
	work = gsl_vector_alloc (6);

	gsl_matrix_transpose_memcpy(U, M_matrix2);
	gsl_linalg_SV_decomp (U, V, S, work);

	/* Get inverse of the M matrix */

	double wv_at;
    gsl_matrix * M_matrix2_inv = gsl_matrix_alloc(ntel*2, nbase);
	for (int j = 0; j < ntel*2; j++) {
		for (int i = 0; i < nbase; i++){
			wv_at = 0;
			for (int ii = 0; ii < nbase; ii++){
				if( gsl_vector_get(S, ii) > 1e-14)
				  wv_at += gsl_matrix_get(V, i, ii) / gsl_vector_get(S, ii) *
					gsl_matrix_get(U, j, ii);
			}
			gsl_matrix_set(M_matrix2_inv, j, i, wv_at);
		}
	}

	gsl_matrix_free (V);
	gsl_matrix_free (U);
	gsl_vector_free (S);
	gsl_vector_free (work);

	CPLCHECK_NUL("Cannot inverse matrix M");

	if (INFO_DEBUG) {
		file = fopen("M_matrix2_inv.txt","w");
		my_gsl_matrix_fprintf (file, M_matrix2_inv, "%g");
		fclose(file);
	}

	/* Check if there is 2 polarization */
	npol = 2;
	if (!gravi_data_get_oi_vis_plist (vis_data, GRAVI_SC, 0, npol) ||
		!gravi_data_get_oi_vis_plist (vis_data, GRAVI_SC, 1, npol) ||
		!gravi_data_get_oi_vis_plist (vis_data, GRAVI_FT, 0, npol) ||
		!gravi_data_get_oi_vis_plist (vis_data, GRAVI_FT, 1, npol) ) {
	  
	  cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
							 "Missing OI_VIS (need split pol for SC and FT)");
	  return NULL;
	}


	/*  (1)
	 *
	 *  Compute the position from FDDL and MET
	 */

    
	cpl_msg_info (cpl_func, "*** 1 ) Compute the real FDDL position from their positions and Metrology ***");
	cpl_msg_info (cpl_func, "Load the FDDL table");

	cpl_table * oi_flux;
	oi_flux = gravi_data_get_oi_flux (vis_data, GRAVI_SC, 0, npol);
	cpl_size nrow = cpl_table_get_nrow (oi_flux) / 4;
	CPLCHECK_NUL ("Cannot get data");

    /* fddl  = [SC_POS,FT_POS] 
     * fddl2 = [SC_POS^2,FT_POS^2] */
	gsl_matrix * fddl  = gsl_matrix_alloc (nrow, ntel*2);
	gsl_matrix * fddl2 = gsl_matrix_alloc (nrow, ntel*2);
	
	for (cpl_size row=0; row<nrow; row++) {
	  for (int tel = 0; tel<ntel; tel++) {
		double value;
		value = cpl_table_get (oi_flux, "SC_POS", row*ntel + tel, NULL);
		gsl_matrix_set (fddl, row,tel,value);
		gsl_matrix_set (fddl2,row,tel,value*value);
		
		value = cpl_table_get (oi_flux, "FT_POS", row*ntel + tel, NULL);
		gsl_matrix_set (fddl, row,ntel+tel,value);
		gsl_matrix_set (fddl2,row,ntel+tel,value*value);
	  }
	}

    /* Rescale fddl */
	gsl_matrix_scale (fddl, 1e-3);
	gsl_matrix_scale (fddl, 1e-6);

	if (INFO_DEBUG) {
		file = fopen("FDDL_matrix.txt","w");
		my_gsl_matrix_fprintf (file, fddl, "%g");
		fclose(file);
	}

	/* MET : allocate and fill the met_data matrix */
	cpl_msg_info (cpl_func, "Load the MET table");

	gsl_matrix * met_data = gsl_matrix_alloc (nrow, 4);

	for (cpl_size row=0; row<nrow; row++) {
	  for (int tel = 0; tel<ntel; tel++) {
		double value = cpl_table_get (oi_flux, "OPD_MET_FC", row*ntel + tel, NULL);
		gsl_matrix_set (met_data,row,tel,value);
	  }
	}

	/* Convert met in [um] */
	gsl_matrix_scale (met_data, 1e6);

	if (INFO_DEBUG) {
		file = fopen("MET_data.txt","w");
		my_gsl_matrix_fprintf (file, met_data, "%g");
		fclose(file);
	}

	/* Compute MET = (met_data*M_matrix.T).T  */
	
	gsl_matrix * MET = gsl_matrix_alloc(nrow, nbase);
	gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1., met_data, M_matrix, 0, MET);

	if (INFO_DEBUG) {
		file = fopen("MET.txt","w");
		my_gsl_matrix_fprintf (file, MET, "%g");
		fclose(file);
	}

	/* Model FDDL=MET+Ki*FT_POSi+Kj*SC_POSj+Kk*FT_POSk^2+Kl*SC_POSl^2 */
	cpl_msg_info(cpl_func, "Compute the FDDL linearity factors");
	int Nmodel = 8*2+6;
	gsl_matrix * model = gsl_matrix_calloc( nrow*nbase, Nmodel);
	for (cpl_size row=0; row<nrow; row++){
		for (int base=0; base<nbase; base++){
			gsl_matrix_set(model, row*nbase+base, base, 1);
			for (int tel=0; tel<ntel*2; tel++) {
				gsl_matrix_set(model, row*nbase+base, nbase+tel,
						gsl_matrix_get(fddl, row, tel) * gsl_matrix_get(M_matrix2, base, tel));
				gsl_matrix_set(model, row*nbase+base, nbase+8+tel,
						gsl_matrix_get(fddl2, row, tel) * gsl_matrix_get(M_matrix2, base, tel));
			}
		}
	}

	if (INFO_DEBUG) {
		file = fopen("model.txt","w");
		my_gsl_matrix_fprintf (file, model, "%g");
		fclose(file);
	}

	/* invert model by SV decomp and solve */
	U = gsl_matrix_alloc (nrow*nbase, Nmodel);
	V = gsl_matrix_alloc (Nmodel, Nmodel);
	S = gsl_vector_alloc (Nmodel);
	work = gsl_vector_alloc (Nmodel);
	gsl_matrix_memcpy(U, model);

	gsl_linalg_SV_decomp (U, V, S, work);

	gsl_vector * K_coeff = gsl_vector_alloc (Nmodel);
	gsl_vector * MET_vector = gsl_vector_alloc(nrow*nbase);
	memcpy (MET_vector->data, MET->data, nbase*nrow*sizeof(double));
	gsl_linalg_SV_solve (U, V, S, MET_vector, K_coeff);

	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_vector_free(MET_vector);

	cpl_msg_info(cpl_func, "Compute the real FDDL positions");
	gsl_matrix * FDDL_matrix=gsl_matrix_alloc(nrow, 8);
	for (cpl_size row=0; row<nrow; row++){
		for (int tel=0; tel<8; tel++){
			gsl_matrix_set(FDDL_matrix, row, tel,
					gsl_vector_get(K_coeff, tel+6)*gsl_matrix_get(fddl, row, tel) +
					gsl_vector_get(K_coeff, tel+6+8)*gsl_matrix_get(fddl2, row, tel));
		}
	}

	/* Correction of met offset	 */
	for (cpl_size row=0; row<nrow; row++){
		for (int base=0; base<nbase; base++){
			gsl_matrix_set(MET, row, base, gsl_matrix_get(MET, row, base)
					- gsl_matrix_get(model, row*nbase+base, base)
					* gsl_vector_get(K_coeff, base));
			//MET-=dot(MODEL[:,:6,:6],K[:6]);
		}
	}

	/* Correction of FDDL
	 * FDDL+=dot(pinv2(M_matrix2),(MET.T-dot(M_matrix2,FDDL.T))).T
	 * or FDDL+=dot(pinv2(M_matrix2),(MET-dot(FDDL,M_matrix2.T)).T).T
	 * or FDDL+=dot((MET-dot(FDDL,M_matrix2.T)), pinv2(M_matrix2).T)
	 * */
	gsl_matrix * FDDLdotM = gsl_matrix_alloc(nrow, nbase);
	gsl_matrix * METsub = gsl_matrix_alloc(nrow, nbase);
	gsl_matrix * FDDLcorrection = gsl_matrix_alloc(nrow, ntel*2);

	memcpy (METsub->data, MET->data, nbase*nrow*sizeof(double));
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., FDDL_matrix, M_matrix2, 0, FDDLdotM);
	gsl_matrix_sub(METsub, FDDLdotM);

	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1., METsub, M_matrix2_inv, 0, FDDLcorrection);
	gsl_matrix_add(FDDL_matrix, FDDLcorrection);

	gsl_matrix_free(FDDLdotM);
	gsl_matrix_free(METsub);
	gsl_matrix_free(FDDLcorrection);


	if (INFO_DEBUG) {
		file = fopen("MET.txt","w");
		my_gsl_matrix_fprintf (file, MET, "%g");
		fclose(file);
		file = fopen("coeff.txt","w");
		for (int Imodel=0; Imodel<Nmodel;Imodel++) fprintf(file, "%e \n", gsl_vector_get(K_coeff, Imodel));
		fclose(file);
		file = fopen("FDDL_matrix.txt","w");
		my_gsl_matrix_fprintf (file, FDDL_matrix, "%g");
		fclose(file);
	}


	/*  (2)
	 *
	 *  Compute the phase
	 */
	/* Load wavelength */
	cpl_msg_info(cpl_func, "*** 2 ) Compute the Phases  ***");
	cpl_msg_info(cpl_func, "Load the phase SC (VISDATA)");
	
	// FIXME: consider same wavelength for both polar
	cpl_table * oi_wavelength = gravi_data_get_oi_wave (vis_data, GRAVI_SC, 0, npol);
	
	cpl_size nwave = cpl_table_get_nrow(oi_wavelength);
	cpl_array * wavenumber = cpl_array_new(nwave, CPL_TYPE_DOUBLE);
	
	/* Compute wavelength in fiber -- FIXME: use the QC parameters to get back */
	for (cpl_size wave=0; wave<nwave; wave++)
		cpl_array_set(wavenumber, wave, (2.*M_PI/(((double)cpl_table_get(oi_wavelength, "EFF_WAVE", wave, NULL)*1.e6-0.021*1.908)/(1.-0.021))));
	if (INFO_DEBUG) {
		file = fopen("wavenumber.txt","w");
		for (cpl_size i=0; i<nwave;i++) fprintf(file, "%5.10g \n", cpl_array_get(wavenumber, i, NULL));
		fclose(file);
	}

	/* load SC phase */
	cpl_table * oi_vis;
	oi_vis = gravi_data_get_oi_vis (vis_data, GRAVI_SC, 0, npol);
	cpl_array ** visdata_p1_orig = cpl_table_get_data_array(oi_vis, "VISDATA");
	cpl_array ** visdata_p1 = cpl_malloc(cpl_table_get_nrow(oi_vis)*sizeof(cpl_array *));
	
	oi_vis = gravi_data_get_oi_vis (vis_data, GRAVI_SC, 1, npol);
	cpl_array ** visdata_p2_orig = cpl_table_get_data_array(oi_vis, "VISDATA");
	cpl_array ** visdata_p2 = cpl_malloc(cpl_table_get_nrow(oi_vis)*sizeof(cpl_array *));
	
	if (INFO_DEBUG) {
		file = fopen("vis_p1.txt","w");
	}

	/* correct phase from metrology
	 * vis_p1 = Vis_p1 * e i(MET*wavenumber)*/
	cpl_msg_info(cpl_func, "Add phase MET to phase SC");
	for (cpl_size row=0; row<nrow; row++){
		for (int base=0; base < nbase; base++){
			visdata_p1[row*nbase+base]=cpl_array_duplicate(visdata_p1_orig[row*nbase+base]);
			gravi_array_multiply_phasor(visdata_p1[row*nbase+base], I*(gsl_matrix_get(MET, row, base)), wavenumber);
			if (INFO_DEBUG) {
				for (cpl_size i=0; i<nwave;i++) fprintf(file, "%f + i%f \t", creal(cpl_array_get_complex(visdata_p1[row*nbase+base], i, NULL)),
						cimag(cpl_array_get_complex(visdata_p1[row*nbase+base], i, NULL)));
			}

			visdata_p2[row*nbase+base]=cpl_array_duplicate(visdata_p2_orig[row*nbase+base]);
			gravi_array_multiply_phasor(visdata_p2[row*nbase+base], I*gsl_matrix_get(MET, row, base), wavenumber);

			if (INFO_DEBUG) {
				fprintf(file, "\n");
			}
		}
	}
	if (INFO_DEBUG) {
		fclose(file);
	}

	/* Compute Group Delay */
	cpl_msg_info(cpl_func, "Compute Group Delay");
	cpl_array * GD1 = cpl_array_new(nrow*nbase, CPL_TYPE_DOUBLE);
	cpl_array * GD2 = cpl_array_new(nrow*nbase, CPL_TYPE_DOUBLE);
	cpl_array * GD1_stat = cpl_array_new(nrow, CPL_TYPE_DOUBLE);
	cpl_array * GD2_stat = cpl_array_new(nrow, CPL_TYPE_DOUBLE);
	cpl_array * GD1_med = cpl_array_new(nbase, CPL_TYPE_DOUBLE);
	cpl_array * GD2_med = cpl_array_new(nbase, CPL_TYPE_DOUBLE);
	double std_gd1 = 0.;
	double std_gd2 = 0.;
	for (int base=0; base < nbase; base++){
		for (cpl_size row=0; row<nrow; row++){
			double complex tmp = 0.0 * I + 0.0;
			for (cpl_size wave=1; wave<nwave; wave++)
			  tmp += cpl_array_get_complex (visdata_p1[row*nbase+base], wave, NULL) * conj(cpl_array_get_complex (visdata_p1[row*nbase+base], wave-1, NULL));
			cpl_array_set(GD1 , row*nbase + base, carg(tmp));
			cpl_array_set(GD1_stat , row, carg(tmp));
			tmp = 0.0 * I + 0.0;
			for (cpl_size wave=1; wave<nwave; wave++)
			  tmp += cpl_array_get_complex (visdata_p2[row*nbase+base], wave, NULL) * conj(cpl_array_get_complex (visdata_p2[row*nbase+base], wave-1, NULL));
			cpl_array_set(GD2 , row*nbase + base, carg(tmp));
			printf("%g ",carg(tmp));
			cpl_array_set(GD2_stat , row, carg(tmp));
		}
		// FIXME: make sure this is /nrow and not /row
		std_gd1+=cpl_array_get_stdev(GD1_stat)/nrow;
		std_gd2+=cpl_array_get_stdev(GD2_stat)/nrow;
		cpl_array_set(GD1_med, base, cpl_array_get_median(GD1_stat));
		cpl_array_set(GD2_med, base, cpl_array_get_median(GD2_stat));
	}
	cpl_array_delete(GD1_stat);
	cpl_array_delete(GD2_stat);

	cpl_msg_info(cpl_func, "GD rms average for pola 1  : %g[radians/element]", std_gd1);
	cpl_msg_info(cpl_func, "GD rms average for pola 2  : %g[radians/element]", std_gd2);

	if (INFO_DEBUG) {
		file = fopen("GD1.txt","w");
		for (cpl_size row=0; row<nrow; row++){
			for (int base=0; base < nbase; base++){
				fprintf(file, "%g \t", cpl_array_get(GD1 , row*nbase + base, NULL));
			}
			fprintf(file, "\n");
		}
		fclose(file);
		file = fopen("GD2.txt","w");
		for (cpl_size row=0; row<nrow; row++){
			for (int base=0; base < nbase; base++){
				fprintf(file, "%g \t", cpl_array_get(GD1 , row*nbase + base, NULL));
			}
			fprintf(file, "\n");
		}
		fclose(file);
	}


	/* slope as function of metrology */
	cpl_msg_info(cpl_func, "Fit the metrology slopes");
	size_t * sort_met[6];
	for (int base=0; base<nbase; base++){
		sort_met[base]=cpl_malloc(nrow*sizeof(size_t));
		gsl_sort_index (sort_met[base], gsl_matrix_ptr(MET, 0,base), nbase, nrow);
	}

	int nspace=1000;
	double step = 0.000001;
	cpl_array * amp_sum=cpl_array_new(nspace, CPL_TYPE_DOUBLE_COMPLEX);
	cpl_array_fill_window_complex(amp_sum, 0, nspace, 0.0 * I + 0.0);
	cpl_array * amp=cpl_array_new(nspace, CPL_TYPE_DOUBLE_COMPLEX);
	gsl_matrix * A1_met = gsl_matrix_alloc (nbase, nwave);
	gsl_matrix * A2_met = gsl_matrix_alloc (nbase, nwave);
	cpl_array * i2step=cpl_array_new(nspace, CPL_TYPE_DOUBLE);
	cpl_size max_pos;
	CPLCHECK_NUL("Error before fitting the metrology slopes");
	for (int i = 0; i < nspace; ++i) cpl_array_set(i2step, i, step*(double)i*2.);
	CPLCHECK_NUL("Error init istep");

	/* POLA 1 */
	for (int base = 0; base < nbase; ++base) {
		for (cpl_size wave = 0; wave < nwave; ++wave) {
			/* for all parameter space compute amp=vis_data*e^i(param*MET) and sum on row*/
			amp_sum=cpl_array_new(nspace, CPL_TYPE_DOUBLE_COMPLEX);
			cpl_array_fill_window_complex(amp_sum, 0, nspace, 0.0 * I + 0.0);
			CPLCHECK_NUL("Error before loop row");
			for (cpl_size row = 0; row < nrow; ++row) {
				cpl_array_fill_window_complex(amp, 0, nspace, (_Complex double)cpl_array_get_float_complex(visdata_p1[sort_met[base][row]*nbase+base], wave, NULL));
				gravi_array_multiply_phasor(amp, I*gsl_matrix_get(MET,  (sort_met[base])[row], base), i2step);
				cpl_array_add(amp_sum, amp);
				CPLCHECK_NUL("Error in loop row");
			}
			cpl_array_abs(amp_sum);
			/* get the best param */
			cpl_array_get_maxpos(amp_sum, &max_pos);
			gsl_matrix_set(A1_met, base, wave, max_pos*step);
			cpl_msg_info_overwritable(cpl_func, "Fit phase sc for polar 1/2 of base %d/6 (wave %lld/%lld)", base+1, wave+1, nwave);
			CPLCHECK_NUL("Error when fitting the metrology slopes");
		}
	}

	/* POLA 2 */
	for (int base = 0; base < nbase; ++base) {
		for (cpl_size wave = 0; wave < nwave; ++wave) {
			/* for all parameter space compute amp=vis_data*e^i(param*MET) and sum on row*/
			amp_sum=cpl_array_new(nspace, CPL_TYPE_DOUBLE_COMPLEX);
			cpl_array_fill_window_complex(amp_sum, 0, nspace, 0.0 * I + 0.0);
			CPLCHECK_NUL("Error before loop row");
			for (cpl_size row = 0; row < nrow; ++row) {
				cpl_array_fill_window_complex(amp, 0, nspace, (_Complex double)cpl_array_get_float_complex(visdata_p2[sort_met[base][row]*nbase+base], wave, NULL));
				gravi_array_multiply_phasor(amp, I*gsl_matrix_get(MET,  (sort_met[base])[row], base), i2step);
				cpl_array_add(amp_sum, amp);
				CPLCHECK_NUL("Error in loop row");
			}
			cpl_array_abs(amp_sum);
			/* get the best param */
			cpl_array_get_maxpos(amp_sum, &max_pos);
			gsl_matrix_set(A2_met, base, wave, max_pos*step);
			cpl_msg_info_overwritable(cpl_func, "Fit phase sc for polar 2/2 of base %d/6 (wave %lld/%lld)", base+1, wave+1, nwave);
			CPLCHECK_NUL("Error when fitting the metrology slopes");
		}
	}

	cpl_array_delete(amp_sum);
	cpl_array_delete(amp);
	cpl_array_delete(i2step);


	if (INFO_DEBUG) {
		file = fopen("A1_met.txt","w");
		my_gsl_matrix_fprintf (file, A1_met, "%g");
		fclose(file);
		file = fopen("A2_met.txt","w");
		my_gsl_matrix_fprintf (file, A2_met, "%g");
		fclose(file);
	}

	/* Reconstruction des phases unwrappés
	*  Slopes=MET[:,:,None]*A1-arange(Nw)*median(GD1,axis=0)[:,None]
	*  SC1b=SC1*exp(1j*Slopes)
	*  K=unwrap(angle(SC1b.mean(axis=0)))
	*
	*/
	cpl_msg_info(cpl_func, "Reconstruct unwrap phase ...");
	cpl_array * SC1b_unwrap;
	cpl_array * SC2b_unwrap;
	cpl_array * slope = cpl_array_new(nwave, CPL_TYPE_DOUBLE_COMPLEX);
	for (int base=0; base < nbase; base++){
		SC1b_unwrap = cpl_array_new (nwave, CPL_TYPE_DOUBLE_COMPLEX);
		cpl_array_fill_window_complex (SC1b_unwrap, 0, nwave, (double complex)(0.0 + 0.0 * I));
		SC2b_unwrap = cpl_array_new (nwave, CPL_TYPE_DOUBLE_COMPLEX);
		cpl_array_fill_window_complex (SC2b_unwrap, 0, nwave, (double complex)(0.0 + 0.0 * I));
		for (cpl_size row=0; row<nrow; row++){
			for (cpl_size wave=0; wave<nwave; wave++){
				cpl_array_set(slope, wave,cpl_array_get(visdata_p1[row*nbase+base], wave, NULL)*
						cexp(I*(
								gsl_matrix_get(MET, row, base)
								*gsl_matrix_get(A1_met, base, wave)
								-wave*cpl_array_get(GD1_med, base, NULL))
								));
			}
			cpl_array_add(SC1b_unwrap, slope);
			for (cpl_size wave=0; wave<nwave; wave++){
				cpl_array_set(slope, wave,cpl_array_get(visdata_p2[row*nbase+base], wave, NULL)*
						cexp(I*(
								gsl_matrix_get(MET, row, base)
								*gsl_matrix_get(A2_met, base, wave)
								-wave*cpl_array_get(GD2_med, base, NULL))
								));
			}
			cpl_array_add(SC2b_unwrap, slope);
		}
		cpl_array_divide_scalar(SC1b_unwrap, nrow);
		cpl_array_divide_scalar(SC2b_unwrap, nrow);
		// unwrap arg(SC1b)
		cpl_array_arg(SC1b_unwrap);
		cpl_array_arg(SC2b_unwrap);
		gravi_array_phase_unwrap (SC1b_unwrap);
		gravi_array_phase_unwrap (SC2b_unwrap);
		for (cpl_size row=0; row<nrow; row++){
			for (cpl_size wave=0; wave<nwave; wave++){
				cpl_array_set(slope, wave,cpl_array_get(visdata_p1[row*nbase+base], wave, NULL)*
						cexp(I*(
								gsl_matrix_get(MET, row, base)
								*gsl_matrix_get(A1_met, base, wave)
								-wave*cpl_array_get(GD1_med, base, NULL))
								));
			}

		}

	}

	cpl_array_delete(SC1b_unwrap);
	cpl_array_delete(SC2b_unwrap);

	disp_data = gravi_data_new (0);
	cpl_propertylist * disp_header = gravi_data_get_header (disp_data);

	/*
	 *  Save QC Parameters
	 */
	 /* Coeff linearity fddl K */
	char * qc_name;
	for (int tel = 0; tel<4; tel++ ) {
		cpl_msg_info (cpl_func, "FDDL linearity  K FT%i = %f, K SC%i = %f   ", tel+1, gsl_vector_get(K_coeff, tel+6),
				tel+1, gsl_vector_get(K_coeff, tel+4+6));
		qc_name=cpl_sprintf("ESO QC FDDL_LIN K_FT%i", tel+1);
		cpl_propertylist_append_double (disp_header, qc_name,gsl_vector_get(K_coeff, tel+6));
		cpl_propertylist_set_comment (disp_header, qc_name, "[-] K fddl linearity factor");
		qc_name=cpl_sprintf("ESO QC FDDL_LIN K_SC%i", tel+1);
		cpl_propertylist_append_double (disp_header, qc_name, gsl_vector_get(K_coeff, tel+4+6));
		cpl_propertylist_set_comment (disp_header, qc_name, "[-] K fddl linearity factor");
	}
	 /* Coeff linearity fddl K2 */
	for (int tel = 0; tel<4; tel++ ) {
		cpl_msg_info (cpl_func, "FDDL linearity  K2 FT%i = %f, K2 SC%i = %f   ", tel,
                      gsl_vector_get (K_coeff, tel+6+8),
                      tel, gsl_vector_get(K_coeff, tel+4+6));
		qc_name=cpl_sprintf("ESO QC FDDL_LIN K2_FT%i", tel+1);
		cpl_propertylist_append_double (disp_header, qc_name,gsl_vector_get(K_coeff, tel+6+8));
		cpl_propertylist_set_comment (disp_header, qc_name, "[-] K2 fddl linearity factor");
		qc_name=cpl_sprintf("ESO QC FDDL_LIN K2_SC%i", tel+1);
		cpl_propertylist_append_double (disp_header, qc_name, gsl_vector_get(K_coeff, tel+4+6+8));
		cpl_propertylist_set_comment (disp_header, qc_name, "[-] K2 fddl linearity factor");
	}
	cpl_free(qc_name);

	/*
	 * Free memory
	 */
	gsl_matrix_free(fddl);
	gsl_matrix_free(fddl2);
	gsl_matrix_free(met_data);
	gsl_matrix_free(M_matrix);
	gsl_matrix_free(M_matrix2);
	gsl_matrix_free(MET);
	gsl_matrix_free(model);
	gsl_vector_free(K_coeff);
	gsl_matrix_free(FDDL_matrix);
	cpl_array_delete(wavenumber);
//	gsl_matrix_complex_free(vis_p1);
//	gsl_matrix_complex_free(vis_p2);
	gsl_matrix_free(A1_met);
	gsl_matrix_free(A2_met);
	cpl_array_delete(GD1_med);
	cpl_array_delete(GD2_med);

	gravi_msg_function_exit(1);
	return disp_data;
}

int my_gsl_matrix_fprintf(FILE *stream, gsl_matrix *m, const char *fmt)
{
        size_t rows=m->size1;
        size_t cols=m->size2;
        size_t row,col,ml;
        int fill;
        char buf[100];
        gsl_vector *maxlen;

        maxlen=gsl_vector_alloc(cols);
        for (col=0;col<cols;++col) {
                ml=0;
                for (row=0;row<rows;++row) {
                        sprintf(buf,fmt,gsl_matrix_get(m,row,col));
                        if (strlen(buf)>ml)
                                ml=strlen(buf);
                }
                gsl_vector_set(maxlen,col,ml);
        }

        for (row=0;row<rows;++row) {
                for (col=0;col<cols;++col) {
                        sprintf(buf,fmt,gsl_matrix_get(m,row,col));
                        fprintf(stream,"%s",buf);
                        fill=gsl_vector_get(maxlen,col)+1-strlen(buf);
                        while (--fill>=0)
                                fprintf(stream," ");
                }
                fprintf(stream,"\n");
        }
        gsl_vector_free(maxlen);
        return 0;
}



/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the opd from metrology for the FC
 * 
 * @param metrology_table The metrology table from a fits file wish all
 * 	  	  	  	  	  	  the shutters are open
 * @param mjd_obs         The Modified Julian Date when metrology stated
 * 
 * @param  p2vm_met       The output p2vm table of the metrology
 * 
 * @return opl_table      This table will contain the opl of each diode
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_opdmet_calibration (cpl_table * metrology_table,
									  double mjd_obs)
{
	gravi_msg_function_start(1);
	cpl_ensure (metrology_table, CPL_ERROR_NULL_INPUT, NULL);
	cpl_vector * vectA, * vectB, *opl_vect;
	const cpl_array * volt;
	char name[90];
	char name2[90];
	int nb_tel = 4 ;
	int n_diode = 80;

	/* Get data */
	cpl_size nb_row = cpl_table_get_nrow (metrology_table);
	CPLCHECK_NUL ("Cannot get data");

	/* Create output table */
	cpl_table * opdmet_table = cpl_table_new (nb_row);
	cpl_table_new_column (opdmet_table, "TIME", CPL_TYPE_INT);
	for (cpl_size row = 0; row < nb_row; row++){
		cpl_table_set_int (opdmet_table, "TIME", row,
			cpl_table_get_int (metrology_table, "TIME",
										row, NULL) + mjd_obs);
	}
    
	/* Loop on SC and FT */
	for (int gravi_type = 0; gravi_type < 2; gravi_type++ ){
		int comb = (gravi_type == GRAVI_SC ? 1 : 2);
		
		/* Loop on beams */
		for (int tel = 0; tel < nb_tel; tel++){

			/* load vectA and vectB from metrology */
			vectA = cpl_vector_new (nb_row);
			vectB = cpl_vector_new (nb_row);

			for (cpl_size row = 0; row < nb_row; row ++){
				volt = cpl_table_get_array ( metrology_table, "VOLT", row);
				cpl_vector_set (vectA, row, cpl_array_get (volt , n_diode - 2*(comb*nb_tel - tel), NULL));
				cpl_vector_set (vectB, row, cpl_array_get (volt , n_diode - 2*(comb*nb_tel - tel) + 1, NULL));
			}

			/* Compute phase from vectA and vectB */
			opl_vect = gravi_vectors_phase_create (vectA, vectB);
			cpl_vector_multiply_scalar (opl_vect, LAMBDA_MET/(2.0*M_PI));

            FREE (cpl_vector_delete, vectA);
            FREE (cpl_vector_delete, vectB);
			CPLCHECK_NUL ("Compute OPD");

			/* Compute OPD, OPL_FT and OPL_SC from the phase */
			if (gravi_type == GRAVI_SC) {
				sprintf (name,"OPL_FC_SC%d", tel + 1);
				cpl_table_wrap_double (opdmet_table, cpl_vector_get_data(opl_vect), name);
				cpl_table_set_column_unit (opdmet_table, name, "m");
				CPLCHECK_NUL("Wrap OPL_FC_SC");

				/* duplicate the OPL_FC_SC into the OPD_FC */
				sprintf (name, "OPD_FC%d", tel + 1);
				sprintf (name2, "OPL_FC_SC%d", tel + 1);
				cpl_table_duplicate_column (opdmet_table, name, opdmet_table, name2);
				cpl_table_set_column_unit (opdmet_table, name, "m");
				CPLCHECK_NUL("Wrap OPD_FC");
			}

			if (gravi_type == GRAVI_FT) {
				sprintf (name,"OPL_FC_FT%d", tel + 1);
				cpl_table_wrap_double (opdmet_table, cpl_vector_get_data(opl_vect), name);
				cpl_table_set_column_unit (opdmet_table, name, "m");
				
				/* compute the OPD as OPL_SC - OPL_FT */
				sprintf (name2, "OPD_FC%d", tel + 1);
				cpl_table_subtract_columns (opdmet_table, name2, name);
				cpl_table_set_column_unit (opdmet_table, name, "m");
				CPLCHECK_NUL("Wrap OPL_FC_FT");
			}
			cpl_vector_unwrap (opl_vect);
			
		} /* End loop beam */
	} /* End loop SC, FT */

	/* Verbose */
	gravi_msg_function_exit(1);
	return opdmet_table;
}


cpl_error_code gravi_phase_correct_closures_new (cpl_table * phase_table)
{
    gravi_msg_function_start(1);
	cpl_ensure_code (phase_table, CPL_ERROR_NULL_INPUT);
    
    cpl_table * input_table = cpl_table_duplicate (phase_table);
    cpl_array ** ref = cpl_table_get_data_array (input_table, "PHASE");
    cpl_array ** out = cpl_table_get_data_array (phase_table, "PHASE");
    
    for (int base = 0; base < nbase; base++) {
        for (int clo = 0; clo < 2; clo ++) {
            for (int b = 0; b < 2; b ++) {
                cpl_array * tmp = ref[ GRAVI_TRI_BASE[base][clo][b] ];
                double sign = GRAVI_TRI_SIGN[base][clo][b];
                cpl_array_multiply_scalar (tmp, sign);
                cpl_array_add (out[base], tmp);
                cpl_array_multiply_scalar (tmp, sign);
                CPLCHECK_MSG ("Cannot correct CP");
            }
        }
        cpl_array_divide_scalar (out[base], 3.0);
    }
    FREE (cpl_table_delete, input_table);
    
    gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief reform the PHASE_FC column of the  metrology for the FC into OPD_FC
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_metrology_reform (cpl_table * vis_met)
{
	gravi_msg_function_start(1);
	cpl_ensure (vis_met, CPL_ERROR_NULL_INPUT, NULL);
	char name[90];
	int ntel = 4;

	/* Get data */
	cpl_size nrow = cpl_table_get_nrow (vis_met)/ntel;
	CPLCHECK_NUL ("Cannot get data");

	/* Create output table */
	/* Create the output table for OPD_MET */
	cpl_table * opd_met = cpl_table_new (nrow);
	cpl_table_new_column (opd_met, "TIME", CPL_TYPE_INT);
	cpl_table_set_column_unit (opd_met, "TIME", "usec");

	/* Create and get the pointer of the table columns */
	double **p_opd_met;
	p_opd_met = cpl_malloc(ntel*sizeof(double*));
	for (int tel=0; tel<ntel; tel++) {
		sprintf (name, "OPD_FC%d", tel + 1);
		cpl_table_new_column (opd_met, name, CPL_TYPE_DOUBLE);
		cpl_table_set_column_unit (opd_met, name, "m");
		p_opd_met[tel]=cpl_table_get_data_double(opd_met, name);
	}
	CPLCHECK_NUL ("Cannot create data");

	/*
	 * Loop on the row to reform the table
	 */
	double * p_vis_met=cpl_table_get_data_double(vis_met, "PHASE_FC");
	for (cpl_size row = 0; row < nrow; row++){

		/* Copy the TIME column */
		cpl_table_set_int (opd_met, "TIME", row,
			cpl_table_get_int (vis_met, "TIME", row*ntel, NULL));

		/* Reform the VIS_MET column */
		for (int tel=0; tel<ntel; tel++) {
			p_opd_met[tel][row] = p_vis_met[row*ntel+tel]*LAMBDA_MET/(2.0*M_PI);
		CPLCHECK_NUL ("Cannot set data");
		}
	}
	CPLCHECK_NUL ("Cannot reform data");

	cpl_free(p_opd_met);

	/* Verbose */
	gravi_msg_function_exit(1);
	return opd_met;
}


            if (cpl_table_get (phase_ft, "PHASE_SC", row*nbase+base, NULL) !=0
                && (row % (nrow/5) == 0) ) {
                double phase_sc  = cpl_table_get (phase_ft, "PHASE_SC", row*nbase+base, NULL);
                double phasec_sc = gsl_matrix_get (A_data, base * nrow + row, 7);
                cpl_msg_info ("TEST", "FT base %i row %lld : phi_sc = %g phase_sc = %g diff = %g [rad]",
                              base, row, phasec_sc, phase_sc, phasec_sc - phase_sc);
                double phase_met  = cpl_table_get (phase_ft, "PHASE_MET_FC", row*nbase+base, NULL);
                double phasec_met = gsl_vector_get (bis_data, base * nrow + row) / LAMBDA_MET * CPL_MATH_2PI;
                cpl_msg_info ("TEST", "FT base %i row %lld : phi_met = %g phase_met = %g diff = %g [rad]",
                              base, row, phasec_met, phase_met, phasec_met - phase_met);
            }



    /* Resample MET and FT phases at the SC period */
    gravi_vis_create_met_ft (phase_ft, vis_met);
    gravi_vis_create_phasesc_ft (phase_ft, phase_sc, dit_sc);
	cpl_table_save (vis_met, NULL, NULL, "vismet.fits", CPL_IO_CREATE);
	cpl_table_save (phase_sc, NULL, NULL, "phase_sc.fits", CPL_IO_CREATE);
	cpl_table_save (phase_ft, NULL, NULL, "phase_ft.fits", CPL_IO_CREATE);



cpl_vector * gravi_phase_fit_opdmet (cpl_table * vis_met,
                                     cpl_table * phase_sc,
                                     cpl_table * phase_ft,
                                     double dit_sc)
{
	int nbase = 6, ntel = 4;
	double opl, time_sc, time_ft, time_met;
	int nv, comp;

    gravi_msg_function_start(1);
	cpl_ensure (vis_met,  CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (phase_sc, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (phase_ft, CPL_ERROR_NULL_INPUT, NULL);

	/* Get the number of acquisitions */
	cpl_size nrow     = cpl_table_get_nrow (phase_ft) / nbase;
	cpl_size nrow_sc  = cpl_table_get_nrow (phase_sc) / nbase;
	cpl_size nrow_met = cpl_table_get_nrow (vis_met) / ntel;
	CPLCHECK_NUL("Cannot get data");

	/* Get time */ 
	int * time_FT = cpl_table_get_data_int (phase_ft, "TIME");
	int * time_SC = cpl_table_get_data_int (phase_sc, "TIME");
    CPLCHECK_NUL ("Cannot get time");
    
    /* Compute the periods of the signals */
    double exptime_ft = time_FT[nbase] - time_FT[0];
    double exptime_sc = time_SC[nbase] - time_SC[0];
    
    /* Get pointer to data */
    double * ft_phase = cpl_table_get_data_double (phase_ft, "PHASE");
    double * sc_phase = cpl_table_get_data_double (phase_sc, "PHASE");
    double * phase_met = cpl_table_get_data_double (vis_met, "PHASE_FC");
    CPLCHECK_NUL ("Cannot get data");

	/* Allocate the matrices */
	gsl_matrix * A_data = gsl_matrix_calloc (nbase * nrow, 8);
	gsl_matrix * U, * A;
	gsl_vector * bis_data = gsl_vector_alloc (nbase * nrow), * bis;
	gsl_matrix * X = gsl_matrix_alloc (8, 8);
    gsl_matrix * V = gsl_matrix_alloc (8, 8);
	gsl_vector * S = gsl_vector_alloc (8);
	gsl_vector * work = gsl_vector_alloc (8);
    gsl_vector * x = gsl_vector_alloc (8);
	CPLCHECK_NUL ("Allocate matrix");
    
	/*
	 *  Extract the OPl metrology and compute the mean of each exposure time
	 */
    int ind_sc;
	for (int base = 0; base < nbase; base ++) {
        int tel0 = GRAVI_BASE_TEL[base][0];
        int tel1 = GRAVI_BASE_TEL[base][1];
        
		ind_sc = 0;
		int im = 0;

		/*
		 * resample the SC and Metrology at FT period
		 */
		for (cpl_size row = 0; row < nrow; row ++){ // loop on row FT

			gsl_matrix_set (A_data, base * nrow + row, 6, - ft_phase[row*nbase+base]);
			CPLCHECK_NUL ("Set FT phase in matrix");

			time_ft = time_FT[row*nbase+base];
			if (ind_sc < nrow_sc) time_sc = time_SC[ind_sc*nbase+base];
			CPLCHECK_NUL ("Get time sc");

			/* increase ind ft if time_ft in the range of the nex DIT */
			while(time_ft  > (time_sc+exptime_sc/2.)){
				ind_sc++;
				if (ind_sc < nrow_sc) time_sc = time_SC[ind_sc*nbase+base];
				else break;
			}
			CPLCHECK_NUL ("Get next sc dit");

			/* get phi sc */
			double phi_sc = 0;
			if (ind_sc > 0){
				if(ind_sc < nrow_sc) phi_sc = sc_phase[ind_sc*nbase+base];
				else                 phi_sc = sc_phase[(nrow_sc - 1)*nbase+base];
			}
			else phi_sc = sc_phase[0*nbase+base];
			CPLCHECK_NUL ("Get phi_sc");


			gsl_matrix_set (A_data, base * nrow + row, 7, phi_sc);
			gsl_matrix_set (A_data, base * nrow + row, base, 1.0);

			/* Metrology case */
			opl = 0;
			comp = 0;

			if (im < nrow_met)
			time_met = cpl_table_get_int (vis_met, "TIME", im*ntel, &nv);

			while ((time_met  < (time_ft + exptime_ft))){

				if (im < nrow_met) {
					opl += (phase_met[im*ntel+tel0] - phase_met[im*ntel+tel1]) * LAMBDA_MET / CPL_MATH_2PI;
					comp ++;
				}

				im++;
				if (im < nrow_met) {
					time_met = cpl_table_get_int (vis_met, "TIME", im*ntel, &nv);
				}
				else {
					break;
				}

				CPLCHECK_NUL("Get time metrology");
			}

			/* average the metrology over the FT DIT */
			if (comp != 0)
				gsl_vector_set (bis_data, base * nrow + row, opl/comp);

			/* if no metrology signal within the FT DIT interpolate the metrology */
			else {
				if (im > 0){
					if(im < nrow_met) {
                        double opd1 = (phase_met[(im-1)*ntel+tel0] - phase_met[(im-1)*ntel+tel1]) * LAMBDA_MET / CPL_MATH_2PI;
                        double opd2 = (phase_met[im*ntel+tel0]     - phase_met[im*ntel+tel1]) * LAMBDA_MET / CPL_MATH_2PI;
						double opd = (opd2 - opd1) *
                            (time_ft  - cpl_table_get_int (vis_met, "TIME", (im-1)*ntel, &nv)) /
							(time_met - cpl_table_get_int (vis_met, "TIME", (im-1)*ntel, &nv)) +
							opd1;
						gsl_vector_set (bis_data, base * nrow + row, opd);
					}
					else {
                        double opd = (phase_met[(nrow_met-1)*ntel+tel0] - phase_met[(nrow_met-1)*ntel+tel1]) * LAMBDA_MET / CPL_MATH_2PI;
						gsl_vector_set (bis_data, base * nrow + row, opd);
					}
				}
				else {
                    double opd = (phase_met[0*ntel+tel0] - phase_met[0*ntel+tel1]) * LAMBDA_MET / CPL_MATH_2PI;
                    gsl_vector_set (bis_data, base * nrow + row, opd);
				}
				CPLCHECK_NUL("Interpolate the metrology");
			}
                        
		} /* end loop on row FT */
	} /* End loop on base */

	/*
	 *  filter the data out of the SC integration time
	 */

	long n_row_A=nrow_sc*(int)(dit_sc/exptime_ft+1);
	cpl_vector *i_A_vector = cpl_vector_new(n_row_A);
	int i_A=0;
	int time_mod;
	int t0_sc   = time_SC[0*nbase];
	int tend_sc = time_SC[(nrow_sc-1)*nbase];

	/* Find the index of the frames within the SC integration 
     *Assume all baseline have the same timing */
	for (cpl_size row=0; row<nrow; row++){
		time_ft = time_FT[row*nbase];
		if ((time_ft >= t0_sc-dit_sc/2) && (time_ft < tend_sc+dit_sc/2)){
			time_mod=((int)(time_ft-(t0_sc-dit_sc/2)) % (int)exptime_sc);
			if ( time_mod >= 0 && time_mod < dit_sc ) {
				cpl_vector_set (i_A_vector, i_A, row);
				i_A++;
			}
		}
	}
	n_row_A=i_A;

	/* copy the frames within SC integration in the matrix A */
	A = gsl_matrix_alloc (nbase * n_row_A, 8);
	bis = gsl_vector_alloc (nbase * n_row_A);
	cpl_vector *time_A = cpl_vector_new(n_row_A);

	for (int base = 0; base < nbase; base ++)
		for (cpl_size i_A = 0; i_A < n_row_A; i_A++){
			cpl_size row = cpl_vector_get(i_A_vector, i_A);
			for (int col = 0; col < 8; col ++)
				gsl_matrix_set (A, base * n_row_A + i_A, col, gsl_matrix_get (A_data, base * nrow + row, col));

			gsl_vector_set (bis, base * n_row_A + i_A, gsl_vector_get (bis_data, base * nrow + row));
			if (base == 0)
				cpl_vector_set (time_A, i_A, time_FT[row*nbase+base]);
		}

	cpl_vector_delete(i_A_vector);
   
	/*
	 * Solve the linear equation Ax = b to get the coefficients to deduce the OPDs
	 */
	nrow=n_row_A;
	U = gsl_matrix_alloc (nbase * nrow, 8);
	gsl_matrix_memcpy (U, A);


	gsl_linalg_SV_decomp (U, V, S, work);
	gsl_linalg_SV_solve (U, V, S, bis, x);
	cpl_vector * opd_coeff = cpl_vector_new (3);
	cpl_vector_set (opd_coeff, 0, gsl_vector_get (x, 7));
	cpl_vector_set (opd_coeff, 1, gsl_vector_get (x, 6));

	cpl_msg_debug (cpl_func, "Wavelength => SC = %g, FT = %g\n",
                   cpl_vector_get(opd_coeff,0)*2*M_PI, cpl_vector_get(opd_coeff,1)*2*M_PI);
    

	/*
	 *  Save the A and B matrix for debug
	 */
	if (INFO_DEBUG) {
	  cpl_table * table_output = cpl_table_new (nrow * nbase);
	  const char * table_name[9] = {"A0","A1","A2","A3","A4","A5","A6","A7","B"};
	  
	  /* Init memory */
	  double ** table_value = cpl_malloc (9 * sizeof(double*));
	  for (int col = 0; col < 9; col++)
		table_value[col] = cpl_malloc (nrow * nbase * sizeof(double));
	  
	  /* Fill values */
	  cpl_msg_info (cpl_func, "fill tmp table");
	  for (int base = 0; base < nbase; base ++)
		for (cpl_size row = 0; row < nrow; row++) {
		  for (int col = 0; col < 8; col ++)
              table_value[col][base * nrow + row] = gsl_matrix_get (A, base * nrow + row, col);
          
		  table_value[8][base * nrow + row] = gsl_vector_get (bis, base * nrow + row);
		}

	  /* Wrap into table */
	  cpl_msg_info (cpl_func, "wrap tmp table");
	  for (int col = 0; col < 9; col++) {
		cpl_table_wrap_double (table_output, table_value[col], table_name[col]);
	  }
	  
	  /* Add best fit in header */
	  cpl_propertylist * table_header = cpl_propertylist_new ();
	  for (int col = 0; col < 8; col ++)
		cpl_propertylist_update_double (table_header, table_name[col], gsl_vector_get (x, col));
	  
	  /* Save tmp table */
	  cpl_msg_info (cpl_func, "save tmp table");
	  cpl_table_save (table_output, NULL, table_header, "matrix_AB.fits",  CPL_IO_CREATE);
	  
	  cpl_msg_info (cpl_func, "done");
	}
	/* --- */
	  
	

	cpl_vector *vect_phase_sc_ft=cpl_vector_new(nrow*nbase);
	cpl_vector *vect_dopd_met=cpl_vector_new(nrow*nbase);
	cpl_vector *vect_diff=cpl_vector_new(nrow*nbase);
	for (cpl_size row = 0; row < nrow*nbase; row ++){
		cpl_vector_set(vect_dopd_met, row, gsl_vector_get (bis, row));
		cpl_vector_set(vect_phase_sc_ft, row,
                       gsl_matrix_get (A, row, 0)*gsl_vector_get (x, 0)+
                       gsl_matrix_get (A, row, 1)*gsl_vector_get (x, 1)+
                       gsl_matrix_get (A, row, 2)*gsl_vector_get (x, 2)+
                       gsl_matrix_get (A, row, 3)*gsl_vector_get (x, 3)+
                       gsl_matrix_get (A, row, 4)*gsl_vector_get (x, 4)+
                       gsl_matrix_get (A, row, 5)*gsl_vector_get (x, 5)+
                       gsl_matrix_get (A, row, 6)*gsl_vector_get (x, 6)+
                       gsl_matrix_get (A, row, 7)*gsl_vector_get (x, 7));
		cpl_vector_set(vect_diff, row, cpl_vector_get(vect_dopd_met, row)-cpl_vector_get(vect_phase_sc_ft, row));
	}

	double rms_fit_phase_met = cpl_vector_get_stdev(vect_diff);
	cpl_msg_info (cpl_func, "RMS of residuals on fit of a.phi_ft+b.phi_sc+c = met :  %g", rms_fit_phase_met);
	cpl_vector_set (opd_coeff, 2, rms_fit_phase_met);
    
    cpl_msg_info ("TEST", "coeff SC = %.20g [um]", cpl_vector_get (opd_coeff, 0) * CPL_MATH_2PI * 1e6);
    cpl_msg_info ("TEST", "coeff FT = %.20g [um]", cpl_vector_get (opd_coeff, 1) * CPL_MATH_2PI * 1e6);
    cpl_msg_info ("TEST", "residual = %.20g [um]", rms_fit_phase_met  * CPL_MATH_2PI * 1e6);
	
	/*
	 * If the residuals are too high print a warning : SC and RMN data may be desynchronized
	 *  */
	if (rms_fit_phase_met > 1e-7 ){
		cpl_msg_info (cpl_func,    "*************************************************");
		cpl_msg_warning (cpl_func, "****  !!! residuals of the fit too high !!!  ****");
		cpl_msg_warning (cpl_func, "****     SC and RMN may be desynchronized    ****");
		cpl_msg_warning (cpl_func, "****     (or out of the envelope in LOW)     ****");
		cpl_msg_info (cpl_func,    "*************************************************");
	}

	if (PLOT_MET_PHASE_FIT)	{
		cpl_msg_info(cpl_func, "Plot fit residuals");

		cpl_errorstate prestate = cpl_errorstate_get();
		const cpl_vector ** vectors=malloc(3 * sizeof(cpl_vector*));
		vectors[0]=NULL;
		vectors[1]=vect_phase_sc_ft;
		vectors[2]=vect_dopd_met;

		cpl_plot_vectors (cpl_sprintf("set title 'Met fit  case a=%g  b=%g'; set xlabel 'time[10-6s]'; set ylabel 'Phase_ft+Phase_sc, dopd_met';", gsl_vector_get (x, 0), gsl_vector_get (x, 1)),
				 "", "", vectors, 3);
		cpl_plot_vector(cpl_sprintf("set title 'Met fit  case a=%g  b=%g'; set xlabel 'time[10-6s]'; set ylabel 'Phase_ft+Phase_sc-dopd_met';", gsl_vector_get (x, 0), gsl_vector_get (x, 1)),
						"", "", vect_diff);
		cpl_free(vectors);
		cpl_errorstate_set (prestate);
	}

	FREE(cpl_vector_delete, vect_phase_sc_ft);
	FREE(cpl_vector_delete, vect_dopd_met);
	FREE(cpl_vector_delete, vect_diff);
	gsl_matrix_free (A);
	gsl_matrix_free (A_data);
	gsl_matrix_free (U);
	gsl_matrix_free (V);
	gsl_vector_free (x);
	gsl_vector_free (bis);
	gsl_vector_free (bis_data);
	cpl_vector_delete(time_A);
	gsl_matrix_free (X);
	gsl_vector_free (S);
	gsl_vector_free (work);

    gravi_msg_function_exit(1);
	return opd_coeff;
}




cpl_table * gravi_vis_create_opdguess_sc (cpl_table * spectrumsc_table,
                                          cpl_table * vis_FT,
                                          cpl_table * vis_MET,
                                          double dit_sc)
{
    gravi_msg_function_start(1);
    cpl_ensure (spectrumsc_table, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (vis_FT,           CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (vis_MET,          CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (dit_sc>0,         CPL_ERROR_ILLEGAL_INPUT, NULL);

    cpl_size nbase = 6, ntel = 4;
    cpl_size nrow_sc  = cpl_table_get_nrow (spectrumsc_table);

    /* Create table */
    cpl_table * vis_SC = cpl_table_new (nrow_sc * nbase);

    /* Create the time colum */
    gravi_table_new_column (vis_SC, "TIME", "us", CPL_TYPE_INT);
    for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++)  {
        double value = cpl_table_get (spectrumsc_table, "TIME", row_sc, NULL);
        for (cpl_size base = 0; base < nbase; base++) 
            cpl_table_set (vis_SC, "TIME", row_sc*nbase+base, value);
    }
    
    /* Create synch columns if not yet existing */
    gravi_signal_create_sync (vis_SC, 6, dit_sc, vis_MET, 4, "MET");
    gravi_signal_create_sync (vis_SC, 6, dit_sc, vis_FT, 6, "FT");

    CPLCHECK_NUL ("Cannot create synch");
  
    /* Get SC data */
    int * first_met = cpl_table_get_data_int (vis_SC, "FIRST_MET");
    int * last_met  = cpl_table_get_data_int (vis_SC, "LAST_MET");
    int * first_ft  = cpl_table_get_data_int (vis_SC, "FIRST_FT");
    int * last_ft   = cpl_table_get_data_int (vis_SC, "LAST_FT");
    
    CPLCHECK_NUL ("Cannot get data");
    
    /* Get MET and data FT data */
    double * phase_met = cpl_table_get_data_double (vis_MET, "PHASE_FC");
    double * phase_ft  = cpl_table_get_data_double (vis_FT, "PHASE");
    
    CPLCHECK_NUL ("Cannot get direct pointer to data");
    
    /* New columns */
    gravi_table_new_column (vis_SC, "OPD", "m", CPL_TYPE_DOUBLE);
    double * opd_sc = cpl_table_get_data_double (vis_SC, "OPD");
    
    /* Loop on base and rows */
    for (cpl_size base = 0; base < nbase; base++) {
        for (cpl_size row_sc = 0; row_sc < nrow_sc; row_sc ++) {
            cpl_size nsc = row_sc * nbase + base;
            
            /* Sum over synch MET frames */
            double opd_met = 0.0;
            for (cpl_size row_met = first_met[nsc] ; row_met < last_met[nsc]; row_met++) {
                cpl_size nmet0 = row_met * ntel + GRAVI_BASE_TEL[base][0];
                cpl_size nmet1 = row_met * ntel + GRAVI_BASE_TEL[base][1];
                opd_met += (phase_met[nmet0] - phase_met[nmet1]) * LAMBDA_MET / CPL_MATH_2PI;
            }
            cpl_size nframe_met = last_met[nsc] - first_met[nsc];
            if (nframe_met != 0 ) opd_met /= nframe_met;
            
            /* Sum over synch FT frames */
            double opd_ft = 0.0;
            for (cpl_size row_ft = first_ft[nsc] ; row_ft < last_ft[nsc]; row_ft++) {
                cpl_size nft = row_ft * nbase + base;
                opd_ft += phase_ft[nft] * 2.2e-6 / CPL_MATH_2PI;
            }
            cpl_size nframe_ft = last_ft[nsc] - first_ft[nsc];
            if (nframe_ft != 0 ) opd_ft /= nframe_ft;
            
            /* Set the guess */
            opd_sc[nsc] = opd_met + opd_ft;
            
            CPLCHECK_NUL ("Cannot build opdguess");
        }  /* End loop on SC rows */
    } /* End loop on base */
    
    gravi_msg_function_exit(1);
    return vis_SC;
}


cpl_error_code gravi_phase_correct_closures (cpl_table * phase_table)
{
    gravi_msg_function_start(1);
	cpl_ensure_code (phase_table, CPL_ERROR_NULL_INPUT);

    int ntel = 4, nbase = 6;
	int tel_1[6] = {0,0,0,1,1,2};
	int tel_2[6] = {1,2,3,2,3,3};
	double ai;
	gsl_matrix * gsl_SC_kernel = NULL;
	gsl_vector * phase_coeff = NULL;
	gsl_matrix * M = NULL;
	gsl_matrix * U = NULL;
	gsl_matrix * V = NULL;
	gsl_vector * S = NULL;
	gsl_vector * work = NULL;
	cpl_matrix * M_inv = NULL;
	cpl_matrix * SC_kernel_mat_final = NULL;
	gsl_vector * Ks = NULL;
	cpl_matrix * kernel = NULL;
	
	/* Compute of the matrix M
	 * 		[1  -1   0   0
	 * 	M =  1   0  -1   0
	 * 		 1   0   0  -1
	 * 		 0   1  -1	 0
	 * 		 0   1   0  -1
	 * 		 0   0   1  -1]     */
	M = gsl_matrix_alloc (nbase,ntel);
	gsl_matrix_set_zero (M);

	for (int j = 0; j < nbase; j++){
		gsl_matrix_set (M, j, tel_1[j], 1);
		gsl_matrix_set (M, j, tel_2[j], -1);
	}

	CPLCHECK_CLEAN ("Cannot set matrix M");

	/*
	 * Compute of the SV dec of M
	 */
	
	U = gsl_matrix_alloc (6, 4);
	V = gsl_matrix_alloc (4, 4);
	S = gsl_vector_alloc (4);
	work = gsl_vector_alloc (4);
	
	gsl_matrix_memcpy (U, M);
	gsl_linalg_SV_decomp (U, V, S, work);

	/*
	 * Get inverse of the M matrix
	 */
	double wv_at;
    double * a_inv_data = cpl_malloc(ntel * nbase * sizeof(double));
	for (int j = 0; j < nbase; j++) {
		for (int i = 0; i < ntel; i++){
			wv_at = 0;
			for (int ii = 0; ii < ntel; ii++){
				if( gsl_vector_get(S, ii) > 1e-14)
				  wv_at += gsl_matrix_get(V, i, ii) / gsl_vector_get(S, ii) *
					gsl_matrix_get(U, j, ii);
			}
			a_inv_data[j + i * nbase] = wv_at;
		}
	}

	gsl_matrix_free (V);
	gsl_matrix_free (U);
	gsl_vector_free (S);
	gsl_vector_free (work);

	M_inv  = cpl_matrix_wrap (ntel, nbase, a_inv_data);

	CPLCHECK_CLEAN ("Cannot invers matrix M");
	
	/*
	 * Compute of the matrix Kernel of M
	 */

	kernel = cpl_matrix_new (nbase, nbase);

	for (int j = 0; j < nbase; j++) {
			for (int k = 0; k < nbase; k++) {
				ai = 0;
				for (int ii = 0; ii < ntel; ii++){
					ai += gsl_matrix_get (M, j, ii)*cpl_matrix_get (M_inv, ii, k);
				}
				if (j == k)
					cpl_matrix_set (kernel, j, k, ai - 1);
				else
					cpl_matrix_set (kernel, j, k, ai);
			}
	}
	
	cpl_matrix_unwrap (M_inv);
	cpl_free (a_inv_data);
	
	CPLCHECK_CLEAN ("Cannot compute the matrix kernel of M");


	/*
	 * Introducing the phase of each base and compute the Ks vector
	 */
	
	int nrow = cpl_table_get_nrow (phase_table) / nbase;
	double * ph_data;
    double * pphase = cpl_table_get_data_double (phase_table, "PHASE");
    cpl_ensure_code (pphase, CPL_ERROR_ILLEGAL_INPUT);

	
	Ks = gsl_vector_alloc (nbase* nrow);
	for (int k = 0; k < nbase; k++) {
		for (cpl_size row = 0;row < nrow; row++) {
			ai = 0;
			for (int ii = 0; ii < nbase; ii++){
 				ai += cpl_matrix_get (kernel, k, ii) * pphase[row*nbase+ii];
			}
			gsl_vector_set (Ks, row + k * nrow, ai);

			CPLCHECK_CLEAN ("Cannot multply phase kernel to compute Ks");
		}
	}

	/* 
	 * Construction of the big matrix T 
	 */

	double * ones_t;
	cpl_array * phase_temp;
	
	double * SC_kernel = cpl_malloc ( nbase*nbase * nrow * sizeof (double));
	double * SC_kernel2 = cpl_malloc ( nbase*nbase * nrow * sizeof (double));
	cpl_array * ones = cpl_array_new (nrow, CPL_TYPE_DOUBLE);

	for (int i = 0; i < nbase; i++) {
		for (int k = 0; k < nbase; k ++){
            
			// phase_temp = cpl_array_duplicate (phase [i]);
            phase_temp = cpl_array_new (nrow, CPL_TYPE_DOUBLE);
            for (cpl_size row = 0; row < nrow; row++)
                cpl_array_set (phase_temp, row, pphase[row*nbase+i]);
            
			cpl_array_fill_window_double (ones, 0, nrow, 1);
			cpl_array_multiply_scalar (ones, cpl_matrix_get (kernel, i, k));
			cpl_array_multiply_scalar (phase_temp, cpl_matrix_get (kernel, i, k));
			ph_data = cpl_array_get_data_double (phase_temp);
			ones_t = cpl_array_get_data_double (ones);

			memcpy (SC_kernel2 +  (i*nbase + k) * nrow , ph_data, nrow * sizeof(double));
			memcpy (SC_kernel + (i*nbase + k) * nrow, ones_t, nrow * sizeof(double));

			CPLCHECK_CLEAN ("Error in loop SC_kernel");
			cpl_array_delete (phase_temp);
		}

	}

    FREE (cpl_array_delete, ones);
	CPLCHECK_CLEAN ("Cannot compute SC_kernel");

	cpl_matrix * temp = cpl_matrix_wrap (nbase, nbase* nrow, SC_kernel);
	cpl_matrix * SC_kernel_mat = cpl_matrix_duplicate (temp);
	cpl_matrix * SC_kernel_mat2 = cpl_matrix_wrap (nbase, nbase* nrow, SC_kernel2);

	cpl_matrix_append (SC_kernel_mat, SC_kernel_mat2, 1);

	SC_kernel_mat_final = cpl_matrix_extract (SC_kernel_mat, 0, 0, 1, 1, 11, nbase*nrow);

	cpl_matrix_unwrap (temp);
	cpl_matrix_unwrap (SC_kernel_mat2);
	cpl_free (SC_kernel2);
	cpl_free (SC_kernel);
	cpl_matrix_delete (SC_kernel_mat);

	gsl_SC_kernel = gsl_matrix_alloc (nbase*nrow, 11);
	
	for (int j = 0; j < cpl_matrix_get_nrow(SC_kernel_mat_final); j++){
		for (int k = 0; k < cpl_matrix_get_ncol(SC_kernel_mat_final); k++)
			gsl_matrix_set (gsl_SC_kernel, k, j, cpl_matrix_get (SC_kernel_mat_final, j, k));
	}

	/*
	 * Resolve the linear equation T * x = Ks
	 */

	U = gsl_matrix_alloc (nbase*nrow, 11);
	V = gsl_matrix_alloc (11, 11);
	S = gsl_vector_alloc (11);
	phase_coeff = gsl_vector_alloc (11);
	work = gsl_vector_alloc (11);
	gsl_matrix_memcpy (U, gsl_SC_kernel);

	gsl_linalg_SV_decomp (U, V, S, work);
	for (int i = 0; i < 11; i++)
		if (gsl_vector_get (S, i) < 1e-10)
			gsl_vector_set (S, i , 0);
	gsl_linalg_SV_solve (U, V, S, Ks, phase_coeff);

    /* 
     * Apply the correction on 5 baselines
     */
    
	for (int base = 0; base < nbase - 1; base ++){
        double f = gsl_vector_get (phase_coeff, nbase + base);
        cpl_msg_info (cpl_func,"correction factor = 1 %+.20f", f);
        for (cpl_size row = 0; row < nrow; row++) {
            pphase[row*nbase+base] *= 1 - f;
        }
	}

	/* 
	 * Cleanup
	 */

 cleanup :
	FREE (gsl_matrix_free, M);
	FREE (gsl_vector_free, Ks);
	FREE (cpl_matrix_delete, SC_kernel_mat_final);
	FREE (cpl_matrix_delete, kernel);
	FREE (gsl_matrix_free, gsl_SC_kernel);
	FREE (gsl_matrix_free, U);
	FREE (gsl_matrix_free, V);
	FREE (gsl_vector_free, S);
	FREE (gsl_vector_free, phase_coeff);
	FREE (gsl_vector_free, work);

    gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}


double gravi_envelope (double lambda, double delta_lambda, double opd)
{
    /* Compute coherent length [m] */
	double coh_len= lambda * lambda / delta_lambda / 3;

	/* Gaussien enveloppe */
    double value = 1.0;
    if (opd != 0.0)
        value = exp (-1 * opd * opd / (coh_len*coh_len/2.) );

    return value;
}



                /* TEST */
                if (wave == nwave/2 && base == 0 & pol == 0) {
                    int step = nwave == 5 ? 100 : 1;
                    
                    cpl_vector ** array = cpl_malloc (3 * sizeof(cpl_vector*));
                    array[0] = gravi_vector_extract (opd_vector, 0, step);
                    cpl_vector_multiply_scalar (array[0], 1e6);
                    
                    int iA = gravi_get_region (detector_table, base, 'A', 0);
                    cpl_vector * tmp = gravi_table_get_vector (spectrum_table, wave, GRAVI_DATA[iA]);
                    cpl_vector_subtract_scalar (tmp, cpl_vector_get_mean (tmp));
                    array[1] = gravi_vector_extract (tmp, 0, step);
                    cpl_vector_divide_scalar (array[1], cpl_vector_get_max (array[1]));
                    cpl_vector_delete (tmp);

                    array[2] = gravi_vector_extract (envelope_vector, 0, step);
                    
                    cpl_plot_vectors (NULL, NULL, NULL, (const cpl_vector **)array, 3);
                    FREELOOP (cpl_vector_delete, array, 3);
                }


                
//                if (!(reg%4) && wave == nwave/2 && iA == 50) {
//                    int step = 100;
//                    cpl_vector ** array = cpl_calloc (4,sizeof(cpl_vector*));
//                    array[0] = gravi_vector_extract (opd_vector, 0, step);
//                    array[1] = gravi_vector_extract (X_vector, 0, step);
//                    array[2] = gravi_vector_extract (R_vector, 0, step);
//                    array[3] = gravi_vector_extract (I_vector, 0, step);
//                    cpl_plot_vectors (NULL, NULL, NULL, (const cpl_vector **)array, 4);
//                    cpl_free (array);
//                    CPLCHECK_NUL ("Cannot plot");
//                }


    // FIXME:
    cpl_msg_info (cpl_func, "Compute WAVE_SCAN for %s", GRAVI_TYPE(type_data));
    cpl_table * wavescan_table;
    wavescan_table = gravi_wave_scan (spectrum_table, detector_table, opd_table);
    
    gravi_data_add_table (wave_map, NULL, type_data == GRAVI_SC ? "WAVE_SCAN_SC" : "WAVE_SCAN_FT",
                          wavescan_table);
    


cpl_table * gravi_wave_scan (cpl_table * spectrum_table,
                             cpl_table * detector_table,
                             cpl_table * opd_table)
{
	gravi_msg_function_start(1);
    cpl_ensure (spectrum_table, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (detector_table, CPL_ERROR_NULL_INPUT, NULL);
    cpl_ensure (opd_table,      CPL_ERROR_NULL_INPUT, NULL);

    int nbase = 6;
    // char name[100];
    
    /* Create the output table */
    cpl_table * wave_table = cpl_table_new (1);

    /* Get the number of wavelength, region, polarisation... */
    cpl_size nwave = gravi_spectrum_get_nwave (spectrum_table);
    cpl_size nregion = cpl_table_get_nrow (detector_table);
    cpl_size nrow = cpl_table_get_nrow (spectrum_table);

    /* To save results */
    cpl_array * wave_array = cpl_array_new (nwave, CPL_TYPE_DOUBLE);

    /*
     * Calibration of each polarization and base
     */

    for (int reg = 0; reg < nregion; reg++) {
        cpl_msg_info (cpl_func, "Compute wave region %i over %lli", reg+1, nregion);

        /* Get the base of this region */
        int base = gravi_region_get_base (detector_table, reg);

        /* Sign of this baseline */
        double phi_sign = gravi_region_get_base_sign (detector_table, base);

        /* Get OPD of this region */
        cpl_vector * opd_vector = cpl_vector_new (nrow);
        for (cpl_size row = 0; row < nrow; row ++ ) {
            double value = cpl_table_get (opd_table, "OPD", row*nbase+base, NULL);
            cpl_vector_set (opd_vector, row, value);
        }
        CPLCHECK_NUL ("Cannot get opd");

        cpl_msg_info (cpl_func,"opd min = %g max = %g [um]",
                      cpl_vector_get_min (opd_vector)*1e6,
                      cpl_vector_get_max (opd_vector)*1e6);

        /* Loop on wave */
        for (cpl_size wave = 0; wave < nwave; wave++) {

            cpl_vector * X_vector;
            X_vector = gravi_table_get_vector (spectrum_table, wave, GRAVI_DATA[reg]);
            cpl_vector_subtract_scalar (X_vector, cpl_vector_get_mean (X_vector));
            cpl_vector_divide_scalar (X_vector, 1.4 * cpl_vector_get_stdev (X_vector));
            CPLCHECK_NUL ("Cannot get data");

            /* Guess wave */
            double lbd0 = 1.97e-6 + (2.48e-6 - 1.97e-6) / nwave * wave;

            /* Compute envelope from OPD for this channel */
            cpl_vector * env_vector = gravi_compute_envelope (opd_vector, wave, nwave);

            /* Search for a better match */
            cpl_size nA = (nwave == 5 ? 20 : 100);
            double searchA = 0.05;
            cpl_vector * V_vector = cpl_vector_new (nA);
            cpl_vector * R_vector = cpl_vector_new (nrow);
            cpl_vector * I_vector = cpl_vector_new (nrow);

            for (cpl_size iA = 0; iA < nA; iA++) {
                double A = phi_sign * (1 - searchA + (2*searchA * iA) / nA ) * CPL_MATH_2PI / lbd0;

                for (cpl_size row = 0; row < nrow; row++) {
                    cpl_vector_set (R_vector, row, cos (A * cpl_vector_get (opd_vector, row)) *
                                    cpl_vector_get (env_vector, row) );
                    cpl_vector_set (I_vector, row, sin (A * cpl_vector_get (opd_vector, row)) *
                                    cpl_vector_get (env_vector, row));
                }

                cpl_vector_multiply (R_vector, X_vector);
                cpl_vector_multiply (I_vector, X_vector);
                double Rvalue = cpl_vector_get_mean (R_vector);
                double Ivalue = cpl_vector_get_mean (I_vector);

                /* Fill chi2 */
                cpl_vector_set (V_vector, iA, Rvalue*Rvalue + Ivalue*Ivalue);
                CPLCHECK_NUL ("Cannot fill V_vector");
            }

            /* Plot chi2 versus iA */
            if (!(reg%8) && (wave == 0 || wave == nwave-1)) {
                cpl_plot_vector (NULL, NULL, NULL, V_vector);
                CPLCHECK_NUL ("Cannot plot");
            }
            
            /* Compute wavelength */
            double lbd = ( (1-searchA) + 2*searchA * gravi_vector_get_maxpos (V_vector) / nA ) * lbd0;
            cpl_array_set (wave_array, wave, lbd);
            
            FREE (cpl_vector_delete, env_vector);
            FREE (cpl_vector_delete, V_vector);
            FREE (cpl_vector_delete, R_vector);
            FREE (cpl_vector_delete, I_vector);
            FREE (cpl_vector_delete, X_vector);
            CPLCHECK_NUL ("Cannot delete");
        }

		/* Add column */
		gravi_table_new_column_array (wave_table, GRAVI_DATA[reg],
                                      "m", CPL_TYPE_DOUBLE, nwave);
		cpl_table_set_array (wave_table, GRAVI_DATA[reg], 0, wave_array);
        
        cpl_vector_delete (opd_vector);
    }

    cpl_array_delete (wave_array);

	gravi_msg_function_exit(1);
	return wave_table;
}

