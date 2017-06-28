/* $Id: gravi_dfs.c,v 1.6 2011/04/31 06:10:40 llundin Exp $
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

/**
 * @defgroup gravi_image  Image reconstruction (MIRA)
 */
/**@{*/

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#define _POSIX_C_SOURCE 200112L

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <math.h>
#include <cpl.h>
#include <time.h>
#include <sys/wait.h>
#include <signal.h>

#include "gravi_image.h"
#include "gravi_utils.h"

const int GRAVI_ERROR_TIMEOUT  = 0 + CPL_ERROR_EOL;

/*-----------------------------------------------------------------------------
                             Functions code
 -----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/**
 * @brief    Compute the image reconstruct from the OIFITS data
 * @param    input_data     the input OIFITS frame
 * @param    input_param    the input list of parameters
 * @return   cpl_image if OK, NULL else
 */
/*----------------------------------------------------------------------------*/

#ifdef YORICK_BIN
cpl_image * gravi_image(const cpl_frame * input_frame,
						const cpl_parameterlist * input_param)
{
    //cpl_errorstate prestate = cpl_errorstate_get();
    cpl_image * reconstruct_image;
    const char *input_filename;
	int status;
	char path[1035];
	int dim=0;
	double * image_buffer=NULL;
	int pid;								// No process
	int cf_pipe[2];							// Child to father pipe
	char* pixelsize_option, *dim_option, *regul_option, *regul_mu_option,
		*maxeval_option;
	const cpl_parameter * param;
	double val;
	int i;
	char *yorick_argv[20];
	int argv_i;
	double timeout;
	char * s;

	/* get the raw file name */
	input_filename=cpl_frame_get_filename(input_frame);
	cpl_msg_info(cpl_func, "Input filename ok %s", input_filename);
	if (input_filename==NULL) {
	        cpl_error_set_message(cpl_func, cpl_error_get_code(),
	                                          "Could not retrieve the input "
	                                          "filename");
	        return(NULL);
	    }

	/* generate the yorick argument table  */
	argv_i=0;
	yorick_argv[argv_i++] = cpl_sprintf("%s",YORICK_BIN);
	yorick_argv[argv_i++] = cpl_sprintf("-batch");
	FILE *file;
	file=fopen(YORICK_MIRA,"r");
	if(file!=NULL) {
		// file exists
		yorick_argv[argv_i++] = cpl_sprintf("%s",YORICK_MIRA);
		cpl_msg_info(cpl_func, "Mira is there: %s", yorick_argv[argv_i-1]);
		fclose(file);
	}
	else {
		yorick_argv[argv_i++] = cpl_sprintf("%s../../mira/mira-script.i",MIRADIR);
		cpl_msg_info(cpl_func, "Mira is there: %s", yorick_argv[argv_i-1]);
	//File not found
	}
	//yorick_argv[argv_i++] = cpl_sprintf("%s",YORICK_MIRA);

	/* RETRIEVE INPUT PARAMETERS */
	cpl_errorstate          prestate = cpl_errorstate_get();

    /* --pixelsize */
    param = cpl_parameterlist_find_const(input_param,
                          "gravi.gravity_image.pixelsize");
    pixelsize_option=cpl_sprintf("-pixelsize=%g",
									cpl_parameter_get_double(param));
    yorick_argv[argv_i++]=pixelsize_option;

    /* --dim */
	param = cpl_parameterlist_find_const(input_param,
						  "gravi.gravity_image.dim");
	dim=cpl_parameter_get_int(param);
	dim_option=cpl_sprintf("-dim=%d", dim);
	yorick_argv[argv_i++]=dim_option;

	/* --regul */
	param = cpl_parameterlist_find_const(input_param,
						  "gravi.gravity_image.regul");
	regul_option=cpl_sprintf("--regul=%s", cpl_parameter_get_string(param));
	yorick_argv[argv_i++]=regul_option;

	/* --regul_mu */
    param = cpl_parameterlist_find_const(input_param,
                          "gravi.gravity_image.regul_mu");
    regul_mu_option=cpl_sprintf("--regul_mu=%g",
									cpl_parameter_get_double(param));
    yorick_argv[argv_i++]=regul_mu_option;

    /* --maxeval */
	param = cpl_parameterlist_find_const(input_param,
						  "gravi.gravity_image.maxeval");
	maxeval_option=cpl_sprintf("-maxeval=%d", cpl_parameter_get_int(param));
	yorick_argv[argv_i++]=maxeval_option;

	if (!cpl_errorstate_is_equal(prestate)) {
	        cpl_error_set_message(cpl_func, cpl_error_get_code(),
	                                          "Could not retrieve the input "
	                                          "parameters");
	        for (i=0; i<argv_i-1; i++) cpl_free(yorick_argv[i]);
	        return(NULL);
	    }

	/* Other arguments */
	yorick_argv[argv_i++]  = cpl_sprintf("-regul_isotropic");
	yorick_argv[argv_i++]  = cpl_sprintf("-ftol=0");
	yorick_argv[argv_i++]  = cpl_sprintf("-gtol=0");
	yorick_argv[argv_i++] = cpl_sprintf("-overwrite");
	yorick_argv[argv_i++] = cpl_sprintf("-normalization=1.0");
	yorick_argv[argv_i++] = cpl_sprintf("-xmin=0.0");
	yorick_argv[argv_i++] = cpl_sprintf("%s",input_filename);
	yorick_argv[argv_i++] = cpl_sprintf("./output_temp.fits"); // not used by mira-script.i
	yorick_argv[argv_i++] = NULL;

	/* Get the timeout parameter */
	param = cpl_parameterlist_find_const(input_param,
                          "gravi.gravity_image.timeout");
    timeout=cpl_parameter_get_double(param);

	/*  Create the pipe */
	if (pipe(cf_pipe) != 0)
	{
		cpl_error_set_message(cpl_func, CPL_ERROR_ASSIGNING_STREAM,
								"Could not create pipe for Yorick process.");
		for (i=0; i<argv_i-1; i++) cpl_free(yorick_argv[i]);
		return NULL;
	}

	/* Create the child process to execute Yorick/Mira */
	switch(pid=fork())
	{
		/* Fork error */
		case (-1): // fork error
			cpl_error_set_message(cpl_func, CPL_ERROR_UNSPECIFIED,
									"Could not create the fork Yorick process.");
			for (i=0; i<argv_i-1; i++) cpl_free(yorick_argv[i]);
			return NULL;

		/* Child process (Yorick mira) : stdout to pipe*/
		case 0:
			cpl_msg_info(cpl_func, "Start Yorick process (pid: %u)",  getpid());

			/* Close the unused side */
			close(cf_pipe[0]);

			/* redirect stdout to the pipe */
			dup2(cf_pipe[1], STDOUT_FILENO);

			/* Execute the mira script */
			if (execv(yorick_argv[0], yorick_argv) == -1 ){
				cpl_msg_error("MIRA/Yorick", "Error in Yorick call (%s)",
						YORICK_BIN);
			}

			// Child process stop
			printf("Child process (%d) - Stop\n", getpid());
			fflush(NULL);
			exit(0);

		/* Parent process : stdin from pipe */
		default:

			/* Close the unused side of the pipe */
			close(cf_pipe[1]);

			/* Redirect the pipe output to STDIN */
			dup2(cf_pipe[0], STDIN_FILENO);

			/* Read the output a line at a time - output it. */
			time_t ticks1, ticks2;


			ticks1=time(NULL);
			while (fgets(path, sizeof(path)-1, stdin) != NULL) {
				/* IMAGE				 */

				/* delete the trailing \n 				 */
			    if ( (s=strrchr(path, '\n') )) s[0]='\0';

				if (strncmp(path, "START_IMAGE", strlen(path)-1) == 0){
					cpl_msg_info(cpl_func, "Receive image");
					image_buffer=cpl_malloc(sizeof(double)*dim*dim);
					/* get the image pixels */
					i=0;
					while (fgets(path, sizeof(path)-1, stdin) != NULL &&
							strncmp(path, "END_IMAGE", strlen(path)-1) != 0){
						sscanf(path, "%lg", &val);
						image_buffer[i++]=val;
					}
					//printf("Image end \n");
				}
				/* MIRA WARNING				 */
				else if (strncmp(path, "WARNING", 7) == 0){
					cpl_msg_warning("MIRA/Yorick", "%s", path);
				}
				/* MIRA WARNING				 */
				else if (strncmp(path, "# warning", 9) == 0){
					cpl_msg_warning("MIRA/Yorick", "%s", path);
				}
				/* MIRA ERROR				 */
				else if (strncmp(path, "ERROR", 5) == 0){
					cpl_msg_error("MIRA/Yorick", "%s", path);
				}
				/* MIRA INFO				 */
				else {
					cpl_msg_info("MIRA/Yorick", "%s", path);
				}

				ticks2=time(NULL);
				if (difftime(ticks2,ticks1) > timeout)
				{
					cpl_error_set_message(cpl_func, GRAVI_ERROR_TIMEOUT,
						"Timeout (%g s) MIRA/Yorick process "
						"killed. Use --timeout option to change it.", timeout);
					cpl_msg_error(cpl_func, "Timeout MIRA/Yorick process Killed");
					kill((pid_t) pid, SIGKILL);
				}
			}
			/* wait end of child process */
			wait(&status);
			cpl_msg_info(cpl_func, "End of Yorick process (pid: %u with satus %d)",
					pid, status);
			/* Create reconstructed image */
			if (image_buffer==NULL){
				cpl_error_set_message(cpl_func, CPL_ERROR_NULL_INPUT,
						"Could not get image from Yorick process.");
				for (i=0; i<argv_i-1; i++) cpl_free(yorick_argv[i]);
				return(NULL);
			}
			reconstruct_image = cpl_image_wrap_double(dim, dim, image_buffer);
	}

	for (i=0; i<argv_i-1; i++) cpl_free(yorick_argv[i]);
    return reconstruct_image;
}
#endif

/**@}*/
