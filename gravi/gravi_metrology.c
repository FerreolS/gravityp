/* $Id: gravi_vis.c,v 1.10 2014/11/12 15:10:40 nazouaoui Exp $
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
 * @defgroup gravi_metrology  Metrology reduction
 */
/**@{*/

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define _XOPEN_SOURCE
#include <cpl.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "gravi_data.h"
#include "gravi_pfits.h"
#include "gravi_dfs.h"

#include "gravi_cpl.h"
#include "gravi_utils.h"
#include "gravi_ellipse.h"
#include "gravi_signal.h"
#include "gravi_eop.h"

#include "gravi_metrology.h"

/*-----------------------------------------------------------------------------
                               Private prototypes
 -----------------------------------------------------------------------------*/


cpl_error_code gravi_metrology_tac (cpl_table * metrology_table,
                                    cpl_table * vismet_table,
                                    cpl_propertylist * header);

cpl_error_code gravi_metrology_acq (cpl_table * visacq_table,
                                    cpl_table * vismet_table,
                                    double delay,
                                    cpl_propertylist * header);

long gravi_round (double number);
long gravi_round (double number)
{
    return (number >= 0) ? (long)(number + 0.5) : (long)(number - 0.5);
}

/*-----------------------------------------------------------------------------
                               TAC definitions
 -----------------------------------------------------------------------------*/

#define DEBUG_LEVEL 2

#define PI 3.14159265359
#define TWOPI 6.28318530718

#define MAX_DATA 1000000
#define MET_MAX_HISTORY 100

#define FLAG_FLUX 1         /* error L on gui */
#define FLAG_SNR 2          /* error K on gui */
#define FLAG_UNWRAP 4       /* error J on gui */
#define FLAG_SPEED 8        /* error I on gui */
#define FLAG_JUMP 16        /* error H on gui */
#define FLAG_OVERLOAD 32    /* error G on gui */
#define FLAG_UNLOCKED 64    /* error F on gui */
#define FLAG_LASER_MV 128   /* error E on gui */
#define FLAG_LASER_MW 256   /* error D on gui */
#define FLAG_LASER_WAVE 512 /* error C on gui */
#define FLAG_VOLT 1024      /* error B on gui */
#define FLAG_PHASOR0 2048   /* error A on gui */
#define FLAG_LOST_PHASE 4096   /* error T on gui */

#define FT 0
#define SC 1
#define BOTH 2

#define SIN 0
#define COS 1

#define NCOL 80

/* Structure definitions */

typedef struct {
    
    double min_allowed_flux_telescope[4][4][2]; /* minimum allowed flux per telescope and diode and side */
    double min_allowed_flux_fiber_coupler[4][2]; /* minimum allowed flux per fiber coupler and side  */
    double min_allowed_snr_telescope[4][4][2]; /* minimum allowed SNR per telescope and diode and side  */
    double min_allowed_snr_fiber_coupler[4][2]; /* minimum allowed SNR per fiber coupler and side  */
    double min_allowed_snr_diff[4][4][2]; /* minimum allowed SNR per telescope and diode and side  */
    long number_to_average; /* how many cycles used for smoothing signal */
    long number_for_rms; /* how many cycles used for calculation of rms */
    long number_for_speed; /* how many cycles used for calculation of speed */
    double nominal_wavelength; /* nominal wavelength of laser in nm */
    double max_allowed_phase_difference; /* in rad */
    double max_allowed_phase_speed; /* in rad per phase step */
    double opl_base; /* proportionality constant between radian and meters */
    double offset_volt_telescope[4][4][2][2]; /* offset voltage per telescope and diode and side and sin/cos */
    double offset_volt_fiber_coupler[4][2][2]; /* offset voltage per fiber coupler and side and sin/cos */
    double rms_flux_telescope[4][4][2]; /* rms of voltage per telescope and diode and side */
    double rms_flux_fiber_coupler[4][2]; /* rms of voltage per fiber coupler and side */
    double lockin_constant_telescope[4][4][2][2]; /* multiplicative constant for lockin voltage per telescope and diode and side and sin/cos */
    double lockin_constant_fiber_coupler[4][2][2]; /* multiplicative constant for lockin voltage per fiber coupler and side and sin/cos */
    double max_allowed_voltage;
    double min_allowed_voltage;
    long number_to_smooth_for_telescope; /* how many cycles are smoothed for difference Tel - FC phasor */
    long sign_of_phase;
    long check_for_jumps;
    long calc_phase_speed;
    double decrement_factor_speed; /* how much the weight of each speed value in buffer is decreased after a cycle */
    double norm_speed; /* normalization factor for weighted speed adding */
    double sigma_clip_speed; /* at what sigma level the speed predictor is clipped */
    
} structTacConfiguration;


typedef struct {
    long sample_number; /* Just a TAC loop counter : 0, 1, 2, ... */
    long buffer_idx_avg;
    
    double volt_lockin_telescope[4][4][2][2]; /* per telescope and diode and side and sin/cos */
    double sum_volt_lockin_telescope[4][4][2][2]; /* per telescope and diode and side and sin/cos */
    double buffer_volt_lockin_telescope[4][4][2][2][MET_MAX_HISTORY]; /* per telescope and diode and side and sin/cos */
    double flux_telescope[4][4][2];
    double sum_flux_telescope[4][4][2]; /* 4 x 4 x FT/SC */
    double sum_sq_flux_telescope[4][4][2]; /* 4 x 4 x FT/SC */
    double buffer_flux_telescope[4][4][2][MET_MAX_HISTORY]; /* 4 x 4 x FT/SC */
    double buffer_sq_flux_telescope[4][4][2][MET_MAX_HISTORY]; /* 4 x 4 x FT/SC */
    double rms_flux_telescope[4][4][2]; /* 4 x 4 x FT/SC */
    double snr_flux_telescope[4][4][2]; /* 4 x 4 x FT/SC */
    double buffer_delta_phasor_telescope[4][4][2][2][MET_MAX_HISTORY];
    double delta_phasor_telescope[4][4][2][2];
    double best_estimate_phasor_telescope[4][4][2][2];
    double best_estimate_flux_telescope[4][4][2];
    
    double buffer_delta_phase_telescope[4][4][2][MET_MAX_HISTORY];
    double sum_delta_phase_telescope[4][4][2];
    double buffer_sq_delta_phase_telescope[4][4][2][MET_MAX_HISTORY];
    double sum_sq_delta_phase_telescope[4][4][2];
    
    double sum_speed_telescope[4][4][2];
    double buffer_speed_telescope[4][4][2][MET_MAX_HISTORY];
    
    int flag_volt_telescope[4][4][2];
    int flag_flux_telescope[4][4][2]; /* 4 x 4 x FT/SC */
    int flag_snr_telescope[4][4][2]; /* 4 x 4 x FT/SC */
    
    double volt_lockin_fiber_coupler[4][2][2]; /* per fiber coupler and side and sin/cos */
    double sum_volt_lockin_fiber_coupler[4][2][2]; /* per fiber coupler and side and sin/cos */
    double buffer_volt_lockin_fiber_coupler[4][2][2][MET_MAX_HISTORY]; /* per fiber coupler and side and sin/cos */
    double flux_fiber_coupler[4][2];
    double sum_flux_fiber_coupler[4][2];
    double sum_sq_flux_fiber_coupler[4][2];
    double buffer_flux_fiber_coupler[4][2][MET_MAX_HISTORY];
    double buffer_sq_flux_fiber_coupler[4][2][MET_MAX_HISTORY];
    double rms_flux_fiber_coupler[4][2];
    double snr_flux_fiber_coupler[4][2];
    double buffer_delta_phasor_fiber_coupler[4][2][2][MET_MAX_HISTORY];
    double delta_phasor_fiber_coupler[4][2][2];
    double best_estimate_phasor_fiber_coupler[4][2][2];
    double best_estimate_flux_fiber_coupler[4][2];
    double buffer_best_estimate_phasor_fiber_coupler[4][2][2][MET_MAX_HISTORY];
    
    double buffer_delta_phase_fiber_coupler[4][2][MET_MAX_HISTORY];
    double sum_delta_phase_fiber_coupler[4][2];
    double buffer_sq_delta_phase_fiber_coupler[4][2][MET_MAX_HISTORY];
    double sum_sq_delta_phase_fiber_coupler[4][2];
    
    double sum_speed_fiber_coupler[4][2];
    double buffer_speed_fiber_coupler[4][2][MET_MAX_HISTORY];
    double buffer_used_speed_fiber_coupler[4][2][MET_MAX_HISTORY];
    
    int flag_volt_fiber_coupler[4][2];
    int flag_flux_fiber_coupler[4][2]; /* 4 fiber couplers, sin/cos */
    int flag_snr_fiber_coupler[4][2]; /* 4 fiber couplers, sin/cos */
    
    double raw_phase_telescope[4][4][2]; /* 4 x 4 (tel, diode), FT/SC */
    double prev_phase_telescope[4][4][2];
    double unwrapped_phase_telescope[4][4][2];
    int flag_unwrap_telescope[4][4][2];
    int flag_jump_telescope[4][4][2];
    int flag_speed_telescope[4][4][2];
    int flag_phasor0_telescope[4][4][2];
    double delta_phase_telescope[4][4]; /* FT - SC */
    double mean_phase_telescope[4]; /* mean over tge 4 diodes */
    
    double raw_phase_fiber_coupler[4][2]; /* 4 fiber couplers , FT/SC */
    double prev_phase_fiber_coupler[4][2];
    double unwrapped_phase_fiber_coupler[4][2];
    double buffer_unwrapped_phase_fiber_coupler[4][2][MET_MAX_HISTORY];
    int flag_unwrap_fiber_coupler[4][2];
    int flag_jump_fiber_coupler[4][2];
    int flag_speed_fiber_coupler[4][2];
    int flag_phasor0_fiber_coupler[4][2];
    double delta_phase_fiber_coupler[4];  /* FT - SC */
    
    double raw_phase_diff[4][4][2]; /* difference: tel 1/2/3/4 minus FC */
    double prev_phase_diff[4][4][2];
    double unwrapped_phase_diff[4][4][2];
    int flag_unwrap_diff[4][4][2];
    int flag_phasor0_diff[4][4][2];
    int flag_snr_diff[4][4][2];
    int flag_flux_diff[4][4][2];
    
    double buffer_phasor_diff[4][4][2][2][MET_MAX_HISTORY]; /* in volt, tel diode, side, sin/cos, index */
    double sum_phasor_diff[4][4][2][2]; /* in volt, tel diode, side, sin/cos */
    double flux_diff[4][4][2]; /* phasor length in volt, tel, diode, side */
    double sum_flux_diff[4][4][2]; /* 4 x 4 x FT/SC */
    double sum_sq_flux_diff[4][4][2]; /* 4 x 4 x FT/SC */
    double buffer_flux_diff[4][4][2][MET_MAX_HISTORY]; /* 4 x 4 x FT/SC */
    double buffer_sq_flux_diff[4][4][2][MET_MAX_HISTORY]; /* 4 x 4 x FT/SC */
    double rms_flux_diff[4][4][2]; /* 4 x 4 x FT/SC */
    double snr_flux_diff[4][4][2]; /* 4 x 4 x FT/SC */
    
    int total_flag_telescope[4][4][2]; /* 4 x 4 (tel, diode), FT/SC */
    int total_flag_fiber_coupler[4][2]; /* 4 fiber couplers , FT/SC */
    
    double opl_telescope[4]; /* diff FT-SC in meter */
    double opl_telescope_diode[4][4]; /* diff FT-SC in meter */
    double opl_fiber_coupler[4]; /* diff FT-SC in meter */
    double opl_zero_telescope[4]; /* zero point in meter */
    double opl_zero_telescope_diode[4][4]; /* zero point in meter */
    double opl_zero_fiber_coupler[4]; /* zero point  meter */
    
    double start_phase_fiber_coupler[4][2]; /* in rad */
    double start_phase_telescope[4][4][2]; /* in rad */
    
    long freeze_count_fiber_coupler[4][2];
    double freeze_phasor_fiber_coupler[4][2][2];
    double freeze_speed_fiber_coupler[4][2];
    double freeze_unwrapped_phase_fiber_coupler[4][2];
    double used_unwrapped_phase_fiber_coupler[4][2];
    
    structTacConfiguration * tacConfiguration;
    
} structTacData;

/*-----------------------------------------------------------------------------
                                 TAC prototypes
 -----------------------------------------------------------------------------*/

double metrology_sq(double x);
double myAtan(double x, double y, int* flag);

structTacConfiguration * metrology_makeDefaultTacConfiguration (double lambda_met);

structTacData * metrology_makeDefaultTacData (double lambda_met);

int metrology_algorithm(structTacData * tacData);

int metrology_unwrap(double raw_phase, double previous_phase, double max_allowed_phase_diff,
                     double *unwrapped_phase);

int metrology_read_voltages(structTacData * tacData, double * volt);


/*-----------------------------------------------------------------------------
                               TAC function code
 -----------------------------------------------------------------------------*/

structTacConfiguration * metrology_makeDefaultTacConfiguration(double lambda_met)
{
    structTacConfiguration *defaultTacConfiguration;
    
    defaultTacConfiguration = cpl_malloc((size_t) sizeof(structTacConfiguration));
    
    long tel, diode, side, comp, i;
    
    for (tel = 0; tel < 4; tel++) {
        for(side = FT; side <= SC; side++) {
            for (diode = 0; diode < 4; diode++) {
                for (comp = SIN; comp <= COS; comp++) {
                    defaultTacConfiguration->offset_volt_telescope[tel][diode][side][comp] = 0.0;
                    defaultTacConfiguration->lockin_constant_telescope[tel][diode][side][comp] = 1.0;
                }
                defaultTacConfiguration->min_allowed_flux_telescope[tel][diode][side] = 0.01;
                defaultTacConfiguration->min_allowed_snr_telescope[tel][diode][side] = 2.0;
                defaultTacConfiguration->min_allowed_snr_diff[tel][diode][side] = 2.0;
                defaultTacConfiguration->rms_flux_telescope[tel][diode][side] = 0.005;
            }
            for (comp = SIN; comp <= COS; comp++) {
                defaultTacConfiguration->offset_volt_fiber_coupler[tel][side][comp] = 0.01;
                defaultTacConfiguration->lockin_constant_fiber_coupler[tel][side][comp] = 1.0;
            }
            defaultTacConfiguration->min_allowed_flux_fiber_coupler[tel][side] = 0.01;
            defaultTacConfiguration->min_allowed_snr_fiber_coupler[tel][side] = 2.0;
            defaultTacConfiguration->rms_flux_fiber_coupler[tel][side] = 0.004;
        }
    }
    
    defaultTacConfiguration->min_allowed_voltage = -11.0;
    defaultTacConfiguration->max_allowed_voltage = 11.0;
    
    defaultTacConfiguration->number_to_average = 20;
    defaultTacConfiguration->number_to_average = 20;
    defaultTacConfiguration->number_for_rms = 50;
    
    defaultTacConfiguration->max_allowed_phase_difference = TWOPI/3;
    defaultTacConfiguration->max_allowed_phase_speed = TWOPI/4;

    defaultTacConfiguration->nominal_wavelength = lambda_met * 1e9;
    defaultTacConfiguration->opl_base = defaultTacConfiguration->nominal_wavelength / 1e9 / TWOPI;
    
    defaultTacConfiguration->number_to_smooth_for_telescope = 50;
    defaultTacConfiguration->sign_of_phase = 1;
    
    defaultTacConfiguration->check_for_jumps = 1;
    defaultTacConfiguration->calc_phase_speed = 1;
    
    defaultTacConfiguration->number_for_speed = 3;
    defaultTacConfiguration->decrement_factor_speed = 0.8;
    defaultTacConfiguration->sigma_clip_speed = 0.0;
    
    defaultTacConfiguration->norm_speed = 0;
    for(i=0; i<defaultTacConfiguration->number_for_speed; i++) {
        defaultTacConfiguration->norm_speed += pow(defaultTacConfiguration->decrement_factor_speed,i);
    }
    
    return defaultTacConfiguration;
}

structTacData * metrology_makeDefaultTacData(double lambda_met)
{
    long tel, diode, idx, side, comp;
    structTacData *defaultTacData;
    
    defaultTacData = cpl_malloc((size_t) sizeof(structTacData));
    
    defaultTacData->tacConfiguration = metrology_makeDefaultTacConfiguration(lambda_met);
    
    defaultTacData->sample_number = 0;
    
    for (tel = 0; tel < 4; tel++) {
        for(side = FT; side <= SC; side++) {
            for (diode = 0; diode < 4; diode++) {
                for (comp = SIN; comp <= COS; comp++) {
                    defaultTacData->volt_lockin_telescope[tel][diode][side][comp] = 0.0;
                    defaultTacData->sum_volt_lockin_telescope[tel][diode][side][comp] = 0.0;
                    defaultTacData->delta_phasor_telescope[tel][diode][side][comp] = 0.0;
                    defaultTacData->best_estimate_phasor_telescope[tel][diode][side][comp] = 0.0;
                    defaultTacData->sum_phasor_diff[tel][diode][side][comp] = 0.0;
                }
                defaultTacData->flux_telescope[tel][diode][side] = 0.0;
                defaultTacData->rms_flux_telescope[tel][diode][side] = 0.0;
                defaultTacData->snr_flux_telescope[tel][diode][side] = 0.0;
                defaultTacData->sum_flux_telescope[tel][diode][side] = 0.0;
                defaultTacData->sum_sq_flux_telescope[tel][diode][side] = 0.0;
                defaultTacData->flag_flux_telescope[tel][diode][side] = 0;
                defaultTacData->flag_snr_telescope[tel][diode][side] = 0;
                defaultTacData->raw_phase_telescope[tel][diode][side] = 0.0;
                defaultTacData->prev_phase_telescope[tel][diode][side] = 0.0;
                defaultTacData->unwrapped_phase_telescope[tel][diode][side] = 0.0;
                defaultTacData->flag_unwrap_telescope[tel][diode][side] = 0;
                defaultTacData->flag_speed_telescope[tel][diode][side] = 0;
                defaultTacData->flag_jump_telescope[tel][diode][side] = 0;
                defaultTacData->flag_phasor0_telescope[tel][diode][side] = 0;
                defaultTacData->total_flag_telescope[tel][diode][side] = 0;
                defaultTacData->flag_volt_telescope[tel][diode][side] = 0;
                defaultTacData->start_phase_telescope[tel][diode][side] = 0.0;
                defaultTacData->best_estimate_flux_telescope[tel][diode][side] = 0.0;
                defaultTacData->sum_speed_telescope[tel][diode][side] = 0.0;
                
                for (idx = 0; idx < MET_MAX_HISTORY; idx++) {
                    for (comp = SIN; comp <= COS; comp++) {
                        defaultTacData->buffer_volt_lockin_telescope[tel][diode][side][comp][idx] = 0.0;
                        defaultTacData->buffer_delta_phasor_telescope[tel][diode][side][comp][idx] = 0.0;
                        defaultTacData->buffer_phasor_diff[tel][diode][side][comp][idx] = 0.0;
                    }
                    defaultTacData->buffer_delta_phase_telescope[tel][diode][side][idx] = 0.0;
                    defaultTacData->buffer_sq_delta_phase_telescope[tel][diode][side][idx] = 0.0;
                    defaultTacData->buffer_flux_telescope[tel][diode][side][idx] = 0.0;
                    defaultTacData->buffer_sq_flux_telescope[tel][diode][side][idx] = 0.0;
                    defaultTacData->buffer_speed_telescope[tel][diode][side][idx] = 0.0;
                    defaultTacData->buffer_flux_diff[tel][diode][side][idx] = 0.0;
                    defaultTacData->buffer_sq_flux_diff[tel][diode][side][idx] = 0.0;
                }
                defaultTacData->flux_diff[tel][diode][side] = 0.0;
                defaultTacData->rms_flux_diff[tel][diode][side] = 0.0;
                defaultTacData->snr_flux_diff[tel][diode][side] = 0.0;
                defaultTacData->sum_flux_diff[tel][diode][side] = 0.0;
                defaultTacData->sum_sq_flux_diff[tel][diode][side] = 0.0;
                
                defaultTacData->sum_delta_phase_telescope[tel][diode][side] = 0.0;
                defaultTacData->sum_sq_delta_phase_telescope[tel][diode][side] = 0.0;
                
                defaultTacData->raw_phase_diff[tel][diode][side] = 0.0;
                defaultTacData->prev_phase_diff[tel][diode][side] = 0.0;
                defaultTacData->unwrapped_phase_diff[tel][diode][side] = 0.0;
                defaultTacData->flag_unwrap_diff[tel][diode][side] = 0;
                defaultTacData->flag_phasor0_diff[tel][diode][side] = 0;
                defaultTacData->flag_snr_diff[tel][diode][side] = 0;
            } /* end loop diode */
            for (comp = SIN; comp <= COS; comp++) {
                defaultTacData->volt_lockin_fiber_coupler[tel][side][comp] = 0.0;
                defaultTacData->sum_volt_lockin_fiber_coupler[tel][side][comp] = 0.0;
                defaultTacData->delta_phasor_fiber_coupler[tel][side][comp] = 0.0;
                defaultTacData->best_estimate_phasor_fiber_coupler[tel][side][comp] = 0.0;
                defaultTacData->freeze_phasor_fiber_coupler[tel][side][comp] = 0.0;
            }
            defaultTacData->flux_fiber_coupler[tel][side] = 0.0;
            defaultTacData->rms_flux_fiber_coupler[tel][side] = 0.0;
            defaultTacData->snr_flux_fiber_coupler[tel][side] = 0.0;
            defaultTacData->sum_flux_fiber_coupler[tel][side] = 0.0;
            defaultTacData->sum_sq_flux_fiber_coupler[tel][side] = 0.0;
            defaultTacData->flag_flux_fiber_coupler[tel][side] = 0;
            defaultTacData->flag_snr_fiber_coupler[tel][side] = 0;
            defaultTacData->raw_phase_fiber_coupler[tel][side] = 0.0;
            defaultTacData->prev_phase_fiber_coupler[tel][side] = 0.0;
            defaultTacData->unwrapped_phase_fiber_coupler[tel][side] = 0.0;
            defaultTacData->flag_unwrap_fiber_coupler[tel][side] = 0;
            defaultTacData->flag_jump_fiber_coupler[tel][side] = 0;
            defaultTacData->flag_speed_fiber_coupler[tel][side] = 0;
            defaultTacData->flag_phasor0_fiber_coupler[tel][side] = 0;
            defaultTacData->total_flag_fiber_coupler[tel][side] = 0;
            defaultTacData->flag_volt_fiber_coupler[tel][side] = 0;
            defaultTacData->start_phase_fiber_coupler[tel][side] = 0.0;
            defaultTacData->best_estimate_flux_fiber_coupler[tel][side] = 0.0;
            defaultTacData->sum_speed_fiber_coupler[tel][side] = 0.0;
            defaultTacData->freeze_count_fiber_coupler[tel][side] = 0;
            defaultTacData->freeze_speed_fiber_coupler[tel][side] = 0.0;
            defaultTacData->freeze_unwrapped_phase_fiber_coupler[tel][side] = 0.0;
            defaultTacData->used_unwrapped_phase_fiber_coupler[tel][side] = 0.0;
            
            for (idx = 0; idx < MET_MAX_HISTORY; idx++) {
                for (comp = SIN; comp <= COS; comp++) {
                    defaultTacData->buffer_volt_lockin_fiber_coupler[tel][side][comp][idx] = 0.0;
                    defaultTacData->buffer_delta_phasor_fiber_coupler[tel][side][comp][idx] = 0.0;
                    defaultTacData->buffer_best_estimate_phasor_fiber_coupler[tel][side][comp][idx] = 0.0;
                }
                defaultTacData->buffer_delta_phase_fiber_coupler[tel][side][idx] = 0.0;
                defaultTacData->buffer_sq_delta_phase_fiber_coupler[tel][side][idx] = 0.0;
                defaultTacData->buffer_flux_fiber_coupler[tel][side][idx] = 0.0;
                defaultTacData->buffer_sq_flux_fiber_coupler[tel][side][idx] = 0.0;
                defaultTacData->buffer_speed_fiber_coupler[tel][side][idx] = 0.0;
                defaultTacData->buffer_used_speed_fiber_coupler[tel][side][idx] = 0.0;
                defaultTacData->buffer_unwrapped_phase_fiber_coupler[tel][side][idx] = 0.0;
            }
            defaultTacData->sum_delta_phase_fiber_coupler[tel][side] = 0.0;
            defaultTacData->sum_sq_delta_phase_fiber_coupler[tel][side] = 0.0;
        }  /* end loop over side */
        defaultTacData->mean_phase_telescope[tel] = 0;
        
        defaultTacData->delta_phase_fiber_coupler[tel] = 0.0;
        for (diode = 0; diode < 4; diode++) {
            defaultTacData->delta_phase_telescope[tel][diode] = 0.0;
            defaultTacData->opl_telescope_diode[tel][diode] = 0.0;
            defaultTacData->opl_zero_telescope_diode[tel][diode] = 0.0;
        }
        
        defaultTacData->opl_fiber_coupler[tel] = 0.0;
        defaultTacData->opl_zero_fiber_coupler[tel] = 0.0;
        defaultTacData->opl_telescope[tel] = 0.0;
        defaultTacData->opl_zero_telescope[tel] = 0.0;
        
    } /* end loop over tel */
    
    return defaultTacData;
    
}

int metrology_read_voltages(structTacData * tacData, double * volt)
{
    int err = 0;
    int tel, diode, side, comp;
    int idx = 0;
    long buffer_idx_avg = tacData->buffer_idx_avg;
    double sqflux, flux;
    structTacConfiguration * tacConfiguration = tacData->tacConfiguration;
    long buffer_idx_rms = ((tacData->sample_number - 1) % tacConfiguration->number_for_rms);

    
    for (side = FT; side <= SC; side++) {
        for (tel = 0; tel < 4; tel++) {
            for (diode = 0; diode < 4; diode++) {
                for (comp = SIN; comp <= COS; comp++) {
                    tacData->volt_lockin_telescope[tel][diode][side][comp] = volt[idx];
                    idx++;
                }
            }
        }
    }
    for (side = FT; side <= SC; side++) {
        for (tel = 0; tel < 4; tel++) {
            for (comp = SIN; comp <= COS; comp++) {
                tacData->volt_lockin_fiber_coupler[tel][side][comp] = volt[idx];
                idx++;
            }
        }
    }
        
    for (tel = 0; tel < 4; tel++) {
        for (side = FT; side <= SC; side++) {
            for (diode = 0; diode < 4; diode++) {
                tacData->flag_volt_telescope[tel][diode][side] = 0;
                if( tacData->volt_lockin_telescope[tel][diode][side][SIN] > tacConfiguration->max_allowed_voltage ||
                   tacData->volt_lockin_telescope[tel][diode][side][SIN] < tacConfiguration->min_allowed_voltage ||
                   tacData->volt_lockin_telescope[tel][diode][side][COS] > tacConfiguration->max_allowed_voltage ||
                   tacData->volt_lockin_telescope[tel][diode][side][COS] < tacConfiguration->min_allowed_voltage) {
                    tacData->flag_volt_telescope[tel][diode][side] = FLAG_VOLT;
                }
                tacData->volt_lockin_telescope[tel][diode][side][SIN] -=
                tacConfiguration->offset_volt_telescope[tel][diode][side][SIN];
                tacData->volt_lockin_telescope[tel][diode][side][SIN] *=
                tacConfiguration->lockin_constant_telescope[tel][diode][side][SIN];
                
                tacData->volt_lockin_telescope[tel][diode][side][COS] -=
                tacConfiguration->offset_volt_telescope[tel][diode][side][COS];
                tacData->volt_lockin_telescope[tel][diode][side][COS] *=
                tacConfiguration->lockin_constant_telescope[tel][diode][side][COS];
                
            }
            tacData->flag_volt_fiber_coupler[tel][side] = 0;
            if( tacData->volt_lockin_fiber_coupler[tel][side][SIN] > tacConfiguration->max_allowed_voltage ||
               tacData->volt_lockin_fiber_coupler[tel][side][SIN] < tacConfiguration->min_allowed_voltage ||
               tacData->volt_lockin_fiber_coupler[tel][side][COS] > tacConfiguration->max_allowed_voltage ||
               tacData->volt_lockin_fiber_coupler[tel][side][COS] < tacConfiguration->min_allowed_voltage) {
                tacData->flag_volt_fiber_coupler[tel][side] = FLAG_VOLT;
            }
            tacData->volt_lockin_fiber_coupler[tel][side][SIN] -=
            tacConfiguration->offset_volt_fiber_coupler[tel][side][SIN];
            tacData->volt_lockin_fiber_coupler[tel][side][SIN] *=
            tacConfiguration->lockin_constant_fiber_coupler[tel][side][SIN];
            tacData->volt_lockin_fiber_coupler[tel][side][COS] -=
            tacConfiguration->offset_volt_fiber_coupler[tel][side][COS];
            tacData->volt_lockin_fiber_coupler[tel][side][COS] *=
            tacConfiguration->lockin_constant_fiber_coupler[tel][side][COS];
        }
    }
    for (tel = 0; tel < 4; tel++) {
        for (side = FT; side <= SC; side++) {
            for (diode = 0; diode < 4; diode++) {
                
                tacData->buffer_volt_lockin_telescope[tel][diode][side][SIN][buffer_idx_avg] = tacData->volt_lockin_telescope[tel][diode][side][SIN];
                tacData->buffer_volt_lockin_telescope[tel][diode][side][COS][buffer_idx_avg] = tacData->volt_lockin_telescope[tel][diode][side][COS];
                sqflux = metrology_sq(tacData->volt_lockin_telescope[tel][diode][side][SIN]) + metrology_sq(tacData->volt_lockin_telescope[tel][diode][side][COS]);
                flux=sqrt(sqflux);
                tacData->flux_telescope[tel][diode][side] = flux;
                tacData->sum_flux_telescope[tel][diode][side] -= tacData->buffer_flux_telescope[tel][diode][side][buffer_idx_rms];
                tacData->sum_sq_flux_telescope[tel][diode][side] -= tacData->buffer_sq_flux_telescope[tel][diode][side][buffer_idx_rms];
                tacData->buffer_flux_telescope[tel][diode][side][buffer_idx_rms] = flux;
                tacData->buffer_sq_flux_telescope[tel][diode][side][buffer_idx_rms] = sqflux;
                tacData->sum_flux_telescope[tel][diode][side] += flux;
                tacData->sum_sq_flux_telescope[tel][diode][side] += sqflux;
                
            } /* end of loop over diodes */
            
            tacData->buffer_volt_lockin_fiber_coupler[tel][side][SIN][buffer_idx_avg] = tacData->volt_lockin_fiber_coupler[tel][side][SIN];
            tacData->buffer_volt_lockin_fiber_coupler[tel][side][COS][buffer_idx_avg] = tacData->volt_lockin_fiber_coupler[tel][side][COS];
            sqflux = metrology_sq(tacData->volt_lockin_fiber_coupler[tel][side][SIN]) + metrology_sq(tacData->volt_lockin_fiber_coupler[tel][side][COS]);
            flux=sqrt(sqflux);
            tacData->flux_fiber_coupler[tel][side] = flux;
            tacData->sum_flux_fiber_coupler[tel][side] -= tacData->buffer_flux_fiber_coupler[tel][side][buffer_idx_rms];
            tacData->sum_sq_flux_fiber_coupler[tel][side] -= tacData->buffer_sq_flux_fiber_coupler[tel][side][buffer_idx_rms];
            tacData->buffer_flux_fiber_coupler[tel][side][buffer_idx_rms] = flux;
            tacData->buffer_sq_flux_fiber_coupler[tel][side][buffer_idx_rms] = sqflux;
            tacData->sum_flux_fiber_coupler[tel][side] += flux;
            tacData->sum_sq_flux_fiber_coupler[tel][side] += sqflux;
            
        }	 /* end side loop */
    } /* end telescope loop */
    
    return err;
}

int metrology_unwrap(double raw_phase, double previous_phase, double max_allowed_phase_diff,
           double *unwrapped_phase) {
    
    double previous_phase_mod, twopicount, diff;
    int flag = 0;
    
    twopicount = floor(previous_phase / TWOPI); /* the number of full circles */
    previous_phase_mod = previous_phase - TWOPI * twopicount; /* subtract the number of full circles to get the actual fraction of circle */
    
    if (previous_phase_mod > PI) { /* if the fraction is bigger the PI subtract 2PI to get the rest */
        previous_phase_mod -= TWOPI;
    }
    
    diff = raw_phase - previous_phase_mod;
    
    if (diff > PI) {
        diff -= TWOPI;
    } else if (diff < (-1.0 * PI)) {
        diff += TWOPI;
    }
    
    if (fabs(diff) > max_allowed_phase_diff) {
        flag = FLAG_UNWRAP;
    }
    
    *unwrapped_phase = (previous_phase + diff);
    
    return flag;
}

double metrology_sq(double x) {
    return x * x;
}

double myAtan(double x, double y, int* flag)
{
    if(x==0 && y == 0) {
        *flag = FLAG_PHASOR0;
        return 0.;
    }
    *flag = 0;
    return atan2(y,x);
}


int metrology_algorithm(structTacData * tacData)
{
    int err = 0;
    long tel, diode, side;
    long i, idx;
    int old_flag_telescope[4][4][2];
    int old_flag_fiber_coupler[4][2];
    double x0, y0, x1, y1, xc, yc, delta_phase, sum_cos, sum_sin, new_phase, sq1, rms, tmp, speed;

    int flag, comb_flag, sum_flag;
    double sqflux, flux;
    double cos_delta_phase, sin_delta_phase;
    double *x0p, *y0p;
    
    structTacConfiguration *tacConfiguration = tacData->tacConfiguration;
    long sample_number = tacData->sample_number;
    
    long buffer_idx_smooth_tel =((sample_number - 1) % tacConfiguration->number_to_smooth_for_telescope);
    long buffer_idx_avg =((sample_number - 1) % tacConfiguration->number_to_average);
    long prev_idx_avg = buffer_idx_avg - 1;
    if(prev_idx_avg < 0) {
        prev_idx_avg += tacConfiguration->number_to_average;
    }
    
    long buffer_idx_rms = ((sample_number - 1) % tacConfiguration->number_for_rms);
    long buffer_idx_speed = ((sample_number - 1) % tacConfiguration->number_for_speed);
    
    long prev_idx_speed = buffer_idx_speed - 1;
    if(prev_idx_speed < 0) {
        prev_idx_speed += tacConfiguration->number_for_speed;
    }
    
    
    /* check if phasor did not move too far */
    
    if(sample_number > tacConfiguration->number_to_average) {
        for (tel = 0; tel < 4; tel++) {
            for (side = FT; side <= SC; side++) {
                tacData->flag_jump_fiber_coupler[tel][side] = 0;
                if((tacData->flux_fiber_coupler[tel][side] > tacConfiguration->min_allowed_flux_fiber_coupler[tel][side]) &&
                   (tacData->flux_fiber_coupler[tel][side] / tacConfiguration->rms_flux_fiber_coupler[tel][side] >
                    tacConfiguration->min_allowed_snr_fiber_coupler[tel][side]) &&
                   (tacData->snr_flux_fiber_coupler[tel][side] > tacConfiguration->min_allowed_snr_fiber_coupler[tel][side]) &&
                   (tacConfiguration->check_for_jumps == 1) ) {
                    x0 = tacData->volt_lockin_fiber_coupler[tel][side][COS];
                    y0 = tacData->volt_lockin_fiber_coupler[tel][side][SIN];
                    x1 = tacData->best_estimate_phasor_fiber_coupler[tel][side][COS];
                    y1 = tacData->best_estimate_phasor_fiber_coupler[tel][side][SIN];
                    xc = x1 * x0 + y1 * y0;
                    yc = x1 * y0 - y1 * x0;
                    if(fabs(myAtan(xc,yc,&flag)) > tacConfiguration->max_allowed_phase_difference) {
                        tacData->flag_jump_fiber_coupler[tel][side] = FLAG_JUMP;
                        x1 = tacData->buffer_volt_lockin_fiber_coupler[tel][side][COS][prev_idx_avg];
                        y1 = tacData->buffer_volt_lockin_fiber_coupler[tel][side][SIN][prev_idx_avg];
                        sq1 = sqrt(x1*x1+y1*y1);
                        if((sq1 > tacConfiguration->min_allowed_flux_fiber_coupler[tel][side]) &&
                           (sq1 / tacConfiguration->rms_flux_fiber_coupler[tel][side] >
                            tacConfiguration->min_allowed_snr_fiber_coupler[tel][side]) &&
                           (sq1 / tacData->rms_flux_fiber_coupler[tel][side] >
                            tacConfiguration->min_allowed_snr_fiber_coupler[tel][side]) ) {
                               xc = x1 * x0 + y1 * y0;
                               yc = x1 * y0 - y1 * x0;
                               if(fabs(myAtan(xc,yc,&flag)) < tacConfiguration->max_allowed_phase_difference) {
                                   tacData->flag_jump_fiber_coupler[tel][side] = 0;
                               }
                           }
                    }
                }
                if(tacConfiguration->number_to_smooth_for_telescope < 0) {
                    for(diode = 0; diode < 4; diode++) {
                        if((tacData->flux_telescope[tel][diode][side] > tacConfiguration->min_allowed_flux_telescope[tel][diode][side]) &&
                           (tacData->flux_telescope[tel][diode][side] / tacConfiguration->rms_flux_telescope[tel][diode][side] >
                            tacConfiguration->min_allowed_snr_telescope[tel][diode][side]) &&
                           (tacData->snr_flux_telescope[tel][diode][side] > tacConfiguration->min_allowed_snr_telescope[tel][diode][side]) &&
                           (tacConfiguration->check_for_jumps == 1) ) {
                            tacData->flag_jump_telescope[tel][diode][side] = 0;
                            x0 = tacData->volt_lockin_telescope[tel][diode][side][COS];
                            y0 = tacData->volt_lockin_telescope[tel][diode][side][SIN];
                            x1 = tacData->best_estimate_phasor_telescope[tel][diode][side][COS];
                            y1 = tacData->best_estimate_phasor_telescope[tel][diode][side][SIN];
                            xc = x1 * x0 + y1 * y0;
                            yc = x1 * y0 - y1 * x0;
                            
                            if(fabs(myAtan(xc,yc,&flag)) > tacConfiguration->max_allowed_phase_difference) {
                                tacData->flag_jump_telescope[tel][diode][side] = FLAG_JUMP;
                                x1 = tacData->buffer_volt_lockin_telescope[tel][diode][side][COS][prev_idx_avg];
                                y1 = tacData->buffer_volt_lockin_telescope[tel][diode][side][SIN][prev_idx_avg];
                                sq1 = sqrt(x1*x1+y1*y1);
                                if((sq1 > tacConfiguration->min_allowed_flux_telescope[tel][diode][side]) &&
                                   (sq1 / tacConfiguration->rms_flux_telescope[tel][diode][side] >
                                    tacConfiguration->min_allowed_snr_telescope[tel][diode][side]) &&
                                   (sq1 / tacData->rms_flux_telescope[tel][diode][side] >
                                    tacConfiguration->min_allowed_snr_telescope[tel][diode][side])) {
                                       xc = x1 * x0 + y1 * y0;
                                       yc = x1 * y0 - y1 * x0;
                                       if(fabs(myAtan(xc,yc,&flag)) < tacConfiguration->max_allowed_phase_difference) {
                                           tacData->flag_jump_telescope[tel][diode][side] = 0;
                                       }
                                   }
                            }
                        }
                    }
                } /* end if num_smooth_tel < 0 */
            }
        } /* end loop tel */
    }
    
    /* raw phase calculation */
    
    for (tel = 0; tel < 4; tel++) {
        for (side = FT; side <= SC; side++) {
            x0 = tacData->volt_lockin_fiber_coupler[tel][side][COS];
            y0 = tacData->volt_lockin_fiber_coupler[tel][side][SIN];
            x0 += tacData->best_estimate_phasor_fiber_coupler[tel][side][COS];
            y0 += tacData->best_estimate_phasor_fiber_coupler[tel][side][SIN];
            tacData->raw_phase_fiber_coupler[tel][side] = myAtan(x0,y0,&(tacData->flag_phasor0_fiber_coupler[tel][side]));
            
            if(tacConfiguration->number_to_smooth_for_telescope < 0) {
                for(diode = 0; diode < 4; diode++) {
                    x0 = tacData->volt_lockin_telescope[tel][diode][side][COS];
                    y0 = tacData->volt_lockin_telescope[tel][diode][side][SIN];
                    x0 += tacData->best_estimate_phasor_telescope[tel][diode][side][COS];
                    y0 += tacData->best_estimate_phasor_telescope[tel][diode][side][SIN];
                    tacData->raw_phase_telescope[tel][diode][side] = myAtan(x0,y0,&(tacData->flag_phasor0_telescope[tel][diode][side]));
                }
            }
        }
    }
    
    /* unwrapping of phase */
    
    for (tel = 0; tel < 4; tel++) {
        for (side = FT; side <= SC; side++) {
            if(sample_number == 1) {
                tacData->prev_phase_fiber_coupler[tel][side] = tacData->raw_phase_fiber_coupler[tel][side];
            } else {
                tacData->prev_phase_fiber_coupler[tel][side] = tacData->unwrapped_phase_fiber_coupler[tel][side];
            }
            
            tacData->flag_unwrap_fiber_coupler[tel][side] = metrology_unwrap(tacData->raw_phase_fiber_coupler[tel][side],
                                                                   tacData->prev_phase_fiber_coupler[tel][side], tacConfiguration->max_allowed_phase_difference,
                                                                   &tacData->unwrapped_phase_fiber_coupler[tel][side]);
            
            tacData->buffer_unwrapped_phase_fiber_coupler[tel][side][buffer_idx_avg] = tacData->unwrapped_phase_fiber_coupler[tel][side];
            
            if(sample_number < tacConfiguration->number_to_average) {
                tacData->flag_unwrap_fiber_coupler[tel][side] = 0;
            }
            
            if(tacConfiguration->number_to_smooth_for_telescope < 0) {
                for(diode = 0; diode < 4; diode++) {
                    if(sample_number == 1) {
                        tacData->prev_phase_telescope[tel][diode][side] = tacData->raw_phase_telescope[tel][diode][side];
                    } else {
                        tacData->prev_phase_telescope[tel][diode][side] = tacData->unwrapped_phase_telescope[tel][diode][side];
                    }
                    
                    tacData->flag_unwrap_telescope[tel][diode][side] = metrology_unwrap(tacData->raw_phase_telescope[tel][diode][side],
                                                                              tacData->prev_phase_telescope[tel][diode][side], tacConfiguration->max_allowed_phase_difference,
                                                                              &tacData->unwrapped_phase_telescope[tel][diode][side]);
                    
                    if(sample_number < tacConfiguration->number_to_average) {
                        tacData->flag_unwrap_telescope[tel][diode][side] = 0;
                    }
                }
            }
        }
    }
    
    /* calculate next best estimate for phasor */
    
    for (tel = 0; tel < 4; tel++) {
        for (side = FT; side <= SC; side++) {
            if(tacConfiguration->number_to_average > 1) {
                tacData->delta_phasor_fiber_coupler[tel][side][COS] -=
                tacData->buffer_delta_phasor_fiber_coupler[tel][side][COS][buffer_idx_avg];
                tacData->delta_phasor_fiber_coupler[tel][side][SIN] -=
                tacData->buffer_delta_phasor_fiber_coupler[tel][side][SIN][buffer_idx_avg];
                
                x0 = tacData->buffer_volt_lockin_fiber_coupler[tel][side][COS][prev_idx_avg];
                y0 = tacData->buffer_volt_lockin_fiber_coupler[tel][side][SIN][prev_idx_avg];
                x1 = tacData->volt_lockin_fiber_coupler[tel][side][COS];
                y1 = tacData->volt_lockin_fiber_coupler[tel][side][SIN];
                xc = x1 * x0 + y1 * y0;
                yc = x1 * y0 - y1 * x0;
                tacData->buffer_delta_phasor_fiber_coupler[tel][side][COS][buffer_idx_avg] = xc;
                tacData->buffer_delta_phasor_fiber_coupler[tel][side][SIN][buffer_idx_avg] = yc;
                
                tacData->delta_phasor_fiber_coupler[tel][side][COS] += xc;
                tacData->delta_phasor_fiber_coupler[tel][side][SIN] += yc;
                
                tacData->sum_delta_phase_fiber_coupler[tel][side] -= tacData->buffer_delta_phase_fiber_coupler[tel][side][buffer_idx_avg];
                tacData->sum_sq_delta_phase_fiber_coupler[tel][side] -= tacData->buffer_sq_delta_phase_fiber_coupler[tel][side][buffer_idx_avg];
                
                tacData->buffer_delta_phase_fiber_coupler[tel][side][buffer_idx_avg] = myAtan(yc, xc, &flag);
                tacData->buffer_sq_delta_phase_fiber_coupler[tel][side][buffer_idx_avg] = metrology_sq(tacData->buffer_delta_phase_fiber_coupler[tel][side][buffer_idx_avg]);
                
                tacData->sum_delta_phase_fiber_coupler[tel][side] += tacData->buffer_delta_phase_fiber_coupler[tel][side][buffer_idx_avg];
                tacData->sum_sq_delta_phase_fiber_coupler[tel][side] += tacData->buffer_sq_delta_phase_fiber_coupler[tel][side][buffer_idx_avg];
                
                rms = sqrt(tacData->sum_sq_delta_phase_fiber_coupler[tel][side]/tacConfiguration->number_to_average - metrology_sq(tacData->sum_delta_phase_fiber_coupler[tel][side]/tacConfiguration->number_to_average));
                
                tacData->sum_speed_fiber_coupler[tel][side] -= tacData->buffer_speed_fiber_coupler[tel][side][buffer_idx_speed] * pow(tacConfiguration->decrement_factor_speed, tacConfiguration->number_for_speed - 1);
                
                if(tacData->delta_phasor_fiber_coupler[tel][side][COS] == 0 && tacData->delta_phasor_fiber_coupler[tel][side][SIN] == 0) {
                    speed = tacData->buffer_speed_fiber_coupler[tel][side][prev_idx_speed];
                } else {
                    speed = -myAtan(tacData->delta_phasor_fiber_coupler[tel][side][COS],tacData->delta_phasor_fiber_coupler[tel][side][SIN],&flag);
                }
                
                if(tacConfiguration->calc_phase_speed == 0) speed = 0.0;
                
                if(tacConfiguration->sigma_clip_speed > 0) {
                    if(fabs(speed) * tacConfiguration->number_to_average < tacConfiguration->sigma_clip_speed * rms) speed = 0.0;
                }
                
                tacData->buffer_speed_fiber_coupler[tel][side][buffer_idx_speed] = speed;
                tacData->sum_speed_fiber_coupler[tel][side] *= tacConfiguration->decrement_factor_speed;
                tacData->sum_speed_fiber_coupler[tel][side] += speed;
                
                delta_phase = tacData->sum_speed_fiber_coupler[tel][side] / tacConfiguration->norm_speed;
                
            } else {
                delta_phase = 0.;
            }
            
            tacData->buffer_used_speed_fiber_coupler[tel][side][buffer_idx_avg] = delta_phase;
            
            tacData->flag_speed_fiber_coupler[tel][side] = 0;
            if(sample_number > tacConfiguration->number_to_average) {
                if(fabs(delta_phase) > tacConfiguration->max_allowed_phase_speed) {
                    tacData->flag_speed_fiber_coupler[tel][side] = FLAG_SPEED;
                }
            }
            
            cos_delta_phase = cos(delta_phase);
            sin_delta_phase = sin(delta_phase);
            
            x1 = 1;
            y1 = 0;
            
            idx = buffer_idx_avg;
            sum_cos = 0;
            sum_sin = 0;
            x0p = tacData->buffer_volt_lockin_fiber_coupler[tel][side][COS];
            y0p = tacData->buffer_volt_lockin_fiber_coupler[tel][side][SIN];
            
            for(i = 0; i < tacConfiguration->number_to_average; i++) {
                x0 = x0p[idx];
                y0 = y0p[idx];
                sum_cos += x0 * x1 - y0 * y1;
                sum_sin += x0 * y1 + x1 * y0;
                idx -= 1;
                if(idx < 0) {
                    idx += tacConfiguration->number_to_average;
                }
                tmp = x1 * cos_delta_phase - y1 * sin_delta_phase;
                y1 = x1 * sin_delta_phase + y1 * cos_delta_phase;
                x1 = tmp;
            }
            
            tacData->best_estimate_flux_fiber_coupler[tel][side] = sqrt(metrology_sq(sum_cos)+metrology_sq(sum_sin));
            new_phase = myAtan(sum_cos, sum_sin, &flag) + delta_phase;
            
            tacData->best_estimate_phasor_fiber_coupler[tel][side][COS] = tacData->best_estimate_flux_fiber_coupler[tel][side] * cos(new_phase);
            tacData->best_estimate_phasor_fiber_coupler[tel][side][SIN] = tacData->best_estimate_flux_fiber_coupler[tel][side] * sin(new_phase);
            tacData->best_estimate_flux_fiber_coupler[tel][side] /= tacConfiguration->number_to_average;
            tacData->buffer_best_estimate_phasor_fiber_coupler[tel][side][COS][buffer_idx_avg] = tacData->best_estimate_phasor_fiber_coupler[tel][side][COS]/ tacConfiguration->number_to_average;
            tacData->buffer_best_estimate_phasor_fiber_coupler[tel][side][SIN][buffer_idx_avg] = tacData->best_estimate_phasor_fiber_coupler[tel][side][SIN]/ tacConfiguration->number_to_average;
            
            tacData->flag_flux_fiber_coupler[tel][side] = 0;
            tacData->flag_snr_fiber_coupler[tel][side] = 0;
            if (sample_number > tacConfiguration->number_to_average) {
                if (tacData->best_estimate_flux_fiber_coupler[tel][side] < tacConfiguration->min_allowed_flux_fiber_coupler[tel][side]) {
                    tacData->flag_flux_fiber_coupler[tel][side] = FLAG_FLUX;
                }
                if (tacData->best_estimate_flux_fiber_coupler[tel][side] / tacConfiguration->rms_flux_fiber_coupler[tel][side] <
                    tacConfiguration->min_allowed_snr_fiber_coupler[tel][side]) {
                    tacData->flag_flux_fiber_coupler[tel][side] = FLAG_FLUX;
                }
            }
            
            if(tacConfiguration->number_to_smooth_for_telescope < 0) {
                for(diode = 0; diode < 4; diode++) {
                    if(tacConfiguration->number_to_average > 1) {
                        tacData->delta_phasor_telescope[tel][diode][side][COS] -=
                        tacData->buffer_delta_phasor_telescope[tel][diode][side][COS][buffer_idx_avg];
                        tacData->delta_phasor_telescope[tel][diode][side][SIN] -=
                        tacData->buffer_delta_phasor_telescope[tel][diode][side][SIN][buffer_idx_avg];
                        
                        x0 = tacData->buffer_volt_lockin_telescope[tel][diode][side][COS][prev_idx_avg];
                        y0 = tacData->buffer_volt_lockin_telescope[tel][diode][side][SIN][prev_idx_avg];
                        x1 = tacData->volt_lockin_telescope[tel][diode][side][COS];
                        y1 = tacData->volt_lockin_telescope[tel][diode][side][SIN];
                        xc = x1 * x0 + y1 * y0;
                        yc = x1 * y0 - y1 * x0;
                        tacData->buffer_delta_phasor_telescope[tel][diode][side][COS][buffer_idx_avg] = xc;
                        tacData->buffer_delta_phasor_telescope[tel][diode][side][SIN][buffer_idx_avg] = yc;
                        
                        tacData->delta_phasor_telescope[tel][diode][side][COS] += xc;
                        tacData->delta_phasor_telescope[tel][diode][side][SIN] += yc;
                        
                        tacData->sum_delta_phase_telescope[tel][diode][side] -= tacData->buffer_delta_phase_telescope[tel][diode][side][buffer_idx_avg];
                        tacData->sum_sq_delta_phase_telescope[tel][diode][side] -= tacData->buffer_sq_delta_phase_telescope[tel][diode][side][buffer_idx_avg];
                        
                        tacData->buffer_delta_phase_telescope[tel][diode][side][buffer_idx_avg] = myAtan(yc, xc, &flag);
                        tacData->buffer_sq_delta_phase_telescope[tel][diode][side][buffer_idx_avg] = metrology_sq(tacData->buffer_delta_phase_telescope[tel][diode][side][buffer_idx_avg]);
                        
                        tacData->sum_delta_phase_telescope[tel][diode][side] += tacData->buffer_delta_phase_telescope[tel][diode][side][buffer_idx_avg];
                        tacData->sum_sq_delta_phase_telescope[tel][diode][side] += tacData->buffer_sq_delta_phase_telescope[tel][diode][side][buffer_idx_avg];
                        
                        rms = sqrt(tacData->sum_sq_delta_phase_telescope[tel][diode][side]/tacConfiguration->number_to_average - metrology_sq(tacData->sum_delta_phase_telescope[tel][diode][side]/tacConfiguration->number_to_average));
                        
                        tacData->sum_speed_telescope[tel][diode][side] -= tacData->buffer_speed_telescope[tel][diode][side][buffer_idx_speed] * pow(tacConfiguration->decrement_factor_speed, tacConfiguration->number_for_speed - 1);
                        
                        if(tacData->delta_phasor_telescope[tel][diode][side][COS] == 0 && tacData->delta_phasor_telescope[tel][diode][side][SIN] == 0) {
                            speed = tacData->buffer_speed_telescope[tel][diode][side][prev_idx_speed];
                        } else {
                            speed = -myAtan(tacData->delta_phasor_telescope[tel][diode][side][COS],tacData->delta_phasor_telescope[tel][diode][side][SIN],&flag);
                        }
                        
                        if(tacConfiguration->calc_phase_speed == 0) speed = 0.0;
                        
                        if(tacConfiguration->sigma_clip_speed > 0) {
                            if(fabs(speed) * tacConfiguration->number_to_average < tacConfiguration->sigma_clip_speed * rms) speed = 0.0;
                        }
                        
                        tacData->buffer_speed_telescope[tel][diode][side][buffer_idx_speed] = speed;
                        tacData->sum_speed_telescope[tel][diode][side] *= tacConfiguration->decrement_factor_speed;
                        tacData->sum_speed_telescope[tel][diode][side] += speed;
                        
                        delta_phase = tacData->sum_speed_telescope[tel][diode][side] / tacConfiguration->norm_speed;
                        
                    } else {
                        delta_phase = 0.;
                    }
                    
                    tacData->flag_speed_telescope[tel][diode][side] = 0;
                    if(sample_number > tacConfiguration->number_to_average) {
                        if(fabs(delta_phase) > tacConfiguration->max_allowed_phase_difference) {
                            tacData->flag_speed_telescope[tel][diode][side] = FLAG_SPEED;
                        }
                    }
                    
                    idx = buffer_idx_avg;
                    sum_cos = 0;
                    sum_sin = 0;
                    
                    cos_delta_phase = cos(delta_phase);
                    sin_delta_phase = sin(delta_phase);
                    
                    x1 = 1;
                    y1 = 0;
                    x0p = tacData->buffer_volt_lockin_telescope[tel][diode][side][COS];
                    y0p = tacData->buffer_volt_lockin_telescope[tel][diode][side][SIN];
                    
                    for(i = 0; i < tacConfiguration->number_to_average; i++) {
                        x0 = x0p[idx];
                        y0 = y0p[idx];
                        sum_cos += x0 * x1 - y0 * y1;
                        sum_sin += x0 * y1 + x1 * y0;
                        idx -= 1;
                        if(idx < 0) {
                            idx += tacConfiguration->number_to_average;
                        }
                        tmp = x1 * cos_delta_phase - y1 * sin_delta_phase;
                        y1 = x1 * sin_delta_phase + y1 * cos_delta_phase;
                        x1 = tmp;
                    }
                    
                    tacData->best_estimate_flux_telescope[tel][diode][side] = sqrt(metrology_sq(sum_cos)+metrology_sq(sum_sin));
                    new_phase = myAtan(sum_cos, sum_sin, &flag) + delta_phase;
                    tacData->best_estimate_phasor_telescope[tel][diode][side][COS] = tacData->best_estimate_flux_telescope[tel][diode][side] * cos(new_phase);
                    tacData->best_estimate_phasor_telescope[tel][diode][side][SIN] = tacData->best_estimate_flux_telescope[tel][diode][side] * sin(new_phase);
                    tacData->best_estimate_flux_telescope[tel][diode][side] /= tacConfiguration->number_to_average;
                    
                    tacData->flag_flux_telescope[tel][diode][side] = 0;
                    tacData->flag_snr_telescope[tel][diode][side] = 0;
                    if (sample_number > tacConfiguration->number_to_average) {
                        if (tacData->best_estimate_flux_telescope[tel][diode][side] < tacConfiguration->min_allowed_flux_telescope[tel][diode][side]) {
                            tacData->flag_flux_telescope[tel][diode][side] = FLAG_FLUX;
                        }
                        if (tacData->best_estimate_flux_telescope[tel][diode][side] / tacConfiguration->rms_flux_telescope[tel][diode][side] <
                            tacConfiguration->min_allowed_snr_telescope[tel][diode][side]) {
                            tacData->flag_flux_telescope[tel][diode][side] = FLAG_FLUX;
                        }
                    }
                }    /* end loop diode */
            }  /* end if nsmooth_tel < 0*/
        } /* end loop side */
    }  /* end loop telescope */
    
    /* calclation of rms of flux and SNR check */
    
    for (tel = 0; tel < 4; tel++) {
        for (side = FT; side <= SC; side++) {
            if(tacConfiguration->number_to_smooth_for_telescope < 0) {
                for(diode = 0; diode < 4; diode++) {
                    tacData->rms_flux_telescope[tel][diode][side] =
                    sqrt(tacData->sum_sq_flux_telescope[tel][diode][side]/tacConfiguration->number_for_rms
                         - metrology_sq(tacData->sum_flux_telescope[tel][diode][side]/tacConfiguration->number_for_rms));
                    tacData->snr_flux_telescope[tel][diode][side] = tacData->best_estimate_flux_telescope[tel][diode][side] /
                    tacData->rms_flux_telescope[tel][diode][side];
                    if (sample_number > tacConfiguration->number_for_rms) {
                        if (tacData->snr_flux_telescope[tel][diode][side]  <
                            tacConfiguration->min_allowed_snr_telescope[tel][diode][side]) {
                            tacData->flag_snr_telescope[tel][diode][side] = FLAG_SNR;
                        }
                    }
                }
            }
            tacData->rms_flux_fiber_coupler[tel][side] =
            sqrt(tacData->sum_sq_flux_fiber_coupler[tel][side]/tacConfiguration->number_for_rms
                 - metrology_sq(tacData->sum_flux_fiber_coupler[tel][side]/tacConfiguration->number_for_rms));
            tacData->snr_flux_fiber_coupler[tel][side] = tacData->best_estimate_flux_fiber_coupler[tel][side] /
            tacData->rms_flux_fiber_coupler[tel][side];
            
            if (sample_number > tacConfiguration->number_for_rms) {
                if (tacData->snr_flux_fiber_coupler[tel][side] <
                    tacConfiguration->min_allowed_snr_fiber_coupler[tel][side]) {
                    tacData->flag_snr_fiber_coupler[tel][side] = FLAG_SNR;
                }
            }
        }
    }
    
    /* phase calculation for difference phasor Tel - FC */
    
    if(tacConfiguration->number_to_smooth_for_telescope > 0) {
        for (tel = 0; tel < 4; tel++) {
            for (side = FT; side <= SC; side++) {
                x0 = tacData->volt_lockin_fiber_coupler[tel][side][COS];
                y0 = tacData->volt_lockin_fiber_coupler[tel][side][SIN];
                for (diode = 0; diode < 4; diode++) {
                    x1 = tacData->volt_lockin_telescope[tel][diode][side][COS];
                    y1 = tacData->volt_lockin_telescope[tel][diode][side][SIN];
                    xc = x0 * x1 + y0 * y1;
                    yc = x0 * y1 - y0 * x1;
                    
                    /* maintain a buffer of x and y values to be able to average */
                    tacData->sum_phasor_diff[tel][diode][side][COS] -= tacData->buffer_phasor_diff[tel][diode][side][COS][buffer_idx_smooth_tel];
                    tacData->sum_phasor_diff[tel][diode][side][SIN] -= tacData->buffer_phasor_diff[tel][diode][side][SIN][buffer_idx_smooth_tel];
                    tacData->buffer_phasor_diff[tel][diode][side][COS][buffer_idx_smooth_tel] = xc;
                    tacData->buffer_phasor_diff[tel][diode][side][SIN][buffer_idx_smooth_tel] = yc;
                    tacData->sum_phasor_diff[tel][diode][side][COS] += xc;
                    tacData->sum_phasor_diff[tel][diode][side][SIN] += yc;
                    /* the flux of the difference is calculated as the sum, not the average ! */
                    sqflux = metrology_sq(tacData->sum_phasor_diff[tel][diode][side][COS])+
                    metrology_sq(tacData->sum_phasor_diff[tel][diode][side][SIN]);
                    flux = sqrt(sqflux);
                    tacData->flux_diff[tel][diode][side] = flux;
                    
                    tacData->raw_phase_diff[tel][diode][side] = myAtan(tacData->sum_phasor_diff[tel][diode][side][COS],
                                                                       tacData->sum_phasor_diff[tel][diode][side][SIN],
                                                                       &(tacData->flag_phasor0_diff[tel][diode][side]));
                    
                    tacData->sum_flux_diff[tel][diode][side] -= tacData->buffer_flux_diff[tel][diode][side][buffer_idx_rms];
                    tacData->sum_sq_flux_diff[tel][diode][side] -= tacData->buffer_sq_flux_diff[tel][diode][side][buffer_idx_rms];
                    tacData->buffer_flux_diff[tel][diode][side][buffer_idx_rms] = flux;
                    tacData->buffer_sq_flux_diff[tel][diode][side][buffer_idx_rms] = sqflux;
                    tacData->sum_flux_diff[tel][diode][side] += flux;
                    tacData->sum_sq_flux_diff[tel][diode][side] += sqflux;
                    
                    tacData->rms_flux_diff[tel][diode][side] =
                    sqrt(tacData->sum_sq_flux_diff[tel][diode][side]/tacConfiguration->number_for_rms
                         - metrology_sq(tacData->sum_flux_diff[tel][diode][side]/tacConfiguration->number_for_rms));
                    tacData->snr_flux_diff[tel][diode][side] = tacData->flux_diff[tel][diode][side] /
                    tacData->rms_flux_diff[tel][diode][side];
                    
                    tacData->flag_flux_diff[tel][diode][side] = 0;
                    tacData->flag_snr_diff[tel][diode][side] = 0;
                    if (sample_number > tacConfiguration->number_to_smooth_for_telescope) {
                        if (tacData->flux_diff[tel][diode][side] < tacConfiguration->min_allowed_flux_telescope[tel][diode][side]) {
                            tacData->flag_flux_diff[tel][diode][side] = FLAG_FLUX;
                        }
                        if (tacData->flux_diff[tel][diode][side] / tacConfiguration->rms_flux_telescope[tel][diode][side] <
                            tacConfiguration->min_allowed_snr_diff[tel][diode][side]) {
                            tacData->flag_flux_diff[tel][diode][side] = FLAG_FLUX;
                        }
                        if(sample_number > tacConfiguration->number_for_rms) {
                            if (tacData->snr_flux_diff[tel][diode][side] <
                                tacConfiguration->min_allowed_snr_diff[tel][diode][side]) {
                                tacData->flag_snr_diff[tel][diode][side] = FLAG_SNR;
                            }
                        }
                    }
                } /* end loop diode */
            } /* end loop side */
        } /* end loop tel */
        
        /* unwrapping of difference phase */
        for (tel = 0; tel < 4; tel++) {
            for (side = FT; side <= SC; side++) {
                for (diode = 0; diode < 4; diode++) {
                    
                    if(sample_number <= tacConfiguration->number_to_smooth_for_telescope) {
                        tacData->prev_phase_diff[tel][diode][side] = tacData->raw_phase_diff[tel][diode][side];
                    } else {
                        tacData->prev_phase_diff[tel][diode][side] = tacData->unwrapped_phase_diff[tel][diode][side];
                    }
                    
                    tacData->flag_unwrap_diff[tel][diode][side] = metrology_unwrap(tacData->raw_phase_diff[tel][diode][side],
                                                                         tacData->prev_phase_diff[tel][diode][side], tacConfiguration->max_allowed_phase_difference,
                                                                         &tacData->unwrapped_phase_diff[tel][diode][side]);
                    
                } /* end loop diode */
            } /* end loop side */
        } /* end loop tel */
    }  /* end if num_smooth_tel > 0 */
    
    
    /* filling the telescope flags and phase from difference phasor */
    if(tacConfiguration->number_to_smooth_for_telescope > 0) {
        for (tel = 0; tel < 4; tel++) {
            for (side = FT; side <= SC; side++) {
                for (diode = 0; diode < 4; diode++) {
                    tacData->flag_unwrap_telescope[tel][diode][side] = tacData->flag_unwrap_diff[tel][diode][side];
                    tacData->flag_flux_telescope[tel][diode][side] = tacData->flag_flux_diff[tel][diode][side];
                    tacData->flag_snr_telescope[tel][diode][side] = tacData->flag_snr_diff[tel][diode][side];
                    tacData->flag_phasor0_telescope[tel][diode][side] = tacData->flag_phasor0_diff[tel][diode][side];
                    tacData->unwrapped_phase_telescope[tel][diode][side] = tacData->unwrapped_phase_fiber_coupler[tel][side] + tacData->unwrapped_phase_diff[tel][diode][side];
                }
            }
        }
    } /* end if num_smooth_tel > 0 */
    
    
    
    /* Calculate the combined flags */
    for (tel = 0; tel < 4; tel++) {
        for (side = FT; side <= SC; side++) {
            comb_flag = 0;
            sum_flag = 0;
            for (diode = 0; diode < 4; diode++) {
                old_flag_telescope[tel][diode][side] = tacData->total_flag_telescope[tel][diode][side];
                tacData->total_flag_telescope[tel][diode][side] =
                tacData->flag_flux_telescope[tel][diode][side]
                + tacData->flag_unwrap_telescope[tel][diode][side]
                + tacData->flag_jump_telescope[tel][diode][side]
                + tacData->flag_speed_telescope[tel][diode][side]
                + tacData->flag_phasor0_telescope[tel][diode][side]
                + tacData->flag_snr_telescope[tel][diode][side]
                + tacData->flag_volt_telescope[tel][diode][side];
                comb_flag |= tacData->total_flag_telescope[tel][diode][side];
                if( tacData->total_flag_telescope[tel][diode][side] != 0) sum_flag++;
            } /* end loop diode */
            
            old_flag_fiber_coupler[tel][side] = tacData->total_flag_fiber_coupler[tel][side];
            tacData->total_flag_fiber_coupler[tel][side] =
            tacData->flag_flux_fiber_coupler[tel][side]
            + tacData->flag_unwrap_fiber_coupler[tel][side]
            + tacData->flag_jump_fiber_coupler[tel][side]
            + tacData->flag_speed_fiber_coupler[tel][side]
            + tacData->flag_phasor0_fiber_coupler[tel][side]
            + tacData->flag_snr_fiber_coupler[tel][side]
            + tacData->flag_volt_fiber_coupler[tel][side];
        } /* end loop side */
    } /* end loop tel */
    
    /* If there was an error, freeze + speed predict */
    for (tel = 0; tel < 4; tel++) {
        for (side = FT; side <= SC; side++) {
            if(tacData->total_flag_fiber_coupler[tel][side] != 0) {
                tacData->freeze_count_fiber_coupler[tel][side] += 1;
                
                if(tacData->freeze_count_fiber_coupler[tel][side] == 1) {
                    
                    
                    tacData->freeze_phasor_fiber_coupler[tel][side][COS] = tacData->buffer_best_estimate_phasor_fiber_coupler[tel][side][COS][prev_idx_avg];
                    tacData->freeze_phasor_fiber_coupler[tel][side][SIN] = tacData->buffer_best_estimate_phasor_fiber_coupler[tel][side][SIN][prev_idx_avg];
                    tacData->freeze_speed_fiber_coupler[tel][side] = tacData->buffer_used_speed_fiber_coupler[tel][side][prev_idx_avg];
                    
                    if(tacConfiguration->calc_phase_speed == 0) tacData->freeze_speed_fiber_coupler[tel][side] = 0.0;
                    if(tacData->flag_speed_fiber_coupler[tel][side] !=0) tacData->freeze_speed_fiber_coupler[tel][side] = 0.0;
                }
                
                new_phase = myAtan(tacData->freeze_phasor_fiber_coupler[tel][side][COS], tacData->freeze_phasor_fiber_coupler[tel][side][SIN], &(tacData->flag_phasor0_fiber_coupler[tel][side])) + tacData->freeze_speed_fiber_coupler[tel][side];
                tmp = sqrt(metrology_sq(tacData->freeze_phasor_fiber_coupler[tel][side][COS])+metrology_sq(tacData->freeze_phasor_fiber_coupler[tel][side][SIN]));
                tacData->freeze_phasor_fiber_coupler[tel][side][COS] = tmp * cos(new_phase);
                tacData->freeze_phasor_fiber_coupler[tel][side][SIN] = tmp * sin(new_phase);
                
                tmp = myAtan(tacData->freeze_phasor_fiber_coupler[tel][side][COS], tacData->freeze_phasor_fiber_coupler[tel][side][SIN], &(tacData->flag_phasor0_fiber_coupler[tel][side]));
                
            } else {
                /* this should be the nominal case, no error occurred and the counter is 0 */
                
                tacData->freeze_count_fiber_coupler[tel][side] = 0;
                
                tmp = myAtan(cos(tacData->unwrapped_phase_fiber_coupler[tel][side]), sin(tacData->unwrapped_phase_fiber_coupler[tel][side]), &(tacData-> flag_phasor0_fiber_coupler[tel][side]));
                
                if(tacData->flag_unwrap_fiber_coupler[tel][side] !=0 ) {
                    printf("After freezing: speed predictor too far from where phase is, pred: %f, phase: %f\n",
                           tacData->freeze_unwrapped_phase_fiber_coupler[tel][side],tacData->unwrapped_phase_fiber_coupler[tel][side]);
                }
            }
            
            if( metrology_unwrap(tmp,tacData->used_unwrapped_phase_fiber_coupler[tel][side],
                       tacConfiguration->max_allowed_phase_difference,
                       &tacData->used_unwrapped_phase_fiber_coupler[tel][side]) != 0) {
                
                if(sample_number >= tacConfiguration->number_to_average) {
                    
                    tacData->total_flag_fiber_coupler[tel][side] += FLAG_LOST_PHASE ;
                    
                }
            }
            
        } /* loop side */
    } /* loop tel */
    
    
    /* difference FT - SC; and calculate mean per telescope */
    
    for (tel = 0; tel < 4; tel++) {
        tacData->mean_phase_telescope[tel] = 0.0;
        for (diode = 0; diode < 4; diode++) {
            tacData->delta_phase_telescope[tel][diode] = tacConfiguration->sign_of_phase * (tacData->unwrapped_phase_telescope[tel][diode][SC] - tacData->unwrapped_phase_telescope[tel][diode][FT]);
            tacData->mean_phase_telescope[tel] += tacData->delta_phase_telescope[tel][diode];
        }
        tacData->delta_phase_fiber_coupler[tel] = tacConfiguration->sign_of_phase * (tacData->used_unwrapped_phase_fiber_coupler[tel][SC] - tacData->used_unwrapped_phase_fiber_coupler[tel][FT]);
        tacData->mean_phase_telescope[tel] /= 4.0;
    }
    
    
    for (tel = 0; tel < 4; tel++) {
        for (diode = 0; diode < 4; diode++) {
           if ( tacData->total_flag_telescope[tel][diode][FT] == 0 && tacData->total_flag_telescope[tel][diode][SC] == 0) {
              tacData->opl_telescope_diode[tel][diode] = tacConfiguration->opl_base * tacData->delta_phase_telescope[tel][diode] - tacData->opl_zero_telescope_diode[tel][diode];
           }
        }
        tacData->opl_fiber_coupler[tel] = tacConfiguration->opl_base * tacData->delta_phase_fiber_coupler[tel]
        - tacData->opl_zero_fiber_coupler[tel];
        
    }
    
    return err;
}

/*-----------------------------------------------------------------------------
                             Functions code
 -----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/**
 * @brief Create the VIS_MET table
 *
 * @param metrology_table: input METROLOGY table
 * @param header:     corresponding HEADER
 *
 * @return a new VIS_MET table
 *
 * Create the VIS_MET table from the METROLOGY table.
 * The table has nsample*ntel rows, so that the
 * measurements for beam 0 are in rows 0,5,10,15... the ones for beam 1
 * are in 1,6,11,16... that is following the same format as OI_FLUX.
 * This function only create the TIME column.
 */
/*----------------------------------------------------------------------------*/

cpl_table * gravi_metrology_create (cpl_table * metrology_table,
                                    cpl_propertylist * header)
{
    gravi_msg_function_start(1);
	cpl_ensure (metrology_table, CPL_ERROR_NULL_INPUT, NULL);
	cpl_ensure (header,  CPL_ERROR_NULL_INPUT, NULL);

    /* Read MJD of PRC */
    double mjd0 = gravi_convert_to_mjd (gravi_pfits_get_start_prcacq (header));
		
	/* Create the output table for VIS_MET */
	int ntel = 4;
	cpl_size nrow = cpl_table_get_nrow (metrology_table);
	cpl_table * vismet_table = cpl_table_new (nrow * ntel);

    /* Create the TIME column */
	cpl_table_new_column (vismet_table, "TIME", CPL_TYPE_INT);
	cpl_table_set_column_unit (vismet_table, "TIME", "usec");
    
    /* Create the MJD column */
	cpl_table_new_column (vismet_table, "MJD", CPL_TYPE_DOUBLE);
	cpl_table_set_column_unit (vismet_table, "MJD", "day");

    /* Fill the TIME and MJD column */
    for (cpl_size row = 0; row < nrow; row++) {
        int time_met = cpl_table_get_int (metrology_table, "TIME", row, NULL);
        double mjd_met = time_met / 86400.E6 + mjd0;
        for (cpl_size tel = 0; tel < ntel; tel++) {
            cpl_table_set_int (vismet_table, "TIME", row*ntel+tel, time_met);
            cpl_table_set_double (vismet_table, "MJD",  row*ntel+tel, mjd_met);
        }
    }

    /* Return */
    gravi_msg_function_exit(1);
    return vismet_table;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Fill the VIS_MET table with the OPD_PUPIL column
 *
 * @param visacq_table: input OI_VIS_ACQ table
 * @param vismet_table: output OI_VIS_MET table
 * @param delay:        delay in [s] between TIME in ACQ and the correction
 *                      seen by the metrology (TIME in OI_VIS_MET)
 * @param header:       corresponding HEADER
 *
 * Fill the OI_VIS_MET table from the OPD_PUPIL column computed from
 * the OPD_PUPIL column of the OI_VIS_ACQ table.
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_metrology_acq (cpl_table * visacq_table,
                                    cpl_table * vismet_table,
                                    double delay,
                                    cpl_propertylist * header)
{
    gravi_msg_function_start(1);
	cpl_ensure_code (visacq_table, CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (vismet_table, CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (header,  CPL_ERROR_NULL_INPUT);

    cpl_msg_info (cpl_func,"Use acq-correction-delay = %.3f [s]",delay);

    /* Get size */
	int ntel = 4;
	cpl_size nrow_met = cpl_table_get_nrow (vismet_table) / ntel;
	cpl_size nrow_acq = cpl_table_get_nrow (visacq_table) / ntel;

    /* Create a temporary table */
    cpl_table * visacq_tmp = cpl_table_new (nrow_acq * ntel);
    cpl_table_duplicate_column (visacq_tmp, "TIME", visacq_table, "TIME");
    cpl_table_duplicate_column (visacq_tmp, "OPD_PUPIL", visacq_table, "OPD_PUPIL");
    cpl_table_duplicate_column (visacq_tmp, "PUPIL_NSPOT", visacq_table, "PUPIL_NSPOT");

    /* Move the time of the temporary table 
     * delta is in [s] */
    cpl_table_add_scalar (visacq_tmp, "TIME", 1e6 * delay);
    
    /* Get the ACQ DIT in [us] */
    double dit_acq = gravi_pfits_get_dit_acqcam (header) * 1e6;
    cpl_msg_info (cpl_func,"dit_acq = %g [us]", dit_acq);

    /* Create SYNC information */
    gravi_signal_create_sync (visacq_tmp, ntel, dit_acq,
                              vismet_table, ntel, "MET");

    /* Create column in output table */
	gravi_table_new_column (vismet_table, "OPD_PUPIL", "m", CPL_TYPE_DOUBLE);
    double * opd_met = cpl_table_get_data_double (vismet_table, "OPD_PUPIL");
    
    /* Get data from input table */
    double * opd_acq = cpl_table_get_data_double (visacq_tmp, "OPD_PUPIL");
    int * nspot = cpl_table_get_data_int (visacq_tmp, "PUPIL_NSPOT");
    int * first = cpl_table_get_data_int (visacq_tmp, "FIRST_MET");
    int * last  = cpl_table_get_data_int (visacq_tmp, "LAST_MET");
    CPLCHECK_MSG ("Cannot load data");

    /* Loop on beam */
    for (cpl_size tel = 0; tel < ntel; tel++) {

        /* Loop on ACQ rows with undetected spot */
        if (gravi_pfits_get_mjd (header) < 57876.5) {
            cpl_msg_info (cpl_func, "Compute OPD_PUPIL for blink ACQ frames");
            for (cpl_size row = 1; row < nrow_acq-1; row++) {
                if (nspot[row*ntel+tel] == 0 &&
                    nspot[(row-1)*ntel+tel] > 0 &&
                    nspot[(row+1)*ntel+tel] > 0) {
                    opd_acq[row*ntel+tel] = 0.5* (opd_acq[(row-1)*ntel+tel] + opd_acq[(row+1)*ntel+tel]);
                }
            }
        }
        
        /* Loop on ACQ rows, fill the corresponding MET rows */
        for (cpl_size row = 0; row < nrow_acq; row++) {            
            for (cpl_size row_met = first[row*ntel+tel]; row_met < last[row*ntel+tel]; row_met++) {
                opd_met[row_met*ntel+tel] = opd_acq[row*ntel+tel];
            }
        }

        /* Loop on MET rows, to fill the empty by the closet futur value */
        double opd = opd_met[nrow_acq*ntel+tel];
        for (cpl_size row = nrow_met-1; row >= 0; row--) {
            if (opd_met[row*ntel+tel] != 0 ) opd = opd_met[row*ntel+tel];
            else opd_met[row*ntel+tel] = opd;
        }

        /* Loop on MET rows, to fill the empty by the closet past value */
        opd = opd_met[0*ntel+tel];
        for (cpl_size row = 0; row < nrow_met; row++) {
            if (opd_met[row*ntel+tel] != 0) opd = opd_met[row*ntel+tel];
            else opd_met[row*ntel+tel] = opd;
        }
        
    }/* End loop on beam */

    /* Free the tmp table */
    FREE (cpl_table_delete, visacq_tmp);
    
    /* Return */
    gravi_msg_function_exit(1);
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Fill the VIS_MET table with the DRS algorithm
 *
 * @param metrology_table: input METROLOGY table
 * @param vismet_table: output OI_VIS_MET table
 * @param header:     corresponding HEADER
 *
 * Fill the VIS_MET table from the METROLOGY table with
 * the pipeline alorithm. The function creates the columns
 * PHASE_MET_FC (scalar) and PHASE_MET_TEL
 * (array of 4 values = diodes).
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_metrology_drs (cpl_table * metrology_table,
                                    cpl_table * vismet_table,
                                    cpl_propertylist * header)
{
    gravi_msg_function_start(1);
	cpl_ensure_code (metrology_table, CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (header,  CPL_ERROR_NULL_INPUT);
	
	int ntel = 4, ndiode = 4;
    char name[100];

    /* Parameters */
	int ind_sintel_FT[4][4]={{0,2,4,6},
	                         	{8,10,12,14},
	                         	{16,18,20,22},
	                         	{24,26,28,30}};
	int ind_costel_FT[4][4]={{1,3,5,7},
	                         	{9,11,13,15},
	                         	{17,19,21,23},
	                         	{25,27,29,31}};
	int ind_sintel_SC[4][4]={{32,34,36,38},
	                         	{40,42,44,46},
	                         	{48,50,52,54},
	                         	{56,58,60,62}};
	int ind_costel_SC[4][4]={{33,35,37,39},
	                         	{41,43,45,47},
	                         	{49,51,53,55},
	                         	{57,59,61,63}};
	int ind_sinfc_FT[4]={64,66,68,70};
	int ind_cosfc_FT[4]={65,67,69,71};
	int ind_sinfc_SC[4]={72,74,76,78};
	int ind_cosfc_SC[4]={73,75,77,79};
    
	cpl_array ** raw_met=cpl_table_get_data_array(metrology_table,"VOLT");
	CPLCHECK_MSG("get data met phase at tel");
    
	
	/* 
	 * Create the vismet_table
	 */

	cpl_msg_info (cpl_func,"Fill the OI_VIS_MET table with the DRS algorithm");
	
	cpl_size nbrow_met = cpl_table_get_nrow (metrology_table);
	 
	double phase_temp[4]={0,0,0,0}, phase_rtc;
	int k_wrap[4]={0,0,0,0};

	/* Create the output table for VIS_MET
	 * PHASE_TEL and PHASE_FC are defined as FT-SC */
	cpl_table_new_column (vismet_table, "PHASE_FC", CPL_TYPE_DOUBLE);
	cpl_table_set_column_unit (vismet_table, "PHASE_FC", "rad");
    cpl_table_fill_column_window (vismet_table, "PHASE_FC", 0, nbrow_met * ntel, 0.0);
	cpl_table_new_column_array (vismet_table, "PHASE_TEL", CPL_TYPE_DOUBLE, ndiode);
	cpl_table_set_column_unit (vismet_table, "PHASE_TEL", "rad");
	
	const char * date = gravi_pfits_get_met_ph(header);
	const char * acq_date = gravi_pfits_get_start_prcacq(header);
    CPLCHECK_MSG ("Cannot get dates");
    
	double met_mjd = 86400*1e6*(gravi_convert_to_mjd(date) -
								gravi_convert_to_mjd(acq_date));
    CPLCHECK_MSG ("Cannot convert dates");
    

    
	/*
	 * Compute the metrology phase at FC, without 
     * any smoothing
	 */
    double k_phase[4]={0,0,0,0};
    
	/* Loop on row of the metrology */
	int met_date_row = -1;
	for (cpl_size row = 0; row < nbrow_met; row ++) {
	  
	  /* Look for the reference row */
	  int time_met = cpl_table_get_int (metrology_table, "TIME", row, NULL);
	  if (row > 0) {
		if ((time_met > met_mjd) && (cpl_table_get_int (metrology_table, "TIME", row-1, NULL) < met_mjd)) {
		  met_date_row = row;
		}
	  }
	  
	  /* Loop on tel */
	  for (int tel = 0; tel < ntel; tel++) {
				
        /* Compute the phase from the VOLTs */
        double phi_ft = atan2 (cpl_array_get (raw_met[row],ind_sinfc_FT[tel], NULL),
                               cpl_array_get (raw_met[row],ind_cosfc_FT[tel], NULL));
        
        double phi_sc = atan2 (cpl_array_get (raw_met[row],ind_sinfc_SC[tel], NULL),
                               cpl_array_get (raw_met[row],ind_cosfc_SC[tel], NULL));
		
		/* Catch errors */
		CPLCHECK_MSG("Error at computation of phase");

		/* unwrap (phi_ft - phi_sc) */
		if ((phi_ft - phi_sc) > M_PI)  phi_ft-=2*M_PI;
		if ((phi_ft - phi_sc) < -M_PI) phi_ft+=2*M_PI;

		/* Keep memory of wraps */
		if      ((phi_ft - phi_sc)-phase_temp[tel] > M_PI)  k_wrap[tel]--;
 		else if ((phi_ft - phi_sc)-phase_temp[tel] < -M_PI) k_wrap[tel]++;

        /* Integrate phase */
		phase_temp[tel] = (phi_ft - phi_sc);
		double phase_d = (phase_temp[tel] + 2 *k_wrap[tel]* M_PI) ;

		/* When we reach the metrology reference date */
		if (met_date_row == row){
            sprintf (name, "ESO OCS MET PH_FC%d_FT", tel+1);
            phase_rtc = cpl_propertylist_get_double (header, name);
            sprintf (name, "ESO OCS MET PH_FC%d_SC", tel+1);
            phase_rtc = phase_rtc - cpl_propertylist_get_double (header, name);
            k_phase[tel] = phase_rtc - (phase_temp[tel] + 2 *k_wrap[tel]* M_PI);
            k_phase[tel] = floor(k_phase[tel]/(2*M_PI)+0.5)*2*M_PI;
		}
		
		/* Set the PHASE_SC and the TIME */
        cpl_table_set (vismet_table, "PHASE_FC", row*ntel+tel, phase_d);
        
        CPLCHECK_MSG("Computing the metrology phase at FC");
	  }
	  /* End loop on tel */
	}
	/* End loop on rows */

	/* IF error of timing */
	if (met_date_row == -1){
	  cpl_error_set_message (cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                             "The metrology phase date is not within"
                             " the boundaries of RMN acquisition");
	  return CPL_ERROR_ILLEGAL_INPUT;
	}
    
    /* Loop on rows to re-apply to k_phase */
    for (cpl_size row = 0; row < nbrow_met; row++) {
        for (int tel = 0; tel < 4 ; tel++) {
            double value = cpl_table_get (vismet_table, "PHASE_FC", row*ntel+tel, NULL);
            cpl_table_set (vismet_table, "PHASE_FC", row*ntel+tel, value + k_phase[tel]);
        }
    } /* End loop on rows */
	
	CPLCHECK_MSG("Computing the metrology phase at FC");

    

	/*
	 * Compute the metrology phase at telescope, with a temporal
     * smoothing over n_filter samples
	 */
    
    cpl_array * tel_phase = cpl_array_new (ndiode, CPL_TYPE_DOUBLE);
    
	int k_wrap_tel[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double phase_temp_tel[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double k_phase_tel[16]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	cpl_array * phase;
	cpl_array * phase_ft = cpl_array_new (16, CPL_TYPE_DOUBLE_COMPLEX);
	cpl_array * phase_sc_conj = cpl_array_new (16, CPL_TYPE_DOUBLE_COMPLEX);
	CPLCHECK_MSG("creat array met phase at tel");

	cpl_array_fill_window_complex(phase_ft, 0, 16, 0.+I*0.);
	cpl_array_fill_window_complex(phase_sc_conj, 0, 16, 0.+I*0.);
	CPLCHECK_MSG("Fill win met phase at tel");
	int n_filter=25;

	/* Compute the sum phase_ft and phase_sc_conj of the first n_filter*2+1 row
	 *  phase_ft=<cos_ft>+i<sin_ft>
	 *  phase_sc_conj=<cos_sc>-i<sin_sc> */
	CPLCHECK_MSG("Init met phase at tel");
	for (cpl_size row=0; row<n_filter*2+1; row++) {
		for (int tel = 0; tel < ntel; tel++){
			for (int diode = 0; diode < ndiode; diode++){
                int ind = tel*ndiode+diode;
				cpl_array_set_complex(phase_ft, ind,
									  cpl_array_get_complex(phase_ft,ind, NULL)
									  +cpl_array_get(raw_met[row],ind_costel_FT[tel][diode], NULL)
									  +I*cpl_array_get(raw_met[row],ind_sintel_FT[tel][diode], NULL));
				cpl_array_set_complex(phase_sc_conj,ind,
									  cpl_array_get_complex(phase_sc_conj,ind, NULL)
									  +cpl_array_get(raw_met[row],ind_costel_SC[tel][diode], NULL)
									  -I*cpl_array_get(raw_met[row],ind_sintel_SC[tel][diode], NULL));
			}
		}
	}
	CPLCHECK_MSG("Computing met phase at tel");

    /* arg{FT - SC} */
	phase = cpl_array_duplicate (phase_ft);
	cpl_array_multiply (phase, phase_sc_conj);
	cpl_array_arg (phase);

	/* case for n_filter first lines */
	for (cpl_size row=0; row<n_filter+1; row++) {
		for (int tel = 0; tel < ntel; tel++) {
            for (int diode = 0; diode < ndiode; diode++) {
                cpl_array_set (tel_phase, diode, cpl_array_get (phase,tel*ndiode+diode, NULL));
            }
            cpl_table_set_array (vismet_table, "PHASE_TEL", row*ntel+tel, tel_phase);
            CPLCHECK_MSG ("Cannot set");
        }
    }

	/* case for all lines */
	for (cpl_size row=n_filter+1; row<nbrow_met-n_filter; row++) {
        
		/* update the mean phase_ft=<cos_ft>+i<sin_ft> and phase_sc_conj=<cos_sc>-i<sin_sc>
		 * 	Remove element row-n_filter-1 and add row+n_filter*/
		for (int tel = 0; tel < ntel; tel++){
			for (int diode = 0; diode < ndiode; diode++){
                int ind = tel*ndiode+diode;
				cpl_array_set_complex(phase_ft, ind,
						cpl_array_get_complex(phase_ft,ind, NULL)
									  -(cpl_array_get(raw_met[row-n_filter-1],ind_costel_FT[tel][diode], NULL)
										+I*cpl_array_get(raw_met[row-n_filter-1],ind_sintel_FT[tel][diode], NULL))
									  +cpl_array_get(raw_met[row+n_filter],ind_costel_FT[tel][diode], NULL)
									  +I*cpl_array_get(raw_met[row+n_filter],ind_sintel_FT[tel][diode], NULL));
				cpl_array_set_complex(phase_sc_conj,ind,
						cpl_array_get_complex(phase_sc_conj,ind, NULL)
									  -(cpl_array_get(raw_met[row-n_filter-1],ind_costel_SC[tel][diode], NULL)
										-I*cpl_array_get(raw_met[row-n_filter-1],ind_sintel_SC[tel][diode], NULL))
									  +cpl_array_get(raw_met[row+n_filter],ind_costel_SC[tel][diode], NULL)
									  -I*cpl_array_get(raw_met[row+n_filter],ind_sintel_SC[tel][diode], NULL));
				CPLCHECK_MSG("loop met phase at tel");
			}
		}

		cpl_array_delete (phase);
        
		/* compute the phase */
		CPLCHECK_MSG("Computing met phase at tel");
		phase=cpl_array_duplicate(phase_ft);
		cpl_array_multiply(phase, phase_sc_conj);
		cpl_array_arg(phase);

		/* Keep memory of wraps */
		for (int ind=0; ind<16; ind++){
			if      (cpl_array_get(phase, ind, NULL)-phase_temp_tel[ind] > M_PI)  k_wrap_tel[ind]--;
			else if (cpl_array_get(phase, ind, NULL)-phase_temp_tel[ind] < -M_PI) k_wrap_tel[ind]++;

			phase_temp_tel[ind]=cpl_array_get(phase, ind, NULL);
			cpl_array_set(phase, ind,(phase_temp_tel[ind] + 2 *k_wrap_tel[ind]* M_PI)) ;
		}
		CPLCHECK_MSG("Unwrap met phase at tel");

		/* When we reach the metrology reference date */
		if (met_date_row == row){
			for (int tel = 0; tel < ntel; tel++){
				for (int diode = 0; diode < 4; diode++){
                    int ind = tel*ndiode+diode;
					sprintf (name, "ESO OCS MET PH_T%d_D%d_FT",tel+1, diode+1);
					phase_rtc = cpl_propertylist_get_double (header, name);
					sprintf (name, "ESO OCS MET PH_T%d_D%d_SC",tel+1, diode+1);
					phase_rtc = phase_rtc - cpl_propertylist_get_double (header, name);
					k_phase_tel[ind] = phase_rtc - (phase_temp_tel[ind] + 2 *k_wrap_tel[ind]* M_PI);
					k_phase_tel[ind] = floor(k_phase_tel[ind]/(2*M_PI)+0.5)*2*M_PI;
				}
			}
			CPLCHECK_MSG("Unwrap met phase at tel");
		}

        /* Descramble and set */
        for (int tel = 0; tel < ntel; tel++) {
            for (int diode = 0; diode < ndiode; diode++) {
                cpl_array_set (tel_phase, diode, cpl_array_get (phase,tel*ndiode+diode, NULL));
            }
            cpl_table_set_array (vismet_table, "PHASE_TEL", row*ntel+tel, tel_phase);
            CPLCHECK_MSG ("Cannot set");
        }
        
	} /* End loop on rows */

    cpl_array_delete (phase_sc_conj);
	cpl_array_delete (phase);
	cpl_array_delete (phase_ft);
    cpl_array_delete (tel_phase);

	/* case for n_filter last lines */
    for (int tel = 0; tel < ntel; tel++) {
        const cpl_array * last_phase = cpl_table_get_array (vismet_table,"PHASE_TEL", (nbrow_met-n_filter-1)*ntel+tel);
        for (cpl_size row = nbrow_met-n_filter; row<nbrow_met; row++) {
            cpl_table_set_array (vismet_table,"PHASE_TEL", row*ntel+tel, last_phase);
            CPLCHECK_MSG ("Cannot set");
        }
    }

	/* IF error of timing */
	if (met_date_row == -1){
        cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
                              "The metrology phase date is not within"
                              " the boundaries of RMN acquisition");
	  return CPL_ERROR_ILLEGAL_INPUT;
	}
    
    /* Loop on rows to re-apply to k_phase */
    cpl_array ** all_phase = cpl_table_get_data_array (vismet_table, "PHASE_TEL");
    
    for (cpl_size row = 0; row < nbrow_met; row++) {
		for (int tel = 0; tel < ntel; tel++) {
			for (int diode = 0; diode < 4; diode++) {
                int ind = tel*ndiode+diode;
                double value = cpl_array_get (all_phase[row*ntel+tel], diode, NULL);
				cpl_array_set (all_phase[row*ntel+tel], diode, value + k_phase_tel[ind]);
            }
        }
    }
    /* End loop on rows */

    
	/*
	 * Compute the TEL-FC complex phasor = Vtel * conj (Vfc)
     * without any smoothing
	 */
    cpl_msg_info (cpl_func,"Compute PHASOR_TELFC");
    
	cpl_table_new_column_array (vismet_table, "PHASOR_TELFC", CPL_TYPE_DOUBLE_COMPLEX, ndiode);
	cpl_table_set_column_unit (vismet_table, "PHASOR_TELFC", "V^4");
    cpl_array ** phasor_telfc = cpl_table_get_data_array (vismet_table, "PHASOR_TELFC");

    for (cpl_size row = 0; row < nbrow_met; row++) {
		for (int tel = 0; tel < ntel; tel++) {
            phasor_telfc[row*ntel+tel] = cpl_array_new (ndiode, CPL_TYPE_DOUBLE_COMPLEX);
            
            double complex V_fc;
            V_fc = (cpl_array_get (raw_met[row],ind_cosfc_FT[tel], NULL) +
                    cpl_array_get (raw_met[row],ind_sinfc_FT[tel], NULL) * I) *
                   (cpl_array_get (raw_met[row],ind_cosfc_SC[tel], NULL) - 
                    cpl_array_get (raw_met[row],ind_sinfc_SC[tel], NULL) * I);
            CPLCHECK_MSG ("Cannot compute V_fc");
            
			for (int diode = 0; diode < 4; diode++) {
                double complex V_tel;
                V_tel = (cpl_array_get (raw_met[row],ind_costel_FT[tel][diode], NULL) +
                         cpl_array_get (raw_met[row],ind_sintel_FT[tel][diode], NULL) * I) *
                        (cpl_array_get (raw_met[row],ind_costel_SC[tel][diode], NULL) - 
                         cpl_array_get (raw_met[row],ind_sintel_SC[tel][diode], NULL) * I);
                cpl_array_set_complex (phasor_telfc[row*ntel+tel], diode,
                                       V_tel * conj (V_fc));
                CPLCHECK_MSG ("Cannot set V_tel * conj (V_fc)");
            } /* End loop on diode */
            
        } /* End loop on tel */
    } /* End loop on rows */ 

    gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}


/*----------------------------------------------------------------------------*/
/**
 * @brief Compute the metrology signal from TAC algorithm
 *  
 * @param metrology_table  The input METROLOGY table (output)
 * @param header      The corresponding HEADER
 * @param vismet_table     An already existing table (output)
 * 
 * Wrapper to interface with the TAC real-time algorithm
 * This fill new columns in the existing VIS_MET table:
 * OPD_MET_FC and OPD_MET_TEL, as well as FLAGs...
 */
/*----------------------------------------------------------------------------*/

cpl_error_code gravi_metrology_tac (cpl_table * metrology_table,
									cpl_table * vismet_table,
									cpl_propertylist * header)
{
    gravi_msg_function_start(1);
	cpl_ensure_code (metrology_table, CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (header,     CPL_ERROR_NULL_INPUT);
	cpl_ensure_code (vismet_table,    CPL_ERROR_NULL_INPUT);

    char card_ft[100];
    char card_sc[100];
	int ndiode = 4, ntel = 4;
	cpl_size nrow_met = cpl_table_get_nrow (metrology_table);
	
	/* 
	 * Copy data to ensure they are of type DOUBLE
	 *    volts[nrow][ndiode] 
	 */
    
	double** volts = cpl_malloc (nrow_met * sizeof(double*));
	cpl_array ** volts_array = cpl_table_get_data_array (metrology_table,"VOLT");
	for (cpl_size row = 0; row < nrow_met; row ++) {
	  volts[row] = cpl_malloc (80 * sizeof(double));
	  for (cpl_size diode = 0; diode < 80; diode ++) {
		volts[row][diode] = cpl_array_get (volts_array[row], diode, NULL);
	  }
	}

	CPLCHECK_MSG ("Cannot load metrology data");
	
	/* 
     * Allocate memory for FC the output
     *     opd_fc[nrow*ntel] 
     */
    
	gravi_table_new_column (vismet_table, "FLAG_FC", NULL, CPL_TYPE_INT);
    int * flag_fc = cpl_table_get_data_int (vismet_table, "FLAG_FC");
    
    gravi_table_new_column (vismet_table, "OPD_FC", "m", CPL_TYPE_DOUBLE);
	double * opd_fc = cpl_table_get_data_double (vismet_table, "OPD_FC");

    gravi_table_new_column (vismet_table, "VAMP_FC_FT", "V", CPL_TYPE_DOUBLE);
	double * coher_fc_ft = cpl_table_get_data_double (vismet_table, "VAMP_FC_FT");

    gravi_table_new_column (vismet_table, "VAMP_FC_SC", "V", CPL_TYPE_DOUBLE);
	double * coher_fc_sc = cpl_table_get_data_double (vismet_table, "VAMP_FC_SC");

	/* 
     * Allocate memory for TEL output
     *     opd_tel[nrow*ntel][ndiode]
     */
	
	gravi_table_new_column_array (vismet_table, "FLAG_TEL", NULL, CPL_TYPE_INT, ndiode);
	cpl_array ** flag_tel_array = cpl_table_get_data_array (vismet_table,"FLAG_TEL");
    
	gravi_table_new_column_array (vismet_table, "OPD_TEL", "m", CPL_TYPE_DOUBLE, ndiode);
	cpl_array ** opd_tel_array = cpl_table_get_data_array (vismet_table,"OPD_TEL");
    
	gravi_table_new_column_array (vismet_table, "VAMP_TEL_FT", "V", CPL_TYPE_DOUBLE, ndiode);
	cpl_array ** coher_tel_ft_array = cpl_table_get_data_array (vismet_table,"VAMP_TEL_FT");
    
	gravi_table_new_column_array (vismet_table, "VAMP_TEL_SC", "V", CPL_TYPE_DOUBLE, ndiode);
	cpl_array ** coher_tel_sc_array = cpl_table_get_data_array (vismet_table,"VAMP_TEL_SC");

    
	/* Wrap the newly computed phi_tel into the cpl_table. Wrapping make
     * the data 'valid' without having to copy them. Thus efficient */
	double ** opd_tel      = cpl_malloc (sizeof(double*) * nrow_met * ntel);
	int ** flag_tel        = cpl_malloc (sizeof(int*) * nrow_met * ntel);
	double ** coher_tel_ft = cpl_malloc (sizeof(double*) * nrow_met * ntel);
	double ** coher_tel_sc = cpl_malloc (sizeof(double*) * nrow_met * ntel);
    
	for (cpl_size row = 0; row < nrow_met * ntel; row ++) {
        flag_tel[row]       = cpl_malloc (sizeof(int) * ndiode);
        flag_tel_array[row] = cpl_array_wrap_int (flag_tel[row], ndiode);
        
        opd_tel[row]        = cpl_malloc (sizeof(double) * ndiode);
        opd_tel_array[row]  = cpl_array_wrap_double (opd_tel[row], ndiode);
        
        coher_tel_ft[row]       = cpl_malloc (sizeof(double) * ndiode);
        coher_tel_ft_array[row] = cpl_array_wrap_double (coher_tel_ft[row], ndiode);
        
        coher_tel_sc[row]       = cpl_malloc (sizeof(double) * ndiode);
        coher_tel_sc_array[row] = cpl_array_wrap_double (coher_tel_sc[row], ndiode);
	}

	CPLCHECK_MSG ("Cannot allocate output memory");

    

	/* Init the TAC algorithm -- this allocate memory */
    structTacData * tacData;
    structTacConfiguration * tacConfiguration;

	/* get the laser wavelength data */
	double lambda_met_mean =  gravi_pfits_get_met_wavelength_mean(header, metrology_table);
	cpl_msg_info (cpl_func,"Lambda met mean :%f nm",  lambda_met_mean*1e9);

    tacData = metrology_makeDefaultTacData(lambda_met_mean);
    tacConfiguration = tacData->tacConfiguration;

    /* Loop on time sample to run the TAC algorithm */
    for(long sample_number = 1; sample_number <= nrow_met; sample_number++) {
        
	  /* Run the TAC algorithm */
	  tacData->sample_number = sample_number;
	  tacData->buffer_idx_avg = ((sample_number - 1) % tacConfiguration->number_to_average);
	  metrology_read_voltages(tacData, volts[sample_number-1]);
	  metrology_algorithm(tacData);
        
	  /* Fill the output arrays as D1T1, D2T1, D3T1, D4T1, D1T2, D2T2... */
	  int idx_ft = 0;
	  int idx_sc = 32;
	  for (int idx = 0, tel = 0; tel < ntel; tel++) {
        cpl_size nmet = (sample_number-1)*ntel + tel;
        
		for (int diode = 0; diode < ndiode; diode++) {
		  /* TAC computation */
		  opd_tel[nmet][diode] = tacData->opl_telescope_diode[tel][diode];
		  flag_tel[nmet][diode] = tacData->total_flag_telescope[tel][diode][FT] | tacData->total_flag_telescope[tel][diode][SC];
          
		  /* Volt amplitude */
		  coher_tel_ft[nmet][diode] = sqrt (volts[sample_number-1][idx_ft+idx*2] * volts[sample_number-1][idx_ft+idx*2] + 
											volts[sample_number-1][idx_ft+idx*2+1] * volts[sample_number-1][idx_ft+idx*2+1]);
		  coher_tel_sc[nmet][diode] = sqrt (volts[sample_number-1][idx_sc+idx*2] * volts[sample_number-1][idx_sc+idx*2] + 
											volts[sample_number-1][idx_sc+idx*2+1] * volts[sample_number-1][idx_sc+idx*2+1]);
		  idx++;
		}
	  }
	  
	  /* Fill the output arrays as T1, T2, T3, T4... */
	  idx_ft = 64;
	  idx_sc = 72;
	  for (int idx = 0, tel = 0; tel < ntel; tel++) {
        cpl_size nmet = (sample_number-1)*ntel + tel;
        
	    /* TAC computation */
	    opd_fc[nmet]  = tacData->opl_fiber_coupler[tel];
		flag_fc[nmet] = tacData->total_flag_fiber_coupler[tel][FT] | tacData->total_flag_fiber_coupler[tel][SC];
	    /* Volt amplitude */
	    coher_fc_ft[nmet] = sqrt (volts[sample_number-1][idx_ft+idx*2] * volts[sample_number-1][idx_ft+idx*2] + 
                                  volts[sample_number-1][idx_ft+idx*2+1] * volts[sample_number-1][idx_ft+idx*2+1]);
	    coher_fc_sc[nmet] = sqrt (volts[sample_number-1][idx_sc+idx*2] * volts[sample_number-1][idx_sc+idx*2] + 
                                  volts[sample_number-1][idx_sc+idx*2+1] * volts[sample_number-1][idx_sc+idx*2+1]);
		idx++;
	  }
    }
    /* End loop on time sample */

	
	/* Get the TIME of the header reference phase */
	const char * date = gravi_pfits_get_met_ph (header);
	const char * acq_date = gravi_pfits_get_start_prcacq (header);
	double time_ref = 86400*1e6*(gravi_convert_to_mjd (date) - gravi_convert_to_mjd (acq_date));
	int * time_met = cpl_table_get_data_int (metrology_table, "TIME");

	cpl_size row_ref = -1;
	for (cpl_size row = 1; row < nrow_met; row++) 
	  if (time_met[row]>time_ref && time_met[row-1]<time_ref) {
	    row_ref = row; break;
	  }
	
	/* Check if we have found the metrology reference date inside the file */
	if (row_ref == -1) {
	  cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
				"The metrology phase date is not within"
				" the boundaries of RMN acquisition");
	  return CPL_ERROR_NULL_INPUT;
	} else {
	  cpl_msg_info (cpl_func,"The metrology phase date is OK (row %lld over %lld).", row_ref, nrow_met);
	}

	/* Re-install the absolute phase for the FC diodes */
	for (int tel = 0; tel < ntel; tel++) {
        sprintf (card_ft, "ESO OCS MET PH_FC%d_FT", tel+1);
        sprintf (card_sc, "ESO OCS MET PH_FC%d_SC", tel+1);
        double opd_ref =  lambda_met_mean / CPL_MATH_2PI *
            (cpl_propertylist_get_double (header, card_sc) - cpl_propertylist_get_double (header, card_ft)) - 
            opd_fc[row_ref*ntel+tel];
        double opd_ref_int = gravi_round (opd_ref / lambda_met_mean) * lambda_met_mean;
        for (cpl_size row = 0; row < nrow_met; row++) {
            opd_fc[row*ntel+tel] += opd_ref_int;
        }
	}
	/* End loop on tel and diodes */

	/* Re-install the absolute phase for the TEL diodes */
	for (int tel = 0; tel < ntel; tel++) {
        for (int diode = 0; diode < ndiode; diode++) {
            sprintf (card_ft,"ESO OCS MET PH_T%d_D%d_FT", tel+1, diode+1);
            sprintf (card_sc,"ESO OCS MET PH_T%d_D%d_SC", tel+1, diode+1);
            double opd_ref = lambda_met_mean / CPL_MATH_2PI *
                (cpl_propertylist_get_double (header, card_sc) - cpl_propertylist_get_double (header, card_ft)) - 
                opd_tel[row_ref*ntel+tel][diode];
            double opd_ref_int = gravi_round (opd_ref / lambda_met_mean) * lambda_met_mean;
            for (cpl_size row = 0; row < nrow_met; row++) {
                opd_tel[row*ntel+tel][diode] += opd_ref_int;
            }
        }
	}
    
	CPLCHECK_MSG ("Cannot fill result of metrology computation");

	/* Free the pointer to pointer to data */
	FREELOOP (cpl_free, volts, nrow_met);
	FREE (cpl_free, opd_tel);
	FREE (cpl_free, flag_tel);
	FREE (cpl_free, coher_tel_ft);
	FREE (cpl_free, coher_tel_sc);

    /* Free the TAC data */
    FREE (cpl_free, tacConfiguration);
    FREE (cpl_free, tacData);

    gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}

/* -------------------------------------------------------------------- */
/**
 * @brief Create the P2VM of the metrology
 *
 * @param metrology_table  A METROLOGY table with modulation on all beams
 * 
 * The METROLOGY volts are fitted with ellipse in order to compute the
 * TRANSMISSION, PHASE and COHERENCE of each diode. The resulting
 * P2VM_MET table is returned.
 */
/* -------------------------------------------------------------------- */
static int sin_lambda(const double x[], const double a[], double *result){

	*result = a[0] + a[1] * sin(x[0] * (2 * M_PI) / a[3]) +
			  a[2] * cos(x[0] * (2 * M_PI) / a[3]);

	return (0);
}

static int dfda_sin(const double x[], const double a[], double result[]){
	result[0] = 1;
	result[1] = sin(x[0] * (2 * M_PI) / a[3]);
	result[2] = cos(x[0] * (2 * M_PI) / a[3]);
	result[3] = a[1] * x[0] * (2 * M_PI) / (a[3] * a[3]) *
				sin(x[0] * (2 * M_PI) / a[3]) - a[1] * x[0] *
				(2 * M_PI) / (a[3] * a[3]) * cos(x[0] * (2 * M_PI) / a[3]);
	return (0);
}

cpl_table * gravi_metrology_compute_p2vm (cpl_table * metrology_table, double wave_met)
{
	gravi_msg_function_start(1);
	cpl_ensure (metrology_table, CPL_ERROR_NULL_INPUT, NULL);
	
	cpl_vector * vectA, * vectB,
			   * y_sigma, * init_val, *opl_vect;
	const cpl_array * volt;
	cpl_array * coherence_fit, * phase_fit;
	cpl_matrix * opd_matrix;
	int infos = 0;
	cpl_vector ** vect;
	int val_to_fit2[] = {1,1,1,0}, nv;
	double mse, red_chisq;
	cpl_array * trans;

	int ntel = 4 ;
	int n_diode = 80;
	
	/* Get data */
	cpl_size nrow = cpl_table_get_nrow (metrology_table);
	CPLCHECK_NUL ("Cannot get data");

	/* Create p2vm_met : (REGNAME, TRANSMISSION, COHERENCE and PHASE) */
	cpl_table * p2vm_met = cpl_table_new (ntel * 2);
	cpl_table_new_column_array (p2vm_met,"TRANSMISSION", CPL_TYPE_DOUBLE, 2);
	cpl_table_new_column_array (p2vm_met,"COHERENCE", CPL_TYPE_DOUBLE, 2);
	cpl_table_new_column_array (p2vm_met,"PHASE", CPL_TYPE_DOUBLE, 2);

	CPLCHECK_NUL ("Allocate the metrology output");

	/* Loop on SC and FT */
	for (int gravi_type = 0; gravi_type < 2; gravi_type++ ){
		int comb = (gravi_type == GRAVI_SC ? 1 : 2);
		
		/* Loop on diodes */
		for (int tel = 0; tel < ntel; tel++){

			/* load vectA and vectB from metrology */
			vectA = cpl_vector_new (nrow);
			vectB = cpl_vector_new (nrow);

			for (cpl_size row = 0; row < nrow; row ++){
				volt = cpl_table_get_array ( metrology_table, "VOLT", row);
				cpl_vector_set (vectA, row, cpl_array_get (volt , n_diode - 2*(comb*ntel - tel), &nv));
				cpl_vector_set (vectB, row, cpl_array_get (volt , n_diode - 2*(comb*ntel - tel) + 1, &nv));
			}

			/* Compute phase from vectA and vectB*/
			opl_vect = gravi_ellipse_phase_create (vectA, vectB, NULL);
			cpl_vector_multiply_scalar (opl_vect, wave_met/(2.0*M_PI));
			CPLCHECK_NUL("Compute OPD");

			/* put OPD into opd_matrix for P2VM_MET calib */
			opd_matrix = cpl_matrix_new(nrow, 1);
			for (cpl_size row = 0; row < nrow; row++){
				cpl_matrix_set (opd_matrix, row, 0, cpl_vector_get (opl_vect, row));
			}
			
			cpl_vector_delete (opl_vect);

			/* fit on a central window of 5*lambda width */

			vect = cpl_malloc (2*sizeof(cpl_vector*));
			vect[0] = vectA;
			vect[1] = vectB;
			phase_fit = cpl_array_new (2, CPL_TYPE_DOUBLE);
			coherence_fit = cpl_array_new (2, CPL_TYPE_DOUBLE);
			trans = cpl_array_new (2, CPL_TYPE_DOUBLE);
			
			for (int i = 0; i < 2; i++){

				/* Get the spectrum of the vector region */
				y_sigma = cpl_vector_new(nrow);
				cpl_vector_fill(y_sigma, 1);

				/* Define and initialize all variables to make a FIT */
				mse = 0;
				red_chisq = 0;
				init_val = cpl_vector_new(4);
				cpl_vector_set(init_val, 0, 1);
				cpl_vector_set(init_val, 1, 1);
				cpl_vector_set(init_val, 2, 1);
				cpl_vector_set(init_val, 3, wave_met);

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

			FREE (cpl_free, vect);
			FREE (cpl_matrix_delete, opd_matrix);
			
			cpl_table_set_array (p2vm_met, "COHERENCE", ntel*(1-gravi_type) + tel, coherence_fit);
			cpl_array_subtract_scalar (phase_fit,
					cpl_array_get_double (phase_fit, 0, &nv));
			cpl_table_set_array (p2vm_met, "PHASE", ntel*(1-gravi_type) + tel, phase_fit);
			cpl_table_set_array (p2vm_met, "TRANSMISSION", ntel*(1-gravi_type) + tel, trans);

			FREE (cpl_array_delete, trans);
			FREE (cpl_array_delete, coherence_fit);
			FREE (cpl_array_delete, phase_fit);
			CPLCHECK_NUL("End loop on tel");
		}
	}
	
	/* Verbose */
	gravi_msg_function_exit(1);
	return p2vm_met;
}

/* -------------------------------------------------------------------- */
/**
 * @brief Reduce the metrology
 *
 * @param data  a gravi_data with a METROLOGY extension (input/output)
 * 
 * The resulting VIS_MET table is created from the METROLOGY table
 * and added this data. It is then updated with the TAC algorithm.
 */
/* -------------------------------------------------------------------- */

cpl_error_code gravi_metrology_reduce (gravi_data * data,
                                       gravi_data * eop_data,
                                       const cpl_parameterlist * parlist)
{
    gravi_msg_function_start(1);
	cpl_ensure_code (data, CPL_ERROR_NULL_INPUT);

	/* Load data */
	cpl_table * metrology_table = gravi_data_get_table (data, GRAVI_METROLOGY_EXT);
	cpl_propertylist * header = gravi_data_get_header (data);
	CPLCHECK_MSG ("Cannot load met extension");

    /* Create the table */
	cpl_table * vismet_table = NULL;
	vismet_table = gravi_metrology_create (metrology_table, header);
	CPLCHECK_MSG ("Cannot create vismet_table");

    /* Compuet the pointing */
    gravi_eop_pointing (vismet_table, header, eop_data);

    /* If VIS_ACQ table exist, we compute the OPD_PUPIL */
    if (gravi_data_has_extension (data, GRAVI_OI_VIS_ACQ_EXT)) {
        cpl_table * visacq_table;        
        visacq_table = gravi_data_get_table (data, GRAVI_OI_VIS_ACQ_EXT);
        double delay = gravi_param_get_double_default (parlist, "gravity.metrology.acq-correction-delay",0.0);
        gravi_metrology_acq (visacq_table, vismet_table, delay, header);
    }

	/* Reduce the metrology with the DRS algorithm, this 
	 * creates the OI_VIS_MET table */
	gravi_metrology_drs (metrology_table, vismet_table, header);
	CPLCHECK_MSG ("Cannot reduce metrology with DRS algo");

	/* Add the columns from TAC algorithm */
	gravi_metrology_tac (metrology_table, vismet_table, header);
	CPLCHECK_MSG ("Cannot reduce metrology with TAC algo");

	/* Add the VISMET_TABLE table to the gravi_data */
	gravi_data_add_table (data, NULL, GRAVI_OI_VIS_MET_EXT, vismet_table);
	CPLCHECK_MSG ("Cannot add OI_VIS_MET in p2vmred_data");
	
    gravi_msg_function_exit(1);
	return CPL_ERROR_NONE;
}

/**@}*/
