
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_image_display.h 153800 2007-07-17 16:45:50Z hlorch $"
 *
 * Tool functions for the display of images
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2007-06-28  created
 */

#ifndef CLIPM_IMAGE_DISPLAY_H
#define CLIPM_IMAGE_DISPLAY_H

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
    Declaration Block
 -----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------
    Types
 -----------------------------------------------------------------------------*/

/**
 * @defgroup    clipm_image_display_fill_options  Fill Pattern Modes
 * @ingroup     clipm_image_display
 * @brief       Fill pattern modes.
 * 
 * These modes define the method of area filling.
 */
/** @{ */

/** @brief Fill mode type */
typedef unsigned int clipm_fill_mode;

/** @brief  Fill with given value. */
extern const clipm_fill_mode CLIPM_FILL_VALUE;
/** @brief  Fill with foreground. */
extern const clipm_fill_mode CLIPM_FILL_FOREGROUND;
/** @brief  Fill with background. */
extern const clipm_fill_mode CLIPM_FILL_BACKGROUND;
/** @brief  Fill with raster of foreground and background. */
extern const clipm_fill_mode CLIPM_FILL_BWRASTER;
/** @brief  Fill with horizontal stripes of foreground and background. */
extern const clipm_fill_mode CLIPM_FILL_HSTRIPES;
/** @brief  Fill with vertical stripes of foreground and background. */
extern const clipm_fill_mode CLIPM_FILL_VSTRIPES;
/** @brief  Fill with caros. */
extern const clipm_fill_mode CLIPM_FILL_CAROS;
/** @brief  Fill with diagonal caros. */
extern const clipm_fill_mode CLIPM_FILL_DCAROS;
/** @brief  Fill with 45` diagonal stripes of foreground over background. */
extern const clipm_fill_mode CLIPM_FILL_DSTRIPES;
/** @brief  Fill with 135` diagonal stripes of foreground over background. */
extern const clipm_fill_mode CLIPM_FILL_DSTRIPES2;

/** @} */

/*-----------------------------------------------------------------------------
    Prototypes
 -----------------------------------------------------------------------------*/

cpl_image       *clipm_image_display_insert_gaps(
                                            const cpl_image *input,
                                            int             x_input_position,
                                            int             y_input_position,
                                            int             x_segmentsize,
                                            int             y_segmentsize,
                                            int             x_gapwidth,
                                            int             y_gapwidth,
                                            int             flag_gaps_bad,
                                            clipm_fill_mode fill_pattern,
                                            double          fill_brightness);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif  /* CLIPM_IMAGE_DISPLAY_H */
