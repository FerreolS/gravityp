
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_align.h 177213 2008-12-05 15:33:51Z hlorch $"
 *
 * Functions for object alignment
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2006-10-04  created
 */

#ifndef CLIPM_ALIGN_H
#define CLIPM_ALIGN_H

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
 * @defgroup    clipm_align_options Alignment Options
 * @ingroup     clipm_align
 * @brief       Available transformation options. These options can be
 *              combined bit-wise to allow respective
 *              transformation options.
 * 
 * @par Supported Alignment Modes:
 * 
 * The following alignment modes are supported by CLIPM aligning functions.
 *  
 * - An alignment mode of type @a clipm_align_opt should be given by a
 *   bit-wise combination of the options below.
 * - CLIPM_ALIGN_FREE overrides scaling, rotation and shift, and allows for
 *   any linear transformation, including deformations and projections.
 * - CLIPM_ALIGN_ROBUST is so far only supported together with (and only with)
 *   CLIPM_ALIGN_SHIFT.
 * - CLIPM_ALIGN_ROTATE is only supported for two-dimensional input (x and y).
 * - The following combinations of modes are supported (shown column-wise):
 *   \n\n
 *   \latexonly
  \begin{tabular}{l|l|l|l|l|l|l|l|l|l}
  Alignment Mode & \multicolumn{9}{l}{Combinations} \\
  \hline
  CLIPM\_ALIGN\_SHIFT  & 1 & 0 & 1 & 0 & 1 & 0 & 1 & X & 1 \\
  CLIPM\_ALIGN\_SCALE  & 0 & 1 & 1 & 0 & 0 & 1 & 1 & X & 0 \\
  CLIPM\_ALIGN\_ROTATE & 0 & 0 & 0 & 1 & 1 & 1 & 1 & X & 0 \\
  CLIPM\_ALIGN\_FREE   & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
  CLIPM\_ALIGN\_ROBUST & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\
  \hline
  Notes:             &   &   &   & \multicolumn{4}{l}{only 2-dim} & & Median
  \end{tabular}
 *   \endlatexonly
 *   \htmlonly
  <table align=center cellpadding=2 cellspacing=1 class=memitem>
  <tr><th><em>Alignment Mode</em></th><th colspan=9><em>Combinations</em></th>
      </tr>
  <tr><td> CLIPM_ALIGN_SHIFT  </td><td> 1 </td><td> 0 </td><td> 1 </td>
      <td> 0 </td><td> 1 </td><td> 0 </td><td> 1 </td><td> X </td><td> 1 </td>
      </tr>
  <tr><td> CLIPM_ALIGN_SCALE  </td><td> 0 </td><td> 1 </td><td> 1 </td>
      <td> 0 </td><td> 0 </td><td> 1 </td><td> 1 </td><td> X </td><td> 0 </td>
      </tr>
  <tr><td> CLIPM_ALIGN_ROTATE </td><td> 0 </td><td> 0 </td><td> 0 </td>
      <td> 1 </td><td> 1 </td><td> 1 </td><td> 1 </td><td> X </td><td> 0 </td>
      </tr>
  <tr><td> CLIPM_ALIGN_FREE   </td><td> 0 </td><td> 0 </td><td> 0 </td>
      <td> 0 </td><td> 0 </td><td> 0 </td><td> 0 </td><td> 1 </td><td> 0 </td>
      </tr>
  <tr><td> CLIPM_ALIGN_ROBUST </td><td> 0 </td><td> 0 </td><td> 0 </td>
      <td> 0 </td><td> 0 </td><td> 0 </td><td> 0 </td><td> 0 </td><td> 1 </td>
      </tr>
  <tr><td><em> Notes: </em></td> <td></td> <td></td> <td></td>
      <td colspan=4 bgcolor="#ffffff">only 2-dim</td> <td></td>
      <td bgcolor="#ffffff">Median</td> </tr>
  </table>
 *   \endhtmlonly
 * 
 * @par Compatibility:
 * 
 * - MIDAS users will find the same functionality as the function ALIGN/IMAGE
 *   with the following mode combinations:
  \htmlonly
  <table align=center cellpadding=2 cellspacing=1 class=memitem>
  <tr><th><em>MIDAS align mode</em></th>
      <th colspan=9><em>CLIPM align mode</em></th></tr>
  <tr><td> SHIFT </td><td> CLIPM_ALIGN_SHIFT  </td></tr>
  <tr><td> UNIT  </td><td> CLIPM_ALIGN_SHIFT | CLIPM_ALIGN_ROTATE  </td></tr>
  <tr><td> EQUAL </td><td> CLIPM_ALIGN_SHIFT | CLIPM_ALIGN_SCALE |
                           CLIPM_ALIGN_ROTATE </td></tr>
  <tr><td> FREE  </td><td> CLIPM_ALIGN_FREE   </td></tr>
  </table>
  \endhtmlonly
 * 
 * @note
 * 
 * - Enabling/disabling rotation and/or scaling influences the @a shift result.
 */
/** @{ */

/* Doxygen requires more than the brief description, so remember to provide it
 * everywhere below. The brief description is interpreted as the first line
 * up to the first point(.). */

/** @brief Alignment option type */
typedef unsigned int clipm_align_opt;

/**
 * @brief  Shifting allowed.
 * 
 * Shifting is allowed as part of the transformation.
 */
extern const clipm_align_opt CLIPM_ALIGN_SHIFT;
/**
 * @brief  Scaling allowed.
 * 
 * Scaling is allowed as part of the transformation.
 */
extern const clipm_align_opt CLIPM_ALIGN_SCALE;
/**
 * @brief  Rotation allowed.
 * 
 * For the two-dimensional case, the rotation will be determined, i.e. the
 * transformation matrix will be a (evt. scaled) rotation matrix. Other
 * numbers of dimensions are not implemented yet.
 */
extern const clipm_align_opt CLIPM_ALIGN_ROTATE;
/**
 * @brief  All parameters are free.
 * 
 * Distortion/shearing is also allowed. CLIPM_ALIGN_SHIFT, CLIPM_ALIGN_SCALE,
 * and CLIPM_ALIGN_ROTATE are overridden.
 */
extern const clipm_align_opt CLIPM_ALIGN_FREE;
/**
 * @brief  Robust computation, ignore outliers.
 * 
 * Outliers are determined and ignored. This is so far only implemented for the
 * mode CLIPM_ALIGN_SHIFT, where it means the median translation.
 */
extern const clipm_align_opt CLIPM_ALIGN_ROBUST;

/** @} */

/*-----------------------------------------------------------------------------
    Prototypes
 -----------------------------------------------------------------------------*/

void            clipm_align_opt_sprint_literal(
                                            char                *str,
                                            clipm_align_opt     opt);

cpl_error_code  clipm_align_points(         const cpl_matrix    *ref_points,
                                            const cpl_matrix    *in_points,
                                            const cpl_matrix    *ref_variances,
                                            const cpl_matrix    *in_variances,
                                            clipm_align_opt align_mode_bitmask,
                                            cpl_matrix  **transform_matrix,
                                            cpl_matrix  **shift,
                                            cpl_matrix  **residuals);

cpl_error_code  clipm_align_correlate(      const cpl_image     *ref_img,
                                            const cpl_image     *test_img,
                                            const cpl_matrix    *ref_locations,
                                            const cpl_matrix    *test_locations,
                                            unsigned int    area_size,
                                            double          max_distance,
                                            clipm_align_opt align_mode_bitmask,
                                            cpl_matrix **pixel_transform_matrix,
                                            cpl_matrix **pixel_transshiftvector,
                                            cpl_matrix      **all_pixelshifts,
                                            cpl_matrix     **all_locationshifts,
                                            cpl_matrix      **all_uncertainties,
                                            cpl_array       **all_error_codes);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif  /* CLIPM_ALIGN_H */
