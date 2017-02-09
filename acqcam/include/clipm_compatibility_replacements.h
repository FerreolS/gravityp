
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_compatibility_replacements.h 169734 2008-07-02 16:20:41Z hlorch $"
 *
 * Define the compatible CPL version
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2007-08-16  created
 */

#ifndef CLIPM_COMPATIBILITY_REPLACEMENTS_H
#define CLIPM_COMPATIBILITY_REPLACEMENTS_H

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "clipm_compatibility_cpl.h"

/*-----------------------------------------------------------------------------
    Function replacements for CPL <= 3
 -----------------------------------------------------------------------------*/
/**
 * @defgroup    clipm_compatibility_replacements_cpl_4  Compatibility Macros for CPL older than 4.0
 * @ingroup     clipm_compatibility_cpl
 * @brief       Macros for automatic replacement of CPL 4.x functions
 *              in compilation with older CPL
 * 
 * @par Synopsis:
 * @code
 *   #include "clipm_compatibility_replacements.h"
 * @endcode
 */
/** @{ */
#if CLIPM_GET_INSTALLED_CPL_VERSION <= 3

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Constant accessor functions
     *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /** @brief  Replaces cpl_array_get_data_int_const() by CPL 3.x function call. */
    #define cpl_array_get_data_int_const(array) \
        cpl_array_get_data_int(array)
    /** @brief  Replaces cpl_array_get_data_float_const() by CPL 3.x function call. */
    #define cpl_array_get_data_float_const(array) \
        cpl_array_get_data_float(array)
    /** @brief  Replaces cpl_array_get_data_double_const() by CPL 3.x function call. */
    #define cpl_array_get_data_double_const(array) \
        cpl_array_get_data_double(array)
    /** @brief  Replaces cpl_array_get_data_string_const() by CPL 3.x function call. */
    #define cpl_array_get_data_string_const(array) \
        (const char**)cpl_array_get_data_string(array) /* whyever necessary */
    /** @brief  Replaces cpl_bivector_get_x_const() by CPL 3.x function call. */
    #define cpl_bivector_get_x_const(bivector) \
        cpl_bivector_get_x(bivector)
    /** @brief  Replaces cpl_bivector_get_y_const() by CPL 3.x function call. */
    #define cpl_bivector_get_y_const(bivector) \
        cpl_bivector_get_y(bivector)
    /** @brief  Replaces cpl_bivector_get_x_data_const() by CPL 3.x function call. */
    #define cpl_bivector_get_x_data_const(bivector) \
        cpl_bivector_get_x_data(bivector)
    /** @brief  Replaces cpl_bivector_get_y_data_const() by CPL 3.x function call. */
    #define cpl_bivector_get_y_data_const(bivector) \
        cpl_bivector_get_y_data(bivector)
    /** @brief  Replaces cpl_frameset_find_const() by CPL 3.x function call. */
    #define cpl_frameset_find_const(frameset) \
        cpl_frameset_find(frameset)
    /** @brief  Replaces cpl_frameset_get_first_const() by CPL 3.x function call. */
    #define cpl_frameset_get_first_const(frameset) \
        cpl_frameset_get_first(frameset)
    /** @brief  Replaces cpl_frameset_get_next_const() by CPL 3.x function call. */
    #define cpl_frameset_get_next_const(frameset) \
        cpl_frameset_get_next(frameset)
    /** @brief  Replaces cpl_frameset_get_frame_const() by CPL 3.x function call. */
    #define cpl_frameset_get_frame_const(frameset) \
        cpl_frameset_get_frame(frameset)
    /** @brief  Replaces cpl_image_get_bpm_const() by CPL 3.x function call. */
    #define cpl_image_get_bpm_const(image) \
        cpl_image_get_bpm(image)
    /** @brief  Replaces cpl_image_get_data_const() by CPL 3.x function call. */
    #define cpl_image_get_data_const(image) \
        cpl_image_get_data(image)
    /** @brief  Replaces cpl_image_get_data_double_const() by CPL 3.x function call. */
    #define cpl_image_get_data_double_const(image) \
        cpl_image_get_data_double(image)
    /** @brief  Replaces cpl_image_get_data_float_const() by CPL 3.x function call. */
    #define cpl_image_get_data_float_const(image) \
        cpl_image_get_data_float(image)
    /** @brief  Replaces cpl_image_get_data_int_const() by CPL 3.x function call. */
    #define cpl_image_get_data_int_const(image) \
        cpl_image_get_data_int(image)
    /** @brief  Replaces cpl_imagelist_get_const() by CPL 3.x function call. */
    #define cpl_imagelist_get_const(imagelist) \
        cpl_imagelist_get(imagelist)
    /** @brief  Replaces cpl_mask_get_data_const() by CPL 3.x function call. */
    #define cpl_mask_get_data_const(mask) \
        cpl_mask_get_data(mask)
    /** @brief  Replaces cpl_matrix_get_data_const() by CPL 3.x function call. */
    #define cpl_matrix_get_data_const(matrix) \
        cpl_matrix_get_data(matrix)
    /** @brief  Replaces cpl_parameterlist_find_const() by CPL 3.x function call. */
    #define cpl_parameterlist_find_const(parameterlist, name) \
        cpl_parameterlist_find(parameterlist, name)
    /** @brief  Replaces cpl_parameterlist_find_context_const() by CPL 3.x function call. */
    #define cpl_parameterlist_find_context_const(parameterlist, context) \
        cpl_parameterlist_find_context(parameterlist, context)
    /** @brief  Replaces cpl_parameterlist_find_tag_const() by CPL 3.x function call. */
    #define cpl_parameterlist_find_tag_const(parameterlist, tag) \
        cpl_parameterlist_find_tag(parameterlist, tag)
    /** @brief  Replaces cpl_parameterlist_find_type_const() by CPL 3.x function call. */
    #define cpl_parameterlist_find_type_const(parameterlist, type) \
        cpl_parameterlist_find_type(parameterlist, type)
    /** @brief  Replaces cpl_parameterlist_get_first_const() by CPL 3.x function call. */
    #define cpl_parameterlist_get_first_const(parameterlist) \
        cpl_parameterlist_get_first(parameterlist)
    /** @brief  Replaces cpl_parameterlist_get_next_const() by CPL 3.x function call. */
    #define cpl_parameterlist_get_next_const(parameterlist) \
        cpl_parameterlist_get_next(parameterlist)
    /** @brief  Replaces cpl_parameterlist_get_last_const() by CPL 3.x function call. */
    #define cpl_parameterlist_get_last_const(parameterlist) \
        cpl_parameterlist_get_last(parameterlist)
    /** @brief  Replaces cpl_propertylist_get_const() by CPL 3.x function call. */
    #define cpl_propertylist_get_const(propertylist) \
        cpl_propertylist_get(propertylist)
    /** @brief  Replaces cpl_propertylist_get_property_const() by CPL 3.x function call. */
    #define cpl_propertylist_get_property_const(propertylist) \
        cpl_propertylist_get_property(propertylist)
    /** @brief  Replaces cpl_propertylist_get_first_const() by CPL 3.x function call. */
    #define cpl_propertylist_get_first_const(propertylist) \
        cpl_propertylist_get_first(propertylist)
    /** @brief  Replaces cpl_propertylist_get_next_const() by CPL 3.x function call. */
    #define cpl_propertylist_get_next_const(propertylist) \
        cpl_propertylist_get_next(propertylist)
    /** @brief  Replaces cpl_table_get_data_int_const() by CPL 3.x function call. */
    #define cpl_table_get_data_int_const(table) \
        cpl_table_get_data_int(table)
    /** @brief  Replaces cpl_table_get_data_float_const() by CPL 3.x function call. */
    #define cpl_table_get_data_float_const(table) \
        cpl_table_get_data_float(table)
    /** @brief  Replaces cpl_table_get_data_double_const() by CPL 3.x function call. */
    #define cpl_table_get_data_double_const(table) \
        cpl_table_get_data_double(table)
    /** @brief  Replaces cpl_table_get_data_string_const() by CPL 3.x function call. */
    #define cpl_table_get_data_string_const(table) \
        cpl_table_get_data_string(table)
    /** @brief  Replaces cpl_table_get_data_array_const() by CPL 3.x function call. */
    #define cpl_table_get_data_array_const(table) \
        cpl_table_get_data_array(table)
    /** @brief  Replaces cpl_vector_get_data_const() by CPL 3.x function call. */
    #define cpl_vector_get_data_const(vector) \
        cpl_vector_get_data(vector)

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Computational functions
     *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    #define cpl_vector_get_median_const(vector) \
        cpl_vector_get_median(vector)
        /**< @brief Replaces cpl_vector_get_median_const()
         *          by CPL 3.x function call. */

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Operational
     *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    #define cpl_init(par) \
        cpl_init()
        /**< @brief Replaces cpl_vector_get_median_const()
         *          by CPL 3.x (parameter-free) function call. */
#endif

/*
 * CPL 4.1 provides mathematical constants. Replace them here if using older
 * CPL.
 */
#ifndef CPL_MATH_CONST_H

    /* The base of the exponential function */
    #define CPL_MATH_E        2.7182818284590452353602874713526624977572470936999595
    
    /* The ratio of a circles circumference to its diameter */
    #define CPL_MATH_PI       3.1415926535897932384626433832795028841971693993751058
    
    /* The natural logarithm of 2 */
    #define CPL_MATH_LN2      0.6931471805599453094172321214581765680755001343602553
    
    /* The natural logarithm of 10 */
    #define CPL_MATH_LN10     2.3025850929940456840179914546843642076011014886287730
    
    /* Derived constants */
    /* 2pi */
    #define CPL_MATH_2PI      6.2831853071795864769252867665590057683943387987502116
    
    /* pi/2 */
    #define CPL_MATH_PI_2     1.5707963267948966192313216916397514420985846996875529
    
    /* pi/4 */
    #define CPL_MATH_PI_4     0.7853981633974483096156608458198757210492923498437765
    
    /* 1/pi */
    #define CPL_MATH_1_PI     0.3183098861837906715377675267450287240689192914809129
    
    /* 2/pi */
    #define CPL_MATH_2_PI     0.6366197723675813430755350534900574481378385829618258
    
    /* 4/pi */
    #define CPL_MATH_4_PI     1.2732395447351626861510701069801148962756771659236516
    
    /* sqrt(2pi) */
    #define CPL_MATH_SQRT2PI  2.5066282746310005024157652848110452530069867406099383
    
    /* 2/sqrt(pi) */
    #define CPL_MATH_2_SQRTPI 1.1283791670955125738961589031215451716881012586579977
    
    /* sqrt(2) */
    #define CPL_MATH_SQRT2    1.4142135623730950488016887242096980785696718753769481
    
    /* sqrt(3) */
    #define CPL_MATH_SQRT3    1.7320508075688772935274463415058723669428052538103806
    
    /* sqrt(1/2) */
    #define CPL_MATH_SQRT1_2  0.7071067811865475244008443621048490392848359376884740
    
    /* log2(e) */
    #define CPL_MATH_LOG2E    1.4426950408889634073599246810018921374266459541529859
    
    /* log10(e) */
    #define CPL_MATH_LOG10E   0.4342944819032518276511289189166050822943970058036666
    
    /* 180/pi */
    #define CPL_MATH_DEG_RAD  57.295779513082320876798154814105170332405472466564322
    
    /* pi/180 */
    #define CPL_MATH_RAD_DEG  0.0174532925199432957692369076848861271344287188854173
    
    /* FWHM per Sigma, 2.0*sqrt(2.0*log(2.0)) */
    #define CPL_MATH_FWHM_SIG 2.3548200450309493820231386529193992754947713787716411
    
    /* Sigma per FWHM, 0.5/sqrt(2.0*log(2.0)) */
    #define CPL_MATH_SIG_FWHM 0.4246609001440095213607514170514448098575705468921770

#endif
/**@}*/

/*----------------------------------------------------------------------------*/
#endif /* CLIPM_COMPATIBILITY_REPLACEMENTS_H */
