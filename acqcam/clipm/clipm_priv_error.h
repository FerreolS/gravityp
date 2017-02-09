
/*********************************************************************
 * E.S.O. - VLT project
 *
 * "@(#) $Id: clipm_priv_error.h 171356 2008-08-08 15:35:24Z hlorch $"
 *
 * PRIVATE functions for clipm internal error handling, requires CPL 4.0
 *
 * who       when        what
 * --------  ----------  ----------------------------------------------
 * hlorch    2006-05-11  created
 */

/**
 * @internal
 * @defgroup clipm_priv_error   Error Handling Macros
 * @ingroup internal_docs
 * 
 * This module provides macros for the handling of CPL errors and the error
 * state.
 * 
 * @par General Rules:
 * 
 * - The macros below take care of setting the right information about the
 *   location where an error happened. Functions like cpl_error_set_where()
 *   don't need to be called by the programmer using them.
 * 
 * @par Synopsis:
 * @code
 *   #include "clipm_priv_error.h"
 * @endcode
 */
/** @{ */

#ifndef CLIPM_PRIV_ERROR_H
#define CLIPM_PRIV_ERROR_H

/*-----------------------------------------------------------------------------
    Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

#include <string.h>

/*-----------------------------------------------------------------------------
    Functional Macros
 -----------------------------------------------------------------------------*/

#define _CLIPM_ERROR_SET_WHERE_() \
do { \
    if (CLIPM_ERROR_IS_SET()) \
    { \
        const char      *msg; \
        int             n = 0; \
        \
        msg = cpl_error_get_message(); \
        /* search the beginning of the last custom message */ \
        while (msg[n] != '\0' && msg[n] != ':') \
            n++; \
        while (msg[n] == ':' || msg[n] == ' ') \
            n++; \
        cpl_error_set_message(              __func__, \
                                            cpl_error_get_code(), \
                                            msg+n); \
        clipm_error_is_set_where = 1; \
    } \
} while (0)

#define _CLIPM_ERROR_SET_MSG_(code, object, msg) \
do { \
    char _clipm_error_msg[256]; \
    _clipm_priv_error_sprint_messages(      _clipm_error_msg, \
                                            object, \
                                            msg, \
                                            255); \
    cpl_error_set_message(                  __func__, \
                                            code, \
                                            _clipm_error_msg); \
    clipm_error_is_set_where = 1; \
} while (0)

/*----------------------------------------------------------------------------*/
/**
 * @brief Beginning of a TRY-block.
 * @hideinitializer
 * 
 * The macro CLIPM_TRY is to be used like a keyword in front of a deeper
 * scope. This scope has to be followed by the macro @ref CLIPM_CATCH. This
 * means that CLIPM_TRY and CLIPM_CATCH build a frame around a code
 * statement or a code scope, called the try-block.
 * 
 * The CLIPM_CATCH macro is to be followed by a statement or scope, which is
 * only executed if a CPL error is set while reaching the CLIPM_CATCH macro,
 * called the catch-block.
 * 
 * The try-block can be exited by using one of the macros below, for
 * example with @ref CLIPM_TRY_EXIT_WITH_ERROR(). In this case, a jump to
 * CLIPM_CATCH is performed, and the catch-block executed if an error is set.
 * 
 * @note
 * 
 * The following constraints have to be fulfilled:
 * - A "return" or "goto" statement inside the try-block is forbidden, because
 *   leaving the try-block without processing the CLIPM_CATCH macro will mess
 *   up the error state information. In the catch-block (which comes after
 *   the CLIPM_CATCH macro), it is allowed.
 * - The macros require some variables, which are declared at the beginning
 *   of the CLIPM_TRY macro. Therefore it is not possible in ANSI-C to have
 *   code statements (except declarations) before the CLIPM_TRY macro. If it is
 *   required, this can be solved by putting a scope around the try-catch
 *   construct.
 * - Only one CLIPM_TRY - CLIPM_CATCH - construct can be inside one function.
 * 
 * @par Example 1:
 * 
 * @code
cpl_error_code my_func()
{
    cpl_object      *obj = NULL;
    
    CLIPM_TRY
    {
        CLIPM_TRY_EXIT_IFN(
            obj = cpl_object_new());
        
        CLIPM_TRY_EXIT_IFN(
            cpl_function(obj) == CPL_ERROR_NONE);
        
        cpl_function(obj);
        CLIPM_TRY_CHECK_ERROR_STATE();
    }
    CLIPM_CATCH
    {
    }

    cpl_object_delete(obj);
    
    return CLIPM_ERROR_GET_NEW_SINCE_TRY;
}
 * @endcode
 * 
 * @par Example 2:
 * 
 * @code
cpl_object  *my_func()
{
    cpl_object      *obj = NULL;
    
    CLIPM_TRY
    {
        CLIPM_TRY_EXIT_IFN(
            obj = cpl_object_new());
        
        CLIPM_TRY_EXIT_IFN(
            cpl_function(obj) == CPL_ERROR_NONE);
        
        cpl_function(obj);
        CLIPM_TRY_CHECK_ERROR_STATE();
    }
    CLIPM_CATCH
    {
        cpl_object_delete(obj);
        obj = NULL;
    }
    
    return obj;
}
 * @endcode
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_TRY \
    int             clipm_error_catch_call_flag, \
                    clipm_error_is_set_where = 0; \
    cpl_errorstate  clipm_error_trystate; \
    \
    clipm_error_trystate = cpl_errorstate_get(); \
    \
    do

/*----------------------------------------------------------------------------*/
/**
 * @brief   End of a TRY-block, beginning of a CATCH-block.
 * @hideinitializer
 * 
 * Please refer to @ref CLIPM_TRY.
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_CATCH \
    while (0); \
    \
    goto _CLIPM_CATCH_LABEL_; /* avoid warning if not used */ \
    _CLIPM_CATCH_LABEL_: \
    \
    if ((clipm_error_catch_call_flag = (CLIPM_ERROR_IS_SET()))) \
    { \
        if (!clipm_error_is_set_where) \
            _CLIPM_ERROR_SET_WHERE_(); \
    } \
    \
    if (clipm_error_catch_call_flag)

/*----------------------------------------------------------------------------*/
/**
 * @brief   Return new CPL error code
 * @hideinitializer
 * 
 * @return  If the CPL error state has changed since CLIPM_TRY, the latest
 *          error code is returned, otherwise CPL_ERROR_NONE.
 * 
 * - May be called outside TRY-block.
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_ERROR_GET_NEW_SINCE_TRY(void) \
    (CLIPM_ERROR_IS_NONE() ? CPL_ERROR_NONE : cpl_error_get_code())

/*----------------------------------------------------------------------------*/
/**
 * @brief   Recover the error state which was present during CLIPM_TRY (at
 *          the beginning of the try-block).
 * @hideinitializer
 * 
 * - May be called outside TRY-block.
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_ERROR_RECOVER_TRYSTATE(void) \
do { \
    cpl_errorstate_set(clipm_error_trystate); \
    clipm_error_is_set_where = 0; \
} while (0)

/*----------------------------------------------------------------------------*/
/**
 * @brief   Set a new error code.
 * @hideinitializer
 * 
 * - @a code must not be CPL_ERROR_NONE.
 * - May be called outside TRY-block.
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_ERROR_SET(code) \
do { \
    CLIPM_TRY_ASSERT(code != CPL_ERROR_NONE); \
    cpl_error_set(__func__, code); \
    clipm_error_is_set_where = 1; \
} while (0)

/*----------------------------------------------------------------------------*/
/**
 * @brief   Set a new error code together with a custom error message.
 * @hideinitializer
 *
 * - @a code must not be CPL_ERROR_NONE.
 * - May be called outside TRY-block.
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_ERROR_SET_MSG(code, object, msg) \
do { \
    CLIPM_TRY_ASSERT(code != CPL_ERROR_NONE); \
    _CLIPM_ERROR_SET_MSG_(code, object, msg); \
} while (0)

/*----------------------------------------------------------------------------*/
/**
 * @brief   Set a new custom error message if a certain error code is already
 *          set.
 * @hideinitializer
 *
 * - @a code must not be CPL_ERROR_NONE.
 * - May be called outside TRY-block.
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_ERROR_SET_MSG_IF_CODE(code, object, msg) \
do { \
    CLIPM_TRY_ASSERT(code != CPL_ERROR_NONE); \
    if (CLIPM_ERROR_GET_NEW_SINCE_TRY() == code) \
        _CLIPM_ERROR_SET_MSG_(code, object, msg); \
} while (0)

/*----------------------------------------------------------------------------*/
/**
 * @brief   Return if a new CPL error is set
 * @hideinitializer
 * 
 * @return  If the CPL error state has changed since CLIPM_TRY, 1 returned,
 *          otherwise 0.
 * 
 * - May be called outside TRY-block.
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_ERROR_IS_SET(void) \
    (!cpl_errorstate_is_equal(clipm_error_trystate))

/*----------------------------------------------------------------------------*/
/**
 * @brief   Return if no new CPL error is set
 * @hideinitializer
 * 
 * @return  If the CPL error state has changed since CLIPM_TRY, 0 returned,
 *          otherwise 1.
 * 
 * - May be called outside TRY-block.
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_ERROR_IS_NONE(void) \
    (cpl_errorstate_is_equal(clipm_error_trystate))

/*----------------------------------------------------------------------------*/
/**
 * @brief   Assure the condition is true, else set the respective error code,
 *          exit the TRY block, and set an error message using the object name
 *          (can be empty string) and a message.
 * @param   condition   Condition
 * @param   code        Error code to be set
 * @param   object      Object name (can be empty or NULL)
 * @param   msg         Message (can be empty or NULL)
 * @hideinitializer
 * 
 * - <b>Only</b> allowed in TRY-block, forbidden in CATCH-block!
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_TRY_CHECK(condition, code, object, msg) \
do { \
    if (!(condition)) \
    { \
        CLIPM_ERROR_SET_MSG( (code), (object), (msg)); \
        CLIPM_TRY_EXIT(); \
    } \
} while (0)

/*----------------------------------------------------------------------------*/
/**
 * @brief   Assure the condition is true, else set the respective error code,
 *          exit the TRY block, and auto-generate an error message (re-using
 *          the condition).
 * @param   condition   Condition
 * @param   code        Error code to be set
 * @hideinitializer
 * 
 * - <b>Only</b> allowed in TRY-block, forbidden in CATCH-block!
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_TRY_CHECK_AUTOMSG(condition, code) \
do { \
    if (!(condition)) \
    { \
        CLIPM_ERROR_SET_MSG( (code), "!("#condition")", ""); \
        CLIPM_TRY_EXIT(); \
    } \
} while (0)

/*----------------------------------------------------------------------------*/
/**
 * @brief   Check the CPL error state, and exit the try-block if not
 *          CPL_ERROR_NONE.
 * @hideinitializer
 * 
 * - Does not overwrite locally set messages.
 * - <b>Only</b> allowed in TRY-block, forbidden in CATCH-block!
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_TRY_CHECK_ERROR_STATE(void) \
do { \
    if (CLIPM_ERROR_IS_SET()) \
    { \
        if (!clipm_error_is_set_where) \
            _CLIPM_ERROR_SET_WHERE_(); \
        CLIPM_TRY_EXIT(); \
    } \
} while (0)

/*----------------------------------------------------------------------------*/
/**
 * @brief   Assert that the given condition is fulfilled, otherwise
 *          set a custom bugreport message and exit the TRY-block with error
 *          code CLIPM_ERROR_UNEXPECTED.
 * @hideinitializer
 * 
 * - <b>Only</b> allowed in TRY-block, forbidden in CATCH-block!
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_TRY_ASSERT(condition) \
do { \
    if (!(condition)) \
    { \
        _CLIPM_ERROR_SET_MSG_(              CLIPM_ERROR_UNEXPECTED, \
                                            "!("#condition")", \
                                            _CLIPM_MSG_ERR_UNEXPECTED); \
        CLIPM_TRY_EXIT(); \
    } \
} while (0)

/*----------------------------------------------------------------------------*/
/**
 * @brief   Assert that the CPL error state is CPL_ERROR_NONE, otherwise
 *          set a custom bugreport message and exit the TRY-block with error
 *          code CLIPM_ERROR_UNEXPECTED.
 * @hideinitializer
 * 
 * - <b>Only</b> allowed in TRY-block, forbidden in CATCH-block!
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_TRY_ASSERT_ERROR_STATE(void) \
do { \
    if (CLIPM_ERROR_IS_SET()) \
    { \
        _CLIPM_ERROR_SET_MSG_(              CLIPM_ERROR_UNEXPECTED, \
                                            NULL, \
                                            _CLIPM_MSG_ERR_UNEXPECTED); \
        CLIPM_TRY_EXIT(); \
    } \
} while (0)

/*----------------------------------------------------------------------------*/
/**
 * @brief   Set a new CPL error, and exit the try-block.
 * @hideinitializer
 * 
 * - @a code must not be CPL_ERROR_NONE!
 * - <b>Only</b> allowed in TRY-block, forbidden in CATCH-block!
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_TRY_EXIT_WITH_ERROR(code) \
do { \
    CLIPM_ERROR_SET((code)); \
    CLIPM_TRY_EXIT(); \
} while (0)

/*----------------------------------------------------------------------------*/
/**
 * @brief   Set a new CPL error together with a custom error message,
 *          and exit the try-block.
 * @hideinitializer
 * 
 * - @a code must not be CPL_ERROR_NONE!
 * - Overwrites any message!
 * - <b>Only</b> allowed in TRY-block, forbidden in CATCH-block!
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_TRY_EXIT_WITH_ERROR_MSG(code, object, msg) \
do { \
    CLIPM_ERROR_SET_MSG((code), (object), (msg)); \
    CLIPM_TRY_EXIT(); \
} while (0)

/*----------------------------------------------------------------------------*/
/**
 * @brief   If @a condition == 0, then
 *          the try-block is exited.
 * @hideinitializer
 * 
 * - It is assumed, that a new error state is already set if @a condition is 0,
 *   otherwise a custom bugreport message is set and the TRY-block is exited
 *   with error code CLIPM_ERROR_UNEXPECTED. If it is desired to instead
 *   eventually set an error here, consider using CLIPM_TRY_CHECK[_AUTOMSG]
 *   instead.
 * - Does not overwrite locally set messages.
 * - <b>Only</b> allowed in TRY-block, forbidden in CATCH-block!
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_TRY_EXIT_IFN(condition) \
do { \
    if (!(condition)) \
    { \
        CLIPM_TRY_ASSERT(CLIPM_ERROR_IS_SET()); \
        if (!clipm_error_is_set_where) \
            _CLIPM_ERROR_SET_WHERE_(); \
        CLIPM_TRY_EXIT(); \
    } \
} while (0)

/*----------------------------------------------------------------------------*/
/**
 * @brief   The try-block is exited.
 * @hideinitializer
 * 
 * - It is not necessary that a new error state is already set. In this case,
 *   the try-block is just exited.
 * - <b>Only</b> allowed in TRY-block, forbidden in CATCH-block!
 */
/*----------------------------------------------------------------------------*/
#define CLIPM_TRY_EXIT(void) \
    goto _CLIPM_CATCH_LABEL_


/*-----------------------------------------------------------------------------
    Declarations
 -----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------------
    Types
 -----------------------------------------------------------------------------*/

/**
 * @brief   Extension to CPL error codes.
 */

typedef enum _clipm_error_code_ {
    CLIPM_ERROR_UNEXPECTED = CPL_ERROR_EOL + 0
    /**< Unexpected internal error */
} clipm_error_code;

/** @brief  Internal error */
extern const char _CLIPM_MSG_ERR_UNEXPECTED[];
/** @brief  Internal error handling bug */
extern const char _CLIPM_MSG_ERR_HANDLING[];

/** @brief  Location matrix must contain 2 rows */
extern const char CLIPM_MSG_ERR_2ROWXY[];
/** @brief  Location matrices differ in size */
extern const char CLIPM_MSG_ERR_DIFFSIZES[];
/** @brief  Location matrices differ in size */
extern const char CLIPM_MSG_ERR_DIFFTYPES[];

/** @} Doxygen group end */
/*-----------------------------------------------------------------------------
    Prototypes
 -----------------------------------------------------------------------------*/

void        _clipm_priv_error_sprint_messages(
                                            char        *outstr,
                                            const char  *msg1,
                                            const char  *msg2,
                                            int         maxlen);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}   /* extern "C" */
#endif

#endif /* CLIPM_PRIV_ERROR_H */
