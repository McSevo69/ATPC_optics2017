/**\file */
#ifndef SLIC_DECLARATIONS_Vectors_H
#define SLIC_DECLARATIONS_Vectors_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define Vectors_maxN (16384)


/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/




/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] param_dataSize Interface Parameter "dataSize".
 * \param [in] instream_inX The stream should be of size (param_dataSize * 4) bytes.
 * \param [in] instream_inX4 The stream should be of size (param_dataSize * 4) bytes.
 * \param [in] instream_inX5 The stream should be of size (param_dataSize * 4) bytes.
 * \param [in] instream_inX6 The stream should be of size (param_dataSize * 4) bytes.
 * \param [in] instream_inX7 The stream should be of size (param_dataSize * 4) bytes.
 * \param [in] instream_inX8 The stream should be of size (param_dataSize * 4) bytes.
 * \param [in] instream_inY The stream should be of size (param_dataSize * 4) bytes.
 * \param [in] instream_inZ The stream should be of size (param_dataSize * 4) bytes.
 * \param [out] outstream_outEuclidian The stream should be of size ((param_dataSize * param_dataSize) * 4) bytes.
 */
void Vectors(
	uint64_t param_dataSize,
	const float *instream_inX,
	const float *instream_inX4,
	const float *instream_inX5,
	const float *instream_inX6,
	const float *instream_inX7,
	const float *instream_inX8,
	const float *instream_inY,
	const float *instream_inZ,
	float *outstream_outEuclidian);

/**
 * \brief Basic static non-blocking function for the interface 'default'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_dataSize Interface Parameter "dataSize".
 * \param [in] instream_inX The stream should be of size (param_dataSize * 4) bytes.
 * \param [in] instream_inX4 The stream should be of size (param_dataSize * 4) bytes.
 * \param [in] instream_inX5 The stream should be of size (param_dataSize * 4) bytes.
 * \param [in] instream_inX6 The stream should be of size (param_dataSize * 4) bytes.
 * \param [in] instream_inX7 The stream should be of size (param_dataSize * 4) bytes.
 * \param [in] instream_inX8 The stream should be of size (param_dataSize * 4) bytes.
 * \param [in] instream_inY The stream should be of size (param_dataSize * 4) bytes.
 * \param [in] instream_inZ The stream should be of size (param_dataSize * 4) bytes.
 * \param [out] outstream_outEuclidian The stream should be of size ((param_dataSize * param_dataSize) * 4) bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Vectors_nonblock(
	uint64_t param_dataSize,
	const float *instream_inX,
	const float *instream_inX4,
	const float *instream_inX5,
	const float *instream_inX6,
	const float *instream_inX7,
	const float *instream_inX8,
	const float *instream_inY,
	const float *instream_inZ,
	float *outstream_outEuclidian);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	uint64_t param_dataSize; /**<  [in] Interface Parameter "dataSize". */
	const float *instream_inX; /**<  [in] The stream should be of size (param_dataSize * 4) bytes. */
	const float *instream_inX4; /**<  [in] The stream should be of size (param_dataSize * 4) bytes. */
	const float *instream_inX5; /**<  [in] The stream should be of size (param_dataSize * 4) bytes. */
	const float *instream_inX6; /**<  [in] The stream should be of size (param_dataSize * 4) bytes. */
	const float *instream_inX7; /**<  [in] The stream should be of size (param_dataSize * 4) bytes. */
	const float *instream_inX8; /**<  [in] The stream should be of size (param_dataSize * 4) bytes. */
	const float *instream_inY; /**<  [in] The stream should be of size (param_dataSize * 4) bytes. */
	const float *instream_inZ; /**<  [in] The stream should be of size (param_dataSize * 4) bytes. */
	float *outstream_outEuclidian; /**<  [out] The stream should be of size ((param_dataSize * param_dataSize) * 4) bytes. */
} Vectors_actions_t;

/**
 * \brief Advanced static function for the interface 'default'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Vectors_run(
	max_engine_t *engine,
	Vectors_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'default'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Vectors_run_nonblock(
	max_engine_t *engine,
	Vectors_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'default'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Vectors_run_group(max_group_t *group, Vectors_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Vectors_run_group_nonblock(max_group_t *group, Vectors_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'default'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Vectors_run_array(max_engarray_t *engarray, Vectors_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Vectors_run_array_nonblock(max_engarray_t *engarray, Vectors_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Vectors_convert(max_file_t *maxfile, Vectors_actions_t *interface_actions);

/**
 * \brief Initialise a maxfile.
 */
max_file_t* Vectors_init(void);

/* Error handling functions */
int Vectors_has_errors(void);
const char* Vectors_get_errors(void);
void Vectors_clear_errors(void);
/* Free statically allocated maxfile data */
void Vectors_free(void);
/* These are dummy functions for hardware builds. */
int Vectors_simulator_start(void);
int Vectors_simulator_stop(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SLIC_DECLARATIONS_Vectors_H */

