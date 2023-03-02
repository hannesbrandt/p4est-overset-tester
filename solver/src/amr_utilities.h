/**
 * \file   amr_utilities.h
 * \author akirby
 */

#ifndef AMR_UTILITIES_H
#define AMR_UTILITIES_H

/* header files */
#include "amr_var_defines.h"
#include "dg4est_solver.hxx"
#include "precision_types.h"

#include <mpi.h>

#ifndef DOXYGEN_IGNORE
#  include <math.h>
#  include <ctype.h>
#  include <unistd.h>
#  include <string.h>
#  include <stdlib.h>
#  include <stdio.h>
#  include <errno.h>
#  include <sys/stat.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Display AMR input options from file to terminal
 *
 * @param [in] ctx        context data related to simulation
 * @param [in] mpi_rank     mpi rank
 */
void amr_utilities_inputs_amr_message(ctx_t *ctx, int mpi_rank);

/** Display AMR grid message to terminal
 *
 * @param [in] ctx        context data related to simulation
 */
void amr_utilities_grid_message(ctx_t *ctx);

/** Display AMR regrid timings
 *
 * @param [in] mpi_rank         mpi rank
 * @param [in] max_time         total regrid time
 * @param [in] balance_time     total amr 2:1 balancing time
 * @param [in] partition_time   total mpi partitioning time
 * @param [in] ship_time        total custom data move time
 * @param [in] max_point_time   total point regridding time
 * @param [in] max_feature_time total feature-based regridding time
 * @param [in] max_spread_time  total feature-based tag spreading time
 * @param [in] max_coarsen_time total coarsening regrid time
 * @param [in] max_buffer_time  total buffer regridding time
 * @param [in] max_level_time   total level-based regrid time
 * @param [in] pdegree_time     total change p-degree time
 * @param [in] log_info         logging threshold
 */
void amr_utilities_write_amr_regrid(ctx_t *ctx,int mpi_rank,Real max_time,
                                    Real balance_time,Real partition_time,
                                    Real ship_time,Real max_point_time,
                                    Real max_feature_time,Real max_spread_time,
                                    Real max_coarsen_time,Real max_buffer_time,
                                    Real max_level_time,Real pdegree_time,
                                    int log_info);

/** MPI wallclock timer wrapper function
 *
 * @param [out] time_out  mpi wall time
 */
Real amr_utilities_timer();

/** MPI reduce Real max function wrapper
 *
 * @param [in] mpicomm      mpi communicator
 * @param [in] Real_in    Real precision number to find global max value
 * @return  global max Real value
 */
Real amr_utilities_mpireducemax_real(int rank,MPI_Comm mpicomm,Real Real_in);

/** Build the output directories for the simulation
 *
 * @param [in] mpi_rank     mpi rank
 * @param [in] mpicomm      mpi communicator
 * @param [in] log_info     log threshold
 */
void amr_utilities_create_directories(int mpi_rank,MPI_Comm mpicomm,int log_info);

/** Input file read line helper function: checks for keyword
 *
 * @param [in]  filename    file to read
 * @param [in]  keyword     keyword to find in file
 * @return                  returns 0 if successfully read variable
 */
char amr_utilities_find_keyword(char *filename,const char *keyword);

/** Input file read line helper function: one integer
 *
 * @param [in]  filename    file to read
 * @param [in]  keyword     keyword to find in file
 * @param [out] integer     integer variable to read
 * @return                  returns 0 if successfully read variable
 */
char amr_utilities_find_keyword_integer(char *filename,const char *keyword,
                                        int *integer);

/** Input file read line helper function: two integers
 *
 * @param [in]  filename    file to read
 * @param [in]  keyword     keyword to find in file
 * @param [out] integer1    integer variable #1 to read
 * @param [out] integer2    integer variable #2 to read
 * @return                  returns 0 if successfully read variable
 */
char amr_utilities_find_keyword_two_integers(char *filename,const char *keyword,
                                             int *integer1,
                                             int *integer2);

/** Input file read line helper function: three integers
 *
 * @param [in]  filename    file to read
 * @param [in]  keyword     keyword to find in file
 * @param [out] integer1    integer variable #1 to read
 * @param [out] integer2    integer variable #2 to read
 * @param [out] integer3    integer variable #3 to read
 * @return                  returns 0 if successfully read variable
 */
char amr_utilities_find_keyword_three_integers(char *filename,const char *keyword,
                                               int *integer1,
                                               int *integer2,
                                               int *integer3);

/** Input file read line helper function: one real
 *
 * @param [in]  filename    file to read
 * @param [in]  keyword     keyword to find in file
 * @param [out] dbl         Real precision variable to read
 * @return                  returns 0 if successfully read variable
 */
char amr_utilities_find_keyword_real(char *filename,const char *keyword,
                                     Real *dbl);

/** Input file read line helper function: two reals
 *
 * @param [in]  filename    file to read
 * @param [in]  keyword     keyword to find in file
 * @param [out] dbl1        Real precision variable #1 to read
 * @param [out] dbl2        Real precision variable #2 to read
 * @return                  returns 0 if successfully read variable
 */
char amr_utilities_find_keyword_two_reals(char *filename,const char *keyword,
                                          Real *dbl1,
                                          Real *dbl2);

/** Input file read line helper function: three reals
 *
 * @param [in]  filename    file to read
 * @param [in]  keyword     keyword to find in file
 * @param [out] dbl1        Real precision variable #1 to read
 * @param [out] dbl2        Real precision variable #2 to read
 * @param [out] dbl3        Real precision variable #3 to read
 * @return                  returns 0 if successfully read variable
 */
char amr_utilities_find_keyword_three_reals(char *filename,const char *keyword,
                                            Real *dbl1,
                                            Real *dbl2,
                                            Real *dbl3);

/** Input file read line helper function: one integer and one real
 *
 * @param [in]  filename    file to read
 * @param [in]  keyword     keyword to find in file
 * @param [out] int1        integer variable #1 to read
 * @param [out] dbl1        Real precision variable #1 to read
 * @return                  returns 0 if successfully read variable
 */
char amr_utilities_find_keyword_int_real(char *filename,const char *keyword,
                                         int *int1,
                                         Real *dbl1);

/** Input file read line helper function: string
 *
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] string      string to read
 * @return                  returns 0 if successfully read variable
 */
char amr_utilities_find_keyword_string(char *filename,const char *keyword,
                                       char *string);

/** Input file read line helper function: string (converted to upper case)
 *
 * @param [in] filename     file to read
 * @param [in] keyword      keyword to find in file
 * @param [out] string      string to read
 * @return                  returns 0 if successfully read variable
 */
char amr_utilities_find_keyword_string2upper(char *filename,const char *keyword,
                                             char *string);

/** Input file read line helper function: string
 *
 * @param [in] filename     file to read
 * @param [in] upper_flag   flag to change to upper case for check
 * @param [in] keyword      keyword to find in file
 * @param [out] string      string to read
 * @return                  returns 0 if successfully read variable
 */
char amr_utilities_find_keyword_string_optional(char *filename,int upper_flag,
                                                const char *keyword,char *string);

/** Input file modification check
 *
 * @param [in] filepath     file to read
 * @param [in] oldMTime     modification time
 * @return modification flag
 */
int file_is_modified(const char *path,time_t oldMTime);

#ifdef __cplusplus
}
#endif
#endif /* AMR_UTILITIES_H */