/**
 * File:   driver.h
 * Author: akirby
 *
 * Created on February 25, 2023, 12:10 PM
 */
#ifndef DRIVER_H
#define DRIVER_H

/* system header files */
#include <malloc.h>

/* header files*/
#include "input.h"
#include "flow_solver.h"
#include "p4est_overset_interface.h"

#ifndef P4_TO_P8
# include <p4est_extended.h>
#else
# include <p8est_extended.h>
#endif

#define LINEBREAK {int e; printf("+"); for(e = 0; e < 60; ++e) printf("="); printf("+\n");}

#ifdef __cplusplus
extern "C" {
#endif

static inline
double memory_usage(int mpi_rank,int timestep,int display){

    /* get malloc info structure */
    struct mallinfo my_mallinfo = mallinfo();

    /*total memory reserved by the system for malloc currently */
    double reserved_mem = my_mallinfo.arena;

    /* get all the memory currently allocated to user by malloc, etc. */
    double used_mem = my_mallinfo.hblkhd
                    + my_mallinfo.usmblks
                    + my_mallinfo.uordblks;

    /* get memory not currently allocated to user but malloc controls */
    double free_mem = my_mallinfo.fsmblks
                    + my_mallinfo.fordblks;

    /* get number of items currently allocated */
    /* double number_allocated = my_mallinfo.ordblks + my_mallinfo.smblks; */

    /* Print out concise malloc info line */
    if(display && mpi_rank == 0){
        printf("Step[%d]: %f MB(%.0f) malloc: %f MB reserved (%.0f unused)\n",
            timestep,
            used_mem / (1024.0 * 1024.0),
            used_mem,
            reserved_mem / (1024.0 * 1024.0),
            free_mem);

 	if(mpi_rank == 0){
            FILE *fp;
	    char filename[] = "new_memusage.dat";
    	    fp=fopen(filename,"a+");
            fprintf(fp,"Step[%d]: %f MB(%.0f) malloc: %f MB reserved (%.0f unused)\n",
	            timestep,used_mem / (1024.0 * 1024.0),used_mem,reserved_mem / (1024.0 * 1024.0), free_mem);
            fclose(fp);
  	}

    }
    return used_mem;
}

static inline
int file_is_modified(const char *path,time_t oldMTime){
    struct stat file_stat;
    int err = stat(path, &file_stat);
    if (err != 0) {
        perror(" [file_is_modified] stat");
        exit(errno);
    }
    return file_stat.st_mtime > oldMTime;
}

static inline
char *trimwhitespace(char *str){
    char *end;

    // Trim leading space
    while(isspace((unsigned char)*str)) str++;

    if(*str == 0) return str;  // All spaces?

    // Trim trailing space
    end = str + strlen(str) - 1;
    while(end > str && isspace((unsigned char)*end)) end--;

    // Write new null terminator character
    end[1] = '\0';
    return str;
}


/* the main global storage struct that holds everything including flow solvers and p4est_overset */
class driver_t {
public:
  char input_file[buff_size]; /** input file name */
  time_t input_file_mod_time; /** input file modification time */

  /** input file variables */
  int ngroups;                /**< total number of groups */
  int nsolvers;               /**< total number of solvers */
  int off_body_group;         /**< identifies the off body group index */

  /** group identification for all processes */
  int *group_solver_id;       /**< solver id for each group */
  int *group_num_procs;       /**< number of procs for each group */
  int *mesh_rank_offsets;     /**< mesh rank offsets */

  /** local group identification */
  int group_index;            /**< my local group number indexing from 0 */

  /** mpi information for all processes */
  MPI_Comm mpicomm;
  int rank;                   /**< mpi rank across all cores */
  int global_leader_flag;     /**< if leader rank, then flag=1,else=0 */
  int num_procs;              /**< total number of processors */

  /** mpi information for group */
  int group_rank;             /**< group rank */
  int group_leader_flag;      /**< if group master rank then flag=1,else=0 */
  MPI_Comm group_comm;        /**< mpi communicator for the group */
  MPI_Comm leader_comm;       /**< communicator between mesh leader ranks

  /** off body code identification */
  int off_body_group_flag;    /**< if a group is an off body group, then = 1, else 0 */
  int off_body_mode_flag;     /**< if any group is an off body group, then mode = 1, else = 0 */

  /** pointers to class data */
  flow_solver *flow; /**< a pointer to an array of flow solver data structs  */
  p4est_overset_interface *p4est_overset; /**< a pointer to an array of p4est_overset data structs */

  int *recvmap;
  int *sendmap;

/* constructors */
  driver_t() :
    recvmap{NULL},
    sendmap{NULL} {}

 ~driver_t(){
     if(recvmap) delete [] recvmap;
     if(sendmap) delete [] sendmap;
 }

public:
  void initialize(int argc,char **argv);
  void overset_test();
  void finalize();

private:
  void initialize_mpi(int argc,char **argv);
  void create_comms();
  void read_input_file(int argc, char **argv);
  void move_into_group_directory(int g){
    char string[12];
    snprintf(string,sizeof(string),"group%d",g);
    (void) chdir(string);
  }
};

#ifdef __cplusplus
}
#endif
#endif /* DRIVER_H */