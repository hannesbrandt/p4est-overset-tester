/**
 * File:   driver.cxx
 * Author: akirby
 *
 * Created on February 25, 2023, 11:40 AM
 */

/* header files */
#include "driver.h"

/* ======================= */
/* PUBLIC CLASS FUNCTIONS */
/* ======================= */
void driver_t::initialize(int argc,char **argv){
    p4est_overset = new p4est_overset_interface();

    /* initialize mpi */
    initialize_mpi(argc,argv);

    /* instantiate flow solver */
    flow = new flow_solver(mpicomm,rank);

    /* read driver inputs */
    read_input_file(argc,argv);

    /* create MPI groups and communicators */
    create_comms();

    /* load the function pointers from the dynamic libraries */
    flow->flow_load_dynamic_library();

    /* move into the group folders, this assumes the folders exist already */
    move_into_group_directory(group_index);

    /* initialize flow solver mpi */
    flow->flow_initialize_group_mpi(group_comm);

    /* initialize flow solver */
    flow->initialize();

    /* get pointers from flow solver */
    flow->flow_set_pointers(group_index);

    /* initialize p4est_overset: p4est bg mesh must be group_index=0 */
    p4est_overset->p4est_overset_init(mpicomm,leader_comm,group_comm,
                                      group_index,ngroups,
                                      mesh_rank_offsets,
                                      flow->p4est,
                                      flow->qpoints);
}

void driver_t::overset_test(){
    int ind = 0;
  //p4est_overset_register_and_connect();
    flow->flow_output_solution(ind);
}

void driver_t::finalize(){
    if(group_solver_id) delete [] group_solver_id;
    if(group_num_procs) delete [] group_num_procs;
    if(mesh_rank_offsets) delete [] mesh_rank_offsets;
    if(recvmap) delete [] recvmap;
    if(sendmap) delete [] sendmap;
    MPI_Finalize();
}

/* ======================= */
/* PRIVATE CLASS FUNCTIONS */
/* ======================= */
void driver_t::initialize_mpi(int argc,char **argv){
    mpicomm = MPI_COMM_WORLD;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(mpicomm,&rank);
    MPI_Comm_size(mpicomm,&num_procs);
}

void driver_t::create_comms(){
    /* split communicator for the different meshes */
    MPI_Comm_split(mpicomm, group_index, rank,&group_comm);
    MPI_Comm_rank(group_comm,&group_rank);

    /* set leader flags */
    global_leader_flag = (rank == 0) ? 1:0;
    group_leader_flag = (group_rank == 0) ? 1:0;

    /* set off-body group flag */
    off_body_group_flag = (off_body_group == group_index) ? 1:0;
    MPI_Allreduce(&off_body_group_flag,
                  &off_body_mode_flag,
                  1,MPI_INT,MPI_MAX,
                  mpicomm);

    /* create new communicator between group leader ranks */
    MPI_Comm_split(mpicomm,
                   group_leader_flag ? 0 : MPI_UNDEFINED,
                   rank,
                  &leader_comm);
}

void driver_t::read_input_file(int argc, char **argv){
    char default_filename[] = "input.pdriver";
    char *filename = (argc < 2) ? default_filename:argv[1];

    char keyword[buff_size];
    char cwd[buff_size];
    struct stat file_stat;
    int sum;
    int err;
    int i;

    /* check if file exists */
    if(access(filename,F_OK) == -1){
        printf("[driver] input file does not exist...stopping\n");
        exit(1);
    }

    if(getcwd(cwd,sizeof(cwd)) != NULL){
        char *pwd = trimwhitespace(cwd);

        /* save input file name with full path*/
        strcpy(input_file,pwd);
        strcat(input_file,"/");
        strcat(input_file,basename(filename));
    }

    /* file state for checking runtime modifications */
    stat(filename,&file_stat);
    input_file_mod_time = file_stat.st_mtime;

    /* mandatory inputs */
    err = find_keyword_integer(filename, "number_of_groups:", &ngroups, 1);
    err+= find_keyword_integer(filename, "number_of_solvers:", &nsolvers, 1);
    if(err) printf("\033[1;31m"
                   "[driver] ERROR input error for number_of_groups or number_of_solvers\n"
                   "\033[0m\n");

    if(ngroups < 1 && rank == 0){
        printf("\033[1;31m"
               "[driver] ERROR input error need at least one group %d\n"
               "\033[0m\n",ngroups);
    }

    if (nsolvers < 1 && rank == 0) {
        printf("\033[1;31m"
               "[driver] ERROR input error need at least one solver %d\n"
               "\033[0m\n",nsolvers);
    }

    off_body_group = 0;
    flow->translation[0] = 0.0;
    flow->translation[1] = 0.0;
    flow->translation[2] = 0.0;

    /* non-mandatory inputs */
    int ob_group;
    find_keyword_integer(filename,"off_body_group:", &ob_group,0);
    if(ob_group != off_body_group){
        printf("\033[1;31m"
               "[driver] ERROR input error off_body_group needs to be group 0\n"
               "\033[0m\n");
        exit(1);
    }

    group_solver_id = new int[ngroups];
    group_num_procs = new int[ngroups];
    mesh_rank_offsets = new int[ngroups+1];

    for(i=0; i<ngroups; i++){
        sprintf(keyword, "group%d:",i);
        find_keyword_two_integers(filename, keyword,
                                  &group_solver_id[i],
                                  &group_num_procs[i],
                                  1);
    }

    mesh_rank_offsets[0] = 0;
    for(i=1; i<=ngroups; ++i) mesh_rank_offsets[i] = mesh_rank_offsets[i-1]+group_num_procs[i-1];
    for(i=0,sum=0; i<ngroups; ++i) sum += group_num_procs[i];

    if(sum!=num_procs && rank==0){
        printf("\033[1;31m"
               "[driver] --WARNING-- Total ranks (%d_ does not match sum of group ranks (%d)\n"
               "\033[0m\n",num_procs,sum);
    }

    for(i=0,sum=0; i<ngroups; i++){
        if(rank >= sum && rank < sum+group_num_procs[i]) group_index = i;
        sum+=group_num_procs[i];
    }

    sprintf(keyword, "translation%d:", group_index);
    find_keyword_three_doubles(filename,keyword,
                               &flow->translation[0],
                               &flow->translation[1],
                               &flow->translation[2],
                               0);

    sprintf(keyword, "so%d:", group_solver_id[group_index]);
    err = find_keyword_string(filename,keyword,flow->solver_so_file,1);
}