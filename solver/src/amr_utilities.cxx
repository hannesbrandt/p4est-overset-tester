/**
 * \file    amr_utilities.c
 * \ingroup amr_group
 * \author  akirby
 *
 * \brief   Utility functions for the AMR code module.
 */

/* header files */
#include "amr_utilities.h"

#define eC "\x1B[0m"
#define GC "\x1B[1;32m"
#define gC "\x1B[0;32m"

void amr_utilities_inputs_amr_message(ctx_t *ctx,int mpi_rank){
    if (mpi_rank == 0 && ctx->log_info < DGLOG_ALL) {
        DGOUT(ctx->log_io,"\n");
        DGOUT(ctx->log_io,"+==========================================+\n");
        DGOUT(ctx->log_io," AMR Inputs:                                \n");
        DGOUT(ctx->log_io," ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾                        \n");

        DGOUT(ctx->log_io," unstructured_flag: %d\n",ctx->d_simulation.unstructured_flag);
        if (ctx->d_simulation.unstructured_flag) {
            DGOUT(ctx->log_io," unstructured_file: %s\n",ctx->d_gridfile.gridfile_name);
        }

        DGOUT(ctx->log_io," periodic_flag: %d %d %d\n",ctx->d_grid.periodic[0],ctx->d_grid.periodic[1],ctx->d_grid.periodic[2]);
        DGOUT(ctx->log_io," domain_lo: %f %f %f\n",ctx->d_grid.xlo[0],ctx->d_grid.xlo[1],ctx->d_grid.xlo[2]);
        DGOUT(ctx->log_io," domain_hi: %f %f %f\n",ctx->d_grid.xhi[0],ctx->d_grid.xhi[1],ctx->d_grid.xhi[2]);
        DGOUT(ctx->log_io," nelem: %d %d %d\n",ctx->d_grid.nelem[0],ctx->d_grid.nelem[1],ctx->d_grid.nelem[2]);
        DGOUT(ctx->log_io," max_amr_level: %d\n",ctx->d_grid.max_level);
        DGOUT(ctx->log_io," min_amr_level: %d\n",ctx->d_grid.min_level);
        DGOUT(ctx->log_io,"+==========================================+\n");
    }
}

void amr_utilities_grid_message(ctx_t *ctx){
    grid_t       *grd   = &ctx->d_grid;
    gridfile_t   *gfile = &ctx->d_gridfile;
    simulation_t *sim   = &ctx->d_simulation;

    if (ctx->rank == 0 && ctx->log_info < DGLOG_ALL) {
        DGOUT(ctx->log_io,"[ amr ]    Levels: %d\n",grd->max_level+1);
        DGOUT(ctx->log_io,"        Max Level: %d\n",grd->max_level);
        DGOUT(ctx->log_io,"        Min Level: %d\n",grd->min_level);

        if(!sim->unstructured_flag){
            DGOUT(ctx->log_io,"[ amr ]      Grid:      dx           dy           dz\n"
                   "         Level[%d]: %e %e %e\n"
                   "         Level[%d]: %e %e %e\n",
                    grd->max_level,grd->min_dx[0],grd->min_dx[1],grd->min_dx[2],
                                 0,grd->max_dx[0],grd->max_dx[1],grd->max_dx[2]);
#ifdef P4_TO_P8
            int nelem = grd->nelem[0]*grd->nelem[1]*grd->nelem[2];
            DGOUT(ctx->log_io,"[ amr ]  Level[0]: %d elements\n"
                   "               nx: %d\n"
                   "               ny: %d\n"
                   "               nz: %d\n"
                   ,nelem,grd->nelem[0],grd->nelem[1],grd->nelem[2]);
            DGOUT(ctx->log_io,"[ amr ]    Domain: \n"
                   "              xlo: %e   xhi: %e \n"
                   "              ylo: %e   yhi: %e \n"
                   "              zlo: %e   zhi: %e \n"
                   "           Volume: %.15e\n"
                   ,grd->xlo[0],grd->xhi[0]
                   ,grd->xlo[1],grd->xhi[1]
                   ,grd->xlo[2],grd->xhi[2]
                   ,grd->domain_volume);
#else
            int nelem = grd->nelem[0]*grd->nelem[1];
            DGOUT(ctx->log_io,"[ amr ]  Level[0]: %d elements\n"
                   "               nx: %d\n"
                   "               ny: %d\n"
                   ,nelem,grd->nelem[0],grd->nelem[1]);
            DGOUT(ctx->log_io,"[ amr ]    Domain: \n"
                   "              xlo: %e   xhi: %e \n"
                   "              ylo: %e   yhi: %e \n"
                   "           Volume: %.15f\n"
                   ,grd->xlo[0],grd->xhi[0]
                   ,grd->xlo[1],grd->xhi[1]
                   ,grd->domain_volume);
#endif
        } else {
            DGOUT(ctx->log_io,"[ amr ] Unstructured Grid file: %s\n",gfile->gridfile_name);
        }
    }
}

void amr_utilities_write_amr_regrid(ctx_t *ctx,int mpi_rank,Real max_time,
                                    Real balance_time,Real partition_time,
                                    Real ship_time,Real max_point_time,
                                    Real max_feature_time,Real max_spread_time,
                                    Real max_coarsen_time,Real max_buffer_time,
                                    Real max_level_time,Real pdegree_time,
                                    int log_info){
    if (mpi_rank == 0 && log_info < DGLOG_MESH_OFF) {
        DGOUT(ctx->log_io,"\n");
        DGOUT(ctx->log_io," Regrid Time (sec):  %1.2e\n",max_time);
        DGOUT(ctx->log_io," ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾\n");
        DGOUT(ctx->log_io,"  -Balance:          %1.2e\n",balance_time);
        DGOUT(ctx->log_io,"  -Partition:        %1.2e\n",partition_time);
        DGOUT(ctx->log_io,"  -Change P-Degree:  %1.2e\n",pdegree_time);
        DGOUT(ctx->log_io,"  -Custom Data Move: %1.2e\n",ship_time);
        DGOUT(ctx->log_io,"   --------------------------\n");
        DGOUT(ctx->log_io,"  -Regrid Points:    %1.2e\n",max_point_time);
        DGOUT(ctx->log_io,"  -Regrid Feature:   %1.2e\n",max_feature_time);
        DGOUT(ctx->log_io,"  -Regrid Spread:    %1.2e\n",max_spread_time);
        DGOUT(ctx->log_io,"  -Regrid Coarsen:   %1.2e\n",max_coarsen_time);
        DGOUT(ctx->log_io,"  -Regrid Buffer:    %1.2e\n",max_buffer_time);
        DGOUT(ctx->log_io,"  -Regrid Level:     %1.2e\n",max_level_time);
        DGOUT(ctx->log_io,"+====================================================+\n");
    }
}

Real amr_utilities_timer(){
    return MPI_Wtime();
}

Real amr_utilities_mpireducemax_real(int rank,MPI_Comm mpicomm,Real Real_in){
    Real max_Real;

    MPI_Reduce(&Real_in,&max_Real,1,MPI_DGREAL,MPI_MAX,0,mpicomm);
    return ((rank == 0) ? max_Real:0.0);
}

char *trimwhitespace(char *str){
    char *end;

    // Trim leading space
    while(isspace((unsigned char)*str)) str++;

    if (*str == 0) {  // All spaces?
        return str;
    }

    // Trim trailing space
    end = str + strlen(str) - 1;
    while(end > str && isspace((unsigned char)*end)) end--;

    // Write new null terminator character
    end[1] = '\0';
    return str;
}

int file_is_modified(const char *path,time_t oldMTime){
    struct stat file_stat;
    int err = stat(path, &file_stat);
    if (err != 0) {
        perror(" [file_is_modified] stat");
        exit(errno);
    }
    return file_stat.st_mtime > oldMTime;
}

void amr_utilities_create_directories(int mpi_rank,MPI_Comm mpicomm,int log_info){
    if (mpi_rank == 0) {
        mkdir("WRK", S_IRWXU | S_IRWXG | S_IRWXO);
    }
    MPI_Barrier(mpicomm);
}

char amr_utilities_find_keyword(char *filename,const char *keyword){
    char buff[BUFF_SIZE] = {'\0'};

    FILE *fp = fopen(filename,"r");
    while (fp != NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if(strstr(buff,keyword)) {fclose(fp); return EXIT_SUCCESS;}
    }

    fclose(fp);
    return EXIT_FAILURE;
}

char amr_utilities_find_keyword_integer(char *filename,const char *keyword,int *integer){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while (fp != NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {

            *integer = atoi(&buff[length]);

            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    DGOUT(stderr,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",keyword);

    fclose(fp);
    return EXIT_FAILURE;
}

char amr_utilities_find_keyword_two_integers(char *filename,const char *keyword,
                                             int *integer1,
                                             int *integer2){
    FILE *fp;
    char buff[BUFF_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while (fp != NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {

          sscanf(&buff[length],"%d %d", integer1,integer2);
          sscanf(&buff[length],"%d, %d",integer1,integer2);

          fclose(fp);
          return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    DGOUT(stderr,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",keyword);

    fclose(fp);
    return EXIT_FAILURE;
}

char amr_utilities_find_keyword_three_integers(char *filename,const char *keyword,
                                               int *integer1,
                                               int *integer2,
                                               int *integer3){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while (fp != NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {

            sscanf(&buff[length],"%d %d %d",  integer1,integer2,integer3);
            sscanf(&buff[length],"%d, %d, %d",integer1,integer2,integer3);
            sscanf(&buff[length],"%d, %d %d", integer1,integer2,integer3);
            sscanf(&buff[length],"%d %d, %d", integer1,integer2,integer3);

            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    DGOUT(stderr,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",keyword);

    fclose(fp);
    return EXIT_FAILURE;
}

char amr_utilities_find_keyword_real(char *filename,const char *keyword,Real *dbl){
    FILE *fp;
    char buff[BUFF_SIZE];
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while (fp != NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {

            sscanf(&buff[length],RealFormat,dbl);

            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    DGOUT(stderr,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",keyword);

    fclose(fp);
    return EXIT_FAILURE;
}

char amr_utilities_find_keyword_two_reals(char *filename,const char *keyword,
                                          Real *dbl1,
                                          Real *dbl2){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    char pref[BUFF_SIZE] = {'\0'};
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while (fp !=NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {
            sprintf(pref,"%s %s",RealFormat,RealFormat);
            sscanf(&buff[length],pref,dbl1,dbl2);

            sprintf(pref,"%s, %s",RealFormat,RealFormat);
            sscanf(&buff[length],pref,dbl1,dbl2);

            sprintf(pref,"%s,%s",RealFormat,RealFormat);
            sscanf(&buff[length],pref,dbl1,dbl2);

            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    DGOUT(stderr,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",keyword);

    fclose(fp);
    return EXIT_FAILURE;
}

char amr_utilities_find_keyword_three_reals(char *filename,const char *keyword,
                                            Real *dbl1,
                                            Real *dbl2,
                                            Real *dbl3){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    char pref[BUFF_SIZE] = {'\0'};
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while (fp !=NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {
            sprintf(pref,"%s %s %s",RealFormat,RealFormat,RealFormat);
            sscanf(&buff[length],pref,  dbl1,dbl2,dbl3);

            sprintf(pref,"%s, %s, %s",RealFormat,RealFormat,RealFormat);
            sscanf(&buff[length],pref,dbl1,dbl2,dbl3);

            sprintf(pref,"%s %s, %s",RealFormat,RealFormat,RealFormat);
            sscanf(&buff[length],pref, dbl1,dbl2,dbl3);

            sprintf(pref,"%s, %s %s",RealFormat,RealFormat,RealFormat);
            sscanf(&buff[length],pref, dbl1,dbl2,dbl3);

            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    DGOUT(stderr,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",keyword);

    fclose(fp);
    return EXIT_FAILURE;
}

char amr_utilities_find_keyword_int_real(char *filename,const char *keyword,
                                         int *int1,
                                         Real *dbl1){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    char pref[BUFF_SIZE] = {'\0'};
    fp = fopen(filename,"r");

    unsigned long length = strlen(keyword);

    while (fp !=NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {
            sprintf(pref,"%s %s",IntFormat,RealFormat);
            sscanf(&buff[length],pref,int1,dbl1);

            sprintf(pref,"%s, %s",IntFormat,RealFormat);
            sscanf(&buff[length],pref,int1,dbl1);

            sprintf(pref,"%s ,%s",IntFormat,RealFormat);
            sscanf(&buff[length],pref,int1,dbl1);

            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    DGOUT(stderr,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",keyword);

    fclose(fp);
    return EXIT_FAILURE;
}

char amr_utilities_find_keyword_string(char *filename,const char *keyword,
                                       char *string){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    unsigned long i,length;

    fp = fopen(filename,"r");
    length = strlen(keyword);

    while(fp != NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {

            int ind = 0;
            for (i = length; i < BUFF_SIZE; ++i) {
                if (buff[i] != ' ') {
                    string[ind] = buff[i];
                    ind = ind+1;
                }
            }

            /* cutoff comments after # sign */
            char *ptr = strchr(string, '#');
            if(ptr != NULL) *ptr = '\0';

            string[strcspn(string, "\n")] = 0;

            if (strlen(string) >= BUFF_SIZE) {
                DGOUT(stderr,"\033[1;31m>>>>> string length is greater than allowable.\033[0m");
                exit(EXIT_FAILURE);
            }

            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    DGOUT(stderr,"\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",keyword);

    fclose(fp);
    return EXIT_FAILURE;
}

char amr_utilities_find_keyword_string2upper(char *filename,const char *keyword,
                                             char *string){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    unsigned long i,length;

    fp = fopen(filename,"r");
    length = strlen(keyword);

    while(fp != NULL && fgets(buff,sizeof(buff),fp) != NULL) {
        if (strstr(buff,keyword)) {

            int ind=0;
            for (i = length; i < BUFF_SIZE; ++i) {
                if (buff[i] != ' ') {
                    string[ind] = toupper(buff[i]);
                    ind = ind+1;
                }
            }

            /* cutoff comments after # sign */
            char *ptr = strchr(string, '#');
            if(ptr != NULL) *ptr = '\0';

            string[strcspn(string, "\n")] = 0;

            if (strlen(string) >= BUFF_SIZE) {
                printf("\033[1;31m>>>>> string length is greater than allowable.\033[0m\n");
                exit(EXIT_FAILURE);
            }

            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    /* otherwise */
    // printf("\033[1;31m>>>>> Keyword [%s] not found! Using default.\033[0m\n",
    //         keyword);

    fclose(fp);
    return EXIT_FAILURE;
}

char amr_utilities_find_keyword_string_optional(char *filename,int upper_flag,
                                                const char *keyword,char *string){
    FILE *fp;
    char buff[BUFF_SIZE] = {'\0'};
    unsigned long i,length;

    fp = fopen(filename,"r");
    length = strlen(keyword);

    while(fp != NULL && fgets(buff,sizeof(buff),fp) != NULL) {

        int keyword_len_match = (buff[length-1] == ':');
        if (strstr(buff,keyword) && keyword_len_match) {

            int ind = 0;
            for (i = length; i < BUFF_SIZE; ++i) {
                // if (buff[i] != ' ') {
                if (upper_flag){
                    string[ind] = toupper(buff[i]);
                } else {
                    string[ind] = buff[i];
                }
                ind = ind+1;
                // }
            }

            /* cutoff comments after # sign */
            char *ptr = strchr(string, '#');
            if(ptr != NULL) *ptr = '\0';

            string[strcspn(string, "\n")] = 0;

            if (strlen(string) >= BUFF_SIZE) {
                DGOUT(stderr,"\033[1;31m>>>>> string length is greater than allowable\033[0m\n");
                exit(EXIT_FAILURE);
            }
            fclose(fp);
            return EXIT_SUCCESS;
        }
    }
    fclose(fp);
    return EXIT_FAILURE;
}