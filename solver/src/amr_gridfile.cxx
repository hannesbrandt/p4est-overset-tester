/**
 * \file    amr_gridfile.cxx
 * \ingroup amr_group
 * \author  mugolott
 *
 * \brief Grid-file IO
 */

/* header files */
#include "amr_gridfile.h"

/*
 * Taken from p4est_connectivity.c, p4est_connectivity_getline_upper
 * Read a line from a file. Obtained from:
 * http://stackoverflow.com/questions/314401/
 * how-to-read-a-line-from-the-console-in-c/314422#314422
 *
 * Using this avoids a dependence on IEEE Std 1003.1-2008 (``POSIX.1'') for the
 * getline function.
 */
char *get_line(FILE *stream,int ln){
    char *line = (char*) malloc(1024*sizeof(char)), *linep = line;
    uint64_t lenmax = 1024, len = lenmax;
    int c;

    for (;;) {
        /* Read character and make it upper case */
        c = fgetc(stream);
        c = toupper(c);

        /* if end of file return NULL */
        if (c == EOF && linep == line) {
            free(linep);
            return NULL;
        }

        /* if filled allocation, reallocate */
        if (--len == 0) {
            char *linen;

            len = lenmax;
            lenmax *= 2;

            linen = (char *) realloc(linep,lenmax);
            if (linen == NULL) {
                DGOUT(stderr,"\x1B[1;31m[GMSH] Error in reading line %d!\x1B[0m\n",ln);
                free(linep);
                return NULL;
            }

            line = linen + (line - linep);
            linep = linen;
        }

        /* if c is a special character do not save it */
        if(c != '"' && c != '\'') *line++ = c;

        /* if end of line, break */
        if(c == '\n') break;
    }
    *line = '\0';
    return linep;
}

/*
 * From StringUtils.cpp in gmsh:
 * https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/Common/StringUtils.cpp
 *
 */
static void SwapBytes(char *array,int size,int n){
    char x[size];
    int i, c;

    for (i = 0; i < n; i++) {
        char *a = &array[i * size];
        memcpy(x, a, size);

        for(c = 0; c < size; c++) a[size - 1 - c] = x[c];
    }
}

static int read_msh4_entity(FILE *stream,int binary, int swap,int version,
                            int dim,int *lines_read,int *lines_free,
                            int grid_dim,size_t ient,geometry_t *geom){
    Real minX,maxX,minY,maxY,minZ,maxZ;
    char *line;
    int tag;

    if (binary) {
        int dataInt;
        if(fread(&dataInt, sizeof(int), 1, stream) != 1) { return 1; }
        if(swap == 1) SwapBytes((char *)&dataInt, sizeof(int), 1);

        Real dataDouble[6];
        const size_t nbb = (dim > 0) ? 6 : 3;

        if(fread(dataDouble, sizeof(Real), nbb, stream) != nbb) { return 1; }
        if(swap == 1) SwapBytes((char *)dataDouble, sizeof(Real), nbb);

        tag = dataInt;
        minX = dataDouble[0];
        minY = dataDouble[1];
        minZ = dataDouble[2];
        maxX = dataDouble[(nbb == 6) ? 3 : 0];
        maxY = dataDouble[(nbb == 6) ? 4 : 1];
        maxZ = dataDouble[(nbb == 6) ? 5 : 2];
    } else {
        /* Read next line */
        line = get_line(stream,++*lines_read);

        if (dim > 0) {
            if (sscanf(line,"%d %lf %lf %lf %lf %lf %lf %[^\n]",&tag,
                            &minX,&minY,&minZ,&maxX,&maxY,&maxZ,line) != 8) {
                return 1;
            }
        } else {
            if (sscanf(line,"%d %lf %lf %lf %[^\n]",&tag,
                            &minX, &minY, &minZ, line) != 5) {
                return 1;
            }
            maxX = minX;
            maxY = minY;
            maxZ = minZ;
        }
    }

    /* store data */
    if (dim == 1 && grid_dim == 2) {
        geom->entity_face_tag[ient]      = tag-1;
        geom->entity_face_patchtag[ient] = NO_PATCH;
    } else
    if (dim == 2 && grid_dim == 2) {
        geom->entity_vol_tag[ient]      = tag-1;
        geom->entity_vol_patchtag[ient] = NO_PATCH;
    } else
    if (dim == 2 && grid_dim == 3) {
        geom->entity_face_tag[ient]      = tag-1;
        geom->entity_face_patchtag[ient] = NO_PATCH;
    } else
    if (dim == 3 && grid_dim == 3) {
        geom->entity_vol_tag[ient]      = tag-1;
        geom->entity_vol_patchtag[ient] = NO_PATCH;
    }

    /* read physical tags */
    size_t numPhy = 0;
    if (binary) {
        if(fread(&numPhy, sizeof(size_t), 1, stream) != 1) { return 1; }
        if(swap == 1) SwapBytes((char *)&numPhy, sizeof(size_t), 1);

        int phyTags[numPhy];
        if(fread(&phyTags[0], sizeof(int), numPhy, stream) != numPhy) { return 1; }
        if(swap == 1) SwapBytes((char *)&phyTags[0], sizeof(int), numPhy);

        /* NOTE: we assume only one tag per entity */
        if (numPhy != 0 && dim == 1 && grid_dim == 2) {
            geom->entity_face_patchtag[ient] = phyTags[0]-1;
        } else
        if (numPhy != 0 && dim == 2 && grid_dim == 2) {
            geom->entity_vol_patchtag[ient] = phyTags[0]-1;
        } else
        if (numPhy != 0 && dim == 2 && grid_dim == 3) {
            geom->entity_face_patchtag[ient] = phyTags[0]-1;
        } else
        if (numPhy != 0 && dim == 3 && grid_dim == 3) {
            geom->entity_vol_patchtag[ient] = phyTags[0]-1;
        }

    } else {
        int phyTags[0];
        sscanf(line, "%lu %d", &numPhy, &phyTags[0]);

        if (numPhy != 0 && dim == 1 && grid_dim == 2) {
            geom->entity_face_patchtag[ient] = phyTags[0]-1;
        } else
        if (numPhy != 0 && dim == 2 && grid_dim == 2) {
            geom->entity_vol_patchtag[ient] = phyTags[0]-1;
        } else
        if (numPhy != 0 && dim == 2 && grid_dim == 3) {
            geom->entity_face_patchtag[ient] = phyTags[0]-1;
        } else
        if (numPhy != 0 && dim == 3 && grid_dim == 3) {
            geom->entity_vol_patchtag[ient] = phyTags[0]-1;
        }
    }

    /* read remaining data on the line: bounding points (not necessary for dg4est) */
    size_t numBrep = 0;
    if (dim > 0 && binary) {
        if(fread(&numBrep, sizeof(size_t), 1, stream) != 1) { return 1; }
        if(swap == 1) SwapBytes((char *)&numBrep, sizeof(size_t), 1);

        int brepTags[numPhy];
        if(fread(&brepTags[0], sizeof(int), numBrep, stream) != numBrep) { return 1; }
        if(swap == 1) SwapBytes((char *)&brepTags[0], sizeof(int), numBrep);
    }

    if (!binary) {
        /* Free line */
        free(line); ++lines_free;
    }
    return 0;
}

static void msh_supported_elem_type(int elem_type,int sim_dim,
                                    char *supported_vol_elem_type,
                                    char *supported_face_elem_type,
                                    int *nvert){
    char supported_1d_elem,supported_2d_elem,supported_3d_elem;

    /* For element type ID, see:
     * https://acdl.mit.edu/csi/userContent/CoverageHTML/Meshing/gmsh/GMSHtoSANS.cpp.gcov.frameset.html
     * https://gmsh.info/doc/texinfo/gmsh.html#High_002dorder-elements
     */
    supported_1d_elem = (elem_type == 1  || // Curve Q=1  (2 nodes)
                         elem_type == 8  || // Curve Q=2  (3 nodes)
                         elem_type == 26 || // Curve Q=3  (4 nodes)
                         elem_type == 27 || // Curve Q=4  (5 nodes)
                         elem_type == 28 || // Curve Q=5  (6 nodes)
                         elem_type == 62 || // Curve Q=6  (7 nodes)
                         elem_type == 63 || // Curve Q=7  (8 nodes)
                         elem_type == 64 || // Curve Q=8  (9 nodes)
                         elem_type == 65 || // Curve Q=9  (10 nodes)
                         elem_type == 66    // Curve Q=10 (11 nodes)
                         );
    supported_2d_elem = (elem_type == 3  || // Quadrangle Q=1  (4 nodes)
                         elem_type == 10 || // Quadrangle Q=2  (9 nodes)
                         elem_type == 36 || // Quadrangle Q=3  (16 nodes)
                         elem_type == 37 || // Quadrangle Q=4  (25 nodes)
                         elem_type == 38 || // Quadrangle Q=5  (36 nodes)
                         elem_type == 47 || // Quadrangle Q=6  (49 nodes)
                         elem_type == 48 || // Quadrangle Q=7  (64 nodes)
                         elem_type == 49 || // Quadrangle Q=8  (81 nodes)
                         elem_type == 50 || // Quadrangle Q=9  (100 nodes)
                         elem_type == 51    // Quadrangle Q=10 (121 nodes)
                         );
    supported_3d_elem = (elem_type == 5  || // Hexahedron Q=1  (8 nodes)
                         elem_type == 12 || // Hexahedron Q=2  (27 nodes)
                         elem_type == 92 || // Hexahedron Q=3  (64 nodes)
                         elem_type == 93 || // Hexahedron Q=4  (125 nodes)
                         elem_type == 94 || // Hexahedron Q=5  (216 nodes)
                         elem_type == 95 || // Hexahedron Q=6  (343 nodes)
                         elem_type == 96 || // Hexahedron Q=7  (512 nodes)
                         elem_type == 97 || // Hexahedron Q=8  (729 nodes)
                         elem_type == 98 || // Hexahedron Q=9  (1000 nodes)
                         elem_type == 99    // Hexahedron Q=10 (1331 nodes) NB: not supported by gmsh!
                         );

    if (sim_dim == 2) {
        *supported_vol_elem_type  = supported_2d_elem;
        *supported_face_elem_type = supported_1d_elem;
    } else
    if (sim_dim == 3) {
        *supported_vol_elem_type  = supported_3d_elem;
        *supported_face_elem_type = supported_2d_elem;
    }

    switch (elem_type) {
        case  1: *nvert =    2; break;
        case  3: *nvert =    4; break;
        case  5: *nvert =    8; break;
        case  8: *nvert =    3; break;
        case 10: *nvert =    9; break;
        case 12: *nvert =   27; break;
        case 15: *nvert =    1; break;
        case 26: *nvert =    4; break;
        case 27: *nvert =    5; break;
        case 28: *nvert =    6; break;
        case 36: *nvert =   16; break;
        case 37: *nvert =   25; break;
        case 38: *nvert =   36; break;
        case 47: *nvert =   49; break;
        case 48: *nvert =   64; break;
        case 49: *nvert =   81; break;
        case 50: *nvert =  100; break;
        case 51: *nvert =  121; break;
        case 62: *nvert =    7; break;
        case 63: *nvert =    8; break;
        case 64: *nvert =    9; break;
        case 65: *nvert =   10; break;
        case 66: *nvert =   11; break;
        case 92: *nvert =   64; break;
        case 93: *nvert =  125; break;
        case 94: *nvert =  216; break;
        case 95: *nvert =  343; break;
        case 96: *nvert =  512; break;
        case 97: *nvert =  729; break;
        case 98: *nvert = 1000; break;
        case 99: *nvert = 1331; break;
        default:
            DGOUT(stderr,"\x1B[1;31m[GMSH] Unsupported element found while reading grid file!\x1B[0m\n");
            exit(1);
    }
}

static void msh_elem_type_to_qdegree(int elem_type,int *qdegree){
    switch (elem_type) {
        case  3: *qdegree = 1; break; // (2D: 4 nodes)
        case  5: *qdegree = 1; break; // (3D: 8 nodes)
        case 10: *qdegree = 2; break; // (2D: 9 nodes)
        case 12: *qdegree = 2; break; // (3D: 27 nodes)
        case 36: *qdegree = 3; break; // (2D: 16 nodes)
        case 92: *qdegree = 3; break; // (3D: 64 nodes)
        case 37: *qdegree = 4; break; // (2D: 25 nodes)
        case 93: *qdegree = 4; break; // (3D: 125 nodes)
        case 38: *qdegree = 5; break; // (2D: 36 nodes)
        case 94: *qdegree = 5; break; // (3D: 216 nodes)
        case 47: *qdegree = 6; break; // (2D: 49 nodes)
        case 95: *qdegree = 6; break; // (3D: 343 nodes)
        case 48: *qdegree = 7; break; // (2D: 64 nodes)
        case 96: *qdegree = 7; break; // (3D: 512 nodes)
        case 49: *qdegree = 8; break; // (2D: 81 nodes)
        case 97: *qdegree = 8; break; // (3D: 729 nodes)
        case 50: *qdegree = 9; break; // (2D: 100 nodes)
        case 98: *qdegree = 9; break; // (3D: 1000 nodes)
        case 51: *qdegree =10; break; // (2D: 121 nodes)
        case 99: *qdegree =10; break; // (3D: 1331 nodes)
        default:
            DGOUT(stderr,"\x1B[1;31m[GMSH] Unsupported element found while reading grid file!\x1B[0m\n");
            exit(1);
    }
}

static void msh_map_node_numbering_flex(int grid_dim,int qorder,Uint *map){
    int inode = 0;
    int fnode = 0;
    int iedge;
    int iface;
    int iloop;
    int floop;
    int iedge_pt;
    int qdeg = qorder;
    int nnodes = DGVDIM(qdeg+1);
    int edge_nnodes = qdeg+1;

    /* 2D */
    if (grid_dim == 2) {
        for (iloop = 1; inode < nnodes; iloop++) {
            /* vertices */
            map[inode++] = (edge_nnodes)*(iloop-1) + iloop - 1;
            map[inode++] = (edge_nnodes)*iloop - (iloop-1) - 1;
            map[inode++] = (edge_nnodes)*(edge_nnodes-(iloop-1)) - (iloop-1) - 1;
            map[inode++] = (edge_nnodes)*(qdeg-(iloop-1)) + iloop - 1;

            /* edges */
            int edge_pts = (qdeg-1) - 2*(iloop-1);
            for (iedge = 1; iedge < 5; iedge++) {
                if (iedge == 1) {
                    /* edge 1 */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes)*(iloop-1)+ iloop + iedge_pt - 1;
                    }
                } else
                if (iedge == 2) {
                    /* edge 2 */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*iloop) + (edge_nnodes*iedge_pt) - (iloop-1) - 1;
                    }
                } else
                if (iedge == 3) {
                    /* edge 3 */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*(edge_nnodes-(iloop-1))) - (iloop-1) - iedge_pt - 1;
                    }
                } else
                if (iedge == 4) {
                    /* edge 4 */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*(edge_nnodes-1-iedge_pt-(iloop-1))) + iloop - 1;
                    }
                }
            }

            /* if have only one center node left (Q is a even number) */
            if (inode == nnodes-1) {
                iloop++;
                map[inode++] = (edge_nnodes)*(iloop-1) + iloop - 1;
            }
            /* if have only 4 vertices left (Q is a odd number) */
            if (inode == nnodes-4) {
                iloop++;
                map[inode++] = (edge_nnodes)*(iloop-1) + iloop - 1;
                map[inode++] = (edge_nnodes)*iloop - (iloop-1) - 1;
                map[inode++] = (edge_nnodes)*(edge_nnodes-(iloop-1)) - (iloop-1) - 1;
                map[inode++] = (edge_nnodes)*(qdeg-(iloop-1)) + iloop - 1;
            }
        }
    }

    /* 3D */
    if (grid_dim == 3) {
        for (iloop = 1; inode < nnodes; iloop++) {
            /* vertices */
            map[inode++] = (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes)*(iloop-1) + iloop - 1;
            map[inode++] = (edge_nnodes*edge_nnodes*(iloop-1)) +(edge_nnodes)*iloop - (iloop-1) - 1;
            map[inode++] = (edge_nnodes*edge_nnodes*(iloop-1)) +(edge_nnodes)*(edge_nnodes-(iloop-1)) - (iloop-1) - 1;
            map[inode++] = (edge_nnodes*edge_nnodes*(iloop-1)) +(edge_nnodes)*(qdeg-(iloop-1)) + iloop - 1;
            map[inode++] = (edge_nnodes*edge_nnodes*qdeg) - (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes)*(iloop-1) + iloop - 1;
            map[inode++] = (edge_nnodes*edge_nnodes*qdeg) - (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes)*iloop - (iloop-1) - 1;
            map[inode++] = (edge_nnodes*edge_nnodes*qdeg) - (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes)*(edge_nnodes-(iloop-1)) - (iloop-1) - 1;
            map[inode++] = (edge_nnodes*edge_nnodes*qdeg) - (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes)*(qdeg-(iloop-1)) +  iloop - 1;

            /* edges */
            for (iedge = 1; iedge < 13; iedge++) {
                int edge_pts = (qdeg-1)-2*(iloop-1);

                if (iedge == 1) {
                    /* edge 1 - 0->1 (x direction) */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes)*(iloop-1)+ iloop + iedge_pt - 1;
                    }
                } else
                if (iedge == 2) {
                    /* edge 2 - 0->3 (y direction) */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes*(iloop-1)) + (edge_nnodes*iedge_pt) + iloop - 1;
                    }
                } else
                if (iedge == 3) {
                    /* edge 3 - 0->4 (z direction) */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes*edge_nnodes*(iedge_pt)) + (edge_nnodes*(iloop-1)) + iloop -  1;
                    }
                } else
                if (iedge == 4) {
                    /* edge 4 - 1->2 (y direction) */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes*(iedge_pt+iloop)) - (iloop-1) - 1;
                    }
                } else
                if (iedge == 5) {
                    /* edge 5 - 1->5 (z direction) */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*edge_nnodes*((iloop-1)+(iedge_pt))) + edge_nnodes*(iloop-1) + edge_nnodes - (iloop-1) - 1;
                    }
                } else
                if (iedge == 6) {
                    /* edge 6 - 2->3 (x direction) */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes*(edge_nnodes-(iloop-1))) - (iloop-1) - iedge_pt - 1;
                    }
                } else
                if (iedge == 7) {
                    /* edge 7 - 2->6 (z direction) */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*edge_nnodes*((iloop-1) +iedge_pt)) + (edge_nnodes*(edge_nnodes-(iloop-1))) - (iloop-1) -1;
                    }
                } else
                if (iedge == 8) {
                    /* edge 8 - 3->7 (z direction) */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*edge_nnodes*((iloop-1) +iedge_pt)) + (edge_nnodes*(edge_nnodes-iloop)) + iloop -1;
                    }
                } else
                if (iedge == 9) {
                    /* edge 9 - 4->5 (x direction) */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*edge_nnodes*(edge_nnodes-iloop)) + (edge_nnodes*(iloop-1)) + iloop + iedge_pt - 1;
                    }
                } else
                if (iedge == 10) {
                    /* edge 10 - 4->7 (y direction) */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*edge_nnodes*(edge_nnodes-iloop)) + (edge_nnodes*(iedge_pt+(iloop-1))) + iloop - 1;
                    }
                } else
                if (iedge == 11) {
                    /* edge 11 - 5->6 (y direction) */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*edge_nnodes*(edge_nnodes-iloop)) + (edge_nnodes*(iedge_pt+iloop)) - (iloop-1) - 1;
                    }
                } else
                if (iedge == 12) {
                    /* edge 12 - 6->7 (x direction) */
                    for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                        map[inode++] = (edge_nnodes*edge_nnodes*(edge_nnodes-iloop)) + (edge_nnodes*(edge_nnodes-(iloop-1))) - (iloop-1) - iedge_pt - 1;
                    }
                }
            }

            /* surfaces */
            for (iface = 1; iface < 7; iface++) {
                int face_pts = ((qdeg-1)-2*(iloop-1))*((qdeg-1)-2*(iloop-1));

                if (iface == 1) {
                    /* face 1 z-min */
                    fnode = 0;
                    for (floop = 1; fnode < face_pts; floop++) {

                        /* if have only one center node left (Q is a even number) */
                        if (fnode == face_pts-1) {
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop-1) + (edge_nnodes)*(iloop) + iloop + (edge_nnodes)*(floop-1) + floop - 1;
                            break;
                        }
                        /* if have only 4 vertices left (Q is a odd number) */
                        if (fnode == face_pts-4) {
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop-1) + (edge_nnodes)*(iloop) + iloop + (edge_nnodes)*(floop-1) + floop - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop-1) + (edge_nnodes)*(qdeg-iloop-(floop-1)) + iloop + floop - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop-1) + (edge_nnodes)*(edge_nnodes-iloop-(floop-1)) - iloop - (floop-1) - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop-1) + (edge_nnodes)*(iloop) - iloop + (edge_nnodes)*floop - (floop-1) - 1;
                            break;
                        }

                        /* vertices */
                        map[inode++] = (edge_nnodes*edge_nnodes)*(iloop-1) + (edge_nnodes)*(iloop) + iloop + (edge_nnodes)*(floop-1) + floop - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(iloop-1) + (edge_nnodes)*(qdeg-iloop-(floop-1)) + iloop + floop - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(iloop-1) + (edge_nnodes)*(edge_nnodes-iloop-(floop-1)) - iloop - (floop-1) - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(iloop-1) + (edge_nnodes)*(iloop) - iloop + (edge_nnodes)*floop - (floop-1) - 1;

                        fnode = fnode + 4;

                        /* edges */
                        for (iedge = 1; iedge < 5; iedge++) {
                            int edge_pts = (qdeg-1)-2*(floop-1)-2*iloop;

                            if (iedge == 1) {
                                /* edge 1 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(iloop-1) + (edge_nnodes)*(iloop) + (edge_nnodes)*(floop-1+iedge_pt) + iloop + floop - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 2) {
                                /* edge 2 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(iloop-1) + (edge_nnodes)*(qdeg-iloop-(floop-1)) + iloop + floop + iedge_pt - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 3) {
                                /* edge 3 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(iloop-1) + (edge_nnodes)*(edge_nnodes-iloop-(floop-1)-iedge_pt) - iloop - (floop-1) - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 4) {
                                /* edge 4 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(iloop-1) + (edge_nnodes)*(iloop) - iloop + (edge_nnodes)*floop - (floop-1) - iedge_pt - 1;
                                    fnode ++;
                                }
                            }
                        }
                    }
                } else
                if (iface == 2) {
                    /* face 2 y-min */
                    fnode = 0;
                    for (floop = 1; fnode < face_pts; floop++) {

                        /* if have only one center node left (Q is a even number) */
                        if (fnode == face_pts-1) {
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(iloop-1) + iloop + floop - 1;
                            break;
                        }
                        /* if have only 4 vertices left (Q is a odd number) */
                        if (fnode == face_pts-4) {
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(iloop-1) + iloop + floop - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(iloop-1) + edge_nnodes - iloop - (floop-1) - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + edge_nnodes*iloop - iloop - (floop-1) - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(iloop-1) + iloop + floop - 1;
                            break;
                        }

                        /* vertices */
                        map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(iloop-1) + iloop + floop - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(iloop-1) + edge_nnodes - iloop - (floop-1) - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + edge_nnodes*iloop - iloop - (floop-1) - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(iloop-1) + iloop + floop - 1;
                        fnode = fnode + 4;

                        /* edges */
                        for (iedge = 1; iedge < 5; iedge++) {
                            int edge_pts = (qdeg-1)-2*(floop-1)-2*iloop;

                            if (iedge == 1) {
                                /* edge 1 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(iloop-1) + iloop + iedge_pt + floop - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 2) {
                                /* edge 2 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)+iedge_pt) + (edge_nnodes)*iloop - (iloop-1) - floop - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 3) {
                                /* edge 3 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes*iloop) - (floop-1) - iloop - iedge_pt - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 4) {
                                /* edge 4 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop-iedge_pt) + edge_nnodes*(iloop-1) + iloop + floop - 1;
                                    fnode ++;
                                }
                            }
                        }
                    }
                } else
                if (iface == 3) {
                    /* face 3 x-min */
                    fnode = 0;
                    for (floop = 1; fnode < face_pts; floop++) {

                        /* if have only one center node left (Q is a even number) */
                        if (fnode == face_pts-1) {
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(iloop+(floop-1)) + iloop - 1;
                            break;
                        }
                        /* if have only 4 vertices left (Q is a odd number) */
                        if (fnode == face_pts-4) {
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(iloop+(floop-1)) + iloop - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(iloop+(floop-1)) + iloop - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(edge_nnodes-1-iloop-(floop-1)) + iloop - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(edge_nnodes-1-iloop-(floop-1)) + iloop - 1;
                            break;
                        }

                        /* vertices */
                        map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(iloop+(floop-1)) + iloop - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(iloop+(floop-1)) + iloop - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(edge_nnodes-1-iloop-(floop-1)) + iloop - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(edge_nnodes-1-iloop-(floop-1)) + iloop - 1;
                        fnode = fnode + 4;

                        /* edges */
                        for (iedge = 1; iedge < 5; iedge++) {
                            int edge_pts = (qdeg-1)-2*(floop-1)-2*iloop;

                            if (iedge == 1) {
                                /* edge 1 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)+iedge_pt) + (edge_nnodes)*(iloop+(floop-1)) + iloop - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 2) {
                                /* edge 2 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(iloop+(floop-1)+iedge_pt) + iloop - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 3) {
                                /* edge 3 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop-iedge_pt) + (edge_nnodes)*(edge_nnodes-1-iloop-(floop-1)) + iloop - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 4) {
                                /* edge 4 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(edge_nnodes-1-iloop-(floop-1)-iedge_pt) + iloop - 1;
                                    fnode ++;
                                }
                            }
                        }
                    }
                } else
                if (iface == 4) {
                    /* face 4 x-max */
                    fnode = 0;
                    for (floop = 1; fnode < face_pts; floop++) {

                        /* if have only one center node left (Q is a even number) */
                        if (fnode == face_pts-1) {
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(iloop+(floop-1)+1) - (iloop-1) - 1;
                            break;
                        }
                        /* if have only 4 vertices left (Q is a odd number) */
                        if (fnode == face_pts-4) {
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(iloop+(floop-1)+1) - (iloop-1) - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(edge_nnodes-iloop-(floop-1)) - (iloop-1) - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(edge_nnodes-iloop-(floop-1)) - (iloop-1) - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(iloop+(floop-1)+1) - (iloop-1) - 1;
                            break;
                        }

                        /* vertices */
                        map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(iloop+(floop-1)+1) - (iloop-1) - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(edge_nnodes-iloop-(floop-1)) - (iloop-1) - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(edge_nnodes-iloop-(floop-1)) - (iloop-1) - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(iloop+(floop-1)+1) - (iloop-1) - 1;
                        fnode = fnode + 4;

                        /* edges */
                        for (iedge = 1; iedge < 5; iedge++) {
                            int edge_pts = (qdeg-1)-2*(floop-1)-2*iloop;

                            if (iedge == 1) {
                                /* edge 1 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(iloop+(floop-1)+iedge_pt+1) - (iloop-1) - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 2) {
                                /* edge 2 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)+iedge_pt) + (edge_nnodes)*(edge_nnodes-iloop-(floop-1)) - (iloop-1) - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 3) {
                                /* edge 3 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(edge_nnodes-iloop-(floop-1)-iedge_pt) - (iloop-1) - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 4) {
                                /* edge 4 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop-iedge_pt) + (edge_nnodes)*(iloop+(floop-1)+1) - (iloop-1) - 1;
                                    fnode ++;
                                }
                            }
                        }
                    }
                } else
                if (iface == 5) {
                    /* face 5 y-max */
                    fnode = 0;
                    for (floop = 1; fnode < face_pts; floop++) {

                        /* if have only one center node left (Q is a even number) */
                        if (fnode == face_pts-1) {
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(edge_nnodes -(iloop-1)) - iloop - (floop-1) - 1;
                            break;
                        }
                        /* if have only 4 vertices left (Q is a odd number) */
                        if (fnode == face_pts-4) {
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(edge_nnodes -(iloop-1)) - iloop - (floop-1) - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(edge_nnodes - iloop) + iloop + floop - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(edge_nnodes - iloop) + iloop + floop - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(edge_nnodes -(iloop-1)) - iloop - (floop-1) - 1;
                            break;
                        }

                        /* vertices */
                        map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(edge_nnodes -(iloop-1)) - iloop - (floop-1) - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(edge_nnodes - iloop) + iloop + floop - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(edge_nnodes - iloop) + iloop + floop - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(edge_nnodes -(iloop-1)) - iloop - (floop-1) - 1;
                        fnode = fnode + 4;

                        /* edges */
                        for (iedge = 1; iedge < 5; iedge++) {
                            int edge_pts = (qdeg-1)-2*(floop-1)-2*iloop;

                            if (iedge == 1) {
                                /* edge 1 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)) + (edge_nnodes)*(edge_nnodes -(iloop-1)) - iloop - (floop-1) - iedge_pt - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 2) {
                                /* edge 2 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(iloop+(floop-1)+iedge_pt) + (edge_nnodes)*(edge_nnodes - iloop) + iloop + floop - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 3) {
                                /* edge 3 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop) + (edge_nnodes)*(edge_nnodes - iloop) + iloop + floop + iedge_pt - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 4) {
                                /* edge 4 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop-floop-iedge_pt) + (edge_nnodes)*(edge_nnodes -(iloop-1)) - iloop - (floop-1) - 1;
                                    fnode ++;
                                }
                            }
                        }
                    }
                } else
                if (iface == 6) {
                    /* face 6 z-max */
                    fnode = 0;
                    for (floop = 1; fnode < face_pts; floop++) {

                        /* if have only one center node left (Q is a even number) */
                        if (fnode == face_pts-1) {
                            map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop) + (edge_nnodes)*(iloop+(floop-1)) + iloop + floop - 1;
                            break;
                        }
                        /* if have only 4 vertices left (Q is a odd number) */
                        if (fnode == face_pts-4) {
                            map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop) + (edge_nnodes)*(iloop+(floop-1)) + iloop + floop - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop) + (edge_nnodes)*(iloop+floop) - iloop - (floop-1) - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop) + (edge_nnodes)*(edge_nnodes-iloop-(floop-1)) - iloop - (floop-1) - 1;
                            map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop) + (edge_nnodes)*(qdeg-iloop-(floop-1)) + iloop + floop - 1;
                            break;
                        }

                        /* vertices */
                        map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop) + (edge_nnodes)*(iloop+(floop-1)) + iloop + floop - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop) + (edge_nnodes)*(iloop+floop) - iloop - (floop-1) - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop) + (edge_nnodes)*(edge_nnodes-iloop-(floop-1)) - iloop - (floop-1) - 1;
                        map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop) + (edge_nnodes)*(qdeg-iloop-(floop-1)) + iloop + floop - 1;
                        fnode = fnode + 4;

                        /* edges */
                        for (iedge = 1; iedge < 5; iedge++) {
                            int edge_pts = (qdeg-1)-2*(floop-1)-2*iloop;

                            if (iedge == 1) {
                                /* edge 1 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop) + (edge_nnodes)*(iloop+(floop-1)) + iloop + floop + iedge_pt - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 2) {
                                /* edge 2 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop) + (edge_nnodes)*(iloop+floop+iedge_pt) - iloop - (floop-1) - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 3) {
                                /* edge 3 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop) + (edge_nnodes)*(edge_nnodes-iloop-(floop-1)) - iloop - (floop-1) - iedge_pt - 1;
                                    fnode ++;
                                }
                            } else
                            if (iedge == 4) {
                                /* edge 4 */
                                for (iedge_pt = 1; iedge_pt < edge_pts+1; iedge_pt++) {
                                    map[inode++] = (edge_nnodes*edge_nnodes)*(edge_nnodes-iloop) + (edge_nnodes)*(qdeg-iloop-(floop-1)-iedge_pt) + iloop + floop - 1;
                                    fnode ++;
                                }
                            }
                        }
                    }
                }
            }

            /* if have only one center node left (Q is a even number) */
            if (inode == nnodes-1) {
                iloop++;
                map[inode++] = (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes)*(iloop-1) + iloop - 1;
            }
            /* if have only 4 vertices left (Q is a odd number) */
            if (inode == nnodes-8) {
                iloop++;
                map[inode++] = (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes)*(iloop-1) + iloop - 1;
                map[inode++] = (edge_nnodes*edge_nnodes*(iloop-1)) +(edge_nnodes)*iloop - (iloop-1) - 1;
                map[inode++] = (edge_nnodes*edge_nnodes*(iloop-1)) +(edge_nnodes)*(edge_nnodes-(iloop-1)) - (iloop-1) - 1;
                map[inode++] = (edge_nnodes*edge_nnodes*(iloop-1)) +(edge_nnodes)*(qdeg-(iloop-1)) + iloop - 1;
                map[inode++] = (edge_nnodes*edge_nnodes*qdeg) - (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes)*(iloop-1) + iloop - 1;
                map[inode++] = (edge_nnodes*edge_nnodes*qdeg) - (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes)*iloop - (iloop-1) - 1;
                map[inode++] = (edge_nnodes*edge_nnodes*qdeg) - (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes)*(edge_nnodes-(iloop-1)) - (iloop-1) - 1;
                map[inode++] = (edge_nnodes*edge_nnodes*qdeg) - (edge_nnodes*edge_nnodes*(iloop-1)) + (edge_nnodes)*(qdeg-(iloop-1)) +  iloop - 1;
            }

        }
    }

}

static void msh_zorder_node_numbering(int grid_dim,int qdegree,Uint *map){
    /* ================================================ */
    /* maps right-hand-rule node ordering to z-ordering */
    /* ================================================ */
    /* map[rhs index] = z index */
    switch (grid_dim) {
        case (2): /* 2D elements: quads */
            if (qdegree == 1) {
                map[0]  = 0;  map[1]  = 1;  map[2]  = 3;   map[3] =  2; // vertices
            } else
            if (qdegree == 2) {
                map[0]  = 0;  map[1]  = 2;  map[2]  = 8;   map[3]  = 6; // vertices
                map[4]  = 1;  map[5]  = 5;  map[6]  = 7;   map[7]  = 3; // edges
                map[8]  = 4; //faces
            } else
            if (qdegree == 3) {
                map[0]  = 0;  map[1]  = 3;  map[2]  = 15;  map[3]  = 12; // vertices
                map[4]  = 1;  map[5]  = 2;  map[6]  = 7;   map[7]  = 11; map[8] = 14; map[9] = 13; map[10] = 8; map[11] = 4; //edges
                map[12] = 5;  map[13] = 6;  map[14] = 10;  map[15] = 9; // faces
            } else
            if (qdegree == 4) {
                map[0]  = 0;  map[1]  = 4;  map[2]  = 24;  map[3]  = 20; // vertices
                map[4]  = 1;  map[5]  = 2;  map[6]  = 3;   //edge 1
                map[7]  = 9;  map[8]  = 14; map[9]  = 19;  //edge 2
                map[10] = 23; map[11] = 22; map[12] = 21;  //edge 3
                map[13] = 15; map[14] = 10; map[15] = 5;   //edge 4
                map[16] = 6;  map[17] = 8;  map[18] = 18; map[19] = 16;// interior vertices
                map[20] = 7;  map[21] = 13; map[22] = 17; map[23] = 11;// interior edges
                map[24] = 12; // interior volume
            } else {
                /* use Q-flexible mapping */
                msh_map_node_numbering_flex(grid_dim,qdegree,map);
            }
            break;
        case (3): /* 3D elements: Hexahedron */
            if (qdegree == 1) {
                map[0]  = 0;  map[1]  = 1;  map[2]  = 3;   map[3]  =  2;   map[4]  = 4;  map[5]  = 5;   map[6]  = 7;   map[7]  =  6; // vertices
            } else
            if (qdegree == 2) {
                map[0]  = 0;  map[1]  = 2;  map[2]  = 8;   map[3]  =  6;   map[4]  = 18; map[5]  = 20;  map[6]  = 26;  map[7]  =  24; // vertices
                map[8]  = 1;  map[9]  = 3;  map[10] = 9;   map[11] =  5;   map[12] = 11; map[13] =  7;  map[14] = 17;  map[15] =  15;  map[16] =  19;  map[17] =  21;  map[18] = 23;   map[19] =  25; // edges
                map[20] = 4;  map[21] = 10; map[22] = 12;  map[23] =  14;  map[24] = 16; map[25] = 22; //faces
                map[26] = 13; //volume
            } else
            if (qdegree == 3) {
                map[0]  = 0;  map[1]  = 3;  map[2]  = 15;  map[3] =  12;   map[4]  = 48; map[5]  = 51;  map[6]  = 63;   map[7] =  60; // vertices
                map[8]  = 1;  map[9]  = 2;  map[10] = 4;   map[11] =  8;   map[12] = 16; map[13] = 32;  map[14] =  7;   map[15] = 11;  map[16] =  19;  map[17] =  35;  map[18] = 14;   map[19] = 13;   map[20] =  31;   map[21] =  47;   map[22] = 28;   map[23] = 44;   map[24] =  49;   map[25] = 50;   map[26] = 52;   map[27] =  56;   map[28] = 55;   map[29] = 59;   map[30] =  62;   map[31] =  61; //edges
                map[32] = 5;  map[33] = 9;  map[34] = 10;  map[35] =  6;   map[36] = 17; map[37] = 18;  map[38] = 34;   map[39] = 33;  map[40] =  20;  map[41] =  36;  map[42] = 40;   map[43] = 24;   map[44] =  23;   map[45] =  27;   map[46] = 43;   map[47] = 39;   map[48] =  30;   map[49] = 29;   map[50] = 45;   map[51] =  46;   map[52] = 53;   map[53] = 54;   map[54] =  58;   map[55] =  57; // faces
                map[56] = 21; map[57] = 22; map[58] = 26;  map[59] = 25;   map[60] = 37; map[61] = 38;  map[62] = 42;   map[63] = 41; // volume (RHS -> kmin,kmax)
            } else
            if (qdegree == 4) {
                map[0]   = 0;   map[1]   = 4;   map[2]   = 24;  map[3]   = 20;  map[4]   = 100; map[5]  = 104;  map[6]  = 124;  map[7]  = 120; // vertices
                map[8]   = 1;   map[9]   = 2;   map[10]  = 3;   //edge 1
                map[11]  = 5;   map[12]  = 10;  map[13]  = 15;  //edge 2
                map[14]  = 25;  map[15]  = 50;  map[16]  = 75;  //edge 3
                map[17]  = 9;   map[18]  = 14;  map[19]  = 19;  //edge 4
                map[20]  = 29;  map[21]  = 54;  map[22]  = 79;  //edge 5
                map[23]  = 23;  map[24]  = 22;  map[25]  = 21;  //edge 6
                map[26]  = 49;  map[27]  = 74;  map[28]  = 99;  //edge 7
                map[29]  = 45;  map[30]  = 70;  map[31]  = 95;  //edge 8
                map[32]  = 101; map[33]  = 102; map[34]  = 103; //edge 9
                map[35]  = 105; map[36]  = 110; map[37]  = 115; //edge 10
                map[38]  = 109; map[39]  = 114; map[40]  = 119; //edge 11
                map[41]  = 123; map[42]  = 122; map[43]  = 121; //edge 12
                map[44]  = 6;   map[45]  = 16;  map[46]  = 18;  map[47]  = 8;   map[48]  = 11;  map[49]  = 17;  map[50]  = 13;  map[51]  = 7;   map[52]  = 12;  // face 1
                map[53]  = 26;  map[54]  = 28;  map[55]  = 78;  map[56]  = 76;  map[57]  = 27;  map[58]  = 53;  map[59]  = 77;  map[60]  = 51;  map[61]  = 52;  // face 2
                map[62]  = 30;  map[63]  = 80;  map[64]  = 90;  map[65]  = 40;  map[66]  = 55;  map[67]  = 85;  map[68]  = 65;  map[69]  = 35;  map[70]  = 60;  // face 3
                map[71]  = 34;  map[72]  = 44;  map[73]  = 94;  map[74]  = 84;  map[75]  = 39;  map[76]  = 69;  map[77]  = 89;  map[78]  = 59;  map[79]  = 64;  // face 4
                map[80]  = 48;  map[81]  = 46;  map[82]  = 96;  map[83]  = 98;  map[84]  = 47;  map[85]  = 71;  map[86]  = 97;  map[87]  = 73;  map[88]  = 72;  // face 5
                map[89]  = 106; map[90]  = 108; map[91]  = 118; map[92]  = 116; map[93]  = 107; map[94]  = 113; map[95]  = 117; map[96]  = 111; map[97]  = 112; // face 6
                map[98]  = 31;  map[99]  = 33;  map[100] = 43;  map[101] = 41;  map[102] = 81;  map[103] = 83;  map[104] = 93;  map[105] = 91; // interior vertices
                map[106] = 32;  map[107] = 36;  map[108] = 56;  map[109] = 38;  map[110] = 58;  map[111] = 42;  map[112] = 68;  map[113] = 66;  map[114] = 82;   map[115] = 86;   map[116] = 88;   map[117] = 92;   // interior edges
                map[118] = 37;  map[119] = 57;  map[120] = 61;  map[121] = 63;  map[122] = 67;  map[123] = 87;     // interior faces
                map[124] = 62;   // interior volume
            } else {
                /* use Q-flexible mapping */
                msh_map_node_numbering_flex(grid_dim,qdegree,map);
            }
            break;
    }
}

p4est_connectivity_t* gridfile_read_failure(p4est_connectivity_t *conn){
    if(conn) p4est_connectivity_destroy(conn);
    return NULL;
}

 /* ===================== */
 /*   GRID FILE READER    */
 /* ===================== */
 /* Supported grid files: */
 /*  - .msh               */
 /* ===================== */
p4est_connectivity_t *gridfile_reader(ctx_t *ctx,const char *filename){
    p4est_topidx_t tree, tv;
    double t1,t2;
    int retval;
    int face;
    int i;

    p4est_connectivity_t *conn = NULL;
    p4est_topidx_t *periodic_tree_to_vertex = NULL;
    p4est_topidx_t *original_tree_to_vertex = NULL;
    p4est_topidx_t *pv2v = NULL;
    FILE *fid = NULL;

    /* ================================= *
     * Rank 0:                           *
     *  1. construct p4est connectivity  *
     *  2. communicate conn and geometry *
     * ================================= */
    if (ctx->rank == 0) {
        DGOUT(ctx->log_io,"[GMSH] Reading connectivity from %s\n",filename);

        /* open grid file */
        fid = fopen(filename, "rb");
        if (fid == NULL) {
            DGOUT(stderr,"\x1B[1;31m[GMSH] Failed to open gridfile %s!\x1B[0m\n",filename);
            exit(EXIT_FAILURE);
        }

        /* read gmsh file and construct connectivity */
        t1 = MPI_Wtime();
        conn = gridfile_reader_msh_stream(fid,ctx,&pv2v);
        t2 = MPI_Wtime();

        /* check user-input periodic conditions */
        if (check_periodic_inputs(ctx)) {
            if(pv2v == NULL) pv2v = P4EST_ALLOC(p4est_topidx_t, conn->num_vertices);
            for(i = 0; i < conn->num_vertices; i++) pv2v[i] = i;

            gridfile_periodic_boundary_construction(ctx,conn,pv2v);
        }

        /* close file */
        retval = fclose(fid);
        if (retval) {
            DGOUT(stderr,"\x1B[1;31m[GMSH] Failed to close gridfile %s!\x1B[0m\n",filename);
            exit(EXIT_FAILURE);
        }

        /* check for errors */
        if (conn == NULL || (conn != NULL && conn->num_vertices == 0 && conn->num_trees == 0)) {
            DGOUT(stderr,"\x1B[1;31m[GMSH] Failed to read gridfile %s\x1B[0m\n",filename);
            exit(EXIT_FAILURE);
        }

        /* display mesh statistics */
        DGOUT(ctx->log_io,"\x1B[1;32m[GMSH] Grid Read Wall Time: %f sec\x1B[0m\n",t2-t1);
    }

    /* ===================================== */
    /* communicate connectivity to all cores */
    /* ===================================== */
    if (ctx->nranks > 1) {
        p4est_topidx_t num_vertices = (conn != NULL) ? conn->num_vertices:0;
        p4est_topidx_t num_trees = (conn != NULL) ? conn->num_trees:0;

        /* communicate number of vertices and trees */
	if(ctx->rank == 0) DGOUT(ctx->log_io,"[GMSH] Broadcasting number of vertices and trees...");
        MPI_Bcast(&num_vertices,1,P4EST_MPI_TOPIDX,0,ctx->comm);
        MPI_Bcast(&num_trees,1,P4EST_MPI_TOPIDX,0,ctx->comm);
        MPI_Barrier(ctx->comm);
	if(ctx->rank == 0) DGOUT(ctx->log_io,"done!\n");

        /* allocate connectivity on all other ranks */
        if (ctx->rank != 0) {
            conn = p4est_connectivity_new(num_vertices,num_trees,
                                    ARG3D(0)
                                    ARG3D(0)
                                          0,
                                          0);
        }

        /* communicate vertices and tree_to_vertex */
	if(ctx->rank == 0) DGOUT(ctx->log_io,"[GMSH] Broadcasting p4est vertices...");
        MPI_Bcast(conn->vertices,3*num_vertices,MPI_DOUBLE,0,ctx->comm);
        MPI_Barrier(ctx->comm);
	if(ctx->rank == 0) DGOUT(ctx->log_io,"done!\n");

	if(ctx->rank == 0) DGOUT(ctx->log_io,"[GMSH] Broadcasting p4est tree_to_vertex...");
        MPI_Bcast(conn->tree_to_vertex,P4EST_CHILDREN*num_trees,P4EST_MPI_TOPIDX,0,ctx->comm);
        MPI_Barrier(ctx->comm);
	if(ctx->rank == 0) DGOUT(ctx->log_io,"done!\n");

        /* ====================================== */
        /* communicate geometry info to all cores */
        /* ====================================== */
	if(ctx->rank == 0) DGOUT(ctx->log_io,"[GMSH] Broadcasting geometry data...");
        geometry_broadcast_data(ctx);
        MPI_Barrier(ctx->comm);
	if(ctx->rank == 0) DGOUT(ctx->log_io,"done!\n");

	if(ctx->rank == 0) DGOUT(ctx->log_io,"[GMSH] Broadcasting gridfile data...");
        gridfile_broadcast_data(ctx,num_trees);
        MPI_Barrier(ctx->comm);
	if(ctx->rank == 0) DGOUT(ctx->log_io,"done!\n");

        /* communicate pv2v to all cores if periodic */
        char periodic = (pv2v != NULL) ? 1:0; // non-zero on root if periodic
        MPI_Bcast(&periodic,1,MPI_CHAR,0,ctx->comm);

        /* allocate and communicate pv2v */
        if (periodic) {
            if(ctx->rank != 0) pv2v = P4EST_ALLOC(p4est_topidx_t, num_vertices);

	    if(ctx->rank == 0) DGOUT(ctx->log_io,"[GMSH] Broadcasting periodic vertices...");
            MPI_Bcast(pv2v,num_vertices,P4EST_MPI_TOPIDX,0,ctx->comm);
            MPI_Barrier(ctx->comm);
            if(ctx->rank == 0) DGOUT(ctx->log_io,"done!\n");
        }
    }

    /* ============================ */
    /* periodic boundary conditions */
    /* ============================ */
    /* replace slave nodes in tree_to_vertex list and store in periodic ttv */
    if (pv2v != NULL) {
        periodic_tree_to_vertex = P4EST_ALLOC (p4est_topidx_t, P4EST_CHILDREN*conn->num_trees);
        for (tv = 0; tv < P4EST_CHILDREN * conn->num_trees; tv++) {
            periodic_tree_to_vertex[tv] = pv2v[conn->tree_to_vertex[tv]];
        }

        /* save pointer to original and point conn to periodic tree_to_vertex */
        original_tree_to_vertex = conn->tree_to_vertex;
        conn->tree_to_vertex = periodic_tree_to_vertex;
    }
    MPI_Barrier(ctx->comm);

    /* ======================================== *
     * REQUIRED:                                *
     * --------                                 *
     *  Fill tree_to_tree and tree_to_face to   *
     *  make sure we have a valid connectivity. *
     * ======================================== */
    if(ctx->rank == 0) DGOUT(ctx->log_io,"[GMSH] Filling p4est tree_to_X data structures...");
    for (tree = 0; tree < conn->num_trees; ++tree) {
        for (face = 0; face < P4EST_FACES; ++face) {
            conn->tree_to_tree[P4EST_FACES*tree + face] = tree;
            conn->tree_to_face[P4EST_FACES*tree + face] = face;
        }
    }
    MPI_Barrier(ctx->comm);
    if(ctx->rank == 0) DGOUT(ctx->log_io,"done!\n");

    retval = p4est_connectivity_is_valid (conn) == 0;
    MPI_Barrier(ctx->comm);
    if (retval) {
        DGOUT(stdout,"\x1B[1;31m[GMSH] Failed to create connectivity structure from gridfile %s!\x1B[0m\n",filename);
	sleep(2);
        exit(EXIT_FAILURE);
    } else {
	if(ctx->rank == 0) DGOUT(ctx->log_io,"[GMSH] p4est_connectivity is valid!\n");
    }

    /* compute real tree_to_* fields and complete (edge and) corner fields. */
    if(ctx->rank == 0) DGOUT(ctx->log_io,"[GMSH] Completing p4est connectivity...");
    p4est_connectivity_complete(conn);
    MPI_Barrier(ctx->comm);
    if(ctx->rank == 0) DGOUT(ctx->log_io,"done!\n");

    /* swap back original tree_to_vertex, destroy pv2v and periodic_tree_to_vertex */
    if (pv2v != NULL) {
        conn->tree_to_vertex = original_tree_to_vertex;
        P4EST_FREE(periodic_tree_to_vertex);
        P4EST_FREE(pv2v);
    }

    /* print grid info */
    if (ctx->rank == 0) {
        DGOUT(ctx->log_io,"[GMSH] New connectivity with %lld trees and %lld vertices\n",
            (long long) conn->num_trees,
            (long long) conn->num_vertices);
    }
    return conn;
}

p4est_connectivity_t* gridfile_reader_msh_stream(FILE *stream,ctx_t *ctx,
                                                 p4est_topidx_t **periodic_v2v){
    simulation_t *sim = &ctx->d_simulation;
    geometry_t *geom = &ctx->d_geometry;
    gridfile_t *gfile = &ctx->d_gridfile;
    int grid_dim = sim->dim;

    char *line;
    int lines_read,lines_free;
    int swap = 0,binary = 0;
    int retval;
    int int_version;
    Real version;

    char buff[BUFF_SIZE] = {'\0'};
    char section_srt[BUFF_SIZE] = { '\0' };
    char section_end[BUFF_SIZE] = { '\0' };
    char gmsh_physicalnames = 0;
    char gmsh_entities = 0;

    p4est_topidx_t num_nodes = 0;
    p4est_topidx_t num_elements = 0;
    p4est_topidx_t num_faces = 0;

    p4est_topidx_t *tree_to_vertex_all = NULL;
    p4est_topidx_t *ttv_map = NULL;
    double *vertices_all = NULL;
    char *corner_nodes = NULL;

    /* return vector */
    p4est_connectivity_t *conn = NULL;
    *periodic_v2v = NULL;

    for (;;) {

        /* Start reading file */
        line = get_line(stream,++lines_read);

        /* if EOF (NULL) exit */
        if(line == NULL) break;

        while (line[0] != '$') {
            free(line); ++lines_free;
            line = get_line(stream,++lines_read);
            if (line == NULL) {
                free(line); ++lines_free;
                break;
            }
        }

        /* get section name and section end name */
        sprintf(section_end,"END");
        strcpy(section_srt,&line[1]);
        strcat(section_end,section_srt);

        /* =========== */
        /* MESH FORMAT */
        /* =========== */
        if (!strncmp(&line[1], "MESHFORMAT", 10)) {
            int data_size;

            /* read next line */
            free(line); ++lines_free;
            line = get_line(stream,++lines_read);

            sprintf(buff,"%s %s %s",RealFormat,IntFormat,IntFormat);
            retval = sscanf(line,buff,&version,&binary,&data_size);
            if (retval != 3) {
                DGOUT(ctx->log_io,"\x1B[1;31m[GMSH] Invalid MeshFormat section. File corrupted!\x1B[0m\n");
                free(line); ++lines_free;
                return gridfile_read_failure(conn);
            }

            /* if binary read extra line (endianness) */
            if (binary) {
                if (ctx->rank == 0 && ctx->log_info < DGLOG_ALL) {
                    DGOUT(ctx->log_io,"[GMSH] Grid file is in binary format %.1f\n",version);
                }

                int one;
                if(fread(&one, sizeof(int), 1, stream) != 1) return gridfile_read_failure(conn);
                if(one != 1) swap = 1;
            }

            /* binary files with obsolete versions are not supported */
            if (binary && version < 4.1) {
                DGOUT(ctx->log_io,"\x1B[1;31m[GMSH] Can only read MSH 4.0 (or lower) format in ASCII mode\n");
                free(line); ++lines_free;
                return gridfile_read_failure(conn);
            }

            /* Make integer version*/
            int_version = (int)(version*10);
        }

        /* ==================== */
        /* GRID FILE FORMAT 4.1 */
        /* ==================== */
        if (int_version == 41) {
            /* ============== */
            /*  PHYSICALNAMES */
            /* ============== */
            if (!strncmp(&line[1], "PHYSICALNAMES", 13)) {
                gmsh_physicalnames = 1;

                int nphys_tag;
                int iphys_tag;

                /* read next line */
                free(line); ++lines_free;
                line = get_line(stream,++lines_read);
                retval = sscanf(line, "%d", &nphys_tag);

                /* record number of physical tags */
                geom->npatch = nphys_tag;
                geometry_allocate_patch_data(geom);

                for (iphys_tag = 0; iphys_tag < nphys_tag; iphys_tag++) {
                    /* read next line */
                    free(line); ++lines_free;
                    line = get_line(stream,++lines_read);

                    /** NOTE: %[^\n] is used to write the remaining line to phys_name.
                     *        This is needed because the physical name could be potentially
                     *        two separate words (i.e. "Farfield Inflow")
                     */
                    int dim,loc;
                    char patch_name[BUFF_SIZE];
                    retval = sscanf(line, "%d %d %[^\n]", &dim, &loc, patch_name);

                    if (loc > nphys_tag) {
                        DGOUT(ctx->log_io,
                              "\x1B[1;31m[GMSH] Physical patch ID is larger than the overall number of patches. "
                              "Make sure the patches are tagged starting from 1 in the grid file!\x1B[0m\n");
                        free(line); ++lines_free;
                        return gridfile_read_failure(conn);
                    }

                    /* store in correct order */
                    geom->patch_dim[loc-1] = dim;
                    strcpy(geom->patch_name[loc-1],patch_name);

                    /* associate wake3d tag to patch */
                    if(strstr(geom->wake3d_obc,patch_name) != NULL) geom->wake3d_bc[loc-1] = W3D_OBC;
                    if(strstr(geom->wake3d_wbc,patch_name) != NULL) geom->wake3d_bc[loc-1] = W3D_WBC;
                }
            } // PHYSICALNAMES

            /* ======== */
            /* ENTITIES */
            /* ======== */
            if (!strncmp(&line[1], "ENTITIES", 8)) {
                gmsh_entities = 1;

                size_t nentities[4];
                size_t ient;
                int idim;

                if (binary) {
                    if (fread(nentities, sizeof(size_t), 4, stream) != 4) {
                        return gridfile_read_failure(conn);
                    }
                    if(swap == 1) SwapBytes((char *)nentities, sizeof(size_t), 4);
                } else {
                    /* read next line */
                    free(line); ++lines_free;
                    line = get_line(stream,++lines_read);
                    if (sscanf(line,"%lu %lu %lu %lu",
                                  &nentities[0],
                                  &nentities[1],
                                  &nentities[2],
                                  &nentities[3]) != 4) {
                        return gridfile_read_failure(conn);
                    }
                }

                /* record number of entities */
                geom->nentity_face = (grid_dim == 2) ? nentities[1]:nentities[2];
                geom->nentity_vol  = (grid_dim == 2) ? nentities[2]:nentities[3];
                geometry_allocate_entity_data(geom);

                /* 0D entities: points */
                for (idim = 0; idim < 4; idim++) {
                    for (ient = 0; ient < nentities[idim]; ient++) {
                        /* read only physical tags */
                        if (read_msh4_entity(stream,binary,swap,int_version,idim,
                                            &lines_read,&lines_free,
                                             grid_dim,ient,geom)) {
                            DGOUT(stderr,
                                  "\x1B[1;31m[GMSH] Problem encountered in reading Entity section for dimension %d!\x1B[0m\n",
                                  idim);
                            return gridfile_read_failure(conn);
                        }
                    }
                }
            } // Entities

            /* ===== */
            /* NODES */
            /* ===== */
            if (!strncmp(&line[1], "NODES", 5)) {
                size_t numBlock = 0, minTag = 0, maxTag = 0;
                size_t nnodes = 0;
                size_t i;

                /* ---------------------------- *
                 * NODES Header:                *
                 *      numEntityBlocks(size_t) *
                 *      numNodes(size_t)        *
                 *      minNodeTag(size_t)      *
                 *      maxNodeTag(size_t)      *
                 * ---------------------------- */
                if (binary) {
                    size_t data[4];
                    if(fread(data, sizeof(size_t), 4, stream) != 4) return gridfile_read_failure(conn);
                    if(swap == 1) SwapBytes((char *)data, sizeof(size_t), 4);

                    numBlock = data[0];
                    nnodes = data[1];
                    minTag = data[2];
                    maxTag = data[3];
                } else {
                    /* read next line */
                    free(line); ++lines_free;
                    line = get_line(stream,++lines_read);
                    if (sscanf(line,"%lu %lu %lu %lu",
                               &numBlock,&nnodes,&minTag,&maxTag) != 4) {
                        DGOUT(stderr,"\x1B[1;31m[GMSH] Invalid Nodes section. Corrupted file!\x1B[0m\n");
                        free(line); ++lines_free;
                        return gridfile_read_failure(conn);
                    }
                }
                DGOUT(ctx->log_io,
                      "\x1B[1;32m[GMSH] Mesh Blocks: %lu, Mesh Nodes: %lu, MaxTag: %lu\x1B[0m\n",
                      numBlock,nnodes,maxTag);

                /* set number of nodes: max of nnodes and maxTag */
                num_nodes = (nnodes > maxTag) ? nnodes:maxTag;

                /* allocate all vertices */
                vertices_all = (double *) malloc(3*num_nodes*sizeof(double));

                /* ---------------- */
                /* loop node blocks */
                /* ---------------- */
                for (i = 0; i < numBlock; i++) {
                    int parametric = 0;
                    int entityTag = 0, entityDim = 0;
                    size_t numNodes = 0;
                    size_t inode;

                    if (binary) {
                        int data[3];
                        if (fread(data, sizeof(int), 3, stream) != 3) {
                            printf("\x1B[1;31m[GMSH] Invalid Nodes section. Corrupted file!\x1B[0m\n");
                            return gridfile_read_failure(conn);
                        }
                        if(swap == 1) SwapBytes((char *)data, sizeof(int), 3);
                        entityDim  = data[0];
                        entityTag  = data[1];
                        parametric = data[2];

                        if (fread(&numNodes, sizeof(size_t), 1, stream) != 1) {
                            DGOUT(stderr,"\x1B[1;31m[GMSH] Invalid Nodes section. Corrupted file!\x1B[0m\n");
                            return gridfile_read_failure(conn);
                        }
                        if(swap == 1) SwapBytes((char *)&numNodes, sizeof(size_t), 1);
                    } else {
                        /* read next line */
                        free(line); ++lines_free;
                        line = get_line(stream,++lines_read);
                        if (sscanf(line,"%d %d %d %lu",
                                        &entityDim,&entityTag,&parametric,&numNodes) != 4) {
                            DGOUT(stderr,"\x1B[1;31m[GMSH] Invalid Nodes section. Corrupted file!\x1B[0m\n");
                            free(line); ++lines_free;
                            return gridfile_read_failure(conn);
                        }
                    }

                    /* loop block node tags */
                    size_t *tags = (size_t *) malloc(numNodes*sizeof(size_t));
                    double *coord = (double *) malloc((3*numNodes)*sizeof(double));

                    /* ---------------- */
                    /* read block nodes */
                    /* ---------------- */
                    if (binary) {
                        /* read block node tags */
                        if (fread(&tags[0], sizeof(size_t), numNodes, stream) != numNodes) {
                            DGOUT(stderr,"\x1B[1;31m[GMSH] Error in reading block node numbers!\x1B[0m\n");
                            free(tags); tags = NULL;
                            free(coord); coord = NULL;
                            return gridfile_read_failure(conn);
                        }
                        if(swap == 1) SwapBytes((char *)&tags[0], sizeof(size_t), numNodes);

                        /* read block node coordinates */
                        if (fread(&coord[0], sizeof(double), 3*numNodes, stream) != 3*numNodes) {
                            DGOUT(stderr,"\x1B[1;31m[GMSH] Error in reading block node coordinates!\x1B[0m\n");
                            free(tags); tags = NULL;
                            free(coord); coord = NULL;
                            return gridfile_read_failure(conn);
                        }
                        if(swap == 1) SwapBytes((char *)&coord[0], sizeof(double), 3*numNodes);

                        /* fill data */
                        for (inode = 0; inode < numNodes; inode++) {
                            vertices_all[3*(tags[inode]-1) + 0] = coord[3*inode+0];
                            vertices_all[3*(tags[inode]-1) + 1] = coord[3*inode+1];
                            vertices_all[3*(tags[inode]-1) + 2] = coord[3*inode+2];
                        }
                    } else {
                        /* loop block nodes tags */
                        for (inode = 0; inode < numNodes; inode++) {
                            /* read next line */
                            free(line); ++lines_free;
                            line = get_line(stream,++lines_read);

                            /* read node tag */
                            retval = sscanf(line, "%lu", &tags[inode]);
                            if (retval != 1) {
                                DGOUT(stderr,"\x1B[1;31m[GMSH] Invalid Nodes section. Corrupted file!\x1B[0m\n");
                                free(line); ++lines_free;
                                free(tags); tags = NULL;
                                free(coord); coord = NULL;
                                free(vertices_all); vertices_all = NULL;
                                return gridfile_read_failure(conn);
                            }
                            if (tags[inode] > num_nodes) {
                                printf("Encountered vertex %lu that will not fit in vertices"
                                       " array of length %ld.  Are the vertices contiguously"
                                       " numbered?\n", tags[inode], (long int) num_nodes);
                                free(line);
                                free(tags); tags = NULL;
                                free(coord); coord = NULL;
                                free(vertices_all); vertices_all = NULL;
                                return gridfile_read_failure(conn);
                            }
                        }

                        /* set ascii coordinate format */
                        sprintf(buff,"%s %s %s",RealFormat,RealFormat,RealFormat);

                        /* loop block nodes coordinates */
                        for (inode = 0; inode < numNodes; inode++) {
                            /* read next line */
                            free(line); ++lines_free;
                            line = get_line(stream,++lines_read);

                            /* read node coordinates */
                            retval = sscanf(line,buff,
                                            &coord[3*inode+0],
                                            &coord[3*inode+1],
                                            &coord[3*inode+2]);
                            if (retval != 3) {
                                DGOUT(stderr,"\x1B[1;31m[GMSH] Invalid Nodes section. Corrupted file!\x1B[0m\n");
                                free(line); ++lines_free;
                                free(tags); tags = NULL;
                                free(coord); coord = NULL;
                                free(vertices_all); vertices_all = NULL;
                                return gridfile_read_failure(conn);
                            }
                            vertices_all[3*(tags[inode]-1) + 0] = coord[3*inode+0];
                            vertices_all[3*(tags[inode]-1) + 1] = coord[3*inode+1];
                            vertices_all[3*(tags[inode]-1) + 2] = coord[3*inode+2];
                        }
                    }

                    /* free memory */
                    free(tags); tags = NULL;
                    free(coord); coord = NULL;
                }
                DGOUT(ctx->log_io,"[GMSH] Done reading Nodes Section.\n");
            } // NODES

            /* ======== */
            /* ELEMENTS */
            /* ======== */
            if (!strncmp(&line[1], "ELEMENTS", 8)) {
                size_t numBlocks = 0, totalNumElements = 0;
                size_t minTag = 0, maxTag = 0;
                size_t element_counter = 0, face_counter = 0;
                size_t nelem_coord_bytes;
                size_t ielem;

                char supported_vol_elem_type, supported_face_elem_type;
                int iblk,nvert,v;

                /* ---------------------------- *
                 * ELEMENTS Header:             *
                 *      numEntityBlocks(size_t) *
                 *      numElements(size_t)     *
                 *      minElementTag(size_t)   *
                 *      maxElementTag(size_t)   *
                 * ---------------------------- */
                if (binary) {
                    size_t data[4];
                    if(fread(data, sizeof(size_t), 4, stream) != 4) return gridfile_read_failure(conn);
                    if(swap == 1) SwapBytes((char *)data, sizeof(size_t), 4);

                    numBlocks        = data[0];
                    totalNumElements = data[1];
                    minTag           = data[2];
                    maxTag           = data[3];
                } else {
                    free(line); ++lines_free;
                    line = get_line(stream,++lines_read);
                    if (sscanf(line,"%lu %lu %lu %lu",
                               &numBlocks,&totalNumElements,&minTag,&maxTag) != 4) {
                        DGOUT(stderr,"\x1B[1;31m[GMSH] Error in reading element header!\x1B[0m\n");
                        free(line); ++lines_free;
                        return gridfile_read_failure(conn);
                    }
                }

                /* allocate gridfile and face node data*/
                gridfile_allocate_data(gfile,totalNumElements);
                geometry_allocate_facenode_data(geom,2*DIM*totalNumElements);

                /* allocate full local data structures */
                tree_to_vertex_all = (p4est_topidx_t *) malloc(totalNumElements*P4EST_CHILDREN*sizeof(p4est_topidx_t));
                ttv_map = (p4est_topidx_t *) malloc(num_nodes*sizeof(p4est_topidx_t));
                corner_nodes = (char *) calloc(num_nodes, sizeof(char)); // zeroed by calloc!

                /* open binary file for dumping the element nodes in order */
                FILE *fp = fopen("QNODES.bin","wb");
                char element_node_header_complete = 0;

                /* loop element blocks */
                for (iblk = 0; iblk < numBlocks; iblk++) {
                    int entityDim = 0, entityTag = 0, elementType = 0;
                    size_t numElements = 0;

                    /* ------------------------------- *
                     * Element Block Data:             *
                     *      entityDim(int)             *
                     *      entityTag(int)             *
                     *      elementType(int)           *
                     *      numElementsInBlock(size_t) *
                     * ------------------------------- */
                    if (binary) {
                        int data[3];
                        if (fread(data, sizeof(int), 3, stream) != 3) {
                            DGOUT(stderr,"\x1B[1;31m[GMSH] Error in reading block element header!\x1B[0m\n");
                            return gridfile_read_failure(conn);
                        }
                        if(swap == 1) SwapBytes((char *)data, sizeof(int), 3);

                        entityDim   = data[0];
                        entityTag   = data[1];
                        elementType = data[2];

                        if(fread(&numElements, sizeof(size_t), 1, stream) != 1) return gridfile_read_failure(conn);
                        if(swap == 1) SwapBytes((char *)&numElements, sizeof(size_t), 1);
                    } else {
                        free(line); ++lines_free;
                        line = get_line(stream,++lines_read);
                        if (sscanf(line,"%d %d %d %lu",
                                   &entityDim,&entityTag,&elementType,&numElements) != 4) {
                            DGOUT(stderr,"\x1B[1;31m[GMSH] Error in reading block element header!\x1B[0m\n");
                            free(line); ++lines_free;
                            return gridfile_read_failure(conn);
                        }
                    }

                    /* If this block contains supported elements, log the number of elements*/
                    msh_supported_elem_type(elementType,grid_dim,
                                           &supported_vol_elem_type,
                                           &supported_face_elem_type,
                                           &nvert);

                    if(supported_vol_elem_type) {num_elements += numElements;}
                    if(supported_face_elem_type){num_faces    += numElements;}

                    /* write element node header: bytes/element, elementType */
                    if (!element_node_header_complete && supported_vol_elem_type) {
                        nelem_coord_bytes = 3*nvert*sizeof(double);
                        size_t header[2] = {nelem_coord_bytes,(size_t)elementType};

                        fwrite(header,sizeof(size_t),2,fp);
                        element_node_header_complete = 1;

                        DGOUT(ctx->log_io,"\x1B[1;32m[GMSH] "
                              "Constructing binary element-node file: elementType-%d, %d nodes/element"
                              "\x1B[0m\n",elementType,nvert);
                    }

                    if (binary) {
                        size_t n = 1 + nvert; // elementTag + Nodes[nvert]
                        size_t *data = (size_t *) malloc((n*numElements)*sizeof(size_t));
                        int v;

                        /* read elementTag(size_t) + nodeTags(size_t) for all elements in block */
                        if (fread(&data[0],sizeof(size_t),n*numElements,stream) != n*numElements) {
                            DGOUT(stderr,"\x1B[1;31m[GMSH] Error in reading block element header!\x1B[0m\n");
                            free(data); data = NULL;
                            return gridfile_read_failure(conn);
                        }
                        if(swap == 1) SwapBytes((char *)&data[0], sizeof(size_t), n*numElements);

                        /* unpack read elements */
                        for (ielem = 0; ielem < numElements; ielem++) {
                            if (supported_vol_elem_type) {
                                /* record volume element entity tag */
                                gfile->elem_tag[element_counter] = entityTag-1; // zero-based

                                /* z-ordered corner nodes */
                                size_t *elem_nodes = &data[n*ielem+1]; // jump to element and skip elementTag

                                size_t sort[P4EST_CHILDREN];
                                sort[0] = elem_nodes[0];
                                sort[1] = elem_nodes[1];
                                sort[2] = elem_nodes[3];
                                sort[3] = elem_nodes[2];
                            E3D(sort[4] = elem_nodes[4])
                            E3D(sort[5] = elem_nodes[5])
                            E3D(sort[6] = elem_nodes[7])
                            E3D(sort[7] = elem_nodes[6])
                                /* hold only corner nodes for p4est */
                                for (v = 0; v < P4EST_CHILDREN; ++v) {
                                    tree_to_vertex_all[P4EST_CHILDREN*element_counter + v] = sort[v]-1; // zero-based
                                }

                                /* flag the corner node ID */
                                for(v = 0; v < P4EST_CHILDREN; ++v) corner_nodes[sort[v]-1] = 1;

                                /* write element vertices */
                                for (v = 0; v < nvert; ++v){
                                    size_t node = elem_nodes[v]-1; // zero-based
                                    double *node_geom = &vertices_all[3*node];

                                    /* write x,y,z to file */
                                    fwrite(node_geom,sizeof(double),3,fp);
                                }

                                /* increment volume element counter */
                                ++element_counter;
                            }

                            /* store entity tag and element face tag */
                            if (supported_face_elem_type) {
                                int  nnodes = P4EST_CHILDREN/2;
                                size_t this_face = (2+nnodes)*face_counter;

                                geom->face_nodes_info[this_face + EFENT_IND] = entityTag-1;     // zero-based
                                geom->face_nodes_info[this_face + EFTAG_IND] = data[n*ielem]-1; // zero-based

                                /* switch from right-hand vertex ordering to z-order
                                 *  NB: .msh order -> right-hand order
                                 *       p4est     -> z-order
                                 */
                                size_t sort[P4EST_CHILDREN/2];
                                sort[0] = data[n*ielem+1+0];
                                sort[1] = data[n*ielem+1+1];
                            E3D(sort[2] = data[n*ielem+1+3])
                            E3D(sort[3] = data[n*ielem+1+2])
                                for (v = 0; v < nnodes; ++v) {
                                    geom->face_nodes_info[this_face + EFNDS_IND + v] = sort[v]-1; // zero-based
                                }

                                /* increment element face counter */
                                ++face_counter;
                            }
                        }

                        /* free memory */
                        free(data); data = NULL;
                    } else {
                        char str[10000]; // q-degree <= 10: 3*10*10*10

                        /* loop block elements */
                        for (ielem = 0; ielem < numElements; ielem++) {
                            if (!supported_vol_elem_type && !supported_face_elem_type) {
                                /* skip line */
                                free(line); ++lines_free;
                                line = get_line(stream,++lines_read);
                                continue;
                            }

                            if (supported_vol_elem_type) {
                                size_t n = nvert; // Nodes[nvert]
                                size_t data[n];
                                size_t elemTag;

                                /* record volume element entity tag */
                                gfile->elem_tag[element_counter] = entityTag-1; // zero-based

                                /* read elementTag */
                                retval = fscanf(stream, "%lu", &elemTag);
                                if (fgets(str, sizeof(str), stream) == NULL) {
                                    printf("Error reading elementTag in gmsh file!\n");
                                    exit(EXIT_FAILURE);
                                }

                                /* read element nodes */
                                for (v = 0; v < nvert; v++) {
                                    if (v != nvert - 1) {
                                        sscanf(str, "%lu %[0-9- ]", &data[v], str);
                                    } else {
                                        sscanf(str, "%lu", &data[v]);
                                    }
                                }

                                /* Note that when we read in the
                                 * vertices we switch from right-hand
                                 * vertex ordering to z-order
                                 *  NB: .msh order -> right-hand order
                                 *       p4est     -> z-order
                                 */
                                size_t sort[P4EST_CHILDREN];
                                sort[0] = data[0];
                                sort[1] = data[1];
                                sort[2] = data[3];
                                sort[3] = data[2];
                            E3D(sort[4] = data[4])
                            E3D(sort[5] = data[5])
                            E3D(sort[6] = data[7])
                            E3D(sort[7] = data[6])
                                /* hold only corner nodes for p4est */
                                for (v = 0; v < P4EST_CHILDREN; ++v) {
                                    tree_to_vertex_all[P4EST_CHILDREN*element_counter + v] = sort[v]-1; // zero-based
                                }

                                /* flag the corner node ID */
                                for(v = 0; v < P4EST_CHILDREN; ++v) corner_nodes[sort[v]-1] = 1;

                                /* write all element vertices */
                                for (v = 0; v < nvert; ++v){
                                    size_t node = data[v]-1; //zero-based
                                    double *node_geom = &vertices_all[3*node];

                                    /* write x,y,z to file */
                                    fwrite(node_geom,sizeof(double),3,fp);
                                }

                                /* increment volume element counter */
                                element_counter++;
                            }

                            if (supported_face_elem_type) {
                                size_t n = nvert; // Nodes[nvert]
                                size_t data[n];
                                size_t elemTag;

                                /* read faceTag */
                                retval = fscanf(stream, "%lu", &elemTag);
                                if (fgets(str, sizeof(str), stream) == NULL) {
                                    printf("Error reading faceTag in gmsh file!\n");
                                    exit(EXIT_FAILURE);
                                }

                                /* read element nodes */
                                for (v = 0; v < nvert; v++) {
                                    if (v != nvert - 1) {
                                        sscanf(str, "%lu %[0-9- ]", &data[v], str);
                                    } else {
                                        sscanf(str, "%lu", &data[v]);
                                    }
                                }

                                /* switch from right-hand vertex ordering to z-order
                                 *  NB: .msh order -> right-hand order
                                 *       p4est     -> z-order
                                 */
                                size_t sort[P4EST_CHILDREN];
                                sort[0] = data[0];
                                sort[1] = data[1];
                            E3D(sort[2] = data[3])
                            E3D(sort[3] = data[2])

                                /* store entity tag and element face tag */
                                int nnodes = P4EST_CHILDREN/2;
                                size_t this_face = (2+nnodes)*face_counter;

                                geom->face_nodes_info[this_face + EFENT_IND] = entityTag-1; // zero-based
                                geom->face_nodes_info[this_face + EFTAG_IND] = elemTag-1;   // zero-based

                                for (v = 0; v < nnodes; ++v) {
                                    geom->face_nodes_info[this_face + EFNDS_IND + v] = sort[v]-1; // zero-based
                                }

                                /* increment element face counter */
                                ++face_counter;
                            }
                        }
                    }
                }
                /* close binary element node file */
                fclose(fp);

                /* display statistics */
                size_t elem_node_bytes = nelem_coord_bytes*num_elements;
                DGOUT(ctx->log_io,
                      "\x1B[1;32m[GMSH] "
                      "Done constructing binary element-node file (%f GB)"
                      "\x1B[0m\n", elem_node_bytes/1024./1024./1024.);

                /* ======================== */
                /* build p4est connectivity */
                /* ======================== */
                /* [Step 1]: trim vertices list */
                p4est_topidx_t node, ncorner, new_id, num_corners = 0;

                /* count corner nodes: flags were set during element nodes reading */
                for(node = 0; node < num_nodes; ++node) num_corners += corner_nodes[node];

                                double conn_bytes_GB = 0.0;
                double B2GB = 1./1024./1024./1024.;

                size_t vertices_bytes       = 3*num_corners*sizeof(double);
                size_t tree_to_vertex_bytes = P4EST_CHILDREN * num_elements * sizeof(p4est_topidx_t);
                size_t tree_to_tree_bytes   = P4EST_FACES * num_elements * sizeof(p4est_topidx_t);
                size_t tree_to_face_bytes   = P4EST_FACES * num_elements * sizeof(int8_t);

                conn_bytes_GB += B2GB*vertices_bytes +
                                 B2GB*tree_to_vertex_bytes +
                                 B2GB*tree_to_tree_bytes +
                                 B2GB*tree_to_face_bytes;

#ifdef P4_TO_P8
                size_t num_edges = P8EST_EDGES * num_elements;
                size_t num_ett = 4*num_edges; // estimate
                size_t tree_to_edge_bytes = P8EST_EDGES * num_elements * sizeof(p4est_topidx_t);
                size_t edge_to_tree_bytes = num_ett * sizeof(p4est_topidx_t);
                size_t edge_to_edge_bytes = num_ett * sizeof(int8_t);
                size_t ett_offset_bytes   = (num_edges + 1) * sizeof(p4est_topidx_t);

                conn_bytes_GB += B2GB*tree_to_edge_bytes +
                                 B2GB*edge_to_tree_bytes +
                                 B2GB*edge_to_edge_bytes +
                                 B2GB*ett_offset_bytes;
#endif

                size_t num_ctt = 8*num_corners; // estimate
                size_t tree_to_corner_bytes = sizeof(p4est_topidx_t) * P4EST_CHILDREN * num_elements;
                size_t corner_to_tree_bytes = sizeof(p4est_topidx_t) * num_ctt;
                size_t corner_to_corner_bytes = sizeof(int8_t) * num_ctt;
                size_t ctt_offset_bytes = sizeof(p4est_topidx_t) * (num_corners + 1);

                conn_bytes_GB += B2GB*tree_to_corner_bytes +
                                 B2GB*corner_to_tree_bytes +
                                 B2GB*corner_to_corner_bytes +
                                 B2GB*ctt_offset_bytes;

                DGOUT(ctx->log_io,
                      "\x1B[1;32m[GMSH] "
                      "Number of corner nodes counted: %d, elements: %d, conn: %f GB/rank"
                      "\x1B[0m\n",num_corners,num_elements,conn_bytes_GB);

                /* allocate connectivity */
                conn = p4est_connectivity_new(num_corners,num_elements,
                                        ARG3D(0)
                                        ARG3D(0)
                                              0,
                                              0);

                /* stack corner nodes */
                ncorner = 0;
                for (node = 0; node < num_nodes; ++node) {
                    if (corner_nodes[node]) {
                        conn->vertices[3*ncorner+0] = vertices_all[3*node+0];
                        conn->vertices[3*ncorner+1] = vertices_all[3*node+1];
                        conn->vertices[3*ncorner+2] = vertices_all[3*node+2];
                        ncorner++;
                    }
                }
                free(vertices_all); vertices_all = NULL;

                /* construct old_node_id to new_node_id map for corner nodes */
                new_id = 0;
                for (node = 0; node < num_nodes; ++node) {
                    if(corner_nodes[node]) ttv_map[node] = new_id++;
                }
                free(corner_nodes); corner_nodes = NULL;

                /* [Step 2]: fill tree_to_vertex list with trimmed vertex indices */
                for (ielem = 0; ielem < num_elements; ielem++) {
                    for (v = 0; v < P4EST_CHILDREN; ++v) {
                        p4est_topidx_t old_node_id = tree_to_vertex_all[P4EST_CHILDREN*ielem + v];
                        conn->tree_to_vertex[P4EST_CHILDREN*ielem + v] = ttv_map[old_node_id];
                    }
                }
                free(tree_to_vertex_all); tree_to_vertex_all = NULL;

                /* [Step 3]: adjust face vertex indices */
                int iface;
                int nfnodes = P4EST_CHILDREN/2;
                for (iface = 0; iface < num_faces; iface++){
                    size_t this_face = 2*DIM*iface;

                    for (v = 0; v < nfnodes; ++v) {
                        p4est_topidx_t old_node_id = geom->face_nodes_info[this_face + EFNDS_IND + v];
                        geom->face_nodes_info[this_face + EFNDS_IND + v] = ttv_map[old_node_id];
                    }
                }
                DGOUT(ctx->log_io,"[GMSH] Done reading Elements Section.\n");
            } // ELEMENTS

            /* ======== */
            /* PERIODIC */
            /* ======== */
            if (!strncmp(&line[1], "PERIODIC", 8)) {
                p4est_topidx_t NumOfVertices = num_nodes;
                size_t numPeriodicLinks, pair_count = 0;
                p4est_topidx_t i,j;

                /* allocate and initialize master-slave nodes */
                if (NumOfVertices > 0) {
                    *periodic_v2v = P4EST_ALLOC (p4est_topidx_t, NumOfVertices);
                    for(i = 0; i < NumOfVertices; i++) (*periodic_v2v)[i] = i;
                }

                if (binary) {
                    if (fread(&numPeriodicLinks, sizeof(size_t), 1, stream) != 1) {
                        return gridfile_read_failure(conn);
                    }
                    if(swap == 1) SwapBytes((char *)&numPeriodicLinks, sizeof(size_t), 1);
                } else {
                    /* read next line */
                    free(line); ++lines_free;
                    line = get_line(stream,++lines_read);
                    if (sscanf(line,"%lu",&numPeriodicLinks) != 1) {
                        DGOUT(stderr,"\x1B[1;31m[GMSH] Invalid Periodic section. Corrupted file!\x1B[0m\n");
                        free(line); ++lines_free;

                        P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                        return gridfile_read_failure(conn);
                    }
                }

                for (i = 0; i < numPeriodicLinks; i++) {
                    int slaveDim = 0, slaveTag = 0, masterTag = 0;

                    /* read and ignore entity dimension and tags */
                    if (binary) {
                        int data[3];
                        if (fread(data, sizeof(int), 3, stream) != 3) {
                            P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                            return gridfile_read_failure(conn);
                        }
                        if(swap == 1) SwapBytes((char *)data, sizeof(int), 3);

                        slaveDim = data[0];
                        slaveTag = data[1];
                        masterTag = data[2];
                    } else {
                        if (fscanf(stream, "%d %d %d", &slaveDim, &slaveTag, &masterTag) != 3) {
                            P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                            return gridfile_read_failure(conn);
                        }
                    }

                    /* read and ignore affine mapping, read number of vertex pairs */
                    size_t correspondingVertexSize = 0;
                    if (binary) {
                        size_t numAffine;

                        if (fread(&numAffine, sizeof(size_t), 1, stream) != 1) {
                            P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                            return gridfile_read_failure(conn);
                        }
                        if(swap == 1) SwapBytes((char *)&numAffine, sizeof(size_t), 1);
                        if (numAffine) {
                            double tfo[numAffine];
                            if (fread(&tfo[0], sizeof(double), numAffine, stream) != numAffine) {
                                P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                                return gridfile_read_failure(conn);
                            }
                        }

                        /* read number of vertex pairs */
                        if (fread(&correspondingVertexSize, sizeof(size_t), 1, stream) != 1) {
                            P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                            return gridfile_read_failure(conn);
                        }
                        if(swap == 1) SwapBytes((char *)&correspondingVertexSize, sizeof(size_t), 1);
                    } else {
                        if (int_version >= 41) {
                            size_t numAffine;

                            if (!fscanf(stream, "%lu", &numAffine)) {
                                P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                                return gridfile_read_failure(conn);
                            }

                            if (numAffine) {
                                double tfo[numAffine];
                                for (j = 0; j < numAffine; j++) {
                                    if (fscanf(stream, "%lf", &tfo[j]) != 1) {
                                        P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                                        return gridfile_read_failure(conn);
                                    }
                                }
                            }
                            if (fscanf(stream, "%lu", &correspondingVertexSize) != 1) {
                                P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                                return gridfile_read_failure(conn);
                            }
                        }
                    }

                    for (j = 0; j < correspondingVertexSize; j++) {
                        size_t slave = 0, master = 0;
                        if (binary) {
                            size_t data[2];
                            if (fread(data, sizeof(size_t), 2, stream) != 2) {
                                P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                                return gridfile_read_failure(conn);
                            }
                            if(swap == 1) SwapBytes((char *)data, sizeof(size_t), 2);
                            slave  = data[0];
                            master = data[1];
                        } else {
                            if (fscanf(stream, "%lu %lu", &slave, &master) != 2) {
                                P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                                return gridfile_read_failure(conn);
                            }
                        }

                        /* store periodic pair map */
                        /* map the old node slave/max nodes to reduced map */
                        (*periodic_v2v)[ttv_map[slave - 1]] = ttv_map[master - 1];
                        pair_count++;
                    }
                }
                if (pair_count == 0) {
                    P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                }
            } // PERIODIC

        /* ==================== */
        /* GRID FILE FORMAT 2.0 */
        /* ==================== */
        } else
        if (int_version < 30) {
            /* ============= */
            /* PHYSICALNAMES */
            /* ============= */
            if (!strncmp(&line[1], "PHYSICALNAMES", 13)) {
                gmsh_physicalnames = 1;
                gmsh_entities = 1; // this section does not exist in this version but we need to check!

                int nphys_tag;
                int iphys_tag;

                /* read next line */
                free(line); ++lines_free;
                line = get_line(stream,++lines_read);
                retval = sscanf(line, "%d", &nphys_tag);

                /* Record number of physical tags */
                /* NOTE: we need the first slot for null bc!!! */
                 /* record number of physical tags */
                geom->npatch = nphys_tag;
                geometry_allocate_patch_data(geom);

                geom->nentity_face = nphys_tag;
                geom->nentity_vol  = nphys_tag;
                geometry_allocate_entity_data(geom);

                for (iphys_tag = 0; iphys_tag < nphys_tag; iphys_tag++) {
                    /* read next line */
                    free(line); ++lines_free;
                    line = get_line(stream,++lines_read);

                    /** NOTE: %[^\n] is used to write the remaining line to phys_name.
                     *        This is needed because the physical name could be potentially
                     *        two separate words (i.e. "Farfield Inflow")
                     */
                    int dim,loc;
                    char patch_name[BUFF_SIZE];
                    retval = sscanf(line, "%d %d %[^\n]", &dim, &loc, patch_name);

                    if (loc >= geom->npatch) {
                        DGOUT(stderr,"\x1B[1;31m[GMSH] Patch location (%d) "
                                     "larger than number of patches (%d). "
                                     "File may be corrupted!\x1B[0m\n",loc,geom->npatch);
                        free(line); ++lines_free;
                        return gridfile_read_failure(conn);
                    }

                    /* Store in correct order (patches are zero-based in version 2.0)*/
                    geom->patch_dim[loc] = dim;
                    strcpy(geom->patch_name[loc],patch_name);

                    geom->entity_face_tag[iphys_tag]      = iphys_tag;
                    geom->entity_vol_tag[iphys_tag]       = iphys_tag;
                    geom->entity_face_patchtag[iphys_tag] = iphys_tag;
                    geom->entity_vol_patchtag[iphys_tag]  = iphys_tag;

                    /* associate wake3d tag to patch */
                    if(strstr(geom->wake3d_obc,patch_name) != NULL) geom->wake3d_bc[loc] = W3D_OBC;
                    if(strstr(geom->wake3d_wbc,patch_name) != NULL) geom->wake3d_bc[loc] = W3D_WBC;
                }
            } // PHYSICALNAMES

            /* ===== */
            /* NODES */
            /* ===== */
            if (!strncmp(&line[1], "NODES", 5)) {
                int nnodes,inode,inode_tag;
                Real xx, yy, zz;

                /* read next line */
                free(line); ++lines_free;
                line = get_line(stream,++lines_read);
                retval = sscanf(line, "%d", &nnodes);
                if (retval != 1) {
                    DGOUT(stderr,"\x1B[1;31m[GMSH] Invalid Nodes section. Corrupted file!\x1B[0m\n");
                    free(line); ++lines_free;
                    return gridfile_read_failure(conn);
                }

                /* set number of nodes */
                num_nodes = nnodes;

                /* allocate all vertices */
                vertices_all = (double *) malloc(3*nnodes*sizeof(double));

                /* set ascii coordinate format */
                sprintf(buff,"%s %s %s %s","%d",RealFormat,RealFormat,RealFormat);

                /* loop block nodes coordinates */
                for (inode = 0; inode < nnodes; inode++) {
                    /* read next line */
                    free(line); ++lines_free;
                    line = get_line(stream,++lines_read);

                    /* read node coordinates */
                    retval = sscanf(line,buff,&inode_tag,&xx,&yy,&zz);
                    if (retval != 4) {
                        DGOUT(stderr,"\x1B[1;31m[GMSH] Invalid Nodes section. Corrupted file!\x1B[0m\n");
                        free(line); ++lines_free;
                        free(vertices_all); vertices_all = NULL;
                        return gridfile_read_failure(conn);
                    }
                    vertices_all[3*(inode_tag-1) + 0] = xx;
                    vertices_all[3*(inode_tag-1) + 1] = yy;
                    vertices_all[3*(inode_tag-1) + 2] = zz;
                }
            } // NODES

            /* ======== */
            /* ELEMENTS */
            /* ======== */
            if (!strncmp(&line[1], "ELEMENTS", 8)) {
                char supported_vol_elem_type, supported_face_elem_type;
                size_t nelem_coord_bytes;
                int ielem, nelems, elem_tag;
                int elem_type;
                int v,j;

                /* read next line */
                free(line); ++lines_free;
                line = get_line(stream,++lines_read);
                retval = sscanf(line, "%d", &nelems);
                if (retval != 1) {
                    printf("\x1B[1;31m[GMSH] Invalid Elements section. Corrupted file!\x1B[0m\n");
                    free(line); ++lines_free;
                    return gridfile_read_failure(conn);
                }

                /* allocate gridfile and face node data*/
                gridfile_allocate_data(gfile,nelems);
                geometry_allocate_facenode_data(geom,2*DIM*nelems);

                /* allocate full local data structures */
                tree_to_vertex_all = (p4est_topidx_t *) malloc(nelems*P4EST_CHILDREN*sizeof(p4est_topidx_t));
                ttv_map = (p4est_topidx_t *) malloc(num_nodes*sizeof(p4est_topidx_t));
                corner_nodes = (char *) calloc(num_nodes, sizeof(char)); // zeroed by calloc!

                /* open binary file for dumping the element nodes in order */
                FILE *fp = fopen("QNODES.bin","wb");
                char element_node_header_complete = 0;

                /* loop block elements */
                for (ielem = 0; ielem < nelems; ielem++) {
                    int num,tag,numTags;
                    int entityTag;

                    /* read: elem_num, elem_type, number_of_tags */
                    if(fscanf(stream, "%d %d %d", &num, &elem_type, &numTags) != 3) {
                        printf("\x1B[1;31m[GMSH] Invalid Elements section. Corrupted file!\x1B[0m\n");
                        free(line); ++lines_free;
                        return gridfile_read_failure(conn);
                    }

                    /* read tags */
                    for (j = 0; j < numTags; j++) {
                        if (fscanf(stream, "%d", &tag) != 1) {
                            printf("\x1B[1;31m[GMSH] Invalid Elements section. Corrupted file!\x1B[0m\n");
                            free(line); ++lines_free;
                            return gridfile_read_failure(conn);
                        }

                        /* save first tag as entityTag */
                        if(j == 0) entityTag = tag;
                    }

                    /* if block contains supported elements, log the number of elements */
                    int nvert;
                    msh_supported_elem_type(elem_type,grid_dim,
                                           &supported_vol_elem_type,
                                           &supported_face_elem_type,
                                           &nvert);

                    if(supported_vol_elem_type)  num_elements += 1;
                    if(supported_face_elem_type) num_faces += 1;

                    /* read vertices */
                    int indices[nvert];
                    for (j = 0; j < nvert; j++) {
                        if(fscanf(stream, "%d", &indices[j]) != 1) {
                            free(line); ++lines_free;
                            return gridfile_read_failure(conn);
                        }
                    }

                    /* write element node header: bytes/element, elementType */
                    if (!element_node_header_complete && supported_vol_elem_type) {
                        nelem_coord_bytes = 3*nvert*sizeof(double);
                        size_t header[2] = {nelem_coord_bytes,(size_t)elem_type};

                        fwrite(header,sizeof(size_t),2,fp);
                        element_node_header_complete = 1;

                        DGOUT(ctx->log_io,"\x1B[1;32m[GMSH] "
                              "Constructing binary element-node file: elementType-%d, %d nodes/element"
                              "\x1B[0m\n",elem_type,nvert);
                    }

                    /* log volume element */
                    if (supported_vol_elem_type) {
                        int element_counter = num_elements-1;

                        /* record volume element entity tag */
                        gfile->elem_tag[element_counter] = entityTag; // already zero-based in format 2.0

                        /* switch from right-hand vertex ordering to z-order
                         *  NB: .msh order -> right-hand order
                         *       p4est     -> z-order
                         */
                        size_t sort[P4EST_CHILDREN];
                        sort[0] = indices[0];
                        sort[1] = indices[1];
                        sort[2] = indices[3];
                        sort[3] = indices[2];
                    E3D(sort[4] = indices[4])
                    E3D(sort[5] = indices[5])
                    E3D(sort[6] = indices[7])
                    E3D(sort[7] = indices[6])
                        /* hold only corner nodes for p4est */
                        for (v = 0; v < P4EST_CHILDREN; ++v) {
                            tree_to_vertex_all[P4EST_CHILDREN*element_counter + v] = sort[v]-1; // zero-based
                        }

                        /* flag the corner node ID */
                        for(v = 0; v < P4EST_CHILDREN; ++v) corner_nodes[sort[v]-1] = 1;

                        /* write all element vertices */
                        for (v = 0; v < nvert; ++v){
                            size_t node = indices[v]-1; //zero-based
                            double *node_geom = &vertices_all[3*node];

                            /* write x,y,z to file */
                            fwrite(node_geom,sizeof(double),3,fp);
                        }
                    }

                    /* log surface element */
                    if (supported_face_elem_type) {
                        int face_counter = num_faces-1;

                        /* switch from right-hand vertex ordering to z-order
                         *  NB: .msh order -> right-hand order
                         *       p4est     -> z-order
                         */
                        size_t sort[P4EST_CHILDREN];
                        sort[0] = indices[0];
                        sort[1] = indices[1];
                    E3D(sort[2] = indices[3])
                    E3D(sort[3] = indices[2])

                        /* store entity tag and element face tag */
                        int nnodes = P4EST_CHILDREN/2;
                        size_t this_face = (2+nnodes)*face_counter;

                        geom->face_nodes_info[this_face + EFENT_IND] = entityTag;  // already zero-based in format 2.0
                        geom->face_nodes_info[this_face + EFTAG_IND] = elem_tag-1; // zero-based

                        for (v = 0; v < nnodes; ++v) {
                            geom->face_nodes_info[this_face + EFNDS_IND + v] = sort[v]-1; // zero-based
                        }
                    }
                }
                /* close binary element node file */
                fclose(fp);

                /* display statistics */
                size_t elem_node_bytes = nelem_coord_bytes*num_elements;
                DGOUT(ctx->log_io,
                      "\x1B[1;32m[GMSH] "
                      "Done constructing binary element-node file (%f GB)"
                      "\x1B[0m\n", elem_node_bytes/1024./1024./1024.);

                /* ======================== */
                /* build p4est connectivity */
                /* ======================== */
                /* [Step 1]: trim vertices list */
                p4est_topidx_t node, ncorner, new_id, num_corners = 0;

                /* count corner nodes: flags were set during element nodes reading */
                for(node = 0; node < num_nodes; ++node) num_corners += corner_nodes[node];

                DGOUT(ctx->log_io,
                      "\x1B[1;32m[GMSH] "
                      "Number of corner nodes counted: %d, elements: %d"
                      "\x1B[0m\n",num_corners,num_elements);

                /* allocate connectivity */
                conn = p4est_connectivity_new(num_corners,num_elements,
                                        ARG3D(0)
                                        ARG3D(0)
                                              0,
                                              0);

                /* stack corner nodes */
                ncorner = 0;
                for (node = 0; node < num_nodes; ++node) {
                    if (corner_nodes[node]) {
                        conn->vertices[3*ncorner+0] = vertices_all[3*node+0];
                        conn->vertices[3*ncorner+1] = vertices_all[3*node+1];
                        conn->vertices[3*ncorner+2] = vertices_all[3*node+2];
                        ncorner++;
                    }
                }
                free(vertices_all); vertices_all = NULL;

                /* construct old_node_id to new_node_id map for corner nodes */
                new_id = 0;
                for (node = 0; node < num_nodes; ++node) {
                    if(corner_nodes[node]) ttv_map[node] = new_id++;
                }
                free(corner_nodes); corner_nodes = NULL;

                /* [Step 2]: fill tree_to_vertex list with trimmed vertex indices */
                for (ielem = 0; ielem < num_elements; ielem++) {
                    for (v = 0; v < P4EST_CHILDREN; ++v) {
                        p4est_topidx_t old_node_id = tree_to_vertex_all[P4EST_CHILDREN*ielem + v];
                        conn->tree_to_vertex[P4EST_CHILDREN*ielem + v] = ttv_map[old_node_id];
                    }
                }
                free(tree_to_vertex_all); tree_to_vertex_all = NULL;

                /* [Step 3]: adjust face vertex indices */
                int iface;
                int nfnodes = P4EST_CHILDREN/2;
                for (iface = 0; iface < num_faces; iface++){
                    size_t this_face = 2*DIM*iface;

                    for (v = 0; v < nfnodes; ++v) {
                        p4est_topidx_t old_node_id = geom->face_nodes_info[this_face + EFNDS_IND + v];
                        geom->face_nodes_info[this_face + EFNDS_IND + v] = ttv_map[old_node_id];
                    }
                }
            } // ELEMENTS

            /* ======== */
            /* PERIODIC */
            /* ======== */
            if (!strncmp(&line[1], "PERIODIC", 8)) {
                p4est_topidx_t NumOfVertices = num_nodes;
                size_t numPeriodicLinks, pair_count = 0;
                p4est_topidx_t i,j;

                /* allocate and initialize master-slave nodes */
                if (NumOfVertices > 0) {
                    *periodic_v2v = P4EST_ALLOC (p4est_topidx_t, NumOfVertices);
                    for(i = 0; i < NumOfVertices; i++) (*periodic_v2v)[i] = i;
                }

                /* read numPeriodicLinks */
                free(line); ++lines_free;
                line = get_line(stream,++lines_read);
                if (sscanf(line,"%lu",&numPeriodicLinks) != 1) {
                    DGOUT(stderr,"\x1B[1;31m[GMSH] Invalid Periodic section. Corrupted file!\x1B[0m\n");
                    free(line); ++lines_free;

                    P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                    return gridfile_read_failure(conn);
                }

                for (i = 0; i < numPeriodicLinks; i++) {
                    int slaveDim = 0, slaveTag = 0, masterTag = 0;

                    /* read and ignore entity dimension and tags */
                    if (fscanf(stream, "%d %d %d", &slaveDim, &slaveTag, &masterTag) != 3) {
                        P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                        return gridfile_read_failure(conn);
                    }

                    /* read and ignore affine mapping, read number of vertex pairs */
                    /* only versions < 4.0 */
                    size_t correspondingVertexSize = 0;
                    char affine[256];
                    if (!fscanf(stream, "%s", affine)) {
                        P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                        return gridfile_read_failure(conn);
                    }
                    if (!strncmp(affine, "Affine", 6)) {
                        if (!fgets(affine, sizeof(affine), stream)) {
                            P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                            return gridfile_read_failure(conn);
                        }
                        double tfo[16];
                        if (sscanf(affine,
                                "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
                                "%lf %lf %lf %lf",
                                &tfo[0], &tfo[1], &tfo[2], &tfo[3], &tfo[4], &tfo[5],
                                &tfo[6], &tfo[7], &tfo[8], &tfo[9], &tfo[10], &tfo[11],
                                &tfo[12], &tfo[13], &tfo[14], &tfo[15]) != 16) {
                            P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                            return gridfile_read_failure(conn);
                        }
                        if (fscanf(stream, "%lu", &correspondingVertexSize) != 1) {
                            P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                            return gridfile_read_failure(conn);
                        }
                    } else {
                        if (sscanf(affine, "%lu", &correspondingVertexSize) != 1) {
                            P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                            return gridfile_read_failure(conn);
                        }
                    }

                    for (j = 0; j < correspondingVertexSize; j++) {
                        size_t slave = 0, master = 0;
                        if (fscanf(stream, "%lu %lu", &slave, &master) != 2) {
                            P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                            return gridfile_read_failure(conn);
                        }

                        /* store periodic pair map */
                        /* map the old node slave/max nodes to reduced map */
                        (*periodic_v2v)[ttv_map[slave - 1]] = ttv_map[master - 1];
                        pair_count++;
                    }
                }
                if (pair_count == 0) {
                    P4EST_FREE(*periodic_v2v); *periodic_v2v = NULL;
                }
            } // PERIODIC

        } else {
            DGOUT(stderr,"\x1B[1;31m[GMSH] Gmsh file format not supported! Available formats 2.0, 2.2 and 4.1.\x1B[0m\n");
            return gridfile_read_failure(conn);
        }

        while (strncmp(&line[1], section_end, strlen(section_end))) {
            free(line); ++lines_free;
            line = get_line(stream,++lines_read);
            if(line == NULL) break;
        }

        /* free line either were reached EOF or we need to move to next line */
        free(line); ++lines_free;
        if(line == NULL) break;
    }

    /* clean up map memory: need at end in case of periodic elements */
    if (ttv_map) {
        free(ttv_map); ttv_map = NULL;
    }

    if (!gmsh_physicalnames || !gmsh_entities ) {
        DGOUT(stderr,"\x1B[1;31m[GMSH] File incomplete! Check existence of all sections.\x1B[0m\n");
        return gridfile_read_failure(conn);
    }

    /* Element faces count*/
    if (num_nodes == 0 || num_elements == 0) {
        DGOUT(stderr,"\x1B[1;31m[GMSH] No elements or nodes found in mesh file.!\x1B[0m\n");
        return gridfile_read_failure(conn);
    }

    /* save element faces count */
    geom->nelem_face = num_faces;

    /* return the periodic vertex to vertex map (may be NULL) */
    return conn;
}

int gridfile_read_highorder(p4est_topidx_t which_tree,Real *xyz){
    const int sim_dim = DIM;
    char supported_vol_elem_type,supported_face_elem_type;
    size_t nelem_coord_bytes,elementType,element_byte_start;
    size_t header[2];
    int nvert,v,z_index;
    int qdegree;
    size_t ret;

    /* open element node file */
    FILE *fp = fopen("QNODES.bin","rb");

    /* read file header:
     *      nelem_coord_bytes(size_t),
     *      elementType(size_t)
     */
    ret = fread(header,sizeof(size_t),2,fp);
    nelem_coord_bytes = header[0]; /* 3*nvert*sizeof(double) */
    elementType = header[1];

    /* get number of vertices from element type */
    msh_supported_elem_type((int)elementType,sim_dim,
                            &supported_vol_elem_type,
                            &supported_face_elem_type,
                            &nvert);

    /* jump to this tree's coordinates */
    element_byte_start = nelem_coord_bytes*which_tree;
    fseek(fp,element_byte_start,SEEK_CUR);

    /* read this elements coordinates in bytes */
    double coordinates[3*nvert];
    ret = fread(coordinates,sizeof(double),3*nvert,fp);

    /* retrieve z-order node map */
    Uint map[nvert];
    msh_elem_type_to_qdegree((int)elementType,&qdegree);
    msh_zorder_node_numbering(sim_dim,qdegree,map);

    /* reorder coordinates for z-ordering */
    for (v = 0; v < nvert; ++v) {
        z_index = map[v];

        double *zco = &coordinates[3*v];
        xyz[3*z_index+0] = (Real) zco[0];
        xyz[3*z_index+1] = (Real) zco[1];
        xyz[3*z_index+2] = (Real) zco[2];
    }

    /* close file */
    fclose(fp);
    return qdegree;
}

void gridfile_allocate_data(gridfile_t *gfile,p4est_topidx_t num_trees){
    const int ntrees = num_trees;
    int i;

    /* allocate memory */
    gfile->elem_tag.malloc(ntrees);
    gfile->p4est_part.malloc(ntrees);

    for(i = 0; i < ntrees; i++) gfile->elem_tag[i] = 0;
    for(i = 0; i < ntrees; i++) gfile->p4est_part[i] = i; // identity map
}

void gridfile_deallocate_data(ctx_t *ctx){
    /* deallocate memory (only for unstructured grids) */
    if (ctx->d_simulation.unstructured_flag) {
        // delete QNODES file
        MPI_Barrier(ctx->comm);
        if(ctx->rank == 0){
            if(remove("QNODES.bin")) printf("Failed to remove QNODES.bin file!\n");
        }
    }
}

void gridfile_broadcast_data(ctx_t *ctx,p4est_topidx_t num_trees){
    gridfile_t *gfile = &ctx->d_gridfile;

    /* allocate gridfile data */
    if(ctx->rank != 0) gridfile_allocate_data(gfile,num_trees);

    MPI_Bcast(gfile->elem_tag.ptr(),  num_trees,MPI_UNSIGNED_LONG_LONG,0,ctx->comm);
    MPI_Bcast(gfile->p4est_part.ptr(),num_trees,MPI_UNSIGNED_LONG_LONG,0,ctx->comm);
}

char check_periodic_inputs(ctx_t *ctx){
    simulation_t *sim = &ctx->d_simulation;
    geometry_t *geom = &ctx->d_geometry;
    int ipatch;
    int i;

    for (i = 0; i < MAX_NBC; i++) {
        for (ipatch = 0; ipatch < geom->npatch; ipatch++) {
            char type_name[BUFF_SIZE] = { '\0' };
            sprintf(type_name,"periodic%d:",i);

            char notfound = amr_utilities_find_keyword(sim->input_file,type_name);
            if(!notfound) return 1;
        }
    }
    /* no periodic inputs found */
    return 0;
}

static
char neareps(const Real x,const Real y){
    return (fabs(x-y) < 1.0E-6);
}

void gridfile_periodic_boundary_construction(ctx_t *ctx,
                                             p4est_connectivity_t *conn,
                                             p4est_topidx_t *pv2v){
    simulation_t *sim = &ctx->d_simulation;
    geometry_t *geom = &ctx->d_geometry;

    int nnodes = P4EST_CHILDREN/2;
    Uint periodic_pair_tag[MAX_NBC][2];
    Uint vi[nnodes];
    Uint vj[nnodes];
    char nperiodic_bc;
    char found_pair;
    char pair_count;
    char notfound;
    int patch_id;
    int pair_id;
    int ipatch;
    int i,j,ni,nj;

    char *node_flag = (char *) calloc(conn->num_vertices,sizeof(char));

    /* count and gather periodic pairs */
    nperiodic_bc = 0;
    for (i = 0; i < MAX_NBC; i++) {
        found_pair = pair_count = 0;
        for (ipatch = 0; ipatch < geom->npatch; ipatch++) {
            char type_name[BUFF_SIZE] = { '\0' };
            sprintf(type_name,"periodic%d:",i);

            char patch_name[BUFF_SIZE]   = { '\0' };
            notfound = amr_utilities_find_keyword_string_optional(sim->input_file,1,type_name,patch_name);
            if (!notfound && strstr(patch_name,geom->patch_name[ipatch])) {
                periodic_pair_tag[nperiodic_bc][pair_count++] = ipatch;
                found_pair = 1;
            }
        }
        if(found_pair) nperiodic_bc++;
    }

    /* assemble periodic pair nodes */
    for (pair_id = 0; pair_id < nperiodic_bc; pair_id++) {
        int nnodes_match = 0;
        for (i = 0; i < geom->nelem_face; i++) {
            Uint *FNI_i = &geom->face_nodes_info[(2 + nnodes)*i];
            Uint entity_tag_i = FNI_i[EFENT_IND];

            /* get entity tag index */
            for (ipatch = 0; ipatch < geom->nentity_face; ipatch++) {
                if (geom->entity_face_tag[ipatch] == entity_tag_i) {
                    patch_id = ipatch;
                    break;
                }
            }

            /* check if main face */
            if(geom->entity_face_patchtag[patch_id] != periodic_pair_tag[pair_id][0]) continue; // skip face

            /* fill i-index nodes */
            for(ni = 0; ni < nnodes; ni++) vi[ni] = FNI_i[EFNDS_IND+ni];

            for (j = 0; j < geom->nelem_face; j++) {
                Uint *FNI_j = &geom->face_nodes_info[(2 + nnodes)*j];
                Uint entity_tag_j = FNI_j[EFENT_IND];

                /* get entity tag index */
                for (ipatch = 0; ipatch < geom->nentity_face; ipatch++) {
                    if (geom->entity_face_tag[ipatch] == entity_tag_j) {
                        patch_id = ipatch;
                        break;
                    }
                }

                /* check if secondary face */
                if(geom->entity_face_patchtag[patch_id] != periodic_pair_tag[pair_id][1]) continue; // skip face

                /* fill j-index nodes */
                for(nj = 0; nj < nnodes; nj++) vj[nj] = FNI_j[EFNDS_IND+nj];

                /* match i-index nodes to j-index nodes */
                char found_match = 0;
                char inode_found = 0;
                for (ni = 0; (ni < nnodes) && (inode_found==0); ni++) {
                    Real xi = conn->vertices[3*vi[ni]+0];
                    Real yi = conn->vertices[3*vi[ni]+1];
                E3D(Real zi = conn->vertices[3*vi[ni]+2])

                    /* search adjacent nodes that haven't been paired */
                    for (nj = 0; (nj < nnodes) && (!node_flag[vj[nj]]); nj++) {
                        Real xj = conn->vertices[3*vj[nj]+0];
                        Real yj = conn->vertices[3*vj[nj]+1];
                    E3D(Real zj = conn->vertices[3*vj[nj]+2])

                        /* compare node coordinates */
                        if(neareps(xi,xj)) found_match++;
                        if(neareps(yi,yj)) found_match++;
                    E3D(if(neareps(zi,zj)) found_match++)

                        /* if (DIM-1) nodes match, then pair face node found */
                        if (found_match  == (DIM-1)) {
                            node_flag[vj[nj]] = 1;  // set flag
                            pv2v[vj[nj]] = vi[ni];  // set paired node
#if 0
                            printf("PER NODE FOUND!      Main Face %d, node(%d)[%d]: %f %f, %s\n"
                                   "                Secondary Face %d, node(%d)[%d]: %f %f, %s\n",
                                i,ni,vi[ni],xi,yi,geom->patch_name[periodic_pair_tag[pair_id][0]],
                                j,nj,vj[nj],xj,yj,geom->patch_name[periodic_pair_tag[pair_id][1]]);
#endif
                            inode_found = 1;
                            nnodes_match++;
                            break; // skip remaining nj nodes
                        }
                        /* reset for post nj loop check */
                        found_match = 0;
                    }

                    /* if no match was found for node ni on face i, skip remaining nj nodes */
                    if(!found_match) break;
                }
            }
        }
        DGOUT(ctx->log_io,"\x1B[1;36m[GMSH] Periodic Face Info: %d matched nodes: [%s] [%s]\x1B[0m\n",
                          nnodes_match,
                          geom->patch_name[periodic_pair_tag[pair_id][0]],
                          geom->patch_name[periodic_pair_tag[pair_id][1]]);
    }

    /* deallocate memory */
    free(node_flag);
}