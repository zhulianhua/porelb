#include <stdlib.h>
#include <stdio.h>

#include "mytype.h"
#include "node.h"
#include "preprocess.h"

void init_v_desc(V_DESC *v_desc_, int num_p_, int num_node_)
{
    v_desc_->num_p = num_p_;
    v_desc_->num_node = num_node_;
    if((v_desc_->V_sub_start = (int* )malloc(sizeof(int)*num_p_))==NULL){
        fprintf(stderr, "Memory allocating failed!\n");
        exit(-3);
    };

    if((v_desc_->node_info = (NODE_INFO*)malloc(sizeof(NODE_INFO)*num_node_))==NULL){
        fprintf(stderr, "Memory allocating failed!\n");
        exit(-3);
    }
}

void free_v_desc(V_DESC *v_desc_)
{
    free(v_desc_->V_sub_start);
    free(v_desc_->node_info);
}

V_DESC build_V(unsigned short num_p_, int nx_, int ny_, int nz_, char* flag_file_name_)
{
    FILE* fp_flag;
    unsigned char* flag;

    int* X_to_n;  //free
    INDEX_NODE* index_node; //free
    int* count_node_of_p; //free

    int x, y, z, k, xp, yp, zp, kp, i;
    int num_node;
    V_DESC v_desc;

    int size;

    size = nx_*ny_*nz_;

    if((flag = (unsigned char*)malloc(size*sizeof(unsigned char))) == NULL) {
        fprintf(stderr, "Memory allocating failed!\n");
        exit(-3);
    }

    if((X_to_n = (int* )calloc(size, sizeof(int))) == NULL){
        fprintf(stderr, "Memory allocating failed!\n");
        exit(-3);
    }

    if((count_node_of_p = (int*)calloc(num_p_, sizeof(int)))==NULL){
        fprintf(stderr, "Memory allocating failed!\n");
        exit(-3);
    }

    if((fp_flag=fopen(flag_file_name_, "rb")) == NULL){
        fprintf(stderr, "flag data file opening error!\n");
        exit(-1);
    }

    if(fread(flag, sizeof(unsigned char), size, fp_flag ) != size){
        fprintf(stderr, "flag data file reading error!\n");
        exit(-2);
    }


    num_node = 0;
    //counting the fluid nodes
    for(z=0; z<nz_; z++)
        for(y=0; y<ny_; y++)
            for(x=0; x<nx_; x++){
                k = GID(z, y, x);
                if(flag[k] == 0){ //fluid nodes
                    X_to_n[k] = num_node;   //save the index from X to  global fluid node id
                    num_node++;
                }
            }

    v_desc.nx = nx_;
    v_desc.ny = ny_;
    v_desc.nz = nz_;
    init_v_desc(&v_desc, num_p_, num_node); //build v_desc, allocating memory....

    if((index_node = (INDEX_NODE*)malloc(sizeof(INDEX_NODE)*num_node))==NULL) {
        fprintf(stderr, "Memory allocating failed!\n");
        exit(-3);
    }

    //index_node: store all of the fluids node info, in the order of their z, y, x index
    i = 0;
    INDEX_NODE index_node_tmp;
    for(z=0; z<nz_; z++)
        for(y=0; y<ny_; y++)
            for(x=0; x<nx_; x++){
                k = GID(z, y, x);
                if(flag[k] == 0){ //fluid nodes
                    index_node_tmp.x = x;
                    index_node_tmp.y = y;
                    index_node_tmp.z = z;
                    index_node_tmp.pid = z/(nz_/num_p_);  //Domain decomposition, here just slicing along z driection
                    index_node_tmp.pos_in_p = count_node_of_p[index_node_tmp.pid];
                    count_node_of_p[index_node_tmp.pid]++;  //counting the nodes in this process
                    index_node[i] = index_node_tmp;  //Assign struct to a struct item of an  struct array
                    i++;
                }
            }

    v_desc.V_sub_start[0] = 0;
    for(i=1; i<num_p_; i++)
        v_desc.V_sub_start[i] = v_desc.V_sub_start[i-1] + count_node_of_p[i-1];


    NODE_INFO node_info_tmp;
    int p;
    int ii = 0;
    int np,  pp; //previous node's node id, process id
    for(z=0; z<nz_; z++)
        for(y=0; y<ny_; y++)
            for(x=0; x<nx_; x++){
                k = GID(z, y, x);
                if(flag[k] == 0){ //fluid nodes
                    node_info_tmp.x = x;
                    node_info_tmp.y = y;
                    node_info_tmp.z = z;
                    p = index_node[ii].pid;
                    for(i=0; i<Q-1; i++){ //big bug here Q -> Q-1
                        //Optimize
                        /*xp = (x  + nx_)%(nx_);  //periodical*/
                        /*yp = (y  + ny_)%(ny_);  //periodical*/
                        /*zp = (z  + nz_)%(nz_);  //periodical*/
                        xp = (x + e[i+1][0] + nx_)%(nx_);  //periodical
                        yp = (y + e[i+1][1] + ny_)%(ny_);  //periodical
                        zp = (z + e[i+1][2] + nz_)%(nz_);  //periodical
                        //end Optimize
                        kp = GID(zp, yp, xp);
                        if(flag[kp] == 1) //previous node is solid node
                            node_info_tmp.nb_info[i].node_type = TYPE_SOLID;
                        else{
                            np = X_to_n[kp]; //previous node's fluid gid
                            pp = index_node[np].pid; //previous node's pid
                            node_info_tmp.nb_info[i].node_id = index_node[np].pos_in_p;
                            node_info_tmp.nb_info[i].node_pid = pp;
                            node_info_tmp.nb_info[i].node_type=(pp == p)?TYPE_FLUID:TYPE_EXTERN_FLUID;
                        }
                    }
                    //flishing build temp node_info, then insert to the appropriate V array
                    v_desc.node_info[v_desc.V_sub_start[p] + index_node[ii].pos_in_p] = node_info_tmp;
                    ii++;
                }
            }

    //free memory;
    free(index_node);
    free(count_node_of_p);
    free(X_to_n);
    return v_desc;
}

void save_V(V_DESC *v_desc_, char* file_name_)
{
    FILE* fp;
    if((fp=fopen(file_name_, "wb"))==NULL){
        fprintf(stderr, "File openning failed\n");
        exit(-1);
    }

    fwrite(&v_desc_->num_node, sizeof(int), 1, fp);
    fwrite(&v_desc_->num_p, sizeof(int), 1, fp);
    fwrite(&v_desc_->nx, sizeof(int), 1, fp);
    fwrite(&v_desc_->ny, sizeof(int), 1, fp);
    fwrite(&v_desc_->nz, sizeof(int), 1, fp);
    if((fwrite(v_desc_->V_sub_start, sizeof(int), v_desc_->num_p, fp)) != v_desc_->num_p) {
        fprintf(stderr, "File writing failed\n");
        fclose(fp);
        exit(-1);
    }

    if((fwrite(v_desc_->node_info, sizeof(NODE_INFO), v_desc_->num_node, fp)) != v_desc_->num_node) {
        fprintf(stderr, "File writing failed\n");
        fclose(fp);
        exit(-1);
    }
    fclose(fp);
}
