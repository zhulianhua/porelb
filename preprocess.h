#ifndef __PREPROCESS_H_
#define __PREPROCESS_H_


#include "porelb.h"

#define GID_(z,y,x) (((z*ny_+y)*nx_ + x))


typedef struct {
    int num_node;  //num. of fluids nodes
    int num_p;     //num of process blocks
    int nx, ny, nz;
    int* V_sub_start;  //Store num. of nodes in each process blocks, length: num_p;
    NODE_INFO* node_info; //Array of fluids node info
} V_DESC;  //Sparse matrix descriptor

typedef struct {
    int pid;      //node's  local process id
    int pos_in_p; //node's index in local process
    int x, y, z;  //node's full matrix index
} INDEX_NODE;      //used in build_V();


//---------------------------------------------------------------------------------------
//Allocating dynamic memory for node map, Called by build_V()
//---------------------------------------------------------------------------------------
void init_v_desc(V_DESC *v_desc_, int num_p_, int num_node_);

//---------------------------------------------------------------------------------------
//Free node map memory after building and saving node map to file
//---------------------------------------------------------------------------------------
void free_v_desc(V_DESC *v_desc_);

//---------------------------------------------------------------------------------------
//build node map from flag file, nodes' order in the node array is raw(not sorted by pid)
//---------------------------------------------------------------------------------------
V_DESC build_V(unsigned short num_p_, int nx_, int ny_, int nz_, char* flag_file_name_);

//---------------------------------------------------------------------------------------
//save node map into binary file, descripted by struct V_DESC
//---------------------------------------------------------------------------------------
void save_V(V_DESC *v_desc_, char* file_name_);

#endif
