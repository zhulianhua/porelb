#ifndef __NODE_H_
#define __NODE_H_

#include "lb.h"

enum NEIGHBOR_TYPE {TYPE_FLUID, TYPE_SOLID, TYPE_EXTERN_FLUID};

/* struct NEIGHBOR_INFO */
typedef struct {
    unsigned short node_pid;
    unsigned char node_type;
    unsigned int node_id;
} NEIGHBOR_INFO;

/* struct NODE_INFO */
typedef struct {
    unsigned short x, y, z;
    NEIGHBOR_INFO nb_info[Q-1];
} NODE_INFO;

#endif
