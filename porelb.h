#include <math.h>
#include <stdlib.h>

#define GID(z,y,x) (((z*ny+y)*nx + x))
#define  Q 19

extern const int e[Q][3];
extern const double tp[Q];
extern const int re[Q];

enum NEIGHBOR_TYPE {TYPE_FLUID, TYPE_SOLID, TYPE_EXTERN};

/* struct NEIGHBOR_INFO */
typedef struct {
    unsigned char node_type;
    unsigned int node_id;
    unsigned short node_pid;
} NEIGHBOR_INFO;

/* struct LINK */
typedef struct {
    unsigned short node_pid;
    int j;
    int nid;
} LINK; /* used only for sorting */

/* struct NODE_INFO */
typedef struct {
    unsigned short x, y, z;
    NEIGHBOR_INFO nb_info[Q-1];
} NODE_INFO;

/* struct PE_DATA */
typedef struct {
    int N_p;  //Num. of neighbor PE;
    int *nei_pid //Store the neighbor PID for each section of send_info, size: N_p
    int *N_c_n; //Num of link with each neighbor PE, integer array, size: num_p

    int N_c;    //Num. of link with neighbor PE, total;

    int *Send_info; //size: 2*N_c
    int *Recv_info; //size: 2*N_c

    double *Send_Dist;
    double *Recv_Dist;

} PE_DATA;

/* struct SIMULATION */
typedef struct {
    // Simulation parameters
    int it_max, it_print, it_save_data, it_save_recovery, it;
    int nx, ny, nz;
    int size;
    int num_f;
    double dt, dx, tau, omega, nu;
    double Xlength, Xforce; 
    char node_map_filename[100];


    //DATA exchange
    int N; //Num of fluid node processed by this processor, read from file
    PE_DATA pe_data;
    NODE_INFO *V;    //subdomain node head pointer, size N

    //MPI 
    int pid;
    int num_p;
    MPI_request *recv_request, *send_request;
    int *tag;

    // LBM data arrays
    double *f0;
    double *f1;
    double *Rho;
    double *Ux;
    double *Uy;
    double *Uz;
} SIMULATION;

/* Init simulatoin data and parameters */

/* set parameters */

/*
 * Executed by all processes
 * Allocating all data arrays
 */
void alloc_memory(SIMULATION* sim_);

/*
 * Executed by all processes
 * Free all data arrays memory
 */
void free_memory(SIMULATION* sim_);

/*
 * exceted by all process
 * Sent out Send_info to neighbor's Recv_info
 */
void build_send_recv_info(SIMULATION* sim_);

void init_para(SIMULATION* sim_);

void init_dist_macro(SIMULATION* sim_);

/* LB evolution function */
void update_send_dist(SIMULATION* sim_);

void send_recv(SIMULATION* sim_);

void stream_inner_nodes(SIMULATION* sim_);

void update_recv_dist(SIMULATION* sim_);

void collision(SIMULATION* sim_);

/* Save and recovery */
void save_state(SIMULATION* sim_);

void recovery_state(SIMULATION* sim_);

void save_data();

int compare (const void * a, const void * b);
