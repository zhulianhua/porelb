#include <math.h>
#include <stdlib.h>

#define GID(z,y,x) (((z*ny+y)*nx + x))
#define  Q 19

const int e[Q][3] = {
    { 0,  0,  0},//0

    { 1,  0,  0},//1 
    {-1,  0,  0},//2
    { 0,  1,  0},//3
    { 0, -1,  0},//4
    { 0,  0,  1},//5
    { 0,  0, -1},//6

    { 0,  1,  1},
    { 0, -1, -1},
    { 0, -1,  1},
    { 0,  1, -1},

    {-1,  0, -1},
    { 1,  0,  1},
    {-1,  0,  1},
    { 1,  0, -1},

    {-1,  1,  0},
    { 1, -1,  0},
    {-1, -1,  0},
    { 1,  1,  0}
};

const float w[Q] = {
    1.f/3, 
    1.f/18, 1.f/18, 1.f/18, 1.f/18, 1.f/18, 1.f/18,
    1.f/36, 1.f/36, 1.f/36, 1.f/36,
    1.f/36, 1.f/36, 1.f/36, 1.f/36, 
    1.f/36, 1.f/36, 1.f/36, 1.f/36
};


enum NEIGHBOR_TYPE {TYPE_FLUID, TYPE_SOLID, TYPE_EXTERN};

/* struct NEIGHBOR_INFO */
typedef struct {
    unsigned char node_type;
    unsigned int node_id;
    unsigned short node_pid;
} NEIGHBOR_INFO;

/* struct NODE_INFO */
typedef struct {
    unsigned short x, y, z;
    NEIGHBOR_INFO[Q-1];
    double fIn[Q];
    double fOut[Q];
} NODE_INFO;

/* struct PE_DATA */
typedef struct {
    NODE_INFO *V_sub;  //subdomain node head pointer
    unsigned short N_p; //Num. of neighbor PE;
    unsigned int N_c; //Num. of link with neighbor PE, total;
    unsigned int *N_c_n; //Num of link with each neighbor PE, integer array, dimension: N_p

    unsigned int **Send_info;
    unsigned int *Send_info_j;
    unsigned int *Send_info_nid;

    unsigned int **Recv_info;
    unsigned int *Recv_info_j;
    unsigned int *Recv_info_nid;

    double *Send_Dist;
    double *Recv_Dist;
} PE_DATA;

/* struct PARAMETERS */
typedef struct {
    unsigned int it_max, it_print, it_save_data, it_save_recovery;
    unsigned int nx, ny, nz;
    unsigned int size;
    unsigned int num_f;
    double dt, dx, tau, omega, nu;
    double L, a, Fx; //new
} PARAMETERS;

/* struct SIMULATION */
typedef struct {
    PE_DATA data;
    unsigned char flag;
    PARAMETERS param;
} SIMULATION;


/* Init simulatoin data and parameters */

/* set parameters */
void set_para(PARAMETER* para_);

/* 
 * compute global physics and LB constants based on inpurt parameters in 
*/
void init_const(PARAMETER* para_);

/*  
 * executed by master process,
 * build global Node_info array from flag data, and saved the array in a file
*/
void build_V(SIMULATION* simulation_);

/* 
 * executed by all process
 * Read precoessing output file, and extract sub node_info array into 
*/
void read_V_sub(SIMULATION* simulation_);


/*
 * exceted by all process
 * Sent out Send_info to neighbor's Recv_info
 */
void build_Send_Recv_info(SIMULATION* simulation_);

/*
 * Init Distributions
 */
void init_Dist(SIMULATION* simulation_);

/* LB evolution function */
void updata_Send_Dist(SIMULATION* simulation_);

void stream_Inner_nodes(SIMULATION* simulation_);

void update_Recv_Dist(SIMULATION* simulation_);

void collision(SIMULATION* simulation_);

/* Save and recovery */
void save_State(SIMULATION* simulation_);

void recovery_State(SIMULATION* simulation_);

void save_Data();
