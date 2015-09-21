#include <mpi.h>
#include "porelb.h"

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

const double tp[Q] = {
    1.0/3, 
    1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/18,
    1.0/36, 1.0/36, 1.0/36, 1.0/36,
    1.0/36, 1.0/36, 1.0/36, 1.0/36, 
    1.0/36, 1.0/36, 1.0/36, 1.0/36
};

const int re[Q] = { 
    0,
    2, 1, 4, 3, 6, 5, 
    8, 7, 10, 9,
    12, 11, 14, 13,
    16, 15, 18, 17
};

double feq(int k, double rho, double u, double v, double w)
{
    double uv,eu,x;
    eu=e[k][0]*u+e[k][1]*v+e[k][2];
    uv=u*u+v*v+w*w;
    x=tp[k]*rho*(1.0+3*eu+4.5*eu*eu-1.5*uv ) ;
    return x;
}

void alloc_memory(SIMULATION *sim_)
{
    FILE *fp_map;
    if((fp_map = fopen(sim_->node_map_filename, "rb")) == NULL)
    {
        fprintf(stderr, "Node map file openning error!\n");
        exit(-1);
    }
    fread(&sim_->num_f, sizeof(int), 1, fp_map);
    if(!sim_->pid) printf("[0]>>> Total num of fluid  nodes = %d\n", sim_->num_f);

    int num_p;
    fread(&num_p,  sizeof(int), 1, fp_map);
    if(!sim_->pid){
        if(!sim_->num_p != num_p) {
            fprintf(stderr, "[0]>>> MPI process num %d differs from num_p %d in V", sim_->num_p, num_p);
            exit(-2);
        }
    }

    fread(&sim_->nx, sizeof(int), 1, fp_map);
    fread(&sim_->ny, sizeof(int), 1, fp_map);
    fread(&sim_->nz, sizeof(int), 1, fp_map);

    //Seek to current 
    fseek(fp_map, sim_->pid*sizeof(int), SEEK_CUR);
    int next_start;
    int this_start;
    //Read the start location of length table processed by me
    fread(&this_start, sizeof(int), 1, fp_map);
    //Read num of fluid process by next process
    fread(&next_start, sizeof(int), 1, fp_map);
    sim_->N = next_start - this_start;
    fprintf(stderr, "[%d]>>> My num of fluid nodes: %d ", sim_->pid, sim_->N);
    //Allocate memory of V
    V = (NODE_INFO *)calloc(sim_->N, sizeof(NODE_INFO));

    //Seek the remaining num_p - pid -1 integers
    fseek(fp_map, (sim_->num_p - sim_->pid -1)*sizeof(int), SEEK_CUR);
    //Seek the to the start of my NODE_INFO section
    fseek(fp_map, this_start*sizeof(NODE_INFO), SEEK_CUR);
    //Read my NODE_INFO array, i.e., V.
    fread(sim_->V, sizeof(NODE_INFO), sim_->N, fp_map);
    fclose(fp_map);

    //Allocate other data arrays
    sim_->fIn  = (double *)calloc(sim_->N*Q, sizeof(double)); 
    sim_->fOut = (double *)calloc(sim_->N*Q, sizeof(double)); 
    sim_->Rho  = (double *)calloc(sim_->N  , sizeof(double)); 
    sim_->Ux   = (double *)calloc(sim_->N  , sizeof(double)); 
    sim_->Uy   = (double *)calloc(sim_->N  , sizeof(double)); 
    sim_->Uz   = (double *)calloc(sim_->N  , sizeof(double)); 
}

void init_para(SIMULATION* sim_)
{

}

void free_memory(SIMULATION *sim_)
{
    free(sim_->V);
    free(sim_->f0);
    free(sim_->f1);
    free(sim_->Rho);
    free(sim_->Ux);
    free(sim_->Uy);
    free(sim_->Uz);
}

void build_send_recv_info(SIMULATION *sim_)
{
    NODE_INFO *V = sim_->V;
    int N = sim_->N;
    LINK *links;
    int count = 0;

    sim_->N_c_n = calloc(sim_->num_p, sizeof(int));

    for(int i=0; i<N; i++)
        for(int j=1; j<Q; j++)
            if(V[i].nb_info[j-1].node_type == TYPE_EXTERN)
            {
                count++;
                sim_->N_c_n[V[i].nb_info[j-1].node_pid]++;
            }

    sim_->N_c = count;

    // Count how many neightour PE
    sim_->N_p = 0;
    for(int i=0; i<sim_->num_p; i++)
        if(!N_c_n[i]) sim_->N_p++;

    //Allocate MPI requsest and tag arrays
    sim_->send_request = (MPI_Request *)calloc(sim_->N_p, sizeof(MPI_Request));
    sim_->recv_request = (MPI_Request *)calloc(sim_->N_p, sizeof(MPI_Request));
    sim_->tag = (int *)calloc(sim_->N_p, sizeof(int));

    //Allocate additional data array
    links = (LINK *)calloc(count, sizeof(LINK));
    Send_info = (int *) calloc(2*count, sizeof(int));
    Recv_info = (int *) calloc(2*count, sizeof(int));
    Send_Dist = (double *) calloc(count, sizeof(double));
    Recv_Dist = (double *) calloc(count, sizeof(double));

    count = 0;
    for(int i=0; i<N; i++)
        for(int j=1; j<Q; j++)
            if(V[i].nb_info[j-1].node_type == TYPE_EXTERN)
            {
                links[count].node_pid = V[i].nb_info[j-1].node_pid;
                links[count].j = j;
                links[count].nid = V[i].nb_info[j-1].node_nid;
            }
    // sorting
    qsort(links, count, sizeof(LINK), compare);

    //count num of link for each neighbour PE

    // Build Send_Recv info
    for(int i=0; i<count; i+=2)
    {
        Send_info[i      ] = links[i].nid;
        Send_info[i+1    ] = links[i].j;
    }
    int start = 0;
    MPI_Status status;
    for(int i=0; i<sim_->num_p; i++)
    {
        if(N_c_n[i]) //have connection with processor i
        {
            MPI_Sendrecv(
                    Send_info+start, N_c_n[i]*2, MPI_INT, i, 1,
                    Recv_info+start, N_c_n[i]*2, MPI_INT, i, 1
                    MPI_COMM_WORLD, &status);
            start+=2*N_c_n[i];
        }
    }
}

void init_dist_macro(SIMULATION* sim_)
{
    for(int i=0; i<sim_->N; i++)
    {
        sim_->rho[i] = 1.0;
        sim_->Ux[i] = 0.0;
        sim_->Uy[i] = 0.0;
        sim_->Uz[i] = 0.0;
        for(int j=0; j<Q; j++)
        {
            f0[i*Q+j] = f1[i*Q+j] = feq(j,sim_->rho[i], 
                    sim_->Ux[i], sim_->Ux[i], sim_->Uz[i])
        }
    }
}

void update_send_dist(SIMULATION* sim_)
{
    //Pack into Send_Dist
    int nid, j;
    for(int i=0; i<sim_->N_c; i++)
    {
        nid = sim_->Recv_info[i*2];
        j   = sim_->Recv_info[i*2+1];
        sim_->Send_Dist[i] = sim_->f0[nid*Q+re[j]];
    }
}

void send_recv(SIMULATION* sim_)
{
    //Send
    int start = 0;
    int ireq = 0;
    for(int i=0; i<sim_->num_p; i++)
    {
        if(!sim_->N_c_n[i])
        {
            MPI_Isend(sim_->Send_Dist+start, sim_->N_c_n[i], MPI_DOUBLE,
                    i, sim_->tag[ireq], MPI_COMM_WORLD, &sim_->send_request[ireq]);
            ireq++;
            start += sim_->N_c_n[i];
        }
    }
    //Recv
    start = 0;
    ireq = 0;
    for(int i=0; i<sim_->num_p; i++)
    {
        if(!sim_->N_c_n[i])
        {
            MPI_Irecv(sim_->Recv_Dist+start, sim_->N_c_n[i], MPI_DOUBLE,
                    i, sim_->tag[ireq], MPI_COMM_WORLD, &sim_->recv_request[ireq]);
            ireq++;
            start += sim_->N_c_n[i];
        }
    }


}

void stream_inner_nodes(SIMULATION* sim_)
{
    int nid, jj;
    for(int i=0; i<N; i++)
    {
        sim_->f1[i*Q] = sim_->f0[i*Q];
        for(int j=1; j<Q; j++)
        {
            if(V[i].nb_info[j-1].node_type == TYPE_FLUID)
                sim_->f1[i*Q+j] = 
                    sim_->f0[V[i].nb_info[j-1].node_id*Q+j]; //normal streaming
            else if(V[i].nb_info[j-1].node_type == TYPE_SOLID)
                sim_->f1[i*Q+j] = sim_->f0[i*Q+re[j]]; //bounce back
        }
    }
}

void update_recv_dist(SIMULATION* sim_)
{
    //Pack from Recv_Dist
    int nid, j;
    for(int i=0; i<sim_->N_c; i++)
    {
        nid = sim_->Recv_info[i*2];
        j   = sim_->Recv_info[i*2+1];
        sim_->f1[nid*Q+j] = sim_->Send_Dist[i];
    }
}

void collision(SIMULATION* sim_)
{
    int rho, rhou, rhov, rhow, u, v, w;
    double omega = sim_->omega;
    for(int i=0; i<N; i++)
    {
        //Update macros
        rho  = 0.0;
        rhou = 0.0;
        rhov = 0.0;
        rhow = 0.0;
        for(int j=0; j<Q; j++)
        {
            rho  += f1[i*Q+j];
            rhou += e[j][0]*f1[i*Q+j];
            rhov += e[j][1]*f1[i*Q+j];
            rhow += e[j][2]*f1[i*Q+j];
        }


        //Correct velocity, only in X direction
        //Assuming dt = 1;
        if(sim_->it>0)
            rhou += 0.5*sim_->Xforce;

        u = rhou/rho;
        v = rhov/rho;
        w = rhow/rho;
        sim_->rho[i] = rho;
        sim_->Ux[i] = u;
        sim_->Uy[i] = v;
        sim_->Uz[i] = w;

        //Colision, assuming dt = 1; cs^2 = 1/3
        for(int j=0; j<Q; j++)
        {
            sim_->f0[i*Q+j] = (1.0-omega)*sim_->f1[i*Q+j]
                + ( omega + (1.0-0.5*omega)*3.0*(e[j][0] - u)*Xforce )
                *feq(j,rho,u, v, w)
        }
    }
}

void save_state(SIMULATION* sim_)
{
}

void recovery_state(SIMULATION* sim_)
{
}

void save_data(SIMULATION* sim_)
{
    double *u, *v, *rho;
    unsigned char *flag;
    if(!sim_->pid)
    {
        u   = (double *)calloc(sim_->num_f, sizeof(double));
        v   = (double *)calloc(sim_->num_f, sizeof(double));
        w   = (double *)calloc(sim_->num_f, sizeof(double));
        rho = (double *)calloc(sim_->num_f, sizeof(double));

    }else
    {

    }



}

int compare(const void * a, const void * b)
{

    LINK *la = (LINK *)a;
    LINK *lb = (LINK *)b;
    return ( la->node_pid - orderA->node_pid);
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    SIMULATION sim;

    MPI_Comm_rank(MPI_COMM_WORLD, &sim.pid);
    MPI_Comm_size(MPI_COMM_WORLD, &sim.num_p);

    alloc_memory(&sim);
    init_para(&sim);
    build_send_recv_info(&sim);
    init_dist_macro(&sim);

    for(int i=0; i<it_max; i++)
    {
        update_dend_Dist(&sim);
        send_recv(&sim);
        stream_inner_nodes(&sim);
        MPI_Wateall(sim.N_p, send_request, MPI_STATUS_IGNORE);
        MPI_Wateall(sim.N_p, recv_request, MPI_STATUS_IGNORE);
        update_recv_dist(&sim);
        collision(&sim);
        if( i%it_print && !sim.pid) printf("Step : %20d \n", i)
    }

    MPI_Finalize();
}
