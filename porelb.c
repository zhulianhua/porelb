#include <mpi.h>
#include <stdio.h> 
#include <stdlib.h>
#include "porelb.h"

double feq(int k, double rho, double u, double v, double w)
{
    double uv,eu,x;
    eu=e[k][0]*u+e[k][1]*v+e[k][2]*w;
    uv=u*u+v*v+w*w;
    x=tp[k]*rho*(1.0+3*eu+4.5*eu*eu-1.5*uv ) ;
    return x;
}

void alloc_memory(SIMULATION *sim_)
{
    FILE* fp_map;
    if((fp_map = fopen("V.bin", "rb")) == NULL)
    {
        fprintf(stderr, "Node map file openning error!\n");
        exit(-1);
    }

    fread(&sim_->num_f, sizeof(int), 1, fp_map);
    if(!sim_->pid) printf("[0]>>> Total num of fluid  nodes = %d\n", sim_->num_f);

    int num_p;
    fread(&num_p,  sizeof(int), 1, fp_map);
    if(!sim_->pid){
        if(sim_->num_p != num_p) {
            fprintf(stderr, "[0]>>> MPI process num %d differs from num_p %d in V\n", sim_->num_p, num_p);
            exit(-2);
        }
    }

    fread(&sim_->nx, sizeof(int), 1, fp_map);
    fread(&sim_->ny, sizeof(int), 1, fp_map);
    fread(&sim_->nz, sizeof(int), 1, fp_map);

    int * start_loc = (int* ) calloc(sim_->num_p, sizeof(int));
    fread(start_loc, sizeof(int), sim_->num_p, fp_map);
    int this_start = start_loc[sim_->pid];
    int next_start = 0;
    if(sim_->pid < sim_->num_p -1) 
        next_start = start_loc[sim_->pid+1];
    else next_start = sim_->num_f;
    free(start_loc);

    //Read the start location of length table processed by me
    sim_->N = next_start - this_start;
    fprintf(stderr, "%d >>> My num of fluid nodes: %d \n", sim_->pid, sim_->N);

    //Allocate memory of V
    sim_->V = (NODE_INFO *)calloc(sim_->N, sizeof(NODE_INFO));
    if( sim_->V == NULL)
    {
        fprintf(stderr, "%d >>> Error calloc memory for V \n", sim_->pid);
        exit(-4);
    }

    //Seek the to the start of my NODE_INFO section
    fseek(fp_map, this_start*sizeof(NODE_INFO), SEEK_CUR);

    //Read my NODE_INFO array, i.e., V.
    if((fread(sim_->V, sizeof(NODE_INFO), sim_->N, fp_map)) != sim_->N)
    {
        fprintf(stderr, "%d >>> Error read V \n", sim_->pid);
        exit(-5);
    }

    /*printf("%d >>> DNOE Read NODE_INFO \n", sim_->pid);*/
    fclose(fp_map);

    //Allocate other data arrays
    sim_->f0   = (double *)calloc(sim_->N*Q, sizeof(double)); 
    sim_->f1   = (double *)calloc(sim_->N*Q, sizeof(double)); 
    sim_->Rho  = (double *)calloc(sim_->N  , sizeof(double)); 
    sim_->Ux   = (double *)calloc(sim_->N  , sizeof(double)); 
    sim_->Uy   = (double *)calloc(sim_->N  , sizeof(double)); 
    sim_->Uz   = (double *)calloc(sim_->N  , sizeof(double)); 
}

void init_para(SIMULATION* sim_)
{
    sim_->dt = 1.0;
    sim_->dx = 1.0;
    sim_->tau = 0.8;
    sim_->omega = 1.0/sim_->tau;
    sim_->nu = 1.0/3.0*(sim_->tau - 0.5);
    sim_->Xlength = sim_->nx;
    sim_->Xforce = 1.0e-5;
    sim_->it_print = 10000;
    sim_->it_save_data = 10000;
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
    free(sim_->N_c_n);
    free(sim_->send_request);
    free(sim_->recv_request);
    free(sim_->tag);
    free(sim_->Send_info);
    free(sim_->Recv_info);
    free(sim_->Send_Dist);
    free(sim_->Recv_Dist);
}

void check_V(SIMULATION *sim_)
{
    printf("%d >>> Start checking V \n", sim_->pid);
    for(int i=0; i<sim_->N; i++)
    {
        if( sim_->V[i].x < 0 || sim_->V[i].x > sim_->nx
          ||sim_->V[i].y < 0 || sim_->V[i].y > sim_->ny
          ||sim_->V[i].z < 0 || sim_->V[i].z > sim_->nz
          ) {
            fprintf(stderr, "%d >>> Check V, bound error!\n", sim_->pid); 
            exit(-3);
        }

        for(int j=0; j<Q-1; j++)
            if(sim_->V[i].nb_info[j].node_type == TYPE_EXTERN)
            {
                if( sim_->V[i].nb_info[j].node_pid > sim_->num_p-1
                  ||sim_->V[i].nb_info[j].node_pid < 0
                  ||sim_->V[i].nb_info[j].node_pid == sim_->pid)
                {
                    fprintf(stderr, "%d >>> Check V, pe error!\n", sim_->pid); 
                    exit(-4);
                }
            }
    }
    printf("%d >>> Check V past\n", sim_->pid);
}

void build_send_recv_info(SIMULATION *sim_)
{
    NODE_INFO *V = sim_->V;
    int N = sim_->N;
    LINK *links;
    int count = 0;

    sim_->N_c_n = (int*) calloc(sim_->num_p, sizeof(int));

    for(int i = 0; i<N; i++)
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
        if(sim_->N_c_n[i] > 0) sim_->N_p++;

    printf("%d >>> N_p = %d \n", sim_->pid, sim_->N_p);
    printf("%d >>> N_c = %d \n", sim_->pid, sim_->N_c);

    //Allocate MPI requsest and tag arrays
    sim_->send_request = (MPI_Request *)calloc(sim_->N_p, sizeof(MPI_Request));
    sim_->recv_request = (MPI_Request *)calloc(sim_->N_p, sizeof(MPI_Request));
    sim_->tag = (int *)calloc(sim_->N_p, sizeof(int));

    //Allocate additional data array
    links = (LINK *) calloc(count, sizeof(LINK));
    sim_->Send_info = (int *) calloc(2*count, sizeof(int));
    sim_->Recv_info = (int *) calloc(2*count, sizeof(int));
    sim_->Send_Dist = (double *) calloc(count, sizeof(double));
    sim_->Recv_Dist = (double *) calloc(count, sizeof(double));

    count = 0;
    for(int i=0; i<N; i++)
        for(int j=1; j<Q; j++)
            if(V[i].nb_info[j-1].node_type == TYPE_EXTERN)
            {
                links[count].node_pid = V[i].nb_info[j-1].node_pid;
                links[count].j = j;
                links[count].nid = V[i].nb_info[j-1].node_id;
                links[count].own_nid = i;
                count++;
            }
    // sorting
    qsort(links, count, sizeof(LINK), compare);

    // Build Send_Recv info
    for(int i=0; i<count; i++)
    {
        sim_->Send_info[2*i  ] = links[i].nid;
        sim_->Send_info[2*i+1] = links[i].j;
    }

    /*if(sim_->pid)*/
    /*{*/
        /*printf("1 >>> Send_info[0], nid = %d, j = %d\n", sim_->Send_info[0], sim_->Send_info[1]);*/
    /*}*/

    int start = 0;
    MPI_Status status;
    for(int i=0; i<sim_->num_p; i++)
    {
        if(sim_->N_c_n[i]) //have connection with processor i
        {
            MPI_Sendrecv(
                    sim_->Send_info+start, sim_->N_c_n[i]*2, MPI_INT, i, 1,
                    sim_->Recv_info+start, sim_->N_c_n[i]*2, MPI_INT, i, 1,
                    MPI_COMM_WORLD, &status);
            start+=2*sim_->N_c_n[i];
        }
    }

    for(int i=0; i<count; i++)
    {
        sim_->Send_info[2*i  ] = links[i].own_nid;
        sim_->Send_info[2*i+1] = links[i].j;
    }
}

void init_dist_macro(SIMULATION* sim_)
{
    for(int i=0; i<sim_->N; i++)
    {
        sim_->Rho[i] = 1.0;
        sim_->Ux[i] = 0.0;
        sim_->Uy[i] = 0.0;
        sim_->Uz[i] = 0.0;
        for(int j=0; j<Q; j++)
        {
            sim_->f0[i*Q+j] = sim_->f1[i*Q+j] = feq(j,sim_->Rho[i], 
                    sim_->Ux[i], sim_->Uy[i], sim_->Uz[i]);
        }
    }
}

void update_send_dist(SIMULATION* sim_)
{
    //Pack into Send_Dist
    int nid, j;
    for(int i=0; i<sim_->N_c; i++)
    {
        nid = sim_->Send_info[2*i];
        j   = sim_->Send_info[2*i+1];
        sim_->Send_Dist[i] = sim_->f0[nid*Q+j];
    }
}

void send_recv(SIMULATION* sim_)
{
    //Send
    int start = 0;
    int ireq = 0;
    for(int i=0; i<sim_->num_p; i++)
    {
        if(sim_->N_c_n[i])
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
        if(sim_->N_c_n[i])
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
    for(int i=0; i<sim_->N; i++)
    {
        sim_->f1[i*Q] = sim_->f0[i*Q];
        for(int j=1; j<Q; j++)
        {
            if(sim_->V[i].nb_info[j-1].node_type == TYPE_FLUID)
                sim_->f1[i*Q+j] = 
                    sim_->f0[sim_->V[i].nb_info[j-1].node_id*Q+j]; //normal streaming
            else if(sim_->V[i].nb_info[j-1].node_type == TYPE_SOLID)
                sim_->f1[i*Q+j] = sim_->f0[i*Q+re[j]]; //bounce back
        }
    }
}

void update_recv_dist(SIMULATION* sim_)
{
    //Unpack from Recv_Dist
    int nid, j;
    for(int i=0; i<sim_->N_c; i++)
    {
        nid = sim_->Recv_info[2*i];
        j   = sim_->Recv_info[2*i+1];
        sim_->f1[nid*Q+j] = sim_->Recv_Dist[i];
    }
}

void check_recv_dist(SIMULATION* sim_)
{
    for(int i=0; i<sim_->N_c; i++)
    {
        printf("N_c = %8d, f = %21.13e\n", i, sim_->Recv_Dist[i]);
        if(sim_->Recv_Dist[i] < 1e-5) {fprintf(stderr, "Err in recv_dist\n"); exit(-8);}
    }
}

void collision(SIMULATION* sim_)
{
    double rho, rhou, rhov, rhow, u, v, w;
    double omega = sim_->omega;

    for(int i=0; i<sim_->N; i++)
    {
        //Update macros
        rho  = 0.0;
        rhou = 0.0;
        rhov = 0.0;
        rhow = 0.0;
        for(int j=0; j<Q; j++)
        {
            rho  += sim_->f1[i*Q+j];
            rhou += e[j][0]*sim_->f1[i*Q+j];
            rhov += e[j][1]*sim_->f1[i*Q+j];
            rhow += e[j][2]*sim_->f1[i*Q+j];
        }

        if(rho > 1.2 | rho < 0.5) 
        {
            fprintf(stderr, "%d >>> at %d (%d, %d, %d),  Macro denstity  = %g , < 0.5 !\n", sim_->pid, i, sim_->V[i].x, sim_->V[i].y, sim_->V[i].z, rho);
            exit(-6);
        }


        //Correct velocity, only in X direction
        //Assuming dt = 1;
        if(sim_->it>0)
            rhou += 0.5*sim_->Xforce;

        u = rhou/rho;
        v = rhov/rho;
        w = rhow/rho;
        sim_->Rho[i] = rho;
        sim_->Ux[i] = u;
        sim_->Uy[i] = v;
        sim_->Uz[i] = w;

        //Colision, assuming dt = 1; cs^2 = 1/3
        for(int j=0; j<Q; j++)
        {
            sim_->f0[i*Q+j] = (1.0-omega)*sim_->f1[i*Q+j]
                + ( omega + (1.0-0.5*omega)*3.0*(e[j][0] - u)*sim_->Xforce )*feq(j,rho,u, v, w);
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
    char filename[100];
    FILE* fp;
    sprintf(filename, "part_%03d_%07d.bin", sim_->pid, sim_->it);
    fp = fopen(filename, "wb");
    if(fp == NULL)
    {
        fprintf(stderr, "Macro binary file open error!\n");
        exit(-3);
    }
    fwrite(sim_->Rho, sizeof(double), sim_->N, fp);
    fwrite(sim_->Ux , sizeof(double), sim_->N, fp);
    fwrite(sim_->Uy , sizeof(double), sim_->N, fp);
    fwrite(sim_->Uz , sizeof(double), sim_->N, fp);
    fclose(fp);
}

int compare(const void * a, const void * b)
{
    LINK *la = (LINK *)a;
    LINK *lb = (LINK *)b;
    if     ( la->node_pid > lb->node_pid ) return 1;
    else if( la->node_pid < lb->node_pid ) return 0;
    else  return (la->nid > lb->nid);
}

int main(int argc, char* argv[])
{
    SIMULATION sim;
    MPI_Init(&argc, &argv);
    if(argc != 2)
    {
        fprintf(stderr, ">>> Error argc, usage ./porelb.x it_max\n");
        exit(-10);
    }
    sim.it_max = atoi(argv[1]);


    MPI_Comm_rank(MPI_COMM_WORLD, &sim.pid);
    MPI_Comm_size(MPI_COMM_WORLD, &sim.num_p);

    alloc_memory(&sim);
    check_V(&sim);
    init_para(&sim);
    build_send_recv_info(&sim);
    init_dist_macro(&sim);

    sim.it_print = 1;
    double t1, t2;
    t1 = MPI_Wtime(); 
    for(sim.it =  0; sim.it<sim.it_max; sim.it++)
    {
        update_send_dist(&sim);
        send_recv(&sim);
        stream_inner_nodes(&sim);
        MPI_Waitall(sim.N_p, sim.send_request, MPI_STATUS_IGNORE);
        MPI_Waitall(sim.N_p, sim.recv_request, MPI_STATUS_IGNORE);
        update_recv_dist(&sim);
        collision(&sim);
        if(!sim.pid) printf(">>> Step : %20d \n", sim.it);
        if(sim.it % sim.it_save_data == 0) save_data(&sim);
    }
    t2 = MPI_Wtime(); 
    if(!sim.pid)printf( "Elapsed time is %f\n", t2 - t1 ); 
    save_data(&sim);
    MPI_Finalize();
    return 0;
}
