#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "simgpu.h"
#include "node.h"
#include "lb.h"
#include "defs.h"
#include "mytype.h"
/*#include "helper_cuda.h"*/
/*#include "helper_functions.h"*/
#include "TECIO.h"

/*Block size*/
#define BX 128
#define BY 1
#define BZ 1

const int e_h[Q][3] = {
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
const MY_DTYPE w_h[Q] = {
    1.0/3, 
    1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/18,
    1.0/36, 1.0/36, 1.0/36, 1.0/36,
    1.0/36, 1.0/36, 1.0/36, 1.0/36, 
    1.0/36, 1.0/36, 1.0/36, 1.0/36
};

__constant__ int e_d[Q][3] = {
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
__constant__ MY_DTYPE w_d[Q] = {
    1.0/3, 
    1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/18,
    1.0/36, 1.0/36, 1.0/36, 1.0/36,
    1.0/36, 1.0/36, 1.0/36, 1.0/36, 
    1.0/36, 1.0/36, 1.0/36, 1.0/36
};


void setPara(SIMGPU *simgpu_,  int argc,  char* argv[])
{
    simgpu_->nuF = 1.0/6;
    simgpu_->nuG = 1.0/6;
    
    simgpu_->tauF = 3.0*simgpu_->nuF + 0.50;
    simgpu_->tauG = 3.0*simgpu_->nuG + 0.50;

    printf("tauF = %f, tauG =%f\n", simgpu_->tauF,  simgpu_->tauG);
    simgpu_->omegaF = 1.0/simgpu_->tauF;
    simgpu_->omegaG = 1.0/simgpu_->tauG;

    simgpu_->rho_f = 1.0;
    simgpu_->rho_g = 1.0;
    simgpu_->rhoS = 1.0;

    simgpu_->G_fg =  0.200;
    simgpu_->G_fs =  0.030;
    simgpu_->G_gs =  -0.030;

    simgpu_->dt = 1.00;
    simgpu_->dx = 1.00;
    simgpu_->Fx = 1.0e-7;

    sprintf(simgpu_->node_map_filename, "%s", "V.bin");

    if(argc !=4){
        fprintf(stderr, "Argc error\n");
        exit(-1);
    }
    simgpu_->itMax = atoi(argv[1]);
    simgpu_->NewOrContinue = argv[2][0];
    //////////////////////////////////////
    printf("G_fg = %f\n", simgpu_->G_fg);
}

void allocMemroy(SIMGPU *simgpu_)
{
    FILE *fp_map;
    if((fp_map = fopen(simgpu_->node_map_filename, "rb")) == NULL)
    {
        fprintf(stderr, "Node map file openning error!\n");
        exit(-1);
    }
    fread(&simgpu_->num_node, sizeof(int), 1, fp_map);
    printf("num_node = %d\n", simgpu_->num_node);
    int num_p;
    fread(&num_p, sizeof(int), 1, fp_map);
    fread(&simgpu_->nx, sizeof(int), 1, fp_map);
    fread(&simgpu_->ny, sizeof(int), 1, fp_map);
    fread(&simgpu_->nz, sizeof(int), 1, fp_map);
    fseek(fp_map, num_p*sizeof(int), SEEK_CUR);
    simgpu_->L = (MY_DTYPE)simgpu_->nx;
    fclose(fp_map);


    //allocating memory at host
    simgpu_->f0_h = (MY_DTYPE *)calloc(Q*simgpu_->num_node, sizeof(MY_DTYPE));
    simgpu_->g0_h = (MY_DTYPE *)calloc(Q*simgpu_->num_node, sizeof(MY_DTYPE));
    simgpu_->rhoF_h = (MY_DTYPE *)calloc(simgpu_->num_node, sizeof(MY_DTYPE));
    simgpu_->rhoG_h = (MY_DTYPE *)calloc(simgpu_->num_node, sizeof(MY_DTYPE));
    simgpu_->vxF_h = (MY_DTYPE *)calloc(simgpu_->num_node, sizeof(MY_DTYPE));
    simgpu_->vyF_h = (MY_DTYPE *)calloc(simgpu_->num_node, sizeof(MY_DTYPE));
    simgpu_->vzF_h = (MY_DTYPE *)calloc(simgpu_->num_node, sizeof(MY_DTYPE));
    simgpu_->vxG_h = (MY_DTYPE *)calloc(simgpu_->num_node, sizeof(MY_DTYPE));
    simgpu_->vyG_h = (MY_DTYPE *)calloc(simgpu_->num_node, sizeof(MY_DTYPE));
    simgpu_->vzG_h = (MY_DTYPE *)calloc(simgpu_->num_node, sizeof(MY_DTYPE));

    simgpu_->vx_h = (MY_DTYPE *)calloc(simgpu_->num_node, sizeof(MY_DTYPE));
    simgpu_->vy_h = (MY_DTYPE *)calloc(simgpu_->num_node, sizeof(MY_DTYPE));
    simgpu_->vz_h = (MY_DTYPE *)calloc(simgpu_->num_node, sizeof(MY_DTYPE));

    simgpu_->node_map_h = (unsigned int *)calloc(Q*simgpu_->num_node, sizeof(unsigned int));
    simgpu_->n_to_XYZ = (unsigned short *)calloc(3*simgpu_->num_node, sizeof(unsigned short));
    if(!( simgpu_->node_map_h && simgpu_->n_to_XYZ 
                && simgpu_->g0_h && simgpu_->f0_h
                && simgpu_->rhoF_h && simgpu_->rhoG_h )) {
        fprintf(stderr, "CPU memory allocating error\n");
        exit(-1);
    }

    //allocating memory at device
    cudaMalloc((void **)&simgpu_->node_map_d, simgpu_->num_node*Q*sizeof(unsigned int));

    cudaMalloc((void **)&simgpu_->f0_d, simgpu_->num_node*Q*sizeof(MY_DTYPE));
    cudaMalloc((void **)&simgpu_->f1_d, simgpu_->num_node*Q*sizeof(MY_DTYPE));
    cudaMalloc((void **)&simgpu_->g0_d, simgpu_->num_node*Q*sizeof(MY_DTYPE));
    cudaMalloc((void **)&simgpu_->g1_d, simgpu_->num_node*Q*sizeof(MY_DTYPE));


    cudaMalloc((void **)&simgpu_->rhoF_d, simgpu_->num_node*sizeof(MY_DTYPE));
    cudaMalloc((void **)&simgpu_->rhoG_d, simgpu_->num_node*sizeof(MY_DTYPE));

    cudaMalloc((void **)&simgpu_->vxF_d, simgpu_->num_node*sizeof(MY_DTYPE));
    cudaMalloc((void **)&simgpu_->vyF_d, simgpu_->num_node*sizeof(MY_DTYPE));
    cudaMalloc((void **)&simgpu_->vzF_d, simgpu_->num_node*sizeof(MY_DTYPE));
    cudaMalloc((void **)&simgpu_->vxG_d, simgpu_->num_node*sizeof(MY_DTYPE));
    cudaMalloc((void **)&simgpu_->vyG_d, simgpu_->num_node*sizeof(MY_DTYPE));
    cudaMalloc((void **)&simgpu_->vzG_d, simgpu_->num_node*sizeof(MY_DTYPE));
}

void freeMemrory(SIMGPU *simgpu_)
{
    free(simgpu_->node_map_h); 
    free(simgpu_->n_to_XYZ); 

    free(simgpu_->f0_h); 
    free(simgpu_->g0_h); 
    free(simgpu_->rhoF_h);
    free(simgpu_->rhoG_h);
    free(simgpu_->vxF_h);
    free(simgpu_->vyF_h);
    free(simgpu_->vzF_h);
    free(simgpu_->vxG_h);
    free(simgpu_->vyG_h);
    free(simgpu_->vzG_h);
    free(simgpu_->vx_h);
    free(simgpu_->vy_h);
    free(simgpu_->vz_h);

    cudaFree(simgpu_->node_map_d); 

    cudaFree(simgpu_->f0_d);
    cudaFree(simgpu_->f1_d);
    cudaFree(simgpu_->g0_d);
    cudaFree(simgpu_->g1_d);

    cudaFree(simgpu_->rhoF_d);
    cudaFree(simgpu_->rhoG_d);

    cudaFree(simgpu_->vxF_d);
    cudaFree(simgpu_->vyF_d);
    cudaFree(simgpu_->vzF_d);
    cudaFree(simgpu_->vxG_d);
    cudaFree(simgpu_->vyG_d);
    cudaFree(simgpu_->vzG_d);
}

void buildDataArray(int *it, SIMGPU *simgpu_)
{
    int num_node = simgpu_->num_node;
    int i,k;
    NODE_INFO node_info_tmp;
    FILE *fp_map;
    if((fp_map = fopen(simgpu_->node_map_filename, "rb")) == NULL) {
        fprintf(stderr, "Node map file openning error!\n");
        exit(-1);
    }
    fseek(fp_map, sizeof(int), SEEK_SET);
    int num_p;
    fread(&num_p, sizeof(int), 1, fp_map);
    fseek(fp_map, 3*sizeof(int), SEEK_CUR);
    fseek(fp_map, num_p*sizeof(int), SEEK_CUR);
    printf("num_p = %d\n", num_p);

    for(k=0; k<num_node; k++) {
        fread(&node_info_tmp,sizeof(NODE_INFO), 1 , fp_map);

        simgpu_->n_to_XYZ[0*num_node + k] = node_info_tmp.x;
        simgpu_->n_to_XYZ[1*num_node + k] = node_info_tmp.y;
        simgpu_->n_to_XYZ[2*num_node + k] = node_info_tmp.z;

        for(i=0;i<Q;i++) {
            if(i>0){
                if(node_info_tmp.nb_info[i-1].node_type == TYPE_SOLID){ //big bug here
                    simgpu_->node_map_h[i*num_node + k] = k + re[i]*num_node;
                }
                else 
                    simgpu_->node_map_h[i*num_node + k] = node_info_tmp.nb_info[i-1].node_id + i*num_node; //big bug here
            }
        }
    }
    fclose(fp_map);

    if(simgpu_->NewOrContinue == 'n')
        initDataArray(simgpu_);
    else
        saveLoadRecovery(it, simgpu_, 'l');

    cudaMemcpy(simgpu_->node_map_d, simgpu_->node_map_h, num_node*Q*sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(simgpu_->f0_d, simgpu_->f0_h, num_node*Q*sizeof(MY_DTYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(simgpu_->g0_d, simgpu_->g0_h, num_node*Q*sizeof(MY_DTYPE), cudaMemcpyHostToDevice);

    cudaMemcpy(simgpu_->rhoF_d, simgpu_->rhoF_h, num_node*sizeof(MY_DTYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(simgpu_->rhoG_d, simgpu_->rhoG_h, num_node*sizeof(MY_DTYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(simgpu_->vxF_d, simgpu_->vxF_h, num_node*sizeof(MY_DTYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(simgpu_->vyF_d, simgpu_->vyF_h, num_node*sizeof(MY_DTYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(simgpu_->vzF_d, simgpu_->vzF_h, num_node*sizeof(MY_DTYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(simgpu_->vxG_d, simgpu_->vxG_h, num_node*sizeof(MY_DTYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(simgpu_->vyG_d, simgpu_->vyG_h, num_node*sizeof(MY_DTYPE), cudaMemcpyHostToDevice);
    cudaMemcpy(simgpu_->vzG_d, simgpu_->vzG_h, num_node*sizeof(MY_DTYPE), cudaMemcpyHostToDevice);
}

void initDataArray(SIMGPU *simgpu_)
{
    FILE *fp_rhoF, *fp_rhoG;
    fp_rhoF = fopen("rhoF.bin", "rb");
    fp_rhoG = fopen("rhoG.bin", "rb");
    if(!(fp_rhoG && fp_rhoF))
    {
        fprintf(stderr, "Rho file openning error!, simgpu.cu\n");
        exit(-1);
    }  
    int nx, ny, nz, num_node;
    unsigned short *n_to_XYZ;
    nx = simgpu_->nx;
    ny = simgpu_->ny;
    nz = simgpu_->nz;
    num_node = simgpu_->num_node;
    n_to_XYZ = simgpu_->n_to_XYZ;

    MY_DTYPE *rhoF, *rhoG;

    rhoF = (MY_DTYPE *)malloc(sizeof(MY_DTYPE)*nx*ny*nz);
    rhoG = (MY_DTYPE *)malloc(sizeof(MY_DTYPE)*nx*ny*nz);
    if(!(rhoF && rhoG)){
        fprintf(stderr, "Memory allocating error, simgpu.c\n");
        exit(-1);
    }

    if((fread(rhoF, sizeof(MY_DTYPE), nx*ny*nz, fp_rhoF)) != nx*ny*nz)
    {
        fprintf(stderr, "Rho file reading error!\n");
        exit(-1);
    }

    if((fread(rhoG, sizeof(MY_DTYPE), nx*ny*nz, fp_rhoG)) != nx*ny*nz)
    {
        fprintf(stderr, "Rho file reading error!\n");
        exit(-1);
    }

    fclose(fp_rhoF);
    fclose(fp_rhoG);

    int x, y, z, k;
    for(k=0; k<num_node; k++)
    {
        x = n_to_XYZ[0*num_node + k];
        y = n_to_XYZ[1*num_node + k];
        z = n_to_XYZ[2*num_node + k];
        simgpu_->rhoF_h[k] = rhoF[z*ny*nx + y*nx + x];
        simgpu_->rhoG_h[k] = rhoG[z*ny*nx + y*nx + x];
        simgpu_->vxF_h[k] = 0.00;
        simgpu_->vyF_h[k] = 0.00;
        simgpu_->vzF_h[k] = 0.00;
        simgpu_->vxG_h[k] = 0.00;
        simgpu_->vyG_h[k] = 0.00;
        simgpu_->vzG_h[k] = 0.00;
    }
    free(rhoF);
    free(rhoG);

    int i;
    for(k=0; k<num_node; k++)
        for(i=0;i<Q;i++) {
            simgpu_->f0_h[i*num_node + k] = w_h[i]*simgpu_->rhoF_h[k];
            simgpu_->g0_h[i*num_node + k] = w_h[i]*simgpu_->rhoG_h[k];
        }
}

void copyBackData(SIMGPU *simgpu_)
{
    cudaMemcpy(simgpu_->f0_h, simgpu_->f0_d, simgpu_->num_node*Q*sizeof(MY_DTYPE), cudaMemcpyDeviceToHost);
    cudaMemcpy(simgpu_->g0_h, simgpu_->g0_d, simgpu_->num_node*Q*sizeof(MY_DTYPE), cudaMemcpyDeviceToHost);

    cudaMemcpy(simgpu_->rhoF_h, simgpu_->rhoF_d, simgpu_->num_node*sizeof(MY_DTYPE), cudaMemcpyDeviceToHost);
    cudaMemcpy(simgpu_->rhoG_h, simgpu_->rhoG_d, simgpu_->num_node*sizeof(MY_DTYPE), cudaMemcpyDeviceToHost);
    cudaMemcpy(simgpu_->vxF_h, simgpu_->vxF_d, simgpu_->num_node*sizeof(MY_DTYPE), cudaMemcpyDeviceToHost);
    cudaMemcpy(simgpu_->vyF_h, simgpu_->vyF_d, simgpu_->num_node*sizeof(MY_DTYPE), cudaMemcpyDeviceToHost);
    cudaMemcpy(simgpu_->vzF_h, simgpu_->vzF_d, simgpu_->num_node*sizeof(MY_DTYPE), cudaMemcpyDeviceToHost);
    cudaMemcpy(simgpu_->vxG_h, simgpu_->vxG_d, simgpu_->num_node*sizeof(MY_DTYPE), cudaMemcpyDeviceToHost);
    cudaMemcpy(simgpu_->vyG_h, simgpu_->vyG_d, simgpu_->num_node*sizeof(MY_DTYPE), cudaMemcpyDeviceToHost);
    cudaMemcpy(simgpu_->vzG_h, simgpu_->vzG_d, simgpu_->num_node*sizeof(MY_DTYPE), cudaMemcpyDeviceToHost);
}

void outputData(SIMGPU *simgpu_, int it)
{
    FILE *fp_rhoF, *fp_rhoG, *fp_vx, *fp_vy, *fp_vz;
    MY_DTYPE *rhoF, *rhoG, *vx, *vy, *vz;
    unsigned short *n_to_XYZ;
    int nz, ny, nx, num_node;
    int size;
    n_to_XYZ = simgpu_->n_to_XYZ;
    nz = simgpu_->nz;
    ny = simgpu_->ny;
    nx = simgpu_->nx;
    size = nx*ny*nz;
    num_node = simgpu_->num_node;
    rhoF = (MY_DTYPE *)calloc(size,sizeof(MY_DTYPE));
    rhoG = (MY_DTYPE *)calloc(size,sizeof(MY_DTYPE));
    vx = (MY_DTYPE *)calloc(size,sizeof(MY_DTYPE));
    vy = (MY_DTYPE *)calloc(size,sizeof(MY_DTYPE));
    vz = (MY_DTYPE *)calloc(size,sizeof(MY_DTYPE));
    if(!(rhoF && rhoG && vx && vy && vz)){
        fprintf(stderr, "Memory allocating error, simgpu.c\n");
        exit(-1);
    }

    int k;
    unsigned short x, y, z;
    /*MY_DTYPE rho_f, rho_g, vx_f, vy_f, vz_f, vx_g, vy_g, vz_g;*/
    for(k=0; k<num_node; k++)
    {
        x = n_to_XYZ[0*num_node + k];
        y = n_to_XYZ[1*num_node + k];
        z = n_to_XYZ[2*num_node + k];
        rhoF[z*nx*ny+y*nx+x] = simgpu_->rhoF_h[k];
        rhoG[z*nx*ny+y*nx+x] = simgpu_->rhoG_h[k];
        vx[z*nx*ny+y*nx+x] = simgpu_->vx_h[k];
        vy[z*nx*ny+y*nx+x] = simgpu_->vy_h[k];
        vz[z*nx*ny+y*nx+x] = simgpu_->vz_h[k];
    }

    char rhoF_filename[100];
    char rhoG_filename[100];
    sprintf(rhoF_filename, "rhoF_%08d.dat", it);
    sprintf(rhoG_filename, "rhoG_%08d.dat", it);

    /*fp_rhoF = fopen("rhoF.dat", "w");*/
    /*fp_rhoG = fopen("rhoG.dat", "w");*/
    fp_rhoF = fopen(rhoF_filename, "w");
    fp_rhoG = fopen(rhoG_filename, "w");

    fp_vx = fopen("vx.dat", "w");
    fp_vy = fopen("vy.dat", "w");
    fp_vz = fopen("vz.dat", "w");
    if(!(fp_rhoF&& fp_rhoG && fp_vx && fp_vy && fp_vz)){
        fprintf(stderr, "rho File openning error, simgpu.cu\n");
        exit(-1);
    }
    /*z = nz/2;*/
    x = nx/2;
    /*for(x=0; x<nx; x++)*/
    for(z=0; z<nz; z++)
    {
        for(y=0; y<ny; y++)
        {
            fprintf(fp_rhoF, "% 16.14e ", rhoF[z*nx*ny + y*nx + x]);
            fprintf(fp_rhoG, "% 16.14e ", rhoG[z*nx*ny + y*nx + x]);
            fprintf(fp_vx, "% 16.14e ", vx[z*nx*ny + y*nx + x]);
            fprintf(fp_vy, "% 16.14e ", vy[z*nx*ny + y*nx + x]);
            fprintf(fp_vz, "% 16.14e ", vz[z*nx*ny + y*nx + x]);
        }
        fprintf(fp_rhoF, "\n");
        fprintf(fp_rhoG, "\n");
        fprintf(fp_vx, "\n");
        fprintf(fp_vy, "\n");
        fprintf(fp_vz, "\n");
    }

    fclose(fp_rhoF);
    fclose(fp_rhoG);
    fclose(fp_vx);
    fclose(fp_vy);
    fclose(fp_vz);
    free(rhoF);
    free(rhoG);
    free(vx);
    free(vy);
    free(vz);
}

void saveTecplot(SIMGPU *simgpu_)
{

    int nx, ny, nz, size, num_node;
    unsigned short *n_to_XYZ;
    nx = simgpu_->nx;
    ny = simgpu_->ny;
    nz = simgpu_->nz;
    size = nx*ny*nz;
    num_node = simgpu_->num_node;
    n_to_XYZ = simgpu_->n_to_XYZ;

    MY_DTYPE rhoff, rhogg,  vxx, vyy, vzz;
    MY_DTYPE *rhoF, *rhoG, *vx, *vy, *vz;
    MY_DTYPE *flag;
    MY_DTYPE *X, *Y, *Z;
    int x, y, z, k ;

    rhoF = (MY_DTYPE *)calloc(size,sizeof(MY_DTYPE));
    rhoG = (MY_DTYPE *)calloc(size,sizeof(MY_DTYPE));
    vx = (MY_DTYPE *)calloc(size,sizeof(MY_DTYPE));
    vy = (MY_DTYPE *)calloc(size,sizeof(MY_DTYPE));
    vz = (MY_DTYPE *)calloc(size,sizeof(MY_DTYPE));
    flag = (MY_DTYPE *)calloc(size,sizeof(MY_DTYPE));
    X = (MY_DTYPE *)calloc(size,sizeof(MY_DTYPE));
    Y = (MY_DTYPE *)calloc(size,sizeof(MY_DTYPE));
    Z = (MY_DTYPE *)calloc(size,sizeof(MY_DTYPE));

    if(!(rhoF && rhoG && vx && vy && vz)){
        fprintf(stderr, "Memory allocating error, simgpu.c\n");
        exit(-1);
    }

    int zi, yi, xi;
    for(zi=0; zi<nz; zi++)
        for(yi=0; yi<ny; yi++)
            for(xi=0; xi<nx; xi++) {
                Z[zi*nx*ny+yi*nx+xi] = (MY_DTYPE)zi;
                Y[zi*nx*ny+yi*nx+xi] = (MY_DTYPE)yi;
                X[zi*nx*ny+yi*nx+xi] = (MY_DTYPE)xi;
                flag[zi*nx*ny+yi*nx+xi] = 1.0;
                rhoF[zi*nx*ny+yi*nx+xi] = -1.0;
                rhoG[zi*nx*ny+yi*nx+xi] = -1.0;
            }

    for(k=0; k<num_node; k++)
    {
        rhoff = simgpu_->rhoF_h[k];
        rhogg = simgpu_->rhoG_h[k];

        vxx = simgpu_->vx_h[k];
        vyy = simgpu_->vy_h[k];
        vzz = simgpu_->vz_h[k];

        x = n_to_XYZ[0*num_node + k];
        y = n_to_XYZ[1*num_node + k];
        z = n_to_XYZ[2*num_node + k];

        rhoF[z*nx*ny+y*nx+x] = rhoff;
        rhoG[z*nx*ny+y*nx+x] = rhogg;
        vx[z*nx*ny+y*nx+x] = vxx;
        vy[z*nx*ny+y*nx+x] = vyy;
        vz[z*nx*ny+y*nx+x] = vzz;
        flag[z*nx*ny+y*nx+x] = 0.0;
    }

    double SolTime;
    INTEGER4 Debug,I, III,DIsDouble,VIsDouble,IMax,JMax,KMax,ZoneType,StrandID,ParentZn,IsBlock;
    INTEGER4 ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn, FileType;

    Debug     = 0;
    VIsDouble = 0;
    DIsDouble = 0;
    IMax      = nx;
    JMax      = ny;
    KMax      = nz;
    ZoneType  = 0;      /* Ordered */
    SolTime   = 0.0;
    StrandID  = 0;     /* StaticZone */
    ParentZn  = 0;      /* No Parent */
    IsBlock   = 1;      /* Block */
    ICellMax  = 0;
    JCellMax  = 0;
    KCellMax  = 0;
    NFConns   = 0;
    FNMode    = 0;
    ShrConn   = 0;
    FileType  = 0;

    I = TECINI112("SIMPLE DATASET",
            "X Y Z rhoF rhoG vx vy vz flag",
            "tec.plt",
            ".",
            &FileType,
            &Debug,
            &VIsDouble);

    I = TECZNE112("Simple Zone",
            &ZoneType,
            &IMax,
            &JMax,
            &KMax,
            &ICellMax,
            &JCellMax,
            &KCellMax,
            &SolTime,
            &StrandID,
            &ParentZn,
            &IsBlock,
            &NFConns,
            &FNMode,
            0,              /* TotalNumFaceNodes */
            0,              /* NumConnectedBoundaryFaces */
            0,              /* TotalNumBoundaryConnections */
            NULL,           /* PassiveVarList */
            NULL,           /* ValueLocation = Nodal */
            NULL,           /* SharVarFromZone */
            &ShrConn);

    III = IMax*JMax*KMax;
    I = TECDAT112(&III,X,&DIsDouble);
    I = TECDAT112(&III,Y,&DIsDouble);
    I = TECDAT112(&III,Z,&DIsDouble);
    I = TECDAT112(&III,rhoF,&DIsDouble);
    I = TECDAT112(&III,rhoG,&DIsDouble);
    I = TECDAT112(&III,vx,&DIsDouble);
    I = TECDAT112(&III,vy,&DIsDouble);
    I = TECDAT112(&III,vz,&DIsDouble);
    I = TECDAT112(&III,flag,&DIsDouble);

    I = TECEND112();
}

void correctVel(SIMGPU *simgpu_)
{
    int num_node, k, i, I, ip;
    MY_DTYPE Fx_F, Fy_F, Fz_F;
    MY_DTYPE Fx_G, Fy_G, Fz_G;
    MY_DTYPE G_F, G_G;
    num_node = simgpu_->num_node;
    MY_DTYPE rhop_F, rhop_G;
    for(k=0; k<num_node; k++)
    {
        Fx_F = Fy_F = Fz_F = 0.0;
        Fx_G = Fy_G = Fz_G = 0.0;
        for(i=0; i<Q; i++)
        {
            I =simgpu_->node_map_h[i*num_node + k]; 
            ip = I/num_node;
            if(ip == i)
            {
                rhop_F = simgpu_->rhoG_h[k];
                rhop_G = simgpu_->rhoF_h[k];
                G_F = G_G = simgpu_->G_fg;
            }else
            {
                rhop_F = rhop_G = simgpu_->rhoS;
                G_F = simgpu_->G_fs;
                G_G = simgpu_->G_gs;
            }
            Fx_F += rhop_F*G_F*w_h[i]*e_h[i][0];
            Fy_F += rhop_F*G_F*w_h[i]*e_h[i][1];
            Fz_F += rhop_F*G_F*w_h[i]*e_h[i][2];
            Fx_G += rhop_G*G_G*w_h[i]*e_h[i][0];
            Fy_G += rhop_G*G_G*w_h[i]*e_h[i][1];
            Fz_G += rhop_G*G_G*w_h[i]*e_h[i][2];
        }

        Fx_F *= -simgpu_->rhoF_h[k]*18.0;
        Fy_F *= -simgpu_->rhoF_h[k]*18.0;
        Fz_F *= -simgpu_->rhoF_h[k]*18.0;
        Fx_G *= -simgpu_->rhoG_h[k]*18.0;
        Fy_G *= -simgpu_->rhoG_h[k]*18.0;
        Fz_G *= -simgpu_->rhoG_h[k]*18.0;

        simgpu_->vx_h[k] = (0.50*simgpu_->dt*Fx_F + simgpu_->vxF_h[k] 
                + 0.50*simgpu_->dt*Fx_G + simgpu_->vxG_h[k])/(
                    simgpu_->rhoF_h[k] + simgpu_->rhoG_h[k]);
        simgpu_->vy_h[k] = (0.50*simgpu_->dt*Fy_F + simgpu_->vyF_h[k] 
                + 0.50*simgpu_->dt*Fy_G + simgpu_->vyG_h[k])/(
                    simgpu_->rhoF_h[k] + simgpu_->rhoG_h[k]);
        simgpu_->vz_h[k] = (0.50*simgpu_->dt*Fz_F + simgpu_->vzF_h[k] 
                + 0.50*simgpu_->dt*Fz_G + simgpu_->vzG_h[k])/(
                    simgpu_->rhoF_h[k] + simgpu_->rhoG_h[k]);
    }
}

MY_DTYPE massError(SIMGPU *simgpu_)
{
    int num_node = simgpu_->num_node;
    int k;
    MY_DTYPE mass=0.0;
    for(k=0; k<num_node; k++)
        mass += simgpu_->rhoF_h[k] + simgpu_->rhoG_h[k];
    return ((mass-num_node)/num_node);
}

void saveLoadRecovery(int *it, SIMGPU *simgpu_, char saveOrLoad)
{
    FILE * fp;
    if(saveOrLoad == 's'){
        if((fp=fopen("recovery.bin", "wb")) == NULL)
        {
            fprintf(stderr, "Recovery data openning error\n");
            exit(-1);
        }

        fwrite(it, sizeof(int), 1, fp);
        if((fwrite(simgpu_->f0_h, sizeof(MY_DTYPE), Q*simgpu_->num_node, fp))!=Q*simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data writing error\n");
            exit(-1);
        }

        if((fwrite(simgpu_->g0_h, sizeof(MY_DTYPE), Q*simgpu_->num_node, fp))!=Q*simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data writing error\n");
            exit(-1);
        }

        if((fwrite(simgpu_->rhoF_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data writing error\n");
            exit(-1);
        }
        if((fwrite(simgpu_->rhoG_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data writing error\n");
            exit(-1);
        }
        if((fwrite(simgpu_->vxF_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data writing error\n");
            exit(-1);
        }
        if((fwrite(simgpu_->vyF_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data writing error\n");
            exit(-1);
        }
        if((fwrite(simgpu_->vzF_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data writing error\n");
            exit(-1);
        }
        if((fwrite(simgpu_->vxG_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data writing error\n");
            exit(-1);
        }
        if((fwrite(simgpu_->vyG_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data writing error\n");
            exit(-1);
        }
        if((fwrite(simgpu_->vzG_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data writing error\n");
            exit(-1);
        }
    }
    else
    {
        if((fp=fopen("recovery.bin", "rb")) == NULL)
        {
            fprintf(stderr, "Recovery data openning error\n");
            exit(-1);
        }
        fread(it, sizeof(int), 1, fp);
        if((fread(simgpu_->f0_h, sizeof(MY_DTYPE), Q*simgpu_->num_node, fp))!=Q*simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data reading error\n");
            exit(-1);
        }
        if((fread(simgpu_->g0_h, sizeof(MY_DTYPE), Q*simgpu_->num_node, fp))!=Q*simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data reading error\n");
            exit(-1);
        }
        if((fread(simgpu_->rhoF_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data reading error\n");
            exit(-1);
        }
        if((fread(simgpu_->rhoG_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data reading error\n");
            exit(-1);
        }
        if((fread(simgpu_->vxF_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data reading error\n");
            exit(-1);
        }
        if((fread(simgpu_->vyF_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data reading error\n");
            exit(-1);
        }
        if((fread(simgpu_->vzF_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data reading error\n");
            exit(-1);
        }
        if((fread(simgpu_->vxG_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data reading error\n");
            exit(-1);
        }
        if((fread(simgpu_->vyG_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data reading error\n");
            exit(-1);
        }
        if((fread(simgpu_->vzG_h, sizeof(MY_DTYPE), simgpu_->num_node, fp))!=simgpu_->num_node)
        {
            fprintf(stderr, "Recovery data reading error\n");
            exit(-1);
        }
    }
    fclose(fp);
}

__global__ void LBCollProp (
        int it,
        int num_node,
        MY_DTYPE omegaF, MY_DTYPE omegaG, MY_DTYPE rhoS,
        MY_DTYPE G_fg, MY_DTYPE G_fs, MY_DTYPE G_gs,
        MY_DTYPE *fInF, MY_DTYPE *fOutF, 
        MY_DTYPE *fInG, MY_DTYPE *fOutG,
        MY_DTYPE *rhoF, MY_DTYPE *rhoG,
        MY_DTYPE *vxF, MY_DTYPE *vxG,
        MY_DTYPE *vyF, MY_DTYPE *vyG,
        MY_DTYPE *vzF, MY_DTYPE *vzG,
        unsigned int *node_map)
{
    int k = blockIdx.y*blockDim.x*gridDim.x + blockIdx.x*blockDim.x + threadIdx.x;
    MY_DTYPE rho_F, vx_F, vy_F, vz_F;
    MY_DTYPE rho_G, vx_G, vy_G, vz_G;
    MY_DTYPE Fx_F, Fy_F, Fz_F;
    MY_DTYPE Fx_G, Fy_G, Fz_G;
    MY_DTYPE F_F, F_G;
    MY_DTYPE vx_eq, vy_eq, vz_eq, vv_eq;
    int  kp;
    int I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I13,I14,I15,I16,I17,I18;
    MY_DTYPE Feq, rhoP_F, rhoP_G;
    MY_DTYPE G_F, G_G;
    int ip;

    if(k < num_node) //valid threads
    {
        rho_F = rhoF[k];
        rho_G = rhoG[k];

        vx_F = vxF[k];
        vy_F = vyF[k];
        vz_F = vzF[k];
        vx_G = vxG[k];
        vy_G = vyG[k];
        vz_G = vzG[k];

        //Optimiaze
        I1  = node_map[ 1*num_node+k]; 
        I2  = node_map[ 2*num_node+k];
        I3  = node_map[ 3*num_node+k];
        I4  = node_map[ 4*num_node+k];
        I5  = node_map[ 5*num_node+k];
        I6  = node_map[ 6*num_node+k];
        I7  = node_map[ 7*num_node+k];
        I8  = node_map[ 8*num_node+k];
        I9  = node_map[ 9*num_node+k];
        I10 = node_map[10*num_node+k];
        I11 = node_map[11*num_node+k];
        I12 = node_map[12*num_node+k];
        I13 = node_map[13*num_node+k];
        I14 = node_map[14*num_node+k];
        I15 = node_map[15*num_node+k];
        I16 = node_map[16*num_node+k];
        I17 = node_map[17*num_node+k];
        I18 = node_map[18*num_node+k];

        Fx_F = Fy_F = Fz_F = 0.00;
        Fx_G = Fy_G = Fz_G = 0.00;

        /////////////////////////////////////////////////////////////////////////////
        ip = I1/num_node;
        kp = I1-ip*num_node;
        if(ip == 1) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F =G_fs; G_G = G_gs; }
        Fx_F += w_d[ 1 ]*rhoP_F*G_F;
        Fx_G += w_d[ 1 ]*rhoP_G*G_G;

        ip = I2/num_node;
        kp = I2-ip*num_node;
        if(ip == 2) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fx_F -= w_d[ 2 ]*rhoP_F*G_F;
        Fx_G -= w_d[ 2 ]*rhoP_G*G_G;

        ip = I3/num_node;
        kp = I3-ip*num_node;
        if(ip == 3) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fy_F += w_d[ 3 ]*rhoP_F*G_F;
        Fy_G += w_d[ 3 ]*rhoP_G*G_G;

        ip = I4/num_node;
        kp = I4-ip*num_node;
        if(ip == 4) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fy_F -= w_d[ 4 ]*rhoP_F*G_F;
        Fy_G -= w_d[ 4 ]*rhoP_G*G_G;

        ip = I5/num_node;
        kp = I5-ip*num_node;
        if(ip == 5) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fz_F += w_d[ 5 ]*rhoP_F*G_F;
        Fz_G += w_d[ 5 ]*rhoP_G*G_G;

        ip = I6/num_node;
        kp = I6-ip*num_node;
        if(ip == 6) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fz_F -= w_d[ 6 ]*rhoP_F*G_F;
        Fz_G -= w_d[ 6 ]*rhoP_G*G_G;

        ip = I7/num_node;
        kp = I7-ip*num_node;
        if(ip == 7) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fy_F += w_d[ 7 ]*rhoP_F*G_F;
        Fy_G += w_d[ 7 ]*rhoP_G*G_G;
        Fz_F += w_d[ 7 ]*rhoP_F*G_F;
        Fz_G += w_d[ 7 ]*rhoP_G*G_G;

        ip = I8/num_node;
        kp = I8-ip*num_node;
        if(ip == 8) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fy_F -= w_d[ 8 ]*rhoP_F*G_F;
        Fy_G -= w_d[ 8 ]*rhoP_G*G_G;
        Fz_F -= w_d[ 8 ]*rhoP_F*G_F;
        Fz_G -= w_d[ 8 ]*rhoP_G*G_G;

        ip = I9/num_node;
        kp = I9-ip*num_node;
        if(ip == 9) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fy_F -= w_d[ 9 ]*rhoP_F*G_F;
        Fy_G -= w_d[ 9 ]*rhoP_G*G_G;
        Fz_F += w_d[ 9 ]*rhoP_F*G_F;
        Fz_G += w_d[ 9 ]*rhoP_G*G_G;

        ip = I10/num_node;
        kp = I10-ip*num_node;
        if(ip == 10) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fy_F += w_d[ 10 ]*rhoP_F*G_F;
        Fy_G += w_d[ 10 ]*rhoP_G*G_G;
        Fz_F -= w_d[ 10 ]*rhoP_F*G_F;
        Fz_G -= w_d[ 10 ]*rhoP_G*G_G;

        ip = I11/num_node;
        kp = I11-ip*num_node;
        if(ip == 11) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fx_F -= w_d[ 11 ]*rhoP_F*G_F;
        Fx_G -= w_d[ 11 ]*rhoP_G*G_G;
        Fz_F -= w_d[ 11 ]*rhoP_F*G_F;
        Fz_G -= w_d[ 11 ]*rhoP_G*G_G;

        ip = I12/num_node;
        kp = I12-ip*num_node;
        if(ip == 12) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fx_F += w_d[ 12 ]*rhoP_F*G_F;
        Fx_G += w_d[ 12 ]*rhoP_G*G_G;
        Fz_F += w_d[ 12 ]*rhoP_F*G_F;
        Fz_G += w_d[ 12 ]*rhoP_G*G_G;

        ip = I13/num_node;
        kp = I13-ip*num_node;
        if(ip == 13) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fx_F -= w_d[ 13 ]*rhoP_F*G_F;
        Fx_G -= w_d[ 13 ]*rhoP_G*G_G;
        Fz_F += w_d[ 13 ]*rhoP_F*G_F;
        Fz_G += w_d[ 13 ]*rhoP_G*G_G;

        ip = I14/num_node;
        kp = I14-ip*num_node;
        if(ip == 14) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fx_F += w_d[ 14 ]*rhoP_F*G_F;
        Fx_G += w_d[ 14 ]*rhoP_G*G_G;
        Fz_F -= w_d[ 14 ]*rhoP_F*G_F;
        Fz_G -= w_d[ 14 ]*rhoP_G*G_G;

        ip = I15/num_node;
        kp = I15-ip*num_node;
        if(ip == 15) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fx_F -= w_d[ 15 ]*rhoP_F*G_F;
        Fx_G -= w_d[ 15 ]*rhoP_G*G_G;
        Fy_F += w_d[ 15 ]*rhoP_F*G_F;
        Fy_G += w_d[ 15 ]*rhoP_G*G_G;

        ip = I16/num_node;
        kp = I16-ip*num_node;
        if(ip == 16) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fx_F += w_d[ 16 ]*rhoP_F*G_F;
        Fx_G += w_d[ 16 ]*rhoP_G*G_G;
        Fy_F -= w_d[ 16 ]*rhoP_F*G_F;
        Fy_G -= w_d[ 16 ]*rhoP_G*G_G;

        ip = I17/num_node;
        kp = I17-ip*num_node;
        if(ip == 17) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fx_F -= w_d[ 17 ]*rhoP_F*G_F;
        Fx_G -= w_d[ 17 ]*rhoP_G*G_G;
        Fy_F -= w_d[ 17 ]*rhoP_F*G_F;
        Fy_G -= w_d[ 17 ]*rhoP_G*G_G;

        ip = I18/num_node;
        kp = I18-ip*num_node;
        if(ip == 18) { rhoP_F = rhoG[kp]; rhoP_G = rhoF[kp]; G_F = G_G = G_fg;
        } else{ rhoP_F = rhoP_G = rhoS; G_F = G_fs; G_G = G_gs; }
        Fx_F += w_d[ 18 ]*rhoP_F*G_F;
        Fx_G += w_d[ 18 ]*rhoP_G*G_G;
        Fy_F += w_d[ 18 ]*rhoP_F*G_F;
        Fy_G += w_d[ 18 ]*rhoP_G*G_G;
        /////////////////////////////////////////////////////////////////////////////

        Fx_F *=-18.00;
        Fy_F *=-18.00;
        Fz_F *=-18.00;
        Fx_G *=-18.00;
        Fy_G *=-18.00;
        Fz_G *=-18.00;

        vx_F += 0.50*Fx_F*rho_F;
        vy_F += 0.50*Fy_F*rho_F;
        vz_F += 0.50*Fz_F*rho_F;
        vx_G += 0.50*Fx_G*rho_G;
        vy_G += 0.50*Fy_G*rho_G;
        vz_G += 0.50*Fz_G*rho_G;

        vx_eq = (vx_F*omegaF + vx_G*omegaG)/(rho_F*omegaF + rho_G*omegaG);
        vy_eq = (vy_F*omegaF + vy_G*omegaG)/(rho_F*omegaF + rho_G*omegaG);
        vz_eq = (vz_F*omegaF + vz_G*omegaG)/(rho_F*omegaF + rho_G*omegaG);
        vv_eq = 1.50*(vx_eq*vx_eq + vy_eq*vy_eq + vz_eq*vz_eq);

        if(it==0)
        {
            vx_eq = 0.00;
            vy_eq = 0.00;
            vz_eq = 0.00;
            vv_eq = 0.00;
        }else{
            vxF[k] = vx_F;
            vyF[k] = vy_F;
            vzF[k] = vz_F;
            vxG[k] = vx_G;
            vyG[k] = vy_G;
            vzG[k] = vz_G;
        }

        //Force
        Fz_F += 1.0e-5;
        Fz_G += 1.0e-5;

        Fx_F *= 3*(1.00 - 0.50*omegaF);
        Fy_F *= 3*(1.00 - 0.50*omegaF);
        Fz_F *= 3*(1.00 - 0.50*omegaF);
        Fx_G *= 3*(1.00 - 0.50*omegaG);
        Fy_G *= 3*(1.00 - 0.50*omegaG);
        Fz_G *= 3*(1.00 - 0.50*omegaG);

        //Optimiaze
        /*Feq = 0.0;*/
        //end Optimiaze

        ////////////////////////////////////////////////////////////////////////////
        F_F = Fx_F*(e_d[0][0] - vx_eq) + Fy_F*(e_d[0][1] - vy_eq) +Fz_F*(e_d[0][2] - vz_eq);
        F_G = Fx_G*(e_d[0][0] - vx_eq) + Fy_G*(e_d[0][1] - vy_eq) +Fz_G*(e_d[0][2] - vz_eq);
        Feq = FEQ_0(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[k ] = (1.0 - omegaF)*fInF[k+0*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[k ] = (1.0 - omegaG)*fInG[k+0*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[1][0] - vx_eq) + Fy_F*(e_d[1][1] - vy_eq) +Fz_F*(e_d[1][2] - vz_eq);
        F_G = Fx_G*(e_d[1][0] - vx_eq) + Fy_G*(e_d[1][1] - vy_eq) +Fz_G*(e_d[1][2] - vz_eq);
        Feq = FEQ_1(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I1] = (1.0 - omegaF)*fInF[k+1*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I1] = (1.0 - omegaG)*fInG[k+1*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[2][0] - vx_eq) + Fy_F*(e_d[2][1] - vy_eq) +Fz_F*(e_d[2][2] - vz_eq);
        F_G = Fx_G*(e_d[2][0] - vx_eq) + Fy_G*(e_d[2][1] - vy_eq) +Fz_G*(e_d[2][2] - vz_eq);
        Feq = FEQ_2(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I2] = (1.0 - omegaF)*fInF[k+2*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I2] = (1.0 - omegaG)*fInG[k+2*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[3][0] - vx_eq) + Fy_F*(e_d[3][1] - vy_eq) +Fz_F*(e_d[3][2] - vz_eq);
        F_G = Fx_G*(e_d[3][0] - vx_eq) + Fy_G*(e_d[3][1] - vy_eq) +Fz_G*(e_d[3][2] - vz_eq);
        Feq = FEQ_3(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I3] = (1.0 - omegaF)*fInF[k+3*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I3] = (1.0 - omegaG)*fInG[k+3*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[4][0] - vx_eq) + Fy_F*(e_d[4][1] - vy_eq) +Fz_F*(e_d[4][2] - vz_eq);
        F_G = Fx_G*(e_d[4][0] - vx_eq) + Fy_G*(e_d[4][1] - vy_eq) +Fz_G*(e_d[4][2] - vz_eq);
        Feq = FEQ_4(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I4] = (1.0 - omegaF)*fInF[k+4*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I4] = (1.0 - omegaG)*fInG[k+4*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[5][0] - vx_eq) + Fy_F*(e_d[5][1] - vy_eq) +Fz_F*(e_d[5][2] - vz_eq);
        F_G = Fx_G*(e_d[5][0] - vx_eq) + Fy_G*(e_d[5][1] - vy_eq) +Fz_G*(e_d[5][2] - vz_eq);
        Feq = FEQ_5(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I5] = (1.0 - omegaF)*fInF[k+5*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I5] = (1.0 - omegaG)*fInG[k+5*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[6][0] - vx_eq) + Fy_F*(e_d[6][1] - vy_eq) +Fz_F*(e_d[6][2] - vz_eq);
        F_G = Fx_G*(e_d[6][0] - vx_eq) + Fy_G*(e_d[6][1] - vy_eq) +Fz_G*(e_d[6][2] - vz_eq);
        Feq = FEQ_6(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I6] = (1.0 - omegaF)*fInF[k+6*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I6] = (1.0 - omegaG)*fInG[k+6*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[7][0] - vx_eq) + Fy_F*(e_d[7][1] - vy_eq) +Fz_F*(e_d[7][2] - vz_eq);
        F_G = Fx_G*(e_d[7][0] - vx_eq) + Fy_G*(e_d[7][1] - vy_eq) +Fz_G*(e_d[7][2] - vz_eq);
        Feq = FEQ_7(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I7] = (1.0 - omegaF)*fInF[k+7*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I7] = (1.0 - omegaG)*fInG[k+7*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[8][0] - vx_eq) + Fy_F*(e_d[8][1] - vy_eq) +Fz_F*(e_d[8][2] - vz_eq);
        F_G = Fx_G*(e_d[8][0] - vx_eq) + Fy_G*(e_d[8][1] - vy_eq) +Fz_G*(e_d[8][2] - vz_eq);
        Feq = FEQ_8(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I8] = (1.0 - omegaF)*fInF[k+8*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I8] = (1.0 - omegaG)*fInG[k+8*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[9][0] - vx_eq) + Fy_F*(e_d[9][1] - vy_eq) +Fz_F*(e_d[9][2] - vz_eq);
        F_G = Fx_G*(e_d[9][0] - vx_eq) + Fy_G*(e_d[9][1] - vy_eq) +Fz_G*(e_d[9][2] - vz_eq);
        Feq = FEQ_9(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I9] = (1.0 - omegaF)*fInF[k+9*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I9] = (1.0 - omegaG)*fInG[k+9*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[10][0] - vx_eq) + Fy_F*(e_d[10][1] - vy_eq) +Fz_F*(e_d[10][2] - vz_eq);
        F_G = Fx_G*(e_d[10][0] - vx_eq) + Fy_G*(e_d[10][1] - vy_eq) +Fz_G*(e_d[10][2] - vz_eq);
        Feq = FEQ_10(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I10] = (1.0 - omegaF)*fInF[k+10*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I10] = (1.0 - omegaG)*fInG[k+10*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[11][0] - vx_eq) + Fy_F*(e_d[11][1] - vy_eq) +Fz_F*(e_d[11][2] - vz_eq);
        F_G = Fx_G*(e_d[11][0] - vx_eq) + Fy_G*(e_d[11][1] - vy_eq) +Fz_G*(e_d[11][2] - vz_eq);
        Feq = FEQ_11(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I11] = (1.0 - omegaF)*fInF[k+11*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I11] = (1.0 - omegaG)*fInG[k+11*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[12][0] - vx_eq) + Fy_F*(e_d[12][1] - vy_eq) +Fz_F*(e_d[12][2] - vz_eq);
        F_G = Fx_G*(e_d[12][0] - vx_eq) + Fy_G*(e_d[12][1] - vy_eq) +Fz_G*(e_d[12][2] - vz_eq);
        Feq = FEQ_12(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I12] = (1.0 - omegaF)*fInF[k+12*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I12] = (1.0 - omegaG)*fInG[k+12*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[13][0] - vx_eq) + Fy_F*(e_d[13][1] - vy_eq) +Fz_F*(e_d[13][2] - vz_eq);
        F_G = Fx_G*(e_d[13][0] - vx_eq) + Fy_G*(e_d[13][1] - vy_eq) +Fz_G*(e_d[13][2] - vz_eq);
        Feq = FEQ_13(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I13] = (1.0 - omegaF)*fInF[k+13*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I13] = (1.0 - omegaG)*fInG[k+13*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[14][0] - vx_eq) + Fy_F*(e_d[14][1] - vy_eq) +Fz_F*(e_d[14][2] - vz_eq);
        F_G = Fx_G*(e_d[14][0] - vx_eq) + Fy_G*(e_d[14][1] - vy_eq) +Fz_G*(e_d[14][2] - vz_eq);
        Feq = FEQ_14(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I14] = (1.0 - omegaF)*fInF[k+14*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I14] = (1.0 - omegaG)*fInG[k+14*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[15][0] - vx_eq) + Fy_F*(e_d[15][1] - vy_eq) +Fz_F*(e_d[15][2] - vz_eq);
        F_G = Fx_G*(e_d[15][0] - vx_eq) + Fy_G*(e_d[15][1] - vy_eq) +Fz_G*(e_d[15][2] - vz_eq);
        Feq = FEQ_15(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I15] = (1.0 - omegaF)*fInF[k+15*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I15] = (1.0 - omegaG)*fInG[k+15*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[16][0] - vx_eq) + Fy_F*(e_d[16][1] - vy_eq) +Fz_F*(e_d[16][2] - vz_eq);
        F_G = Fx_G*(e_d[16][0] - vx_eq) + Fy_G*(e_d[16][1] - vy_eq) +Fz_G*(e_d[16][2] - vz_eq);
        Feq = FEQ_16(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I16] = (1.0 - omegaF)*fInF[k+16*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I16] = (1.0 - omegaG)*fInG[k+16*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[17][0] - vx_eq) + Fy_F*(e_d[17][1] - vy_eq) +Fz_F*(e_d[17][2] - vz_eq);
        F_G = Fx_G*(e_d[17][0] - vx_eq) + Fy_G*(e_d[17][1] - vy_eq) +Fz_G*(e_d[17][2] - vz_eq);
        Feq = FEQ_17(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I17] = (1.0 - omegaF)*fInF[k+17*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I17] = (1.0 - omegaG)*fInG[k+17*num_node] + (omegaG + F_G)*Feq*rho_G;

        F_F = Fx_F*(e_d[18][0] - vx_eq) + Fy_F*(e_d[18][1] - vy_eq) +Fz_F*(e_d[18][2] - vz_eq);
        F_G = Fx_G*(e_d[18][0] - vx_eq) + Fy_G*(e_d[18][1] - vy_eq) +Fz_G*(e_d[18][2] - vz_eq);
        Feq = FEQ_18(vx_eq, vy_eq, vz_eq, vv_eq);
        fOutF[I18] = (1.0 - omegaF)*fInF[k+18*num_node] + (omegaF + F_F)*Feq*rho_F;
        fOutG[I18] = (1.0 - omegaG)*fInG[k+18*num_node] + (omegaG + F_G)*Feq*rho_G;
    }
}


__global__ void LBUpdateMacro(
        int num_node,
        MY_DTYPE *rhoF, MY_DTYPE *rhoG,
        MY_DTYPE *vxF, MY_DTYPE *vxG,
        MY_DTYPE *vyF, MY_DTYPE *vyG,
        MY_DTYPE *vzF, MY_DTYPE *vzG,
        MY_DTYPE *fInF, MY_DTYPE *fInG)
{
    int k = blockIdx.y*blockDim.x*gridDim.x + blockIdx.x*blockDim.x + threadIdx.x;
    MY_DTYPE fi, rho, vx, vy, vz;
    int i;
    if(k < num_node) //valid threads
    {

        rho = vx = vy = vz = 0.00;
#pragma unroll
        for(i=0; i<Q; i++){
            fi = fInF[i*num_node + k];
            rho += fi;
            vx += fi*e_d[i][0];
            vy += fi*e_d[i][1];
            vz += fi*e_d[i][2];
        }
        rhoF[k] = rho;
        vxF[k] = vx;
        vyF[k] = vy;
        vzF[k] = vz;

        rho = vx = vy = vz = 0.00;
#pragma unroll
        for(i=0; i<Q; i++){
            fi = fInG[i*num_node + k];
            rho += fi;
            vx += fi*e_d[i][0];
            vy += fi*e_d[i][1];
            vz += fi*e_d[i][2];
        }
        rhoG[k] = rho;
        vxG[k] = vx;
        vyG[k] = vy;
        vzG[k] = vz;
    }
}
