#include <stdio.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <time.h>
#include <string.h>
#include "cuda_runtime.h"
#include "mytype.h"
#include "simgpu.cu"

int main(int argc, char* argv[])
{
    float elpsd_time;
    clock_t time_begin, time_end; 

    SIMGPU sim;
    int it=0;
    if(argc != 4){
        fprintf(stderr, "Useage: ./simgpu itMax n|c GPU_ID\n");
        exit(-2);
    }

    setPara(&sim, argc, argv);
    allocMemroy(&sim);
    buildDataArray(&it, &sim);
    /*printf("It = %d \n", it);*/
    sim.itMax += it;

    int GX;
    int GXmax = 256;
    GX = (sim.num_node + BX-1)/BX;
    GX = (GX>GXmax) ? GXmax : GX;
    dim3 grid(GX, (sim.num_node + GX*BX-1)/(GX*BX), 1);  //grid layout for LBCollProp kernel

    cudaSetDevice(atoi(argv[3]));
    printf("Using GPU %d\n", atoi(argv[3]));
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start,0);

    time_begin = clock();
    //----------------------8<-----------------------------------------------
    /*K_old = 1.0f;*/
    MY_DTYPE *fOld, *fNew;
    MY_DTYPE *gOld, *gNew;
    for(; it<sim.itMax; it++) {
        if(it%2==0) {
            fOld = sim.f0_d;
            fNew = sim.f1_d;
            gOld = sim.g0_d;
            gNew = sim.g1_d;
        } else {
            fOld = sim.f1_d;
            fNew = sim.f0_d;
            gOld = sim.g1_d;
            gNew = sim.g0_d;
        }

        LBCollProp<<<grid,  BX>>>(
                it,
                sim.num_node,
                sim.omegaF, sim.omegaG, sim.rhoS,
                sim.G_fg, sim.G_fs, sim.G_gs,
                fOld, fNew, gOld, gNew,
                sim.rhoF_d, sim.rhoG_d,
                sim.vxF_d, sim.vxG_d, sim.vyF_d, sim.vyG_d, sim.vzF_d, sim.vzG_d,
                sim.node_map_d);
        /*getLastCudaError("LBCollProp Launch failed!\n");*/

        LBUpdateMacro<<<grid, BX>>>(
                sim.num_node,
                sim.rhoF_d, sim.rhoG_d,
                sim.vxF_d, sim.vxG_d, sim.vyF_d, sim.vyG_d, sim.vzF_d, sim.vzG_d,
                fNew, gNew);
        /*getLastCudaError("LBUpdateMacro Launch failed!\n");*/

        if(it%1000==0) {
            copyBackData(&sim);
            /*outputData(&sim, it);*/
            printf("it : %8d, f11=% e, massError=% e\n", it, sim.f0_h[sim.num_node*12 + 1000], massError(&sim));
        }
    }
    time_end = clock();
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elpsd_time, start, stop);
    printf("Time          : %f seconds \n", elpsd_time/1000.0);
    printf("Speed(CPU API): %f MLUPS\n", (float)(sim.num_node)*(float)(sim.itMax)/1000000.0/((float)(time_end-time_begin)/CLOCKS_PER_SEC));
    printf("Speed(CUDA API): %fMLUPS\n", (float)(sim.num_node)*(float)(sim.itMax)/1000000.0/elpsd_time*1000.0);

    copyBackData(&sim);
    correctVel(&sim);
    outputData(&sim, it);
    saveLoadRecovery(&it, &sim, 's');
    saveTecplot(&sim);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    freeMemrory(&sim);
    return 0;
}
