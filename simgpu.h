#ifndef __SIMGPU_H_
#define __SIMGPU_H_

typedef struct {
    //parameters
    int nx, ny, nz;
    int num_node;
    MY_DTYPE L,  dt, dx;
    MY_DTYPE nuF, omegaF, tauF;
    MY_DTYPE nuG, omegaG, tauG;
    MY_DTYPE G_fg, G_fs, G_gs;
    MY_DTYPE rho_f, rho_g, rhoS;

    MY_DTYPE Fx;
    int itMax, itPrint, itSave;
    char node_map_filename[100];

    char NewOrContinue;

    //data arrays
    MY_DTYPE *f0_h, *f0_d, *f1_d; //0 and 1 are two swap lattice 
    MY_DTYPE *g0_h, *g0_d, *g1_d;
    MY_DTYPE *rhoF_h, *rhoG_h;
    MY_DTYPE *vxF_h, *vyF_h, *vzF_h;
    MY_DTYPE *vxG_h, *vyG_h, *vzG_h;
    MY_DTYPE *vx_h, *vy_h, *vz_h;

    //MY_DTYPE *rho0_h, *rho1_h, *vx_h, *vy_h, *vz_h;
    MY_DTYPE *rhoF_d, *rhoG_d;

    MY_DTYPE *vxF_d, *vyF_d, *vzF_d;
    MY_DTYPE *vxG_d, *vyG_d, *vzG_d;

    unsigned int *node_map_h, *node_map_d;
    unsigned short *n_to_XYZ;
} SIMGPU;


void setPara(SIMGPU *simgpu_,  int argc,  char* argv[]);

void allocMemroy(SIMGPU *simgpu_);

void freeMemrory(SIMGPU *simgpu_);

void initDataArray(SIMGPU *simgpu_);

void buildDataArray(int *it, SIMGPU *simgpu_);


void copyBackData(SIMGPU *simgpu_);

void outPutData(SIMGPU *simgpu_, int it);

MY_DTYPE massError(SIMGPU *simgpu_);

void saveTecplot(SIMGPU *simgpu_);

void correctVel(SIMGPU *simgpu_);

void saveLoadRecovery(int *it, SIMGPU *simgpu_, char saveOrLoad);

__global__ void LBCollProp (
        MY_DTYPE omega, MY_DTYPE FxDt_omega, 
        MY_DTYPE *fIn, MY_DTYPE *fOut, int *node_map);

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
        unsigned int *node_map);

__global__ void LBUpdateMacro(
        int num_node,
        MY_DTYPE *rhoF, MY_DTYPE *rhoG,
        MY_DTYPE *vxF, MY_DTYPE *vxG,
        MY_DTYPE *vyF, MY_DTYPE *vyG,
        MY_DTYPE *vzF, MY_DTYPE *vzG,
        MY_DTYPE *fInF, MY_DTYPE *fInG);
#endif
