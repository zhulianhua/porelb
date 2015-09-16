#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mytype.h"
#include "preprocess.h"
#include "node.h"

#define N 256
#define NZ 64
#define NY 64
#define NX 64

unsigned char flag[NZ][NY][NX];
MY_DTYPE rhoF[NZ][NY][NX];
MY_DTYPE rhoG[NZ][NY][NX];

void buildSaveFlag(MY_DTYPE R);

void buildSaveFlag(MY_DTYPE R)
{
    MY_DTYPE L = (MY_DTYPE)(N-1);
    int x, y, z;
    MY_DTYPE xp, yp, zp;
    MY_DTYPE xc, yc, zc;
    MY_DTYPE ll;
    FILE *fp, *fp_rhoF, *fp_rhoG;
    MY_DTYPE rho_f, rho_g;

    rho_f = 1.0;
    rho_g = 1.0;

    fp = fopen("flag.bin", "rb");
    if(fp==NULL)
    {
        fprintf(stderr, "Flag file openning error!\n");
        exit(-2);
    }
    if((fread(flag, sizeof(unsigned char), NZ*NY*NX, fp)) != NZ*NY*NX)
    {
        fprintf(stderr, "Reading flag file error!\n");
        exit(-1);
    }
    fclose(fp);

    int count=0;
    MY_DTYPE rand_tmp;
    srand(time(NULL));
    for(z=0; z<NZ; z++)
        for(y=0; y<NY; y++)
            for(x=0; x<NX; x++)
            {
                /*flag[z][y][x] == 0;*/
                rhoF[z][y][x] = 1.0;
                rhoG[z][y][x] = 1.0;
                if(flag[z][y][x] == 1){
                    rhoG[z][y][x] = -1;
                    rhoF[z][y][x] = -1;
                }
                if(flag[z][y][x] == 0){
                    rand_tmp = rand()/(MY_DTYPE)RAND_MAX;
                    if(rand_tmp < 0.5)
                    {
                        rhoF[z][y][x] = 0;
                        rhoG[z][y][x] = 1;
                    }else
                    {
                        rhoF[z][y][x] = 1;
                        rhoG[z][y][x] = 0;
                    }
                    count++;
                }
            }

    printf("count: %d\n", count);
    /*printf()*/
    /*fp = fopen("flag.bin", "wb");*/

    fp_rhoF = fopen("rhoF.bin", "wb");
    fp_rhoG = fopen("rhoG.bin", "wb");

    fp = fopen("flag.bin", "wb");
    if(fp==NULL)
    {
        fprintf(stderr, "Flag file openning error!\n");
        exit(-2);
    }
    if((fwrite(flag, sizeof(unsigned char), NZ*NY*NX, fp)) != NZ*NY*NX)
    {
        fprintf(stderr, "Writing flag file error!\n");
        exit(-1);
    }
    fclose(fp);

    if(!(fp_rhoF && fp_rhoG))
    {
        fprintf(stderr, "Openning file error!, testpre.c\n");
        exit(-1);
    }


    if((fwrite(rhoF, sizeof(MY_DTYPE), NZ*NY*NX, fp_rhoF)) != NZ*NY*NX)
    {
        fprintf(stderr, "Writing flag file error!\n");
        exit(-1);
    }
    fclose(fp_rhoF);

    if((fwrite(rhoG, sizeof(MY_DTYPE), NZ*NY*NX, fp_rhoG)) != NZ*NY*NX)
    {
        fprintf(stderr, "Writing flag file error!\n");
        exit(-1);
    }
    fclose(fp_rhoG);
}

void saveNode()
{
    FILE *fp;
    int num_p;
    int num_node;
    fp = fopen("V.bin", "rb");
    NODE_INFO node_info_tmp;
    fread(&num_node, sizeof(int), 1,fp);
    fread(&num_p, sizeof(int), 1,fp);
    fseek(fp, 3*sizeof(int),SEEK_CUR);
    fseek(fp, num_p*sizeof(int),SEEK_CUR);
    int i;
    for(i=0; i<num_node; i++)
    {
        if((fread(&node_info_tmp,sizeof(NODE_INFO), 1, fp)) != 1)
        {
            fprintf(stderr, "Node map read error!\n");
            exit(-1);
        }
    }
    fclose(fp);
}

int main(int argc, char* argv[])
{

    MY_DTYPE R;
    if(argc != 2)
    {
        fprintf(stderr, "Argc error: ./a.out R\n");
        exit(-2);
    }

    R = atof(argv[1]);
    V_DESC v_desc;
    buildSaveFlag(R);
    v_desc = build_V(1, NX,NY,NZ, "flag.bin");
    save_V(&v_desc, "V.bin");
    free_v_desc(&v_desc);
    saveNode();
    return 0;
}
