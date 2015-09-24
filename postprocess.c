#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "porelb.h"
#include "TECIO.h"

double *u, *v, *w, *rho, *X, *Y, *Z;
unsigned char * flag;
int * time_list;
int nx, ny, nz, nt;

void print_help()
{
    printf("Usage: ./postprocess.x  [A|X|Y|Z] [n] time_list\n");
    printf("Output 3D full data     : ./postprocess.x A time_list\n");
    printf("Output slice at  X = 20 : ./postprocess.x X 20  time_list\n");
}

void output_slice_X(int n)
{
    FILE* fp;
    char filename[100];
    for(int it = 0; it<nt; it++)
    {
        sprintf(filename, "Slice_X_%07d.dat", time_list[it]);
        fp = fopen(filename, "w");
        fprintf(fp, "TITLE = \" X slice at %04d\"\n", n);
        fprintf(fp, "Variables = \"X\", \"Y\", \"Rho\", \"Ux\", \"Uy\", \"Uz\", \"flag\"\n");
        fprintf(fp, "ZONE I = %4d, J = %4d, F=POINT\n", ny, nz);
        for(int k=0; k<nz; k++)
            for(int j=0; j<ny; j++)
            {
                long int loc = k*nx*ny + j*nx + n;
                fprintf(fp,  "%f %f % 21.14e  % 21.14E  % 21.14E  % 21.14E  % 8.1e\n" ,
                        Z[loc], Y[loc], rho[loc], u[loc],  v[loc],   w[loc],   (int)flag[loc]);
            }

        fclose(fp);
    }
}

int main(int argc, char* argv[])
{
    int num_p, num_f;
    int n_slice;

    if(argc < 2) {print_help(); exit(-2);}
    /*if(    strcmp(argv[1], "A")*/
            /*|| strcmp(argv[1], "X")*/
            /*|| strcmp(argv[1], "Y")*/
            /*|| strcmp(argv[1], "Z") ){print_help(); exit(-2);}*/
    if (!strcmp(argv[1], "A"))
    {
        nt = argc - 2;
        time_list = (int*)calloc(nt, sizeof(int));
        for(int i=0; i<nt; i++)
            time_list[i] = atoi(argv[i+2]);
    }
    else 
    {
        printf("argc = %d\n", argc);
        nt = argc - 3;
        time_list = (int*)calloc(nt, sizeof(int));
        n_slice = atoi(argv[2]);
        for(int i=0; i<nt; i++)
            time_list[i] = atoi(argv[i+3]);
    }

    printf("Total %d time seriese data\n", nt);


    FILE *fp_map;
    if((fp_map = fopen("V.bin", "rb")) == NULL)
    {
        fprintf(stderr, "Node map file openning error!\n");
        exit(-1);
    }
    fread(&num_f, sizeof(int), 1, fp_map);
    printf(">>> Total num of fluid  nodes = %d\n", num_f);

    fread(&num_p,  sizeof(int), 1, fp_map);
    printf(">>> num_p = %d\n", num_p);
    int *start_loc = (int*) calloc(num_p, sizeof(int));

    fread(&nx,  sizeof(int), 1, fp_map);
    fread(&ny,  sizeof(int), 1, fp_map);
    fread(&nz,  sizeof(int), 1, fp_map);

    long size = nx*ny*nz;

    rho = (double*) calloc(size, sizeof(double));
    u   = (double*) calloc(size, sizeof(double));
    v   = (double*) calloc(size, sizeof(double));
    w   = (double*) calloc(size, sizeof(double));
    X   = (double*) calloc(size, sizeof(double));
    Y   = (double*) calloc(size, sizeof(double));
    Z   = (double*) calloc(size, sizeof(double));

    flag= (unsigned char*) calloc(size, sizeof(unsigned char));

    FILE* fp_flag;
    if((fp_flag = fopen("flag.bin", "rb")) == NULL)
    {
        fprintf(stderr, "Flag file openning error!\n");
        exit(-1);
    }
    fread(flag, sizeof(unsigned char), size, fp_flag);
    fclose(fp_flag);

    fread(start_loc, sizeof(int), num_p, fp_map);
    NODE_INFO *V = (NODE_INFO *)calloc(num_f, sizeof(NODE_INFO));
    fread(V, sizeof(NODE_INFO), num_f, fp_map);
    fclose(fp_map);

    for(int zi=0; zi<nz; zi++)
        for(int yi=0; yi<ny; yi++)
            for(int xi=0; xi<nx; xi++) {
                Z[zi*nx*ny+yi*nx+xi] = (double)zi;
                Y[zi*nx*ny+yi*nx+xi] = (double)yi;
                X[zi*nx*ny+yi*nx+xi] = (double)xi;
            }

    char partname[100];
    char fullname[100];
    FILE* fp_part;

    double *rho_p, *u_p, *v_p, *w_p;
    int N;

    for(int it = 0; it<nt; it++)
    {
        sprintf(fullname, "full_%07d.dat", time_list[it]);
        for(int i=0; i<num_p; i++)
        {
            sprintf(partname, "part_%03d_%07d.bin", i, time_list[it]);
            if((fp_part = fopen(partname, "rb")) == NULL)
            {
                fprintf(stderr, "Part file openning error!\n");
                exit(-1);
            }
            if(i < num_p -1) 
                N = start_loc[i+1] - start_loc[i];
            else 
                N = num_f - start_loc[i];

            rho_p = (double*) calloc(N, sizeof(double));
            u_p   = (double*) calloc(N, sizeof(double));
            v_p   = (double*) calloc(N, sizeof(double));
            w_p   = (double*) calloc(N, sizeof(double));

            fread(rho_p, sizeof(double), N, fp_part);
            fread(u_p,   sizeof(double), N, fp_part);
            fread(v_p,   sizeof(double), N, fp_part);
            fread(w_p,   sizeof(double), N, fp_part);
            fclose(fp_part);

            int x, y, z, loc;
            for(int n=0; n<N; n++)
            {
                x = V[start_loc[i] + n].x;
                y = V[start_loc[i] + n].y;
                z = V[start_loc[i] + n].z;
                loc = z*nx*ny + y*nx + x;
                rho[loc] = rho_p[n];
                u[loc] = u_p[n];
                v[loc] = v_p[n];
                w[loc] = w_p[n];
            }
        } //processors


        if (!strcmp(argv[1], "A"))
        {
            double SolTime;
            INTEGER4 Debug,I, III,DIsDouble,VIsDouble,IMax,JMax,KMax,ZoneType,StrandID,ParentZn,IsBlock;
            INTEGER4 ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn, FileType;

            Debug     = 0;
            VIsDouble = 1;
            DIsDouble = 1;
            IMax      = nx;
            JMax      = ny;
            KMax      = nz;
            ZoneType  = 0;      /* Ordered */
            SolTime   = time_list[it];
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
                    "X Y Z rho  vx vy vz flag",
                    fullname,
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
            I = TECDAT112(&III,rho,&DIsDouble);
            I = TECDAT112(&III,u,&DIsDouble);
            I = TECDAT112(&III,v,&DIsDouble);
            I = TECDAT112(&III,w,&DIsDouble);
            I = TECDAT112(&III,flag,&DIsDouble);
            I = TECEND112();
        }
        else if (!strcmp(argv[1], "X"))
        {
            output_slice_X(n_slice);
        }
        else
        {
            fprintf(stderr, "Not implemented yet!\n");
            exit(-3);
        }

        free(u_p);
        free(v_p);
        free(w_p);
        free(rho_p);

    } //Time list

    free(time_list);
    free(flag);
    free(start_loc);
    free(u);
    free(v);
    free(w);
    free(X);
    free(Y);
    free(Z);
    free(rho);
    free(V);
}
