#include <cstdlib>
#include <cstdio>
#include "porelb.h"
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkImageData.h>

int main(int argc, char* argv[])
{
    int num_p, num_f, nx, ny, nz, nt;
    int * time_list;
    nt = argc-1;
    printf("Total %d time seriese data\n", nt);
    time_list = (int*)calloc(nt, sizeof(int));

    for(int i=0; i<nt; i++)
        time_list[i] = atoi(argv[i+1]);


    FILE *fp_map;
    if((fp_map = fopen("V.bin", "rb")) == NULL)
    {
        fprintf(stderr, "Node map file openning error!\n");
        exit(-1);
    }
    fread(&num_f, sizeof(int), 1, fp_map);
    printf("[0]>>> Total num of fluid  nodes = %d\n", num_f);


    fread(&num_p,  sizeof(int), 1, fp_map);
    int *start_loc = (int*) calloc(num_p, sizeof(int));

    fread(&nx,  sizeof(int), 1, fp_map);
    fread(&ny,  sizeof(int), 1, fp_map);
    fread(&nz,  sizeof(int), 1, fp_map);

    long size = nx*ny*nz;

    double *u, *v, *w, *rho;
    unsigned char * flag;
    rho = (double*) calloc(size, sizeof(double));
    u   = (double*) calloc(size, sizeof(double));
    v   = (double*) calloc(size, sizeof(double));
    w   = (double*) calloc(size, sizeof(double));
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

    char partname[100];
    char fullname[100];
    FILE* fp_part;

    double *rho_p, *u_p, *v_p, *w_p;
    int N;

    //VTK stuff
    vtkSmartPointer<vtkImageData> imageData =
            vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(nx,ny,nz);
    imageData->AllocateScalars(VTK_DOUBLE, 5);
    vtkSmartPointer<vtkXMLImageDataWriter> writer =
        vtkSmartPointer<vtkXMLImageDataWriter>::New();


    for(int it = 0; it<nt; it++)
    {
        sprintf(fullname, "full_%07d.vti", time_list[it]);
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
        }

        int count = 0;
        for (int z = 0; z < nz; z++)
            for (int y = 0; y < ny; y++)
                for (int x = 0; x < nx; x++)
                {
                    double* pixel = static_cast<double*>(imageData->GetScalarPointer(x,y,z));
                    pixel[0] = rho[count];
                    pixel[1] = u[count];
                    pixel[2] = v[count];
                    pixel[3] = w[count];
                    pixel[4] = flag[count];
                    count++;
                }

        writer->SetFileName(fullname);
        writer->SetInputData(imageData);
        writer->Write();

        free(u_p);
        free(v_p);
        free(w_p);
        free(rho_p);
    }

    free(time_list);
    free(flag);
    free(start_loc);
    free(u);
    free(v);
    free(w);
    free(rho);
    free(V);
}
