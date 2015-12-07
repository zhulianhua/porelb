#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double qs_BCC[31]={
0.1000000e+01,
0.1575834e+01,
0.2483254e+01,
0.3233022e+01,
0.4022864e+01,
0.4650320e+01,
0.5281412e+01,
0.5826374e+01,
0.6258376e+01,
0.6544504e+01,
0.6878396e+01,
0.7190839e+01,
0.7268068e+01,
0.7304025e+01,
0.7301217e+01,
0.7236410e+01,
0.7298014e+01,
0.7369849e+01,
0.7109497e+01,
0.6228418e+01,
0.5235796e+01,
0.4476874e+01,
0.3541982e+01,
0.2939353e+01,
0.3935484e+01,
0.5179097e+01,
0.3959872e+01,
0.2227627e+01,
0.3393390e+01,
0.4491369e+01,
0.2200686e+01
};

int main(int argc, char *argv[])
{
    double Chi, K;
    int i, nx;
    double a, L;
    double  epsilon;
    double C, S;
    double d_star;
    double K1, K2;
    int ii;
           
    if(argc !=3)
    {
        fprintf(stderr, "argc err\n");
        printf("Usage: ./testpre Chi nx\n")
        printf("0<Chi<1 controls the size of the sphere, and nx is the resolution\n")
        exit(-1);
    }
    Chi = atof(argv[1]);
    nx = atoi(argv[2]);
    L = (double)(nx-1);
    /*L = 127.0;*/
    a = sqrt(3.0)*L*Chi/4.0; //sphere radius in lattice unit
    epsilon = 1.0 - 8.0/3*M_PI*a*a*a/L/L/L;
    C = 5;
    S = 2*4*M_PI*a*a/L/L/L;
    /*S = 2*M_PI/a;*/

    d_star = 0.0;
    for(i=0; i<31; i++)
        d_star += qs_BCC[i]*pow(Chi, i);
    K1 = 2.0/9*a*a/d_star/(1-epsilon);
    printf("K1=%e\n", K1);
    printf("K2=%e\n", epsilon*epsilon*epsilon/C/S/S);

#ifdef TABLE
    printf("Chi\tepsilon\tK1\tK2\n");
    for(ii=0; ii<31; ii++)
    {
        Chi = ii*(1.0/30);
        a = sqrt(3.0)*L*Chi/4.0; //sphere radius in lattice unit
        S = 2*4*M_PI*a*a/L/L/L;
        epsilon = 1.0 - 8.0/3*M_PI*a*a*a/L/L/L;


        d_star = 0.0;
        for(i=0; i<31; i++)
            d_star += qs_BCC[i]*pow(Chi, i);
        K1 = 2.0/9*a*a*epsilon/d_star/(1-epsilon);
        /*K1 = 1.0/d_star*epsilon*L*L*L/12/M_PI/a;*/
        K2 = epsilon*epsilon*epsilon/C/S/S;

        printf("%e\t%e\t%e\t%e\n", Chi, epsilon, K1, K2);
    }
#endif
    return 0;
}
