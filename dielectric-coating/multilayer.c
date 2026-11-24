#include "clib.h"
#include <math.h>
#include <stdio.h>

#define PI 3.14159265359f

int main() {
    /* Refractive indices of air, dielectric layers, and glass substrate */
    float n0 = 1.0;
    float nL = 1.38;
    float ns = 1.5;
    float nH = 2.1;

    /* Physical thickness of layers */
    float dL = 100e-9;
    float dH = 66e-9;

    /* Optical thicknesses of layers */
    float hL = nL * dL;
    float hH = nH * dH;

    /* Optical admittances (Y = n at normal incidence),
        the constant factor sqrt(epsilon_0/mu_0) cancels in the
        reflectivity expression, so only refractive indices are used. */
    float Y0 = n0;
    float YL = nL;
    float YH = nH;
    float Ys = ns;

    /* User input of number of bilayers */
    int N;
    printf("Number of bilayers: ");
    scanf("%d", &N);

    /* Generates an output file name depending on the number of bilayers */
    char filename[50];
    sprintf(filename, "multi_spectrum_%d.dat", N);

    /* Opens a file for writing the reflectance
    spectrum for multi-layer dielectrics */
    FILE *fptr;
    fptr = fopen(filename, "w");

    if (!fptr) {
        printf("File error\n");
        return 1;
    }

    /* Loop over wavelength range (in nm) */
    for (int wl_nm = 300; wl_nm <= 800; wl_nm += 10) {
        float wl = wl_nm * 1e-9; /* Converting wavelength from nm to m */

        float k0 = 2.0 * PI / wl; /* Free space wavenumber */

        /* Phase thickness of each dielectric layer */
        float cL = k0 * hL;
        float cH = k0 * hH;

        /* Characteristic matrices for the individual layers */
        fcomplex ML[2][2], MH[2][2];

        /* Matrices for total system and intermediate multiplication */
        fcomplex M[2][2], Mtemp[2][2], pair[2][2];

        /* Low index layer matrix */
        ML[0][0] = complex(cos(cL), 0);
        ML[0][1] = cdiv(complex(0, sin(cL)), complex(YL, 0));
        ML[1][0] = cmult(complex(0, YL), complex(sin(cL), 0));
        ML[1][1] = complex(cos(cL), 0);

        /* High index layer matrix */
        MH[0][0] = complex(cos(cH), 0);
        MH[0][1] = cdiv(complex(0, sin(cH)), complex(YH, 0));
        MH[1][0] = cmult(complex(0, YH), complex(sin(cH), 0));
        MH[1][1] = complex(cos(cH), 0);

        /* Combining one bilayer pair (pair = M_L * M_H) */
        matmult(ML, MH, pair);

        /* Make the initial total matrix an identity matrix */
        M[0][0] = complex(1, 0);
        M[0][1] = complex(0, 0);
        M[1][0] = complex(0, 0);
        M[1][1] = complex(1, 0);

        /* Multiplying the bilayer matrix N times (M = (M_L * M_H)^N) */
        for (int i = 0; i < N; i++) {
            matmult(M, pair, Mtemp);
            matcopy(Mtemp, M);
        }

        /* Construct terms for reflection coefficient calculation */
        fcomplex a1 = rcmult(Y0, M[0][0]);
        fcomplex a2 = rcmult(Y0 * Ys, M[0][1]);
        fcomplex a3 = M[1][0];
        fcomplex a4 = rcmult(Ys, M[1][1]);

        /* Numerator and denominator of complex reflectivity r */
        fcomplex r1 = csub(cadd(a1, a2), cadd(a3, a4));
        fcomplex r2 = cadd(cadd(a1, a2), cadd(a3, a4));

        fcomplex r = cdiv(r1, r2); /* Complex reflection coefficient */

        /* Reflectance: R = |r|^2 */
        float R = ccabs(r);
        R = R * R;
        /* Writes the wavelength (nm) and reflectance of the
        multi-layer bilayer to the file */
        fprintf(fptr, "%d %f\n", wl_nm, R);
    }

    fclose(fptr); /* Closes output file */

    printf("multi_spectrum_N.dat file has been made\n");
    return (0);
}
