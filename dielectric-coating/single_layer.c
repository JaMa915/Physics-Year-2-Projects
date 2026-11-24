#include "clib.h"
#include <math.h>
#include <stdio.h>

#define PI 3.14159265359f

int main() {
    /* Refractive indices of air, dielectric layer, and glass substrate */
    float n0 = 1.0;
    float nL = 1.38;
    float ns = 1.5;

    /* Physical and optical thickness of the layer */
    float dL = 100e-9;
    float hL = nL * dL;

    /* Optical admittances (Y = n at normal incidence),
        the constant factor sqrt(epsilon_0/mu_0) cancels in the
        reflectivity expression, so only refractive indices are used. */
    float Y0 = n0;
    float YL = nL;
    float Ys = ns;

    /* Opens file to store computed reflectance spectrum */
    FILE *fptr;
    fptr = fopen("spectrum.dat", "w");

    if (!fptr) {
        printf("File error\n");
        return 1;
    }

    /* Loop over wavelength range (in nm) */
    for (int wl_nm = 300; wl_nm <= 800; wl_nm += 10) {
        float wl = wl_nm * 1e-9; /* Converting wavelength from nm to m */

        float k0 = 2.0 * PI / wl; /* Free space wavenumber */

        float c = k0 * hL; /* Phase thickness of the dielectric layer */

        /* Characteristic matrix of the dielectric layer,
relating electric and magnetic fields across the layer */
        fcomplex M[2][2];
        M[0][0] = complex(cos(c), 0);
        M[0][1] = cdiv(complex(0, sin(c)), complex(YL, 0));
        M[1][0] = cmult(complex(0, YL), complex(sin(c), 0));
        M[1][1] = complex(cos(c), 0);

        /* Construct terms for reflection coefficient
    using transfer matrix elements */
        fcomplex a1 = rcmult(Y0, M[0][0]);
        fcomplex a2 = rcmult(Y0 * Ys, M[0][1]);
        fcomplex a3 = M[1][0];
        fcomplex a4 = rcmult(Ys, M[1][1]);

        /* Numerator and denominator of the complex reflectivity r */
        fcomplex r1 = csub(cadd(a1, a2), cadd(a3, a4));
        fcomplex r2 = cadd(cadd(a1, a2), cadd(a3, a4));

        fcomplex r = cdiv(r1, r2); /* Complex reflection coefficient */

        /* Reflectance: R = |r|^2 */
        float R = ccabs(r);
        R = R * R;

        /* Writes the wavelength (nm) and reflectance to the file */
        fprintf(fptr, "%d %f\n", wl_nm, R);
    }

    fclose(fptr); /* Closes the output file */

    printf("spectrum.dat file has been made\n");
    return (0);
}
