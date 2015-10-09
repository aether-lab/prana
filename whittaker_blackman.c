// whittaker_blackman7.c
#include <math.h>
#include <stdio.h>
#include "mex.h" /* Always include this */
void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
int nrhs, const mxArray *prhs[]) /* Input variables */
{
	
	// INPUTS
	
	// nlhs is the number of output (left-side) arguments, or the size of the plhs array.
	
	// plhs is the array of output arguments.
	
	// nrhs is the number of input (right hand side) arguments, or the size of the prhs array.
	
	// prhs is the array of input arguments.

    #define B_OUT plhs[0]

    #define A_IN prhs[0]
	
    #define X_IN prhs[1]

    #define Y_IN prhs[2]
	
	// Kernel radius
    #define rad_IN prhs[3]

	// Window type. If this variable is true, window the interpolation function
	// using the Blackman filter. If zero, don't window the function.
	#define window_flag prhs[4]

    double *B, *A, *X, *Y;
    int M, N, m, n, rad, h, k, L;
	int do_blackman;
    double eps = 1.0E-15;
    double pi = 3.141592653589793;
    double *Xfl, *Yfl;
	double apodization_norm_factor;


	// Get data from the matlab structures.
    M = mxGetM(A_IN); /* Get the dimensions of A */
    N = mxGetN(A_IN);
    A = mxGetPr(A_IN); /* Get the pointer to the data of A */
    X = mxGetPr(X_IN); /* Get the pointer to the data of X */
    Y = mxGetPr(Y_IN); /* Get the pointer to the data of Y */
    rad = mxGetScalar(rad_IN); /* Get the pointer to the data of rad */
	do_blackman = mxGetScalar(window_flag);
    B_OUT = mxCreateDoubleMatrix(M, N, mxREAL); /* Create the output matrix */
    B = mxGetPr(B_OUT); /* Get the pointer to the data of B */

	// Length of each side of the stencil
	L = 2 * rad + 1;

	// Normalization factor for the apodization function. 
	// If the Blackman window is used, then the apodization function
	// needs to be re-scaled by exactly 2.9584 so that its max value is one.
	if(do_blackman > 0) {
		 apodization_norm_factor = 2.9584;
	}
	else {
		 apodization_norm_factor = 1;
	}

    for(n = 0; n < N; n++) /* Compute a matrix with normalized columns */
    {
        for(m = 0; m < M; m++) {

            int rad1, rad2, rad3, rad4, xn_a, yn_a;
            double xn, yn, dx_a, dy_a, sincx_a, sincy_a;
			double apodization_x = 1;
			double apodization_y = 1;
            double buf = 0.0;

            xn = X[m + M*n]-1;
            yn = Y[m + M*n]-1;
            xn_a = floor(xn);
            yn_a = floor(yn);

            rad1 = xn_a-rad; if (rad1 < 0)     { rad1 = 0;   }
            rad2 = xn_a+rad; if (rad2 > N-1)   { rad2 = N-1; }
            rad3 = yn_a-rad; if (rad3 < 0)     { rad3 = 0;   }
            rad4 = yn_a+rad; if (rad4 > M-1)   { rad4 = M-1; }

            for (k = rad1; k <= rad2; k++) {
            dx_a        = (k-xn)*pi;
            sincx_a = sin(dx_a+eps)/(dx_a+eps);
			
			// Blackman apodization function (rows)
			// Only calculate this if apodization
			// has been specified.
			if(do_blackman > 0){
				
				// Blackman apodization function (columns)
				apodization_x = 0.42 + 0.50 * cos(dx_a / L) +
					0.80 * cos(2.0 * dx_a / L);
			}
			
            for (h = rad3; h<= rad4; h++) {
                dy_a   = (h-yn)*pi;
                sincy_a = sin(dy_a+eps)/(dy_a+eps);
				
				// Blackman apodization function (rows)
				// Only calculate this if apodization
				// has been specified.
				if(do_blackman > 0){
					apodization_y = 0.42 + 0.50 * cos(dy_a / L) + 
						0.80 * cos(2.0 * dy_a / L);
				}

				// Convolution of the kernel with the image data.
				buf += (sincx_a * apodization_x) * 
					(sincy_a * apodization_y) * (A[h + M * k]) / 					apodization_norm_factor;
				
                } // for h
            } // for k

            B[m + M*n] = buf;
        }
    }
    return;
}
