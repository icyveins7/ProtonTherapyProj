#include "mex.h"
#include<math.h>
#include<gsl/gsl_integration.h>
#define _USE_MATH_DEFINES

/*
 * e1_C.c
 *
 * e1 part of Dielectric function
 *
 * Input with Evalue in eV
 *
 * The calling syntax is:
 *
 *		out = e2(Evalue)
 *
 * This is a MEX-file for MATLAB, compiled with version 8.4.0.150421(R2014b),
 * using Microsoft Visual Studio 2013.
 *
 * If unable to load, download Visual C++ Redistributable Packages for Visual Studio 2013
 * http://www.microsoft.com/en-us/download/details.aspx?id=40784
 *
 * Command for build as follows:
 *
 * mex e1_C.c -v -largeArrayDims -I"F:\gsl\include" -L"F:\gsl\lib" -lgsl -lcblas -L"E:\MATLAB\R2014b\extern\lib\win64\microsoft" -llibmwblas -llibmwlapack LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:libcmt.lib"
 *
 * NOTE: Requires 64-Bit GSL libraries.
*/
//using the analytic integral from paper

struct my_f_params { double E; int j; };

double e1_exc(double E){
	double Ek[5] = { 8.17, 10.13, 11.31, 12.91, 14.5 };
	double yk[5] = { 1.62, 2.2, 2.1, 3.1, 3.9 };
	double fk[5] = { 0.0118, 0.023, 0.01675, 0.0285, 0.028 };

	double value = 0;
	for (int k = 0; k <= 4; k++){
		value += (fk[k] * (Ek[k] * Ek[k] - E * E)*((Ek[k] * Ek[k] - E*E)*(Ek[k] * Ek[k] - E*E) + 3 * yk[k] * yk[k] * E*E)) / (((Ek[k] * Ek[k] - E*E)*(Ek[k] * Ek[k] - E*E) + yk[k] * yk[k] * E*E)*((Ek[k] * Ek[k] - E*E)*(Ek[k] * Ek[k] - E*E) + yk[k] * yk[k] * E*E));
	}
	return 460.5316*value;
}

double e1_ion_integrand(double omega, void *p){
	struct my_f_params * params = (struct my_f_params *)p;
	double E = (params->E);
	int j = (params->j);

	double Ej[5] = { 11.95, 14.7, 16.6, 33.3, 540 };
	double yj[5] = { 12.5, 16.1, 19.4, 95, 220 };
	double trij[5] = { 1.16, 1.31, 0.55, 1, 1 };
	double fj[5] = { .225, .206, .161, .121, .179 };
	//double fj[5] = { .34, .31, .2394, .1595, .311 };



	double value = fj[j] * yj[j] / (2 * M_PI*((omega * omega - E*E)*(omega * omega - E*E)+yj[j]*yj[j]*E*E));
	if (yj[j] < (2 * omega)){
		value = value*(2 * E*log(fabs((omega + E) / (omega - E))) + M_PI * (omega * omega - E*E) / yj[j] - (omega * omega + E*E)*log((2 * omega + sqrt(4 * omega * omega - yj[j] * yj[j])) / (2 * omega - sqrt(4 * omega * omega - yj[j] * yj[j]))) / sqrt(4 * omega * omega - yj[j] * yj[j]));
	}
	else{
		value = value*(2 * E*log(fabs((omega + E) / (omega - E))) + M_PI * (omega * omega - E*E) / yj[j] - 2 * (omega * omega + E*E)*atan(sqrt(yj[j] * yj[j] - 4 * omega * omega) / (2 * omega)) / sqrt(yj[j] * yj[j] - 4 * omega * omega));
	}
	double gaussfunc = exp(-(omega - Ej[j])*(omega - Ej[j]) / (2 * trij[j] * trij[j]));
	return 460.5316*value*gaussfunc;

	
}

void e1(double Evalue, double *out){
    gsl_integration_workspace * w
			= gsl_integration_workspace_alloc(1000);
    
	*out = 0;
	double result, error;
	double trij[5] = { 1.16, 1.31, 0.55, 1, 1 };
	double Ej[5] = { 11.95, 14.7, 16.6, 33.3, 540 };

	for (int j = 0; j <= 4; j++){


		struct my_f_params alpha = { Evalue, j };

		gsl_function F;
		F.function = &e1_ion_integrand;
		F.params = &alpha;

		double pts[3] = { Ej[alpha.j] - trij[alpha.j], Ej[j], Ej[alpha.j] + trij[alpha.j] };

		//gsl_integration_qags(&F, Ej[alpha.j] - trij[alpha.j], Ej[alpha.j] + trij[alpha.j], 0, 100, 1000,
		//	w, &result, &error);

		gsl_integration_qagp(&F, &pts, 3, 0, 1, 1000, w, &result, &error);
		*out += result;
		
	}
	*out += e1_exc(Evalue);
	*out = 1 + *out;
    
    gsl_integration_workspace_free(w);
}

/* The gateway function */
void mexFunction( mwSize nlhs, mxArray *plhs[],
                  mwSize nrhs, const mxArray *prhs[])
{
    double Evalue;                                             
    double *out;            
    
    /* check for proper number of arguments */
    if(nrhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","1 inputs required.");
    }
//     if(nlhs!=1) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
//     }
    /* make sure the first input argument is scalar */
    if(  !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input a must be type double.");
    }
    
//     /* make sure the second input argument is type double */
//     if( !mxIsDouble(prhs[1]) || 
//          mxIsComplex(prhs[1])) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input b must be type double.");
//     }
//     
//     /* make sure the 3rd input argument is type double */
//     if( !mxIsDouble(prhs[2]) || 
//          mxIsComplex(prhs[2])) {
//         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input z must be type double.");
//     }
    
    /* get the value of the scalar inputs  */
    Evalue = mxGetScalar(prhs[0]);

    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleScalar(0);
    
    /* get a pointer to the real data in the output matrix */
    out = mxGetPr(plhs[0]);

    /* call the computational routine */
    e1(Evalue,out);
}