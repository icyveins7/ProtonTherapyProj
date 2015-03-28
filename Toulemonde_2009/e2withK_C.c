#include "mex.h"
#include<math.h>
#include<gsl/gsl_integration.h>

/*
 * e2_C.c
 *
 * e2 part of Dielectric function
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
 * mex e2_C.c -v -largeArrayDims -I"F:\gsl\include" -L"F:\gsl\lib" -lgsl -lcblas -L"E:\MATLAB\R2014b\extern\lib\win64\microsoft" -llibmwblas -llibmwlapack LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:libcmt.lib"
 *
 * NOTE: Requires 64-Bit GSL libraries.
*/
struct my_f_params { double E; int j; double Kvalue};

double e2_ion_integrand(double omega, void *p){

	struct my_f_params * params = (struct my_f_params *)p;
	double E = (params->E);
	int j = (params->j);
    double Kvalue = (params->Kvalue);

	double Ej[5] = { 11.95, 14.7, 16.6, 33.3, 540 };
	double yj[5] = { 12.5, 16.1, 19.4, 95, 220 };
	double trij[5] = { 1.16, 1.31, 0.55, 1, 1 };
	double fj[5] = {.225, .206, .161, .121, .179};
    
    double fk[5] = { 0.0118, 0.023, 0.01675, 0.0285, 0.028 };
    double aj[5] = {3.82, 2.47, 2.47, 3.01, 2.44};
    double bj[5] = {0.0272, 0.0295, 0.0311, 0.0111, 0.0633};
    double cj[5] = {0.098, 0.075, 0.074, 0.765, 0.425};
    
    double fk_K[5];
    double fk_Ksum=0;
    double fksum=0;
    
    for (int i=0;i<=4;i++){
        fk_K[i]=fk[i]*(exp(-ak[i]*Kvalue*Kvalue) + bk[i]*Kvalue*Kvalue*exp(-ck[i]*Kvalue*Kvalue) );
        fk_Ksum+=fk_K[i];
        fksum+=fk[i];
    }
    
    for (int l=0;l<=4;l++){
        fj[l] = fj * (10 - fk_ksum) / (10 - fksum) ;
    }
    

	if (omega <= E){
		double Dlist = fj[j] * yj[j] * E / ((omega*omega - E*E)*(omega*omega - E*E) + yj[j] * yj[j] * E*E);
		double gaussfunc = exp(-(omega - Ej[j])*(omega - Ej[j]) / (2 * trij[j] * trij[j]));
		return Dlist*gaussfunc;
	}
	else{
		return 0;
	}
}

double e2_exc(double E, double Kvalue){
	double Ek[5] = { 8.17, 10.13, 11.31, 12.91, 14.5 };
	double yk[5] = { 1.62, 2.2, 2.1, 3.1, 3.9 };
	double fk[5] = { 0.0118, 0.023, 0.01675, 0.0285, 0.028 };
    double ak[5] = {3.82, 2.47, 2.47, 3.01, 2.44};
    double bk[5] = {0.0272, 0.0295, 0.0311, 0.0111, 0.0633};
    double ck[5] = {0.098, 0.075, 0.074, 0.765, 0.425};
    
    for (int i=0;i<=4;i++){fk[i]=fk[i]*(exp(-ak[i]*Kvalue*Kvalue) + bk[i]*Kvalue*Kvalue*exp(-ck[i]*Kvalue*Kvalue) );}

	double Dstar = 0;
	for (int k = 0; k <= 4; k++){
		Dstar += 2.*fk[k] * yk[k] * yk[k] * yk[k] * E*E*E / pow((Ek[k] * Ek[k] - E*E)*(Ek[k] * Ek[k] - E*E) + yk[k] * yk[k] * E*E, 2);
	}

	return Dstar;
}


void e2(double Evalue, double Kvalue, double *out){
	//current term
	gsl_integration_workspace * w
		= gsl_integration_workspace_alloc(1000);

	double result, error;
	double trij[5] = { 1.16, 1.31, 0.55, 1, 1 };
	double Ej[5] = { 11.95, 14.7, 16.6, 33.3, 540 };
    
    for (int j = 0;j<=4;j++){
        struct my_f_params alpha = { Evalue, j, Kvalue };

        gsl_function F;
        F.function = &e2_ion_integrand;
        F.params = &alpha;

        gsl_integration_qags(&F, Ej[alpha.j]-trij[alpha.j], Ej[alpha.j]+trij[alpha.j], 0, 1e-7, 1000,
            w, &result, &error);

        *out+=result;
    }
    *out+=e2_exc(Evalue);
    *out=*out*460.5316;
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
    e2(Evalue,out);
}