#include "rotation.h"
#include <time.h>

int rotation_cl_at_l(
	struct rotation * pro,
	int l,
	double * cl_rotated
	) {
}

int rotation_init(
	struct precision * ppr,
	struct perturbations * ppt,
	struct harmoic * phr,
	struct fourier * pfo,
	struct rotation * pro
	) {

	double * mu; /* mu[index_mu]: discretized values of mu
                  between -1 and 1, roots of Legendre polynomial */
	double * w8; /* Corresponding Gauss-Legendre quadrature weights */
	double theta,delta_theta;

	double ** d02; /* dmn[index_mu][index_l] */
	double ** d22;
	double ** dm22;
	double * buf_dxx; /* buffer */

	double * C_a; /* C_a[index_mu] */
	double * W_Ea; /* W_Ea[index_mu] */

	double * ksip = NULL; /* ksip[index_mu] */
	double * ksim = NULL; /* ksim[index_mu] */
	double * ksiX = NULL; /* ksiX[index_mu] */

	double fac;

	int num_nu,index_mu,icount;
	int l;
	double ll;
	double * cl_lensed; /* cl_lensed[index_ct] */
	double * cl_ee = NULL; /* lensed  cl, to be filled to avoid repeated calls to harmonic_cl_at_l */
	double * cl_bb = NULL; /* lensed  cl, to be filled to avoid repeated calls to harmonic_cl_at_l */
	double * cl_aa; /* potential cl_aa, to be filled to avoid repeated calls to harmonic_cl_at_l */

	double res, resX, rot;
	double resp, resm, rotp, rotm;

	double ** cl_md_ic; /* array with argument
						   cl_md_ic[index_md][index_ic1_ic2*phr->ct_size+index_ct] */

	double ** cl_md;    /* array with argument
						   cl_md[index_md][index_ct] */

	int index_md;

	/* Timing */
	//double debut, fin;
	//double cpu_time;

	/** - check that we really want to compute at least one spectrum */
	

}


int rotation_free(
	struct rotation * pro
	) {

}

int rotaion_indices(
	struct precision * ppr,
	struct harmonic * phr,
	struct rotation * pro
	){


}


int rotation_rotated_cl_te(
	double *)
