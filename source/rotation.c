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

	double * C_alpha; /* C_alpha[index_mu] */
	double * W_Ealpha; /* W_Ealpha[index_mu] */

	double * ksip = NULL; /* ksip[index_mu] */
	double * ksim = NULL; /* ksim[index_mu] */
	double * ksiX = NULL; /* ksiX[index_mu] */

	int num_nu,index_mu,icount;
	int l;
	double ll;




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
