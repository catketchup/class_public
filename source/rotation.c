#include "rotation.h"
#include "lensing.h"
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
	struct lensing * ple,
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

	double * Ca; /* C_a[index_mu] */
	double * W_Ea; /* W_Ea[index_mu] */

	double * ksip = NULL; /* ksip[index_mu] */
	double * ksim = NULL; /* ksim[index_mu] */
	double * ksiX = NULL; /* ksiX[index_mu] */

	double fac;

	int num_nu,index_mu,icount;
	int l;
	double ll;
	double * cl_unrotated; /* cl_unrotated[index_ct] */
	double * cl_te = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or lensing_cl_at_l */
	double * cl_ee = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or lensing_cl_at_l */
	double * cl_bb = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or lensing_cl_at_l */
	double * cl_aa; /* potential cl_aa */

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

	if (pro->has_rotated_cls == _FALSE_) {
		if (pro->rotation_verbose > 0)
			printf("No rotation requested. Rotation module skipped.\n");
		return _SUCCESS_;
	}
	else {
		if (pro->rotation_verbose > 0) {
			printf("Computing rotated spectra ");
			if (ppr->accurate_rotation==_TRUE_)
				printf("(accurate mode)\n");
			else
				printf("(fast mode)\n");
		}
	}

	/** - initialize indices and allocate some of the arrays in the
		rotation structure */
	class_call(rotation_indices(ppr,phr,pro),
			   pro->error_message,
			   pro->error_message);



	/***************/
	/* configuring */
	/***************/

	/** - Compute \f$ d^l_{mm'} (\mu) \f$*/

	icount = 0;
	if(pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {
		class_alloc(d22,
				num_nu*sizeof(double*),
				pro->error_message);

		class_alloc(dm22,
				num_nu*sizeof(double*),
				pro->error_message);
		icount += 2*num_nu*(pro->l_unrotated_max+1);
	}

	/** - Allocate main contiguous buffer **/
	class_alloc(buf_dxx,
				icount * sizeof(double),
				pro->error_message);

	icount = 0;
	if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {
		for (index_mu=0; index_mu<num_nu; index_mu++){
			d22[index_mu] = &(buf_dxx[icount+index_mu            * (pro->l_unrotated_max+1)]);
			dm22[index_mu] = &(buf_dxx[icount+(index_mu+num_mu)  * (pro->l_unrotated_max+1)]);
		}
		icount +=2*num_nu*(pro->l_unrotated_max+1);
	}

	//debut = omp_get_wtime();
	if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {
		class_call(lensing_d22(mu,num_mu,pro->l_unrotated_max,d22),
               pro->error_message,
               pro->error_message);
		class_call(lensing_d2m2(mu,num_mu,pro->l_unrotated_max,dm22),
             pro->error_message,
             pro->error_message);
	}

	/** - compute \f$ Ca(\mu)\f$ */
	class_alloc(Ca,
				num_nu*sizeof(double),
				pro->error_message);

	/** - Locally store unrotated temperature \f$ cl\f$ and potential \f$ cl_{aa}\f$ spectra **/
	if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {
		class_alloc(cl_ee,
					  (pro->l_unrotated_max+1)*sizeof(double),
					  pro->error_message);
		class_alloc(cl_bb,
					(pro->l_unrotated_max+1)*sizeof(double),
					pro->error_message);
	}
	class_alloc(cl_aa,
				(pro->l_unrotated_max+1)*sizeof(double),
				pro->error_message);

	class_alloc(cl_md_ic,
				phr->md_size*sizeof(double *),
				pro->error_message);

	class_alloc(cl_md,
				phr->md_size*sizeof(double *),
				pro->error_message);

	for (index_md = 0; index_md < phr->md_size; index_md++) {

		if (phr->md_size > 1)

			class_alloc(cl_md[index_md],
						phr->ct_size*sizeof(double),
						pro->error_message);

		if (phr->ic_size[index_md] > 1)

			class_alloc(cl_md_ic[index_md],
						phr->ic_ic_size[index_md]*phr->ct_size*sizeof(double),
						pro->error_message);
	}

  for (l=2; l<=pro->l_unrotated_max; l++) {
	  class_call(harmonic_cl_at_l(phr,l,cl_unrotated,cl_md,cl_md_ic),
				 phr->error_message,
				 pro->error_message);
	  cl_aa[l] = cl_unrotated[pro->index_lt_aa];

	  if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {
		  cl_ee[l] = cl_unrotated[pro->index_lt_ee];
		  cl_bb[l] = cl_unrotated[pro->index_lt_bb];
	  }
  }

  for (index_md = 0; index_md < phr->md_size; index_md++) {

	  if (phr->md_size > 1)
		  free(cl_md[index_md]);

	  if (phr->ic_size[index_md] > 1)
		  free(cl_md_ic[index_md]);

  }

  free(cl_md_ic);
  free(cl_md);

  /** - Compute Ca\f$(\mu)\f$ **/

  //debut = omp_get_wtime();
#pragma omp parallel for                        \
	private (index_mu,l)						\
	schedule (static)
  for




}


int rotation_free(
	struct rotation * pro
	) {

}

int rotaion_indices(
	struct precision * ppr,
	struct harmonic * phr,
	struct lensing * ple,
	struct rotation * pro
	){

	int index_l;

	double ** cl_md_ic; /* array with argument
						   cl_md_ic[index_md][index_ic1_ic2*phr->ct_size+index_ct] */

	double ** cl_md;    /* array with argument
						   cl_md[index_md][index_ct] */

	int index_md;
	int index_lt;

	/* indices of all Cl types (lensed and unlensed) */

	if (phr->has_te == _TRUE_) {
		pro->has_te = _TRUE_;
		pro->index_lt_te=phr->index_ct_te;
	}
	else {
		pro->has_te==_FALSE_;
	}

	if (phr->has_ee == _TRUE_) {
		pro->has_ee = _TRUE_;
		pro->index_lt_ee=phr->index_ct_ee;
	}
	else {
		pro->has_ee = _FALSE_;
	}

	if (phr->has_bb == _TRUE_) {
		pro->has_bb = _TRUE_;
		pro->index_lt_bb=phr->index_ct_bb;
	}
	else {
		pro->has_bb = _FALSE_;
	}

	if (phr->has_aa == _TRUE_) {
		pro->has_aa = _TRUE_;
		pro->index_lt_aa=phr->index_ct_aa;
	}
	else {
		pro->has_aa = _FALSE_;
	}

	if (phr->has_ea == _TRUE_) {
		pro->has_ea = _TRUE_;
		pro->index_lt_ea=phr->index_ct_ea;
	}
	else {
		pro->has_ea = _FALSE_;
	}

	pro->lt_size = phr->l_max_tot;

	pro->l_rotated_max = pro->


}


int rotation_rotated_cl_te(
	double *)
