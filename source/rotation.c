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
	struct fourier * pfo,
	struct rotation * pro
	) {

	double * mu; /* mu[index_mu]: discretized values of mu
                  between -1 and 1, roots of Legendre polynomial */
	double * w8; /* Corresponding Gauss-Legendre quadrature weights */
	double theta,delta_theta;

	double ** d00; /* dmn[index_mu][index_l] */
	double ** d02;
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
	double * cl_tt = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or roting_cl_at_l */
	double * cl_te = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or roting_cl_at_l */
	double * cl_tb = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or roting_cl_at_l */
	double * cl_ee = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or roting_cl_at_l */
	double * cl_bb = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or roting_cl_at_l */
	double * cl_eb = NULL; /* unrotated cl, to be filled to avoid repeated calls to harmonic_cl_at_l or roting_cl_at_l */
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

	/** - put all precision variables hare; will be stored later in precision structure */
	if (ppr->accurate_rotation == _TRUE_) {
		num_mu=(pro->l_unrotated_max+ppr->num_mu_minus_lmax); /* Must be even ?? CHECK */
		num_mu += num_mu%2; /* Force it to be even */
	} else {
		/* Integrate correlation function difference on [0,pi/16] */
		num_mu = (pro->l_unrotated_max * 2 )/16;
	}

	/** - allocate array of \f$ \mu \f$ values, as well as quadrature weights */
	mu[num_mu-1] = 1.0;

	class_alloc(w8,
				(num_mu-1)*sizeof(double),
				pro->error_message);

	if (ppr->accurate_rotation == _TRUE_) {

		//debut = omp_get_wtime();
		class_call(quadrature_gauss_legendre(mu,
											 w8,
											 num_mu-1,
											 ppr->tol_gauss_legendre,
											 pro->error_message),
				   pro->error_message,
				   pro->error_message);
		//fin = omp_get_wtime();
		//cpu_time = (fin-debut);
		//printf("time in quadrature_gauss_legendre=%4.3f s\n",cpu_time);
	} else { /* Crude integration on [0,pi/16]: Riemann sum on theta */
		delta_theta = _PI_/16. / (double)(num_mu-1);
		for (index_mu=0;index_mu<num_mu-1;index_mu++) {
			theta = (index_mu+1)*delta_theta;
			mu[index_mu] = cos(theta);
			w8[index_mu] = sin(theta)*delta_theta; /* We integrate on mu */
		}
	}

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
	if(pro->has_eb==_TRUE_) {
		class_alloc(dm22,
				num_nu*sizeof(double*),
				pro->error_message);
		icount += num_nu*(pro->l_unrotated_max+1);
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
	if (pro->has_eb==_TRUE_) {
		for (index_mu=0; index_mu<num_nu; index_mu++){
			dm22[index_mu] = &(buf_dxx[icount+(index_mu+num_mu)  * (pro->l_unrotated_max+1)]);
		}
		icount +=num_nu*(pro->l_unrotated_max+1);
	}


	//debut = omp_get_wtime();
	class_call(lensing_d00(mu,num_mu,pro->l_unrotated_max,d00),
             pro->error_message,
             pro->error_message);

	if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {
		class_call(lensing_d22(mu,num_mu,pro->l_unrotated_max,d22),
               pro->error_message,
               pro->error_message);
		class_call(lensing_d2m2(mu,num_mu,pro->l_unrotated_max,dm22),
             pro->error_message,
             pro->error_message);
	}

	if (pro->has_eb==_TRUE_) {
		class_call(lensing_d2m2(mu,num_mu,pro->l_unrotated_max,dm22),
				   pro->error_message,
				   pro->error_message);
	}

	/** - compute \f$ Ca(\mu)\f$ */
	class_alloc(Ca,
				num_nu*sizeof(double),
				pro->error_message);

	/** - Locally store unrotated temperature \f$ cl\f$ and potential \f$ cl_{aa}\f$ spectra **/
	class_alloc(cl_tt,
				(pro->l_unrotated_max+1)*sizeof(double),
				pro->error_message);
	if (pro->has_te==_TRUE_ || pro->has_tb==_TRUE_) {
		class_alloc(cl_te,
					  (pro->l_unrotated_max+1)*sizeof(double),
					  pro->error_message);
	}
	if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {
		class_alloc(cl_ee,
					  (pro->l_unrotated_max+1)*sizeof(double),
					  pro->error_message);
		class_alloc(cl_bb,
					(pro->l_unrotated_max+1)*sizeof(double),
					pro->error_message);
	}
	if (pro->has_eb==_TRUE_) {
		class_alloc(cl_eb,
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
	  cl_tt[l] = cl_unrotated[pro->index_lt_tt];
	  /* generate cl_aa */
	  cl_aa[l] = 2.*_PI_*pro->A_cb*10E-5/(l*(l+1));
	  if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_ || pro->has_eb==_TRUE_) {
		  cl_ee[l] = cl_unrotated[pro->index_lt_ee];
		  cl_bb[l] = cl_unrotated[pro->index_lt_bb];
	  }

	  if (pro->has_te==_TRUE_ || pro->has_tb==_TRUE_) {
		  cl_te[l] = cl_unrotated[pro->index_lt_te];
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
  for (index_mu=0; index_mu<num_mu; index_mu++) {
	  Ca[index_mu] = 0;

	  for (l=2; l<=pro->l_unrotated_max; l++) {
		  Ca[index_mu] += (2.*l+1.)*cl_aa[l]*d00[index_mu];
	  }

	  Ca[index_mu] /= 4.*_PI_;
  }

  /** - compute ksi+, ksi-, ksiX for rotation*/

  /** - --> ksip, ksim for EE, BB **/
  if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {

	  class_calloc(ksi,
				   (num_mu-1),
				   sizeof(double),
				   pro->error_message);
  }

  /** - --> ksiX is for EB **/
  if (pro->has_eb==_TRUE_) {

	  class_calloc(ksiX,
				   (num_mu-1),
				   sizeof(double),
				   pro->error_message);
  }
  //debut = omp_get_wtime();
#pragma omp parallel for
  private (index_mu, l, ll, resp, resm, resX)  \
	  schedule (static)

	  for (index_mu;index_mu<num_mu-1;index_mu++) {
		  for (l=2;l<=pro->l_unrotated_max;l++) {

			  ll = (double)l;

			  fac = 2*ll+1;

			  if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {

				  resp = fac*d22[index_mu][l]*(cl_ee[l]+cl_bb[l]);
				  resm = fac*d22[index_mu][l]*(cl_ee[l]-cl_bb[l]);

				  ksip[index_mu] += resp;
				  ksim[index_mu] += resm;
			  }

			  if (pro->has_eb==_TRUE_) {
				  resX = fac*d22[index_mu][l]*(cl_ee[l]-cl_bb[l]);

				  ksiX[index_mu] += resX;
			  }
		  }

		  ksip[index_mu] *= exp(4*Ca[index_mu]);
		  ksim[index_mu] *= exp(-4*Ca[index_mu]);
		  ksiX[index_mu] *= exp(-4*Ca[index_mu]);

	  }
  //fin = omp_get_wtime();
  //cpu_time = (fin-debut);
  //printf("time in ksi=%4.3f s\n",cpu_time);

  /** - compute rotated \f$ C_l\f$'s by multiplation or integration */
  //debut = omp_get_wtime();
  if (pro->has_tt==_TRUE_) {
	  class_call(rotation_rotated_cl_tt(cl_tt, pro),
				 pro->error_message,
				 pro->error_message);

  }

  if (pro->has_te==_TRUE_) {
	  class_call(rotation_rotated_cl_te(cl_te, pro, pro->alpha, Ca[0]),
				 pro->error_message,
				 pro->error_message);
  }

  if (pro->has_tb==_TRUE_) {
	  class_call(rotation_rotated_cl_tb(cl_te, pro, pro->alpha, Ca[0]),
				 pro->error_message,
				 pro->error_message);
  }

  if (pro->has_ee==_TRUE_ || pro->has_bb==_TRUE_) {


  }


}

/**
 * This routine frees all the memory space allocated by rotation_init().
 *
 * To be called at the end of each run, only when no further calls to
 * rotation_cl_at_l() are needed.
 *
 * @param ple Input: pointer to rotation structure (which fields must be freed)
 * @return the error status
 */

int rotation_free(
                 struct rotation * pro
                 ) {

  if (pro->has_rotated_cls == _TRUE_) {

    free(pro->l);
    free(pro->cl_rot);
    free(pro->ddcl_rot);
    free(pro->l_max_lt);

  }

  return _SUCCESS_;

}

/* rotation does not change Cl_TT */
int rotation_rotated_cl_tt(double *cl_tt,
						   struct rotation * pro){
	int index_l;
	for(index_l=0; index_l<pro->l_size; index_l++){
		pro->cl_rotated[index_l*pro->lt_size+pro->index_lt_tt] = cl_tt[index_l];
	}

	return _SUCCESS_;
}

int rotation_rotated_cl_te(double *cl_te,
						   struct rotation * pro,
						   double alpha,
						   double Ca0){
	int index_l;
	for(index_l=0; index_l<pro->l_size; index_l++){
		pro->cl_rotated[index_l*pro->lt_size+pro->index_lt_te] = cl_te[index_l]*cos(2*pro->alpha)*exp(-2*Ca0);
	}

	return _SUCCESS_;
}

int rotation_rotated_cl_tb(double *cl_te,
						   struct rotation * pro,
						   double alpha,
						   double Ca0){
	int index_l;
	for(index_l=0; index_l<pro->l_size; index_l++){
		pro->cl_rotated[index_l*pro->lt_size+pro->index_lt_te] = cl_te[index_l]*sin(2*pro->alpha)*exp(-2*Ca0);
	}

	return _SUCCESS_;
}

/**
 * This routine computes the rotated power spectra by Gaussian quadrature
 *
 * @param ksip  Input: rotated correlation function (ksip[index_mu])
 * @param d22  Input: Legendre polynomials (\f$ d^l_{22}\f$[l][index_mu])
 * @param w8   Input: Legendre quadrature weights (w8[index_mu])
 * @param nmu  Input: Number of quadrature points (0<=index_mu<=nmu)
 * @param pro  Input/output: Pointer to the rotation structure
 * @return the error status
 */

int rotation_rotated_cl_ee(
	double *ksip,
	double **d22,
	double *w8,
	int nmu,
	struct rotation *pro
	){

}

int rotation_rotated_cl_bb(){

}

int rotation_rotated_cl_eb(){

}


int rotation_indices(
	struct precision * ppr,
	struct harmonic * phr,
	struct rotation * pro
	){

	int index_l;

	double ** cl_md_ic; /* array with argument
						   cl_md_ic[index_md][index_ic1_ic2*phr->ct_size+index_ct] */

	double ** cl_md;    /* array with argument
						   cl_md[index_md][index_ct] */

	int index_md;
	int index_lt;

	/* indices of all Cl types (rotated and unrotated) */
	index_ct=0;

	if (phr->has_tt == _TRUE_) {
		pro->has_tt = _TRUE_;
		pro->index_lt_tt=index_ct;
		index_ct++;
	}
	else {
		pro->has_tt==_FALSE_;
	}

	if ((phr->has_ee == _TRUE_) && (phr->has_bb == _TRUE_)) {
		pro->has_ee = _TRUE_;
		pro->has_bb = _TRUE_;
		pro->has_eb = _TRUE_;
		pro->index_lt_ee=index_ct;
		index_ct++;
		pro->index_lt_bb=index_ct;
		index_ct++;
		pro->index_lt_eb=index_ct+1;


	}
	else {
		pro->has_ee = _FALSE_;
		pro->has_bb = _FALSE_;
		pro->has_eb = _FALSE_;
	}

	if (phr->has_te == _TRUE_) {
		pro->has_te = _TRUE_;
		pro->has_tb = _TRUE_;
		pro->index_lt_te=index_ct;
		index_ct++;
		pro->index_lt_tb=index_ct+1;
	}
	else {
		pro->has_te = _FALSE_;
		pro->has_tb = _FALSE_;
	}

	/* if (phr->has_aa == _TRUE_) { */
	/* 	pro->has_aa = _TRUE_; */
	/* 	pro->index_lt_aa=phr->index_ct_aa; */
	/* } */
	/* else { */
	/* 	pro->has_aa = _FALSE_; */
	/* } */

	/* if (phr->has_ea == _TRUE_) { */
	/* 	pro->has_ea = _TRUE_; */
	/* 	pro->index_lt_ea=phr->index_ct_ea; */
	/* } */
	/* else { */
	/* 	pro->has_ea = _FALSE_; */
	/* } */

	pro->lt_size = phr->ct_size;

	/* number of multipoles */

	pro->l_unrotated_max = phr->l_max_tot;

	pro->l_rotated_max = pro->l_unrotated_max - ppr->delta_l_max;

	for (index_l=0; (index_l < phr->l_size_max) && (phr->l[index_l] <= pro->l_rotated_max);index_l++);

	if (index_l < phr->l_size_max) index_l++; /* one more point in order to be able to interpolate till pro->l_rotated_max */

	pro->l_size = index_l+1;

	class_alloc(pro->, pro->l_size*sizeof(double),pro->error_message);

	for (index_l=0; index_l < pro->l_size; index_l++) {

		pro->l[index_l] = phr->l[index_l];

	}

	/* allocate table where results will be stored */
	class_alloc(pro->cl_rot,
				pro->l_size*pro->lt_size*sizeof(double),
				pro->error_message);

	class_alloc(pro->ddcl_rot,
				pro->l_size*pro->lt_size*sizeof(double),
				pro->error_message);

	/* fill with unrotated cls */

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

	for (index_l=0; index_l<pro->l_size; index_l++) {

		class_call(harmonic_cl_at_l(phr,pro->l[index_l],&(pro->cl_rot[index_l*pro->lt_size]),cl_md,cl_md_ic),
				   phr->error_message,
				   pro->error_message);

	}

	for (index_md = 0; index_md < phr->md_size; index_md++) {

		if (phr->md_size > 1)
			free(cl_md[index_md]);

		if (phr->ic_size[index_md] > 1)
			free(cl_md_ic[index_md]);

	}

	free(cl_md_ic);
	free(cl_md);

	/* we want to output Cl_rotated up to the same l_max as Cl_unroated
	   (even if a number delta_l_max of extra values of l have been used
	   internally for more accurate results). Notable exception to the
	   above rule: ClBB_rotated(scalars) must be outputed at least up to the same l_max as
	   ClEE_unrotated(scalars) (since ClBB_unroated is null for scalars)
	*/

	class_alloc(pro->l_max_lt,pro->lt_size*sizeof(double),pro->error_message);
	for (index_lt = 0; index_lt < pro->lt_size; index_lt++) {
		pro->l_max_lt[index_lt]=0.;
		for (index_md = 0; index_md < phr->md_size; index_md++) {
			pro->l_max_lt[index_lt]=MAX(pro->l_max_lt[index_lt],phr->l_max_ct[index_md][index_lt]);

			if ((pro->has_bb == _TRUE_) && (pro->has_ee == _TRUE_) && (index_lt == pro->index_lt_bb)) {
				pro->l_max_lt[index_lt]=MAX(pro->l_max_lt[index_lt],phr->l_max_ct[index_md][pro->index_lt_ee]);
			}

		}
	}

	return _SUCCESS_;

}


int rotation_rotated_cl_te(
	double *)
