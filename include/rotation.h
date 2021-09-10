/** @file rotation.h Documented includes for harmonic module */

#ifndef __ROTATION__
#define __ROTATION__

#include "harmonic.h"

struct rotation {
	short has_rotated_cls;

	int has_tt;
	int has_te;
	int has_tb;
	int has_ee;
	int has_bb;
	int has_eb;
	int has_aa; /**< do we want \f$ C_l^{\alpha\alpha}\f$? (\f$ \alpha \f$ = rotation field) */

	int l_unrotated_max;
	int l_rotated_max;





};


#ifdef __cplusplus
extern "C" {
#endif

	int rotation_cl_at_l(
		struct rotation * pro,
		int l,
		double * cl_rotated
		);

	int rotation_init(
		struct precision * ppr,
		struct perturbations * ppt,
		struct harmonic * phr,
		struct fourier * pfo,
		struct rotation * pro
		);

	int rotation_free(struct rotation * pro
		);

	int rotation_indices(
		struct precision * ppr,
		struct harmonic * phr,
		struct rotation * pro
		);

	/* int rotation_rotated_cl_tt( */
	/* 	); */






#ifdef __cplusplus
}
#endif

#endif
