/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_EWALD_PME_GATHER_H
#define GMX_EWALD_PME_GATHER_H

#include "gromacs/utility/real.h"

#include "pme-internal.h"
#include "se.h"
#include "se_util.h"

#if GMX_DOUBLE==1
#if GMX_SIMD_X86_AVX_256
#include "se_int_sse_double.h"
#include "se_int_avx_256_double.h"
#else
#include "se_int_sse_double.h"
#endif  //AVX

#else  //SINGLE

#if GMX_SIMD_X86_AVX_256
#include "se_int_sse_single.h"
#include "se_int_avx_256_single.h"
#else
#include "se_int_sse_single.h"
#endif //AVX
#endif //DOUBLE


void SE_int_dispatch(rvec *force, real *grid, real *q,
                     splinedata_t *spline,
                     const SE_params *params, real scale,
                     gmx_bool bClearF);


void
gather_f_bsplines(struct gmx_pme_t *pme, real *grid,
                  gmx_bool bClearF, pme_atomcomm_t *atc,
                  splinedata_t *spline,
                  real scale,
		  SE_params *se_params,
		  gmx_bool       se_set
		  );

real
gather_energy_bsplines(struct gmx_pme_t *pme, real *grid,
                       pme_atomcomm_t *atc);

#endif
