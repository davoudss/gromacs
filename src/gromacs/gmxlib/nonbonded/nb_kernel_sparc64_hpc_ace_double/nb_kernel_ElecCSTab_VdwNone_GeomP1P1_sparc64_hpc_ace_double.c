/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2017, by the GROMACS development team, led by
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
/*
 * Note: this file was generated by the GROMACS sparc64_hpc_ace_double kernel generator.
 */
#include "gmxpre.h"

#include "config.h"

#include <math.h>

#include "../nb_kernel.h"
#include "gromacs/gmxlib/nrnb.h"

#include "kernelutil_sparc64_hpc_ace_double.h"

/*
 * Gromacs nonbonded kernel:   nb_kernel_ElecCSTab_VdwNone_GeomP1P1_VF_sparc64_hpc_ace_double
 * Electrostatics interaction: CubicSplineTable
 * VdW interaction:            None
 * Geometry:                   Particle-Particle
 * Calculate force/pot:        PotentialAndForce
 */
void
nb_kernel_ElecCSTab_VdwNone_GeomP1P1_VF_sparc64_hpc_ace_double
                    (t_nblist                    * gmx_restrict       nlist,
                     rvec                        * gmx_restrict          xx,
                     rvec                        * gmx_restrict          ff,
                     struct t_forcerec           * gmx_restrict          fr,
                     t_mdatoms                   * gmx_restrict     mdatoms,
                     nb_kernel_data_t gmx_unused * gmx_restrict kernel_data,
                     t_nrnb                      * gmx_restrict        nrnb)
{
    /* Suffixes 0,1,2,3 refer to particle indices for waters in the inner or outer loop, or
     * just 0 for non-waters.
     * Suffixes A,B refer to j loop unrolling done with double precision SIMD, e.g. for the two different
     * jnr indices corresponding to data put in the four positions in the SIMD register.
     */
    int              i_shift_offset,i_coord_offset,outeriter,inneriter;
    int              j_index_start,j_index_end,jidx,nri,inr,ggid,iidx;
    int              jnrA,jnrB;
    int              j_coord_offsetA,j_coord_offsetB;
    int              *iinr,*jindex,*jjnr,*shiftidx,*gid;
    real             rcutoff_scalar;
    real             *shiftvec,*fshift,*x,*f;
    _fjsp_v2r8       tx,ty,tz,fscal,rcutoff,rcutoff2,jidxall;
    int              vdwioffset0;
    _fjsp_v2r8       ix0,iy0,iz0,fix0,fiy0,fiz0,iq0,isai0;
    int              vdwjidx0A,vdwjidx0B;
    _fjsp_v2r8       jx0,jy0,jz0,fjx0,fjy0,fjz0,jq0,isaj0;
    _fjsp_v2r8       dx00,dy00,dz00,rsq00,rinv00,rinvsq00,r00,qq00,c6_00,c12_00;
    _fjsp_v2r8       velec,felec,velecsum,facel,crf,krf,krf2;
    real             *charge;
    _fjsp_v2r8       rt,vfeps,vftabscale,Y,F,G,H,Heps,Fp,VV,FF,twovfeps;
    real             *vftab;
    _fjsp_v2r8       itab_tmp;
    _fjsp_v2r8       dummy_mask,cutoff_mask;
    _fjsp_v2r8       one     = gmx_fjsp_set1_v2r8(1.0);
    _fjsp_v2r8       two     = gmx_fjsp_set1_v2r8(2.0);
    union { _fjsp_v2r8 simd; long long int i[2]; } vfconv,gbconv,ewconv;

    x                = xx[0];
    f                = ff[0];

    nri              = nlist->nri;
    iinr             = nlist->iinr;
    jindex           = nlist->jindex;
    jjnr             = nlist->jjnr;
    shiftidx         = nlist->shift;
    gid              = nlist->gid;
    shiftvec         = fr->shift_vec[0];
    fshift           = fr->fshift[0];
    facel            = gmx_fjsp_set1_v2r8(fr->epsfac);
    charge           = mdatoms->chargeA;

    vftab            = kernel_data->table_elec->data;
    vftabscale       = gmx_fjsp_set1_v2r8(kernel_data->table_elec->scale);

    /* Avoid stupid compiler warnings */
    jnrA = jnrB = 0;
    j_coord_offsetA = 0;
    j_coord_offsetB = 0;

    outeriter        = 0;
    inneriter        = 0;

    /* Start outer loop over neighborlists */
    for(iidx=0; iidx<nri; iidx++)
    {
        /* Load shift vector for this list */
        i_shift_offset   = DIM*shiftidx[iidx];

        /* Load limits for loop over neighbors */
        j_index_start    = jindex[iidx];
        j_index_end      = jindex[iidx+1];

        /* Get outer coordinate index */
        inr              = iinr[iidx];
        i_coord_offset   = DIM*inr;

        /* Load i particle coords and add shift vector */
        gmx_fjsp_load_shift_and_1rvec_broadcast_v2r8(shiftvec+i_shift_offset,x+i_coord_offset,&ix0,&iy0,&iz0);

        fix0             = _fjsp_setzero_v2r8();
        fiy0             = _fjsp_setzero_v2r8();
        fiz0             = _fjsp_setzero_v2r8();

        /* Load parameters for i particles */
        iq0              = _fjsp_mul_v2r8(facel,gmx_fjsp_load1_v2r8(charge+inr+0));

        /* Reset potential sums */
        velecsum         = _fjsp_setzero_v2r8();

        /* Start inner kernel loop */
        for(jidx=j_index_start; jidx<j_index_end-1; jidx+=2)
        {

            /* Get j neighbor index, and coordinate index */
            jnrA             = jjnr[jidx];
            jnrB             = jjnr[jidx+1];
            j_coord_offsetA  = DIM*jnrA;
            j_coord_offsetB  = DIM*jnrB;

            /* load j atom coordinates */
            gmx_fjsp_load_1rvec_2ptr_swizzle_v2r8(x+j_coord_offsetA,x+j_coord_offsetB,
                                              &jx0,&jy0,&jz0);

            /* Calculate displacement vector */
            dx00             = _fjsp_sub_v2r8(ix0,jx0);
            dy00             = _fjsp_sub_v2r8(iy0,jy0);
            dz00             = _fjsp_sub_v2r8(iz0,jz0);

            /* Calculate squared distance and things based on it */
            rsq00            = gmx_fjsp_calc_rsq_v2r8(dx00,dy00,dz00);

            rinv00           = gmx_fjsp_invsqrt_v2r8(rsq00);

            /* Load parameters for j particles */
            jq0              = gmx_fjsp_load_2real_swizzle_v2r8(charge+jnrA+0,charge+jnrB+0);

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            r00              = _fjsp_mul_v2r8(rsq00,rinv00);

            /* Compute parameters for interactions between i and j atoms */
            qq00             = _fjsp_mul_v2r8(iq0,jq0);

            /* Calculate table index by multiplying r with table scale and truncate to integer */
            rt               = _fjsp_mul_v2r8(r00,vftabscale);
            itab_tmp         = _fjsp_dtox_v2r8(rt);
            vfeps            = _fjsp_sub_v2r8(rt, _fjsp_xtod_v2r8(itab_tmp));
            twovfeps         = _fjsp_add_v2r8(vfeps,vfeps);
            _fjsp_store_v2r8(&vfconv.simd,itab_tmp);

            vfconv.i[0]     *= 4;
            vfconv.i[1]     *= 4;

            /* CUBIC SPLINE TABLE ELECTROSTATICS */
            Y                = _fjsp_load_v2r8( vftab + vfconv.i[0] );
            F                = _fjsp_load_v2r8( vftab + vfconv.i[1] );
            GMX_FJSP_TRANSPOSE2_V2R8(Y,F);
            G                = _fjsp_load_v2r8( vftab + vfconv.i[0] +2);
            H                = _fjsp_load_v2r8( vftab + vfconv.i[1] +2);
            GMX_FJSP_TRANSPOSE2_V2R8(G,H);
            Fp               = _fjsp_madd_v2r8(vfeps,_fjsp_madd_v2r8(vfeps,H,G),F);
            VV               = _fjsp_madd_v2r8(vfeps,Fp,Y);
            velec            = _fjsp_mul_v2r8(qq00,VV);
            FF               = _fjsp_madd_v2r8(_fjsp_madd_v2r8(twovfeps,H,G),vfeps,Fp);
            felec            = _fjsp_neg_v2r8(_fjsp_mul_v2r8(_fjsp_mul_v2r8(qq00,FF),_fjsp_mul_v2r8(vftabscale,rinv00)));

            /* Update potential sum for this i atom from the interaction with this j atom. */
            velecsum         = _fjsp_add_v2r8(velecsum,velec);

            fscal            = felec;

            /* Update vectorial force */
            fix0             = _fjsp_madd_v2r8(dx00,fscal,fix0);
            fiy0             = _fjsp_madd_v2r8(dy00,fscal,fiy0);
            fiz0             = _fjsp_madd_v2r8(dz00,fscal,fiz0);
            
            gmx_fjsp_decrement_fma_1rvec_2ptr_swizzle_v2r8(f+j_coord_offsetA,f+j_coord_offsetB,fscal,dx00,dy00,dz00);

            /* Inner loop uses 46 flops */
        }

        if(jidx<j_index_end)
        {

            jnrA             = jjnr[jidx];
            j_coord_offsetA  = DIM*jnrA;

            /* load j atom coordinates */
            gmx_fjsp_load_1rvec_1ptr_swizzle_v2r8(x+j_coord_offsetA,
                                              &jx0,&jy0,&jz0);

            /* Calculate displacement vector */
            dx00             = _fjsp_sub_v2r8(ix0,jx0);
            dy00             = _fjsp_sub_v2r8(iy0,jy0);
            dz00             = _fjsp_sub_v2r8(iz0,jz0);

            /* Calculate squared distance and things based on it */
            rsq00            = gmx_fjsp_calc_rsq_v2r8(dx00,dy00,dz00);

            rinv00           = gmx_fjsp_invsqrt_v2r8(rsq00);

            /* Load parameters for j particles */
            jq0              = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(),charge+jnrA+0);

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            r00              = _fjsp_mul_v2r8(rsq00,rinv00);

            /* Compute parameters for interactions between i and j atoms */
            qq00             = _fjsp_mul_v2r8(iq0,jq0);

            /* Calculate table index by multiplying r with table scale and truncate to integer */
            rt               = _fjsp_mul_v2r8(r00,vftabscale);
            itab_tmp         = _fjsp_dtox_v2r8(rt);
            vfeps            = _fjsp_sub_v2r8(rt, _fjsp_xtod_v2r8(itab_tmp));
            twovfeps         = _fjsp_add_v2r8(vfeps,vfeps);
            _fjsp_store_v2r8(&vfconv.simd,itab_tmp);

            vfconv.i[0]     *= 4;
            vfconv.i[1]     *= 4;

            /* CUBIC SPLINE TABLE ELECTROSTATICS */
            Y                = _fjsp_load_v2r8( vftab + vfconv.i[0] );
            F                = _fjsp_setzero_v2r8();
            GMX_FJSP_TRANSPOSE2_V2R8(Y,F);
            G                = _fjsp_load_v2r8( vftab + vfconv.i[0] +2);
            H                = _fjsp_setzero_v2r8();
            GMX_FJSP_TRANSPOSE2_V2R8(G,H);
            Fp               = _fjsp_madd_v2r8(vfeps,_fjsp_madd_v2r8(vfeps,H,G),F);
            VV               = _fjsp_madd_v2r8(vfeps,Fp,Y);
            velec            = _fjsp_mul_v2r8(qq00,VV);
            FF               = _fjsp_madd_v2r8(_fjsp_madd_v2r8(twovfeps,H,G),vfeps,Fp);
            felec            = _fjsp_neg_v2r8(_fjsp_mul_v2r8(_fjsp_mul_v2r8(qq00,FF),_fjsp_mul_v2r8(vftabscale,rinv00)));

            /* Update potential sum for this i atom from the interaction with this j atom. */
            velec            = _fjsp_unpacklo_v2r8(velec,_fjsp_setzero_v2r8());
            velecsum         = _fjsp_add_v2r8(velecsum,velec);

            fscal            = felec;

            fscal            = _fjsp_unpacklo_v2r8(fscal,_fjsp_setzero_v2r8());

            /* Update vectorial force */
            fix0             = _fjsp_madd_v2r8(dx00,fscal,fix0);
            fiy0             = _fjsp_madd_v2r8(dy00,fscal,fiy0);
            fiz0             = _fjsp_madd_v2r8(dz00,fscal,fiz0);
            
            gmx_fjsp_decrement_fma_1rvec_1ptr_swizzle_v2r8(f+j_coord_offsetA,fscal,dx00,dy00,dz00);

            /* Inner loop uses 46 flops */
        }

        /* End of innermost loop */

        gmx_fjsp_update_iforce_1atom_swizzle_v2r8(fix0,fiy0,fiz0,
                                              f+i_coord_offset,fshift+i_shift_offset);

        ggid                        = gid[iidx];
        /* Update potential energies */
        gmx_fjsp_update_1pot_v2r8(velecsum,kernel_data->energygrp_elec+ggid);

        /* Increment number of inner iterations */
        inneriter                  += j_index_end - j_index_start;

        /* Outer loop uses 8 flops */
    }

    /* Increment number of outer iterations */
    outeriter        += nri;

    /* Update outer/inner flops */

    inc_nrnb(nrnb,eNR_NBKERNEL_ELEC_VF,outeriter*8 + inneriter*46);
}
/*
 * Gromacs nonbonded kernel:   nb_kernel_ElecCSTab_VdwNone_GeomP1P1_F_sparc64_hpc_ace_double
 * Electrostatics interaction: CubicSplineTable
 * VdW interaction:            None
 * Geometry:                   Particle-Particle
 * Calculate force/pot:        Force
 */
void
nb_kernel_ElecCSTab_VdwNone_GeomP1P1_F_sparc64_hpc_ace_double
                    (t_nblist                    * gmx_restrict       nlist,
                     rvec                        * gmx_restrict          xx,
                     rvec                        * gmx_restrict          ff,
                     struct t_forcerec           * gmx_restrict          fr,
                     t_mdatoms                   * gmx_restrict     mdatoms,
                     nb_kernel_data_t gmx_unused * gmx_restrict kernel_data,
                     t_nrnb                      * gmx_restrict        nrnb)
{
    /* Suffixes 0,1,2,3 refer to particle indices for waters in the inner or outer loop, or
     * just 0 for non-waters.
     * Suffixes A,B refer to j loop unrolling done with double precision SIMD, e.g. for the two different
     * jnr indices corresponding to data put in the four positions in the SIMD register.
     */
    int              i_shift_offset,i_coord_offset,outeriter,inneriter;
    int              j_index_start,j_index_end,jidx,nri,inr,ggid,iidx;
    int              jnrA,jnrB;
    int              j_coord_offsetA,j_coord_offsetB;
    int              *iinr,*jindex,*jjnr,*shiftidx,*gid;
    real             rcutoff_scalar;
    real             *shiftvec,*fshift,*x,*f;
    _fjsp_v2r8       tx,ty,tz,fscal,rcutoff,rcutoff2,jidxall;
    int              vdwioffset0;
    _fjsp_v2r8       ix0,iy0,iz0,fix0,fiy0,fiz0,iq0,isai0;
    int              vdwjidx0A,vdwjidx0B;
    _fjsp_v2r8       jx0,jy0,jz0,fjx0,fjy0,fjz0,jq0,isaj0;
    _fjsp_v2r8       dx00,dy00,dz00,rsq00,rinv00,rinvsq00,r00,qq00,c6_00,c12_00;
    _fjsp_v2r8       velec,felec,velecsum,facel,crf,krf,krf2;
    real             *charge;
    _fjsp_v2r8       rt,vfeps,vftabscale,Y,F,G,H,Heps,Fp,VV,FF,twovfeps;
    real             *vftab;
    _fjsp_v2r8       itab_tmp;
    _fjsp_v2r8       dummy_mask,cutoff_mask;
    _fjsp_v2r8       one     = gmx_fjsp_set1_v2r8(1.0);
    _fjsp_v2r8       two     = gmx_fjsp_set1_v2r8(2.0);
    union { _fjsp_v2r8 simd; long long int i[2]; } vfconv,gbconv,ewconv;

    x                = xx[0];
    f                = ff[0];

    nri              = nlist->nri;
    iinr             = nlist->iinr;
    jindex           = nlist->jindex;
    jjnr             = nlist->jjnr;
    shiftidx         = nlist->shift;
    gid              = nlist->gid;
    shiftvec         = fr->shift_vec[0];
    fshift           = fr->fshift[0];
    facel            = gmx_fjsp_set1_v2r8(fr->epsfac);
    charge           = mdatoms->chargeA;

    vftab            = kernel_data->table_elec->data;
    vftabscale       = gmx_fjsp_set1_v2r8(kernel_data->table_elec->scale);

    /* Avoid stupid compiler warnings */
    jnrA = jnrB = 0;
    j_coord_offsetA = 0;
    j_coord_offsetB = 0;

    outeriter        = 0;
    inneriter        = 0;

    /* Start outer loop over neighborlists */
    for(iidx=0; iidx<nri; iidx++)
    {
        /* Load shift vector for this list */
        i_shift_offset   = DIM*shiftidx[iidx];

        /* Load limits for loop over neighbors */
        j_index_start    = jindex[iidx];
        j_index_end      = jindex[iidx+1];

        /* Get outer coordinate index */
        inr              = iinr[iidx];
        i_coord_offset   = DIM*inr;

        /* Load i particle coords and add shift vector */
        gmx_fjsp_load_shift_and_1rvec_broadcast_v2r8(shiftvec+i_shift_offset,x+i_coord_offset,&ix0,&iy0,&iz0);

        fix0             = _fjsp_setzero_v2r8();
        fiy0             = _fjsp_setzero_v2r8();
        fiz0             = _fjsp_setzero_v2r8();

        /* Load parameters for i particles */
        iq0              = _fjsp_mul_v2r8(facel,gmx_fjsp_load1_v2r8(charge+inr+0));

        /* Start inner kernel loop */
        for(jidx=j_index_start; jidx<j_index_end-1; jidx+=2)
        {

            /* Get j neighbor index, and coordinate index */
            jnrA             = jjnr[jidx];
            jnrB             = jjnr[jidx+1];
            j_coord_offsetA  = DIM*jnrA;
            j_coord_offsetB  = DIM*jnrB;

            /* load j atom coordinates */
            gmx_fjsp_load_1rvec_2ptr_swizzle_v2r8(x+j_coord_offsetA,x+j_coord_offsetB,
                                              &jx0,&jy0,&jz0);

            /* Calculate displacement vector */
            dx00             = _fjsp_sub_v2r8(ix0,jx0);
            dy00             = _fjsp_sub_v2r8(iy0,jy0);
            dz00             = _fjsp_sub_v2r8(iz0,jz0);

            /* Calculate squared distance and things based on it */
            rsq00            = gmx_fjsp_calc_rsq_v2r8(dx00,dy00,dz00);

            rinv00           = gmx_fjsp_invsqrt_v2r8(rsq00);

            /* Load parameters for j particles */
            jq0              = gmx_fjsp_load_2real_swizzle_v2r8(charge+jnrA+0,charge+jnrB+0);

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            r00              = _fjsp_mul_v2r8(rsq00,rinv00);

            /* Compute parameters for interactions between i and j atoms */
            qq00             = _fjsp_mul_v2r8(iq0,jq0);

            /* Calculate table index by multiplying r with table scale and truncate to integer */
            rt               = _fjsp_mul_v2r8(r00,vftabscale);
            itab_tmp         = _fjsp_dtox_v2r8(rt);
            vfeps            = _fjsp_sub_v2r8(rt, _fjsp_xtod_v2r8(itab_tmp));
            twovfeps         = _fjsp_add_v2r8(vfeps,vfeps);
            _fjsp_store_v2r8(&vfconv.simd,itab_tmp);

            vfconv.i[0]     *= 4;
            vfconv.i[1]     *= 4;

            /* CUBIC SPLINE TABLE ELECTROSTATICS */
            Y                = _fjsp_load_v2r8( vftab + vfconv.i[0] );
            F                = _fjsp_load_v2r8( vftab + vfconv.i[1] );
            GMX_FJSP_TRANSPOSE2_V2R8(Y,F);
            G                = _fjsp_load_v2r8( vftab + vfconv.i[0] +2);
            H                = _fjsp_load_v2r8( vftab + vfconv.i[1] +2);
            GMX_FJSP_TRANSPOSE2_V2R8(G,H);
            Fp               = _fjsp_madd_v2r8(vfeps,_fjsp_madd_v2r8(vfeps,H,G),F);
            FF               = _fjsp_madd_v2r8(_fjsp_madd_v2r8(twovfeps,H,G),vfeps,Fp);
            felec            = _fjsp_neg_v2r8(_fjsp_mul_v2r8(_fjsp_mul_v2r8(qq00,FF),_fjsp_mul_v2r8(vftabscale,rinv00)));

            fscal            = felec;

            /* Update vectorial force */
            fix0             = _fjsp_madd_v2r8(dx00,fscal,fix0);
            fiy0             = _fjsp_madd_v2r8(dy00,fscal,fiy0);
            fiz0             = _fjsp_madd_v2r8(dz00,fscal,fiz0);
            
            gmx_fjsp_decrement_fma_1rvec_2ptr_swizzle_v2r8(f+j_coord_offsetA,f+j_coord_offsetB,fscal,dx00,dy00,dz00);

            /* Inner loop uses 42 flops */
        }

        if(jidx<j_index_end)
        {

            jnrA             = jjnr[jidx];
            j_coord_offsetA  = DIM*jnrA;

            /* load j atom coordinates */
            gmx_fjsp_load_1rvec_1ptr_swizzle_v2r8(x+j_coord_offsetA,
                                              &jx0,&jy0,&jz0);

            /* Calculate displacement vector */
            dx00             = _fjsp_sub_v2r8(ix0,jx0);
            dy00             = _fjsp_sub_v2r8(iy0,jy0);
            dz00             = _fjsp_sub_v2r8(iz0,jz0);

            /* Calculate squared distance and things based on it */
            rsq00            = gmx_fjsp_calc_rsq_v2r8(dx00,dy00,dz00);

            rinv00           = gmx_fjsp_invsqrt_v2r8(rsq00);

            /* Load parameters for j particles */
            jq0              = _fjsp_loadl_v2r8(_fjsp_setzero_v2r8(),charge+jnrA+0);

            /**************************
             * CALCULATE INTERACTIONS *
             **************************/

            r00              = _fjsp_mul_v2r8(rsq00,rinv00);

            /* Compute parameters for interactions between i and j atoms */
            qq00             = _fjsp_mul_v2r8(iq0,jq0);

            /* Calculate table index by multiplying r with table scale and truncate to integer */
            rt               = _fjsp_mul_v2r8(r00,vftabscale);
            itab_tmp         = _fjsp_dtox_v2r8(rt);
            vfeps            = _fjsp_sub_v2r8(rt, _fjsp_xtod_v2r8(itab_tmp));
            twovfeps         = _fjsp_add_v2r8(vfeps,vfeps);
            _fjsp_store_v2r8(&vfconv.simd,itab_tmp);

            vfconv.i[0]     *= 4;
            vfconv.i[1]     *= 4;

            /* CUBIC SPLINE TABLE ELECTROSTATICS */
            Y                = _fjsp_load_v2r8( vftab + vfconv.i[0] );
            F                = _fjsp_setzero_v2r8();
            GMX_FJSP_TRANSPOSE2_V2R8(Y,F);
            G                = _fjsp_load_v2r8( vftab + vfconv.i[0] +2);
            H                = _fjsp_setzero_v2r8();
            GMX_FJSP_TRANSPOSE2_V2R8(G,H);
            Fp               = _fjsp_madd_v2r8(vfeps,_fjsp_madd_v2r8(vfeps,H,G),F);
            FF               = _fjsp_madd_v2r8(_fjsp_madd_v2r8(twovfeps,H,G),vfeps,Fp);
            felec            = _fjsp_neg_v2r8(_fjsp_mul_v2r8(_fjsp_mul_v2r8(qq00,FF),_fjsp_mul_v2r8(vftabscale,rinv00)));

            fscal            = felec;

            fscal            = _fjsp_unpacklo_v2r8(fscal,_fjsp_setzero_v2r8());

            /* Update vectorial force */
            fix0             = _fjsp_madd_v2r8(dx00,fscal,fix0);
            fiy0             = _fjsp_madd_v2r8(dy00,fscal,fiy0);
            fiz0             = _fjsp_madd_v2r8(dz00,fscal,fiz0);
            
            gmx_fjsp_decrement_fma_1rvec_1ptr_swizzle_v2r8(f+j_coord_offsetA,fscal,dx00,dy00,dz00);

            /* Inner loop uses 42 flops */
        }

        /* End of innermost loop */

        gmx_fjsp_update_iforce_1atom_swizzle_v2r8(fix0,fiy0,fiz0,
                                              f+i_coord_offset,fshift+i_shift_offset);

        /* Increment number of inner iterations */
        inneriter                  += j_index_end - j_index_start;

        /* Outer loop uses 7 flops */
    }

    /* Increment number of outer iterations */
    outeriter        += nri;

    /* Update outer/inner flops */

    inc_nrnb(nrnb,eNR_NBKERNEL_ELEC_F,outeriter*7 + inneriter*42);
}
