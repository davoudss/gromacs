/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifndef _pme_h
#define _pme_h

#include <stdio.h>
#include "visibility.h"
#include "typedefs.h"
#include "gmxcomplex.h"
#include "gmx_wallcycle.h"
#include "gmx_parallel_3dfft.h"
#include "../src/mdlib/SE_fgg.h"


#ifdef __cplusplus
extern "C" {
#endif

typedef real *splinevec[DIM];

enum {
    GMX_SUM_QGRID_FORWARD, GMX_SUM_QGRID_BACKWARD
};

GMX_LIBMD_EXPORT
int gmx_pme_init(gmx_pme_t *pmedata, t_commrec *cr,
                 int nnodes_major, int nnodes_minor,
                 t_inputrec *ir, int homenr,
                 gmx_bool bFreeEnergy, gmx_bool bReproducible, int nthread);
/* Initialize the pme data structures resepectively.
 * Return value 0 indicates all well, non zero is an error code.
 */

GMX_LIBMD_EXPORT
int gmx_pme_reinit(gmx_pme_t *         pmedata,
                   t_commrec *         cr,
                   gmx_pme_t           pme_src,
                   const t_inputrec *  ir,
                   ivec                grid_size);
/* As gmx_pme_init, but takes most settings, except the grid, from pme_src */

int gmx_pme_destroy(FILE *log, gmx_pme_t *pmedata);
/* Destroy the pme data structures resepectively.
 * Return value 0 indicates all well, non zero is an error code.
 */

#define GMX_PME_SPREAD_Q      (1<<0)
#define GMX_PME_SOLVE         (1<<1)
#define GMX_PME_CALC_F        (1<<2)
#define GMX_PME_CALC_ENER_VIR (1<<3)
/* This forces the grid to be backtransformed even without GMX_PME_CALC_F */
#define GMX_PME_CALC_POT      (1<<4)
#define GMX_PME_DO_ALL_F  (GMX_PME_SPREAD_Q | GMX_PME_SOLVE | GMX_PME_CALC_F)

int gmx_pme_do(gmx_pme_t pme,
               int start,       int homenr,
               rvec x[],        rvec f[],
               real chargeA[],  real chargeB[],
               matrix box,      t_commrec *cr,
               int  maxshift_x, int maxshift_y,
               t_nrnb *nrnb,    gmx_wallcycle_t wcycle,
               matrix lrvir,    real ewaldcoeff,
               real *energy,    real lambda,
               real *dvdlambda, int flags);
/* Do a PME calculation for the long range electrostatics.
 * flags, defined above, determine which parts of the calculation are performed.
 * Return value 0 indicates all well, non zero is an error code.
 */

GMX_LIBMD_EXPORT
int gmx_pmeonly(gmx_pme_t pme,
                t_commrec *cr,     t_nrnb *mynrnb,
                gmx_wallcycle_t wcycle,
                real ewaldcoeff,   gmx_bool bGatherOnly,
                t_inputrec *ir);
/* Called on the nodes that do PME exclusively (as slaves)
 */

void gmx_pme_calc_energy(gmx_pme_t pme, int n, rvec *x, real *q, real *V);
/* Calculate the PME grid energy V for n charges with a potential
 * in the pme struct determined before with a call to gmx_pme_do
 * with at least GMX_PME_SPREAD_Q and GMX_PME_SOLVE specified.
 * Note that the charges are not spread on the grid in the pme struct.
 * Currently does not work in parallel or with free energy.
 */

/* The following three routines are for PME/PP node splitting in pme_pp.c */

/* Abstract type for PME <-> PP communication */
typedef struct gmx_pme_pp *gmx_pme_pp_t;

gmx_pme_pp_t gmx_pme_pp_init(t_commrec *cr);
/* Initialize the PME-only side of the PME <-> PP communication */

void gmx_pme_send_q(t_commrec *cr,
                    gmx_bool bFreeEnergy, real *chargeA, real *chargeB,
                    int maxshift_x, int maxshift_y);
/* Send the charges and maxshift to out PME-only node. */

void gmx_pme_send_x(t_commrec *cr, matrix box, rvec *x,
                    gmx_bool bFreeEnergy, real lambda,
                    gmx_bool bEnerVir,
                    gmx_large_int_t step);
/* Send the coordinates to our PME-only node and request a PME calculation */

GMX_LIBMD_EXPORT
void gmx_pme_send_finish(t_commrec *cr);
/* Tell our PME-only node to finish */

GMX_LIBMD_EXPORT
void gmx_pme_send_switchgrid(t_commrec *cr, ivec grid_size, real ewaldcoeff);
/* Tell our PME-only node to switch to a new grid size */

GMX_LIBMD_EXPORT
void gmx_pme_send_resetcounters(t_commrec *cr, gmx_large_int_t step);
/* Tell our PME-only node to reset all cycle and flop counters */

void gmx_pme_receive_f(t_commrec *cr,
                       rvec f[], matrix vir,
                       real *energy, real *dvdlambda,
                       float *pme_cycles);
/* PP nodes receive the long range forces from the PME nodes */

/* Return values for gmx_pme_recv_q_x */
enum {
    pmerecvqxX,            /* calculate PME mesh interactions for new x    */
    pmerecvqxFINISH,       /* the simulation should finish, we should quit */
    pmerecvqxSWITCHGRID,   /* change the PME grid size                     */
    pmerecvqxRESETCOUNTERS /* reset the cycle and flop counters            */
};

int gmx_pme_recv_q_x(gmx_pme_pp_t pme_pp,
                     int *natoms,
                     real **chargeA, real **chargeB,
                     matrix box, rvec **x, rvec **f,
                     int *maxshift_x, int *maxshift_y,
                     gmx_bool *bFreeEnergy, real *lambda,
                     gmx_bool *bEnerVir,
                     gmx_large_int_t *step,
                     ivec grid_size, real *ewaldcoeff);
;
/* With return value:
 * pmerecvqxX:             all parameters set, chargeA and chargeB can be NULL
 * pmerecvqxFINISH:        no parameters set
 * pmerecvqxSWITCHGRID:    only grid_size and *ewaldcoeff are set
 * pmerecvqxRESETCOUNTERS: *step is set
 */

void gmx_pme_send_force_vir_ener(gmx_pme_pp_t pme_pp,
                                 rvec *f, matrix vir,
                                 real energy, real dvdlambda,
                                 float cycles);
/* Send the PME mesh force, virial and energy to the PP-only nodes */

#ifdef __cplusplus
}
#endif

#endif





// MOVED FROM PME.C TO SIMPLIFY THE CODE:
typedef struct {
    int send_index0;
    int send_nindex;
    int recv_index0;
    int recv_nindex;
    int recv_size;   // Receive buffer width, used with OpenMP 
} pme_grid_comm_t;


typedef struct {
#ifdef GMX_MPI
    MPI_Comm         mpi_comm;
#endif
    int              nnodes, nodeid;
    int             *s2g0;
    int             *s2g1;
    int              noverlap_nodes;
    int             *send_id, *recv_id;
    int              send_size; /* Send buffer width, used with OpenMP */
    pme_grid_comm_t *comm_data;
    real            *sendbuf;
    real            *recvbuf;
} pme_overlap_t;

typedef struct {
    int *n;      /* Cumulative counts of the number of particles per thread */
    int  nalloc; /* Allocation size of i */
    int *i;      /* Particle indices ordered on thread index (n) */
} thread_plist_t;

typedef struct {
    int      *thread_one;
    int       n;
    int      *ind;
    splinevec theta;
    real     *ptr_theta_z;
    splinevec dtheta;
    real     *ptr_dtheta_z;
    // DDDD
    real     *zs;
    int      *idx;
} splinedata_t;

typedef struct {
    int      dimind;        /* The index of the dimension, 0=x, 1=y */
    int      nslab;
    int      nodeid;
#ifdef GMX_MPI
    MPI_Comm mpi_comm;
#endif

    int     *node_dest;     /* The nodes to send x and q to with DD */
    int     *node_src;      /* The nodes to receive x and q from with DD */
    int     *buf_index;     /* Index for commnode into the buffers */

    int      maxshift;

    int      npd;
    int      pd_nalloc;
    int     *pd;
    int     *count;         /* The number of atoms to send to each node */
    int    **count_thread;
    int     *rcount;        /* The number of atoms to receive */

    int      n;
    int      nalloc;
    rvec    *x;
    real    *q;
    rvec    *f;
    gmx_bool bSpread;       /* These coordinates are used for spreading */
    int      pme_order;
    ivec    *idx;
    rvec    *fractx;            /* Fractional coordinate relative to the
                                 * lower cell boundary
                                 */
    int             nthread;
    int            *thread_idx; /* Which thread should spread which charge */
    thread_plist_t *thread_plist;
    splinedata_t   *spline;
} pme_atomcomm_t;

#define FLBS  3
#define FLBSZ 4

typedef struct {
    ivec  ci;     /* The spatial location of this grid         */
    ivec  n;      /* The used size of *grid, including order-1 */
    ivec  offset; /* The grid offset from the full node grid   */
    int   order;  /* PME spreading order                       */
    ivec  s;      /* The allocated size of *grid, s >= n       */
    real *grid;   /* The grid local thread, size n             */
} pmegrid_t;

typedef struct {
  pmegrid_t  grid;         /* The full node grid (non thread-local)            */
  int        nthread;      /* The number of threads operating on this grid     */
  ivec       nc;           /* The local spatial decomposition over the threads */
  pmegrid_t *grid_th;      /* Array of grids for each thread                   */
  real      *grid_all;     /* Allocated array for the grids in *grid_th        */
  int      **g2t;          /* The grid to thread index                         */
  ivec       nthread_comm; /* The number of threads to communicate with        */
} pmegrids_t;


typedef struct {
#ifdef PME_SSE
    /* Masks for SSE aligned spreading and gathering */
    __m128 mask_SSE0[6], mask_SSE1[6];
#else
    int    dummy; /* C89 requires that struct has at least one member */
#endif

} pme_spline_work_t;

typedef struct {
    /* work data for solve_pme */
    int      nalloc;
    real *   mhx;
    real *   mhy;
    real *   mhz;
    real *   m2;
    real *   denom;
    real *   tmp1_alloc;
    real *   tmp1;
    real *   eterm;
    real *   m2inv;

    real     energy;
    matrix   vir;
} pme_work_t;

typedef struct gmx_pme {
    int           ndecompdim; /* The number of decomposition dimensions */
    int           nodeid;     /* Our nodeid in mpi->mpi_comm */
    int           nodeid_major;
    int           nodeid_minor;
    int           nnodes;    /* The number of nodes doing PME */
    int           nnodes_major;
    int           nnodes_minor;

    MPI_Comm      mpi_comm;
    MPI_Comm      mpi_comm_d[2]; /* Indexed on dimension, 0=x, 1=y */
#ifdef GMX_MPI
    MPI_Datatype  rvec_mpi;      /* the pme vector's MPI type */
#endif

    int        nthread;       /* The number of threads doing PME */

    gmx_bool   bPPnode;       /* Node also does particle-particle forces */
    gmx_bool   bFEP;          /* Compute Free energy contribution */
    int        nkx, nky, nkz; /* Grid dimensions */
    gmx_bool   bP3M;          /* Do P3M: optimize the influence function */
    int        pme_order;
    real       epsilon_r;

    pmegrids_t pmegridA;  /* Grids on which we do spreading/interpolation, includes overlap */
    pmegrids_t pmegridB;
    /* The PME charge spreading grid sizes/strides, includes pme_order-1 */
    int        pmegrid_nx, pmegrid_ny, pmegrid_nz;
    /* pmegrid_nz might be larger than strictly necessary to ensure
     * memory alignment, pmegrid_nz_base gives the real base size.
     */
    int     pmegrid_nz_base;
    /* The local PME grid starting indices */
    int     pmegrid_start_ix, pmegrid_start_iy, pmegrid_start_iz;

    /* Work data for spreading and gathering */
    pme_spline_work_t    *spline_work;

    real                 *fftgridA; /* Grids for FFT. With 1D FFT decomposition this can be a pointer */
    real                 *fftgridB; /* inside the interpolation grid, but separate for 2D PME decomp. */
    int                   fftgrid_nx, fftgrid_ny, fftgrid_nz;

    t_complex            *cfftgridA;  /* Grids for complex FFT data */
    t_complex            *cfftgridB;
    int                   cfftgrid_nx, cfftgrid_ny, cfftgrid_nz;

    gmx_parallel_3dfft_t  pfft_setupA;
    gmx_parallel_3dfft_t  pfft_setupB;

    int                  *nnx, *nny, *nnz;
    real                 *fshx, *fshy, *fshz;

    pme_atomcomm_t        atc[2]; /* Indexed on decomposition index */
    matrix                recipbox;
    splinevec             bsp_mod;

    pme_overlap_t         overlap[2]; /* Indexed on dimension, 0=x, 1=y */

    pme_atomcomm_t        atc_energy; /* Only for gmx_pme_calc_energy */

    rvec                 *bufv;       /* Communication buffer */
    real                 *bufr;       /* Communication buffer */
    int                   buf_nalloc; /* The communication buffer size */

    /* thread local work data for solve_pme */
    pme_work_t *work;

    /* Work data for PME_redist */
    gmx_bool redist_init;
    int *    scounts;
    int *    rcounts;
    int *    sdispls;
    int *    rdispls;
    int *    sidx;
    int *    idxa;
    real *   redist_buf;
    int      redist_buf_nalloc;

    /* Work data for sum_qgrid */
    real *   sum_qgrid_tmp;
    real *   sum_qgrid_dd_tmp;

} t_gmx_pme;
