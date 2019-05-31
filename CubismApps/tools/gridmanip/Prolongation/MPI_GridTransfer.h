/* File:   MPI_GridTransfer.h */
/* Date:   August 2015 */
/* Author: Ursula Rasthofer */
/* Tag:    Transfer solution fields from coarse to fine grid */
/* Copyright 2015 ETH Zurich. All Rights Reserved. */
#ifndef MPI_GRIDTRANSFER_H_JCMHPZQU
#define MPI_GRIDTRANSFER_H_JCMHPZQU

#include <iostream>
#include <vector>

#include <Cubism/BlockInfo.h>
#include "Types.h"
#include "BlockProcessor_MPI.h"
using namespace cubism;

struct grid_smoother
{
    StencilInfo stencil;
    int stencil_start[3];
    int stencil_end[3];

    grid_smoother():
        stencil(-1,-1,-1,2,2,2, false, 8, 0,1,2,3,4,5,6,7)
    {
        stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
    }

    grid_smoother(const grid_smoother& c):
        stencil(-1,-1,-1,2,2,2, false, 8, 0,1,2,3,4,5,6,7)
    {
        stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;
    }

    template<typename TLab, typename TBlock>
    inline void operator()(TLab& lab, const BlockInfo& info, TBlock& o) const
    {
        typedef typename TBlock::ElementType TElement;
        for(int iz=0; iz<TBlock::sizeZ; iz++)
            for(int iy=0; iy<TBlock::sizeY; iy++)
                for(int ix=0; ix<TBlock::sizeX; ix++)
                {
                    // myself
                    TElement sum = lab(ix,iy,iz);
                    // neighbor faces
                    sum = sum + 0.5*(
                            lab(ix-1,iy,iz) + lab(ix+1,iy,iz) +
                            lab(ix,iy-1,iz) + lab(ix,iy+1,iz) +
                            lab(ix,iy,iz-1) + lab(ix,iy,iz+1)
                            );

                    // // tensorial
                    // // neighbor edges
                    // sum = sum + 0.25*(
                    //         lab(ix,iy-1,iz-1) + lab(ix,iy-1,iz+1) + lab(ix-1,iy-1,iz) + lab(ix+1,iy-1,iz) +
                    //         lab(ix,iy+1,iz-1) + lab(ix,iy+1,iz+1) + lab(ix-1,iy+1,iz) + lab(ix+1,iy+1,iz) +
                    //         lab(ix-1,iy,iz-1) + lab(ix-1,iy,iz+1) + lab(ix+1,iy,iz-1) + lab(ix+1,iy,iz+1)
                    //         );
                    // // neighbor corners
                    // sum = sum + 0.125*(
                    //         lab(ix-1,iy-1,iz-1) + lab(ix-1,iy-1,iz+1) + lab(ix-1,iy+1,iz-1) + lab(ix-1,iy+1,iz+1) +
                    //         lab(ix+1,iy-1,iz-1) + lab(ix+1,iy-1,iz+1) + lab(ix+1,iy+1,iz-1) + lab(ix+1,iy+1,iz+1)
                    //         );
                    // o(ix,iy,iz) = 0.125*sum;

                    o(ix,iy,iz) = 0.25*sum;
                }
    }
};

template <typename TGrid>
struct grid_transfer
{
    // define variables

    // number of neighboring (coarse) cells that are
    // required to compute the current values at a fine cell
    // (corresponds to stencil when discrtizing pdes)
    StencilInfo stencil;
    int stencil_start[3], stencil_end[3];

    // point to the fine mesh which should be filled
    typedef typename TGrid::BlockType TBlock;
    TGrid& fine_grid;

    // number of coarse blocks in each spatial dimension
    const int cNX;
    const int cNY;
    const int cNZ;

    // number of fine blocks in each spatial dimension
    const int fNX;
    const int fNY;
    const int fNZ;

    // contructor
    grid_transfer(TGrid& fgrid,
            const int cnx, const int cny, const int cnz,
            const int fnx, const int fny, const int fnz):
        stencil(-2,-2,-2,3,3,3, true, 8, 0,1,2,3,4,5,6,7),
        fine_grid(fgrid),
        cNX(cnx),cNY(cny),cNZ(cnz),
        fNX(fnx),fNY(fny),fNZ(fnz)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -2;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 3;
    }

    // copy constructor
    grid_transfer(const grid_transfer& gt):
        stencil(-2,-2,-2,3,3,3, true, 8, 0,1,2,3,4,5,6,7),
        fine_grid(gt.fine_grid),
        cNX(gt.cNX),cNY(gt.cNY),cNZ(gt.cNZ),
        fNX(gt.fNX),fNY(gt.fNY),fNZ(gt.fNZ)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -2;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 3;
    }

    // function to tranfer x,y,z-index space into block id for fine grid
    int calc_ID_fine(const int ix, const int iy, const int iz)
    {
        return ix + fNX*(iy + fNY*iz);
    }

    // function to calculate offset for index space of cell for coarse grid
    std::vector<int> calc_idx_offset(const int ix, const int iy, const int iz)
    {
        std::vector<int> idx_off(3,0);
        idx_off[0] = (ix%2) * TBlock::sizeX / 2;
        idx_off[1] = (iy%2) * TBlock::sizeY / 2;
        idx_off[2] = (iz%2) * TBlock::sizeZ / 2;
        return idx_off;
    }

    // references:
    // - A. Harten, Multiresolution algorithms for the numerical solution of hyperbolic conservation laws,
    //   Comm. Pur Appl. Math. 48 (1995) 1305-1342.
    // - B. L. Bihari & A. Harten, Multiresolution schemes for the numerical solution of of 2-d conservation laws I.
    //   SIAM J. Sci. Comput. 18 (1997) 315-354.
    // - O. Roussel et al., A conservative fully adaptive multiresolution algorithm for parabolic PDEs,
    //   Journal Comput. Phys. 188 (2003) 493-523.
    // Be careful with the signs given in those references!.
    template<typename LabType, typename ReturnType>
    ReturnType _interpolate_HARTEN(const int f_cell_ix,     // first/x index of fine cell
            const int f_cell_iy,     // second/y index of fine cell
            const int f_cell_iz,     // third/z index of fine cell
            const int c_cell_ix_off, // first/x index offset of coarse cell
            const int c_cell_iy_off, // second/y index offset of coarse cell
            const int c_cell_iz_off, // third/z index offset of coarse cell
            LabType &lab,            // the lab corresponding to the coarse cell
            const int s,             // stencil width related to the order r as r=2s+1
            const Real gamma[2])     // weights
    {
        // to be filled
        ReturnType out;

        // index to distinguish first and second fine cell in coarse cell in one spatial direction
        const int n = f_cell_ix%2;
        const int p = f_cell_iy%2;
        const int q = f_cell_iz%2;

        // get coarse cell indices containing the current fine cell
        const int c_cell_ix = c_cell_ix_off + (f_cell_ix - n)/2;
        const int c_cell_iy = c_cell_iy_off + (f_cell_iy - p)/2;
        const int c_cell_iz = c_cell_iz_off + (f_cell_iz - q)/2;

        // compute all terms
        ReturnType Qsx;
        Qsx.clear();
        ReturnType Qsy;
        Qsy.clear();
        ReturnType Qsz;
        Qsz.clear();
        ReturnType Qsxy;
        Qsxy.clear();
        ReturnType Qsxz;
        Qsxz.clear();
        ReturnType Qsyz;
        Qsyz.clear();
        ReturnType Qsxyz;
        Qsxyz.clear();

        for (int rr=1; rr<=s; rr++)
        {
            Qsx = Qsx + gamma[rr-1] * (lab(c_cell_ix+rr, c_cell_iy, c_cell_iz) - lab(c_cell_ix-rr, c_cell_iy, c_cell_iz));
            Qsy = Qsy + gamma[rr-1] * (lab(c_cell_ix, c_cell_iy+rr, c_cell_iz) - lab(c_cell_ix, c_cell_iy-rr, c_cell_iz));
            Qsz = Qsz + gamma[rr-1] * (lab(c_cell_ix, c_cell_iy, c_cell_iz+rr) - lab(c_cell_ix, c_cell_iy, c_cell_iz-rr));

            for (int ss=1; ss<=s; ss++)
            {
                Qsxy = Qsxy + gamma[rr-1] * gamma[ss-1] * (lab(c_cell_ix+rr, c_cell_iy+ss, c_cell_iz) - lab(c_cell_ix+rr, c_cell_iy-ss, c_cell_iz)
                        - lab(c_cell_ix-rr, c_cell_iy+ss, c_cell_iz) + lab(c_cell_ix-rr, c_cell_iy-ss, c_cell_iz));
                Qsxz = Qsxz + gamma[rr-1] * gamma[ss-1] * (lab(c_cell_ix+rr, c_cell_iy, c_cell_iz+ss) - lab(c_cell_ix+rr, c_cell_iy, c_cell_iz-ss)
                        - lab(c_cell_ix-rr, c_cell_iy, c_cell_iz+ss) + lab(c_cell_ix-rr, c_cell_iy, c_cell_iz-ss));
                Qsyz = Qsyz + gamma[rr-1] * gamma[ss-1] * (lab(c_cell_ix, c_cell_iy+rr, c_cell_iz+ss) - lab(c_cell_ix, c_cell_iy+rr, c_cell_iz-ss)
                        - lab(c_cell_ix, c_cell_iy-rr, c_cell_iz+ss) + lab(c_cell_ix, c_cell_iy-rr, c_cell_iz-ss));

                for (int tt=1; tt<=s; tt++)
                    Qsxyz = Qsxyz + gamma[rr-1] * gamma[ss-1] * gamma[tt-1] * (lab(c_cell_ix+rr, c_cell_iy+ss, c_cell_iz+tt) - lab(c_cell_ix+rr, c_cell_iy+ss, c_cell_iz-tt)
                            - lab(c_cell_ix+rr, c_cell_iy-ss, c_cell_iz+tt) - lab(c_cell_ix-rr, c_cell_iy+ss, c_cell_iz+tt)
                            + lab(c_cell_ix+rr, c_cell_iy-ss, c_cell_iz-tt) + lab(c_cell_ix-rr, c_cell_iy+ss, c_cell_iz-tt)
                            + lab(c_cell_ix-rr, c_cell_iy-ss, c_cell_iz+tt) - lab(c_cell_ix-rr, c_cell_iy-ss, c_cell_iz-tt));
            }
        }

        // sum up
        out = lab(c_cell_ix, c_cell_iy, c_cell_iz)
        + static_cast<Real>(pow(-1.0,(double) n)) * Qsx + static_cast<Real>(pow(-1.0,(double) p)) * Qsy + static_cast<Real>(pow(-1.0,(double) q)) * Qsz
        + static_cast<Real>(pow(-1.0,(double) (n+p))) * Qsxy + static_cast<Real>(pow(-1.0,(double) (n+q))) * Qsxz + static_cast<Real>(pow(-1.0,(double) (p+q))) * Qsyz
        + static_cast<Real>(pow(-1.0,(double) (n+p+q))) * Qsxyz;

        return out;
    }


    template<typename LabType, typename BlockType>
    void operator()(LabType& lab, const BlockInfo& info, BlockType& o)
    {
        // get id of current coarse grid block
        const long long cID = info.blockID;

        // get coarse block in x,y,z-index space
        const int cbiz = (int) floor(((double) cID) / ((double) (cNX * cNY)));
        const int tmp = cID % (cNX * cNY);
        const int cbiy = (int) floor(((double) tmp) / ((double) (cNX)));
        const int cbix = tmp % cNX;

        // get first fine block in x,y,z-index space
        const int fbix = 2 * cbix;
        const int fbiy = 2 * cbiy;
        const int fbiz = 2 * cbiz;
        // get vector of all fine blocks in current coarse block
        std::vector< std::vector<int> > fidx(8);
        std::vector<int> vtmp(3);
        vtmp[0] = fbix; vtmp[1] = fbiy; vtmp[2] = fbiz;
        fidx[0] = vtmp;
        vtmp[0] = fbix+1; vtmp[1] = fbiy; vtmp[2] = fbiz;
        fidx[1] = vtmp;
        vtmp[0] = fbix; vtmp[1] = fbiy+1; vtmp[2] = fbiz;
        fidx[2] = vtmp;
        vtmp[0] = fbix; vtmp[1] = fbiy; vtmp[2] = fbiz+1;
        fidx[3] = vtmp;
        vtmp[0] = fbix+1; vtmp[1] = fbiy+1; vtmp[2] = fbiz;
        fidx[4] = vtmp;
        vtmp[0] = fbix+1; vtmp[1] = fbiy; vtmp[2] = fbiz+1;
        fidx[5] = vtmp;
        vtmp[0] = fbix; vtmp[1] = fbiy+1; vtmp[2] = fbiz+1;
        fidx[6] = vtmp;
        vtmp[0] = fbix+1; vtmp[1] = fbiy+1; vtmp[2] = fbiz+1;
        fidx[7] = vtmp;

        // get all blocks of fine grid
        std::vector<BlockInfo> f_info = fine_grid.getBlocksInfo();

        typedef typename BlockType::ElementType TElement;

        // loop all blocks of fine grid within the current block of the coarse grid
        for(int i=0; i<8; i++)
        {
            // get fine-grid block via its id calculated from indices
            const int fID =  calc_ID_fine(fidx[i][0], fidx[i][1], fidx[i][2]);
            BlockInfo cf_info = f_info[fID];
            BlockType& b = *(BlockType*)cf_info.ptrBlock;

            // calculate starting index for coarse cells comprising current fine grid block
            std::vector<int> c_idx_off = calc_idx_offset(fidx[i][0], fidx[i][1], fidx[i][2]);

            // loop all cells of this fine grid block
            for (int iz=0; iz<BlockType::sizeZ; iz++)
                for (int iy=0; iy<BlockType::sizeY; iy++)
                    for (int ix=0; ix<BlockType::sizeX; ix++)
                    {
                        // 3rd order accurate interpolation
#ifndef _HARTEN5_
                        const Real gamma[2] = {-0.125,0.0};
                        const TElement actcell = _interpolate_HARTEN<LabType, TElement>(ix, iy, iz, c_idx_off[0], c_idx_off[1], c_idx_off[2], lab, 1, gamma);
                        // 5th order accurate interpolation
#else
                        const Real gamma1 = -22.0/128.0;
                        const Real gamma2 = 3.0/128.0;
                        const Real gamma[2] = {gamma1, gamma2};
                        const TElement actcell = _interpolate_HARTEN<LabType, TElement>(ix, iy, iz, c_idx_off[0], c_idx_off[1], c_idx_off[2], lab, 2, gamma);
#endif
                        b(ix,iy,iz) = actcell;
                    }
        }

        return;
    }

};


template<typename TGridIn, typename TGridOut, typename LabMPI>
void interpolate_from_coarse(const TGridIn& coarse_grid, TGridOut& fine_grid,
        const int in_bpdx, const int in_bpdy, const int in_bpdz,
        const int out_bpdx, const int out_bpdy, const int out_bpdz,
        const int smooth_iter=0, const bool verbose=true)
{
    // construct (previous) coarser grid
    // we assume that the resolution increases by a factor of 2
    assert(out_bpdx%2 == 0 && out_bpdz%2 == 0 && out_bpdz%2 == 0);
    assert(2*in_bpdx == out_bpdx);
    assert(2*in_bpdy == out_bpdy);
    assert(2*in_bpdz == out_bpdz);

    // do the interpolation from the coarse to the fine mesh
    // we fill the fine grid on each rank
    // therefore, we have to set up a halo/ghost cell layer for the
    // coarse grid (hence, the coarse grid goes into the lab)
    grid_transfer<TGridOut> coarse_to_fine(fine_grid,
            in_bpdx, in_bpdy, in_bpdz,
            out_bpdx, out_bpdy, out_bpdz);

    // This requires the that coarse block on THIS process contains the fine
    // blocks
    //
    // o------x------o
    // |      |      |
    // |      |      |
    // x------x------x
    // |      |      |
    // |      |      |
    // o------x------o
    // Coarse block corners "o", fine blocks corners "x", both must be in the
    // memory space of this process
    process<LabMPI>(coarse_to_fine, const_cast<TGridIn&>(coarse_grid), 0, 0);

    // smooth grid
    for (int i = 0; i < smooth_iter; ++i)
    {
        if (verbose)
            std::cout << "smoothing grid: iteration " << i+1 << std::endl;
        grid_smoother smoother;
        process<LabMPI>(smoother, fine_grid, 0, 0);
    }
};

#endif /* MPI_GRIDTRANSFER_H_JCMHPZQU */
