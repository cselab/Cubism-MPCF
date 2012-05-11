/*
 *  Test_CVT.h
 *  MPCFnode
 *
 *  Created by Babak Hejazialhosseini  on 2/10/12.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include "Test_SteadyState.h"
#include <HDF5Dumper.h>

struct EvaluateVorticity_CPP
{   
    StencilInfo stencil;
    int stencil_start[3], stencil_end[3];
    
    EvaluateVorticity_CPP(): stencil(-1,-1,-1,2,2,2, false, 4, 0,1,2,3)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -1;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;   
    }
    
    EvaluateVorticity_CPP(const EvaluateVorticity_CPP& c): stencil(-1,-1,-1,2,2,2, false, 4, 0,1,2,3)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -1;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;    
    }
    
    template<typename LabType, typename BlockType>
    void operator()(LabType& lab, const BlockInfo& info, BlockType& o) const
    {
        typedef BlockType B;
        
        const double h = info.h_gridpoint;
        
        for(int iz=0; iz<BlockType::sizeZ; iz++)
            for(int iy=0; iy<BlockType::sizeY; iy++)
                for(int ix=0; ix<BlockType::sizeX; ix++)
                {
                    const Real uy = lab(ix,iy+1,iz).u/lab(ix,iy+1,iz).rho-lab(ix,iy-1,iz).u/lab(ix,iy-1,iz).rho;
                    const Real vx = lab(ix+1,iy,iz).v/lab(ix+1,iy,iz).rho-lab(ix-1,iy,iz).v/lab(ix-1,iy,iz).rho;
                    const Real uz = lab(ix,iy,iz+1).u/lab(ix,iy,iz+1).rho-lab(ix,iy,iz-1).u/lab(ix,iy,iz-1).rho;
                    const Real wx = lab(ix+1,iy,iz).w/lab(ix+1,iy,iz).rho-lab(ix-1,iy,iz).w/lab(ix-1,iy,iz).rho;
                    const Real vz = lab(ix,iy,iz+1).v/lab(ix,iy,iz+1).rho-lab(ix,iy,iz-1).v/lab(ix,iy,iz-1).rho;
                    const Real wy = lab(ix,iy+1,iz).w/lab(ix,iy+1,iz).rho-lab(ix,iy-1,iz).w/lab(ix,iy-1,iz).rho;
                    
                    o.tmp[iz][iy][ix][0] = 0.5*(wy-vz)/h;
                    o.tmp[iz][iy][ix][1] = 0.5*(uz-wx)/h;
                    o.tmp[iz][iy][ix][2] = 0.5*(vx-uy)/h;
                }
    }
};

template <typename TGrid, typename Kvort, typename TLab>
struct EvaluateVorticity
{
    TGrid & grid;
    
    EvaluateVorticity(TGrid & _grid): grid(_grid) {}
    
    void execute()
    {        
        Kvort eval_vort;
		vector<BlockInfo> vInfo = grid.getBlocksInfo();
        BlockProcessing::process< TLab >(vInfo, eval_vort, (TGrid&)grid, 0);
    }
};


template<typename TGrid, typename TLab>
void DumpHDF5_vorticity(TGrid &grid, const int iCounter, const string f_name, const string dump_path)
{
#ifdef _USE_HDF_
	EvaluateVorticity< TGrid, EvaluateVorticity_CPP, TLab > eval_vort(grid);
	
	eval_vort.execute();
	
	typedef typename TGrid::BlockType B;
	
	char filename[256];
	herr_t status;
	hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;
	
	static const int NCHANNELS = 3;
	const int NX = grid.getBlocksPerDimension(0)*B::sizeX;
	const int NY = grid.getBlocksPerDimension(1)*B::sizeY;
	const int NZ = grid.getBlocksPerDimension(2)*B::sizeZ;
	
	Real * array_all = new Real[NX * NY * NZ * NCHANNELS];
	
	vector<BlockInfo> vInfo_local = grid.getBlocksInfo();
	
	const int sX = 0;
	const int sY = 0;
	const int sZ = 0;
	
	const int eX = B::sizeX;
	const int eY = B::sizeY;
	const int eZ = B::sizeZ;
	
	hsize_t count[4] = {
		grid.getBlocksPerDimension(0)*B::sizeX,
		grid.getBlocksPerDimension(1)*B::sizeY,
		grid.getBlocksPerDimension(2)*B::sizeZ, NCHANNELS};
	
	hsize_t dims[4] = {
		grid.getBlocksPerDimension(0)*B::sizeX,
		grid.getBlocksPerDimension(1)*B::sizeY,
		grid.getBlocksPerDimension(2)*B::sizeZ, NCHANNELS};
	
	hsize_t offset[4] = {0, 0, 0, 0};
	
	sprintf(filename, "%s/%s.h5", dump_path.c_str(), f_name.c_str());
	
	H5open();
	fapl_id = H5Pcreate(H5P_FILE_ACCESS);
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
	status = H5Pclose(fapl_id);
	
#pragma omp parallel for   
	for(int i=0; i<vInfo_local.size(); i++)
	{
		BlockInfo& info = vInfo_local[i];
		const int idx[3] = {info.index[0], info.index[1], info.index[2]};
		B & b = *(B*)info.ptrBlock;
		
		for(int iz=sZ; iz<eZ; iz++)
			for(int iy=sY; iy<eY; iy++)
				for(int ix=sX; ix<eX; ix++)
				{     
					const int gx = idx[0]*B::sizeX + ix;
					const int gy = idx[1]*B::sizeY + iy;
					const int gz = idx[2]*B::sizeZ + iz;
					
					Real * const ptr = array_all + NCHANNELS*(gz + NZ * (gy + NY * gx));
					
					ptr[0] = b.tmp[iz][iy][ix][0];
					ptr[1] = b.tmp[iz][iy][ix][1];
					ptr[2] = b.tmp[iz][iy][ix][2];
				}
	}
	
	fapl_id = H5Pcreate(H5P_DATASET_XFER);
	
	fspace_id = H5Screate_simple(4, dims, NULL);
	dataset_id = H5Dcreate(file_id, "data", HDF_REAL, fspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	fspace_id = H5Dget_space(dataset_id);
	H5Sselect_hyperslab(fspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
	mspace_id = H5Screate_simple(4, count, NULL);        
	status = H5Dwrite(dataset_id, HDF_REAL, mspace_id, fspace_id, fapl_id, array_all);
	
	status = H5Sclose(mspace_id);
	status = H5Sclose(fspace_id);
	status = H5Dclose(dataset_id);
	status = H5Pclose(fapl_id);
	status = H5Fclose(file_id);
	H5close();
	
	{
		char wrapper[256];
		sprintf(wrapper, "%s/%s.xmf", dump_path.c_str(), f_name.c_str());
		FILE *xmf = 0;
		xmf = fopen(wrapper, "w");
		fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
		fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
		fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
		fprintf(xmf, " <Domain>\n");
		fprintf(xmf, "   <Grid GridType=\"Uniform\">\n");
		fprintf(xmf, "     <Time Value=\"%05d\"/>\n", iCounter);
		fprintf(xmf, "     <Topology TopologyType=\"3DCORECTMesh\" Dimensions=\"%d %d %d\"/>\n", (int)dims[0], (int)dims[1], (int)dims[2]);
		fprintf(xmf, "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
		fprintf(xmf, "       <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
		fprintf(xmf, "        %e %e %e\n", 0.,0.,0.);
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "       <DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
		fprintf(xmf, "        %e %e %e\n", 1./(Real)dims[0],1./(Real)dims[0],1./(Real)dims[0]);
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "     </Geometry>\n");
		
		fprintf(xmf, "     <Attribute Name=\"data\" AttributeType=\"Vector\" Center=\"Node\">\n");
		fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (int)dims[0], (int)dims[1], (int)dims[2], (int)dims[3]);
		fprintf(xmf, "        %s:/data\n",(f_name+".h5").c_str());
		fprintf(xmf, "       </DataItem>\n");
		fprintf(xmf, "     </Attribute>\n");
		
		fprintf(xmf, "   </Grid>\n");
		fprintf(xmf, " </Domain>\n");
		fprintf(xmf, "</Xdmf>\n");
		fclose(xmf);
	}
#else
#warning USE OF HDF WAS DISABLED AT COMPILE TIME
#endif
}

class Test_CVT: public Test_SteadyState
{
    friend class Test_CVTMPI;
    
	Real radius, bubble_pos[3], nu1, re;
    Real gamma_max;
    
	void _ic(FluidGrid& grid);
	void _setup_constants();
    void _dumpStatistics(FluidGrid& grid, const int counter, const Real t, const Real dt);
	void _restart();
	

public:	
	Test_CVT(const int argc, const char ** argv): Test_SteadyState(argc, argv) { }
    
	void run();
	void setup();
};

template<typename BlockType, template<typename X> class allocator=std::allocator>
class BlockLab_CVT: public BlockLab<BlockType,allocator>
{		
	typedef typename BlockType::ElementType ElementTypeBlock;

protected:
	bool is_xperiodic() {return false;}
	bool is_yperiodic() {return true;}
	bool is_zperiodic() {return false;}

public:
	BlockLab_CVT(): BlockLab<BlockType,allocator>(){}
	
	void _apply_bc(const BlockInfo& info, const Real t=0)
	{
	  //BlockLab<BlockType, allocator>::load(info);
		
		BoundaryCondition<BlockType,ElementTypeBlock,allocator> bc(this->m_stencilStart, this->m_stencilEnd, this->m_cacheBlock);
		
		if (info.index[0]==0)               bc.template applyBC_reflecting<0,0>();
		if (info.index[0]==this->NX-1)  bc.template applyBC_reflecting<0,1>();
		if (info.index[2]==0)				bc.template applyBC_reflecting<2,0>();
		if (info.index[2]==this->NZ-1)	bc.template applyBC_reflecting<2,1>();
	}
};



