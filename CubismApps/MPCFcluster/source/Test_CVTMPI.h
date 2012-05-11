/*
 *  Test_CVTMPI.h
 *  MPCFcluster
 *
 *  Created by Babak Hejazialhosseini on 2/10/12.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */
#pragma once

#include <limits>

#include <Test_CVT.h>
#include "PoissonSolver_MPI.h"

template<typename TGrid>
void DumpVP_VorticityMagnitudeMPI(TGrid& grid, string fileName, const string dump_path=".");

#ifdef _USE_HDF_
namespace MPCFcluster
{
	template<typename TGrid>
	void DumpHDF5_vorticity(TGrid &grid, const bool bIC, const int iCounter, const string f_name, const string dump_path=".");
}
#endif

struct Evaluate_VelocityDivergence_CPP
{   
    StencilInfo stencil;
    int stencil_start[3], stencil_end[3];
    
    Evaluate_VelocityDivergence_CPP(): stencil(-1,-1,-1,2,2,2, false, 4, 0,1,2,3)
    {
        stencil_start[0] = stencil_start[1] =  stencil_start[2] = -1;
        stencil_end[0] = stencil_end[1] = stencil_end[2] = 2;   
    }
    
    Evaluate_VelocityDivergence_CPP(const Evaluate_VelocityDivergence_CPP& c): stencil(-1,-1,-1,2,2,2, false, 4, 0,1,2,3)
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
                    const Real ux= lab(ix+1,iy,iz).u/lab(ix+1,iy,iz).rho-lab(ix-1,iy,iz).u/lab(ix-1,iy,iz).rho;
                    const Real vy= lab(ix,iy+1,iz).v/lab(ix,iy+1,iz).rho-lab(ix,iy-1,iz).v/lab(ix,iy-1,iz).rho;
                    const Real wz= lab(ix,iy,iz+1).w/lab(ix,iy,iz+1).rho-lab(ix,iy,iz-1).w/lab(ix,iy,iz-1).rho;
                    
                    o.tmp[iz][iy][ix][3] = 0.5*(ux+vy+wz)/h;
                }
    }
};

template <typename TGrid, typename Kvort, typename TLab>
struct EvaluateVorticityMPI
{
    TGrid & grid;
    
    EvaluateVorticityMPI(TGrid & _grid): grid(_grid) {}
    
    void execute()
    {        
        Kvort eval_vort;
        SynchronizerMPI& synch = ((TGrid&)grid).sync(eval_vort);
        while (!synch.done())
        {                
            vector<BlockInfo> avail = synch.avail(2);
            BlockProcessingMPI::process< TLab >(avail, eval_vort, (TGrid&)grid);
        }
    }
};

template <typename TGrid, typename Kveldiv, typename TLab>
struct EvaluateVelocityDivergence
{
    TGrid & grid;
    
    EvaluateVelocityDivergence(TGrid & _grid): grid(_grid) {}
    
    void execute()
    {        
        Kveldiv eval_veldiv;
        SynchronizerMPI& synch = ((TGrid&)grid).sync(eval_veldiv);
        while (!synch.done())
        {                
            vector<BlockInfo> avail = synch.avail(2);
            BlockProcessingMPI::process< TLab >(avail, eval_veldiv, (TGrid&)grid);
        }
    }
};

class Test_CVTMPI: public Test_CVT
{
    Test_SteadyStateMPI * t_ssmpi;
    Real gamma_max;
    
protected:
    int XPESIZE, YPESIZE, ZPESIZE;
	
    G * grid;
	FlowStep_LSRK3MPI<G> * stepper;
	
	void _dump_vorticity(bool bIC)
	{
		EvaluateVorticityMPI< G, EvaluateVorticity_CPP, LabMPI > eval_vort(*grid);
		
		if (bIC)
			EvaluateVorticityMPI< G, EvaluateVorticity_CPP, LabMPI_PossionSolver > eval_vort(*grid);
		
		eval_vort.execute();
		
		std::stringstream streamer_vort;
        streamer_vort<<"vorticity-"<<step_id;
		
		const string path = parser("-fpath").asString(".");

		DumpHDF5_MPI<G, StreamerFromTemp_HDF5>(*grid, step_id, streamer_vort.str(), path);		
	}
    
public:
    bool isroot;
    
	Test_CVTMPI(const bool isroot, const int argc, const char ** argv):
    Test_CVT(argc, argv), isroot(isroot)
    {
		t_ssmpi = new Test_SteadyStateMPI(isroot, argc, argv);
    }
    
	void setup()
	{
		if (isroot && VERBOSITY)
		{
			printf("////////////////////////////////////////////////////////////\n");
			printf("////     TEST Compressible Vortex Reconnection MPI   ///////\n");
			printf("////////////////////////////////////////////////////////////\n");
		}
		
        _setup_constants();
        t_ssmpi->setup_mpi_constants(XPESIZE, YPESIZE, ZPESIZE);
        
		if (!isroot)
			VERBOSITY = 0;
		
		grid = new G(XPESIZE, YPESIZE, ZPESIZE, BPDX, BPDY, BPDZ);
		
		assert(grid != NULL);
		
        stepper = new FlowStep_LSRK3MPI<G>(*grid, CFL, Simulation_Environment::GAMMA1, Simulation_Environment::GAMMA2, parser, VERBOSITY, &profiler);
                
		if(bRESTART)
		{
            t_ssmpi->restart(*grid);
			const string path = parser("-fpath").asString(".");
			std::stringstream streamer_vort;
			streamer_vort<<"vorticity_restarted-"<<step_id;
#ifdef _USE_HDF_
			MPCFcluster::DumpHDF5_vorticity< G > (*grid, true, step_id, streamer_vort.str(), path);
#endif
			t = t_ssmpi->get_time();
            step_id = t_ssmpi->get_stepid();
		}
		else
            _ic(*grid);
        
        if (parser("-ic").asBool(0)==1 && !bRESTART) 
            solve_possion();
	}	
    
	void solve_possion()
    {
		PoissonSolver_MPI< G > ps_mpi(*grid, Simulation_Environment::GAMMA1);
        ps_mpi.solve_poisson();
        t_ssmpi->save(*grid, 0, t);
		
		if (bVP)
			vp(*grid, 0, 1);
		else
		{
			const string path = parser("-fpath").asString(".");
			std::stringstream streamer_vort;
			streamer_vort<<"vorticityIC-"<<step_id;
#ifdef _USE_HDF_
			MPCFcluster::DumpHDF5_vorticity< G > (*grid, true, step_id, streamer_vort.str(), path);
#endif
			
			if (isroot)
			{
				if (bASCIIFILES)
					SerializerIO<FluidGrid, StreamerGridPointASCII>().Write(*grid, "restart");
				else 
					SerializerIO<FluidGrid, StreamerGridPoint>().Write(*grid, "restart");
			}
        }
		
		integrals(*grid, 0, 0, true);
		
		if (isroot) cout << "Initial condition is ready for a restart. Make sure you have used a 1D MPI topology for FFTW and change it to a better MPI topology when restaring." << endl;
        MPI_Finalize();
		exit(0);
    }
    
	void run()
	{
        if (isroot) printf("HELLO RUN\n");
        bool bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);       
        
        while(bLoop)
	    {
            if (isroot) printf("Step id %d,Time %f\n", step_id, t);
            
            profiler.push_start("STATISTICS");
            if (step_id%10==0) integrals(*grid, step_id, t);
            profiler.pop_stop();
            
            if (step_id%DUMPPERIOD==0)
            {
				if (bVP)
					vp(*grid, step_id, bVP);
				else
				{
					const string path = parser("-fpath").asString(".");
					std::stringstream streamer_vort;
					streamer_vort<<"vorticity-"<<step_id;
#ifdef _USE_HDF_
					MPCFcluster::DumpHDF5_vorticity< G > (*grid, false, step_id, streamer_vort.str(), path);
#endif
				}
				
                std::stringstream streamer;
                streamer<<"data-"<<step_id;
                t_ssmpi->dump(*grid, step_id, streamer.str());
                t_ssmpi->vp(*grid, step_id, bVP);
            }
            
            if (step_id%SAVEPERIOD==0) t_ssmpi->save(*grid, step_id, t);
            
            const Real dt = (*stepper)(TEND-t);
            
            if(step_id%REPORT_FREQ == 0 && isroot && step_id > 0)
                profiler.printSummary();
            
            t+=dt;
            step_id++;
            
            bLoop = (NSTEPS>0) ? (step_id<NSTEPS) : (fabs(t-TEND) > std::numeric_limits<Real>::epsilon()*1e1);
		}
		
		if (bVP)
			vp(*grid, step_id, bVP);
		else
		{
			const string path = parser("-fpath").asString(".");
			std::stringstream streamer_vort;
			streamer_vort<<"vorticity-"<<step_id;
#ifdef _USE_HDF_
			MPCFcluster::DumpHDF5_vorticity< G > (*grid, true, step_id, streamer_vort.str(), path);
#endif
        }
        std::stringstream streamer;
        streamer<<"data-"<<step_id;;
        t_ssmpi->dump(*grid, step_id, streamer.str());
        
        if (isroot) printf("Finishing RUN\n");
        MPI_Finalize();
        exit(0);
	}
    
	void integrals(G& grid, int iCounter, Real time, bool bIC=false)
	{    
        const string path = parser("-fpath").asString(".");
        vector<BlockInfo> vInfo = grid.getResidentBlocksInfo();
        vector<BlockInfo> vInfo_global = grid.getBlocksInfo();
        double rInt=0, uInt=0, vInt=0, wInt=0, sInt=0, omegaMin=HUGE_VAL, omegaMax=-HUGE_VAL, densityMin=HUGE_VAL, densityMax=-HUGE_VAL, omegaMaxPiS=-HUGE_VAL, omegaMaxPiD=-HUGE_VAL;
        double gammaPiS=0, gammaPiD=0, kenergy=0, kenstrophy=0, veldiv=-HUGE_VAL;
        double x[3];        
        
        int peidx[3];
        grid.peindex(peidx);
        
        EvaluateVorticityMPI< G, EvaluateVorticity_CPP, LabMPI > eval_vort(grid);        
        if (bIC)
            EvaluateVorticityMPI< G, EvaluateVorticity_CPP, LabMPI_PossionSolver > eval_vort(grid);
        eval_vort.execute();
        
        EvaluateVelocityDivergence<G, Evaluate_VelocityDivergence_CPP, LabMPI > eval_veldiv(grid);
        eval_veldiv.execute();
        
        const double h  = vInfo_global[0].h_gridpoint;
        const double h2 = h*h;
        const double h3 = h2*h;
        
        for(int i=0; i<vInfo.size(); i++)
	    {
            const BlockInfo& info = vInfo[i];
            const BlockInfo& info_global = vInfo_global[i];
            FluidBlock& block = *(FluidBlock *)info.ptrBlock;
            
            for(int iz=0; iz<FluidBlock::sizeZ; iz++)
                for(int iy=0; iy<FluidBlock::sizeY; iy++)
                    for(int ix=0; ix<FluidBlock::sizeX; ix++)
                    {
                        info_global.pos(x, ix, iy, iz);
                        
                        rInt +=        block(ix,iy,iz).rho*h3;
                        uInt +=        block(ix,iy,iz).u*h3;
                        vInt +=        block(ix,iy,iz).v*h3;
                        wInt +=        block(ix,iy,iz).w*h3;
                        sInt +=        block(ix,iy,iz).energy*h3;
                        
                        const Real  vorticity_magnitude =  sqrt(pow(block.tmp[iz][iy][ix][0],2)+pow(block.tmp[iz][iy][ix][1],2)+pow(block.tmp[iz][iy][ix][2],2));
                        omegaMin = min((Real)omegaMin,(Real)vorticity_magnitude);
                        omegaMax = max((Real)omegaMax, (Real)vorticity_magnitude);
                        
                        const bool bLastZ = info_global.index[2]==ZPESIZE*BPDZ-1 && iz==FluidBlock::sizeZ-1;
                        const bool bLastX = info_global.index[0]==XPESIZE*BPDX-1 && ix==FluidBlock::sizeX-1;
                                                
                        if (bLastX)
                        {
                            omegaMaxPiS = max((Real)omegaMaxPiS, (Real)vorticity_magnitude);
                            gammaPiS += block.tmp[iz][iy][ix][0];
                        }
                        
                        if (bLastZ)
                        {
                            omegaMaxPiD = max((Real)omegaMaxPiD, (Real)vorticity_magnitude);
                            gammaPiD += block.tmp[iz][iy][ix][2];
                        }
                        
                        densityMin = min((Real)densityMin, (Real)block(ix,iy,iz).rho);
                        densityMax = max((Real)densityMax, (Real)block(ix,iy,iz).rho);
                        
                        kenergy    += (pow(block(ix,iy,iz).u,2)+pow(block(ix,iy,iz).v,2)+pow(block(ix,iy,iz).w,2))/pow(block(ix,iy,iz).rho,2) * h3;
                        kenstrophy += pow(vorticity_magnitude,2) * h3;
                        
                        veldiv = max((Real)veldiv, (Real)block.tmp[iz][iy][ix][3]);
                        
                        block.tmp[iz][iy][ix][0]=0;
                        block.tmp[iz][iy][ix][1]=0;
                        block.tmp[iz][iy][ix][2]=0;
                        block.tmp[iz][iy][ix][3]=0;
                    }
	    }
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        double rIntGlobal, uIntGlobal, vIntGlobal, wIntGlobal, sIntGlobal, omegaMinGlobal, omegaMaxGlobal, densityMinGlobal, densityMaxGlobal, omegaMaxPiSGlobal, omegaMaxPiDGlobal;
        double gammaPiSGlobal, gammaPiDGlobal, kenergyGlobal, kenstrophyGlobal, veldivGlobal;
		
#ifdef _USE_FFTW_
        MPI_Allreduce(&rInt, &rIntGlobal, 1, MPI_FP, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&uInt, &uIntGlobal, 1, MPI_FP, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&vInt, &vIntGlobal, 1, MPI_FP, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&wInt, &wIntGlobal, 1, MPI_FP, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&sInt, &sIntGlobal, 1, MPI_FP, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&omegaMin, &omegaMinGlobal, 1, MPI_FP, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&omegaMax, &omegaMaxGlobal, 1, MPI_FP, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&densityMin, &densityMinGlobal, 1, MPI_FP, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&densityMax, &densityMaxGlobal, 1, MPI_FP, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&omegaMaxPiS, &omegaMaxPiSGlobal, 1, MPI_FP, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&omegaMaxPiD, &omegaMaxPiDGlobal, 1, MPI_FP, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&gammaPiS, &gammaPiSGlobal, 1, MPI_FP, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&gammaPiD, &gammaPiDGlobal, 1, MPI_FP, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&kenergy, &kenergyGlobal, 1, MPI_FP, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&kenstrophy, &kenstrophyGlobal, 1, MPI_FP, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&veldiv, &veldivGlobal, 1, MPI_FP, MPI_MAX, MPI_COMM_WORLD);
#endif
        if (step_id==0)
        {
            gamma_max = gammaPiDGlobal;
            
            const string restart_status = path+"/restart_omega.status";
            ofstream status(restart_status.c_str());
            assert(status.good());
            status << gamma_max;
        }
        else if (bRESTART)
        {
            const string restart_status = path+"/restart_omega.status";
            ifstream status(restart_status.c_str());
            assert(status.good());
            status >> gamma_max;
        }
        else
            cout << "**** Maximum inital vorticity is not set. Either it is not t=0 or no data was available to the restarted run" << endl;

        if (isroot) 
	    {
            FILE * f = fopen("integrals.dat", "a");
            FILE * f2 = fopen("ranges.dat", "a");
            
            fprintf(f,"%d %e %e %e %e %e %e %e %e %e %e\n", iCounter, time, rIntGlobal, uIntGlobal, vIntGlobal, wIntGlobal, sIntGlobal, gammaPiSGlobal/gamma_max, gammaPiDGlobal/gamma_max, kenergyGlobal, kenstrophyGlobal);
            fprintf(f2,"%d %e %e %e %e %e %e %e %e\n", iCounter, time, densityMinGlobal, densityMaxGlobal, omegaMinGlobal, omegaMaxGlobal, omegaMaxPiSGlobal, omegaMaxPiDGlobal, veldivGlobal);
            
            fclose(f);
            fclose(f2);
	    }
	}
	
	void vp(G& grid, const int step_id, const bool bVP)
	{
		if (bVP)
		{
			if (isroot) cout << "dumping MPI VP ..." ;
			
			const string path = parser("-fpath").asString(".");
			std::stringstream streamer;
			streamer<<"vorticity";
			streamer.setf(ios::dec | ios::right);
			streamer.width(5);
			streamer.fill('0');
			streamer<<MPI::COMM_WORLD.Get_rank();
			streamer<<"_";
			streamer.setf(ios::dec | ios::right);
			streamer.width(5);
			streamer.fill('0');
			streamer<<step_id;
			
			DumpVP_VorticityMagnitudeMPI<G>(grid,streamer.str(), path);
			
			if (isroot) cout << "done" << endl;
		}
	}
};

template<typename TGrid>
void DumpVP_VorticityMagnitudeMPI(TGrid& grid, string fileName, const string dump_path)
{
    EvaluateVorticityMPI< TGrid, EvaluateVorticity_CPP, LabMPI > eval_vort(grid);
    eval_vort.execute();
	
	typedef typename GridType::BlockType TBlock;
	
	vector<BlockInfo> vInfo = grid.getResidentBlocksInfo();
	vector<BlockInfo> vInfo_global = grid.getBlocksInfo();
    
	char fullnametxt[256];
	char fullnamegrid[256];
	sprintf(fullnametxt, "%s/%s.txt", dump_path.c_str(), fileName.c_str());
	sprintf(fullnamegrid, "%s/%s.grid", dump_path.c_str(), fileName.c_str());
	
	{
		FILE * file = fopen(fullnametxt, "w");
		assert(file!=NULL);
		fprintf(file, "sizeof(Matrix4D): %d\n",(int)sizeof(Matrix4D<float,true,std::allocator>)); 
		fprintf(file, "Blocks: %d\n", (int)vInfo.size());
		fprintf(file, "Block size: %d %d %d\n", TBlock::sizeX, TBlock::sizeY, TBlock::sizeZ);
		fprintf(file, "Blocks per dim per cpu: %d %d %d\n", grid.getResidentBlocksPerDimension(0), grid.getResidentBlocksPerDimension(1), grid.getResidentBlocksPerDimension(2));
		
		for(int i=0; i<vInfo.size(); i++)
		{
			BlockInfo& info = vInfo[i];
			BlockInfo& info_global = vInfo_global[i];
			fprintf(file, "Block %d: Tree Index: %d %d %d %d %d %d\n", i, info_global.index[0], info_global.index[1], info_global.index[2], info.index[0], info.index[1], info.index[2]); 
		}
		
		fclose(file);
	}
	
	{
		Matrix4D<float, true, std::allocator> * matData = new Matrix4D<float,true,std::allocator>(1, TBlock::sizeX, TBlock::sizeY, TBlock::sizeZ);
		
		FILE * file = fopen(fullnamegrid, "wb");
		assert(file!=NULL);
		
		for(int i=0; i<vInfo.size(); i++)
		{
			BlockInfo& info = vInfo[i];
			TBlock & b = *(TBlock *)info.ptrBlock;
			
			for(int iz=0; iz<TBlock::sizeZ; iz++)
				for(int iy=0; iy<TBlock::sizeY; iy++)
					for(int ix=0; ix<TBlock::sizeX; ix++)
						matData->Access(0, ix, iy, iz) = sqrt(pow(b.tmp[iz][iy][ix][0],2)+pow(b.tmp[iz][iy][ix][1],2)+pow(b.tmp[iz][iy][ix][2],2));
			
			matData->Serialize(file);
		}
		
		fclose(file);
		delete matData;
	}
}

#ifdef _USE_HDF_
template<typename TGrid>
void MPCFcluster::DumpHDF5_vorticity(TGrid &grid, const bool bIC, const int iCounter, const string f_name, const string dump_path)
{
    EvaluateVorticityMPI< TGrid, EvaluateVorticity_CPP, LabMPI > eval_vort(grid);
    
    if (bIC)
        EvaluateVorticityMPI< TGrid, EvaluateVorticity_CPP, LabMPI_PossionSolver > eval_vort(grid);
    
    eval_vort.execute();
    
    typedef typename TGrid::BlockType B;
    
    int rank;
    char filename[256];
    herr_t status;
    hid_t file_id, dataset_id, fspace_id, fapl_id, mspace_id;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int coords[3];
    grid.peindex(coords);
    
    Real array_all[grid.getResidentBlocksPerDimension(0)*_BLOCKSIZEX_][grid.getResidentBlocksPerDimension(1)*_BLOCKSIZEY_][grid.getResidentBlocksPerDimension(2)*_BLOCKSIZEZ_][3];
    
    vector<BlockInfo> vInfo_local = grid.getResidentBlocksInfo();
    
    const int sX = 0;
    const int sY = 0;
    const int sZ = 0;
    
    const int eX = _BLOCKSIZEX_;
    const int eY = _BLOCKSIZEY_;
    const int eZ = _BLOCKSIZEZ_;
    
    hsize_t count[4] = {
        grid.getResidentBlocksPerDimension(0)*_BLOCKSIZEX_,
        grid.getResidentBlocksPerDimension(1)*_BLOCKSIZEY_,
        grid.getResidentBlocksPerDimension(2)*_BLOCKSIZEZ_,3};
    
    hsize_t dims[4] = {
        grid.getBlocksPerDimension(0)*_BLOCKSIZEX_,
        grid.getBlocksPerDimension(1)*_BLOCKSIZEY_,
        grid.getBlocksPerDimension(2)*_BLOCKSIZEZ_, 3};
    
    hsize_t offset[4] = {
        coords[0]*grid.getResidentBlocksPerDimension(0)*_BLOCKSIZEX_,
        coords[1]*grid.getResidentBlocksPerDimension(1)*_BLOCKSIZEY_,
        coords[2]*grid.getResidentBlocksPerDimension(2)*_BLOCKSIZEZ_, 0};
    
    sprintf(filename, "%s/%s.h5", dump_path.c_str(), f_name.c_str());
    
    H5open();
    fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    status = H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL);
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
                    array_all[idx[0]*B::sizeX+ix][idx[1]*B::sizeY+iy][idx[2]*B::sizeZ+iz][0] = b.tmp[iz][iy][ix][0];		
                    array_all[idx[0]*B::sizeX+ix][idx[1]*B::sizeY+iy][idx[2]*B::sizeZ+iz][1] = b.tmp[iz][iy][ix][1];		
                    array_all[idx[0]*B::sizeX+ix][idx[1]*B::sizeY+iy][idx[2]*B::sizeZ+iz][2] = b.tmp[iz][iy][ix][2];
                }
    }
    
    fapl_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);
    
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
    
    if (rank==0)
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
        fprintf(xmf, "        %e %e %e\n", 1./(double)max(dims[0],max(dims[1],dims[2])),1./(double)max(dims[0],max(dims[1],dims[2])),1./(double)max(dims[0],max(dims[1],dims[2])));
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
}
#endif