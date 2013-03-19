#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

namespace WaveletsOnInterval 
{
	template<bool lifting>
	struct WI4
	{
		const Real interp_first(const Real f0, const Real f1, const Real f2, const Real f3) 
		{
			return 5./16 * f0 + 15./16 * f1 - 5./16 * f2 + 1./16 * f3;
		};
		
		const Real interp_middle(const Real f0, const Real f1, const Real f2, const Real f3) 
		{
			return -1./16 * f0 + 9./16 * f1 + 9./16 * f2 - 1./16 * f3;
		};
		
		const Real interp_onetolast(const Real f0, const Real f1, const Real f2, const Real f3) 
		{
			return 1./16 * f0 - 5./16 * f1 + 15./16 * f2 + 5./16 * f3;
		};
		
		const Real interp_last(const Real f0, const Real f1, const Real f2, const Real f3) 
		{
			return -5./16*f0 + 21./16*f1 -35./16*f2 + 35./16*f3;
		};
		
		vector<Real> forward_sweep(const Real * const data, const int N)
		{			
			assert(N >= 8);
			assert(N%2==0);
			
			const int Nhalf = N / 2;
			vector<Real> details(Nhalf);
			
			// compute first detail
			details.front() = data[1] - interp_first(data[0], data[2], data[4], data[6]);
			
			// compute middle details
			for(int i = 1; i < Nhalf - 2; ++i)
			{
				const int s = 2 * i;
				details[i] = data[2 * i + 1] - interp_middle(data[s - 2], data[s], data[s + 2], data[s + 4]); 
			}
			
			//compue one to the last
			details[Nhalf-2] = data[N-3]- interp_onetolast(data[N-8],data[N-6],data[N-4],data[N-2]);
			
			// compute last detail
			details.back() = data[N-1] - interp_last(data[N-8], data[N-6], data[N-4], data[N-2]);
			
			vector<Real> scalings(Nhalf);
			for(int i=0; i<Nhalf; i++)
				scalings[i] = data[2*i];
			
			if (lifting)
				for(int i=0; i<Nhalf; i++)
					scalings[i] += 0.5 * details[i];
			
			vector<Real> retval = scalings;
			retval.insert(retval.end(), details.begin(), details.end());
			
			return retval;
		}
		
		vector<Real> inverse_sweep(const Real * const fwt_data, const int N)
		{			
			assert(N >= 8);
			assert(N%2==0);
			
			const int Nhalf = N / 2;
			
			vector<Real> scalings(fwt_data,fwt_data+Nhalf);
			vector<Real> details(fwt_data+Nhalf, fwt_data+N);
			
			vector<Real> refined(2*Nhalf);
			
			if (lifting)
				for(int i=0; i<Nhalf; i++)
					scalings[i] -= 0.5 * details[i];
			
			for(int i=0; i<Nhalf; i++)
				refined[2*i] = scalings[i];// - 0.5 * details[i];
			
			refined[1] = interp_first(scalings[0], scalings[1], scalings[2], scalings[3]) + details[0];
			
			for(int i=1; i<Nhalf-2; i++)
				refined[2*i+1] = interp_middle(scalings[i-1],scalings[i],scalings[i+1],scalings[i+2]) + details[i];
			
			refined[N-3] = interp_onetolast(scalings[Nhalf-4], scalings[Nhalf-3], scalings[Nhalf-2], scalings[Nhalf-1]) + details[Nhalf-2];
			
			refined[N-1] = interp_last(scalings[Nhalf-4], scalings[Nhalf-3], scalings[Nhalf-2], scalings[Nhalf-1]) + details[Nhalf-1];
			
			return refined;
		}
		
		template<int BS, bool bForward>
		void sweep2D(Real data[BS][BS])
		{
//			auto sweep1D = [&] {
#define sweep1D()	do {	\
				for(int iy=0; iy<BS; iy++) \
				{ \
					vector<Real> temp = forward_sweep(&data[iy][0], BS); \
					copy(temp.begin(), temp.end(), &data[iy][0]);  \
				} \
			} while(0)
			
//			auto sweep1D_inv = [&] {
#define sweep1D_inv()	do { \
				for(int iy=0; iy<BS; iy++) \
				{ \
					vector<Real> temp = inverse_sweep(&data[iy][0], BS); \
					copy(temp.begin(), temp.end(), &data[iy][0]); \
				} \
			} while(0)
			
//			auto transpose = [&] {
#define transpose()	do { \
				for(int iy=0; iy<BS; iy++) \
					for(int ix=iy+1; ix<BS; ix++) \
					{ \
						const Real temp = data[iy][ix]; \
						data[iy][ix] = data[ix][iy]; \
						data[ix][iy] = temp; \
					} \
			} while(0)
			
//			auto sweep = [&] {
#define	sweep()		do { \ 
				if (bForward) \
					sweep1D(); \
				else \
					sweep1D_inv(); \
			} while(0)
			
			sweep();
			transpose();
			sweep();
		}
		
		template<int BS, bool bForward>
		void sweep3D(Real data[BS][BS][BS])
		{			
//			auto xz_transpose = [&] {
#define xz_transpose()	do { \
				for(int iy=0; iy<BS; iy++)	 \
					for(int iz=0; iz<BS; iz++) \
						for(int ix=iz+1; ix<BS; ix++) \
						{ \
							const Real temp = data[iz][iy][ix]; \
							data[iz][iy][ix] = data[ix][iy][iz]; \
							data[ix][iy][iz] = temp; \
						} \
			} while(0)
			
			if(bForward)
			{
				for(int iz=0; iz<BS; iz++)
					sweep2D<BS, true>(data[iz]);
				
				xz_transpose();
				
				for(int iz=0; iz<BS; iz++)
					for(int iy=0; iy<BS; iy++)
					{
						vector<Real> temp = forward_sweep(&data[iz][iy][0], BS);
						copy(temp.begin(), temp.end(), &data[iz][iy][0]);
					}
			}
			else
			{
				for(int iz=0; iz<BS; iz++)
					for(int iy=0; iy<BS; iy++)
					{
						vector<Real> temp = inverse_sweep(&data[iz][iy][0], BS);
						copy(temp.begin(), temp.end(), &data[iz][iy][0]);
					}	
				
				xz_transpose();
				
				for(int iz=0; iz<BS; iz++)
					sweep2D<BS, false>(data[iz]);
			}
		}
		
		template<int BS>
		Real getMaxDetail(const Real data[BS][BS][BS]) 
		{	
			Real mymax = 0;
			const int Nhalf = BS/2;
			
			for(int iz=0; iz<BS; iz++)
				for(int iy=0; iy<BS; iy++)
					for(int ix=0; ix<BS; ix++)
					{
						const Real outside = (Real)(!(ix<Nhalf && iy<Nhalf && iz<Nhalf));
						
						mymax = max(mymax, outside*abs(data[iz][iy][ix]));
					}
			
			return mymax;
		}
	};
	
	template<int BS, bool lifting>
	struct FullTransform : WI4<lifting>
	{
		Real data[BS][BS][BS];
		
		FullTransform<BS/2, lifting> child;
		
		void fwt()
		{
			this->template sweep3D<BS, true>(data);
			
			static const int BSH = BS/2;
			for(int iz = 0; iz < BSH; ++iz)
				for(int iy = 0; iy < BSH; ++iy)
					for(int ix = 0; ix < BSH; ++ix)
					{
						child.data[iz][iy][ix] = data[iz][iy][ix];
						//data[iz][iy][ix] = 0;
					}
			
			child.fwt();
		}
		
		void inverse()
		{
			child.inverse();
			
			static const int BSH = BS/2;
			for(int iz = 0; iz < BSH; ++iz)
				for(int iy = 0; iy < BSH; ++iy)
					for(int ix = 0; ix < BSH; ++ix)
						data[iz][iy][ix] = child.data[iz][iy][ix];
			
			this->template sweep3D<BS, false>(data);
		}
		
		void print()
		{
			printf("**** PRINTING AT BS %d ****\n", BS);
			
			for(int iz=0; iz<BS; iz++)
			{
				printf("Slice ID %d\n", iz);
				
				for(int iy=0; iy<BS; iy++)
				{
					for(int ix=0; ix<BS; ix++)
						printf("%1.2e ", data[iz][iy][ix]);
					
					printf("\n");
				}
			}
			
			child.print();
		}
		
		vector<vector<Real> > collect()
		{
			static const int BSH = BS / 2;
			
			vector<vector<Real> > other = child.collect();
			
			vector<Real> retval;
			
			for(int code = 1; code < 8; ++code)
			{
				const int xs = BSH * (code & 1);
				const int ys = BSH * (code / 2 & 1);
				const int zs = BSH * (code / 4 & 1);
				
				for(int iz = 0; iz < BSH; ++iz)
					for(int iy = 0; iy < BSH; ++iy)
						for(int ix = 0; ix < BSH; ++ix)
							retval.push_back(data[zs+iz][ys+iy][xs+ix]);
				
			}
			
			other.push_back(retval);
			
			return other;
		}
		
		void load(vector<Real>& datastream, bitset<BS * BS * BS> mask)
		{
			static const int BSH = BS / 2;
			
			for(int code = 7; code >= 1; --code)
			{
				const int xstart = BSH * (code & 1);
				const int ystart = BSH * (code / 2 & 1);
				const int zstart = BSH * (code / 4 & 1);
				
				for(int iz = BSH - 1; iz >= 0; --iz)
					for(int iy = BSH - 1; iy >= 0; --iy)
						for(int ix = BSH - 1; ix >= 0; --ix)
						{
							const int myx = xstart + ix;
							const int myy = ystart + iy;
							const int myz = zstart + iz;
							
							const int srcidx = myx + BS * (myy + BS * myz);
							
							assert(srcidx >= 0);
							assert(srcidx < BS * BS * BS);
							
							const bool eat = mask[srcidx];
							
							const Real mydata = eat ? datastream.back() : 0;
							
							data[myz][myy][myx] = mydata;
							
							if (eat) 
								datastream.pop_back();						
						}
			}
			
			//code 0
			{
				bitset<BSH * BSH * BSH> childmask;
				
				for(int iz = 0; iz < BSH; ++iz)
					for(int iy = 0; iy < BSH; ++iy)
						for(int ix = 0; ix < BSH; ++ix)
						{
							const int src = ix + BS * (iy + BS * iz);
							const int dst = ix + BSH * (iy + BSH * iz);
							
							childmask[dst] = mask[src];
						}
				
				child.load(datastream, childmask);
			}
		}
		
		pair<vector<Real>, bitset<BS * BS * BS> > threshold(const Real eps)
		{
			static const int BSH = BS / 2;
			
			pair<vector<Real>, bitset<BS * BS * BS> > retval;
			
			//code 0
			{
				pair<vector<Real>, bitset<BSH * BSH * BSH> > childretval = child.threshold(eps);
				
				retval.first = childretval.first;
				
				for(int iz = 0; iz < BSH; ++iz)
					for(int iy = 0; iy < BSH; ++iy)
						for(int ix = 0; ix < BSH; ++ix)
						{
							const int dst = ix + BS * (iy + BS * iz);
							const int src = ix + BSH * (iy + BSH * iz);
							
							retval.second[dst] = childretval.second[src];
						}
			}
			
			for(int code = 1; code < 8; ++code)
			{
				const int xstart = BSH * (code & 1);
				const int ystart = BSH * (code / 2 & 1);
				const int zstart = BSH * (code / 4 & 1);
				
				for(int iz = 0; iz < BSH; ++iz)
					for(int iy = 0; iy < BSH; ++iy)
						for(int ix = 0; ix < BSH; ++ix)
						{
							const int xsrc = xstart + ix;
							const int ysrc = ystart + iy;
							const int zsrc = zstart + iz;
							
							const int dst = xsrc + BS * (ysrc + BS * zsrc);
							
							const Real mydata = data[zsrc][ysrc][xsrc];
							
							const bool accepted = fabs(mydata) > eps; 
							
							retval.second[dst] = accepted;
							
							if (accepted)
								retval.first.push_back(mydata);
							else 
								data[zsrc][ysrc][xsrc] = 0;
						}
			}
			
			return retval;
		}
	};
	
	template<bool lifting>
	struct FullTransform<4, lifting>
	{
		enum { BS = 4 } ;
		
		Real data[BS][BS][BS];
		
		void fwt() { }
		
		void inverse() { }
		
		pair<vector<Real>, bitset<BS * BS * BS> > threshold(const Real eps) 
		{ 
			enum { N = BS * BS * BS };
			
			const Real * const e = &data[0][0][0];
			
			vector<Real> v(N);
			
			for(int i = 0; i < N; ++i)
				v[i] = e[i];
			
			pair<vector<Real>, bitset<BS * BS * BS> > retval;
			retval.first = v;
			retval.second.set();
			
			return retval; 
		}
		
		void load(vector<Real>& datastream, bitset<BS * BS * BS> mask)
		{
			assert(datastream.size() == BS * BS * BS);
			
			for(int iz = 0, c = 0; iz < BS; ++iz)
				for(int iy = 0; iy < BS; ++iy)
					for(int ix = 0; ix < BS; ++ix, ++c)
						data[iz][iy][ix] = datastream[c];
			
			datastream.clear();
		}
		
		vector<vector<Real> > collect()
		{
			vector<Real> retval;
			
			for(int iz = 0; iz < BS; ++iz)
				for(int iy = 0; iy < BS; ++iy)
					for(int ix = 0; ix < BS; ++ix)
						retval.push_back(data[iz][iy][ix]);
			
			return vector<vector<Real> >(1,retval);
		}
		
		void print()
		{
			printf("**** LAST LEVEL! ****\n");
			
			for(int iz=0; iz<BS; iz++)
			{
				printf("Slice ID %d\n", iz);
				
				for(int iy=0; iy<BS; iy++)
				{
					for(int ix=0; ix<BS; ix++)
						printf("%1.2e ", data[iz][iy][ix]);
					
					printf("\n");
				}
			}
		}
	};
}
