/*
 *  SliceSOA.h
 *  
 *
 *  Created by Diego Rossinelli on 5/15/12.
 *  Copyright 2012 ETH Zurich. All rights reserved.
 *
 */

template < int _SX, int _EX, int _SY, int _EY, typename TReal=Real > 
struct SOA2D
{
	static const int _CPERALIGNBYTES = _ALIGNBYTES_/sizeof(TReal);
	
	static const int SX = _CPERALIGNBYTES*((_SX - (_CPERALIGNBYTES-1))/_CPERALIGNBYTES);
	static const int EX = _CPERALIGNBYTES*((_EX + (_CPERALIGNBYTES-1))/_CPERALIGNBYTES);
	static const int NX = _EX - _SX;
	static const int NY = _EY - _SY;
	
	static const int PITCH = EX - SX;
	
	TReal __attribute__((aligned(_ALIGNBYTES_))) data[NY][PITCH];
	
	SOA2D()
	{
		assert(((size_t)(&data[0][0]) % _ALIGNBYTES_) == 0);
	}
	
	inline TReal operator()(const int ix, const int iy) const
	{
		assert(ix >= _SX); assert(ix < _EX);
		assert(iy >= _SY); assert(iy < _EY);
		
		return data[iy-_SY][ix-SX];
	}
	
	inline const TReal * ptr(const int ix, const int iy) const
	{
		assert(ix >= SX); assert(ix < EX);
		assert(iy >= _SY); assert(iy < _EY);
		
		return &data[iy-_SY][ix-SX];
	}
	
	inline TReal& ref(const int ix, const int iy)
	{
		assert(ix >= SX); assert(ix < EX);
		assert(iy >= _SY); assert(iy < _EY);
		
		return data[iy-_SY][ix-SX];
	}
	
	static float kB(int nobjects=1)
	{
		return nobjects*sizeof(SOA2D)/1024.;
	}
	
	inline SOA2D& operator= (const SOA2D& c)
	{
		for(int y=0; y<NY; ++y)
			for(int x=0; x<PITCH; ++x)
				data[y][x] = c.data[y][x];
		
		return *this;
	}
};

template< int _SX, int _EX, int _SY, int _EY, int _NSLICES, typename TReal=Real>
struct RingSOA2D
{
	int currslice;
	
	SOA2D<_SX, _EX, _SY, _EY, TReal> slices[_NSLICES];
	
	RingSOA2D(): currslice(0){ }
	
	inline const SOA2D<_SX, _EX, _SY, _EY, TReal>& operator()(const int relativeid=0) const
	{
		return slices[(relativeid + currslice + _NSLICES) % _NSLICES];
	}
	
	inline SOA2D<_SX, _EX, _SY, _EY, TReal>& ref(const int relativeid=0)
	{
		return slices[(relativeid + currslice + _NSLICES) % _NSLICES];
	}
	
	void next(){ currslice = (currslice + 1) % _NSLICES; } 
	
	static float kB(int nobjects=1)
	{
		return nobjects*sizeof(RingSOA2D)/1024.;
	}
};

typedef SOA2D<0, _BLOCKSIZE_, 0, _BLOCKSIZE_> OutputSOA;
