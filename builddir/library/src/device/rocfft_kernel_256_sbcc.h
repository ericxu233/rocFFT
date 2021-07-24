#pragma once
#include "rocfft_butterfly_template.h"

////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb, bool TwdLarge>
__device__ void
FwdPass0_len256_sbcc(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	FwdRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 1 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 2 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 3 ) ] = (*R3);
	}

}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void
FwdPass1_len256_sbcc(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	{
		T W = twiddles[3 + 3*((1*me + 0)%4) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R1).x; ry = (*R1).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R1).x = TR;
		(*R1).y = TI;
	}

	{
		T W = twiddles[3 + 3*((1*me + 0)%4) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R2).x; ry = (*R2).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R2).x = TR;
		(*R2).y = TI;
	}

	{
		T W = twiddles[3 + 3*((1*me + 0)%4) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	FwdRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 4 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 8 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 12 ) ] = (*R3);
	}

}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void
FwdPass2_len256_sbcc(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	{
		T W = twiddles[15 + 3*((1*me + 0)%16) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R1).x; ry = (*R1).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R1).x = TR;
		(*R1).y = TI;
	}

	{
		T W = twiddles[15 + 3*((1*me + 0)%16) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R2).x; ry = (*R2).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R2).x = TR;
		(*R2).y = TI;
	}

	{
		T W = twiddles[15 + 3*((1*me + 0)%16) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	FwdRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 16 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 32 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 48 ) ] = (*R3);
	}

}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void
FwdPass3_len256_sbcc(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	{
		T W = twiddles[63 + 3*((1*me + 0)%64) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R1).x; ry = (*R1).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R1).x = TR;
		(*R1).y = TI;
	}

	{
		T W = twiddles[63 + 3*((1*me + 0)%64) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R2).x; ry = (*R2).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R2).x = TR;
		(*R2).y = TI;
	}

	{
		T W = twiddles[63 + 3*((1*me + 0)%64) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	FwdRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((1*me + 0)%64 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R0).x; ry = (*R0).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R0).x = TR;
		(*R0).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((1*me + 0)%64 + 64) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R1).x; ry = (*R1).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R1).x = TR;
		(*R1).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((1*me + 0)%64 + 128) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R2).x; ry = (*R2).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R2).x = TR;
		(*R2).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((1*me + 0)%64 + 192) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	if(rw)
	{
	bufOut[outOffset + ( 1*me + 0 + 0 ) ] = (*R0);
	bufOut[outOffset + ( 1*me + 0 + 64 ) ] = (*R1);
	bufOut[outOffset + ( 1*me + 0 + 128 ) ] = (*R2);
	bufOut[outOffset + ( 1*me + 0 + 192 ) ] = (*R3);
	}

}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void
InvPass0_len256_sbcc(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	InvRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 1 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 2 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 3 ) ] = (*R3);
	}

}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void
InvPass1_len256_sbcc(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	{
		T W = twiddles[3 + 3*((1*me + 0)%4) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R1).x; ry = (*R1).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R1).x = TR;
		(*R1).y = TI;
	}

	{
		T W = twiddles[3 + 3*((1*me + 0)%4) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R2).x; ry = (*R2).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R2).x = TR;
		(*R2).y = TI;
	}

	{
		T W = twiddles[3 + 3*((1*me + 0)%4) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	InvRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 4 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 8 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 12 ) ] = (*R3);
	}

}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void
InvPass2_len256_sbcc(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	{
		T W = twiddles[15 + 3*((1*me + 0)%16) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R1).x; ry = (*R1).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R1).x = TR;
		(*R1).y = TI;
	}

	{
		T W = twiddles[15 + 3*((1*me + 0)%16) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R2).x; ry = (*R2).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R2).x = TR;
		(*R2).y = TI;
	}

	{
		T W = twiddles[15 + 3*((1*me + 0)%16) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	InvRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 16 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 32 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 48 ) ] = (*R3);
	}

}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void
InvPass3_len256_sbcc(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	{
		T W = twiddles[63 + 3*((1*me + 0)%64) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R1).x; ry = (*R1).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R1).x = TR;
		(*R1).y = TI;
	}

	{
		T W = twiddles[63 + 3*((1*me + 0)%64) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R2).x; ry = (*R2).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R2).x = TR;
		(*R2).y = TI;
	}

	{
		T W = twiddles[63 + 3*((1*me + 0)%64) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	InvRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((1*me + 0)%64 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R0).x; ry = (*R0).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R0).x = TR;
		(*R0).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((1*me + 0)%64 + 64) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R1).x; ry = (*R1).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R1).x = TR;
		(*R1).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((1*me + 0)%64 + 128) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R2).x; ry = (*R2).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R2).x = TR;
		(*R2).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((1*me + 0)%64 + 192) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	if(rw)
	{
	bufOut[outOffset + ( 1*me + 0 + 0 ) ] = (*R0);
	bufOut[outOffset + ( 1*me + 0 + 64 ) ] = (*R1);
	bufOut[outOffset + ( 1*me + 0 + 128 ) ] = (*R2);
	bufOut[outOffset + ( 1*me + 0 + 192 ) ] = (*R3);
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
fwd_len256_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, T *lds)
{
	T R0, R1, R2, R3;
	FwdPass0_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  lwbIn, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass1_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass2_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass3_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  lwbOut, &R0, &R1, &R2, &R3);
	__syncthreads();
}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
back_len256_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, T *lds)
{
	T R0, R1, R2, R3;
	InvPass0_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  lwbIn, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass1_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass2_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass3_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  lwbOut, &R0, &R1, &R2, &R3);
	__syncthreads();
}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
fwd_len256_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3;
	FwdPass0_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  lwbIn, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass1_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass2_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass3_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3);
	__syncthreads();
}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
back_len256_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3;
	InvPass0_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  lwbIn, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass1_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass2_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass3_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3);
	__syncthreads();
}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
fwd_len256_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, T *lds)
{
	T R0, R1, R2, R3;
	FwdPass0_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass1_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass2_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass3_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  lwbOut, &R0, &R1, &R2, &R3);
	__syncthreads();
}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
back_len256_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, T *lds)
{
	T R0, R1, R2, R3;
	InvPass0_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass1_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass2_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass3_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  lwbOut, &R0, &R1, &R2, &R3);
	__syncthreads();
}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
fwd_len256_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3;
	FwdPass0_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass1_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass2_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass3_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3);
	__syncthreads();
}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
back_len256_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3;
	InvPass0_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass1_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass2_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass3_len256_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3);
	__syncthreads();
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(256)
fft_fwd_ip_len256_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2056), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len256_sbcc>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int ioOffset = 0;
	T *lwb;

	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/8);

			ioOffset += (counter_mod/currentLength)*stride_in[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8;

		ioOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			ioOffset *= (stride_in[1]);
		ioOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

	}
	lwb = gb + ioOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0 = lwb[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*32 + (me%8)*256 + (me/8)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/8))*8 + t*4 + (me/64);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len256_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_in[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		lwb[(me%8) + (me/8)*stride_in[0] + t*stride_in[0]*32] = R0;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(256)
fft_back_ip_len256_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2056), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len256_sbcc>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int ioOffset = 0;
	T *lwb;

	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/8);

			ioOffset += (counter_mod/currentLength)*stride_in[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8;

		ioOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			ioOffset *= (stride_in[1]);
		ioOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

	}
	lwb = gb + ioOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0 = lwb[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*32 + (me%8)*256 + (me/8)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/8))*8 + t*4 + (me/64);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len256_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_in[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		lwb[(me%8) + (me/8)*stride_in[0] + t*stride_in[0]*32] = R0;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(256)
fft_fwd_ip_len256_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2056), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len256_sbcc>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int ioOffset = 0;
	real_type_t<T> *lwbRe;
	real_type_t<T> *lwbIm;

	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/8);

			ioOffset += (counter_mod/currentLength)*stride_in[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8;

		ioOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			ioOffset *= (stride_in[1]);
		ioOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

	}
	lwbRe = gbRe + ioOffset;
	lwbIm = gbIm + ioOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0.x = lwbRe[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		R0.y = lwbIm[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*32 + (me%8)*256 + (me/8)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/8))*8 + t*4 + (me/64);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len256_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_in[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		lwbRe[(me%8) + (me/8)*stride_in[0] + t*stride_in[0]*32] = R0.x;
		lwbIm[(me%8) + (me/8)*stride_in[0] + t*stride_in[0]*32] = R0.y;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(256)
fft_back_ip_len256_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2056), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len256_sbcc>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int ioOffset = 0;
	real_type_t<T> *lwbRe;
	real_type_t<T> *lwbIm;

	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/8);

			ioOffset += (counter_mod/currentLength)*stride_in[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8;

		ioOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			ioOffset *= (stride_in[1]);
		ioOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

	}
	lwbRe = gbRe + ioOffset;
	lwbIm = gbIm + ioOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0.x = lwbRe[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		R0.y = lwbIm[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*32 + (me%8)*256 + (me/8)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/8))*8 + t*4 + (me/64);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len256_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_in[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		lwbRe[(me%8) + (me/8)*stride_in[0] + t*stride_in[0]*32] = R0.x;
		lwbIm[(me%8) + (me/8)*stride_in[0] + t*stride_in[0]*32] = R0.y;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(256)
fft_fwd_op_len256_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2056), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len256_sbcc>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	T *lwbOut;

	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbIn = gbIn + iOffset;
	lwbOut = gbOut + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0 = lwbIn[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*32 + (me%8)*256 + (me/8)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/8))*8 + t*4 + (me/64);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len256_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		lwbOut[(me%8) + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(256)
fft_back_op_len256_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2056), cgh);
	cgh.parallel_for<class kern_fft_back_op_len256_sbcc>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	T *lwbOut;

	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbIn = gbIn + iOffset;
	lwbOut = gbOut + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0 = lwbIn[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*32 + (me%8)*256 + (me/8)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/8))*8 + t*4 + (me/64);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len256_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		lwbOut[(me%8) + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(256)
fft_fwd_op_len256_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutRe = gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutIm = gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2056), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len256_sbcc>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	real_type_t<T> *lwbOutRe;
	real_type_t<T> *lwbOutIm;

	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbIn = gbIn + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0 = lwbIn[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*32 + (me%8)*256 + (me/8)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/8))*8 + t*4 + (me/64);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len256_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		lwbOutRe[(me%8) + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.x;
		lwbOutIm[(me%8) + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.y;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(256)
fft_back_op_len256_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutRe = gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutIm = gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2056), cgh);
	cgh.parallel_for<class kern_fft_back_op_len256_sbcc>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	real_type_t<T> *lwbOutRe;
	real_type_t<T> *lwbOutIm;

	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbIn = gbIn + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0 = lwbIn[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*32 + (me%8)*256 + (me/8)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/8))*8 + t*4 + (me/64);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len256_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		lwbOutRe[(me%8) + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.x;
		lwbOutIm[(me%8) + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.y;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(256)
fft_fwd_op_len256_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbInRe = gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbInIm = gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2056), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len256_sbcc>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	real_type_t<T> *lwbInRe;
	real_type_t<T> *lwbInIm;
	T *lwbOut;

	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOut = gbOut + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0.x = lwbInRe[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		R0.y = lwbInIm[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*32 + (me%8)*256 + (me/8)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/8))*8 + t*4 + (me/64);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len256_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		lwbOut[(me%8) + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(256)
fft_back_op_len256_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbInRe = gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbInIm = gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2056), cgh);
	cgh.parallel_for<class kern_fft_back_op_len256_sbcc>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	real_type_t<T> *lwbInRe;
	real_type_t<T> *lwbInIm;
	T *lwbOut;

	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOut = gbOut + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0.x = lwbInRe[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		R0.y = lwbInIm[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*32 + (me%8)*256 + (me/8)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/8))*8 + t*4 + (me/64);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len256_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		lwbOut[(me%8) + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(256)
fft_fwd_op_len256_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbInRe = gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbInIm = gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutRe = gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutIm = gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2056), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len256_sbcc>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	real_type_t<T> *lwbInRe;
	real_type_t<T> *lwbInIm;
	real_type_t<T> *lwbOutRe;
	real_type_t<T> *lwbOutIm;

	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0.x = lwbInRe[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		R0.y = lwbInIm[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*32 + (me%8)*256 + (me/8)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/8))*8 + t*4 + (me/64);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len256_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		lwbOutRe[(me%8) + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.x;
		lwbOutIm[(me%8) + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.y;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(256)
fft_back_op_len256_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbInRe = gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbInIm = gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutRe = gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutIm = gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2056), cgh);
	cgh.parallel_for<class kern_fft_back_op_len256_sbcc>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	real_type_t<T> *lwbInRe;
	real_type_t<T> *lwbInIm;
	real_type_t<T> *lwbOutRe;
	real_type_t<T> *lwbOutIm;

	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0.x = lwbInRe[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		R0.y = lwbInIm[(me%8) * stride_in[1] + (me/8)*stride_in[0] + t*stride_in[0]*32];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*32 + (me%8)*256 + (me/8)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/8))*8 + t*4 + (me/64);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len256_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		lwbOutRe[(me%8) + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.x;
		lwbOutIm[(me%8) + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.y;
	}

});
});
}

