#pragma once
#include "rocfft_butterfly_template.h"
#include "real2complex.h"


////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb>
__device__ void
FwdPass0_len128_sbrc(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 16 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 32 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 48 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*1 + 0 + 80 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*1 + 0 + 96 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*1 + 0 + 112 ) ];
	}



	FwdRad8B1(R0, R1, R2, R3, R4, R5, R6, R7);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 1 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 2 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 3 ) ] = (*R3);
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 4 ) ] = (*R4);
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 5 ) ] = (*R5);
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 6 ) ] = (*R6);
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 7 ) ] = (*R7);
	}

}

template <typename T, StrideBin sb>
__device__ void
FwdPass1_len128_sbrc(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*2 + 0 + 32 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*2 + 1 + 32 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*2 + 0 + 64 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*2 + 1 + 64 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*2 + 0 + 96 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*2 + 1 + 96 ) ];
	}



	{
		T W = twiddles[7 + 3*((2*me + 0)%8) + 0];
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
		T W = twiddles[7 + 3*((2*me + 0)%8) + 1];
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
		T W = twiddles[7 + 3*((2*me + 0)%8) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	{
		T W = twiddles[7 + 3*((2*me + 1)%8) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R5).x; ry = (*R5).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R5).x = TR;
		(*R5).y = TI;
	}

	{
		T W = twiddles[7 + 3*((2*me + 1)%8) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R6).x; ry = (*R6).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R6).x = TR;
		(*R6).y = TI;
	}

	{
		T W = twiddles[7 + 3*((2*me + 1)%8) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R7).x; ry = (*R7).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R7).x = TR;
		(*R7).y = TI;
	}

	FwdRad4B1(R0, R1, R2, R3);
	FwdRad4B1(R4, R5, R6, R7);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((2*me + 0)/8)*32 + (2*me + 0)%8 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((2*me + 0)/8)*32 + (2*me + 0)%8 + 8 ) ] = (*R1);
	bufOut[outOffset + ( ((2*me + 0)/8)*32 + (2*me + 0)%8 + 16 ) ] = (*R2);
	bufOut[outOffset + ( ((2*me + 0)/8)*32 + (2*me + 0)%8 + 24 ) ] = (*R3);
	bufOut[outOffset + ( ((2*me + 1)/8)*32 + (2*me + 1)%8 + 0 ) ] = (*R4);
	bufOut[outOffset + ( ((2*me + 1)/8)*32 + (2*me + 1)%8 + 8 ) ] = (*R5);
	bufOut[outOffset + ( ((2*me + 1)/8)*32 + (2*me + 1)%8 + 16 ) ] = (*R6);
	bufOut[outOffset + ( ((2*me + 1)/8)*32 + (2*me + 1)%8 + 24 ) ] = (*R7);
	}

}

template <typename T, StrideBin sb>
__device__ void
FwdPass2_len128_sbrc(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*2 + 0 + 32 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*2 + 1 + 32 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*2 + 0 + 64 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*2 + 1 + 64 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*2 + 0 + 96 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*2 + 1 + 96 ) ];
	}



	{
		T W = twiddles[31 + 3*((2*me + 0)%32) + 0];
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
		T W = twiddles[31 + 3*((2*me + 0)%32) + 1];
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
		T W = twiddles[31 + 3*((2*me + 0)%32) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	{
		T W = twiddles[31 + 3*((2*me + 1)%32) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R5).x; ry = (*R5).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R5).x = TR;
		(*R5).y = TI;
	}

	{
		T W = twiddles[31 + 3*((2*me + 1)%32) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R6).x; ry = (*R6).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R6).x = TR;
		(*R6).y = TI;
	}

	{
		T W = twiddles[31 + 3*((2*me + 1)%32) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R7).x; ry = (*R7).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R7).x = TR;
		(*R7).y = TI;
	}

	FwdRad4B1(R0, R1, R2, R3);
	FwdRad4B1(R4, R5, R6, R7);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( 2*me + 0 + 0 ) ] = (*R0);
	bufOut[outOffset + ( 2*me + 1 + 0 ) ] = (*R4);
	bufOut[outOffset + ( 2*me + 0 + 32 ) ] = (*R1);
	bufOut[outOffset + ( 2*me + 1 + 32 ) ] = (*R5);
	bufOut[outOffset + ( 2*me + 0 + 64 ) ] = (*R2);
	bufOut[outOffset + ( 2*me + 1 + 64 ) ] = (*R6);
	bufOut[outOffset + ( 2*me + 0 + 96 ) ] = (*R3);
	bufOut[outOffset + ( 2*me + 1 + 96 ) ] = (*R7);
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass0_len128_sbrc(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 16 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 32 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 48 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*1 + 0 + 80 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*1 + 0 + 96 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*1 + 0 + 112 ) ];
	}



	InvRad8B1(R0, R1, R2, R3, R4, R5, R6, R7);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 1 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 2 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 3 ) ] = (*R3);
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 4 ) ] = (*R4);
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 5 ) ] = (*R5);
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 6 ) ] = (*R6);
	bufOut[outOffset + ( ((1*me + 0)/1)*8 + (1*me + 0)%1 + 7 ) ] = (*R7);
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass1_len128_sbrc(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*2 + 0 + 32 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*2 + 1 + 32 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*2 + 0 + 64 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*2 + 1 + 64 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*2 + 0 + 96 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*2 + 1 + 96 ) ];
	}



	{
		T W = twiddles[7 + 3*((2*me + 0)%8) + 0];
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
		T W = twiddles[7 + 3*((2*me + 0)%8) + 1];
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
		T W = twiddles[7 + 3*((2*me + 0)%8) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	{
		T W = twiddles[7 + 3*((2*me + 1)%8) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R5).x; ry = (*R5).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R5).x = TR;
		(*R5).y = TI;
	}

	{
		T W = twiddles[7 + 3*((2*me + 1)%8) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R6).x; ry = (*R6).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R6).x = TR;
		(*R6).y = TI;
	}

	{
		T W = twiddles[7 + 3*((2*me + 1)%8) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R7).x; ry = (*R7).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R7).x = TR;
		(*R7).y = TI;
	}

	InvRad4B1(R0, R1, R2, R3);
	InvRad4B1(R4, R5, R6, R7);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((2*me + 0)/8)*32 + (2*me + 0)%8 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((2*me + 0)/8)*32 + (2*me + 0)%8 + 8 ) ] = (*R1);
	bufOut[outOffset + ( ((2*me + 0)/8)*32 + (2*me + 0)%8 + 16 ) ] = (*R2);
	bufOut[outOffset + ( ((2*me + 0)/8)*32 + (2*me + 0)%8 + 24 ) ] = (*R3);
	bufOut[outOffset + ( ((2*me + 1)/8)*32 + (2*me + 1)%8 + 0 ) ] = (*R4);
	bufOut[outOffset + ( ((2*me + 1)/8)*32 + (2*me + 1)%8 + 8 ) ] = (*R5);
	bufOut[outOffset + ( ((2*me + 1)/8)*32 + (2*me + 1)%8 + 16 ) ] = (*R6);
	bufOut[outOffset + ( ((2*me + 1)/8)*32 + (2*me + 1)%8 + 24 ) ] = (*R7);
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass2_len128_sbrc(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*2 + 0 + 32 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*2 + 1 + 32 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*2 + 0 + 64 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*2 + 1 + 64 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*2 + 0 + 96 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*2 + 1 + 96 ) ];
	}



	{
		T W = twiddles[31 + 3*((2*me + 0)%32) + 0];
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
		T W = twiddles[31 + 3*((2*me + 0)%32) + 1];
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
		T W = twiddles[31 + 3*((2*me + 0)%32) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	{
		T W = twiddles[31 + 3*((2*me + 1)%32) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R5).x; ry = (*R5).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R5).x = TR;
		(*R5).y = TI;
	}

	{
		T W = twiddles[31 + 3*((2*me + 1)%32) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R6).x; ry = (*R6).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R6).x = TR;
		(*R6).y = TI;
	}

	{
		T W = twiddles[31 + 3*((2*me + 1)%32) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R7).x; ry = (*R7).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R7).x = TR;
		(*R7).y = TI;
	}

	InvRad4B1(R0, R1, R2, R3);
	InvRad4B1(R4, R5, R6, R7);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( 2*me + 0 + 0 ) ] = (*R0);
	bufOut[outOffset + ( 2*me + 1 + 0 ) ] = (*R4);
	bufOut[outOffset + ( 2*me + 0 + 32 ) ] = (*R1);
	bufOut[outOffset + ( 2*me + 1 + 32 ) ] = (*R5);
	bufOut[outOffset + ( 2*me + 0 + 64 ) ] = (*R2);
	bufOut[outOffset + ( 2*me + 1 + 64 ) ] = (*R6);
	bufOut[outOffset + ( 2*me + 0 + 96 ) ] = (*R3);
	bufOut[outOffset + ( 2*me + 1 + 96 ) ] = (*R7);
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb>
__device__ void 
fwd_len128_sbrc_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7;
	FwdPass0_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  lwbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	FwdPass1_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	FwdPass2_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
}

template <typename T, StrideBin sb>
__device__ void 
back_len128_sbrc_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7;
	InvPass0_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  lwbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	InvPass1_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	InvPass2_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len128_sbrc_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7;
	FwdPass0_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  lwbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	FwdPass1_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	FwdPass2_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
}

template <typename T, StrideBin sb>
__device__ void 
back_len128_sbrc_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7;
	InvPass0_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  lwbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	InvPass1_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	InvPass2_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len128_sbrc_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7;
	FwdPass0_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	FwdPass1_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	FwdPass2_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
}

template <typename T, StrideBin sb>
__device__ void 
back_len128_sbrc_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7;
	InvPass0_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	InvPass1_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	InvPass2_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len128_sbrc_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7;
	FwdPass0_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	FwdPass1_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	FwdPass2_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
}

template <typename T, StrideBin sb>
__device__ void 
back_len128_sbrc_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7;
	InvPass0_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	InvPass1_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
	InvPass2_len128_sbrc<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7);
	__syncthreads();
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 128, transforms: 4, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose>
__global__ void
__launch_bounds__(128)
fft_fwd_op_len128_sbrc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1032), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len128_sbrc>(
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
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
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset *= (stride_out[1]);
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 128 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 16;
		unsigned int tileBlockIdx_x = ((bid / 16) + tileBlockIdx_y) % 128;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}
	lwbIn = gbIn + iOffset;
	lwbOut = gbOut + oOffset;


	unsigned int lds_row_padding = 0;
	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 128 * 8 * 128 / 128 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = lwbIn[me + t*128];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = lwbIn[me % 128 * stride_in[0] + ((me /128 * 16) + t % 8)*stride_in[2] + t / 8 * 128 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 128 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *128 + t / 8 * 128 + me % 128 + me / 128 * 1024] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len128_sbrc_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, t * (128 + lds_row_padding) * 8 + (me/16)*(128+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true>(me, 128 - me, 128, 64, &lds[r * (128 + lds_row_padding)], &twiddles[128]);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			lwbOut[(me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*16] = R0;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			lwbOut[(me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*16] = R0;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			lwbOut[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

			lwbOut[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0;
		}
	}
});
});
}

//Kernel configuration: number of threads per thread block: 128, transforms: 4, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose>
__global__ void
__launch_bounds__(128)
fft_back_op_len128_sbrc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1032), cgh);
	cgh.parallel_for<class kern_fft_back_op_len128_sbrc>(
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
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
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset *= (stride_out[1]);
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 128 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 16;
		unsigned int tileBlockIdx_x = ((bid / 16) + tileBlockIdx_y) % 128;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}
	lwbIn = gbIn + iOffset;
	lwbOut = gbOut + oOffset;


	unsigned int lds_row_padding = 0;
	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 128 * 8 * 128 / 128 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = lwbIn[me + t*128];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = lwbIn[me % 128 * stride_in[0] + ((me /128 * 16) + t % 8)*stride_in[2] + t / 8 * 128 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 128 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *128 + t / 8 * 128 + me % 128 + me / 128 * 1024] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len128_sbrc_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, t * (128 + lds_row_padding) * 8 + (me/16)*(128+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true>(me, 128 - me, 128, 64, &lds[r * (128 + lds_row_padding)], &twiddles[128]);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			lwbOut[(me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*16] = R0;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			lwbOut[(me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*16] = R0;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			lwbOut[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

			lwbOut[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0;
		}
	}
});
});
}

//Kernel configuration: number of threads per thread block: 128, transforms: 4, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose>
__global__ void
__launch_bounds__(128)
fft_fwd_op_len128_sbrc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutRe = gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutIm = gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1032), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len128_sbrc>(
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
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
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset *= (stride_out[1]);
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 128 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 16;
		unsigned int tileBlockIdx_x = ((bid / 16) + tileBlockIdx_y) % 128;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}
	lwbIn = gbIn + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;


	unsigned int lds_row_padding = 0;
	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 128 * 8 * 128 / 128 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = lwbIn[me + t*128];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = lwbIn[me % 128 * stride_in[0] + ((me /128 * 16) + t % 8)*stride_in[2] + t / 8 * 128 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 128 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *128 + t / 8 * 128 + me % 128 + me / 128 * 1024] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len128_sbrc_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, t * (128 + lds_row_padding) * 8 + (me/16)*(128+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true>(me, 128 - me, 128, 64, &lds[r * (128 + lds_row_padding)], &twiddles[128]);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			lwbOutRe[(me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*16] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			lwbOutRe[(me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*16] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			lwbOutRe[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

			lwbOutRe[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.y;
		}
	}
});
});
}

//Kernel configuration: number of threads per thread block: 128, transforms: 4, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose>
__global__ void
__launch_bounds__(128)
fft_back_op_len128_sbrc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutRe = gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutIm = gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1032), cgh);
	cgh.parallel_for<class kern_fft_back_op_len128_sbrc>(
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
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
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset *= (stride_out[1]);
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 128 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 16;
		unsigned int tileBlockIdx_x = ((bid / 16) + tileBlockIdx_y) % 128;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}
	lwbIn = gbIn + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;


	unsigned int lds_row_padding = 0;
	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 128 * 8 * 128 / 128 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = lwbIn[me + t*128];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = lwbIn[me % 128 * stride_in[0] + ((me /128 * 16) + t % 8)*stride_in[2] + t / 8 * 128 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 128 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *128 + t / 8 * 128 + me % 128 + me / 128 * 1024] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len128_sbrc_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, t * (128 + lds_row_padding) * 8 + (me/16)*(128+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true>(me, 128 - me, 128, 64, &lds[r * (128 + lds_row_padding)], &twiddles[128]);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			lwbOutRe[(me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*16] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			lwbOutRe[(me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*16] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			lwbOutRe[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

			lwbOutRe[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.y;
		}
	}
});
});
}

//Kernel configuration: number of threads per thread block: 128, transforms: 4, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose>
__global__ void
__launch_bounds__(128)
fft_fwd_op_len128_sbrc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbInRe = gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbInIm = gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1032), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len128_sbrc>(
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
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
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset *= (stride_out[1]);
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 128 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 16;
		unsigned int tileBlockIdx_x = ((bid / 16) + tileBlockIdx_y) % 128;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOut = gbOut + oOffset;


	unsigned int lds_row_padding = 0;
	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 128 * 8 * 128 / 128 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = lwbInRe[me + t*128];
			R0.y = lwbInIm[me + t*128];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = lwbInRe[me % 128 * stride_in[0] + ((me /128 * 16) + t % 8)*stride_in[2] + t / 8 * 128 * stride_in[0]];
			R0.y = lwbInIm[me % 128 * stride_in[0] + ((me /128 * 16) + t % 8)*stride_in[2] + t / 8 * 128 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 128 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *128 + t / 8 * 128 + me % 128 + me / 128 * 1024] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len128_sbrc_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, t * (128 + lds_row_padding) * 8 + (me/16)*(128+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true>(me, 128 - me, 128, 64, &lds[r * (128 + lds_row_padding)], &twiddles[128]);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			lwbOut[(me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*16] = R0;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			lwbOut[(me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*16] = R0;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			lwbOut[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

			lwbOut[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0;
		}
	}
});
});
}

//Kernel configuration: number of threads per thread block: 128, transforms: 4, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose>
__global__ void
__launch_bounds__(128)
fft_back_op_len128_sbrc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbInRe = gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbInIm = gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1032), cgh);
	cgh.parallel_for<class kern_fft_back_op_len128_sbrc>(
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
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
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset *= (stride_out[1]);
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 128 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 16;
		unsigned int tileBlockIdx_x = ((bid / 16) + tileBlockIdx_y) % 128;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOut = gbOut + oOffset;


	unsigned int lds_row_padding = 0;
	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 128 * 8 * 128 / 128 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = lwbInRe[me + t*128];
			R0.y = lwbInIm[me + t*128];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = lwbInRe[me % 128 * stride_in[0] + ((me /128 * 16) + t % 8)*stride_in[2] + t / 8 * 128 * stride_in[0]];
			R0.y = lwbInIm[me % 128 * stride_in[0] + ((me /128 * 16) + t % 8)*stride_in[2] + t / 8 * 128 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 128 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *128 + t / 8 * 128 + me % 128 + me / 128 * 1024] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len128_sbrc_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, t * (128 + lds_row_padding) * 8 + (me/16)*(128+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true>(me, 128 - me, 128, 64, &lds[r * (128 + lds_row_padding)], &twiddles[128]);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			lwbOut[(me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*16] = R0;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			lwbOut[(me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*16] = R0;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			lwbOut[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

			lwbOut[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0;
		}
	}
});
});
}

//Kernel configuration: number of threads per thread block: 128, transforms: 4, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose>
__global__ void
__launch_bounds__(128)
fft_fwd_op_len128_sbrc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbInRe = gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbInIm = gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutRe = gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutIm = gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1032), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len128_sbrc>(
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
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
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset *= (stride_out[1]);
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 128 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 16;
		unsigned int tileBlockIdx_x = ((bid / 16) + tileBlockIdx_y) % 128;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;


	unsigned int lds_row_padding = 0;
	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 128 * 8 * 128 / 128 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = lwbInRe[me + t*128];
			R0.y = lwbInIm[me + t*128];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = lwbInRe[me % 128 * stride_in[0] + ((me /128 * 16) + t % 8)*stride_in[2] + t / 8 * 128 * stride_in[0]];
			R0.y = lwbInIm[me % 128 * stride_in[0] + ((me /128 * 16) + t % 8)*stride_in[2] + t / 8 * 128 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 128 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *128 + t / 8 * 128 + me % 128 + me / 128 * 1024] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len128_sbrc_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, t * (128 + lds_row_padding) * 8 + (me/16)*(128+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true>(me, 128 - me, 128, 64, &lds[r * (128 + lds_row_padding)], &twiddles[128]);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			lwbOutRe[(me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*16] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			lwbOutRe[(me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*16] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			lwbOutRe[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

			lwbOutRe[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.y;
		}
	}
});
});
}

//Kernel configuration: number of threads per thread block: 128, transforms: 4, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose>
__global__ void
__launch_bounds__(128)
fft_back_op_len128_sbrc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbInRe = gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbInIm = gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutRe = gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutIm = gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1032), cgh);
	cgh.parallel_for<class kern_fft_back_op_len128_sbrc>(
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
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
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset *= (stride_out[1]);
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 128 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 16;
		unsigned int tileBlockIdx_x = ((bid / 16) + tileBlockIdx_y) % 128;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;


	unsigned int lds_row_padding = 0;
	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1657
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 128 * 8 * 128 / 128 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = lwbInRe[me + t*128];
			R0.y = lwbInIm[me + t*128];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = lwbInRe[me % 128 * stride_in[0] + ((me /128 * 16) + t % 8)*stride_in[2] + t / 8 * 128 * stride_in[0]];
			R0.y = lwbInIm[me % 128 * stride_in[0] + ((me /128 * 16) + t % 8)*stride_in[2] + t / 8 * 128 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 128 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *128 + t / 8 * 128 + me % 128 + me / 128 * 1024] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len128_sbrc_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, t * (128 + lds_row_padding) * 8 + (me/16)*(128+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true>(me, 128 - me, 128, 64, &lds[r * (128 + lds_row_padding)], &twiddles[128]);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1970
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			lwbOutRe[(me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*16] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			lwbOutRe[(me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*16] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			lwbOutRe[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*16 + (me%8)*(128 + lds_row_padding) + (me/8)];

			lwbOutRe[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.x;
			lwbOutIm[(me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*16] = R0.y;
		}
	}
});
});
}

