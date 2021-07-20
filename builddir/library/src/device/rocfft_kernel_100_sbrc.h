#pragma once
#include "rocfft_butterfly_template.h"
#include "real2complex.h"


////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass0_len100_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 10 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 20 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 30 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*1 + 0 + 40 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*1 + 0 + 50 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*1 + 0 + 60 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*1 + 0 + 70 ) ];
	(*R8) = bufIn[inOffset + ( 0 + me*1 + 0 + 80 ) ];
	(*R9) = bufIn[inOffset + ( 0 + me*1 + 0 + 90 ) ];
	}



	FwdRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 1 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 2 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 3 ) ] = (*R3);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 4 ) ] = (*R4);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 5 ) ] = (*R5);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 6 ) ] = (*R6);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 7 ) ] = (*R7);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 8 ) ] = (*R8);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 9 ) ] = (*R9);
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass1_len100_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 10 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 20 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 30 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*1 + 0 + 40 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*1 + 0 + 50 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*1 + 0 + 60 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*1 + 0 + 70 ) ];
	(*R8) = bufIn[inOffset + ( 0 + me*1 + 0 + 80 ) ];
	(*R9) = bufIn[inOffset + ( 0 + me*1 + 0 + 90 ) ];
	}



	{
		T W = twiddles[9 + 9*((1*me + 0)%10) + 0];
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
		T W = twiddles[9 + 9*((1*me + 0)%10) + 1];
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
		T W = twiddles[9 + 9*((1*me + 0)%10) + 2];
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
		T W = twiddles[9 + 9*((1*me + 0)%10) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R4).x; ry = (*R4).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R4).x = TR;
		(*R4).y = TI;
	}

	{
		T W = twiddles[9 + 9*((1*me + 0)%10) + 4];
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
		T W = twiddles[9 + 9*((1*me + 0)%10) + 5];
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
		T W = twiddles[9 + 9*((1*me + 0)%10) + 6];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R7).x; ry = (*R7).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R7).x = TR;
		(*R7).y = TI;
	}

	{
		T W = twiddles[9 + 9*((1*me + 0)%10) + 7];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	{
		T W = twiddles[9 + 9*((1*me + 0)%10) + 8];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R9).x; ry = (*R9).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R9).x = TR;
		(*R9).y = TI;
	}

	FwdRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( 1*me + 0 + 0 ) ] = (*R0);
	bufOut[outOffset + ( 1*me + 0 + 10 ) ] = (*R1);
	bufOut[outOffset + ( 1*me + 0 + 20 ) ] = (*R2);
	bufOut[outOffset + ( 1*me + 0 + 30 ) ] = (*R3);
	bufOut[outOffset + ( 1*me + 0 + 40 ) ] = (*R4);
	bufOut[outOffset + ( 1*me + 0 + 50 ) ] = (*R5);
	bufOut[outOffset + ( 1*me + 0 + 60 ) ] = (*R6);
	bufOut[outOffset + ( 1*me + 0 + 70 ) ] = (*R7);
	bufOut[outOffset + ( 1*me + 0 + 80 ) ] = (*R8);
	bufOut[outOffset + ( 1*me + 0 + 90 ) ] = (*R9);
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass0_len100_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 10 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 20 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 30 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*1 + 0 + 40 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*1 + 0 + 50 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*1 + 0 + 60 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*1 + 0 + 70 ) ];
	(*R8) = bufIn[inOffset + ( 0 + me*1 + 0 + 80 ) ];
	(*R9) = bufIn[inOffset + ( 0 + me*1 + 0 + 90 ) ];
	}



	InvRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 1 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 2 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 3 ) ] = (*R3);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 4 ) ] = (*R4);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 5 ) ] = (*R5);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 6 ) ] = (*R6);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 7 ) ] = (*R7);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 8 ) ] = (*R8);
	bufOut[outOffset + ( ((1*me + 0)/1)*10 + (1*me + 0)%1 + 9 ) ] = (*R9);
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass1_len100_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 10 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 20 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 30 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*1 + 0 + 40 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*1 + 0 + 50 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*1 + 0 + 60 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*1 + 0 + 70 ) ];
	(*R8) = bufIn[inOffset + ( 0 + me*1 + 0 + 80 ) ];
	(*R9) = bufIn[inOffset + ( 0 + me*1 + 0 + 90 ) ];
	}



	{
		T W = twiddles[9 + 9*((1*me + 0)%10) + 0];
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
		T W = twiddles[9 + 9*((1*me + 0)%10) + 1];
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
		T W = twiddles[9 + 9*((1*me + 0)%10) + 2];
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
		T W = twiddles[9 + 9*((1*me + 0)%10) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R4).x; ry = (*R4).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R4).x = TR;
		(*R4).y = TI;
	}

	{
		T W = twiddles[9 + 9*((1*me + 0)%10) + 4];
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
		T W = twiddles[9 + 9*((1*me + 0)%10) + 5];
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
		T W = twiddles[9 + 9*((1*me + 0)%10) + 6];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R7).x; ry = (*R7).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R7).x = TR;
		(*R7).y = TI;
	}

	{
		T W = twiddles[9 + 9*((1*me + 0)%10) + 7];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	{
		T W = twiddles[9 + 9*((1*me + 0)%10) + 8];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R9).x; ry = (*R9).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R9).x = TR;
		(*R9).y = TI;
	}

	InvRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( 1*me + 0 + 0 ) ] = (*R0);
	bufOut[outOffset + ( 1*me + 0 + 10 ) ] = (*R1);
	bufOut[outOffset + ( 1*me + 0 + 20 ) ] = (*R2);
	bufOut[outOffset + ( 1*me + 0 + 30 ) ] = (*R3);
	bufOut[outOffset + ( 1*me + 0 + 40 ) ] = (*R4);
	bufOut[outOffset + ( 1*me + 0 + 50 ) ] = (*R5);
	bufOut[outOffset + ( 1*me + 0 + 60 ) ] = (*R6);
	bufOut[outOffset + ( 1*me + 0 + 70 ) ] = (*R7);
	bufOut[outOffset + ( 1*me + 0 + 80 ) ] = (*R8);
	bufOut[outOffset + ( 1*me + 0 + 90 ) ] = (*R9);
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len100_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, T *gbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9;
	FwdPass0_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  gbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
	FwdPass1_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len100_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, T *gbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9;
	InvPass0_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  gbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
	InvPass1_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len100_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9;
	FwdPass0_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  gbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
	FwdPass1_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len100_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9;
	InvPass0_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  gbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
	InvPass1_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len100_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *gbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9;
	FwdPass0_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
	FwdPass1_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len100_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *gbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9;
	InvPass0_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
	InvPass1_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len100_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9;
	FwdPass0_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
	FwdPass1_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len100_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9;
	InvPass0_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
	InvPass1_len100_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9);
	__syncthreads();
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 50, transforms: 12, Passes: 2
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(50)
fft_fwd_op_len100_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
{

	extern __shared__ T lds[];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 5 - 1) / 5);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 5 - 1) / 5);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/5);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/5));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/5));
		tileOffset_x	= 5*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 5 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 100 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 100;
		iOffset += tileBlockIdx_y * (5 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (5 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 5;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (5 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 5 + t % 5 + me / 100 * 5 * 100 / 50 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = load_cb(gbIn, iOffset + me + t*50, load_cb_data, nullptr);
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = load_cb(gbIn, iOffset + me % 100 * stride_in[0] + ((me /100 * 10) + t % 5)*stride_in[2] + t / 5 * 50 * stride_in[0], load_cb_data, nullptr);
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 50 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 5 *100 + t / 5 * 50 + me % 100 + me / 100 * 1000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		fwd_len100_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (100 + lds_row_padding) * 5 + (me/10)*(100+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 5; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 100 - me, 100, 50, &lds[r * (100 + lds_row_padding)], 0, &twiddles[100], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 5 + me % 5 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			store_cb(gbOut, oOffset + (me%5) * stride_out[1] + (me/5)*stride_out[0] + t*stride_out[0]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			store_cb(gbOut, oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[1] + t*stride_out[1]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			store_cb(gbOut, oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 5)
		{
			unsigned int t = 10;
			T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

			store_cb(gbOut, oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}
}

//Kernel configuration: number of threads per thread block: 50, transforms: 12, Passes: 2
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(50)
fft_back_op_len100_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
{

	extern __shared__ T lds[];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 5 - 1) / 5);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 5 - 1) / 5);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/5);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/5));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/5));
		tileOffset_x	= 5*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 5 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 100 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 100;
		iOffset += tileBlockIdx_y * (5 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (5 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 5;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (5 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 5 + t % 5 + me / 100 * 5 * 100 / 50 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = load_cb(gbIn, iOffset + me + t*50, load_cb_data, nullptr);
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = load_cb(gbIn, iOffset + me % 100 * stride_in[0] + ((me /100 * 10) + t % 5)*stride_in[2] + t / 5 * 50 * stride_in[0], load_cb_data, nullptr);
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 50 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 5 *100 + t / 5 * 50 + me % 100 + me / 100 * 1000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		back_len100_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (100 + lds_row_padding) * 5 + (me/10)*(100+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 5; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 100 - me, 100, 50, &lds[r * (100 + lds_row_padding)], 0, &twiddles[100], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 5 + me % 5 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			store_cb(gbOut, oOffset + (me%5) * stride_out[1] + (me/5)*stride_out[0] + t*stride_out[0]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			store_cb(gbOut, oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[1] + t*stride_out[1]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			store_cb(gbOut, oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 5)
		{
			unsigned int t = 10;
			T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

			store_cb(gbOut, oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}
}

//Kernel configuration: number of threads per thread block: 50, transforms: 12, Passes: 2
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(50)
fft_fwd_op_len100_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	extern __shared__ T lds[];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 5 - 1) / 5);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 5 - 1) / 5);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/5);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/5));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/5));
		tileOffset_x	= 5*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 5 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 100 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 100;
		iOffset += tileBlockIdx_y * (5 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (5 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 5;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (5 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 5 + t % 5 + me / 100 * 5 * 100 / 50 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = load_cb(gbIn, iOffset + me + t*50, load_cb_data, nullptr);
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = load_cb(gbIn, iOffset + me % 100 * stride_in[0] + ((me /100 * 10) + t % 5)*stride_in[2] + t / 5 * 50 * stride_in[0], load_cb_data, nullptr);
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 50 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 5 *100 + t / 5 * 50 + me % 100 + me / 100 * 1000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		fwd_len100_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (100 + lds_row_padding) * 5 + (me/10)*(100+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 5; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 100 - me, 100, 50, &lds[r * (100 + lds_row_padding)], 0, &twiddles[100], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 5 + me % 5 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			gbOutRe[oOffset + (me%5) * stride_out[1] + (me/5)*stride_out[0] + t*stride_out[0]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[1] + (me/5)*stride_out[0] + t*stride_out[0]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			gbOutRe[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[1] + t*stride_out[1]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[1] + t*stride_out[1]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			gbOutRe[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 5)
		{
			unsigned int t = 10;
			T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

			gbOutRe[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}
}

//Kernel configuration: number of threads per thread block: 50, transforms: 12, Passes: 2
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(50)
fft_back_op_len100_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	extern __shared__ T lds[];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 5 - 1) / 5);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 5 - 1) / 5);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/5);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/5));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/5));
		tileOffset_x	= 5*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 5 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 100 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 100;
		iOffset += tileBlockIdx_y * (5 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (5 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 5;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (5 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 5 + t % 5 + me / 100 * 5 * 100 / 50 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = load_cb(gbIn, iOffset + me + t*50, load_cb_data, nullptr);
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = load_cb(gbIn, iOffset + me % 100 * stride_in[0] + ((me /100 * 10) + t % 5)*stride_in[2] + t / 5 * 50 * stride_in[0], load_cb_data, nullptr);
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 50 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 5 *100 + t / 5 * 50 + me % 100 + me / 100 * 1000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		back_len100_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (100 + lds_row_padding) * 5 + (me/10)*(100+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 5; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 100 - me, 100, 50, &lds[r * (100 + lds_row_padding)], 0, &twiddles[100], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 5 + me % 5 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			gbOutRe[oOffset + (me%5) * stride_out[1] + (me/5)*stride_out[0] + t*stride_out[0]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[1] + (me/5)*stride_out[0] + t*stride_out[0]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			gbOutRe[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[1] + t*stride_out[1]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[1] + t*stride_out[1]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			gbOutRe[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 5)
		{
			unsigned int t = 10;
			T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

			gbOutRe[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}
}

//Kernel configuration: number of threads per thread block: 50, transforms: 12, Passes: 2
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(50)
fft_fwd_op_len100_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
{

	extern __shared__ T lds[];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 5 - 1) / 5);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 5 - 1) / 5);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/5);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/5));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/5));
		tileOffset_x	= 5*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 5 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 100 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 100;
		iOffset += tileBlockIdx_y * (5 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (5 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 5;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (5 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 5 + t % 5 + me / 100 * 5 * 100 / 50 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = gbInRe[iOffset + me + t*50];
			R0.y = gbInIm[iOffset + me + t*50];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = gbInRe[iOffset + me % 100 * stride_in[0] + ((me /100 * 10) + t % 5)*stride_in[2] + t / 5 * 50 * stride_in[0]];
			R0.y = gbInIm[iOffset + me % 100 * stride_in[0] + ((me /100 * 10) + t % 5)*stride_in[2] + t / 5 * 50 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 50 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 5 *100 + t / 5 * 50 + me % 100 + me / 100 * 1000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		fwd_len100_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (100 + lds_row_padding) * 5 + (me/10)*(100+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 5; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 100 - me, 100, 50, &lds[r * (100 + lds_row_padding)], 0, &twiddles[100], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 5 + me % 5 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			store_cb(gbOut, oOffset + (me%5) * stride_out[1] + (me/5)*stride_out[0] + t*stride_out[0]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			store_cb(gbOut, oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[1] + t*stride_out[1]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			store_cb(gbOut, oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 5)
		{
			unsigned int t = 10;
			T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

			store_cb(gbOut, oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}
}

//Kernel configuration: number of threads per thread block: 50, transforms: 12, Passes: 2
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(50)
fft_back_op_len100_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
{

	extern __shared__ T lds[];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 5 - 1) / 5);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 5 - 1) / 5);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/5);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/5));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/5));
		tileOffset_x	= 5*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 5 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 100 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 100;
		iOffset += tileBlockIdx_y * (5 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (5 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 5;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (5 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 5 + t % 5 + me / 100 * 5 * 100 / 50 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = gbInRe[iOffset + me + t*50];
			R0.y = gbInIm[iOffset + me + t*50];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = gbInRe[iOffset + me % 100 * stride_in[0] + ((me /100 * 10) + t % 5)*stride_in[2] + t / 5 * 50 * stride_in[0]];
			R0.y = gbInIm[iOffset + me % 100 * stride_in[0] + ((me /100 * 10) + t % 5)*stride_in[2] + t / 5 * 50 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 50 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 5 *100 + t / 5 * 50 + me % 100 + me / 100 * 1000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		back_len100_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (100 + lds_row_padding) * 5 + (me/10)*(100+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 5; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 100 - me, 100, 50, &lds[r * (100 + lds_row_padding)], 0, &twiddles[100], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 5 + me % 5 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			store_cb(gbOut, oOffset + (me%5) * stride_out[1] + (me/5)*stride_out[0] + t*stride_out[0]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			store_cb(gbOut, oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[1] + t*stride_out[1]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			store_cb(gbOut, oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 5)
		{
			unsigned int t = 10;
			T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

			store_cb(gbOut, oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}
}

//Kernel configuration: number of threads per thread block: 50, transforms: 12, Passes: 2
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(50)
fft_fwd_op_len100_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	extern __shared__ T lds[];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 5 - 1) / 5);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 5 - 1) / 5);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/5);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/5));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/5));
		tileOffset_x	= 5*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 5 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 100 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 100;
		iOffset += tileBlockIdx_y * (5 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (5 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 5;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (5 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 5 + t % 5 + me / 100 * 5 * 100 / 50 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = gbInRe[iOffset + me + t*50];
			R0.y = gbInIm[iOffset + me + t*50];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = gbInRe[iOffset + me % 100 * stride_in[0] + ((me /100 * 10) + t % 5)*stride_in[2] + t / 5 * 50 * stride_in[0]];
			R0.y = gbInIm[iOffset + me % 100 * stride_in[0] + ((me /100 * 10) + t % 5)*stride_in[2] + t / 5 * 50 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 50 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 5 *100 + t / 5 * 50 + me % 100 + me / 100 * 1000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		fwd_len100_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (100 + lds_row_padding) * 5 + (me/10)*(100+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 5; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 100 - me, 100, 50, &lds[r * (100 + lds_row_padding)], 0, &twiddles[100], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 5 + me % 5 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			gbOutRe[oOffset + (me%5) * stride_out[1] + (me/5)*stride_out[0] + t*stride_out[0]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[1] + (me/5)*stride_out[0] + t*stride_out[0]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			gbOutRe[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[1] + t*stride_out[1]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[1] + t*stride_out[1]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			gbOutRe[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 5)
		{
			unsigned int t = 10;
			T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

			gbOutRe[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}
}

//Kernel configuration: number of threads per thread block: 50, transforms: 12, Passes: 2
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(50)
fft_back_op_len100_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	extern __shared__ T lds[];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 5 - 1) / 5);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 5 - 1) / 5);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/5);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/5));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/5));
		tileOffset_x	= 5*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 5 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 100 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 100;
		iOffset += tileBlockIdx_y * (5 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (5 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 5;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (5 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 5 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 5 + t % 5 + me / 100 * 5 * 100 / 50 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = gbInRe[iOffset + me + t*50];
			R0.y = gbInIm[iOffset + me + t*50];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = gbInRe[iOffset + me % 100 * stride_in[0] + ((me /100 * 10) + t % 5)*stride_in[2] + t / 5 * 50 * stride_in[0]];
			R0.y = gbInIm[iOffset + me % 100 * stride_in[0] + ((me /100 * 10) + t % 5)*stride_in[2] + t / 5 * 50 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 50 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 5 *100 + t / 5 * 50 + me % 100 + me / 100 * 1000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		back_len100_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (100 + lds_row_padding) * 5 + (me/10)*(100+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 5; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 100 - me, 100, 50, &lds[r * (100 + lds_row_padding)], 0, &twiddles[100], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<10; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 5 + me % 5 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			gbOutRe[oOffset + (me%5) * stride_out[1] + (me/5)*stride_out[0] + t*stride_out[0]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[1] + (me/5)*stride_out[0] + t*stride_out[0]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			gbOutRe[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[1] + t*stride_out[1]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[1] + t*stride_out[1]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			gbOutRe[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 5)
		{
			unsigned int t = 10;
			T R0 = lds[t*10 + (me%5)*(100 + lds_row_padding) + (me/5)];

			gbOutRe[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%5) * stride_out[0] + (me/5)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}
}

