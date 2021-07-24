#pragma once
#include "rocfft_butterfly_template.h"

////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb>
__device__ void
FwdPass0_len169(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 )*stride_in];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 13 )*stride_in];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 26 )*stride_in];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 39 )*stride_in];
	(*R4) = bufIn[inOffset + ( 0 + me*1 + 0 + 52 )*stride_in];
	(*R5) = bufIn[inOffset + ( 0 + me*1 + 0 + 65 )*stride_in];
	(*R6) = bufIn[inOffset + ( 0 + me*1 + 0 + 78 )*stride_in];
	(*R7) = bufIn[inOffset + ( 0 + me*1 + 0 + 91 )*stride_in];
	(*R8) = bufIn[inOffset + ( 0 + me*1 + 0 + 104 )*stride_in];
	(*R9) = bufIn[inOffset + ( 0 + me*1 + 0 + 117 )*stride_in];
	(*R10) = bufIn[inOffset + ( 0 + me*1 + 0 + 130 )*stride_in];
	(*R11) = bufIn[inOffset + ( 0 + me*1 + 0 + 143 )*stride_in];
	(*R12) = bufIn[inOffset + ( 0 + me*1 + 0 + 156 )*stride_in];
	}



	FwdRad13B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12);


	if(rw)
	{
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 10 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 11 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 12 ) ] = (*R12).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 13 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 26 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 39 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 52 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 65 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 78 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 91 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 104 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 117 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 130 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 143 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 156 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 10 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 11 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 12 ) ] = (*R12).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 13 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 26 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 39 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 52 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 65 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 78 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 91 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 104 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 117 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 130 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 143 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 156 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass0_len169(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 13 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 13 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 26 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 26 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 39 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 39 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 52 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 52 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 65 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 65 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 78 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 78 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 91 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 91 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 104 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 104 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 117 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 117 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 130 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 130 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 143 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 143 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 156 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 156 )*stride_in];
	}



	FwdRad13B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12);


	if(rw)
	{
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 10 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 11 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 12 ) ] = (*R12).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 13 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 26 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 39 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 52 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 65 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 78 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 91 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 104 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 117 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 130 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 143 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 156 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 10 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 11 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 12 ) ] = (*R12).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 13 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 26 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 39 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 52 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 65 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 78 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 91 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 104 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 117 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 130 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 143 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 156 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass1_len169(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12)
{




	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 0];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 1];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 2];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 3];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 4];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 5];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 6];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 7];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 8];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R9).x; ry = (*R9).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R9).x = TR;
		(*R9).y = TI;
	}

	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 9];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R10).x; ry = (*R10).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R10).x = TR;
		(*R10).y = TI;
	}

	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 10];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R11).x; ry = (*R11).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R11).x = TR;
		(*R11).y = TI;
	}

	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 11];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R12).x; ry = (*R12).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R12).x = TR;
		(*R12).y = TI;
	}

	FwdRad13B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12);


	if(rw)
	{
	bufOut[outOffset + ( 1*me + 0 + 0 )*stride_out] = (*R0);
	bufOut[outOffset + ( 1*me + 0 + 13 )*stride_out] = (*R1);
	bufOut[outOffset + ( 1*me + 0 + 26 )*stride_out] = (*R2);
	bufOut[outOffset + ( 1*me + 0 + 39 )*stride_out] = (*R3);
	bufOut[outOffset + ( 1*me + 0 + 52 )*stride_out] = (*R4);
	bufOut[outOffset + ( 1*me + 0 + 65 )*stride_out] = (*R5);
	bufOut[outOffset + ( 1*me + 0 + 78 )*stride_out] = (*R6);
	bufOut[outOffset + ( 1*me + 0 + 91 )*stride_out] = (*R7);
	bufOut[outOffset + ( 1*me + 0 + 104 )*stride_out] = (*R8);
	bufOut[outOffset + ( 1*me + 0 + 117 )*stride_out] = (*R9);
	bufOut[outOffset + ( 1*me + 0 + 130 )*stride_out] = (*R10);
	bufOut[outOffset + ( 1*me + 0 + 143 )*stride_out] = (*R11);
	bufOut[outOffset + ( 1*me + 0 + 156 )*stride_out] = (*R12);
	}

}

template <typename T, StrideBin sb>
__device__ void
FwdPass1_len169(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12)
{




	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 0];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 1];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 2];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 3];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 4];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 5];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 6];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 7];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 8];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R9).x; ry = (*R9).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R9).x = TR;
		(*R9).y = TI;
	}

	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 9];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R10).x; ry = (*R10).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R10).x = TR;
		(*R10).y = TI;
	}

	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 10];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R11).x; ry = (*R11).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R11).x = TR;
		(*R11).y = TI;
	}

	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 11];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R12).x; ry = (*R12).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R12).x = TR;
		(*R12).y = TI;
	}

	FwdRad13B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12);


	if(rw)
	{
	bufOutRe[outOffset + ( 1*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 1*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 1*me + 0 + 13 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 1*me + 0 + 13 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 1*me + 0 + 26 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 1*me + 0 + 26 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 1*me + 0 + 39 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 1*me + 0 + 39 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 1*me + 0 + 52 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 1*me + 0 + 52 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 1*me + 0 + 65 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 1*me + 0 + 65 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 1*me + 0 + 78 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 1*me + 0 + 78 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 1*me + 0 + 91 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 1*me + 0 + 91 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 1*me + 0 + 104 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 1*me + 0 + 104 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 1*me + 0 + 117 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 1*me + 0 + 117 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 1*me + 0 + 130 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 1*me + 0 + 130 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 1*me + 0 + 143 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 1*me + 0 + 143 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 1*me + 0 + 156 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 1*me + 0 + 156 )*stride_out] = (*R12).y;
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass0_len169(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 )*stride_in];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 13 )*stride_in];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 26 )*stride_in];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 39 )*stride_in];
	(*R4) = bufIn[inOffset + ( 0 + me*1 + 0 + 52 )*stride_in];
	(*R5) = bufIn[inOffset + ( 0 + me*1 + 0 + 65 )*stride_in];
	(*R6) = bufIn[inOffset + ( 0 + me*1 + 0 + 78 )*stride_in];
	(*R7) = bufIn[inOffset + ( 0 + me*1 + 0 + 91 )*stride_in];
	(*R8) = bufIn[inOffset + ( 0 + me*1 + 0 + 104 )*stride_in];
	(*R9) = bufIn[inOffset + ( 0 + me*1 + 0 + 117 )*stride_in];
	(*R10) = bufIn[inOffset + ( 0 + me*1 + 0 + 130 )*stride_in];
	(*R11) = bufIn[inOffset + ( 0 + me*1 + 0 + 143 )*stride_in];
	(*R12) = bufIn[inOffset + ( 0 + me*1 + 0 + 156 )*stride_in];
	}



	InvRad13B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12);


	if(rw)
	{
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 10 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 11 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 12 ) ] = (*R12).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 13 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 26 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 39 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 52 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 65 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 78 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 91 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 104 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 117 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 130 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 143 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 156 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 10 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 11 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 12 ) ] = (*R12).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 13 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 26 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 39 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 52 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 65 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 78 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 91 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 104 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 117 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 130 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 143 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 156 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass0_len169(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 13 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 13 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 26 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 26 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 39 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 39 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 52 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 52 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 65 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 65 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 78 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 78 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 91 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 91 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 104 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 104 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 117 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 117 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 130 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 130 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 143 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 143 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 156 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 156 )*stride_in];
	}



	InvRad13B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12);


	if(rw)
	{
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 10 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 11 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 12 ) ] = (*R12).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 13 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 26 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 39 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 52 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 65 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 78 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 91 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 104 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 117 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 130 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 143 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 156 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 10 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 11 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*13 + (1*me + 0)%1 + 12 ) ] = (*R12).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 13 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 26 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 39 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 52 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 65 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 78 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 91 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 104 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 117 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 130 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 143 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 156 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass1_len169(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12)
{




	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 0];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 1];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 2];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 3];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 4];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 5];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 6];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 7];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 8];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R9).x; ry = (*R9).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R9).x = TR;
		(*R9).y = TI;
	}

	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 9];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R10).x; ry = (*R10).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R10).x = TR;
		(*R10).y = TI;
	}

	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 10];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R11).x; ry = (*R11).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R11).x = TR;
		(*R11).y = TI;
	}

	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 11];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R12).x; ry = (*R12).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R12).x = TR;
		(*R12).y = TI;
	}

	InvRad13B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12);


	if(rw)
	{
	bufOut[outOffset + ( 1*me + 0 + 0 )*stride_out] = (*R0);
	bufOut[outOffset + ( 1*me + 0 + 13 )*stride_out] = (*R1);
	bufOut[outOffset + ( 1*me + 0 + 26 )*stride_out] = (*R2);
	bufOut[outOffset + ( 1*me + 0 + 39 )*stride_out] = (*R3);
	bufOut[outOffset + ( 1*me + 0 + 52 )*stride_out] = (*R4);
	bufOut[outOffset + ( 1*me + 0 + 65 )*stride_out] = (*R5);
	bufOut[outOffset + ( 1*me + 0 + 78 )*stride_out] = (*R6);
	bufOut[outOffset + ( 1*me + 0 + 91 )*stride_out] = (*R7);
	bufOut[outOffset + ( 1*me + 0 + 104 )*stride_out] = (*R8);
	bufOut[outOffset + ( 1*me + 0 + 117 )*stride_out] = (*R9);
	bufOut[outOffset + ( 1*me + 0 + 130 )*stride_out] = (*R10);
	bufOut[outOffset + ( 1*me + 0 + 143 )*stride_out] = (*R11);
	bufOut[outOffset + ( 1*me + 0 + 156 )*stride_out] = (*R12);
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass1_len169(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12)
{




	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 0];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 1];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 2];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 3];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 4];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 5];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 6];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 7];
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
		T W = twiddles[12 + 12*((1*me + 0)%13) + 8];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R9).x; ry = (*R9).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R9).x = TR;
		(*R9).y = TI;
	}

	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 9];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R10).x; ry = (*R10).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R10).x = TR;
		(*R10).y = TI;
	}

	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 10];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R11).x; ry = (*R11).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R11).x = TR;
		(*R11).y = TI;
	}

	{
		T W = twiddles[12 + 12*((1*me + 0)%13) + 11];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R12).x; ry = (*R12).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R12).x = TR;
		(*R12).y = TI;
	}

	InvRad13B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12);


	if(rw)
	{
	bufOutRe[outOffset + ( 1*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 1*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 1*me + 0 + 13 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 1*me + 0 + 13 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 1*me + 0 + 26 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 1*me + 0 + 26 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 1*me + 0 + 39 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 1*me + 0 + 39 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 1*me + 0 + 52 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 1*me + 0 + 52 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 1*me + 0 + 65 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 1*me + 0 + 65 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 1*me + 0 + 78 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 1*me + 0 + 78 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 1*me + 0 + 91 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 1*me + 0 + 91 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 1*me + 0 + 104 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 1*me + 0 + 104 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 1*me + 0 + 117 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 1*me + 0 + 117 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 1*me + 0 + 130 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 1*me + 0 + 130 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 1*me + 0 + 143 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 1*me + 0 + 143 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 1*me + 0 + 156 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 1*me + 0 + 156 )*stride_out] = (*R12).y;
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb>
__device__ void 
fwd_len169_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12;
	FwdPass0_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
	FwdPass1_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
}

template <typename T, StrideBin sb>
__device__ void 
back_len169_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12;
	InvPass0_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
	InvPass1_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len169_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12;
	FwdPass0_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
	FwdPass1_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
}

template <typename T, StrideBin sb>
__device__ void 
back_len169_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12;
	InvPass0_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
	InvPass1_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len169_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12;
	FwdPass0_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
	FwdPass1_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
}

template <typename T, StrideBin sb>
__device__ void 
back_len169_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12;
	InvPass0_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
	InvPass1_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len169_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12;
	FwdPass0_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
	FwdPass1_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
}

template <typename T, StrideBin sb>
__device__ void 
back_len169_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12;
	InvPass0_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
	InvPass1_len169<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12);
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 169, maximum transforms: 13, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(169)
fft_fwd_ip_len169( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2197), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len169>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int ioOffset = 0;
	T *lwb;

	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*13)*13) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*13 + (me/13));
	if(dim == 1){
		ioOffset += counter_mod*stride_in[1];
	}
	else if(dim == 2){
		int counter_1 = counter_mod / lengths[1];
		int counter_mod_1 = counter_mod % lengths[1];
		ioOffset += counter_1*stride_in[2] + counter_mod_1*stride_in[1];
	}
	else if(dim == 3){
		int counter_2 = counter_mod / (lengths[1] * lengths[2]);
		int counter_mod_2 = counter_mod % (lengths[1] * lengths[2]);
		int counter_1 = counter_mod_2 / lengths[1];
		int counter_mod_1 = counter_mod_2 % lengths[1];
		ioOffset += counter_2*stride_in[3] + counter_1*stride_in[2] + counter_mod_1*stride_in[1];
	}
	else{
		for(int i = dim; i>1; i--){
			int currentLength = 1;
			for(int j=1; j<i; j++){
				currentLength *= lengths[j];
			}

			ioOffset += (counter_mod / currentLength)*stride_in[i];
			counter_mod = counter_mod % currentLength;
		}
		ioOffset+= counter_mod * stride_in[1];
	}
	lwb = gb + ioOffset;

	// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, lwb, lds
	fwd_len169_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%13, (me/13)*169, lwb, lwb, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 169, maximum transforms: 13, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(169)
fft_back_ip_len169( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2197), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len169>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int ioOffset = 0;
	T *lwb;

	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*13)*13) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*13 + (me/13));
	if(dim == 1){
		ioOffset += counter_mod*stride_in[1];
	}
	else if(dim == 2){
		int counter_1 = counter_mod / lengths[1];
		int counter_mod_1 = counter_mod % lengths[1];
		ioOffset += counter_1*stride_in[2] + counter_mod_1*stride_in[1];
	}
	else if(dim == 3){
		int counter_2 = counter_mod / (lengths[1] * lengths[2]);
		int counter_mod_2 = counter_mod % (lengths[1] * lengths[2]);
		int counter_1 = counter_mod_2 / lengths[1];
		int counter_mod_1 = counter_mod_2 % lengths[1];
		ioOffset += counter_2*stride_in[3] + counter_1*stride_in[2] + counter_mod_1*stride_in[1];
	}
	else{
		for(int i = dim; i>1; i--){
			int currentLength = 1;
			for(int j=1; j<i; j++){
				currentLength *= lengths[j];
			}

			ioOffset += (counter_mod / currentLength)*stride_in[i];
			counter_mod = counter_mod % currentLength;
		}
		ioOffset+= counter_mod * stride_in[1];
	}
	lwb = gb + ioOffset;

	// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, lwb, lds
	back_len169_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%13, (me/13)*169, lwb, lwb, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 169, maximum transforms: 13, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(169)
fft_fwd_ip_len169( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2197), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len169>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int ioOffset = 0;
	real_type_t<T> *lwbRe;
	real_type_t<T> *lwbIm;

	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*13)*13) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*13 + (me/13));
	if(dim == 1){
		ioOffset += counter_mod*stride_in[1];
	}
	else if(dim == 2){
		int counter_1 = counter_mod / lengths[1];
		int counter_mod_1 = counter_mod % lengths[1];
		ioOffset += counter_1*stride_in[2] + counter_mod_1*stride_in[1];
	}
	else if(dim == 3){
		int counter_2 = counter_mod / (lengths[1] * lengths[2]);
		int counter_mod_2 = counter_mod % (lengths[1] * lengths[2]);
		int counter_1 = counter_mod_2 / lengths[1];
		int counter_mod_1 = counter_mod_2 % lengths[1];
		ioOffset += counter_2*stride_in[3] + counter_1*stride_in[2] + counter_mod_1*stride_in[1];
	}
	else{
		for(int i = dim; i>1; i--){
			int currentLength = 1;
			for(int j=1; j<i; j++){
				currentLength *= lengths[j];
			}

			ioOffset += (counter_mod / currentLength)*stride_in[i];
			counter_mod = counter_mod % currentLength;
		}
		ioOffset+= counter_mod * stride_in[1];
	}
	lwbRe = gbRe + ioOffset;
	lwbIm = gbIm + ioOffset;

	// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, lwb, lds
	fwd_len169_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%13, (me/13)*169, lwbRe, lwbIm, lwbRe, lwbIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 169, maximum transforms: 13, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(169)
fft_back_ip_len169( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2197), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len169>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int ioOffset = 0;
	real_type_t<T> *lwbRe;
	real_type_t<T> *lwbIm;

	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*13)*13) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*13 + (me/13));
	if(dim == 1){
		ioOffset += counter_mod*stride_in[1];
	}
	else if(dim == 2){
		int counter_1 = counter_mod / lengths[1];
		int counter_mod_1 = counter_mod % lengths[1];
		ioOffset += counter_1*stride_in[2] + counter_mod_1*stride_in[1];
	}
	else if(dim == 3){
		int counter_2 = counter_mod / (lengths[1] * lengths[2]);
		int counter_mod_2 = counter_mod % (lengths[1] * lengths[2]);
		int counter_1 = counter_mod_2 / lengths[1];
		int counter_mod_1 = counter_mod_2 % lengths[1];
		ioOffset += counter_2*stride_in[3] + counter_1*stride_in[2] + counter_mod_1*stride_in[1];
	}
	else{
		for(int i = dim; i>1; i--){
			int currentLength = 1;
			for(int j=1; j<i; j++){
				currentLength *= lengths[j];
			}

			ioOffset += (counter_mod / currentLength)*stride_in[i];
			counter_mod = counter_mod % currentLength;
		}
		ioOffset+= counter_mod * stride_in[1];
	}
	lwbRe = gbRe + ioOffset;
	lwbIm = gbIm + ioOffset;

	// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, lwb, lds
	back_len169_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%13, (me/13)*169, lwbRe, lwbIm, lwbRe, lwbIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 169, maximum transforms: 13, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(169)
fft_fwd_op_len169( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2197), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len169>(
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

	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*13)*13) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*13 + (me/13));
	if(dim == 1){
		iOffset += counter_mod*stride_in[1];
		oOffset += counter_mod*stride_out[1];
	}
	else if(dim == 2){
		int counter_1 = counter_mod / lengths[1];
		int counter_mod_1 = counter_mod % lengths[1];
		iOffset += counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else if(dim == 3){
		int counter_2 = counter_mod / (lengths[1] * lengths[2]);
		int counter_mod_2 = counter_mod % (lengths[1] * lengths[2]);
		int counter_1 = counter_mod_2 / lengths[1];
		int counter_mod_1 = counter_mod_2 % lengths[1];
		iOffset += counter_2*stride_in[3] + counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_2*stride_out[3] + counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else{
		for(int i = dim; i>1; i--){
			int currentLength = 1;
			for(int j=1; j<i; j++){
				currentLength *= lengths[j];
			}

			iOffset += (counter_mod / currentLength)*stride_in[i];
			oOffset += (counter_mod / currentLength)*stride_out[i];
			counter_mod = counter_mod % currentLength;
		}
		iOffset+= counter_mod * stride_in[1];
		oOffset+= counter_mod * stride_out[1];
	}
	lwbIn = gbIn + iOffset;
	lwbOut = gbOut + oOffset;

	// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, lwb, lds
	fwd_len169_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%13, (me/13)*169, lwbIn, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 169, maximum transforms: 13, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(169)
fft_back_op_len169( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2197), cgh);
	cgh.parallel_for<class kern_fft_back_op_len169>(
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

	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*13)*13) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*13 + (me/13));
	if(dim == 1){
		iOffset += counter_mod*stride_in[1];
		oOffset += counter_mod*stride_out[1];
	}
	else if(dim == 2){
		int counter_1 = counter_mod / lengths[1];
		int counter_mod_1 = counter_mod % lengths[1];
		iOffset += counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else if(dim == 3){
		int counter_2 = counter_mod / (lengths[1] * lengths[2]);
		int counter_mod_2 = counter_mod % (lengths[1] * lengths[2]);
		int counter_1 = counter_mod_2 / lengths[1];
		int counter_mod_1 = counter_mod_2 % lengths[1];
		iOffset += counter_2*stride_in[3] + counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_2*stride_out[3] + counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else{
		for(int i = dim; i>1; i--){
			int currentLength = 1;
			for(int j=1; j<i; j++){
				currentLength *= lengths[j];
			}

			iOffset += (counter_mod / currentLength)*stride_in[i];
			oOffset += (counter_mod / currentLength)*stride_out[i];
			counter_mod = counter_mod % currentLength;
		}
		iOffset+= counter_mod * stride_in[1];
		oOffset+= counter_mod * stride_out[1];
	}
	lwbIn = gbIn + iOffset;
	lwbOut = gbOut + oOffset;

	// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, lwb, lds
	back_len169_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%13, (me/13)*169, lwbIn, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 169, maximum transforms: 13, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(169)
fft_fwd_op_len169( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2197), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len169>(
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

	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*13)*13) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*13 + (me/13));
	if(dim == 1){
		iOffset += counter_mod*stride_in[1];
		oOffset += counter_mod*stride_out[1];
	}
	else if(dim == 2){
		int counter_1 = counter_mod / lengths[1];
		int counter_mod_1 = counter_mod % lengths[1];
		iOffset += counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else if(dim == 3){
		int counter_2 = counter_mod / (lengths[1] * lengths[2]);
		int counter_mod_2 = counter_mod % (lengths[1] * lengths[2]);
		int counter_1 = counter_mod_2 / lengths[1];
		int counter_mod_1 = counter_mod_2 % lengths[1];
		iOffset += counter_2*stride_in[3] + counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_2*stride_out[3] + counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else{
		for(int i = dim; i>1; i--){
			int currentLength = 1;
			for(int j=1; j<i; j++){
				currentLength *= lengths[j];
			}

			iOffset += (counter_mod / currentLength)*stride_in[i];
			oOffset += (counter_mod / currentLength)*stride_out[i];
			counter_mod = counter_mod % currentLength;
		}
		iOffset+= counter_mod * stride_in[1];
		oOffset+= counter_mod * stride_out[1];
	}
	lwbIn = gbIn + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;

	// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, lwb, lds
	fwd_len169_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%13, (me/13)*169, lwbIn, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 169, maximum transforms: 13, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(169)
fft_back_op_len169( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2197), cgh);
	cgh.parallel_for<class kern_fft_back_op_len169>(
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

	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*13)*13) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*13 + (me/13));
	if(dim == 1){
		iOffset += counter_mod*stride_in[1];
		oOffset += counter_mod*stride_out[1];
	}
	else if(dim == 2){
		int counter_1 = counter_mod / lengths[1];
		int counter_mod_1 = counter_mod % lengths[1];
		iOffset += counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else if(dim == 3){
		int counter_2 = counter_mod / (lengths[1] * lengths[2]);
		int counter_mod_2 = counter_mod % (lengths[1] * lengths[2]);
		int counter_1 = counter_mod_2 / lengths[1];
		int counter_mod_1 = counter_mod_2 % lengths[1];
		iOffset += counter_2*stride_in[3] + counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_2*stride_out[3] + counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else{
		for(int i = dim; i>1; i--){
			int currentLength = 1;
			for(int j=1; j<i; j++){
				currentLength *= lengths[j];
			}

			iOffset += (counter_mod / currentLength)*stride_in[i];
			oOffset += (counter_mod / currentLength)*stride_out[i];
			counter_mod = counter_mod % currentLength;
		}
		iOffset+= counter_mod * stride_in[1];
		oOffset+= counter_mod * stride_out[1];
	}
	lwbIn = gbIn + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;

	// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, lwb, lds
	back_len169_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%13, (me/13)*169, lwbIn, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 169, maximum transforms: 13, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(169)
fft_fwd_op_len169( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2197), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len169>(
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

	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*13)*13) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*13 + (me/13));
	if(dim == 1){
		iOffset += counter_mod*stride_in[1];
		oOffset += counter_mod*stride_out[1];
	}
	else if(dim == 2){
		int counter_1 = counter_mod / lengths[1];
		int counter_mod_1 = counter_mod % lengths[1];
		iOffset += counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else if(dim == 3){
		int counter_2 = counter_mod / (lengths[1] * lengths[2]);
		int counter_mod_2 = counter_mod % (lengths[1] * lengths[2]);
		int counter_1 = counter_mod_2 / lengths[1];
		int counter_mod_1 = counter_mod_2 % lengths[1];
		iOffset += counter_2*stride_in[3] + counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_2*stride_out[3] + counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else{
		for(int i = dim; i>1; i--){
			int currentLength = 1;
			for(int j=1; j<i; j++){
				currentLength *= lengths[j];
			}

			iOffset += (counter_mod / currentLength)*stride_in[i];
			oOffset += (counter_mod / currentLength)*stride_out[i];
			counter_mod = counter_mod % currentLength;
		}
		iOffset+= counter_mod * stride_in[1];
		oOffset+= counter_mod * stride_out[1];
	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOut = gbOut + oOffset;

	// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, lwb, lds
	fwd_len169_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%13, (me/13)*169, lwbInRe, lwbInIm, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 169, maximum transforms: 13, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(169)
fft_back_op_len169( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2197), cgh);
	cgh.parallel_for<class kern_fft_back_op_len169>(
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

	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*13)*13) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*13 + (me/13));
	if(dim == 1){
		iOffset += counter_mod*stride_in[1];
		oOffset += counter_mod*stride_out[1];
	}
	else if(dim == 2){
		int counter_1 = counter_mod / lengths[1];
		int counter_mod_1 = counter_mod % lengths[1];
		iOffset += counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else if(dim == 3){
		int counter_2 = counter_mod / (lengths[1] * lengths[2]);
		int counter_mod_2 = counter_mod % (lengths[1] * lengths[2]);
		int counter_1 = counter_mod_2 / lengths[1];
		int counter_mod_1 = counter_mod_2 % lengths[1];
		iOffset += counter_2*stride_in[3] + counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_2*stride_out[3] + counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else{
		for(int i = dim; i>1; i--){
			int currentLength = 1;
			for(int j=1; j<i; j++){
				currentLength *= lengths[j];
			}

			iOffset += (counter_mod / currentLength)*stride_in[i];
			oOffset += (counter_mod / currentLength)*stride_out[i];
			counter_mod = counter_mod % currentLength;
		}
		iOffset+= counter_mod * stride_in[1];
		oOffset+= counter_mod * stride_out[1];
	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOut = gbOut + oOffset;

	// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, lwb, lds
	back_len169_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%13, (me/13)*169, lwbInRe, lwbInIm, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 169, maximum transforms: 13, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(169)
fft_fwd_op_len169( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2197), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len169>(
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

	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*13)*13) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*13 + (me/13));
	if(dim == 1){
		iOffset += counter_mod*stride_in[1];
		oOffset += counter_mod*stride_out[1];
	}
	else if(dim == 2){
		int counter_1 = counter_mod / lengths[1];
		int counter_mod_1 = counter_mod % lengths[1];
		iOffset += counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else if(dim == 3){
		int counter_2 = counter_mod / (lengths[1] * lengths[2]);
		int counter_mod_2 = counter_mod % (lengths[1] * lengths[2]);
		int counter_1 = counter_mod_2 / lengths[1];
		int counter_mod_1 = counter_mod_2 % lengths[1];
		iOffset += counter_2*stride_in[3] + counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_2*stride_out[3] + counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else{
		for(int i = dim; i>1; i--){
			int currentLength = 1;
			for(int j=1; j<i; j++){
				currentLength *= lengths[j];
			}

			iOffset += (counter_mod / currentLength)*stride_in[i];
			oOffset += (counter_mod / currentLength)*stride_out[i];
			counter_mod = counter_mod % currentLength;
		}
		iOffset+= counter_mod * stride_in[1];
		oOffset+= counter_mod * stride_out[1];
	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;

	// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, lwb, lds
	fwd_len169_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%13, (me/13)*169, lwbInRe, lwbInIm, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 169, maximum transforms: 13, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(169)
fft_back_op_len169( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2197), cgh);
	cgh.parallel_for<class kern_fft_back_op_len169>(
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

	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*13)*13) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*13 + (me/13));
	if(dim == 1){
		iOffset += counter_mod*stride_in[1];
		oOffset += counter_mod*stride_out[1];
	}
	else if(dim == 2){
		int counter_1 = counter_mod / lengths[1];
		int counter_mod_1 = counter_mod % lengths[1];
		iOffset += counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else if(dim == 3){
		int counter_2 = counter_mod / (lengths[1] * lengths[2]);
		int counter_mod_2 = counter_mod % (lengths[1] * lengths[2]);
		int counter_1 = counter_mod_2 / lengths[1];
		int counter_mod_1 = counter_mod_2 % lengths[1];
		iOffset += counter_2*stride_in[3] + counter_1*stride_in[2] + counter_mod_1*stride_in[1];
		oOffset += counter_2*stride_out[3] + counter_1*stride_out[2] + counter_mod_1*stride_out[1];
	}
	else{
		for(int i = dim; i>1; i--){
			int currentLength = 1;
			for(int j=1; j<i; j++){
				currentLength *= lengths[j];
			}

			iOffset += (counter_mod / currentLength)*stride_in[i];
			oOffset += (counter_mod / currentLength)*stride_out[i];
			counter_mod = counter_mod % currentLength;
		}
		iOffset+= counter_mod * stride_in[1];
		oOffset+= counter_mod * stride_out[1];
	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;

	// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, lwb, lds
	back_len169_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%13, (me/13)*169, lwbInRe, lwbInIm, lwbOutRe, lwbOutIm, lds);
});
});
}

