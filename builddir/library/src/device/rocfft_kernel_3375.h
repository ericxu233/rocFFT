#pragma once
#include "rocfft_butterfly_template.h"

////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb>
__device__ void
FwdPass0_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R5) = bufIn[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R10) = bufIn[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R1) = bufIn[inOffset + ( 0 + me*3 + 0 + 675 )*stride_in];
	(*R6) = bufIn[inOffset + ( 0 + me*3 + 1 + 675 )*stride_in];
	(*R11) = bufIn[inOffset + ( 0 + me*3 + 2 + 675 )*stride_in];
	(*R2) = bufIn[inOffset + ( 0 + me*3 + 0 + 1350 )*stride_in];
	(*R7) = bufIn[inOffset + ( 0 + me*3 + 1 + 1350 )*stride_in];
	(*R12) = bufIn[inOffset + ( 0 + me*3 + 2 + 1350 )*stride_in];
	(*R3) = bufIn[inOffset + ( 0 + me*3 + 0 + 2025 )*stride_in];
	(*R8) = bufIn[inOffset + ( 0 + me*3 + 1 + 2025 )*stride_in];
	(*R13) = bufIn[inOffset + ( 0 + me*3 + 2 + 2025 )*stride_in];
	(*R4) = bufIn[inOffset + ( 0 + me*3 + 0 + 2700 )*stride_in];
	(*R9) = bufIn[inOffset + ( 0 + me*3 + 1 + 2700 )*stride_in];
	(*R14) = bufIn[inOffset + ( 0 + me*3 + 2 + 2700 )*stride_in];
	}



	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 1 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 2 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 3 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 4 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 4 ) ] = (*R14).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 675 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 675 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 675 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1350 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1350 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1350 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2025 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2025 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2025 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 1 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 2 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 3 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 4 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 4 ) ] = (*R14).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 675 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 675 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 675 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1350 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1350 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1350 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2025 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2025 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2025 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass0_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 675 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 675 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 675 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 675 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 675 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 675 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1350 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1350 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1350 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1350 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1350 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1350 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2025 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2025 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2025 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2025 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2025 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2025 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2700 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2700 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2700 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2700 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2700 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2700 )*stride_in];
	}



	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 1 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 2 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 3 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 4 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 4 ) ] = (*R14).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 675 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 675 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 675 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1350 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1350 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1350 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2025 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2025 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2025 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 1 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 2 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 3 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 4 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 4 ) ] = (*R14).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 675 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 675 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 675 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1350 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1350 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1350 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2025 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2025 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2025 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass1_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{




	{
		T W = twiddles[4 + 4*((3*me + 0)%5) + 0];
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
		T W = twiddles[4 + 4*((3*me + 0)%5) + 1];
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
		T W = twiddles[4 + 4*((3*me + 0)%5) + 2];
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
		T W = twiddles[4 + 4*((3*me + 0)%5) + 3];
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
		T W = twiddles[4 + 4*((3*me + 1)%5) + 0];
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
		T W = twiddles[4 + 4*((3*me + 1)%5) + 1];
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
		T W = twiddles[4 + 4*((3*me + 1)%5) + 2];
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
		T W = twiddles[4 + 4*((3*me + 1)%5) + 3];
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
		T W = twiddles[4 + 4*((3*me + 2)%5) + 0];
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
		T W = twiddles[4 + 4*((3*me + 2)%5) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R12).x; ry = (*R12).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R12).x = TR;
		(*R12).y = TI;
	}

	{
		T W = twiddles[4 + 4*((3*me + 2)%5) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		T W = twiddles[4 + 4*((3*me + 2)%5) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 5 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 10 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 15 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 20 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 5 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 10 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 15 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 20 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 5 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 10 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 15 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 20 ) ] = (*R14).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 675 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 675 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 675 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1350 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1350 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1350 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2025 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2025 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2025 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 5 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 10 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 15 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 20 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 5 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 10 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 15 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 20 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 5 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 10 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 15 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 20 ) ] = (*R14).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 675 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 675 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 675 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1350 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1350 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1350 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2025 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2025 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2025 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass2_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{




	{
		T W = twiddles[24 + 4*((3*me + 0)%25) + 0];
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
		T W = twiddles[24 + 4*((3*me + 0)%25) + 1];
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
		T W = twiddles[24 + 4*((3*me + 0)%25) + 2];
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
		T W = twiddles[24 + 4*((3*me + 0)%25) + 3];
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
		T W = twiddles[24 + 4*((3*me + 1)%25) + 0];
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
		T W = twiddles[24 + 4*((3*me + 1)%25) + 1];
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
		T W = twiddles[24 + 4*((3*me + 1)%25) + 2];
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
		T W = twiddles[24 + 4*((3*me + 1)%25) + 3];
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
		T W = twiddles[24 + 4*((3*me + 2)%25) + 0];
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
		T W = twiddles[24 + 4*((3*me + 2)%25) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R12).x; ry = (*R12).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R12).x = TR;
		(*R12).y = TI;
	}

	{
		T W = twiddles[24 + 4*((3*me + 2)%25) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		T W = twiddles[24 + 4*((3*me + 2)%25) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 25 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 50 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 75 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 100 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 25 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 50 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 75 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 100 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 25 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 50 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 75 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 100 ) ] = (*R14).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1125 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1125 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1125 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1125 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1125 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2250 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2250 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2250 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2250 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2250 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 25 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 50 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 75 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 100 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 25 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 50 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 75 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 100 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 25 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 50 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 75 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 100 ) ] = (*R14).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1125 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1125 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1125 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1125 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1125 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2250 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2250 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2250 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2250 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2250 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass3_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{




	{
		T W = twiddles[124 + 2*((5*me + 0)%125) + 0];
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
		T W = twiddles[124 + 2*((5*me + 0)%125) + 1];
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
		T W = twiddles[124 + 2*((5*me + 1)%125) + 0];
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
		T W = twiddles[124 + 2*((5*me + 1)%125) + 1];
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
		T W = twiddles[124 + 2*((5*me + 2)%125) + 0];
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
		T W = twiddles[124 + 2*((5*me + 2)%125) + 1];
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
		T W = twiddles[124 + 2*((5*me + 3)%125) + 0];
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
		T W = twiddles[124 + 2*((5*me + 3)%125) + 1];
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
		T W = twiddles[124 + 2*((5*me + 4)%125) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		T W = twiddles[124 + 2*((5*me + 4)%125) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);
	FwdRad3B1(R9, R10, R11);
	FwdRad3B1(R12, R13, R14);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/125)*375 + (5*me + 0)%125 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*375 + (5*me + 0)%125 + 125 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*375 + (5*me + 0)%125 + 250 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*375 + (5*me + 1)%125 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*375 + (5*me + 1)%125 + 125 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*375 + (5*me + 1)%125 + 250 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*375 + (5*me + 2)%125 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*375 + (5*me + 2)%125 + 125 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*375 + (5*me + 2)%125 + 250 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*375 + (5*me + 3)%125 + 0 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*375 + (5*me + 3)%125 + 125 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*375 + (5*me + 3)%125 + 250 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*375 + (5*me + 4)%125 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*375 + (5*me + 4)%125 + 125 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*375 + (5*me + 4)%125 + 250 ) ] = (*R14).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1125 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1125 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1125 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1125 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1125 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2250 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2250 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2250 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2250 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2250 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/125)*375 + (5*me + 0)%125 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*375 + (5*me + 0)%125 + 125 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*375 + (5*me + 0)%125 + 250 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*375 + (5*me + 1)%125 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*375 + (5*me + 1)%125 + 125 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*375 + (5*me + 1)%125 + 250 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*375 + (5*me + 2)%125 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*375 + (5*me + 2)%125 + 125 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*375 + (5*me + 2)%125 + 250 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*375 + (5*me + 3)%125 + 0 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*375 + (5*me + 3)%125 + 125 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*375 + (5*me + 3)%125 + 250 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*375 + (5*me + 4)%125 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*375 + (5*me + 4)%125 + 125 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*375 + (5*me + 4)%125 + 250 ) ] = (*R14).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1125 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1125 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1125 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1125 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1125 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2250 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2250 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2250 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2250 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2250 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass4_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{




	{
		T W = twiddles[374 + 2*((5*me + 0)%375) + 0];
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
		T W = twiddles[374 + 2*((5*me + 0)%375) + 1];
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
		T W = twiddles[374 + 2*((5*me + 1)%375) + 0];
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
		T W = twiddles[374 + 2*((5*me + 1)%375) + 1];
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
		T W = twiddles[374 + 2*((5*me + 2)%375) + 0];
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
		T W = twiddles[374 + 2*((5*me + 2)%375) + 1];
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
		T W = twiddles[374 + 2*((5*me + 3)%375) + 0];
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
		T W = twiddles[374 + 2*((5*me + 3)%375) + 1];
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
		T W = twiddles[374 + 2*((5*me + 4)%375) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		T W = twiddles[374 + 2*((5*me + 4)%375) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);
	FwdRad3B1(R9, R10, R11);
	FwdRad3B1(R12, R13, R14);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/375)*1125 + (5*me + 0)%375 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/375)*1125 + (5*me + 0)%375 + 375 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/375)*1125 + (5*me + 0)%375 + 750 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 1)/375)*1125 + (5*me + 1)%375 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 1)/375)*1125 + (5*me + 1)%375 + 375 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/375)*1125 + (5*me + 1)%375 + 750 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 2)/375)*1125 + (5*me + 2)%375 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 2)/375)*1125 + (5*me + 2)%375 + 375 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 2)/375)*1125 + (5*me + 2)%375 + 750 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 3)/375)*1125 + (5*me + 3)%375 + 0 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 3)/375)*1125 + (5*me + 3)%375 + 375 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 3)/375)*1125 + (5*me + 3)%375 + 750 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 4)/375)*1125 + (5*me + 4)%375 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 4)/375)*1125 + (5*me + 4)%375 + 375 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 4)/375)*1125 + (5*me + 4)%375 + 750 ) ] = (*R14).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1125 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1125 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1125 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1125 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1125 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2250 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2250 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2250 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2250 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2250 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/375)*1125 + (5*me + 0)%375 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/375)*1125 + (5*me + 0)%375 + 375 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/375)*1125 + (5*me + 0)%375 + 750 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 1)/375)*1125 + (5*me + 1)%375 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 1)/375)*1125 + (5*me + 1)%375 + 375 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/375)*1125 + (5*me + 1)%375 + 750 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 2)/375)*1125 + (5*me + 2)%375 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 2)/375)*1125 + (5*me + 2)%375 + 375 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 2)/375)*1125 + (5*me + 2)%375 + 750 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 3)/375)*1125 + (5*me + 3)%375 + 0 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 3)/375)*1125 + (5*me + 3)%375 + 375 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 3)/375)*1125 + (5*me + 3)%375 + 750 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 4)/375)*1125 + (5*me + 4)%375 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 4)/375)*1125 + (5*me + 4)%375 + 375 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 4)/375)*1125 + (5*me + 4)%375 + 750 ) ] = (*R14).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1125 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1125 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1125 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1125 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1125 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2250 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2250 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2250 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2250 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2250 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass5_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{




	{
		T W = twiddles[1124 + 2*((5*me + 0)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 0)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 1)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 1)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 2)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 2)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 3)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 3)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 4)%1125) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		T W = twiddles[1124 + 2*((5*me + 4)%1125) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);
	FwdRad3B1(R9, R10, R11);
	FwdRad3B1(R12, R13, R14);


	if(rw)
	{
	bufOut[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0);
	bufOut[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R3);
	bufOut[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R6);
	bufOut[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R9);
	bufOut[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R12);
	bufOut[outOffset + ( 5*me + 0 + 1125 )*stride_out] = (*R1);
	bufOut[outOffset + ( 5*me + 1 + 1125 )*stride_out] = (*R4);
	bufOut[outOffset + ( 5*me + 2 + 1125 )*stride_out] = (*R7);
	bufOut[outOffset + ( 5*me + 3 + 1125 )*stride_out] = (*R10);
	bufOut[outOffset + ( 5*me + 4 + 1125 )*stride_out] = (*R13);
	bufOut[outOffset + ( 5*me + 0 + 2250 )*stride_out] = (*R2);
	bufOut[outOffset + ( 5*me + 1 + 2250 )*stride_out] = (*R5);
	bufOut[outOffset + ( 5*me + 2 + 2250 )*stride_out] = (*R8);
	bufOut[outOffset + ( 5*me + 3 + 2250 )*stride_out] = (*R11);
	bufOut[outOffset + ( 5*me + 4 + 2250 )*stride_out] = (*R14);
	}

}

template <typename T, StrideBin sb>
__device__ void
FwdPass5_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{




	{
		T W = twiddles[1124 + 2*((5*me + 0)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 0)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 1)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 1)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 2)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 2)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 3)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 3)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 4)%1125) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		T W = twiddles[1124 + 2*((5*me + 4)%1125) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);
	FwdRad3B1(R9, R10, R11);
	FwdRad3B1(R12, R13, R14);


	if(rw)
	{
	bufOutRe[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R12).y;
	bufOutRe[outOffset + ( 5*me + 0 + 1125 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 5*me + 0 + 1125 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 5*me + 1 + 1125 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 5*me + 1 + 1125 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 5*me + 2 + 1125 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 5*me + 2 + 1125 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 5*me + 3 + 1125 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 5*me + 3 + 1125 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 5*me + 4 + 1125 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 5*me + 4 + 1125 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 5*me + 0 + 2250 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 5*me + 0 + 2250 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 5*me + 1 + 2250 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 5*me + 1 + 2250 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 5*me + 2 + 2250 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 5*me + 2 + 2250 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 5*me + 3 + 2250 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 5*me + 3 + 2250 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 5*me + 4 + 2250 )*stride_out] = (*R14).x;
	bufOutIm[outOffset + ( 5*me + 4 + 2250 )*stride_out] = (*R14).y;
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass0_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R5) = bufIn[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R10) = bufIn[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R1) = bufIn[inOffset + ( 0 + me*3 + 0 + 675 )*stride_in];
	(*R6) = bufIn[inOffset + ( 0 + me*3 + 1 + 675 )*stride_in];
	(*R11) = bufIn[inOffset + ( 0 + me*3 + 2 + 675 )*stride_in];
	(*R2) = bufIn[inOffset + ( 0 + me*3 + 0 + 1350 )*stride_in];
	(*R7) = bufIn[inOffset + ( 0 + me*3 + 1 + 1350 )*stride_in];
	(*R12) = bufIn[inOffset + ( 0 + me*3 + 2 + 1350 )*stride_in];
	(*R3) = bufIn[inOffset + ( 0 + me*3 + 0 + 2025 )*stride_in];
	(*R8) = bufIn[inOffset + ( 0 + me*3 + 1 + 2025 )*stride_in];
	(*R13) = bufIn[inOffset + ( 0 + me*3 + 2 + 2025 )*stride_in];
	(*R4) = bufIn[inOffset + ( 0 + me*3 + 0 + 2700 )*stride_in];
	(*R9) = bufIn[inOffset + ( 0 + me*3 + 1 + 2700 )*stride_in];
	(*R14) = bufIn[inOffset + ( 0 + me*3 + 2 + 2700 )*stride_in];
	}



	InvRad5B1(R0, R1, R2, R3, R4);
	InvRad5B1(R5, R6, R7, R8, R9);
	InvRad5B1(R10, R11, R12, R13, R14);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 1 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 2 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 3 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 4 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 4 ) ] = (*R14).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 675 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 675 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 675 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1350 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1350 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1350 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2025 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2025 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2025 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 1 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 2 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 3 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 4 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 4 ) ] = (*R14).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 675 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 675 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 675 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1350 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1350 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1350 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2025 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2025 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2025 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass0_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 675 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 675 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 675 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 675 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 675 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 675 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1350 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1350 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1350 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1350 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1350 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1350 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2025 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2025 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2025 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2025 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2025 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2025 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2700 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2700 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2700 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2700 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2700 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2700 )*stride_in];
	}



	InvRad5B1(R0, R1, R2, R3, R4);
	InvRad5B1(R5, R6, R7, R8, R9);
	InvRad5B1(R10, R11, R12, R13, R14);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 1 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 2 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 3 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 4 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 4 ) ] = (*R14).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 675 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 675 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 675 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1350 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1350 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1350 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2025 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2025 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2025 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*5 + (3*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 1 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 2 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 3 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*5 + (3*me + 1)%1 + 4 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*5 + (3*me + 2)%1 + 4 ) ] = (*R14).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 675 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 675 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 675 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1350 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1350 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1350 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2025 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2025 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2025 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass1_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{




	{
		T W = twiddles[4 + 4*((3*me + 0)%5) + 0];
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
		T W = twiddles[4 + 4*((3*me + 0)%5) + 1];
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
		T W = twiddles[4 + 4*((3*me + 0)%5) + 2];
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
		T W = twiddles[4 + 4*((3*me + 0)%5) + 3];
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
		T W = twiddles[4 + 4*((3*me + 1)%5) + 0];
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
		T W = twiddles[4 + 4*((3*me + 1)%5) + 1];
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
		T W = twiddles[4 + 4*((3*me + 1)%5) + 2];
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
		T W = twiddles[4 + 4*((3*me + 1)%5) + 3];
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
		T W = twiddles[4 + 4*((3*me + 2)%5) + 0];
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
		T W = twiddles[4 + 4*((3*me + 2)%5) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R12).x; ry = (*R12).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R12).x = TR;
		(*R12).y = TI;
	}

	{
		T W = twiddles[4 + 4*((3*me + 2)%5) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		T W = twiddles[4 + 4*((3*me + 2)%5) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	InvRad5B1(R0, R1, R2, R3, R4);
	InvRad5B1(R5, R6, R7, R8, R9);
	InvRad5B1(R10, R11, R12, R13, R14);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 5 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 10 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 15 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 20 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 5 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 10 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 15 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 20 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 5 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 10 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 15 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 20 ) ] = (*R14).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 675 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 675 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 675 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1350 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1350 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1350 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2025 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2025 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2025 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 5 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 10 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 15 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/5)*25 + (3*me + 0)%5 + 20 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 5 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 10 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 15 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 1)/5)*25 + (3*me + 1)%5 + 20 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 5 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 10 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 15 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 2)/5)*25 + (3*me + 2)%5 + 20 ) ] = (*R14).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 675 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 675 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 675 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1350 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1350 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1350 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2025 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2025 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2025 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass2_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{




	{
		T W = twiddles[24 + 4*((3*me + 0)%25) + 0];
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
		T W = twiddles[24 + 4*((3*me + 0)%25) + 1];
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
		T W = twiddles[24 + 4*((3*me + 0)%25) + 2];
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
		T W = twiddles[24 + 4*((3*me + 0)%25) + 3];
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
		T W = twiddles[24 + 4*((3*me + 1)%25) + 0];
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
		T W = twiddles[24 + 4*((3*me + 1)%25) + 1];
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
		T W = twiddles[24 + 4*((3*me + 1)%25) + 2];
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
		T W = twiddles[24 + 4*((3*me + 1)%25) + 3];
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
		T W = twiddles[24 + 4*((3*me + 2)%25) + 0];
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
		T W = twiddles[24 + 4*((3*me + 2)%25) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R12).x; ry = (*R12).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R12).x = TR;
		(*R12).y = TI;
	}

	{
		T W = twiddles[24 + 4*((3*me + 2)%25) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		T W = twiddles[24 + 4*((3*me + 2)%25) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	InvRad5B1(R0, R1, R2, R3, R4);
	InvRad5B1(R5, R6, R7, R8, R9);
	InvRad5B1(R10, R11, R12, R13, R14);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 25 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 50 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 75 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 100 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 25 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 50 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 75 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 100 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 25 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 50 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 75 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 100 ) ] = (*R14).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1125 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1125 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1125 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1125 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1125 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2250 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2250 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2250 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2250 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2250 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 25 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 50 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 75 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/25)*125 + (3*me + 0)%25 + 100 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 25 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 50 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 75 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 1)/25)*125 + (3*me + 1)%25 + 100 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 25 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 50 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 75 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 2)/25)*125 + (3*me + 2)%25 + 100 ) ] = (*R14).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1125 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1125 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1125 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1125 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1125 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2250 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2250 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2250 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2250 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2250 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass3_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{




	{
		T W = twiddles[124 + 2*((5*me + 0)%125) + 0];
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
		T W = twiddles[124 + 2*((5*me + 0)%125) + 1];
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
		T W = twiddles[124 + 2*((5*me + 1)%125) + 0];
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
		T W = twiddles[124 + 2*((5*me + 1)%125) + 1];
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
		T W = twiddles[124 + 2*((5*me + 2)%125) + 0];
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
		T W = twiddles[124 + 2*((5*me + 2)%125) + 1];
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
		T W = twiddles[124 + 2*((5*me + 3)%125) + 0];
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
		T W = twiddles[124 + 2*((5*me + 3)%125) + 1];
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
		T W = twiddles[124 + 2*((5*me + 4)%125) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		T W = twiddles[124 + 2*((5*me + 4)%125) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);
	InvRad3B1(R9, R10, R11);
	InvRad3B1(R12, R13, R14);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/125)*375 + (5*me + 0)%125 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*375 + (5*me + 0)%125 + 125 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*375 + (5*me + 0)%125 + 250 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*375 + (5*me + 1)%125 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*375 + (5*me + 1)%125 + 125 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*375 + (5*me + 1)%125 + 250 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*375 + (5*me + 2)%125 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*375 + (5*me + 2)%125 + 125 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*375 + (5*me + 2)%125 + 250 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*375 + (5*me + 3)%125 + 0 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*375 + (5*me + 3)%125 + 125 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*375 + (5*me + 3)%125 + 250 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*375 + (5*me + 4)%125 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*375 + (5*me + 4)%125 + 125 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*375 + (5*me + 4)%125 + 250 ) ] = (*R14).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1125 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1125 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1125 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1125 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1125 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2250 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2250 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2250 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2250 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2250 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/125)*375 + (5*me + 0)%125 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*375 + (5*me + 0)%125 + 125 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*375 + (5*me + 0)%125 + 250 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*375 + (5*me + 1)%125 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*375 + (5*me + 1)%125 + 125 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*375 + (5*me + 1)%125 + 250 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*375 + (5*me + 2)%125 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*375 + (5*me + 2)%125 + 125 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*375 + (5*me + 2)%125 + 250 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*375 + (5*me + 3)%125 + 0 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*375 + (5*me + 3)%125 + 125 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*375 + (5*me + 3)%125 + 250 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*375 + (5*me + 4)%125 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*375 + (5*me + 4)%125 + 125 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*375 + (5*me + 4)%125 + 250 ) ] = (*R14).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1125 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1125 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1125 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1125 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1125 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2250 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2250 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2250 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2250 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2250 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass4_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{




	{
		T W = twiddles[374 + 2*((5*me + 0)%375) + 0];
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
		T W = twiddles[374 + 2*((5*me + 0)%375) + 1];
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
		T W = twiddles[374 + 2*((5*me + 1)%375) + 0];
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
		T W = twiddles[374 + 2*((5*me + 1)%375) + 1];
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
		T W = twiddles[374 + 2*((5*me + 2)%375) + 0];
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
		T W = twiddles[374 + 2*((5*me + 2)%375) + 1];
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
		T W = twiddles[374 + 2*((5*me + 3)%375) + 0];
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
		T W = twiddles[374 + 2*((5*me + 3)%375) + 1];
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
		T W = twiddles[374 + 2*((5*me + 4)%375) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		T W = twiddles[374 + 2*((5*me + 4)%375) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);
	InvRad3B1(R9, R10, R11);
	InvRad3B1(R12, R13, R14);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/375)*1125 + (5*me + 0)%375 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/375)*1125 + (5*me + 0)%375 + 375 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/375)*1125 + (5*me + 0)%375 + 750 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 1)/375)*1125 + (5*me + 1)%375 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 1)/375)*1125 + (5*me + 1)%375 + 375 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/375)*1125 + (5*me + 1)%375 + 750 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 2)/375)*1125 + (5*me + 2)%375 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 2)/375)*1125 + (5*me + 2)%375 + 375 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 2)/375)*1125 + (5*me + 2)%375 + 750 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 3)/375)*1125 + (5*me + 3)%375 + 0 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 3)/375)*1125 + (5*me + 3)%375 + 375 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 3)/375)*1125 + (5*me + 3)%375 + 750 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 4)/375)*1125 + (5*me + 4)%375 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 4)/375)*1125 + (5*me + 4)%375 + 375 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 4)/375)*1125 + (5*me + 4)%375 + 750 ) ] = (*R14).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1125 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1125 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1125 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1125 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1125 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2250 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2250 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2250 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2250 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2250 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/375)*1125 + (5*me + 0)%375 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/375)*1125 + (5*me + 0)%375 + 375 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/375)*1125 + (5*me + 0)%375 + 750 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 1)/375)*1125 + (5*me + 1)%375 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 1)/375)*1125 + (5*me + 1)%375 + 375 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/375)*1125 + (5*me + 1)%375 + 750 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 2)/375)*1125 + (5*me + 2)%375 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 2)/375)*1125 + (5*me + 2)%375 + 375 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 2)/375)*1125 + (5*me + 2)%375 + 750 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 3)/375)*1125 + (5*me + 3)%375 + 0 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 3)/375)*1125 + (5*me + 3)%375 + 375 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 3)/375)*1125 + (5*me + 3)%375 + 750 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 4)/375)*1125 + (5*me + 4)%375 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 4)/375)*1125 + (5*me + 4)%375 + 375 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 4)/375)*1125 + (5*me + 4)%375 + 750 ) ] = (*R14).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1125 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1125 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1125 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1125 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1125 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2250 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2250 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2250 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2250 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2250 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass5_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{




	{
		T W = twiddles[1124 + 2*((5*me + 0)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 0)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 1)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 1)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 2)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 2)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 3)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 3)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 4)%1125) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		T W = twiddles[1124 + 2*((5*me + 4)%1125) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);
	InvRad3B1(R9, R10, R11);
	InvRad3B1(R12, R13, R14);


	if(rw)
	{
	bufOut[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0);
	bufOut[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R3);
	bufOut[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R6);
	bufOut[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R9);
	bufOut[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R12);
	bufOut[outOffset + ( 5*me + 0 + 1125 )*stride_out] = (*R1);
	bufOut[outOffset + ( 5*me + 1 + 1125 )*stride_out] = (*R4);
	bufOut[outOffset + ( 5*me + 2 + 1125 )*stride_out] = (*R7);
	bufOut[outOffset + ( 5*me + 3 + 1125 )*stride_out] = (*R10);
	bufOut[outOffset + ( 5*me + 4 + 1125 )*stride_out] = (*R13);
	bufOut[outOffset + ( 5*me + 0 + 2250 )*stride_out] = (*R2);
	bufOut[outOffset + ( 5*me + 1 + 2250 )*stride_out] = (*R5);
	bufOut[outOffset + ( 5*me + 2 + 2250 )*stride_out] = (*R8);
	bufOut[outOffset + ( 5*me + 3 + 2250 )*stride_out] = (*R11);
	bufOut[outOffset + ( 5*me + 4 + 2250 )*stride_out] = (*R14);
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass5_len3375(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14)
{




	{
		T W = twiddles[1124 + 2*((5*me + 0)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 0)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 1)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 1)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 2)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 2)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 3)%1125) + 0];
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
		T W = twiddles[1124 + 2*((5*me + 3)%1125) + 1];
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
		T W = twiddles[1124 + 2*((5*me + 4)%1125) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	{
		T W = twiddles[1124 + 2*((5*me + 4)%1125) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);
	InvRad3B1(R9, R10, R11);
	InvRad3B1(R12, R13, R14);


	if(rw)
	{
	bufOutRe[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R12).y;
	bufOutRe[outOffset + ( 5*me + 0 + 1125 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 5*me + 0 + 1125 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 5*me + 1 + 1125 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 5*me + 1 + 1125 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 5*me + 2 + 1125 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 5*me + 2 + 1125 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 5*me + 3 + 1125 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 5*me + 3 + 1125 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 5*me + 4 + 1125 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 5*me + 4 + 1125 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 5*me + 0 + 2250 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 5*me + 0 + 2250 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 5*me + 1 + 2250 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 5*me + 1 + 2250 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 5*me + 2 + 2250 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 5*me + 2 + 2250 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 5*me + 3 + 2250 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 5*me + 3 + 2250 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 5*me + 4 + 2250 )*stride_out] = (*R14).x;
	bufOutIm[outOffset + ( 5*me + 4 + 2250 )*stride_out] = (*R14).y;
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb>
__device__ void 
fwd_len3375_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14;
	FwdPass0_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass1_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass2_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass3_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass4_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass5_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
}

template <typename T, StrideBin sb>
__device__ void 
back_len3375_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14;
	InvPass0_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass1_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass2_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass3_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass4_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass5_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len3375_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14;
	FwdPass0_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass1_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass2_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass3_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass4_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass5_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
}

template <typename T, StrideBin sb>
__device__ void 
back_len3375_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14;
	InvPass0_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass1_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass2_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass3_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass4_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass5_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len3375_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14;
	FwdPass0_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass1_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass2_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass3_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass4_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass5_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
}

template <typename T, StrideBin sb>
__device__ void 
back_len3375_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14;
	InvPass0_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass1_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass2_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass3_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass4_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass5_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len3375_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14;
	FwdPass0_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass1_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass2_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass3_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass4_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	FwdPass5_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
}

template <typename T, StrideBin sb>
__device__ void 
back_len3375_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14;
	InvPass0_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass1_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass2_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass3_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass4_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
	InvPass5_len3375<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14);
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 225, maximum transforms: 1, Passes: 6
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(225)
fft_fwd_ip_len3375( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3375), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len3375>(
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
	size_t counter_mod = batch;
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
	fwd_len3375_device<T, sb>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, lwb, lwb, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 225, maximum transforms: 1, Passes: 6
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(225)
fft_back_ip_len3375( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3375), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len3375>(
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
	size_t counter_mod = batch;
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
	back_len3375_device<T, sb>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, lwb, lwb, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 225, maximum transforms: 1, Passes: 6
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(225)
fft_fwd_ip_len3375( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3375), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len3375>(
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
	size_t counter_mod = batch;
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
	fwd_len3375_device<T, sb>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, lwbRe, lwbIm, lwbRe, lwbIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 225, maximum transforms: 1, Passes: 6
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(225)
fft_back_ip_len3375( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3375), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len3375>(
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
	size_t counter_mod = batch;
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
	back_len3375_device<T, sb>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, lwbRe, lwbIm, lwbRe, lwbIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 225, maximum transforms: 1, Passes: 6
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(225)
fft_fwd_op_len3375( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3375), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len3375>(
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
	size_t counter_mod = batch;
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
	fwd_len3375_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbIn, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 225, maximum transforms: 1, Passes: 6
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(225)
fft_back_op_len3375( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3375), cgh);
	cgh.parallel_for<class kern_fft_back_op_len3375>(
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
	size_t counter_mod = batch;
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
	back_len3375_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbIn, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 225, maximum transforms: 1, Passes: 6
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(225)
fft_fwd_op_len3375( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3375), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len3375>(
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
	size_t counter_mod = batch;
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
	fwd_len3375_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbIn, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 225, maximum transforms: 1, Passes: 6
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(225)
fft_back_op_len3375( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3375), cgh);
	cgh.parallel_for<class kern_fft_back_op_len3375>(
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
	size_t counter_mod = batch;
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
	back_len3375_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbIn, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 225, maximum transforms: 1, Passes: 6
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(225)
fft_fwd_op_len3375( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3375), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len3375>(
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
	size_t counter_mod = batch;
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
	fwd_len3375_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbInRe, lwbInIm, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 225, maximum transforms: 1, Passes: 6
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(225)
fft_back_op_len3375( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3375), cgh);
	cgh.parallel_for<class kern_fft_back_op_len3375>(
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
	size_t counter_mod = batch;
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
	back_len3375_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbInRe, lwbInIm, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 225, maximum transforms: 1, Passes: 6
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(225)
fft_fwd_op_len3375( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3375), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len3375>(
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
	size_t counter_mod = batch;
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
	fwd_len3375_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbInRe, lwbInIm, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 225, maximum transforms: 1, Passes: 6
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(225)
fft_back_op_len3375( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3375), cgh);
	cgh.parallel_for<class kern_fft_back_op_len3375>(
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
	size_t counter_mod = batch;
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
	back_len3375_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbInRe, lwbInIm, lwbOutRe, lwbOutIm, lds);
});
});
}

