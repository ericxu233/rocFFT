#pragma once
#include "rocfft_butterfly_template.h"

////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass0_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	(*R0) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 0 )*stride_in, load_cb_data, nullptr);

	(*R3) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 0 )*stride_in, load_cb_data, nullptr);

	(*R6) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 0 )*stride_in, load_cb_data, nullptr);

	(*R1) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 729 )*stride_in, load_cb_data, nullptr);

	(*R4) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 729 )*stride_in, load_cb_data, nullptr);

	(*R7) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 729 )*stride_in, load_cb_data, nullptr);

	(*R2) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 1458 )*stride_in, load_cb_data, nullptr);

	(*R5) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 1458 )*stride_in, load_cb_data, nullptr);

	(*R8) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 1458 )*stride_in, load_cb_data, nullptr);

	}



	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 1 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 2 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 1 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 2 ) ] = (*R8).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 1 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 2 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 1 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 2 ) ] = (*R8).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass0_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 729 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 729 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 729 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 729 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 729 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 729 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1458 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1458 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1458 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1458 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1458 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1458 )*stride_in];
	}



	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 1 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 2 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 1 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 2 ) ] = (*R8).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 1 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 2 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 1 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 2 ) ] = (*R8).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass1_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8)
{




	{
		T W = twiddles[2 + 2*((3*me + 0)%3) + 0];
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
		T W = twiddles[2 + 2*((3*me + 0)%3) + 1];
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
		T W = twiddles[2 + 2*((3*me + 1)%3) + 0];
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
		T W = twiddles[2 + 2*((3*me + 1)%3) + 1];
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
		T W = twiddles[2 + 2*((3*me + 2)%3) + 0];
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
		T W = twiddles[2 + 2*((3*me + 2)%3) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/3)*9 + (3*me + 0)%3 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/3)*9 + (3*me + 0)%3 + 3 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/3)*9 + (3*me + 0)%3 + 6 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 1)/3)*9 + (3*me + 1)%3 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 1)/3)*9 + (3*me + 1)%3 + 3 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/3)*9 + (3*me + 1)%3 + 6 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 2)/3)*9 + (3*me + 2)%3 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 2)/3)*9 + (3*me + 2)%3 + 3 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 2)/3)*9 + (3*me + 2)%3 + 6 ) ] = (*R8).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/3)*9 + (3*me + 0)%3 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/3)*9 + (3*me + 0)%3 + 3 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/3)*9 + (3*me + 0)%3 + 6 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 1)/3)*9 + (3*me + 1)%3 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 1)/3)*9 + (3*me + 1)%3 + 3 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/3)*9 + (3*me + 1)%3 + 6 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 2)/3)*9 + (3*me + 2)%3 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 2)/3)*9 + (3*me + 2)%3 + 3 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 2)/3)*9 + (3*me + 2)%3 + 6 ) ] = (*R8).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass2_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8)
{




	{
		T W = twiddles[8 + 2*((3*me + 0)%9) + 0];
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
		T W = twiddles[8 + 2*((3*me + 0)%9) + 1];
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
		T W = twiddles[8 + 2*((3*me + 1)%9) + 0];
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
		T W = twiddles[8 + 2*((3*me + 1)%9) + 1];
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
		T W = twiddles[8 + 2*((3*me + 2)%9) + 0];
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
		T W = twiddles[8 + 2*((3*me + 2)%9) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/9)*27 + (3*me + 0)%9 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/9)*27 + (3*me + 0)%9 + 9 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/9)*27 + (3*me + 0)%9 + 18 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 1)/9)*27 + (3*me + 1)%9 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 1)/9)*27 + (3*me + 1)%9 + 9 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/9)*27 + (3*me + 1)%9 + 18 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 2)/9)*27 + (3*me + 2)%9 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 2)/9)*27 + (3*me + 2)%9 + 9 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 2)/9)*27 + (3*me + 2)%9 + 18 ) ] = (*R8).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/9)*27 + (3*me + 0)%9 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/9)*27 + (3*me + 0)%9 + 9 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/9)*27 + (3*me + 0)%9 + 18 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 1)/9)*27 + (3*me + 1)%9 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 1)/9)*27 + (3*me + 1)%9 + 9 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/9)*27 + (3*me + 1)%9 + 18 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 2)/9)*27 + (3*me + 2)%9 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 2)/9)*27 + (3*me + 2)%9 + 9 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 2)/9)*27 + (3*me + 2)%9 + 18 ) ] = (*R8).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass3_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8)
{




	{
		T W = twiddles[26 + 2*((3*me + 0)%27) + 0];
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
		T W = twiddles[26 + 2*((3*me + 0)%27) + 1];
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
		T W = twiddles[26 + 2*((3*me + 1)%27) + 0];
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
		T W = twiddles[26 + 2*((3*me + 1)%27) + 1];
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
		T W = twiddles[26 + 2*((3*me + 2)%27) + 0];
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
		T W = twiddles[26 + 2*((3*me + 2)%27) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/27)*81 + (3*me + 0)%27 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/27)*81 + (3*me + 0)%27 + 27 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/27)*81 + (3*me + 0)%27 + 54 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 1)/27)*81 + (3*me + 1)%27 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 1)/27)*81 + (3*me + 1)%27 + 27 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/27)*81 + (3*me + 1)%27 + 54 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 2)/27)*81 + (3*me + 2)%27 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 2)/27)*81 + (3*me + 2)%27 + 27 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 2)/27)*81 + (3*me + 2)%27 + 54 ) ] = (*R8).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/27)*81 + (3*me + 0)%27 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/27)*81 + (3*me + 0)%27 + 27 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/27)*81 + (3*me + 0)%27 + 54 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 1)/27)*81 + (3*me + 1)%27 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 1)/27)*81 + (3*me + 1)%27 + 27 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/27)*81 + (3*me + 1)%27 + 54 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 2)/27)*81 + (3*me + 2)%27 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 2)/27)*81 + (3*me + 2)%27 + 27 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 2)/27)*81 + (3*me + 2)%27 + 54 ) ] = (*R8).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass4_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8)
{




	{
		T W = twiddles[80 + 2*((3*me + 0)%81) + 0];
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
		T W = twiddles[80 + 2*((3*me + 0)%81) + 1];
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
		T W = twiddles[80 + 2*((3*me + 1)%81) + 0];
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
		T W = twiddles[80 + 2*((3*me + 1)%81) + 1];
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
		T W = twiddles[80 + 2*((3*me + 2)%81) + 0];
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
		T W = twiddles[80 + 2*((3*me + 2)%81) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/81)*243 + (3*me + 0)%81 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/81)*243 + (3*me + 0)%81 + 81 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/81)*243 + (3*me + 0)%81 + 162 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 1)/81)*243 + (3*me + 1)%81 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 1)/81)*243 + (3*me + 1)%81 + 81 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/81)*243 + (3*me + 1)%81 + 162 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 2)/81)*243 + (3*me + 2)%81 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 2)/81)*243 + (3*me + 2)%81 + 81 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 2)/81)*243 + (3*me + 2)%81 + 162 ) ] = (*R8).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/81)*243 + (3*me + 0)%81 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/81)*243 + (3*me + 0)%81 + 81 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/81)*243 + (3*me + 0)%81 + 162 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 1)/81)*243 + (3*me + 1)%81 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 1)/81)*243 + (3*me + 1)%81 + 81 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/81)*243 + (3*me + 1)%81 + 162 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 2)/81)*243 + (3*me + 2)%81 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 2)/81)*243 + (3*me + 2)%81 + 81 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 2)/81)*243 + (3*me + 2)%81 + 162 ) ] = (*R8).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass5_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8)
{




	{
		T W = twiddles[242 + 2*((3*me + 0)%243) + 0];
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
		T W = twiddles[242 + 2*((3*me + 0)%243) + 1];
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
		T W = twiddles[242 + 2*((3*me + 1)%243) + 0];
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
		T W = twiddles[242 + 2*((3*me + 1)%243) + 1];
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
		T W = twiddles[242 + 2*((3*me + 2)%243) + 0];
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
		T W = twiddles[242 + 2*((3*me + 2)%243) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/243)*729 + (3*me + 0)%243 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/243)*729 + (3*me + 0)%243 + 243 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/243)*729 + (3*me + 0)%243 + 486 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 1)/243)*729 + (3*me + 1)%243 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 1)/243)*729 + (3*me + 1)%243 + 243 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/243)*729 + (3*me + 1)%243 + 486 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 2)/243)*729 + (3*me + 2)%243 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 2)/243)*729 + (3*me + 2)%243 + 243 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 2)/243)*729 + (3*me + 2)%243 + 486 ) ] = (*R8).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/243)*729 + (3*me + 0)%243 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/243)*729 + (3*me + 0)%243 + 243 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/243)*729 + (3*me + 0)%243 + 486 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 1)/243)*729 + (3*me + 1)%243 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 1)/243)*729 + (3*me + 1)%243 + 243 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/243)*729 + (3*me + 1)%243 + 486 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 2)/243)*729 + (3*me + 2)%243 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 2)/243)*729 + (3*me + 2)%243 + 243 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 2)/243)*729 + (3*me + 2)%243 + 486 ) ] = (*R8).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass6_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[728 + 2*((3*me + 0)%729) + 0];
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
		T W = twiddles[728 + 2*((3*me + 0)%729) + 1];
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
		T W = twiddles[728 + 2*((3*me + 1)%729) + 0];
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
		T W = twiddles[728 + 2*((3*me + 1)%729) + 1];
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
		T W = twiddles[728 + 2*((3*me + 2)%729) + 0];
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
		T W = twiddles[728 + 2*((3*me + 2)%729) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);


	if(rw)
	{
	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);

	store_cb(bufOut, outOffset + ( 3*me + 0 + 0 )*stride_out, (*R0), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 1 + 0 )*stride_out, (*R3), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 2 + 0 )*stride_out, (*R6), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 0 + 729 )*stride_out, (*R1), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 1 + 729 )*stride_out, (*R4), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 2 + 729 )*stride_out, (*R7), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 0 + 1458 )*stride_out, (*R2), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 1 + 1458 )*stride_out, (*R5), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 2 + 1458 )*stride_out, (*R8), store_cb_data, nullptr);

	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass6_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[728 + 2*((3*me + 0)%729) + 0];
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
		T W = twiddles[728 + 2*((3*me + 0)%729) + 1];
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
		T W = twiddles[728 + 2*((3*me + 1)%729) + 0];
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
		T W = twiddles[728 + 2*((3*me + 1)%729) + 1];
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
		T W = twiddles[728 + 2*((3*me + 2)%729) + 0];
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
		T W = twiddles[728 + 2*((3*me + 2)%729) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( 3*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 3*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 3*me + 1 + 0 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 3*me + 1 + 0 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 3*me + 2 + 0 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 3*me + 2 + 0 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 3*me + 0 + 729 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 3*me + 0 + 729 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 3*me + 1 + 729 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 3*me + 1 + 729 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 3*me + 2 + 729 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 3*me + 2 + 729 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 3*me + 0 + 1458 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 3*me + 0 + 1458 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 3*me + 1 + 1458 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 3*me + 1 + 1458 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 3*me + 2 + 1458 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 3*me + 2 + 1458 )*stride_out] = (*R8).y;
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass0_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	(*R0) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 0 )*stride_in, load_cb_data, nullptr);

	(*R3) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 0 )*stride_in, load_cb_data, nullptr);

	(*R6) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 0 )*stride_in, load_cb_data, nullptr);

	(*R1) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 729 )*stride_in, load_cb_data, nullptr);

	(*R4) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 729 )*stride_in, load_cb_data, nullptr);

	(*R7) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 729 )*stride_in, load_cb_data, nullptr);

	(*R2) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 1458 )*stride_in, load_cb_data, nullptr);

	(*R5) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 1458 )*stride_in, load_cb_data, nullptr);

	(*R8) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 1458 )*stride_in, load_cb_data, nullptr);

	}



	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 1 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 2 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 1 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 2 ) ] = (*R8).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 1 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 2 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 1 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 2 ) ] = (*R8).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass0_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 729 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 729 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 729 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 729 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 729 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 729 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1458 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1458 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1458 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1458 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1458 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1458 )*stride_in];
	}



	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 1 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 2 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 1 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 2 ) ] = (*R8).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*3 + (3*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 1 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*3 + (3*me + 1)%1 + 2 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 1 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*3 + (3*me + 2)%1 + 2 ) ] = (*R8).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass1_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8)
{




	{
		T W = twiddles[2 + 2*((3*me + 0)%3) + 0];
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
		T W = twiddles[2 + 2*((3*me + 0)%3) + 1];
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
		T W = twiddles[2 + 2*((3*me + 1)%3) + 0];
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
		T W = twiddles[2 + 2*((3*me + 1)%3) + 1];
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
		T W = twiddles[2 + 2*((3*me + 2)%3) + 0];
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
		T W = twiddles[2 + 2*((3*me + 2)%3) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/3)*9 + (3*me + 0)%3 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/3)*9 + (3*me + 0)%3 + 3 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/3)*9 + (3*me + 0)%3 + 6 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 1)/3)*9 + (3*me + 1)%3 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 1)/3)*9 + (3*me + 1)%3 + 3 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/3)*9 + (3*me + 1)%3 + 6 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 2)/3)*9 + (3*me + 2)%3 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 2)/3)*9 + (3*me + 2)%3 + 3 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 2)/3)*9 + (3*me + 2)%3 + 6 ) ] = (*R8).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/3)*9 + (3*me + 0)%3 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/3)*9 + (3*me + 0)%3 + 3 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/3)*9 + (3*me + 0)%3 + 6 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 1)/3)*9 + (3*me + 1)%3 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 1)/3)*9 + (3*me + 1)%3 + 3 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/3)*9 + (3*me + 1)%3 + 6 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 2)/3)*9 + (3*me + 2)%3 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 2)/3)*9 + (3*me + 2)%3 + 3 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 2)/3)*9 + (3*me + 2)%3 + 6 ) ] = (*R8).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass2_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8)
{




	{
		T W = twiddles[8 + 2*((3*me + 0)%9) + 0];
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
		T W = twiddles[8 + 2*((3*me + 0)%9) + 1];
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
		T W = twiddles[8 + 2*((3*me + 1)%9) + 0];
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
		T W = twiddles[8 + 2*((3*me + 1)%9) + 1];
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
		T W = twiddles[8 + 2*((3*me + 2)%9) + 0];
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
		T W = twiddles[8 + 2*((3*me + 2)%9) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/9)*27 + (3*me + 0)%9 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/9)*27 + (3*me + 0)%9 + 9 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/9)*27 + (3*me + 0)%9 + 18 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 1)/9)*27 + (3*me + 1)%9 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 1)/9)*27 + (3*me + 1)%9 + 9 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/9)*27 + (3*me + 1)%9 + 18 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 2)/9)*27 + (3*me + 2)%9 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 2)/9)*27 + (3*me + 2)%9 + 9 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 2)/9)*27 + (3*me + 2)%9 + 18 ) ] = (*R8).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/9)*27 + (3*me + 0)%9 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/9)*27 + (3*me + 0)%9 + 9 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/9)*27 + (3*me + 0)%9 + 18 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 1)/9)*27 + (3*me + 1)%9 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 1)/9)*27 + (3*me + 1)%9 + 9 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/9)*27 + (3*me + 1)%9 + 18 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 2)/9)*27 + (3*me + 2)%9 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 2)/9)*27 + (3*me + 2)%9 + 9 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 2)/9)*27 + (3*me + 2)%9 + 18 ) ] = (*R8).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass3_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8)
{




	{
		T W = twiddles[26 + 2*((3*me + 0)%27) + 0];
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
		T W = twiddles[26 + 2*((3*me + 0)%27) + 1];
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
		T W = twiddles[26 + 2*((3*me + 1)%27) + 0];
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
		T W = twiddles[26 + 2*((3*me + 1)%27) + 1];
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
		T W = twiddles[26 + 2*((3*me + 2)%27) + 0];
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
		T W = twiddles[26 + 2*((3*me + 2)%27) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/27)*81 + (3*me + 0)%27 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/27)*81 + (3*me + 0)%27 + 27 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/27)*81 + (3*me + 0)%27 + 54 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 1)/27)*81 + (3*me + 1)%27 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 1)/27)*81 + (3*me + 1)%27 + 27 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/27)*81 + (3*me + 1)%27 + 54 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 2)/27)*81 + (3*me + 2)%27 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 2)/27)*81 + (3*me + 2)%27 + 27 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 2)/27)*81 + (3*me + 2)%27 + 54 ) ] = (*R8).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/27)*81 + (3*me + 0)%27 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/27)*81 + (3*me + 0)%27 + 27 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/27)*81 + (3*me + 0)%27 + 54 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 1)/27)*81 + (3*me + 1)%27 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 1)/27)*81 + (3*me + 1)%27 + 27 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/27)*81 + (3*me + 1)%27 + 54 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 2)/27)*81 + (3*me + 2)%27 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 2)/27)*81 + (3*me + 2)%27 + 27 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 2)/27)*81 + (3*me + 2)%27 + 54 ) ] = (*R8).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass4_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8)
{




	{
		T W = twiddles[80 + 2*((3*me + 0)%81) + 0];
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
		T W = twiddles[80 + 2*((3*me + 0)%81) + 1];
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
		T W = twiddles[80 + 2*((3*me + 1)%81) + 0];
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
		T W = twiddles[80 + 2*((3*me + 1)%81) + 1];
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
		T W = twiddles[80 + 2*((3*me + 2)%81) + 0];
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
		T W = twiddles[80 + 2*((3*me + 2)%81) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/81)*243 + (3*me + 0)%81 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/81)*243 + (3*me + 0)%81 + 81 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/81)*243 + (3*me + 0)%81 + 162 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 1)/81)*243 + (3*me + 1)%81 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 1)/81)*243 + (3*me + 1)%81 + 81 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/81)*243 + (3*me + 1)%81 + 162 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 2)/81)*243 + (3*me + 2)%81 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 2)/81)*243 + (3*me + 2)%81 + 81 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 2)/81)*243 + (3*me + 2)%81 + 162 ) ] = (*R8).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/81)*243 + (3*me + 0)%81 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/81)*243 + (3*me + 0)%81 + 81 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/81)*243 + (3*me + 0)%81 + 162 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 1)/81)*243 + (3*me + 1)%81 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 1)/81)*243 + (3*me + 1)%81 + 81 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/81)*243 + (3*me + 1)%81 + 162 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 2)/81)*243 + (3*me + 2)%81 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 2)/81)*243 + (3*me + 2)%81 + 81 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 2)/81)*243 + (3*me + 2)%81 + 162 ) ] = (*R8).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass5_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8)
{




	{
		T W = twiddles[242 + 2*((3*me + 0)%243) + 0];
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
		T W = twiddles[242 + 2*((3*me + 0)%243) + 1];
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
		T W = twiddles[242 + 2*((3*me + 1)%243) + 0];
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
		T W = twiddles[242 + 2*((3*me + 1)%243) + 1];
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
		T W = twiddles[242 + 2*((3*me + 2)%243) + 0];
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
		T W = twiddles[242 + 2*((3*me + 2)%243) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/243)*729 + (3*me + 0)%243 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/243)*729 + (3*me + 0)%243 + 243 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/243)*729 + (3*me + 0)%243 + 486 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 1)/243)*729 + (3*me + 1)%243 + 0 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 1)/243)*729 + (3*me + 1)%243 + 243 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 1)/243)*729 + (3*me + 1)%243 + 486 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 2)/243)*729 + (3*me + 2)%243 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 2)/243)*729 + (3*me + 2)%243 + 243 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 2)/243)*729 + (3*me + 2)%243 + 486 ) ] = (*R8).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/243)*729 + (3*me + 0)%243 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/243)*729 + (3*me + 0)%243 + 243 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/243)*729 + (3*me + 0)%243 + 486 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 1)/243)*729 + (3*me + 1)%243 + 0 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 1)/243)*729 + (3*me + 1)%243 + 243 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 1)/243)*729 + (3*me + 1)%243 + 486 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 2)/243)*729 + (3*me + 2)%243 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 2)/243)*729 + (3*me + 2)%243 + 243 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 2)/243)*729 + (3*me + 2)%243 + 486 ) ] = (*R8).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 729 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 729 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 729 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1458 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1458 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1458 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass6_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[728 + 2*((3*me + 0)%729) + 0];
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
		T W = twiddles[728 + 2*((3*me + 0)%729) + 1];
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
		T W = twiddles[728 + 2*((3*me + 1)%729) + 0];
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
		T W = twiddles[728 + 2*((3*me + 1)%729) + 1];
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
		T W = twiddles[728 + 2*((3*me + 2)%729) + 0];
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
		T W = twiddles[728 + 2*((3*me + 2)%729) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);


	if(rw)
	{
	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);

	store_cb(bufOut, outOffset + ( 3*me + 0 + 0 )*stride_out, (*R0), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 1 + 0 )*stride_out, (*R3), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 2 + 0 )*stride_out, (*R6), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 0 + 729 )*stride_out, (*R1), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 1 + 729 )*stride_out, (*R4), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 2 + 729 )*stride_out, (*R7), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 0 + 1458 )*stride_out, (*R2), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 1 + 1458 )*stride_out, (*R5), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 3*me + 2 + 1458 )*stride_out, (*R8), store_cb_data, nullptr);

	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass6_len2187(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[728 + 2*((3*me + 0)%729) + 0];
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
		T W = twiddles[728 + 2*((3*me + 0)%729) + 1];
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
		T W = twiddles[728 + 2*((3*me + 1)%729) + 0];
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
		T W = twiddles[728 + 2*((3*me + 1)%729) + 1];
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
		T W = twiddles[728 + 2*((3*me + 2)%729) + 0];
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
		T W = twiddles[728 + 2*((3*me + 2)%729) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);


	if(rw)
	{
	bufOutRe[outOffset + ( 3*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 3*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 3*me + 1 + 0 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 3*me + 1 + 0 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 3*me + 2 + 0 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 3*me + 2 + 0 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 3*me + 0 + 729 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 3*me + 0 + 729 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 3*me + 1 + 729 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 3*me + 1 + 729 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 3*me + 2 + 729 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 3*me + 2 + 729 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 3*me + 0 + 1458 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 3*me + 0 + 1458 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 3*me + 1 + 1458 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 3*me + 1 + 1458 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 3*me + 2 + 1458 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 3*me + 2 + 1458 )*stride_out] = (*R8).y;
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len2187_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8;
	FwdPass0_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_fn, load_cb_data);
	FwdPass1_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass2_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass3_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass4_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass5_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass6_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len2187_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8;
	InvPass0_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_fn, load_cb_data);
	InvPass1_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass2_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass3_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass4_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass5_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass6_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len2187_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8;
	FwdPass0_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_fn, load_cb_data);
	FwdPass1_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass2_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass3_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass4_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass5_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass6_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len2187_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8;
	InvPass0_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_fn, load_cb_data);
	InvPass1_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass2_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass3_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass4_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass5_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass6_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len2187_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8;
	FwdPass0_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_fn, load_cb_data);
	FwdPass1_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass2_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass3_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass4_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass5_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass6_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len2187_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8;
	InvPass0_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_fn, load_cb_data);
	InvPass1_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass2_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass3_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass4_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass5_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass6_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len2187_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8;
	FwdPass0_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_fn, load_cb_data);
	FwdPass1_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass2_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass3_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass4_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass5_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	FwdPass6_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len2187_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8;
	InvPass0_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_fn, load_cb_data);
	InvPass1_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass2_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass3_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass4_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass5_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8);
	InvPass6_len2187<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 243, maximum transforms: 1, Passes: 7
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(243)
fft_fwd_ip_len2187( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gb)
{

	__shared__ real_type_t<T> lds[2187];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int ioOffset = 0;
	
	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
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
	// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, gb, lds
	fwd_len2187_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, gb, ioOffset, gb, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 243, maximum transforms: 1, Passes: 7
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(243)
fft_back_ip_len2187( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gb)
{

	__shared__ real_type_t<T> lds[2187];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int ioOffset = 0;
	
	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
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
	// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, gb, lds
	back_len2187_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, gb, ioOffset, gb, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 243, maximum transforms: 1, Passes: 7
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(243)
fft_fwd_ip_len2187( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbRe, real_type_t<T> * __restrict__ gbIm)
{

	__shared__ real_type_t<T> lds[2187];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int ioOffset = 0;
	
	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
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
	// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, gb, lds
	fwd_len2187_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, gbRe, gbIm, ioOffset, gbRe, gbIm, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 243, maximum transforms: 1, Passes: 7
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(243)
fft_back_ip_len2187( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbRe, real_type_t<T> * __restrict__ gbIm)
{

	__shared__ real_type_t<T> lds[2187];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int ioOffset = 0;
	
	unsigned int rw = 1;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
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
	// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, gb, lds
	back_len2187_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, gbRe, gbIm, ioOffset, gbRe, gbIm, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 243, maximum transforms: 1, Passes: 7
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(243)
fft_fwd_op_len2187( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[2187];
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
	// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, gb, lds
	fwd_len2187_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 243, maximum transforms: 1, Passes: 7
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(243)
fft_back_op_len2187( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[2187];
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
	// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, gb, lds
	back_len2187_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 243, maximum transforms: 1, Passes: 7
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(243)
fft_fwd_op_len2187( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[2187];
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
	// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, gb, lds
	fwd_len2187_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 243, maximum transforms: 1, Passes: 7
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(243)
fft_back_op_len2187( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[2187];
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
	// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, gb, lds
	back_len2187_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 243, maximum transforms: 1, Passes: 7
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(243)
fft_fwd_op_len2187( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[2187];
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
	// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, gb, lds
	fwd_len2187_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 243, maximum transforms: 1, Passes: 7
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(243)
fft_back_op_len2187( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[2187];
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
	// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, gb, lds
	back_len2187_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 243, maximum transforms: 1, Passes: 7
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(243)
fft_fwd_op_len2187( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[2187];
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
	// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, gb, lds
	fwd_len2187_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbInRe, gbInIm, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 243, maximum transforms: 1, Passes: 7
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(243)
fft_back_op_len2187( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[2187];
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
	// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
	// rw, b, me% control read/write; then ldsOffset, gb, lds
	back_len2187_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbInRe, gbInIm, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

