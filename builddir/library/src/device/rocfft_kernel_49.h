#pragma once
#include "rocfft_butterfly_template.h"

////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass0_len49(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	(*R0) = load_cb(bufIn, inOffset + ( 0 + me*1 + 0 + 0 )*stride_in, load_cb_data, nullptr);

	(*R1) = load_cb(bufIn, inOffset + ( 0 + me*1 + 0 + 7 )*stride_in, load_cb_data, nullptr);

	(*R2) = load_cb(bufIn, inOffset + ( 0 + me*1 + 0 + 14 )*stride_in, load_cb_data, nullptr);

	(*R3) = load_cb(bufIn, inOffset + ( 0 + me*1 + 0 + 21 )*stride_in, load_cb_data, nullptr);

	(*R4) = load_cb(bufIn, inOffset + ( 0 + me*1 + 0 + 28 )*stride_in, load_cb_data, nullptr);

	(*R5) = load_cb(bufIn, inOffset + ( 0 + me*1 + 0 + 35 )*stride_in, load_cb_data, nullptr);

	(*R6) = load_cb(bufIn, inOffset + ( 0 + me*1 + 0 + 42 )*stride_in, load_cb_data, nullptr);

	}



	FwdRad7B1(R0, R1, R2, R3, R4, R5, R6);


	if(rw)
	{
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 6 ) ] = (*R6).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 7 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 14 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 21 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 28 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 35 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 42 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 6 ) ] = (*R6).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 7 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 14 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 21 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 28 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 35 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 42 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass0_len49(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 7 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 7 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 14 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 14 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 21 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 21 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 28 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 28 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 35 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 35 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 42 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 42 )*stride_in];
	}



	FwdRad7B1(R0, R1, R2, R3, R4, R5, R6);


	if(rw)
	{
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 6 ) ] = (*R6).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 7 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 14 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 21 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 28 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 35 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 42 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 6 ) ] = (*R6).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 7 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 14 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 21 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 28 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 35 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 42 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass1_len49(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[6 + 6*((1*me + 0)%7) + 0];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 1];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 2];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 3];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 4];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 5];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R6).x; ry = (*R6).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R6).x = TR;
		(*R6).y = TI;
	}

	FwdRad7B1(R0, R1, R2, R3, R4, R5, R6);


	if(rw)
	{
	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);

	store_cb(bufOut, outOffset + ( 1*me + 0 + 0 )*stride_out, (*R0), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 1*me + 0 + 7 )*stride_out, (*R1), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 1*me + 0 + 14 )*stride_out, (*R2), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 1*me + 0 + 21 )*stride_out, (*R3), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 1*me + 0 + 28 )*stride_out, (*R4), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 1*me + 0 + 35 )*stride_out, (*R5), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 1*me + 0 + 42 )*stride_out, (*R6), store_cb_data, nullptr);

	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass1_len49(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[6 + 6*((1*me + 0)%7) + 0];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 1];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 2];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 3];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 4];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 5];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R6).x; ry = (*R6).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R6).x = TR;
		(*R6).y = TI;
	}

	FwdRad7B1(R0, R1, R2, R3, R4, R5, R6);


	if(rw)
	{
	bufOutRe[outOffset + ( 1*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 1*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 1*me + 0 + 7 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 1*me + 0 + 7 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 1*me + 0 + 14 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 1*me + 0 + 14 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 1*me + 0 + 21 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 1*me + 0 + 21 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 1*me + 0 + 28 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 1*me + 0 + 28 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 1*me + 0 + 35 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 1*me + 0 + 35 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 1*me + 0 + 42 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 1*me + 0 + 42 )*stride_out] = (*R6).y;
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass0_len49(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	(*R0) = load_cb(bufIn, inOffset + ( 0 + me*1 + 0 + 0 )*stride_in, load_cb_data, nullptr);

	(*R1) = load_cb(bufIn, inOffset + ( 0 + me*1 + 0 + 7 )*stride_in, load_cb_data, nullptr);

	(*R2) = load_cb(bufIn, inOffset + ( 0 + me*1 + 0 + 14 )*stride_in, load_cb_data, nullptr);

	(*R3) = load_cb(bufIn, inOffset + ( 0 + me*1 + 0 + 21 )*stride_in, load_cb_data, nullptr);

	(*R4) = load_cb(bufIn, inOffset + ( 0 + me*1 + 0 + 28 )*stride_in, load_cb_data, nullptr);

	(*R5) = load_cb(bufIn, inOffset + ( 0 + me*1 + 0 + 35 )*stride_in, load_cb_data, nullptr);

	(*R6) = load_cb(bufIn, inOffset + ( 0 + me*1 + 0 + 42 )*stride_in, load_cb_data, nullptr);

	}



	InvRad7B1(R0, R1, R2, R3, R4, R5, R6);


	if(rw)
	{
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 6 ) ] = (*R6).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 7 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 14 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 21 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 28 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 35 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 42 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 6 ) ] = (*R6).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 7 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 14 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 21 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 28 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 35 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 42 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass0_len49(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 7 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 7 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 14 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 14 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 21 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 21 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 28 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 28 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 35 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 35 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*1 + 0 + 42 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*1 + 0 + 42 )*stride_in];
	}



	InvRad7B1(R0, R1, R2, R3, R4, R5, R6);


	if(rw)
	{
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 6 ) ] = (*R6).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 7 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 14 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 21 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 28 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 35 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*1 + 0 + 42 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((1*me + 0)/1)*7 + (1*me + 0)%1 + 6 ) ] = (*R6).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 7 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 14 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 21 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 28 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 35 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*1 + 0 + 42 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass1_len49(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[6 + 6*((1*me + 0)%7) + 0];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 1];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 2];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 3];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 4];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 5];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R6).x; ry = (*R6).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R6).x = TR;
		(*R6).y = TI;
	}

	InvRad7B1(R0, R1, R2, R3, R4, R5, R6);


	if(rw)
	{
	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);

	store_cb(bufOut, outOffset + ( 1*me + 0 + 0 )*stride_out, (*R0), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 1*me + 0 + 7 )*stride_out, (*R1), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 1*me + 0 + 14 )*stride_out, (*R2), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 1*me + 0 + 21 )*stride_out, (*R3), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 1*me + 0 + 28 )*stride_out, (*R4), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 1*me + 0 + 35 )*stride_out, (*R5), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 1*me + 0 + 42 )*stride_out, (*R6), store_cb_data, nullptr);

	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass1_len49(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[6 + 6*((1*me + 0)%7) + 0];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 1];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 2];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 3];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 4];
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
		T W = twiddles[6 + 6*((1*me + 0)%7) + 5];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R6).x; ry = (*R6).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R6).x = TR;
		(*R6).y = TI;
	}

	InvRad7B1(R0, R1, R2, R3, R4, R5, R6);


	if(rw)
	{
	bufOutRe[outOffset + ( 1*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 1*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 1*me + 0 + 7 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 1*me + 0 + 7 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 1*me + 0 + 14 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 1*me + 0 + 14 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 1*me + 0 + 21 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 1*me + 0 + 21 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 1*me + 0 + 28 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 1*me + 0 + 28 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 1*me + 0 + 35 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 1*me + 0 + 35 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 1*me + 0 + 42 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 1*me + 0 + 42 )*stride_out] = (*R6).y;
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len49_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6;
	FwdPass0_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_fn, load_cb_data);
	FwdPass1_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len49_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6;
	InvPass0_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_fn, load_cb_data);
	InvPass1_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len49_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6;
	FwdPass0_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_fn, load_cb_data);
	FwdPass1_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len49_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6;
	InvPass0_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_fn, load_cb_data);
	InvPass1_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len49_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6;
	FwdPass0_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_fn, load_cb_data);
	FwdPass1_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len49_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6;
	InvPass0_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_fn, load_cb_data);
	InvPass1_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len49_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6;
	FwdPass0_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_fn, load_cb_data);
	FwdPass1_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len49_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6;
	InvPass0_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_fn, load_cb_data);
	InvPass1_len49<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 49, maximum transforms: 7, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(49)
fft_fwd_ip_len49( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gb)
{

	__shared__ real_type_t<T> lds[343];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int ioOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*7)*7) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*7 + (me/7));
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
	fwd_len49_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  rw, b, me%7, (me/7)*49, gb, ioOffset, gb, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 49, maximum transforms: 7, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(49)
fft_back_ip_len49( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gb)
{

	__shared__ real_type_t<T> lds[343];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int ioOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*7)*7) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*7 + (me/7));
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
	back_len49_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  rw, b, me%7, (me/7)*49, gb, ioOffset, gb, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 49, maximum transforms: 7, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(49)
fft_fwd_ip_len49( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbRe, real_type_t<T> * __restrict__ gbIm)
{

	__shared__ real_type_t<T> lds[343];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int ioOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*7)*7) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*7 + (me/7));
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
	fwd_len49_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  rw, b, me%7, (me/7)*49, gbRe, gbIm, ioOffset, gbRe, gbIm, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 49, maximum transforms: 7, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(49)
fft_back_ip_len49( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbRe, real_type_t<T> * __restrict__ gbIm)
{

	__shared__ real_type_t<T> lds[343];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int ioOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*7)*7) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*7 + (me/7));
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
	back_len49_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  rw, b, me%7, (me/7)*49, gbRe, gbIm, ioOffset, gbRe, gbIm, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 49, maximum transforms: 7, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(49)
fft_fwd_op_len49( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[343];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*7)*7) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*7 + (me/7));
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
	fwd_len49_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%7, (me/7)*49, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 49, maximum transforms: 7, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(49)
fft_back_op_len49( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[343];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*7)*7) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*7 + (me/7));
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
	back_len49_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%7, (me/7)*49, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 49, maximum transforms: 7, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(49)
fft_fwd_op_len49( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[343];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*7)*7) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*7 + (me/7));
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
	fwd_len49_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%7, (me/7)*49, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 49, maximum transforms: 7, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(49)
fft_back_op_len49( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[343];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*7)*7) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*7 + (me/7));
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
	back_len49_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%7, (me/7)*49, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 49, maximum transforms: 7, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(49)
fft_fwd_op_len49( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[343];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*7)*7) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*7 + (me/7));
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
	fwd_len49_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%7, (me/7)*49, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 49, maximum transforms: 7, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(49)
fft_back_op_len49( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[343];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*7)*7) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*7 + (me/7));
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
	back_len49_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%7, (me/7)*49, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 49, maximum transforms: 7, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(49)
fft_fwd_op_len49( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[343];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*7)*7) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*7 + (me/7));
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
	fwd_len49_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%7, (me/7)*49, gbInRe, gbInIm, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 49, maximum transforms: 7, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(49)
fft_back_op_len49( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[343];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*7)*7) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*7 + (me/7));
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
	back_len49_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%7, (me/7)*49, gbInRe, gbInIm, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

