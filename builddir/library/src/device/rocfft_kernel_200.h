#pragma once
#include "rocfft_butterfly_template.h"

////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass0_len200(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	(*R0) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 0 )*stride_in, load_cb_data, nullptr);

	(*R10) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 0 )*stride_in, load_cb_data, nullptr);

	(*R1) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 20 )*stride_in, load_cb_data, nullptr);

	(*R11) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 20 )*stride_in, load_cb_data, nullptr);

	(*R2) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 40 )*stride_in, load_cb_data, nullptr);

	(*R12) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 40 )*stride_in, load_cb_data, nullptr);

	(*R3) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 60 )*stride_in, load_cb_data, nullptr);

	(*R13) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 60 )*stride_in, load_cb_data, nullptr);

	(*R4) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 80 )*stride_in, load_cb_data, nullptr);

	(*R14) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 80 )*stride_in, load_cb_data, nullptr);

	(*R5) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 100 )*stride_in, load_cb_data, nullptr);

	(*R15) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 100 )*stride_in, load_cb_data, nullptr);

	(*R6) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 120 )*stride_in, load_cb_data, nullptr);

	(*R16) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 120 )*stride_in, load_cb_data, nullptr);

	(*R7) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 140 )*stride_in, load_cb_data, nullptr);

	(*R17) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 140 )*stride_in, load_cb_data, nullptr);

	(*R8) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 160 )*stride_in, load_cb_data, nullptr);

	(*R18) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 160 )*stride_in, load_cb_data, nullptr);

	(*R9) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 180 )*stride_in, load_cb_data, nullptr);

	(*R19) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 180 )*stride_in, load_cb_data, nullptr);

	}



	FwdRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	FwdRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);


	if(rw)
	{
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 4 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 5 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 6 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 7 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 8 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 9 ) ] = (*R19).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 20 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 20 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 40 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 40 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 60 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 60 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 80 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 80 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 100 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 100 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 120 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 120 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 140 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 140 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 160 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 160 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 180 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 180 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 4 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 5 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 6 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 7 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 8 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 9 ) ] = (*R19).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 20 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 20 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 40 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 40 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 60 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 60 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 80 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 80 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 100 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 100 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 120 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 120 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 140 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 140 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 160 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 160 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 180 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 180 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass0_len200(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 20 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 20 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 20 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 20 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 40 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 40 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 40 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 40 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 60 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 60 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 60 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 60 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 80 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 80 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 80 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 80 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 100 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 100 )*stride_in];
	(*R15).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 100 )*stride_in];
	(*R15).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 100 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 120 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 120 )*stride_in];
	(*R16).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 120 )*stride_in];
	(*R16).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 120 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 140 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 140 )*stride_in];
	(*R17).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 140 )*stride_in];
	(*R17).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 140 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 160 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 160 )*stride_in];
	(*R18).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 160 )*stride_in];
	(*R18).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 160 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 180 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 180 )*stride_in];
	(*R19).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 180 )*stride_in];
	(*R19).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 180 )*stride_in];
	}



	FwdRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	FwdRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);


	if(rw)
	{
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 4 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 5 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 6 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 7 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 8 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 9 ) ] = (*R19).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 20 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 20 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 40 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 40 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 60 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 60 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 80 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 80 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 100 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 100 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 120 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 120 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 140 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 140 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 160 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 160 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 180 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 180 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 4 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 5 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 6 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 7 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 8 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 9 ) ] = (*R19).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 20 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 20 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 40 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 40 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 60 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 60 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 80 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 80 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 100 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 100 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 120 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 120 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 140 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 140 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 160 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 160 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 180 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 180 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass1_len200(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{




	{
		T W = twiddles[9 + 9*((2*me + 0)%10) + 0];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 1];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 2];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 3];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 4];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 5];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 6];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 7];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 8];
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
		T W = twiddles[9 + 9*((2*me + 1)%10) + 0];
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
		T W = twiddles[9 + 9*((2*me + 1)%10) + 1];
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
		T W = twiddles[9 + 9*((2*me + 1)%10) + 2];
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
		T W = twiddles[9 + 9*((2*me + 1)%10) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	{
		T W = twiddles[9 + 9*((2*me + 1)%10) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R15).x; ry = (*R15).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R15).x = TR;
		(*R15).y = TI;
	}

	{
		T W = twiddles[9 + 9*((2*me + 1)%10) + 5];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R16).x; ry = (*R16).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R16).x = TR;
		(*R16).y = TI;
	}

	{
		T W = twiddles[9 + 9*((2*me + 1)%10) + 6];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R17).x; ry = (*R17).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R17).x = TR;
		(*R17).y = TI;
	}

	{
		T W = twiddles[9 + 9*((2*me + 1)%10) + 7];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R18).x; ry = (*R18).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R18).x = TR;
		(*R18).y = TI;
	}

	{
		T W = twiddles[9 + 9*((2*me + 1)%10) + 8];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	FwdRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	FwdRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);


	if(rw)
	{
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 10 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 20 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 30 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 40 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 50 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 60 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 70 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 80 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 90 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 10 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 20 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 30 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 40 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 50 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 60 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 70 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 80 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 90 ) ] = (*R19).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*10 + 0 + 0 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*10 + 1 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*10 + 2 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*10 + 3 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*10 + 4 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*10 + 5 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*10 + 6 + 0 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*10 + 7 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*10 + 8 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*10 + 9 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*10 + 0 + 100 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*10 + 1 + 100 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*10 + 2 + 100 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*10 + 3 + 100 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*10 + 4 + 100 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*10 + 5 + 100 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*10 + 6 + 100 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*10 + 7 + 100 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*10 + 8 + 100 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*10 + 9 + 100 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 10 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 20 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 30 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 40 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 50 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 60 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 70 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 80 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 90 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 10 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 20 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 30 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 40 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 50 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 60 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 70 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 80 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 90 ) ] = (*R19).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*10 + 0 + 0 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*10 + 1 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*10 + 2 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*10 + 3 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*10 + 4 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*10 + 5 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*10 + 6 + 0 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*10 + 7 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*10 + 8 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*10 + 9 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*10 + 0 + 100 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*10 + 1 + 100 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*10 + 2 + 100 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*10 + 3 + 100 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*10 + 4 + 100 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*10 + 5 + 100 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*10 + 6 + 100 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*10 + 7 + 100 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*10 + 8 + 100 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*10 + 9 + 100 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass2_len200(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[99 + 1*((10*me + 0)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 1)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 2)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 3)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 4)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 5)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 6)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 7)%100) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R15).x; ry = (*R15).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R15).x = TR;
		(*R15).y = TI;
	}

	{
		T W = twiddles[99 + 1*((10*me + 8)%100) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R17).x; ry = (*R17).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R17).x = TR;
		(*R17).y = TI;
	}

	{
		T W = twiddles[99 + 1*((10*me + 9)%100) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	FwdRad2B1(R0, R1);
	FwdRad2B1(R2, R3);
	FwdRad2B1(R4, R5);
	FwdRad2B1(R6, R7);
	FwdRad2B1(R8, R9);
	FwdRad2B1(R10, R11);
	FwdRad2B1(R12, R13);
	FwdRad2B1(R14, R15);
	FwdRad2B1(R16, R17);
	FwdRad2B1(R18, R19);


	if(rw)
	{
	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);

	 //Optimization: coalescing into float4/double4 write
	if(sb == SB_UNIT && !store_cb_fn) {
	vector4_type_t<T> *buff4g = (vector4_type_t<T>*)(bufOut+outOffset);
	
	buff4g[ 5*me + 0 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R0).x, (*R0).y, (*R2).x, (*R2).y) ;
	buff4g[ 5*me + 1 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R4).x, (*R4).y, (*R6).x, (*R6).y) ;
	buff4g[ 5*me + 2 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R8).x, (*R8).y, (*R10).x, (*R10).y) ;
	buff4g[ 5*me + 3 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R12).x, (*R12).y, (*R14).x, (*R14).y) ;
	buff4g[ 5*me + 4 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R16).x, (*R16).y, (*R18).x, (*R18).y) ;
	buff4g[ 5*me + 0 + 50 ] = lib_make_vector4< vector4_type_t<T> >((*R1).x, (*R1).y, (*R3).x, (*R3).y) ;
	buff4g[ 5*me + 1 + 50 ] = lib_make_vector4< vector4_type_t<T> >((*R5).x, (*R5).y, (*R7).x, (*R7).y) ;
	buff4g[ 5*me + 2 + 50 ] = lib_make_vector4< vector4_type_t<T> >((*R9).x, (*R9).y, (*R11).x, (*R11).y) ;
	buff4g[ 5*me + 3 + 50 ] = lib_make_vector4< vector4_type_t<T> >((*R13).x, (*R13).y, (*R15).x, (*R15).y) ;
	buff4g[ 5*me + 4 + 50 ] = lib_make_vector4< vector4_type_t<T> >((*R17).x, (*R17).y, (*R19).x, (*R19).y) ;
	}
	else{ // such optimization is not possible 
	store_cb(bufOut, outOffset + ( 10*me + 0 + 0 )*stride_out, (*R0), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 1 + 0 )*stride_out, (*R2), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 2 + 0 )*stride_out, (*R4), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 3 + 0 )*stride_out, (*R6), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 4 + 0 )*stride_out, (*R8), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 5 + 0 )*stride_out, (*R10), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 6 + 0 )*stride_out, (*R12), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 7 + 0 )*stride_out, (*R14), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 8 + 0 )*stride_out, (*R16), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 9 + 0 )*stride_out, (*R18), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 0 + 100 )*stride_out, (*R1), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 1 + 100 )*stride_out, (*R3), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 2 + 100 )*stride_out, (*R5), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 3 + 100 )*stride_out, (*R7), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 4 + 100 )*stride_out, (*R9), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 5 + 100 )*stride_out, (*R11), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 6 + 100 )*stride_out, (*R13), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 7 + 100 )*stride_out, (*R15), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 8 + 100 )*stride_out, (*R17), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 9 + 100 )*stride_out, (*R19), store_cb_data, nullptr);

	}
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass2_len200(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[99 + 1*((10*me + 0)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 1)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 2)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 3)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 4)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 5)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 6)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 7)%100) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R15).x; ry = (*R15).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R15).x = TR;
		(*R15).y = TI;
	}

	{
		T W = twiddles[99 + 1*((10*me + 8)%100) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R17).x; ry = (*R17).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R17).x = TR;
		(*R17).y = TI;
	}

	{
		T W = twiddles[99 + 1*((10*me + 9)%100) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	FwdRad2B1(R0, R1);
	FwdRad2B1(R2, R3);
	FwdRad2B1(R4, R5);
	FwdRad2B1(R6, R7);
	FwdRad2B1(R8, R9);
	FwdRad2B1(R10, R11);
	FwdRad2B1(R12, R13);
	FwdRad2B1(R14, R15);
	FwdRad2B1(R16, R17);
	FwdRad2B1(R18, R19);


	if(rw)
	{
	bufOutRe[outOffset + ( 10*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 10*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 10*me + 1 + 0 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 10*me + 1 + 0 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 10*me + 2 + 0 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 10*me + 2 + 0 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 10*me + 3 + 0 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 10*me + 3 + 0 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 10*me + 4 + 0 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 10*me + 4 + 0 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 10*me + 5 + 0 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 10*me + 5 + 0 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 10*me + 6 + 0 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 10*me + 6 + 0 )*stride_out] = (*R12).y;
	bufOutRe[outOffset + ( 10*me + 7 + 0 )*stride_out] = (*R14).x;
	bufOutIm[outOffset + ( 10*me + 7 + 0 )*stride_out] = (*R14).y;
	bufOutRe[outOffset + ( 10*me + 8 + 0 )*stride_out] = (*R16).x;
	bufOutIm[outOffset + ( 10*me + 8 + 0 )*stride_out] = (*R16).y;
	bufOutRe[outOffset + ( 10*me + 9 + 0 )*stride_out] = (*R18).x;
	bufOutIm[outOffset + ( 10*me + 9 + 0 )*stride_out] = (*R18).y;
	bufOutRe[outOffset + ( 10*me + 0 + 100 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 10*me + 0 + 100 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 10*me + 1 + 100 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 10*me + 1 + 100 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 10*me + 2 + 100 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 10*me + 2 + 100 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 10*me + 3 + 100 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 10*me + 3 + 100 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 10*me + 4 + 100 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 10*me + 4 + 100 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 10*me + 5 + 100 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 10*me + 5 + 100 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 10*me + 6 + 100 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 10*me + 6 + 100 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 10*me + 7 + 100 )*stride_out] = (*R15).x;
	bufOutIm[outOffset + ( 10*me + 7 + 100 )*stride_out] = (*R15).y;
	bufOutRe[outOffset + ( 10*me + 8 + 100 )*stride_out] = (*R17).x;
	bufOutIm[outOffset + ( 10*me + 8 + 100 )*stride_out] = (*R17).y;
	bufOutRe[outOffset + ( 10*me + 9 + 100 )*stride_out] = (*R19).x;
	bufOutIm[outOffset + ( 10*me + 9 + 100 )*stride_out] = (*R19).y;
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass0_len200(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	(*R0) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 0 )*stride_in, load_cb_data, nullptr);

	(*R10) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 0 )*stride_in, load_cb_data, nullptr);

	(*R1) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 20 )*stride_in, load_cb_data, nullptr);

	(*R11) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 20 )*stride_in, load_cb_data, nullptr);

	(*R2) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 40 )*stride_in, load_cb_data, nullptr);

	(*R12) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 40 )*stride_in, load_cb_data, nullptr);

	(*R3) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 60 )*stride_in, load_cb_data, nullptr);

	(*R13) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 60 )*stride_in, load_cb_data, nullptr);

	(*R4) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 80 )*stride_in, load_cb_data, nullptr);

	(*R14) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 80 )*stride_in, load_cb_data, nullptr);

	(*R5) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 100 )*stride_in, load_cb_data, nullptr);

	(*R15) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 100 )*stride_in, load_cb_data, nullptr);

	(*R6) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 120 )*stride_in, load_cb_data, nullptr);

	(*R16) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 120 )*stride_in, load_cb_data, nullptr);

	(*R7) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 140 )*stride_in, load_cb_data, nullptr);

	(*R17) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 140 )*stride_in, load_cb_data, nullptr);

	(*R8) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 160 )*stride_in, load_cb_data, nullptr);

	(*R18) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 160 )*stride_in, load_cb_data, nullptr);

	(*R9) = load_cb(bufIn, inOffset + ( 0 + me*2 + 0 + 180 )*stride_in, load_cb_data, nullptr);

	(*R19) = load_cb(bufIn, inOffset + ( 0 + me*2 + 1 + 180 )*stride_in, load_cb_data, nullptr);

	}



	InvRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	InvRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);


	if(rw)
	{
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 4 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 5 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 6 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 7 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 8 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 9 ) ] = (*R19).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 20 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 20 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 40 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 40 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 60 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 60 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 80 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 80 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 100 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 100 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 120 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 120 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 140 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 140 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 160 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 160 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 180 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 180 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 4 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 5 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 6 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 7 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 8 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 9 ) ] = (*R19).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 20 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 20 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 40 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 40 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 60 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 60 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 80 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 80 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 100 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 100 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 120 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 120 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 140 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 140 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 160 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 160 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 180 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 180 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass0_len200(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 20 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 20 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 20 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 20 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 40 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 40 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 40 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 40 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 60 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 60 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 60 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 60 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 80 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 80 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 80 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 80 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 100 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 100 )*stride_in];
	(*R15).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 100 )*stride_in];
	(*R15).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 100 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 120 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 120 )*stride_in];
	(*R16).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 120 )*stride_in];
	(*R16).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 120 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 140 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 140 )*stride_in];
	(*R17).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 140 )*stride_in];
	(*R17).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 140 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 160 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 160 )*stride_in];
	(*R18).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 160 )*stride_in];
	(*R18).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 160 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 180 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 180 )*stride_in];
	(*R19).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 180 )*stride_in];
	(*R19).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 180 )*stride_in];
	}



	InvRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	InvRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);


	if(rw)
	{
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 4 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 5 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 6 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 7 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 8 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 9 ) ] = (*R19).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 20 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 20 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 40 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 40 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 60 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 60 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 80 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 80 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 100 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 100 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 120 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 120 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 140 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 140 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 160 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 160 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*2 + 0 + 180 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*2 + 1 + 180 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 4 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 5 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 6 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 7 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 8 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 9 ) ] = (*R19).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 20 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 20 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 40 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 40 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 60 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 60 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 80 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 80 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 100 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 100 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 120 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 120 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 140 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 140 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 160 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 160 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*2 + 0 + 180 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*2 + 1 + 180 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass1_len200(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{




	{
		T W = twiddles[9 + 9*((2*me + 0)%10) + 0];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 1];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 2];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 3];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 4];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 5];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 6];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 7];
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
		T W = twiddles[9 + 9*((2*me + 0)%10) + 8];
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
		T W = twiddles[9 + 9*((2*me + 1)%10) + 0];
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
		T W = twiddles[9 + 9*((2*me + 1)%10) + 1];
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
		T W = twiddles[9 + 9*((2*me + 1)%10) + 2];
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
		T W = twiddles[9 + 9*((2*me + 1)%10) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	{
		T W = twiddles[9 + 9*((2*me + 1)%10) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R15).x; ry = (*R15).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R15).x = TR;
		(*R15).y = TI;
	}

	{
		T W = twiddles[9 + 9*((2*me + 1)%10) + 5];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R16).x; ry = (*R16).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R16).x = TR;
		(*R16).y = TI;
	}

	{
		T W = twiddles[9 + 9*((2*me + 1)%10) + 6];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R17).x; ry = (*R17).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R17).x = TR;
		(*R17).y = TI;
	}

	{
		T W = twiddles[9 + 9*((2*me + 1)%10) + 7];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R18).x; ry = (*R18).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R18).x = TR;
		(*R18).y = TI;
	}

	{
		T W = twiddles[9 + 9*((2*me + 1)%10) + 8];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	InvRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	InvRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);


	if(rw)
	{
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 10 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 20 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 30 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 40 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 50 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 60 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 70 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 80 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 90 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 10 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 20 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 30 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 40 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 50 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 60 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 70 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 80 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 90 ) ] = (*R19).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*10 + 0 + 0 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*10 + 1 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*10 + 2 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*10 + 3 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*10 + 4 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*10 + 5 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*10 + 6 + 0 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*10 + 7 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*10 + 8 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*10 + 9 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*10 + 0 + 100 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*10 + 1 + 100 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*10 + 2 + 100 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*10 + 3 + 100 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*10 + 4 + 100 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*10 + 5 + 100 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*10 + 6 + 100 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*10 + 7 + 100 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*10 + 8 + 100 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*10 + 9 + 100 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 10 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 20 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 30 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 40 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 50 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 60 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 70 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 80 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 90 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 10 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 20 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 30 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 40 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 50 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 60 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 70 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 80 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 90 ) ] = (*R19).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*10 + 0 + 0 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*10 + 1 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*10 + 2 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*10 + 3 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*10 + 4 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*10 + 5 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*10 + 6 + 0 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*10 + 7 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*10 + 8 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*10 + 9 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*10 + 0 + 100 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*10 + 1 + 100 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*10 + 2 + 100 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*10 + 3 + 100 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*10 + 4 + 100 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*10 + 5 + 100 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*10 + 6 + 100 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*10 + 7 + 100 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*10 + 8 + 100 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*10 + 9 + 100 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass2_len200(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[99 + 1*((10*me + 0)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 1)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 2)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 3)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 4)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 5)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 6)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 7)%100) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R15).x; ry = (*R15).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R15).x = TR;
		(*R15).y = TI;
	}

	{
		T W = twiddles[99 + 1*((10*me + 8)%100) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R17).x; ry = (*R17).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R17).x = TR;
		(*R17).y = TI;
	}

	{
		T W = twiddles[99 + 1*((10*me + 9)%100) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	InvRad2B1(R0, R1);
	InvRad2B1(R2, R3);
	InvRad2B1(R4, R5);
	InvRad2B1(R6, R7);
	InvRad2B1(R8, R9);
	InvRad2B1(R10, R11);
	InvRad2B1(R12, R13);
	InvRad2B1(R14, R15);
	InvRad2B1(R16, R17);
	InvRad2B1(R18, R19);


	if(rw)
	{
	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);

	 //Optimization: coalescing into float4/double4 write
	if(sb == SB_UNIT && !store_cb_fn) {
	vector4_type_t<T> *buff4g = (vector4_type_t<T>*)(bufOut+outOffset);
	
	buff4g[ 5*me + 0 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R0).x, (*R0).y, (*R2).x, (*R2).y) ;
	buff4g[ 5*me + 1 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R4).x, (*R4).y, (*R6).x, (*R6).y) ;
	buff4g[ 5*me + 2 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R8).x, (*R8).y, (*R10).x, (*R10).y) ;
	buff4g[ 5*me + 3 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R12).x, (*R12).y, (*R14).x, (*R14).y) ;
	buff4g[ 5*me + 4 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R16).x, (*R16).y, (*R18).x, (*R18).y) ;
	buff4g[ 5*me + 0 + 50 ] = lib_make_vector4< vector4_type_t<T> >((*R1).x, (*R1).y, (*R3).x, (*R3).y) ;
	buff4g[ 5*me + 1 + 50 ] = lib_make_vector4< vector4_type_t<T> >((*R5).x, (*R5).y, (*R7).x, (*R7).y) ;
	buff4g[ 5*me + 2 + 50 ] = lib_make_vector4< vector4_type_t<T> >((*R9).x, (*R9).y, (*R11).x, (*R11).y) ;
	buff4g[ 5*me + 3 + 50 ] = lib_make_vector4< vector4_type_t<T> >((*R13).x, (*R13).y, (*R15).x, (*R15).y) ;
	buff4g[ 5*me + 4 + 50 ] = lib_make_vector4< vector4_type_t<T> >((*R17).x, (*R17).y, (*R19).x, (*R19).y) ;
	}
	else{ // such optimization is not possible 
	store_cb(bufOut, outOffset + ( 10*me + 0 + 0 )*stride_out, (*R0), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 1 + 0 )*stride_out, (*R2), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 2 + 0 )*stride_out, (*R4), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 3 + 0 )*stride_out, (*R6), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 4 + 0 )*stride_out, (*R8), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 5 + 0 )*stride_out, (*R10), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 6 + 0 )*stride_out, (*R12), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 7 + 0 )*stride_out, (*R14), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 8 + 0 )*stride_out, (*R16), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 9 + 0 )*stride_out, (*R18), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 0 + 100 )*stride_out, (*R1), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 1 + 100 )*stride_out, (*R3), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 2 + 100 )*stride_out, (*R5), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 3 + 100 )*stride_out, (*R7), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 4 + 100 )*stride_out, (*R9), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 5 + 100 )*stride_out, (*R11), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 6 + 100 )*stride_out, (*R13), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 7 + 100 )*stride_out, (*R15), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 8 + 100 )*stride_out, (*R17), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 9 + 100 )*stride_out, (*R19), store_cb_data, nullptr);

	}
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass2_len200(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[99 + 1*((10*me + 0)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 1)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 2)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 3)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 4)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 5)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 6)%100) + 0];
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
		T W = twiddles[99 + 1*((10*me + 7)%100) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R15).x; ry = (*R15).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R15).x = TR;
		(*R15).y = TI;
	}

	{
		T W = twiddles[99 + 1*((10*me + 8)%100) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R17).x; ry = (*R17).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R17).x = TR;
		(*R17).y = TI;
	}

	{
		T W = twiddles[99 + 1*((10*me + 9)%100) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	InvRad2B1(R0, R1);
	InvRad2B1(R2, R3);
	InvRad2B1(R4, R5);
	InvRad2B1(R6, R7);
	InvRad2B1(R8, R9);
	InvRad2B1(R10, R11);
	InvRad2B1(R12, R13);
	InvRad2B1(R14, R15);
	InvRad2B1(R16, R17);
	InvRad2B1(R18, R19);


	if(rw)
	{
	bufOutRe[outOffset + ( 10*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 10*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 10*me + 1 + 0 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 10*me + 1 + 0 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 10*me + 2 + 0 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 10*me + 2 + 0 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 10*me + 3 + 0 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 10*me + 3 + 0 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 10*me + 4 + 0 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 10*me + 4 + 0 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 10*me + 5 + 0 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 10*me + 5 + 0 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 10*me + 6 + 0 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 10*me + 6 + 0 )*stride_out] = (*R12).y;
	bufOutRe[outOffset + ( 10*me + 7 + 0 )*stride_out] = (*R14).x;
	bufOutIm[outOffset + ( 10*me + 7 + 0 )*stride_out] = (*R14).y;
	bufOutRe[outOffset + ( 10*me + 8 + 0 )*stride_out] = (*R16).x;
	bufOutIm[outOffset + ( 10*me + 8 + 0 )*stride_out] = (*R16).y;
	bufOutRe[outOffset + ( 10*me + 9 + 0 )*stride_out] = (*R18).x;
	bufOutIm[outOffset + ( 10*me + 9 + 0 )*stride_out] = (*R18).y;
	bufOutRe[outOffset + ( 10*me + 0 + 100 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 10*me + 0 + 100 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 10*me + 1 + 100 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 10*me + 1 + 100 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 10*me + 2 + 100 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 10*me + 2 + 100 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 10*me + 3 + 100 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 10*me + 3 + 100 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 10*me + 4 + 100 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 10*me + 4 + 100 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 10*me + 5 + 100 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 10*me + 5 + 100 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 10*me + 6 + 100 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 10*me + 6 + 100 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 10*me + 7 + 100 )*stride_out] = (*R15).x;
	bufOutIm[outOffset + ( 10*me + 7 + 100 )*stride_out] = (*R15).y;
	bufOutRe[outOffset + ( 10*me + 8 + 100 )*stride_out] = (*R17).x;
	bufOutIm[outOffset + ( 10*me + 8 + 100 )*stride_out] = (*R17).y;
	bufOutRe[outOffset + ( 10*me + 9 + 100 )*stride_out] = (*R19).x;
	bufOutIm[outOffset + ( 10*me + 9 + 100 )*stride_out] = (*R19).y;
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len200_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_fn, load_cb_data);
	FwdPass1_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass2_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len200_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_fn, load_cb_data);
	InvPass1_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass2_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len200_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_fn, load_cb_data);
	FwdPass1_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass2_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len200_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_fn, load_cb_data);
	InvPass1_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass2_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len200_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_fn, load_cb_data);
	FwdPass1_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass2_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len200_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_fn, load_cb_data);
	InvPass1_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass2_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len200_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_fn, load_cb_data);
	FwdPass1_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass2_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len200_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_fn, load_cb_data);
	InvPass1_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass2_len200<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 60, maximum transforms: 6, Passes: 3
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(60)
fft_fwd_ip_len200( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gb)
{

	__shared__ real_type_t<T> lds[1200];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int ioOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*6)*10) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*6 + (me/10));
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
	fwd_len200_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  rw, b, me%10, (me/10)*200, gb, ioOffset, gb, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 6, Passes: 3
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(60)
fft_back_ip_len200( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gb)
{

	__shared__ real_type_t<T> lds[1200];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int ioOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*6)*10) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*6 + (me/10));
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
	back_len200_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  rw, b, me%10, (me/10)*200, gb, ioOffset, gb, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 6, Passes: 3
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(60)
fft_fwd_ip_len200( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbRe, real_type_t<T> * __restrict__ gbIm)
{

	__shared__ real_type_t<T> lds[1200];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int ioOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*6)*10) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*6 + (me/10));
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
	fwd_len200_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  rw, b, me%10, (me/10)*200, gbRe, gbIm, ioOffset, gbRe, gbIm, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 6, Passes: 3
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(60)
fft_back_ip_len200( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbRe, real_type_t<T> * __restrict__ gbIm)
{

	__shared__ real_type_t<T> lds[1200];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int ioOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*6)*10) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*6 + (me/10));
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
	back_len200_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  rw, b, me%10, (me/10)*200, gbRe, gbIm, ioOffset, gbRe, gbIm, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 6, Passes: 3
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(60)
fft_fwd_op_len200( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[1200];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*6)*10) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*6 + (me/10));
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
	fwd_len200_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, (me/10)*200, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 6, Passes: 3
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(60)
fft_back_op_len200( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[1200];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*6)*10) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*6 + (me/10));
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
	back_len200_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, (me/10)*200, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 6, Passes: 3
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(60)
fft_fwd_op_len200( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[1200];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*6)*10) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*6 + (me/10));
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
	fwd_len200_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, (me/10)*200, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 6, Passes: 3
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(60)
fft_back_op_len200( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[1200];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*6)*10) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*6 + (me/10));
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
	back_len200_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, (me/10)*200, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 6, Passes: 3
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(60)
fft_fwd_op_len200( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[1200];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*6)*10) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*6 + (me/10));
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
	fwd_len200_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, (me/10)*200, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 6, Passes: 3
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(60)
fft_back_op_len200( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[1200];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*6)*10) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*6 + (me/10));
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
	back_len200_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, (me/10)*200, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 6, Passes: 3
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(60)
fft_fwd_op_len200( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[1200];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*6)*10) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*6 + (me/10));
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
	fwd_len200_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, (me/10)*200, gbInRe, gbInIm, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 6, Passes: 3
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(60)
fft_back_op_len200( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[1200];
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	unsigned int upper_count = batch_count;
	for(int i=1; i<dim; i++){
		upper_count *= lengths[i];
	}
	// do signed math to guard against underflow
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*6)*10) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	size_t counter_mod = (batch*6 + (me/10));
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
	back_len200_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, (me/10)*200, gbInRe, gbInIm, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

