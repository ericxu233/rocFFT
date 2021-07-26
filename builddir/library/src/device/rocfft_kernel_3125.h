#pragma once
#include "rocfft_butterfly_template.h"

////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass0_len3125(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	(*R0) = load_cb(bufIn, inOffset + ( 0 + me*5 + 0 + 0 )*stride_in, load_cb_data, nullptr);

	(*R5) = load_cb(bufIn, inOffset + ( 0 + me*5 + 1 + 0 )*stride_in, load_cb_data, nullptr);

	(*R10) = load_cb(bufIn, inOffset + ( 0 + me*5 + 2 + 0 )*stride_in, load_cb_data, nullptr);

	(*R15) = load_cb(bufIn, inOffset + ( 0 + me*5 + 3 + 0 )*stride_in, load_cb_data, nullptr);

	(*R20) = load_cb(bufIn, inOffset + ( 0 + me*5 + 4 + 0 )*stride_in, load_cb_data, nullptr);

	(*R1) = load_cb(bufIn, inOffset + ( 0 + me*5 + 0 + 625 )*stride_in, load_cb_data, nullptr);

	(*R6) = load_cb(bufIn, inOffset + ( 0 + me*5 + 1 + 625 )*stride_in, load_cb_data, nullptr);

	(*R11) = load_cb(bufIn, inOffset + ( 0 + me*5 + 2 + 625 )*stride_in, load_cb_data, nullptr);

	(*R16) = load_cb(bufIn, inOffset + ( 0 + me*5 + 3 + 625 )*stride_in, load_cb_data, nullptr);

	(*R21) = load_cb(bufIn, inOffset + ( 0 + me*5 + 4 + 625 )*stride_in, load_cb_data, nullptr);

	(*R2) = load_cb(bufIn, inOffset + ( 0 + me*5 + 0 + 1250 )*stride_in, load_cb_data, nullptr);

	(*R7) = load_cb(bufIn, inOffset + ( 0 + me*5 + 1 + 1250 )*stride_in, load_cb_data, nullptr);

	(*R12) = load_cb(bufIn, inOffset + ( 0 + me*5 + 2 + 1250 )*stride_in, load_cb_data, nullptr);

	(*R17) = load_cb(bufIn, inOffset + ( 0 + me*5 + 3 + 1250 )*stride_in, load_cb_data, nullptr);

	(*R22) = load_cb(bufIn, inOffset + ( 0 + me*5 + 4 + 1250 )*stride_in, load_cb_data, nullptr);

	(*R3) = load_cb(bufIn, inOffset + ( 0 + me*5 + 0 + 1875 )*stride_in, load_cb_data, nullptr);

	(*R8) = load_cb(bufIn, inOffset + ( 0 + me*5 + 1 + 1875 )*stride_in, load_cb_data, nullptr);

	(*R13) = load_cb(bufIn, inOffset + ( 0 + me*5 + 2 + 1875 )*stride_in, load_cb_data, nullptr);

	(*R18) = load_cb(bufIn, inOffset + ( 0 + me*5 + 3 + 1875 )*stride_in, load_cb_data, nullptr);

	(*R23) = load_cb(bufIn, inOffset + ( 0 + me*5 + 4 + 1875 )*stride_in, load_cb_data, nullptr);

	(*R4) = load_cb(bufIn, inOffset + ( 0 + me*5 + 0 + 2500 )*stride_in, load_cb_data, nullptr);

	(*R9) = load_cb(bufIn, inOffset + ( 0 + me*5 + 1 + 2500 )*stride_in, load_cb_data, nullptr);

	(*R14) = load_cb(bufIn, inOffset + ( 0 + me*5 + 2 + 2500 )*stride_in, load_cb_data, nullptr);

	(*R19) = load_cb(bufIn, inOffset + ( 0 + me*5 + 3 + 2500 )*stride_in, load_cb_data, nullptr);

	(*R24) = load_cb(bufIn, inOffset + ( 0 + me*5 + 4 + 2500 )*stride_in, load_cb_data, nullptr);

	}



	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);
	FwdRad5B1(R15, R16, R17, R18, R19);
	FwdRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 1 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 2 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 3 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 4 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 4 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 0 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 1 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 2 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 3 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 4 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 1 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 2 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 3 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 4 ) ] = (*R24).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 1 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 2 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 3 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 4 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 4 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 0 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 1 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 2 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 3 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 4 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 1 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 2 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 3 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 4 ) ] = (*R24).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass0_len3125(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*5 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*5 + 0 + 0 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*5 + 1 + 0 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*5 + 1 + 0 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*5 + 2 + 0 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*5 + 2 + 0 )*stride_in];
	(*R15).x = bufInRe[inOffset + ( 0 + me*5 + 3 + 0 )*stride_in];
	(*R15).y = bufInIm[inOffset + ( 0 + me*5 + 3 + 0 )*stride_in];
	(*R20).x = bufInRe[inOffset + ( 0 + me*5 + 4 + 0 )*stride_in];
	(*R20).y = bufInIm[inOffset + ( 0 + me*5 + 4 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*5 + 0 + 625 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*5 + 0 + 625 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*5 + 1 + 625 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*5 + 1 + 625 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*5 + 2 + 625 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*5 + 2 + 625 )*stride_in];
	(*R16).x = bufInRe[inOffset + ( 0 + me*5 + 3 + 625 )*stride_in];
	(*R16).y = bufInIm[inOffset + ( 0 + me*5 + 3 + 625 )*stride_in];
	(*R21).x = bufInRe[inOffset + ( 0 + me*5 + 4 + 625 )*stride_in];
	(*R21).y = bufInIm[inOffset + ( 0 + me*5 + 4 + 625 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*5 + 0 + 1250 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*5 + 0 + 1250 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*5 + 1 + 1250 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*5 + 1 + 1250 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*5 + 2 + 1250 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*5 + 2 + 1250 )*stride_in];
	(*R17).x = bufInRe[inOffset + ( 0 + me*5 + 3 + 1250 )*stride_in];
	(*R17).y = bufInIm[inOffset + ( 0 + me*5 + 3 + 1250 )*stride_in];
	(*R22).x = bufInRe[inOffset + ( 0 + me*5 + 4 + 1250 )*stride_in];
	(*R22).y = bufInIm[inOffset + ( 0 + me*5 + 4 + 1250 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*5 + 0 + 1875 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*5 + 0 + 1875 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*5 + 1 + 1875 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*5 + 1 + 1875 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*5 + 2 + 1875 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*5 + 2 + 1875 )*stride_in];
	(*R18).x = bufInRe[inOffset + ( 0 + me*5 + 3 + 1875 )*stride_in];
	(*R18).y = bufInIm[inOffset + ( 0 + me*5 + 3 + 1875 )*stride_in];
	(*R23).x = bufInRe[inOffset + ( 0 + me*5 + 4 + 1875 )*stride_in];
	(*R23).y = bufInIm[inOffset + ( 0 + me*5 + 4 + 1875 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*5 + 0 + 2500 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*5 + 0 + 2500 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*5 + 1 + 2500 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*5 + 1 + 2500 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*5 + 2 + 2500 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*5 + 2 + 2500 )*stride_in];
	(*R19).x = bufInRe[inOffset + ( 0 + me*5 + 3 + 2500 )*stride_in];
	(*R19).y = bufInIm[inOffset + ( 0 + me*5 + 3 + 2500 )*stride_in];
	(*R24).x = bufInRe[inOffset + ( 0 + me*5 + 4 + 2500 )*stride_in];
	(*R24).y = bufInIm[inOffset + ( 0 + me*5 + 4 + 2500 )*stride_in];
	}



	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);
	FwdRad5B1(R15, R16, R17, R18, R19);
	FwdRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 1 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 2 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 3 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 4 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 4 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 0 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 1 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 2 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 3 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 4 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 1 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 2 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 3 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 4 ) ] = (*R24).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 1 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 2 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 3 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 4 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 4 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 0 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 1 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 2 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 3 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 4 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 1 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 2 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 3 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 4 ) ] = (*R24).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass1_len3125(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24)
{




	{
		T W = twiddles[4 + 4*((5*me + 0)%5) + 0];
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
		T W = twiddles[4 + 4*((5*me + 0)%5) + 1];
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
		T W = twiddles[4 + 4*((5*me + 0)%5) + 2];
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
		T W = twiddles[4 + 4*((5*me + 0)%5) + 3];
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
		T W = twiddles[4 + 4*((5*me + 1)%5) + 0];
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
		T W = twiddles[4 + 4*((5*me + 1)%5) + 1];
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
		T W = twiddles[4 + 4*((5*me + 1)%5) + 2];
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
		T W = twiddles[4 + 4*((5*me + 1)%5) + 3];
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
		T W = twiddles[4 + 4*((5*me + 2)%5) + 0];
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
		T W = twiddles[4 + 4*((5*me + 2)%5) + 1];
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
		T W = twiddles[4 + 4*((5*me + 2)%5) + 2];
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
		T W = twiddles[4 + 4*((5*me + 2)%5) + 3];
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
		T W = twiddles[4 + 4*((5*me + 3)%5) + 0];
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
		T W = twiddles[4 + 4*((5*me + 3)%5) + 1];
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
		T W = twiddles[4 + 4*((5*me + 3)%5) + 2];
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
		T W = twiddles[4 + 4*((5*me + 3)%5) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	{
		T W = twiddles[4 + 4*((5*me + 4)%5) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R21).x; ry = (*R21).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R21).x = TR;
		(*R21).y = TI;
	}

	{
		T W = twiddles[4 + 4*((5*me + 4)%5) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R22).x; ry = (*R22).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R22).x = TR;
		(*R22).y = TI;
	}

	{
		T W = twiddles[4 + 4*((5*me + 4)%5) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R23).x; ry = (*R23).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R23).x = TR;
		(*R23).y = TI;
	}

	{
		T W = twiddles[4 + 4*((5*me + 4)%5) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);
	FwdRad5B1(R15, R16, R17, R18, R19);
	FwdRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 5 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 10 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 15 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 20 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 5 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 10 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 15 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 20 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 5 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 10 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 15 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 20 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 0 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 5 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 10 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 15 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 20 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 5 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 10 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 15 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 20 ) ] = (*R24).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 5 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 10 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 15 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 20 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 5 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 10 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 15 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 20 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 5 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 10 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 15 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 20 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 0 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 5 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 10 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 15 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 20 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 5 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 10 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 15 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 20 ) ] = (*R24).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass2_len3125(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24)
{




	{
		T W = twiddles[24 + 4*((5*me + 0)%25) + 0];
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
		T W = twiddles[24 + 4*((5*me + 0)%25) + 1];
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
		T W = twiddles[24 + 4*((5*me + 0)%25) + 2];
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
		T W = twiddles[24 + 4*((5*me + 0)%25) + 3];
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
		T W = twiddles[24 + 4*((5*me + 1)%25) + 0];
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
		T W = twiddles[24 + 4*((5*me + 1)%25) + 1];
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
		T W = twiddles[24 + 4*((5*me + 1)%25) + 2];
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
		T W = twiddles[24 + 4*((5*me + 1)%25) + 3];
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
		T W = twiddles[24 + 4*((5*me + 2)%25) + 0];
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
		T W = twiddles[24 + 4*((5*me + 2)%25) + 1];
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
		T W = twiddles[24 + 4*((5*me + 2)%25) + 2];
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
		T W = twiddles[24 + 4*((5*me + 2)%25) + 3];
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
		T W = twiddles[24 + 4*((5*me + 3)%25) + 0];
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
		T W = twiddles[24 + 4*((5*me + 3)%25) + 1];
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
		T W = twiddles[24 + 4*((5*me + 3)%25) + 2];
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
		T W = twiddles[24 + 4*((5*me + 3)%25) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	{
		T W = twiddles[24 + 4*((5*me + 4)%25) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R21).x; ry = (*R21).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R21).x = TR;
		(*R21).y = TI;
	}

	{
		T W = twiddles[24 + 4*((5*me + 4)%25) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R22).x; ry = (*R22).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R22).x = TR;
		(*R22).y = TI;
	}

	{
		T W = twiddles[24 + 4*((5*me + 4)%25) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R23).x; ry = (*R23).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R23).x = TR;
		(*R23).y = TI;
	}

	{
		T W = twiddles[24 + 4*((5*me + 4)%25) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);
	FwdRad5B1(R15, R16, R17, R18, R19);
	FwdRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 25 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 50 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 75 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 100 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 25 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 50 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 75 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 100 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 25 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 50 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 75 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 100 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 0 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 25 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 50 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 75 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 100 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 25 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 50 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 75 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 100 ) ] = (*R24).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 25 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 50 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 75 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 100 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 25 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 50 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 75 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 100 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 25 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 50 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 75 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 100 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 0 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 25 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 50 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 75 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 100 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 25 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 50 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 75 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 100 ) ] = (*R24).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass3_len3125(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24)
{




	{
		T W = twiddles[124 + 4*((5*me + 0)%125) + 0];
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
		T W = twiddles[124 + 4*((5*me + 0)%125) + 1];
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
		T W = twiddles[124 + 4*((5*me + 0)%125) + 2];
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
		T W = twiddles[124 + 4*((5*me + 0)%125) + 3];
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
		T W = twiddles[124 + 4*((5*me + 1)%125) + 0];
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
		T W = twiddles[124 + 4*((5*me + 1)%125) + 1];
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
		T W = twiddles[124 + 4*((5*me + 1)%125) + 2];
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
		T W = twiddles[124 + 4*((5*me + 1)%125) + 3];
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
		T W = twiddles[124 + 4*((5*me + 2)%125) + 0];
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
		T W = twiddles[124 + 4*((5*me + 2)%125) + 1];
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
		T W = twiddles[124 + 4*((5*me + 2)%125) + 2];
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
		T W = twiddles[124 + 4*((5*me + 2)%125) + 3];
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
		T W = twiddles[124 + 4*((5*me + 3)%125) + 0];
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
		T W = twiddles[124 + 4*((5*me + 3)%125) + 1];
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
		T W = twiddles[124 + 4*((5*me + 3)%125) + 2];
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
		T W = twiddles[124 + 4*((5*me + 3)%125) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	{
		T W = twiddles[124 + 4*((5*me + 4)%125) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R21).x; ry = (*R21).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R21).x = TR;
		(*R21).y = TI;
	}

	{
		T W = twiddles[124 + 4*((5*me + 4)%125) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R22).x; ry = (*R22).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R22).x = TR;
		(*R22).y = TI;
	}

	{
		T W = twiddles[124 + 4*((5*me + 4)%125) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R23).x; ry = (*R23).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R23).x = TR;
		(*R23).y = TI;
	}

	{
		T W = twiddles[124 + 4*((5*me + 4)%125) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);
	FwdRad5B1(R15, R16, R17, R18, R19);
	FwdRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 125 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 250 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 375 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 500 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 125 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 250 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 375 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 500 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 125 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 250 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 375 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 500 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 0 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 125 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 250 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 375 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 500 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 125 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 250 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 375 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 500 ) ] = (*R24).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 125 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 250 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 375 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 500 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 125 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 250 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 375 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 500 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 125 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 250 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 375 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 500 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 0 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 125 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 250 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 375 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 500 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 125 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 250 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 375 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 500 ) ] = (*R24).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass4_len3125(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[624 + 4*((5*me + 0)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 0)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 0)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 0)%625) + 3];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 3];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 3];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R21).x; ry = (*R21).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R21).x = TR;
		(*R21).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R22).x; ry = (*R22).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R22).x = TR;
		(*R22).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R23).x; ry = (*R23).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R23).x = TR;
		(*R23).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);
	FwdRad5B1(R15, R16, R17, R18, R19);
	FwdRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);

	store_cb(bufOut, outOffset + ( 5*me + 0 + 0 )*stride_out, (*R0), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 1 + 0 )*stride_out, (*R5), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 2 + 0 )*stride_out, (*R10), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 3 + 0 )*stride_out, (*R15), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 4 + 0 )*stride_out, (*R20), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 0 + 625 )*stride_out, (*R1), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 1 + 625 )*stride_out, (*R6), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 2 + 625 )*stride_out, (*R11), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 3 + 625 )*stride_out, (*R16), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 4 + 625 )*stride_out, (*R21), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 0 + 1250 )*stride_out, (*R2), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 1 + 1250 )*stride_out, (*R7), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 2 + 1250 )*stride_out, (*R12), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 3 + 1250 )*stride_out, (*R17), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 4 + 1250 )*stride_out, (*R22), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 0 + 1875 )*stride_out, (*R3), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 1 + 1875 )*stride_out, (*R8), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 2 + 1875 )*stride_out, (*R13), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 3 + 1875 )*stride_out, (*R18), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 4 + 1875 )*stride_out, (*R23), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 0 + 2500 )*stride_out, (*R4), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 1 + 2500 )*stride_out, (*R9), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 2 + 2500 )*stride_out, (*R14), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 3 + 2500 )*stride_out, (*R19), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 4 + 2500 )*stride_out, (*R24), store_cb_data, nullptr);

	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass4_len3125(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[624 + 4*((5*me + 0)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 0)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 0)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 0)%625) + 3];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 3];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 3];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R21).x; ry = (*R21).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R21).x = TR;
		(*R21).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R22).x; ry = (*R22).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R22).x = TR;
		(*R22).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R23).x; ry = (*R23).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R23).x = TR;
		(*R23).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	FwdRad5B1(R0, R1, R2, R3, R4);
	FwdRad5B1(R5, R6, R7, R8, R9);
	FwdRad5B1(R10, R11, R12, R13, R14);
	FwdRad5B1(R15, R16, R17, R18, R19);
	FwdRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R15).x;
	bufOutIm[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R15).y;
	bufOutRe[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R20).x;
	bufOutIm[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R20).y;
	bufOutRe[outOffset + ( 5*me + 0 + 625 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 5*me + 0 + 625 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 5*me + 1 + 625 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 5*me + 1 + 625 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 5*me + 2 + 625 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 5*me + 2 + 625 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 5*me + 3 + 625 )*stride_out] = (*R16).x;
	bufOutIm[outOffset + ( 5*me + 3 + 625 )*stride_out] = (*R16).y;
	bufOutRe[outOffset + ( 5*me + 4 + 625 )*stride_out] = (*R21).x;
	bufOutIm[outOffset + ( 5*me + 4 + 625 )*stride_out] = (*R21).y;
	bufOutRe[outOffset + ( 5*me + 0 + 1250 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 5*me + 0 + 1250 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 5*me + 1 + 1250 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 5*me + 1 + 1250 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 5*me + 2 + 1250 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 5*me + 2 + 1250 )*stride_out] = (*R12).y;
	bufOutRe[outOffset + ( 5*me + 3 + 1250 )*stride_out] = (*R17).x;
	bufOutIm[outOffset + ( 5*me + 3 + 1250 )*stride_out] = (*R17).y;
	bufOutRe[outOffset + ( 5*me + 4 + 1250 )*stride_out] = (*R22).x;
	bufOutIm[outOffset + ( 5*me + 4 + 1250 )*stride_out] = (*R22).y;
	bufOutRe[outOffset + ( 5*me + 0 + 1875 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 5*me + 0 + 1875 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 5*me + 1 + 1875 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 5*me + 1 + 1875 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 5*me + 2 + 1875 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 5*me + 2 + 1875 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 5*me + 3 + 1875 )*stride_out] = (*R18).x;
	bufOutIm[outOffset + ( 5*me + 3 + 1875 )*stride_out] = (*R18).y;
	bufOutRe[outOffset + ( 5*me + 4 + 1875 )*stride_out] = (*R23).x;
	bufOutIm[outOffset + ( 5*me + 4 + 1875 )*stride_out] = (*R23).y;
	bufOutRe[outOffset + ( 5*me + 0 + 2500 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 5*me + 0 + 2500 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 5*me + 1 + 2500 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 5*me + 1 + 2500 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 5*me + 2 + 2500 )*stride_out] = (*R14).x;
	bufOutIm[outOffset + ( 5*me + 2 + 2500 )*stride_out] = (*R14).y;
	bufOutRe[outOffset + ( 5*me + 3 + 2500 )*stride_out] = (*R19).x;
	bufOutIm[outOffset + ( 5*me + 3 + 2500 )*stride_out] = (*R19).y;
	bufOutRe[outOffset + ( 5*me + 4 + 2500 )*stride_out] = (*R24).x;
	bufOutIm[outOffset + ( 5*me + 4 + 2500 )*stride_out] = (*R24).y;
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass0_len3125(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	(*R0) = load_cb(bufIn, inOffset + ( 0 + me*5 + 0 + 0 )*stride_in, load_cb_data, nullptr);

	(*R5) = load_cb(bufIn, inOffset + ( 0 + me*5 + 1 + 0 )*stride_in, load_cb_data, nullptr);

	(*R10) = load_cb(bufIn, inOffset + ( 0 + me*5 + 2 + 0 )*stride_in, load_cb_data, nullptr);

	(*R15) = load_cb(bufIn, inOffset + ( 0 + me*5 + 3 + 0 )*stride_in, load_cb_data, nullptr);

	(*R20) = load_cb(bufIn, inOffset + ( 0 + me*5 + 4 + 0 )*stride_in, load_cb_data, nullptr);

	(*R1) = load_cb(bufIn, inOffset + ( 0 + me*5 + 0 + 625 )*stride_in, load_cb_data, nullptr);

	(*R6) = load_cb(bufIn, inOffset + ( 0 + me*5 + 1 + 625 )*stride_in, load_cb_data, nullptr);

	(*R11) = load_cb(bufIn, inOffset + ( 0 + me*5 + 2 + 625 )*stride_in, load_cb_data, nullptr);

	(*R16) = load_cb(bufIn, inOffset + ( 0 + me*5 + 3 + 625 )*stride_in, load_cb_data, nullptr);

	(*R21) = load_cb(bufIn, inOffset + ( 0 + me*5 + 4 + 625 )*stride_in, load_cb_data, nullptr);

	(*R2) = load_cb(bufIn, inOffset + ( 0 + me*5 + 0 + 1250 )*stride_in, load_cb_data, nullptr);

	(*R7) = load_cb(bufIn, inOffset + ( 0 + me*5 + 1 + 1250 )*stride_in, load_cb_data, nullptr);

	(*R12) = load_cb(bufIn, inOffset + ( 0 + me*5 + 2 + 1250 )*stride_in, load_cb_data, nullptr);

	(*R17) = load_cb(bufIn, inOffset + ( 0 + me*5 + 3 + 1250 )*stride_in, load_cb_data, nullptr);

	(*R22) = load_cb(bufIn, inOffset + ( 0 + me*5 + 4 + 1250 )*stride_in, load_cb_data, nullptr);

	(*R3) = load_cb(bufIn, inOffset + ( 0 + me*5 + 0 + 1875 )*stride_in, load_cb_data, nullptr);

	(*R8) = load_cb(bufIn, inOffset + ( 0 + me*5 + 1 + 1875 )*stride_in, load_cb_data, nullptr);

	(*R13) = load_cb(bufIn, inOffset + ( 0 + me*5 + 2 + 1875 )*stride_in, load_cb_data, nullptr);

	(*R18) = load_cb(bufIn, inOffset + ( 0 + me*5 + 3 + 1875 )*stride_in, load_cb_data, nullptr);

	(*R23) = load_cb(bufIn, inOffset + ( 0 + me*5 + 4 + 1875 )*stride_in, load_cb_data, nullptr);

	(*R4) = load_cb(bufIn, inOffset + ( 0 + me*5 + 0 + 2500 )*stride_in, load_cb_data, nullptr);

	(*R9) = load_cb(bufIn, inOffset + ( 0 + me*5 + 1 + 2500 )*stride_in, load_cb_data, nullptr);

	(*R14) = load_cb(bufIn, inOffset + ( 0 + me*5 + 2 + 2500 )*stride_in, load_cb_data, nullptr);

	(*R19) = load_cb(bufIn, inOffset + ( 0 + me*5 + 3 + 2500 )*stride_in, load_cb_data, nullptr);

	(*R24) = load_cb(bufIn, inOffset + ( 0 + me*5 + 4 + 2500 )*stride_in, load_cb_data, nullptr);

	}



	InvRad5B1(R0, R1, R2, R3, R4);
	InvRad5B1(R5, R6, R7, R8, R9);
	InvRad5B1(R10, R11, R12, R13, R14);
	InvRad5B1(R15, R16, R17, R18, R19);
	InvRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 1 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 2 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 3 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 4 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 4 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 0 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 1 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 2 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 3 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 4 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 1 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 2 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 3 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 4 ) ] = (*R24).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 1 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 2 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 3 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 4 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 4 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 0 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 1 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 2 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 3 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 4 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 1 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 2 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 3 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 4 ) ] = (*R24).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass0_len3125(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*5 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*5 + 0 + 0 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*5 + 1 + 0 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*5 + 1 + 0 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*5 + 2 + 0 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*5 + 2 + 0 )*stride_in];
	(*R15).x = bufInRe[inOffset + ( 0 + me*5 + 3 + 0 )*stride_in];
	(*R15).y = bufInIm[inOffset + ( 0 + me*5 + 3 + 0 )*stride_in];
	(*R20).x = bufInRe[inOffset + ( 0 + me*5 + 4 + 0 )*stride_in];
	(*R20).y = bufInIm[inOffset + ( 0 + me*5 + 4 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*5 + 0 + 625 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*5 + 0 + 625 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*5 + 1 + 625 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*5 + 1 + 625 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*5 + 2 + 625 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*5 + 2 + 625 )*stride_in];
	(*R16).x = bufInRe[inOffset + ( 0 + me*5 + 3 + 625 )*stride_in];
	(*R16).y = bufInIm[inOffset + ( 0 + me*5 + 3 + 625 )*stride_in];
	(*R21).x = bufInRe[inOffset + ( 0 + me*5 + 4 + 625 )*stride_in];
	(*R21).y = bufInIm[inOffset + ( 0 + me*5 + 4 + 625 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*5 + 0 + 1250 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*5 + 0 + 1250 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*5 + 1 + 1250 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*5 + 1 + 1250 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*5 + 2 + 1250 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*5 + 2 + 1250 )*stride_in];
	(*R17).x = bufInRe[inOffset + ( 0 + me*5 + 3 + 1250 )*stride_in];
	(*R17).y = bufInIm[inOffset + ( 0 + me*5 + 3 + 1250 )*stride_in];
	(*R22).x = bufInRe[inOffset + ( 0 + me*5 + 4 + 1250 )*stride_in];
	(*R22).y = bufInIm[inOffset + ( 0 + me*5 + 4 + 1250 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*5 + 0 + 1875 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*5 + 0 + 1875 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*5 + 1 + 1875 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*5 + 1 + 1875 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*5 + 2 + 1875 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*5 + 2 + 1875 )*stride_in];
	(*R18).x = bufInRe[inOffset + ( 0 + me*5 + 3 + 1875 )*stride_in];
	(*R18).y = bufInIm[inOffset + ( 0 + me*5 + 3 + 1875 )*stride_in];
	(*R23).x = bufInRe[inOffset + ( 0 + me*5 + 4 + 1875 )*stride_in];
	(*R23).y = bufInIm[inOffset + ( 0 + me*5 + 4 + 1875 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*5 + 0 + 2500 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*5 + 0 + 2500 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*5 + 1 + 2500 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*5 + 1 + 2500 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*5 + 2 + 2500 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*5 + 2 + 2500 )*stride_in];
	(*R19).x = bufInRe[inOffset + ( 0 + me*5 + 3 + 2500 )*stride_in];
	(*R19).y = bufInIm[inOffset + ( 0 + me*5 + 3 + 2500 )*stride_in];
	(*R24).x = bufInRe[inOffset + ( 0 + me*5 + 4 + 2500 )*stride_in];
	(*R24).y = bufInIm[inOffset + ( 0 + me*5 + 4 + 2500 )*stride_in];
	}



	InvRad5B1(R0, R1, R2, R3, R4);
	InvRad5B1(R5, R6, R7, R8, R9);
	InvRad5B1(R10, R11, R12, R13, R14);
	InvRad5B1(R15, R16, R17, R18, R19);
	InvRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 1 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 2 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 3 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 4 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 4 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 0 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 1 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 2 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 3 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 4 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 1 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 2 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 3 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 4 ) ] = (*R24).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/1)*5 + (5*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 1 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 2 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 3 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/1)*5 + (5*me + 1)%1 + 4 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/1)*5 + (5*me + 2)%1 + 4 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 0 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 1 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 2 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 3 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/1)*5 + (5*me + 3)%1 + 4 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 1 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 2 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 3 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/1)*5 + (5*me + 4)%1 + 4 ) ] = (*R24).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass1_len3125(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24)
{




	{
		T W = twiddles[4 + 4*((5*me + 0)%5) + 0];
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
		T W = twiddles[4 + 4*((5*me + 0)%5) + 1];
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
		T W = twiddles[4 + 4*((5*me + 0)%5) + 2];
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
		T W = twiddles[4 + 4*((5*me + 0)%5) + 3];
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
		T W = twiddles[4 + 4*((5*me + 1)%5) + 0];
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
		T W = twiddles[4 + 4*((5*me + 1)%5) + 1];
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
		T W = twiddles[4 + 4*((5*me + 1)%5) + 2];
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
		T W = twiddles[4 + 4*((5*me + 1)%5) + 3];
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
		T W = twiddles[4 + 4*((5*me + 2)%5) + 0];
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
		T W = twiddles[4 + 4*((5*me + 2)%5) + 1];
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
		T W = twiddles[4 + 4*((5*me + 2)%5) + 2];
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
		T W = twiddles[4 + 4*((5*me + 2)%5) + 3];
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
		T W = twiddles[4 + 4*((5*me + 3)%5) + 0];
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
		T W = twiddles[4 + 4*((5*me + 3)%5) + 1];
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
		T W = twiddles[4 + 4*((5*me + 3)%5) + 2];
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
		T W = twiddles[4 + 4*((5*me + 3)%5) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	{
		T W = twiddles[4 + 4*((5*me + 4)%5) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R21).x; ry = (*R21).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R21).x = TR;
		(*R21).y = TI;
	}

	{
		T W = twiddles[4 + 4*((5*me + 4)%5) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R22).x; ry = (*R22).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R22).x = TR;
		(*R22).y = TI;
	}

	{
		T W = twiddles[4 + 4*((5*me + 4)%5) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R23).x; ry = (*R23).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R23).x = TR;
		(*R23).y = TI;
	}

	{
		T W = twiddles[4 + 4*((5*me + 4)%5) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	InvRad5B1(R0, R1, R2, R3, R4);
	InvRad5B1(R5, R6, R7, R8, R9);
	InvRad5B1(R10, R11, R12, R13, R14);
	InvRad5B1(R15, R16, R17, R18, R19);
	InvRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 5 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 10 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 15 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 20 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 5 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 10 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 15 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 20 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 5 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 10 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 15 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 20 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 0 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 5 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 10 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 15 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 20 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 5 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 10 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 15 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 20 ) ] = (*R24).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 5 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 10 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 15 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/5)*25 + (5*me + 0)%5 + 20 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 5 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 10 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 15 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/5)*25 + (5*me + 1)%5 + 20 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 5 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 10 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 15 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/5)*25 + (5*me + 2)%5 + 20 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 0 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 5 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 10 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 15 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/5)*25 + (5*me + 3)%5 + 20 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 5 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 10 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 15 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/5)*25 + (5*me + 4)%5 + 20 ) ] = (*R24).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass2_len3125(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24)
{




	{
		T W = twiddles[24 + 4*((5*me + 0)%25) + 0];
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
		T W = twiddles[24 + 4*((5*me + 0)%25) + 1];
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
		T W = twiddles[24 + 4*((5*me + 0)%25) + 2];
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
		T W = twiddles[24 + 4*((5*me + 0)%25) + 3];
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
		T W = twiddles[24 + 4*((5*me + 1)%25) + 0];
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
		T W = twiddles[24 + 4*((5*me + 1)%25) + 1];
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
		T W = twiddles[24 + 4*((5*me + 1)%25) + 2];
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
		T W = twiddles[24 + 4*((5*me + 1)%25) + 3];
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
		T W = twiddles[24 + 4*((5*me + 2)%25) + 0];
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
		T W = twiddles[24 + 4*((5*me + 2)%25) + 1];
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
		T W = twiddles[24 + 4*((5*me + 2)%25) + 2];
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
		T W = twiddles[24 + 4*((5*me + 2)%25) + 3];
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
		T W = twiddles[24 + 4*((5*me + 3)%25) + 0];
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
		T W = twiddles[24 + 4*((5*me + 3)%25) + 1];
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
		T W = twiddles[24 + 4*((5*me + 3)%25) + 2];
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
		T W = twiddles[24 + 4*((5*me + 3)%25) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	{
		T W = twiddles[24 + 4*((5*me + 4)%25) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R21).x; ry = (*R21).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R21).x = TR;
		(*R21).y = TI;
	}

	{
		T W = twiddles[24 + 4*((5*me + 4)%25) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R22).x; ry = (*R22).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R22).x = TR;
		(*R22).y = TI;
	}

	{
		T W = twiddles[24 + 4*((5*me + 4)%25) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R23).x; ry = (*R23).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R23).x = TR;
		(*R23).y = TI;
	}

	{
		T W = twiddles[24 + 4*((5*me + 4)%25) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	InvRad5B1(R0, R1, R2, R3, R4);
	InvRad5B1(R5, R6, R7, R8, R9);
	InvRad5B1(R10, R11, R12, R13, R14);
	InvRad5B1(R15, R16, R17, R18, R19);
	InvRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 25 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 50 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 75 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 100 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 25 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 50 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 75 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 100 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 25 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 50 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 75 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 100 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 0 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 25 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 50 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 75 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 100 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 25 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 50 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 75 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 100 ) ] = (*R24).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 25 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 50 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 75 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/25)*125 + (5*me + 0)%25 + 100 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 25 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 50 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 75 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/25)*125 + (5*me + 1)%25 + 100 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 25 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 50 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 75 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/25)*125 + (5*me + 2)%25 + 100 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 0 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 25 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 50 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 75 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/25)*125 + (5*me + 3)%25 + 100 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 25 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 50 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 75 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/25)*125 + (5*me + 4)%25 + 100 ) ] = (*R24).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass3_len3125(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24)
{




	{
		T W = twiddles[124 + 4*((5*me + 0)%125) + 0];
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
		T W = twiddles[124 + 4*((5*me + 0)%125) + 1];
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
		T W = twiddles[124 + 4*((5*me + 0)%125) + 2];
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
		T W = twiddles[124 + 4*((5*me + 0)%125) + 3];
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
		T W = twiddles[124 + 4*((5*me + 1)%125) + 0];
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
		T W = twiddles[124 + 4*((5*me + 1)%125) + 1];
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
		T W = twiddles[124 + 4*((5*me + 1)%125) + 2];
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
		T W = twiddles[124 + 4*((5*me + 1)%125) + 3];
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
		T W = twiddles[124 + 4*((5*me + 2)%125) + 0];
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
		T W = twiddles[124 + 4*((5*me + 2)%125) + 1];
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
		T W = twiddles[124 + 4*((5*me + 2)%125) + 2];
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
		T W = twiddles[124 + 4*((5*me + 2)%125) + 3];
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
		T W = twiddles[124 + 4*((5*me + 3)%125) + 0];
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
		T W = twiddles[124 + 4*((5*me + 3)%125) + 1];
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
		T W = twiddles[124 + 4*((5*me + 3)%125) + 2];
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
		T W = twiddles[124 + 4*((5*me + 3)%125) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	{
		T W = twiddles[124 + 4*((5*me + 4)%125) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R21).x; ry = (*R21).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R21).x = TR;
		(*R21).y = TI;
	}

	{
		T W = twiddles[124 + 4*((5*me + 4)%125) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R22).x; ry = (*R22).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R22).x = TR;
		(*R22).y = TI;
	}

	{
		T W = twiddles[124 + 4*((5*me + 4)%125) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R23).x; ry = (*R23).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R23).x = TR;
		(*R23).y = TI;
	}

	{
		T W = twiddles[124 + 4*((5*me + 4)%125) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	InvRad5B1(R0, R1, R2, R3, R4);
	InvRad5B1(R5, R6, R7, R8, R9);
	InvRad5B1(R10, R11, R12, R13, R14);
	InvRad5B1(R15, R16, R17, R18, R19);
	InvRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 125 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 250 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 375 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 500 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 0 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 125 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 250 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 375 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 500 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 125 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 250 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 375 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 500 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 0 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 125 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 250 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 375 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 500 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 125 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 250 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 375 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 500 ) ] = (*R24).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 125 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 250 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 375 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/125)*625 + (5*me + 0)%125 + 500 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 0 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 125 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 250 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 375 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/125)*625 + (5*me + 1)%125 + 500 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 125 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 250 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 375 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/125)*625 + (5*me + 2)%125 + 500 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 0 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 125 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 250 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 375 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/125)*625 + (5*me + 3)%125 + 500 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 125 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 250 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 375 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/125)*625 + (5*me + 4)%125 + 500 ) ] = (*R24).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 625 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 625 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 625 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 625 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 625 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1250 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1250 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1250 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1250 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1250 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1875 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1875 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1875 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1875 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1875 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2500 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2500 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2500 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2500 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2500 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass4_len3125(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[624 + 4*((5*me + 0)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 0)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 0)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 0)%625) + 3];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 3];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 3];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R21).x; ry = (*R21).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R21).x = TR;
		(*R21).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R22).x; ry = (*R22).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R22).x = TR;
		(*R22).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R23).x; ry = (*R23).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R23).x = TR;
		(*R23).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	InvRad5B1(R0, R1, R2, R3, R4);
	InvRad5B1(R5, R6, R7, R8, R9);
	InvRad5B1(R10, R11, R12, R13, R14);
	InvRad5B1(R15, R16, R17, R18, R19);
	InvRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);

	store_cb(bufOut, outOffset + ( 5*me + 0 + 0 )*stride_out, (*R0), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 1 + 0 )*stride_out, (*R5), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 2 + 0 )*stride_out, (*R10), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 3 + 0 )*stride_out, (*R15), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 4 + 0 )*stride_out, (*R20), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 0 + 625 )*stride_out, (*R1), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 1 + 625 )*stride_out, (*R6), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 2 + 625 )*stride_out, (*R11), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 3 + 625 )*stride_out, (*R16), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 4 + 625 )*stride_out, (*R21), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 0 + 1250 )*stride_out, (*R2), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 1 + 1250 )*stride_out, (*R7), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 2 + 1250 )*stride_out, (*R12), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 3 + 1250 )*stride_out, (*R17), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 4 + 1250 )*stride_out, (*R22), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 0 + 1875 )*stride_out, (*R3), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 1 + 1875 )*stride_out, (*R8), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 2 + 1875 )*stride_out, (*R13), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 3 + 1875 )*stride_out, (*R18), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 4 + 1875 )*stride_out, (*R23), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 0 + 2500 )*stride_out, (*R4), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 1 + 2500 )*stride_out, (*R9), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 2 + 2500 )*stride_out, (*R14), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 3 + 2500 )*stride_out, (*R19), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 5*me + 4 + 2500 )*stride_out, (*R24), store_cb_data, nullptr);

	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass4_len3125(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[624 + 4*((5*me + 0)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 0)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 0)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 0)%625) + 3];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 1)%625) + 3];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 2)%625) + 3];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 0];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 1];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 2];
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
		T W = twiddles[624 + 4*((5*me + 3)%625) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R21).x; ry = (*R21).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R21).x = TR;
		(*R21).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R22).x; ry = (*R22).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R22).x = TR;
		(*R22).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R23).x; ry = (*R23).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R23).x = TR;
		(*R23).y = TI;
	}

	{
		T W = twiddles[624 + 4*((5*me + 4)%625) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	InvRad5B1(R0, R1, R2, R3, R4);
	InvRad5B1(R5, R6, R7, R8, R9);
	InvRad5B1(R10, R11, R12, R13, R14);
	InvRad5B1(R15, R16, R17, R18, R19);
	InvRad5B1(R20, R21, R22, R23, R24);


	if(rw)
	{
	bufOutRe[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R15).x;
	bufOutIm[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R15).y;
	bufOutRe[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R20).x;
	bufOutIm[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R20).y;
	bufOutRe[outOffset + ( 5*me + 0 + 625 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 5*me + 0 + 625 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 5*me + 1 + 625 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 5*me + 1 + 625 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 5*me + 2 + 625 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 5*me + 2 + 625 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 5*me + 3 + 625 )*stride_out] = (*R16).x;
	bufOutIm[outOffset + ( 5*me + 3 + 625 )*stride_out] = (*R16).y;
	bufOutRe[outOffset + ( 5*me + 4 + 625 )*stride_out] = (*R21).x;
	bufOutIm[outOffset + ( 5*me + 4 + 625 )*stride_out] = (*R21).y;
	bufOutRe[outOffset + ( 5*me + 0 + 1250 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 5*me + 0 + 1250 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 5*me + 1 + 1250 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 5*me + 1 + 1250 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 5*me + 2 + 1250 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 5*me + 2 + 1250 )*stride_out] = (*R12).y;
	bufOutRe[outOffset + ( 5*me + 3 + 1250 )*stride_out] = (*R17).x;
	bufOutIm[outOffset + ( 5*me + 3 + 1250 )*stride_out] = (*R17).y;
	bufOutRe[outOffset + ( 5*me + 4 + 1250 )*stride_out] = (*R22).x;
	bufOutIm[outOffset + ( 5*me + 4 + 1250 )*stride_out] = (*R22).y;
	bufOutRe[outOffset + ( 5*me + 0 + 1875 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 5*me + 0 + 1875 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 5*me + 1 + 1875 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 5*me + 1 + 1875 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 5*me + 2 + 1875 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 5*me + 2 + 1875 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 5*me + 3 + 1875 )*stride_out] = (*R18).x;
	bufOutIm[outOffset + ( 5*me + 3 + 1875 )*stride_out] = (*R18).y;
	bufOutRe[outOffset + ( 5*me + 4 + 1875 )*stride_out] = (*R23).x;
	bufOutIm[outOffset + ( 5*me + 4 + 1875 )*stride_out] = (*R23).y;
	bufOutRe[outOffset + ( 5*me + 0 + 2500 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 5*me + 0 + 2500 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 5*me + 1 + 2500 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 5*me + 1 + 2500 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 5*me + 2 + 2500 )*stride_out] = (*R14).x;
	bufOutIm[outOffset + ( 5*me + 2 + 2500 )*stride_out] = (*R14).y;
	bufOutRe[outOffset + ( 5*me + 3 + 2500 )*stride_out] = (*R19).x;
	bufOutIm[outOffset + ( 5*me + 3 + 2500 )*stride_out] = (*R19).y;
	bufOutRe[outOffset + ( 5*me + 4 + 2500 )*stride_out] = (*R24).x;
	bufOutIm[outOffset + ( 5*me + 4 + 2500 )*stride_out] = (*R24).y;
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len3125_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24;
	FwdPass0_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_fn, load_cb_data);
	FwdPass1_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass2_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass3_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass4_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len3125_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24;
	InvPass0_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_fn, load_cb_data);
	InvPass1_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	InvPass2_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	InvPass3_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	InvPass4_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len3125_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24;
	FwdPass0_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_fn, load_cb_data);
	FwdPass1_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass2_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass3_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass4_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len3125_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24;
	InvPass0_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_fn, load_cb_data);
	InvPass1_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	InvPass2_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	InvPass3_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	InvPass4_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len3125_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24;
	FwdPass0_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_fn, load_cb_data);
	FwdPass1_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass2_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass3_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass4_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len3125_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24;
	InvPass0_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_fn, load_cb_data);
	InvPass1_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	InvPass2_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	InvPass3_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	InvPass4_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len3125_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24;
	FwdPass0_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_fn, load_cb_data);
	FwdPass1_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass2_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass3_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	FwdPass4_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len3125_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24;
	InvPass0_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_fn, load_cb_data);
	InvPass1_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	InvPass2_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	InvPass3_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24);
	InvPass4_len3125<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 125, maximum transforms: 1, Passes: 5
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(125)
fft_fwd_ip_len3125( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gb)
{

	__shared__ real_type_t<T> lds[3125];
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
	fwd_len3125_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, gb, ioOffset, gb, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 125, maximum transforms: 1, Passes: 5
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(125)
fft_back_ip_len3125( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gb)
{

	__shared__ real_type_t<T> lds[3125];
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
	back_len3125_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, gb, ioOffset, gb, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 125, maximum transforms: 1, Passes: 5
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(125)
fft_fwd_ip_len3125( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbRe, real_type_t<T> * __restrict__ gbIm)
{

	__shared__ real_type_t<T> lds[3125];
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
	fwd_len3125_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, gbRe, gbIm, ioOffset, gbRe, gbIm, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 125, maximum transforms: 1, Passes: 5
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(125)
fft_back_ip_len3125( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbRe, real_type_t<T> * __restrict__ gbIm)
{

	__shared__ real_type_t<T> lds[3125];
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
	back_len3125_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, gbRe, gbIm, ioOffset, gbRe, gbIm, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 125, maximum transforms: 1, Passes: 5
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(125)
fft_fwd_op_len3125( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[3125];
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
	fwd_len3125_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 125, maximum transforms: 1, Passes: 5
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(125)
fft_back_op_len3125( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[3125];
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
	back_len3125_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 125, maximum transforms: 1, Passes: 5
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(125)
fft_fwd_op_len3125( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[3125];
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
	fwd_len3125_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 125, maximum transforms: 1, Passes: 5
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(125)
fft_back_op_len3125( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[3125];
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
	back_len3125_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 125, maximum transforms: 1, Passes: 5
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(125)
fft_fwd_op_len3125( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[3125];
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
	fwd_len3125_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 125, maximum transforms: 1, Passes: 5
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(125)
fft_back_op_len3125( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[3125];
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
	back_len3125_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 125, maximum transforms: 1, Passes: 5
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(125)
fft_fwd_op_len3125( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[3125];
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
	fwd_len3125_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbInRe, gbInIm, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 125, maximum transforms: 1, Passes: 5
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(125)
fft_back_op_len3125( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[3125];
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
	back_len3125_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbInRe, gbInIm, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}
