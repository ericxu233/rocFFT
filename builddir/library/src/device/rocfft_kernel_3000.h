#pragma once
#include "rocfft_butterfly_template.h"

////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass0_len3000(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	(*R0) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 0 )*stride_in, load_cb_data, nullptr);

	(*R10) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 0 )*stride_in, load_cb_data, nullptr);

	(*R20) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 0 )*stride_in, load_cb_data, nullptr);

	(*R1) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 300 )*stride_in, load_cb_data, nullptr);

	(*R11) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 300 )*stride_in, load_cb_data, nullptr);

	(*R21) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 300 )*stride_in, load_cb_data, nullptr);

	(*R2) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 600 )*stride_in, load_cb_data, nullptr);

	(*R12) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 600 )*stride_in, load_cb_data, nullptr);

	(*R22) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 600 )*stride_in, load_cb_data, nullptr);

	(*R3) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 900 )*stride_in, load_cb_data, nullptr);

	(*R13) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 900 )*stride_in, load_cb_data, nullptr);

	(*R23) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 900 )*stride_in, load_cb_data, nullptr);

	(*R4) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 1200 )*stride_in, load_cb_data, nullptr);

	(*R14) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 1200 )*stride_in, load_cb_data, nullptr);

	(*R24) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 1200 )*stride_in, load_cb_data, nullptr);

	(*R5) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 1500 )*stride_in, load_cb_data, nullptr);

	(*R15) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 1500 )*stride_in, load_cb_data, nullptr);

	(*R25) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 1500 )*stride_in, load_cb_data, nullptr);

	(*R6) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 1800 )*stride_in, load_cb_data, nullptr);

	(*R16) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 1800 )*stride_in, load_cb_data, nullptr);

	(*R26) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 1800 )*stride_in, load_cb_data, nullptr);

	(*R7) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 2100 )*stride_in, load_cb_data, nullptr);

	(*R17) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 2100 )*stride_in, load_cb_data, nullptr);

	(*R27) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 2100 )*stride_in, load_cb_data, nullptr);

	(*R8) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 2400 )*stride_in, load_cb_data, nullptr);

	(*R18) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 2400 )*stride_in, load_cb_data, nullptr);

	(*R28) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 2400 )*stride_in, load_cb_data, nullptr);

	(*R9) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 2700 )*stride_in, load_cb_data, nullptr);

	(*R19) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 2700 )*stride_in, load_cb_data, nullptr);

	(*R29) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 2700 )*stride_in, load_cb_data, nullptr);

	}



	FwdRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	FwdRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);
	FwdRad10B1(R20, R21, R22, R23, R24, R25, R26, R27, R28, R29);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 4 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 5 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 6 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 7 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 8 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 9 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 1 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 2 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 3 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 4 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 5 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 6 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 7 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 8 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 9 ) ] = (*R29).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 300 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 300 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 300 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 600 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 600 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 600 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 900 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 900 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 900 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1200 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1200 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1200 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1500 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1500 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1500 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2100 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2100 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2100 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2400 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2400 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2400 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 4 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 5 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 6 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 7 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 8 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 9 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 1 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 2 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 3 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 4 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 5 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 6 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 7 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 8 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 9 ) ] = (*R29).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 300 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 300 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 300 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 600 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 600 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 600 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 900 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 900 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 900 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1200 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1200 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1200 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1500 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1500 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1500 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2100 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2100 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2100 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2400 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2400 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2400 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass0_len3000(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R20).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R20).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 300 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 300 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 300 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 300 )*stride_in];
	(*R21).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 300 )*stride_in];
	(*R21).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 300 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 600 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 600 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 600 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 600 )*stride_in];
	(*R22).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 600 )*stride_in];
	(*R22).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 600 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 900 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 900 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 900 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 900 )*stride_in];
	(*R23).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 900 )*stride_in];
	(*R23).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 900 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1200 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1200 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1200 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1200 )*stride_in];
	(*R24).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1200 )*stride_in];
	(*R24).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1200 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1500 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1500 )*stride_in];
	(*R15).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1500 )*stride_in];
	(*R15).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1500 )*stride_in];
	(*R25).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1500 )*stride_in];
	(*R25).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1500 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1800 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1800 )*stride_in];
	(*R16).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1800 )*stride_in];
	(*R16).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1800 )*stride_in];
	(*R26).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1800 )*stride_in];
	(*R26).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1800 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2100 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2100 )*stride_in];
	(*R17).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2100 )*stride_in];
	(*R17).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2100 )*stride_in];
	(*R27).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2100 )*stride_in];
	(*R27).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2100 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2400 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2400 )*stride_in];
	(*R18).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2400 )*stride_in];
	(*R18).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2400 )*stride_in];
	(*R28).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2400 )*stride_in];
	(*R28).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2400 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2700 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2700 )*stride_in];
	(*R19).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2700 )*stride_in];
	(*R19).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2700 )*stride_in];
	(*R29).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2700 )*stride_in];
	(*R29).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2700 )*stride_in];
	}



	FwdRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	FwdRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);
	FwdRad10B1(R20, R21, R22, R23, R24, R25, R26, R27, R28, R29);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 4 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 5 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 6 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 7 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 8 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 9 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 1 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 2 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 3 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 4 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 5 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 6 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 7 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 8 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 9 ) ] = (*R29).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 300 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 300 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 300 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 600 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 600 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 600 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 900 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 900 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 900 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1200 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1200 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1200 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1500 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1500 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1500 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2100 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2100 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2100 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2400 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2400 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2400 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 4 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 5 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 6 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 7 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 8 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 9 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 1 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 2 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 3 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 4 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 5 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 6 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 7 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 8 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 9 ) ] = (*R29).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 300 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 300 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 300 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 600 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 600 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 600 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 900 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 900 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 900 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1200 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1200 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1200 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1500 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1500 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1500 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2100 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2100 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2100 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2400 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2400 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2400 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass1_len3000(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
{




	{
		T W = twiddles[9 + 9*((3*me + 0)%10) + 0];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 1];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 2];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 3];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 4];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 5];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 6];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 7];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 8];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 0];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 1];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 2];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 3];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 4];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 5];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 6];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 7];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 8];
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
		T W = twiddles[9 + 9*((3*me + 2)%10) + 0];
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
		T W = twiddles[9 + 9*((3*me + 2)%10) + 1];
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
		T W = twiddles[9 + 9*((3*me + 2)%10) + 2];
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
		T W = twiddles[9 + 9*((3*me + 2)%10) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	{
		T W = twiddles[9 + 9*((3*me + 2)%10) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
	}

	{
		T W = twiddles[9 + 9*((3*me + 2)%10) + 5];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R26).x; ry = (*R26).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R26).x = TR;
		(*R26).y = TI;
	}

	{
		T W = twiddles[9 + 9*((3*me + 2)%10) + 6];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R27).x; ry = (*R27).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R27).x = TR;
		(*R27).y = TI;
	}

	{
		T W = twiddles[9 + 9*((3*me + 2)%10) + 7];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R28).x; ry = (*R28).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R28).x = TR;
		(*R28).y = TI;
	}

	{
		T W = twiddles[9 + 9*((3*me + 2)%10) + 8];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	FwdRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	FwdRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);
	FwdRad10B1(R20, R21, R22, R23, R24, R25, R26, R27, R28, R29);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 10 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 20 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 30 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 40 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 50 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 60 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 70 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 80 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 90 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 10 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 20 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 30 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 40 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 50 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 60 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 70 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 80 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 90 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 10 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 20 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 30 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 40 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 50 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 60 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 70 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 80 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 90 ) ] = (*R29).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 300 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 300 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 300 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 600 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 600 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 600 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 900 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 900 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 900 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1200 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1200 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1200 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1500 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1500 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1500 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2100 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2100 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2100 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2400 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2400 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2400 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 10 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 20 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 30 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 40 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 50 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 60 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 70 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 80 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 90 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 10 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 20 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 30 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 40 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 50 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 60 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 70 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 80 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 90 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 10 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 20 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 30 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 40 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 50 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 60 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 70 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 80 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 90 ) ] = (*R29).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 300 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 300 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 300 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 600 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 600 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 600 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 900 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 900 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 900 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1200 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1200 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1200 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1500 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1500 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1500 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2100 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2100 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2100 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2400 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2400 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2400 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass2_len3000(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
{




	{
		T W = twiddles[99 + 9*((3*me + 0)%100) + 0];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 1];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 2];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 3];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 4];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 5];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 6];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 7];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 8];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 0];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 1];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 2];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 3];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 4];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 5];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 6];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 7];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 8];
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
		T W = twiddles[99 + 9*((3*me + 2)%100) + 0];
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
		T W = twiddles[99 + 9*((3*me + 2)%100) + 1];
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
		T W = twiddles[99 + 9*((3*me + 2)%100) + 2];
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
		T W = twiddles[99 + 9*((3*me + 2)%100) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	{
		T W = twiddles[99 + 9*((3*me + 2)%100) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
	}

	{
		T W = twiddles[99 + 9*((3*me + 2)%100) + 5];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R26).x; ry = (*R26).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R26).x = TR;
		(*R26).y = TI;
	}

	{
		T W = twiddles[99 + 9*((3*me + 2)%100) + 6];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R27).x; ry = (*R27).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R27).x = TR;
		(*R27).y = TI;
	}

	{
		T W = twiddles[99 + 9*((3*me + 2)%100) + 7];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R28).x; ry = (*R28).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R28).x = TR;
		(*R28).y = TI;
	}

	{
		T W = twiddles[99 + 9*((3*me + 2)%100) + 8];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	FwdRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	FwdRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);
	FwdRad10B1(R20, R21, R22, R23, R24, R25, R26, R27, R28, R29);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 100 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 200 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 300 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 400 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 500 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 600 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 700 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 800 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 900 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 100 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 200 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 300 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 400 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 500 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 600 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 700 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 800 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 900 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 100 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 200 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 300 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 400 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 500 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 600 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 700 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 800 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 900 ) ] = (*R29).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*10 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*10 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*10 + 2 + 0 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*10 + 3 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*10 + 4 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*10 + 5 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*10 + 6 + 0 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*10 + 7 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*10 + 8 + 0 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*10 + 9 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*10 + 0 + 1000 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*10 + 1 + 1000 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*10 + 2 + 1000 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*10 + 3 + 1000 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*10 + 4 + 1000 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*10 + 5 + 1000 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*10 + 6 + 1000 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*10 + 7 + 1000 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*10 + 8 + 1000 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*10 + 9 + 1000 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*10 + 0 + 2000 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*10 + 1 + 2000 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*10 + 2 + 2000 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*10 + 3 + 2000 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*10 + 4 + 2000 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*10 + 5 + 2000 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*10 + 6 + 2000 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*10 + 7 + 2000 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*10 + 8 + 2000 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*10 + 9 + 2000 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 100 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 200 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 300 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 400 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 500 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 600 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 700 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 800 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 900 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 100 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 200 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 300 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 400 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 500 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 600 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 700 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 800 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 900 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 100 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 200 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 300 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 400 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 500 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 600 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 700 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 800 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 900 ) ] = (*R29).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*10 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*10 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*10 + 2 + 0 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*10 + 3 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*10 + 4 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*10 + 5 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*10 + 6 + 0 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*10 + 7 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*10 + 8 + 0 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*10 + 9 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*10 + 0 + 1000 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*10 + 1 + 1000 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*10 + 2 + 1000 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*10 + 3 + 1000 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*10 + 4 + 1000 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*10 + 5 + 1000 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*10 + 6 + 1000 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*10 + 7 + 1000 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*10 + 8 + 1000 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*10 + 9 + 1000 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*10 + 0 + 2000 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*10 + 1 + 2000 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*10 + 2 + 2000 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*10 + 3 + 2000 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*10 + 4 + 2000 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*10 + 5 + 2000 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*10 + 6 + 2000 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*10 + 7 + 2000 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*10 + 8 + 2000 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*10 + 9 + 2000 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass3_len3000(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[999 + 2*((10*me + 0)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 0)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 1)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 1)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 2)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 2)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 3)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 3)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 4)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 4)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 5)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 5)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 6)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 6)%1000) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R20).x; ry = (*R20).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R20).x = TR;
		(*R20).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 7)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 7)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 8)%1000) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 8)%1000) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R26).x; ry = (*R26).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R26).x = TR;
		(*R26).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 9)%1000) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R28).x; ry = (*R28).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R28).x = TR;
		(*R28).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 9)%1000) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);
	FwdRad3B1(R9, R10, R11);
	FwdRad3B1(R12, R13, R14);
	FwdRad3B1(R15, R16, R17);
	FwdRad3B1(R18, R19, R20);
	FwdRad3B1(R21, R22, R23);
	FwdRad3B1(R24, R25, R26);
	FwdRad3B1(R27, R28, R29);


	if(rw)
	{
	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);

	 //Optimization: coalescing into float4/double4 write
	if(sb == SB_UNIT && !store_cb_fn) {
	vector4_type_t<T> *buff4g = (vector4_type_t<T>*)(bufOut+outOffset);
	
	buff4g[ 5*me + 0 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R0).x, (*R0).y, (*R3).x, (*R3).y) ;
	buff4g[ 5*me + 1 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R6).x, (*R6).y, (*R9).x, (*R9).y) ;
	buff4g[ 5*me + 2 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R12).x, (*R12).y, (*R15).x, (*R15).y) ;
	buff4g[ 5*me + 3 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R18).x, (*R18).y, (*R21).x, (*R21).y) ;
	buff4g[ 5*me + 4 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R24).x, (*R24).y, (*R27).x, (*R27).y) ;
	buff4g[ 5*me + 0 + 500 ] = lib_make_vector4< vector4_type_t<T> >((*R1).x, (*R1).y, (*R4).x, (*R4).y) ;
	buff4g[ 5*me + 1 + 500 ] = lib_make_vector4< vector4_type_t<T> >((*R7).x, (*R7).y, (*R10).x, (*R10).y) ;
	buff4g[ 5*me + 2 + 500 ] = lib_make_vector4< vector4_type_t<T> >((*R13).x, (*R13).y, (*R16).x, (*R16).y) ;
	buff4g[ 5*me + 3 + 500 ] = lib_make_vector4< vector4_type_t<T> >((*R19).x, (*R19).y, (*R22).x, (*R22).y) ;
	buff4g[ 5*me + 4 + 500 ] = lib_make_vector4< vector4_type_t<T> >((*R25).x, (*R25).y, (*R28).x, (*R28).y) ;
	buff4g[ 5*me + 0 + 1000 ] = lib_make_vector4< vector4_type_t<T> >((*R2).x, (*R2).y, (*R5).x, (*R5).y) ;
	buff4g[ 5*me + 1 + 1000 ] = lib_make_vector4< vector4_type_t<T> >((*R8).x, (*R8).y, (*R11).x, (*R11).y) ;
	buff4g[ 5*me + 2 + 1000 ] = lib_make_vector4< vector4_type_t<T> >((*R14).x, (*R14).y, (*R17).x, (*R17).y) ;
	buff4g[ 5*me + 3 + 1000 ] = lib_make_vector4< vector4_type_t<T> >((*R20).x, (*R20).y, (*R23).x, (*R23).y) ;
	buff4g[ 5*me + 4 + 1000 ] = lib_make_vector4< vector4_type_t<T> >((*R26).x, (*R26).y, (*R29).x, (*R29).y) ;
	}
	else{ // such optimization is not possible 
	store_cb(bufOut, outOffset + ( 10*me + 0 + 0 )*stride_out, (*R0), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 1 + 0 )*stride_out, (*R3), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 2 + 0 )*stride_out, (*R6), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 3 + 0 )*stride_out, (*R9), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 4 + 0 )*stride_out, (*R12), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 5 + 0 )*stride_out, (*R15), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 6 + 0 )*stride_out, (*R18), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 7 + 0 )*stride_out, (*R21), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 8 + 0 )*stride_out, (*R24), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 9 + 0 )*stride_out, (*R27), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 0 + 1000 )*stride_out, (*R1), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 1 + 1000 )*stride_out, (*R4), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 2 + 1000 )*stride_out, (*R7), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 3 + 1000 )*stride_out, (*R10), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 4 + 1000 )*stride_out, (*R13), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 5 + 1000 )*stride_out, (*R16), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 6 + 1000 )*stride_out, (*R19), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 7 + 1000 )*stride_out, (*R22), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 8 + 1000 )*stride_out, (*R25), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 9 + 1000 )*stride_out, (*R28), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 0 + 2000 )*stride_out, (*R2), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 1 + 2000 )*stride_out, (*R5), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 2 + 2000 )*stride_out, (*R8), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 3 + 2000 )*stride_out, (*R11), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 4 + 2000 )*stride_out, (*R14), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 5 + 2000 )*stride_out, (*R17), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 6 + 2000 )*stride_out, (*R20), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 7 + 2000 )*stride_out, (*R23), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 8 + 2000 )*stride_out, (*R26), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 9 + 2000 )*stride_out, (*R29), store_cb_data, nullptr);

	}
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass3_len3000(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[999 + 2*((10*me + 0)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 0)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 1)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 1)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 2)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 2)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 3)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 3)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 4)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 4)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 5)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 5)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 6)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 6)%1000) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R20).x; ry = (*R20).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R20).x = TR;
		(*R20).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 7)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 7)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 8)%1000) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 8)%1000) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R26).x; ry = (*R26).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R26).x = TR;
		(*R26).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 9)%1000) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R28).x; ry = (*R28).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R28).x = TR;
		(*R28).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 9)%1000) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	FwdRad3B1(R0, R1, R2);
	FwdRad3B1(R3, R4, R5);
	FwdRad3B1(R6, R7, R8);
	FwdRad3B1(R9, R10, R11);
	FwdRad3B1(R12, R13, R14);
	FwdRad3B1(R15, R16, R17);
	FwdRad3B1(R18, R19, R20);
	FwdRad3B1(R21, R22, R23);
	FwdRad3B1(R24, R25, R26);
	FwdRad3B1(R27, R28, R29);


	if(rw)
	{
	bufOutRe[outOffset + ( 10*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 10*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 10*me + 1 + 0 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 10*me + 1 + 0 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 10*me + 2 + 0 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 10*me + 2 + 0 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 10*me + 3 + 0 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 10*me + 3 + 0 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 10*me + 4 + 0 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 10*me + 4 + 0 )*stride_out] = (*R12).y;
	bufOutRe[outOffset + ( 10*me + 5 + 0 )*stride_out] = (*R15).x;
	bufOutIm[outOffset + ( 10*me + 5 + 0 )*stride_out] = (*R15).y;
	bufOutRe[outOffset + ( 10*me + 6 + 0 )*stride_out] = (*R18).x;
	bufOutIm[outOffset + ( 10*me + 6 + 0 )*stride_out] = (*R18).y;
	bufOutRe[outOffset + ( 10*me + 7 + 0 )*stride_out] = (*R21).x;
	bufOutIm[outOffset + ( 10*me + 7 + 0 )*stride_out] = (*R21).y;
	bufOutRe[outOffset + ( 10*me + 8 + 0 )*stride_out] = (*R24).x;
	bufOutIm[outOffset + ( 10*me + 8 + 0 )*stride_out] = (*R24).y;
	bufOutRe[outOffset + ( 10*me + 9 + 0 )*stride_out] = (*R27).x;
	bufOutIm[outOffset + ( 10*me + 9 + 0 )*stride_out] = (*R27).y;
	bufOutRe[outOffset + ( 10*me + 0 + 1000 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 10*me + 0 + 1000 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 10*me + 1 + 1000 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 10*me + 1 + 1000 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 10*me + 2 + 1000 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 10*me + 2 + 1000 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 10*me + 3 + 1000 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 10*me + 3 + 1000 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 10*me + 4 + 1000 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 10*me + 4 + 1000 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 10*me + 5 + 1000 )*stride_out] = (*R16).x;
	bufOutIm[outOffset + ( 10*me + 5 + 1000 )*stride_out] = (*R16).y;
	bufOutRe[outOffset + ( 10*me + 6 + 1000 )*stride_out] = (*R19).x;
	bufOutIm[outOffset + ( 10*me + 6 + 1000 )*stride_out] = (*R19).y;
	bufOutRe[outOffset + ( 10*me + 7 + 1000 )*stride_out] = (*R22).x;
	bufOutIm[outOffset + ( 10*me + 7 + 1000 )*stride_out] = (*R22).y;
	bufOutRe[outOffset + ( 10*me + 8 + 1000 )*stride_out] = (*R25).x;
	bufOutIm[outOffset + ( 10*me + 8 + 1000 )*stride_out] = (*R25).y;
	bufOutRe[outOffset + ( 10*me + 9 + 1000 )*stride_out] = (*R28).x;
	bufOutIm[outOffset + ( 10*me + 9 + 1000 )*stride_out] = (*R28).y;
	bufOutRe[outOffset + ( 10*me + 0 + 2000 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 10*me + 0 + 2000 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 10*me + 1 + 2000 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 10*me + 1 + 2000 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 10*me + 2 + 2000 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 10*me + 2 + 2000 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 10*me + 3 + 2000 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 10*me + 3 + 2000 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 10*me + 4 + 2000 )*stride_out] = (*R14).x;
	bufOutIm[outOffset + ( 10*me + 4 + 2000 )*stride_out] = (*R14).y;
	bufOutRe[outOffset + ( 10*me + 5 + 2000 )*stride_out] = (*R17).x;
	bufOutIm[outOffset + ( 10*me + 5 + 2000 )*stride_out] = (*R17).y;
	bufOutRe[outOffset + ( 10*me + 6 + 2000 )*stride_out] = (*R20).x;
	bufOutIm[outOffset + ( 10*me + 6 + 2000 )*stride_out] = (*R20).y;
	bufOutRe[outOffset + ( 10*me + 7 + 2000 )*stride_out] = (*R23).x;
	bufOutIm[outOffset + ( 10*me + 7 + 2000 )*stride_out] = (*R23).y;
	bufOutRe[outOffset + ( 10*me + 8 + 2000 )*stride_out] = (*R26).x;
	bufOutIm[outOffset + ( 10*me + 8 + 2000 )*stride_out] = (*R26).y;
	bufOutRe[outOffset + ( 10*me + 9 + 2000 )*stride_out] = (*R29).x;
	bufOutIm[outOffset + ( 10*me + 9 + 2000 )*stride_out] = (*R29).y;
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass0_len3000(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	(*R0) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 0 )*stride_in, load_cb_data, nullptr);

	(*R10) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 0 )*stride_in, load_cb_data, nullptr);

	(*R20) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 0 )*stride_in, load_cb_data, nullptr);

	(*R1) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 300 )*stride_in, load_cb_data, nullptr);

	(*R11) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 300 )*stride_in, load_cb_data, nullptr);

	(*R21) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 300 )*stride_in, load_cb_data, nullptr);

	(*R2) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 600 )*stride_in, load_cb_data, nullptr);

	(*R12) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 600 )*stride_in, load_cb_data, nullptr);

	(*R22) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 600 )*stride_in, load_cb_data, nullptr);

	(*R3) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 900 )*stride_in, load_cb_data, nullptr);

	(*R13) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 900 )*stride_in, load_cb_data, nullptr);

	(*R23) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 900 )*stride_in, load_cb_data, nullptr);

	(*R4) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 1200 )*stride_in, load_cb_data, nullptr);

	(*R14) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 1200 )*stride_in, load_cb_data, nullptr);

	(*R24) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 1200 )*stride_in, load_cb_data, nullptr);

	(*R5) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 1500 )*stride_in, load_cb_data, nullptr);

	(*R15) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 1500 )*stride_in, load_cb_data, nullptr);

	(*R25) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 1500 )*stride_in, load_cb_data, nullptr);

	(*R6) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 1800 )*stride_in, load_cb_data, nullptr);

	(*R16) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 1800 )*stride_in, load_cb_data, nullptr);

	(*R26) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 1800 )*stride_in, load_cb_data, nullptr);

	(*R7) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 2100 )*stride_in, load_cb_data, nullptr);

	(*R17) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 2100 )*stride_in, load_cb_data, nullptr);

	(*R27) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 2100 )*stride_in, load_cb_data, nullptr);

	(*R8) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 2400 )*stride_in, load_cb_data, nullptr);

	(*R18) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 2400 )*stride_in, load_cb_data, nullptr);

	(*R28) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 2400 )*stride_in, load_cb_data, nullptr);

	(*R9) = load_cb(bufIn, inOffset + ( 0 + me*3 + 0 + 2700 )*stride_in, load_cb_data, nullptr);

	(*R19) = load_cb(bufIn, inOffset + ( 0 + me*3 + 1 + 2700 )*stride_in, load_cb_data, nullptr);

	(*R29) = load_cb(bufIn, inOffset + ( 0 + me*3 + 2 + 2700 )*stride_in, load_cb_data, nullptr);

	}



	InvRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	InvRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);
	InvRad10B1(R20, R21, R22, R23, R24, R25, R26, R27, R28, R29);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 4 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 5 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 6 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 7 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 8 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 9 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 1 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 2 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 3 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 4 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 5 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 6 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 7 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 8 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 9 ) ] = (*R29).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 300 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 300 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 300 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 600 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 600 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 600 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 900 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 900 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 900 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1200 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1200 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1200 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1500 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1500 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1500 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2100 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2100 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2100 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2400 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2400 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2400 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 4 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 5 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 6 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 7 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 8 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 9 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 1 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 2 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 3 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 4 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 5 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 6 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 7 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 8 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 9 ) ] = (*R29).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 300 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 300 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 300 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 600 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 600 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 600 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 900 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 900 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 900 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1200 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1200 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1200 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1500 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1500 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1500 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2100 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2100 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2100 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2400 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2400 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2400 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass0_len3000(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R20).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R20).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 300 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 300 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 300 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 300 )*stride_in];
	(*R21).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 300 )*stride_in];
	(*R21).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 300 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 600 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 600 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 600 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 600 )*stride_in];
	(*R22).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 600 )*stride_in];
	(*R22).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 600 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 900 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 900 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 900 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 900 )*stride_in];
	(*R23).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 900 )*stride_in];
	(*R23).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 900 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1200 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1200 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1200 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1200 )*stride_in];
	(*R24).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1200 )*stride_in];
	(*R24).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1200 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1500 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1500 )*stride_in];
	(*R15).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1500 )*stride_in];
	(*R15).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1500 )*stride_in];
	(*R25).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1500 )*stride_in];
	(*R25).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1500 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1800 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1800 )*stride_in];
	(*R16).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1800 )*stride_in];
	(*R16).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1800 )*stride_in];
	(*R26).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1800 )*stride_in];
	(*R26).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1800 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2100 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2100 )*stride_in];
	(*R17).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2100 )*stride_in];
	(*R17).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2100 )*stride_in];
	(*R27).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2100 )*stride_in];
	(*R27).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2100 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2400 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2400 )*stride_in];
	(*R18).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2400 )*stride_in];
	(*R18).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2400 )*stride_in];
	(*R28).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2400 )*stride_in];
	(*R28).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2400 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2700 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2700 )*stride_in];
	(*R19).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2700 )*stride_in];
	(*R19).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2700 )*stride_in];
	(*R29).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2700 )*stride_in];
	(*R29).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2700 )*stride_in];
	}



	InvRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	InvRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);
	InvRad10B1(R20, R21, R22, R23, R24, R25, R26, R27, R28, R29);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 1 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 2 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 3 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 4 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 5 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 6 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 7 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 8 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 9 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 1 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 2 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 3 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 4 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 5 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 6 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 7 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 8 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 9 ) ] = (*R29).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 300 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 300 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 300 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 600 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 600 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 600 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 900 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 900 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 900 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1200 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1200 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1200 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1500 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1500 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1500 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2100 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2100 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2100 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2400 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2400 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2400 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 0)/1)*10 + (3*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 1 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 2 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 3 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 4 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 5 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 6 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 7 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 8 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((3*me + 1)/1)*10 + (3*me + 1)%1 + 9 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 1 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 2 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 3 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 4 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 5 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 6 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 7 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 8 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((3*me + 2)/1)*10 + (3*me + 2)%1 + 9 ) ] = (*R29).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 300 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 300 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 300 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 600 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 600 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 600 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 900 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 900 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 900 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1200 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1200 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1200 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1500 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1500 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1500 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2100 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2100 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2100 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2400 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2400 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2400 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass1_len3000(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
{




	{
		T W = twiddles[9 + 9*((3*me + 0)%10) + 0];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 1];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 2];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 3];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 4];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 5];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 6];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 7];
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
		T W = twiddles[9 + 9*((3*me + 0)%10) + 8];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 0];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 1];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 2];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 3];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 4];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 5];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 6];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 7];
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
		T W = twiddles[9 + 9*((3*me + 1)%10) + 8];
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
		T W = twiddles[9 + 9*((3*me + 2)%10) + 0];
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
		T W = twiddles[9 + 9*((3*me + 2)%10) + 1];
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
		T W = twiddles[9 + 9*((3*me + 2)%10) + 2];
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
		T W = twiddles[9 + 9*((3*me + 2)%10) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	{
		T W = twiddles[9 + 9*((3*me + 2)%10) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
	}

	{
		T W = twiddles[9 + 9*((3*me + 2)%10) + 5];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R26).x; ry = (*R26).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R26).x = TR;
		(*R26).y = TI;
	}

	{
		T W = twiddles[9 + 9*((3*me + 2)%10) + 6];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R27).x; ry = (*R27).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R27).x = TR;
		(*R27).y = TI;
	}

	{
		T W = twiddles[9 + 9*((3*me + 2)%10) + 7];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R28).x; ry = (*R28).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R28).x = TR;
		(*R28).y = TI;
	}

	{
		T W = twiddles[9 + 9*((3*me + 2)%10) + 8];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	InvRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	InvRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);
	InvRad10B1(R20, R21, R22, R23, R24, R25, R26, R27, R28, R29);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 10 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 20 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 30 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 40 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 50 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 60 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 70 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 80 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 90 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 10 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 20 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 30 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 40 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 50 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 60 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 70 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 80 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 90 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 10 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 20 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 30 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 40 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 50 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 60 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 70 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 80 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 90 ) ] = (*R29).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 300 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 300 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 300 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 600 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 600 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 600 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 900 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 900 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 900 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1200 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1200 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1200 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1500 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1500 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1500 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2100 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2100 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2100 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2400 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2400 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2400 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 10 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 20 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 30 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 40 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 50 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 60 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 70 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 80 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 0)/10)*100 + (3*me + 0)%10 + 90 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 10 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 20 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 30 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 40 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 50 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 60 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 70 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 80 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((3*me + 1)/10)*100 + (3*me + 1)%10 + 90 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 10 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 20 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 30 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 40 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 50 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 60 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 70 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 80 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((3*me + 2)/10)*100 + (3*me + 2)%10 + 90 ) ] = (*R29).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 300 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 300 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 300 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 600 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 600 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 600 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 900 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 900 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 900 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1200 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1200 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1200 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1500 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1500 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1500 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2100 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2100 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2100 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2400 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2400 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2400 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2700 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2700 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2700 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass2_len3000(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
{




	{
		T W = twiddles[99 + 9*((3*me + 0)%100) + 0];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 1];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 2];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 3];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 4];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 5];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 6];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 7];
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
		T W = twiddles[99 + 9*((3*me + 0)%100) + 8];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 0];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 1];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 2];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 3];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 4];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 5];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 6];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 7];
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
		T W = twiddles[99 + 9*((3*me + 1)%100) + 8];
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
		T W = twiddles[99 + 9*((3*me + 2)%100) + 0];
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
		T W = twiddles[99 + 9*((3*me + 2)%100) + 1];
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
		T W = twiddles[99 + 9*((3*me + 2)%100) + 2];
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
		T W = twiddles[99 + 9*((3*me + 2)%100) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R24).x; ry = (*R24).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R24).x = TR;
		(*R24).y = TI;
	}

	{
		T W = twiddles[99 + 9*((3*me + 2)%100) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
	}

	{
		T W = twiddles[99 + 9*((3*me + 2)%100) + 5];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R26).x; ry = (*R26).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R26).x = TR;
		(*R26).y = TI;
	}

	{
		T W = twiddles[99 + 9*((3*me + 2)%100) + 6];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R27).x; ry = (*R27).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R27).x = TR;
		(*R27).y = TI;
	}

	{
		T W = twiddles[99 + 9*((3*me + 2)%100) + 7];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R28).x; ry = (*R28).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R28).x = TR;
		(*R28).y = TI;
	}

	{
		T W = twiddles[99 + 9*((3*me + 2)%100) + 8];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	InvRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	InvRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);
	InvRad10B1(R20, R21, R22, R23, R24, R25, R26, R27, R28, R29);


	if(rw)
	{
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 100 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 200 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 300 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 400 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 500 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 600 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 700 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 800 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 900 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 100 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 200 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 300 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 400 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 500 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 600 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 700 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 800 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 900 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 100 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 200 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 300 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 400 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 500 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 600 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 700 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 800 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 900 ) ] = (*R29).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*10 + 0 + 0 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*10 + 1 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*10 + 2 + 0 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*10 + 3 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*10 + 4 + 0 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*10 + 5 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*10 + 6 + 0 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*10 + 7 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*10 + 8 + 0 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*10 + 9 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*10 + 0 + 1000 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*10 + 1 + 1000 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*10 + 2 + 1000 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*10 + 3 + 1000 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*10 + 4 + 1000 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*10 + 5 + 1000 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*10 + 6 + 1000 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*10 + 7 + 1000 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*10 + 8 + 1000 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*10 + 9 + 1000 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*10 + 0 + 2000 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*10 + 1 + 2000 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*10 + 2 + 2000 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*10 + 3 + 2000 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*10 + 4 + 2000 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*10 + 5 + 2000 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*10 + 6 + 2000 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*10 + 7 + 2000 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*10 + 8 + 2000 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*10 + 9 + 2000 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 100 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 200 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 300 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 400 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 500 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 600 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 700 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 800 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((3*me + 0)/100)*1000 + (3*me + 0)%100 + 900 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 100 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 200 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 300 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 400 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 500 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 600 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 700 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 800 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((3*me + 1)/100)*1000 + (3*me + 1)%100 + 900 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 100 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 200 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 300 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 400 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 500 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 600 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 700 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 800 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((3*me + 2)/100)*1000 + (3*me + 2)%100 + 900 ) ] = (*R29).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*10 + 0 + 0 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*10 + 1 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*10 + 2 + 0 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*10 + 3 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*10 + 4 + 0 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*10 + 5 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*10 + 6 + 0 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*10 + 7 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*10 + 8 + 0 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*10 + 9 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*10 + 0 + 1000 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*10 + 1 + 1000 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*10 + 2 + 1000 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*10 + 3 + 1000 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*10 + 4 + 1000 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*10 + 5 + 1000 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*10 + 6 + 1000 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*10 + 7 + 1000 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*10 + 8 + 1000 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*10 + 9 + 1000 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*10 + 0 + 2000 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*10 + 1 + 2000 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*10 + 2 + 2000 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*10 + 3 + 2000 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*10 + 4 + 2000 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*10 + 5 + 2000 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*10 + 6 + 2000 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*10 + 7 + 2000 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*10 + 8 + 2000 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*10 + 9 + 2000 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass3_len3000(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[999 + 2*((10*me + 0)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 0)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 1)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 1)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 2)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 2)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 3)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 3)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 4)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 4)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 5)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 5)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 6)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 6)%1000) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R20).x; ry = (*R20).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R20).x = TR;
		(*R20).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 7)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 7)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 8)%1000) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 8)%1000) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R26).x; ry = (*R26).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R26).x = TR;
		(*R26).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 9)%1000) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R28).x; ry = (*R28).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R28).x = TR;
		(*R28).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 9)%1000) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);
	InvRad3B1(R9, R10, R11);
	InvRad3B1(R12, R13, R14);
	InvRad3B1(R15, R16, R17);
	InvRad3B1(R18, R19, R20);
	InvRad3B1(R21, R22, R23);
	InvRad3B1(R24, R25, R26);
	InvRad3B1(R27, R28, R29);


	if(rw)
	{
	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);

	 //Optimization: coalescing into float4/double4 write
	if(sb == SB_UNIT && !store_cb_fn) {
	vector4_type_t<T> *buff4g = (vector4_type_t<T>*)(bufOut+outOffset);
	
	buff4g[ 5*me + 0 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R0).x, (*R0).y, (*R3).x, (*R3).y) ;
	buff4g[ 5*me + 1 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R6).x, (*R6).y, (*R9).x, (*R9).y) ;
	buff4g[ 5*me + 2 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R12).x, (*R12).y, (*R15).x, (*R15).y) ;
	buff4g[ 5*me + 3 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R18).x, (*R18).y, (*R21).x, (*R21).y) ;
	buff4g[ 5*me + 4 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R24).x, (*R24).y, (*R27).x, (*R27).y) ;
	buff4g[ 5*me + 0 + 500 ] = lib_make_vector4< vector4_type_t<T> >((*R1).x, (*R1).y, (*R4).x, (*R4).y) ;
	buff4g[ 5*me + 1 + 500 ] = lib_make_vector4< vector4_type_t<T> >((*R7).x, (*R7).y, (*R10).x, (*R10).y) ;
	buff4g[ 5*me + 2 + 500 ] = lib_make_vector4< vector4_type_t<T> >((*R13).x, (*R13).y, (*R16).x, (*R16).y) ;
	buff4g[ 5*me + 3 + 500 ] = lib_make_vector4< vector4_type_t<T> >((*R19).x, (*R19).y, (*R22).x, (*R22).y) ;
	buff4g[ 5*me + 4 + 500 ] = lib_make_vector4< vector4_type_t<T> >((*R25).x, (*R25).y, (*R28).x, (*R28).y) ;
	buff4g[ 5*me + 0 + 1000 ] = lib_make_vector4< vector4_type_t<T> >((*R2).x, (*R2).y, (*R5).x, (*R5).y) ;
	buff4g[ 5*me + 1 + 1000 ] = lib_make_vector4< vector4_type_t<T> >((*R8).x, (*R8).y, (*R11).x, (*R11).y) ;
	buff4g[ 5*me + 2 + 1000 ] = lib_make_vector4< vector4_type_t<T> >((*R14).x, (*R14).y, (*R17).x, (*R17).y) ;
	buff4g[ 5*me + 3 + 1000 ] = lib_make_vector4< vector4_type_t<T> >((*R20).x, (*R20).y, (*R23).x, (*R23).y) ;
	buff4g[ 5*me + 4 + 1000 ] = lib_make_vector4< vector4_type_t<T> >((*R26).x, (*R26).y, (*R29).x, (*R29).y) ;
	}
	else{ // such optimization is not possible 
	store_cb(bufOut, outOffset + ( 10*me + 0 + 0 )*stride_out, (*R0), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 1 + 0 )*stride_out, (*R3), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 2 + 0 )*stride_out, (*R6), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 3 + 0 )*stride_out, (*R9), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 4 + 0 )*stride_out, (*R12), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 5 + 0 )*stride_out, (*R15), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 6 + 0 )*stride_out, (*R18), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 7 + 0 )*stride_out, (*R21), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 8 + 0 )*stride_out, (*R24), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 9 + 0 )*stride_out, (*R27), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 0 + 1000 )*stride_out, (*R1), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 1 + 1000 )*stride_out, (*R4), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 2 + 1000 )*stride_out, (*R7), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 3 + 1000 )*stride_out, (*R10), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 4 + 1000 )*stride_out, (*R13), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 5 + 1000 )*stride_out, (*R16), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 6 + 1000 )*stride_out, (*R19), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 7 + 1000 )*stride_out, (*R22), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 8 + 1000 )*stride_out, (*R25), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 9 + 1000 )*stride_out, (*R28), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 0 + 2000 )*stride_out, (*R2), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 1 + 2000 )*stride_out, (*R5), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 2 + 2000 )*stride_out, (*R8), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 3 + 2000 )*stride_out, (*R11), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 4 + 2000 )*stride_out, (*R14), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 5 + 2000 )*stride_out, (*R17), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 6 + 2000 )*stride_out, (*R20), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 7 + 2000 )*stride_out, (*R23), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 8 + 2000 )*stride_out, (*R26), store_cb_data, nullptr);

	store_cb(bufOut, outOffset + ( 10*me + 9 + 2000 )*stride_out, (*R29), store_cb_data, nullptr);

	}
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass3_len3000(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{




	{
		T W = twiddles[999 + 2*((10*me + 0)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 0)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 1)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 1)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 2)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 2)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 3)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 3)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 4)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 4)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 5)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 5)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 6)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 6)%1000) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R20).x; ry = (*R20).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R20).x = TR;
		(*R20).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 7)%1000) + 0];
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
		T W = twiddles[999 + 2*((10*me + 7)%1000) + 1];
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
		T W = twiddles[999 + 2*((10*me + 8)%1000) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 8)%1000) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R26).x; ry = (*R26).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R26).x = TR;
		(*R26).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 9)%1000) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R28).x; ry = (*R28).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R28).x = TR;
		(*R28).y = TI;
	}

	{
		T W = twiddles[999 + 2*((10*me + 9)%1000) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	InvRad3B1(R0, R1, R2);
	InvRad3B1(R3, R4, R5);
	InvRad3B1(R6, R7, R8);
	InvRad3B1(R9, R10, R11);
	InvRad3B1(R12, R13, R14);
	InvRad3B1(R15, R16, R17);
	InvRad3B1(R18, R19, R20);
	InvRad3B1(R21, R22, R23);
	InvRad3B1(R24, R25, R26);
	InvRad3B1(R27, R28, R29);


	if(rw)
	{
	bufOutRe[outOffset + ( 10*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 10*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 10*me + 1 + 0 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 10*me + 1 + 0 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 10*me + 2 + 0 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 10*me + 2 + 0 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 10*me + 3 + 0 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 10*me + 3 + 0 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 10*me + 4 + 0 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 10*me + 4 + 0 )*stride_out] = (*R12).y;
	bufOutRe[outOffset + ( 10*me + 5 + 0 )*stride_out] = (*R15).x;
	bufOutIm[outOffset + ( 10*me + 5 + 0 )*stride_out] = (*R15).y;
	bufOutRe[outOffset + ( 10*me + 6 + 0 )*stride_out] = (*R18).x;
	bufOutIm[outOffset + ( 10*me + 6 + 0 )*stride_out] = (*R18).y;
	bufOutRe[outOffset + ( 10*me + 7 + 0 )*stride_out] = (*R21).x;
	bufOutIm[outOffset + ( 10*me + 7 + 0 )*stride_out] = (*R21).y;
	bufOutRe[outOffset + ( 10*me + 8 + 0 )*stride_out] = (*R24).x;
	bufOutIm[outOffset + ( 10*me + 8 + 0 )*stride_out] = (*R24).y;
	bufOutRe[outOffset + ( 10*me + 9 + 0 )*stride_out] = (*R27).x;
	bufOutIm[outOffset + ( 10*me + 9 + 0 )*stride_out] = (*R27).y;
	bufOutRe[outOffset + ( 10*me + 0 + 1000 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 10*me + 0 + 1000 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 10*me + 1 + 1000 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 10*me + 1 + 1000 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 10*me + 2 + 1000 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 10*me + 2 + 1000 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 10*me + 3 + 1000 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 10*me + 3 + 1000 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 10*me + 4 + 1000 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 10*me + 4 + 1000 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 10*me + 5 + 1000 )*stride_out] = (*R16).x;
	bufOutIm[outOffset + ( 10*me + 5 + 1000 )*stride_out] = (*R16).y;
	bufOutRe[outOffset + ( 10*me + 6 + 1000 )*stride_out] = (*R19).x;
	bufOutIm[outOffset + ( 10*me + 6 + 1000 )*stride_out] = (*R19).y;
	bufOutRe[outOffset + ( 10*me + 7 + 1000 )*stride_out] = (*R22).x;
	bufOutIm[outOffset + ( 10*me + 7 + 1000 )*stride_out] = (*R22).y;
	bufOutRe[outOffset + ( 10*me + 8 + 1000 )*stride_out] = (*R25).x;
	bufOutIm[outOffset + ( 10*me + 8 + 1000 )*stride_out] = (*R25).y;
	bufOutRe[outOffset + ( 10*me + 9 + 1000 )*stride_out] = (*R28).x;
	bufOutIm[outOffset + ( 10*me + 9 + 1000 )*stride_out] = (*R28).y;
	bufOutRe[outOffset + ( 10*me + 0 + 2000 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 10*me + 0 + 2000 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 10*me + 1 + 2000 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 10*me + 1 + 2000 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 10*me + 2 + 2000 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 10*me + 2 + 2000 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 10*me + 3 + 2000 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 10*me + 3 + 2000 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 10*me + 4 + 2000 )*stride_out] = (*R14).x;
	bufOutIm[outOffset + ( 10*me + 4 + 2000 )*stride_out] = (*R14).y;
	bufOutRe[outOffset + ( 10*me + 5 + 2000 )*stride_out] = (*R17).x;
	bufOutIm[outOffset + ( 10*me + 5 + 2000 )*stride_out] = (*R17).y;
	bufOutRe[outOffset + ( 10*me + 6 + 2000 )*stride_out] = (*R20).x;
	bufOutIm[outOffset + ( 10*me + 6 + 2000 )*stride_out] = (*R20).y;
	bufOutRe[outOffset + ( 10*me + 7 + 2000 )*stride_out] = (*R23).x;
	bufOutIm[outOffset + ( 10*me + 7 + 2000 )*stride_out] = (*R23).y;
	bufOutRe[outOffset + ( 10*me + 8 + 2000 )*stride_out] = (*R26).x;
	bufOutIm[outOffset + ( 10*me + 8 + 2000 )*stride_out] = (*R26).y;
	bufOutRe[outOffset + ( 10*me + 9 + 2000 )*stride_out] = (*R29).x;
	bufOutIm[outOffset + ( 10*me + 9 + 2000 )*stride_out] = (*R29).y;
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len3000_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	FwdPass0_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_fn, load_cb_data);
	FwdPass1_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass2_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass3_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len3000_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	InvPass0_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_fn, load_cb_data);
	InvPass1_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass2_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass3_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len3000_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	FwdPass0_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_fn, load_cb_data);
	FwdPass1_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass2_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass3_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len3000_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	InvPass0_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  gbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_fn, load_cb_data);
	InvPass1_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass2_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass3_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len3000_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	FwdPass0_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_fn, load_cb_data);
	FwdPass1_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass2_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass3_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len3000_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, T *gbOut, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	InvPass0_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_fn, load_cb_data);
	InvPass1_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass2_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass3_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len3000_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	FwdPass0_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_fn, load_cb_data);
	FwdPass1_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass2_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass3_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len3000_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, unsigned int iOffset, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, unsigned int oOffset, real_type_t<T> *lds, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	InvPass0_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, iOffset, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_fn, load_cb_data);
	InvPass1_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass2_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass3_len3000<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, oOffset, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 100, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_fwd_ip_len3000( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gb)
{

	__shared__ real_type_t<T> lds[3000];
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
	fwd_len3000_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, gb, ioOffset, gb, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 100, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_back_ip_len3000( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gb)
{

	__shared__ real_type_t<T> lds[3000];
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
	back_len3000_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, gb, ioOffset, gb, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 100, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_fwd_ip_len3000( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbRe, real_type_t<T> * __restrict__ gbIm)
{

	__shared__ real_type_t<T> lds[3000];
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
	fwd_len3000_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, gbRe, gbIm, ioOffset, gbRe, gbIm, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 100, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_back_ip_len3000( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbRe, real_type_t<T> * __restrict__ gbIm)
{

	__shared__ real_type_t<T> lds[3000];
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
	back_len3000_device<T, sb, cbtype>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, gbRe, gbIm, ioOffset, gbRe, gbIm, ioOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 100, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_fwd_op_len3000( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[3000];
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
	fwd_len3000_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 100, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_back_op_len3000( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[3000];
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
	back_len3000_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 100, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_fwd_op_len3000( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[3000];
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
	fwd_len3000_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 100, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_back_op_len3000( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[3000];
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
	back_len3000_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 100, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_fwd_op_len3000( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[3000];
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
	fwd_len3000_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 100, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_back_op_len3000( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
{

	__shared__ real_type_t<T> lds[3000];
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
	back_len3000_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 100, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_fwd_op_len3000( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[3000];
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
	fwd_len3000_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbInRe, gbInIm, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

//Kernel configuration: number of threads per thread block: 100, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_back_op_len3000( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	__shared__ real_type_t<T> lds[3000];
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
	back_len3000_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, gbInRe, gbInIm, iOffset, gbOutRe, gbOutIm, oOffset, lds, load_cb_fn, load_cb_data, load_cb_lds_bytes, store_cb_fn, store_cb_data);
}

