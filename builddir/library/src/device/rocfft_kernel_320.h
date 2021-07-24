#pragma once
#include "rocfft_butterfly_template.h"

////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb>
__device__ void
FwdPass0_len320(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R10) = bufIn[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R1) = bufIn[inOffset + ( 0 + me*2 + 0 + 32 )*stride_in];
	(*R11) = bufIn[inOffset + ( 0 + me*2 + 1 + 32 )*stride_in];
	(*R2) = bufIn[inOffset + ( 0 + me*2 + 0 + 64 )*stride_in];
	(*R12) = bufIn[inOffset + ( 0 + me*2 + 1 + 64 )*stride_in];
	(*R3) = bufIn[inOffset + ( 0 + me*2 + 0 + 96 )*stride_in];
	(*R13) = bufIn[inOffset + ( 0 + me*2 + 1 + 96 )*stride_in];
	(*R4) = bufIn[inOffset + ( 0 + me*2 + 0 + 128 )*stride_in];
	(*R14) = bufIn[inOffset + ( 0 + me*2 + 1 + 128 )*stride_in];
	(*R5) = bufIn[inOffset + ( 0 + me*2 + 0 + 160 )*stride_in];
	(*R15) = bufIn[inOffset + ( 0 + me*2 + 1 + 160 )*stride_in];
	(*R6) = bufIn[inOffset + ( 0 + me*2 + 0 + 192 )*stride_in];
	(*R16) = bufIn[inOffset + ( 0 + me*2 + 1 + 192 )*stride_in];
	(*R7) = bufIn[inOffset + ( 0 + me*2 + 0 + 224 )*stride_in];
	(*R17) = bufIn[inOffset + ( 0 + me*2 + 1 + 224 )*stride_in];
	(*R8) = bufIn[inOffset + ( 0 + me*2 + 0 + 256 )*stride_in];
	(*R18) = bufIn[inOffset + ( 0 + me*2 + 1 + 256 )*stride_in];
	(*R9) = bufIn[inOffset + ( 0 + me*2 + 0 + 288 )*stride_in];
	(*R19) = bufIn[inOffset + ( 0 + me*2 + 1 + 288 )*stride_in];
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
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 80 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 80 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 80 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 80 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 80 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 160 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 160 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 160 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 160 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 160 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 240 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 240 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 240 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 240 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 240 ) ];
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
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 80 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 80 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 80 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 80 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 80 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 160 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 160 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 160 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 160 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 160 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 240 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 240 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 240 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 240 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 240 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass0_len320(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 32 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 32 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 32 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 32 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 64 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 64 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 64 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 64 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 96 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 96 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 96 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 96 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 128 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 128 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 128 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 128 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 160 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 160 )*stride_in];
	(*R15).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 160 )*stride_in];
	(*R15).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 160 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 192 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 192 )*stride_in];
	(*R16).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 192 )*stride_in];
	(*R16).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 192 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 224 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 224 )*stride_in];
	(*R17).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 224 )*stride_in];
	(*R17).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 224 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 256 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 256 )*stride_in];
	(*R18).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 256 )*stride_in];
	(*R18).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 256 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 288 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 288 )*stride_in];
	(*R19).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 288 )*stride_in];
	(*R19).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 288 )*stride_in];
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
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 80 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 80 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 80 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 80 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 80 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 160 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 160 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 160 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 160 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 160 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 240 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 240 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 240 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 240 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 240 ) ];
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
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 80 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 80 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 80 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 80 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 80 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 160 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 160 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 160 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 160 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 160 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 240 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 240 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 240 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 240 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 240 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass1_len320(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{




	{
		T W = twiddles[9 + 3*((5*me + 0)%10) + 0];
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
		T W = twiddles[9 + 3*((5*me + 0)%10) + 1];
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
		T W = twiddles[9 + 3*((5*me + 0)%10) + 2];
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
		T W = twiddles[9 + 3*((5*me + 1)%10) + 0];
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
		T W = twiddles[9 + 3*((5*me + 1)%10) + 1];
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
		T W = twiddles[9 + 3*((5*me + 1)%10) + 2];
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
		T W = twiddles[9 + 3*((5*me + 2)%10) + 0];
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
		T W = twiddles[9 + 3*((5*me + 2)%10) + 1];
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
		T W = twiddles[9 + 3*((5*me + 2)%10) + 2];
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
		T W = twiddles[9 + 3*((5*me + 3)%10) + 0];
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
		T W = twiddles[9 + 3*((5*me + 3)%10) + 1];
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
		T W = twiddles[9 + 3*((5*me + 3)%10) + 2];
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
		T W = twiddles[9 + 3*((5*me + 4)%10) + 0];
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
		T W = twiddles[9 + 3*((5*me + 4)%10) + 1];
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
		T W = twiddles[9 + 3*((5*me + 4)%10) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	FwdRad4B1(R0, R1, R2, R3);
	FwdRad4B1(R4, R5, R6, R7);
	FwdRad4B1(R8, R9, R10, R11);
	FwdRad4B1(R12, R13, R14, R15);
	FwdRad4B1(R16, R17, R18, R19);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 10 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 20 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 30 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 0 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 10 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 20 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 30 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 0 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 10 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 20 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 30 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 10 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 20 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 30 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 0 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 10 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 20 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 30 ) ] = (*R19).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 80 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 80 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 80 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 80 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 80 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 160 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 160 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 160 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 160 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 160 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 240 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 240 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 240 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 240 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 240 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 10 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 20 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 30 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 0 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 10 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 20 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 30 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 0 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 10 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 20 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 30 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 10 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 20 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 30 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 0 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 10 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 20 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 30 ) ] = (*R19).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 80 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 80 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 80 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 80 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 80 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 160 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 160 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 160 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 160 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 160 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 240 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 240 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 240 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 240 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 240 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass2_len320(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{




	{
		T W = twiddles[39 + 3*((5*me + 0)%40) + 0];
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
		T W = twiddles[39 + 3*((5*me + 0)%40) + 1];
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
		T W = twiddles[39 + 3*((5*me + 0)%40) + 2];
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
		T W = twiddles[39 + 3*((5*me + 1)%40) + 0];
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
		T W = twiddles[39 + 3*((5*me + 1)%40) + 1];
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
		T W = twiddles[39 + 3*((5*me + 1)%40) + 2];
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
		T W = twiddles[39 + 3*((5*me + 2)%40) + 0];
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
		T W = twiddles[39 + 3*((5*me + 2)%40) + 1];
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
		T W = twiddles[39 + 3*((5*me + 2)%40) + 2];
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
		T W = twiddles[39 + 3*((5*me + 3)%40) + 0];
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
		T W = twiddles[39 + 3*((5*me + 3)%40) + 1];
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
		T W = twiddles[39 + 3*((5*me + 3)%40) + 2];
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
		T W = twiddles[39 + 3*((5*me + 4)%40) + 0];
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
		T W = twiddles[39 + 3*((5*me + 4)%40) + 1];
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
		T W = twiddles[39 + 3*((5*me + 4)%40) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	FwdRad4B1(R0, R1, R2, R3);
	FwdRad4B1(R4, R5, R6, R7);
	FwdRad4B1(R8, R9, R10, R11);
	FwdRad4B1(R12, R13, R14, R15);
	FwdRad4B1(R16, R17, R18, R19);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 40 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 80 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 120 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 0 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 40 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 80 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 120 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 0 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 40 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 80 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 120 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 40 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 80 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 120 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 0 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 40 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 80 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 120 ) ] = (*R19).x;
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
	(*R1).x = bufOutRe[outOffset + ( 0 + me*10 + 0 + 160 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*10 + 1 + 160 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*10 + 2 + 160 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*10 + 3 + 160 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*10 + 4 + 160 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*10 + 5 + 160 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*10 + 6 + 160 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*10 + 7 + 160 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*10 + 8 + 160 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*10 + 9 + 160 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 40 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 80 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 120 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 0 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 40 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 80 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 120 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 0 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 40 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 80 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 120 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 40 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 80 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 120 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 0 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 40 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 80 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 120 ) ] = (*R19).y;
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
	(*R1).y = bufOutIm[outOffset + ( 0 + me*10 + 0 + 160 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*10 + 1 + 160 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*10 + 2 + 160 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*10 + 3 + 160 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*10 + 4 + 160 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*10 + 5 + 160 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*10 + 6 + 160 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*10 + 7 + 160 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*10 + 8 + 160 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*10 + 9 + 160 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass3_len320(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{




	{
		T W = twiddles[159 + 1*((10*me + 0)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 1)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 2)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 3)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 4)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 5)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 6)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 7)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 8)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 9)%160) + 0];
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
	 //Optimization: coalescing into float4/double4 write
	if(sb == SB_UNIT) {
	vector4_type_t<T> *buff4g = (vector4_type_t<T>*)bufOut;
	
	buff4g[ 5*me + 0 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R0).x, (*R0).y, (*R2).x, (*R2).y) ;
	buff4g[ 5*me + 1 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R4).x, (*R4).y, (*R6).x, (*R6).y) ;
	buff4g[ 5*me + 2 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R8).x, (*R8).y, (*R10).x, (*R10).y) ;
	buff4g[ 5*me + 3 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R12).x, (*R12).y, (*R14).x, (*R14).y) ;
	buff4g[ 5*me + 4 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R16).x, (*R16).y, (*R18).x, (*R18).y) ;
	buff4g[ 5*me + 0 + 80 ] = lib_make_vector4< vector4_type_t<T> >((*R1).x, (*R1).y, (*R3).x, (*R3).y) ;
	buff4g[ 5*me + 1 + 80 ] = lib_make_vector4< vector4_type_t<T> >((*R5).x, (*R5).y, (*R7).x, (*R7).y) ;
	buff4g[ 5*me + 2 + 80 ] = lib_make_vector4< vector4_type_t<T> >((*R9).x, (*R9).y, (*R11).x, (*R11).y) ;
	buff4g[ 5*me + 3 + 80 ] = lib_make_vector4< vector4_type_t<T> >((*R13).x, (*R13).y, (*R15).x, (*R15).y) ;
	buff4g[ 5*me + 4 + 80 ] = lib_make_vector4< vector4_type_t<T> >((*R17).x, (*R17).y, (*R19).x, (*R19).y) ;
	}
	else{ // such optimization is not possible 
	bufOut[outOffset + ( 10*me + 0 + 0 )*stride_out] = (*R0);
	bufOut[outOffset + ( 10*me + 1 + 0 )*stride_out] = (*R2);
	bufOut[outOffset + ( 10*me + 2 + 0 )*stride_out] = (*R4);
	bufOut[outOffset + ( 10*me + 3 + 0 )*stride_out] = (*R6);
	bufOut[outOffset + ( 10*me + 4 + 0 )*stride_out] = (*R8);
	bufOut[outOffset + ( 10*me + 5 + 0 )*stride_out] = (*R10);
	bufOut[outOffset + ( 10*me + 6 + 0 )*stride_out] = (*R12);
	bufOut[outOffset + ( 10*me + 7 + 0 )*stride_out] = (*R14);
	bufOut[outOffset + ( 10*me + 8 + 0 )*stride_out] = (*R16);
	bufOut[outOffset + ( 10*me + 9 + 0 )*stride_out] = (*R18);
	bufOut[outOffset + ( 10*me + 0 + 160 )*stride_out] = (*R1);
	bufOut[outOffset + ( 10*me + 1 + 160 )*stride_out] = (*R3);
	bufOut[outOffset + ( 10*me + 2 + 160 )*stride_out] = (*R5);
	bufOut[outOffset + ( 10*me + 3 + 160 )*stride_out] = (*R7);
	bufOut[outOffset + ( 10*me + 4 + 160 )*stride_out] = (*R9);
	bufOut[outOffset + ( 10*me + 5 + 160 )*stride_out] = (*R11);
	bufOut[outOffset + ( 10*me + 6 + 160 )*stride_out] = (*R13);
	bufOut[outOffset + ( 10*me + 7 + 160 )*stride_out] = (*R15);
	bufOut[outOffset + ( 10*me + 8 + 160 )*stride_out] = (*R17);
	bufOut[outOffset + ( 10*me + 9 + 160 )*stride_out] = (*R19);
	}
	}

}

template <typename T, StrideBin sb>
__device__ void
FwdPass3_len320(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{




	{
		T W = twiddles[159 + 1*((10*me + 0)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 1)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 2)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 3)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 4)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 5)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 6)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 7)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 8)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 9)%160) + 0];
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
	bufOutRe[outOffset + ( 10*me + 0 + 160 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 10*me + 0 + 160 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 10*me + 1 + 160 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 10*me + 1 + 160 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 10*me + 2 + 160 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 10*me + 2 + 160 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 10*me + 3 + 160 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 10*me + 3 + 160 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 10*me + 4 + 160 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 10*me + 4 + 160 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 10*me + 5 + 160 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 10*me + 5 + 160 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 10*me + 6 + 160 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 10*me + 6 + 160 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 10*me + 7 + 160 )*stride_out] = (*R15).x;
	bufOutIm[outOffset + ( 10*me + 7 + 160 )*stride_out] = (*R15).y;
	bufOutRe[outOffset + ( 10*me + 8 + 160 )*stride_out] = (*R17).x;
	bufOutIm[outOffset + ( 10*me + 8 + 160 )*stride_out] = (*R17).y;
	bufOutRe[outOffset + ( 10*me + 9 + 160 )*stride_out] = (*R19).x;
	bufOutIm[outOffset + ( 10*me + 9 + 160 )*stride_out] = (*R19).y;
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass0_len320(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R10) = bufIn[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R1) = bufIn[inOffset + ( 0 + me*2 + 0 + 32 )*stride_in];
	(*R11) = bufIn[inOffset + ( 0 + me*2 + 1 + 32 )*stride_in];
	(*R2) = bufIn[inOffset + ( 0 + me*2 + 0 + 64 )*stride_in];
	(*R12) = bufIn[inOffset + ( 0 + me*2 + 1 + 64 )*stride_in];
	(*R3) = bufIn[inOffset + ( 0 + me*2 + 0 + 96 )*stride_in];
	(*R13) = bufIn[inOffset + ( 0 + me*2 + 1 + 96 )*stride_in];
	(*R4) = bufIn[inOffset + ( 0 + me*2 + 0 + 128 )*stride_in];
	(*R14) = bufIn[inOffset + ( 0 + me*2 + 1 + 128 )*stride_in];
	(*R5) = bufIn[inOffset + ( 0 + me*2 + 0 + 160 )*stride_in];
	(*R15) = bufIn[inOffset + ( 0 + me*2 + 1 + 160 )*stride_in];
	(*R6) = bufIn[inOffset + ( 0 + me*2 + 0 + 192 )*stride_in];
	(*R16) = bufIn[inOffset + ( 0 + me*2 + 1 + 192 )*stride_in];
	(*R7) = bufIn[inOffset + ( 0 + me*2 + 0 + 224 )*stride_in];
	(*R17) = bufIn[inOffset + ( 0 + me*2 + 1 + 224 )*stride_in];
	(*R8) = bufIn[inOffset + ( 0 + me*2 + 0 + 256 )*stride_in];
	(*R18) = bufIn[inOffset + ( 0 + me*2 + 1 + 256 )*stride_in];
	(*R9) = bufIn[inOffset + ( 0 + me*2 + 0 + 288 )*stride_in];
	(*R19) = bufIn[inOffset + ( 0 + me*2 + 1 + 288 )*stride_in];
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
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 80 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 80 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 80 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 80 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 80 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 160 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 160 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 160 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 160 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 160 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 240 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 240 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 240 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 240 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 240 ) ];
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
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 80 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 80 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 80 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 80 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 80 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 160 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 160 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 160 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 160 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 160 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 240 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 240 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 240 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 240 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 240 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass0_len320(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 32 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 32 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 32 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 32 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 64 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 64 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 64 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 64 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 96 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 96 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 96 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 96 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 128 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 128 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 128 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 128 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 160 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 160 )*stride_in];
	(*R15).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 160 )*stride_in];
	(*R15).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 160 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 192 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 192 )*stride_in];
	(*R16).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 192 )*stride_in];
	(*R16).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 192 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 224 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 224 )*stride_in];
	(*R17).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 224 )*stride_in];
	(*R17).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 224 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 256 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 256 )*stride_in];
	(*R18).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 256 )*stride_in];
	(*R18).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 256 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 288 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 288 )*stride_in];
	(*R19).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 288 )*stride_in];
	(*R19).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 288 )*stride_in];
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
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 80 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 80 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 80 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 80 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 80 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 160 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 160 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 160 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 160 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 160 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 240 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 240 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 240 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 240 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 240 ) ];
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
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 80 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 80 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 80 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 80 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 80 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 160 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 160 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 160 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 160 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 160 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 240 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 240 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 240 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 240 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 240 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass1_len320(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{




	{
		T W = twiddles[9 + 3*((5*me + 0)%10) + 0];
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
		T W = twiddles[9 + 3*((5*me + 0)%10) + 1];
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
		T W = twiddles[9 + 3*((5*me + 0)%10) + 2];
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
		T W = twiddles[9 + 3*((5*me + 1)%10) + 0];
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
		T W = twiddles[9 + 3*((5*me + 1)%10) + 1];
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
		T W = twiddles[9 + 3*((5*me + 1)%10) + 2];
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
		T W = twiddles[9 + 3*((5*me + 2)%10) + 0];
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
		T W = twiddles[9 + 3*((5*me + 2)%10) + 1];
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
		T W = twiddles[9 + 3*((5*me + 2)%10) + 2];
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
		T W = twiddles[9 + 3*((5*me + 3)%10) + 0];
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
		T W = twiddles[9 + 3*((5*me + 3)%10) + 1];
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
		T W = twiddles[9 + 3*((5*me + 3)%10) + 2];
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
		T W = twiddles[9 + 3*((5*me + 4)%10) + 0];
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
		T W = twiddles[9 + 3*((5*me + 4)%10) + 1];
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
		T W = twiddles[9 + 3*((5*me + 4)%10) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	InvRad4B1(R0, R1, R2, R3);
	InvRad4B1(R4, R5, R6, R7);
	InvRad4B1(R8, R9, R10, R11);
	InvRad4B1(R12, R13, R14, R15);
	InvRad4B1(R16, R17, R18, R19);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 10 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 20 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 30 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 0 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 10 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 20 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 30 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 0 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 10 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 20 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 30 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 10 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 20 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 30 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 0 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 10 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 20 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 30 ) ] = (*R19).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 80 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 80 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 80 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 80 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 80 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 160 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 160 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 160 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 160 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 160 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 240 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 240 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 240 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 240 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 240 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 10 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 20 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/10)*40 + (5*me + 0)%10 + 30 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 0 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 10 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 20 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/10)*40 + (5*me + 1)%10 + 30 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 0 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 10 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 20 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/10)*40 + (5*me + 2)%10 + 30 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 10 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 20 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/10)*40 + (5*me + 3)%10 + 30 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 0 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 10 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 20 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 4)/10)*40 + (5*me + 4)%10 + 30 ) ] = (*R19).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 80 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 80 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 80 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 80 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 80 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 160 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 160 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 160 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 160 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 160 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 240 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 240 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 240 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 240 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 240 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass2_len320(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{




	{
		T W = twiddles[39 + 3*((5*me + 0)%40) + 0];
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
		T W = twiddles[39 + 3*((5*me + 0)%40) + 1];
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
		T W = twiddles[39 + 3*((5*me + 0)%40) + 2];
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
		T W = twiddles[39 + 3*((5*me + 1)%40) + 0];
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
		T W = twiddles[39 + 3*((5*me + 1)%40) + 1];
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
		T W = twiddles[39 + 3*((5*me + 1)%40) + 2];
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
		T W = twiddles[39 + 3*((5*me + 2)%40) + 0];
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
		T W = twiddles[39 + 3*((5*me + 2)%40) + 1];
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
		T W = twiddles[39 + 3*((5*me + 2)%40) + 2];
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
		T W = twiddles[39 + 3*((5*me + 3)%40) + 0];
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
		T W = twiddles[39 + 3*((5*me + 3)%40) + 1];
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
		T W = twiddles[39 + 3*((5*me + 3)%40) + 2];
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
		T W = twiddles[39 + 3*((5*me + 4)%40) + 0];
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
		T W = twiddles[39 + 3*((5*me + 4)%40) + 1];
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
		T W = twiddles[39 + 3*((5*me + 4)%40) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	InvRad4B1(R0, R1, R2, R3);
	InvRad4B1(R4, R5, R6, R7);
	InvRad4B1(R8, R9, R10, R11);
	InvRad4B1(R12, R13, R14, R15);
	InvRad4B1(R16, R17, R18, R19);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 40 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 80 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 120 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 0 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 40 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 80 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 120 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 0 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 40 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 80 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 120 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 40 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 80 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 120 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 0 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 40 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 80 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 120 ) ] = (*R19).x;
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
	(*R1).x = bufOutRe[outOffset + ( 0 + me*10 + 0 + 160 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*10 + 1 + 160 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*10 + 2 + 160 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*10 + 3 + 160 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*10 + 4 + 160 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*10 + 5 + 160 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*10 + 6 + 160 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*10 + 7 + 160 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*10 + 8 + 160 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*10 + 9 + 160 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 40 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 80 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/40)*160 + (5*me + 0)%40 + 120 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 0 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 40 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 80 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/40)*160 + (5*me + 1)%40 + 120 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 0 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 40 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 80 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 2)/40)*160 + (5*me + 2)%40 + 120 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 40 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 80 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 3)/40)*160 + (5*me + 3)%40 + 120 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 0 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 40 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 80 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 4)/40)*160 + (5*me + 4)%40 + 120 ) ] = (*R19).y;
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
	(*R1).y = bufOutIm[outOffset + ( 0 + me*10 + 0 + 160 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*10 + 1 + 160 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*10 + 2 + 160 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*10 + 3 + 160 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*10 + 4 + 160 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*10 + 5 + 160 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*10 + 6 + 160 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*10 + 7 + 160 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*10 + 8 + 160 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*10 + 9 + 160 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass3_len320(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{




	{
		T W = twiddles[159 + 1*((10*me + 0)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 1)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 2)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 3)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 4)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 5)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 6)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 7)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 8)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 9)%160) + 0];
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
	 //Optimization: coalescing into float4/double4 write
	if(sb == SB_UNIT) {
	vector4_type_t<T> *buff4g = (vector4_type_t<T>*)bufOut;
	
	buff4g[ 5*me + 0 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R0).x, (*R0).y, (*R2).x, (*R2).y) ;
	buff4g[ 5*me + 1 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R4).x, (*R4).y, (*R6).x, (*R6).y) ;
	buff4g[ 5*me + 2 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R8).x, (*R8).y, (*R10).x, (*R10).y) ;
	buff4g[ 5*me + 3 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R12).x, (*R12).y, (*R14).x, (*R14).y) ;
	buff4g[ 5*me + 4 + 0 ] = lib_make_vector4< vector4_type_t<T> >((*R16).x, (*R16).y, (*R18).x, (*R18).y) ;
	buff4g[ 5*me + 0 + 80 ] = lib_make_vector4< vector4_type_t<T> >((*R1).x, (*R1).y, (*R3).x, (*R3).y) ;
	buff4g[ 5*me + 1 + 80 ] = lib_make_vector4< vector4_type_t<T> >((*R5).x, (*R5).y, (*R7).x, (*R7).y) ;
	buff4g[ 5*me + 2 + 80 ] = lib_make_vector4< vector4_type_t<T> >((*R9).x, (*R9).y, (*R11).x, (*R11).y) ;
	buff4g[ 5*me + 3 + 80 ] = lib_make_vector4< vector4_type_t<T> >((*R13).x, (*R13).y, (*R15).x, (*R15).y) ;
	buff4g[ 5*me + 4 + 80 ] = lib_make_vector4< vector4_type_t<T> >((*R17).x, (*R17).y, (*R19).x, (*R19).y) ;
	}
	else{ // such optimization is not possible 
	bufOut[outOffset + ( 10*me + 0 + 0 )*stride_out] = (*R0);
	bufOut[outOffset + ( 10*me + 1 + 0 )*stride_out] = (*R2);
	bufOut[outOffset + ( 10*me + 2 + 0 )*stride_out] = (*R4);
	bufOut[outOffset + ( 10*me + 3 + 0 )*stride_out] = (*R6);
	bufOut[outOffset + ( 10*me + 4 + 0 )*stride_out] = (*R8);
	bufOut[outOffset + ( 10*me + 5 + 0 )*stride_out] = (*R10);
	bufOut[outOffset + ( 10*me + 6 + 0 )*stride_out] = (*R12);
	bufOut[outOffset + ( 10*me + 7 + 0 )*stride_out] = (*R14);
	bufOut[outOffset + ( 10*me + 8 + 0 )*stride_out] = (*R16);
	bufOut[outOffset + ( 10*me + 9 + 0 )*stride_out] = (*R18);
	bufOut[outOffset + ( 10*me + 0 + 160 )*stride_out] = (*R1);
	bufOut[outOffset + ( 10*me + 1 + 160 )*stride_out] = (*R3);
	bufOut[outOffset + ( 10*me + 2 + 160 )*stride_out] = (*R5);
	bufOut[outOffset + ( 10*me + 3 + 160 )*stride_out] = (*R7);
	bufOut[outOffset + ( 10*me + 4 + 160 )*stride_out] = (*R9);
	bufOut[outOffset + ( 10*me + 5 + 160 )*stride_out] = (*R11);
	bufOut[outOffset + ( 10*me + 6 + 160 )*stride_out] = (*R13);
	bufOut[outOffset + ( 10*me + 7 + 160 )*stride_out] = (*R15);
	bufOut[outOffset + ( 10*me + 8 + 160 )*stride_out] = (*R17);
	bufOut[outOffset + ( 10*me + 9 + 160 )*stride_out] = (*R19);
	}
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass3_len320(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{




	{
		T W = twiddles[159 + 1*((10*me + 0)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 1)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 2)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 3)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 4)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 5)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 6)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 7)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 8)%160) + 0];
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
		T W = twiddles[159 + 1*((10*me + 9)%160) + 0];
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
	bufOutRe[outOffset + ( 10*me + 0 + 160 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 10*me + 0 + 160 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 10*me + 1 + 160 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 10*me + 1 + 160 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 10*me + 2 + 160 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 10*me + 2 + 160 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 10*me + 3 + 160 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 10*me + 3 + 160 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 10*me + 4 + 160 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 10*me + 4 + 160 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 10*me + 5 + 160 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 10*me + 5 + 160 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 10*me + 6 + 160 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 10*me + 6 + 160 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 10*me + 7 + 160 )*stride_out] = (*R15).x;
	bufOutIm[outOffset + ( 10*me + 7 + 160 )*stride_out] = (*R15).y;
	bufOutRe[outOffset + ( 10*me + 8 + 160 )*stride_out] = (*R17).x;
	bufOutIm[outOffset + ( 10*me + 8 + 160 )*stride_out] = (*R17).y;
	bufOutRe[outOffset + ( 10*me + 9 + 160 )*stride_out] = (*R19).x;
	bufOutIm[outOffset + ( 10*me + 9 + 160 )*stride_out] = (*R19).y;
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb>
__device__ void 
fwd_len320_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass1_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass2_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass3_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
}

template <typename T, StrideBin sb>
__device__ void 
back_len320_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass1_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass2_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass3_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len320_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass1_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass2_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass3_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
}

template <typename T, StrideBin sb>
__device__ void 
back_len320_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass1_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass2_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass3_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len320_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass1_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass2_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass3_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
}

template <typename T, StrideBin sb>
__device__ void 
back_len320_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass1_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass2_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass3_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len320_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass1_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass2_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	FwdPass3_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
}

template <typename T, StrideBin sb>
__device__ void 
back_len320_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass1_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass2_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	InvPass3_len320<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 64, maximum transforms: 4, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(64)
fft_fwd_ip_len320( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1280), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len320>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*4)*16) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*4 + (me/16));
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
	fwd_len320_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%16, (me/16)*320, lwb, lwb, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 4, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(64)
fft_back_ip_len320( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1280), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len320>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*4)*16) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*4 + (me/16));
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
	back_len320_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%16, (me/16)*320, lwb, lwb, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 4, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(64)
fft_fwd_ip_len320( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1280), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len320>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*4)*16) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*4 + (me/16));
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
	fwd_len320_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%16, (me/16)*320, lwbRe, lwbIm, lwbRe, lwbIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 4, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(64)
fft_back_ip_len320( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1280), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len320>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*4)*16) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*4 + (me/16));
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
	back_len320_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%16, (me/16)*320, lwbRe, lwbIm, lwbRe, lwbIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 4, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(64)
fft_fwd_op_len320( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1280), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len320>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*4)*16) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*4 + (me/16));
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
	fwd_len320_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*320, lwbIn, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 4, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(64)
fft_back_op_len320( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1280), cgh);
	cgh.parallel_for<class kern_fft_back_op_len320>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*4)*16) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*4 + (me/16));
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
	back_len320_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*320, lwbIn, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 4, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(64)
fft_fwd_op_len320( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1280), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len320>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*4)*16) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*4 + (me/16));
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
	fwd_len320_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*320, lwbIn, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 4, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(64)
fft_back_op_len320( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1280), cgh);
	cgh.parallel_for<class kern_fft_back_op_len320>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*4)*16) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*4 + (me/16));
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
	back_len320_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*320, lwbIn, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 4, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(64)
fft_fwd_op_len320( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1280), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len320>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*4)*16) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*4 + (me/16));
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
	fwd_len320_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*320, lwbInRe, lwbInIm, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 4, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(64)
fft_back_op_len320( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1280), cgh);
	cgh.parallel_for<class kern_fft_back_op_len320>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*4)*16) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*4 + (me/16));
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
	back_len320_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*320, lwbInRe, lwbInIm, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 4, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(64)
fft_fwd_op_len320( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1280), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len320>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*4)*16) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*4 + (me/16));
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
	fwd_len320_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*320, lwbInRe, lwbInIm, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 4, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(64)
fft_back_op_len320( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(1280), cgh);
	cgh.parallel_for<class kern_fft_back_op_len320>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*4)*16) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*4 + (me/16));
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
	back_len320_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*320, lwbInRe, lwbInIm, lwbOutRe, lwbOutIm, lds);
});
});
}

