#pragma once
#include "rocfft_butterfly_template.h"

////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb>
__device__ void
FwdPass0_len104(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R13) = bufIn[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R1) = bufIn[inOffset + ( 0 + me*2 + 0 + 8 )*stride_in];
	(*R14) = bufIn[inOffset + ( 0 + me*2 + 1 + 8 )*stride_in];
	(*R2) = bufIn[inOffset + ( 0 + me*2 + 0 + 16 )*stride_in];
	(*R15) = bufIn[inOffset + ( 0 + me*2 + 1 + 16 )*stride_in];
	(*R3) = bufIn[inOffset + ( 0 + me*2 + 0 + 24 )*stride_in];
	(*R16) = bufIn[inOffset + ( 0 + me*2 + 1 + 24 )*stride_in];
	(*R4) = bufIn[inOffset + ( 0 + me*2 + 0 + 32 )*stride_in];
	(*R17) = bufIn[inOffset + ( 0 + me*2 + 1 + 32 )*stride_in];
	(*R5) = bufIn[inOffset + ( 0 + me*2 + 0 + 40 )*stride_in];
	(*R18) = bufIn[inOffset + ( 0 + me*2 + 1 + 40 )*stride_in];
	(*R6) = bufIn[inOffset + ( 0 + me*2 + 0 + 48 )*stride_in];
	(*R19) = bufIn[inOffset + ( 0 + me*2 + 1 + 48 )*stride_in];
	(*R7) = bufIn[inOffset + ( 0 + me*2 + 0 + 56 )*stride_in];
	(*R20) = bufIn[inOffset + ( 0 + me*2 + 1 + 56 )*stride_in];
	(*R8) = bufIn[inOffset + ( 0 + me*2 + 0 + 64 )*stride_in];
	(*R21) = bufIn[inOffset + ( 0 + me*2 + 1 + 64 )*stride_in];
	(*R9) = bufIn[inOffset + ( 0 + me*2 + 0 + 72 )*stride_in];
	(*R22) = bufIn[inOffset + ( 0 + me*2 + 1 + 72 )*stride_in];
	(*R10) = bufIn[inOffset + ( 0 + me*2 + 0 + 80 )*stride_in];
	(*R23) = bufIn[inOffset + ( 0 + me*2 + 1 + 80 )*stride_in];
	(*R11) = bufIn[inOffset + ( 0 + me*2 + 0 + 88 )*stride_in];
	(*R24) = bufIn[inOffset + ( 0 + me*2 + 1 + 88 )*stride_in];
	(*R12) = bufIn[inOffset + ( 0 + me*2 + 0 + 96 )*stride_in];
	(*R25) = bufIn[inOffset + ( 0 + me*2 + 1 + 96 )*stride_in];
	}



	FwdRad13B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12);
	FwdRad13B1(R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25);


	if(rw)
	{
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 10 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 11 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 12 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 0 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 1 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 2 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 3 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 4 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 5 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 6 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 7 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 8 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 9 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 10 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 11 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 12 ) ] = (*R25).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 10 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 11 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 12 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 0 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 1 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 2 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 3 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 4 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 5 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 6 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 7 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 8 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 9 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 10 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 11 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 12 ) ] = (*R25).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass0_len104(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 8 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 8 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 8 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 8 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 16 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 16 )*stride_in];
	(*R15).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 16 )*stride_in];
	(*R15).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 16 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 24 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 24 )*stride_in];
	(*R16).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 24 )*stride_in];
	(*R16).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 24 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 32 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 32 )*stride_in];
	(*R17).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 32 )*stride_in];
	(*R17).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 32 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 40 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 40 )*stride_in];
	(*R18).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 40 )*stride_in];
	(*R18).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 40 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 48 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 48 )*stride_in];
	(*R19).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 48 )*stride_in];
	(*R19).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 48 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 56 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 56 )*stride_in];
	(*R20).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 56 )*stride_in];
	(*R20).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 56 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 64 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 64 )*stride_in];
	(*R21).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 64 )*stride_in];
	(*R21).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 64 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 72 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 72 )*stride_in];
	(*R22).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 72 )*stride_in];
	(*R22).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 72 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 80 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 80 )*stride_in];
	(*R23).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 80 )*stride_in];
	(*R23).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 80 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 88 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 88 )*stride_in];
	(*R24).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 88 )*stride_in];
	(*R24).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 88 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 96 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 96 )*stride_in];
	(*R25).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 96 )*stride_in];
	(*R25).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 96 )*stride_in];
	}



	FwdRad13B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12);
	FwdRad13B1(R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25);


	if(rw)
	{
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 10 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 11 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 12 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 0 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 1 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 2 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 3 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 4 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 5 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 6 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 7 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 8 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 9 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 10 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 11 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 12 ) ] = (*R25).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 10 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 11 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 12 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 0 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 1 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 2 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 3 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 4 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 5 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 6 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 7 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 8 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 9 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 10 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 11 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 12 ) ] = (*R25).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass1_len104(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25)
{




	{
		T W = twiddles[12 + 1*((13*me + 0)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 1)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 2)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 3)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 4)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 5)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 6)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 7)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 8)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 9)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 10)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 11)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 12)%13) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
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
	FwdRad2B1(R20, R21);
	FwdRad2B1(R22, R23);
	FwdRad2B1(R24, R25);


	if(rw)
	{
	bufOutRe[outOffset + ( ((13*me + 0)/13)*26 + (13*me + 0)%13 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((13*me + 0)/13)*26 + (13*me + 0)%13 + 13 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((13*me + 1)/13)*26 + (13*me + 1)%13 + 0 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((13*me + 1)/13)*26 + (13*me + 1)%13 + 13 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((13*me + 2)/13)*26 + (13*me + 2)%13 + 0 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((13*me + 2)/13)*26 + (13*me + 2)%13 + 13 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((13*me + 3)/13)*26 + (13*me + 3)%13 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((13*me + 3)/13)*26 + (13*me + 3)%13 + 13 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((13*me + 4)/13)*26 + (13*me + 4)%13 + 0 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((13*me + 4)/13)*26 + (13*me + 4)%13 + 13 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((13*me + 5)/13)*26 + (13*me + 5)%13 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((13*me + 5)/13)*26 + (13*me + 5)%13 + 13 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((13*me + 6)/13)*26 + (13*me + 6)%13 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((13*me + 6)/13)*26 + (13*me + 6)%13 + 13 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((13*me + 7)/13)*26 + (13*me + 7)%13 + 0 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((13*me + 7)/13)*26 + (13*me + 7)%13 + 13 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((13*me + 8)/13)*26 + (13*me + 8)%13 + 0 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((13*me + 8)/13)*26 + (13*me + 8)%13 + 13 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((13*me + 9)/13)*26 + (13*me + 9)%13 + 0 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((13*me + 9)/13)*26 + (13*me + 9)%13 + 13 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((13*me + 10)/13)*26 + (13*me + 10)%13 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((13*me + 10)/13)*26 + (13*me + 10)%13 + 13 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((13*me + 11)/13)*26 + (13*me + 11)%13 + 0 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((13*me + 11)/13)*26 + (13*me + 11)%13 + 13 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((13*me + 12)/13)*26 + (13*me + 12)%13 + 0 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((13*me + 12)/13)*26 + (13*me + 12)%13 + 13 ) ] = (*R25).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((13*me + 0)/13)*26 + (13*me + 0)%13 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((13*me + 0)/13)*26 + (13*me + 0)%13 + 13 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((13*me + 1)/13)*26 + (13*me + 1)%13 + 0 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((13*me + 1)/13)*26 + (13*me + 1)%13 + 13 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((13*me + 2)/13)*26 + (13*me + 2)%13 + 0 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((13*me + 2)/13)*26 + (13*me + 2)%13 + 13 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((13*me + 3)/13)*26 + (13*me + 3)%13 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((13*me + 3)/13)*26 + (13*me + 3)%13 + 13 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((13*me + 4)/13)*26 + (13*me + 4)%13 + 0 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((13*me + 4)/13)*26 + (13*me + 4)%13 + 13 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((13*me + 5)/13)*26 + (13*me + 5)%13 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((13*me + 5)/13)*26 + (13*me + 5)%13 + 13 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((13*me + 6)/13)*26 + (13*me + 6)%13 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((13*me + 6)/13)*26 + (13*me + 6)%13 + 13 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((13*me + 7)/13)*26 + (13*me + 7)%13 + 0 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((13*me + 7)/13)*26 + (13*me + 7)%13 + 13 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((13*me + 8)/13)*26 + (13*me + 8)%13 + 0 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((13*me + 8)/13)*26 + (13*me + 8)%13 + 13 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((13*me + 9)/13)*26 + (13*me + 9)%13 + 0 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((13*me + 9)/13)*26 + (13*me + 9)%13 + 13 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((13*me + 10)/13)*26 + (13*me + 10)%13 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((13*me + 10)/13)*26 + (13*me + 10)%13 + 13 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((13*me + 11)/13)*26 + (13*me + 11)%13 + 0 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((13*me + 11)/13)*26 + (13*me + 11)%13 + 13 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((13*me + 12)/13)*26 + (13*me + 12)%13 + 0 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((13*me + 12)/13)*26 + (13*me + 12)%13 + 13 ) ] = (*R25).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass2_len104(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25)
{




	{
		T W = twiddles[25 + 1*((13*me + 0)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 1)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 2)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 3)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 4)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 5)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 6)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 7)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 8)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 9)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 10)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 11)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 12)%26) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
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
	FwdRad2B1(R20, R21);
	FwdRad2B1(R22, R23);
	FwdRad2B1(R24, R25);


	if(rw)
	{
	bufOutRe[outOffset + ( ((13*me + 0)/26)*52 + (13*me + 0)%26 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((13*me + 0)/26)*52 + (13*me + 0)%26 + 26 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((13*me + 1)/26)*52 + (13*me + 1)%26 + 0 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((13*me + 1)/26)*52 + (13*me + 1)%26 + 26 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((13*me + 2)/26)*52 + (13*me + 2)%26 + 0 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((13*me + 2)/26)*52 + (13*me + 2)%26 + 26 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((13*me + 3)/26)*52 + (13*me + 3)%26 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((13*me + 3)/26)*52 + (13*me + 3)%26 + 26 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((13*me + 4)/26)*52 + (13*me + 4)%26 + 0 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((13*me + 4)/26)*52 + (13*me + 4)%26 + 26 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((13*me + 5)/26)*52 + (13*me + 5)%26 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((13*me + 5)/26)*52 + (13*me + 5)%26 + 26 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((13*me + 6)/26)*52 + (13*me + 6)%26 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((13*me + 6)/26)*52 + (13*me + 6)%26 + 26 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((13*me + 7)/26)*52 + (13*me + 7)%26 + 0 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((13*me + 7)/26)*52 + (13*me + 7)%26 + 26 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((13*me + 8)/26)*52 + (13*me + 8)%26 + 0 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((13*me + 8)/26)*52 + (13*me + 8)%26 + 26 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((13*me + 9)/26)*52 + (13*me + 9)%26 + 0 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((13*me + 9)/26)*52 + (13*me + 9)%26 + 26 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((13*me + 10)/26)*52 + (13*me + 10)%26 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((13*me + 10)/26)*52 + (13*me + 10)%26 + 26 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((13*me + 11)/26)*52 + (13*me + 11)%26 + 0 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((13*me + 11)/26)*52 + (13*me + 11)%26 + 26 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((13*me + 12)/26)*52 + (13*me + 12)%26 + 0 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((13*me + 12)/26)*52 + (13*me + 12)%26 + 26 ) ] = (*R25).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((13*me + 0)/26)*52 + (13*me + 0)%26 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((13*me + 0)/26)*52 + (13*me + 0)%26 + 26 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((13*me + 1)/26)*52 + (13*me + 1)%26 + 0 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((13*me + 1)/26)*52 + (13*me + 1)%26 + 26 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((13*me + 2)/26)*52 + (13*me + 2)%26 + 0 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((13*me + 2)/26)*52 + (13*me + 2)%26 + 26 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((13*me + 3)/26)*52 + (13*me + 3)%26 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((13*me + 3)/26)*52 + (13*me + 3)%26 + 26 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((13*me + 4)/26)*52 + (13*me + 4)%26 + 0 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((13*me + 4)/26)*52 + (13*me + 4)%26 + 26 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((13*me + 5)/26)*52 + (13*me + 5)%26 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((13*me + 5)/26)*52 + (13*me + 5)%26 + 26 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((13*me + 6)/26)*52 + (13*me + 6)%26 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((13*me + 6)/26)*52 + (13*me + 6)%26 + 26 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((13*me + 7)/26)*52 + (13*me + 7)%26 + 0 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((13*me + 7)/26)*52 + (13*me + 7)%26 + 26 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((13*me + 8)/26)*52 + (13*me + 8)%26 + 0 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((13*me + 8)/26)*52 + (13*me + 8)%26 + 26 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((13*me + 9)/26)*52 + (13*me + 9)%26 + 0 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((13*me + 9)/26)*52 + (13*me + 9)%26 + 26 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((13*me + 10)/26)*52 + (13*me + 10)%26 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((13*me + 10)/26)*52 + (13*me + 10)%26 + 26 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((13*me + 11)/26)*52 + (13*me + 11)%26 + 0 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((13*me + 11)/26)*52 + (13*me + 11)%26 + 26 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((13*me + 12)/26)*52 + (13*me + 12)%26 + 0 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((13*me + 12)/26)*52 + (13*me + 12)%26 + 26 ) ] = (*R25).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass3_len104(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25)
{




	{
		T W = twiddles[51 + 1*((13*me + 0)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 1)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 2)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 3)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 4)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 5)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 6)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 7)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 8)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 9)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 10)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 11)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 12)%52) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
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
	FwdRad2B1(R20, R21);
	FwdRad2B1(R22, R23);
	FwdRad2B1(R24, R25);


	if(rw)
	{
	bufOut[outOffset + ( 13*me + 0 + 0 )*stride_out] = (*R0);
	bufOut[outOffset + ( 13*me + 1 + 0 )*stride_out] = (*R2);
	bufOut[outOffset + ( 13*me + 2 + 0 )*stride_out] = (*R4);
	bufOut[outOffset + ( 13*me + 3 + 0 )*stride_out] = (*R6);
	bufOut[outOffset + ( 13*me + 4 + 0 )*stride_out] = (*R8);
	bufOut[outOffset + ( 13*me + 5 + 0 )*stride_out] = (*R10);
	bufOut[outOffset + ( 13*me + 6 + 0 )*stride_out] = (*R12);
	bufOut[outOffset + ( 13*me + 7 + 0 )*stride_out] = (*R14);
	bufOut[outOffset + ( 13*me + 8 + 0 )*stride_out] = (*R16);
	bufOut[outOffset + ( 13*me + 9 + 0 )*stride_out] = (*R18);
	bufOut[outOffset + ( 13*me + 10 + 0 )*stride_out] = (*R20);
	bufOut[outOffset + ( 13*me + 11 + 0 )*stride_out] = (*R22);
	bufOut[outOffset + ( 13*me + 12 + 0 )*stride_out] = (*R24);
	bufOut[outOffset + ( 13*me + 0 + 52 )*stride_out] = (*R1);
	bufOut[outOffset + ( 13*me + 1 + 52 )*stride_out] = (*R3);
	bufOut[outOffset + ( 13*me + 2 + 52 )*stride_out] = (*R5);
	bufOut[outOffset + ( 13*me + 3 + 52 )*stride_out] = (*R7);
	bufOut[outOffset + ( 13*me + 4 + 52 )*stride_out] = (*R9);
	bufOut[outOffset + ( 13*me + 5 + 52 )*stride_out] = (*R11);
	bufOut[outOffset + ( 13*me + 6 + 52 )*stride_out] = (*R13);
	bufOut[outOffset + ( 13*me + 7 + 52 )*stride_out] = (*R15);
	bufOut[outOffset + ( 13*me + 8 + 52 )*stride_out] = (*R17);
	bufOut[outOffset + ( 13*me + 9 + 52 )*stride_out] = (*R19);
	bufOut[outOffset + ( 13*me + 10 + 52 )*stride_out] = (*R21);
	bufOut[outOffset + ( 13*me + 11 + 52 )*stride_out] = (*R23);
	bufOut[outOffset + ( 13*me + 12 + 52 )*stride_out] = (*R25);
	}

}

template <typename T, StrideBin sb>
__device__ void
FwdPass3_len104(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25)
{




	{
		T W = twiddles[51 + 1*((13*me + 0)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 1)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 2)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 3)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 4)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 5)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 6)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 7)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 8)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 9)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 10)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 11)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 12)%52) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
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
	FwdRad2B1(R20, R21);
	FwdRad2B1(R22, R23);
	FwdRad2B1(R24, R25);


	if(rw)
	{
	bufOutRe[outOffset + ( 13*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 13*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 13*me + 1 + 0 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 13*me + 1 + 0 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 13*me + 2 + 0 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 13*me + 2 + 0 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 13*me + 3 + 0 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 13*me + 3 + 0 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 13*me + 4 + 0 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 13*me + 4 + 0 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 13*me + 5 + 0 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 13*me + 5 + 0 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 13*me + 6 + 0 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 13*me + 6 + 0 )*stride_out] = (*R12).y;
	bufOutRe[outOffset + ( 13*me + 7 + 0 )*stride_out] = (*R14).x;
	bufOutIm[outOffset + ( 13*me + 7 + 0 )*stride_out] = (*R14).y;
	bufOutRe[outOffset + ( 13*me + 8 + 0 )*stride_out] = (*R16).x;
	bufOutIm[outOffset + ( 13*me + 8 + 0 )*stride_out] = (*R16).y;
	bufOutRe[outOffset + ( 13*me + 9 + 0 )*stride_out] = (*R18).x;
	bufOutIm[outOffset + ( 13*me + 9 + 0 )*stride_out] = (*R18).y;
	bufOutRe[outOffset + ( 13*me + 10 + 0 )*stride_out] = (*R20).x;
	bufOutIm[outOffset + ( 13*me + 10 + 0 )*stride_out] = (*R20).y;
	bufOutRe[outOffset + ( 13*me + 11 + 0 )*stride_out] = (*R22).x;
	bufOutIm[outOffset + ( 13*me + 11 + 0 )*stride_out] = (*R22).y;
	bufOutRe[outOffset + ( 13*me + 12 + 0 )*stride_out] = (*R24).x;
	bufOutIm[outOffset + ( 13*me + 12 + 0 )*stride_out] = (*R24).y;
	bufOutRe[outOffset + ( 13*me + 0 + 52 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 13*me + 0 + 52 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 13*me + 1 + 52 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 13*me + 1 + 52 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 13*me + 2 + 52 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 13*me + 2 + 52 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 13*me + 3 + 52 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 13*me + 3 + 52 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 13*me + 4 + 52 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 13*me + 4 + 52 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 13*me + 5 + 52 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 13*me + 5 + 52 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 13*me + 6 + 52 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 13*me + 6 + 52 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 13*me + 7 + 52 )*stride_out] = (*R15).x;
	bufOutIm[outOffset + ( 13*me + 7 + 52 )*stride_out] = (*R15).y;
	bufOutRe[outOffset + ( 13*me + 8 + 52 )*stride_out] = (*R17).x;
	bufOutIm[outOffset + ( 13*me + 8 + 52 )*stride_out] = (*R17).y;
	bufOutRe[outOffset + ( 13*me + 9 + 52 )*stride_out] = (*R19).x;
	bufOutIm[outOffset + ( 13*me + 9 + 52 )*stride_out] = (*R19).y;
	bufOutRe[outOffset + ( 13*me + 10 + 52 )*stride_out] = (*R21).x;
	bufOutIm[outOffset + ( 13*me + 10 + 52 )*stride_out] = (*R21).y;
	bufOutRe[outOffset + ( 13*me + 11 + 52 )*stride_out] = (*R23).x;
	bufOutIm[outOffset + ( 13*me + 11 + 52 )*stride_out] = (*R23).y;
	bufOutRe[outOffset + ( 13*me + 12 + 52 )*stride_out] = (*R25).x;
	bufOutIm[outOffset + ( 13*me + 12 + 52 )*stride_out] = (*R25).y;
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass0_len104(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R13) = bufIn[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R1) = bufIn[inOffset + ( 0 + me*2 + 0 + 8 )*stride_in];
	(*R14) = bufIn[inOffset + ( 0 + me*2 + 1 + 8 )*stride_in];
	(*R2) = bufIn[inOffset + ( 0 + me*2 + 0 + 16 )*stride_in];
	(*R15) = bufIn[inOffset + ( 0 + me*2 + 1 + 16 )*stride_in];
	(*R3) = bufIn[inOffset + ( 0 + me*2 + 0 + 24 )*stride_in];
	(*R16) = bufIn[inOffset + ( 0 + me*2 + 1 + 24 )*stride_in];
	(*R4) = bufIn[inOffset + ( 0 + me*2 + 0 + 32 )*stride_in];
	(*R17) = bufIn[inOffset + ( 0 + me*2 + 1 + 32 )*stride_in];
	(*R5) = bufIn[inOffset + ( 0 + me*2 + 0 + 40 )*stride_in];
	(*R18) = bufIn[inOffset + ( 0 + me*2 + 1 + 40 )*stride_in];
	(*R6) = bufIn[inOffset + ( 0 + me*2 + 0 + 48 )*stride_in];
	(*R19) = bufIn[inOffset + ( 0 + me*2 + 1 + 48 )*stride_in];
	(*R7) = bufIn[inOffset + ( 0 + me*2 + 0 + 56 )*stride_in];
	(*R20) = bufIn[inOffset + ( 0 + me*2 + 1 + 56 )*stride_in];
	(*R8) = bufIn[inOffset + ( 0 + me*2 + 0 + 64 )*stride_in];
	(*R21) = bufIn[inOffset + ( 0 + me*2 + 1 + 64 )*stride_in];
	(*R9) = bufIn[inOffset + ( 0 + me*2 + 0 + 72 )*stride_in];
	(*R22) = bufIn[inOffset + ( 0 + me*2 + 1 + 72 )*stride_in];
	(*R10) = bufIn[inOffset + ( 0 + me*2 + 0 + 80 )*stride_in];
	(*R23) = bufIn[inOffset + ( 0 + me*2 + 1 + 80 )*stride_in];
	(*R11) = bufIn[inOffset + ( 0 + me*2 + 0 + 88 )*stride_in];
	(*R24) = bufIn[inOffset + ( 0 + me*2 + 1 + 88 )*stride_in];
	(*R12) = bufIn[inOffset + ( 0 + me*2 + 0 + 96 )*stride_in];
	(*R25) = bufIn[inOffset + ( 0 + me*2 + 1 + 96 )*stride_in];
	}



	InvRad13B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12);
	InvRad13B1(R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25);


	if(rw)
	{
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 10 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 11 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 12 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 0 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 1 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 2 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 3 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 4 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 5 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 6 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 7 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 8 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 9 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 10 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 11 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 12 ) ] = (*R25).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 10 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 11 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 12 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 0 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 1 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 2 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 3 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 4 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 5 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 6 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 7 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 8 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 9 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 10 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 11 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 12 ) ] = (*R25).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass0_len104(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 0 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 8 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 8 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 8 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 8 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 16 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 16 )*stride_in];
	(*R15).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 16 )*stride_in];
	(*R15).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 16 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 24 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 24 )*stride_in];
	(*R16).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 24 )*stride_in];
	(*R16).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 24 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 32 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 32 )*stride_in];
	(*R17).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 32 )*stride_in];
	(*R17).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 32 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 40 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 40 )*stride_in];
	(*R18).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 40 )*stride_in];
	(*R18).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 40 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 48 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 48 )*stride_in];
	(*R19).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 48 )*stride_in];
	(*R19).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 48 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 56 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 56 )*stride_in];
	(*R20).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 56 )*stride_in];
	(*R20).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 56 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 64 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 64 )*stride_in];
	(*R21).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 64 )*stride_in];
	(*R21).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 64 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 72 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 72 )*stride_in];
	(*R22).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 72 )*stride_in];
	(*R22).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 72 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 80 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 80 )*stride_in];
	(*R23).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 80 )*stride_in];
	(*R23).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 80 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 88 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 88 )*stride_in];
	(*R24).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 88 )*stride_in];
	(*R24).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 88 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*2 + 0 + 96 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*2 + 0 + 96 )*stride_in];
	(*R25).x = bufInRe[inOffset + ( 0 + me*2 + 1 + 96 )*stride_in];
	(*R25).y = bufInIm[inOffset + ( 0 + me*2 + 1 + 96 )*stride_in];
	}



	InvRad13B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12);
	InvRad13B1(R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25);


	if(rw)
	{
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 8 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 9 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 10 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 11 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 12 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 0 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 1 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 2 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 3 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 4 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 5 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 6 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 7 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 8 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 9 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 10 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 11 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 12 ) ] = (*R25).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 8 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 9 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 10 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 11 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((2*me + 0)/1)*13 + (2*me + 0)%1 + 12 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 0 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 1 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 2 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 3 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 4 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 5 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 6 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 7 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 8 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 9 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 10 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 11 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((2*me + 1)/1)*13 + (2*me + 1)%1 + 12 ) ] = (*R25).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass1_len104(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25)
{




	{
		T W = twiddles[12 + 1*((13*me + 0)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 1)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 2)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 3)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 4)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 5)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 6)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 7)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 8)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 9)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 10)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 11)%13) + 0];
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
		T W = twiddles[12 + 1*((13*me + 12)%13) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
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
	InvRad2B1(R20, R21);
	InvRad2B1(R22, R23);
	InvRad2B1(R24, R25);


	if(rw)
	{
	bufOutRe[outOffset + ( ((13*me + 0)/13)*26 + (13*me + 0)%13 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((13*me + 0)/13)*26 + (13*me + 0)%13 + 13 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((13*me + 1)/13)*26 + (13*me + 1)%13 + 0 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((13*me + 1)/13)*26 + (13*me + 1)%13 + 13 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((13*me + 2)/13)*26 + (13*me + 2)%13 + 0 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((13*me + 2)/13)*26 + (13*me + 2)%13 + 13 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((13*me + 3)/13)*26 + (13*me + 3)%13 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((13*me + 3)/13)*26 + (13*me + 3)%13 + 13 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((13*me + 4)/13)*26 + (13*me + 4)%13 + 0 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((13*me + 4)/13)*26 + (13*me + 4)%13 + 13 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((13*me + 5)/13)*26 + (13*me + 5)%13 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((13*me + 5)/13)*26 + (13*me + 5)%13 + 13 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((13*me + 6)/13)*26 + (13*me + 6)%13 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((13*me + 6)/13)*26 + (13*me + 6)%13 + 13 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((13*me + 7)/13)*26 + (13*me + 7)%13 + 0 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((13*me + 7)/13)*26 + (13*me + 7)%13 + 13 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((13*me + 8)/13)*26 + (13*me + 8)%13 + 0 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((13*me + 8)/13)*26 + (13*me + 8)%13 + 13 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((13*me + 9)/13)*26 + (13*me + 9)%13 + 0 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((13*me + 9)/13)*26 + (13*me + 9)%13 + 13 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((13*me + 10)/13)*26 + (13*me + 10)%13 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((13*me + 10)/13)*26 + (13*me + 10)%13 + 13 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((13*me + 11)/13)*26 + (13*me + 11)%13 + 0 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((13*me + 11)/13)*26 + (13*me + 11)%13 + 13 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((13*me + 12)/13)*26 + (13*me + 12)%13 + 0 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((13*me + 12)/13)*26 + (13*me + 12)%13 + 13 ) ] = (*R25).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((13*me + 0)/13)*26 + (13*me + 0)%13 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((13*me + 0)/13)*26 + (13*me + 0)%13 + 13 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((13*me + 1)/13)*26 + (13*me + 1)%13 + 0 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((13*me + 1)/13)*26 + (13*me + 1)%13 + 13 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((13*me + 2)/13)*26 + (13*me + 2)%13 + 0 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((13*me + 2)/13)*26 + (13*me + 2)%13 + 13 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((13*me + 3)/13)*26 + (13*me + 3)%13 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((13*me + 3)/13)*26 + (13*me + 3)%13 + 13 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((13*me + 4)/13)*26 + (13*me + 4)%13 + 0 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((13*me + 4)/13)*26 + (13*me + 4)%13 + 13 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((13*me + 5)/13)*26 + (13*me + 5)%13 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((13*me + 5)/13)*26 + (13*me + 5)%13 + 13 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((13*me + 6)/13)*26 + (13*me + 6)%13 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((13*me + 6)/13)*26 + (13*me + 6)%13 + 13 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((13*me + 7)/13)*26 + (13*me + 7)%13 + 0 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((13*me + 7)/13)*26 + (13*me + 7)%13 + 13 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((13*me + 8)/13)*26 + (13*me + 8)%13 + 0 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((13*me + 8)/13)*26 + (13*me + 8)%13 + 13 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((13*me + 9)/13)*26 + (13*me + 9)%13 + 0 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((13*me + 9)/13)*26 + (13*me + 9)%13 + 13 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((13*me + 10)/13)*26 + (13*me + 10)%13 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((13*me + 10)/13)*26 + (13*me + 10)%13 + 13 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((13*me + 11)/13)*26 + (13*me + 11)%13 + 0 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((13*me + 11)/13)*26 + (13*me + 11)%13 + 13 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((13*me + 12)/13)*26 + (13*me + 12)%13 + 0 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((13*me + 12)/13)*26 + (13*me + 12)%13 + 13 ) ] = (*R25).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass2_len104(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25)
{




	{
		T W = twiddles[25 + 1*((13*me + 0)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 1)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 2)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 3)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 4)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 5)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 6)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 7)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 8)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 9)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 10)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 11)%26) + 0];
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
		T W = twiddles[25 + 1*((13*me + 12)%26) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
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
	InvRad2B1(R20, R21);
	InvRad2B1(R22, R23);
	InvRad2B1(R24, R25);


	if(rw)
	{
	bufOutRe[outOffset + ( ((13*me + 0)/26)*52 + (13*me + 0)%26 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((13*me + 0)/26)*52 + (13*me + 0)%26 + 26 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((13*me + 1)/26)*52 + (13*me + 1)%26 + 0 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((13*me + 1)/26)*52 + (13*me + 1)%26 + 26 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((13*me + 2)/26)*52 + (13*me + 2)%26 + 0 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((13*me + 2)/26)*52 + (13*me + 2)%26 + 26 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((13*me + 3)/26)*52 + (13*me + 3)%26 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((13*me + 3)/26)*52 + (13*me + 3)%26 + 26 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((13*me + 4)/26)*52 + (13*me + 4)%26 + 0 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((13*me + 4)/26)*52 + (13*me + 4)%26 + 26 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((13*me + 5)/26)*52 + (13*me + 5)%26 + 0 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((13*me + 5)/26)*52 + (13*me + 5)%26 + 26 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((13*me + 6)/26)*52 + (13*me + 6)%26 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((13*me + 6)/26)*52 + (13*me + 6)%26 + 26 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((13*me + 7)/26)*52 + (13*me + 7)%26 + 0 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((13*me + 7)/26)*52 + (13*me + 7)%26 + 26 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((13*me + 8)/26)*52 + (13*me + 8)%26 + 0 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((13*me + 8)/26)*52 + (13*me + 8)%26 + 26 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((13*me + 9)/26)*52 + (13*me + 9)%26 + 0 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((13*me + 9)/26)*52 + (13*me + 9)%26 + 26 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((13*me + 10)/26)*52 + (13*me + 10)%26 + 0 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((13*me + 10)/26)*52 + (13*me + 10)%26 + 26 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((13*me + 11)/26)*52 + (13*me + 11)%26 + 0 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((13*me + 11)/26)*52 + (13*me + 11)%26 + 26 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((13*me + 12)/26)*52 + (13*me + 12)%26 + 0 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((13*me + 12)/26)*52 + (13*me + 12)%26 + 26 ) ] = (*R25).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((13*me + 0)/26)*52 + (13*me + 0)%26 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((13*me + 0)/26)*52 + (13*me + 0)%26 + 26 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((13*me + 1)/26)*52 + (13*me + 1)%26 + 0 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((13*me + 1)/26)*52 + (13*me + 1)%26 + 26 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((13*me + 2)/26)*52 + (13*me + 2)%26 + 0 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((13*me + 2)/26)*52 + (13*me + 2)%26 + 26 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((13*me + 3)/26)*52 + (13*me + 3)%26 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((13*me + 3)/26)*52 + (13*me + 3)%26 + 26 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((13*me + 4)/26)*52 + (13*me + 4)%26 + 0 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((13*me + 4)/26)*52 + (13*me + 4)%26 + 26 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((13*me + 5)/26)*52 + (13*me + 5)%26 + 0 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((13*me + 5)/26)*52 + (13*me + 5)%26 + 26 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((13*me + 6)/26)*52 + (13*me + 6)%26 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((13*me + 6)/26)*52 + (13*me + 6)%26 + 26 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((13*me + 7)/26)*52 + (13*me + 7)%26 + 0 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((13*me + 7)/26)*52 + (13*me + 7)%26 + 26 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((13*me + 8)/26)*52 + (13*me + 8)%26 + 0 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((13*me + 8)/26)*52 + (13*me + 8)%26 + 26 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((13*me + 9)/26)*52 + (13*me + 9)%26 + 0 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((13*me + 9)/26)*52 + (13*me + 9)%26 + 26 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((13*me + 10)/26)*52 + (13*me + 10)%26 + 0 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((13*me + 10)/26)*52 + (13*me + 10)%26 + 26 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((13*me + 11)/26)*52 + (13*me + 11)%26 + 0 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((13*me + 11)/26)*52 + (13*me + 11)%26 + 26 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((13*me + 12)/26)*52 + (13*me + 12)%26 + 0 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((13*me + 12)/26)*52 + (13*me + 12)%26 + 26 ) ] = (*R25).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 0 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 0 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 0 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*13 + 0 + 52 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*13 + 1 + 52 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*13 + 2 + 52 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*13 + 3 + 52 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*13 + 4 + 52 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*13 + 5 + 52 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*13 + 6 + 52 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*13 + 7 + 52 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*13 + 8 + 52 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*13 + 9 + 52 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*13 + 10 + 52 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*13 + 11 + 52 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*13 + 12 + 52 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass3_len104(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25)
{




	{
		T W = twiddles[51 + 1*((13*me + 0)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 1)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 2)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 3)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 4)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 5)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 6)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 7)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 8)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 9)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 10)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 11)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 12)%52) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
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
	InvRad2B1(R20, R21);
	InvRad2B1(R22, R23);
	InvRad2B1(R24, R25);


	if(rw)
	{
	bufOut[outOffset + ( 13*me + 0 + 0 )*stride_out] = (*R0);
	bufOut[outOffset + ( 13*me + 1 + 0 )*stride_out] = (*R2);
	bufOut[outOffset + ( 13*me + 2 + 0 )*stride_out] = (*R4);
	bufOut[outOffset + ( 13*me + 3 + 0 )*stride_out] = (*R6);
	bufOut[outOffset + ( 13*me + 4 + 0 )*stride_out] = (*R8);
	bufOut[outOffset + ( 13*me + 5 + 0 )*stride_out] = (*R10);
	bufOut[outOffset + ( 13*me + 6 + 0 )*stride_out] = (*R12);
	bufOut[outOffset + ( 13*me + 7 + 0 )*stride_out] = (*R14);
	bufOut[outOffset + ( 13*me + 8 + 0 )*stride_out] = (*R16);
	bufOut[outOffset + ( 13*me + 9 + 0 )*stride_out] = (*R18);
	bufOut[outOffset + ( 13*me + 10 + 0 )*stride_out] = (*R20);
	bufOut[outOffset + ( 13*me + 11 + 0 )*stride_out] = (*R22);
	bufOut[outOffset + ( 13*me + 12 + 0 )*stride_out] = (*R24);
	bufOut[outOffset + ( 13*me + 0 + 52 )*stride_out] = (*R1);
	bufOut[outOffset + ( 13*me + 1 + 52 )*stride_out] = (*R3);
	bufOut[outOffset + ( 13*me + 2 + 52 )*stride_out] = (*R5);
	bufOut[outOffset + ( 13*me + 3 + 52 )*stride_out] = (*R7);
	bufOut[outOffset + ( 13*me + 4 + 52 )*stride_out] = (*R9);
	bufOut[outOffset + ( 13*me + 5 + 52 )*stride_out] = (*R11);
	bufOut[outOffset + ( 13*me + 6 + 52 )*stride_out] = (*R13);
	bufOut[outOffset + ( 13*me + 7 + 52 )*stride_out] = (*R15);
	bufOut[outOffset + ( 13*me + 8 + 52 )*stride_out] = (*R17);
	bufOut[outOffset + ( 13*me + 9 + 52 )*stride_out] = (*R19);
	bufOut[outOffset + ( 13*me + 10 + 52 )*stride_out] = (*R21);
	bufOut[outOffset + ( 13*me + 11 + 52 )*stride_out] = (*R23);
	bufOut[outOffset + ( 13*me + 12 + 52 )*stride_out] = (*R25);
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass3_len104(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25)
{




	{
		T W = twiddles[51 + 1*((13*me + 0)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 1)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 2)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 3)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 4)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 5)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 6)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 7)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 8)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 9)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 10)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 11)%52) + 0];
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
		T W = twiddles[51 + 1*((13*me + 12)%52) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R25).x; ry = (*R25).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R25).x = TR;
		(*R25).y = TI;
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
	InvRad2B1(R20, R21);
	InvRad2B1(R22, R23);
	InvRad2B1(R24, R25);


	if(rw)
	{
	bufOutRe[outOffset + ( 13*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 13*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 13*me + 1 + 0 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 13*me + 1 + 0 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 13*me + 2 + 0 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 13*me + 2 + 0 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 13*me + 3 + 0 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 13*me + 3 + 0 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 13*me + 4 + 0 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 13*me + 4 + 0 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 13*me + 5 + 0 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 13*me + 5 + 0 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 13*me + 6 + 0 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 13*me + 6 + 0 )*stride_out] = (*R12).y;
	bufOutRe[outOffset + ( 13*me + 7 + 0 )*stride_out] = (*R14).x;
	bufOutIm[outOffset + ( 13*me + 7 + 0 )*stride_out] = (*R14).y;
	bufOutRe[outOffset + ( 13*me + 8 + 0 )*stride_out] = (*R16).x;
	bufOutIm[outOffset + ( 13*me + 8 + 0 )*stride_out] = (*R16).y;
	bufOutRe[outOffset + ( 13*me + 9 + 0 )*stride_out] = (*R18).x;
	bufOutIm[outOffset + ( 13*me + 9 + 0 )*stride_out] = (*R18).y;
	bufOutRe[outOffset + ( 13*me + 10 + 0 )*stride_out] = (*R20).x;
	bufOutIm[outOffset + ( 13*me + 10 + 0 )*stride_out] = (*R20).y;
	bufOutRe[outOffset + ( 13*me + 11 + 0 )*stride_out] = (*R22).x;
	bufOutIm[outOffset + ( 13*me + 11 + 0 )*stride_out] = (*R22).y;
	bufOutRe[outOffset + ( 13*me + 12 + 0 )*stride_out] = (*R24).x;
	bufOutIm[outOffset + ( 13*me + 12 + 0 )*stride_out] = (*R24).y;
	bufOutRe[outOffset + ( 13*me + 0 + 52 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 13*me + 0 + 52 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 13*me + 1 + 52 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 13*me + 1 + 52 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 13*me + 2 + 52 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 13*me + 2 + 52 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 13*me + 3 + 52 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 13*me + 3 + 52 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 13*me + 4 + 52 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 13*me + 4 + 52 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 13*me + 5 + 52 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 13*me + 5 + 52 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 13*me + 6 + 52 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 13*me + 6 + 52 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 13*me + 7 + 52 )*stride_out] = (*R15).x;
	bufOutIm[outOffset + ( 13*me + 7 + 52 )*stride_out] = (*R15).y;
	bufOutRe[outOffset + ( 13*me + 8 + 52 )*stride_out] = (*R17).x;
	bufOutIm[outOffset + ( 13*me + 8 + 52 )*stride_out] = (*R17).y;
	bufOutRe[outOffset + ( 13*me + 9 + 52 )*stride_out] = (*R19).x;
	bufOutIm[outOffset + ( 13*me + 9 + 52 )*stride_out] = (*R19).y;
	bufOutRe[outOffset + ( 13*me + 10 + 52 )*stride_out] = (*R21).x;
	bufOutIm[outOffset + ( 13*me + 10 + 52 )*stride_out] = (*R21).y;
	bufOutRe[outOffset + ( 13*me + 11 + 52 )*stride_out] = (*R23).x;
	bufOutIm[outOffset + ( 13*me + 11 + 52 )*stride_out] = (*R23).y;
	bufOutRe[outOffset + ( 13*me + 12 + 52 )*stride_out] = (*R25).x;
	bufOutIm[outOffset + ( 13*me + 12 + 52 )*stride_out] = (*R25).y;
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb>
__device__ void 
fwd_len104_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25;
	FwdPass0_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	FwdPass1_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	FwdPass2_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	FwdPass3_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
}

template <typename T, StrideBin sb>
__device__ void 
back_len104_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25;
	InvPass0_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	InvPass1_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	InvPass2_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	InvPass3_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len104_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25;
	FwdPass0_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	FwdPass1_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	FwdPass2_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	FwdPass3_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
}

template <typename T, StrideBin sb>
__device__ void 
back_len104_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25;
	InvPass0_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	InvPass1_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	InvPass2_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	InvPass3_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len104_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25;
	FwdPass0_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	FwdPass1_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	FwdPass2_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	FwdPass3_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
}

template <typename T, StrideBin sb>
__device__ void 
back_len104_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25;
	InvPass0_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	InvPass1_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	InvPass2_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	InvPass3_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len104_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25;
	FwdPass0_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	FwdPass1_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	FwdPass2_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	FwdPass3_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
}

template <typename T, StrideBin sb>
__device__ void 
back_len104_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25;
	InvPass0_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	InvPass1_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	InvPass2_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
	InvPass3_len104<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25);
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 128, maximum transforms: 32, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(128)
fft_fwd_ip_len104( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3328), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len104>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*32)*4) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*32 + (me/4));
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
	fwd_len104_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%4, (me/4)*104, lwb, lwb, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 128, maximum transforms: 32, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(128)
fft_back_ip_len104( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3328), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len104>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*32)*4) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*32 + (me/4));
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
	back_len104_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%4, (me/4)*104, lwb, lwb, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 128, maximum transforms: 32, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(128)
fft_fwd_ip_len104( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3328), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len104>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*32)*4) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*32 + (me/4));
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
	fwd_len104_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%4, (me/4)*104, lwbRe, lwbIm, lwbRe, lwbIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 128, maximum transforms: 32, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(128)
fft_back_ip_len104( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3328), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len104>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*32)*4) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*32 + (me/4));
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
	back_len104_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%4, (me/4)*104, lwbRe, lwbIm, lwbRe, lwbIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 128, maximum transforms: 32, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(128)
fft_fwd_op_len104( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3328), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len104>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*32)*4) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*32 + (me/4));
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
	fwd_len104_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*104, lwbIn, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 128, maximum transforms: 32, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(128)
fft_back_op_len104( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3328), cgh);
	cgh.parallel_for<class kern_fft_back_op_len104>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*32)*4) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*32 + (me/4));
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
	back_len104_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*104, lwbIn, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 128, maximum transforms: 32, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(128)
fft_fwd_op_len104( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3328), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len104>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*32)*4) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*32 + (me/4));
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
	fwd_len104_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*104, lwbIn, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 128, maximum transforms: 32, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(128)
fft_back_op_len104( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3328), cgh);
	cgh.parallel_for<class kern_fft_back_op_len104>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*32)*4) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*32 + (me/4));
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
	back_len104_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*104, lwbIn, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 128, maximum transforms: 32, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(128)
fft_fwd_op_len104( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3328), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len104>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*32)*4) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*32 + (me/4));
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
	fwd_len104_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*104, lwbInRe, lwbInIm, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 128, maximum transforms: 32, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(128)
fft_back_op_len104( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3328), cgh);
	cgh.parallel_for<class kern_fft_back_op_len104>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*32)*4) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*32 + (me/4));
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
	back_len104_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*104, lwbInRe, lwbInIm, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 128, maximum transforms: 32, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(128)
fft_fwd_op_len104( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3328), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len104>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*32)*4) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*32 + (me/4));
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
	fwd_len104_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*104, lwbInRe, lwbInIm, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 128, maximum transforms: 32, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(128)
fft_back_op_len104( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3328), cgh);
	cgh.parallel_for<class kern_fft_back_op_len104>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*32)*4) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*32 + (me/4));
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
	back_len104_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*104, lwbInRe, lwbInIm, lwbOutRe, lwbOutIm, lds);
});
});
}

