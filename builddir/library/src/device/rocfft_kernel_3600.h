#pragma once
#include "rocfft_butterfly_template.h"

////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb>
__device__ void
FwdPass0_len3600(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R10) = bufIn[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R20) = bufIn[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R1) = bufIn[inOffset + ( 0 + me*3 + 0 + 360 )*stride_in];
	(*R11) = bufIn[inOffset + ( 0 + me*3 + 1 + 360 )*stride_in];
	(*R21) = bufIn[inOffset + ( 0 + me*3 + 2 + 360 )*stride_in];
	(*R2) = bufIn[inOffset + ( 0 + me*3 + 0 + 720 )*stride_in];
	(*R12) = bufIn[inOffset + ( 0 + me*3 + 1 + 720 )*stride_in];
	(*R22) = bufIn[inOffset + ( 0 + me*3 + 2 + 720 )*stride_in];
	(*R3) = bufIn[inOffset + ( 0 + me*3 + 0 + 1080 )*stride_in];
	(*R13) = bufIn[inOffset + ( 0 + me*3 + 1 + 1080 )*stride_in];
	(*R23) = bufIn[inOffset + ( 0 + me*3 + 2 + 1080 )*stride_in];
	(*R4) = bufIn[inOffset + ( 0 + me*3 + 0 + 1440 )*stride_in];
	(*R14) = bufIn[inOffset + ( 0 + me*3 + 1 + 1440 )*stride_in];
	(*R24) = bufIn[inOffset + ( 0 + me*3 + 2 + 1440 )*stride_in];
	(*R5) = bufIn[inOffset + ( 0 + me*3 + 0 + 1800 )*stride_in];
	(*R15) = bufIn[inOffset + ( 0 + me*3 + 1 + 1800 )*stride_in];
	(*R25) = bufIn[inOffset + ( 0 + me*3 + 2 + 1800 )*stride_in];
	(*R6) = bufIn[inOffset + ( 0 + me*3 + 0 + 2160 )*stride_in];
	(*R16) = bufIn[inOffset + ( 0 + me*3 + 1 + 2160 )*stride_in];
	(*R26) = bufIn[inOffset + ( 0 + me*3 + 2 + 2160 )*stride_in];
	(*R7) = bufIn[inOffset + ( 0 + me*3 + 0 + 2520 )*stride_in];
	(*R17) = bufIn[inOffset + ( 0 + me*3 + 1 + 2520 )*stride_in];
	(*R27) = bufIn[inOffset + ( 0 + me*3 + 2 + 2520 )*stride_in];
	(*R8) = bufIn[inOffset + ( 0 + me*3 + 0 + 2880 )*stride_in];
	(*R18) = bufIn[inOffset + ( 0 + me*3 + 1 + 2880 )*stride_in];
	(*R28) = bufIn[inOffset + ( 0 + me*3 + 2 + 2880 )*stride_in];
	(*R9) = bufIn[inOffset + ( 0 + me*3 + 0 + 3240 )*stride_in];
	(*R19) = bufIn[inOffset + ( 0 + me*3 + 1 + 3240 )*stride_in];
	(*R29) = bufIn[inOffset + ( 0 + me*3 + 2 + 3240 )*stride_in];
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
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 360 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 360 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 360 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 720 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 720 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 720 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1080 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1080 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1080 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1440 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1440 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1440 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2160 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2160 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2160 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2520 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2520 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2520 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2880 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2880 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2880 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 3240 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 3240 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 3240 ) ];
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
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 360 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 360 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 360 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 720 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 720 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 720 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1080 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1080 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1080 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1440 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1440 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1440 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2160 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2160 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2160 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2520 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2520 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2520 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2880 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2880 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2880 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 3240 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 3240 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 3240 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass0_len3600(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R20).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R20).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 360 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 360 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 360 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 360 )*stride_in];
	(*R21).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 360 )*stride_in];
	(*R21).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 360 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 720 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 720 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 720 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 720 )*stride_in];
	(*R22).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 720 )*stride_in];
	(*R22).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 720 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1080 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1080 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1080 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1080 )*stride_in];
	(*R23).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1080 )*stride_in];
	(*R23).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1080 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1440 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1440 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1440 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1440 )*stride_in];
	(*R24).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1440 )*stride_in];
	(*R24).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1440 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1800 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1800 )*stride_in];
	(*R15).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1800 )*stride_in];
	(*R15).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1800 )*stride_in];
	(*R25).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1800 )*stride_in];
	(*R25).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1800 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2160 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2160 )*stride_in];
	(*R16).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2160 )*stride_in];
	(*R16).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2160 )*stride_in];
	(*R26).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2160 )*stride_in];
	(*R26).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2160 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2520 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2520 )*stride_in];
	(*R17).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2520 )*stride_in];
	(*R17).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2520 )*stride_in];
	(*R27).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2520 )*stride_in];
	(*R27).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2520 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2880 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2880 )*stride_in];
	(*R18).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2880 )*stride_in];
	(*R18).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2880 )*stride_in];
	(*R28).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2880 )*stride_in];
	(*R28).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2880 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 3240 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 3240 )*stride_in];
	(*R19).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 3240 )*stride_in];
	(*R19).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 3240 )*stride_in];
	(*R29).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 3240 )*stride_in];
	(*R29).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 3240 )*stride_in];
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
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 360 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 360 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 360 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 720 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 720 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 720 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1080 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1080 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1080 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1440 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1440 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1440 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2160 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2160 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2160 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2520 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2520 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2520 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2880 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2880 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2880 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 3240 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 3240 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 3240 ) ];
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
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 360 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 360 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 360 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 720 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 720 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 720 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1080 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1080 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1080 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1440 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1440 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1440 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2160 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2160 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2160 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2520 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2520 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2520 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2880 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2880 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2880 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 3240 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 3240 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 3240 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass1_len3600(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
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
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 600 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 600 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 600 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 600 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 600 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1200 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1200 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1200 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1200 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1200 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1800 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1800 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1800 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1800 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1800 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2400 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2400 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2400 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2400 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2400 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 3000 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 3000 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 3000 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 3000 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 3000 ) ];
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
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 600 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 600 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 600 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 600 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 600 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1200 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1200 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1200 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1200 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1200 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1800 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1800 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1800 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1800 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1800 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2400 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2400 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2400 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2400 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2400 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 3000 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 3000 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 3000 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 3000 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 3000 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass2_len3600(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
{




	{
		T W = twiddles[99 + 5*((5*me + 0)%100) + 0];
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
		T W = twiddles[99 + 5*((5*me + 0)%100) + 1];
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
		T W = twiddles[99 + 5*((5*me + 0)%100) + 2];
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
		T W = twiddles[99 + 5*((5*me + 0)%100) + 3];
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
		T W = twiddles[99 + 5*((5*me + 0)%100) + 4];
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
		T W = twiddles[99 + 5*((5*me + 1)%100) + 0];
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
		T W = twiddles[99 + 5*((5*me + 1)%100) + 1];
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
		T W = twiddles[99 + 5*((5*me + 1)%100) + 2];
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
		T W = twiddles[99 + 5*((5*me + 1)%100) + 3];
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
		T W = twiddles[99 + 5*((5*me + 1)%100) + 4];
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
		T W = twiddles[99 + 5*((5*me + 2)%100) + 0];
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
		T W = twiddles[99 + 5*((5*me + 2)%100) + 1];
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
		T W = twiddles[99 + 5*((5*me + 2)%100) + 2];
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
		T W = twiddles[99 + 5*((5*me + 2)%100) + 3];
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
		T W = twiddles[99 + 5*((5*me + 2)%100) + 4];
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
		T W = twiddles[99 + 5*((5*me + 3)%100) + 0];
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
		T W = twiddles[99 + 5*((5*me + 3)%100) + 1];
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
		T W = twiddles[99 + 5*((5*me + 3)%100) + 2];
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
		T W = twiddles[99 + 5*((5*me + 3)%100) + 3];
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
		T W = twiddles[99 + 5*((5*me + 3)%100) + 4];
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
		T W = twiddles[99 + 5*((5*me + 4)%100) + 0];
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
		T W = twiddles[99 + 5*((5*me + 4)%100) + 1];
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
		T W = twiddles[99 + 5*((5*me + 4)%100) + 2];
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
		T W = twiddles[99 + 5*((5*me + 4)%100) + 3];
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
		T W = twiddles[99 + 5*((5*me + 4)%100) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	FwdRad6B1(R0, R1, R2, R3, R4, R5);
	FwdRad6B1(R6, R7, R8, R9, R10, R11);
	FwdRad6B1(R12, R13, R14, R15, R16, R17);
	FwdRad6B1(R18, R19, R20, R21, R22, R23);
	FwdRad6B1(R24, R25, R26, R27, R28, R29);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 100 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 200 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 300 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 400 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 500 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 100 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 200 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 300 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 400 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 500 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 100 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 200 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 300 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 400 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 500 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 0 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 100 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 200 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 300 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 400 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 500 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 0 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 100 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 200 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 300 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 400 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 500 ) ] = (*R29).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 600 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 600 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 600 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 600 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 600 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1200 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1200 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1200 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1200 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1200 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1800 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1800 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1800 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1800 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1800 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2400 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2400 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2400 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2400 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2400 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 3000 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 3000 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 3000 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 3000 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 3000 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 100 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 200 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 300 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 400 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 500 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 100 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 200 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 300 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 400 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 500 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 100 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 200 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 300 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 400 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 500 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 0 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 100 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 200 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 300 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 400 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 500 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 0 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 100 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 200 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 300 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 400 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 500 ) ] = (*R29).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 600 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 600 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 600 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 600 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 600 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1200 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1200 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1200 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1200 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1200 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1800 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1800 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1800 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1800 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1800 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2400 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2400 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2400 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2400 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2400 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 3000 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 3000 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 3000 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 3000 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 3000 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass3_len3600(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
{




	{
		T W = twiddles[599 + 5*((5*me + 0)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	FwdRad6B1(R0, R1, R2, R3, R4, R5);
	FwdRad6B1(R6, R7, R8, R9, R10, R11);
	FwdRad6B1(R12, R13, R14, R15, R16, R17);
	FwdRad6B1(R18, R19, R20, R21, R22, R23);
	FwdRad6B1(R24, R25, R26, R27, R28, R29);


	if(rw)
	{
	bufOut[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0);
	bufOut[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R6);
	bufOut[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R12);
	bufOut[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R18);
	bufOut[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R24);
	bufOut[outOffset + ( 5*me + 0 + 600 )*stride_out] = (*R1);
	bufOut[outOffset + ( 5*me + 1 + 600 )*stride_out] = (*R7);
	bufOut[outOffset + ( 5*me + 2 + 600 )*stride_out] = (*R13);
	bufOut[outOffset + ( 5*me + 3 + 600 )*stride_out] = (*R19);
	bufOut[outOffset + ( 5*me + 4 + 600 )*stride_out] = (*R25);
	bufOut[outOffset + ( 5*me + 0 + 1200 )*stride_out] = (*R2);
	bufOut[outOffset + ( 5*me + 1 + 1200 )*stride_out] = (*R8);
	bufOut[outOffset + ( 5*me + 2 + 1200 )*stride_out] = (*R14);
	bufOut[outOffset + ( 5*me + 3 + 1200 )*stride_out] = (*R20);
	bufOut[outOffset + ( 5*me + 4 + 1200 )*stride_out] = (*R26);
	bufOut[outOffset + ( 5*me + 0 + 1800 )*stride_out] = (*R3);
	bufOut[outOffset + ( 5*me + 1 + 1800 )*stride_out] = (*R9);
	bufOut[outOffset + ( 5*me + 2 + 1800 )*stride_out] = (*R15);
	bufOut[outOffset + ( 5*me + 3 + 1800 )*stride_out] = (*R21);
	bufOut[outOffset + ( 5*me + 4 + 1800 )*stride_out] = (*R27);
	bufOut[outOffset + ( 5*me + 0 + 2400 )*stride_out] = (*R4);
	bufOut[outOffset + ( 5*me + 1 + 2400 )*stride_out] = (*R10);
	bufOut[outOffset + ( 5*me + 2 + 2400 )*stride_out] = (*R16);
	bufOut[outOffset + ( 5*me + 3 + 2400 )*stride_out] = (*R22);
	bufOut[outOffset + ( 5*me + 4 + 2400 )*stride_out] = (*R28);
	bufOut[outOffset + ( 5*me + 0 + 3000 )*stride_out] = (*R5);
	bufOut[outOffset + ( 5*me + 1 + 3000 )*stride_out] = (*R11);
	bufOut[outOffset + ( 5*me + 2 + 3000 )*stride_out] = (*R17);
	bufOut[outOffset + ( 5*me + 3 + 3000 )*stride_out] = (*R23);
	bufOut[outOffset + ( 5*me + 4 + 3000 )*stride_out] = (*R29);
	}

}

template <typename T, StrideBin sb>
__device__ void
FwdPass3_len3600(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
{




	{
		T W = twiddles[599 + 5*((5*me + 0)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	FwdRad6B1(R0, R1, R2, R3, R4, R5);
	FwdRad6B1(R6, R7, R8, R9, R10, R11);
	FwdRad6B1(R12, R13, R14, R15, R16, R17);
	FwdRad6B1(R18, R19, R20, R21, R22, R23);
	FwdRad6B1(R24, R25, R26, R27, R28, R29);


	if(rw)
	{
	bufOutRe[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R12).y;
	bufOutRe[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R18).x;
	bufOutIm[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R18).y;
	bufOutRe[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R24).x;
	bufOutIm[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R24).y;
	bufOutRe[outOffset + ( 5*me + 0 + 600 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 5*me + 0 + 600 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 5*me + 1 + 600 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 5*me + 1 + 600 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 5*me + 2 + 600 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 5*me + 2 + 600 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 5*me + 3 + 600 )*stride_out] = (*R19).x;
	bufOutIm[outOffset + ( 5*me + 3 + 600 )*stride_out] = (*R19).y;
	bufOutRe[outOffset + ( 5*me + 4 + 600 )*stride_out] = (*R25).x;
	bufOutIm[outOffset + ( 5*me + 4 + 600 )*stride_out] = (*R25).y;
	bufOutRe[outOffset + ( 5*me + 0 + 1200 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 5*me + 0 + 1200 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 5*me + 1 + 1200 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 5*me + 1 + 1200 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 5*me + 2 + 1200 )*stride_out] = (*R14).x;
	bufOutIm[outOffset + ( 5*me + 2 + 1200 )*stride_out] = (*R14).y;
	bufOutRe[outOffset + ( 5*me + 3 + 1200 )*stride_out] = (*R20).x;
	bufOutIm[outOffset + ( 5*me + 3 + 1200 )*stride_out] = (*R20).y;
	bufOutRe[outOffset + ( 5*me + 4 + 1200 )*stride_out] = (*R26).x;
	bufOutIm[outOffset + ( 5*me + 4 + 1200 )*stride_out] = (*R26).y;
	bufOutRe[outOffset + ( 5*me + 0 + 1800 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 5*me + 0 + 1800 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 5*me + 1 + 1800 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 5*me + 1 + 1800 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 5*me + 2 + 1800 )*stride_out] = (*R15).x;
	bufOutIm[outOffset + ( 5*me + 2 + 1800 )*stride_out] = (*R15).y;
	bufOutRe[outOffset + ( 5*me + 3 + 1800 )*stride_out] = (*R21).x;
	bufOutIm[outOffset + ( 5*me + 3 + 1800 )*stride_out] = (*R21).y;
	bufOutRe[outOffset + ( 5*me + 4 + 1800 )*stride_out] = (*R27).x;
	bufOutIm[outOffset + ( 5*me + 4 + 1800 )*stride_out] = (*R27).y;
	bufOutRe[outOffset + ( 5*me + 0 + 2400 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 5*me + 0 + 2400 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 5*me + 1 + 2400 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 5*me + 1 + 2400 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 5*me + 2 + 2400 )*stride_out] = (*R16).x;
	bufOutIm[outOffset + ( 5*me + 2 + 2400 )*stride_out] = (*R16).y;
	bufOutRe[outOffset + ( 5*me + 3 + 2400 )*stride_out] = (*R22).x;
	bufOutIm[outOffset + ( 5*me + 3 + 2400 )*stride_out] = (*R22).y;
	bufOutRe[outOffset + ( 5*me + 4 + 2400 )*stride_out] = (*R28).x;
	bufOutIm[outOffset + ( 5*me + 4 + 2400 )*stride_out] = (*R28).y;
	bufOutRe[outOffset + ( 5*me + 0 + 3000 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 5*me + 0 + 3000 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 5*me + 1 + 3000 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 5*me + 1 + 3000 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 5*me + 2 + 3000 )*stride_out] = (*R17).x;
	bufOutIm[outOffset + ( 5*me + 2 + 3000 )*stride_out] = (*R17).y;
	bufOutRe[outOffset + ( 5*me + 3 + 3000 )*stride_out] = (*R23).x;
	bufOutIm[outOffset + ( 5*me + 3 + 3000 )*stride_out] = (*R23).y;
	bufOutRe[outOffset + ( 5*me + 4 + 3000 )*stride_out] = (*R29).x;
	bufOutIm[outOffset + ( 5*me + 4 + 3000 )*stride_out] = (*R29).y;
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass0_len3600(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R10) = bufIn[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R20) = bufIn[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R1) = bufIn[inOffset + ( 0 + me*3 + 0 + 360 )*stride_in];
	(*R11) = bufIn[inOffset + ( 0 + me*3 + 1 + 360 )*stride_in];
	(*R21) = bufIn[inOffset + ( 0 + me*3 + 2 + 360 )*stride_in];
	(*R2) = bufIn[inOffset + ( 0 + me*3 + 0 + 720 )*stride_in];
	(*R12) = bufIn[inOffset + ( 0 + me*3 + 1 + 720 )*stride_in];
	(*R22) = bufIn[inOffset + ( 0 + me*3 + 2 + 720 )*stride_in];
	(*R3) = bufIn[inOffset + ( 0 + me*3 + 0 + 1080 )*stride_in];
	(*R13) = bufIn[inOffset + ( 0 + me*3 + 1 + 1080 )*stride_in];
	(*R23) = bufIn[inOffset + ( 0 + me*3 + 2 + 1080 )*stride_in];
	(*R4) = bufIn[inOffset + ( 0 + me*3 + 0 + 1440 )*stride_in];
	(*R14) = bufIn[inOffset + ( 0 + me*3 + 1 + 1440 )*stride_in];
	(*R24) = bufIn[inOffset + ( 0 + me*3 + 2 + 1440 )*stride_in];
	(*R5) = bufIn[inOffset + ( 0 + me*3 + 0 + 1800 )*stride_in];
	(*R15) = bufIn[inOffset + ( 0 + me*3 + 1 + 1800 )*stride_in];
	(*R25) = bufIn[inOffset + ( 0 + me*3 + 2 + 1800 )*stride_in];
	(*R6) = bufIn[inOffset + ( 0 + me*3 + 0 + 2160 )*stride_in];
	(*R16) = bufIn[inOffset + ( 0 + me*3 + 1 + 2160 )*stride_in];
	(*R26) = bufIn[inOffset + ( 0 + me*3 + 2 + 2160 )*stride_in];
	(*R7) = bufIn[inOffset + ( 0 + me*3 + 0 + 2520 )*stride_in];
	(*R17) = bufIn[inOffset + ( 0 + me*3 + 1 + 2520 )*stride_in];
	(*R27) = bufIn[inOffset + ( 0 + me*3 + 2 + 2520 )*stride_in];
	(*R8) = bufIn[inOffset + ( 0 + me*3 + 0 + 2880 )*stride_in];
	(*R18) = bufIn[inOffset + ( 0 + me*3 + 1 + 2880 )*stride_in];
	(*R28) = bufIn[inOffset + ( 0 + me*3 + 2 + 2880 )*stride_in];
	(*R9) = bufIn[inOffset + ( 0 + me*3 + 0 + 3240 )*stride_in];
	(*R19) = bufIn[inOffset + ( 0 + me*3 + 1 + 3240 )*stride_in];
	(*R29) = bufIn[inOffset + ( 0 + me*3 + 2 + 3240 )*stride_in];
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
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 360 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 360 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 360 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 720 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 720 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 720 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1080 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1080 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1080 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1440 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1440 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1440 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2160 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2160 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2160 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2520 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2520 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2520 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2880 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2880 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2880 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 3240 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 3240 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 3240 ) ];
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
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 360 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 360 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 360 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 720 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 720 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 720 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1080 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1080 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1080 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1440 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1440 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1440 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2160 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2160 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2160 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2520 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2520 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2520 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2880 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2880 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2880 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 3240 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 3240 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 3240 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass0_len3600(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 0 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 0 )*stride_in];
	(*R20).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R20).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 360 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 360 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 360 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 360 )*stride_in];
	(*R21).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 360 )*stride_in];
	(*R21).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 360 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 720 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 720 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 720 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 720 )*stride_in];
	(*R22).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 720 )*stride_in];
	(*R22).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 720 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1080 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1080 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1080 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1080 )*stride_in];
	(*R23).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1080 )*stride_in];
	(*R23).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1080 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1440 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1440 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1440 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1440 )*stride_in];
	(*R24).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1440 )*stride_in];
	(*R24).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1440 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 1800 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 1800 )*stride_in];
	(*R15).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 1800 )*stride_in];
	(*R15).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 1800 )*stride_in];
	(*R25).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 1800 )*stride_in];
	(*R25).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 1800 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2160 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2160 )*stride_in];
	(*R16).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2160 )*stride_in];
	(*R16).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2160 )*stride_in];
	(*R26).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2160 )*stride_in];
	(*R26).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2160 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2520 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2520 )*stride_in];
	(*R17).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2520 )*stride_in];
	(*R17).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2520 )*stride_in];
	(*R27).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2520 )*stride_in];
	(*R27).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2520 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 2880 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 2880 )*stride_in];
	(*R18).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 2880 )*stride_in];
	(*R18).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 2880 )*stride_in];
	(*R28).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 2880 )*stride_in];
	(*R28).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 2880 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*3 + 0 + 3240 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*3 + 0 + 3240 )*stride_in];
	(*R19).x = bufInRe[inOffset + ( 0 + me*3 + 1 + 3240 )*stride_in];
	(*R19).y = bufInIm[inOffset + ( 0 + me*3 + 1 + 3240 )*stride_in];
	(*R29).x = bufInRe[inOffset + ( 0 + me*3 + 2 + 3240 )*stride_in];
	(*R29).y = bufInIm[inOffset + ( 0 + me*3 + 2 + 3240 )*stride_in];
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
	(*R1).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 360 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 360 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 360 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 720 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 720 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 720 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1080 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1080 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1080 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1440 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1440 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1440 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2160 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2160 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2160 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2520 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2520 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2520 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 2880 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 2880 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 2880 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*3 + 0 + 3240 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*3 + 1 + 3240 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*3 + 2 + 3240 ) ];
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
	(*R1).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 360 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 360 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 360 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 720 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 720 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 720 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1080 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1080 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1080 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1440 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1440 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1440 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 1800 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 1800 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 1800 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2160 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2160 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2160 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2520 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2520 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2520 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 2880 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 2880 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 2880 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*3 + 0 + 3240 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*3 + 1 + 3240 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*3 + 2 + 3240 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass1_len3600(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
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
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 600 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 600 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 600 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 600 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 600 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1200 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1200 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1200 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1200 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1200 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1800 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1800 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1800 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1800 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1800 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2400 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2400 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2400 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2400 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2400 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 3000 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 3000 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 3000 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 3000 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 3000 ) ];
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
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 600 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 600 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 600 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 600 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 600 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1200 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1200 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1200 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1200 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1200 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1800 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1800 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1800 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1800 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1800 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2400 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2400 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2400 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2400 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2400 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 3000 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 3000 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 3000 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 3000 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 3000 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass2_len3600(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
{




	{
		T W = twiddles[99 + 5*((5*me + 0)%100) + 0];
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
		T W = twiddles[99 + 5*((5*me + 0)%100) + 1];
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
		T W = twiddles[99 + 5*((5*me + 0)%100) + 2];
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
		T W = twiddles[99 + 5*((5*me + 0)%100) + 3];
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
		T W = twiddles[99 + 5*((5*me + 0)%100) + 4];
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
		T W = twiddles[99 + 5*((5*me + 1)%100) + 0];
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
		T W = twiddles[99 + 5*((5*me + 1)%100) + 1];
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
		T W = twiddles[99 + 5*((5*me + 1)%100) + 2];
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
		T W = twiddles[99 + 5*((5*me + 1)%100) + 3];
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
		T W = twiddles[99 + 5*((5*me + 1)%100) + 4];
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
		T W = twiddles[99 + 5*((5*me + 2)%100) + 0];
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
		T W = twiddles[99 + 5*((5*me + 2)%100) + 1];
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
		T W = twiddles[99 + 5*((5*me + 2)%100) + 2];
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
		T W = twiddles[99 + 5*((5*me + 2)%100) + 3];
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
		T W = twiddles[99 + 5*((5*me + 2)%100) + 4];
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
		T W = twiddles[99 + 5*((5*me + 3)%100) + 0];
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
		T W = twiddles[99 + 5*((5*me + 3)%100) + 1];
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
		T W = twiddles[99 + 5*((5*me + 3)%100) + 2];
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
		T W = twiddles[99 + 5*((5*me + 3)%100) + 3];
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
		T W = twiddles[99 + 5*((5*me + 3)%100) + 4];
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
		T W = twiddles[99 + 5*((5*me + 4)%100) + 0];
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
		T W = twiddles[99 + 5*((5*me + 4)%100) + 1];
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
		T W = twiddles[99 + 5*((5*me + 4)%100) + 2];
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
		T W = twiddles[99 + 5*((5*me + 4)%100) + 3];
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
		T W = twiddles[99 + 5*((5*me + 4)%100) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	InvRad6B1(R0, R1, R2, R3, R4, R5);
	InvRad6B1(R6, R7, R8, R9, R10, R11);
	InvRad6B1(R12, R13, R14, R15, R16, R17);
	InvRad6B1(R18, R19, R20, R21, R22, R23);
	InvRad6B1(R24, R25, R26, R27, R28, R29);


	if(rw)
	{
	bufOutRe[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 100 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 200 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 300 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 400 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 500 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 100 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 200 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 300 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 400 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 500 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 100 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 200 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 300 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 400 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 500 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 0 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 100 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 200 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 300 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 400 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 500 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 0 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 100 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 200 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 300 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 400 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 500 ) ] = (*R29).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 600 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 600 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 600 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 600 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 600 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1200 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1200 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1200 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1200 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1200 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 1800 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 1800 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 1800 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 1800 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 1800 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 2400 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 2400 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 2400 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 2400 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 2400 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*5 + 0 + 3000 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*5 + 1 + 3000 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*5 + 2 + 3000 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*5 + 3 + 3000 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*5 + 4 + 3000 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 100 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 200 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 300 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 400 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((5*me + 0)/100)*600 + (5*me + 0)%100 + 500 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 100 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 200 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 300 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 400 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((5*me + 1)/100)*600 + (5*me + 1)%100 + 500 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 100 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 200 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 300 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 400 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((5*me + 2)/100)*600 + (5*me + 2)%100 + 500 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 0 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 100 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 200 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 300 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 400 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((5*me + 3)/100)*600 + (5*me + 3)%100 + 500 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 0 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 100 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 200 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 300 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 400 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((5*me + 4)/100)*600 + (5*me + 4)%100 + 500 ) ] = (*R29).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 600 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 600 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 600 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 600 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 600 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1200 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1200 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1200 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1200 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1200 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 1800 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 1800 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 1800 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 1800 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 1800 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 2400 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 2400 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 2400 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 2400 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 2400 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*5 + 0 + 3000 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*5 + 1 + 3000 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*5 + 2 + 3000 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*5 + 3 + 3000 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*5 + 4 + 3000 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass3_len3600(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
{




	{
		T W = twiddles[599 + 5*((5*me + 0)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	InvRad6B1(R0, R1, R2, R3, R4, R5);
	InvRad6B1(R6, R7, R8, R9, R10, R11);
	InvRad6B1(R12, R13, R14, R15, R16, R17);
	InvRad6B1(R18, R19, R20, R21, R22, R23);
	InvRad6B1(R24, R25, R26, R27, R28, R29);


	if(rw)
	{
	bufOut[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0);
	bufOut[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R6);
	bufOut[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R12);
	bufOut[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R18);
	bufOut[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R24);
	bufOut[outOffset + ( 5*me + 0 + 600 )*stride_out] = (*R1);
	bufOut[outOffset + ( 5*me + 1 + 600 )*stride_out] = (*R7);
	bufOut[outOffset + ( 5*me + 2 + 600 )*stride_out] = (*R13);
	bufOut[outOffset + ( 5*me + 3 + 600 )*stride_out] = (*R19);
	bufOut[outOffset + ( 5*me + 4 + 600 )*stride_out] = (*R25);
	bufOut[outOffset + ( 5*me + 0 + 1200 )*stride_out] = (*R2);
	bufOut[outOffset + ( 5*me + 1 + 1200 )*stride_out] = (*R8);
	bufOut[outOffset + ( 5*me + 2 + 1200 )*stride_out] = (*R14);
	bufOut[outOffset + ( 5*me + 3 + 1200 )*stride_out] = (*R20);
	bufOut[outOffset + ( 5*me + 4 + 1200 )*stride_out] = (*R26);
	bufOut[outOffset + ( 5*me + 0 + 1800 )*stride_out] = (*R3);
	bufOut[outOffset + ( 5*me + 1 + 1800 )*stride_out] = (*R9);
	bufOut[outOffset + ( 5*me + 2 + 1800 )*stride_out] = (*R15);
	bufOut[outOffset + ( 5*me + 3 + 1800 )*stride_out] = (*R21);
	bufOut[outOffset + ( 5*me + 4 + 1800 )*stride_out] = (*R27);
	bufOut[outOffset + ( 5*me + 0 + 2400 )*stride_out] = (*R4);
	bufOut[outOffset + ( 5*me + 1 + 2400 )*stride_out] = (*R10);
	bufOut[outOffset + ( 5*me + 2 + 2400 )*stride_out] = (*R16);
	bufOut[outOffset + ( 5*me + 3 + 2400 )*stride_out] = (*R22);
	bufOut[outOffset + ( 5*me + 4 + 2400 )*stride_out] = (*R28);
	bufOut[outOffset + ( 5*me + 0 + 3000 )*stride_out] = (*R5);
	bufOut[outOffset + ( 5*me + 1 + 3000 )*stride_out] = (*R11);
	bufOut[outOffset + ( 5*me + 2 + 3000 )*stride_out] = (*R17);
	bufOut[outOffset + ( 5*me + 3 + 3000 )*stride_out] = (*R23);
	bufOut[outOffset + ( 5*me + 4 + 3000 )*stride_out] = (*R29);
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass3_len3600(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29)
{




	{
		T W = twiddles[599 + 5*((5*me + 0)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 0)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 1)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 2)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 3)%600) + 4];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 0];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 1];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 2];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 3];
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
		T W = twiddles[599 + 5*((5*me + 4)%600) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	InvRad6B1(R0, R1, R2, R3, R4, R5);
	InvRad6B1(R6, R7, R8, R9, R10, R11);
	InvRad6B1(R12, R13, R14, R15, R16, R17);
	InvRad6B1(R18, R19, R20, R21, R22, R23);
	InvRad6B1(R24, R25, R26, R27, R28, R29);


	if(rw)
	{
	bufOutRe[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 5*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 5*me + 1 + 0 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 5*me + 2 + 0 )*stride_out] = (*R12).y;
	bufOutRe[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R18).x;
	bufOutIm[outOffset + ( 5*me + 3 + 0 )*stride_out] = (*R18).y;
	bufOutRe[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R24).x;
	bufOutIm[outOffset + ( 5*me + 4 + 0 )*stride_out] = (*R24).y;
	bufOutRe[outOffset + ( 5*me + 0 + 600 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 5*me + 0 + 600 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 5*me + 1 + 600 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 5*me + 1 + 600 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 5*me + 2 + 600 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 5*me + 2 + 600 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 5*me + 3 + 600 )*stride_out] = (*R19).x;
	bufOutIm[outOffset + ( 5*me + 3 + 600 )*stride_out] = (*R19).y;
	bufOutRe[outOffset + ( 5*me + 4 + 600 )*stride_out] = (*R25).x;
	bufOutIm[outOffset + ( 5*me + 4 + 600 )*stride_out] = (*R25).y;
	bufOutRe[outOffset + ( 5*me + 0 + 1200 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 5*me + 0 + 1200 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 5*me + 1 + 1200 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 5*me + 1 + 1200 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 5*me + 2 + 1200 )*stride_out] = (*R14).x;
	bufOutIm[outOffset + ( 5*me + 2 + 1200 )*stride_out] = (*R14).y;
	bufOutRe[outOffset + ( 5*me + 3 + 1200 )*stride_out] = (*R20).x;
	bufOutIm[outOffset + ( 5*me + 3 + 1200 )*stride_out] = (*R20).y;
	bufOutRe[outOffset + ( 5*me + 4 + 1200 )*stride_out] = (*R26).x;
	bufOutIm[outOffset + ( 5*me + 4 + 1200 )*stride_out] = (*R26).y;
	bufOutRe[outOffset + ( 5*me + 0 + 1800 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 5*me + 0 + 1800 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 5*me + 1 + 1800 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 5*me + 1 + 1800 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 5*me + 2 + 1800 )*stride_out] = (*R15).x;
	bufOutIm[outOffset + ( 5*me + 2 + 1800 )*stride_out] = (*R15).y;
	bufOutRe[outOffset + ( 5*me + 3 + 1800 )*stride_out] = (*R21).x;
	bufOutIm[outOffset + ( 5*me + 3 + 1800 )*stride_out] = (*R21).y;
	bufOutRe[outOffset + ( 5*me + 4 + 1800 )*stride_out] = (*R27).x;
	bufOutIm[outOffset + ( 5*me + 4 + 1800 )*stride_out] = (*R27).y;
	bufOutRe[outOffset + ( 5*me + 0 + 2400 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 5*me + 0 + 2400 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 5*me + 1 + 2400 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 5*me + 1 + 2400 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 5*me + 2 + 2400 )*stride_out] = (*R16).x;
	bufOutIm[outOffset + ( 5*me + 2 + 2400 )*stride_out] = (*R16).y;
	bufOutRe[outOffset + ( 5*me + 3 + 2400 )*stride_out] = (*R22).x;
	bufOutIm[outOffset + ( 5*me + 3 + 2400 )*stride_out] = (*R22).y;
	bufOutRe[outOffset + ( 5*me + 4 + 2400 )*stride_out] = (*R28).x;
	bufOutIm[outOffset + ( 5*me + 4 + 2400 )*stride_out] = (*R28).y;
	bufOutRe[outOffset + ( 5*me + 0 + 3000 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 5*me + 0 + 3000 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 5*me + 1 + 3000 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 5*me + 1 + 3000 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 5*me + 2 + 3000 )*stride_out] = (*R17).x;
	bufOutIm[outOffset + ( 5*me + 2 + 3000 )*stride_out] = (*R17).y;
	bufOutRe[outOffset + ( 5*me + 3 + 3000 )*stride_out] = (*R23).x;
	bufOutIm[outOffset + ( 5*me + 3 + 3000 )*stride_out] = (*R23).y;
	bufOutRe[outOffset + ( 5*me + 4 + 3000 )*stride_out] = (*R29).x;
	bufOutIm[outOffset + ( 5*me + 4 + 3000 )*stride_out] = (*R29).y;
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb>
__device__ void 
fwd_len3600_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	FwdPass0_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass1_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass2_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass3_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
}

template <typename T, StrideBin sb>
__device__ void 
back_len3600_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	InvPass0_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass1_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass2_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass3_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len3600_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	FwdPass0_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass1_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass2_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass3_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
}

template <typename T, StrideBin sb>
__device__ void 
back_len3600_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	InvPass0_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass1_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass2_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass3_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len3600_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	FwdPass0_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass1_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass2_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass3_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
}

template <typename T, StrideBin sb>
__device__ void 
back_len3600_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	InvPass0_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass1_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass2_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass3_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len3600_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	FwdPass0_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass1_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass2_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	FwdPass3_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
}

template <typename T, StrideBin sb>
__device__ void 
back_len3600_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29;
	InvPass0_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass1_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass2_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
	InvPass3_len3600<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29);
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 120, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(120)
fft_fwd_ip_len3600( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3600), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len3600>(
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
	fwd_len3600_device<T, sb>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, lwb, lwb, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 120, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(120)
fft_back_ip_len3600( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3600), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len3600>(
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
	back_len3600_device<T, sb>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, lwb, lwb, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 120, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(120)
fft_fwd_ip_len3600( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3600), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len3600>(
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
	fwd_len3600_device<T, sb>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, lwbRe, lwbIm, lwbRe, lwbIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 120, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(120)
fft_back_ip_len3600( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3600), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len3600>(
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
	back_len3600_device<T, sb>(twiddles, stride_in[0], stride_in[0],  1, b, me, 0, lwbRe, lwbIm, lwbRe, lwbIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 120, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(120)
fft_fwd_op_len3600( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3600), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len3600>(
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
	fwd_len3600_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbIn, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 120, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(120)
fft_back_op_len3600( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3600), cgh);
	cgh.parallel_for<class kern_fft_back_op_len3600>(
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
	back_len3600_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbIn, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 120, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(120)
fft_fwd_op_len3600( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3600), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len3600>(
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
	fwd_len3600_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbIn, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 120, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(120)
fft_back_op_len3600( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3600), cgh);
	cgh.parallel_for<class kern_fft_back_op_len3600>(
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
	back_len3600_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbIn, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 120, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(120)
fft_fwd_op_len3600( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3600), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len3600>(
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
	fwd_len3600_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbInRe, lwbInIm, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 120, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(120)
fft_back_op_len3600( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3600), cgh);
	cgh.parallel_for<class kern_fft_back_op_len3600>(
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
	back_len3600_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbInRe, lwbInIm, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 120, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(120)
fft_fwd_op_len3600( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3600), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len3600>(
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
	fwd_len3600_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbInRe, lwbInIm, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 120, maximum transforms: 1, Passes: 4
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(120)
fft_back_op_len3600( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(3600), cgh);
	cgh.parallel_for<class kern_fft_back_op_len3600>(
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
	back_len3600_device<T, sb>(twiddles, stride_in[0], stride_out[0],  1, b, me, 0, lwbInRe, lwbInIm, lwbOutRe, lwbOutIm, lds);
});
});
}

