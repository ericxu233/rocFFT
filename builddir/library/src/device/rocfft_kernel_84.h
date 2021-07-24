#pragma once
#include "rocfft_butterfly_template.h"

////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb>
__device__ void
FwdPass0_len84(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, T *R30, T *R31, T *R32, T *R33, T *R34, T *R35, T *R36, T *R37, T *R38, T *R39, T *R40, T *R41)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*6 + 0 + 0 )*stride_in];
	(*R7) = bufIn[inOffset + ( 0 + me*6 + 1 + 0 )*stride_in];
	(*R14) = bufIn[inOffset + ( 0 + me*6 + 2 + 0 )*stride_in];
	(*R21) = bufIn[inOffset + ( 0 + me*6 + 3 + 0 )*stride_in];
	(*R28) = bufIn[inOffset + ( 0 + me*6 + 4 + 0 )*stride_in];
	(*R35) = bufIn[inOffset + ( 0 + me*6 + 5 + 0 )*stride_in];
	(*R1) = bufIn[inOffset + ( 0 + me*6 + 0 + 12 )*stride_in];
	(*R8) = bufIn[inOffset + ( 0 + me*6 + 1 + 12 )*stride_in];
	(*R15) = bufIn[inOffset + ( 0 + me*6 + 2 + 12 )*stride_in];
	(*R22) = bufIn[inOffset + ( 0 + me*6 + 3 + 12 )*stride_in];
	(*R29) = bufIn[inOffset + ( 0 + me*6 + 4 + 12 )*stride_in];
	(*R36) = bufIn[inOffset + ( 0 + me*6 + 5 + 12 )*stride_in];
	(*R2) = bufIn[inOffset + ( 0 + me*6 + 0 + 24 )*stride_in];
	(*R9) = bufIn[inOffset + ( 0 + me*6 + 1 + 24 )*stride_in];
	(*R16) = bufIn[inOffset + ( 0 + me*6 + 2 + 24 )*stride_in];
	(*R23) = bufIn[inOffset + ( 0 + me*6 + 3 + 24 )*stride_in];
	(*R30) = bufIn[inOffset + ( 0 + me*6 + 4 + 24 )*stride_in];
	(*R37) = bufIn[inOffset + ( 0 + me*6 + 5 + 24 )*stride_in];
	(*R3) = bufIn[inOffset + ( 0 + me*6 + 0 + 36 )*stride_in];
	(*R10) = bufIn[inOffset + ( 0 + me*6 + 1 + 36 )*stride_in];
	(*R17) = bufIn[inOffset + ( 0 + me*6 + 2 + 36 )*stride_in];
	(*R24) = bufIn[inOffset + ( 0 + me*6 + 3 + 36 )*stride_in];
	(*R31) = bufIn[inOffset + ( 0 + me*6 + 4 + 36 )*stride_in];
	(*R38) = bufIn[inOffset + ( 0 + me*6 + 5 + 36 )*stride_in];
	(*R4) = bufIn[inOffset + ( 0 + me*6 + 0 + 48 )*stride_in];
	(*R11) = bufIn[inOffset + ( 0 + me*6 + 1 + 48 )*stride_in];
	(*R18) = bufIn[inOffset + ( 0 + me*6 + 2 + 48 )*stride_in];
	(*R25) = bufIn[inOffset + ( 0 + me*6 + 3 + 48 )*stride_in];
	(*R32) = bufIn[inOffset + ( 0 + me*6 + 4 + 48 )*stride_in];
	(*R39) = bufIn[inOffset + ( 0 + me*6 + 5 + 48 )*stride_in];
	(*R5) = bufIn[inOffset + ( 0 + me*6 + 0 + 60 )*stride_in];
	(*R12) = bufIn[inOffset + ( 0 + me*6 + 1 + 60 )*stride_in];
	(*R19) = bufIn[inOffset + ( 0 + me*6 + 2 + 60 )*stride_in];
	(*R26) = bufIn[inOffset + ( 0 + me*6 + 3 + 60 )*stride_in];
	(*R33) = bufIn[inOffset + ( 0 + me*6 + 4 + 60 )*stride_in];
	(*R40) = bufIn[inOffset + ( 0 + me*6 + 5 + 60 )*stride_in];
	(*R6) = bufIn[inOffset + ( 0 + me*6 + 0 + 72 )*stride_in];
	(*R13) = bufIn[inOffset + ( 0 + me*6 + 1 + 72 )*stride_in];
	(*R20) = bufIn[inOffset + ( 0 + me*6 + 2 + 72 )*stride_in];
	(*R27) = bufIn[inOffset + ( 0 + me*6 + 3 + 72 )*stride_in];
	(*R34) = bufIn[inOffset + ( 0 + me*6 + 4 + 72 )*stride_in];
	(*R41) = bufIn[inOffset + ( 0 + me*6 + 5 + 72 )*stride_in];
	}



	FwdRad7B1(R0, R1, R2, R3, R4, R5, R6);
	FwdRad7B1(R7, R8, R9, R10, R11, R12, R13);
	FwdRad7B1(R14, R15, R16, R17, R18, R19, R20);
	FwdRad7B1(R21, R22, R23, R24, R25, R26, R27);
	FwdRad7B1(R28, R29, R30, R31, R32, R33, R34);
	FwdRad7B1(R35, R36, R37, R38, R39, R40, R41);


	if(rw)
	{
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 0 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 1 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 2 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 3 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 4 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 5 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 6 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 0 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 1 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 2 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 3 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 4 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 5 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 6 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 0 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 1 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 2 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 3 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 4 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 5 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 6 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 0 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 1 ) ] = (*R29).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 2 ) ] = (*R30).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 3 ) ] = (*R31).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 4 ) ] = (*R32).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 5 ) ] = (*R33).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 6 ) ] = (*R34).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 0 ) ] = (*R35).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 1 ) ] = (*R36).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 2 ) ] = (*R37).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 3 ) ] = (*R38).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 4 ) ] = (*R39).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 5 ) ] = (*R40).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 6 ) ] = (*R41).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 0 ) ];
	(*R30).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 0 ) ];
	(*R36).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 14 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 14 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 14 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 14 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 14 ) ];
	(*R31).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 14 ) ];
	(*R37).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 14 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 28 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 28 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 28 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 28 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 28 ) ];
	(*R32).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 28 ) ];
	(*R38).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 28 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 42 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 42 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 42 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 42 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 42 ) ];
	(*R33).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 42 ) ];
	(*R39).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 42 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 56 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 56 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 56 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 56 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 56 ) ];
	(*R34).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 56 ) ];
	(*R40).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 56 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 70 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 70 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 70 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 70 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 70 ) ];
	(*R35).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 70 ) ];
	(*R41).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 70 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 0 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 1 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 2 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 3 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 4 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 5 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 6 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 0 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 1 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 2 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 3 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 4 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 5 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 6 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 0 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 1 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 2 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 3 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 4 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 5 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 6 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 0 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 1 ) ] = (*R29).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 2 ) ] = (*R30).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 3 ) ] = (*R31).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 4 ) ] = (*R32).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 5 ) ] = (*R33).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 6 ) ] = (*R34).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 0 ) ] = (*R35).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 1 ) ] = (*R36).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 2 ) ] = (*R37).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 3 ) ] = (*R38).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 4 ) ] = (*R39).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 5 ) ] = (*R40).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 6 ) ] = (*R41).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 0 ) ];
	(*R30).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 0 ) ];
	(*R36).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 14 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 14 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 14 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 14 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 14 ) ];
	(*R31).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 14 ) ];
	(*R37).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 14 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 28 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 28 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 28 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 28 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 28 ) ];
	(*R32).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 28 ) ];
	(*R38).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 28 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 42 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 42 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 42 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 42 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 42 ) ];
	(*R33).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 42 ) ];
	(*R39).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 42 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 56 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 56 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 56 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 56 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 56 ) ];
	(*R34).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 56 ) ];
	(*R40).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 56 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 70 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 70 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 70 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 70 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 70 ) ];
	(*R35).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 70 ) ];
	(*R41).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 70 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass0_len84(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, T *R30, T *R31, T *R32, T *R33, T *R34, T *R35, T *R36, T *R37, T *R38, T *R39, T *R40, T *R41)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*6 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*6 + 0 + 0 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*6 + 1 + 0 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*6 + 1 + 0 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*6 + 2 + 0 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*6 + 2 + 0 )*stride_in];
	(*R21).x = bufInRe[inOffset + ( 0 + me*6 + 3 + 0 )*stride_in];
	(*R21).y = bufInIm[inOffset + ( 0 + me*6 + 3 + 0 )*stride_in];
	(*R28).x = bufInRe[inOffset + ( 0 + me*6 + 4 + 0 )*stride_in];
	(*R28).y = bufInIm[inOffset + ( 0 + me*6 + 4 + 0 )*stride_in];
	(*R35).x = bufInRe[inOffset + ( 0 + me*6 + 5 + 0 )*stride_in];
	(*R35).y = bufInIm[inOffset + ( 0 + me*6 + 5 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*6 + 0 + 12 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*6 + 0 + 12 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*6 + 1 + 12 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*6 + 1 + 12 )*stride_in];
	(*R15).x = bufInRe[inOffset + ( 0 + me*6 + 2 + 12 )*stride_in];
	(*R15).y = bufInIm[inOffset + ( 0 + me*6 + 2 + 12 )*stride_in];
	(*R22).x = bufInRe[inOffset + ( 0 + me*6 + 3 + 12 )*stride_in];
	(*R22).y = bufInIm[inOffset + ( 0 + me*6 + 3 + 12 )*stride_in];
	(*R29).x = bufInRe[inOffset + ( 0 + me*6 + 4 + 12 )*stride_in];
	(*R29).y = bufInIm[inOffset + ( 0 + me*6 + 4 + 12 )*stride_in];
	(*R36).x = bufInRe[inOffset + ( 0 + me*6 + 5 + 12 )*stride_in];
	(*R36).y = bufInIm[inOffset + ( 0 + me*6 + 5 + 12 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*6 + 0 + 24 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*6 + 0 + 24 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*6 + 1 + 24 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*6 + 1 + 24 )*stride_in];
	(*R16).x = bufInRe[inOffset + ( 0 + me*6 + 2 + 24 )*stride_in];
	(*R16).y = bufInIm[inOffset + ( 0 + me*6 + 2 + 24 )*stride_in];
	(*R23).x = bufInRe[inOffset + ( 0 + me*6 + 3 + 24 )*stride_in];
	(*R23).y = bufInIm[inOffset + ( 0 + me*6 + 3 + 24 )*stride_in];
	(*R30).x = bufInRe[inOffset + ( 0 + me*6 + 4 + 24 )*stride_in];
	(*R30).y = bufInIm[inOffset + ( 0 + me*6 + 4 + 24 )*stride_in];
	(*R37).x = bufInRe[inOffset + ( 0 + me*6 + 5 + 24 )*stride_in];
	(*R37).y = bufInIm[inOffset + ( 0 + me*6 + 5 + 24 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*6 + 0 + 36 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*6 + 0 + 36 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*6 + 1 + 36 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*6 + 1 + 36 )*stride_in];
	(*R17).x = bufInRe[inOffset + ( 0 + me*6 + 2 + 36 )*stride_in];
	(*R17).y = bufInIm[inOffset + ( 0 + me*6 + 2 + 36 )*stride_in];
	(*R24).x = bufInRe[inOffset + ( 0 + me*6 + 3 + 36 )*stride_in];
	(*R24).y = bufInIm[inOffset + ( 0 + me*6 + 3 + 36 )*stride_in];
	(*R31).x = bufInRe[inOffset + ( 0 + me*6 + 4 + 36 )*stride_in];
	(*R31).y = bufInIm[inOffset + ( 0 + me*6 + 4 + 36 )*stride_in];
	(*R38).x = bufInRe[inOffset + ( 0 + me*6 + 5 + 36 )*stride_in];
	(*R38).y = bufInIm[inOffset + ( 0 + me*6 + 5 + 36 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*6 + 0 + 48 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*6 + 0 + 48 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*6 + 1 + 48 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*6 + 1 + 48 )*stride_in];
	(*R18).x = bufInRe[inOffset + ( 0 + me*6 + 2 + 48 )*stride_in];
	(*R18).y = bufInIm[inOffset + ( 0 + me*6 + 2 + 48 )*stride_in];
	(*R25).x = bufInRe[inOffset + ( 0 + me*6 + 3 + 48 )*stride_in];
	(*R25).y = bufInIm[inOffset + ( 0 + me*6 + 3 + 48 )*stride_in];
	(*R32).x = bufInRe[inOffset + ( 0 + me*6 + 4 + 48 )*stride_in];
	(*R32).y = bufInIm[inOffset + ( 0 + me*6 + 4 + 48 )*stride_in];
	(*R39).x = bufInRe[inOffset + ( 0 + me*6 + 5 + 48 )*stride_in];
	(*R39).y = bufInIm[inOffset + ( 0 + me*6 + 5 + 48 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*6 + 0 + 60 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*6 + 0 + 60 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*6 + 1 + 60 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*6 + 1 + 60 )*stride_in];
	(*R19).x = bufInRe[inOffset + ( 0 + me*6 + 2 + 60 )*stride_in];
	(*R19).y = bufInIm[inOffset + ( 0 + me*6 + 2 + 60 )*stride_in];
	(*R26).x = bufInRe[inOffset + ( 0 + me*6 + 3 + 60 )*stride_in];
	(*R26).y = bufInIm[inOffset + ( 0 + me*6 + 3 + 60 )*stride_in];
	(*R33).x = bufInRe[inOffset + ( 0 + me*6 + 4 + 60 )*stride_in];
	(*R33).y = bufInIm[inOffset + ( 0 + me*6 + 4 + 60 )*stride_in];
	(*R40).x = bufInRe[inOffset + ( 0 + me*6 + 5 + 60 )*stride_in];
	(*R40).y = bufInIm[inOffset + ( 0 + me*6 + 5 + 60 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*6 + 0 + 72 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*6 + 0 + 72 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*6 + 1 + 72 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*6 + 1 + 72 )*stride_in];
	(*R20).x = bufInRe[inOffset + ( 0 + me*6 + 2 + 72 )*stride_in];
	(*R20).y = bufInIm[inOffset + ( 0 + me*6 + 2 + 72 )*stride_in];
	(*R27).x = bufInRe[inOffset + ( 0 + me*6 + 3 + 72 )*stride_in];
	(*R27).y = bufInIm[inOffset + ( 0 + me*6 + 3 + 72 )*stride_in];
	(*R34).x = bufInRe[inOffset + ( 0 + me*6 + 4 + 72 )*stride_in];
	(*R34).y = bufInIm[inOffset + ( 0 + me*6 + 4 + 72 )*stride_in];
	(*R41).x = bufInRe[inOffset + ( 0 + me*6 + 5 + 72 )*stride_in];
	(*R41).y = bufInIm[inOffset + ( 0 + me*6 + 5 + 72 )*stride_in];
	}



	FwdRad7B1(R0, R1, R2, R3, R4, R5, R6);
	FwdRad7B1(R7, R8, R9, R10, R11, R12, R13);
	FwdRad7B1(R14, R15, R16, R17, R18, R19, R20);
	FwdRad7B1(R21, R22, R23, R24, R25, R26, R27);
	FwdRad7B1(R28, R29, R30, R31, R32, R33, R34);
	FwdRad7B1(R35, R36, R37, R38, R39, R40, R41);


	if(rw)
	{
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 0 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 1 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 2 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 3 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 4 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 5 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 6 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 0 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 1 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 2 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 3 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 4 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 5 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 6 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 0 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 1 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 2 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 3 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 4 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 5 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 6 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 0 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 1 ) ] = (*R29).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 2 ) ] = (*R30).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 3 ) ] = (*R31).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 4 ) ] = (*R32).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 5 ) ] = (*R33).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 6 ) ] = (*R34).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 0 ) ] = (*R35).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 1 ) ] = (*R36).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 2 ) ] = (*R37).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 3 ) ] = (*R38).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 4 ) ] = (*R39).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 5 ) ] = (*R40).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 6 ) ] = (*R41).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 0 ) ];
	(*R30).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 0 ) ];
	(*R36).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 14 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 14 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 14 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 14 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 14 ) ];
	(*R31).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 14 ) ];
	(*R37).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 14 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 28 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 28 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 28 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 28 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 28 ) ];
	(*R32).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 28 ) ];
	(*R38).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 28 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 42 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 42 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 42 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 42 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 42 ) ];
	(*R33).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 42 ) ];
	(*R39).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 42 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 56 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 56 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 56 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 56 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 56 ) ];
	(*R34).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 56 ) ];
	(*R40).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 56 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 70 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 70 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 70 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 70 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 70 ) ];
	(*R35).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 70 ) ];
	(*R41).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 70 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 0 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 1 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 2 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 3 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 4 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 5 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 6 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 0 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 1 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 2 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 3 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 4 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 5 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 6 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 0 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 1 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 2 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 3 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 4 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 5 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 6 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 0 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 1 ) ] = (*R29).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 2 ) ] = (*R30).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 3 ) ] = (*R31).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 4 ) ] = (*R32).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 5 ) ] = (*R33).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 6 ) ] = (*R34).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 0 ) ] = (*R35).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 1 ) ] = (*R36).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 2 ) ] = (*R37).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 3 ) ] = (*R38).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 4 ) ] = (*R39).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 5 ) ] = (*R40).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 6 ) ] = (*R41).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 0 ) ];
	(*R30).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 0 ) ];
	(*R36).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 14 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 14 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 14 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 14 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 14 ) ];
	(*R31).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 14 ) ];
	(*R37).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 14 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 28 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 28 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 28 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 28 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 28 ) ];
	(*R32).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 28 ) ];
	(*R38).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 28 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 42 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 42 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 42 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 42 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 42 ) ];
	(*R33).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 42 ) ];
	(*R39).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 42 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 56 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 56 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 56 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 56 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 56 ) ];
	(*R34).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 56 ) ];
	(*R40).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 56 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 70 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 70 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 70 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 70 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 70 ) ];
	(*R35).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 70 ) ];
	(*R41).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 70 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass1_len84(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, T *R30, T *R31, T *R32, T *R33, T *R34, T *R35, T *R36, T *R37, T *R38, T *R39, T *R40, T *R41)
{




	{
		T W = twiddles[6 + 5*((7*me + 0)%7) + 0];
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
		T W = twiddles[6 + 5*((7*me + 0)%7) + 1];
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
		T W = twiddles[6 + 5*((7*me + 0)%7) + 2];
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
		T W = twiddles[6 + 5*((7*me + 0)%7) + 3];
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
		T W = twiddles[6 + 5*((7*me + 0)%7) + 4];
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
		T W = twiddles[6 + 5*((7*me + 1)%7) + 0];
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
		T W = twiddles[6 + 5*((7*me + 1)%7) + 1];
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
		T W = twiddles[6 + 5*((7*me + 1)%7) + 2];
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
		T W = twiddles[6 + 5*((7*me + 1)%7) + 3];
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
		T W = twiddles[6 + 5*((7*me + 1)%7) + 4];
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
		T W = twiddles[6 + 5*((7*me + 2)%7) + 0];
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
		T W = twiddles[6 + 5*((7*me + 2)%7) + 1];
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
		T W = twiddles[6 + 5*((7*me + 2)%7) + 2];
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
		T W = twiddles[6 + 5*((7*me + 2)%7) + 3];
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
		T W = twiddles[6 + 5*((7*me + 2)%7) + 4];
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
		T W = twiddles[6 + 5*((7*me + 3)%7) + 0];
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
		T W = twiddles[6 + 5*((7*me + 3)%7) + 1];
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
		T W = twiddles[6 + 5*((7*me + 3)%7) + 2];
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
		T W = twiddles[6 + 5*((7*me + 3)%7) + 3];
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
		T W = twiddles[6 + 5*((7*me + 3)%7) + 4];
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
		T W = twiddles[6 + 5*((7*me + 4)%7) + 0];
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
		T W = twiddles[6 + 5*((7*me + 4)%7) + 1];
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
		T W = twiddles[6 + 5*((7*me + 4)%7) + 2];
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
		T W = twiddles[6 + 5*((7*me + 4)%7) + 3];
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
		T W = twiddles[6 + 5*((7*me + 4)%7) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 5)%7) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R31).x; ry = (*R31).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R31).x = TR;
		(*R31).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 5)%7) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R32).x; ry = (*R32).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R32).x = TR;
		(*R32).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 5)%7) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R33).x; ry = (*R33).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R33).x = TR;
		(*R33).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 5)%7) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R34).x; ry = (*R34).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R34).x = TR;
		(*R34).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 5)%7) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R35).x; ry = (*R35).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R35).x = TR;
		(*R35).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 6)%7) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R37).x; ry = (*R37).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R37).x = TR;
		(*R37).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 6)%7) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R38).x; ry = (*R38).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R38).x = TR;
		(*R38).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 6)%7) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R39).x; ry = (*R39).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R39).x = TR;
		(*R39).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 6)%7) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R40).x; ry = (*R40).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R40).x = TR;
		(*R40).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 6)%7) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R41).x; ry = (*R41).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R41).x = TR;
		(*R41).y = TI;
	}

	FwdRad6B1(R0, R1, R2, R3, R4, R5);
	FwdRad6B1(R6, R7, R8, R9, R10, R11);
	FwdRad6B1(R12, R13, R14, R15, R16, R17);
	FwdRad6B1(R18, R19, R20, R21, R22, R23);
	FwdRad6B1(R24, R25, R26, R27, R28, R29);
	FwdRad6B1(R30, R31, R32, R33, R34, R35);
	FwdRad6B1(R36, R37, R38, R39, R40, R41);


	if(rw)
	{
	bufOutRe[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 7 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 14 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 21 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 28 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 35 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 14 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 21 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 28 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 35 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 7 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 14 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 21 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 28 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 35 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 0 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 7 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 14 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 21 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 28 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 35 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 0 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 7 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 14 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 21 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 28 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 35 ) ] = (*R29).x;
	bufOutRe[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 0 ) ] = (*R30).x;
	bufOutRe[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 7 ) ] = (*R31).x;
	bufOutRe[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 14 ) ] = (*R32).x;
	bufOutRe[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 21 ) ] = (*R33).x;
	bufOutRe[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 28 ) ] = (*R34).x;
	bufOutRe[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 35 ) ] = (*R35).x;
	bufOutRe[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 0 ) ] = (*R36).x;
	bufOutRe[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 7 ) ] = (*R37).x;
	bufOutRe[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 14 ) ] = (*R38).x;
	bufOutRe[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 21 ) ] = (*R39).x;
	bufOutRe[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 28 ) ] = (*R40).x;
	bufOutRe[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 35 ) ] = (*R41).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*21 + 0 + 0 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*21 + 1 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*21 + 2 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*21 + 3 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*21 + 4 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*21 + 5 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*21 + 6 + 0 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*21 + 7 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*21 + 8 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*21 + 9 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*21 + 10 + 0 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*21 + 11 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*21 + 12 + 0 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*21 + 13 + 0 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*21 + 14 + 0 ) ];
	(*R30).x = bufOutRe[outOffset + ( 0 + me*21 + 15 + 0 ) ];
	(*R32).x = bufOutRe[outOffset + ( 0 + me*21 + 16 + 0 ) ];
	(*R34).x = bufOutRe[outOffset + ( 0 + me*21 + 17 + 0 ) ];
	(*R36).x = bufOutRe[outOffset + ( 0 + me*21 + 18 + 0 ) ];
	(*R38).x = bufOutRe[outOffset + ( 0 + me*21 + 19 + 0 ) ];
	(*R40).x = bufOutRe[outOffset + ( 0 + me*21 + 20 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*21 + 0 + 42 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*21 + 1 + 42 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*21 + 2 + 42 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*21 + 3 + 42 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*21 + 4 + 42 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*21 + 5 + 42 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*21 + 6 + 42 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*21 + 7 + 42 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*21 + 8 + 42 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*21 + 9 + 42 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*21 + 10 + 42 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*21 + 11 + 42 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*21 + 12 + 42 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*21 + 13 + 42 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*21 + 14 + 42 ) ];
	(*R31).x = bufOutRe[outOffset + ( 0 + me*21 + 15 + 42 ) ];
	(*R33).x = bufOutRe[outOffset + ( 0 + me*21 + 16 + 42 ) ];
	(*R35).x = bufOutRe[outOffset + ( 0 + me*21 + 17 + 42 ) ];
	(*R37).x = bufOutRe[outOffset + ( 0 + me*21 + 18 + 42 ) ];
	(*R39).x = bufOutRe[outOffset + ( 0 + me*21 + 19 + 42 ) ];
	(*R41).x = bufOutRe[outOffset + ( 0 + me*21 + 20 + 42 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 7 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 14 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 21 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 28 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 35 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 14 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 21 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 28 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 35 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 7 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 14 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 21 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 28 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 35 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 0 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 7 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 14 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 21 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 28 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 35 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 0 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 7 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 14 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 21 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 28 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 35 ) ] = (*R29).y;
	bufOutIm[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 0 ) ] = (*R30).y;
	bufOutIm[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 7 ) ] = (*R31).y;
	bufOutIm[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 14 ) ] = (*R32).y;
	bufOutIm[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 21 ) ] = (*R33).y;
	bufOutIm[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 28 ) ] = (*R34).y;
	bufOutIm[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 35 ) ] = (*R35).y;
	bufOutIm[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 0 ) ] = (*R36).y;
	bufOutIm[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 7 ) ] = (*R37).y;
	bufOutIm[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 14 ) ] = (*R38).y;
	bufOutIm[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 21 ) ] = (*R39).y;
	bufOutIm[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 28 ) ] = (*R40).y;
	bufOutIm[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 35 ) ] = (*R41).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*21 + 0 + 0 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*21 + 1 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*21 + 2 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*21 + 3 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*21 + 4 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*21 + 5 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*21 + 6 + 0 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*21 + 7 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*21 + 8 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*21 + 9 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*21 + 10 + 0 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*21 + 11 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*21 + 12 + 0 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*21 + 13 + 0 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*21 + 14 + 0 ) ];
	(*R30).y = bufOutIm[outOffset + ( 0 + me*21 + 15 + 0 ) ];
	(*R32).y = bufOutIm[outOffset + ( 0 + me*21 + 16 + 0 ) ];
	(*R34).y = bufOutIm[outOffset + ( 0 + me*21 + 17 + 0 ) ];
	(*R36).y = bufOutIm[outOffset + ( 0 + me*21 + 18 + 0 ) ];
	(*R38).y = bufOutIm[outOffset + ( 0 + me*21 + 19 + 0 ) ];
	(*R40).y = bufOutIm[outOffset + ( 0 + me*21 + 20 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*21 + 0 + 42 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*21 + 1 + 42 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*21 + 2 + 42 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*21 + 3 + 42 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*21 + 4 + 42 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*21 + 5 + 42 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*21 + 6 + 42 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*21 + 7 + 42 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*21 + 8 + 42 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*21 + 9 + 42 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*21 + 10 + 42 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*21 + 11 + 42 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*21 + 12 + 42 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*21 + 13 + 42 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*21 + 14 + 42 ) ];
	(*R31).y = bufOutIm[outOffset + ( 0 + me*21 + 15 + 42 ) ];
	(*R33).y = bufOutIm[outOffset + ( 0 + me*21 + 16 + 42 ) ];
	(*R35).y = bufOutIm[outOffset + ( 0 + me*21 + 17 + 42 ) ];
	(*R37).y = bufOutIm[outOffset + ( 0 + me*21 + 18 + 42 ) ];
	(*R39).y = bufOutIm[outOffset + ( 0 + me*21 + 19 + 42 ) ];
	(*R41).y = bufOutIm[outOffset + ( 0 + me*21 + 20 + 42 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
FwdPass2_len84(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, T *R30, T *R31, T *R32, T *R33, T *R34, T *R35, T *R36, T *R37, T *R38, T *R39, T *R40, T *R41)
{




	{
		T W = twiddles[41 + 1*((21*me + 0)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 1)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 2)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 3)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 4)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 5)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 6)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 7)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 8)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 9)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 10)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 11)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 12)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 13)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 14)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 15)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R31).x; ry = (*R31).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R31).x = TR;
		(*R31).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 16)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R33).x; ry = (*R33).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R33).x = TR;
		(*R33).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 17)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R35).x; ry = (*R35).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R35).x = TR;
		(*R35).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 18)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R37).x; ry = (*R37).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R37).x = TR;
		(*R37).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 19)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R39).x; ry = (*R39).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R39).x = TR;
		(*R39).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 20)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R41).x; ry = (*R41).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R41).x = TR;
		(*R41).y = TI;
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
	FwdRad2B1(R26, R27);
	FwdRad2B1(R28, R29);
	FwdRad2B1(R30, R31);
	FwdRad2B1(R32, R33);
	FwdRad2B1(R34, R35);
	FwdRad2B1(R36, R37);
	FwdRad2B1(R38, R39);
	FwdRad2B1(R40, R41);


	if(rw)
	{
	bufOut[outOffset + ( 21*me + 0 + 0 )*stride_out] = (*R0);
	bufOut[outOffset + ( 21*me + 1 + 0 )*stride_out] = (*R2);
	bufOut[outOffset + ( 21*me + 2 + 0 )*stride_out] = (*R4);
	bufOut[outOffset + ( 21*me + 3 + 0 )*stride_out] = (*R6);
	bufOut[outOffset + ( 21*me + 4 + 0 )*stride_out] = (*R8);
	bufOut[outOffset + ( 21*me + 5 + 0 )*stride_out] = (*R10);
	bufOut[outOffset + ( 21*me + 6 + 0 )*stride_out] = (*R12);
	bufOut[outOffset + ( 21*me + 7 + 0 )*stride_out] = (*R14);
	bufOut[outOffset + ( 21*me + 8 + 0 )*stride_out] = (*R16);
	bufOut[outOffset + ( 21*me + 9 + 0 )*stride_out] = (*R18);
	bufOut[outOffset + ( 21*me + 10 + 0 )*stride_out] = (*R20);
	bufOut[outOffset + ( 21*me + 11 + 0 )*stride_out] = (*R22);
	bufOut[outOffset + ( 21*me + 12 + 0 )*stride_out] = (*R24);
	bufOut[outOffset + ( 21*me + 13 + 0 )*stride_out] = (*R26);
	bufOut[outOffset + ( 21*me + 14 + 0 )*stride_out] = (*R28);
	bufOut[outOffset + ( 21*me + 15 + 0 )*stride_out] = (*R30);
	bufOut[outOffset + ( 21*me + 16 + 0 )*stride_out] = (*R32);
	bufOut[outOffset + ( 21*me + 17 + 0 )*stride_out] = (*R34);
	bufOut[outOffset + ( 21*me + 18 + 0 )*stride_out] = (*R36);
	bufOut[outOffset + ( 21*me + 19 + 0 )*stride_out] = (*R38);
	bufOut[outOffset + ( 21*me + 20 + 0 )*stride_out] = (*R40);
	bufOut[outOffset + ( 21*me + 0 + 42 )*stride_out] = (*R1);
	bufOut[outOffset + ( 21*me + 1 + 42 )*stride_out] = (*R3);
	bufOut[outOffset + ( 21*me + 2 + 42 )*stride_out] = (*R5);
	bufOut[outOffset + ( 21*me + 3 + 42 )*stride_out] = (*R7);
	bufOut[outOffset + ( 21*me + 4 + 42 )*stride_out] = (*R9);
	bufOut[outOffset + ( 21*me + 5 + 42 )*stride_out] = (*R11);
	bufOut[outOffset + ( 21*me + 6 + 42 )*stride_out] = (*R13);
	bufOut[outOffset + ( 21*me + 7 + 42 )*stride_out] = (*R15);
	bufOut[outOffset + ( 21*me + 8 + 42 )*stride_out] = (*R17);
	bufOut[outOffset + ( 21*me + 9 + 42 )*stride_out] = (*R19);
	bufOut[outOffset + ( 21*me + 10 + 42 )*stride_out] = (*R21);
	bufOut[outOffset + ( 21*me + 11 + 42 )*stride_out] = (*R23);
	bufOut[outOffset + ( 21*me + 12 + 42 )*stride_out] = (*R25);
	bufOut[outOffset + ( 21*me + 13 + 42 )*stride_out] = (*R27);
	bufOut[outOffset + ( 21*me + 14 + 42 )*stride_out] = (*R29);
	bufOut[outOffset + ( 21*me + 15 + 42 )*stride_out] = (*R31);
	bufOut[outOffset + ( 21*me + 16 + 42 )*stride_out] = (*R33);
	bufOut[outOffset + ( 21*me + 17 + 42 )*stride_out] = (*R35);
	bufOut[outOffset + ( 21*me + 18 + 42 )*stride_out] = (*R37);
	bufOut[outOffset + ( 21*me + 19 + 42 )*stride_out] = (*R39);
	bufOut[outOffset + ( 21*me + 20 + 42 )*stride_out] = (*R41);
	}

}

template <typename T, StrideBin sb>
__device__ void
FwdPass2_len84(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, T *R30, T *R31, T *R32, T *R33, T *R34, T *R35, T *R36, T *R37, T *R38, T *R39, T *R40, T *R41)
{




	{
		T W = twiddles[41 + 1*((21*me + 0)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 1)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 2)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 3)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 4)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 5)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 6)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 7)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 8)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 9)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 10)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 11)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 12)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 13)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 14)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 15)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R31).x; ry = (*R31).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R31).x = TR;
		(*R31).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 16)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R33).x; ry = (*R33).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R33).x = TR;
		(*R33).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 17)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R35).x; ry = (*R35).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R35).x = TR;
		(*R35).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 18)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R37).x; ry = (*R37).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R37).x = TR;
		(*R37).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 19)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R39).x; ry = (*R39).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R39).x = TR;
		(*R39).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 20)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R41).x; ry = (*R41).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R41).x = TR;
		(*R41).y = TI;
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
	FwdRad2B1(R26, R27);
	FwdRad2B1(R28, R29);
	FwdRad2B1(R30, R31);
	FwdRad2B1(R32, R33);
	FwdRad2B1(R34, R35);
	FwdRad2B1(R36, R37);
	FwdRad2B1(R38, R39);
	FwdRad2B1(R40, R41);


	if(rw)
	{
	bufOutRe[outOffset + ( 21*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 21*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 21*me + 1 + 0 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 21*me + 1 + 0 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 21*me + 2 + 0 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 21*me + 2 + 0 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 21*me + 3 + 0 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 21*me + 3 + 0 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 21*me + 4 + 0 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 21*me + 4 + 0 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 21*me + 5 + 0 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 21*me + 5 + 0 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 21*me + 6 + 0 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 21*me + 6 + 0 )*stride_out] = (*R12).y;
	bufOutRe[outOffset + ( 21*me + 7 + 0 )*stride_out] = (*R14).x;
	bufOutIm[outOffset + ( 21*me + 7 + 0 )*stride_out] = (*R14).y;
	bufOutRe[outOffset + ( 21*me + 8 + 0 )*stride_out] = (*R16).x;
	bufOutIm[outOffset + ( 21*me + 8 + 0 )*stride_out] = (*R16).y;
	bufOutRe[outOffset + ( 21*me + 9 + 0 )*stride_out] = (*R18).x;
	bufOutIm[outOffset + ( 21*me + 9 + 0 )*stride_out] = (*R18).y;
	bufOutRe[outOffset + ( 21*me + 10 + 0 )*stride_out] = (*R20).x;
	bufOutIm[outOffset + ( 21*me + 10 + 0 )*stride_out] = (*R20).y;
	bufOutRe[outOffset + ( 21*me + 11 + 0 )*stride_out] = (*R22).x;
	bufOutIm[outOffset + ( 21*me + 11 + 0 )*stride_out] = (*R22).y;
	bufOutRe[outOffset + ( 21*me + 12 + 0 )*stride_out] = (*R24).x;
	bufOutIm[outOffset + ( 21*me + 12 + 0 )*stride_out] = (*R24).y;
	bufOutRe[outOffset + ( 21*me + 13 + 0 )*stride_out] = (*R26).x;
	bufOutIm[outOffset + ( 21*me + 13 + 0 )*stride_out] = (*R26).y;
	bufOutRe[outOffset + ( 21*me + 14 + 0 )*stride_out] = (*R28).x;
	bufOutIm[outOffset + ( 21*me + 14 + 0 )*stride_out] = (*R28).y;
	bufOutRe[outOffset + ( 21*me + 15 + 0 )*stride_out] = (*R30).x;
	bufOutIm[outOffset + ( 21*me + 15 + 0 )*stride_out] = (*R30).y;
	bufOutRe[outOffset + ( 21*me + 16 + 0 )*stride_out] = (*R32).x;
	bufOutIm[outOffset + ( 21*me + 16 + 0 )*stride_out] = (*R32).y;
	bufOutRe[outOffset + ( 21*me + 17 + 0 )*stride_out] = (*R34).x;
	bufOutIm[outOffset + ( 21*me + 17 + 0 )*stride_out] = (*R34).y;
	bufOutRe[outOffset + ( 21*me + 18 + 0 )*stride_out] = (*R36).x;
	bufOutIm[outOffset + ( 21*me + 18 + 0 )*stride_out] = (*R36).y;
	bufOutRe[outOffset + ( 21*me + 19 + 0 )*stride_out] = (*R38).x;
	bufOutIm[outOffset + ( 21*me + 19 + 0 )*stride_out] = (*R38).y;
	bufOutRe[outOffset + ( 21*me + 20 + 0 )*stride_out] = (*R40).x;
	bufOutIm[outOffset + ( 21*me + 20 + 0 )*stride_out] = (*R40).y;
	bufOutRe[outOffset + ( 21*me + 0 + 42 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 21*me + 0 + 42 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 21*me + 1 + 42 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 21*me + 1 + 42 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 21*me + 2 + 42 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 21*me + 2 + 42 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 21*me + 3 + 42 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 21*me + 3 + 42 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 21*me + 4 + 42 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 21*me + 4 + 42 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 21*me + 5 + 42 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 21*me + 5 + 42 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 21*me + 6 + 42 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 21*me + 6 + 42 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 21*me + 7 + 42 )*stride_out] = (*R15).x;
	bufOutIm[outOffset + ( 21*me + 7 + 42 )*stride_out] = (*R15).y;
	bufOutRe[outOffset + ( 21*me + 8 + 42 )*stride_out] = (*R17).x;
	bufOutIm[outOffset + ( 21*me + 8 + 42 )*stride_out] = (*R17).y;
	bufOutRe[outOffset + ( 21*me + 9 + 42 )*stride_out] = (*R19).x;
	bufOutIm[outOffset + ( 21*me + 9 + 42 )*stride_out] = (*R19).y;
	bufOutRe[outOffset + ( 21*me + 10 + 42 )*stride_out] = (*R21).x;
	bufOutIm[outOffset + ( 21*me + 10 + 42 )*stride_out] = (*R21).y;
	bufOutRe[outOffset + ( 21*me + 11 + 42 )*stride_out] = (*R23).x;
	bufOutIm[outOffset + ( 21*me + 11 + 42 )*stride_out] = (*R23).y;
	bufOutRe[outOffset + ( 21*me + 12 + 42 )*stride_out] = (*R25).x;
	bufOutIm[outOffset + ( 21*me + 12 + 42 )*stride_out] = (*R25).y;
	bufOutRe[outOffset + ( 21*me + 13 + 42 )*stride_out] = (*R27).x;
	bufOutIm[outOffset + ( 21*me + 13 + 42 )*stride_out] = (*R27).y;
	bufOutRe[outOffset + ( 21*me + 14 + 42 )*stride_out] = (*R29).x;
	bufOutIm[outOffset + ( 21*me + 14 + 42 )*stride_out] = (*R29).y;
	bufOutRe[outOffset + ( 21*me + 15 + 42 )*stride_out] = (*R31).x;
	bufOutIm[outOffset + ( 21*me + 15 + 42 )*stride_out] = (*R31).y;
	bufOutRe[outOffset + ( 21*me + 16 + 42 )*stride_out] = (*R33).x;
	bufOutIm[outOffset + ( 21*me + 16 + 42 )*stride_out] = (*R33).y;
	bufOutRe[outOffset + ( 21*me + 17 + 42 )*stride_out] = (*R35).x;
	bufOutIm[outOffset + ( 21*me + 17 + 42 )*stride_out] = (*R35).y;
	bufOutRe[outOffset + ( 21*me + 18 + 42 )*stride_out] = (*R37).x;
	bufOutIm[outOffset + ( 21*me + 18 + 42 )*stride_out] = (*R37).y;
	bufOutRe[outOffset + ( 21*me + 19 + 42 )*stride_out] = (*R39).x;
	bufOutIm[outOffset + ( 21*me + 19 + 42 )*stride_out] = (*R39).y;
	bufOutRe[outOffset + ( 21*me + 20 + 42 )*stride_out] = (*R41).x;
	bufOutIm[outOffset + ( 21*me + 20 + 42 )*stride_out] = (*R41).y;
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass0_len84(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, T *R30, T *R31, T *R32, T *R33, T *R34, T *R35, T *R36, T *R37, T *R38, T *R39, T *R40, T *R41)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*6 + 0 + 0 )*stride_in];
	(*R7) = bufIn[inOffset + ( 0 + me*6 + 1 + 0 )*stride_in];
	(*R14) = bufIn[inOffset + ( 0 + me*6 + 2 + 0 )*stride_in];
	(*R21) = bufIn[inOffset + ( 0 + me*6 + 3 + 0 )*stride_in];
	(*R28) = bufIn[inOffset + ( 0 + me*6 + 4 + 0 )*stride_in];
	(*R35) = bufIn[inOffset + ( 0 + me*6 + 5 + 0 )*stride_in];
	(*R1) = bufIn[inOffset + ( 0 + me*6 + 0 + 12 )*stride_in];
	(*R8) = bufIn[inOffset + ( 0 + me*6 + 1 + 12 )*stride_in];
	(*R15) = bufIn[inOffset + ( 0 + me*6 + 2 + 12 )*stride_in];
	(*R22) = bufIn[inOffset + ( 0 + me*6 + 3 + 12 )*stride_in];
	(*R29) = bufIn[inOffset + ( 0 + me*6 + 4 + 12 )*stride_in];
	(*R36) = bufIn[inOffset + ( 0 + me*6 + 5 + 12 )*stride_in];
	(*R2) = bufIn[inOffset + ( 0 + me*6 + 0 + 24 )*stride_in];
	(*R9) = bufIn[inOffset + ( 0 + me*6 + 1 + 24 )*stride_in];
	(*R16) = bufIn[inOffset + ( 0 + me*6 + 2 + 24 )*stride_in];
	(*R23) = bufIn[inOffset + ( 0 + me*6 + 3 + 24 )*stride_in];
	(*R30) = bufIn[inOffset + ( 0 + me*6 + 4 + 24 )*stride_in];
	(*R37) = bufIn[inOffset + ( 0 + me*6 + 5 + 24 )*stride_in];
	(*R3) = bufIn[inOffset + ( 0 + me*6 + 0 + 36 )*stride_in];
	(*R10) = bufIn[inOffset + ( 0 + me*6 + 1 + 36 )*stride_in];
	(*R17) = bufIn[inOffset + ( 0 + me*6 + 2 + 36 )*stride_in];
	(*R24) = bufIn[inOffset + ( 0 + me*6 + 3 + 36 )*stride_in];
	(*R31) = bufIn[inOffset + ( 0 + me*6 + 4 + 36 )*stride_in];
	(*R38) = bufIn[inOffset + ( 0 + me*6 + 5 + 36 )*stride_in];
	(*R4) = bufIn[inOffset + ( 0 + me*6 + 0 + 48 )*stride_in];
	(*R11) = bufIn[inOffset + ( 0 + me*6 + 1 + 48 )*stride_in];
	(*R18) = bufIn[inOffset + ( 0 + me*6 + 2 + 48 )*stride_in];
	(*R25) = bufIn[inOffset + ( 0 + me*6 + 3 + 48 )*stride_in];
	(*R32) = bufIn[inOffset + ( 0 + me*6 + 4 + 48 )*stride_in];
	(*R39) = bufIn[inOffset + ( 0 + me*6 + 5 + 48 )*stride_in];
	(*R5) = bufIn[inOffset + ( 0 + me*6 + 0 + 60 )*stride_in];
	(*R12) = bufIn[inOffset + ( 0 + me*6 + 1 + 60 )*stride_in];
	(*R19) = bufIn[inOffset + ( 0 + me*6 + 2 + 60 )*stride_in];
	(*R26) = bufIn[inOffset + ( 0 + me*6 + 3 + 60 )*stride_in];
	(*R33) = bufIn[inOffset + ( 0 + me*6 + 4 + 60 )*stride_in];
	(*R40) = bufIn[inOffset + ( 0 + me*6 + 5 + 60 )*stride_in];
	(*R6) = bufIn[inOffset + ( 0 + me*6 + 0 + 72 )*stride_in];
	(*R13) = bufIn[inOffset + ( 0 + me*6 + 1 + 72 )*stride_in];
	(*R20) = bufIn[inOffset + ( 0 + me*6 + 2 + 72 )*stride_in];
	(*R27) = bufIn[inOffset + ( 0 + me*6 + 3 + 72 )*stride_in];
	(*R34) = bufIn[inOffset + ( 0 + me*6 + 4 + 72 )*stride_in];
	(*R41) = bufIn[inOffset + ( 0 + me*6 + 5 + 72 )*stride_in];
	}



	InvRad7B1(R0, R1, R2, R3, R4, R5, R6);
	InvRad7B1(R7, R8, R9, R10, R11, R12, R13);
	InvRad7B1(R14, R15, R16, R17, R18, R19, R20);
	InvRad7B1(R21, R22, R23, R24, R25, R26, R27);
	InvRad7B1(R28, R29, R30, R31, R32, R33, R34);
	InvRad7B1(R35, R36, R37, R38, R39, R40, R41);


	if(rw)
	{
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 0 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 1 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 2 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 3 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 4 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 5 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 6 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 0 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 1 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 2 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 3 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 4 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 5 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 6 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 0 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 1 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 2 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 3 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 4 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 5 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 6 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 0 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 1 ) ] = (*R29).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 2 ) ] = (*R30).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 3 ) ] = (*R31).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 4 ) ] = (*R32).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 5 ) ] = (*R33).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 6 ) ] = (*R34).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 0 ) ] = (*R35).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 1 ) ] = (*R36).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 2 ) ] = (*R37).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 3 ) ] = (*R38).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 4 ) ] = (*R39).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 5 ) ] = (*R40).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 6 ) ] = (*R41).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 0 ) ];
	(*R30).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 0 ) ];
	(*R36).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 14 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 14 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 14 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 14 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 14 ) ];
	(*R31).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 14 ) ];
	(*R37).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 14 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 28 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 28 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 28 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 28 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 28 ) ];
	(*R32).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 28 ) ];
	(*R38).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 28 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 42 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 42 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 42 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 42 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 42 ) ];
	(*R33).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 42 ) ];
	(*R39).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 42 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 56 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 56 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 56 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 56 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 56 ) ];
	(*R34).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 56 ) ];
	(*R40).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 56 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 70 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 70 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 70 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 70 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 70 ) ];
	(*R35).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 70 ) ];
	(*R41).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 70 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 0 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 1 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 2 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 3 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 4 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 5 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 6 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 0 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 1 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 2 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 3 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 4 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 5 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 6 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 0 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 1 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 2 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 3 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 4 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 5 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 6 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 0 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 1 ) ] = (*R29).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 2 ) ] = (*R30).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 3 ) ] = (*R31).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 4 ) ] = (*R32).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 5 ) ] = (*R33).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 6 ) ] = (*R34).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 0 ) ] = (*R35).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 1 ) ] = (*R36).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 2 ) ] = (*R37).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 3 ) ] = (*R38).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 4 ) ] = (*R39).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 5 ) ] = (*R40).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 6 ) ] = (*R41).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 0 ) ];
	(*R30).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 0 ) ];
	(*R36).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 14 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 14 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 14 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 14 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 14 ) ];
	(*R31).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 14 ) ];
	(*R37).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 14 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 28 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 28 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 28 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 28 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 28 ) ];
	(*R32).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 28 ) ];
	(*R38).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 28 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 42 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 42 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 42 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 42 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 42 ) ];
	(*R33).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 42 ) ];
	(*R39).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 42 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 56 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 56 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 56 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 56 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 56 ) ];
	(*R34).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 56 ) ];
	(*R40).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 56 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 70 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 70 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 70 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 70 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 70 ) ];
	(*R35).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 70 ) ];
	(*R41).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 70 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass0_len84(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, T *R30, T *R31, T *R32, T *R33, T *R34, T *R35, T *R36, T *R37, T *R38, T *R39, T *R40, T *R41)
{


	if(rw)
	{
	(*R0).x = bufInRe[inOffset + ( 0 + me*6 + 0 + 0 )*stride_in];
	(*R0).y = bufInIm[inOffset + ( 0 + me*6 + 0 + 0 )*stride_in];
	(*R7).x = bufInRe[inOffset + ( 0 + me*6 + 1 + 0 )*stride_in];
	(*R7).y = bufInIm[inOffset + ( 0 + me*6 + 1 + 0 )*stride_in];
	(*R14).x = bufInRe[inOffset + ( 0 + me*6 + 2 + 0 )*stride_in];
	(*R14).y = bufInIm[inOffset + ( 0 + me*6 + 2 + 0 )*stride_in];
	(*R21).x = bufInRe[inOffset + ( 0 + me*6 + 3 + 0 )*stride_in];
	(*R21).y = bufInIm[inOffset + ( 0 + me*6 + 3 + 0 )*stride_in];
	(*R28).x = bufInRe[inOffset + ( 0 + me*6 + 4 + 0 )*stride_in];
	(*R28).y = bufInIm[inOffset + ( 0 + me*6 + 4 + 0 )*stride_in];
	(*R35).x = bufInRe[inOffset + ( 0 + me*6 + 5 + 0 )*stride_in];
	(*R35).y = bufInIm[inOffset + ( 0 + me*6 + 5 + 0 )*stride_in];
	(*R1).x = bufInRe[inOffset + ( 0 + me*6 + 0 + 12 )*stride_in];
	(*R1).y = bufInIm[inOffset + ( 0 + me*6 + 0 + 12 )*stride_in];
	(*R8).x = bufInRe[inOffset + ( 0 + me*6 + 1 + 12 )*stride_in];
	(*R8).y = bufInIm[inOffset + ( 0 + me*6 + 1 + 12 )*stride_in];
	(*R15).x = bufInRe[inOffset + ( 0 + me*6 + 2 + 12 )*stride_in];
	(*R15).y = bufInIm[inOffset + ( 0 + me*6 + 2 + 12 )*stride_in];
	(*R22).x = bufInRe[inOffset + ( 0 + me*6 + 3 + 12 )*stride_in];
	(*R22).y = bufInIm[inOffset + ( 0 + me*6 + 3 + 12 )*stride_in];
	(*R29).x = bufInRe[inOffset + ( 0 + me*6 + 4 + 12 )*stride_in];
	(*R29).y = bufInIm[inOffset + ( 0 + me*6 + 4 + 12 )*stride_in];
	(*R36).x = bufInRe[inOffset + ( 0 + me*6 + 5 + 12 )*stride_in];
	(*R36).y = bufInIm[inOffset + ( 0 + me*6 + 5 + 12 )*stride_in];
	(*R2).x = bufInRe[inOffset + ( 0 + me*6 + 0 + 24 )*stride_in];
	(*R2).y = bufInIm[inOffset + ( 0 + me*6 + 0 + 24 )*stride_in];
	(*R9).x = bufInRe[inOffset + ( 0 + me*6 + 1 + 24 )*stride_in];
	(*R9).y = bufInIm[inOffset + ( 0 + me*6 + 1 + 24 )*stride_in];
	(*R16).x = bufInRe[inOffset + ( 0 + me*6 + 2 + 24 )*stride_in];
	(*R16).y = bufInIm[inOffset + ( 0 + me*6 + 2 + 24 )*stride_in];
	(*R23).x = bufInRe[inOffset + ( 0 + me*6 + 3 + 24 )*stride_in];
	(*R23).y = bufInIm[inOffset + ( 0 + me*6 + 3 + 24 )*stride_in];
	(*R30).x = bufInRe[inOffset + ( 0 + me*6 + 4 + 24 )*stride_in];
	(*R30).y = bufInIm[inOffset + ( 0 + me*6 + 4 + 24 )*stride_in];
	(*R37).x = bufInRe[inOffset + ( 0 + me*6 + 5 + 24 )*stride_in];
	(*R37).y = bufInIm[inOffset + ( 0 + me*6 + 5 + 24 )*stride_in];
	(*R3).x = bufInRe[inOffset + ( 0 + me*6 + 0 + 36 )*stride_in];
	(*R3).y = bufInIm[inOffset + ( 0 + me*6 + 0 + 36 )*stride_in];
	(*R10).x = bufInRe[inOffset + ( 0 + me*6 + 1 + 36 )*stride_in];
	(*R10).y = bufInIm[inOffset + ( 0 + me*6 + 1 + 36 )*stride_in];
	(*R17).x = bufInRe[inOffset + ( 0 + me*6 + 2 + 36 )*stride_in];
	(*R17).y = bufInIm[inOffset + ( 0 + me*6 + 2 + 36 )*stride_in];
	(*R24).x = bufInRe[inOffset + ( 0 + me*6 + 3 + 36 )*stride_in];
	(*R24).y = bufInIm[inOffset + ( 0 + me*6 + 3 + 36 )*stride_in];
	(*R31).x = bufInRe[inOffset + ( 0 + me*6 + 4 + 36 )*stride_in];
	(*R31).y = bufInIm[inOffset + ( 0 + me*6 + 4 + 36 )*stride_in];
	(*R38).x = bufInRe[inOffset + ( 0 + me*6 + 5 + 36 )*stride_in];
	(*R38).y = bufInIm[inOffset + ( 0 + me*6 + 5 + 36 )*stride_in];
	(*R4).x = bufInRe[inOffset + ( 0 + me*6 + 0 + 48 )*stride_in];
	(*R4).y = bufInIm[inOffset + ( 0 + me*6 + 0 + 48 )*stride_in];
	(*R11).x = bufInRe[inOffset + ( 0 + me*6 + 1 + 48 )*stride_in];
	(*R11).y = bufInIm[inOffset + ( 0 + me*6 + 1 + 48 )*stride_in];
	(*R18).x = bufInRe[inOffset + ( 0 + me*6 + 2 + 48 )*stride_in];
	(*R18).y = bufInIm[inOffset + ( 0 + me*6 + 2 + 48 )*stride_in];
	(*R25).x = bufInRe[inOffset + ( 0 + me*6 + 3 + 48 )*stride_in];
	(*R25).y = bufInIm[inOffset + ( 0 + me*6 + 3 + 48 )*stride_in];
	(*R32).x = bufInRe[inOffset + ( 0 + me*6 + 4 + 48 )*stride_in];
	(*R32).y = bufInIm[inOffset + ( 0 + me*6 + 4 + 48 )*stride_in];
	(*R39).x = bufInRe[inOffset + ( 0 + me*6 + 5 + 48 )*stride_in];
	(*R39).y = bufInIm[inOffset + ( 0 + me*6 + 5 + 48 )*stride_in];
	(*R5).x = bufInRe[inOffset + ( 0 + me*6 + 0 + 60 )*stride_in];
	(*R5).y = bufInIm[inOffset + ( 0 + me*6 + 0 + 60 )*stride_in];
	(*R12).x = bufInRe[inOffset + ( 0 + me*6 + 1 + 60 )*stride_in];
	(*R12).y = bufInIm[inOffset + ( 0 + me*6 + 1 + 60 )*stride_in];
	(*R19).x = bufInRe[inOffset + ( 0 + me*6 + 2 + 60 )*stride_in];
	(*R19).y = bufInIm[inOffset + ( 0 + me*6 + 2 + 60 )*stride_in];
	(*R26).x = bufInRe[inOffset + ( 0 + me*6 + 3 + 60 )*stride_in];
	(*R26).y = bufInIm[inOffset + ( 0 + me*6 + 3 + 60 )*stride_in];
	(*R33).x = bufInRe[inOffset + ( 0 + me*6 + 4 + 60 )*stride_in];
	(*R33).y = bufInIm[inOffset + ( 0 + me*6 + 4 + 60 )*stride_in];
	(*R40).x = bufInRe[inOffset + ( 0 + me*6 + 5 + 60 )*stride_in];
	(*R40).y = bufInIm[inOffset + ( 0 + me*6 + 5 + 60 )*stride_in];
	(*R6).x = bufInRe[inOffset + ( 0 + me*6 + 0 + 72 )*stride_in];
	(*R6).y = bufInIm[inOffset + ( 0 + me*6 + 0 + 72 )*stride_in];
	(*R13).x = bufInRe[inOffset + ( 0 + me*6 + 1 + 72 )*stride_in];
	(*R13).y = bufInIm[inOffset + ( 0 + me*6 + 1 + 72 )*stride_in];
	(*R20).x = bufInRe[inOffset + ( 0 + me*6 + 2 + 72 )*stride_in];
	(*R20).y = bufInIm[inOffset + ( 0 + me*6 + 2 + 72 )*stride_in];
	(*R27).x = bufInRe[inOffset + ( 0 + me*6 + 3 + 72 )*stride_in];
	(*R27).y = bufInIm[inOffset + ( 0 + me*6 + 3 + 72 )*stride_in];
	(*R34).x = bufInRe[inOffset + ( 0 + me*6 + 4 + 72 )*stride_in];
	(*R34).y = bufInIm[inOffset + ( 0 + me*6 + 4 + 72 )*stride_in];
	(*R41).x = bufInRe[inOffset + ( 0 + me*6 + 5 + 72 )*stride_in];
	(*R41).y = bufInIm[inOffset + ( 0 + me*6 + 5 + 72 )*stride_in];
	}



	InvRad7B1(R0, R1, R2, R3, R4, R5, R6);
	InvRad7B1(R7, R8, R9, R10, R11, R12, R13);
	InvRad7B1(R14, R15, R16, R17, R18, R19, R20);
	InvRad7B1(R21, R22, R23, R24, R25, R26, R27);
	InvRad7B1(R28, R29, R30, R31, R32, R33, R34);
	InvRad7B1(R35, R36, R37, R38, R39, R40, R41);


	if(rw)
	{
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 1 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 2 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 3 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 4 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 5 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 6 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 0 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 1 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 2 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 3 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 4 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 5 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 6 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 0 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 1 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 2 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 3 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 4 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 5 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 6 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 0 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 1 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 2 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 3 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 4 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 5 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 6 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 0 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 1 ) ] = (*R29).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 2 ) ] = (*R30).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 3 ) ] = (*R31).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 4 ) ] = (*R32).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 5 ) ] = (*R33).x;
	bufOutRe[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 6 ) ] = (*R34).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 0 ) ] = (*R35).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 1 ) ] = (*R36).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 2 ) ] = (*R37).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 3 ) ] = (*R38).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 4 ) ] = (*R39).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 5 ) ] = (*R40).x;
	bufOutRe[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 6 ) ] = (*R41).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 0 ) ];
	(*R30).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 0 ) ];
	(*R36).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 14 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 14 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 14 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 14 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 14 ) ];
	(*R31).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 14 ) ];
	(*R37).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 14 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 28 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 28 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 28 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 28 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 28 ) ];
	(*R32).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 28 ) ];
	(*R38).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 28 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 42 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 42 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 42 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 42 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 42 ) ];
	(*R33).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 42 ) ];
	(*R39).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 42 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 56 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 56 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 56 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 56 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 56 ) ];
	(*R34).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 56 ) ];
	(*R40).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 56 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*7 + 0 + 70 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*7 + 1 + 70 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*7 + 2 + 70 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*7 + 3 + 70 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*7 + 4 + 70 ) ];
	(*R35).x = bufOutRe[outOffset + ( 0 + me*7 + 5 + 70 ) ];
	(*R41).x = bufOutRe[outOffset + ( 0 + me*7 + 6 + 70 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 1 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 2 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 3 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 4 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 5 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((6*me + 0)/1)*7 + (6*me + 0)%1 + 6 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 0 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 1 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 2 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 3 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 4 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 5 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((6*me + 1)/1)*7 + (6*me + 1)%1 + 6 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 0 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 1 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 2 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 3 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 4 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 5 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((6*me + 2)/1)*7 + (6*me + 2)%1 + 6 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 0 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 1 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 2 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 3 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 4 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 5 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((6*me + 3)/1)*7 + (6*me + 3)%1 + 6 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 0 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 1 ) ] = (*R29).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 2 ) ] = (*R30).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 3 ) ] = (*R31).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 4 ) ] = (*R32).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 5 ) ] = (*R33).y;
	bufOutIm[outOffset + ( ((6*me + 4)/1)*7 + (6*me + 4)%1 + 6 ) ] = (*R34).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 0 ) ] = (*R35).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 1 ) ] = (*R36).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 2 ) ] = (*R37).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 3 ) ] = (*R38).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 4 ) ] = (*R39).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 5 ) ] = (*R40).y;
	bufOutIm[outOffset + ( ((6*me + 5)/1)*7 + (6*me + 5)%1 + 6 ) ] = (*R41).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 0 ) ];
	(*R30).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 0 ) ];
	(*R36).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 14 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 14 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 14 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 14 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 14 ) ];
	(*R31).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 14 ) ];
	(*R37).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 14 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 28 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 28 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 28 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 28 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 28 ) ];
	(*R32).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 28 ) ];
	(*R38).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 28 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 42 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 42 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 42 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 42 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 42 ) ];
	(*R33).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 42 ) ];
	(*R39).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 42 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 56 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 56 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 56 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 56 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 56 ) ];
	(*R34).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 56 ) ];
	(*R40).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 56 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*7 + 0 + 70 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*7 + 1 + 70 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*7 + 2 + 70 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*7 + 3 + 70 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*7 + 4 + 70 ) ];
	(*R35).y = bufOutIm[outOffset + ( 0 + me*7 + 5 + 70 ) ];
	(*R41).y = bufOutIm[outOffset + ( 0 + me*7 + 6 + 70 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass1_len84(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, T *R30, T *R31, T *R32, T *R33, T *R34, T *R35, T *R36, T *R37, T *R38, T *R39, T *R40, T *R41)
{




	{
		T W = twiddles[6 + 5*((7*me + 0)%7) + 0];
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
		T W = twiddles[6 + 5*((7*me + 0)%7) + 1];
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
		T W = twiddles[6 + 5*((7*me + 0)%7) + 2];
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
		T W = twiddles[6 + 5*((7*me + 0)%7) + 3];
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
		T W = twiddles[6 + 5*((7*me + 0)%7) + 4];
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
		T W = twiddles[6 + 5*((7*me + 1)%7) + 0];
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
		T W = twiddles[6 + 5*((7*me + 1)%7) + 1];
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
		T W = twiddles[6 + 5*((7*me + 1)%7) + 2];
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
		T W = twiddles[6 + 5*((7*me + 1)%7) + 3];
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
		T W = twiddles[6 + 5*((7*me + 1)%7) + 4];
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
		T W = twiddles[6 + 5*((7*me + 2)%7) + 0];
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
		T W = twiddles[6 + 5*((7*me + 2)%7) + 1];
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
		T W = twiddles[6 + 5*((7*me + 2)%7) + 2];
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
		T W = twiddles[6 + 5*((7*me + 2)%7) + 3];
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
		T W = twiddles[6 + 5*((7*me + 2)%7) + 4];
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
		T W = twiddles[6 + 5*((7*me + 3)%7) + 0];
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
		T W = twiddles[6 + 5*((7*me + 3)%7) + 1];
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
		T W = twiddles[6 + 5*((7*me + 3)%7) + 2];
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
		T W = twiddles[6 + 5*((7*me + 3)%7) + 3];
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
		T W = twiddles[6 + 5*((7*me + 3)%7) + 4];
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
		T W = twiddles[6 + 5*((7*me + 4)%7) + 0];
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
		T W = twiddles[6 + 5*((7*me + 4)%7) + 1];
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
		T W = twiddles[6 + 5*((7*me + 4)%7) + 2];
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
		T W = twiddles[6 + 5*((7*me + 4)%7) + 3];
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
		T W = twiddles[6 + 5*((7*me + 4)%7) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 5)%7) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R31).x; ry = (*R31).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R31).x = TR;
		(*R31).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 5)%7) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R32).x; ry = (*R32).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R32).x = TR;
		(*R32).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 5)%7) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R33).x; ry = (*R33).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R33).x = TR;
		(*R33).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 5)%7) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R34).x; ry = (*R34).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R34).x = TR;
		(*R34).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 5)%7) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R35).x; ry = (*R35).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R35).x = TR;
		(*R35).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 6)%7) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R37).x; ry = (*R37).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R37).x = TR;
		(*R37).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 6)%7) + 1];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R38).x; ry = (*R38).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R38).x = TR;
		(*R38).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 6)%7) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R39).x; ry = (*R39).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R39).x = TR;
		(*R39).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 6)%7) + 3];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R40).x; ry = (*R40).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R40).x = TR;
		(*R40).y = TI;
	}

	{
		T W = twiddles[6 + 5*((7*me + 6)%7) + 4];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R41).x; ry = (*R41).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R41).x = TR;
		(*R41).y = TI;
	}

	InvRad6B1(R0, R1, R2, R3, R4, R5);
	InvRad6B1(R6, R7, R8, R9, R10, R11);
	InvRad6B1(R12, R13, R14, R15, R16, R17);
	InvRad6B1(R18, R19, R20, R21, R22, R23);
	InvRad6B1(R24, R25, R26, R27, R28, R29);
	InvRad6B1(R30, R31, R32, R33, R34, R35);
	InvRad6B1(R36, R37, R38, R39, R40, R41);


	if(rw)
	{
	bufOutRe[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 0 ) ] = (*R0).x;
	bufOutRe[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 7 ) ] = (*R1).x;
	bufOutRe[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 14 ) ] = (*R2).x;
	bufOutRe[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 21 ) ] = (*R3).x;
	bufOutRe[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 28 ) ] = (*R4).x;
	bufOutRe[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 35 ) ] = (*R5).x;
	bufOutRe[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 0 ) ] = (*R6).x;
	bufOutRe[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 7 ) ] = (*R7).x;
	bufOutRe[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 14 ) ] = (*R8).x;
	bufOutRe[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 21 ) ] = (*R9).x;
	bufOutRe[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 28 ) ] = (*R10).x;
	bufOutRe[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 35 ) ] = (*R11).x;
	bufOutRe[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 0 ) ] = (*R12).x;
	bufOutRe[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 7 ) ] = (*R13).x;
	bufOutRe[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 14 ) ] = (*R14).x;
	bufOutRe[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 21 ) ] = (*R15).x;
	bufOutRe[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 28 ) ] = (*R16).x;
	bufOutRe[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 35 ) ] = (*R17).x;
	bufOutRe[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 0 ) ] = (*R18).x;
	bufOutRe[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 7 ) ] = (*R19).x;
	bufOutRe[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 14 ) ] = (*R20).x;
	bufOutRe[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 21 ) ] = (*R21).x;
	bufOutRe[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 28 ) ] = (*R22).x;
	bufOutRe[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 35 ) ] = (*R23).x;
	bufOutRe[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 0 ) ] = (*R24).x;
	bufOutRe[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 7 ) ] = (*R25).x;
	bufOutRe[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 14 ) ] = (*R26).x;
	bufOutRe[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 21 ) ] = (*R27).x;
	bufOutRe[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 28 ) ] = (*R28).x;
	bufOutRe[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 35 ) ] = (*R29).x;
	bufOutRe[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 0 ) ] = (*R30).x;
	bufOutRe[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 7 ) ] = (*R31).x;
	bufOutRe[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 14 ) ] = (*R32).x;
	bufOutRe[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 21 ) ] = (*R33).x;
	bufOutRe[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 28 ) ] = (*R34).x;
	bufOutRe[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 35 ) ] = (*R35).x;
	bufOutRe[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 0 ) ] = (*R36).x;
	bufOutRe[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 7 ) ] = (*R37).x;
	bufOutRe[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 14 ) ] = (*R38).x;
	bufOutRe[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 21 ) ] = (*R39).x;
	bufOutRe[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 28 ) ] = (*R40).x;
	bufOutRe[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 35 ) ] = (*R41).x;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).x = bufOutRe[outOffset + ( 0 + me*21 + 0 + 0 ) ];
	(*R2).x = bufOutRe[outOffset + ( 0 + me*21 + 1 + 0 ) ];
	(*R4).x = bufOutRe[outOffset + ( 0 + me*21 + 2 + 0 ) ];
	(*R6).x = bufOutRe[outOffset + ( 0 + me*21 + 3 + 0 ) ];
	(*R8).x = bufOutRe[outOffset + ( 0 + me*21 + 4 + 0 ) ];
	(*R10).x = bufOutRe[outOffset + ( 0 + me*21 + 5 + 0 ) ];
	(*R12).x = bufOutRe[outOffset + ( 0 + me*21 + 6 + 0 ) ];
	(*R14).x = bufOutRe[outOffset + ( 0 + me*21 + 7 + 0 ) ];
	(*R16).x = bufOutRe[outOffset + ( 0 + me*21 + 8 + 0 ) ];
	(*R18).x = bufOutRe[outOffset + ( 0 + me*21 + 9 + 0 ) ];
	(*R20).x = bufOutRe[outOffset + ( 0 + me*21 + 10 + 0 ) ];
	(*R22).x = bufOutRe[outOffset + ( 0 + me*21 + 11 + 0 ) ];
	(*R24).x = bufOutRe[outOffset + ( 0 + me*21 + 12 + 0 ) ];
	(*R26).x = bufOutRe[outOffset + ( 0 + me*21 + 13 + 0 ) ];
	(*R28).x = bufOutRe[outOffset + ( 0 + me*21 + 14 + 0 ) ];
	(*R30).x = bufOutRe[outOffset + ( 0 + me*21 + 15 + 0 ) ];
	(*R32).x = bufOutRe[outOffset + ( 0 + me*21 + 16 + 0 ) ];
	(*R34).x = bufOutRe[outOffset + ( 0 + me*21 + 17 + 0 ) ];
	(*R36).x = bufOutRe[outOffset + ( 0 + me*21 + 18 + 0 ) ];
	(*R38).x = bufOutRe[outOffset + ( 0 + me*21 + 19 + 0 ) ];
	(*R40).x = bufOutRe[outOffset + ( 0 + me*21 + 20 + 0 ) ];
	(*R1).x = bufOutRe[outOffset + ( 0 + me*21 + 0 + 42 ) ];
	(*R3).x = bufOutRe[outOffset + ( 0 + me*21 + 1 + 42 ) ];
	(*R5).x = bufOutRe[outOffset + ( 0 + me*21 + 2 + 42 ) ];
	(*R7).x = bufOutRe[outOffset + ( 0 + me*21 + 3 + 42 ) ];
	(*R9).x = bufOutRe[outOffset + ( 0 + me*21 + 4 + 42 ) ];
	(*R11).x = bufOutRe[outOffset + ( 0 + me*21 + 5 + 42 ) ];
	(*R13).x = bufOutRe[outOffset + ( 0 + me*21 + 6 + 42 ) ];
	(*R15).x = bufOutRe[outOffset + ( 0 + me*21 + 7 + 42 ) ];
	(*R17).x = bufOutRe[outOffset + ( 0 + me*21 + 8 + 42 ) ];
	(*R19).x = bufOutRe[outOffset + ( 0 + me*21 + 9 + 42 ) ];
	(*R21).x = bufOutRe[outOffset + ( 0 + me*21 + 10 + 42 ) ];
	(*R23).x = bufOutRe[outOffset + ( 0 + me*21 + 11 + 42 ) ];
	(*R25).x = bufOutRe[outOffset + ( 0 + me*21 + 12 + 42 ) ];
	(*R27).x = bufOutRe[outOffset + ( 0 + me*21 + 13 + 42 ) ];
	(*R29).x = bufOutRe[outOffset + ( 0 + me*21 + 14 + 42 ) ];
	(*R31).x = bufOutRe[outOffset + ( 0 + me*21 + 15 + 42 ) ];
	(*R33).x = bufOutRe[outOffset + ( 0 + me*21 + 16 + 42 ) ];
	(*R35).x = bufOutRe[outOffset + ( 0 + me*21 + 17 + 42 ) ];
	(*R37).x = bufOutRe[outOffset + ( 0 + me*21 + 18 + 42 ) ];
	(*R39).x = bufOutRe[outOffset + ( 0 + me*21 + 19 + 42 ) ];
	(*R41).x = bufOutRe[outOffset + ( 0 + me*21 + 20 + 42 ) ];
	}


	__syncthreads();

	if(rw)
	{
	bufOutIm[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 0 ) ] = (*R0).y;
	bufOutIm[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 7 ) ] = (*R1).y;
	bufOutIm[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 14 ) ] = (*R2).y;
	bufOutIm[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 21 ) ] = (*R3).y;
	bufOutIm[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 28 ) ] = (*R4).y;
	bufOutIm[outOffset + ( ((7*me + 0)/7)*42 + (7*me + 0)%7 + 35 ) ] = (*R5).y;
	bufOutIm[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 0 ) ] = (*R6).y;
	bufOutIm[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 7 ) ] = (*R7).y;
	bufOutIm[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 14 ) ] = (*R8).y;
	bufOutIm[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 21 ) ] = (*R9).y;
	bufOutIm[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 28 ) ] = (*R10).y;
	bufOutIm[outOffset + ( ((7*me + 1)/7)*42 + (7*me + 1)%7 + 35 ) ] = (*R11).y;
	bufOutIm[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 0 ) ] = (*R12).y;
	bufOutIm[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 7 ) ] = (*R13).y;
	bufOutIm[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 14 ) ] = (*R14).y;
	bufOutIm[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 21 ) ] = (*R15).y;
	bufOutIm[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 28 ) ] = (*R16).y;
	bufOutIm[outOffset + ( ((7*me + 2)/7)*42 + (7*me + 2)%7 + 35 ) ] = (*R17).y;
	bufOutIm[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 0 ) ] = (*R18).y;
	bufOutIm[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 7 ) ] = (*R19).y;
	bufOutIm[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 14 ) ] = (*R20).y;
	bufOutIm[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 21 ) ] = (*R21).y;
	bufOutIm[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 28 ) ] = (*R22).y;
	bufOutIm[outOffset + ( ((7*me + 3)/7)*42 + (7*me + 3)%7 + 35 ) ] = (*R23).y;
	bufOutIm[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 0 ) ] = (*R24).y;
	bufOutIm[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 7 ) ] = (*R25).y;
	bufOutIm[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 14 ) ] = (*R26).y;
	bufOutIm[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 21 ) ] = (*R27).y;
	bufOutIm[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 28 ) ] = (*R28).y;
	bufOutIm[outOffset + ( ((7*me + 4)/7)*42 + (7*me + 4)%7 + 35 ) ] = (*R29).y;
	bufOutIm[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 0 ) ] = (*R30).y;
	bufOutIm[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 7 ) ] = (*R31).y;
	bufOutIm[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 14 ) ] = (*R32).y;
	bufOutIm[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 21 ) ] = (*R33).y;
	bufOutIm[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 28 ) ] = (*R34).y;
	bufOutIm[outOffset + ( ((7*me + 5)/7)*42 + (7*me + 5)%7 + 35 ) ] = (*R35).y;
	bufOutIm[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 0 ) ] = (*R36).y;
	bufOutIm[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 7 ) ] = (*R37).y;
	bufOutIm[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 14 ) ] = (*R38).y;
	bufOutIm[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 21 ) ] = (*R39).y;
	bufOutIm[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 28 ) ] = (*R40).y;
	bufOutIm[outOffset + ( ((7*me + 6)/7)*42 + (7*me + 6)%7 + 35 ) ] = (*R41).y;
	}


	__syncthreads();

	if(rw)
	{
	(*R0).y = bufOutIm[outOffset + ( 0 + me*21 + 0 + 0 ) ];
	(*R2).y = bufOutIm[outOffset + ( 0 + me*21 + 1 + 0 ) ];
	(*R4).y = bufOutIm[outOffset + ( 0 + me*21 + 2 + 0 ) ];
	(*R6).y = bufOutIm[outOffset + ( 0 + me*21 + 3 + 0 ) ];
	(*R8).y = bufOutIm[outOffset + ( 0 + me*21 + 4 + 0 ) ];
	(*R10).y = bufOutIm[outOffset + ( 0 + me*21 + 5 + 0 ) ];
	(*R12).y = bufOutIm[outOffset + ( 0 + me*21 + 6 + 0 ) ];
	(*R14).y = bufOutIm[outOffset + ( 0 + me*21 + 7 + 0 ) ];
	(*R16).y = bufOutIm[outOffset + ( 0 + me*21 + 8 + 0 ) ];
	(*R18).y = bufOutIm[outOffset + ( 0 + me*21 + 9 + 0 ) ];
	(*R20).y = bufOutIm[outOffset + ( 0 + me*21 + 10 + 0 ) ];
	(*R22).y = bufOutIm[outOffset + ( 0 + me*21 + 11 + 0 ) ];
	(*R24).y = bufOutIm[outOffset + ( 0 + me*21 + 12 + 0 ) ];
	(*R26).y = bufOutIm[outOffset + ( 0 + me*21 + 13 + 0 ) ];
	(*R28).y = bufOutIm[outOffset + ( 0 + me*21 + 14 + 0 ) ];
	(*R30).y = bufOutIm[outOffset + ( 0 + me*21 + 15 + 0 ) ];
	(*R32).y = bufOutIm[outOffset + ( 0 + me*21 + 16 + 0 ) ];
	(*R34).y = bufOutIm[outOffset + ( 0 + me*21 + 17 + 0 ) ];
	(*R36).y = bufOutIm[outOffset + ( 0 + me*21 + 18 + 0 ) ];
	(*R38).y = bufOutIm[outOffset + ( 0 + me*21 + 19 + 0 ) ];
	(*R40).y = bufOutIm[outOffset + ( 0 + me*21 + 20 + 0 ) ];
	(*R1).y = bufOutIm[outOffset + ( 0 + me*21 + 0 + 42 ) ];
	(*R3).y = bufOutIm[outOffset + ( 0 + me*21 + 1 + 42 ) ];
	(*R5).y = bufOutIm[outOffset + ( 0 + me*21 + 2 + 42 ) ];
	(*R7).y = bufOutIm[outOffset + ( 0 + me*21 + 3 + 42 ) ];
	(*R9).y = bufOutIm[outOffset + ( 0 + me*21 + 4 + 42 ) ];
	(*R11).y = bufOutIm[outOffset + ( 0 + me*21 + 5 + 42 ) ];
	(*R13).y = bufOutIm[outOffset + ( 0 + me*21 + 6 + 42 ) ];
	(*R15).y = bufOutIm[outOffset + ( 0 + me*21 + 7 + 42 ) ];
	(*R17).y = bufOutIm[outOffset + ( 0 + me*21 + 8 + 42 ) ];
	(*R19).y = bufOutIm[outOffset + ( 0 + me*21 + 9 + 42 ) ];
	(*R21).y = bufOutIm[outOffset + ( 0 + me*21 + 10 + 42 ) ];
	(*R23).y = bufOutIm[outOffset + ( 0 + me*21 + 11 + 42 ) ];
	(*R25).y = bufOutIm[outOffset + ( 0 + me*21 + 12 + 42 ) ];
	(*R27).y = bufOutIm[outOffset + ( 0 + me*21 + 13 + 42 ) ];
	(*R29).y = bufOutIm[outOffset + ( 0 + me*21 + 14 + 42 ) ];
	(*R31).y = bufOutIm[outOffset + ( 0 + me*21 + 15 + 42 ) ];
	(*R33).y = bufOutIm[outOffset + ( 0 + me*21 + 16 + 42 ) ];
	(*R35).y = bufOutIm[outOffset + ( 0 + me*21 + 17 + 42 ) ];
	(*R37).y = bufOutIm[outOffset + ( 0 + me*21 + 18 + 42 ) ];
	(*R39).y = bufOutIm[outOffset + ( 0 + me*21 + 19 + 42 ) ];
	(*R41).y = bufOutIm[outOffset + ( 0 + me*21 + 20 + 42 ) ];
	}


	__syncthreads();

}

template <typename T, StrideBin sb>
__device__ void
InvPass2_len84(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, T *R30, T *R31, T *R32, T *R33, T *R34, T *R35, T *R36, T *R37, T *R38, T *R39, T *R40, T *R41)
{




	{
		T W = twiddles[41 + 1*((21*me + 0)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 1)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 2)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 3)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 4)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 5)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 6)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 7)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 8)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 9)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 10)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 11)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 12)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 13)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 14)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 15)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R31).x; ry = (*R31).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R31).x = TR;
		(*R31).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 16)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R33).x; ry = (*R33).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R33).x = TR;
		(*R33).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 17)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R35).x; ry = (*R35).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R35).x = TR;
		(*R35).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 18)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R37).x; ry = (*R37).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R37).x = TR;
		(*R37).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 19)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R39).x; ry = (*R39).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R39).x = TR;
		(*R39).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 20)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R41).x; ry = (*R41).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R41).x = TR;
		(*R41).y = TI;
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
	InvRad2B1(R26, R27);
	InvRad2B1(R28, R29);
	InvRad2B1(R30, R31);
	InvRad2B1(R32, R33);
	InvRad2B1(R34, R35);
	InvRad2B1(R36, R37);
	InvRad2B1(R38, R39);
	InvRad2B1(R40, R41);


	if(rw)
	{
	bufOut[outOffset + ( 21*me + 0 + 0 )*stride_out] = (*R0);
	bufOut[outOffset + ( 21*me + 1 + 0 )*stride_out] = (*R2);
	bufOut[outOffset + ( 21*me + 2 + 0 )*stride_out] = (*R4);
	bufOut[outOffset + ( 21*me + 3 + 0 )*stride_out] = (*R6);
	bufOut[outOffset + ( 21*me + 4 + 0 )*stride_out] = (*R8);
	bufOut[outOffset + ( 21*me + 5 + 0 )*stride_out] = (*R10);
	bufOut[outOffset + ( 21*me + 6 + 0 )*stride_out] = (*R12);
	bufOut[outOffset + ( 21*me + 7 + 0 )*stride_out] = (*R14);
	bufOut[outOffset + ( 21*me + 8 + 0 )*stride_out] = (*R16);
	bufOut[outOffset + ( 21*me + 9 + 0 )*stride_out] = (*R18);
	bufOut[outOffset + ( 21*me + 10 + 0 )*stride_out] = (*R20);
	bufOut[outOffset + ( 21*me + 11 + 0 )*stride_out] = (*R22);
	bufOut[outOffset + ( 21*me + 12 + 0 )*stride_out] = (*R24);
	bufOut[outOffset + ( 21*me + 13 + 0 )*stride_out] = (*R26);
	bufOut[outOffset + ( 21*me + 14 + 0 )*stride_out] = (*R28);
	bufOut[outOffset + ( 21*me + 15 + 0 )*stride_out] = (*R30);
	bufOut[outOffset + ( 21*me + 16 + 0 )*stride_out] = (*R32);
	bufOut[outOffset + ( 21*me + 17 + 0 )*stride_out] = (*R34);
	bufOut[outOffset + ( 21*me + 18 + 0 )*stride_out] = (*R36);
	bufOut[outOffset + ( 21*me + 19 + 0 )*stride_out] = (*R38);
	bufOut[outOffset + ( 21*me + 20 + 0 )*stride_out] = (*R40);
	bufOut[outOffset + ( 21*me + 0 + 42 )*stride_out] = (*R1);
	bufOut[outOffset + ( 21*me + 1 + 42 )*stride_out] = (*R3);
	bufOut[outOffset + ( 21*me + 2 + 42 )*stride_out] = (*R5);
	bufOut[outOffset + ( 21*me + 3 + 42 )*stride_out] = (*R7);
	bufOut[outOffset + ( 21*me + 4 + 42 )*stride_out] = (*R9);
	bufOut[outOffset + ( 21*me + 5 + 42 )*stride_out] = (*R11);
	bufOut[outOffset + ( 21*me + 6 + 42 )*stride_out] = (*R13);
	bufOut[outOffset + ( 21*me + 7 + 42 )*stride_out] = (*R15);
	bufOut[outOffset + ( 21*me + 8 + 42 )*stride_out] = (*R17);
	bufOut[outOffset + ( 21*me + 9 + 42 )*stride_out] = (*R19);
	bufOut[outOffset + ( 21*me + 10 + 42 )*stride_out] = (*R21);
	bufOut[outOffset + ( 21*me + 11 + 42 )*stride_out] = (*R23);
	bufOut[outOffset + ( 21*me + 12 + 42 )*stride_out] = (*R25);
	bufOut[outOffset + ( 21*me + 13 + 42 )*stride_out] = (*R27);
	bufOut[outOffset + ( 21*me + 14 + 42 )*stride_out] = (*R29);
	bufOut[outOffset + ( 21*me + 15 + 42 )*stride_out] = (*R31);
	bufOut[outOffset + ( 21*me + 16 + 42 )*stride_out] = (*R33);
	bufOut[outOffset + ( 21*me + 17 + 42 )*stride_out] = (*R35);
	bufOut[outOffset + ( 21*me + 18 + 42 )*stride_out] = (*R37);
	bufOut[outOffset + ( 21*me + 19 + 42 )*stride_out] = (*R39);
	bufOut[outOffset + ( 21*me + 20 + 42 )*stride_out] = (*R41);
	}

}

template <typename T, StrideBin sb>
__device__ void
InvPass2_len84(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19, T *R20, T *R21, T *R22, T *R23, T *R24, T *R25, T *R26, T *R27, T *R28, T *R29, T *R30, T *R31, T *R32, T *R33, T *R34, T *R35, T *R36, T *R37, T *R38, T *R39, T *R40, T *R41)
{




	{
		T W = twiddles[41 + 1*((21*me + 0)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 1)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 2)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 3)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 4)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 5)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 6)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 7)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 8)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 9)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 10)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 11)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 12)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 13)%42) + 0];
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
		T W = twiddles[41 + 1*((21*me + 14)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R29).x; ry = (*R29).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R29).x = TR;
		(*R29).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 15)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R31).x; ry = (*R31).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R31).x = TR;
		(*R31).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 16)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R33).x; ry = (*R33).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R33).x = TR;
		(*R33).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 17)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R35).x; ry = (*R35).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R35).x = TR;
		(*R35).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 18)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R37).x; ry = (*R37).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R37).x = TR;
		(*R37).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 19)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R39).x; ry = (*R39).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R39).x = TR;
		(*R39).y = TI;
	}

	{
		T W = twiddles[41 + 1*((21*me + 20)%42) + 0];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R41).x; ry = (*R41).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R41).x = TR;
		(*R41).y = TI;
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
	InvRad2B1(R26, R27);
	InvRad2B1(R28, R29);
	InvRad2B1(R30, R31);
	InvRad2B1(R32, R33);
	InvRad2B1(R34, R35);
	InvRad2B1(R36, R37);
	InvRad2B1(R38, R39);
	InvRad2B1(R40, R41);


	if(rw)
	{
	bufOutRe[outOffset + ( 21*me + 0 + 0 )*stride_out] = (*R0).x;
	bufOutIm[outOffset + ( 21*me + 0 + 0 )*stride_out] = (*R0).y;
	bufOutRe[outOffset + ( 21*me + 1 + 0 )*stride_out] = (*R2).x;
	bufOutIm[outOffset + ( 21*me + 1 + 0 )*stride_out] = (*R2).y;
	bufOutRe[outOffset + ( 21*me + 2 + 0 )*stride_out] = (*R4).x;
	bufOutIm[outOffset + ( 21*me + 2 + 0 )*stride_out] = (*R4).y;
	bufOutRe[outOffset + ( 21*me + 3 + 0 )*stride_out] = (*R6).x;
	bufOutIm[outOffset + ( 21*me + 3 + 0 )*stride_out] = (*R6).y;
	bufOutRe[outOffset + ( 21*me + 4 + 0 )*stride_out] = (*R8).x;
	bufOutIm[outOffset + ( 21*me + 4 + 0 )*stride_out] = (*R8).y;
	bufOutRe[outOffset + ( 21*me + 5 + 0 )*stride_out] = (*R10).x;
	bufOutIm[outOffset + ( 21*me + 5 + 0 )*stride_out] = (*R10).y;
	bufOutRe[outOffset + ( 21*me + 6 + 0 )*stride_out] = (*R12).x;
	bufOutIm[outOffset + ( 21*me + 6 + 0 )*stride_out] = (*R12).y;
	bufOutRe[outOffset + ( 21*me + 7 + 0 )*stride_out] = (*R14).x;
	bufOutIm[outOffset + ( 21*me + 7 + 0 )*stride_out] = (*R14).y;
	bufOutRe[outOffset + ( 21*me + 8 + 0 )*stride_out] = (*R16).x;
	bufOutIm[outOffset + ( 21*me + 8 + 0 )*stride_out] = (*R16).y;
	bufOutRe[outOffset + ( 21*me + 9 + 0 )*stride_out] = (*R18).x;
	bufOutIm[outOffset + ( 21*me + 9 + 0 )*stride_out] = (*R18).y;
	bufOutRe[outOffset + ( 21*me + 10 + 0 )*stride_out] = (*R20).x;
	bufOutIm[outOffset + ( 21*me + 10 + 0 )*stride_out] = (*R20).y;
	bufOutRe[outOffset + ( 21*me + 11 + 0 )*stride_out] = (*R22).x;
	bufOutIm[outOffset + ( 21*me + 11 + 0 )*stride_out] = (*R22).y;
	bufOutRe[outOffset + ( 21*me + 12 + 0 )*stride_out] = (*R24).x;
	bufOutIm[outOffset + ( 21*me + 12 + 0 )*stride_out] = (*R24).y;
	bufOutRe[outOffset + ( 21*me + 13 + 0 )*stride_out] = (*R26).x;
	bufOutIm[outOffset + ( 21*me + 13 + 0 )*stride_out] = (*R26).y;
	bufOutRe[outOffset + ( 21*me + 14 + 0 )*stride_out] = (*R28).x;
	bufOutIm[outOffset + ( 21*me + 14 + 0 )*stride_out] = (*R28).y;
	bufOutRe[outOffset + ( 21*me + 15 + 0 )*stride_out] = (*R30).x;
	bufOutIm[outOffset + ( 21*me + 15 + 0 )*stride_out] = (*R30).y;
	bufOutRe[outOffset + ( 21*me + 16 + 0 )*stride_out] = (*R32).x;
	bufOutIm[outOffset + ( 21*me + 16 + 0 )*stride_out] = (*R32).y;
	bufOutRe[outOffset + ( 21*me + 17 + 0 )*stride_out] = (*R34).x;
	bufOutIm[outOffset + ( 21*me + 17 + 0 )*stride_out] = (*R34).y;
	bufOutRe[outOffset + ( 21*me + 18 + 0 )*stride_out] = (*R36).x;
	bufOutIm[outOffset + ( 21*me + 18 + 0 )*stride_out] = (*R36).y;
	bufOutRe[outOffset + ( 21*me + 19 + 0 )*stride_out] = (*R38).x;
	bufOutIm[outOffset + ( 21*me + 19 + 0 )*stride_out] = (*R38).y;
	bufOutRe[outOffset + ( 21*me + 20 + 0 )*stride_out] = (*R40).x;
	bufOutIm[outOffset + ( 21*me + 20 + 0 )*stride_out] = (*R40).y;
	bufOutRe[outOffset + ( 21*me + 0 + 42 )*stride_out] = (*R1).x;
	bufOutIm[outOffset + ( 21*me + 0 + 42 )*stride_out] = (*R1).y;
	bufOutRe[outOffset + ( 21*me + 1 + 42 )*stride_out] = (*R3).x;
	bufOutIm[outOffset + ( 21*me + 1 + 42 )*stride_out] = (*R3).y;
	bufOutRe[outOffset + ( 21*me + 2 + 42 )*stride_out] = (*R5).x;
	bufOutIm[outOffset + ( 21*me + 2 + 42 )*stride_out] = (*R5).y;
	bufOutRe[outOffset + ( 21*me + 3 + 42 )*stride_out] = (*R7).x;
	bufOutIm[outOffset + ( 21*me + 3 + 42 )*stride_out] = (*R7).y;
	bufOutRe[outOffset + ( 21*me + 4 + 42 )*stride_out] = (*R9).x;
	bufOutIm[outOffset + ( 21*me + 4 + 42 )*stride_out] = (*R9).y;
	bufOutRe[outOffset + ( 21*me + 5 + 42 )*stride_out] = (*R11).x;
	bufOutIm[outOffset + ( 21*me + 5 + 42 )*stride_out] = (*R11).y;
	bufOutRe[outOffset + ( 21*me + 6 + 42 )*stride_out] = (*R13).x;
	bufOutIm[outOffset + ( 21*me + 6 + 42 )*stride_out] = (*R13).y;
	bufOutRe[outOffset + ( 21*me + 7 + 42 )*stride_out] = (*R15).x;
	bufOutIm[outOffset + ( 21*me + 7 + 42 )*stride_out] = (*R15).y;
	bufOutRe[outOffset + ( 21*me + 8 + 42 )*stride_out] = (*R17).x;
	bufOutIm[outOffset + ( 21*me + 8 + 42 )*stride_out] = (*R17).y;
	bufOutRe[outOffset + ( 21*me + 9 + 42 )*stride_out] = (*R19).x;
	bufOutIm[outOffset + ( 21*me + 9 + 42 )*stride_out] = (*R19).y;
	bufOutRe[outOffset + ( 21*me + 10 + 42 )*stride_out] = (*R21).x;
	bufOutIm[outOffset + ( 21*me + 10 + 42 )*stride_out] = (*R21).y;
	bufOutRe[outOffset + ( 21*me + 11 + 42 )*stride_out] = (*R23).x;
	bufOutIm[outOffset + ( 21*me + 11 + 42 )*stride_out] = (*R23).y;
	bufOutRe[outOffset + ( 21*me + 12 + 42 )*stride_out] = (*R25).x;
	bufOutIm[outOffset + ( 21*me + 12 + 42 )*stride_out] = (*R25).y;
	bufOutRe[outOffset + ( 21*me + 13 + 42 )*stride_out] = (*R27).x;
	bufOutIm[outOffset + ( 21*me + 13 + 42 )*stride_out] = (*R27).y;
	bufOutRe[outOffset + ( 21*me + 14 + 42 )*stride_out] = (*R29).x;
	bufOutIm[outOffset + ( 21*me + 14 + 42 )*stride_out] = (*R29).y;
	bufOutRe[outOffset + ( 21*me + 15 + 42 )*stride_out] = (*R31).x;
	bufOutIm[outOffset + ( 21*me + 15 + 42 )*stride_out] = (*R31).y;
	bufOutRe[outOffset + ( 21*me + 16 + 42 )*stride_out] = (*R33).x;
	bufOutIm[outOffset + ( 21*me + 16 + 42 )*stride_out] = (*R33).y;
	bufOutRe[outOffset + ( 21*me + 17 + 42 )*stride_out] = (*R35).x;
	bufOutIm[outOffset + ( 21*me + 17 + 42 )*stride_out] = (*R35).y;
	bufOutRe[outOffset + ( 21*me + 18 + 42 )*stride_out] = (*R37).x;
	bufOutIm[outOffset + ( 21*me + 18 + 42 )*stride_out] = (*R37).y;
	bufOutRe[outOffset + ( 21*me + 19 + 42 )*stride_out] = (*R39).x;
	bufOutIm[outOffset + ( 21*me + 19 + 42 )*stride_out] = (*R39).y;
	bufOutRe[outOffset + ( 21*me + 20 + 42 )*stride_out] = (*R41).x;
	bufOutIm[outOffset + ( 21*me + 20 + 42 )*stride_out] = (*R41).y;
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb>
__device__ void 
fwd_len84_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29, R30, R31, R32, R33, R34, R35, R36, R37, R38, R39, R40, R41;
	FwdPass0_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	FwdPass1_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	FwdPass2_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
}

template <typename T, StrideBin sb>
__device__ void 
back_len84_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29, R30, R31, R32, R33, R34, R35, R36, R37, R38, R39, R40, R41;
	InvPass0_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	InvPass1_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	InvPass2_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len84_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29, R30, R31, R32, R33, R34, R35, R36, R37, R38, R39, R40, R41;
	FwdPass0_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	FwdPass1_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	FwdPass2_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
}

template <typename T, StrideBin sb>
__device__ void 
back_len84_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29, R30, R31, R32, R33, R34, R35, R36, R37, R38, R39, R40, R41;
	InvPass0_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  lwbIn, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	InvPass1_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	InvPass2_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len84_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29, R30, R31, R32, R33, R34, R35, R36, R37, R38, R39, R40, R41;
	FwdPass0_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	FwdPass1_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	FwdPass2_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
}

template <typename T, StrideBin sb>
__device__ void 
back_len84_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29, R30, R31, R32, R33, R34, R35, R36, R37, R38, R39, R40, R41;
	InvPass0_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	InvPass1_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	InvPass2_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
}

template <typename T, StrideBin sb>
__device__ void 
fwd_len84_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29, R30, R31, R32, R33, R34, R35, R36, R37, R38, R39, R40, R41;
	FwdPass0_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	FwdPass1_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	FwdPass2_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
}

template <typename T, StrideBin sb>
__device__ void 
back_len84_device(const T *twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, real_type_t<T> *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29, R30, R31, R32, R33, R34, R35, R36, R37, R38, R39, R40, R41;
	InvPass0_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, 0, ldsOffset,  bufInRe, bufInIm, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	InvPass1_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
	InvPass2_len84<T, sb>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, 0, lds, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19, &R20, &R21, &R22, &R23, &R24, &R25, &R26, &R27, &R28, &R29, &R30, &R31, &R32, &R33, &R34, &R35, &R36, &R37, &R38, &R39, &R40, &R41);
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 60, maximum transforms: 30, Passes: 3
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(60)
fft_fwd_ip_len84( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2520), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len84>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*30)*2) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*30 + (me/2));
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
	fwd_len84_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%2, (me/2)*84, lwb, lwb, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 30, Passes: 3
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(60)
fft_back_ip_len84( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2520), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len84>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*30)*2) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*30 + (me/2));
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
	back_len84_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%2, (me/2)*84, lwb, lwb, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 30, Passes: 3
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(60)
fft_fwd_ip_len84( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2520), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len84>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*30)*2) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*30 + (me/2));
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
	fwd_len84_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%2, (me/2)*84, lwbRe, lwbIm, lwbRe, lwbIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 30, Passes: 3
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(60)
fft_back_ip_len84( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2520), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len84>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*30)*2) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*30 + (me/2));
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
	back_len84_device<T, sb>(twiddles, stride_in[0], stride_in[0],  rw, b, me%2, (me/2)*84, lwbRe, lwbIm, lwbRe, lwbIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 30, Passes: 3
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(60)
fft_fwd_op_len84( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2520), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len84>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*30)*2) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*30 + (me/2));
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
	fwd_len84_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*84, lwbIn, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 30, Passes: 3
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(60)
fft_back_op_len84( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2520), cgh);
	cgh.parallel_for<class kern_fft_back_op_len84>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*30)*2) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*30 + (me/2));
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
	back_len84_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*84, lwbIn, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 30, Passes: 3
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(60)
fft_fwd_op_len84( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2520), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len84>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*30)*2) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*30 + (me/2));
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
	fwd_len84_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*84, lwbIn, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 30, Passes: 3
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(60)
fft_back_op_len84( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2520), cgh);
	cgh.parallel_for<class kern_fft_back_op_len84>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*30)*2) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*30 + (me/2));
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
	back_len84_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*84, lwbIn, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 30, Passes: 3
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(60)
fft_fwd_op_len84( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2520), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len84>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*30)*2) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*30 + (me/2));
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
	fwd_len84_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*84, lwbInRe, lwbInIm, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 30, Passes: 3
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(60)
fft_back_op_len84( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2520), cgh);
	cgh.parallel_for<class kern_fft_back_op_len84>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*30)*2) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*30 + (me/2));
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
	back_len84_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*84, lwbInRe, lwbInIm, lwbOut, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 30, Passes: 3
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(60)
fft_fwd_op_len84( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2520), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len84>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*30)*2) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*30 + (me/2));
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
	fwd_len84_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*84, lwbInRe, lwbInIm, lwbOutRe, lwbOutIm, lds);
});
});
}

//Kernel configuration: number of threads per thread block: 60, maximum transforms: 30, Passes: 3
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(60)
fft_back_op_len84( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
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
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2520), cgh);
	cgh.parallel_for<class kern_fft_back_op_len84>(
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
	unsigned int rw = (static_cast<int>(me) < (static_cast<int>(upper_count)  - static_cast<int>(batch)*30)*2) ? 1 : 0;

	//suppress warning
	#ifdef __NVCC__
		(void)(rw == rw);
	#else
		(void)rw;
	#endif
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	size_t counter_mod = (batch*30 + (me/2));
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
	back_len84_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*84, lwbInRe, lwbInIm, lwbOutRe, lwbOutIm, lds);
});
});
}

