#pragma once
#include "rocfft_butterfly_template.h"

////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb, bool TwdLarge>
__device__ void
FwdPass0_len200_sbcc(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R10) = bufIn[inOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*2 + 0 + 20 ) ];
	(*R11) = bufIn[inOffset + ( 0 + me*2 + 1 + 20 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*2 + 0 + 40 ) ];
	(*R12) = bufIn[inOffset + ( 0 + me*2 + 1 + 40 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*2 + 0 + 60 ) ];
	(*R13) = bufIn[inOffset + ( 0 + me*2 + 1 + 60 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*2 + 0 + 80 ) ];
	(*R14) = bufIn[inOffset + ( 0 + me*2 + 1 + 80 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*2 + 0 + 100 ) ];
	(*R15) = bufIn[inOffset + ( 0 + me*2 + 1 + 100 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*2 + 0 + 120 ) ];
	(*R16) = bufIn[inOffset + ( 0 + me*2 + 1 + 120 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*2 + 0 + 140 ) ];
	(*R17) = bufIn[inOffset + ( 0 + me*2 + 1 + 140 ) ];
	(*R8) = bufIn[inOffset + ( 0 + me*2 + 0 + 160 ) ];
	(*R18) = bufIn[inOffset + ( 0 + me*2 + 1 + 160 ) ];
	(*R9) = bufIn[inOffset + ( 0 + me*2 + 0 + 180 ) ];
	(*R19) = bufIn[inOffset + ( 0 + me*2 + 1 + 180 ) ];
	}



	FwdRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	FwdRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 1 ) ] = (*R1);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 2 ) ] = (*R2);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 3 ) ] = (*R3);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 4 ) ] = (*R4);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 5 ) ] = (*R5);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 6 ) ] = (*R6);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 7 ) ] = (*R7);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 8 ) ] = (*R8);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 9 ) ] = (*R9);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 0 ) ] = (*R10);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 1 ) ] = (*R11);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 2 ) ] = (*R12);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 3 ) ] = (*R13);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 4 ) ] = (*R14);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 5 ) ] = (*R15);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 6 ) ] = (*R16);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 7 ) ] = (*R17);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 8 ) ] = (*R18);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 9 ) ] = (*R19);
	}

}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void
FwdPass1_len200_sbcc(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R10) = bufIn[inOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*2 + 0 + 20 ) ];
	(*R11) = bufIn[inOffset + ( 0 + me*2 + 1 + 20 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*2 + 0 + 40 ) ];
	(*R12) = bufIn[inOffset + ( 0 + me*2 + 1 + 40 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*2 + 0 + 60 ) ];
	(*R13) = bufIn[inOffset + ( 0 + me*2 + 1 + 60 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*2 + 0 + 80 ) ];
	(*R14) = bufIn[inOffset + ( 0 + me*2 + 1 + 80 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*2 + 0 + 100 ) ];
	(*R15) = bufIn[inOffset + ( 0 + me*2 + 1 + 100 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*2 + 0 + 120 ) ];
	(*R16) = bufIn[inOffset + ( 0 + me*2 + 1 + 120 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*2 + 0 + 140 ) ];
	(*R17) = bufIn[inOffset + ( 0 + me*2 + 1 + 140 ) ];
	(*R8) = bufIn[inOffset + ( 0 + me*2 + 0 + 160 ) ];
	(*R18) = bufIn[inOffset + ( 0 + me*2 + 1 + 160 ) ];
	(*R9) = bufIn[inOffset + ( 0 + me*2 + 0 + 180 ) ];
	(*R19) = bufIn[inOffset + ( 0 + me*2 + 1 + 180 ) ];
	}



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

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 10 ) ] = (*R1);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 20 ) ] = (*R2);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 30 ) ] = (*R3);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 40 ) ] = (*R4);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 50 ) ] = (*R5);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 60 ) ] = (*R6);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 70 ) ] = (*R7);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 80 ) ] = (*R8);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 90 ) ] = (*R9);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 0 ) ] = (*R10);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 10 ) ] = (*R11);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 20 ) ] = (*R12);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 30 ) ] = (*R13);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 40 ) ] = (*R14);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 50 ) ] = (*R15);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 60 ) ] = (*R16);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 70 ) ] = (*R17);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 80 ) ] = (*R18);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 90 ) ] = (*R19);
	}

}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void
FwdPass2_len200_sbcc(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*10 + 0 + 0 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*10 + 1 + 0 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*10 + 2 + 0 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*10 + 3 + 0 ) ];
	(*R8) = bufIn[inOffset + ( 0 + me*10 + 4 + 0 ) ];
	(*R10) = bufIn[inOffset + ( 0 + me*10 + 5 + 0 ) ];
	(*R12) = bufIn[inOffset + ( 0 + me*10 + 6 + 0 ) ];
	(*R14) = bufIn[inOffset + ( 0 + me*10 + 7 + 0 ) ];
	(*R16) = bufIn[inOffset + ( 0 + me*10 + 8 + 0 ) ];
	(*R18) = bufIn[inOffset + ( 0 + me*10 + 9 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*10 + 0 + 100 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*10 + 1 + 100 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*10 + 2 + 100 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*10 + 3 + 100 ) ];
	(*R9) = bufIn[inOffset + ( 0 + me*10 + 4 + 100 ) ];
	(*R11) = bufIn[inOffset + ( 0 + me*10 + 5 + 100 ) ];
	(*R13) = bufIn[inOffset + ( 0 + me*10 + 6 + 100 ) ];
	(*R15) = bufIn[inOffset + ( 0 + me*10 + 7 + 100 ) ];
	(*R17) = bufIn[inOffset + ( 0 + me*10 + 8 + 100 ) ];
	(*R19) = bufIn[inOffset + ( 0 + me*10 + 9 + 100 ) ];
	}



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

	__syncthreads();



	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 0)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R0).x; ry = (*R0).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R0).x = TR;
		(*R0).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 0)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R1).x; ry = (*R1).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R1).x = TR;
		(*R1).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 1)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R2).x; ry = (*R2).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R2).x = TR;
		(*R2).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 1)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 2)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R4).x; ry = (*R4).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R4).x = TR;
		(*R4).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 2)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R5).x; ry = (*R5).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R5).x = TR;
		(*R5).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 3)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R6).x; ry = (*R6).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R6).x = TR;
		(*R6).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 3)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R7).x; ry = (*R7).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R7).x = TR;
		(*R7).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 4)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 4)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R9).x; ry = (*R9).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R9).x = TR;
		(*R9).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 5)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R10).x; ry = (*R10).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R10).x = TR;
		(*R10).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 5)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R11).x; ry = (*R11).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R11).x = TR;
		(*R11).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 6)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R12).x; ry = (*R12).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R12).x = TR;
		(*R12).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 6)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 7)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 7)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R15).x; ry = (*R15).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R15).x = TR;
		(*R15).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 8)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R16).x; ry = (*R16).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R16).x = TR;
		(*R16).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 8)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R17).x; ry = (*R17).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R17).x = TR;
		(*R17).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 9)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R18).x; ry = (*R18).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R18).x = TR;
		(*R18).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 9)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	if(rw)
	{
	bufOut[outOffset + ( 10*me + 0 + 0 ) ] = (*R0);
	bufOut[outOffset + ( 10*me + 1 + 0 ) ] = (*R2);
	bufOut[outOffset + ( 10*me + 2 + 0 ) ] = (*R4);
	bufOut[outOffset + ( 10*me + 3 + 0 ) ] = (*R6);
	bufOut[outOffset + ( 10*me + 4 + 0 ) ] = (*R8);
	bufOut[outOffset + ( 10*me + 5 + 0 ) ] = (*R10);
	bufOut[outOffset + ( 10*me + 6 + 0 ) ] = (*R12);
	bufOut[outOffset + ( 10*me + 7 + 0 ) ] = (*R14);
	bufOut[outOffset + ( 10*me + 8 + 0 ) ] = (*R16);
	bufOut[outOffset + ( 10*me + 9 + 0 ) ] = (*R18);
	bufOut[outOffset + ( 10*me + 0 + 100 ) ] = (*R1);
	bufOut[outOffset + ( 10*me + 1 + 100 ) ] = (*R3);
	bufOut[outOffset + ( 10*me + 2 + 100 ) ] = (*R5);
	bufOut[outOffset + ( 10*me + 3 + 100 ) ] = (*R7);
	bufOut[outOffset + ( 10*me + 4 + 100 ) ] = (*R9);
	bufOut[outOffset + ( 10*me + 5 + 100 ) ] = (*R11);
	bufOut[outOffset + ( 10*me + 6 + 100 ) ] = (*R13);
	bufOut[outOffset + ( 10*me + 7 + 100 ) ] = (*R15);
	bufOut[outOffset + ( 10*me + 8 + 100 ) ] = (*R17);
	bufOut[outOffset + ( 10*me + 9 + 100 ) ] = (*R19);
	}

}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void
InvPass0_len200_sbcc(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R10) = bufIn[inOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*2 + 0 + 20 ) ];
	(*R11) = bufIn[inOffset + ( 0 + me*2 + 1 + 20 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*2 + 0 + 40 ) ];
	(*R12) = bufIn[inOffset + ( 0 + me*2 + 1 + 40 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*2 + 0 + 60 ) ];
	(*R13) = bufIn[inOffset + ( 0 + me*2 + 1 + 60 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*2 + 0 + 80 ) ];
	(*R14) = bufIn[inOffset + ( 0 + me*2 + 1 + 80 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*2 + 0 + 100 ) ];
	(*R15) = bufIn[inOffset + ( 0 + me*2 + 1 + 100 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*2 + 0 + 120 ) ];
	(*R16) = bufIn[inOffset + ( 0 + me*2 + 1 + 120 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*2 + 0 + 140 ) ];
	(*R17) = bufIn[inOffset + ( 0 + me*2 + 1 + 140 ) ];
	(*R8) = bufIn[inOffset + ( 0 + me*2 + 0 + 160 ) ];
	(*R18) = bufIn[inOffset + ( 0 + me*2 + 1 + 160 ) ];
	(*R9) = bufIn[inOffset + ( 0 + me*2 + 0 + 180 ) ];
	(*R19) = bufIn[inOffset + ( 0 + me*2 + 1 + 180 ) ];
	}



	InvRad10B1(R0, R1, R2, R3, R4, R5, R6, R7, R8, R9);
	InvRad10B1(R10, R11, R12, R13, R14, R15, R16, R17, R18, R19);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 1 ) ] = (*R1);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 2 ) ] = (*R2);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 3 ) ] = (*R3);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 4 ) ] = (*R4);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 5 ) ] = (*R5);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 6 ) ] = (*R6);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 7 ) ] = (*R7);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 8 ) ] = (*R8);
	bufOut[outOffset + ( ((2*me + 0)/1)*10 + (2*me + 0)%1 + 9 ) ] = (*R9);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 0 ) ] = (*R10);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 1 ) ] = (*R11);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 2 ) ] = (*R12);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 3 ) ] = (*R13);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 4 ) ] = (*R14);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 5 ) ] = (*R15);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 6 ) ] = (*R16);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 7 ) ] = (*R17);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 8 ) ] = (*R18);
	bufOut[outOffset + ( ((2*me + 1)/1)*10 + (2*me + 1)%1 + 9 ) ] = (*R19);
	}

}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void
InvPass1_len200_sbcc(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*2 + 0 + 0 ) ];
	(*R10) = bufIn[inOffset + ( 0 + me*2 + 1 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*2 + 0 + 20 ) ];
	(*R11) = bufIn[inOffset + ( 0 + me*2 + 1 + 20 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*2 + 0 + 40 ) ];
	(*R12) = bufIn[inOffset + ( 0 + me*2 + 1 + 40 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*2 + 0 + 60 ) ];
	(*R13) = bufIn[inOffset + ( 0 + me*2 + 1 + 60 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*2 + 0 + 80 ) ];
	(*R14) = bufIn[inOffset + ( 0 + me*2 + 1 + 80 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*2 + 0 + 100 ) ];
	(*R15) = bufIn[inOffset + ( 0 + me*2 + 1 + 100 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*2 + 0 + 120 ) ];
	(*R16) = bufIn[inOffset + ( 0 + me*2 + 1 + 120 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*2 + 0 + 140 ) ];
	(*R17) = bufIn[inOffset + ( 0 + me*2 + 1 + 140 ) ];
	(*R8) = bufIn[inOffset + ( 0 + me*2 + 0 + 160 ) ];
	(*R18) = bufIn[inOffset + ( 0 + me*2 + 1 + 160 ) ];
	(*R9) = bufIn[inOffset + ( 0 + me*2 + 0 + 180 ) ];
	(*R19) = bufIn[inOffset + ( 0 + me*2 + 1 + 180 ) ];
	}



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

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 10 ) ] = (*R1);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 20 ) ] = (*R2);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 30 ) ] = (*R3);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 40 ) ] = (*R4);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 50 ) ] = (*R5);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 60 ) ] = (*R6);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 70 ) ] = (*R7);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 80 ) ] = (*R8);
	bufOut[outOffset + ( ((2*me + 0)/10)*100 + (2*me + 0)%10 + 90 ) ] = (*R9);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 0 ) ] = (*R10);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 10 ) ] = (*R11);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 20 ) ] = (*R12);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 30 ) ] = (*R13);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 40 ) ] = (*R14);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 50 ) ] = (*R15);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 60 ) ] = (*R16);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 70 ) ] = (*R17);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 80 ) ] = (*R18);
	bufOut[outOffset + ( ((2*me + 1)/10)*100 + (2*me + 1)%10 + 90 ) ] = (*R19);
	}

}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void
InvPass2_len200_sbcc(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*10 + 0 + 0 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*10 + 1 + 0 ) ];
	(*R4) = bufIn[inOffset + ( 0 + me*10 + 2 + 0 ) ];
	(*R6) = bufIn[inOffset + ( 0 + me*10 + 3 + 0 ) ];
	(*R8) = bufIn[inOffset + ( 0 + me*10 + 4 + 0 ) ];
	(*R10) = bufIn[inOffset + ( 0 + me*10 + 5 + 0 ) ];
	(*R12) = bufIn[inOffset + ( 0 + me*10 + 6 + 0 ) ];
	(*R14) = bufIn[inOffset + ( 0 + me*10 + 7 + 0 ) ];
	(*R16) = bufIn[inOffset + ( 0 + me*10 + 8 + 0 ) ];
	(*R18) = bufIn[inOffset + ( 0 + me*10 + 9 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*10 + 0 + 100 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*10 + 1 + 100 ) ];
	(*R5) = bufIn[inOffset + ( 0 + me*10 + 2 + 100 ) ];
	(*R7) = bufIn[inOffset + ( 0 + me*10 + 3 + 100 ) ];
	(*R9) = bufIn[inOffset + ( 0 + me*10 + 4 + 100 ) ];
	(*R11) = bufIn[inOffset + ( 0 + me*10 + 5 + 100 ) ];
	(*R13) = bufIn[inOffset + ( 0 + me*10 + 6 + 100 ) ];
	(*R15) = bufIn[inOffset + ( 0 + me*10 + 7 + 100 ) ];
	(*R17) = bufIn[inOffset + ( 0 + me*10 + 8 + 100 ) ];
	(*R19) = bufIn[inOffset + ( 0 + me*10 + 9 + 100 ) ];
	}



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

	__syncthreads();



	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 0)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R0).x; ry = (*R0).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R0).x = TR;
		(*R0).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 0)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R1).x; ry = (*R1).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R1).x = TR;
		(*R1).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 1)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R2).x; ry = (*R2).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R2).x = TR;
		(*R2).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 1)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 2)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R4).x; ry = (*R4).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R4).x = TR;
		(*R4).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 2)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R5).x; ry = (*R5).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R5).x = TR;
		(*R5).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 3)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R6).x; ry = (*R6).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R6).x = TR;
		(*R6).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 3)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R7).x; ry = (*R7).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R7).x = TR;
		(*R7).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 4)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R8).x; ry = (*R8).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R8).x = TR;
		(*R8).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 4)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R9).x; ry = (*R9).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R9).x = TR;
		(*R9).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 5)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R10).x; ry = (*R10).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R10).x = TR;
		(*R10).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 5)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R11).x; ry = (*R11).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R11).x = TR;
		(*R11).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 6)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R12).x; ry = (*R12).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R12).x = TR;
		(*R12).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 6)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R13).x; ry = (*R13).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R13).x = TR;
		(*R13).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 7)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R14).x; ry = (*R14).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R14).x = TR;
		(*R14).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 7)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R15).x; ry = (*R15).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R15).x = TR;
		(*R15).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 8)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R16).x; ry = (*R16).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R16).x = TR;
		(*R16).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 8)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R17).x; ry = (*R17).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R17).x = TR;
		(*R17).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 9)%100 + 0) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R18).x; ry = (*R18).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R18).x = TR;
		(*R18).y = TI;
	}

	if(TwdLarge)
	{
		T W = TW2step<T>(twiddles_large, ((10*me + 9)%100 + 100) * b );
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R19).x; ry = (*R19).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R19).x = TR;
		(*R19).y = TI;
	}

	if(rw)
	{
	bufOut[outOffset + ( 10*me + 0 + 0 ) ] = (*R0);
	bufOut[outOffset + ( 10*me + 1 + 0 ) ] = (*R2);
	bufOut[outOffset + ( 10*me + 2 + 0 ) ] = (*R4);
	bufOut[outOffset + ( 10*me + 3 + 0 ) ] = (*R6);
	bufOut[outOffset + ( 10*me + 4 + 0 ) ] = (*R8);
	bufOut[outOffset + ( 10*me + 5 + 0 ) ] = (*R10);
	bufOut[outOffset + ( 10*me + 6 + 0 ) ] = (*R12);
	bufOut[outOffset + ( 10*me + 7 + 0 ) ] = (*R14);
	bufOut[outOffset + ( 10*me + 8 + 0 ) ] = (*R16);
	bufOut[outOffset + ( 10*me + 9 + 0 ) ] = (*R18);
	bufOut[outOffset + ( 10*me + 0 + 100 ) ] = (*R1);
	bufOut[outOffset + ( 10*me + 1 + 100 ) ] = (*R3);
	bufOut[outOffset + ( 10*me + 2 + 100 ) ] = (*R5);
	bufOut[outOffset + ( 10*me + 3 + 100 ) ] = (*R7);
	bufOut[outOffset + ( 10*me + 4 + 100 ) ] = (*R9);
	bufOut[outOffset + ( 10*me + 5 + 100 ) ] = (*R11);
	bufOut[outOffset + ( 10*me + 6 + 100 ) ] = (*R13);
	bufOut[outOffset + ( 10*me + 7 + 100 ) ] = (*R15);
	bufOut[outOffset + ( 10*me + 8 + 100 ) ] = (*R17);
	bufOut[outOffset + ( 10*me + 9 + 100 ) ] = (*R19);
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
fwd_len200_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  lwbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass1_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass2_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
back_len200_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, T *lwbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  lwbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass1_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass2_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
fwd_len200_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  lwbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass1_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass2_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
back_len200_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *lwbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  lwbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass1_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass2_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
fwd_len200_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass1_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass2_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
back_len200_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *lwbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass1_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass2_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  lwbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
fwd_len200_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass1_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass2_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}

template <typename T, StrideBin sb, bool TwdLarge>
__device__ void 
back_len200_sbcc_device(const T *twiddles, const T *twiddles_large, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass1_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass2_len200_sbcc<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(100)
fft_fwd_ip_len200_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2010), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len200_sbcc>(
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
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/10);

			ioOffset += (counter_mod/currentLength)*stride_in[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/10));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/10));
		tileOffset_x	= 10;

		ioOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			ioOffset *= (stride_in[1]);
		ioOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

	}
	lwb = gb + ioOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0 = lwb[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*10 + (me%10)*200 + (me/10)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/10))*10 + t*10 + (me/10);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len200_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_in[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		lwb[(me%10) + (me/10)*stride_in[0] + t*stride_in[0]*10] = R0;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(100)
fft_back_ip_len200_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gb = gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2010), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len200_sbcc>(
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
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/10);

			ioOffset += (counter_mod/currentLength)*stride_in[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/10));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/10));
		tileOffset_x	= 10;

		ioOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			ioOffset *= (stride_in[1]);
		ioOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

	}
	lwb = gb + ioOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0 = lwb[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*10 + (me%10)*200 + (me/10)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/10))*10 + t*10 + (me/10);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len200_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_in[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		lwb[(me%10) + (me/10)*stride_in[0] + t*stride_in[0]*10] = R0;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(100)
fft_fwd_ip_len200_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2010), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_len200_sbcc>(
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
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/10);

			ioOffset += (counter_mod/currentLength)*stride_in[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/10));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/10));
		tileOffset_x	= 10;

		ioOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			ioOffset *= (stride_in[1]);
		ioOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

	}
	lwbRe = gbRe + ioOffset;
	lwbIm = gbIm + ioOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0.x = lwbRe[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		R0.y = lwbIm[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*10 + (me%10)*200 + (me/10)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/10))*10 + t*10 + (me/10);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len200_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_in[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		lwbRe[(me%10) + (me/10)*stride_in[0] + t*stride_in[0]*10] = R0.x;
		lwbIm[(me%10) + (me/10)*stride_in[0] + t*stride_in[0]*10] = R0.y;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(100)
fft_back_ip_len200_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbRe = gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbIm = gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2010), cgh);
	cgh.parallel_for<class kern_fft_back_ip_len200_sbcc>(
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
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/10);

			ioOffset += (counter_mod/currentLength)*stride_in[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/10));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/10));
		tileOffset_x	= 10;

		ioOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			ioOffset *= (stride_in[1]);
		ioOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

	}
	lwbRe = gbRe + ioOffset;
	lwbIm = gbIm + ioOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0.x = lwbRe[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		R0.y = lwbIm[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*10 + (me%10)*200 + (me/10)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/10))*10 + t*10 + (me/10);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len200_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_in[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		lwbRe[(me%10) + (me/10)*stride_in[0] + t*stride_in[0]*10] = R0.x;
		lwbIm[(me%10) + (me/10)*stride_in[0] + t*stride_in[0]*10] = R0.y;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(100)
fft_fwd_op_len200_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2010), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len200_sbcc>(
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
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/10);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/10));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/10));
		tileOffset_x	= 10;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbIn = gbIn + iOffset;
	lwbOut = gbOut + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0 = lwbIn[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*10 + (me%10)*200 + (me/10)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/10))*10 + t*10 + (me/10);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len200_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		lwbOut[(me%10) + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(100)
fft_back_op_len200_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2010), cgh);
	cgh.parallel_for<class kern_fft_back_op_len200_sbcc>(
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
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/10);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/10));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/10));
		tileOffset_x	= 10;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbIn = gbIn + iOffset;
	lwbOut = gbOut + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0 = lwbIn[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*10 + (me%10)*200 + (me/10)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/10))*10 + t*10 + (me/10);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len200_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		lwbOut[(me%10) + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(100)
fft_fwd_op_len200_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutRe = gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutIm = gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2010), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len200_sbcc>(
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
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/10);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/10));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/10));
		tileOffset_x	= 10;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbIn = gbIn + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0 = lwbIn[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*10 + (me%10)*200 + (me/10)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/10))*10 + t*10 + (me/10);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len200_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		lwbOutRe[(me%10) + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.x;
		lwbOutIm[(me%10) + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.y;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(100)
fft_back_op_len200_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbIn = gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutRe = gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutIm = gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2010), cgh);
	cgh.parallel_for<class kern_fft_back_op_len200_sbcc>(
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
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/10);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/10));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/10));
		tileOffset_x	= 10;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbIn = gbIn + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0 = lwbIn[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*10 + (me%10)*200 + (me/10)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/10))*10 + t*10 + (me/10);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len200_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		lwbOutRe[(me%10) + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.x;
		lwbOutIm[(me%10) + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.y;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(100)
fft_fwd_op_len200_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbInRe = gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbInIm = gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2010), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len200_sbcc>(
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
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/10);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/10));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/10));
		tileOffset_x	= 10;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOut = gbOut + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0.x = lwbInRe[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		R0.y = lwbInIm[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*10 + (me%10)*200 + (me/10)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/10))*10 + t*10 + (me/10);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len200_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		lwbOut[(me%10) + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(100)
fft_back_op_len200_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<T, 1> gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbInRe = gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbInIm = gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOut = gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2010), cgh);
	cgh.parallel_for<class kern_fft_back_op_len200_sbcc>(
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
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/10);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/10));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/10));
		tileOffset_x	= 10;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOut = gbOut + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0.x = lwbInRe[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		R0.y = lwbInIm[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*10 + (me%10)*200 + (me/10)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/10))*10 + t*10 + (me/10);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len200_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		lwbOut[(me%10) + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(100)
fft_fwd_op_len200_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbInRe = gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbInIm = gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutRe = gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutIm = gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2010), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_len200_sbcc>(
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
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/10);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/10));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/10));
		tileOffset_x	= 10;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0.x = lwbInRe[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		R0.y = lwbInIm[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*10 + (me%10)*200 + (me/10)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/10))*10 + t*10 + (me/10);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		fwd_len200_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		lwbOutRe[(me%10) + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.x;
		lwbOutIm[(me%10) + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.y;
	}

});
});
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, bool TwdLarge>
__global__ void
__launch_bounds__(100)
fft_back_op_len200_sbcc( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, cl::sycl::buffer<T, 1> twiddles_large_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *lengths_GB, cl::sycl::buffer<size_t, 1> *stride_in_GB, cl::sycl::buffer<size_t, 1> *stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto twiddles_large = twiddles_large_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto lengths = lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_in = stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto stride_out = stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto gbInRe = gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbInIm = gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutRe = gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto gbOutIm = gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<T,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2010), cgh);
	cgh.parallel_for<class kern_fft_back_op_len200_sbcc>(
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
	{
		// SBCC+SBRC fold higher dimensions into the batch_count, so we need
		// extra math to work out how many 'true' batches we really have
		unsigned int batch_block_size = hipGridDim_x / batch_count; //To opt: it can be calc on host.
		unsigned int counter_mod = batch % batch_block_size;
		unsigned int batch_local_count = batch / batch_block_size; //To check: technically it should be done in one instruction.
		for(int i = dim; i>2; i--){
			int currentLength = 1;
			for(int j=2; j<i; j++){
				currentLength *= lengths[j];
			}
			currentLength *= (lengths[1]/10);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/10));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/10));
		tileOffset_x	= 10;

		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
			iOffset *= (stride_in[1]);
		iOffset += (batch_local_count * stride_in[dim >= 3 ? dim-1 : 2]);

		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10;

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		oOffset += (batch_local_count * stride_out[dim >= 3 ? dim-1 : 2]);

	}
	lwbInRe = gbInRe + iOffset;
	lwbInIm = gbInIm + iOffset;
	lwbOutRe = gbOutRe + oOffset;
	lwbOutIm = gbOutIm + oOffset;


	unsigned int lds_row_padding = 0;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1657
	{
		T R0;
		// Calc global offset within a tile and read
		R0.x = lwbInRe[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		R0.y = lwbInIm[(me%10) * stride_in[1] + (me/10)*stride_in[0] + t*stride_in[0]*10];
		// Write into lds in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		lds[t*10 + (me%10)*200 + (me/10)] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1847
	{

		b = (batch % (lengths[1]/10))*10 + t*10 + (me/10);

		// Perform FFT input: lwb(In) ; output: lwb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, lwb, lds
		back_len200_sbcc_device<T, sb, TwdLarge>(twiddles, twiddles_large, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1970
	{
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		lwbOutRe[(me%10) + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.x;
		lwbOutIm[(me%10) + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.y;
	}

});
});
}

