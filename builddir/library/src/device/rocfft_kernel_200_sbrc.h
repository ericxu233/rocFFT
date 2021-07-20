#pragma once
#include "rocfft_butterfly_template.h"
#include "real2complex.h"


////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass0_len200_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
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

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass1_len200_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
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

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass2_len200_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
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

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass0_len200_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
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

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass1_len200_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
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

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass2_len200_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3, T *R4, T *R5, T *R6, T *R7, T *R8, T *R9, T *R10, T *R11, T *R12, T *R13, T *R14, T *R15, T *R16, T *R17, T *R18, T *R19)
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
template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len200_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, T *gbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  gbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass1_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass2_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len200_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, T *gbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  gbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass1_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass2_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len200_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  gbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass1_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass2_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len200_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  gbIn, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass1_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass2_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len200_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *gbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass1_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass2_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len200_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *gbOut, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass1_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass2_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  gbOut, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len200_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	FwdPass0_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass1_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	FwdPass2_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len200_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3, R4, R5, R6, R7, R8, R9, R10, R11, R12, R13, R14, R15, R16, R17, R18, R19;
	InvPass0_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass1_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
	InvPass2_len200_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3, &R4, &R5, &R6, &R7, &R8, &R9, &R10, &R11, &R12, &R13, &R14, &R15, &R16, &R17, &R18, &R19);
	__syncthreads();
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_fwd_op_len200_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
{

	extern __shared__ T lds[];
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 10 - 1) / 10);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 10 - 1) / 10);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
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
		tileOffset_x	= 10*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 200 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 200;
		iOffset += tileBlockIdx_y * (10 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (10 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 10;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (10 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 10 + t % 10 + me / 200 * 10 * 200 / 100 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = load_cb(gbIn, iOffset + me + t*100, load_cb_data, nullptr);
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = load_cb(gbIn, iOffset + me % 200 * stride_in[0] + ((me /200 * 10) + t % 10)*stride_in[2] + t / 10 * 100 * stride_in[0], load_cb_data, nullptr);
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 100 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 10 *200 + t / 10 * 100 + me % 200 + me / 200 * 4000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		fwd_len200_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 10; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 200 - me, 200, 100, &lds[r * (200 + lds_row_padding)], 0, &twiddles[200], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 10 + me % 10 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			store_cb(gbOut, oOffset + (me%10) * stride_out[1] + (me/10)*stride_out[0] + t*stride_out[0]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			store_cb(gbOut, oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[1] + t*stride_out[1]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			store_cb(gbOut, oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 10)
		{
			unsigned int t = 20;
			T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

			store_cb(gbOut, oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_back_op_len200_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
{

	extern __shared__ T lds[];
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 10 - 1) / 10);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 10 - 1) / 10);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
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
		tileOffset_x	= 10*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 200 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 200;
		iOffset += tileBlockIdx_y * (10 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (10 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 10;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (10 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 10 + t % 10 + me / 200 * 10 * 200 / 100 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = load_cb(gbIn, iOffset + me + t*100, load_cb_data, nullptr);
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = load_cb(gbIn, iOffset + me % 200 * stride_in[0] + ((me /200 * 10) + t % 10)*stride_in[2] + t / 10 * 100 * stride_in[0], load_cb_data, nullptr);
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 100 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 10 *200 + t / 10 * 100 + me % 200 + me / 200 * 4000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		back_len200_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 10; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 200 - me, 200, 100, &lds[r * (200 + lds_row_padding)], 0, &twiddles[200], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 10 + me % 10 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			store_cb(gbOut, oOffset + (me%10) * stride_out[1] + (me/10)*stride_out[0] + t*stride_out[0]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			store_cb(gbOut, oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[1] + t*stride_out[1]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			store_cb(gbOut, oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 10)
		{
			unsigned int t = 20;
			T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

			store_cb(gbOut, oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_fwd_op_len200_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	extern __shared__ T lds[];
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 10 - 1) / 10);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 10 - 1) / 10);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
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
		tileOffset_x	= 10*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 200 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 200;
		iOffset += tileBlockIdx_y * (10 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (10 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 10;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (10 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 10 + t % 10 + me / 200 * 10 * 200 / 100 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = load_cb(gbIn, iOffset + me + t*100, load_cb_data, nullptr);
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = load_cb(gbIn, iOffset + me % 200 * stride_in[0] + ((me /200 * 10) + t % 10)*stride_in[2] + t / 10 * 100 * stride_in[0], load_cb_data, nullptr);
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 100 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 10 *200 + t / 10 * 100 + me % 200 + me / 200 * 4000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		fwd_len200_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 10; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 200 - me, 200, 100, &lds[r * (200 + lds_row_padding)], 0, &twiddles[200], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 10 + me % 10 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			gbOutRe[oOffset + (me%10) * stride_out[1] + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[1] + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			gbOutRe[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[1] + t*stride_out[1]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[1] + t*stride_out[1]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			gbOutRe[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 10)
		{
			unsigned int t = 20;
			T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

			gbOutRe[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_back_op_len200_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	extern __shared__ T lds[];
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 10 - 1) / 10);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 10 - 1) / 10);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
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
		tileOffset_x	= 10*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 200 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 200;
		iOffset += tileBlockIdx_y * (10 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (10 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 10;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (10 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 10 + t % 10 + me / 200 * 10 * 200 / 100 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = load_cb(gbIn, iOffset + me + t*100, load_cb_data, nullptr);
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = load_cb(gbIn, iOffset + me % 200 * stride_in[0] + ((me /200 * 10) + t % 10)*stride_in[2] + t / 10 * 100 * stride_in[0], load_cb_data, nullptr);
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 100 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 10 *200 + t / 10 * 100 + me % 200 + me / 200 * 4000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		back_len200_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 10; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 200 - me, 200, 100, &lds[r * (200 + lds_row_padding)], 0, &twiddles[200], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 10 + me % 10 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			gbOutRe[oOffset + (me%10) * stride_out[1] + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[1] + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			gbOutRe[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[1] + t*stride_out[1]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[1] + t*stride_out[1]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			gbOutRe[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 10)
		{
			unsigned int t = 20;
			T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

			gbOutRe[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_fwd_op_len200_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
{

	extern __shared__ T lds[];
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 10 - 1) / 10);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 10 - 1) / 10);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
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
		tileOffset_x	= 10*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 200 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 200;
		iOffset += tileBlockIdx_y * (10 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (10 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 10;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (10 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 10 + t % 10 + me / 200 * 10 * 200 / 100 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = gbInRe[iOffset + me + t*100];
			R0.y = gbInIm[iOffset + me + t*100];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = gbInRe[iOffset + me % 200 * stride_in[0] + ((me /200 * 10) + t % 10)*stride_in[2] + t / 10 * 100 * stride_in[0]];
			R0.y = gbInIm[iOffset + me % 200 * stride_in[0] + ((me /200 * 10) + t % 10)*stride_in[2] + t / 10 * 100 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 100 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 10 *200 + t / 10 * 100 + me % 200 + me / 200 * 4000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		fwd_len200_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 10; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 200 - me, 200, 100, &lds[r * (200 + lds_row_padding)], 0, &twiddles[200], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 10 + me % 10 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			store_cb(gbOut, oOffset + (me%10) * stride_out[1] + (me/10)*stride_out[0] + t*stride_out[0]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			store_cb(gbOut, oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[1] + t*stride_out[1]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			store_cb(gbOut, oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 10)
		{
			unsigned int t = 20;
			T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

			store_cb(gbOut, oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_back_op_len200_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
{

	extern __shared__ T lds[];
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 10 - 1) / 10);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 10 - 1) / 10);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
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
		tileOffset_x	= 10*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 200 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 200;
		iOffset += tileBlockIdx_y * (10 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (10 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 10;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (10 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 10 + t % 10 + me / 200 * 10 * 200 / 100 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = gbInRe[iOffset + me + t*100];
			R0.y = gbInIm[iOffset + me + t*100];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = gbInRe[iOffset + me % 200 * stride_in[0] + ((me /200 * 10) + t % 10)*stride_in[2] + t / 10 * 100 * stride_in[0]];
			R0.y = gbInIm[iOffset + me % 200 * stride_in[0] + ((me /200 * 10) + t % 10)*stride_in[2] + t / 10 * 100 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 100 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 10 *200 + t / 10 * 100 + me % 200 + me / 200 * 4000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		back_len200_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 10; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 200 - me, 200, 100, &lds[r * (200 + lds_row_padding)], 0, &twiddles[200], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 10 + me % 10 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			store_cb(gbOut, oOffset + (me%10) * stride_out[1] + (me/10)*stride_out[0] + t*stride_out[0]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			store_cb(gbOut, oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[1] + t*stride_out[1]*10, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			store_cb(gbOut, oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 10)
		{
			unsigned int t = 20;
			T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

			store_cb(gbOut, oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10, R0, store_cb_data, nullptr );
		}
	}
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_fwd_op_len200_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	extern __shared__ T lds[];
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 10 - 1) / 10);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 10 - 1) / 10);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
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
		tileOffset_x	= 10*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 200 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 200;
		iOffset += tileBlockIdx_y * (10 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (10 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 10;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (10 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 10 + t % 10 + me / 200 * 10 * 200 / 100 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = gbInRe[iOffset + me + t*100];
			R0.y = gbInIm[iOffset + me + t*100];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = gbInRe[iOffset + me % 200 * stride_in[0] + ((me /200 * 10) + t % 10)*stride_in[2] + t / 10 * 100 * stride_in[0]];
			R0.y = gbInIm[iOffset + me % 200 * stride_in[0] + ((me /200 * 10) + t % 10)*stride_in[2] + t / 10 * 100 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 100 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 10 *200 + t / 10 * 100 + me % 200 + me / 200 * 4000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		fwd_len200_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 10; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 200 - me, 200, 100, &lds[r * (200 + lds_row_padding)], 0, &twiddles[200], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 10 + me % 10 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			gbOutRe[oOffset + (me%10) * stride_out[1] + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[1] + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			gbOutRe[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[1] + t*stride_out[1]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[1] + t*stride_out[1]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			gbOutRe[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 10)
		{
			unsigned int t = 20;
			T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

			gbOutRe[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}
}

//Kernel configuration: number of threads per thread block: 100, transforms: 6, Passes: 3
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(100)
fft_back_op_len200_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
{

	extern __shared__ T lds[];
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
	// for SBRC 3D kernels: number of blocks needed to deal with one 3D batch
	unsigned int blocks_per_batch;
	if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		blocks_per_batch = lengths[1] * ((lengths[2] + 10 - 1) / 10);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 10 - 1) / 10);
	// index of this block in the current 3D batch
	unsigned int block_in_batch = batch % blocks_per_batch;

	if(Tsbrc == SBRC_2D)
	{
		unsigned int counter_mod = batch;
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
		tileOffset_x	= 10*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 10 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 200 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 10;
		unsigned int tileBlockIdx_x = ((bid / 10) + tileBlockIdx_y) % 200;
		iOffset += tileBlockIdx_y * (10 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (10 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 10;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (10 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 10 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 10 + t % 10 + me / 200 * 10 * 200 / 100 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = gbInRe[iOffset + me + t*100];
			R0.y = gbInIm[iOffset + me + t*100];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = gbInRe[iOffset + me % 200 * stride_in[0] + ((me /200 * 10) + t % 10)*stride_in[2] + t / 10 * 100 * stride_in[0]];
			R0.y = gbInIm[iOffset + me % 200 * stride_in[0] + ((me /200 * 10) + t % 10)*stride_in[2] + t / 10 * 100 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 100 + t / 2 * lds_row_padding + me] = R0;
		else
			lds[t % 10 *200 + t / 10 * 100 + me % 200 + me / 200 * 4000] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<1; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		back_len200_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%10, t * (200 + lds_row_padding) * 10 + (me/10)*(200+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 10; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 200 - me, 200, 100, &lds[r * (200 + lds_row_padding)], 0, &twiddles[200], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<20; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 10 + me % 10 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			gbOutRe[oOffset + (me%10) * stride_out[1] + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[1] + (me/10)*stride_out[0] + t*stride_out[0]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			gbOutRe[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[1] + t*stride_out[1]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[1] + t*stride_out[1]*10] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			gbOutRe[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 10)
		{
			unsigned int t = 20;
			T R0 = lds[t*10 + (me%10)*(200 + lds_row_padding) + (me/10)];

			gbOutRe[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.x;
			gbOutIm[oOffset + (me%10) * stride_out[0] + (me/10)*stride_out[2] + t*stride_out[2]*10] = R0.y;
		}
	}
}

