#pragma once
#include "rocfft_butterfly_template.h"
#include "real2complex.h"


////////////////////////////////////////Passes kernels
template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass0_len256_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	FwdRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 1 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 2 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 3 ) ] = (*R3);
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass1_len256_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	{
		T W = twiddles[3 + 3*((1*me + 0)%4) + 0];
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
		T W = twiddles[3 + 3*((1*me + 0)%4) + 1];
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
		T W = twiddles[3 + 3*((1*me + 0)%4) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	FwdRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 4 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 8 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 12 ) ] = (*R3);
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass2_len256_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	{
		T W = twiddles[15 + 3*((1*me + 0)%16) + 0];
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
		T W = twiddles[15 + 3*((1*me + 0)%16) + 1];
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
		T W = twiddles[15 + 3*((1*me + 0)%16) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	FwdRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 16 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 32 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 48 ) ] = (*R3);
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
FwdPass3_len256_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	{
		T W = twiddles[63 + 3*((1*me + 0)%64) + 0];
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
		T W = twiddles[63 + 3*((1*me + 0)%64) + 1];
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
		T W = twiddles[63 + 3*((1*me + 0)%64) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx - wy * ry;
		TI = wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	FwdRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( 1*me + 0 + 0 ) ] = (*R0);
	bufOut[outOffset + ( 1*me + 0 + 64 ) ] = (*R1);
	bufOut[outOffset + ( 1*me + 0 + 128 ) ] = (*R2);
	bufOut[outOffset + ( 1*me + 0 + 192 ) ] = (*R3);
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass0_len256_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	InvRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 1 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 2 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/1)*4 + (1*me + 0)%1 + 3 ) ] = (*R3);
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass1_len256_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	{
		T W = twiddles[3 + 3*((1*me + 0)%4) + 0];
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
		T W = twiddles[3 + 3*((1*me + 0)%4) + 1];
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
		T W = twiddles[3 + 3*((1*me + 0)%4) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	InvRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 4 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 8 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/4)*16 + (1*me + 0)%4 + 12 ) ] = (*R3);
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass2_len256_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	{
		T W = twiddles[15 + 3*((1*me + 0)%16) + 0];
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
		T W = twiddles[15 + 3*((1*me + 0)%16) + 1];
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
		T W = twiddles[15 + 3*((1*me + 0)%16) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	InvRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 0 ) ] = (*R0);
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 16 ) ] = (*R1);
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 32 ) ] = (*R2);
	bufOut[outOffset + ( ((1*me + 0)/16)*64 + (1*me + 0)%16 + 48 ) ] = (*R3);
	}

}

template <typename T, StrideBin sb, CallbackType cbtype>
__device__ void
InvPass3_len256_sbrc(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int inOffset, unsigned int outOffset, T *bufIn, T *bufOut, T *R0, T *R1, T *R2, T *R3)
{


	if(rw)
	{
	(*R0) = bufIn[inOffset + ( 0 + me*1 + 0 + 0 ) ];
	(*R1) = bufIn[inOffset + ( 0 + me*1 + 0 + 64 ) ];
	(*R2) = bufIn[inOffset + ( 0 + me*1 + 0 + 128 ) ];
	(*R3) = bufIn[inOffset + ( 0 + me*1 + 0 + 192 ) ];
	}



	{
		T W = twiddles[63 + 3*((1*me + 0)%64) + 0];
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
		T W = twiddles[63 + 3*((1*me + 0)%64) + 1];
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
		T W = twiddles[63 + 3*((1*me + 0)%64) + 2];
		real_type_t<T> TR, TI;
		real_type_t<T>wx, wy, rx, ry;
		wx = W.x; wy = W.y;
		rx = (*R3).x; ry = (*R3).y;
		TR = wx * rx + wy * ry;
		TI = - wy * rx + wx * ry;
		(*R3).x = TR;
		(*R3).y = TI;
	}

	InvRad4B1(R0, R1, R2, R3);

	__syncthreads();



	if(rw)
	{
	bufOut[outOffset + ( 1*me + 0 + 0 ) ] = (*R0);
	bufOut[outOffset + ( 1*me + 0 + 64 ) ] = (*R1);
	bufOut[outOffset + ( 1*me + 0 + 128 ) ] = (*R2);
	bufOut[outOffset + ( 1*me + 0 + 192 ) ] = (*R3);
	}

}


////////////////////////////////////////Encapsulated passes kernels
template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len256_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, T *gbOut, T *lds)
{
	T R0, R1, R2, R3;
	FwdPass0_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  gbIn, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass1_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass2_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass3_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  gbOut, &R0, &R1, &R2, &R3);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len256_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, T *gbOut, T *lds)
{
	T R0, R1, R2, R3;
	InvPass0_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  gbIn, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass1_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass2_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass3_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  gbOut, &R0, &R1, &R2, &R3);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len256_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3;
	FwdPass0_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  gbIn, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass1_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass2_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass3_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len256_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, T *gbIn, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3;
	InvPass0_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  gbIn, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass1_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass2_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass3_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len256_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *gbOut, T *lds)
{
	T R0, R1, R2, R3;
	FwdPass0_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass1_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass2_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass3_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  gbOut, &R0, &R1, &R2, &R3);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len256_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, T *gbOut, T *lds)
{
	T R0, R1, R2, R3;
	InvPass0_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass1_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass2_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass3_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds,  gbOut, &R0, &R1, &R2, &R3);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
fwd_len256_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3;
	FwdPass0_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass1_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass2_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	FwdPass3_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3);
	__syncthreads();
}

template <typename T, StrideBin sb,CallbackType cbtype>
__device__ void 
back_len256_sbrc_device(const T * const twiddles, const size_t stride_in, const size_t stride_out, unsigned int rw, unsigned int b, unsigned int me, unsigned int ldsOffset, real_type_t<T> *bufInRe, real_type_t<T> *bufInIm, real_type_t<T> *bufOutRe, real_type_t<T> *bufOutIm, T *lds)
{
	T R0, R1, R2, R3;
	InvPass0_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset,  bufInRe, bufInIm, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass1_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass2_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, lds, &R0, &R1, &R2, &R3);
	__syncthreads();
	InvPass3_len256_sbrc<T, sb, cbtype>(twiddles, stride_in, stride_out, rw, b, me, ldsOffset, ldsOffset, lds, bufOutRe, bufOutIm, &R0, &R1, &R2, &R3);
	__syncthreads();
}


////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(256)
fft_fwd_op_len256_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
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
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
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
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 256 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 32;
		unsigned int tileBlockIdx_x = ((bid / 32) + tileBlockIdx_y) % 256;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 256 * 8 * 256 / 256 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = load_cb(gbIn, iOffset + me + t*256, load_cb_data, nullptr);
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = load_cb(gbIn, iOffset + me % 256 * stride_in[0] + ((me /256 * 32) + t % 8)*stride_in[2] + t / 8 * 256 * stride_in[0], load_cb_data, nullptr);
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 256 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *256 + t / 8 * 256 + me % 256 + me / 256 * 2048] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		fwd_len256_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 256 - me, 256, 128, &lds[r * (256 + lds_row_padding)], 0, &twiddles[256], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			store_cb(gbOut, oOffset + (me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*32, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			store_cb(gbOut, oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*32, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			store_cb(gbOut, oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32, R0, store_cb_data, nullptr );
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

			store_cb(gbOut, oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32, R0, store_cb_data, nullptr );
		}
	}
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(256)
fft_back_op_len256_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, T * __restrict__ gbOut)
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
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
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
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 256 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 32;
		unsigned int tileBlockIdx_x = ((bid / 32) + tileBlockIdx_y) % 256;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 256 * 8 * 256 / 256 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = load_cb(gbIn, iOffset + me + t*256, load_cb_data, nullptr);
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = load_cb(gbIn, iOffset + me % 256 * stride_in[0] + ((me /256 * 32) + t % 8)*stride_in[2] + t / 8 * 256 * stride_in[0], load_cb_data, nullptr);
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 256 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *256 + t / 8 * 256 + me % 256 + me / 256 * 2048] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		back_len256_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 256 - me, 256, 128, &lds[r * (256 + lds_row_padding)], 0, &twiddles[256], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			store_cb(gbOut, oOffset + (me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*32, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			store_cb(gbOut, oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*32, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			store_cb(gbOut, oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32, R0, store_cb_data, nullptr );
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

			store_cb(gbOut, oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32, R0, store_cb_data, nullptr );
		}
	}
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(256)
fft_fwd_op_len256_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
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
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
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
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 256 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 32;
		unsigned int tileBlockIdx_x = ((bid / 32) + tileBlockIdx_y) % 256;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 256 * 8 * 256 / 256 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = load_cb(gbIn, iOffset + me + t*256, load_cb_data, nullptr);
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = load_cb(gbIn, iOffset + me % 256 * stride_in[0] + ((me /256 * 32) + t % 8)*stride_in[2] + t / 8 * 256 * stride_in[0], load_cb_data, nullptr);
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 256 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *256 + t / 8 * 256 + me % 256 + me / 256 * 2048] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		fwd_len256_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 256 - me, 256, 128, &lds[r * (256 + lds_row_padding)], 0, &twiddles[256], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			gbOutRe[oOffset + (me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			gbOutRe[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*32] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			gbOutRe[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

			gbOutRe[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.y;
		}
	}
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(256)
fft_back_op_len256_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ gbIn, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
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
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
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
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 256 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 32;
		unsigned int tileBlockIdx_x = ((bid / 32) + tileBlockIdx_y) % 256;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 256 * 8 * 256 / 256 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0 = load_cb(gbIn, iOffset + me + t*256, load_cb_data, nullptr);
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0 = load_cb(gbIn, iOffset + me % 256 * stride_in[0] + ((me /256 * 32) + t % 8)*stride_in[2] + t / 8 * 256 * stride_in[0], load_cb_data, nullptr);
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 256 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *256 + t / 8 * 256 + me % 256 + me / 256 * 2048] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		back_len256_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 256 - me, 256, 128, &lds[r * (256 + lds_row_padding)], 0, &twiddles[256], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			gbOutRe[oOffset + (me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			gbOutRe[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*32] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			gbOutRe[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

			gbOutRe[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.y;
		}
	}
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(256)
fft_fwd_op_len256_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
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
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
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
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 256 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 32;
		unsigned int tileBlockIdx_x = ((bid / 32) + tileBlockIdx_y) % 256;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 256 * 8 * 256 / 256 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = gbInRe[iOffset + me + t*256];
			R0.y = gbInIm[iOffset + me + t*256];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = gbInRe[iOffset + me % 256 * stride_in[0] + ((me /256 * 32) + t % 8)*stride_in[2] + t / 8 * 256 * stride_in[0]];
			R0.y = gbInIm[iOffset + me % 256 * stride_in[0] + ((me /256 * 32) + t % 8)*stride_in[2] + t / 8 * 256 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 256 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *256 + t / 8 * 256 + me % 256 + me / 256 * 2048] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		fwd_len256_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 256 - me, 256, 128, &lds[r * (256 + lds_row_padding)], 0, &twiddles[256], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			store_cb(gbOut, oOffset + (me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*32, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			store_cb(gbOut, oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*32, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			store_cb(gbOut, oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32, R0, store_cb_data, nullptr );
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

			store_cb(gbOut, oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32, R0, store_cb_data, nullptr );
		}
	}
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(256)
fft_back_op_len256_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, T * __restrict__ gbOut)
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
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
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
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 256 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 32;
		unsigned int tileBlockIdx_x = ((bid / 32) + tileBlockIdx_y) % 256;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 256 * 8 * 256 / 256 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = gbInRe[iOffset + me + t*256];
			R0.y = gbInIm[iOffset + me + t*256];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = gbInRe[iOffset + me % 256 * stride_in[0] + ((me /256 * 32) + t % 8)*stride_in[2] + t / 8 * 256 * stride_in[0]];
			R0.y = gbInIm[iOffset + me % 256 * stride_in[0] + ((me /256 * 32) + t % 8)*stride_in[2] + t / 8 * 256 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 256 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *256 + t / 8 * 256 + me % 256 + me / 256 * 2048] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		back_len256_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 256 - me, 256, 128, &lds[r * (256 + lds_row_padding)], 0, &twiddles[256], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			store_cb(gbOut, oOffset + (me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*32, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			store_cb(gbOut, oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*32, R0, store_cb_data, nullptr );
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			store_cb(gbOut, oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32, R0, store_cb_data, nullptr );
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

			store_cb(gbOut, oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32, R0, store_cb_data, nullptr );
		}
	}
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(256)
fft_fwd_op_len256_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
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
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
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
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 256 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 32;
		unsigned int tileBlockIdx_x = ((bid / 32) + tileBlockIdx_y) % 256;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 256 * 8 * 256 / 256 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = gbInRe[iOffset + me + t*256];
			R0.y = gbInIm[iOffset + me + t*256];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = gbInRe[iOffset + me % 256 * stride_in[0] + ((me /256 * 32) + t % 8)*stride_in[2] + t / 8 * 256 * stride_in[0]];
			R0.y = gbInIm[iOffset + me % 256 * stride_in[0] + ((me /256 * 32) + t % 8)*stride_in[2] + t / 8 * 256 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 256 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *256 + t / 8 * 256 + me % 256 + me / 256 * 2048] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		fwd_len256_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 256 - me, 256, 128, &lds[r * (256 + lds_row_padding)], 0, &twiddles[256], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			gbOutRe[oOffset + (me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			gbOutRe[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*32] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			gbOutRe[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

			gbOutRe[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.y;
		}
	}
}

//Kernel configuration: number of threads per thread block: 256, transforms: 1, Passes: 4
template <typename T, StrideBin sb, SBRC_TYPE Tsbrc, SBRC_TRANSPOSE_TYPE Ttranspose, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(256)
fft_back_op_len256_sbrc( const T * __restrict__ twiddles, const size_t dim, const size_t *lengths, const size_t *stride_in, const size_t *stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ gbInRe, real_type_t<T> * __restrict__ gbInIm, real_type_t<T> * __restrict__ gbOutRe, real_type_t<T> * __restrict__ gbOutIm)
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
		blocks_per_batch = lengths[1] * ((lengths[2] + 8 - 1) / 8);
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		blocks_per_batch = lengths[2] * ((lengths[1] + 8 - 1) / 8);
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
			currentLength *= (lengths[1]/8);

			iOffset += (counter_mod/currentLength)*stride_in[i];
			oOffset += (counter_mod/currentLength)*stride_out[i];
			counter_mod = (counter_mod % currentLength); 
		}

		// We handle a 2D tile block with one work-group threads.
		// In the below, '_x' moves along the fast dimension of the tile.
		unsigned int tileIdx_x, tileIdx_y, tileOffset_x, tileOffset_y;

		// Calc input tile start offset
		tileIdx_y		= (counter_mod / (lengths[1]/8));
		tileOffset_y	= stride_in[2];
		tileIdx_x		= (counter_mod % (lengths[1]/8));
		tileOffset_x	= 8*lengths[0];
		iOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
		// Calc output tile start offset
		tileOffset_y	= stride_out[2];
		tileOffset_x	= 8 * stride_out[1];

		oOffset += tileIdx_y * tileOffset_y + tileIdx_x * tileOffset_x;
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
	{
		unsigned int readTileIdx_x = batch % lengths[1];
		unsigned int readTileIdx_y = batch % blocks_per_batch / lengths[1];
		if(Ttranspose == DIAGONAL)
		{
		// diagonal transpose for power of 2 length
		unsigned int bid = readTileIdx_x + 256 * readTileIdx_y;
		unsigned int tileBlockIdx_y = bid % 32;
		unsigned int tileBlockIdx_x = ((bid / 32) + tileBlockIdx_y) % 256;
		iOffset += tileBlockIdx_y * (8 * stride_in[2]) + tileBlockIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];
		unsigned int writeTileIdx_x = tileBlockIdx_y;
		unsigned int writeTileIdx_y = tileBlockIdx_x;
		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
		else
		{
		iOffset += readTileIdx_y * (8 * stride_in[2]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[2] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
		}
	}
	else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		dim3 tgs; // tile grid size
		tgs.x = 1;
		tgs.y = lengths[1] * lengths[2] / 8;
		unsigned int blocks_per_batch = tgs.x * tgs.y;
		unsigned int readTileIdx_x = 0; // batch % tgs.x;
		unsigned int readTileIdx_y = (batch % blocks_per_batch) / tgs.x;
		iOffset += readTileIdx_y * (8 * stride_in[1]) + readTileIdx_x  * stride_in[1] + batch / blocks_per_batch * stride_in[3];

		unsigned int writeTileIdx_x = readTileIdx_y;
		unsigned int writeTileIdx_y = readTileIdx_x;

		oOffset += writeTileIdx_y * stride_out[3] + writeTileIdx_x * 8 * stride_out[0] + batch / blocks_per_batch * stride_out[3];
	}

	unsigned int lds_row_padding = 0;
	auto load_cb = get_load_cb<T,cbtype>(load_cb_fn);

	if(Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		lds_row_padding = 1;
	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1477
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't read more rows than are there
			if(block_in_batch / lengths[1] * 8 + t % 8 + me / 256 * 8 * 256 / 256 >= lengths[2])
				continue;
		}
		T R0;
		// Calc global offset within a tile and read
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			R0.x = gbInRe[iOffset + me + t*256];
			R0.y = gbInIm[iOffset + me + t*256];
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			R0.x = gbInRe[iOffset + me % 256 * stride_in[0] + ((me /256 * 32) + t % 8)*stride_in[2] + t / 8 * 256 * stride_in[0]];
			R0.y = gbInIm[iOffset + me % 256 * stride_in[0] + ((me /256 * 32) + t % 8)*stride_in[2] + t / 8 * 256 * stride_in[0]];
		}

		// Write into lds in row-major
		if(Tsbrc == SBRC_2D || Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
			lds[t * 256 + t / 1 * lds_row_padding + me] = R0;
		else
			lds[t % 8 *256 + t / 8 * 256 + me % 256 + me / 256 * 2048] = R0;
	}

	__syncthreads();


	for(unsigned int t=0; t<2; t++) // generator.kernel.hpp:1669
	{

		// Perform FFT input: gb(In) ; output: gb(Out); working space: lds 
		// rw, b, me% control read/write; then ldsOffset, gb, lds
		back_len256_sbrc_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  1, b, me%64, t * (256 + lds_row_padding) * 4 + (me/64)*(256+ lds_row_padding), lds, lds, lds);

	}

	__syncthreads();


	auto store_cb = get_store_cb<T,cbtype>(store_cb_fn);
	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		for(unsigned int r = 0; r < 8; r++)
		{
			post_process_interleaved_inplace<T, true,CallbackType::NONE>(me, 256 - me, 256, 128, &lds[r * (256 + lds_row_padding)], 0, &twiddles[256], nullptr, nullptr, 0, nullptr, nullptr);
		}
		__syncthreads();
	}

	for(unsigned int t=0; t<8; t++) // generator.kernel.hpp:1800
	{
		if(Ttranspose == TILE_UNALIGNED && Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			// for rectangular cases, make sure we don't write more rows than are there
			if(block_in_batch / lengths[1] * 8 + me % 8 >= lengths[2])
				continue;
		}
		// Read from lds and write to global mem in column-major
		// In lds, the offset = blockIdx * blockOffset + threadIdx_x * threadOffset_x + threadIdx_y * 1
		// which is    R0 = lds[   t     *  (wgs/bwd)  +  (me%bwd)   *  length[0]     + (me/bwd)    * 1]
		T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

		// Calc global offset within a tile and write
		if(Tsbrc == SBRC_2D)
		{
			gbOutRe[oOffset + (me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[1] + (me/8)*stride_out[0] + t*stride_out[0]*32] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_XY_Z)
		{
			gbOutRe[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[1] + t*stride_out[1]*32] = R0.y;
		}
		else if(Tsbrc == SBRC_3D_FFT_TRANS_Z_XY || Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
		{
			gbOutRe[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.y;
		}
	}

	if (Tsbrc == SBRC_3D_FFT_ERC_TRANS_Z_XY)
	{
		if(me < 8)
		{
			unsigned int t = 8;
			T R0 = lds[t*32 + (me%8)*(256 + lds_row_padding) + (me/8)];

			gbOutRe[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.x;
			gbOutIm[oOffset + (me%8) * stride_out[0] + (me/8)*stride_out[2] + t*stride_out[2]*32] = R0.y;
		}
	}
}

