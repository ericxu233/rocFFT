#pragma once
#include "rocfft_kernel_4.h"
#include "rocfft_kernel_125.h"

////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 250, maximum transforms: 32, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(250)
fft_fwd_ip_2D_4_125( const T * __restrict__ twiddles, const size_t dim, const size_t *_lengths, const size_t *_stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ _gb)
{

	__shared__ real_type_t<T> lds[500];
	__shared__ T lds_data[4*125];
	// use supplied input stride for row transform
	size_t stride_in[4];
	stride_in[0] = _stride_in[0];
	stride_in[1] = _stride_in[1];
	stride_in[2] = _stride_in[2];
	stride_in[3] = _stride_in[3];
	// set unit output stride, since we're writing to LDS
	size_t stride_out[4];
	stride_out[0] = 1;
	stride_out[1] = _lengths[0];
	stride_out[2] = _lengths[1];
	stride_out[3] = _lengths[2];
	// use supplied lengths for row transform
	size_t lengths[3];
	lengths[0] = _lengths[0];
	lengths[1] = _lengths[1];
	lengths[2] = _lengths[2];
	// declare input/output pointers
	T* gbIn = _gb;
	// write to LDS
	T* gbOut = lds_data;
	// transform each row
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 2);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/2));
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
	fwd_len4_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*4, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, 0, nullptr, nullptr);
	}
	// we have two twiddle tables back to back in device
	// memory - move to the second table (if nonsquare)
	twiddles = twiddles + lengths[0];
	// write output to original in-place destination
	stride_out[0] = _stride_in[1];
	stride_out[1] = _stride_in[0];
	stride_out[2] = _stride_in[2];
	stride_out[3] = _stride_in[3];
	// get unit stride input from LDS
	stride_in[0] = _lengths[0];
	stride_in[1] = 1;
	stride_in[2] = _lengths[2];
	stride_in[3] = _lengths[3];
	
	// flip dimensions and transform each column

	auto temp = lengths[0];
	lengths[0] = lengths[1];
	lengths[1] = temp;
	// Let the row transform finish before starting column transform
	__syncthreads();
	// declare input/output pointers for column transform
	gbIn = lds_data;
	gbOut = _gb;
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 5);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/25));
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
	fwd_len125_device<T, SB_NONUNIT, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%25, (me/25)*125, gbIn, iOffset, gbOut, oOffset, lds, nullptr, nullptr, load_cb_lds_bytes, store_cb_fn, store_cb_data);
	}
}

//Kernel configuration: number of threads per thread block: 250, maximum transforms: 32, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(250)
fft_back_ip_2D_4_125( const T * __restrict__ twiddles, const size_t dim, const size_t *_lengths, const size_t *_stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ _gb)
{

	__shared__ real_type_t<T> lds[500];
	__shared__ T lds_data[4*125];
	// use supplied input stride for row transform
	size_t stride_in[4];
	stride_in[0] = _stride_in[0];
	stride_in[1] = _stride_in[1];
	stride_in[2] = _stride_in[2];
	stride_in[3] = _stride_in[3];
	// set unit output stride, since we're writing to LDS
	size_t stride_out[4];
	stride_out[0] = 1;
	stride_out[1] = _lengths[0];
	stride_out[2] = _lengths[1];
	stride_out[3] = _lengths[2];
	// use supplied lengths for row transform
	size_t lengths[3];
	lengths[0] = _lengths[0];
	lengths[1] = _lengths[1];
	lengths[2] = _lengths[2];
	// declare input/output pointers
	T* gbIn = _gb;
	// write to LDS
	T* gbOut = lds_data;
	// transform each row
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 2);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/2));
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
	back_len4_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*4, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, 0, nullptr, nullptr);
	}
	// we have two twiddle tables back to back in device
	// memory - move to the second table (if nonsquare)
	twiddles = twiddles + lengths[0];
	// write output to original in-place destination
	stride_out[0] = _stride_in[1];
	stride_out[1] = _stride_in[0];
	stride_out[2] = _stride_in[2];
	stride_out[3] = _stride_in[3];
	// get unit stride input from LDS
	stride_in[0] = _lengths[0];
	stride_in[1] = 1;
	stride_in[2] = _lengths[2];
	stride_in[3] = _lengths[3];
	
	// flip dimensions and transform each column

	auto temp = lengths[0];
	lengths[0] = lengths[1];
	lengths[1] = temp;
	// Let the row transform finish before starting column transform
	__syncthreads();
	// declare input/output pointers for column transform
	gbIn = lds_data;
	gbOut = _gb;
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 5);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/25));
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
	back_len125_device<T, SB_NONUNIT, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%25, (me/25)*125, gbIn, iOffset, gbOut, oOffset, lds, nullptr, nullptr, load_cb_lds_bytes, store_cb_fn, store_cb_data);
	}
}

//Kernel configuration: number of threads per thread block: 250, maximum transforms: 32, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(250)
fft_fwd_ip_2D_4_125( const T * __restrict__ twiddles, const size_t dim, const size_t *_lengths, const size_t *_stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ _gbRe, real_type_t<T> * __restrict__ _gbIm)
{

	__shared__ real_type_t<T> lds[500];
	__shared__ T lds_data[4*125];
	// use supplied input stride for row transform
	size_t stride_in[4];
	stride_in[0] = _stride_in[0];
	stride_in[1] = _stride_in[1];
	stride_in[2] = _stride_in[2];
	stride_in[3] = _stride_in[3];
	// set unit output stride, since we're writing to LDS
	size_t stride_out[4];
	stride_out[0] = 1;
	stride_out[1] = _lengths[0];
	stride_out[2] = _lengths[1];
	stride_out[3] = _lengths[2];
	// use supplied lengths for row transform
	size_t lengths[3];
	lengths[0] = _lengths[0];
	lengths[1] = _lengths[1];
	lengths[2] = _lengths[2];
	// declare input/output pointers
	real_type_t<T>* gbInRe = _gbRe;
	real_type_t<T>* gbInIm = _gbIm;
	// write to LDS
	T* gbOut = lds_data;
	// transform each row
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 2);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/2));
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
	fwd_len4_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*4, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, 0, nullptr, nullptr);
	}
	// we have two twiddle tables back to back in device
	// memory - move to the second table (if nonsquare)
	twiddles = twiddles + lengths[0];
	// write output to original in-place destination
	stride_out[0] = _stride_in[1];
	stride_out[1] = _stride_in[0];
	stride_out[2] = _stride_in[2];
	stride_out[3] = _stride_in[3];
	// get unit stride input from LDS
	stride_in[0] = _lengths[0];
	stride_in[1] = 1;
	stride_in[2] = _lengths[2];
	stride_in[3] = _lengths[3];
	
	// flip dimensions and transform each column

	auto temp = lengths[0];
	lengths[0] = lengths[1];
	lengths[1] = temp;
	// Let the row transform finish before starting column transform
	__syncthreads();
	// declare input/output pointers for column transform
	T* gbIn = lds_data;
	real_type_t<T>* gbOutRe = _gbRe;
	real_type_t<T>* gbOutIm = _gbIm;
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 5);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/25));
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
	fwd_len125_device<T, SB_NONUNIT, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%25, (me/25)*125, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, nullptr, nullptr, load_cb_lds_bytes, store_cb_fn, store_cb_data);
	}
}

//Kernel configuration: number of threads per thread block: 250, maximum transforms: 32, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(250)
fft_back_ip_2D_4_125( const T * __restrict__ twiddles, const size_t dim, const size_t *_lengths, const size_t *_stride_in, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ _gbRe, real_type_t<T> * __restrict__ _gbIm)
{

	__shared__ real_type_t<T> lds[500];
	__shared__ T lds_data[4*125];
	// use supplied input stride for row transform
	size_t stride_in[4];
	stride_in[0] = _stride_in[0];
	stride_in[1] = _stride_in[1];
	stride_in[2] = _stride_in[2];
	stride_in[3] = _stride_in[3];
	// set unit output stride, since we're writing to LDS
	size_t stride_out[4];
	stride_out[0] = 1;
	stride_out[1] = _lengths[0];
	stride_out[2] = _lengths[1];
	stride_out[3] = _lengths[2];
	// use supplied lengths for row transform
	size_t lengths[3];
	lengths[0] = _lengths[0];
	lengths[1] = _lengths[1];
	lengths[2] = _lengths[2];
	// declare input/output pointers
	real_type_t<T>* gbInRe = _gbRe;
	real_type_t<T>* gbInIm = _gbIm;
	// write to LDS
	T* gbOut = lds_data;
	// transform each row
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 2);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/2));
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
	back_len4_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*4, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, 0, nullptr, nullptr);
	}
	// we have two twiddle tables back to back in device
	// memory - move to the second table (if nonsquare)
	twiddles = twiddles + lengths[0];
	// write output to original in-place destination
	stride_out[0] = _stride_in[1];
	stride_out[1] = _stride_in[0];
	stride_out[2] = _stride_in[2];
	stride_out[3] = _stride_in[3];
	// get unit stride input from LDS
	stride_in[0] = _lengths[0];
	stride_in[1] = 1;
	stride_in[2] = _lengths[2];
	stride_in[3] = _lengths[3];
	
	// flip dimensions and transform each column

	auto temp = lengths[0];
	lengths[0] = lengths[1];
	lengths[1] = temp;
	// Let the row transform finish before starting column transform
	__syncthreads();
	// declare input/output pointers for column transform
	T* gbIn = lds_data;
	real_type_t<T>* gbOutRe = _gbRe;
	real_type_t<T>* gbOutIm = _gbIm;
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 5);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/25));
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
	back_len125_device<T, SB_NONUNIT, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%25, (me/25)*125, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, nullptr, nullptr, load_cb_lds_bytes, store_cb_fn, store_cb_data);
	}
}

//Kernel configuration: number of threads per thread block: 250, maximum transforms: 32, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(250)
fft_fwd_op_2D_4_125( const T * __restrict__ twiddles, const size_t dim, const size_t *_lengths, const size_t *_stride_in, const size_t *_stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ _gbIn, T * __restrict__ _gbOut)
{

	__shared__ real_type_t<T> lds[500];
	__shared__ T lds_data[4*125];
	// use supplied input stride for row transform
	size_t stride_in[4];
	stride_in[0] = _stride_in[0];
	stride_in[1] = _stride_in[1];
	stride_in[2] = _stride_in[2];
	stride_in[3] = _stride_in[3];
	// set unit output stride, since we're writing to LDS
	size_t stride_out[4];
	stride_out[0] = 1;
	stride_out[1] = _lengths[0];
	stride_out[2] = _lengths[1];
	stride_out[3] = _lengths[2];
	// use supplied lengths for row transform
	size_t lengths[3];
	lengths[0] = _lengths[0];
	lengths[1] = _lengths[1];
	lengths[2] = _lengths[2];
	// declare input/output pointers
	T* gbIn = _gbIn;
	// write to LDS
	T* gbOut = lds_data;
	// transform each row
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 2);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/2));
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
	fwd_len4_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*4, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, 0, nullptr, nullptr);
	}
	// we have two twiddle tables back to back in device
	// memory - move to the second table (if nonsquare)
	twiddles = twiddles + lengths[0];
	// write output to original out-of-place destination
	stride_out[0] = _stride_out[1];
	stride_out[1] = _stride_out[0];
	stride_out[2] = _stride_out[2];
	stride_out[3] = _stride_out[3];
	// get unit stride input from LDS
	stride_in[0] = _lengths[0];
	stride_in[1] = 1;
	stride_in[2] = _lengths[2];
	stride_in[3] = _lengths[3];
	
	// flip dimensions and transform each column

	auto temp = lengths[0];
	lengths[0] = lengths[1];
	lengths[1] = temp;
	// Let the row transform finish before starting column transform
	__syncthreads();
	// declare input/output pointers for column transform
	gbIn = lds_data;
	gbOut = _gbOut;
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 5);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/25));
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
	fwd_len125_device<T, SB_NONUNIT, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%25, (me/25)*125, gbIn, iOffset, gbOut, oOffset, lds, nullptr, nullptr, load_cb_lds_bytes, store_cb_fn, store_cb_data);
	}
}

//Kernel configuration: number of threads per thread block: 250, maximum transforms: 32, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(250)
fft_back_op_2D_4_125( const T * __restrict__ twiddles, const size_t dim, const size_t *_lengths, const size_t *_stride_in, const size_t *_stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ _gbIn, T * __restrict__ _gbOut)
{

	__shared__ real_type_t<T> lds[500];
	__shared__ T lds_data[4*125];
	// use supplied input stride for row transform
	size_t stride_in[4];
	stride_in[0] = _stride_in[0];
	stride_in[1] = _stride_in[1];
	stride_in[2] = _stride_in[2];
	stride_in[3] = _stride_in[3];
	// set unit output stride, since we're writing to LDS
	size_t stride_out[4];
	stride_out[0] = 1;
	stride_out[1] = _lengths[0];
	stride_out[2] = _lengths[1];
	stride_out[3] = _lengths[2];
	// use supplied lengths for row transform
	size_t lengths[3];
	lengths[0] = _lengths[0];
	lengths[1] = _lengths[1];
	lengths[2] = _lengths[2];
	// declare input/output pointers
	T* gbIn = _gbIn;
	// write to LDS
	T* gbOut = lds_data;
	// transform each row
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 2);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/2));
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
	back_len4_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*4, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, 0, nullptr, nullptr);
	}
	// we have two twiddle tables back to back in device
	// memory - move to the second table (if nonsquare)
	twiddles = twiddles + lengths[0];
	// write output to original out-of-place destination
	stride_out[0] = _stride_out[1];
	stride_out[1] = _stride_out[0];
	stride_out[2] = _stride_out[2];
	stride_out[3] = _stride_out[3];
	// get unit stride input from LDS
	stride_in[0] = _lengths[0];
	stride_in[1] = 1;
	stride_in[2] = _lengths[2];
	stride_in[3] = _lengths[3];
	
	// flip dimensions and transform each column

	auto temp = lengths[0];
	lengths[0] = lengths[1];
	lengths[1] = temp;
	// Let the row transform finish before starting column transform
	__syncthreads();
	// declare input/output pointers for column transform
	gbIn = lds_data;
	gbOut = _gbOut;
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 5);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/25));
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
	back_len125_device<T, SB_NONUNIT, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%25, (me/25)*125, gbIn, iOffset, gbOut, oOffset, lds, nullptr, nullptr, load_cb_lds_bytes, store_cb_fn, store_cb_data);
	}
}

//Kernel configuration: number of threads per thread block: 250, maximum transforms: 32, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(250)
fft_fwd_op_2D_4_125( const T * __restrict__ twiddles, const size_t dim, const size_t *_lengths, const size_t *_stride_in, const size_t *_stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ _gbIn, real_type_t<T> * __restrict__ _gbOutRe, real_type_t<T> * __restrict__ _gbOutIm)
{

	__shared__ real_type_t<T> lds[500];
	__shared__ T lds_data[4*125];
	// use supplied input stride for row transform
	size_t stride_in[4];
	stride_in[0] = _stride_in[0];
	stride_in[1] = _stride_in[1];
	stride_in[2] = _stride_in[2];
	stride_in[3] = _stride_in[3];
	// set unit output stride, since we're writing to LDS
	size_t stride_out[4];
	stride_out[0] = 1;
	stride_out[1] = _lengths[0];
	stride_out[2] = _lengths[1];
	stride_out[3] = _lengths[2];
	// use supplied lengths for row transform
	size_t lengths[3];
	lengths[0] = _lengths[0];
	lengths[1] = _lengths[1];
	lengths[2] = _lengths[2];
	// declare input/output pointers
	T* gbIn = _gbIn;
	// write to LDS
	T* gbOut = lds_data;
	// transform each row
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 2);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/2));
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
	fwd_len4_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*4, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, 0, nullptr, nullptr);
	}
	// we have two twiddle tables back to back in device
	// memory - move to the second table (if nonsquare)
	twiddles = twiddles + lengths[0];
	// write output to original out-of-place destination
	stride_out[0] = _stride_out[1];
	stride_out[1] = _stride_out[0];
	stride_out[2] = _stride_out[2];
	stride_out[3] = _stride_out[3];
	// get unit stride input from LDS
	stride_in[0] = _lengths[0];
	stride_in[1] = 1;
	stride_in[2] = _lengths[2];
	stride_in[3] = _lengths[3];
	
	// flip dimensions and transform each column

	auto temp = lengths[0];
	lengths[0] = lengths[1];
	lengths[1] = temp;
	// Let the row transform finish before starting column transform
	__syncthreads();
	// declare input/output pointers for column transform
	gbIn = lds_data;
	real_type_t<T>* gbOutRe = _gbOutRe;
	real_type_t<T>* gbOutIm = _gbOutIm;
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 5);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/25));
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
	fwd_len125_device<T, SB_NONUNIT, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%25, (me/25)*125, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, nullptr, nullptr, load_cb_lds_bytes, store_cb_fn, store_cb_data);
	}
}

//Kernel configuration: number of threads per thread block: 250, maximum transforms: 32, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(250)
fft_back_op_2D_4_125( const T * __restrict__ twiddles, const size_t dim, const size_t *_lengths, const size_t *_stride_in, const size_t *_stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, T * __restrict__ _gbIn, real_type_t<T> * __restrict__ _gbOutRe, real_type_t<T> * __restrict__ _gbOutIm)
{

	__shared__ real_type_t<T> lds[500];
	__shared__ T lds_data[4*125];
	// use supplied input stride for row transform
	size_t stride_in[4];
	stride_in[0] = _stride_in[0];
	stride_in[1] = _stride_in[1];
	stride_in[2] = _stride_in[2];
	stride_in[3] = _stride_in[3];
	// set unit output stride, since we're writing to LDS
	size_t stride_out[4];
	stride_out[0] = 1;
	stride_out[1] = _lengths[0];
	stride_out[2] = _lengths[1];
	stride_out[3] = _lengths[2];
	// use supplied lengths for row transform
	size_t lengths[3];
	lengths[0] = _lengths[0];
	lengths[1] = _lengths[1];
	lengths[2] = _lengths[2];
	// declare input/output pointers
	T* gbIn = _gbIn;
	// write to LDS
	T* gbOut = lds_data;
	// transform each row
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 2);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/2));
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
	back_len4_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*4, gbIn, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, 0, nullptr, nullptr);
	}
	// we have two twiddle tables back to back in device
	// memory - move to the second table (if nonsquare)
	twiddles = twiddles + lengths[0];
	// write output to original out-of-place destination
	stride_out[0] = _stride_out[1];
	stride_out[1] = _stride_out[0];
	stride_out[2] = _stride_out[2];
	stride_out[3] = _stride_out[3];
	// get unit stride input from LDS
	stride_in[0] = _lengths[0];
	stride_in[1] = 1;
	stride_in[2] = _lengths[2];
	stride_in[3] = _lengths[3];
	
	// flip dimensions and transform each column

	auto temp = lengths[0];
	lengths[0] = lengths[1];
	lengths[1] = temp;
	// Let the row transform finish before starting column transform
	__syncthreads();
	// declare input/output pointers for column transform
	gbIn = lds_data;
	real_type_t<T>* gbOutRe = _gbOutRe;
	real_type_t<T>* gbOutIm = _gbOutIm;
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 5);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/25));
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
	back_len125_device<T, SB_NONUNIT, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%25, (me/25)*125, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, nullptr, nullptr, load_cb_lds_bytes, store_cb_fn, store_cb_data);
	}
}

//Kernel configuration: number of threads per thread block: 250, maximum transforms: 32, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(250)
fft_fwd_op_2D_4_125( const T * __restrict__ twiddles, const size_t dim, const size_t *_lengths, const size_t *_stride_in, const size_t *_stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ _gbInRe, real_type_t<T> * __restrict__ _gbInIm, T * __restrict__ _gbOut)
{

	__shared__ real_type_t<T> lds[500];
	__shared__ T lds_data[4*125];
	// use supplied input stride for row transform
	size_t stride_in[4];
	stride_in[0] = _stride_in[0];
	stride_in[1] = _stride_in[1];
	stride_in[2] = _stride_in[2];
	stride_in[3] = _stride_in[3];
	// set unit output stride, since we're writing to LDS
	size_t stride_out[4];
	stride_out[0] = 1;
	stride_out[1] = _lengths[0];
	stride_out[2] = _lengths[1];
	stride_out[3] = _lengths[2];
	// use supplied lengths for row transform
	size_t lengths[3];
	lengths[0] = _lengths[0];
	lengths[1] = _lengths[1];
	lengths[2] = _lengths[2];
	// declare input/output pointers
	real_type_t<T>* gbInRe = _gbInRe;
	real_type_t<T>* gbInIm = _gbInIm;
	// write to LDS
	T* gbOut = lds_data;
	// transform each row
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 2);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/2));
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
	fwd_len4_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*4, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, 0, nullptr, nullptr);
	}
	// we have two twiddle tables back to back in device
	// memory - move to the second table (if nonsquare)
	twiddles = twiddles + lengths[0];
	// write output to original out-of-place destination
	stride_out[0] = _stride_out[1];
	stride_out[1] = _stride_out[0];
	stride_out[2] = _stride_out[2];
	stride_out[3] = _stride_out[3];
	// get unit stride input from LDS
	stride_in[0] = _lengths[0];
	stride_in[1] = 1;
	stride_in[2] = _lengths[2];
	stride_in[3] = _lengths[3];
	
	// flip dimensions and transform each column

	auto temp = lengths[0];
	lengths[0] = lengths[1];
	lengths[1] = temp;
	// Let the row transform finish before starting column transform
	__syncthreads();
	// declare input/output pointers for column transform
	T* gbIn = lds_data;
	gbOut = _gbOut;
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 5);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/25));
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
	fwd_len125_device<T, SB_NONUNIT, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%25, (me/25)*125, gbIn, iOffset, gbOut, oOffset, lds, nullptr, nullptr, load_cb_lds_bytes, store_cb_fn, store_cb_data);
	}
}

//Kernel configuration: number of threads per thread block: 250, maximum transforms: 32, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(250)
fft_back_op_2D_4_125( const T * __restrict__ twiddles, const size_t dim, const size_t *_lengths, const size_t *_stride_in, const size_t *_stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ _gbInRe, real_type_t<T> * __restrict__ _gbInIm, T * __restrict__ _gbOut)
{

	__shared__ real_type_t<T> lds[500];
	__shared__ T lds_data[4*125];
	// use supplied input stride for row transform
	size_t stride_in[4];
	stride_in[0] = _stride_in[0];
	stride_in[1] = _stride_in[1];
	stride_in[2] = _stride_in[2];
	stride_in[3] = _stride_in[3];
	// set unit output stride, since we're writing to LDS
	size_t stride_out[4];
	stride_out[0] = 1;
	stride_out[1] = _lengths[0];
	stride_out[2] = _lengths[1];
	stride_out[3] = _lengths[2];
	// use supplied lengths for row transform
	size_t lengths[3];
	lengths[0] = _lengths[0];
	lengths[1] = _lengths[1];
	lengths[2] = _lengths[2];
	// declare input/output pointers
	real_type_t<T>* gbInRe = _gbInRe;
	real_type_t<T>* gbInIm = _gbInIm;
	// write to LDS
	T* gbOut = lds_data;
	// transform each row
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 2);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/2));
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
	back_len4_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*4, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, 0, nullptr, nullptr);
	}
	// we have two twiddle tables back to back in device
	// memory - move to the second table (if nonsquare)
	twiddles = twiddles + lengths[0];
	// write output to original out-of-place destination
	stride_out[0] = _stride_out[1];
	stride_out[1] = _stride_out[0];
	stride_out[2] = _stride_out[2];
	stride_out[3] = _stride_out[3];
	// get unit stride input from LDS
	stride_in[0] = _lengths[0];
	stride_in[1] = 1;
	stride_in[2] = _lengths[2];
	stride_in[3] = _lengths[3];
	
	// flip dimensions and transform each column

	auto temp = lengths[0];
	lengths[0] = lengths[1];
	lengths[1] = temp;
	// Let the row transform finish before starting column transform
	__syncthreads();
	// declare input/output pointers for column transform
	T* gbIn = lds_data;
	gbOut = _gbOut;
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 5);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/25));
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
	back_len125_device<T, SB_NONUNIT, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%25, (me/25)*125, gbIn, iOffset, gbOut, oOffset, lds, nullptr, nullptr, load_cb_lds_bytes, store_cb_fn, store_cb_data);
	}
}

//Kernel configuration: number of threads per thread block: 250, maximum transforms: 32, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(250)
fft_fwd_op_2D_4_125( const T * __restrict__ twiddles, const size_t dim, const size_t *_lengths, const size_t *_stride_in, const size_t *_stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ _gbInRe, real_type_t<T> * __restrict__ _gbInIm, real_type_t<T> * __restrict__ _gbOutRe, real_type_t<T> * __restrict__ _gbOutIm)
{

	__shared__ real_type_t<T> lds[500];
	__shared__ T lds_data[4*125];
	// use supplied input stride for row transform
	size_t stride_in[4];
	stride_in[0] = _stride_in[0];
	stride_in[1] = _stride_in[1];
	stride_in[2] = _stride_in[2];
	stride_in[3] = _stride_in[3];
	// set unit output stride, since we're writing to LDS
	size_t stride_out[4];
	stride_out[0] = 1;
	stride_out[1] = _lengths[0];
	stride_out[2] = _lengths[1];
	stride_out[3] = _lengths[2];
	// use supplied lengths for row transform
	size_t lengths[3];
	lengths[0] = _lengths[0];
	lengths[1] = _lengths[1];
	lengths[2] = _lengths[2];
	// declare input/output pointers
	real_type_t<T>* gbInRe = _gbInRe;
	real_type_t<T>* gbInIm = _gbInIm;
	// write to LDS
	T* gbOut = lds_data;
	// transform each row
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 2);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/2));
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
	fwd_len4_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*4, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, 0, nullptr, nullptr);
	}
	// we have two twiddle tables back to back in device
	// memory - move to the second table (if nonsquare)
	twiddles = twiddles + lengths[0];
	// write output to original out-of-place destination
	stride_out[0] = _stride_out[1];
	stride_out[1] = _stride_out[0];
	stride_out[2] = _stride_out[2];
	stride_out[3] = _stride_out[3];
	// get unit stride input from LDS
	stride_in[0] = _lengths[0];
	stride_in[1] = 1;
	stride_in[2] = _lengths[2];
	stride_in[3] = _lengths[3];
	
	// flip dimensions and transform each column

	auto temp = lengths[0];
	lengths[0] = lengths[1];
	lengths[1] = temp;
	// Let the row transform finish before starting column transform
	__syncthreads();
	// declare input/output pointers for column transform
	T* gbIn = lds_data;
	real_type_t<T>* gbOutRe = _gbOutRe;
	real_type_t<T>* gbOutIm = _gbOutIm;
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 5);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/25));
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
	fwd_len125_device<T, SB_NONUNIT, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%25, (me/25)*125, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, nullptr, nullptr, load_cb_lds_bytes, store_cb_fn, store_cb_data);
	}
}

//Kernel configuration: number of threads per thread block: 250, maximum transforms: 32, Passes: 2
template <typename T, StrideBin sb, EmbeddedType ebtype, CallbackType cbtype>
__global__ void
__launch_bounds__(250)
fft_back_op_2D_4_125( const T * __restrict__ twiddles, const size_t dim, const size_t *_lengths, const size_t *_stride_in, const size_t *_stride_out, const size_t batch_count, const unsigned int lds_padding, void* __restrict__ load_cb_fn, void* __restrict__ load_cb_data, uint32_t load_cb_lds_bytes, void* __restrict__ store_cb_fn, void* __restrict__ store_cb_data, real_type_t<T> * __restrict__ _gbInRe, real_type_t<T> * __restrict__ _gbInIm, real_type_t<T> * __restrict__ _gbOutRe, real_type_t<T> * __restrict__ _gbOutIm)
{

	__shared__ real_type_t<T> lds[500];
	__shared__ T lds_data[4*125];
	// use supplied input stride for row transform
	size_t stride_in[4];
	stride_in[0] = _stride_in[0];
	stride_in[1] = _stride_in[1];
	stride_in[2] = _stride_in[2];
	stride_in[3] = _stride_in[3];
	// set unit output stride, since we're writing to LDS
	size_t stride_out[4];
	stride_out[0] = 1;
	stride_out[1] = _lengths[0];
	stride_out[2] = _lengths[1];
	stride_out[3] = _lengths[2];
	// use supplied lengths for row transform
	size_t lengths[3];
	lengths[0] = _lengths[0];
	lengths[1] = _lengths[1];
	lengths[2] = _lengths[2];
	// declare input/output pointers
	real_type_t<T>* gbInRe = _gbInRe;
	real_type_t<T>* gbInIm = _gbInIm;
	// write to LDS
	T* gbOut = lds_data;
	// transform each row
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 2);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/2));
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
	back_len4_device<T, sb, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%2, (me/2)*4, gbInRe, gbInIm, iOffset, gbOut, oOffset, lds, load_cb_fn, load_cb_data, 0, nullptr, nullptr);
	}
	// we have two twiddle tables back to back in device
	// memory - move to the second table (if nonsquare)
	twiddles = twiddles + lengths[0];
	// write output to original out-of-place destination
	stride_out[0] = _stride_out[1];
	stride_out[1] = _stride_out[0];
	stride_out[2] = _stride_out[2];
	stride_out[3] = _stride_out[3];
	// get unit stride input from LDS
	stride_in[0] = _lengths[0];
	stride_in[1] = 1;
	stride_in[2] = _lengths[2];
	stride_in[3] = _lengths[3];
	
	// flip dimensions and transform each column

	auto temp = lengths[0];
	lengths[0] = lengths[1];
	lengths[1] = temp;
	// Let the row transform finish before starting column transform
	__syncthreads();
	// declare input/output pointers for column transform
	T* gbIn = lds_data;
	real_type_t<T>* gbOutRe = _gbOutRe;
	real_type_t<T>* gbOutIm = _gbOutIm;
	{
	unsigned int me = (unsigned int)hipThreadIdx_x;
	unsigned int batch = (unsigned int)hipBlockIdx_x;

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	
	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 5);
	unsigned int b = 0;

    // generator.kernel.hpp:1453
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/25));
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
	back_len125_device<T, SB_NONUNIT, cbtype>(twiddles, stride_in[0], stride_out[0],  rw, b, me%25, (me/25)*125, gbIn, iOffset, gbOutRe, gbOutIm, oOffset, lds, nullptr, nullptr, load_cb_lds_bytes, store_cb_fn, store_cb_data);
	}
}

