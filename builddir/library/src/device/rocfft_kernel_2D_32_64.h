#pragma once
#include "rocfft_kernel_32.h"
#include "rocfft_kernel_64.h"

////////////////////////////////////////Global kernels
//Kernel configuration: number of threads per thread block: 64, maximum transforms: 16, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(256)
fft_fwd_ip_2D_32_64( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *_lengths_GB, cl::sycl::buffer<size_t, 1> *_stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> _gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _lengths = _lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_in = _stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _gb = _gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2048), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_2D_32_64>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 8);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/4));
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
	fwd_len32_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*32, lwbIn, lwbOut, lds);
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 4);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/16));
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
	fwd_len64_device<T, SB_NONUNIT>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*64, lwbIn, lwbOut, lds);
	}
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 16, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(256)
fft_back_ip_2D_32_64( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *_lengths_GB, cl::sycl::buffer<size_t, 1> *_stride_in_GB, const size_t batch_count, cl::sycl::buffer<T, 1> _gb_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _lengths = _lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_in = _stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _gb = _gb_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2048), cgh);
	cgh.parallel_for<class kern_fft_back_ip_2D_32_64>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 8);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/4));
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
	back_len32_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*32, lwbIn, lwbOut, lds);
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 4);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/16));
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
	back_len64_device<T, SB_NONUNIT>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*64, lwbIn, lwbOut, lds);
	}
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 16, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(256)
fft_fwd_ip_2D_32_64( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *_lengths_GB, cl::sycl::buffer<size_t, 1> *_stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> _gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> _gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _lengths = _lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_in = _stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _gbRe = _gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbIm = _gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2048), cgh);
	cgh.parallel_for<class kern_fft_fwd_ip_2D_32_64>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	real_type_t<T> *lwbInRe;
	real_type_t<T> *lwbInIm;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 8);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/4));
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
	fwd_len32_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*32, lwbInRe, lwbInIm, lwbOut, lds);
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	real_type_t<T> *lwbOutRe;
	real_type_t<T> *lwbOutIm;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 4);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/16));
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
	fwd_len64_device<T, SB_NONUNIT>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*64, lwbIn, lwbOutRe, lwbOutIm, lds);
	}
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 16, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(256)
fft_back_ip_2D_32_64( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *_lengths_GB, cl::sycl::buffer<size_t, 1> *_stride_in_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> _gbRe_GB, cl::sycl::buffer<real_type_t<T>, 1> _gbIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _lengths = _lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_in = _stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _gbRe = _gbRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbIm = _gbIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2048), cgh);
	cgh.parallel_for<class kern_fft_back_ip_2D_32_64>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	real_type_t<T> *lwbInRe;
	real_type_t<T> *lwbInIm;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 8);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/4));
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
	back_len32_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*32, lwbInRe, lwbInIm, lwbOut, lds);
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	real_type_t<T> *lwbOutRe;
	real_type_t<T> *lwbOutIm;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 4);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/16));
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
	back_len64_device<T, SB_NONUNIT>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*64, lwbIn, lwbOutRe, lwbOutIm, lds);
	}
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 16, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(256)
fft_fwd_op_2D_32_64( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *_lengths_GB, cl::sycl::buffer<size_t, 1> *_stride_in_GB, cl::sycl::buffer<size_t, 1> *_stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> _gbIn_GB, cl::sycl::buffer<T, 1> _gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _lengths = _lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_in = _stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_out = _stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _gbIn = _gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbOut = _gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2048), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_2D_32_64>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 8);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/4));
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
	fwd_len32_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*32, lwbIn, lwbOut, lds);
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 4);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/16));
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
	fwd_len64_device<T, SB_NONUNIT>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*64, lwbIn, lwbOut, lds);
	}
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 16, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(256)
fft_back_op_2D_32_64( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *_lengths_GB, cl::sycl::buffer<size_t, 1> *_stride_in_GB, cl::sycl::buffer<size_t, 1> *_stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> _gbIn_GB, cl::sycl::buffer<T, 1> _gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _lengths = _lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_in = _stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_out = _stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _gbIn = _gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbOut = _gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2048), cgh);
	cgh.parallel_for<class kern_fft_back_op_2D_32_64>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 8);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/4));
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
	back_len32_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*32, lwbIn, lwbOut, lds);
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 4);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/16));
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
	back_len64_device<T, SB_NONUNIT>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*64, lwbIn, lwbOut, lds);
	}
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 16, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(256)
fft_fwd_op_2D_32_64( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *_lengths_GB, cl::sycl::buffer<size_t, 1> *_stride_in_GB, cl::sycl::buffer<size_t, 1> *_stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> _gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> _gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> _gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _lengths = _lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_in = _stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_out = _stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _gbIn = _gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbOutRe = _gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbOutIm = _gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2048), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_2D_32_64>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 8);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/4));
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
	fwd_len32_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*32, lwbIn, lwbOut, lds);
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	real_type_t<T> *lwbOutRe;
	real_type_t<T> *lwbOutIm;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 4);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/16));
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
	fwd_len64_device<T, SB_NONUNIT>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*64, lwbIn, lwbOutRe, lwbOutIm, lds);
	}
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 16, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(256)
fft_back_op_2D_32_64( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *_lengths_GB, cl::sycl::buffer<size_t, 1> *_stride_in_GB, cl::sycl::buffer<size_t, 1> *_stride_out_GB, const size_t batch_count, cl::sycl::buffer<T, 1> _gbIn_GB, cl::sycl::buffer<real_type_t<T>, 1> _gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> _gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _lengths = _lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_in = _stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_out = _stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _gbIn = _gbIn_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbOutRe = _gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbOutIm = _gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2048), cgh);
	cgh.parallel_for<class kern_fft_back_op_2D_32_64>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 8);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/4));
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
	back_len32_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*32, lwbIn, lwbOut, lds);
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	real_type_t<T> *lwbOutRe;
	real_type_t<T> *lwbOutIm;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 4);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/16));
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
	back_len64_device<T, SB_NONUNIT>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*64, lwbIn, lwbOutRe, lwbOutIm, lds);
	}
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 16, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(256)
fft_fwd_op_2D_32_64( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *_lengths_GB, cl::sycl::buffer<size_t, 1> *_stride_in_GB, cl::sycl::buffer<size_t, 1> *_stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> _gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> _gbInIm_GB, cl::sycl::buffer<T, 1> _gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _lengths = _lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_in = _stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_out = _stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _gbInRe = _gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbInIm = _gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbOut = _gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2048), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_2D_32_64>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	real_type_t<T> *lwbInRe;
	real_type_t<T> *lwbInIm;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 8);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/4));
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
	fwd_len32_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*32, lwbInRe, lwbInIm, lwbOut, lds);
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 4);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/16));
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
	fwd_len64_device<T, SB_NONUNIT>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*64, lwbIn, lwbOut, lds);
	}
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 16, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(256)
fft_back_op_2D_32_64( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *_lengths_GB, cl::sycl::buffer<size_t, 1> *_stride_in_GB, cl::sycl::buffer<size_t, 1> *_stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> _gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> _gbInIm_GB, cl::sycl::buffer<T, 1> _gbOut_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _lengths = _lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_in = _stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_out = _stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _gbInRe = _gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbInIm = _gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbOut = _gbOut_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2048), cgh);
	cgh.parallel_for<class kern_fft_back_op_2D_32_64>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	real_type_t<T> *lwbInRe;
	real_type_t<T> *lwbInIm;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 8);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/4));
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
	back_len32_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*32, lwbInRe, lwbInIm, lwbOut, lds);
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 4);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/16));
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
	back_len64_device<T, SB_NONUNIT>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*64, lwbIn, lwbOut, lds);
	}
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 16, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(256)
fft_fwd_op_2D_32_64( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *_lengths_GB, cl::sycl::buffer<size_t, 1> *_stride_in_GB, cl::sycl::buffer<size_t, 1> *_stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> _gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> _gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> _gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> _gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _lengths = _lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_in = _stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_out = _stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _gbInRe = _gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbInIm = _gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbOutRe = _gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbOutIm = _gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2048), cgh);
	cgh.parallel_for<class kern_fft_fwd_op_2D_32_64>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	real_type_t<T> *lwbInRe;
	real_type_t<T> *lwbInIm;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 8);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/4));
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
	fwd_len32_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*32, lwbInRe, lwbInIm, lwbOut, lds);
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	real_type_t<T> *lwbOutRe;
	real_type_t<T> *lwbOutIm;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 4);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/16));
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
	fwd_len64_device<T, SB_NONUNIT>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*64, lwbIn, lwbOutRe, lwbOutIm, lds);
	}
});
});
}

//Kernel configuration: number of threads per thread block: 64, maximum transforms: 16, Passes: 2
template <typename T, StrideBin sb>
__global__ void
__launch_bounds__(256)
fft_back_op_2D_32_64( cl::sycl::range<3> blocks, cl::sycl::range<3> threads, cl::sycl::queue rocfftQueue, cl::sycl::buffer<T, 1> twiddles_GB, const size_t dim, cl::sycl::buffer<size_t, 1> *_lengths_GB, cl::sycl::buffer<size_t, 1> *_stride_in_GB, cl::sycl::buffer<size_t, 1> *_stride_out_GB, const size_t batch_count, cl::sycl::buffer<real_type_t<T>, 1> _gbInRe_GB, cl::sycl::buffer<real_type_t<T>, 1> _gbInIm_GB, cl::sycl::buffer<real_type_t<T>, 1> _gbOutRe_GB, cl::sycl::buffer<real_type_t<T>, 1> _gbOutIm_GB)
{
rocfftQueue.submit([&](cl::sycl::handler &cgh)
{
	auto twiddles = twiddles_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _lengths = _lengths_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_in = _stride_in_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _stride_out = _stride_out_GB.get_access<cl::sycl::access::mode::read>(cgh);
	auto _gbInRe = _gbInRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbInIm = _gbInIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbOutRe = _gbOutRe_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	auto _gbOutIm = _gbOutIm_GB.get_access<cl::sycl::access::mode::read_write>(cgh);
	cl::sycl::accessor<real_type_t<T>,1,sycl::access::mode::read_write,sycl::access::target::local> lds(cl::sycl::range<1>(2048), cgh);
	cgh.parallel_for<class kern_fft_back_op_2D_32_64>(
	                   cl::sycl::nd_range<3>(blocks, threads),
	                   [=](cl::sycl::nd_item<3> wItem)
{
//// Print kernel code here (after function prototype)
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	real_type_t<T> *lwbInRe;
	real_type_t<T> *lwbInIm;
	T *lwbOut;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 8);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// row transform writes to LDS, so respect non-unit strides for input
	// and assume unit stride for output
	iOffset = dim == 2 ?
		batch * _stride_in[2] :
		batch / lengths[2] * _stride_in[3] + batch % lengths[2] * _stride_in[2];
	size_t counter_mod = (batch*0 + (me/4));
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
	back_len32_device<T, sb>(twiddles, stride_in[0], stride_out[0],  rw, b, me%4, (me/4)*32, lwbInRe, lwbInIm, lwbOut, lds);
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
	unsigned int me = (unsigned int)wItem.get_local_id(0);
	unsigned int batch = (unsigned int)wItem.get_group(0);

	unsigned int iOffset = 0;
	unsigned int oOffset = 0;
	T *lwbIn;
	real_type_t<T> *lwbOutRe;
	real_type_t<T> *lwbOutIm;

	// set rw for enough threads to cover total number of 2D elements
	unsigned int rw = me < (lengths[0] * lengths[1] / 4);
	unsigned int b = 0;

    // generator.kernel.hpp:1581
	// col transform reads from LDS, so respect non-unit strides for output
	// and assume unit stride for input
	oOffset = dim == 2 ?
		batch * stride_out[2] :
		batch / lengths[2] * stride_out[3] + batch % lengths[2] * stride_out[2];
	size_t counter_mod = (batch*0 + (me/16));
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
	back_len64_device<T, SB_NONUNIT>(twiddles, stride_in[0], stride_out[0],  rw, b, me%16, (me/16)*64, lwbIn, lwbOutRe, lwbOutIm, lds);
	}
});
});
}

