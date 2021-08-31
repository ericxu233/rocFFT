// Copyright (c) 2016 - present Advanced Micro Devices, Inc. All rights reserved.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "transpose.h"
#include "kernel_launch.h"
#include "rocfft_hip.h"
#include <CL/sycl.hpp>
#include <iostream>

namespace sycl = cl::sycl;

// chain of macros to iterate over transpose kernel template parameters, to
// set a function pointer 'kernel_func'
#define TRANSPOSE_KERNEL2_CALLBACK(TWL, DIR, ALL, STRIDE, DIAG, CALLBACK)                        \
    else if(twl == TWL && dir == DIR && all == ALL && unit_stride0 == STRIDE && diagonal == DIAG \
            && cbtype == CALLBACK) kernel_func                                                   \
        = transpose_kernel2<T,                                                                   \
                            TA,                                                                  \
                            TB,                                                                  \
                            TRANSPOSE_DIM_X,                                                     \
                            TRANSPOSE_DIM_Y,                                                     \
                            true,                                                                \
                            TWL,                                                                 \
                            DIR,                                                                 \
                            ALL,                                                                 \
                            STRIDE,                                                              \
                            DIAG,                                                                \
                            CALLBACK>;

#define TRANSPOSE_KERNEL2_DIAG(TWL, DIR, ALL, STRIDE, DIAG)                                \
    TRANSPOSE_KERNEL2_CALLBACK(TWL, DIR, ALL, STRIDE, DIAG, CallbackType::USER_LOAD_STORE) \
    TRANSPOSE_KERNEL2_CALLBACK(TWL, DIR, ALL, STRIDE, DIAG, CallbackType::NONE)

#define TRANSPOSE_KERNEL2_STRIDE(TWL, DIR, ALL, STRIDE)  \
    TRANSPOSE_KERNEL2_DIAG(TWL, DIR, ALL, STRIDE, false) \
    TRANSPOSE_KERNEL2_DIAG(TWL, DIR, ALL, STRIDE, true)

#define TRANSPOSE_KERNEL2_ALL(TWL, DIR, ALL)       \
    TRANSPOSE_KERNEL2_STRIDE(TWL, DIR, ALL, false) \
    TRANSPOSE_KERNEL2_STRIDE(TWL, DIR, ALL, true)

#define TRANSPOSE_KERNEL2_DIR(TWL, DIR)    \
    TRANSPOSE_KERNEL2_ALL(TWL, DIR, false) \
    TRANSPOSE_KERNEL2_ALL(TWL, DIR, true)

#define TRANSPOSE_KERNEL2_TWL(TWL) \
    TRANSPOSE_KERNEL2_DIR(TWL, 1)  \
    TRANSPOSE_KERNEL2_DIR(TWL, -1)

// every function pointer assignment is generated as an "else if", so
// start with "if(false)" to get the ball rolling
#define GET_TRANSPOSE_KERNEL2_FUNC() \
    if(false)                        \
        ;                            \
    TRANSPOSE_KERNEL2_TWL(0)         \
    TRANSPOSE_KERNEL2_TWL(1)         \
    TRANSPOSE_KERNEL2_TWL(2)         \
    TRANSPOSE_KERNEL2_TWL(3)         \
    TRANSPOSE_KERNEL2_TWL(4)

#define TRANSPOSE_KERNEL2_SCHEME_CALLBACK(ALL, STRIDE, DIAG, CALLBACK)                      \
    else if(all == ALL && unit_stride0 == STRIDE && diagonal == DIAG && cbtype == CALLBACK) \
        kernel_func                                                                         \
        = transpose_kernel2_scheme<T,                                                       \
                                   TA,                                                      \
                                   TB,                                                      \
                                   TRANSPOSE_DIM_X,                                         \
                                   TRANSPOSE_DIM_Y,                                         \
                                   ALL,                                                     \
                                   STRIDE,                                                  \
                                   DIAG,                                                    \
                                   CALLBACK>;

#define TRANSPOSE_KERNEL2_SCHEME_DIAG(ALL, STRIDE, DIAG)                                \
    TRANSPOSE_KERNEL2_SCHEME_CALLBACK(ALL, STRIDE, DIAG, CallbackType::USER_LOAD_STORE) \
    TRANSPOSE_KERNEL2_SCHEME_CALLBACK(ALL, STRIDE, DIAG, CallbackType::NONE)
#define TRANSPOSE_KERNEL2_SCHEME_STRIDE(ALL, STRIDE)  \
    TRANSPOSE_KERNEL2_SCHEME_DIAG(ALL, STRIDE, false) \
    TRANSPOSE_KERNEL2_SCHEME_DIAG(ALL, STRIDE, true)

#define TRANSPOSE_KERNEL2_SCHEME_ALL(ALL)       \
    TRANSPOSE_KERNEL2_SCHEME_STRIDE(ALL, false) \
    TRANSPOSE_KERNEL2_SCHEME_STRIDE(ALL, true)

// every function pointer assignment is generated as an "else if", so
// start with "if(false)" to get the ball rolling
#define GET_TRANSPOSE_KERNEL2_SCHEME_FUNC() \
    if(false)                               \
        ;                                   \
    TRANSPOSE_KERNEL2_SCHEME_ALL(false)     \
    TRANSPOSE_KERNEL2_SCHEME_ALL(true)

/// \brief FFT Transpose out-of-place API
/// \details transpose matrix A of size (m row by n cols) to matrix B (n row by m cols)
///    both A and B are in row major
///
/// @param[in]    m size_t.
/// @param[in]    n size_t.
/// @param[in]    A pointer storing batch_count of A matrix on the GPU.
/// @param[inout] B pointer storing batch_count of B matrix on the GPU.
/// @param[in]    count size_t number of matrices processed
template <typename T, typename TA, typename TB, int TRANSPOSE_DIM_X, int TRANSPOSE_DIM_Y>
rocfft_status rocfft_transpose_outofplace_template(size_t       m,
                                                   size_t       n,
                                                   const TA*    A,             //not sure if conversion requires
                                                   TB*          B,             //not sure if conversion requires
                                                   sycl::buffer<T, 1> twiddles_large,
                                                   size_t       count,
                                                   size_t*      lengths,
                                                   size_t*      stride_in,
                                                   size_t*      stride_out,
                                                   int          twl,
                                                   int          dir,
                                                   int          scheme,
                                                   bool         unit_stride0,
                                                   bool         diagonal,
                                                   size_t       ld_in,
                                                   size_t       ld_out,
                                                   sycl::queue  rocfft_queue,
                                                   CallbackType cbtype,
                                                   void* __restrict__ load_cb_fn,
                                                   void* __restrict__ load_cb_data,
                                                   uint32_t load_cb_lds_bytes,
                                                   void* __restrict__ store_cb_fn,
                                                   void* __restrict__ store_cb_data)
{
    /*
        @@ Changes start here, by Eric Xu
    */
    sycl::range<3> grid((n - 1) / TRANSPOSE_DIM_X + 1, ((m - 1) / TRANSPOSE_DIM_X + 1), count);
    sycl::range<3> threads(TRANSPOSE_DIM_X, TRANSPOSE_DIM_Y, 1);

    // working threads match problem sizes, no partial cases
    const bool all = (n % TRANSPOSE_DIM_X == 0) && (m % TRANSPOSE_DIM_X == 0);

    if(scheme == 0)
    {
        void (*kernel_func)(const TA*,
                            TB*,
                            sycl::buffer<T, 1>,
                            size_t*,
                            size_t*,
                            size_t*,
                            void* __restrict__,
                            void* __restrict__,
                            uint32_t,
                            void* __restrict__,
                            void* __restrict__)
            = nullptr;
        GET_TRANSPOSE_KERNEL2_FUNC();

        if(kernel_func)
        //where are the template arguments? are they processed by the macros from kernel_func? <<<>>>>
            rocFFTLaunchKernelGGL(kernel_func, 
                               sycl::range<3>(grid),
                               sycl::range<3>(threads),
                               0,
                               rocfft_queue,
                               A,
                               B,
                               twiddles_large,
                               lengths,
                               stride_in,
                               stride_out,
                               load_cb_fn,
                               load_cb_data,
                               load_cb_lds_bytes,
                               store_cb_fn,
                               store_cb_data);
        else
        {
            rocfft_cout << "scheme: " << scheme << std::endl;
            rocfft_cout << "twl: " << twl << std::endl;
            rocfft_cout << "dir: " << dir << std::endl;
            rocfft_cout << "all: " << all << std::endl;
            rocfft_cout << "diagonal: " << diagonal << std::endl;
        }
    }
    else
    {
        void (*kernel_func)(const TA*,
                            TB*,
                            sycl::buffer<T, 1>,
                            sycl::buffer<size_t, 1>*,
                            sycl::buffer<size_t, 1>*,
                            sycl::buffer<size_t, 1>*,
                            size_t,
                            size_t,
                            size_t,
                            size_t,
                            void* __restrict__,
                            void* __restrict__,
                            uint32_t,
                            void* __restrict__,
                            void* __restrict__)
            = nullptr;
        GET_TRANSPOSE_KERNEL2_SCHEME_FUNC();

        if(kernel_func)
        {
            rocFFTLaunchKernelGGL(kernel_func,
                               sycl::range<3>(grid),
                               sycl::range<3>(threads),
                               0,
                               rocfft_queue,
                               A,
                               B,
                               twiddles_large,
                               lengths,
                               stride_in,
                               stride_out,
                               ld_in,
                               ld_out,
                               m,
                               n,
                               load_cb_fn,
                               load_cb_data,
                               load_cb_lds_bytes,
                               store_cb_fn,
                               store_cb_data);
        }
        else
        {
            rocfft_cout << "scheme: " << scheme << std::endl;
            rocfft_cout << "twl: " << twl << std::endl;
            rocfft_cout << "dir: " << dir << std::endl;
            rocfft_cout << "all: " << all << std::endl;
            rocfft_cout << "diagonal: " << diagonal << std::endl;
        }
    }

    return rocfft_status_success;
}

void rocfft_internal_transpose_var2(const void* data_p, void* back_p)
{
    DeviceCallIn* data = (DeviceCallIn*)data_p;

    size_t m = data->node->length[1];
    size_t n = data->node->length[0];

    int scheme = 0;
    if(data->node->scheme == CS_KERNEL_TRANSPOSE_XY_Z)
    {
        scheme = 1;
        m      = data->node->length[2];
        n      = data->node->length[0] * data->node->length[1];
    }
    else if(data->node->scheme == CS_KERNEL_TRANSPOSE_Z_XY)
    {
        scheme = 2;
        m      = data->node->length[1] * data->node->length[2];
        n      = data->node->length[0];
    }

    size_t ld_in  = scheme == 1 ? data->node->inStride[2] : data->node->inStride[1];
    size_t ld_out = scheme == 1 ? data->node->outStride[1] : data->node->outStride[2];

    // TODO:
    //   - Might open this option to upstream
    //   - Check it for non-unit stride and other cases
    //   - It seems diagonal transpose is not friendly to XY_Z pow-of-2 cases. Need more investigation.
    bool diagonal = false;

    // NB:
    //    Need better arch comparison other than strings.
    //    Need to check it for Navi and other oncoming new Archs.
    std::string arch(data->deviceProp.gcnArchName); //Eric_critical: what should i do about this?
    if(arch.compare(0, 6, "gfx908") == 0 || arch.compare(0, 6, "gfx90a") == 0)
    {
        diagonal = (m % 256) == 0 && (data->node->outStride[1] % 256 == 0)
                   && (data->node->scheme != CS_KERNEL_TRANSPOSE_XY_Z);
    }
    else
    {
        diagonal = (m % 256) == 0 && (data->node->outStride[1] % 256 == 0);
    }

    // size_t ld_in = data->node->inStride[1];
    // size_t ld_out = data->node->outStride[1];

    // if (ld_in < m )
    //     return rocfft_status_invalid_dimensions;
    // else if (ld_out < n )
    //     return rocfft_status_invalid_dimensions;

    // if(m == 0 || n == 0 ) return rocfft_status_success;

    int twl = 0;

    if(data->node->large1D > (size_t)256 * 256 * 256 * 256)
        printf("large1D twiddle size too large error");
    else if(data->node->large1D > (size_t)256 * 256 * 256)
        twl = 4;
    else if(data->node->large1D > (size_t)256 * 256)
        twl = 3;
    // TODO- possibly using a smaller LargeTwdBase for transpose by large_twiddle_base
    else if(data->node->large1D > (size_t)256)
        twl = 2;
    else if(data->node->large1D > 0)
        twl = 1;
    else
        twl = 0;

    int dir = data->node->direction;

    size_t count = data->node->batch;

    size_t extraDimStart = 2;
    if(scheme != 0)
        extraDimStart = 3;

    sycl::queue rocfft_queue = data->rocfft_stream;

    bool unit_stride0
        = (data->node->inStride[0] == 1 && data->node->outStride[0] == 1) ? true : false;

    for(size_t i = extraDimStart; i < data->node->length.size(); i++)
        count *= data->node->length[i];

    CallbackType cbtype = data->get_callback_type();

    // double2 must use 32 otherwise exceed the shared memory (LDS) size

    // FIXME: push planar ptr on device in better way!!!
    if((data->node->inArrayType == rocfft_array_type_complex_planar
        || data->node->inArrayType == rocfft_array_type_hermitian_planar)
       && (data->node->outArrayType == rocfft_array_type_complex_interleaved
           || data->node->outArrayType == rocfft_array_type_hermitian_interleaved))
    {
        if(data->node->precision == rocfft_precision_single)
        {
            cmplx_float_planar in_planar;
            in_planar.R = (real_type_t<float2>*)data->bufIn[0];
            in_planar.I = (real_type_t<float2>*)data->bufIn[1];

            void* d_in_planar;
            //hipMalloc(&d_in_planar, sizeof(cmplx_float_planar));
            //hipMemcpy(d_in_planar, &in_planar, sizeof(cmplx_float_planar), hipMemcpyHostToDevice);
            d_in_planar = static_cast<void *>(sycl::malloc_device(sizeof(cmplx_float_planar), rocfft_queue));
            rocfft_queue.memcpy(&in_planar, d_in_planar, sizeof(cmplx_float_planar));

            //sycl::buffer<int, 1> d_in_planar(range<1>(sizeof(cmplx_float_planar))); //line that potentially replaces the previous three lines

            //pinned memory

            rocfft_transpose_outofplace_template<cmplx_float,
                                                 cmplx_float_planar,
                                                 cmplx_float,
                                                 64,
                                                 16>(
                m,
                n,
                (const cmplx_float_planar*)d_in_planar,
                (cmplx_float*)data->bufOut[0],
                data->node->twiddles_large.data(),
                count,
                data->node->devKernArg.data(),
                data->node->devKernArg.data() + 1 * KERN_ARGS_ARRAY_WIDTH,
                data->node->devKernArg.data() + 2 * KERN_ARGS_ARRAY_WIDTH,
                twl,
                dir,
                scheme,
                unit_stride0,
                diagonal,
                ld_in,
                ld_out,
                rocfft_queue,
                cbtype,
                data->callbacks.load_cb_fn,
                data->callbacks.load_cb_data,
                data->callbacks.load_cb_lds_bytes,
                data->callbacks.store_cb_fn,
                data->callbacks.store_cb_data);

            sycl::free(d_in_planar, rocfft_queue); // needs to free this if implementation to sycl is made
                                  // sycl buffers are not required to be freed since the scoping will automatically destruct it
        }
        else
        {
            cmplx_double_planar in_planar;
            in_planar.R = (real_type_t<double2>*)data->bufIn[0];
            in_planar.I = (real_type_t<double2>*)data->bufIn[1];

            void* d_in_planar;
            //hipMalloc(&d_in_planar, sizeof(cmplx_double_planar));
            //hipMemcpy(d_in_planar, &in_planar, sizeof(cmplx_double_planar), hipMemcpyHostToDevice);
            d_in_planar = static_cast<void *>(sycl::malloc_device(sizeof(cmplx_double_planar), rocfft_queue));
            rocfft_queue.memcpy(&in_planar, d_in_planar, sizeof(cmplx_double_planar));

            //sycl::buffer<int, 1> d_in_planar(range<1>(sizeof(cmplx_double_planar))); //line that potentially replaces the previous three lines

            rocfft_transpose_outofplace_template<cmplx_double,
                                                 cmplx_double_planar,
                                                 cmplx_double,
                                                 32,
                                                 32>(
                m,
                n,
                (const cmplx_double_planar*)d_in_planar,
                (double2*)data->bufOut[0],
                data->node->twiddles_large.data(),
                count,
                data->node->devKernArg.data(),
                data->node->devKernArg.data() + 1 * KERN_ARGS_ARRAY_WIDTH,
                data->node->devKernArg.data() + 2 * KERN_ARGS_ARRAY_WIDTH,
                twl,
                dir,
                scheme,
                unit_stride0,
                diagonal,
                ld_in,
                ld_out,
                rocfft_queue,
                cbtype,
                data->callbacks.load_cb_fn,
                data->callbacks.load_cb_data,
                data->callbacks.load_cb_lds_bytes,
                data->callbacks.store_cb_fn,
                data->callbacks.store_cb_data);

            sycl::free(d_in_planar, rocfft_queue);
        }
    }
    else if((data->node->inArrayType == rocfft_array_type_complex_interleaved
             || data->node->inArrayType == rocfft_array_type_hermitian_interleaved)
            && (data->node->outArrayType == rocfft_array_type_complex_planar
                || data->node->outArrayType == rocfft_array_type_hermitian_planar))
    {
        if(data->node->precision == rocfft_precision_single)
        {
            cmplx_float_planar out_planar;
            out_planar.R = (real_type_t<float2>*)data->bufOut[0];
            out_planar.I = (real_type_t<float2>*)data->bufOut[1];

            void* d_out_planar;
            //hipMalloc(&d_out_planar, sizeof(cmplx_float_planar));
            //hipMemcpy(d_out_planar, &out_planar, sizeof(cmplx_float_planar), hipMemcpyHostToDevice);
            d_out_planar = static_cast<void *>(sycl::malloc_device(sizeof(cmplx_float_planar), rocfft_queue));
            rocfft_queue.memcpy(&out_planar, d_out_planar, sizeof(cmplx_float_planar));
            

            //sycl::buffer<int, 1> d_out_planar(range<1>(sizeof(cmplx_float_planar))); //line that potentially replaces the previous three lines

            rocfft_transpose_outofplace_template<cmplx_float,
                                                 cmplx_float,
                                                 cmplx_float_planar,
                                                 64,
                                                 16>(
                m,
                n,
                (const cmplx_float*)data->bufIn[0],
                (cmplx_float_planar*)d_out_planar,
                data->node->twiddles_large.data(),
                count,
                data->node->devKernArg.data(),
                data->node->devKernArg.data() + 1 * KERN_ARGS_ARRAY_WIDTH,
                data->node->devKernArg.data() + 2 * KERN_ARGS_ARRAY_WIDTH,
                twl,
                dir,
                scheme,
                unit_stride0,
                diagonal,
                ld_in,
                ld_out,
                rocfft_queue,
                cbtype,
                data->callbacks.load_cb_fn,
                data->callbacks.load_cb_data,
                data->callbacks.load_cb_lds_bytes,
                data->callbacks.store_cb_fn,
                data->callbacks.store_cb_data);

            sycl::free(d_out_planar, rocfft_queue);
        }
        else
        {
            cmplx_double_planar out_planar;
            out_planar.R = (real_type_t<double2>*)data->bufOut[0];
            out_planar.I = (real_type_t<double2>*)data->bufOut[1];

            void* d_out_planar;
            //hipMalloc(&d_out_planar, sizeof(cmplx_double_planar));
            //hipMemcpy(d_out_planar, &out_planar, sizeof(cmplx_double_planar), hipMemcpyHostToDevice);
            d_out_planar = static_cast<void *>(sycl::malloc_device(sizeof(cmplx_double_planar), rocfft_queue));
            rocfft_queue.memcpy(&out_planar, d_out_planar, sizeof(cmplx_double_planar));

            rocfft_transpose_outofplace_template<cmplx_double,
                                                 cmplx_double,
                                                 cmplx_double_planar,
                                                 32,
                                                 32>(
                m,
                n,
                (const cmplx_double*)data->bufIn[0],
                (cmplx_double_planar*)d_out_planar,
                data->node->twiddles_large.data(),
                count,
                data->node->devKernArg.data(),
                data->node->devKernArg.data() + 1 * KERN_ARGS_ARRAY_WIDTH,
                data->node->devKernArg.data() + 2 * KERN_ARGS_ARRAY_WIDTH,
                twl,
                dir,
                scheme,
                unit_stride0,
                diagonal,
                ld_in,
                ld_out,
                rocfft_queue,
                cbtype,
                data->callbacks.load_cb_fn,
                data->callbacks.load_cb_data,
                data->callbacks.load_cb_lds_bytes,
                data->callbacks.store_cb_fn,
                data->callbacks.store_cb_data);

            sycl::free(d_out_planar, rocfft_queue);
        }
    }
    else if((data->node->inArrayType == rocfft_array_type_complex_planar
             || data->node->inArrayType == rocfft_array_type_hermitian_planar)
            && (data->node->outArrayType == rocfft_array_type_complex_planar
                || data->node->outArrayType == rocfft_array_type_hermitian_planar))
    {
        if(data->node->precision == rocfft_precision_single)
        {
            cmplx_float_planar in_planar;
            in_planar.R = (real_type_t<float2>*)data->bufIn[0];
            in_planar.I = (real_type_t<float2>*)data->bufIn[1];
            cmplx_float_planar out_planar;
            out_planar.R = (real_type_t<float2>*)data->bufOut[0];
            out_planar.I = (real_type_t<float2>*)data->bufOut[1];

            void* d_in_planar;
            //hipMalloc(&d_in_planar, sizeof(cmplx_float_planar));
            //hipMemcpy(d_in_planar, &in_planar, sizeof(cmplx_float_planar), hipMemcpyHostToDevice);
            d_in_planar = static_cast<void *>(sycl::malloc_device(sizeof(cmplx_float_planar), rocfft_queue));
            rocfft_queue.memcpy(&out_planar, d_in_planar, sizeof(cmplx_float_planar));


            //sycl::buffer<int, 1> d_in_planar(range<1>(sizeof(cmplx_float_planar))); //line that potentially replaces the previous three lines

            void* d_out_planar;
            //hipMalloc(&d_out_planar, sizeof(cmplx_float_planar));
            //hipMemcpy(d_out_planar, &out_planar, sizeof(cmplx_float_planar), hipMemcpyHostToDevice);
            d_out_planar = static_cast<void *>(sycl::malloc_device(sizeof(cmplx_float_planar), rocfft_queue));
            rocfft_queue.memcpy(&out_planar, d_out_planar, sizeof(cmplx_float_planar));

            //sycl::buffer<int, 1> d_out_planar(range<1>(sizeof(cmplx_float_planar))); //line that potentially replaces the previous three lines

            rocfft_transpose_outofplace_template<cmplx_float,
                                                 cmplx_float_planar,
                                                 cmplx_float_planar,
                                                 64,
                                                 16>(
                m,
                n,
                (const cmplx_float_planar*)d_in_planar,
                (cmplx_float_planar*)d_out_planar,
                data->node->twiddles_large.data(),
                count,
                data->node->devKernArg.data(),
                data->node->devKernArg.data() + 1 * KERN_ARGS_ARRAY_WIDTH,
                data->node->devKernArg.data() + 2 * KERN_ARGS_ARRAY_WIDTH,
                twl,
                dir,
                scheme,
                unit_stride0,
                diagonal,
                ld_in,
                ld_out,
                rocfft_queue,
                cbtype,
                data->callbacks.load_cb_fn,
                data->callbacks.load_cb_data,
                data->callbacks.load_cb_lds_bytes,
                data->callbacks.store_cb_fn,
                data->callbacks.store_cb_data);

            sycl::free(d_in_planar, rocfft_queue);
            sycl::free(d_out_planar, rocfft_queue);
        }
        else
        {
            cmplx_double_planar in_planar;
            in_planar.R = (real_type_t<double2>*)data->bufIn[0];
            in_planar.I = (real_type_t<double2>*)data->bufIn[1];
            cmplx_double_planar out_planar;
            out_planar.R = (real_type_t<double2>*)data->bufOut[0];
            out_planar.I = (real_type_t<double2>*)data->bufOut[1];

            void* d_in_planar;
            //hipMalloc(&d_in_planar, sizeof(cmplx_double_planar));
            //hipMemcpy(d_in_planar, &in_planar, sizeof(cmplx_double_planar), hipMemcpyHostToDevice);
            d_in_planar = static_cast<void *>(sycl::malloc_device(sizeof(cmplx_double_planar), rocfft_queue));
            rocfft_queue.memcpy(&out_planar, d_in_planar, sizeof(cmplx_double_planar));
            //sycl::buffer<int, 1> d_in_planar(range<1>(sizeof(cmplx_double_planar))); //line that potentially replaces the previous three lines

            void* d_out_planar;
            //hipMalloc(&d_out_planar, sizeof(cmplx_double_planar));
            //hipMemcpy(d_out_planar, &out_planar, sizeof(cmplx_double_planar), hipMemcpyHostToDevice);
            d_out_planar = static_cast<void *>(sycl::malloc_device(sizeof(cmplx_double_planar), rocfft_queue));
            rocfft_queue.memcpy(&out_planar, d_out_planar, sizeof(cmplx_double_planar));
            //sycl::buffer<int, 1> d_out_planar(range<1>(sizeof(cmplx_double_planar))); //line that potentially replaces the previous three lines

            rocfft_transpose_outofplace_template<cmplx_double,
                                                 cmplx_double_planar,
                                                 cmplx_double_planar,
                                                 32,
                                                 32>(
                m,
                n,
                (const cmplx_double_planar*)d_in_planar,
                (cmplx_double_planar*)d_out_planar,
                data->node->twiddles_large.data(),
                count,
                data->node->devKernArg.data(),
                data->node->devKernArg.data() + 1 * KERN_ARGS_ARRAY_WIDTH,
                data->node->devKernArg.data() + 2 * KERN_ARGS_ARRAY_WIDTH,
                twl,
                dir,
                scheme,
                unit_stride0,
                diagonal,
                ld_in,
                ld_out,
                rocfft_queue,
                cbtype,
                data->callbacks.load_cb_fn,
                data->callbacks.load_cb_data,
                data->callbacks.load_cb_lds_bytes,
                data->callbacks.store_cb_fn,
                data->callbacks.store_cb_data);

            sycl::free(d_in_planar, rocfft_queue);
            sycl::free(d_out_planar, rocfft_queue);
        }
    }
    else
    {
        //FIXME:
        //  there are more cases than
        //      if(data->node->inArrayType == rocfft_array_type_complex_interleaved
        //      && data->node->outArrayType == rocfft_array_type_complex_interleaved)
        //  fall into this default case which might to correct
        if(data->node->precision == rocfft_precision_single)
            rocfft_transpose_outofplace_template<cmplx_float, cmplx_float, cmplx_float, 64, 16>(
                m,
                n,
                (const cmplx_float*)data->bufIn[0],
                (cmplx_float*)data->bufOut[0],
                data->node->twiddles_large.data(),
                count,
                data->node->devKernArg.data(),
                data->node->devKernArg.data() + 1 * KERN_ARGS_ARRAY_WIDTH,
                data->node->devKernArg.data() + 2 * KERN_ARGS_ARRAY_WIDTH,
                twl,
                dir,
                scheme,
                unit_stride0,
                diagonal,
                ld_in,
                ld_out,
                rocfft_queue,
                cbtype,
                data->callbacks.load_cb_fn,
                data->callbacks.load_cb_data,
                data->callbacks.load_cb_lds_bytes,
                data->callbacks.store_cb_fn,
                data->callbacks.store_cb_data);
        else
            rocfft_transpose_outofplace_template<cmplx_double, cmplx_double, cmplx_double, 32, 32>(
                m,
                n,
                (const cmplx_double*)data->bufIn[0],
                (cmplx_double*)data->bufOut[0],
                data->node->twiddles_large.data(),
                count,
                data->node->devKernArg.data(),
                data->node->devKernArg.data() + 1 * KERN_ARGS_ARRAY_WIDTH,
                data->node->devKernArg.data() + 2 * KERN_ARGS_ARRAY_WIDTH,
                twl,
                dir,
                scheme,
                unit_stride0,
                diagonal,
                ld_in,
                ld_out,
                rocfft_queue,
                cbtype,
                data->callbacks.load_cb_fn,
                data->callbacks.load_cb_data,
                data->callbacks.load_cb_lds_bytes,
                data->callbacks.store_cb_fn,
                data->callbacks.store_cb_data);
    }
}
