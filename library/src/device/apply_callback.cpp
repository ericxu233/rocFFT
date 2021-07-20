#include "kernel_launch.h"
#include "kernels/callback.h"

static const size_t APPLY_REAL_CALLBACK_THREADS = 64;
template <typename Treal>
__global__ void __launch_bounds__(APPLY_REAL_CALLBACK_THREADS)
    apply_real_callback_kernel(const size_t input_size,
                               const size_t dist1D,
                               Treal* __restrict__ input0,
                               const size_t dist,
                               void* __restrict__ load_cb_fn,
                               void* __restrict__ load_cb_data,
                               uint32_t load_cb_lds_bytes,
                               void* __restrict__ store_cb_fn,
                               void* __restrict__ store_cb_data)
{
    const size_t tid = hipBlockIdx_x * hipBlockDim_x + hipThreadIdx_x;
    if(tid < input_size)
    {
        auto load_cb  = get_load_cb<Treal, CallbackType::USER_LOAD_STORE>(load_cb_fn);
        auto store_cb = get_store_cb<Treal, CallbackType::USER_LOAD_STORE>(store_cb_fn);

        // blockIdx.y gives the multi-dimensional offset
        // blockIdx.z gives the batch offset
        const auto idx = blockIdx.y * dist1D + blockIdx.z * dist;

        auto elem = load_cb(input0, idx + tid, load_cb_data, nullptr);
        store_cb(input0, idx + tid, elem, store_cb_data, nullptr);
    }
}

void apply_real_callback(const void* data_p, void* back)
{
    DeviceCallIn* data = (DeviceCallIn*)data_p;

    size_t input_size = data->node->length[0];

    size_t input_distance = data->node->iDist;

    size_t input_stride
        = (data->node->length.size() > 1) ? data->node->inStride[1] : input_distance;

    void* input_buffer = data->bufIn[0];

    size_t batch          = data->node->batch;
    size_t high_dimension = 1;
    if(data->node->length.size() > 1)
    {
        for(int i = 1; i < data->node->length.size(); i++)
        {
            high_dimension *= data->node->length[i];
        }
    }
    rocfft_precision precision = data->node->precision;

    size_t blocks = (input_size - 1) / APPLY_REAL_CALLBACK_THREADS + 1;

    dim3 grid(blocks, high_dimension, batch);
    dim3 threads(APPLY_REAL_CALLBACK_THREADS, 1, 1);

    switch(precision)
    {
    case rocfft_precision_single:
        rocFFTLaunchKernelGGL(HIP_KERNEL_NAME(apply_real_callback_kernel<float>),
                           grid,
                           threads,
                           0,
                           data->rocfft_stream,
                           input_size,
                           input_stride,
                           static_cast<float*>(input_buffer),
                           input_distance,
                           data->callbacks.load_cb_fn,
                           data->callbacks.load_cb_data,
                           data->callbacks.load_cb_lds_bytes,
                           data->callbacks.store_cb_fn,
                           data->callbacks.store_cb_data);
        break;
    case rocfft_precision_double:
        rocFFTLaunchKernelGGL(HIP_KERNEL_NAME(apply_real_callback_kernel<double>),
                           grid,
                           threads,
                           0,
                           data->rocfft_stream,
                           input_size,
                           input_stride,
                           static_cast<double*>(input_buffer),
                           input_distance,
                           data->callbacks.load_cb_fn,
                           data->callbacks.load_cb_data,
                           data->callbacks.load_cb_lds_bytes,
                           data->callbacks.store_cb_fn,
                           data->callbacks.store_cb_data);
        break;
    }
}
