#include <CL/sycl.hpp>

namespace sycl = cl::sycl;

template <typename T,
          typename T_I,
          typename T_O,
          size_t       DIM_X,
          size_t       DIM_Y,
          bool         ALL,
          bool         UNIT_STRIDE_0,
          CallbackType cbtype>
void transpose_tile_device_scheme(const T_I*   input,
                                             T_O*         output,
                                             size_t       in_offset,
                                             size_t       out_offset,
                                             const size_t m,
                                             const size_t n,
                                             size_t       ld_in,
                                             size_t       ld_out,
                                             size_t       stride_0_in,
                                             size_t       stride_0_out,
                                             void* __restrict__ load_cb_fn,
                                             void* __restrict__ load_cb_data,
                                             uint32_t load_cb_lds_bytes,
                                             void* __restrict__ store_cb_fn,
                                             void* __restrict__ store_cb_data,
                                             sycl::nd_item<3> wItem,
                                             sycl::accessor <T, 2, sycl::access::mode::read_write, sycl::access::target::local> shared)
{
    //__shared__ T shared[DIM_X][DIM_X];

    size_t tid = wItem.get_local_id(0)/*hipThreadIdx_x*/ + wItem.get_local_id(1)/*hipThreadIdx_y*/ * wItem.get_group_range(0)/*hipBlockDim_x*/;
    size_t tx1 = tid % DIM_X;
    size_t ty1 = tid / DIM_X;
    // what is program unroll?
    if(ALL)
    {
#pragma unroll
        for(int i = 0; i < DIM_X; i += DIM_Y)
        {
            T tmp;
            if(UNIT_STRIDE_0)
            {
                //tmp = Handler<T_I, cbtype>::read(
                    //input, in_offset + tx1 + (ty1 + i) * ld_in, load_cb_fn, load_cb_data);
            }
            else
            {
                // tmp = Handler<T_I, cbtype>::read(input,
                //                                  in_offset + tx1 * stride_0_in + (ty1 + i) * ld_in,
                //                                  load_cb_fn,
                //                                  load_cb_data);
            }
            shared[tx1][ty1 + i] = tmp; // the transpose taking place here
        }

        wItem.barrier(sycl::access::fence_space::local_space);// __syncthreads();
        T val[DIM_X / DIM_Y];
#pragma unroll
        for(int i = 0, j = 0; i < DIM_X; i += DIM_Y, j++)
        {
            val[j] = shared[ty1 + i][tx1];
        }
#pragma unroll
        for(int i = 0, j = 0; i < DIM_X; i += DIM_Y, j++)
        {
            // reconfigure the threads
            if(UNIT_STRIDE_0)
            {
                // Handler<T_O, cbtype>::write(output,
                //                             out_offset + tx1 + (i + ty1) * ld_out,
                //                             val[j],
                //                             store_cb_fn,
                //                             store_cb_data);
            }
            else
            {
                // Handler<T_O, cbtype>::write(output,
                //                             out_offset + tx1 * stride_0_out + (i + ty1) * ld_out,
                //                             val[j],
                //                             store_cb_fn,
                //                             store_cb_data);
            }
        }
    }
    else
    {
        T val[DIM_X / DIM_Y];
#pragma unroll
        for(int i = 0, j = 0; i < DIM_X; i += DIM_Y, j++)
        {
            if(tx1 < n && (ty1 + i) < m && i < m)
            {
                if(UNIT_STRIDE_0)
                {
                    // val[j] = Handler<T_I, cbtype>::read(
                    //     input, in_offset + tx1 + (ty1 + i) * ld_in, load_cb_fn, load_cb_data);
                }
                else
                {
                    // val[j] = Handler<T_I, cbtype>::read(input,
                    //                                     in_offset + tx1 * stride_0_in
                    //                                         + (ty1 + i) * ld_in,
                    //                                     load_cb_fn,
                    //                                     load_cb_data);
                }
            }
        }
#pragma unroll
        for(int i = 0, j = 0; i < DIM_X; i += DIM_Y, j++)
        {
            if(tx1 < n && (ty1 + i) < m && i < m)
            {
                shared[tx1][ty1 + i] = val[j]; // the transpose taking place here
            }
        }
        wItem.barrier(sycl::access::fence_space::local_space);// __syncthreads();
#pragma unroll
        for(int i = 0, j = 0; i < DIM_X; i += DIM_Y, j++)
        {
            if(tx1 < m && (ty1 + i) < n && i < n)
            {
                val[j] = shared[ty1 + i][tx1];
            }
        }
#pragma unroll
        for(int i = 0, j = 0; i < DIM_X; i += DIM_Y, j++)
        {
            // reconfigure the threads
            if(tx1 < m && (ty1 + i) < n && i < n)
            {
                if(UNIT_STRIDE_0)
                {
                    // Handler<T_O, cbtype>::write(output,
                    //                             out_offset + tx1 + (i + ty1) * ld_out,
                    //                             val[j],
                    //                             store_cb_fn,
                    //                             store_cb_data);
                }
                else
                {
                    // Handler<T_O, cbtype>::write(output,
                    //                             out_offset + tx1 * stride_0_out
                    //                                 + (i + ty1) * ld_out,
                    //                             val[j],
                    //                             store_cb_fn,
                    //                             store_cb_data);
                }
            }
        }
    }
}

// global function for transpose scheme
template <typename T,
          typename T_I,
          typename T_O,
          size_t       DIM_X,
          size_t       DIM_Y,
          bool         ALL,
          bool         UNIT_STRIDE_0,
          bool         DIAGONAL,
          CallbackType cbtype>
void //__launch_bounds__(DIM_X* DIM_Y) //
    transpose_kernel2_scheme(
                             sycl::range<3> grid,
                             sycl::rang<3> threads,
                             size_t shared,           
                             sycl::queue rocfftQueue,
                             const T_I* input,
                             T_O*       output,
                             T*         twiddles_large,
                             size_t*    lengths,
                             size_t*    stride_in,
                             size_t*    stride_out,
                             size_t     ld_in,
                             size_t     ld_out,
                             size_t     m,
                             size_t     n,
                             void* __restrict__ load_cb_fn,
                             void* __restrict__ load_cb_data,
                             uint32_t load_cb_lds_bytes,
                             void* __restrict__ store_cb_fn,
                             void* __restrict__ store_cb_data)
{
    
    rocfftQueue.submit([&](cl::sycl::handler &cgh) {
    //missing accessors
    sycl::accessor <T, 2, sycl::access::mode::read_write, sycl::access::target::local>
                        shared(sycl::range<2>(DIM_X, DIM_X), cgh);
    cgh.parallel_for<class transpose_kernel2_scheme>(sycl::nd_range<3>(grid, threads),
	                   [=](sycl::nd_item<3> wItem)) {
    size_t iOffset = 0;
    size_t oOffset = 0;

    size_t counter_mod =  wItem.get_group(2);//hipBlockIdx_z;

    iOffset += counter_mod * stride_in[3];
    oOffset += counter_mod * stride_out[3];

    size_t tileBlockIdx_x, tileBlockIdx_y;
    if(DIAGONAL) // diagonal reordering
    {
        //TODO: template and simplify index calc for square case if necessary
        size_t bid     = wItem.get_group(0) /*hipBlockIdx_x*/ + wItem.get_group_range(0) /*gridDim.x*/ * wItem.get_group(1)/*hipBlockIdx_y*/;
        tileBlockIdx_y = bid % wItem.get_group_range(1) /*hipGridDim_y*/;
        tileBlockIdx_x = ((bid / wItem.get_group_range(1) /*hipGridDim_y*/) + tileBlockIdx_y) % wItem.get_group_range(0) /*hipGridDim_x*/;
    }
    else
    {
        tileBlockIdx_x = wItem.get_group(0) /*hipBlockIdx_x*/;
        tileBlockIdx_y = wItem.get_group(1) /*hipBlockIdx_y*/;
    }

    iOffset += tileBlockIdx_x * DIM_X * stride_in[0] + tileBlockIdx_y * DIM_X * ld_in;  //what is DIM_X
    oOffset += tileBlockIdx_x * DIM_X * ld_out + tileBlockIdx_y * DIM_X * stride_out[0];

    if(ALL)
    {
        transpose_tile_device_scheme<T, T_I, T_O, DIM_X, DIM_Y, ALL, UNIT_STRIDE_0, cbtype>(
            input,
            output,
            iOffset,
            oOffset,
            DIM_X,
            DIM_X,
            ld_in,
            ld_out,
            stride_in[0],
            stride_out[0],
            load_cb_fn,
            load_cb_data,
            load_cb_lds_bytes,
            store_cb_fn,
            store_cb_data,
            wItem,
            shared);
    }
    else
    {
        size_t mm = min(m - tileBlockIdx_y * DIM_X, DIM_X); // the partial case along m
        size_t nn = min(n - tileBlockIdx_x * DIM_X, DIM_X); // the partial case along n
        transpose_tile_device_scheme<T, T_I, T_O, DIM_X, DIM_Y, ALL, UNIT_STRIDE_0, cbtype>(
            input,
            output,
            iOffset,
            oOffset,
            mm,
            nn,
            ld_in,
            ld_out,
            stride_in[0],
            stride_out[0],
            load_cb_fn,
            load_cb_data,
            load_cb_lds_bytes,
            store_cb_fn,
            store_cb_data,
            wItem,
            shared);
    }
                       });
}


using namespace std;

int main() {
    sycl::device device = sycl::default_selector{}.select_device();
    
    //sycl::queue queue(sycl::default_selector{});

    
    sycl::queue queue(device, [] (sycl::exception_list el) {
       for (auto ex : el) { std::rethrow_exception(ex); }
    });

    sycl::range<3> pp(1, 1, 1);
    sycl::range<3> pp1(1, 1, 1);
    sycl::range<3> pp2(1, 1, 1);

    int* a1 = nullptr;
    int* a2 = nullptr;
    int* a3 = nullptr;
    transpose_kernel2_scheme<int, int, int, 2, 2, flase, false, false, void*>(pp, pp1, 2, queue, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 2, 2, 2, 2, nullptr, nullptr, 3, nullptr, nullptr);
}


