
#include "kernel_launch.h" 

#include "rocfft_kernel_384.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_384, fft_fwd_ip_len384, fft_back_ip_len384, fft_fwd_op_len384, fft_back_op_len384, float2)
