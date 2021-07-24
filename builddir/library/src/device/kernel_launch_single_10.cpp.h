
#include "kernel_launch.h" 

#include "rocfft_kernel_24.h" 
#include "rocfft_kernel_25.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_24, fft_fwd_ip_len24, fft_back_ip_len24, fft_fwd_op_len24, fft_back_op_len24, float2)
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_25, fft_fwd_ip_len25, fft_back_ip_len25, fft_fwd_op_len25, fft_back_op_len25, float2)
