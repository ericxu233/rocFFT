
#include "kernel_launch.h" 

#include "rocfft_kernel_540.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_540, fft_fwd_ip_len540, fft_back_ip_len540, fft_fwd_op_len540, fft_back_op_len540, float2)
