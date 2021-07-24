
#include "kernel_launch.h" 

#include "rocfft_kernel_50.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_50, fft_fwd_ip_len50, fft_back_ip_len50, fft_fwd_op_len50, fft_back_op_len50, float2)
