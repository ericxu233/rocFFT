
#include "kernel_launch.h" 

#include "rocfft_kernel_300.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_300, fft_fwd_ip_len300, fft_back_ip_len300, fft_fwd_op_len300, fft_back_op_len300, float2)
