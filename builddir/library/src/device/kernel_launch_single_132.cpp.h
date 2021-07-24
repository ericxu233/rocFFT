
#include "kernel_launch.h" 

#include "rocfft_kernel_2880.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_2880, fft_fwd_ip_len2880, fft_back_ip_len2880, fft_fwd_op_len2880, fft_back_op_len2880, float2)
