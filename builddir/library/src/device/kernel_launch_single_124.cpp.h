
#include "kernel_launch.h" 

#include "rocfft_kernel_2250.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_2250, fft_fwd_ip_len2250, fft_back_ip_len2250, fft_fwd_op_len2250, fft_back_op_len2250, float2)
