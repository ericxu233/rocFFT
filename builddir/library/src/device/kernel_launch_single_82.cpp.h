
#include "kernel_launch.h" 

#include "rocfft_kernel_1800.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_1800, fft_fwd_ip_len1800, fft_back_ip_len1800, fft_fwd_op_len1800, fft_back_op_len1800, float2)
