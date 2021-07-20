
#include "kernel_launch.h" 

#include "rocfft_kernel_1024.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_1024, fft_fwd_ip_len1024, fft_back_ip_len1024, fft_fwd_op_len1024, fft_back_op_len1024, float2)
