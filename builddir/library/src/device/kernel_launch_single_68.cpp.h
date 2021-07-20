
#include "kernel_launch.h" 

#include "rocfft_kernel_1152.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_1152, fft_fwd_ip_len1152, fft_back_ip_len1152, fft_fwd_op_len1152, fft_back_op_len1152, float2)
