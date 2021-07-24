
#include "kernel_launch.h" 

#include "rocfft_kernel_1620.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_1620, fft_fwd_ip_len1620, fft_back_ip_len1620, fft_fwd_op_len1620, fft_back_op_len1620, float2)
