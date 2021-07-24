
#include "kernel_launch.h" 

#include "rocfft_kernel_1296.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_1296, fft_fwd_ip_len1296, fft_back_ip_len1296, fft_fwd_op_len1296, fft_back_op_len1296, float2)
