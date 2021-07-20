
#include "kernel_launch.h" 

#include "rocfft_kernel_3645.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_3645, fft_fwd_ip_len3645, fft_back_ip_len3645, fft_fwd_op_len3645, fft_back_op_len3645, float2)
