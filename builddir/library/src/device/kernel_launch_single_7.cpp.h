
#include "kernel_launch.h" 

#include "rocfft_kernel_8.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_8, fft_fwd_ip_len8, fft_back_ip_len8, fft_fwd_op_len8, fft_back_op_len8, float2)
