
#include "kernel_launch.h" 

#include "rocfft_kernel_176.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_176, fft_fwd_ip_len176, fft_back_ip_len176, fft_fwd_op_len176, fft_back_op_len176, float2)
