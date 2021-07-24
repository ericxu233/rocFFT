
#include "kernel_launch.h" 

#include "rocfft_kernel_3000.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_3000, fft_fwd_ip_len3000, fft_back_ip_len3000, fft_fwd_op_len3000, fft_back_op_len3000, float2)
