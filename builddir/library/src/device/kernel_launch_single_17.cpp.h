
#include "kernel_launch.h" 

#include "rocfft_kernel_64.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_64, fft_fwd_ip_len64, fft_back_ip_len64, fft_fwd_op_len64, fft_back_op_len64, float2)
