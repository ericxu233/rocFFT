
#include "kernel_launch.h" 

#include "rocfft_kernel_16.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_16, fft_fwd_ip_len16, fft_back_ip_len16, fft_fwd_op_len16, fft_back_op_len16, float2)
