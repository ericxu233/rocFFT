
#include "kernel_launch.h" 

#include "rocfft_kernel_36.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_36, fft_fwd_ip_len36, fft_back_ip_len36, fft_fwd_op_len36, fft_back_op_len36, float2)
