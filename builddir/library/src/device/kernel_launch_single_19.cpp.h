
#include "kernel_launch.h" 

#include "rocfft_kernel_81.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_81, fft_fwd_ip_len81, fft_back_ip_len81, fft_fwd_op_len81, fft_back_op_len81, float2)
