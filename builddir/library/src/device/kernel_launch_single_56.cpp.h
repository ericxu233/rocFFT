
#include "kernel_launch.h" 

#include "rocfft_kernel_224.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_224, fft_fwd_ip_len224, fft_back_ip_len224, fft_fwd_op_len224, fft_back_op_len224, float2)
