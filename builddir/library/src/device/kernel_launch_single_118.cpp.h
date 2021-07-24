
#include "kernel_launch.h" 

#include "rocfft_kernel_1944.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_1944, fft_fwd_ip_len1944, fft_back_ip_len1944, fft_fwd_op_len1944, fft_back_op_len1944, float2)
