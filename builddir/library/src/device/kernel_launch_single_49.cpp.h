
#include "kernel_launch.h" 

#include "rocfft_kernel_169.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_169, fft_fwd_ip_len169, fft_back_ip_len169, fft_fwd_op_len169, fft_back_op_len169, float2)
