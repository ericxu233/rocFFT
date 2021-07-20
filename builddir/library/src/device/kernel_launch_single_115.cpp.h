
#include "kernel_launch.h" 

#include "rocfft_kernel_4096.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_4096, fft_fwd_ip_len4096, fft_back_ip_len4096, fft_fwd_op_len4096, fft_back_op_len4096, float2)
