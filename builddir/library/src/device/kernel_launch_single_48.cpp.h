
#include "kernel_launch.h" 

#include "rocfft_kernel_168.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_168, fft_fwd_ip_len168, fft_back_ip_len168, fft_fwd_op_len168, fft_back_op_len168, float2)
