
#include "kernel_launch.h" 

#include "rocfft_kernel_52.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_52, fft_fwd_ip_len52, fft_back_ip_len52, fft_fwd_op_len52, fft_back_op_len52, float2)
