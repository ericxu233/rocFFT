
#include "kernel_launch.h" 

#include "rocfft_kernel_18.h" 
#include "rocfft_kernel_20.h" 

//single precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_18, fft_fwd_ip_len18, fft_back_ip_len18, fft_fwd_op_len18, fft_back_op_len18, float2)
POWX_SMALL_GENERATOR( rocfft_internal_dfn_sp_ci_ci_stoc_20, fft_fwd_ip_len20, fft_back_ip_len20, fft_fwd_op_len20, fft_back_op_len20, float2)
