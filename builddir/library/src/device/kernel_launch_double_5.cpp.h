
#include "kernel_launch.h" 

#include "rocfft_kernel_11.h" 
#include "rocfft_kernel_12.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_11, fft_fwd_ip_len11, fft_back_ip_len11, fft_fwd_op_len11, fft_back_op_len11, double2)
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_12, fft_fwd_ip_len12, fft_back_ip_len12, fft_fwd_op_len12, fft_back_op_len12, double2)
