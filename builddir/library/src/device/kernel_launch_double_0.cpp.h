
#include "kernel_launch.h" 

#include "rocfft_kernel_1.h" 
#include "rocfft_kernel_2.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_1, fft_fwd_ip_len1, fft_back_ip_len1, fft_fwd_op_len1, fft_back_op_len1, double2)
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_2, fft_fwd_ip_len2, fft_back_ip_len2, fft_fwd_op_len2, fft_back_op_len2, double2)
