
#include "kernel_launch.h" 

#include "rocfft_kernel_2160.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_2160, fft_fwd_ip_len2160, fft_back_ip_len2160, fft_fwd_op_len2160, fft_back_op_len2160, double2)
