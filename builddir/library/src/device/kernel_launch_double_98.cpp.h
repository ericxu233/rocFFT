
#include "kernel_launch.h" 

#include "rocfft_kernel_2700.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_2700, fft_fwd_ip_len2700, fft_back_ip_len2700, fft_fwd_op_len2700, fft_back_op_len2700, double2)
