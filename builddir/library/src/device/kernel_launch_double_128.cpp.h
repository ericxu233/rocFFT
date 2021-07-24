
#include "kernel_launch.h" 

#include "rocfft_kernel_2500.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_2500, fft_fwd_ip_len2500, fft_back_ip_len2500, fft_fwd_op_len2500, fft_back_op_len2500, double2)
