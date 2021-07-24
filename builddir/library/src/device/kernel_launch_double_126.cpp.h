
#include "kernel_launch.h" 

#include "rocfft_kernel_2400.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_2400, fft_fwd_ip_len2400, fft_back_ip_len2400, fft_fwd_op_len2400, fft_back_op_len2400, double2)
