
#include "kernel_launch.h" 

#include "rocfft_kernel_42.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_42, fft_fwd_ip_len42, fft_back_ip_len42, fft_fwd_op_len42, fft_back_op_len42, double2)
