
#include "kernel_launch.h" 

#include "rocfft_kernel_128.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_128, fft_fwd_ip_len128, fft_back_ip_len128, fft_fwd_op_len128, fft_back_op_len128, double2)