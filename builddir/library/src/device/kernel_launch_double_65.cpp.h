
#include "kernel_launch.h" 

#include "rocfft_kernel_320.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_320, fft_fwd_ip_len320, fft_back_ip_len320, fft_fwd_op_len320, fft_back_op_len320, double2)
