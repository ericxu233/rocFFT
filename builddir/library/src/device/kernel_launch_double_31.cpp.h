
#include "kernel_launch.h" 

#include "rocfft_kernel_84.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_84, fft_fwd_ip_len84, fft_back_ip_len84, fft_fwd_op_len84, fft_back_op_len84, double2)
