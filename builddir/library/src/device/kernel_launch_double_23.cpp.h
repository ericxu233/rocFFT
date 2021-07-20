
#include "kernel_launch.h" 

#include "rocfft_kernel_121.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_121, fft_fwd_ip_len121, fft_back_ip_len121, fft_fwd_op_len121, fft_back_op_len121, double2)
