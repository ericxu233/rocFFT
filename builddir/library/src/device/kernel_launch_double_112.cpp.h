
#include "kernel_launch.h" 

#include "rocfft_kernel_1600.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_1600, fft_fwd_ip_len1600, fft_back_ip_len1600, fft_fwd_op_len1600, fft_back_op_len1600, double2)
