
#include "kernel_launch.h" 

#include "rocfft_kernel_1440.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_1440, fft_fwd_ip_len1440, fft_back_ip_len1440, fft_fwd_op_len1440, fft_back_op_len1440, double2)
