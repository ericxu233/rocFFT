
#include "kernel_launch.h" 

#include "rocfft_kernel_3200.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_3200, fft_fwd_ip_len3200, fft_back_ip_len3200, fft_fwd_op_len3200, fft_back_op_len3200, double2)
