
#include "kernel_launch.h" 

#include "rocfft_kernel_240.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_240, fft_fwd_ip_len240, fft_back_ip_len240, fft_fwd_op_len240, fft_back_op_len240, double2)
