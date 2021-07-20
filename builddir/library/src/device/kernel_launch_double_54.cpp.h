
#include "kernel_launch.h" 

#include "rocfft_kernel_640.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_640, fft_fwd_ip_len640, fft_back_ip_len640, fft_fwd_op_len640, fft_back_op_len640, double2)
