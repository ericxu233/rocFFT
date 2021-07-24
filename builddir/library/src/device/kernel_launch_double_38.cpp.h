
#include "kernel_launch.h" 

#include "rocfft_kernel_112.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_112, fft_fwd_ip_len112, fft_back_ip_len112, fft_fwd_op_len112, fft_back_op_len112, double2)
