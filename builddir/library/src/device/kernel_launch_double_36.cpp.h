
#include "kernel_launch.h" 

#include "rocfft_kernel_250.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_250, fft_fwd_ip_len250, fft_back_ip_len250, fft_fwd_op_len250, fft_back_op_len250, double2)
