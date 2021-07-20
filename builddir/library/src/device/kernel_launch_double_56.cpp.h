
#include "kernel_launch.h" 

#include "rocfft_kernel_729.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_729, fft_fwd_ip_len729, fft_back_ip_len729, fft_fwd_op_len729, fft_back_op_len729, double2)
