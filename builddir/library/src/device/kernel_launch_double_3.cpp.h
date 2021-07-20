
#include "kernel_launch.h" 

#include "rocfft_kernel_4.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_4, fft_fwd_ip_len4, fft_back_ip_len4, fft_fwd_op_len4, fft_back_op_len4, double2)
