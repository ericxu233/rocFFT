
#include "kernel_launch.h" 

#include "rocfft_kernel_512.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_512, fft_fwd_ip_len512, fft_back_ip_len512, fft_fwd_op_len512, fft_back_op_len512, double2)
