
#include "kernel_launch.h" 

#include "rocfft_kernel_40.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_40, fft_fwd_ip_len40, fft_back_ip_len40, fft_fwd_op_len40, fft_back_op_len40, double2)
