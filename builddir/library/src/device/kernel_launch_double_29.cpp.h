
#include "kernel_launch.h" 

#include "rocfft_kernel_80.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_80, fft_fwd_ip_len80, fft_back_ip_len80, fft_fwd_op_len80, fft_back_op_len80, double2)
