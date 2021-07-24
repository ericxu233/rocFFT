
#include "kernel_launch.h" 

#include "rocfft_kernel_13.h" 
#include "rocfft_kernel_14.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_13, fft_fwd_ip_len13, fft_back_ip_len13, fft_fwd_op_len13, fft_back_op_len13, double2)
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_14, fft_fwd_ip_len14, fft_back_ip_len14, fft_fwd_op_len14, fft_back_op_len14, double2)
