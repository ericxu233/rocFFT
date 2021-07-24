
#include "kernel_launch.h" 

#include "rocfft_kernel_15.h" 
#include "rocfft_kernel_16.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_15, fft_fwd_ip_len15, fft_back_ip_len15, fft_fwd_op_len15, fft_back_op_len15, double2)
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_16, fft_fwd_ip_len16, fft_back_ip_len16, fft_fwd_op_len16, fft_back_op_len16, double2)
