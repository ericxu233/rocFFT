
#include "kernel_launch.h" 

#include "rocfft_kernel_28.h" 
#include "rocfft_kernel_30.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_28, fft_fwd_ip_len28, fft_back_ip_len28, fft_fwd_op_len28, fft_back_op_len28, double2)
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_30, fft_fwd_ip_len30, fft_back_ip_len30, fft_fwd_op_len30, fft_back_op_len30, double2)
