
#include "kernel_launch.h" 

#include "rocfft_kernel_21.h" 
#include "rocfft_kernel_22.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_21, fft_fwd_ip_len21, fft_back_ip_len21, fft_fwd_op_len21, fft_back_op_len21, double2)
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_22, fft_fwd_ip_len22, fft_back_ip_len22, fft_fwd_op_len22, fft_back_op_len22, double2)
