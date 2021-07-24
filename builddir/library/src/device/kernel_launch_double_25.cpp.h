
#include "kernel_launch.h" 

#include "rocfft_kernel_60.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_60, fft_fwd_ip_len60, fft_back_ip_len60, fft_fwd_op_len60, fft_back_op_len60, double2)
