
#include "kernel_launch.h" 

#include "rocfft_kernel_2025.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_2025, fft_fwd_ip_len2025, fft_back_ip_len2025, fft_fwd_op_len2025, fft_back_op_len2025, double2)
