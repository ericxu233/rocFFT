
#include "kernel_launch.h" 

#include "rocfft_kernel_180.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_180, fft_fwd_ip_len180, fft_back_ip_len180, fft_fwd_op_len180, fft_back_op_len180, double2)
