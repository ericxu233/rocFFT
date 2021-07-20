
#include "kernel_launch.h" 

#include "rocfft_kernel_750.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_750, fft_fwd_ip_len750, fft_back_ip_len750, fft_fwd_op_len750, fft_back_op_len750, double2)
