
#include "kernel_launch.h" 

#include "rocfft_kernel_270.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_270, fft_fwd_ip_len270, fft_back_ip_len270, fft_fwd_op_len270, fft_back_op_len270, double2)
