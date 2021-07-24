
#include "kernel_launch.h" 

#include "rocfft_kernel_960.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_960, fft_fwd_ip_len960, fft_back_ip_len960, fft_fwd_op_len960, fft_back_op_len960, double2)
