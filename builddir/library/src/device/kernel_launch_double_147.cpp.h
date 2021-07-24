
#include "kernel_launch.h" 

#include "rocfft_kernel_4050.h" 

//double precision 
POWX_SMALL_GENERATOR( rocfft_internal_dfn_dp_ci_ci_stoc_4050, fft_fwd_ip_len4050, fft_back_ip_len4050, fft_fwd_op_len4050, fft_back_op_len4050, double2)
