#include "kernel_launch.h"
#include "rocfft_kernel_2D_25_25.h"
POWX_SMALL_GENERATOR(rocfft_internal_dfn_sp_ci_ci_2D_25_25, fft_fwd_ip_2D_25_25, fft_back_ip_2D_25_25, fft_fwd_op_2D_25_25, fft_back_op_2D_25_25, float2)
#include "rocfft_kernel_2D_25_125.h"
POWX_SMALL_GENERATOR(rocfft_internal_dfn_sp_ci_ci_2D_25_125, fft_fwd_ip_2D_25_125, fft_back_ip_2D_25_125, fft_fwd_op_2D_25_125, fft_back_op_2D_25_125, float2)
#include "rocfft_kernel_2D_125_25.h"
POWX_SMALL_GENERATOR(rocfft_internal_dfn_sp_ci_ci_2D_125_25, fft_fwd_ip_2D_125_25, fft_back_ip_2D_125_25, fft_fwd_op_2D_125_25, fft_back_op_2D_125_25, float2)
