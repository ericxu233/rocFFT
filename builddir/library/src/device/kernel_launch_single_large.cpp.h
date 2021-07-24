
#include "kernel_launch.h" 

//single precision 

#include "rocfft_kernel_50_sbcc.h" 
POWX_LARGE_SBCC_GENERATOR( rocfft_internal_dfn_sp_ci_ci_sbcc_50, fft_fwd_ip_len50_sbcc, fft_back_ip_len50_sbcc, fft_fwd_op_len50_sbcc, fft_back_op_len50_sbcc, float2)
#include "rocfft_kernel_50_sbrc.h" 
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc_50, fft_fwd_op_len50_sbrc, fft_back_op_len50_sbrc, float2, SBRC_2D, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_xy_z_tile_aligned_50, fft_fwd_op_len50_sbrc, fft_back_op_len50_sbrc, float2, SBRC_3D_FFT_TRANS_XY_Z, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_z_xy_tile_aligned_50, fft_fwd_op_len50_sbrc, fft_back_op_len50_sbrc, float2, SBRC_3D_FFT_TRANS_Z_XY, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_erc_trans_z_xy_tile_aligned_50, fft_fwd_op_len50_sbrc, fft_back_op_len50_sbrc, float2, SBRC_3D_FFT_ERC_TRANS_Z_XY, TILE_ALIGNED)
#include "rocfft_kernel_64_sbcc.h" 
POWX_LARGE_SBCC_GENERATOR( rocfft_internal_dfn_sp_ci_ci_sbcc_64, fft_fwd_ip_len64_sbcc, fft_back_ip_len64_sbcc, fft_fwd_op_len64_sbcc, fft_back_op_len64_sbcc, float2)
#include "rocfft_kernel_64_sbrc.h" 
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc_64, fft_fwd_op_len64_sbrc, fft_back_op_len64_sbrc, float2, SBRC_2D, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_xy_z_tile_aligned_64, fft_fwd_op_len64_sbrc, fft_back_op_len64_sbrc, float2, SBRC_3D_FFT_TRANS_XY_Z, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_z_xy_tile_aligned_64, fft_fwd_op_len64_sbrc, fft_back_op_len64_sbrc, float2, SBRC_3D_FFT_TRANS_Z_XY, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_erc_trans_z_xy_tile_aligned_64, fft_fwd_op_len64_sbrc, fft_back_op_len64_sbrc, float2, SBRC_3D_FFT_ERC_TRANS_Z_XY, TILE_ALIGNED)
#include "rocfft_kernel_81_sbcc.h" 
POWX_LARGE_SBCC_GENERATOR( rocfft_internal_dfn_sp_ci_ci_sbcc_81, fft_fwd_ip_len81_sbcc, fft_back_ip_len81_sbcc, fft_fwd_op_len81_sbcc, fft_back_op_len81_sbcc, float2)
#include "rocfft_kernel_81_sbrc.h" 
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc_81, fft_fwd_op_len81_sbrc, fft_back_op_len81_sbrc, float2, SBRC_2D, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_xy_z_tile_aligned_81, fft_fwd_op_len81_sbrc, fft_back_op_len81_sbrc, float2, SBRC_3D_FFT_TRANS_XY_Z, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_z_xy_tile_aligned_81, fft_fwd_op_len81_sbrc, fft_back_op_len81_sbrc, float2, SBRC_3D_FFT_TRANS_Z_XY, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_erc_trans_z_xy_tile_aligned_81, fft_fwd_op_len81_sbrc, fft_back_op_len81_sbrc, float2, SBRC_3D_FFT_ERC_TRANS_Z_XY, TILE_ALIGNED)
#include "rocfft_kernel_100_sbcc.h" 
POWX_LARGE_SBCC_GENERATOR( rocfft_internal_dfn_sp_ci_ci_sbcc_100, fft_fwd_ip_len100_sbcc, fft_back_ip_len100_sbcc, fft_fwd_op_len100_sbcc, fft_back_op_len100_sbcc, float2)
#include "rocfft_kernel_100_sbrc.h" 
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc_100, fft_fwd_op_len100_sbrc, fft_back_op_len100_sbrc, float2, SBRC_2D, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_xy_z_tile_aligned_100, fft_fwd_op_len100_sbrc, fft_back_op_len100_sbrc, float2, SBRC_3D_FFT_TRANS_XY_Z, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_z_xy_tile_aligned_100, fft_fwd_op_len100_sbrc, fft_back_op_len100_sbrc, float2, SBRC_3D_FFT_TRANS_Z_XY, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_erc_trans_z_xy_tile_aligned_100, fft_fwd_op_len100_sbrc, fft_back_op_len100_sbrc, float2, SBRC_3D_FFT_ERC_TRANS_Z_XY, TILE_ALIGNED)
#include "rocfft_kernel_128_sbcc.h" 
POWX_LARGE_SBCC_GENERATOR( rocfft_internal_dfn_sp_ci_ci_sbcc_128, fft_fwd_ip_len128_sbcc, fft_back_ip_len128_sbcc, fft_fwd_op_len128_sbcc, fft_back_op_len128_sbcc, float2)
#include "rocfft_kernel_128_sbrc.h" 
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc_128, fft_fwd_op_len128_sbrc, fft_back_op_len128_sbrc, float2, SBRC_2D, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_xy_z_tile_aligned_128, fft_fwd_op_len128_sbrc, fft_back_op_len128_sbrc, float2, SBRC_3D_FFT_TRANS_XY_Z, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_z_xy_tile_aligned_128, fft_fwd_op_len128_sbrc, fft_back_op_len128_sbrc, float2, SBRC_3D_FFT_TRANS_Z_XY, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_erc_trans_z_xy_tile_aligned_128, fft_fwd_op_len128_sbrc, fft_back_op_len128_sbrc, float2, SBRC_3D_FFT_ERC_TRANS_Z_XY, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_xy_z_diagonal_128, fft_fwd_op_len128_sbrc, fft_back_op_len128_sbrc, float2, SBRC_3D_FFT_TRANS_XY_Z, DIAGONAL)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_z_xy_diagonal_128, fft_fwd_op_len128_sbrc, fft_back_op_len128_sbrc, float2, SBRC_3D_FFT_TRANS_Z_XY, DIAGONAL)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_erc_trans_z_xy_diagonal_128, fft_fwd_op_len128_sbrc, fft_back_op_len128_sbrc, float2, SBRC_3D_FFT_ERC_TRANS_Z_XY, DIAGONAL)
#include "rocfft_kernel_200_sbcc.h" 
POWX_LARGE_SBCC_GENERATOR( rocfft_internal_dfn_sp_ci_ci_sbcc_200, fft_fwd_ip_len200_sbcc, fft_back_ip_len200_sbcc, fft_fwd_op_len200_sbcc, fft_back_op_len200_sbcc, float2)
#include "rocfft_kernel_200_sbrc.h" 
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc_200, fft_fwd_op_len200_sbrc, fft_back_op_len200_sbrc, float2, SBRC_2D, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_xy_z_tile_aligned_200, fft_fwd_op_len200_sbrc, fft_back_op_len200_sbrc, float2, SBRC_3D_FFT_TRANS_XY_Z, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_z_xy_tile_aligned_200, fft_fwd_op_len200_sbrc, fft_back_op_len200_sbrc, float2, SBRC_3D_FFT_TRANS_Z_XY, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_erc_trans_z_xy_tile_aligned_200, fft_fwd_op_len200_sbrc, fft_back_op_len200_sbrc, float2, SBRC_3D_FFT_ERC_TRANS_Z_XY, TILE_ALIGNED)
#include "rocfft_kernel_256_sbcc.h" 
POWX_LARGE_SBCC_GENERATOR( rocfft_internal_dfn_sp_ci_ci_sbcc_256, fft_fwd_ip_len256_sbcc, fft_back_ip_len256_sbcc, fft_fwd_op_len256_sbcc, fft_back_op_len256_sbcc, float2)
#include "rocfft_kernel_256_sbrc.h" 
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc_256, fft_fwd_op_len256_sbrc, fft_back_op_len256_sbrc, float2, SBRC_2D, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_xy_z_tile_aligned_256, fft_fwd_op_len256_sbrc, fft_back_op_len256_sbrc, float2, SBRC_3D_FFT_TRANS_XY_Z, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_z_xy_tile_aligned_256, fft_fwd_op_len256_sbrc, fft_back_op_len256_sbrc, float2, SBRC_3D_FFT_TRANS_Z_XY, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_erc_trans_z_xy_tile_aligned_256, fft_fwd_op_len256_sbrc, fft_back_op_len256_sbrc, float2, SBRC_3D_FFT_ERC_TRANS_Z_XY, TILE_ALIGNED)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_xy_z_diagonal_256, fft_fwd_op_len256_sbrc, fft_back_op_len256_sbrc, float2, SBRC_3D_FFT_TRANS_XY_Z, DIAGONAL)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_trans_z_xy_diagonal_256, fft_fwd_op_len256_sbrc, fft_back_op_len256_sbrc, float2, SBRC_3D_FFT_TRANS_Z_XY, DIAGONAL)
POWX_LARGE_SBRC_GENERATOR( rocfft_internal_dfn_sp_op_ci_ci_sbrc3d_fft_erc_trans_z_xy_diagonal_256, fft_fwd_op_len256_sbrc, fft_back_op_len256_sbrc, float2, SBRC_3D_FFT_ERC_TRANS_Z_XY, DIAGONAL)
