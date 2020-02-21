// gtclang (0.0.4-5fc8d35e-x86_64-GNU-9.2.1)
// based on LLVM/Clang (9.0.0), Dawn (0.0.2)
// Generated on 2020-02-20  14:11:44

#define DAWN_GENERATED 1
#define DAWN_BACKEND_T CUDA
#ifndef BOOST_RESULT_OF_USE_TR1
 #define BOOST_RESULT_OF_USE_TR1 1
#endif
#ifndef BOOST_NO_CXX11_DECLTYPE
 #define BOOST_NO_CXX11_DECLTYPE 1
#endif
#ifndef GRIDTOOLS_DAWN_HALO_EXTENT
 #define GRIDTOOLS_DAWN_HALO_EXTENT 3
#endif
#ifndef BOOST_PP_VARIADICS
 #define BOOST_PP_VARIADICS 1
#endif
#ifndef BOOST_FUSION_DONT_USE_PREPROCESSED_FILES
 #define BOOST_FUSION_DONT_USE_PREPROCESSED_FILES 1
#endif
#ifndef BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
 #define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS 1
#endif
#ifndef GT_VECTOR_LIMIT_SIZE
 #define GT_VECTOR_LIMIT_SIZE 30
#endif
#ifndef BOOST_FUSION_INVOKE_MAX_ARITY
 #define BOOST_FUSION_INVOKE_MAX_ARITY GT_VECTOR_LIMIT_SIZE
#endif
#ifndef FUSION_MAX_VECTOR_SIZE
 #define FUSION_MAX_VECTOR_SIZE GT_VECTOR_LIMIT_SIZE
#endif
#ifndef FUSION_MAX_MAP_SIZE
 #define FUSION_MAX_MAP_SIZE GT_VECTOR_LIMIT_SIZE
#endif
#ifndef BOOST_MPL_LIMIT_VECTOR_SIZE
 #define BOOST_MPL_LIMIT_VECTOR_SIZE GT_VECTOR_LIMIT_SIZE
#endif
#include <driver-includes/gridtools_includes.hpp>
using namespace gridtools::dawn;
//===--------------------------------------------------------------------------------*- C++ -*-===//
//                         _       _
//                        | |     | |
//                    __ _| |_ ___| | __ _ _ __   __ _
//                   / _` | __/ __| |/ _` | '_ \ / _` |
//                  | (_| | || (__| | (_| | | | | (_| |
//                   \__, |\__\___|_|\__,_|_| |_|\__, | - GridTools Clang DSL
//                    __/ |                       __/ |
//                   |___/                       |___/
//
//
//  This file is distributed under the MIT License (MIT).
//  See LICENSE.txt for details.
//
//===------------------------------------------------------------------------------------------===//

#include "gtclang_dsl_defs/gtclang_dsl.hpp"

using namespace gtclang::dsl;

namespace dawn_generated {
namespace cuda {
__global__ void __launch_bounds__(128)
    coriolis_stencil_stencil50_ms104_kernel(const int isize, const int jsize, const int ksize, const int stride_111_1,
                                            const int stride_111_2, ::dawn::float_type* const u_nnow,
                                            ::dawn::float_type* const v_nnow, ::dawn::float_type* const fc,
                                            ::dawn::float_type* const u_tens, ::dawn::float_type* const v_tens) {
  // Start kernel
  const unsigned int nx = isize;
  const unsigned int ny = jsize;
  const int block_size_i = (blockIdx.x + 1) * 32 < nx ? 32 : nx - blockIdx.x * 32;
  const int block_size_j = (blockIdx.y + 1) * 4 < ny ? 4 : ny - blockIdx.y * 4;

  // computing the global position in the physical domain

  // In a typical cuda block we have the following regions

  // aa bbbbbbbb cc

  // aa bbbbbbbb cc

  // hh dddddddd ii

  // hh dddddddd ii

  // hh dddddddd ii

  // hh dddddddd ii

  // ee ffffffff gg

  // ee ffffffff gg

  // Regions b,d,f have warp (or multiple of warp size)

  // Size of regions a, c, h, i, e, g are determined by max_extent_t

  // Regions b,d,f are easily executed by dedicated warps (one warp for each line)

  // Regions (a,h,e) and (c,i,g) are executed by two specialized warp
  int iblock = 0 - 1;
  int jblock = 0 - 1;
  if (threadIdx.y < +4) {
    iblock = threadIdx.x;
    jblock = (int)threadIdx.y + 0;
  }
  // initialized iterators
  int idx111 = (blockIdx.x * 32 + iblock) * 1 + (blockIdx.y * 4 + jblock) * stride_111_1;

  // jump iterators to match the intersection of beginning of next interval and the parallel execution block
  idx111 += max(0, blockIdx.z * 4) * stride_111_2;
  int kleg_lower_bound = max(0, blockIdx.z * 4);
  int kleg_upper_bound = min(ksize - 1 + 0, (blockIdx.z + 1) * 4 - 1);
  ;
  for (int k = kleg_lower_bound + 0; k <= kleg_upper_bound + 0; ++k) {
    if (iblock >= 0 && iblock <= block_size_i - 1 + 0 && jblock >= 0 && jblock <= block_size_j - 1 + 0) {
      ::dawn::float_type __local_z_fv_north_98 =
          (__ldg(&(fc[idx111])) * (__ldg(&(v_nnow[idx111])) + __ldg(&(v_nnow[idx111 + 1 * 1]))));
      ::dawn::float_type __local_z_fv_south_99 =
          (__ldg(&(fc[idx111 + stride_111_1 * -1])) *
           (__ldg(&(v_nnow[idx111 + stride_111_1 * -1])) + __ldg(&(v_nnow[idx111 + 1 * 1 + stride_111_1 * -1]))));
      u_tens[idx111] += ((::dawn::float_type)0.25 * (__local_z_fv_north_98 + __local_z_fv_south_99));
      ::dawn::float_type __local_z_fu_east_101 =
          (__ldg(&(fc[idx111])) * (__ldg(&(u_nnow[idx111])) + __ldg(&(u_nnow[idx111 + stride_111_1 * 1]))));
      ::dawn::float_type __local_z_fu_west_102 =
          (__ldg(&(fc[idx111 + 1 * -1])) *
           (__ldg(&(u_nnow[idx111 + 1 * -1])) + __ldg(&(u_nnow[idx111 + 1 * -1 + stride_111_1 * 1]))));
      v_tens[idx111] -= ((::dawn::float_type)0.25 * (__local_z_fu_east_101 + __local_z_fu_west_102));
    }
    // Slide kcaches

    // increment iterators
    idx111 += stride_111_2;
  }
}

class coriolis_stencil {
 public:
  struct sbase : public timer_cuda {
    sbase(std::string name) : timer_cuda(name) {}

    double get_time() { return total_time(); }
  };

  struct stencil_50 : public sbase {
    // Members

    // Temporary storage typedefs
    using tmp_halo_t = gridtools::halo<0, 0, 0, 0, 0>;
    using tmp_meta_data_t = storage_traits_t::storage_info_t<0, 5, tmp_halo_t>;
    using tmp_storage_t = storage_traits_t::data_store_t<::dawn::float_type, tmp_meta_data_t>;
    const gridtools::dawn::domain m_dom;

   public:
    stencil_50(const gridtools::dawn::domain& dom_, int rank, int xcols, int ycols)
        : sbase("stencil_50"), m_dom(dom_) {}

    void run(storage_ijk_t u_nnow_ds, storage_ijk_t v_nnow_ds, storage_ijk_t fc_ds,
             storage_ijk_t u_tens_ds, storage_ijk_t v_tens_ds) {
      // starting timers
      start();
      {
        ;
        gridtools::data_view<storage_ijk_t> u_nnow = gridtools::make_device_view(u_nnow_ds);
        gridtools::data_view<storage_ijk_t> v_nnow = gridtools::make_device_view(v_nnow_ds);
        gridtools::data_view<storage_ijk_t> fc = gridtools::make_device_view(fc_ds);
        gridtools::data_view<storage_ijk_t> u_tens = gridtools::make_device_view(u_tens_ds);
        gridtools::data_view<storage_ijk_t> v_tens = gridtools::make_device_view(v_tens_ds);
        const unsigned int nx = m_dom.isize() - m_dom.iminus() - m_dom.iplus();
        const unsigned int ny = m_dom.jsize() - m_dom.jminus() - m_dom.jplus();
        const unsigned int nz = m_dom.ksize() - m_dom.kminus() - m_dom.kplus();
        dim3 threads(32, 4 + 0, 1);
        const unsigned int nbx = (nx + 32 - 1) / 32;
        const unsigned int nby = (ny + 4 - 1) / 4;
        const unsigned int nbz = (m_dom.ksize() + 4 - 1) / 4;
        dim3 blocks(nbx, nby, nbz);
        coriolis_stencil_stencil50_ms104_kernel<<<blocks, threads>>>(
            nx, ny, nz, u_tens_ds.strides()[1], u_tens_ds.strides()[2],
            (u_nnow.data() + u_nnow_ds.get_storage_info_ptr()->index(u_nnow.begin<0>(), u_nnow.begin<1>(), 0)),
            (v_nnow.data() + v_nnow_ds.get_storage_info_ptr()->index(v_nnow.begin<0>(), v_nnow.begin<1>(), 0)),
            (fc.data() + fc_ds.get_storage_info_ptr()->index(fc.begin<0>(), fc.begin<1>(), 0)),
            (u_tens.data() + u_tens_ds.get_storage_info_ptr()->index(u_tens.begin<0>(), u_tens.begin<1>(), 0)),
            (v_tens.data() + v_tens_ds.get_storage_info_ptr()->index(v_tens.begin<0>(), v_tens.begin<1>(), 0)));
      };

      // stopping timers
      pause();
    }
  };
  static constexpr const char* s_name = "coriolis_stencil";
  stencil_50 m_stencil_50;

 public:
  coriolis_stencil(const coriolis_stencil&) = delete;

  // Members

  // Stencil-Data

  coriolis_stencil(const gridtools::dawn::domain& dom, int rank = 1, int xcols = 1, int ycols = 1)
      : m_stencil_50(dom, rank, xcols, ycols) {}

  template <typename S>
  void sync_storages(S field) {
    field.sync();
  }

  template <typename S0, typename... S>
  void sync_storages(S0 f0, S... fields) {
    f0.sync();
    sync_storages(fields...);
  }

  void run(storage_ijk_t u_nnow, storage_ijk_t v_nnow, storage_ijk_t fc, storage_ijk_t u_tens, storage_ijk_t v_tens) {
    sync_storages(u_nnow, v_nnow, fc, u_tens, v_tens);
    m_stencil_50.run(u_nnow, v_nnow, fc, u_tens, v_tens);
    sync_storages(u_nnow, v_nnow, fc, u_tens, v_tens);
  }

  std::string get_name() const { return std::string(s_name); }

  void reset_meters() { m_stencil_50.reset(); }

  double get_total_time() {
    double res = 0;
    res += m_stencil_50.get_time();
    return res;
  }
};
}  // namespace cuda
}  // namespace dawn_generated

#include "driver-includes/verify.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>

void print(const domain& dom, const gridtools::data_view<storage_ijk_t>& view) {
  for(int i = dom.iminus(); i < std::min(int(dom.isize() - dom.iplus()), view.total_length<0>()); ++i)
    for(int j = dom.jminus(); j < std::min(int(dom.jsize() - dom.jplus()), view.total_length<1>()); ++j)
      for(int k = dom.kminus(); k < std::min(int(dom.ksize() - dom.kplus()), view.total_length<2>()); ++k)
        std::cout << std::setprecision(9) << view(i, j, k) << ' ';
  std::cout << std::endl;
}

int main(int argc, const char** argv) {
  int isize, jsize, ksize, halo;
  if(argc > 1)
    isize = jsize = ksize = atoi(argv[1]);
  else
    isize = jsize = ksize = 12;
  if(argc > 2)
    halo = atoi(argv[2]);
  else
    halo = GRIDTOOLS_DAWN_HALO_EXTENT;

  float time;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  domain dom(isize, jsize, ksize);
  dom.set_halos(halo, halo, halo, halo, 0, 0);
  meta_data_t meta(isize, jsize, ksize+1);
  storage_t u_nnow(meta, "u_nnow"),v_nnow(meta, "v_nnow"),fc(meta, "fc"),u_tens(meta, "u_tens"),v_tens(meta, "v_tens");
  verifier verif(dom);
  verif.fillMath(8,2,1.5,1.5,2,4,u_nnow);
  verif.fillMath(5,1.2,1.3,1.7,2.2,3.5,v_nnow);
  verif.fillMath(2,1.3,1.4,1.6,2.1,3,fc);
  verif.fill(-1,u_tens,v_tens);

  dawn_generated::cuda::coriolis_stencil stencil(dom);
  cudaEventRecord(start, 0);
  stencil.run(u_nnow,v_nnow,fc,u_tens,v_tens);
  cudaEventRecord(stop, 0);

  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);
  std::cerr << "cuda_time = " << (time * 1E-3)  << std::endl;

  //print(dom, make_host_view(u_tens));
  //print(dom, make_host_view(v_tens));
}
