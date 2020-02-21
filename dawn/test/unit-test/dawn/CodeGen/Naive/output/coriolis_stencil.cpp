#define DAWN_GENERATED 1
#define DAWN_BACKEND_T CXXNAIVE
#ifndef BOOST_RESULT_OF_USE_TR1
 #define BOOST_RESULT_OF_USE_TR1 1
#endif
#ifndef BOOST_NO_CXX11_DECLTYPE
 #define BOOST_NO_CXX11_DECLTYPE 1
#endif
#ifndef GRIDTOOLS_DAWN_HALO_EXTENT
 #define GRIDTOOLS_DAWN_HALO_EXTENT 0
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
namespace dawn_generated{
namespace cxxnaive{

class coriolis_stencil {
private:

  struct stencil_50 {

    // Members

    // Temporary storages
    using tmp_halo_t = gridtools::halo< GRIDTOOLS_DAWN_HALO_EXTENT, GRIDTOOLS_DAWN_HALO_EXTENT, 0>;
    using tmp_meta_data_t = storage_traits_t::storage_info_t< 0, 3, tmp_halo_t >;
    using tmp_storage_t = storage_traits_t::data_store_t< ::dawn::float_type, tmp_meta_data_t>;
    const gridtools::dawn::domain m_dom;

    // Input/Output storages
  public:

    stencil_50(const gridtools::dawn::domain& dom_, int rank, int xcols, int ycols) : m_dom(dom_){}

    void run(storage_ijk_t& u_nnow_, storage_ijk_t& v_nnow_, storage_ijk_t& fc_, storage_ijk_t& u_tens_, storage_ijk_t& v_tens_) {
      int iMin = m_dom.iminus();
      int iMax = m_dom.isize() - m_dom.iplus() - 1;
      int jMin = m_dom.jminus();
      int jMax = m_dom.jsize() - m_dom.jplus() - 1;
      int kMin = m_dom.kminus();
      int kMax = m_dom.ksize() - m_dom.kplus() - 1;
      u_nnow_.sync();
      v_nnow_.sync();
      fc_.sync();
      u_tens_.sync();
      v_tens_.sync();
{      gridtools::data_view<storage_ijk_t> u_nnow= gridtools::make_host_view(u_nnow_);
      std::array<int,3> u_nnow_offsets{0,0,0};
      gridtools::data_view<storage_ijk_t> v_nnow= gridtools::make_host_view(v_nnow_);
      std::array<int,3> v_nnow_offsets{0,0,0};
      gridtools::data_view<storage_ijk_t> fc= gridtools::make_host_view(fc_);
      std::array<int,3> fc_offsets{0,0,0};
      gridtools::data_view<storage_ijk_t> u_tens= gridtools::make_host_view(u_tens_);
      std::array<int,3> u_tens_offsets{0,0,0};
      gridtools::data_view<storage_ijk_t> v_tens= gridtools::make_host_view(v_tens_);
      std::array<int,3> v_tens_offsets{0,0,0};
    for(int k = kMin + 0+0; k <= kMax + 0+0; ++k) {
      for(int i = iMin+0; i  <=  iMax+0; ++i) {
        for(int j = jMin+0; j  <=  jMax+0; ++j) {
::dawn::float_type __local_z_fv_north_98 = (fc(i+0, j+0, k+0) * (v_nnow(i+0, j+0, k+0) + v_nnow(i+1, j+0, k+0)));
::dawn::float_type __local_z_fv_south_99 = (fc(i+0, j+-1, k+0) * (v_nnow(i+0, j+-1, k+0) + v_nnow(i+1, j+-1, k+0)));
u_tens(i+0, j+0, k+0) += ((::dawn::float_type) 0.25 * (__local_z_fv_north_98 + __local_z_fv_south_99));
::dawn::float_type __local_z_fu_east_101 = (fc(i+0, j+0, k+0) * (u_nnow(i+0, j+0, k+0) + u_nnow(i+0, j+1, k+0)));
::dawn::float_type __local_z_fu_west_102 = (fc(i+-1, j+0, k+0) * (u_nnow(i+-1, j+0, k+0) + u_nnow(i+-1, j+1, k+0)));
v_tens(i+0, j+0, k+0) -= ((::dawn::float_type) 0.25 * (__local_z_fu_east_101 + __local_z_fu_west_102));
        }      }    }}      u_nnow_.sync();
      v_nnow_.sync();
      fc_.sync();
      u_tens_.sync();
      v_tens_.sync();
    }
  };
  static constexpr const char* s_name = "coriolis_stencil";
  stencil_50 m_stencil_50;
public:

  coriolis_stencil(const coriolis_stencil&) = delete;

  coriolis_stencil(const gridtools::dawn::domain& dom, int rank = 1, int xcols = 1, int ycols = 1) : m_stencil_50(dom, rank, xcols, ycols){
    assert(dom.isize() >= dom.iminus() + dom.iplus());
    assert(dom.jsize() >= dom.jminus() + dom.jplus());
    assert(dom.ksize() >= dom.kminus() + dom.kplus());
    assert(dom.ksize() >= 1);
  }

  void run(storage_ijk_t u_nnow, storage_ijk_t v_nnow, storage_ijk_t fc, storage_ijk_t u_tens, storage_ijk_t v_tens) {
    m_stencil_50.run(u_nnow,v_nnow,fc,u_tens,v_tens);
  }
};
} // namespace cxxnaive
} // namespace dawn_generated

#include "driver-includes/verify.hpp"
#include <iostream>
#include <iomanip>
#include <omp.h>

void print(const domain& dom, const gridtools::data_view<storage_ijk_t>& view) {
  for(int i = dom.iminus(); i < std::min(int(dom.isize() - dom.iplus()), view.total_length<0>()); ++i)
    for(int j = dom.jminus(); j < std::min(int(dom.jsize() - dom.jplus()), view.total_length<1>()); ++j)
      for(int k = dom.kminus(); k < std::min(int(dom.ksize() - dom.kplus()), view.total_length<2>()); ++k)
        std::cout << std::setprecision(9) << view(i, j, k) << ' ';
  std::cout << std::endl;
}

int main(int argc, const char** argv) {
  int isize = atoi(argv[1]);
  int jsize = isize;
  int ksize = atoi(argv[2]);
  int halo = 3;
  domain dom(isize, jsize, ksize);
  dom.set_halos(halo, halo, halo, halo, 0, 0);
  meta_data_t meta(isize, jsize, ksize + 1);
  storage_t u_nnow(meta, "u_nnow"),v_nnow(meta, "v_nnow"),fc(meta, "fc"),u_tens(meta, "u_tens"),v_tens(meta, "v_tens");

//  #pragma omp parallel
//  {
//    nthreads = omp_get_num_threads();
//  }
//  printf("nthreads=%d\n", nthreads);


  verifier verif(dom);
  verif.fillMath(8,2,1.5,1.5,2,4,u_nnow);
  verif.fillMath(5,1.2,1.3,1.7,2.2,3.5,v_nnow);
  verif.fillMath(2,1.3,1.4,1.6,2.1,3,fc);
  verif.fill(-1,u_tens,v_tens);

  dawn_generated::cxxnaive::coriolis_stencil stencil(dom);
  float time = omp_get_wtime();
  stencil.run(u_nnow,v_nnow,fc,u_tens,v_tens);
  time = omp_get_wtime() - time;
  std::cerr << "omp_time = " << std::setprecision(9) << time  << std::endl;

//  print(dom, make_host_view(u_tens));
//  print(dom, make_host_view(v_tens));
}
