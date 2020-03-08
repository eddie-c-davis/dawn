//===--------------------------------------------------------------------------------*- C++ -*-===//
//                          _
//                         | |
//                       __| | __ ___      ___ ___
//                      / _` |/ _` \ \ /\ / / '_  |
//                     | (_| | (_| |\ V  V /| | | |
//                      \__,_|\__,_| \_/\_/ |_| |_| - Compiler Toolchain
//
//
//  This file is distributed under the MIT License (MIT).
//  See LICENSE.txt for details.
//
//===------------------------------------------------------------------------------------------===//

#include "dawn/CodeGen/CodeGen.h"
#include "dawn/Compiler/DawnCompiler.h"
#include "dawn/Compiler/Options.h"
#include "dawn/Optimizer/OptimizerContext.h"
#include "dawn/SIR/SIR.h"
#include "dawn/Support/DiagnosticsEngine.h"
#include "dawn/Support/FileSystem.h"
#include "dawn/Support/FileUtil.h"
#include "dawn/Support/StringUtil.h"
#include "dawn/Unittest/CompilerUtil.h"
#include "dawn/Unittest/IIRBuilder.h"

#include <fstream>
#include <gtest/gtest.h>

#define GRID_SIZE 12
#ifndef GRIDTOOLS_DAWN_HALO_EXTENT
#define GRIDTOOLS_DAWN_HALO_EXTENT 3
#endif

using namespace dawn;

namespace {

class TestCodeGen : public ::testing::Test {
protected:
  dawn::OptimizerContext::OptimizerContextOptions options_;
  std::unique_ptr<OptimizerContext> context_;

  virtual void SetUp() {
    dawn::UIDGenerator::getInstance()->reset();
    // CompilerUtil::Verbose = true;
  }

  template <unsigned M, unsigned N = 1, unsigned P = 1, unsigned D = 3>
  void fillMath(const std::vector<double>& args, std::array<double, M * N * P>& data) {
    double a = args[0];
    double b = args[1];
    double c = args[2];
    double d = args[3];
    double e = args[4];
    double f = args[5];
    double pi = std::atan(1.) * 4.;

    for(int i = 0; i < M; ++i) {
      for(int j = 0; j < N; ++j) {
        double x = i / (double)M;
        double y = j / (double)N;

        for(int k = 0; k < P; ++k) {
          double val = k * 10e-3 + a * (b + cos(pi * (x + c * y)) + sin(d * pi * (x + e * y))) / f;
          data[k + (j + i * M) * N] = val;
        }
      }
    }
  }

  template <unsigned M, unsigned N = 1, unsigned P = 1, unsigned D = 3>
  void fillValue(std::array<double, M * N * P>& storage, double val = -1.0) {
    for(int i = 0; i < M; ++i) {
      for(int j = 0; j < N; ++j) {
        for(int k = 0; k < P; ++k) {
          storage[k + (j + i * M) * N] = val;
        }
      }
    }
  }

  template <unsigned M, unsigned N = 1, unsigned P = 1, unsigned D = 3>
  void fillValue(std::vector<std::array<double, M * N * P>>& storages, double val = 0.0) {
    for(const auto& storage : storages) {
      fillValue<M, N, P, D>(storage, val);
    }
  }

  template <unsigned M, unsigned N = 1, unsigned P = 1>
  void verify(std::array<double, M * N * P>& refData, std::vector<double>& testData,
              const unsigned halo = GRIDTOOLS_DAWN_HALO_EXTENT, const double eps = 1e-6) {
    unsigned n = 0;
    for(int i = halo; i < M - halo; ++i) {
      for(int j = halo; j < N - halo; ++j) {
        for(int k = 0; k < P - 1; ++k) {
          double diff = fabs(refData[k + (j + i * M) * N] - testData[n]);
          ASSERT_LT(diff, eps) << "Test data does not match reference data at n=" << n;
          n += 1;
        }
      }
    }
  }

  void writeHeader(std::ofstream& ofs) {
    ofs << "#include \"driver-includes/verify.hpp\"\n";
    ofs << "#include <iostream>\n";
    ofs << "#include <iomanip>\n";

    // Add print function...
    ofs << "\nvoid print(const domain& dom, const gridtools::data_view<storage_ijk_t>& view) {\n";
    ofs << "  for(int i = dom.iminus(); i < std::min(int(dom.isize() - dom.iplus()), "
           "view.total_length<0>()); ++i)\n";
    ofs << "    for(int j = dom.jminus(); j < std::min(int(dom.jsize() - dom.jplus()), "
           "view.total_length<1>()); ++j)\n";
    ofs << "      for(int k = dom.kminus(); k < std::min(int(dom.ksize() - dom.kplus()), "
           "view.total_length<2>()); ++k)\n";
    ofs << "        std::cout << std::setprecision(9) << view(i, j, k) << ' ';\n";
    ofs << "  std::cout << std::endl;\n";
    ofs << "}\n";
    ofs << "\nint main() {\n";
  }

  void writeVarDecls(unsigned M, unsigned N, unsigned P, unsigned halo,
                     const std::vector<std::string>& inputNames,
                     const std::vector<std::string>& outputNames, std::ofstream& ofs) {
    ofs << "  domain dom(" << M << ", " << N << ", " << P - 1 << ");\n";
    ofs << "  dom.set_halos(" << halo << ", " << halo << ", " << halo << ", " << halo
        << ", 0, 0);\n";
    ofs << "  meta_data_t meta(" << M << ", " << N << ", " << P << ");\n";
    ofs << "  storage_t";

    char delim = ' ';
    for(const auto& input : inputNames) {
      ofs << delim << input << "(meta, \"" << input << "\")";
      delim = ',';
    }
    for(const auto& output : outputNames) {
      ofs << delim << output << "(meta, \"" << output << "\")";
      delim = ',';
    }
    ofs << ";\n";
  }

  void writeFieldInits(const std::vector<std::vector<double>>& inFillValues,
                       const double outFillValue, const std::vector<std::string>& inputNames,
                       const std::vector<std::string>& outputNames, std::ofstream& ofs) {
    ofs << "  verifier verif(dom);\n";
    char delim = '(';
    if(inFillValues.size() == 1) {
      ofs << "  verif.fillMath";
      for(const double fillValue : inFillValues[0]) {
        ofs << delim << fillValue;
        delim = ',';
      }
      for(const auto& input : inputNames) {
        ofs << delim << input;
      }
      ofs << ");\n";
    } else {
      unsigned n = 0;
      for(const auto& fillValues : inFillValues) {
        ofs << "  verif.fillMath";
        delim = '(';
        for(const double fillValue : fillValues) {
          ofs << delim << fillValue;
          delim = ',';
        }
        ofs << delim << inputNames[n] << ");\n";
        n += 1;
      }
    }

    ofs << "  verif.fill(" << outFillValue;
    delim = ',';
    for(const auto& output : outputNames) {
      ofs << delim << output;
    }
    ofs << ");\n";
  }

  void writeStencilCall(const BackendType& backend, const std::string& stencilName,
                        const std::vector<std::string>& inputNames,
                        const std::vector<std::string>& outputNames, std::ofstream& ofs) {
    ofs << "  dawn_generated::" << backend << "::" << stencilName << " stencil(dom);\n";
    ofs << "  stencil.run";

    char delim = '(';
    for(const auto& input : inputNames) {
      ofs << delim << input;
      delim = ',';
    }
    for(const auto& output : outputNames) {
      ofs << delim << output;
      delim = ',';
    }
    ofs << ");\n";

    for(const auto& output : outputNames) {
      ofs << "  print(dom, make_host_view(" << output << "));\n";
    }
    ofs << "}\n";
  }

  template <unsigned M, unsigned N = 1, unsigned P = 1, unsigned D = 3>
  void runTest(std::shared_ptr<iir::StencilInstantiation>& instantiation,
               std::unordered_map<std::string, std::vector<double>>& outData,
               std::unordered_map<std::string, std::array<double, M * N * P>>& refData,
               const std::vector<std::vector<double>>& inFillValues,
               const double outFillValue = -1.0, const unsigned halo = GRIDTOOLS_DAWN_HALO_EXTENT,
               const std::vector<std::string>& inputNames = {"in"},
               const std::vector<std::string>& outputNames = {"out"},
               const std::string& srcFile = "",
               const BackendType& backend = BackendType::CXXNaive) {
    std::string stencilName = instantiation->getMetaData().getStencilName();

    std::string gen = CompilerUtil::generate(instantiation, srcFile, backend);
    ASSERT_GT(gen.size(), 0);

    std::string genFile = srcFile;
    if(genFile.empty()) {
      genFile = fs::temp_directory_path().string() + "/" + stencilName + ".cpp";
    }

    // Create wrapper
    std::ofstream ofs(genFile);
    ofs << gen << std::endl;
    writeHeader(ofs);

    // Add main function...
    writeVarDecls(M, N, P, halo, inputNames, outputNames, ofs);
    writeFieldInits(inFillValues, outFillValue, inputNames, outputNames, ofs);
    writeStencilCall(backend, stencilName, inputNames, outputNames, ofs);

    ofs.close();

    std::string outFile = genFile;
    unsigned pos = outFile.rfind('.');
    if(pos != std::string::npos) {
      outFile = outFile.substr(0, pos);
    }

    std::string buildOut = CompilerUtil::build(genFile, outFile, backend);
    ASSERT_TRUE(buildOut.empty());
    ASSERT_TRUE(fs::exists(outFile));

    std::string output = readPipe(outFile);
    ASSERT_FALSE(output.empty());
    ASSERT_TRUE(output.find("error") == std::string::npos) << output;

    std::vector<std::string> lines;
    tokenize(output, '\n', lines);
    for(int i = 0; i < outputNames.size(); ++i) {
      outData[outputNames[i]].clear();
      tokenize(lines[i], ' ', outData[outputNames[i]]);
    }

    // Verify data
    for(const std::string& outputName : outputNames) {
      verify<M, N, P>(refData[outputName], outData[outputName]);
    }
  }
};

#define ndx(i, j, k) (k) + ((j) + (i)*N) * N

TEST_F(TestCodeGen, Asymmetric) {
  auto instantiation = CompilerUtil::load("input/asymmetric.iir", options_, context_);

  const unsigned N = GRID_SIZE;
  const unsigned halo = GRIDTOOLS_DAWN_HALO_EXTENT;
  const unsigned size = N * N * (N + 1);

  const int iMax = N - halo - 1;
  const int jMax = N - halo - 1;
  const int kMax = N - 1;

  std::array<double, size> inData;
  std::array<double, size> outRef;
  std::array<double, size> tmpData;
  std::unordered_map<std::string, std::vector<double>> outData({{"out", std::vector<double>()}});
  std::unordered_map<std::string, std::array<double, size>> refData;

  // Initialize data
  std::vector<double> inFills = {8.0, 2.0, 1.5, 1.5, 2.0, 4.0};
  fillMath<N, N, N + 1>(inFills, inData);
  fillValue<N, N, N + 1>(outRef, -1.0);
  fillValue<N, N, N + 1>(tmpData, 0.0);

  // Populate reference data...
  int k = 0;
  for(int i = halo; i <= iMax; ++i) {
    for(int j = halo; j <= jMax; ++j) {
      tmpData[ndx(i, j, k)] = inData[ndx(i + 1, j, k)] + inData[ndx(i, j - 2, k)];
      outRef[ndx(i, j, k)] = 1;
    }
  }

  for(k = 1; k <= kMax; ++k) {
    for(int i = halo; i <= iMax; ++i) {
      for(int j = halo; j <= jMax; ++j) {
        tmpData[ndx(i, j, k)] = inData[ndx(i - 3, j, k)] + inData[ndx(i, j + 1, k)];
        outRef[ndx(i, j, k)] = tmpData[ndx(i, j, k - 1)];
      }
    }
  }

  refData["out"] = outRef;

  // Run the generated code
  runTest<N, N, N + 1>(instantiation, outData, refData, {inFills});

  if(CompilerUtil::hasCudaGPU()) {
    runTest<N, N, N + 1>(instantiation, outData, refData, {inFills}, -1.0, halo, {"in"}, {"out"},
                         "", BackendType::CUDA);
  }
}

TEST_F(TestCodeGen, ConditionalStencil) {
  auto instantiation = CompilerUtil::load("input/conditional_stencil.iir", options_, context_);

  const unsigned N = GRID_SIZE;
  const unsigned halo = GRIDTOOLS_DAWN_HALO_EXTENT;
  const unsigned size = N * N * (N + 1);

  std::array<double, size> inData;
  std::array<double, size> outRef;
  std::unordered_map<std::string, std::vector<double>> outData({{"out", std::vector<double>()}});
  std::unordered_map<std::string, std::array<double, size>> refData;

  // Initialize data
  std::vector<double> inFills = {8.0, 2.0, 1.5, 1.5, 2.0, 4.0};
  fillMath<N, N, N + 1>(inFills, inData);
  fillValue<N, N, N + 1>(outRef);

  // Populate reference data...
  for(int k = 0; k <= N - 1; ++k) {
    for(int i = halo; i <= N - halo - 1; ++i) {
      for(int j = halo; j <= N - halo - 1; ++j) {
        outRef[ndx(i, j, k)] = inData[ndx(i, j + 1, k)];
      }
    }
  }

  refData["out"] = outRef;

  // Per issue #724 (MeteoSwiss-APN/dawn/issues/724), globals are not
  //   deserialized properly, so need to correct before running the test...
  instantiation->getIIR()->insertGlobalVariable("var1", sir::Global(int(1)));
  instantiation->getIIR()->insertGlobalVariable("var2", sir::Global(true));

  // Run the generated code
  runTest<N, N, N + 1>(instantiation, outData, refData, {inFills});

  if(CompilerUtil::hasCudaGPU()) {
    runTest<N, N, N + 1>(instantiation, outData, refData, {inFills}, -1.0, halo, {"in"}, {"out"},
                         "", BackendType::CUDA);
  }
}

TEST_F(TestCodeGen, CoriolisStencil) {
  const unsigned N = GRID_SIZE;
  const unsigned halo = GRIDTOOLS_DAWN_HALO_EXTENT;
  const unsigned size = N * N * (N + 1);

  std::array<double, size> u_nnow, v_nnow, fc;
  std::array<double, size> u_ref, v_ref;
  std::vector<double> u_tens, v_tens;
  std::unordered_map<std::string, std::vector<double>> outData(
      {{"u_tens", u_tens}, {"v_tens", v_tens}});
  std::unordered_map<std::string, std::array<double, size>> refData;

  // Initialize data
  fillValue<N, N, N + 1>(u_ref);
  fillValue<N, N, N + 1>(v_ref);
  std::vector<double> u_fill = {8.0, 2.0, 1.5, 1.5, 2.0, 4.0};
  fillMath<N, N, N + 1>(u_fill, u_nnow);
  std::vector<double> v_fill = {5.0, 1.2, 1.3, 1.7, 2.2, 3.5};
  fillMath<N, N, N + 1>(v_fill, v_nnow);
  std::vector<double> fc_fill = {2.0, 1.3, 1.4, 1.6, 2.1, 3.0};
  fillMath<N, N, N + 1>(fc_fill, fc);

  // Populate reference data...
  for(int k = 0; k <= N - 1; ++k) {
    for(int i = halo; i <= N - halo - 1; ++i) {
      for(int j = halo; j <= N - halo - 1; ++j) {
        double z_fv_north = (fc[ndx(i, j, k)] * (v_nnow[ndx(i, j, k)] + v_nnow[ndx(i + 1, j, k)]));
        double z_fv_south =
            fc[ndx(i, j - 1, k)] * (v_nnow[ndx(i, j - 1, k)] + v_nnow[ndx(i + 1, j - 1, k)]);
        u_ref[ndx(i, j, k)] += 0.25 * (z_fv_north + z_fv_south);
        double z_fu_east = fc[ndx(i, j, k)] * (u_nnow[ndx(i, j, k)] + u_nnow[ndx(i, j + 1, k)]);
        double z_fu_west =
            (fc[ndx(i - 1, j, k)] * (u_nnow[ndx(i - 1, j, k)] + u_nnow[ndx(i - 1, j + 1, k)]));
        v_ref[ndx(i, j, k)] -= 0.25 * (z_fu_east + z_fu_west);
      }
    }
  }

  refData["u_tens"] = u_ref;
  refData["v_tens"] = v_ref;

  // Deserialize the IIR
  auto instantiation = CompilerUtil::load("input/coriolis_stencil.iir", options_, context_);

  // Run the generated code
  runTest<N, N, N + 1>(instantiation, outData, refData, {u_fill, v_fill, fc_fill}, -1.0, halo,
                       {"u_nnow", "v_nnow", "fc"}, {"u_tens", "v_tens"});

  if(CompilerUtil::hasCudaGPU()) {
    runTest<N, N, N + 1>(instantiation, outData, refData, {u_fill, v_fill, fc_fill}, -1.0, halo,
                         {"u_nnow", "v_nnow", "fc"}, {"u_tens", "v_tens"}, "", BackendType::CUDA);
  }
}

inline bool checkOffset(unsigned int min, unsigned int max, unsigned int val) {
  return (min <= val && val < max);
}

TEST_F(TestCodeGen, YPPMStencil) {
  const unsigned N = GRID_SIZE;
  const unsigned halo = GRIDTOOLS_DAWN_HALO_EXTENT;
  const unsigned size = N * N * (N + 1);

  const double c1 = 0.142857;
  const double c2 = 0.785714;
  const double c3 = 0.357143;
  const double p1 = 0.583333;
  const double p2 = 0.083333;

  std::array<double, size> q, dya;
  std::array<double, size> al_ref;
  std::vector<double> al_out;
  std::unordered_map<std::string, std::vector<double>> outData({{"al", al_out}});
  std::unordered_map<std::string, std::array<double, size>> refData;

  // Initialize data
  fillValue<N, N, N + 1>(al_ref);
  std::vector<double> q_fill = {8.0, 2.0, 1.5, 1.5, 2.0, 4.0};
  fillMath<N, N, N + 1>(q_fill, q);
  std::vector<double> dy_fill = {5.0, 1.2, 1.3, 1.7, 2.2, 3.5};
  fillMath<N, N, N + 1>(dy_fill, dya);

  std::array<unsigned, 2> stage1GlobalJIndices = {3, 4};
  std::array<unsigned, 2> stage2GlobalJIndices = {8, 9};

  // Populate reference data...
  for(int k = 0; k <= N - 1; ++k) {
    for(int i = halo; i <= N - halo - 1; ++i) {
      for(int j = halo; j <= N - halo - 1; ++j) {
        al_ref[ndx(i, j, k)] = (p1 * (q[ndx(i, j - 1, k)] + q[ndx(i, j, k)])) -
                               (p2 * (q[ndx(i, j - 2, k)] + q[ndx(i, j + 1, k)]));
      }
    }
    for(int i = halo; i <= N - halo - 1; ++i) {
      for(int j = halo; j <= N - halo - 1; ++j) {
        if(checkOffset(stage1GlobalJIndices[0], stage1GlobalJIndices[1], j)) {
          al_ref[ndx(i, j, k)] =
              -c1 * q[ndx(i, j - 2, k)] + c2 * q[ndx(i, j - 1, k)] + c3 * q[ndx(i, j, k)];
          al_ref[ndx(i, j + 1, k)] =
              (0.5 *
               (((((2.0 * dya[ndx(i, j, k)] + dya[ndx(i, j - 1, k)]) * q[ndx(i, j, k)]) -
                  dya[ndx(i, j, k)] * q[ndx(i, j - 1, k)]) /
                 (dya[ndx(i, j - 1, k)] + dya[ndx(i, j, k)])) +
                ((((2.0 * dya[ndx(i, j + 1, k)] + dya[ndx(i, j + 2, k)]) * q[ndx(i, j + 1, k)]) -
                  dya[ndx(i, j + 1, k)] * q[ndx(i, j + 2, k)]) /
                 (dya[ndx(i, j + 1, k)] + dya[ndx(i, j + 2, k)]))));
          al_ref[ndx(i, j + 2, k)] =
              ((c3 * q[ndx(i, j + 1, k)] + c2 * q[ndx(i, j + 2, k)]) - c1 * q[ndx(i, j + 3, k)]);
        }
      }
    }
    for(int i = halo; i <= N - halo - 1; ++i) {
      for(int j = halo; j <= N - halo - 1; ++j) {
        if(checkOffset(stage2GlobalJIndices[0], stage2GlobalJIndices[1], j)) {
          al_ref[ndx(i, j, k)] =
              (-c1 * q[ndx(i, j - 2, k)] + c2 * q[ndx(i, j - 1, k)]) + (c3 * q[ndx(i, j, k)]);
          al_ref[ndx(i, j + 1, k)] =
              0.5 *
              (((((2.0 * dya[ndx(i, j, k)] + dya[ndx(i, j - 1, k)]) * q[ndx(i, j, k)]) -
                 dya[ndx(i, j, k)] * q[ndx(i, j - 1, k)]) /
                (dya[ndx(i, j - 1, k)] + dya[ndx(i, j, k)])) +
               ((((2.0 * dya[ndx(i, j + 1, k)] + dya[ndx(i, j + 2, k)]) * q[ndx(i, j + 1, k)]) -
                 (dya[ndx(i, j + 1, k)] * q[ndx(i, j + 2, k)])) /
                (dya[ndx(i, j + 1, k)] + dya[ndx(i, j + 2, k)])));
          al_ref[ndx(i, j + 2, k)] =
              (c3 * q[ndx(i, j + 1, k)] + c2 * q[ndx(i, j + 2, k)]) - (c1 * q[ndx(i, j + 3, k)]);
        }
      }
    }
  }

  refData["al"] = al_ref;

  // Deserialize the IIR
  // TOOD: Appears global indices (iteration spaces) are not serialized in the IIR, not sure
  //       if in the SIR for that matter, need to follow it through the tool chain...
  auto instantiation = CompilerUtil::load("input/yppm.iir", options_, context_);

  // Run the generated code
  runTest<N, N, N + 1>(instantiation, outData, refData, {q_fill, dy_fill}, -1.0, halo, {"q", "dya"},
                       {"al"});

  if(CompilerUtil::hasCudaGPU()) {
    runTest<N, N, N + 1>(instantiation, outData, refData, {q_fill, dy_fill}, -1.0, halo,
                         {"q", "dya"}, {"al"}, "", BackendType::CUDA);
  }
}

} // anonymous namespace
