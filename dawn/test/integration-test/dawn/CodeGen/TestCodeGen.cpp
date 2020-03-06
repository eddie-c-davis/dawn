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
  void fillMath(std::array<double, M * N * P>& data, double a, double b, double c, double d,
                double e, double f) {
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
  void fillValue(std::array<double, M * N * P>& storage, double val = 0.0) {
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
               std::unordered_map<std::string, std::vector<double>>& outData, const unsigned halo,
               const std::vector<std::vector<double>>& inFillValues, const double outFillValue,
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
      tokenize(lines[i], ' ', outData[outputNames[i]]);
    }
  }

  template <unsigned M, unsigned N = 1, unsigned P = 1, unsigned H = 3>
  void verify(std::array<double, M * N * P>& refData, std::vector<double>& testData,
              const double eps = 1e-6) {
    unsigned n = 0;
    for(int i = H; i < M - H; ++i) {
      for(int j = H; j < N - H; ++j) {
        for(int k = 0; k < P - 1; ++k) {
          double diff = fabs(refData[k + (j + i * M) * N] - testData[n]);
          ASSERT_LT(diff, eps) << "Test data does not match reference data at n=" << n;
          n += 1;
        }
      }
    }
  }
};

#define ndx(i, j, k) (k) + ((j) + (i)*N) * N

TEST_F(TestCodeGen, Asymmetric) {
  auto instantiation = CompilerUtil::load("input/asymmetric.iir", options_, context_);

  const unsigned N = 12;
  const unsigned halo = 3;
  const unsigned size = N * N * (N + 1);

  const int iMax = N - halo - 1;
  const int jMax = N - halo - 1;
  const int kMax = N - 1;

  std::array<double, size> inData;
  std::array<double, size> refData;
  std::array<double, size> tmpData;
  std::unordered_map<std::string, std::vector<double>> outData({{"out", std::vector<double>()}});

  // Initialize data
  fillMath<N, N, N + 1>(inData, 8.0, 2.0, 1.5, 1.5, 2.0, 4.0);
  fillValue<N, N, N + 1>(refData, -1.0);
  fillValue<N, N, N + 1>(tmpData, 0.0);

  // Populate reference data...
  int k = 0;
  for(int i = halo; i <= iMax; ++i) {
    for(int j = halo; j <= jMax; ++j) {
      tmpData[ndx(i, j, k)] = inData[ndx(i + 1, j, k)] + inData[ndx(i, j - 2, k)];
      refData[ndx(i, j, k)] = 1;
    }
  }

  for(k = 1; k <= kMax; ++k) {
    for(int i = halo; i <= iMax; ++i) {
      for(int j = halo; j <= jMax; ++j) {
        tmpData[ndx(i, j, k)] = inData[ndx(i - 3, j, k)] + inData[ndx(i, j + 1, k)];
        refData[ndx(i, j, k)] = tmpData[ndx(i, j, k - 1)];
      }
    }
  }

  // Run the generated code
  runTest<N, N, N + 1>(instantiation, outData, halo, {{8.0, 2.0, 1.5, 1.5, 2.0, 4.0}}, -1.0);

  // Verify data
  verify<N, N, N + 1, halo>(refData, outData["out"]);
}

TEST_F(TestCodeGen, ConditionalStencil) {
  auto instantiation = CompilerUtil::load("input/conditional_stencil.iir", options_, context_);

  const unsigned N = 12;
  const unsigned halo = 3;
  const unsigned size = N * N * (N + 1);

  std::array<double, size> inData;
  std::array<double, size> refData;
  std::unordered_map<std::string, std::vector<double>> outData({{"out", std::vector<double>()}});

  // Initialize data
  fillMath<N, N, N + 1>(inData, 8.0, 2.0, 1.5, 1.5, 2.0, 4.0);
  fillValue<N, N, N + 1>(refData, -1.0);

  // Populate reference data...
  for(int k = 0; k <= N - 1; ++k) {
    for(int i = halo; i <= N - halo - 1; ++i) {
      for(int j = halo; j <= N - halo - 1; ++j) {
        refData[ndx(i, j, k)] = inData[ndx(i, j + 1, k)];
      }
    }
  }

  // Per issue #724 (MeteoSwiss-APN/dawn/issues/724), globals are not
  //   deserialized properly, so need to correct before running the test...
  instantiation->getIIR()->insertGlobalVariable("var1", sir::Global(int(1)));
  instantiation->getIIR()->insertGlobalVariable("var2", sir::Global(true));

  // Run the generated code
  runTest<N, N, N + 1>(instantiation, outData, halo, {{8.0, 2.0, 1.5, 1.5, 2.0, 4.0}}, -1.0);

  // Verify data
  verify<N, N, N + 1, halo>(refData, outData["out"]);
}

TEST_F(TestCodeGen, CoriolisStencil) {
  const unsigned N = 12;
  const unsigned halo = 3;
  const unsigned size = N * N * (N + 1);

  std::array<double, size> u_nnow, v_nnow, fc;
  std::array<double, size> u_ref, v_ref;
  std::vector<double> u_tens, v_tens;
  std::unordered_map<std::string, std::vector<double>> outData(
      {{"u_tens", u_tens}, {"v_tens", v_tens}});

  // Initialize data
  fillValue<N, N, N + 1>(u_ref, -1.0);
  fillValue<N, N, N + 1>(v_ref, -1.0);
  fillMath<N, N, N + 1>(u_nnow, 8.0, 2.0, 1.5, 1.5, 2.0, 4.0);
  fillMath<N, N, N + 1>(v_nnow, 5.0, 1.2, 1.3, 1.7, 2.2, 3.5);
  fillMath<N, N, N + 1>(fc, 2.0, 1.3, 1.4, 1.6, 2.1, 3.0);

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

  // Deserialize the IIR
  auto instantiation = CompilerUtil::load("input/coriolis_stencil.iir", options_, context_);

  // Run the generated code
  runTest<N, N, N + 1>(instantiation, outData, halo,
                       {{8.0, 2.0, 1.5, 1.5, 2.0, 4.0},
                        {5.0, 1.2, 1.3, 1.7, 2.2, 3.5},
                        {2.0, 1.3, 1.4, 1.6, 2.1, 3.0}},
                       -1.0, {"u_nnow", "v_nnow", "fc"}, {"u_tens", "v_tens"});

  // Verify data
  verify<N, N, N + 1, halo>(u_ref, outData["u_tens"]);
  verify<N, N, N + 1, halo>(v_ref, outData["v_tens"]);

  if(CompilerUtil::hasCudaGPU()) {
    runTest<N, N, N + 1>(instantiation, outData, halo,
                         {{8.0, 2.0, 1.5, 1.5, 2.0, 4.0},
                          {5.0, 1.2, 1.3, 1.7, 2.2, 3.5},
                          {2.0, 1.3, 1.4, 1.6, 2.1, 3.0}},
                         -1.0, {"u_nnow", "v_nnow", "fc"}, {"u_tens", "v_tens"}, "",
                         BackendType::CUDA);
  }
}

} // anonymous namespace
