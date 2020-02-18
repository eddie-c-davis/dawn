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
#include "dawn/Unittest/CompilerUtil.h"
#include "dawn/Unittest/IIRBuilder.h"
#include "dawn/Unittest/UnittestLogger.h"

#include <fstream>
#include <gtest/gtest.h>

using namespace dawn;

namespace {

class TestCodeGenNaive : public ::testing::Test {
protected:
  dawn::OptimizerContext::OptimizerContextOptions options_;
  std::unique_ptr<OptimizerContext> context_;

  virtual void SetUp() { dawn::UIDGenerator::getInstance()->reset(); }

  template<unsigned M, unsigned N = 1, unsigned P = 1, unsigned D = 3>
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

  template<unsigned M, unsigned N = 1, unsigned P = 1, unsigned D = 3>
  void fillValue(std::array<double, M * N * P>& data, double val = 0.0) {
    for(int i = 0; i < M; ++i) {
      for(int j = 0; j < N; ++j) {
        for(int k = 0; k < P; ++k) {
          data[k + (j + i * M) * N] = val;
        }
      }
    }
  }

  void genTest(std::shared_ptr<iir::StencilInstantiation>& instantiation,
               const std::string& ref_file = "") {

    // Run prerequisite groups
    //    ASSERT_TRUE(CompilerUtil::runGroup(PassGroup::Parallel, context_, instantiation));
    //    ASSERT_TRUE(CompilerUtil::runGroup(PassGroup::ReorderStages, context_, instantiation));
    //    ASSERT_TRUE(CompilerUtil::runGroup(PassGroup::MergeStages, context_, instantiation));

    // Expect codegen to succeed...
    std::string gen = CompilerUtil::generate(instantiation);
    ASSERT_GT(gen.size(), 0);

    if(!ref_file.empty()) {
      std::string ref = dawn::readFile(ref_file);
      ASSERT_EQ(gen, ref) << "Generated code does not match reference code";
    }
  }

  template <unsigned M, unsigned N = 1, unsigned P = 1, unsigned D = 3>
  void runTest(std::shared_ptr<iir::StencilInstantiation>& instantiation,
               std::array<double, M * N * P>& data, const std::string& gen_file) {
    std::string gen = CompilerUtil::generate(instantiation, gen_file);
    ASSERT_GT(gen.size(), 0);

    std::string out_file = gen_file;
    unsigned pos = out_file.rfind('.');
    if(pos != std::string::npos) {
      out_file = out_file.substr(0, pos);
    }

    std::string build_out = CompilerUtil::build(gen_file, out_file);
    ASSERT_TRUE(build_out.empty());
    ASSERT_FALSE(out_file.empty());
    ASSERT_TRUE(fs::exists(out_file));
  }

  template <unsigned M, unsigned N = 1, unsigned P = 1, unsigned D = 3>
  void verify(std::array<double, M * N * P>& refData, std::array<double, M * N * P>& testData,
              const double eps = 10e-6) {
    for(int n = 0; n < M * N * P; ++n) {
      double diff = fabs(refData[n] - testData[n]);
      ASSERT_LT(diff, eps) << "Test data does not match reference data at n=" << n;
    }
  }
};

TEST_F(TestCodeGenNaive, Asymmetric) {
  std::string filename = "input/asymmetric.iir";
  auto instantiation = CompilerUtil::load(filename, options_, context_);

  const unsigned N = 12;
  const unsigned halo = N / 4;
  const unsigned size = N * N * (N + 1);

  const int iMax = N - halo - 1;
  const int jMax = N - halo - 1;
  const int kMax = N - 1;

  std::array<double, size> inData;
  std::array<double, size> outData;
  std::array<double, size> refData;
  std::array<double, size> tmpData;

  // Initialize data
  fillMath<N, N, N + 1>(inData, 8.0, 2.0, 1.5, 1.5, 2.0, 4.0);
  fillValue<N, N, N + 1>(outData, -1.0);
  fillValue<N, N, N + 1>(refData, -1.0);
  fillValue<N, N, N + 1>(tmpData, 0.0);

  // Populate reference data...
  int k = 0;
  for(int i = halo; i <= iMax; ++i) {
    for(int j = halo; j <= jMax; ++j) {
      tmpData[k + (j + i * N) * N] =
          inData[k + (j + (i + 1) * N) * N] + inData[k + ((j - 2) + i * N) * N];
      refData[k + (j + i * N) * N] = 1;
    }
  }

  for(k = 1; k <= kMax; ++k) {
    for(int i = halo; i <= iMax; ++i) {
      for(int j = halo; j <= jMax; ++j) {
        tmpData[k + (j + i * N) * N] =
            inData[k + (j + (i - 3) * N) * N] + inData[k + ((j + 1) + i * N) * N];
        refData[k + (j + i * N) * N] = tmpData[(k - 1) + (j + i * N) * N];
      }
    }
  }

  // Run the generated code
  runTest<N, N, N + 1>(instantiation, outData, "output/asymmetric.cpp");
}

TEST_F(TestCodeGenNaive, GlobalIndexStencil) {
  using namespace dawn::iir;

  CartesianIIRBuilder b;
  auto in_f = b.field("in_field", FieldType::ijk);
  auto out_f = b.field("out_field", FieldType::ijk);

  auto instantiation =
      b.build("generated",
              b.stencil(b.multistage(
                  LoopOrderKind::Parallel,
                  b.stage(b.doMethod(dawn::sir::Interval::Start, dawn::sir::Interval::End,
                                     b.block(b.stmt(b.assignExpr(b.at(out_f), b.at(in_f)))))),
                  b.stage(1, {0, 2},
                          b.doMethod(dawn::sir::Interval::Start, dawn::sir::Interval::End,
                                     b.block(b.stmt(b.assignExpr(b.at(out_f), b.lit(10)))))))));

  genTest(instantiation, "reference/global_indexing.cpp");
}

TEST_F(TestCodeGenNaive, NonOverlappingInterval) {
  using namespace dawn::iir;
  using SInterval = dawn::sir::Interval;

  CartesianIIRBuilder b;
  auto in = b.field("in", FieldType::ijk);
  auto out = b.field("out", FieldType::ijk);
  auto dx = b.localvar("dx", dawn::BuiltinTypeID::Double);

  auto instantiation = b.build(
      "generated",
      b.stencil(b.multistage(
          LoopOrderKind::Parallel,
          b.stage(b.doMethod(
              SInterval(SInterval::Start, 10), b.declareVar(dx),
              b.block(b.stmt(b.assignExpr(
                  b.at(out),
                  b.binaryExpr(
                      b.binaryExpr(
                          b.lit(-4),
                          b.binaryExpr(
                              b.at(in),
                              b.binaryExpr(b.at(in, {1, 0, 0}),
                                           b.binaryExpr(b.at(in, {-1, 0, 0}),
                                                        b.binaryExpr(b.at(in, {0, -1, 0}),
                                                                     b.at(in, {0, 1, 0}))))),
                          Op::multiply),
                      b.binaryExpr(b.at(dx), b.at(dx), Op::multiply), Op::divide)))))),
          b.stage(b.doMethod(SInterval(15, SInterval::End),
                             b.block(b.stmt(b.assignExpr(b.at(out), b.lit(10)))))))));

  genTest(instantiation, "reference/nonoverlapping_stencil.cpp");
}

TEST_F(TestCodeGenNaive, LaplacianStencil) {
  using namespace dawn::iir;
  using SInterval = dawn::sir::Interval;

  CartesianIIRBuilder b;
  auto in = b.field("in", FieldType::ijk);
  auto out = b.field("out", FieldType::ijk);
  auto dx = b.localvar("dx", dawn::BuiltinTypeID::Double);

  auto instantiation = b.build(
      "generated",
      b.stencil(b.multistage(
          LoopOrderKind::Parallel,
          b.stage(b.doMethod(
              SInterval::Start, SInterval::End, b.declareVar(dx),
              b.block(b.stmt(b.assignExpr(
                  b.at(out),
                  b.binaryExpr(
                      b.binaryExpr(
                          b.lit(-4),
                          b.binaryExpr(
                              b.at(in),
                              b.binaryExpr(b.at(in, {1, 0, 0}),
                                           b.binaryExpr(b.at(in, {-1, 0, 0}),
                                                        b.binaryExpr(b.at(in, {0, -1, 0}),
                                                                     b.at(in, {0, 1, 0}))))),
                          Op::multiply),
                      b.binaryExpr(b.at(dx), b.at(dx), Op::multiply), Op::divide)))))))));

  genTest(instantiation);
}

} // anonymous namespace
