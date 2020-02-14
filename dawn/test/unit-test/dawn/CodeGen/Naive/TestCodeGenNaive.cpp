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

  void genTest(std::shared_ptr<iir::StencilInstantiation>& instantiation,
               const std::string& ref_file = "") {
//    std::shared_ptr<iir::StencilInstantiation> instantiation =
//        CompilerUtil::load(filename, options_, context_);

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
};

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
