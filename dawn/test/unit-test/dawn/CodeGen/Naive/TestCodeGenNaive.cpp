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

//#include "driver-includes/gridtools_includes.hpp"
//#include "driver-includes/storage_runtime.hpp"

#include <fstream>
#include <gtest/gtest.h>

using namespace dawn;
//using namespace gridtools::dawn;

namespace {

class TestCodeGenNaive : public ::testing::Test {
protected:
  dawn::OptimizerContext::OptimizerContextOptions options_;
  std::unique_ptr<OptimizerContext> context_;
  const int domainSize_ = 12;

  virtual void SetUp() { dawn::UIDGenerator::getInstance()->reset(); }

  void genTest(std::shared_ptr<iir::StencilInstantiation>& instantiation,
               const std::string& gen_file = "", const std::string& ref_file = "") {

    // Run prerequisite groups
//    ASSERT_TRUE(CompilerUtil::runGroup(PassGroup::Parallel, context_, instantiation));
//    ASSERT_TRUE(CompilerUtil::runGroup(PassGroup::ReorderStages, context_, instantiation));
//    ASSERT_TRUE(CompilerUtil::runGroup(PassGroup::MergeStages, context_, instantiation));

    // Expect codegen to succeed...
    std::string gen = CompilerUtil::generate(instantiation, gen_file);
    ASSERT_GT(gen.size(), 0);

    if(!ref_file.empty()) {
      std::string ref = dawn::readFile(ref_file);
      ASSERT_EQ(gen, ref) << "Generated code does not match reference code";
    }

    if(!gen_file.empty()) {
      std::string build_out = CompilerUtil::build(gen_file);
    }
  }
};

TEST_F(TestCodeGenNaive, Asymmetric) {
  std::string filename = "input/asymmetric.iir";
  auto instantiation = CompilerUtil::load(filename, options_, context_);
  genTest(instantiation, "output/asymmetric.cpp");
  std::array<double, 432> inData{2.46927,2.47927,2.48927,2.49927,2.50927,2.51927,2.52927,2.53927,2.54927,2.55927,2.56927,2.57927,0.738027,0.748027,0.758027,0.768027,0.778027,0.788027,0.798027,0.808027,0.818027,0.828027,0.838027,0.848027,0.304482,0.314482,0.324482,0.334482,0.344482,0.354482,0.364482,0.374482,0.384482,0.394482,0.404482,0.414482,1.23463,1.24463,1.25463,1.26463,1.27463,1.28463,1.29463,1.30463,1.31463,1.32463,1.33463,1.34463,2.91761,2.92761,2.93761,2.94761,2.95761,2.96761,2.97761,2.98761,2.99761,3.00761,3.01761,3.02761,4.43355,4.44355,4.45355,4.46355,4.47355,4.48355,4.49355,4.50355,4.51355,4.52355,4.53355,4.54355,1.36826,1.37826,1.38826,1.39826,1.40826,1.41826,1.42826,1.43826,1.44826,1.45826,1.46826,1.47826,0.267949,0.277949,0.287949,0.297949,0.307949,0.317949,0.327949,0.337949,0.347949,0.357949,0.367949,0.377949,0.602897,0.612897,0.622897,0.632897,0.642897,0.652897,0.662897,0.672897,0.682897,0.692897,0.702897,0.712897,2.06815,2.07815,2.08815,2.09815,2.10815,2.11815,2.12815,2.13815,2.14815,2.15815,2.16815,2.17815,3.82751,3.83751,3.84751,3.85751,3.86751,3.87751,3.88751,3.89751,3.90751,3.91751,3.92751,3.93751,5,5.01,5.02,5.03,5.04,5.05,5.06,5.07,5.08,5.09,5.1,5.11,0.565534,0.575534,0.585534,0.595534,0.605534,0.615534,0.625534,0.635534,0.645534,0.655534,0.665534,0.675534,0.220389,0.230389,0.240389,0.250389,0.260389,0.270389,0.280389,0.290389,0.300389,0.310389,0.320389,0.330389,1.25174,1.26174,1.27174,1.28174,1.29174,1.30174,1.31174,1.32174,1.33174,1.34174,1.35174,1.36174,3.03332,3.04332,3.05332,3.06332,3.07332,3.08332,3.09332,3.10332,3.11332,3.12332,3.13332,3.14332,4.63024,4.64024,4.65024,4.66024,4.67024,4.68024,4.69024,4.70024,4.71024,4.72024,4.73024,4.74024,5.33012,5.34012,5.35012,5.36012,5.37012,5.38012,5.39012,5.40012,5.41012,5.42012,5.43012,5.44012,0.152241,0.162241,0.172241,0.182241,0.192241,0.202241,0.212241,0.222241,0.232241,0.242241,0.252241,0.262241,0.585786,0.595786,0.605786,0.615786,0.625786,0.635786,0.645786,0.655786,0.665786,0.675786,0.685786,0.695786,2.15224,2.16224,2.17224,2.18224,2.19224,2.20224,2.21224,2.22224,2.23224,2.24224,2.25224,2.26224,4,4.01,4.02,4.03,4.04,4.05,4.06,4.07,4.08,4.09,4.1,4.11,5.23463,5.24463,5.25463,5.26463,5.27463,5.28463,5.29463,5.30463,5.31463,5.32463,5.33463,5.34463,5.41421,5.42421,5.43421,5.44421,5.45421,5.46421,5.47421,5.48421,5.49421,5.50421,5.51421,5.52421,0.169351,0.179351,0.189351,0.199351,0.209351,0.219351,0.229351,0.239351,0.249351,0.259351,0.269351,0.279351,1.30278,1.31278,1.32278,1.33278,1.34278,1.35278,1.36278,1.37278,1.38278,1.39278,1.40278,1.41278,3.17866,3.18866,3.19866,3.20866,3.21866,3.22866,3.23866,3.24866,3.25866,3.26866,3.27866,3.28866,4.84776,4.85776,4.86776,4.87776,4.88776,4.89776,4.90776,4.91776,4.92776,4.93776,4.94776,4.95776,5.58671,5.59671,5.60671,5.61671,5.62671,5.63671,5.64671,5.65671,5.66671,5.67671,5.68671,5.69671,5.283,5.293,5.303,5.313,5.323,5.333,5.343,5.353,5.363,5.373,5.383,5.393,0.602897,0.612897,0.622897,0.632897,0.642897,0.652897,0.662897,0.672897,0.682897,0.692897,0.702897,0.712897,2.26795,2.27795,2.28795,2.29795,2.30795,2.31795,2.32795,2.33795,2.34795,2.35795,2.36795,2.37795,4.19669,4.20669,4.21669,4.22669,4.23669,4.24669,4.25669,4.26669,4.27669,4.28669,4.29669,4.30669,5.48236,5.49236,5.50236,5.51236,5.52236,5.53236,5.54236,5.55236,5.56236,5.57236,5.58236,5.59236,5.67527,5.68527,5.69527,5.70527,5.71527,5.72527,5.73527,5.74527,5.75527,5.76527,5.77527,5.78527,5,5.01,5.02,5.03,5.04,5.05,5.06,5.07,5.08,5.09,5.1,5.11};


//  domain dom(domainSize_, domainSize_, domainSize_);
//  dom.set_halos(halo::value, halo::value, halo::value, halo::value, 0, 0);
//
//  meta_data_t meta_data(dom.isize(), dom.jsize(), dom.ksize() + 1);
//  storage_t in(meta_data, "in"), out(meta_data, "out");


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

  genTest(instantiation, "", "reference/global_indexing.cpp");
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

  genTest(instantiation, "", "reference/nonoverlapping_stencil.cpp");
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
