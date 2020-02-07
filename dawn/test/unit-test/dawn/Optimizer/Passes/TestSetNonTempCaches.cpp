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

#include "dawn/Compiler/DawnCompiler.h"
#include "dawn/Compiler/Options.h"
#include "dawn/IIR/IIR.h"
#include "dawn/IIR/StencilInstantiation.h"
#include "dawn/Optimizer/PassSetNonTempCaches.h"
#include "dawn/Serialization/IIRSerializer.h"
#include "dawn/Unittest/CompilerUtil.h"
#include "test/unit-test/dawn/Optimizer/TestEnvironment.h"

#include <fstream>
#include <gtest/gtest.h>

using namespace dawn;

namespace {

class TestSetNonTempCaches : public ::testing::Test {
protected:
  dawn::OptimizerContext::OptimizerContextOptions options_;
  std::unique_ptr<OptimizerContext> context_;

  virtual void SetUp() { options_.UseNonTempCaches = true; }

  void runTest(const std::string& filename, const std::vector<std::string>& cacheNames) {
    std::shared_ptr<iir::StencilInstantiation> instantiation =
        CompilerUtil::load(filename, options_, context_, TestEnvironment::path_);

    // Expect pass to succeed...
    PassSetNonTempCaches pass(*context_);
    ASSERT_TRUE(pass.run(instantiation));
    ASSERT_EQ(cacheNames, pass.getCachedFieldNames());
  }
};

TEST_F(TestSetNonTempCaches, NoCache1) { runTest("NoCache_1.sir", {}); }

TEST_F(TestSetNonTempCaches, NoCache2) { runTest("NoCache_2.sir", {}); }

TEST_F(TestSetNonTempCaches, MultipleCaches1) {
  runTest("MultipleCaches_1.sir", {"field_b", "field_c"});
}

} // anonymous namespace