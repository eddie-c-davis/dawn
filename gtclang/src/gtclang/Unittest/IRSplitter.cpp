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

#include "gtclang/Unittest/IRSplitter.h"
#include "dawn/CodeGen/CXXNaive/CXXNaiveCodeGen.h"
#include "dawn/CodeGen/Cuda/CudaCodeGen.h"
#include "dawn/Optimizer/PassComputeStageExtents.h"
#include "dawn/Optimizer/PassDataLocalityMetric.h"
#include "dawn/Optimizer/PassFieldVersioning.h"
#include "dawn/Optimizer/PassFixVersionedInputFields.h"
#include "dawn/Optimizer/PassInlining.h"
#include "dawn/Optimizer/PassIntervalPartitioner.h"
#include "dawn/Optimizer/PassMultiStageSplitter.h"
#include "dawn/Optimizer/PassSSA.h"
#include "dawn/Optimizer/PassSetBlockSize.h"
#include "dawn/Optimizer/PassSetCaches.h"
#include "dawn/Optimizer/PassSetNonTempCaches.h"
#include "dawn/Optimizer/PassSetStageGraph.h"
#include "dawn/Optimizer/PassSetStageName.h"
#include "dawn/Optimizer/PassSetSyncStage.h"
#include "dawn/Optimizer/PassStageMerger.h"
#include "dawn/Optimizer/PassStageReordering.h"
#include "dawn/Optimizer/PassStageSplitter.h"
#include "dawn/Optimizer/PassTemporaryMerger.h"
#include "dawn/Optimizer/PassTemporaryToStencilFunction.h"
#include "dawn/Optimizer/PassTemporaryType.h"
#include "dawn/SIR/SIR.h"
#include "dawn/Serialization/IIRSerializer.h"
#include "dawn/Serialization/SIRSerializer.h"
#include "dawn/Support/DiagnosticsEngine.h"
#include "dawn/Support/FileSystem.h"
#include "dawn/Support/UIDGenerator.h"
#include "gtclang/Unittest/Config.h"
#include "gtclang/Unittest/GTClang.h"
#include <fstream>

namespace gtclang {

IRSplitter::IRSplitter(const std::string& destDir, unsigned maxLevel)
    : filePrefix_(destDir), maxLevel_(maxLevel) {}

void IRSplitter::split(const std::string& dslFile, const std::vector<std::string>& args) {
  fs::path filePath(dslFile);
  if(filePrefix_.empty())
    filePrefix_ = filePath.root_directory().string();
  filePrefix_ += "/" + filePath.stem().string();

  std::vector<std::string> flags = {"-std=c++11",
                                    std::string{"-I"} + std::string{GTCLANG_UNITTEST_INCLUDES}};
  for(const auto& arg : args) {
    flags.emplace_back(arg);
  }

  auto [success, sir] = GTClang::run({dslFile, "-fno-codegen"}, flags);

  // Serialize the SIR
  writeSIR(sir);

  // Use SIR to create context
  createContext(sir);

  // Lower to unoptimized IIR and serialize
  writeIIR();

  if(success && maxLevel_ > 0) {
    // Run parallelization passes
    parallelize();
    writeIIR(1);

    // Run optimization passes
    optimize();
  }
}

void IRSplitter::generate(const std::string& outFile) {
  std::unique_ptr<dawn::codegen::TranslationUnit> tu;
  dawn::DiagnosticsEngine diagnostics;
  auto& ctx = context_->getStencilInstantiationMap();

  if(outFile.find(".cu") != std::string::npos) {
    dawn::codegen::cuda::CudaCodeGen generator(ctx, diagnostics, 0, 0, 0, {0, 0, 0});
    tu = generator.generateCode();
  } else {
    dawn::codegen::cxxnaive::CXXNaiveCodeGen generator(ctx, diagnostics, 0);
    tu = generator.generateCode();
  }

  std::ostringstream ss;
  for(auto const& macroDefine : tu->getPPDefines())
    ss << macroDefine << "\n";

  ss << tu->getGlobals();
  for(auto const& s : tu->getStencils())
    ss << s.second;

  if(outFile.empty()) {
    std::cerr << ss.str();
  } else {
    std::ofstream ofs(outFile.c_str());
    ofs << ss.str();
  }
}

void IRSplitter::createContext(const std::shared_ptr<dawn::SIR>& sir) {
  dawn::OptimizerContext::OptimizerContextOptions options;
  context_ = std::make_unique<dawn::OptimizerContext>(diag_, options, sir);
}

void IRSplitter::parallelize() {
  using MultistageSplitStrategy = dawn::PassMultiStageSplitter::MultiStageSplittingStrategy;
  MultistageSplitStrategy mssSplitStrategy = MultistageSplitStrategy::Optimized;

  for(auto& [name, instantiation] : context_->getStencilInstantiationMap()) {
    runPass<dawn::PassInlining>(name, instantiation, true,
                                dawn::PassInlining::InlineStrategy::InlineProcedures);
    runPass<dawn::PassFieldVersioning>(name, instantiation);
    runPass<dawn::PassSSA>(name, instantiation);
    runPass<dawn::PassMultiStageSplitter>(name, instantiation, mssSplitStrategy);
    runPass<dawn::PassStageSplitter>(name, instantiation);
    runPass<dawn::PassTemporaryType>(name, instantiation);
    runPass<dawn::PassFixVersionedInputFields>(name, instantiation);
    runPass<dawn::PassComputeStageExtents>(name, instantiation);
    runPass<dawn::PassSetSyncStage>(name, instantiation);
  }
}

void IRSplitter::optimize() {
  unsigned level = 1;

  // Reorder stages
  reorderStages();
  level += 1;
  writeIIR(level);

  // Merge stages
  mergeStages();
  level += 1;
  writeIIR(level);

  // Merge temporaries
  mergeTemporaries();
  level += 1;
  writeIIR(level);

  // Next inlining step
  inlining();
  level += 1;
  writeIIR(level);

  // Interval partitioning...
  partitionIntervals();
  level += 1;
  writeIIR(level);

  // Pass temporaries to functions
  // OFF by default (dawn/Optimizer/OptimizerOptions.inc)
  //  passTmpToFunction();
  //  level += 1;
  //  writeIIR(level);

  // OFF by default (dawn/Optimizer/OptimizerOptions.inc)
  //  setNonTempCaches();
  //  level += 1;
  //  writeIIR(level);

  // OFF by default (dawn/Optimizer/OptimizerOptions.inc)
  //  setCaches();
  //  level += 1;
  //  writeIIR(level);

  // Unsure whether this is ON by default -- probably only for CudaCodeGen
  //  setBlockSize();
  //  level += 1;
  //  writeIIR(level);

  // Unsure whether this is ON by default -- diagnostics only
  //  dataLocalityMetric();
  //  level += 1;
  //  writeIIR(level);
}

void IRSplitter::reorderStages() {
  dawn::ReorderStrategy::Kind reorderStrategy = dawn::ReorderStrategy::Kind::Greedy;
  for(auto& [name, instantiation] : context_->getStencilInstantiationMap()) {
    runPass<dawn::PassSetStageGraph>(name, instantiation);
    runPass<dawn::PassStageReordering>(name, instantiation, reorderStrategy);
    runPass<dawn::PassSetSyncStage>(name, instantiation);
    runPass<dawn::PassSetStageName>(name, instantiation);
  }
}

void IRSplitter::mergeStages() {
  for(auto& [name, instantiation] : context_->getStencilInstantiationMap()) {
    runPass<dawn::PassStageMerger>(name, instantiation);
    // since this can change the scope of temporaries ...
    runPass<dawn::PassTemporaryType>(name, instantiation);
    runPass<dawn::PassFixVersionedInputFields>(name, instantiation);
    // modify stages and their extents ...
    runPass<dawn::PassComputeStageExtents>(name, instantiation);
    // and changes their dependencies
    runPass<dawn::PassSetSyncStage>(name, instantiation);
  }
}

void IRSplitter::mergeTemporaries() {
  for(auto& [name, instantiation] : context_->getStencilInstantiationMap()) {
    runPass<dawn::PassTemporaryMerger>(name, instantiation);
    // this should not affect the temporaries but since we're touching them it would probably be a
    // safe idea
    runPass<dawn::PassTemporaryType>(name, instantiation);
  }
}

void IRSplitter::inlining() {
  for(auto& [name, instantiation] : context_->getStencilInstantiationMap()) {
    runPass<dawn::PassInlining>(name, instantiation, false,
                                dawn::PassInlining::InlineStrategy::ComputationsOnTheFly);
  }
}

void IRSplitter::partitionIntervals() {
  for(auto& [name, instantiation] : context_->getStencilInstantiationMap()) {
    runPass<dawn::PassIntervalPartitioner>(name, instantiation);
    // since this can change the scope of temporaries ...
    runPass<dawn::PassTemporaryType>(name, instantiation);
    runPass<dawn::PassFixVersionedInputFields>(name, instantiation);
  }
}

void IRSplitter::passTmpToFunction() {
  for(auto& [name, instantiation] : context_->getStencilInstantiationMap()) {
    runPass<dawn::PassTemporaryToStencilFunction>(name, instantiation);
  }
}

void IRSplitter::setNonTempCaches() {
  for(auto& [name, instantiation] : context_->getStencilInstantiationMap()) {
    runPass<dawn::PassSetNonTempCaches>(name, instantiation);
    // this should not affect the temporaries but since we're touching them it would probably be a
    // safe idea
    runPass<dawn::PassTemporaryType>(name, instantiation);
  }
}

void IRSplitter::setCaches() {
  for(auto& [name, instantiation] : context_->getStencilInstantiationMap()) {
    runPass<dawn::PassSetCaches>(name, instantiation);
  }
}

void IRSplitter::setBlockSize() {
  for(auto& [name, instantiation] : context_->getStencilInstantiationMap()) {
    runPass<dawn::PassSetBlockSize>(name, instantiation);
  }
}

void IRSplitter::dataLocalityMetric() {
  for(auto& [name, instantiation] : context_->getStencilInstantiationMap()) {
    // Plain diagnostics, should not even be a pass but is independent
    runPass<dawn::PassDataLocalityMetric>(name, instantiation);
  }
}

template <class T, typename... Args>
bool IRSplitter::runPass(const std::string& name,
                         std::shared_ptr<dawn::iir::StencilInstantiation>& instantiation,
                         Args&&... args) {
  T pass(*context_, std::forward<Args>(args)...);
  return pass.run(instantiation);
}

void IRSplitter::writeSIR(const std::shared_ptr<dawn::SIR>& sir) {
  dawn::UIDGenerator::getInstance()->reset();
  std::cerr << "Writing SIR file '" << filePrefix_ << ".sir'" << std::endl;
  dawn::SIRSerializer::serialize(filePrefix_ + ".sir", sir.get());
}

void IRSplitter::writeIIR(const unsigned level) {
  if(level > maxLevel_)
    return;

  unsigned nstencils = context_->getStencilInstantiationMap().size();
  unsigned stencil_id = 0;
  for(auto& [name, instantiation] : context_->getStencilInstantiationMap()) {
    std::string iirFile = filePrefix_;
    if(nstencils > 1)
      iirFile += "." + std::to_string(stencil_id);
    if(maxLevel_ > 0)
      iirFile += ".O" + std::to_string(level);
    iirFile += ".iir";

    std::cerr << "Writing IIR file '" << iirFile << "'" << std::endl;
    dawn::IIRSerializer::serialize(iirFile, instantiation);
    stencil_id += 1;
  }
}

} // namespace gtclang

// TODO: Refactor this before PR!!!
#ifdef MAIN_ENABLED
int main(int argc, char* argv[]) {
  if(argc < 2) {
    std::cerr << "usage: " << argv[0] << " <DSL File> [Dest Dir] [Max Opt. Level]" << std::endl;
    return 1;
  }

  std::string filename{argv[1]};
  std::string dest_dir;
  if(argc > 2)
    dest_dir = std::string{argv[2]};
  unsigned max_level = 1000;
  if(argc > 3)
    max_level = atoi(argv[3]);

  gtclang::IRSplitter splitter(dest_dir, max_level);
  splitter.split(filename);

  return 0;
}
#endif
