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

#ifndef DAWN_OPTIMIZER_INTEGRITYCHECKER_H
#define DAWN_OPTIMIZER_INTEGRITYCHECKER_H

#include "dawn/IIR/ASTUtil.h"
#include "dawn/IIR/StencilInstantiation.h"
#include "dawn/Optimizer/OptimizerContext.h"

namespace dawn {

//===------------------------------------------------------------------------------------------===//
//     IntegrityChecker
//===------------------------------------------------------------------------------------------===//
/// @brief Perform basic integrity checks on the AST (e.g., ensure constant expressions
/// are not reassigned).
class IntegrityChecker : public iir::ASTVisitor {
  iir::StencilInstantiation* instantiation_;
  iir::StencilMetaInformation& metadata_;
  //OptimizerContext& context_;

public:
  IntegrityChecker(iir::StencilInstantiation* instantiation); //, OptimizerContext& context);

  void visit(const std::shared_ptr<iir::BlockStmt>& stmt) override {}

  void visit(const std::shared_ptr<iir::ExprStmt>& stmt) override {}

  void visit(const std::shared_ptr<iir::ReturnStmt>& stmt) override {}

  void visit(const std::shared_ptr<iir::IfStmt>& stmt) override {}

  void visit(const std::shared_ptr<iir::VarDeclStmt>& stmt) override {}

  void visit(const std::shared_ptr<iir::VerticalRegionDeclStmt>& stmt) override {}

  void visit(const std::shared_ptr<iir::StencilCallDeclStmt>& stmt) override {}

  void visit(const std::shared_ptr<iir::BoundaryConditionDeclStmt>& stmt) override {}

  void visit(const std::shared_ptr<iir::AssignmentExpr>& expr) override;

  void visit(const std::shared_ptr<iir::UnaryOperator>& expr) override {}

  void visit(const std::shared_ptr<iir::BinaryOperator>& expr) override {}

  void visit(const std::shared_ptr<iir::TernaryOperator>& expr) override {}

  void visit(const std::shared_ptr<iir::FunCallExpr>& expr) override {}

  void visit(const std::shared_ptr<iir::StencilFunCallExpr>& expr) override {}

  virtual void visit(const std::shared_ptr<iir::StencilFunArgExpr>& expr) override {}

  void visit(const std::shared_ptr<iir::VarAccessExpr>& expr) override {}

  void visit(const std::shared_ptr<iir::LiteralAccessExpr>& expr) override {}

  void visit(const std::shared_ptr<iir::FieldAccessExpr>& expr) override {}
  void visit(const std::shared_ptr<iir::ReductionOverNeighborExpr>& expr) override {}
};

} // namespace dawn

#endif  // DAWN_OPTIMIZER_INTEGRITYCHECKER_H