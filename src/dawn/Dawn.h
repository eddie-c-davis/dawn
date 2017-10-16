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

#ifndef DAWN_DAWN_H
#define DAWN_DAWN_H

#ifdef DAWN_DOXYGEN_INVOKED
/**

@defgroup codegen CodeGen
@brief Various code generation routines

@defgroup compiler Compiler
@brief Compiler driver API of Dawn

@defgroup optimizer Optimizer
@brief Optimizer infrastructure (including static analysis passes)

@defgroup sir SIR
@brief Implementation of the Stencil Intermediate Representation

@defgroup support Support
@brief Utility functionality used by the Dawn project

@defgroup unittest Unittest
@brief Unittest utility functionality

**/
#endif DAWN_DOXYGEN_INVOKED

#include "dawn/Compiler/Compiler.h"
#include "dawn/SIR/SIR.h"

#endif
