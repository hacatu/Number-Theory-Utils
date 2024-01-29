#pragma once

/// @file
/// @author hacatu
/// @version 0.2.0
/// @section LICENSE
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
/// @section DESCRIPTION
/// OpenCL implementations of functions in dirichlet.h and dirichlet_powerful.h
/// Generally, these functions require an OpenCL manager ({ @link nut_ClMgr })
/// and if mgr->flags does not have die or verbose set, these functions will
/// usually return 0 on error so if 0 is returned and neither die nor verbose is set,
/// check mgr->err before relying on the result (eg by calling { @link nut_cl_check_err }).

#include <nut/cl.h>

uint64_t nut_cl_dirichlet_D(nut_ClMgr *mgr, uint64_t n, uint64_t m);

