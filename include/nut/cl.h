#pragma once

/// @file
/// @author hacatu
/// @version 0.2.0
/// @section LICENSE
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this
/// file, You can obtain one at http://mozilla.org/MPL/2.0/.
/// @section DESCRIPTION
/// Common OpenCL wrapper functions, used by all programs that want to use
/// OpenCL implementations of any algorithm supported by Number Theory Utils

#define CL_TARGET_OPENCL_VERSION 220

#include <CL/cl.h>
#include <pthread.h>

/// Set on { @link nut_ClMgr } to cause { @link nut_cl_check_err } and other
/// functions to terminate the program if an error is encountered
#define NUT_CL_FLAG_DIE 0x1
/// Set on { @link nut_ClMgr } to cause { @link nut_cl_check_err } and other
/// functions to display more information on stderr
#define NUT_CL_FLAG_VERBOSE 0x2
/// Set on { @link nut_ClMgr } to cause functions, especially the
/// error reporting callback passed to clCreateContext, to access stdio under
/// a lock
#define NUT_CL_FLAG_MT_SAFE 0x4

typedef struct{
	cl_kernel kernel;
	size_t preferred_work_group_size_multiple;
	size_t work_group_size;
	size_t work_item_sizes[];
} nut_ClKernel;

/// Manages "global" OpenCL state including programs, kernels, the context, and the command queue.
/// Multiple can be created, but this would require manually setting up all but the first to use
/// distinct devices, otherwise there would be multiple for the same device.
typedef struct{
	int flags;
	cl_int err;
	cl_uint num_platforms;
	cl_platform_id *platform_ids;
	cl_uint num_devices;
	cl_device_id *device_ids;
	cl_uint max_compute_units;
	cl_uint max_work_item_dimensions;
	size_t max_work_group_size;
	size_t *max_work_item_sizes;
	cl_context context;
	cl_command_queue queue;
	size_t programs_len, programs_cap;
	cl_program *programs;
	size_t kernels_len, kernels_cap;
	nut_ClKernel *kernels;
	pthread_mutex_t log_lock;
} nut_ClMgr;

/// Check and report/resolve OpenCL errors after manually calling OpenCL functions
/// Will print errors and/or terminate the program if there are errors, based on the flags
/// set in mgr.  Wrapper functions already do this, so only call this function after
/// manually calling OpenCL functions.  In this case, ensure that mgr->err gets set to the
/// OpenCL error from the function, ie by doing `mgr->err = clErrRetFn(...)` or `clErrOutArgFn(&mgr->err)`
/// as appropriate.
/// @param[in] mgr: mgr->flags influences whether or not this function will print,
/// lock stdio, and/or kill the program on errors
/// @param[in] title: printed before the description of the error if printing occurs
/// @param[in] subtitle: printed after the description of the error if printing occurs
/// @return true if mgr->err == CL_SUCCESS, false otherwise (unless mgr->flags has
/// NUT_CL_FLAG_DIE set, in which case the program terminates so this function never returns)
bool nut_cl_check_err(nut_ClMgr *mgr, const char *title, const char *subtitle);

/// Initialize OpenCL and store the relevant state in mgr
/// The first available platform and device will be chosen, this might not be
/// right for systems with multiple gpus
/// The settings in flags will be passed on to mgr and govern
/// how errors will be handled during this function call, but mgr->flags can be
/// changed afterwards (eg you can set/unset verbose/die to help clear up clutter and avoid handling
/// unresolvable errors).  However, mt_safe should not be changed, as this changes the
/// state expected of the log lock.
/// @return true on success, false on failure (unless flags has die set, in which case the program terminates)
bool nut_cl_setup(nut_ClMgr *mgr, int flags);

/// Shut down OpenCL and clean up all resources associated with mgr
/// This cleans up all programs, kernels, and other opencl objects associated with mgr
/// and then shuts down the context and command queue.
/// However, OpenCL still has a lot of "still reachable" memory leaks and a couple other
/// spurious issues under valgrind/msan, so not only is it not necessary to call this as
/// long as you don't leak mgr or need the memory back before the program ends,
/// opencl will actually leave some memory still allocated even after this function is called.
void nut_cl_close(nut_ClMgr *mgr);

/// Read an OpenCL source file (.cl) into a string
/// The filename can be any relative path, and is appended to
/// the base path, which defaults to "../../src/opencl" but can be
/// overridden by setting the environment variable "NUT_CL_BASEDIR"
/// The returned string is malloc'd and should be free'd after it is
/// used (once it has been passed to { @link nut_cl_make_program_from_source }).
const char *nut_cl_read_source(nut_ClMgr *mgr, const char *filename);

/// Compile OpenCL source read into an OpenCL program
/// The source should come from { @link nut_cl_read_source }.
/// To determine the program, note mgr->programs_len before this call;
/// if p is its value before calling this function, on success mgr->programs[p]
/// will hold the compiled program.
/// mgr->programs is subject to reallocation, so do not store pointers to its elements;
/// rather store copies of them (`cl_program` is a typedef of a pointer so this is ok)
/// or their indicies.
/// Besides { @link nut_cl_read source }, the source code could also be embedded
/// as a string within the program code or stored any other way you see fit.
/// Unfortunately, loading programs from spirv IL is not currently supported,
/// it seems there are severe limitations to compiling OpenCL programs to this form
/// as any kernels that use 64 bit integer arithmetic will fail to translate from
/// llvm bytecode to spriv IL
bool nut_cl_make_program_from_source(nut_ClMgr *mgr, const char *source);

/// Create one kernel with a given name
/// program_idx should be the index of the program in mgr.  This will be
/// 0 for the first program created with { @link nut_cl_make_program_from_source },
/// 1 for the second, and so on.
/// All kernels from all programs are currently stored in one array, mgr->kernels.
/// To find where the kernel will show up, take note of mgr->num_kernels before calling
/// this function.
/// mgr->kernels is subject to reallocation, so do not store pointers to its elements;
/// rather store copies of them (`cl_kernel` is a typedef of a pointer so this is ok)
/// or their indicies.
bool nut_cl_create_kernel(nut_ClMgr *mgr, size_t program_idx, const char *name);

/// Create all kernels within a given program
/// Automatically extracts kernel names from the program and gets references to all of them.
/// If mgr->num_kernels equals p before calling this, on success it increases by n the number
/// of kernels, and the kernels are stored at `mgr->kernels[p]` through `mgr->kernels[p + n - 1]`
/// inclusive.  The kernels should be created in the order they appear in the source,
/// but this might not be guaranteed by `clGetProgramInfo`.
bool nut_cl_create_all_kernels(nut_ClMgr *mgr, size_t program_idx);

