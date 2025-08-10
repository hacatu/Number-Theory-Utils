#define _GNU_SOURCE
#include <inttypes.h>
#include <stddef.h>
#include <stdalign.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <nut/cl.h>
#include <nut/debug.h>
#include <nut/factorization.h>

// need to double include so vscode doesn't automagically insert it above feature test macro defs
// and break the code
#include <CL/cl.h>

#define CLKIDX_DIRICHLET_D 0

static size_t nut_cl_get_nck_size(nut_ClMgr *mgr){
	size_t min_size = offsetof(nut_ClKernel, work_item_sizes) + mgr->max_work_item_dimensions*sizeof(size_t);
	size_t align = alignof(nut_ClKernel);
	// round up min_size so it is a multiple of align
	min_size += align - 1;
	return min_size - min_size%align;
}

static void nut_cl_log_error(nut_ClMgr *mgr, const char *title, const char *subtitle, const char *evar, const char *what){
	if(mgr->flags & NUT_CL_FLAG_VERBOSE){
		bool need_lock = mgr->flags & NUT_CL_FLAG_MT_SAFE;
		if(need_lock){
			pthread_mutex_lock(&mgr->log_lock);
		}
		fprintf(stderr, "\e[1;31mERROR: %s (%s=%d): %s%s\e[0m\n", title, evar, mgr->err, what, subtitle);
		if(need_lock){
			pthread_mutex_unlock(&mgr->log_lock);
		}
	}
	if(mgr->flags & NUT_CL_FLAG_DIE){
		exit(EXIT_FAILURE);
	}
}

bool nut_cl_check_err(nut_ClMgr *mgr, const char *title, const char *subtitle){
	// https://bashbaug.github.io/OpenCL-Docs/html/OpenCL_API.html#error_codes
	if(mgr->err != CL_SUCCESS){
		char *what = "UNKNOWN";
		switch(mgr->err){
#define MAKE_CASE(name) case name: what = #name; break;
MAKE_CASE(CL_SUCCESS)
MAKE_CASE(CL_BUILD_PROGRAM_FAILURE)
MAKE_CASE(CL_COMPILE_PROGRAM_FAILURE)
MAKE_CASE(CL_COMPILER_NOT_AVAILABLE)
MAKE_CASE(CL_DEVICE_NOT_FOUND)
MAKE_CASE(CL_DEVICE_NOT_AVAILABLE)
MAKE_CASE(CL_DEVICE_PARTITION_FAILED)
MAKE_CASE(CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST)
MAKE_CASE(CL_IMAGE_FORMAT_MISMATCH)
MAKE_CASE(CL_IMAGE_FORMAT_NOT_SUPPORTED)
MAKE_CASE(CL_INVALID_ARG_INDEX)
MAKE_CASE(CL_INVALID_ARG_SIZE)
MAKE_CASE(CL_INVALID_ARG_VALUE)
MAKE_CASE(CL_INVALID_BINARY)
MAKE_CASE(CL_INVALID_BUFFER_SIZE)
MAKE_CASE(CL_INVALID_BUILD_OPTIONS)
MAKE_CASE(CL_INVALID_COMMAND_QUEUE)
MAKE_CASE(CL_INVALID_COMPILER_OPTIONS)
MAKE_CASE(CL_INVALID_CONTEXT)
MAKE_CASE(CL_INVALID_DEVICE)
MAKE_CASE(CL_INVALID_DEVICE_PARTITION_COUNT)
MAKE_CASE(CL_INVALID_DEVICE_QUEUE)
MAKE_CASE(CL_INVALID_DEVICE_TYPE)
MAKE_CASE(CL_INVALID_EVENT)
MAKE_CASE(CL_INVALID_EVENT_WAIT_LIST)
MAKE_CASE(CL_INVALID_GLOBAL_OFFSET)
MAKE_CASE(CL_INVALID_GLOBAL_WORK_SIZE)
MAKE_CASE(CL_INVALID_HOST_PTR)
MAKE_CASE(CL_INVALID_IMAGE_DESCRIPTOR)
MAKE_CASE(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR)
MAKE_CASE(CL_INVALID_IMAGE_SIZE)
MAKE_CASE(CL_INVALID_KERNEL)
MAKE_CASE(CL_INVALID_KERNEL_ARGS)
MAKE_CASE(CL_INVALID_KERNEL_DEFINITION)
MAKE_CASE(CL_INVALID_KERNEL_NAME)
MAKE_CASE(CL_INVALID_LINKER_OPTIONS)
MAKE_CASE(CL_INVALID_MEM_OBJECT)
MAKE_CASE(CL_INVALID_OPERATION)
MAKE_CASE(CL_INVALID_PIPE_SIZE)
MAKE_CASE(CL_INVALID_PLATFORM)
MAKE_CASE(CL_INVALID_PROGRAM)
MAKE_CASE(CL_INVALID_PROGRAM_EXECUTABLE)
MAKE_CASE(CL_INVALID_PROPERTY)
MAKE_CASE(CL_INVALID_QUEUE_PROPERTIES)
MAKE_CASE(CL_INVALID_SAMPLER)
MAKE_CASE(CL_INVALID_SPEC_ID)
MAKE_CASE(CL_INVALID_VALUE)
MAKE_CASE(CL_INVALID_WORK_DIMENSION)
MAKE_CASE(CL_INVALID_WORK_GROUP_SIZE)
MAKE_CASE(CL_INVALID_WORK_ITEM_SIZE)
MAKE_CASE(CL_KERNEL_ARG_INFO_NOT_AVAILABLE)
MAKE_CASE(CL_LINK_PROGRAM_FAILURE)
MAKE_CASE(CL_LINKER_NOT_AVAILABLE)
MAKE_CASE(CL_MAP_FAILURE)
MAKE_CASE(CL_MEM_COPY_OVERLAP)
MAKE_CASE(CL_MEM_OBJECT_ALLOCATION_FAILURE)
MAKE_CASE(CL_MISALIGNED_SUB_BUFFER_OFFSET)
MAKE_CASE(CL_OUT_OF_HOST_MEMORY)
MAKE_CASE(CL_OUT_OF_RESOURCES)
MAKE_CASE(CL_MAX_SIZE_RESTRICTION_EXCEEDED)
MAKE_CASE(CL_PROFILING_INFO_NOT_AVAILABLE)
#undef MAKE_CASE
		}
		nut_cl_log_error(mgr, title, subtitle, "err", what);
	}
	return mgr->err == CL_SUCCESS;
}

void nut_cl_check_errno(nut_ClMgr *mgr, const char *title, const char *subtitle){
	mgr->err = errno;
	nut_cl_log_error(mgr, title, subtitle, "errno", strerror(errno));
}

static void cl_ctx_err_fn(const char *errinfo, const void *private_info, size_t cb, void *user_data){
	nut_ClMgr *mgr = user_data;
	nut_cl_log_error(mgr, "", "", "err", errinfo);
}

bool nut_cl_setup(nut_ClMgr *mgr, int flags){
	*mgr = (nut_ClMgr){};
	mgr->flags = flags;
	if(mgr->flags & NUT_CL_FLAG_MT_SAFE){
		pthread_mutex_init(&mgr->log_lock, NULL);
	}

	bool verbose = flags & NUT_CL_FLAG_VERBOSE;
	if(verbose){
		fprintf(stderr, ">>> Starting OpenCL...\n");
	}
	
	mgr->err = clGetPlatformIDs(0, NULL, &mgr->num_platforms);
	if(!nut_cl_check_err(mgr, "OpenCL failed to get num platform ids", "")){
		return false;
	}
	mgr->platform_ids = malloc(mgr->num_platforms*sizeof(cl_platform_id));
	if(!mgr->platform_ids){
		mgr->err = CL_OUT_OF_HOST_MEMORY;
		nut_cl_check_err(mgr, "Could not allocate platform_ids buffer", "");
		return false;
	}
	mgr->err = clGetPlatformIDs(mgr->num_platforms, mgr->platform_ids, NULL);
	if(!nut_cl_check_err(mgr, "OpenCL failed to get platform ids list", "")){
		return false;
	}
	char name_buf[1024] = {};
	if(verbose){
		fprintf(stderr, "\e[1;33mFound %u platform%s:\e[0m\n", mgr->num_platforms, mgr->num_platforms == 1 ? "" : "s");
	}
	bool found_gpu = false;
	for(size_t i = 0; i < mgr->num_platforms; ++i){
		mgr->err = clGetPlatformInfo(mgr->platform_ids[i], CL_PLATFORM_PROFILE, 1024, name_buf, NULL);
		if(!nut_cl_check_err(mgr, "OpenCL failed to get next platform profile", "")){
			return false;
		}
		if(verbose){
			fprintf(stderr, "%zu: CL_PLATFORM_PROFILE: %s\n", i, name_buf);
		}
		mgr->err = clGetPlatformInfo(mgr->platform_ids[i], CL_PLATFORM_VERSION, 1024, name_buf, NULL);
		if(!nut_cl_check_err(mgr, "OpenCL failed to get next platform version", "")){
			return false;
		}
		if(verbose){
			fprintf(stderr, "%zu: CL_PLATFORM_VERSION: %s\n", i, name_buf);
		}
		
		if(found_gpu){
			if(verbose){
				fprintf(stderr, "%zu: (skipping GPU search since we already found one)\n", i);
			}
			continue;
		}
		mgr->err = clGetDeviceIDs(mgr->platform_ids[i], CL_DEVICE_TYPE_GPU, 0, NULL, &mgr->num_devices);
		if(mgr->err == CL_DEVICE_NOT_FOUND){
			if(verbose){
				fprintf(stderr, "%zu: (no GPU found for this platform)\n", i);
			}
			continue;
		}
		if(!nut_cl_check_err(mgr, "OpenCL failed to get num device ids", "")){
			return false;
		}
		if(verbose){
			fprintf(stderr, "%zu: Found %u GPUs:\n", i, mgr->num_devices);
		}
		found_gpu = true;
		mgr->device_ids = malloc(mgr->num_devices*sizeof(cl_device_id));
		if(!mgr->device_ids){
			mgr->err = CL_OUT_OF_HOST_MEMORY;
			nut_cl_check_err(mgr, "Could not allocate device_ids buffer", "");
			return false;
		}
		mgr->err = clGetDeviceIDs(mgr->platform_ids[i], CL_DEVICE_TYPE_GPU, mgr->num_devices, mgr->device_ids, NULL);
		if(!nut_cl_check_err(mgr, "OpenCL failed to get device ids", "")){
			return false;
		}
		if(verbose){
			fprintf(stderr, "\e[1;33mFound %u device%s:\e[0m\n", mgr->num_devices, mgr->num_devices == 1 ? "" : "s");
		}
		for(size_t i = 0; i < mgr->num_devices; ++i){
			mgr->err = clGetDeviceInfo(mgr->device_ids[i], CL_DEVICE_NAME, 1024, name_buf, NULL);
			if(!nut_cl_check_err(mgr, "OpenCL failed to get next device name", "")){
				return false;
			}
			if(verbose){
				fprintf(stderr, "%zu: CL_DEVICE_NAME: %s\n", i, name_buf);
			}
		}
		if(verbose){
			fprintf(stderr, "----------\n");
		}

	}
	if(verbose){
		fprintf(stderr, "----------\n");
	}

	mgr->context = clCreateContext(NULL, 1, &mgr->device_ids[0], cl_ctx_err_fn, &mgr, &mgr->err);
	if(!nut_cl_check_err(mgr, "OpenCL failed to create context", "")){
		return false;
	}

	mgr->err = clGetDeviceInfo(mgr->device_ids[0], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &mgr->max_compute_units, NULL);
	if(!nut_cl_check_err(mgr, "OpenCL failed to get max compute units", "")){
		return false;
	}else if(verbose){
		fprintf(stderr, "   MAX_COMPUTE_UNITS: %u\n", mgr->max_compute_units);
	}

	mgr->err = clGetDeviceInfo(mgr->device_ids[0], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &mgr->max_work_group_size, NULL);
	if(!nut_cl_check_err(mgr, "OpenCL failed to get max work group size", "")){
		return false;
	}else if(verbose){
		fprintf(stderr, "   MAX_WORK_GROUP_SIZE: %zu\n", mgr->max_work_group_size);
	}

	mgr->err = clGetDeviceInfo(mgr->device_ids[0], CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(cl_uint), &mgr->max_work_item_dimensions, NULL);
	if(!nut_cl_check_err(mgr, "OpenCL failed to get max work item dimensions", "")){
		return false;
	}else if(verbose){
		fprintf(stderr, "   MAX_WORK_ITEM_DIMENSIONS: %u\n", mgr->max_compute_units);
	}

	mgr->max_work_item_sizes = malloc(mgr->max_work_item_dimensions*sizeof(size_t));
	if(!mgr->max_work_item_sizes){
		mgr->err = CL_OUT_OF_HOST_MEMORY;
		nut_cl_check_err(mgr, "Could not allocate mgr->max_work_item_sizes", "");
		return false;
	}
	mgr->err = clGetDeviceInfo(mgr->device_ids[0], CL_DEVICE_MAX_WORK_ITEM_SIZES, mgr->max_work_item_dimensions*sizeof(size_t), mgr->max_work_item_sizes, NULL);
	if(!nut_cl_check_err(mgr, "OpenCL failed to get max work item sizes", "")){
		return false;
	}else if(verbose){
		fprintf(stderr, "   MAX_WORK_ITEM_SIZES: %zu", mgr->max_work_item_sizes[0]);
		for(size_t i = 1; i < mgr->max_work_item_dimensions; ++i){
			fprintf(stderr, ", %zu", mgr->max_work_item_sizes[i]);
		}
		fprintf(stderr, "\n");
	}

	mgr->queue = clCreateCommandQueueWithProperties(mgr->context, mgr->device_ids[0], NULL, &mgr->err);
	if(!nut_cl_check_err(mgr, "OpenCL failed to create command queue", "")){
		return false;
	}

	mgr->programs = malloc(1*sizeof(cl_program));
	if(!mgr->programs){
		mgr->err = CL_OUT_OF_HOST_MEMORY;
		nut_cl_check_err(mgr, "Could not allocate programs buffer", "");
		return false;
	}
	mgr->kernels = malloc(1*nut_cl_get_nck_size(mgr));
	if(!mgr->kernels){
		mgr->err = CL_OUT_OF_HOST_MEMORY;
		nut_cl_check_err(mgr, "Could not allocate kernels buffer", "");
		return false;
	}

	mgr->programs_cap = 1;
	mgr->kernels_cap = 1;
	
	return true;
}

void nut_cl_close(nut_ClMgr *mgr){
	for(size_t i = 0; i < mgr->kernels_len; ++i){
		clReleaseKernel(mgr->kernels[i].kernel);
	}
	for(size_t i = 0; i < mgr->programs_len; ++i){
		clReleaseProgram(mgr->programs[i]);
	}
	free(mgr->max_work_item_sizes);
	free(mgr->kernels);
	free(mgr->programs);
	clReleaseCommandQueue(mgr->queue);
	clReleaseContext(mgr->context);
	free(mgr->device_ids);
	free(mgr->platform_ids);
	if(mgr->flags & NUT_CL_FLAG_MT_SAFE){
		pthread_mutex_destroy(&mgr->log_lock);
	}
	*mgr = (nut_ClMgr){};
}

const char *nut_cl_read_source(nut_ClMgr *mgr, const char *filename){
	struct stat stat_buf;
	char *kernel_source = NULL;
	const char *base_dir = getenv("NUT_CL_BASEDIR") ?: "../../src/opencl";
	char *full_path [[gnu::cleanup(cleanup_free)]] = NULL;
	if(asprintf(&full_path, "%s/%s", base_dir, filename) < 0){
		nut_cl_check_errno(mgr, "Could not concatenate path", filename);
		goto CLEANUP;
	}
	FILE *file = fopen(full_path, "r");
	if(!file){
		nut_cl_check_errno(mgr, "Could not open cl source file", full_path);
		goto CLEANUP;
	}
	if(fstat(fileno(file), &stat_buf) < 0){
		nut_cl_check_errno(mgr, "Could not stat cl source file", full_path);
		goto CLEANUP;
	}
	kernel_source = malloc(stat_buf.st_size + 1);
	if(!kernel_source){
		nut_cl_check_errno(mgr, "Could not allocate buffer for cl source file", full_path);
		goto CLEANUP;
	}
	(void)!fread(kernel_source, 1, stat_buf.st_size, file);
	if(ferror(file)){
		nut_cl_check_errno(mgr, "Could not read cl source file", full_path);
		goto CLEANUP;
	}
	if(file){
		fclose(file);
	}
	kernel_source[stat_buf.st_size] = '\0';
	return kernel_source;
	CLEANUP:;
	free(kernel_source);
	return NULL;
}

static bool nut_cl_try_ensure_buffer_cap(nut_ClMgr *mgr, void **buf, size_t *cap, size_t min_len, size_t elem_size, const char *err_subtitle){
	if(min_len > *cap){
		void *tmp = NULL;
		size_t t = 64 - __builtin_clzll(min_len - 1);
		size_t new_cap = 1;
		if(t != 64){
			new_cap <<= t;
			if(*cap){
				tmp = realloc(*buf, new_cap*elem_size);
			}else{
				tmp = malloc(new_cap*elem_size);
			}
		}
		if(!tmp){
			mgr->err = CL_OUT_OF_HOST_MEMORY;
			nut_cl_check_err(mgr, "Could not extend buffer", err_subtitle);
			return false;
		}
		*buf = tmp;
		*cap = new_cap;
	}
	return true;
}

bool nut_cl_make_program_from_source(nut_ClMgr *mgr, const char *source){
	if(!nut_cl_try_ensure_buffer_cap(mgr, (void**)&mgr->programs, &mgr->programs_cap, mgr->programs_len + 1, sizeof(cl_program), "programs")){
		return false;
	}
	cl_program program = clCreateProgramWithSource(mgr->context, 1, &source, NULL, &mgr->err);
	if(!nut_cl_check_err(mgr, "OpenCL failed to compile program", "")){
		return false;
	}
	mgr->err = clBuildProgram(program, 0, NULL, "", NULL, NULL);
	int orig_die = mgr->flags & NUT_CL_FLAG_DIE;
	mgr->flags &= ~NUT_CL_FLAG_DIE;
	bool built = nut_cl_check_err(mgr, "OpenCL failed to build program", "\nContinuing in order to show errors...");
	mgr->flags ^= orig_die;
	
	size_t log_size;
	mgr->err = clGetProgramBuildInfo(program, mgr->device_ids[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
	if(!nut_cl_check_err(mgr, "clGetProgramBuildInfo failed getting size", "")){
		return false;
	}
	char *messages = malloc((log_size + 1)*sizeof(char));
	if(!messages){
		mgr->err = CL_OUT_OF_HOST_MEMORY;
		nut_cl_check_err(mgr, "Could not allocate messages buffer", "");
		return false;
	}
	mgr->err = clGetProgramBuildInfo(program, mgr->device_ids[0], CL_PROGRAM_BUILD_LOG, log_size, messages, NULL);
	if(!nut_cl_check_err(mgr, "clGetProgramBuildInfo failed getting log", "")){
		return false;
	}
	if(log_size > 0){
		messages[log_size] = '\0';
		printf(">>> OpenCL Compiler message: %s\n", messages);
	}
	if(!built){
		puts("End of messages.  Build failed.");
		if(orig_die){
			exit(EXIT_FAILURE);
		}
		return false;
	}
	free(messages);
	mgr->programs[mgr->programs_len++] = program;
	return true;
}

static bool nut_cl_create_kernel_single(nut_ClMgr *mgr, size_t program_idx, const char *name){
	cl_kernel kernel = clCreateKernel(mgr->programs[program_idx], name, &mgr->err);
	if(!nut_cl_check_err(mgr, "OpenCL failed to get kernel ", name)){
		return false;
	}
	mgr->kernels[mgr->kernels_len].kernel = kernel;
	mgr->err = clGetKernelWorkGroupInfo(kernel, mgr->device_ids[0], CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t),
		&mgr->kernels[mgr->kernels_len].preferred_work_group_size_multiple, NULL);
	if(!nut_cl_check_err(mgr, "OpenCL failed to get kernel preferred work group size multiple for ", name)){
		return false;
	}
	mgr->err = clGetKernelWorkGroupInfo(kernel, mgr->device_ids[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t),
		&mgr->kernels[mgr->kernels_len].work_group_size, NULL);
	if(!nut_cl_check_err(mgr, "OpenCL failed to get kernel work group size for ", name)){
		return false;
	}
	mgr->err = clGetKernelWorkGroupInfo(kernel, mgr->device_ids[0], CL_KERNEL_WORK_GROUP_SIZE, 3*sizeof(size_t),
		&mgr->kernels[mgr->kernels_len].work_item_sizes, NULL);
	++mgr->kernels_len;
	return true;
}

bool nut_cl_create_kernel(nut_ClMgr *mgr, size_t program_idx, const char *name){
	if(!nut_cl_try_ensure_buffer_cap(mgr, (void**)&mgr->kernels, &mgr->kernels_cap, mgr->kernels_len, nut_cl_get_nck_size(mgr), "kernels")){
		return false;
	}
	return nut_cl_create_kernel_single(mgr, program_idx, name);
}

bool nut_cl_create_all_kernels(nut_ClMgr *mgr, size_t program_idx){
	size_t num_kernels;
	mgr->err = clGetProgramInfo(mgr->programs[program_idx], CL_PROGRAM_NUM_KERNELS, sizeof(size_t), &num_kernels, NULL);
	if(!nut_cl_check_err(mgr, "OpenCL failed to get number of kernels in program", "")){
		return false;
	}
	if(!nut_cl_try_ensure_buffer_cap(mgr, (void**)&mgr->kernels, &mgr->kernels_cap, mgr->kernels_len + num_kernels, nut_cl_get_nck_size(mgr), "kernels")){
		return false;
	}
	size_t names_cap = 64*num_kernels;
	size_t names_len;
	char *names = malloc(names_cap);
	// Try to get the concatenated kernel names using clGetProgramInfo, doubling the allocation size of names as necessary
	while(1){
		if(!names){
			mgr->err = CL_OUT_OF_HOST_MEMORY;
			nut_cl_check_err(mgr, "Could not allocate buffer for kernel names", "");
			return false;
		}
		mgr->err = clGetProgramInfo(mgr->programs[program_idx], CL_PROGRAM_KERNEL_NAMES, names_cap, names, &names_len);
		if(mgr->err == CL_SUCCESS){
			names = realloc(names, names_len + 1) ?: names;
			break;
		}else if(mgr->err != CL_INVALID_VALUE){
			nut_cl_check_err(mgr, "OpenCL failed to get kernel names in program", "");
			return false;
		}
		free(names);
		if(__builtin_mul_overflow(names_cap, 2, &names_cap)){
			names = NULL;
		}else{
			names = malloc(names_cap);
		}
	}
	for(char *str1 = names, *saveptr1;; str1 = NULL){
		const char *name = strtok_r(str1, ";", &saveptr1);
		if(!name){
			break;
		}
		nut_cl_create_kernel_single(mgr, program_idx, name);
	}
	free(names);
	return true;
}

uint64_t nut_cl_dirichlet_D(nut_ClMgr *mgr, uint64_t n, uint64_t m){
	bool verbose = mgr->flags & NUT_CL_FLAG_VERBOSE;

	uint64_t y = nut_u64_nth_root(n, 2);
	cl_mem buf_sums = clCreateBuffer(mgr->context, CL_MEM_READ_WRITE, y*sizeof(uint64_t), NULL, &mgr->err);
	if(!nut_cl_check_err(mgr, "OpenCL could not create work buffer for kernel", "")){
		return 0;
	}
	//clEnqueueWriteBuffer(mgr->queue, buf_sums, CL_TRUE, 0, y*sizeof(uint64_t), sums, 0, NULL, NULL);
	mgr->err = clSetKernelArg(mgr->kernels[CLKIDX_DIRICHLET_D].kernel, 0, sizeof(uint64_t), &n);
	if(!nut_cl_check_err(mgr, "OpenCL failed to set kernel argument", "")){
		return 0;
	}
	mgr->err = clSetKernelArg(mgr->kernels[CLKIDX_DIRICHLET_D].kernel, 1, sizeof(uint64_t), &y);
	if(!nut_cl_check_err(mgr, "OpenCL failed to set kernel argument", "")){
		return 0;
	}
	mgr->err = clSetKernelArg(mgr->kernels[CLKIDX_DIRICHLET_D].kernel, 2, sizeof(uint64_t), &m);
	if(!nut_cl_check_err(mgr, "OpenCL failed to set kernel argument", "")){
		return 0;
	}
	mgr->err = clSetKernelArg(mgr->kernels[CLKIDX_DIRICHLET_D].kernel, 3, sizeof(cl_mem), &buf_sums);
	if(!nut_cl_check_err(mgr, "OpenCL failed to set kernel argument", "")){
		return 0;
	}

	size_t max_wg_size = mgr->kernels[CLKIDX_DIRICHLET_D].work_group_size;
	size_t pref_wg_size_mult = mgr->kernels[CLKIDX_DIRICHLET_D].preferred_work_group_size_multiple;
	if(max_wg_size > pref_wg_size_mult){
		max_wg_size -= max_wg_size%pref_wg_size_mult;
	}
	size_t max_wg_count = mgr->max_compute_units;
	size_t local_dims[1] = {max_wg_size};
	size_t global_dims[1] = {max_wg_size*max_wg_count};
	if(max_wg_size*max_wg_count > y){
		size_t trunc_y = y - y%pref_wg_size_mult;
		if(!trunc_y){
			local_dims[0] = y;
			global_dims[0] = y;
		}else if(trunc_y <= max_wg_size){
			local_dims[0] = trunc_y;
			global_dims[0] = trunc_y;
		}else{
			max_wg_count = trunc_y/max_wg_size;
			global_dims[0] = max_wg_count;
		}
	}

	cl_event event = NULL;
	if(verbose){
		fprintf(stderr, ">>> Ready to launch OpenCL kernel\n");
	}
	mgr->err = clEnqueueNDRangeKernel(mgr->queue, mgr->kernels[CLKIDX_DIRICHLET_D].kernel, 1, NULL, global_dims, local_dims, 0, NULL, &event);
	if(!nut_cl_check_err(mgr, "OpenCL failed to run kernel", "")){
		return 0;
	}
	mgr->err = clWaitForEvents(1, &event);
	if(!nut_cl_check_err(mgr, "OpenCL kernel running failed", "")){
		return 0;
	}
	if(verbose){
		fprintf(stderr, ">>> Finished running OpenCL kernel; reading results\n");
	}
	uint64_t res;
	mgr->err = clEnqueueReadBuffer(mgr->queue, buf_sums, CL_TRUE, 0, sizeof(uint64_t), &res, 0, NULL, NULL);
	if(!nut_cl_check_err(mgr, "OpenCL failed to read result buffer", "")){
		return 0;
	}
	mgr->err = clReleaseMemObject(buf_sums);
	if(!nut_cl_check_err(mgr, "OpenCL failed to release result buffer", "")){
		return 0;
	}
	return res;
}

