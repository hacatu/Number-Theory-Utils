#include <inttypes.h>
#include <stdio.h>

#include <nut/cl.h>
#include <nut/dirichlet_cl.h>

int main(){
	nut_ClMgr mgr;
	nut_cl_setup(&mgr, NUT_CL_FLAG_DIE | NUT_CL_FLAG_VERBOSE);

	const char *kernel_source = nut_cl_read_source(&mgr, "../../src/opencl/dirichlet.cl");
	nut_cl_make_program_from_source(&mgr, kernel_source);
	free((char*)kernel_source);
	nut_cl_create_all_kernels(&mgr, 0);
	
	uint64_t n = 1000, m = 0;
	uint64_t Dn = nut_cl_dirichlet_D(&mgr, n, m);

	fprintf(stderr, "cl_dircichlet_D(%"PRIu64", %"PRIu64") = %"PRIu64"\n", n, m, Dn);
	
	nut_cl_close(&mgr);
}

