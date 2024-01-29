__kernel void dirichlet_D(const unsigned long n, const unsigned long y, const unsigned long m, __global unsigned long *sums){
	unsigned long x = get_global_id(0);
	sums[x] = n/(x + 1);
	barrier(CLK_GLOBAL_MEM_FENCE);
	if(m){
		for(unsigned long t = 1; t < y; t <<= 1){
			if(!(x & ((t << 1) - 1)) && x + t < y){
				sums[x] = (sums[x] + sums[x + t])%m;
			}
			barrier(CLK_GLOBAL_MEM_FENCE);
			if(t == (1ull << 63)){
				break;
			}
		}
		if(!x){
			sums[x] = (m + 2*sums[x] - y*y%m)%m;
		}
	}else{
		for(unsigned long t = 1; t < y; t <<= 1){
			if(!(x & ((t << 1) - 1)) && x + t < y){
				sums[x] += sums[x + t];
			}
			barrier(CLK_GLOBAL_MEM_FENCE);
			if(t == (1ull << 63)){
				break;
			}
		}
		if(!x){
			sums[x] = 2*sums[x] - y*y;
		}
	}
}

