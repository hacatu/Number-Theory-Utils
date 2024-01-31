/*
Compute the divisor summatory function D(n) using the dirichlet hyperbola method
This is the sum of the divisor count for every integer 1 <= x <= n
The basic formula is D(n) = 2*sum(n/x for 1 <= x <= y) - y*y, where y = floor(sqrt(n))
This function breaks the sum for 1 <= x <= y into w sums, where w is the number of global work items.
Each work item sums n/x for all x equal to its global id modulo w, and stores this result in a global array.
Then, these partial sums are combined using the classic parallel sum algorithm: for every 2^t
 up to w, each work item whose id is an exact multiple of 2^(t+1) adds the partial sum offset + 2^t to its own
 if that partial sum is not out of range.  This computes the overall sum in only log(w) time

 TODO: having m be a variable is insanely expensive and most programs will not need many values of m
 We should have it be a constant, but this would require fully separating the m == 0 and m != 0 cases
 into different functions and more problematically having some kind of templating/preprocessing step for opencl code
*/
__kernel void dirichlet_D(const unsigned long n, const unsigned long y, const unsigned long m, __global unsigned long *sums){
	const unsigned long x = get_global_id(0);
	const unsigned long w = get_global_size(0);
	unsigned long partial_sum = n/(x + 1);
	if(m){
		for(unsigned long xx = x + 1; !__builtin_add_overflow(xx, w, &xx) && xx <= y;){
			partial_sum = (partial_sum + n/xx)%m;
		}
	}else{
		for(unsigned long xx = x + 1; !__builtin_add_overflow(xx, w, &xx) && xx <= y;){
			partial_sum += n/xx;
		}
	}
	sums[x] = partial_sum;
	barrier(CLK_GLOBAL_MEM_FENCE);
	if(m){
		for(unsigned long t = 1; t < w; t <<= 1){
			if(!(x & ((t << 1) - 1)) && x + t < w){
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
		for(unsigned long t = 1; t < w; t <<= 1){
			if(!(x & ((t << 1) - 1)) && x + t < w){
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

