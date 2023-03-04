import sys
import itertools as itt

def print_whole_graph(n):
	print("digraph {")
	for i in range(n):
		print(f"\t{i} -> {(i*i + 1)%n}")
	print("}")

def get_cycles(seeds, p):
	frontier = set(seeds)
	cycles = {}
	while frontier:
		x = frontier.pop()
		visited = {x}
		while (x := (x*x + 1)%p) not in visited:
			visited.add(x)
			frontier.discard(x)
		it = x
		def advance_it():
			nonlocal it
			return (it := (it*it + 1)%p)
		cyc = {x} | set(iter(advance_it, x))
		cycles[min(cyc)] = cyc
	return cycles.values()

def get_p_cycles(p):
	return get_cycles(range((p+1)//2), p)

def get_p2_cycles(p):
	p_cyclants = itt.chain.from_iterable(get_p_cycles(p))
	return get_cycles((a + b for (a, b) in itt.product(range(0, p*p, p), p_cyclants)), p*p)

if __name__ == "__main__":
	if len(sys.argv) == 2:
		print_whole_graph(int(sys.argv[1]))

