from math import gcd
from contextlib import contextmanager
import sys

class CodegenIsCompositeMark:
	def __init__(self, wheel_size):
		self.names = {
			"function.mark_is_composite.name": "mark_is_composite",
			"function.mark_is_composite.decltype": "static void ",
			"param.is_composite.name": "is_composite",
			"param.is_composite.decltype": "uint8_t *",
			"param.q.name": "q",
			"param.q.decltype": "uint64_t ",
			"param.r.name": "r",
			"param.r.decltype": "uint64_t ",
			"param.q_ub.name": "q_ub",
			"param.q_ub.decltype": "uint64_t ",
			"var.mq.name": "mq",
			"var.mq.decltype": "uint64_t ",
		}
		self.wheel_size = wheel_size
		self.wheel_spokes = [d for d in range(wheel_size) if gcd(wheel_size, d) == 1]
		self.lines = []
		self.indent = 0
	
	def addLine(self, text):
		self.lines.append("\t"*self.indent + text)
	
	def addMultilineComment(self, *lines):
		for i, line in enumerate(lines):
			prefix = "/* " if i == 0 else " * "
			self.addLine(prefix + line)
		self.addLine(" */")
	
	def addFunctionDef(self, fn_name, *param_names):
		self.addLine(self.names[f"function.{fn_name}.decltype"] + self.names[f"function.{fn_name}.name"] + "(" + ", ".join(self.names[f"param.{n}.decltype"] + self.names[f"param.{n}.name"] for n in param_names) + "){")
	
	def addVarDecl(self, var_name, ini = None):
		self.addLine(self.names[f"var.{var_name}.decltype"] + self.names[f"var.{var_name}.name"] + (";" if ini is None else ini + ";"))
	
	@contextmanager
	def makeIndentedBlock(self, end="}"):
		try:
			self.indent += 1
			yield
		finally:
			self.indent -= 1
			if end is not None:
				self.addLine(end)
	
	def codegen(self):
		self.addFunctionDef("mark_is_composite", "is_composite", "q", "r", "q_ub")
		with self.makeIndentedBlock():
			self.addMultilineComment(
				"We need to find the first multiple of the prime 30*q + r which could possibly be prime and is at least the prime squared.",
				"Since the prime squared is still coprime to the wheel size, we start at the prime squared.",
				"(30*q + r)**2 = 30*30*q**2 + 30*2*q*r + r**2 = 30*(30*q**2 + 2*q*r) + r**2")
			self.addVarDecl("mq");
			self.addLine("switch(" + self.names["param.r.name"] + "){")
			with self.makeIndentedBlock():
				for wheel_spoke in self.wheel_spokes:
					self.codegenWheelCase(wheel_spoke)
				self.addLine("default:")
				with self.makeIndentedBlock(None):
					self.addLine("__builtin_unreachable();")
	
	def codegenWheelCase(self, r):		
		self.addLine(f"case {r}:")
		with self.makeIndentedBlock(None):
			w = self.wheel_size
			i = self.wheel_spokes.index(r)
			# (w*q + r)*(w*q + r) = w*w*q*q + w*q*2*r + r*r = w*(q*(w*q + 2*r) + r*r//w) + r*r%w
			mq = self.names["var.mq.name"]
			q = self.names["param.q.name"]
			q_ub = self.names["param.q_ub.name"]
			self.addLine(f"{mq} = {q}*({w}*{q} + {2*r}) + {r*r//w};// mr = {r*r%w}")
			l = len(self.wheel_spokes)
			mq_mults = []
			for j in range(l):
				s = self.wheel_spokes[(i + j)%l]
				if s < r:
					s += w
				q_coeff = (s - r)%w
				q_shift = r*s//w - r*r//w;
				mr = r*s%w
				mq_mults.append((q_coeff, q_shift, mr))
			self.addLine(f"for(; {mq} + {mq_mults[-1][0]}*{q} + {mq_mults[-1][1]} < {q_ub}; {mq} += {w}*{q} + {r}){{")
			is_composite = self.names["param.is_composite.name"]
			with self.makeIndentedBlock():
				mr = mq_mults[0][2]
				mask = hex(1 << self.wheel_spokes.index(mr))
				self.addLine(f"{is_composite}[{mq}] |= {mask};// mr = {mr}")
				for q_coeff, q_shift, mr in mq_mults[1:]:
					mask = hex(1 << self.wheel_spokes.index(mr))
					self.addLine(f"{is_composite}[{mq} + {q_coeff}*{q} + {q_shift}] |= {mask};// mr = {mr}")
			mr = mq_mults[0][2]
			mask = hex(1 << self.wheel_spokes.index(mr))
			self.addLine(f"if({mq} >= {q_ub}){{break;}}")
			self.addLine(f"{is_composite}[{mq}] |= {mask};// mr = {mr}")
			for i, (q_coeff, q_shift, mr) in enumerate(mq_mults[1:-1]):
				mask = hex(1 << self.wheel_spokes.index(mr))
				self.addLine(f"if(({mq} += {q_coeff - mq_mults[i][0]}*{q} + {q_shift - mq_mults[i][1]}) >= {q_ub}){{break;}}")
				self.addLine(f"{is_composite}[{mq}] |= {mask};// mr = {mr}")
			self.addLine(f"{mq} += {mq_mults[-1][0] - mq_mults[-2][0]}*{q} + {mq_mults[-1][1] - mq_mults[-2][1]}; break; // If we get here, we have {mq} >= {q_ub} by the fact that we made it out of the for loop")

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print("Please specify snippet to codegen!")
		sys.exit(1)
	match sys.argv[1]:
		case "mark_is_composite":
			if len(sys.argv) != 3:
				print("Please specify wheel size!")
				sys.exit(1)
			g = CodegenIsCompositeMark(int(sys.argv[2]))
			g.codegen()
			for line in g.lines:
				print(line)
		case _:
			print("Unrecognized snippet name!")
			sys.exit(1)

