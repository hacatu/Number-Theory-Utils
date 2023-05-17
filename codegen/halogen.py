from dataclasses import dataclass
from enum import Enum, auto
from typing import Callable
from numbers import Integral
import string
import itertools as itt
from collections import deque
from datetime import datetime
import re
import sys
from pathlib import Path
import io

def iter_len(it):
	c = itt.count()
	deque(zip(it, c), 0)
	return next(c)

class PPTokenKind(Enum):
	header_name = auto()
	identifier = auto()
	pp_number = auto()
	character_constant = auto()
	string_literal = auto()
	punctuator = auto()
	universal_character_name = auto()
	non_whitespace_character = auto()
	placemarker = auto()

@dataclass
class PPLine:
	text: str
	physical_lines: int
	logical_line: int
	logical_file: str

@dataclass
class PPToken:
	spelling: str
	kind: PPTokenKind
	line: PPLine
	leading_ws: str

class ParseError(Exception):
	pass

class CachedTokenStream:
	def __init__(self, gen):
		self._gen = gen
		self._cache = []
	
	def _extend(self, i):
		if i < len(self._cache):
			return True
		self._cache.extend(itt.islice(self._gen, i - len(self._cache) + 1))
		return i < len(self._cache)
	
	def __bool__(self):
		return self._cache or self._extend(0)
	
	def __len__(self):
		return iter_len(self)
	
	def _sliceFromStart(self, i):
		a, b, s = i.start, i.stop, i.step
		if s is None:
			s = 1
		elif not s:
			raise IndexError("Slice must have a nonzero step size")
		if s < 0:
			if a is None:
				a = 0
			elif a < 0:
				a = max(len(self._cache) + a, 0)
			if b is not None and b < 0:
				b = max(len(self._cache) + b, 0)
			return itt.islice(self, a, b, s)
		if a is None:
			a = len(self._cache) - 1
		elif a < 0:
			a = len(self._cache) + a
			if a < 0:
				return ()
		if b is not None and b < 0:
			b = max(len(self._cache) + b, -1)
		self._extend(a)
		return self._cache[a, b, s]
	
	def _sliceFromEnd(self, i):
		a, b, s = i.start, i.stop, i.step
		if s is None:
			s = 1
		elif not s:
			raise IndexError("Slice must have a nonzero step size")
		if s < 0:
			if a is None:
				a = len(self._cache)
			elif a < 0:
				a = max(len(self._cache) + a, 0)
			else:
				a += len(self._cache)
			if b is not None:
				b = max(len(self._cache) + b, -1)
			return itt.islice(self, a, b, s)
		if a is None:
			a = len(self._cache) - 1
		else:
			a += len(self._cache)
			if a < 0:
				return ()
		if b is not None:
			b += len(self._cache)
		self._extend(a)
		return self._cache[a, b, s]
	
	def __getitem__(self, i):
		if isinstance(i, Integral):
			if i < 0:
				return self._cache[i]
			self._extend(i)
			return self._cache[i]
		elif isinstance(i, slice):
			return self._sliceFromStart(i))
		raise TypeError("Index must be a slice or Integral number")
	
	def __call__(self, i):
		if isinstance(i, Integral):
			return self[len(self._cache) + i]
		elif isinstance(i, slice):
			return self._sliceFromEnd(i))
		raise TypeError("Index must be a slice or Integral number")
	
	def __iter__(self, j=0):
		while True:
			while j < len(self._cache):
				yield self._cache[j]
				j += 1
			tok = next(self._gen)
			self._cache.append(tok)
			yield tok
			j += 1
	
	def __reversed__(self):
		return reversed(self._cache)

class PPParser:
	def __init__(self, ini):
		self.ini = ini
	
	def combineLines(self):
		self.lines = []
		line_builder = []
		if isinstance(self.ini, str):
			lines = self.ini.split("\n")
			self.logical_file = "<text>"
		elif isinstance(self.ini, io.TextIOBase):
			lines = self.ini.readlines()
			self.logical_file = self.ini.name
		else:
			raise TypeError("Unsupported ini argument type")
		logical_line = 1
		for line in lines:
			if line.endswith("\\"):
				line_builder.append(line[:-2])
			else:
				line_builder.append(line)
				physical_lines = len(line_builder)
				self.lines.append(PPLine("".join(line_builder), physical_lines, logical_line, self.logical_file))
				logical_line += physical_lines
				line_builder = []
		if line_builder:
			line_builder[-1] += "\\"
			self.lines.append(PPLine("".join(line_builder), len(line_builder), logical_line, self.logical_file))
		self.col = 0
		self.row = 0
		self.presumed_line = 1
		self.curr_directive_prefix = "" # used to test if we are in a #include, #embed, #if, #elif, #ifdef, #elifdef, #ifndef, #elifndef
		self.in_has_include = False # used to test if we are in a __has_include or __has_embed in an #if or #elif
	
	def getCurLine(self):
		if self.row >= len(self.lines):
			return None
		return self.lines[self.row]
	
	def getCurLineText(self):
		if self.row >= len(self.lines):
			return ""
		return self.lines[self.row].text
	
	def getRestOfLine(self):
		if self.row >= len(self.lines):
			return ""
		line = self.lines[self.row].text
		if self.col >= len(line):
			return ""
		return line[self.col:]
	
	def _scanSingleComment(self):
		text = self.getRestOfLine()
		if text.startswith("//"):
			self.row += 1
			self.col = 0
			self.curr_directive_prefix = " "
			return text
		return ""
	
	def _scanMultiComment(self):
		text = self.getRestOfLine()
		if text.startswith("/*"):
			res_builder = []
			while (end_idx := text.find("*/")) == -1:
				self.row += 1
				if self.row == len(self.lines):
					logical_line = self.lines[-1-len(res_builder)].logical_line
					if res_builder:
						col = len(self.lines[-1-len(res_builder)].text) - len(res_builder[0]) + 1
					else:
						col = self.col + 1
					raise ParseError(f"Unterminated multiline comment on line {logical_line} col {col}")
				self.col = 0
				self.curr_directive_prefix = " "
				res_builder.append(text)
				text = self.getCurLineText()
			res_builder.append(text[:end_idx + 2])
			if end_idx + 2 == len(text):
				self.row += 1
				self.col = 0
				self.curr_directive_prefix = " "
			else:
				self.col += end_idx + 2
			return "\n".join(res_builder)
		return ""
	
	def _scanWhitespace(self):
		res_builder = []
		while self.row != len(self.lines):
			text = self.getRestOfLine()
			if not text:
				res.append("\n")
				self.col = 0
				self.row += 1
				continue
			idx = 0
			while text[idx].isspace():
				idx += 1
				self.col += 1
				if idx == len(text):
					self.row += 1
					self.col = 0
					self.curr_directive_prefix = " "
					break
			if idx:
				res_builder.append(text[:idx])
			if space := self._scanMultiComment():
				res_builder.append(space)
			elif space := self._scanSingleComment():
				res_builder.append(space)
			else:
				break
		return "".join(res_builder)
	
	def _guessTokenKind(self, text):
		"""
		pp_token:
			header_name | identifier | pp_number | character_constant | string_literal | punctuator | universal_character_name | non_whitespace_character
		"""
		if self.in_has_include:
			if text.startswith("<"):
				return PPTokenKind.header_name
			elif text.startswith("\""):
				return PPTokenKind.header_name
		if text.startswith("u8"):
			if len(text) >= 3:
				if text[2] == "'":
					return PPTokenKind.character_constant
				elif text[2] == "\"":
					return PPTokenKind.string_literal
			return PPTokenKind.identifier
		elif text[0] in "uUL":
			if len(text) != 1:
				if text[1] == "'":
					return PPTokenKind.character_constant
				elif text[1] == "\"":
					return PPTokenKind.string_literal
			return PPTokenKind.identifier
		elif text[0] in string.ascii_letters or text[0] == "_":
			return PPTokenKind.identifier
		elif text[0] in string.digits or (text[0] == "." and len(text) != 1 and text[1] in string.digits):
			return PPTokenKind.pp_number
		elif text[0] == "'":
			return PPTokenKind.character_constant
		elif text[0] == "\"":
			return PPTokenKind.string_literal
		elif text[0] in "[](){}.-+&*+!/%<>=^|&?:;,#":
			return PPTokenKind.punctuator
		elif text[0] == "\\" and len(text) != 1 and text[1] in "uU":
			return PPTokenKind.universal_character_name
		else:
			return PPTokenKind.non_whitespace_character
	
	def _scanHeaderName(self, text, space, end):
		"""
		header_name:
			"<" h_char+ ">"
			"\"" q_char+ "\""
		h_char:
			[^>\n]
		q_char:
			[^"\n]
		"""
		for end_idx, c in enumerate(text[1:], 1):
			if c == end:
				break
			elif c in "'\\\"" or (end_idx + 1 != len(text) and c == "/" and text[end_idx + 1] in "/*"):
				raise ParseError("Illegal character sequence in header name")
		else:
			raise ParseError("Illegal new line in header name")
		res = PPToken(text[:end_idx + 1], PPTokenKind.header_name, self.getCurLine(), space)
		if end_idx + 1 == len(text):
			self.row += 1
			self.col = 0
		else:
			self.col += end_idx + 1
		return res
	
	def _scanIdentifier(self, text, space):
		"""
		identifier:
			nondigit (nondigit | digit)*
		nondigit:
			[_a-zA-Z]
		digit:
			[0-9]
		"""
		for end_idx, c in enumerate(text[1:], 1):
			if c not in string.ascii_letters and c not in string.digits and c != "_":
				break
		else:
			res = PPToken(text, PPTokenKind.identifier, self.getCurLine(), space)
			self.row += 1
			self.col = 0
			return res
		res = PPToken(text[:end_idx], PPTokenKind.identifier, self.getCurLine(), space)
		self.col += end_idx
		return res
	
	def _scanPpNumber(self, text, space):
		"""
		pp_number:
			"."? [0-9] (("'"? (digit | nondigit)) | ([eEpP] [+-]) | ".")*
		"""
		idx = 0
		if text[idx] == ".":
			idx += 2
		else:
			idx += 1
		while idx + 1 < len(text):
			if text[idx] == "'":
				if text[idx + 1] in string.ascii_letters or text[idx + 1] in string.digits or text[idx + 1] == "_":
					idx += 2
				else:
					break
			elif text[idx] in "eEpP":
				if text[idx + 1] in "+-":
					idx += 2
				else:
					idx += 1
			elif text[idx] in string.ascii_letters or text[idx] in string.digits or text[idx] == "_":
				idx += 1
			else:
				break
		if idx < len(text) and (text[idx] in string.ascii_letters or text[idx] in string.digits or text[idx] == "_"):
			idx += 1
		res = PPToken(text[:idx], PPTokenKind.pp_number, self.getCurLine(), space)
		if idx == len(text):
			self.row += 1
			self.col = 0
		else:
			self.col += idx
		return res
	
	RE_ESCAPE_SEQ = re.compile(r"""^\['"?\abfnrtv]|[0-7]{1,3}|(?:x[0-9a-fA-F]+)|(?:u[0-9a-fA-F]{4})|(?:U[0-9a-fA-F]{8})""")
	
	def _scanStringLiteral(self, text, space, end):
		r"""
		character_constant:
			([uUL] | "u8")? "'" c_char+ "'"
		c_char:
			[^'\\\n] | escape_sequence
		escape_sequence:
			simple_escape_sequence | octal_escape_sequence | hex_escape_sequence | universal_escape_sequence
		simple_escape_sequence:
			"\\'" | "\\\"" | "\\?" | "\\\\" | "\\a" | "\\b" | "\\f" | "\\n" | "\\r" | "\\t" | "\\v"
		octal_escape_sequence:
			"\\" [0-7]{1,3}
		hex_escape_sequence:
			"\\x" [0-9a-fA-F]+
		universal_character_name:
			"\u" hex_digit{4}
			"\U" hex_digit{8}
		string_literal:
			([uUL] | "u8") "\"" s_char* "\""
		s_char:
			[^"\\\n] | escape_sequence
		"""
		idx = text.find(end) + 1
		while True:
			if idx == len(text):
				raise ParseError("Unexpected new line in char/string literal")
			elif text[idx] == end:
				idx += 1
				break
			elif text[idx] == "\\":
				if (m := RE_ESCAPE_SEQ.match(text[idx])) is not None:
					idx += m.end
				else:
					raise ParseError("Invalid escape sequence")
		res = PPToken(text[:idx], PPTokenKind.string_literal, self.getCurLine(), space)
		if idx == len(text):
			self.row += 1
			self.col = 0
		else:
			self.col += idx
		return res
	
	def _scanPunctuator(self, text, space):
		"""
		punctuator:
			"[ ] ( ) { } . ->
			++ -- & * + - ~ !
			/ % << >> < > <= >= == != ^ | && ||
			? : :: ; ...
			= *= /= %= += -= <<= >>= &= ^= |=
			, # ## <: :> <% %> %: %:%:"
		"""
		if text.startswith("%:%:"):
			idx = 4
		elif any(map(text.startswith, "... <<= >>=".split(" "))):
			idx = 3
		elif any(map(text.startswith, "-> ++ -- << >> <= >= ++ != && || :: *= /= %= += -= &= ## <: :> <% %> %:".split(" "))):
			idx = 2
		else:
			idx = 1
		spelling = text[:idx]
		res = PPToken(spelling, PPTokenKind.punctuator, self.getCurLine(), space)
		if idx == len(text):
			self.row += 1
			self.col = 0
		else:
			self.col += idx
		return res
	
	RE_UNIVERSAL_ESCAPE_SEQ = re.compile(r"""^\\(?:u[0-9a-fA-F]{4})|(?:U[0-9a-fA-F]{8})""")
	
	def _scanUniversalCharacterName(self, text, space):
		if (m := RE_UNIVERSAL_ESCAPE_SEQ.match(text)) is None:
			raise ParseError("Invalid universal character name")
		idx = m.end
		res = PPToken(text[:idx], PPTokenKind.universal_character_name, self.getCurLine(), space)
		if idx == len(text):
			self.row += 1
			self.col = 0
		else:
			self.col += idx
		return res
	
	def _scanNonWhitespaceCharacter(self, text, space):
		res = PPToken(text[:1], PPTokenKind.non_whitespace_character, self.getCurLine(), space)
		if len(text) == 1:
			self.row += 1
			self.col = 0
		else:
			self.col += 1
		return res
	
	def _updateDirective(self):
		if self.curr_directive_prefix in (" ", " #", " #include"):
			self.curr_directive_prefix = ""
			self.in_has_include = False
		elif self.curr_directive_prefix.endswith("__has_include"):
			self.curr_directive_prefix = " #if"
			self.in_has_include = False
	
	def getToken(self):
		space = self._scanWhitespace()
		if self.row == len(self.lines):
			return PPToken("", PPTokenKind.placemarker, self.getCurLine(), space)
		text = self.getRestOfLine()
		k = self._guessTokenKind(text)
		if k is PPTokenKind.header_name:
			self.in_has_include = False
			end = "\"" if text.startswith("\"") else ">"
			return self._scanHeaderName(text, space, end)
		elif k is PPTokenKind.identifier:
			res = self._scanIdentifier(text, space)
			if self.curr_directive_prefix.startswith(" #if "):
				if self.curr_directive_prefix == " #if " and res.spelling in ("__has_include", "__has_embed"):
					self.curr_direcive_prefix += "__has_include"
			elif self.curr_directive_prefix == " #":
				if res.spelling in ("if", "elif"):
					self.curr_directive_prefix += "if "
				elif res.spelling in ("include", "embed"):
					self.curr_directive_prefix += "include"
					self.in_has_include = True
				else:
					self.curr_directive_prefix = ""
					self.in_has_include = False
			else:
				self.curr_directive_prefix = ""
				self.in_has_include = False
			return res
		elif k is PPTokenKind.pp_number:
			self._updateDirective()
			return self._scanPpNumber(text, space)
		elif k is PPTokenKind.character_constant:
			self._updateDirective()
			return self._scanStringLiteral(text, space, "'")
		elif k is PPTokenKind.string_literal:
			self._updateDirective()
			return self._scanStringLiteral(text, space, "\"")
		elif k is PPTokenKind.punctuator:
			res = self._scanPunctuator(text, space)
			if self.curr_directive_prefix.endswith("__has_include") or self.curr_directive_prefix.endswith("__has_embed"):
				if res.spelling == "(":
					self.curr_directive_prefix = " #if ("
					self.in_has_include = True
			elif self.curr_directive_prefix.startswith(" #if "):
				if res.spelling == "(":
					self.curr_directive_prefix += "("
				elif res.spelling == ")" and self.curr_directive_prefix.endswith(")"):
					self.curr_directive_prefix = self.curr_directive_prefix[:-1]
			return res
		elif k is PPTokenKind.universal_character_name:
			self._updateDirective()
			return self._scanUniversalCharacterName(text, space)
		else:
			self._updateDirective()
			return self._scanNonWhitespaceCharacter(text, space)
	
	def genTokens(self):
		while True:
			tok = self.getToken()
			yield tok
			if tok.kind == PPTokenKind.placemarker:
				return

class PPMacroKind(Enum):
	builtin_object = auto()
	builtin_function = auto()
	object_like = auto()
	function_like = auto()

@dataclass
class PPMacroDef:
	spelling: str
	kind: PPMacroKind
	param_list: tuple[PPToken]
	replacement: tuple[PPToken] | Callable[[tuple[PPToken]], tuple[PPToken]]
	line: PPLine

builtin_line = PPLine("", 0, 0, "<builtin>")
	
class PPExecutor:
	RESERVED_MACRO_NAMES = set("__has_c_attribute __has_include __has_embed __VA_ARGS__ __VA_OPT__ __DATE__ __FILE__ __LINE__ __STDC__ __STDC_HOSTED__ __STDC_UTF_16__ __STDC_UTF_32__ __STDC_VERSION__ __TIME__ __STDC_ISO_10646__ __STDC_MB_MIGHT_NEQ_WC__ __STDC_ANALYZABLE__ __STDC_IEC_60559__ __STDC_IEC_559__ __STDC_IEC_60559_DFP__ __STDC_IEC_60559_COMPLEX__ __STDC_IEC_60559_TYPES__ __STDC_IEC_559_COMPLEX__ __STDC_LIB_EXT1__ __STDC_NO_ATOMICS__ __STDC_NO_COMPLEX__ __STDC_NO_THREADS__ __STDC_NO_VLA__".split())
	SYSTEM_INCLUDE_DIRS = ("/usr/local/include", "/usr/include")
	
	def __init__(self, tok_stream, namespace=None):
		now = datetime.now()
		self.tok_stream = tok_stream
		if namespace is None:
			self.makeNamespace()
		self.date_str = now.strftime("%b %d %Y")
		if self.date_str[len("Mmm ")] == '0':
			self.date_str = self.date_str[:len("Mmm ")] + " " + self.date_str[len("Mmm 0"):]
		self.time_str = now.strftime("%H:%M:%S")
		self.curr_line = None
	
	def _validateMacro(self, mdef):
		if self.curr_line is not None and mdef.name in PPExecutor.RESERVED_MACRO_NAMES:
			raise ParseError(f"{mdef.line.logical_file}:{mdef.line.logical_line} Macro name \"{mdef.spelling}\" is reserved")
		if isinstance(mdef.replacement, Callable):
			if mdef.kind == PPMacroKind.builtin_object and mdef.param_list:
				raise ParseError(f"{mdef.line.logical_file}:{mdef.line.logical_line} Builtin object-like macro \"{mdef.spelling}\" should not have parameter list")
			elif mdef.kind not in (PPMacroKind.builtin_object, PPMacroKind.builtin_function):
				raise ParseError(f"{mdef.line.logical_file}:{mdef.line.logical_line} Macro \"{mdef.spelling}\" with replacement callback should be builtin")
		if mdef.kind in (PPMacroKind.builtin_object, PPMacroKind.object_like):
			if mdef.param_list:
				raise ParseError(f"{mdef.line.logical_file}:{mdef.line.logical_line} Object-like macro \"{mdef.spelling}\" should not have parameter list")
		elif mdef.kind in (PPMacroKind.builtin_function, PPMacroKind.function_like):
			if len(set(tok.spelling for tok in self.param_list)) != len(self.param_list):
				raise ParseError(f"{mdef.line.logical_file}:{mdef.line.logical_line} Function-like macro \"{mdef.spelling}\" has duplicate parameter name")
		elif (old := self.definitions.get(self.spelling, None)) is not None:
			if (mdef.kind != old.kind or
			    len(mdef.param_list) != len(old.param_list) or
			    any(p.spelling != q.spelling or p.kind != q.kind for (p, q) in zip(mdef.param_list, old.param_list)) or
			    isinstance(mdef.replacement, Callable) or
			    isinstance(old.replacement, Callable) or
			    len(mdef.replacement) != len(mdef.replacement) or
			    any(p.spelling != q.spelling or p.kind != q.kind for (p, q) in zip(mdef.replacement, old.replacement))
			):
				raise ParseError(f"{mdef.line.logical_file}:{mdef.line.logical_line} Illegal redefinition of macro \"{mdef.spelling}\"\noriginally defined at {old.line.logical_file}:{old.line.logical_line}")
		return True
	
	def _registerMacro(self, *mdef_args):
		mdef = PPMacroDef(*mdef_args)
		self._validateMacro(mdef)
		self.definitions[mdef.spelling] = mdef
	
	def makeDateStrToken(self):
		return PPToken(self.date_str, PPTokenKind.string_literal, self.curr_line, "")
	
	def makeFileStrToken(self):
		if self.curr_line is not None:
			logical_file = self.curr_line.logical_file
		else:
			logical_file = "<text>"
		return PPToken(logical_file, PPTokenKind.string_literal, self.curr_line, "")
	
	def makeLineNumToken(self):
		if self.curr_line is not None:
			logical_line = self.curr_line.logical_line
		else:
			logical_line = 0
		return PPToken(str(logical_line), PPTokenKind.pp_number, self.curr_line, "")
	
	def makeTimeStrToken(self):
		return PPToken(self.time_str, PPTokenKind.string_literal, self.curr_line, "")
	
	def makeNamespace(self):
		self.definitions = {}
		self._registerMacro("__DATE__", PPMacroKind.builtin_object, (), self.makeDateStrToken, builtin_line)
		self._registerMacro("__FILE__", PPMacroKind.builtin_object, (), self.makeFileStrToken, builtin_line)
		self._registerMacro("__LINE__", PPMacroKind.builtin_object, (), self.makeLineNumToken, builtin_line)
		self._registerMacro("__STDC__", PPMacroKind.builtin_object, (), (PPToken("1", PPTokenKind.pp_number, None, ""),))
		# __STDC_HOSTED__ may also be 0
		self._registerMacro("__STDC_HOSTED__", PPMacroKind.builtin_object, (), (PPToken("1", PPTokenKind.pp_number, None, ""),))
		self._registerMacro("__STDC_UTF_16__", PPMacroKind.builtin_object, (), (PPToken("1", PPTokenKind.pp_number, None, ""),))
		self._registerMacro("__STDC_UTF_32__", PPMacroKind.builtin_object, (), (PPToken("1", PPTokenKind.pp_number, None, ""),))
		self._registerMacro("__STDC_VERSION__", PPMacroKind.builtin_object, (), (PPToken("202311L", PPTokenKind.pp_number, None, ""),))
		self._registerMacro("__TIME__", PPMacroKind.builtin_object, (), self.makeTimeStrToken, builtin_line)
		# begin OPTIONAL
		self._registerMacro("__STDC_ISO_10646__", PPMacroKind.builtin_object, (), (PPToken("1", PPTokenKind.pp_number, None, ""),))
		#self._registerMacro("__STDC_MB_MIGHT_NEQ_WC__", PPMacroKind.builtin_object, (), (PPToken("1", PPTokenKind.pp_number, None, ""),))
		self._registerMacro("__STDC_ANALYZABLE__", PPMacroKind.builtin_object, (), (PPToken("1", PPTokenKind.pp_number, None, ""),))
		self._registerMacro("__STDC_IEC_60559__", PPMacroKind.builtin_object, (), (PPToken("202311L", PPTokenKind.pp_number, None, ""),))
		self._registerMacro("__STDC_IEC_559__", PPMacroKind.builtin_object, (), (PPToken("1", PPTokenKind.pp_number, None, ""),))
		self._registerMacro("__STDC_IEC_60559_DFP__", PPMacroKind.builtin_object, (), (PPToken("202311L", PPTokenKind.pp_number, None, ""),))
		self._registerMacro("__STDC_IEC_60559_COMPLEX__", PPMacroKind.builtin_object, (), (PPToken("202311L", PPTokenKind.pp_number, None, ""),))
		self._registerMacro("__STDC_IEC_60559_TYPES__", PPMacroKind.builtin_object, (), (PPToken("202311L", PPTokenKind.pp_number, None, ""),))
		self._registerMacro("__STDC_IEC_559_COMPLEX__", PPMacroKind.builtin_object, (), (PPToken("1", PPTokenKind.pp_number, None, ""),))
		self._registerMacro("__STDC_LIB_EXT1__", PPMacroKind.builtin_object, (), (PPToken("202311L", PPTokenKind.pp_number, None, ""),))
		#self._registerMacro("__STDC_NO_ATOMICS__", PPMacroKind.builtin_object, (), (PPToken("1", PPTokenKind.pp_number, None, ""),))
		#self._registerMacro("__STDC_NO_COMPLEX__", PPMacroKind.builtin_object, (), (PPToken("1", PPTokenKind.pp_number, None, ""),))
		self._registerMacro("__STDC_NO_THREADS__", PPMacroKind.builtin_object, (), (PPToken("1", PPTokenKind.pp_number, None, ""),))
		#self._registerMacro("__STDC_NO_VLA__", PPMacroKind.builtin_object, (), (PPToken("1", PPTokenKind.pp_number, None, ""),))
	
	def isMacroDefined(self, spelling):
		return spelling in ("__has_include", "__has_embed", "__has_c_attribute") or spelling in self.definitions
	
	def _macroSubSingle(self, line_tokens):
		...
	
	def _applyPpExprs(self, line_tokens):
		
	
	def _processIf(self, line_tokens):
		if self.if_layers and self.if_layers[-1] != 2:
			self.if_layers.append(1)
			return
		line_tokens = list(self._applyPpExprs(line_tokens))
	
	def _processIfdef(self, line_tokens, invert=False):
		if self.if_layers and self.if_layers[-1] != 2:
			self.if_layers.append(1)
		elif len(line_tokens) != 3 or line_tokens[2].kind != PPTokenKind.identifier:
			raise ParseError(f"{self.curr_line.logical_file}:{self.curr_line.logical_line} Invalid argument to #ifdef")
		else:
			spelling = line_tokens[2].spelling
			status = self.isMacroDefined(spelling) == invert
			self.if_layers.append(2*int(status))
	
	def _processElif(self, line_tokens):
		...
	
	def _processElifdef(self, line_tokens, invert=False):
		if not self.if_layers:
			raise ParseError(f"{self.curr_line.logical_file}:{self.curr_line.logical_line} #elifdef outside of #if")
		elif self.if_layers[-1]:
			self.if_layers[-1] = 1
		elif len(line_tokens) != 3 or line_tokens[2].kind != PPTokenKind.identifier:
			raise ParseError(f"{self.curr_line.logical_file}:{self.curr_line.logical_line} Invalid argument to #elifdef")
		else:
			spelling = line_tokens[2].spelling
			status = (spelling in ("__has_include", "__has_embed", "__has_c_attribute") or spelling in self.definitions) == invert
			self.if_layers[-1] = 2*int(status)
	
	def _processElse(self, line_tokens):
		if len(line_tokens) != 2:
			raise ParseError(f"{self.curr_line.logical_file}:{self.curr_line.logical_line} #else must be followed by newline")
		if not self.if_layers:
			raise ParseError(f"{self.curr_line.logical_file}:{self.curr_line.logical_line} #else outside of #if")
		if self.if_layers[-1]:
			self.if_layers[-1] = 1
		else:
			self.if_layers[-1] = 2
	
	def _processEndif(self, line_tokens):
		if len(line_tokens) != 2:
			raise ParseError(f"{self.curr_line.logical_file}:{self.curr_line.logical_line} #endif must be followed by newline")
		if not self.if_layers:
			raise ParseError(f"{self.curr_line.logical_file}:{self.curr_line.logical_line} #endif outside of #if")
		self.if_layers.pop()
	
	def _processInclude(self, line_tokens):
		if self.if_layers and self.if_layers[-1] != 2:
			return #in inactive #if group
		line_tokens = self._macroSubSingle(line_tokens[2:])
		if len(line_tokens) > 2:
			if line_tokens[0].spelling == "<" and line_tokens[-1].spelling == ">":
				header_name = "".join(tok.spelling for tok in line_tokens)
			elif line_tokens[0].spelling == "\"" and line_tokens[-1].spelling == "\"":
				header_name = "".join(tok.spelling for tok in line_tokens)
			else:
				raise ParseError(f"{self.curr_line.logical_file}:{self.curr_line.logical_line} #include did not resolve to \"header_name\" or <header_name>")
		elif len(line_tokens) == 1:
			if line_tokens[0].kind == PPTokenKind.string_literal:
				if line_tokens[0].spelling[::len(line_tokens[0].spelling)-1] != "\"\"":
					raise ParseError(f"{self.curr_line.logical_file}:{self.curr_line.logical_line} #include cannot accept a prefixed/suffixed string literal")
			elif line_tokens[0].kind != PPTokenKind.header_name:
				raise ParseError(f"{self.curr_line.logical_file}:{self.curr_line.logical_line} Invalid #include argument")
			header_name = line_tokens[0].spelling
		else:
			raise ParseError(f"{self.curr_line.logical_file}:{self.curr_line.logical_line} #include needs exactly one argument")
		if header_name[0] == "\"":
			search_dirs = itt.chain((Path.cwd(),), PPExecutor.SYSTEM_INCLUDE_DIRS)
		else:
			search_dirs = PPExecutor.SYSTEM_INCLUDE_DIRS
		header_path = Path(header_name[1:-1])
		if not header_path.is_absolute():
			for d in search_dirs:
				p = d/header_path
				if p.is_file():
					header_path = p
					break
			else:
				raise ParseError(f"Unknown header {header_name}")
		elif not header_path.is_file():
			raise ParseError(f"Unknown header {header_name}")
		with open(header_path) as f:
			ppp = PPParser(f)
			ppp.combineLines()
		token_gen = ppp.genTokens()
		ppx = PPExecutor(token_gen, self.namespace)
		return ppx.process()
	
	def _processEmbed(self, line_tokens):
		...
	
	def _processDefine(self, line_tokens):
		...
	
	def _processUndef(self, line_tokens):
		...
	
	def _processLine(self, line_tokens):
		...
	
	def _processDiagnostic(self, line_tokens):
		...
	
	def _processPragma(self, line_tokens):
		...
	
	def _processEmpty(self, line_tokens):
		...
	
	def _processNonDirective(self, line_tokens):
		...
	
	def _processTextLine(self, line_tokens):
		...
	
	def process(self):
		self.if_layers = []
		def tokLine(tok):
			if tok.line is None:
				return self.curr_line.logical_line
			return tok.line.logical_line
		for _, line_tokens in itt.groupby(self.tok_stream, tokLine):
			line_tokens = list(line_tokens)
			if not line_tokens:
				continue
			if line_tokens[0].line is not None:
				self.curr_line = line_tokens[0].line
			if line_tokens[0].spelling == "#":
				if len(line_tokens) == 1:
					self._processEmpty(line_tokens)
				elif line_tokens[1].spelling == "if":
					self._processIf(line_tokens)
				elif line_tokens[1].spelling == "ifdef":
					self._processIfdef(line_tokens)
				elif line_tokens[1].spelling == "ifndef":
					self._processIfdef(line_tokens, invert=True)
				elif line_tokens[1].spelling == "elif":
					cond = self._processIf(line_tokens)
				elif line_tokens[1].spelling == "elifdef":
					self._processIfdef(line_tokens)
				elif line_tokens[1].spelling == "elifndef":
					self._processIfdef(line_tokens, invert=True)
				elif line_tokens[1].spelling == "else":
					self._processElse(line_tokens)
				elif line_tokens[1].spelling == "endif":
					self._processEndif(line_tokens)
				elif line_tokens[1].spelling == "include":
					self._processInclude(line_tokens)
				elif line_tokens[1].spelling == "embed":
					self._processEmbed(line_tokens)
				elif line_tokens[1].spelling == "define":
					self._processDefine(line_tokens)
				elif line_tokens[1].spelling == "undef":
					self._processUndef(line_tokens)
				elif line_tokens[1].spelling == "line":
					self._processLine(line_tokens)
				elif line_tokens[1].spelling == "error":
					self._processDiagnostic(line_tokens, error=True)
				elif line_tokens[1].spelling == "warning":
					self._processDiagnostic(line_tokens)
				elif line_tokens[1].spelling == "pragma":
					self._processPragma(line_tokens)
				else:
					self._processNonDirective(line_tokens)
			else:
				self._processTextLine(line_tokens)

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("Please specify file to parse!")
		sys.exit(1)
	with open(sys.argv[1]) as f:
		ppp = PPParser(f)
		ppp.combineLines()
	pp_tokens = list(ppp.genTokens())
	curr_line = None
	for tok in pp_tokens:
		if tok.line is not curr_line:
			curr_line = tok.line
			if curr_line is not None:
				print(f"LINE #{curr_line.logical_line}")
		#if tok.leading_ws:
		#	print(f"  WHITESPACE: {tok.leading_ws}")
		print(f"  {tok.kind}: {tok.spelling}")

