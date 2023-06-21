import functools as func
import itertools as itt
from collections import defaultdict, deque
import urllib.request
import urllib.parse
import hashlib
import json
import os
import shutil
import tempfile
import subprocess
import multiprocessing as mp
import re
from zipfile import ZipFile
from pathlib import Path
from typing import Any, Optional
from collections.abc import Iterable, MutableSequence

SDK_VERSION = "10.0.19041"
VS_VERSION = "14.30.17.0"
VS_MAJOR = "17"
VS_CHANNEL = "release"

packages = [
	f"Win10SDK_{SDK_VERSION}",
	f"Microsoft.VisualStudio.Component.VC.{VS_VERSION}.x86.x64",
	"Microsoft.VisualStudio.Workload.VCTools"]

def getPkgOrd(pkg: Any) -> tuple[bool, bool]:
	# prioritize packages first by whether or not they are x86_64, then by whether or not they are english
	return (pkg.get("chip", "").lower() == "x64", pkg.get("language", "").lower().startswith("en-"))

class Downloader:
	BUILD_TOOLS_VERSION_RE = re.compile(r"^Microsoft\.VisualStudio\.Component\.VC\.(.*)\.x86\.x64$", re.IGNORECASE)
	WIN_SDK_VERSION_RE = re.compile(r"^Win(\d+)SDK_(.*)", re.IGNORECASE)
	
	def __init__(self, manifest_path: Optional[os.PathLike] = None, warn: bool = True):
		self.manifest_path = manifest_path
		self.warn = warn
		self.manifest = None
		self.package_db = None
		self._initManifest()
		self._initPackageDb()
	
	def _initManifest(self):
		if self.manifest_path is None:
			with urllib.request.urlopen(f"https://aka.ms/vs/{VS_MAJOR}/{VS_CHANNEL}/channel") as f:
				manifest = json.load(f)
			# find the packages manifest entry in the channelItems manifest
			url = next(it["payloads"][0]["url"] for it in manifest["channelItems"] if it["type"] == "Manifest")
			with urllib.request.urlopen(url) as f:
				data = f.read()
				self.manifest_path = f"{manifest['info']['productSemanticVersion']}.vsman"
				with open(self.manifest_path, "wb") as o:
					o.write(data)
				self.manifest = json.loads(data)
		else:
			with open(self.manifest_path) as f:
				self.manifest = json.load(f)

	def _initPackageDb(self):
		res = defaultdict(list)
		for pkg in self.manifest.get("packages", ()):
			res[pkg["id"].lower()].append(pkg)
		for k in res.keys():
			res[k].sort(key=getPkgOrd, reverse=True)
		self.package_db = res

	def findPackage(self, pkg_id: str, chip: Optional[str]) -> Any:
		if not (pkgs := self.package_db.get(pkg_id.lower(), ())):
			if self.warn:
				print(f"WARNING: {pkg_id} not found!")
			return None
		if chip is not None:
			chip = chip.lower()
			return next((pkg for pkg in pkgs if pkg.get("chip", "").lower() == chip), pkgs[0])
		return pkgs[0]

	def _tryGetFirstPayload(self, pkg_id: str, out: MutableSequence[Any], chip: Optional[str], warn: bool = True):
		if (pkg := self.findPackage(pkg_id, chip)) is not None and (payloads := pkg.get("payloads", ())):
			out.append(payloads[0])
	
	def _getCrtPayloads(self, chips: set[str], variants: set[str], include_atl: bool, out: MutableSequence[Any]):
		if (build_tools := self.findPackage("Microsoft.VisualStudio.Product.BuildTools")) is None:
			return
		crt_version = max(
			(Version(version) for k in build_tools.get("dependencies", {})
			if (version := self.BUILD_TOOLS_VERSION_RE.fullmatch(k)) is not None),
			default=None)
		if crt_version is None:
			return
		self._tryGetFirstPayload(f"Microsoft.VC.{crt_version}.CRT.Headers.base", out)
		variants.add("store")
		for chip, variant in itt.product(chips, variants):
			pkg_id = f"Microsoft.VC.{crt_version}.CRT.{chip}.{variant}"
			if variant != "store" and "spectre" in variants:
				pkg_id += ".spectre"
			pkg_id += ".base"
			self._tryGetFirstPayload(pkg_id, out)
		if include_atl:
			self._tryGetFirstPayload(f"Microsoft.VC.{crt_version}.ATL.Headers.base", out)
			for chip, want_spectre in itt.product(chips, (False, True) if "spectre" in variants else (False,)):
				pkg_id = f"Microsoft.VC.{crt_version}.ATL.{chip}"
				if want_spectre:
					pkg_id += ".spectre"
				pkg_id += ".base"
				self._tryGetFirstPayload(pkg_id, out)
	
	def _getWsdkPayloads(self, chips: set[str], variants: set[str], out: MutableSequence[Any]):
		win_sdk_version_major, win_sdk_version_full = max(
			((m.group(1), Version(m.group(2))) for k in package_db
			if (m := self.WIN_SDK_VERSION_RE.match(k)) is not None),
			default=(None, None))
		if win_sdk_version_full is None:
			return
		win_sdk_version = f"Win{win_sdk_version_major}SDK_{win_sdk_version_full}"
		if (pkg := self.findPackage(win_sdk_version)) is None:
			return
		pl = next(
			(pl for pl in pkg.get("payloads", ())
			if pl.get("fileName", "").endswith("Windows SDK Desktop Headers x86-x86_en-us.msi")),
			None)
		if pl is None:
			return
		out.append({
			"fileName": f"{pkg['id']}_headers.msi",
			"sha256": pl["sha256"],
			"url": pl["url"],
			# install_size
			"type": "SdkHeaders",
			"variant": None,
			"chip": None
		})

	def getMinimalPayloads(self, chips: Iterable[str], variants: Iterable[str], include_atl: bool) -> list[Any]:
		chips = set(map(chips, str.lower))
		variants = set(map(variants, str.lower))
		res = []
		self._getCrtPayloads(chips, variants, include_atl, res)
		self._getWsdkPayloads(chips, variants, res)

	def getDependentTree(package_db: dict[str, Any], pkg_ids: Iterable[str], chip: Optional[str] = None) -> list[Any]:
		visited = {*zip(pkg_ids, itt.repeat(chip))}
		count = len(visited)
		targets = deque(pkg for (pkg_id, chip) in visited if (pkg := findPackage(package_db, pkg_id, chip)) is not None)
		res = list(targets)
		while targets:
			pkg = targets.pop()
			for dep_id, dep in pkg.get("dependencies", {}).items():
				ty = dep.get("type", None) if isinstance(dep, dict) else None
				if ty in ("Optional", "Recommended"):
					continue
				chip = dep.get("chip", None) if isinstance(dep, dict) else None
				visited.add((dep_id, chip))
				if len(visited) != count and (pkg := findPackage(package_db, dep_id, chip)) is not None:
					count += 1
					targets.appendleft(pkg)
					res.append(pkg)
		return res

	def sumInstalledSize(pkgs: Iterable[Any]) -> int:
		return sum(sum(pkg.get("installSizes", {}).values()) for pkg in pkgs)

	def sumDownloadSize(pkgs: Iterable[Any]) -> int:
		return sum(sum(pl.get("size", 0) for pl in pkg.get("payloads", ())) for pkg in pkgs)

	def formatSize(n: int) -> str:
		BYTE_UNITS = {1: 'bytes', 1024: 'KB', 1024**2: 'MB', 1024**3: 'GB', 1024**4: 'TB'}
		for i, (unit_count, unit_name) in enumerate(BYTE_UNITS.items()):
			if 1024*unit_count > n or i + 1 == len(BYTE_UNITS):
				return f"{n/unit_count:.1f} {unit_name}"

	def sha256File(filename: os.PathLike) -> str:
		hasher = hashlib.sha256()
		with open(file, "rb") as f:
			while len(b := f.read(4096)):
				hasher.update(b)
			return hasher.hexdigest()

	def getPkgDirname(pkg: Any) -> str:
		return f"{pkg['id']}-{pkg.get('version','')}-{pkg.get('chip','')}"

	def downloadPackages(pkgs: Iterable[Any], path: Path, check_hash: bool = True, pool: Optional[mp.Pool] = None) -> int:
		tasks = []
		os.makedirs(path, exist_ok=True)
		for pkg in pkgs:
			if "payloads" not in pkg:
				continue
			pkg_path = path / getPkgDirname(pkg)
			os.makedirs(pkg_path, exist_ok=True)
			for pl in pkg["payloads"]:
				name = Path(pl["fileName"]).name
				dest = pkg_path / name
				if pool is None:
					tasks.append(downloadPayload(pl, dest, check_hash))
				else:
					tasks.append(pool.apply_async(downloadPayload, (pl, dest, check_hash)))
		return sum(tasks) if pool is None else sum(task.get() for task in tasks)

	def downloadPayload(pl: Any, dest: Path, check_hash: bool, attempts: int = 5) -> int:
		path_tail = Path(dest.parts[-2:])
		for attempt in range(attempts):
			try:
				if dest.is_file():
					if not check_hash or (sha := pl.get("sha256", None)) is None or sha256File(dest).lower() == sha.lower():
						return 0
					print(f"Existing file {path_tail} has invalid hash, attempting to replace...", flush=True)
				size = pl.get("size", 0)
				print(f"Downloading {path_tail} ({formatSize(size)})", flush=True)
				with urllib.request.urlopen(pl["url"]) as f, open(dest, "wb") as o:
					o.write(f.read())
				if (sha := pl.get("sha256", None)) is not None and sha256File(dest).lower() != sha.lower():
					if check_hash:
						raise OSError(f"Newly downloaded file {path_tail} has invalid hash, aborting")
					else:
						print("WARNING: Newly downloaded file {path_tail} has invalid hash", flush=True)
				return size
			except OSError as e:
				if attempt + 1 == attempts:
					raise
				print(e, flush=True)

	def mergeTrees(src: Path, dest: Path):
		if not src.is_dir():
			return
		elif not dest.is_dir():
			shutil.move(src, dest)
			return
		dest_subs = {(n := dest_sub.name).lower(): n for dest_sub in dest.iterdir()}
		for src_sub in src.iterdir():
			sub_name = src_sub.name
			if src_sub.is_dir():
				if (dest_name := dest_subs.get(sub_name, None)) is None:
					dest_subs[sub_name] = (dest_name := sub_name)
				mergeTrees(src_sub, dest / dest_name)
			else:
				shutil.move(src_sub, dest / sub_name)

	def unzipFiltered(z: ZipFile, dest: Path):
		with tempfile.TemporaryDirectory() as temp:
			for f in z.infolist():
				name = urllib.parse.unquote(f.filename)
				if (idx := name.rfind("/")) != -1:
					os.makedirs(dest / name[:idx], exist_ok=True)
				shutil.move(z.extract(f, temp), dest / name)

	def unpackVsix(filename: Path, dest: Path, listing: Path):
		with tempfile.TemporaryDirectory() as temp, ZipFile(file) as z:
			unzipFiltered(z, temp)
			with open(listing, "w") as o:
				for n in z.namelist():
					print(n, file=o)
			if (contents := temp / "Contents").exists():
				mergeTrees(contents, dest)

	def unpackWin10SDK(src: Path, payloads: Iterable[Any], dest: Path):
		for pl in payloads:
			if (name := Path(pl["fileName"]).name).endswith(".msi"):
				msifile = src / name
				cmd = ("msiextract", "-C", dest, msifile)
				with open(dest / f"WinSDK-{name}.log") as o:
					subprocess.check_call(cmd, stdout=o)

	def extractPackages(pkgs: Iterable[Any], src: Path, dest: Path):
		os.makedirs(dest, exist_ok=True)
		for pkg in pkgs:
			pkg_dirname = getPkgDirname(pkg)
			pkg_src = src / pkg_dirname
			match (ty := pkg.get("type", None)):
				case "Component" | "Workload" | "Group": continue
				case "Vsix":
					print(f"Unpacking {p['id']}")
					for pl in pkg.get("payloads", ()):
						unpackVsix(pkg_src / Path(pl["fileName"]).name, dest, pkg / f"{pkg_dirname}.listing.txt")
				case _ if pkg["id"].startswith("Win10SDK") or pkg["id"].startswith("Win11SDK"):
					print(f"Unpacking {p['id']}")
					unpackWin10SDK(pkg_src, pkg.get("payloads", ()), dest)
				case _:
					print(f"Skipping unpacking {pkg['id']} of type {ty}")

	def moveWsdk(src: Path, dest: Path):
		os.makedirs(dest / "kits")
		mergeTrees(src / "VC", dest / "VC")
		kits_path = src / "Program Files" / "Windows Kits" / "10"
		mergeTrees(kits_path, dest / "kits" / "10")
		mergeTrees(src / "DIA SDK", dest / "DIA SDK")

if __name__ == "__main__":
	downloader = Downloader("17.6.2+33723.286.vsman")

