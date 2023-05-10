#!/usr/bin/env python3
import os, subprocess, sys, json, argparse

def runNoRedTest(test_path, test_argv, log_name):
	pipe_out, pipe_in = os.pipe()
	if args.valgrind:
		popen_args = ["valgrind", "--leak-check=full", "--track-origins=yes", "--read-inline-info=yes", "--read-var-info=yes", "--enable-debuginfod=yes", test_path] + test_argv
	else:
		popen_args = [test_path] + test_argv
	child = subprocess.Popen(popen_args, stdout=pipe_in, stderr=pipe_in)
	tee = subprocess.Popen(["tee", log_name], stdin=pipe_out)
	child.wait()
	os.close(pipe_in)
	tee.wait()
	os.close(pipe_out)
	status = subprocess.call(["grep", "-qF", "\033[1;31m", log_name])
	os.remove(log_name)
	return status == 1

def runLogDiffTests(test_path, i, test_argv, log_name):
	pipe_out, pipe_in = os.pipe()
	if args.valgrind:
		popen_args = ["valgrind", "--leak-check=full", "--track-origins=yes", "--read-inline-info=yes", "--read-var-info=yes", "--enable-debuginfod=yes", test_path] + test_argv
	else:
		popen_args = [test_path] + test_argv
	child = subprocess.Popen(popen_args, stdout=pipe_in, stderr=pipe_in)
	tee = subprocess.Popen(["tee", log_name], stdin=pipe_out)
	child.wait()
	os.close(pipe_in)
	tee.wait()
	os.close(pipe_out)
	expected_log_name = os.path.join(args.cfg_dir, f"{test_name}.{i}.log")
	status = subprocess.call(["diff", "-q", expected_log_name, log_name])
	os.remove(log_name)
	return status == 0

if __name__ == "__main__":
	arg_parser = argparse.ArgumentParser(description="Automatically run tests stored in a json file")
	arg_parser.add_argument("cfg_dir", help="directory containing tests.json and expected output logs if used")
	arg_parser.add_argument("bin_dir", help="directory containing test binaries")
	arg_parser.add_argument("-g", "--valgrind", action="store_true", help="run tests under valgrind")
	args = arg_parser.parse_args()

	with open(os.path.join(args.cfg_dir, "tests.json"), "r") as test_json:
		test_cfg = json.load(test_json)

	passed = 0
	tested = 0
	for test_name, test_spec in test_cfg.items():
		no_red_tests = test_spec.get("no_red_tests", ())
		test_path = os.path.join(args.bin_dir, test_name)
		log_name = os.path.join(args.cfg_dir, "tmp.log")
		for test_argv in no_red_tests:
			if runNoRedTest(test_path, test_argv, log_name):
				passed += 1
			tested += 1
		log_diff_tests = test_spec.get("log_diff_tests", ())
		for i, test_argv in enumerate(log_diff_tests):
			if runLogDiffTests(test_path, i, test_argv, log_name):
				passed += 1
			tested += 1

	print(f"{passed}/{tested} tests passed.")
	if passed == tested:
		print("\033[1;32mtest.py suite passed all tests!\033[0m")
	sys.exit(passed != tested)

