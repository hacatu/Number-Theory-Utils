exe:
	mkdir -p obj/test bin/test
	clang -target x86_64-pc-windows-wsdk -m64 -isystem $(HOME)/wsdk/10/Include/10.0.22621.0/ucrt -isystem $(HOME)/wsdk/10/Include/10.0.22621.0/um -isystem $(HOME)/wsdk/10/Include/10.0.22621.0/shared -fmsc-version=1936 -fms-extensions -fms-compatibility -fdelayed-template-parsing -DWIN32 -c ../../src/test/test_windows.c -o obj/test/test_windows.o
	lld -flavor link /libpath:$(HOME)/wsdk/10/Lib/10.0.22621.0/ucrt/x64 /libpath:$(HOME)/wsdk/10/Lib/10.0.22621.0/um/x64 libcmt.lib /subsystem:console /out:bin/test/test_windows.exe obj/test/test_windows.o

