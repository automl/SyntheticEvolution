CC=clang++
CFLAGS=-O3 -ffast-math

USalign: x86_app arm_app
	lipo -create -output $@ x86_app arm_app

x86_app: USalign.cpp
	${CC} ${CFLAGS} USalign.cpp -o x86_app -target x86_64-apple-macos10.12

arm_app: USalign.cpp
	${CC} ${CFLAGS} USalign.cpp -o arm_app -target arm64-apple-macos11
