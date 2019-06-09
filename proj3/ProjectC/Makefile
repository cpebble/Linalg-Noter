.PHONY: clean

FILES=*.cs Core/*.cs

build:	ProjectC.csproj
	msbuild ProjectC.csproj

run:
	mono ./bin/Debug/ProjectC.exe

c:	${FILES}
	csc ${FILES}

r:
	mono *.exe

clean:
	rm -rf bin/ obj/ .idea/ *.exe
