.PHONY: clean

FILES=*.cs Core/*.cs

build:	ProjectB.csproj
	msbuild ProjectB.csproj

run:
	mono ./bin/Debug/ProjectB.exe

mcs:	${FILES}
	mcs ${FILES}

clean:
	rm -rf bin/ obj/ *.exe
