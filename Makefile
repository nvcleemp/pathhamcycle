
SOURCES = pathhamcycle.c bitset.h boolean.h\
          Makefile LICENSE README.md

all: build/pathhamcycle

clean:
	rm -rf build
	rm -rf dist

build/pathhamcycle: pathhamcycle.c
	mkdir -p build
	cc -o $@ -O4 -Wall $^

sources: dist/pathhamcycle-sources.zip dist/pathhamcycle-sources.tar.gz

dist/pathhamcycle-sources.zip: $(SOURCES)
	mkdir -p dist
	zip dist/pathhamcycle-sources $(SOURCES)

dist/pathhamcycle-sources.tar.gz: $(SOURCES)
	mkdir -p dist
	tar czf dist/pathhamcycle-sources.tar.gz $(SOURCES)
