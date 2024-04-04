#master makefile for OpenExtrap

all: mkdirs
	cd emmod		; $(MAKE) install

mkdirs:
	-mkdir bin

clean:
	cd emmod		; $(MAKE) $@

realclean:
	cd emmod		; $(MAKE) $@
	rm -f bin/*
