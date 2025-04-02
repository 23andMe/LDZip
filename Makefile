.PHONY: all cpp r-package clean

all: cpp r-package

cpp:
	$(MAKE) -C cpp build

r-package:
	$(MAKE) -C R install

clean:
	$(MAKE) -C cpp clean || true
	$(MAKE) -C R clean || true
