.PHONY: all roxygenise install clean

all: roxygenise install

roxygenise:
	R -e 'roxygen2::roxygenise(".")'
install: 
	R CMD INSTALL .
build_tgz:
	R CMD build .
clean:
	rm -rf man NAMESPACE
	



