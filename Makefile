.PHONY: all roxygenise install clean

all: roxygenise install

roxygenise:
	R -e 'roxygen2::roxygenise(".")'
install: 
	R CMD INSTALL . 
clean:
	rm -rf man NAMESPACE
	



