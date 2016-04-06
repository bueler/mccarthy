# Generic make file for LaTeX: requires GNU make

TEXFILE	= ideas2014.tex

all: $(TEXFILE:.tex=.pdf)

# link ice-bib.bib needed for bibtex to work  (to that file in pism-dev/doc/)

%.pdf: %.tex ice-bib.bib
	pdflatex $<
	bibtex $(<:.tex=.aux)
	pdflatex $<
	pdflatex $<

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.synctex.gz *.dvi *.toc *.nav *.snm *~

