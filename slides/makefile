all: slides.pdf

slides.pdf: slides.tex slides.bib
	xelatex slides
	bibtex slides
	xelatex slides
	xelatex slides

.PHONY: clean

clean:
	@rm -f *~ *.out *.aux *.log *.bbl *.blg *.snm *.toc *.nav *.vrb *.synctex.gz
	@rm -rf tmp/
