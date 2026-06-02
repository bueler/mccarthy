all: slides.pdf

slides.pdf: slides.tex slides.bib
	xelatex slides
	biber slides
	xelatex slides
	xelatex slides

.PHONY: clean

clean:
	@rm -f *~ *.out *.aux *.log *.bbl *.blg *.bcf *.snm *.toc *.nav *.vrb *.synctex.gz *.run.xml
	@rm -rf tmp/
