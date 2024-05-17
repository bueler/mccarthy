all: slides.pdf

slides.pdf:
	pdflatex slides
	bibtex slides
	pdflatex slides
	pdflatex slides

.PHONY: clean

clean:
	@rm -f *~ *.out *.aux *.log *.bbl *.blg *.snm *.toc *.nav *.vrb *.synctex.gz
	@rm -rf tmp/
