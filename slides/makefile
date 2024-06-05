all: slides.pdf exercises.pdf

slides.pdf: slides.tex slides.bib
	xelatex slides
	bibtex slides
	xelatex slides
	xelatex slides

exercises.pdf: exercises.tex
	pdflatex exercises

.PHONY: clean

clean:
	@rm -f *~ *.out *.aux *.log *.bbl *.blg *.snm *.toc *.nav *.vrb *.synctex.gz
	@rm -rf tmp/
