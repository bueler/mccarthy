all: stokes.pdf

stokes.pdf: stokes.tex
	pdflatex stokes
	bibtex stokes
	pdflatex stokes
	pdflatex stokes

.PHONY: clean

clean:
	@rm -f *~ *.out *.aux *.log *.blg *.snm *.toc *.nav *.vrb *.synctex.gz *.bbl

