all: syntheticglacier.pdf

syntheticglacier.pdf: syntheticglacier.tex
	pdflatex syntheticglacier
	pdflatex syntheticglacier

.PHONY: clean

clean:
	@rm -f *~ *.out *.aux *.log *.blg *.snm *.toc *.nav *.vrb *.synctex.gz *.bbl

