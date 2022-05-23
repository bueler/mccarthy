all: synglac.pdf

synglac.pdf: synglac.tex synglac.bib
	pdflatex synglac
	bibtex synglac
	pdflatex synglac
	pdflatex synglac

.PHONY: clean

clean:
	@rm -f *~ *.out *.aux *.log *.blg *.snm *.toc *.nav *.vrb *.synctex.gz *.bbl

