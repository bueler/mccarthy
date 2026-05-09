all: doc.pdf

doc.pdf: doc.tex doc.bib
	pdflatex doc
	bibtex doc
	pdflatex doc
	pdflatex doc

.PHONY: clean

clean:
	@rm -f *~ *.out *.aux *.log *.blg *.snm *.toc *.nav *.vrb *.synctex.gz *.bbl

