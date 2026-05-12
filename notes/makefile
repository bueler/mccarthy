all: notes.pdf

notes.pdf: notes.tex
	pdflatex notes
	bibtex notes
	pdflatex notes
	pdflatex notes

.PHONY: clean

clean:
	@rm -f *~ *.out *.aux *.log *.blg *.snm *.toc *.nav *.vrb *.synctex.gz
