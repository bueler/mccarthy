all: exercises.pdf papersolns.pdf

exercises.pdf: exercises.tex
	pdflatex exercises.tex
	pdflatex exercises.tex

papersolns.pdf: papersolns.tex
	pdflatex papersolns.tex
	pdflatex papersolns.tex

.PHONY: clean

clean:
	@rm -f *~ *.out *.aux *.log *.bbl *.blg *.snm *.toc *.nav *.vrb *.synctex.gz
	@rm -rf tmp/
