all: lecture.pdf

# FIXME: much to do here
# FIXME need to do "make" in karthaus/, and add links to all graphics, for this too work
lecture.pdf:  lecture.tex style.tex preamble.tex sia.tex masscont.tex ssa.tex advanced.tex
	pdflatex lecture
	pdflatex lecture

# a good way to get a copy for print with multiple pages per sheet is:
printable.pdf:
	pdfnup --nup '1x2' --scale 0.9 --frame true --no-landscape -- lecture.pdf -
	pdftk lecture-nup.pdf ../extras/refs_stable.pdf cat output printable.pdf

.PHONY: clean

clean:
	@rm -f *.out *.aux *.log *.bbl *.blg *.snm *.toc *.nav *.vrb *.synctex.gz

