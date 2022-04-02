all: notes.pdf

minputs := heat diffstag siaflat flowline ssaflowline
minputs := $(addsuffix .m,$(minputs))

mslim := $(subst .m,.slim.m,$(minputs))
msource := $(addprefix tmp/,$(mslim))

notes.pdf: $(msource) notes.tex
	pdflatex notes
	pdflatex notes

# the purpose of this nonsense is to put only the "meat" of codes into the
# PDF lectures.  thus we remove comments and print statements.  the flag
# "STRIPFROMHERE" can be used in a comment to remove everything after the flag
tmp/%.slim.m: ../mfiles/%.m
	mkdir -p tmp/
	sed -e '/STRIPFROMHERE/,$$ d' $< > $@.tmp
	sed -e 's/ *fprintf/%/' -e 's/ *%/%/' -e '/^%/ d' -e 's/%.*//' $@.tmp > $@
	rm -rf $@.tmp

.PHONY: clean

clean:
	@rm -f *~ *.out *.aux *.log *.blg *.snm *.toc *.nav *.vrb *.synctex.gz *.dvi
	@rm -rf tmp/
