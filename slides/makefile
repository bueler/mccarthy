all: slides.pdf

minputs := heat heatadapt diffstag siaflat flowline ssaflowline
minputs := $(addsuffix .m,$(minputs))

mslim := $(subst .m,.slim.m,$(minputs))
msource := $(addprefix tmp/,$(mslim))

slides.pdf: $(msource) anim/halfar-0.png slides.tex preamble.tex sia.tex masscont.tex ssa.tex
	pdflatex slides
	pdflatex slides

# the purpose of this nonsense is to put only the "meat" of codes into the
# PDF slides.  thus we remove comments and print statements.  the flag
# "STRIPFROMHERE" can be used in a comment to remove everything after the flag
tmp/%.slim.m: ../mfiles/%.m
	mkdir -p tmp/
	sed -e '/STRIPFROMHERE/,$$ d' $< > $@.tmp
	sed -e 's/ *fprintf/%/' -e 's/ *%/%/' -e '/^%/ d' -e 's/%.*//' $@.tmp > $@
	rm -rf $@.tmp

# un-tarring puts files halfar-?.png in directory which are older than the tar
# ... which explains using "touch" after creation
anim/halfar-0.png: anim/halfarfigs.tar.gz
	(cd anim/ && tar -xzvf halfarfigs.tar.gz && touch *.png)

.PHONY: clean

clean:
	@rm -f *~ *.out *.aux *.log *.bbl *.blg *.snm *.toc *.nav *.vrb *.synctex.gz
	@rm -rf tmp/ anim/halfar*png
