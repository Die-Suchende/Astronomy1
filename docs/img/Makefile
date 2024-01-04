all: $(addsuffix .pdf, $(basename $(wildcard SP2_lecture*.tex)))

clean: tidy
	latexmk -pdf -C SP2_lecture*.tex

tidy:
	latexmk -pdf -c SP2_lecture*.tex
	- rm SP2_lecture*.{snm,nav}

%.pdf: %.tex
	latexmk -pdf $<
