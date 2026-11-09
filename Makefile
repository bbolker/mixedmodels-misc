ALL: docs/glmmFAQ.html docs/glmmbib.html

docs/glmmbib.html: glmmbib.html
	mv glmmbib.html docs/

glmmbib.html: glmm.bib
	. ./mkbibhtml

docs/glmmFAQ.html: glmmFAQ.html
	mv glmmFAQ.html docs/

%.html: %.rmd
	Rscript -e "library(\"rmarkdown\"); render(\"$<\")"

%.html: %.qmd
	quarto render $<

%.md: %.rmd
	Rscript -e "library(\"knitr\"); knit(\"$<\")"

%.tex: %.md
	pandoc -s -S -t latex -V documentclass=tufte-handout $*.md -o $*.tex

%.pdf: %.tex
	pdflatex --interaction=nonstopmode $*

clean:
	rm -f *.log *.aux *.md *.out texput.log

