SRCDIR=../

%.html: ${SRCDIR}/%.[Rr]md
	Rscript -e "rmarkdown::render(\"$<\",output_dir='.')"

%.pdf: ${SRCDIR}/%.rmd
	Rscript -e "rmarkdown::render(\"$<\",output_dir='.',output_format=\"pdf_document\")"

glmmbib.html: ../glmm.bib
	cd ..; ./mkbibhtml; cp -r glmmbib.html gh-pages
clean:
	rm -f *.log *.aux *.md *.out texput.log


