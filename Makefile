SRCDIR=../

%.html: ${SRCDIR}/%.rmd
	Rscript -e "rmarkdown::render(\"$<\",output_dir='.')"

%.pdf: ${SRCDIR}/%.rmd
	Rscript -e "rmarkdown::render(\"$<\",output_dir='.',output_format=\"pdf_document\")"

clean:
	rm -f *.log *.aux *.md *.out texput.log


