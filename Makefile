spa2.pdf: spa2.Rnw
	R --vanilla -e "library(knitr); knit2pdf('spa2.Rnw');"
	pdflatex spa2.tex

clean:
	rm -f spa2.aux spa2.log spa2.out spa2.tex
	rm -rf figure
	rm -f *~
