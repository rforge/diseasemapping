all: bym.pdf

clean:
	rm -rf *.tex *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz

bym.pdf: bym.tex
	xelatex bym;
	xelatex bym;
	'$(R_HOME)/bin/Rscript' -e "tools::compactPDF('bym.pdf', gs_quality='ebook')" 


