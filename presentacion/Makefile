ROBLEMS:=${shell ls *.tex}

TARGET_FORMAT=pdf

.PHONY: all

all:			${PROBLEMS:.tex=.${TARGET_FORMAT}}

%.pdf:			%.tex force
				latexmk -xelatex $<

force:			;

.PHONY: clean

clean:
			rm -f *.bbl *.xml *.snm *.nav *.out *.fdb_latexmk *.vrb *.table *.gnuplot *-blx.bib *.log *.toc *.fls *.aux
