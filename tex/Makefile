all:minigraph.pdf

minigraph.bbl:minigraph.bib minigraph.tex
	pdflatex minigraph; bibtex minigraph

minigraph.pdf:minigraph.tex minigraph.bbl
	pdflatex minigraph; pdflatex minigraph

lite:minigraph.tex
	pdflatex minigraph

clean:
	rm -f minigraph.bbl minigraph.log minigraph.aux
