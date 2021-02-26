doc:
	asciidoctor-latex -b latex doc.adoc
	latexmk -pdf doc
	latexmk -c
	rm doc.tex

clean:
	latexmk -C doc
