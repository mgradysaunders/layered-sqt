TITLE=Layered-SQT
$(TITLE).pdf: $(TITLE).tex References.bib
	xelatex $(TITLE).tex && biber $(TITLE) && \
	xelatex $(TITLE).tex && xelatex $(TITLE).tex

.PHONY: clean
clean:
	rm -f $(TITLE).aux
	rm -f $(TITLE).bbl
	rm -f $(TITLE).bcf
	rm -f $(TITLE).blg
	rm -f $(TITLE).log
	rm -f $(TITLE).out
	rm -f $(TITLE).run.xml
