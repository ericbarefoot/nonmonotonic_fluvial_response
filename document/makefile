all: barefoot_manu_flood_regime.pdf clean

clean:
	rm -f *.aux *.aux *.lof *.log *.lot *.fls *.out *.toc *.fmt
	rm -f *.fot *.cb *.cb2 .*.lb *.bcf *.run.xml *.bbl *.blg *.tex.bak
	# rm -f *.tex

figures:
	cp ../../../figures/outputs/manuscript_figures/* figures/

.PHONY: all clean fix_biber figures

barefoot_manu_flood_regime.pdf: barefoot_manu_flood_regime.rnw  flood_regime_experiments.bib makefile
	Rscript --vanilla knitit.r barefoot_manu_flood_regime.rnw barefoot_manu_flood_regime.tex
	# make figures
	pdflatex	barefoot_manu_flood_regime.tex
	# biber 		barefoot_manu_flood_regime.bcf
	# biber 		barefoot_manu_flood_regime.bcf
	bibtex	barefoot_manu_flood_regime.aux
	pdflatex	barefoot_manu_flood_regime.tex
	pdflatex	barefoot_manu_flood_regime.tex
	# make 			barefoot_manu_flood_regime.docx
	make clean

barefoot_manu_flood_regime.docx:  barefoot_manu_flood_regime.tex #flood_regime_experiments.bib
	pandoc --reference-doc=pandoc_reference_doc.docx -o barefoot_manu_flood_regime.docx barefoot_manu_flood_regime.tex

fix_biber:
	rm -rf `biber --cache`

backup:
	zip -9y versions/version-`date +"%Y%m%d"`.zip *
