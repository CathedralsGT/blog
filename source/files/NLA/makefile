BIB = /home/aaa/Desktop/math-soft/Mandelbrot/doc

all : iterative_methods.tex $(BIB)
	xelatex iterative_methods.tex
	bibtex iterative_methods
	xelatex iterative_methods.tex
	xelatex iterative_methods.tex

clean :
	rm -f *.aux *.bbl *.blg *.log *.out 

.PHONY : all clean
