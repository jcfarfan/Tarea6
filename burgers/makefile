all: burger.gif matrizdata500.dat graficar.py ./solver_burger.x

burger.gif: graficar.py
	python graficar.py

graficar.py: matrizdata500.dat
	python graficar.py

matrizdata500.dat: solver_burger.x
	./solver_burger.x

./solver_burger.x: solver_burger.c
	cc solver_burger.c -lm -o solver_burger.x
clean:
	rm -f *.dat *.gif *.png
