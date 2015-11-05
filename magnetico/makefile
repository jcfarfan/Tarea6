
trayectoria_1.000000_30.000000.pdf:trayectoria_1.000000_30.000000.dat
	python graficar_v.py trayectoria_1.000000_30.000000.dat
trayectoria_1.000000_30.000000.dat: solver_v.x
	./solver_v.x 1.0 30.0
./solver_v.x: solver_v.c
	cc solver_v.c -lm -o solver_v.x

clean:
	rm *.dat *.pdf
