obj = main.o func.o dc.o

exec: ${obj}
	g++ -o exec ${obj} -fcilkplus -O3

main.o : func.h
	g++ -c main.cpp -fcilkplus -O3
func.o : 
	g++ -c func.cpp -fcilkplus -O3
dc.o   : func.h
	g++ -c dc.cpp   -fcilkplus -O3

.PHONY : clean
clean:
	rm exec ${obj}
