obj = main.o func.o dc.o

exec: ${obj}
	g++ -o exec ${obj}  -O3

main.o : func.h
	g++ -c main.cpp  -O3
func.o : 
	g++ -c func.cpp  -O3
dc.o   : func.h
	g++ -c dc.cpp    -O3

.PHONY : clean
clean:
	rm exec ${obj}
