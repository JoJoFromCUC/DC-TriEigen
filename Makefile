obj = main.o func.o dc.o

exec: ${obj}
	g++ -o exec ${obj}  -static -DNDEBUG -O3

main.o : func.h
	g++ -c main.cpp  -static -DNDEBUG -O3
func.o : 
	g++ -c func.cpp  -static -DNDEBUG -O3
dc.o   : func.h
	g++ -c dc.cpp  -static -DNDEBUG   -O3

.PHONY : clean
clean:
	rm exec ${obj}
