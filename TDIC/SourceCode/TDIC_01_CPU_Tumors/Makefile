OBJS = GTMatrix.o PanCanGTMatrix.o PanCanTDIC.o TDIC.o TDIMatrix.o TDIC_Main.o 
DEPS = GTMatrix.h PanCanGTMatrix.h PanCanTDIC.h TDIC.h TDIMatrix.h

%.o: %.cpp $(DEPS)
	g++  -O3 -fopenmp -c $< -o $@

PanCanTDICexeOOMP: $(OBJS)
	g++  -O3 -fopenmp $(OBJS) -o $@

clean:
	rm -f *.o PanCanTDICexeOOMP


	