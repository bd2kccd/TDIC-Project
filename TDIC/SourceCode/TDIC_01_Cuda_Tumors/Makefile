OBJS = GTMatrix.o PanCanGTMatrix.o PanCanTDIC.o PanCanTDIC_GeLoop_Cuda.o TDIC.o TDIC_GeLoop_Cuda.o TDIMatrix.o TDIC_Main.o 
DEPS = GTMatrix.h PanCanGTMatrix.h PanCanTDIC.h PanCanTDIC_GeLoop_Cuda.h TDIC.h TDIC_GeLoop_Cuda.h TDIMatrix.h

%.o: %.cpp $(DEPS)
	nvcc -O3 -x cu -arch=sm_20 -I. -dc $< -o $@

PanCanTDICexeGPU: $(OBJS)
	nvcc -arch=sm_20 $(OBJS) -o $@

clean:
	rm -f *.o PanCanTDICexeGPU



	