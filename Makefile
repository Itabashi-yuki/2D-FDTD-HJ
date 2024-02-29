OBJS = fdtd2d.o current_source.o update_E.o update_E_PML.o update_H.o update_H_PML.o output_Ez.o initialize_PML.o
HEADERS = fdtd2d.h

OPTS = -O3 -Wall

all: main

main: $(OBJS)
	g++ -o $@ $(OBJS)

%.o: %.cpp $(HEADERS)
	g++ -c $< $(OPTS)

clean:
	rm -rf main *.o

.PHONY: all clean