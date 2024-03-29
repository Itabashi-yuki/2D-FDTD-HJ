OBJS = fdtd2d.o current_source.o update_E.o update_E_PML.o update_H.o update_H_PML.o update_J.o output_Ez.o initialize_PML.o \
		 initialize_Plasma.o output_pal.o cal_obs_n0.o make_dir.o judge_Ez_div.o allocate.o output_obs.o
HEADERS = fdtd2d.h

OPTS = -O3 -Wall -std=gnu++17

all: main

main: $(OBJS)
	g++ -o $@ $(OBJS)

%.o: %.cpp $(HEADERS)
	g++ -c $< $(OPTS)

clean:
	rm -rf main *.o

.PHONY: all clean