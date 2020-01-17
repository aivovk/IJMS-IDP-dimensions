# \todo check which uncommented flags are not needed
# \todo write non MKL version
# \todo write cuda/GPU version

CC=g++
DEBUG = -g -O0 -fno-omit-frame-pointer #-pg
OPTIMIZE = -O3

CFLAGS = -Wall $(OPTIMIZE) #-Wl,--no-as-needed

GSL=-lgsl -lm #-lgslcblas

MKL=-L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core #-lpthread -ldl

# for visualization
SFML=
# moved to SFML below
#FL = -lGL -lGL #-lrt

#SRCS := $(filter-out src/main.cpp, $(wildcard src/*.cpp))
#HDRS := $(wildcard src/*.h)
OBJS := $(filter-out obj/main.o, $(patsubst src/%.cpp,obj/%.o,$(wildcard src/*.cpp)))

# \todo does not work, can compile one or the other, not both
all: nup_model_temp nup_model_sfml

# compile without SFML
nup_model_temp: obj/main.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ obj/main.o $(OBJS) $(SFML) $(GSL) $(MKL)

# SFML for visualization
nup_model_sfml: SFML=-lsfml-graphics -lsfml-window -lsfml-system -D SFML -lGL -lGLU -D GL
nup_model_sfml: obj/main.o.sfml $(OBJS)
	$(CC) $(CFLAGS) -o $@ obj/main.o.sfml $(OBJS) $(SFML) $(GSL) $(MKL)

obj/%.o : src/%.cpp src/%.h
	$(CC) $(CFLAGS) $(SFML) -c $< -o $@

obj/main.o: src/main.cpp src/*.h
	$(CC) $(CFLAGS) $(SFML) -c src/main.cpp -o $@

obj/main.o.sfml:
	$(CC) $(CFLAGS) $(SFML) -c src/main.cpp -o $@

clean:
	rm *~ obj/*.o obj/main.o.sfml src/*~ nup_model_temp nup_model_sfml


