OBJECTS = system.o main.o render.o measure.o
WARNINGS = -pedantic -Wextra -Wall -Wwrite-strings -Wshadow -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes
PROFILE = 
OPTIM = -O4 -flto -DNDEBUG -fexcess-precision=fast -ffast-math -finline-limit=2000
CFLAGS = $(WARNINGS) -std=c99 -pipe -march=native -ggdb $(OPTIM) $(PROFILE)

all: main
profile:
	make -B PROFILE="-pg -DNDEBUG"
debug:
	make -B OPTIM="-O0 -DDEBUG"
O0:
	make -B OPTIM="-O0 -DNDEBUG"
O1:
	make -B OPTIM="-O1 -DNDEBUG"
O2:
	make -B OPTIM="-O2 -DNDEBUG"
O3:
	make -B OPTIM="-O3 -DNDEBUG"

main: $(OBJECTS)
	@echo "LD	main"
	@gcc -o main $(CFLAGS) $(OBJECTS) -lm -lSDL -lGL -ggdb $(PROFILE)

.c.o:
	@echo "CC	$@"
	@gcc -o $@ $(CFLAGS) -c $<

tex:
	cd latex; make

plots:
	./runEnergies.sh
	./runPairCorrelation.sh
	./runPressure.sh

clean:
	rm -f *.o
