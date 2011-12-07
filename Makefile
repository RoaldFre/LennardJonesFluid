OBJECTS = system.o main.o render.o measure.o
WARNINGS = -pedantic -Wextra -Wall -Wwrite-strings -Wshadow -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes
PROFILE = 
OPTIM = -O4 -flto -DNDEBUG -fexcess-precision=fast
CFLAGS = $(WARNINGS) -std=c99 -march=native -ggdb $(OPTIM) $(PROFILE) $(BROWN)

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
	@gcc -o main $(OBJECTS) -lm -lSDL -lGL -ggdb $(PROFILE)

.c.o:
	@echo "CC	$@"
	@gcc -o $@ $(CFLAGS) -c $<

clean:
	rm -f *.o

