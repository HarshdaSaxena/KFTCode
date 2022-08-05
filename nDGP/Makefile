CPP	:= g++
CCFLAGS	:= -O -g -Wall -I/usr/local/include -pedantic -std=c++11
LDFLAGS	:= /usr/local/lib
LIBS	:= -lKFT -lastrocpp -lgsl -lgslcblas

files	:= \
	main.cpp \
        forceTerm.cpp \
	gScreening.cpp
objects	:= $(files:.cpp=.o)

%.o: %.cpp
	$(CPP) $(CCFLAGS) -c $<

default: main.x

main.x: $(objects)
	$(CPP) $(CCFLAGS) -o $@ $(objects) -L$(LDFLAGS) $(LIBS)

clean:
	$(RM) $(objects)
distclean: clean
	$(RM) main.x
