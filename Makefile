# simple Makefile

.PHONY: clean realclean

# optimization flags, can be changed
OPT_FLAGS = -O3
INC=-I/opt/local/include
##################################################################################
## no need to edit anything beyond this point, unless you know what you are doing 
##

## sources
HDRS=   $(wildcard *.h)
OBJS=   $(patsubst %.C,%.o, $(wildcard *.C)) $(patsubst %.c,%.o, $(wildcard *.c))

# flags
CFLAGS=-Wall ${OPT_FLAGS} ${INC}
CXXFLAGS=-Wall -std=c++0x ${OPT_FLAGS} ${INC}
LINK_FLAGS=-L/opt/local/lib -lnetcdf

TARGET = Pre

$(TARGET): $(OBJS)
	$(CXX) $^ -o ${TARGET} ${OPT_FLAGS} ${LINK_FLAGS}

# dependencies, blanket
$(OBJS):$(HDRS)

## other options
clean:
	rm -rf $(OBJS)

realclean:
	rm -rf $(OBJS) $(TARGET)

## end

