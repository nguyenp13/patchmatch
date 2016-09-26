CXX=g++
CXXFLAGS= -fmax-errors=3 -fopenmp

all: CXXFLAGS += -O3 -funroll-loops -DBUILD_DEBUG=0 -std=c++0x 
debug: CXXFLAGS += -g -DBUILD_DEBUG=1 -std=c++11
all: LDFLAGS = -lpng -fopenmp
debug: LDFLAGS = -lpng -fopenmp

SRCS=main.cpp
OBJS=$(SRCS:.cpp=.o)
EXE=main

all: $(SRCS) $(EXE)
debug: $(SRCS) $(EXE)

$(EXE): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS)

-include $(OBJS:.o=.d)

.cpp.o: 
	$(CXX) -c $(CXXFLAGS) $< -o $@
	$(CXX) -M $(CXXFLAGS) $< > $*.d

clean:
	rm -f $(OBJS) $(EXE) $(SRCS:.cpp=.d)

.PHONY: clean all debug
