CC       = g++
CXX      = g++
CXXFLAGS = -fopenmp -Wall -O3 -pedantic -g
CCFLAGS = -fopenmp -Wall -O3 -pedantic -g
LDLIBS   = -lblas -lboost_program_options -fopenmp
TARGET   = main
TARGETNEXT  = ShallowWater
HDRS     = ShallowWater.h

default: $(TARGET)

$(TARGET): $(TARGET).o $(TARGETNEXT).o

$(TARGET).o: $(TARGET).cpp 

$(TARGETNEXT).o: $(TARGETNEXT).cpp $(HDRS)

.PHONY: test1 test2 test3 test4 clean git plot

test1: $(TARGET)
	time ./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 1 --option1 1

test2: $(TARGET)
	time ./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 2 --option1 1

test3: $(TARGET)
	time ./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3 --option1 1

test4: $(TARGET)
	time ./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 4 --option1 1

test5: $(TARGET)
	time ./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 1 --option1 2

test6: $(TARGET)
	time ./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 2 --option1 2

test7: $(TARGET)
	time ./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3 --option1 2

test8: $(TARGET)
	time ./$(TARGET) --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 4 --option1 2

clean:
	rm -f $(TARGET) *.o
	rm -f $(TARGETNEXT) *.o
	rm -f $(TARGET)
	rm -f repository.log

git:
	git log --name-status > repository.log