
SDIR = ../../src/cpp/
IDIR = $(SDIR)
CXX = g++
#OPT = -O3
OPT = -g

# Common objects
OBJ = $(SDIR)creep_law.o $(SDIR)stiffness_matrix.o $(SDIR)fault.o $(SDIR)earthquake.o $(SDIR)temperature_profile.o $(SDIR)strength_profile.o $(SDIR)fault_cooling.o $(SDIR)romberg_integration.o $(SDIR)activ_energy.o

all: test_model

test_model: test_model.cpp $(OBJ) $(SDIR)side_creep.o input_parameters.h
	$(CXX) $(filter %.cpp, $^) $(filter %.o, $^) -fopenmp -o $@ $(OPT) -I$(IDIR)


%.o:%.cpp %.h
	$(CXX) $< -c $(OPT) -I$(IDIR) -fopenmp -o $@

clear: clean
	rm -f test_model

clean:
	rm -f $(SDIR)*.o