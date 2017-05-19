
CXX       = g++ -pipe -std=c++14
CXX_FLAGS = -mtune=native -march=native -m64 -O3 -fPIC -fopenmp
CXX_INC   = -I/usr/local/include -I.


all: harmonic


harmonic:
	$(CXX) -o harmonic.exe harmonic.cpp -g $(CXX_FLAGS) $(CXX_INC)

roots:
	$(CXX) -o bisection.exe root_finding/bisection.cpp -g $(CXX_FLAGS) $(CXX_INC)
	$(CXX) -o brents.exe root_finding/brents.cpp -g $(CXX_FLAGS) $(CXX_INC)
	$(CXX) -o gsl1.exe root_finding/gsl1.cpp -g $(CXX_FLAGS) $(CXX_INC) -I./root_finding/ -lgsl -lopenblas
	$(CXX) -o gsl2.exe root_finding/gsl2.cpp -g $(CXX_FLAGS) $(CXX_INC) -I./root_finding/ -lgsl -lopenblas
	#$(CXX) -o roots_boost.exe root_finding.cpp -g $(CXX_FLAGS) $(CXX_INC) -lrt


mini:
	$(CXX) -o minimization/gsl_minimization.exe  minimization/gsl_minimization.c -g $(CXX_FLAGS) $(CXX_INC) -lgsl -lopenblas

lab6:
	$(CXX) -o lab6.exe lab6.cpp $(CXX_FLAGS) $(CXX_INC) -lgsl -lopenblas #-DDEBUG

lab61:
	$(CXX) -o lab6_1.exe lab6_1.cpp $(CXX_FLAGS) $(CXX_INC) -lgsl -lopenblas #-DDEBUG

clean:
	@rm *.exe