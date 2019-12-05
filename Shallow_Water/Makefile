CXXFLAGS=-std=c++11 -g
INCLUDE = -I../Common


all: sw_riemann_solver_test


#dambreak_test: dambreak_test.cpp sw_riemann_problem.o
#	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@

#borewave_test: borewave_test.cpp sw_riemann_problem.o
#	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@

sw_riemann_solver_test: sw_riemann_solver_test.cpp sw_riemann_problem.o ../Common/util.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@

#shallow_water1d: shallow_water1d.cpp riemann_problem.o ../Common/util.o
#	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@



sw_riemann_problem.o: sw_riemann_problem.cpp sw_riemann_problem.h ../Common/secant_method.h
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

../Common/util.o: ../Common/util.h


tags:
	ctags *.cpp *.h ../Common/*.cpp ../Common/*.h

clean:
	$(RM) *.o
	$(RM) dambreak_test
	$(RM) borewave_test
	$(RM) sw_riemann_solver_test
	$(RM) shallow_water1d
