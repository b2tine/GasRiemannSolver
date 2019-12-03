CXXFLAGS=-std=c++11 -g
INCLUDE = -I../Common


all: riemann_solver_test


#dambreak_test: dambreak_test.cpp riemann_problem.o
#	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@

#borewave_test: borewave_test.cpp riemann_problem.o
#	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@

riemann_solver_test: riemann_solver_test.cpp riemann_problem.o ../Common/util.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@

#shallow1d: shallow1d.cpp riemann_problem.o ../Common/util.o
#	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@



riemann_problem.o: riemann_problem.cpp riemann_problem.h ../Common/secant_method.h
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

../Common/util.o: ../Common/util.h


tags:
	ctags *.cpp *.h ../Common/*.cpp ../Common/*.h

clean:
	$(RM) *.o
	$(RM) dambreak_test
	$(RM) borewave_test
	$(RM) riemann_solver_test
	$(RM) shallow1d
