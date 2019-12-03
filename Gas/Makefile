CXXFLAGS=-std=c++11 -g


all: withdraw_test compress_test \
	secant_method_test riemann_solver_test \
	gas1d


withdraw_test: withdraw_test.cpp riemann_problem.o
	$(CXX) $(CXXFLAGS) $^ -o $@

compress_test: compress_test.cpp riemann_problem.o
	$(CXX) $(CXXFLAGS) $^ -o $@

secant_method_test: secant_method_test.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

riemann_solver_test: riemann_solver_test.cpp riemann_problem.o util.o
	$(CXX) $(CXXFLAGS) $^ -o $@

gas1d: gas1d.cpp riemann_problem.o util.o
	$(CXX) $(CXXFLAGS) $^ -o $@



riemann_problem.o: riemann_problem.h secant_method.h

util.o: util.h


tags:
	ctags *.cpp *.h

clean:
	$(RM) *.o
	$(RM) withdraw_test
	$(RM) compress_test
	$(RM) secant_method_test
	$(RM) riemann_solver_test
	$(RM) gas1d
