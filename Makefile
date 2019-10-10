CXXFLAGS=-std=c++11


all: driver withdraw_test secantmethod_test


driver: driver.cpp riemann_problem.o
	$(CXX) $(CXXFLAGS) $^ -o $@

withdraw_test: withdraw_test.cpp riemann_problem.o
	$(CXX) $(CXXFLAGS) $^ -o $@

secantmethod_test: secantmethod_test.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@


clean:
	$(RM) *.o
	$(RM) driver
	$(RM) withdraw_test
	$(RM) secantmethod_test
