CXXFLAGS=-std=c++11 -g

secant_method_test: secant_method_test.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

tags:
	ctags *.cpp *.h

clean:
	$(RM) *.o
	$(RM) secant_method_test
