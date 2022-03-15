run_tests: matching_tests biclique_cover_tests run_matching_tests run_biclique_cover_tests
	
functions.o: functions.cpp functions.h
	g++ -std=c++14 -g -c functions.cpp
matching_tests.o: matching_tests.cpp 
	g++ -std=c++14 -g -c matching_tests.cpp 
biclique_cover_tests.o: biclique_cover_tests.cpp
	g++ -std=c++14 -g -c biclique_cover_tests.cpp
matching_tests: functions.o matching_tests.o
	g++ -g functions.o matching_tests.o -lm -o matching_tests
biclique_cover_tests: functions.o biclique_cover_tests.o
	g++ -g functions.o biclique_cover_tests.o -lm -o biclique_cover_tests
clear:
	rm -f biclique_cover_tests
	rm -f matching_tests
	rm -f *.o
run_matching_tests: matching_tests
	./matching_tests
run_biclique_cover_tests:  
	./biclique_cover_tests

