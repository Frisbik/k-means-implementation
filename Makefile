k_means: k_means.o datapoint.o center_initialization.o
	g++ -Wall -std=c++14 -stdlib=libc++ -g -I boost -o k_means k_means.o datapoint.o center_initialization.o

k_means.o: k_means.cpp datapoint.hpp k_means.hpp center_initialization.hpp
	g++ -Wall -std=c++14 -stdlib=libc++ -g -I boost -c k_means.cpp

center_initialization.o: datapoint.hpp center_initialization.cpp center_initialization.hpp
	g++ -Wall -std=c++14 -stdlib=libc++ -g -I boost -c center_initialization.cpp

datapoint.o: datapoint.cpp datapoint.hpp
	g++ -Wall -std=c++14 -stdlib=libc++ -g -I boost -c datapoint.cpp


clean :
	rm −f ∗˜ *.o
# DO NOT DELETE
