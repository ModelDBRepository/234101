CPPFLAGS = -std=c++11
BOOSTINC = -I/usr/include 
LINKFLAGS = -lnetcdf_c++4


all: Network.o main.o icpg


icpg: Network.o main.o 
	g++ Network.o main.o $(BOOSTINC) $(LINKFLAGS) -o icpg  $(CPPFLAGS)

Network.o: ./src/Network.cpp ./src/Network.h
	g++ -c ./src/Network.cpp $(BOOSTINC) -o Network.o $(CPPFLAGS) -O3

main.o: ./src/main.cpp ./src/Network.h
	g++ -c ./src/main.cpp $(BOOSTINC) -o main.o $(CPPFLAGS) -O3


clean:
	rm *.o
	rm icpg
