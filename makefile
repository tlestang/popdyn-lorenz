popdyn_lorenz: randNormal.o lorenz.o lorenz_transient.o popdyn.o
	g++ -o popdyn_lorenz randNormal.o lorenz.o lorenz_transient.o popdyn.o
popdyn.o: popdyn.cpp lorenz.h
	g++ -c popdyn.cpp
lorenz.o: lorenz.cpp randNormal.h
	g++ -c lorenz.cpp
lorenz_transient.o: lorenz_transient.cpp randNormal.h
	g++ -c lorenz_transient.cpp
randNormal.o: randNormal.cpp
	g++ -c randNormal.cpp
