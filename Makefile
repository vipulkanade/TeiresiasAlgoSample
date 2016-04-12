#Teiresias make file

teiresias: alphabet.o config.o convolution.o pattern.o seqs.o main.o
	g++ -o teiresias alphabet.o config.o convolution.o pattern.o seqs.o main.o

alphabet.o: alphabet.cpp alphabet.h
	g++ -c alphabet.cpp

config.o: config.cpp config.h
	g++ -c config.cpp

convolution.o: convolution.cpp convolution.h
	g++ -c convolution.cpp

pattern.o: pattern.cpp pattern.h
	g++ -c pattern.cpp

seqs.o: seqs.cpp seqs.h
	g++ -c seqs.cpp

main.o: main.cpp
	g++ -c main.cpp

clean:
	rm teiresias alphabet.o config.o convolution.o pattern.o seqs.o main.o
