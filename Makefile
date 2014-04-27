
all: obj/sequencer.o obj/filters.o
	g++ -o sequencer124 obj/sequencer.o obj/filters.o

obj/sequencer.o: src/sequencer.cpp
	g++ -c src/sequencer.cpp -o obj/sequencer.o

obj/filters.o: src/filters.cpp
	g++ -c src/filters.cpp -o obj/filters.o


clean:
	rm -f obj/*.o sequencer124
