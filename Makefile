CC = g++
CFLAGS = --std=c++11 -DNDEBU -O3
target = fast-wclq
all : $(target)

$(target) : fastwclq.cpp basic.h upper_bound.h
	$(CC) $(CFLAGS) fastwclq.cpp -o $(target)

clean :
	-rm $(target)
