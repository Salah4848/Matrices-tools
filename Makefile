CC = g++
CFLAGS = -std=c++11 -Wall

TARGET = myprogram
SRCS = main.cpp overloads.cpp
OBJS = $(SRCS:.cpp=.o)
DEPS = Matrix.h

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

