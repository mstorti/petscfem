CPPFLAGS=-g -frepo 
LDFLAGS=-lg++

user: user.o stack.o
	g++ -o user $^

clean:
	rm -f *.o *~ user

stack.o user.o: stack.h

