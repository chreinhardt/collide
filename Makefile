
CXXFLAGS  = -I../tipsy_wrapper
LDFLAGS = -L../tipsy_wrapper
LIBS    = -ltipsy

collide: collide.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f collide *.o
