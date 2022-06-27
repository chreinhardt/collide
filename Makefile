
CXXFLAGS  = -I../tipsy_wrapper
LDFLAGS = -L../tipsy_wrapper
LIBS    = -ltipsy

# tirpc library (needed if glibc >= 2.32)
RPC_LIB = -ltirpc

LIBS += $(RPC_LIB)

collide: collide.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f collide *.o
