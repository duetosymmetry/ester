SRC=matrix.cpp mat_math.cpp mat_spec.cpp matrix_block_diag.cpp matrix_map.cpp
SRC_H=matrix.h
DEPS=

LIB=lib$(DIR).a
OBJS=$(SRC:.cpp=.o)
HEADERS=$(addprefix $(ESTER_ROOT)/include/,$(SRC_H))
DEPS:=$(addprefix $(ESTER_ROOT)/include/,$(DEPS))

.PHONY: clean distclean

$(LIB): $(HEADERS) $(OBJS)
	ar rvu $(LIB) $(OBJS)
	ranlib $(LIB)

%.o: %.cpp $(HEADERS) $(DEPS)
	$(CC) -c $< $(CCFLAGS) $(VARS)

$(ESTER_ROOT)/include/%.h: %.h $(SRC_H) $(DEPS)
	cp $< $(ESTER_ROOT)/include

clean:
	rm -f *.o
	rm -f $(LIB)

distclean:
	rm -f $(HEADERS)
