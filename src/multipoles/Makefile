CFLAGS += -I../../include
CXXFLAGS := $(CFLAGS)
LDFLAGS += -L../../lib -lester

LD := $(CXX)

SRCS := multipoles.cpp
TARGETS := multipoles

.PHONY: all clean depend

all: $(TARGETS)

multipoles: multipoles.o
	$(LD) -o $@ $^ $(LDFLAGS)

clean:
	rm -f *.o *~ core $(TARGETS)


###
# makedepends related settings
# -Y		Don't search standard library paths, this is excessive
# -f		output dependencies to this file
# -s		use this token string
###

DEPFILE 	= .depends
DEPTOKEN	= '\# MAKEDEPENDS'
DEPFLAGS	= -Y -f $(DEPFILE) -s $(DEPTOKEN)

depend:
	rm -f $(DEPFILE)
	make $(DEPFILE)

$(DEPFILE):
	@echo $(DEPTOKEN) > $(DEPFILE)
	makedepend $(DEPFLAGS) -- $(CFLAGS) -- $(SRCS) >&/dev/null

# put this file in the last line of your Makefile
sinclude $(DEPFILE)
# DO NOT DELETE
