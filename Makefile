CXX = g++
CXXFLAGS  += -O2 -Wall -Werror -Wextra
RCXXFLAGS := `root-config --cflags --libs`
LDFLAGS   += -lconf -L./git/config/lib
RLDFLAGS  := -lEG

BINDIR = ./bin
BLDDIR = ./build
SRCDIR = ./src

RPPSRCS = $(wildcard $(SRCDIR)/*.C)
RPPEXES = $(patsubst $(SRCDIR)/%.C,$(BINDIR)/%,$(RPPSRCS))
RPPDEPS = $(patsubst $(SRCDIR)/%.C,$(BLDDIR)/%.d,$(RPPSRCS))

CPPSRCS = $(wildcard $(SRCDIR)/*.cpp)
CPPEXES = $(patsubst $(SRCDIR)/%.cpp,$(BINDIR)/%,$(CPPSRCS))
CPPDEPS = $(patsubst $(SRCDIR)/%.cpp,$(BLDDIR)/%.d,$(CPPSRCS))

EXES = $(RPPEXES) $(CPPEXES)
DEPS = $(RPPDEPS) $(CPPDEPS)

.PHONY: all clean

all: $(RPPEXES) $(CPPEXES)

$(BINDIR)/% : $(SRCDIR)/%.C
	@mkdir -p $(BINDIR)
	@mkdir -p $(BLDDIR)/$(@D)
	$(CXX) $(CXXFLAGS) $(RCXXFLAGS) -MMD -MF $(BLDDIR)/$(@D)/$(*F).d $< -o $@ \
		$(LDFLAGS) $(RLDFLAGS)

$(BINDIR)/% : $(SRCDIR)/%.cpp
	@mkdir -p $(BINDIR)
	@mkdir -p $(BLDDIR)/$(@D)
	$(CXX) $(CXXFLAGS) -MMD -MF $(BLDDIR)/$(@D)/$(*F).d $< -o $@ \
		$(LDFLAGS)

clean:
	@$(RM) $(EXES) $(DEPS)
	@rm -f $(BINDIR)/*
	@rm -rf $(BLDDIR)/*

-include $(DEPS)
