CXX = g++
CXXFLAGS = -std=c++17 -g

INCLUDES = include
SOURCE = src
OUT = out
MAIN = main.cc

EXENAME = raytracer.exe

FILES = $(patsubst $(SOURCE)/%.cc,$(OUT)/%.o,$(wildcard $(SOURCE)/*.cc))

all: $(EXENAME)

$(OUT)/%.o: $(SOURCE)/%.cc $(INCLUDES)/%.h
	mkdir -p $(OUT)
	$(CXX) $(CXXFLAGS) -I$(INCLUDES) -c $< -o $@

$(EXENAME): $(MAIN) $(FILES)
	$(CXX) $(CXXFLAGS) -I$(INCLUDES) $< $(FILES) -o $@

clean:
	rm -rf $(EXENAME)
	rm -rf $(OUT)
	rm -rf *.o
