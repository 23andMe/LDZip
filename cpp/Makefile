CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -Iinclude -Iexternal
DEBUG_CXXFLAGS = -std=c++17 -Wall -Wextra -O0 -g3 -fno-inline -Iinclude -Iexternal

SRC_DIR = src
OBJ_DIR = obj
DEBUG_OBJ_DIR = obj_debug
BIN = bin/ldzip
DEBUG_BIN = bin/ldzip_debug
.PHONY: build debug clean test

SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
DEBUG_OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(DEBUG_OBJ_DIR)/%.o,$(SOURCES))

build: $(BIN)
	@echo "✅ Build complete."
	@echo "👉 To use this binary, either:"
	@echo "   1) Add '$(abspath bin)' to your \$$PATH, e.g."
	@echo "        export PATH=\"$(abspath bin):\$$PATH\""
	@echo "   2) Create a symlink in a directory already on your PATH, e.g."
	@echo "        ln -s $(abspath $(BIN)) /usr/local/bin/"

$(BIN): $(OBJECTS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $^


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

debug: $(DEBUG_BIN)

$(DEBUG_BIN): $(DEBUG_OBJECTS)
	@mkdir -p bin
	$(CXX) $(DEBUG_CXXFLAGS) -o $@ $^

$(DEBUG_OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(DEBUG_OBJ_DIR)
	$(CXX) $(DEBUG_CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(DEBUG_OBJ_DIR) bin/

quick-test: create-dir test-fast

test-fast:
	@rm -rf tests/output/*
	@cd tests && ./test_square_unit_lossless.sh

test-square:  
	@rm -rf tests/output/*
	@cd tests && ./test_square_unit_lossless.sh && ./test_square_plink_lossless.sh && ./test_square_combination.sh

test-tabular:
	@rm -rf tests/output/*
	@cd tests && ./test_tabular_plink_lossless.sh && ./test_tabular_plink_combination.sh

test-concat:
	@rm -rf tests/output/*
	@cd tests && ./test_tabular_concat.sh && ./test_square_concat.sh && ./test_tabular_concat_overlap.sh && ./test_concat_sequential.sh && ./test_concat_fail.sh

test-overflow:
	@rm -rf tests/output/*
	@cd tests && ./test_overflow_tabular.sh && ./test_overflow_concat.sh && ./test_overflow_concat_overlap.sh

test-filter:
	@rm -rf tests/output/*
	@cd tests && ./test_filter_square_unit_lossless.sh && ./test_filter_tabular_plink_lossless.sh && ./test_filter_tabular_plink_combination.sh 


test: create-dir test-square test-tabular test-concat test-overflow test-filter

create-dir:
	@cd tests && mkdir -p output