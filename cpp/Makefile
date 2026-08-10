CXX = g++
ZSTD_PREFIX ?= $(shell brew --prefix zstd 2>/dev/null || pkg-config --variable=prefix zstd 2>/dev/null || echo /usr/local)
CXXFLAGS = -std=c++17 -Wall -Wextra -Wno-deprecated-declarations -O2 -Iinclude -Iexternal -I$(ZSTD_PREFIX)/include
DEBUG_CXXFLAGS = -std=c++17 -Wall -Wextra -Wno-deprecated-declarations -O0 -g3 -fno-inline -Iinclude -Iexternal -I$(ZSTD_PREFIX)/include
LDFLAGS = -L$(ZSTD_PREFIX)/lib -lzstd

SRC_DIR = src
OBJ_DIR = obj
DEBUG_OBJ_DIR = obj_debug
DEV_OBJ_DIR = obj_dev
BIN = bin/ldzip
DEBUG_BIN = bin/ldzip_debug
DEV_BIN = bin/ldzip-dev
.PHONY: build debug dev clean test

SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SOURCES))
DEBUG_OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(DEBUG_OBJ_DIR)/%.o,$(SOURCES))
DEV_OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(DEV_OBJ_DIR)/%.o,$(SOURCES))

build: $(BIN)
	@echo "✅ Build complete."
	@echo "👉 To use this binary, either:"
	@echo "   1) Add '$(abspath bin)' to your \$$PATH, e.g."
	@echo "        export PATH=\"$(abspath bin):\$$PATH\""
	@echo "   2) Create a symlink in a directory already on your PATH, e.g."
	@echo "        ln -s $(abspath $(BIN)) /usr/local/bin/"

$(BIN): $(OBJECTS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

debug: $(DEBUG_BIN)

$(DEBUG_BIN): $(DEBUG_OBJECTS)
	@mkdir -p bin
	$(CXX) $(DEBUG_CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(DEBUG_OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(DEBUG_OBJ_DIR)
	$(CXX) $(DEBUG_CXXFLAGS) -c $< -o $@

dev: $(DEV_BIN)
	@echo "✅ Development build complete: $(DEV_BIN)"
	@echo "👉 Use './$(DEV_BIN)' for testing new features"

$(DEV_BIN): $(DEV_OBJECTS)
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(DEV_OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(DEV_OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(DEBUG_OBJ_DIR) $(DEV_OBJ_DIR) bin/

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
	@cd tests && ./test_tabular_concat_non_overlap.sh && ./test_square_concat_non_overlap.sh && ./test_tabular_concat_overlap.sh && ./test_concat_sequential.sh

test-filter:
	@rm -rf tests/output/*
	@cd tests && ./test_filter_square_unit_lossless.sh && ./test_filter_tabular_plink_lossless.sh && ./test_filter_tabular_plink_combination.sh

test-find-tag-variants:
	@rm -rf tests/output/*
	@cd tests && ./test_find_tag_variants.sh

test-ld-pruning:
	@rm -rf tests/output/*
	@cd tests && ./test_ld_pruning.sh

test: create-dir test-square test-tabular test-concat test-filter test-find-tag-variants test-ld-pruning

create-dir:
	@cd tests && mkdir -p output