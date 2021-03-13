PROCESS_NUM = 2
CC = g++
LINT_FLAGS = -Wall -Wextra -pedantic -Werror -Wfatal-errors
MY_PYTHON_PATH = python
EXECUTION_PREFIX =
ifeq ($(MODE),debug)
	OPTIMIZATION_FLAGS := -O0 -g -pg
	LFLAGS := -pg -lclipper
else ifeq ($(MODE),parallel)
	CC := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI -O3
	EXECUTION_PREFIX := mpirun -n $(PROCESS_NUM)
	LINT_FLAGS = -Wall -Wextra -pedantic -Wfatal-errors 
	LFLAGS = -lclipper
else ifeq ($(MODE),debug_parallel)
	CC := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI -O0 -g -pg -frecord-gcc-switches
	LFLAGS := -pg
	EXECUTION_PREFIX := mpirun -n $(PROCESS_NUM)
	LINT_FLAGS :=
	LFLAGS = -lboost_mpi -lboost_serialization
else ifeq ($(MODE),intel)
	CC := icpc
	OPTIMIZATION_FLAGS := -O3 -ipo -xHost -fp-model precise
	LINT_FLAGS :=
else ifeq ($(MODE),clang)
	CC := clang++
	OPTIMIZATION_FLAGS := -Weverything -Werror -ferror-limit=1 -Wno-error=padded
	LINT_FLAGS :=
else
	MODE = production
	OPTIMIZATION_FLAGS := -O3 -march=native
	LFLAGS = -lclipper
endif
CFLAGS = $(LINT_FLAGS) $(OPTIMIZATION_FLAGS)
SOURCE_LIST := $(shell find . -name 'test.cpp')
EMPTY :=
TEMP1 := $(subst /test.cpp,$(EMPTY),$(SOURCE_LIST))
TEMP2 := $(subst ./,/,$(TEMP1))
FOLDER_LIST := $(addprefix $(MODE),$(TEMP2))
BINARY_LIST := $(addsuffix /test.exe,$(FOLDER_LIST))
OBJECT_LIST := $(addsuffix /test.o,$(FOLDER_LIST))
RESULT_LIST := $(addsuffix /test_passed.res,$(FOLDER_LIST))
NORMAL_TERMINATION_LIST := $(addsuffix /normal_termination.res,$(FOLDER_LIST))
TEST_ROOT := $(shell pwd)
SOURCE_DEPENDENCIES = $(shell find ../source/ -name '*.*pp')

.PRECIOUS: $(NORMAL_TERMINATION_LIST) $(RESULT_LIST)

all: $(RESULT_LIST)
	@echo "done"

$(RESULT_LIST): $(MODE)/%/test_passed.res: $(MODE)/%/normal_termination.res %/test.py
	cd $(MODE)/$* && $(MY_PYTHON_PATH) $(TEST_ROOT)/$*/test.py && tail ./test_passed.res

$(NORMAL_TERMINATION_LIST): $(MODE)/%/normal_termination.res: $(MODE)/%/test.exe
	cd $(MODE)/$* && $(EXECUTION_PREFIX) ./test.exe && touch normal_termination.res

$(BINARY_LIST): $(MODE)/%/test.exe: $(MODE)/%/test.o
	$(MAKE) -C .. MODE=$(MODE)
	$(CC) $(MODE)/$*/test.o ../library_$(MODE)/librich.a -L $(HDF5_LIB_PATH) -lhdf5 -lhdf5_cpp -lz $(LFLAGS) -o $(MODE)/$*/test.exe

$(OBJECT_LIST): $(MODE)/%/test.o: %/test.cpp $(SOURCE_DEPENDENCIES) | $(MODE)/%
	$(CC) -c $(CFLAGS) $*/test.cpp -o $(MODE)/$*/test.o -I ..
	$(CC) -MM $*/test.cpp -o $(MODE)/$*/test.d -I ..
	cp $*/*.txt $(MODE)/$* 2> /dev/null || :

$(FOLDER_LIST): %:
	mkdir -p $*

test:
	@echo $(FOLDER_LIST)


