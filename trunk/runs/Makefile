CC := g++
DUMP_FOLDER := temp_$(TARGET)
ifeq ($(MODE),debug)
	OPTIMIZATION_FLAGS := -O0 -g -pg
	LFLAGS := -pg
else ifeq ($(MODE),intel)
	CC := icpc
	OPTIMIZATION_FLAGS := -O3 -ipo -xHost -fp-model precise
	LINT_FLAGS := 
	ARCHIVER_FUNC := xiar
else ifeq ($(MODE),parallel)
	CC := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI
	LINT_FLAGS = -Wextra -pedantic -Wfatal-errors -Weffc++ -Wshadow -Wmissing-declarations
else
	MODE := production
	OPTIMIZATION_FLAGS := -O2
	LFLAGS :=
endif
LINT_FLAGS := -Wall -Wextra -pedantic -Werror -Wfatal-errors
CFLAGS := $(LINT_FLAGS) $(OPTIMIZATION_FLAGS)

$*DUMP_FOLDER/normal_termination.res: $(DUMP_FOLDER)/test.exe
	$(MAKE) -C .. set_environ_vars.sh && \
	`cat ../set_environ_vars.sh` && \
	cd $(DUMP_FOLDER) && ./test.exe

$(DUMP_FOLDER)/test.exe: $(DUMP_FOLDER)/test.o | ../library_$(MODE)/librich.a
	$(MAKE) -C .. set_environ_vars.sh
	`cat ../set_environ_vars.sh` && \
	$(CC) $^ $| -o $@ -lhdf5 -lhdf5_cpp -L $$HDF5_LIB_PATH $$LFLAGS

../library_$(MODE)/librich.a:
	$(MAKE) -C .. set_environ_vars.sh
	`cat ../set_environ_vars.sh` && \
	$(MAKE) -C .. MODE=$(MODE)

$(DUMP_FOLDER)/test.o: $(TARGET)/test.cpp
	mkdir -p $(DUMP_FOLDER)
	$(MAKE) -C .. set_environ_vars.sh
	`cat ../set_environ_vars.sh` && \
	$(CC) $(CFLAGS) -c $< -o $@ -I ..
