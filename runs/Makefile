CC := g++
DUMP_FOLDER := temp_$(TARGET)
LINT_FLAGS := -Wall -Wextra -pedantic -Werror -Wfatal-errors
NUM_PROCESSES := 2
PREXECUTE := 
ENDEMIC_SOURCES := $(shell find $(TARGET) -name '*.cpp')
ENDEMIC_OBJECTS := $(patsubst $(TARGET)/%.cpp,$(DUMP_FOLDER)/%.o,$(ENDEMIC_SOURCES))
ifeq ($(MODE),debug)
	OPTIMIZATION_FLAGS := -O0 -g -pg -std=c++11
	LFLAGS := -pg -lclipper -lr3d
else ifeq ($(MODE),intel)
	CC := icpc
	OPTIMIZATION_FLAGS := -O3 -ipo -xHost -fp-model precise -std=c++11
	LINT_FLAGS := 
	LFLAGS := -lclipper -lr3d
	ARCHIVER_FUNC := xiar
else ifeq ($(MODE),parallel_intel)
	PREXECUTE := mpirun -np $(NUM_PROCESSES)
	CC := mpiicpc
	OPTIMIZATION_FLAGS := -O3 -DRICH_MPI -ipo -xHost -fp-model precise  -std=c++11 -qopenmp
	LINT_FLAGS := 
	LFLAGS := -lclipper -lr3d
	ARCHIVER_FUNC := xiar
else ifeq ($(MODE),icpl)
	CC := icpc
	OPTIMIZATION_FLAGS := -O3 -march=core-avx2 -I $(HDF5_DIR)/include -I $(Boost_DIR)/include -fp-model precise  -std=c++11 -qopenmp
	LINT_FLAGS := 
	LFLAGS := -L $(HDF5_DIR)/lib -lclipper -lr3d -qopenmp
else ifeq ($(MODE),icpl_mpi)
	CC = mpic++
	OPTIMIZATION_FLAGS := -O3 -DRICH_MPI -march=core-avx2 -I $(HDF5_DIR)/include -I $(Boost_DIR)/include -fp-model precise  -std=c++11 -qopenmp
	LINT_FLAGS := 
	LFLAGS := -L $(HDF5_DIR)/lib -lclipper -lr3d -qopenmp
else ifeq ($(MODE),nasa_intel)
	PREXECUTE := mpirun -np $(NUM_PROCESSES)
	CC := icpc
	OPTIMIZATION_FLAGS := -O3 -DRICH_MPI -ipo -xAVX -fp-model precise  -std=c++11 -DMPI_NO_CPPBIND
	LINT_FLAGS := 
	LFLAGS := -lclipper -lr3d -lmpi -lmpi++
	ARCHIVER_FUNC := xiar
else ifeq ($(MODE),TACC)
	CCC := icc
	CC := mpicxx
	OPTIMIZATION_FLAGS := -O3 -DRICH_MPI -ipo -xSKYLAKE -fp-model precise  -std=c++11 -DMPI_NO_CPPBIND
	LINT_FLAGS := 
	LFLAGS := -lclipper -lr3d 
	ARCHIVER_FUNC := xiar
else ifeq ($(MODE),parallel_intel_check)
	CCC := icc
	CC := mpiicpc
	OPTIMIZATION_FLAGS := -DRICH_MPI -O2 -xHost -fp-model precise -std=c++11 -DOMPI_SKIP_MPICXX -D_GLIBCXX_DEBUG -fp-trap=common
	LINT_FLAGS =  -lclipper -lr3d
	ARCHIVER_FUNC := xiar
else ifeq ($(MODE),parallel)
	PREXECUTE := mpirun -np $(NUM_PROCESSES)
	CC := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI -O3 -DOMPI_SKIP_MPICXX -std=c++11
	LINT_FLAGS = -Werror -Wall -Wextra -pedantic -Wfatal-errors -Wmissing-declarations -Wno-long-long
	LFLAGS := -lclipper -lr3d
else ifeq ($(MODE),parallel_time)
	PREXECUTE := mpirun -np $(NUM_PROCESSES)
	CC := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI -O3 -DOMPI_SKIP_MPICXX -Dtiming -std=c++11
	LINT_FLAGS = -Werror -Wall -Wextra -pedantic -Wfatal-errors -Wmissing-declarations -Wno-long-long
	LFLAGS := -lclipper -lr3d
else ifeq ($(MODE),parallel_check)
	PREXECUTE := mpirun -np $(NUM_PROCESSES)
	CC := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI -Og -g -D_GLIBCXX_DEBUG -std=c++11
	LINT_FLAGS = -Werror -Wall -Wextra -pedantic -Wfatal-errors -Wmissing-declarations -Wno-long-long
	LFLAGS := -lclipper -lr3d
else ifeq ($(MODE),parallel_profile)
	PREXECUTE := mpirun -np $(NUM_PROCESSES)
	CC := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI -g -pg -O3 -std=c++11
	LINT_FLAGS = -Werror -Wall -Wextra -pedantic -Wfatal-errors -Wmissing-declarations -Wno-long-long
	LFLAGS := -lclipper -lr3d
else ifeq ($(MODE),debug_parallel)
	PREXECUTE := mpirun -np $(NUM_PROCESSES)
	CC := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI -O0 -g -pg -frecord-gcc-switches -DOMPI_SKIP_MPICXX -std=c++11
	LINT_FLAGS := -Wall -Wextra -pedantic -Wfatal-errors -Weffc++ -Wmissing-declarations
	LFLAGS := -pg -lclipper  -lr3d
else ifeq ($(MODE),parallel_check)
	PREXECUTE := mpirun -np $(NUM_PROCESSES)
	CC := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI -O3 -D_GLIBCXX_DEBUG -DOMPI_SKIP_MPICXX -std=c++11
	LINT_FLAGS = -Werror -Wall -Wextra -pedantic -Wfatal-errors -Wmissing-declarations -Wno-long-long
	LFLAGS := -lclipper -lr3d
else
	MODE := production
	OPTIMIZATION_FLAGS := -O3 -std=c++11
	LFLAGS := -lclipper -lr3d
endif
ifeq ($(DEBUG),on)
	OPTIMIZATION_FLAGS += -DRICH_DEBUG
endif
CFLAGS := $(LINT_FLAGS) $(OPTIMIZATION_FLAGS)

$(DUMP_FOLDER)/test.exe: $(ENDEMIC_OBJECTS) | ../library_$(MODE)/librich.a
	echo "$(ENDEMIC_SOURCES)"
	$(CC) $^ $| -o $@ $(LFLAGS) -lhdf5 -lhdf5_cpp -L $(HDF5_LIB_PATH)

../library_$(MODE)/librich.a:
ifeq ($(MODE),$(filter $(MODE),icpl icpl_mpi))
	$(MAKE) -C .. set_environ_vars_icpl.sh
	`cat ../set_environ_vars_icpl.sh` && \
	$(MAKE) -C .. MODE=$(MODE)
else
	$(MAKE) -C .. set_environ_vars.sh
	`cat ../set_environ_vars.sh` && \
	$(MAKE) -C .. MODE=$(MODE)
endif

$(ENDEMIC_OBJECTS): $(DUMP_FOLDER)/%.o: $(TARGET)/%.cpp
ifeq ($(MODE),$(filter $(MODE),icpl icpl_mpi))
	mkdir -p $(DUMP_FOLDER)
	$(MAKE) -C .. set_environ_vars_icpl.sh
	`cat ../set_environ_vars_icpl.sh` && \
	$(CC) $(CFLAGS) -c $< -o $@ -I ..
else
	mkdir -p $(DUMP_FOLDER)
	$(MAKE) -C .. set_environ_vars.sh
	`cat ../set_environ_vars.sh` && \
	$(CC) $(CFLAGS) -c $< -o $@ -I ..
endif
	
