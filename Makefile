SOURCE_DIR := source
RAW_SOURCES := $(shell find $(SOURCE_DIR) -name '*.cpp')
SOURCES := $(RAW_SOURCES)
LIB_FILE = librich.a
CC := g++
LINT_FLAGS = -Werror -Wall -Wextra -pedantic -Wno-long-long -Wfatal-errors -Weffc++ -Wshadow -Wmissing-declarations -Wconversion
ARCHIVER_FUNC := ar
ifeq ($(MODE),debug)
	OPTIMIZATION_FLAGS := -O0 -g -pg 
	LINT_FLAGS :=
else ifeq ($(MODE),parallel)
	CC := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI -O3
	LINT_FLAGS = -Werror -Wall -Wextra -pedantic -Wfatal-errors -Weffc++ -Wshadow -Wmissing-declarations -Wno-long-long -Wno-effc++ -Wno-parentheses -Wno-reorder -Wno-shadow -Wconversion
else ifeq ($(MODE),debug_parallel)
	CC := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI -O0 -g -pg -frecord-gcc-switches
	LINT_FLAGS := -Werror -Wall -Wextra -pedantic -Wfatal-errors -Weffc++ -Wshadow -Wmissing-declarations
else ifeq ($(MODE),intel)
	CC := icpc
	OPTIMIZATION_FLAGS := -O3 -ipo -xHost -fp-model precise
	LINT_FLAGS := 
	ARCHIVER_FUNC := xiar
else ifeq ($(MODE),clang)
	CC := clang++
	OPTIMIZATION_FLAGS := -Weverything -Werror -ferror-limit=1 -Wno-error=padded
	LINT_FLAGS := 
else
	MODE = production
	OPTIMIZATION_FLAGS := -O3 -march=native
endif
LIBRARY_FOLDER := library_$(MODE)
OBJECTS := $(patsubst $(SOURCE_DIR)/%.cpp,$(LIBRARY_FOLDER)/%.o,$(SOURCES))
TREECODE_OBJECTS := $(patsubst $(SOURCE_DIR)/%.cpp,$(LIBRARY_FOLDER)/%.o,$(TREECODE_SOURCES))

$(LIBRARY_FOLDER)/$(LIB_FILE): $(OBJECTS)
	export RICH_ROOT=`pwd`
	$(ARCHIVER_FUNC) cr $@ $^

$(OBJECTS): $(LIBRARY_FOLDER)/%.o: $(SOURCE_DIR)/%.cpp
	mkdir -p `dirname $@`
	$(CC) -c $(OPTIMIZATION_FLAGS) $(LINT_FLAGS) $< -o $@
	$(CC) -MM $(LINT_FLAGS) $< -o $(LIBRARY_FOLDER)/$*.d
	@sed 's,\(\w*\)\.o,$@,g' -i $(LIBRARY_FOLDER)/$*.d

$(TREECODE_OBJECTS): $(LIBRARY_FOLDER)/%.o: $(SOURCE_DIR)/%.cpp
	mkdir -p `dirname $@`
	$(CC) -c $(OPTIMIZATION_FLAGS) $< -o $@
	$(CC) -MM $< -o $(LIBRARY_FOLDER)/$*.d
	@sed 's,\(\w*\)\.o,$(LIBRARY_FOLDER)/\1.o,g' -i $(LIBRARY_FOLDER)/$*.d

-include $(OBJECTS:.o=.d)

.PHONY: clean set_environ_vars.sh

clean:
	rm -rf ./$(LIBRARY_FOLDER)

set_environ_vars.sh: | external_libraries/include/H5Cpp.h external_libraries/boost_dump/boost_1_59_0/boost/container/static_vector.hpp external_libraries/ann_tree_dump/ann_1.1.2/lib/libANN.a
	$(eval MY_BOOST_PATH=`pwd`/external_libraries/boost_dump/boost_1_59_0)
	$(eval MY_HDF5_PATH=`pwd`/external_libraries/include)
	$(eval MY_ANN_PATH=`pwd`/external_libraries/ann_tree_dump/ann_1.1.2/include)
	echo export\ CPLUS_INCLUDE_PATH=$(CPLUS_INCLUDE_PATH):$(MY_BOOST_PATH):$(MY_HDF5_PATH):$(MY_ANN_PATH) > set_environ_vars.sh
	echo export\ HDF5_LIB_PATH=`pwd`/external_libraries/lib >> set_environ_vars.sh
	echo export\ LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):`pwd`/external_libraries/lib:`pwd`/external_libraries/ann_tree_dump/ann_1.1.2/lib >> set_environ_vars.sh
	echo export\ LD_PATH=$(LD_PATH):`pwd`/external_libraries/lib:`pwd`/external_libraries/ann_tree_dump/ann_1.1.2/lib >> set_environ_vars.sh
	echo export\ RICH_ROOT=`pwd` >> set_environ_vars.sh

external_libraries/include/H5Cpp.h: external_libraries/hdf5_dump/hdf5-1.8.15-patch1/c++/src/H5Cpp.h
	cd external_libraries/hdf5_dump/hdf5-1.8.15-patch1 && \
	./configure --enable-cxx --prefix=`cd ../.. && pwd`
	cd external_libraries/hdf5_dump/hdf5-1.8.15-patch1 && make
	cd external_libraries/hdf5_dump/hdf5-1.8.15-patch1 && make install

external_libraries/hdf5_dump/hdf5-1.8.15-patch1/c++/src/H5Cpp.h: | external_libraries/hdf5_dump/hdf5-1.8.15-patch1.tar.gz
	cd external_libraries/hdf5_dump/ && tar xvf ./hdf5-1.8.15-patch1.tar.gz

external_libraries/boost_dump/boost_1_59_0/boost/container/static_vector.hpp: | external_libraries/boost_dump/boost_1_59_0.tar.bz2
	cd external_libraries/boost_dump/ && tar xvf ./boost_1_59_0.tar.bz2

external_libraries/ann_tree_dump/ann_1.1.2/lib/libANN.a: | external_libraries/ann_tree_dump/ann_1.1.2/include/ANN/ANN.h
	cd external_libraries/ann_tree_dump/ann_1.1.2 && make linux-g++

external_libraries/ann_tree_dump/ann_1.1.2/include/ANN/ANN.h: | external_libraries/ann_tree_dump/ann_1.1.2.tar.gz
	cd external_libraries/ann_tree_dump/ && tar xvf ./ann_1.1.2.tar.gz

external_libraries/hdf5_dump/hdf5-1.8.15-patch1.tar.gz:
	mkdir -p external_libraries/hdf5_dump
	cd external_libraries/hdf5_dump && \
	wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.15-patch1.tar.gz

external_libraries/boost_dump/boost_1_59_0.tar.bz2:
	mkdir -p external_libraries/boost_dump
	cd external_libraries/boost_dump && \
	wget 'http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.tar.bz2/download'
	cd external_libraries/boost_dump && mv download boost_1_59_0.tar.bz2

external_libraries/ann_tree_dump/ann_1.1.2.tar.gz:
	mkdir -p external_libraries/ann_tree_dump
	cd external_libraries/ann_tree_dump && \
	wget http://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz

