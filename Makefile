SOURCE_DIR := source
RAW_SOURCES := $(shell find $(SOURCE_DIR) -name '*.cpp')
TREECODE_SOURCES := $(shell find $(SOURCE_DIR)/treecode -name '*.cpp')
SOURCES := $(filter-out $(TREECODE_SOURCES),$(RAW_SOURCES))
LIB_FILE = librich.a
CC := g++
LINT_FLAGS = -Werror -Wextra -pedantic -Wfatal-errors -Weffc++ -Wshadow -Wmissing-declarations
ARCHIVER_FUNC := ar
ifeq ($(MODE),debug)
	OPTIMIZATION_FLAGS := -O0 -g -pg 
else ifeq ($(MODE),parallel)
	CC := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI
	LINT_FLAGS = -Wextra -pedantic -Wfatal-errors -Weffc++ -Wshadow -Wmissing-declarations
else ifeq ($(MODE),intel)
	CC := icpc
	OPTIMIZATION_FLAGS := -O3 -ipo -xHost -fp-model precise
	LINT_FLAGS := 
	ARCHIVER_FUNC := xiar
else
	MODE = production
	OPTIMIZATION_FLAGS := -O2
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
	$(CC) -MM $(CFLAGS) $< -o $(LIBRARY_FOLDER)/$*.d
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

set_environ_vars.sh: | external_libraries/include/H5Cpp.h external_libraries/boost_dump/boost_1_55_0/boost/container/static_vector.hpp external_libraries/boost_dump/boost_1_55_0/stage/lib/libboost_mpi.a
	$(eval MY_BOOST_PATH=`pwd`/external_libraries/boost_dump/boost_1_55_0)
	$(eval MY_HDF5_PATH=`pwd`/external_libraries/include)
	echo export\ CPLUS_INCLUDE_PATH=$(CPLUS_INCLUDE_PATH):$(MY_BOOST_PATH):$(MY_HDF5_PATH) > set_environ_vars.sh
	echo export\ HDF5_LIB_PATH=`pwd`/external_libraries/lib >> set_environ_vars.sh
	echo export\ BOOST_LIB_PATH=`pwd`/external_libraries/boost_dump/boost_1_55_0/stage/lib >> set_environ_vars.sh
	echo export\ LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):`pwd`/external_libraries/lib:`pwd`/external_libraries/boost_dump/boost_1_55_0/stage/lib >> set_environ_vars.sh
	echo export\ LD_PATH=$(LD_PATH):`pwd`/external_libraries/lib:`pwd`/external_libraries/boost_dump/boost_1_55_0/stage/lib >> set_environ_vars.sh
	echo export\ RICH_ROOT=`pwd` >> set_environ_vars.sh

external_libraries/include/H5Cpp.h: external_libraries/hdf5_dump/hdf5-1.8.13/c++/src/H5Cpp.h
	cd external_libraries/hdf5_dump/hdf5-1.8.13 && \
	./configure --enable-cxx --prefix=`cd ../.. && pwd`
	cd external_libraries/hdf5_dump/hdf5-1.8.13 && make
	cd external_libraries/hdf5_dump/hdf5-1.8.13 && make install

external_libraries/hdf5_dump/hdf5-1.8.13/c++/src/H5Cpp.h: | external_libraries/hdf5_dump/hdf5-1.8.13.tar.gz
	cd external_libraries/hdf5_dump/ && tar xvf ./hdf5-1.8.13.tar.gz

external_libraries/boost_dump/boost_1_55_0/stage/lib/libboost_mpi.a: external_libraries/boost_dump/boost_1_55_0/boost/container/static_vector.hpp
	cd external_libraries/boost_dump/ && \
	./bootstrap.sh --prefix=`cd ../.. && pwd` && \
	if ! grep -q "using mpi ;" tools/build/v2/user-config.jam
	then
		echo "using mpi ;" >> tools/build/v2/user-config.jam
	fi && \
	./b2 --with-mpi && bjam --with-mpi install

external_libraries/boost_dump/boost_1_55_0/boost/container/static_vector.hpp: | external_libraries/boost_dump/boost_1_55_0.tar.bz2
	cd external_libraries/boost_dump/ && tar xvf ./boost_1_55_0.tar.bz2

external_libraries/hdf5_dump/hdf5-1.8.13.tar.gz:
	mkdir -p external_libraries/hdf5_dump
	cd external_libraries/hdf5_dump && \
	wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.13.tar.gz

external_libraries/boost_dump/boost_1_55_0.tar.bz2:
	mkdir -p external_libraries/boost_dump
	cd external_libraries/boost_dump && \
	wget 'http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.bz2/download'
	cd external_libraries/boost_dump && mv download boost_1_55_0.tar.bz2

