SOURCE_DIR := source
RAW_SOURCES := $(shell find $(SOURCE_DIR) -name '*.cpp')
SOURCES := $(RAW_SOURCES)
LIB_FILE = librich.a
CC := gcc
CXX := g++
LINT_FLAGS = -Werror -Wall -Wextra -pedantic -Wno-long-long -Wfatal-errors -Weffc++ -Wshadow -Wmissing-declarations -Wconversion
ARCHIVER_FUNC := ar
ifeq ($(MODE),debug)
	OPTIMIZATION_FLAGS := -O0 -g -pg -std=c++11
	LINT_FLAGS :=
else ifeq ($(MODE),parallel)
	CXX := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI -O3 -std=c++11
	LINT_FLAGS = -Werror -Wall -Wextra -pedantic -Wfatal-errors -Weffc++ -Wshadow -Wmissing-declarations -Wno-long-long -Wno-effc++ -Wno-parentheses -Wno-reorder -Wno-shadow -Wconversion
else ifeq ($(MODE),parallel_profile)
	CXX := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI -O3 -g -pg -std=c++11
	LINT_FLAGS = -Werror -Wall -Wextra -pedantic -Wfatal-errors -Weffc++ -Wshadow -Wmissing-declarations -Wno-long-long -Wno-effc++ -Wno-parentheses -Wno-reorder -Wno-shadow -Wconversion
else ifeq ($(MODE),debug_parallel)
	CXX := mpiCC
	OPTIMIZATION_FLAGS := -DRICH_MPI -O0 -g -pg -frecord-gcc-switches -std=c++11
	LINT_FLAGS := 
else ifeq ($(MODE),intel)
	CXX := icpc
	OPTIMIZATION_FLAGS := -O3 -ipo -xHost -fp-model precise -std=c++11
	LINT_FLAGS := 
	ARCHIVER_FUNC := xiar
else ifeq ($(MODE),parallel_intel)
	CXX := mpiicpc
	OPTIMIZATION_FLAGS := -O3 -ipo -xHost -fp-model precise -DRICH_MPI -std=c++11
	LINT_FLAGS := 
	ARCHIVER_FUNC := xiar
else ifeq ($(MODE),clang)
	CXX := clang++
	OPTIMIZATION_FLAGS := -Weverything -Werror -ferror-limit=1 -Wno-error=padded -std=c++11 -Wno-zero-as-null-pointer-constant
	LINT_FLAGS := -std=c++11
else ifeq ($(MODE),TACC)
	CC := icc
	CXX := mpicxx
	OPTIMIZATION_FLAGS := -O3 -DRICH_MPI -ipo -xSKYLAKE -fp-model precise  -std=c++11 -DMPI_NO_CPPBIND
	LINT_FLAGS := 
	ARCHIVER_FUNC := xiar
else
	MODE = production
	OPTIMIZATION_FLAGS := -O3 -march=native -std=c++11
endif
LIBRARY_FOLDER := library_$(MODE)
OBJECTS := $(patsubst $(SOURCE_DIR)/%.cpp,$(LIBRARY_FOLDER)/%.o,$(SOURCES))
TREECODE_OBJECTS := $(patsubst $(SOURCE_DIR)/%.cpp,$(LIBRARY_FOLDER)/%.o,$(TREECODE_SOURCES))

$(LIBRARY_FOLDER)/$(LIB_FILE): $(OBJECTS)
	export RICH_ROOT=`pwd`
	$(ARCHIVER_FUNC) cr $@ $^

$(OBJECTS): $(LIBRARY_FOLDER)/%.o: $(SOURCE_DIR)/%.cpp
	mkdir -p `dirname $@`
	$(CXX) -c $(OPTIMIZATION_FLAGS) $(LINT_FLAGS) $< -o $@
	$(CXX) -MM $(OPTIMIZATION_FLAGS) $(LINT_FLAGS) $< -o $(LIBRARY_FOLDER)/$*.d
	@sed 's,\(\w*\)\.o,$@,g' -i $(LIBRARY_FOLDER)/$*.d

$(TREECODE_OBJECTS): $(LIBRARY_FOLDER)/%.o: $(SOURCE_DIR)/%.cpp
	mkdir -p `dirname $@`
	$(CXX) -c $(OPTIMIZATION_FLAGS) $(LINT_FLAGS) $< -o $@
	$(CXX) -MM $(OPTIMIZATION_FLAGS) $(LINT_FLAGS) $< -o $(LIBRARY_FOLDER)/$*.d
	@sed 's,\(\w*\)\.o,$(LIBRARY_FOLDER)/\1.o,g' -i $(LIBRARY_FOLDER)/$*.d

-include $(OBJECTS:.o=.d)

.PHONY: clean set_environ_vars.sh

clean:
	rm -rf ./$(LIBRARY_FOLDER)

set_environ_vars.sh: | external_libraries/include/H5Cpp.h external_libraries/boost_dump/boost_1_66_0/boost/container/static_vector.hpp external_libraries/ann_tree_dump/ann_1.1.2/lib/libANN.a external_libraries/lib/libclipper.a external_libraries/lib/libdclipper.a
	$(eval MY_BOOST_PATH=`pwd`/external_libraries/boost_dump/boost_1_66_0)
	$(eval MY_HDF5_PATH=`pwd`/external_libraries/include)
	$(eval MY_ANN_PATH=`pwd`/external_libraries/ann_tree_dump/ann_1.1.2/include)
	echo export\ CPLUS_INCLUDE_PATH=$(CPLUS_INCLUDE_PATH):$(MY_BOOST_PATH):$(MY_HDF5_PATH):$(MY_ANN_PATH) > set_environ_vars.sh
	echo export\ HDF5_LIB_PATH=`pwd`/external_libraries/lib >> set_environ_vars.sh
	echo export\ LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):`pwd`/external_libraries/lib:`pwd`/external_libraries/ann_tree_dump/ann_1.1.2/lib >> set_environ_vars.sh
	echo export\ LD_PATH=$(LD_PATH):`pwd`/external_libraries/lib:`pwd`/external_libraries/ann_tree_dump/ann_1.1.2/lib >> set_environ_vars.sh
	echo export\ RICH_ROOT=`pwd` >> set_environ_vars.sh

external_libraries/include/H5Cpp.h: external_libraries/hdf5_dump/hdf5-1.10.4/c++/src/H5Cpp.h
	cd external_libraries/hdf5_dump/hdf5-1.10.4 && \
	./configure --enable-cxx --prefix=`cd ../.. && pwd`
	cd external_libraries/hdf5_dump/hdf5-1.10.4 && make
	cd external_libraries/hdf5_dump/hdf5-1.10.4 && make check
	cd external_libraries/hdf5_dump/hdf5-1.10.4 && make install

external_libraries/include/clipper.hpp:
	mkdir -p external_libraries/dump_clipper
	cd external_libraries/dump_clipper && wget https://sourceforge.net/projects/polyclipping/files/latest/download && mv download clipper_ver6.4.2.zip && unzip clipper_ver6.4.2.zip && cp cpp/clipper.hpp ../include
#	cd external_libraries/dump_clipper && wget http://sourceforge.net/projects/polyclipping/files/latest/download?source=files && mv download?source=files clipper.zip && unzip clipper.zip && cp cpp/clipper.hpp ../include

external_libraries/dump_clipper/clipper.o: external_libraries/include/clipper.hpp
	cd external_libraries/dump_clipper && $(CXX) -c -O3 cpp/clipper.cpp -o clipper.o

external_libraries/lib/libclipper.a: external_libraries/dump_clipper/clipper.o
	$(ARCHIVER_FUNC) cr $@ $^ 

external_libraries/dump_clipper/dclipper.o: external_libraries/include/clipper.hpp
	cd external_libraries/dump_clipper && $(CXX) -c -O0 -g -pg -D_GLIBCXX_DEBUG cpp/clipper.cpp -o dclipper.o

external_libraries/lib/libdclipper.a: external_libraries/dump_clipper/dclipper.o
	$(ARCHIVER_FUNC) cr $@ $^ 

external_libraries/hdf5_dump/hdf5-1.10.4/c++/src/H5Cpp.h: | external_libraries/hdf5_dump/hdf5-1.10.4.tar.gz
	cd external_libraries/hdf5_dump/ && tar xf ./hdf5-1.10.4.tar.gz

external_libraries/boost_dump/boost_1_66_0/boost/container/static_vector.hpp: | external_libraries/boost_dump/boost_1_66_0.tar.bz2
	cd external_libraries/boost_dump/ && tar xf ./boost_1_66_0.tar.bz2

external_libraries/ann_tree_dump/ann_1.1.2/lib/libANN.a: | external_libraries/ann_tree_dump/ann_1.1.2/include/ANN/ANN.h
	cd external_libraries/ann_tree_dump/ann_1.1.2 && make linux-g++

external_libraries/ann_tree_dump/ann_1.1.2/include/ANN/ANN.h: | external_libraries/ann_tree_dump/ann_1.1.2.tar.gz
	cd external_libraries/ann_tree_dump/ && tar xf ./ann_1.1.2.tar.gz

external_libraries/hdf5_dump/hdf5-1.10.4.tar.gz:
	mkdir -p external_libraries/hdf5_dump
	cd external_libraries/hdf5_dump && \
	wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.gz

external_libraries/boost_dump/boost_1_66_0.tar.bz2:
	mkdir -p external_libraries/boost_dump
	cd external_libraries/boost_dump && \
	wget 'http://sourceforge.net/projects/boost/files/boost/1.66.0/boost_1_66_0.tar.bz2/download'
	cd external_libraries/boost_dump && mv download boost_1_66_0.tar.bz2

external_libraries/ann_tree_dump/ann_1.1.2.tar.gz:
	mkdir -p external_libraries/ann_tree_dump
	cd external_libraries/ann_tree_dump && \
	wget http://www.cs.umd.edu/~mount/ANN/Files/1.1.2/ann_1.1.2.tar.gz

