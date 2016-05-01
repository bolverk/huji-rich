def find_files(directory, pattern):

    import os, fnmatch

    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename

import os

debug = ARGUMENTS.get('debug',0)
compiler = ARGUMENTS.get('compiler','g++')

linkflags = ''
if compiler=='g++':
    cflags = '-Wfatal-errors'
    if int(debug):
        cflags +=' -O0 -g -pg'
        linkflags = ' -g -pg'
    else:
        cflags +=' -O3'
elif compiler=='clang++':
    cflags = '-Weverything -Werror -ferror-limit=1 -Wno-error=padded'
    if int(debug):
        cflags += ' -O0 -g -pg'
        linkflags = ' -g -pg'
    else:
		cflags += ' -O3 -march=native'
elif compiler=='profile':
    cflags = ' -Wfatal-errors '
    cflags += ' -O3 -g'
    linkflags = ' -g '
elif compiler=='icc':
    cflags = ' -O3 -ipo -xHost -fp-model precise -g '
    linkflags = ' -g '
elif compiler=='mpiCC':
    cflags = '-DRICH_MPI -Wfatal-errors '
    if int(debug):
        cflags += ' -O0 -g -pg'
    else:
        cflags += ' -O3 -march=native'
else:
    raise NameError('unsupported compiler')
if compiler=='profile':
	compiler = 'g++'

build_dir = 'build/'+compiler
if int(debug):
    build_dir += '/debug'
else:
    build_dir += '/release'

env = Environment(ENV = os.environ,
                  CXX=compiler,
                  CPPPATH=[os.environ['RICH_ROOT']+'/source',
                           os.environ['RICH_ROOT']],
                  LIBPATH=['.',os.environ['HDF5_LIB_PATH']],
                  LIBS=['rich','hdf5','hdf5_cpp'],
                  LINKFLAGS=linkflags,
                  CXXFLAGS=cflags)
env.VariantDir(build_dir,'source')
source_file_list = [build_dir+'/'+fname.replace('source/','') for fname in find_files('source','*.cpp')]
lib_file = env.Library(build_dir+'/rich',source_file_list)
