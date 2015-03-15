def find_files(directory, pattern):

    import os, fnmatch

    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename

import os

debug = ARGUMENTS.get('debug',0)
source = ARGUMENTS.get('source',None)
target = ARGUMENTS.get('target',None)
compiler = ARGUMENTS.get('compiler','g++')

if compiler=='g++':
    cflags = '-Wfatal-errors'
    if int(debug):
        cflags +=' -O0 -g -pg'
    else:
        cflags +=' -O3'
elif compiler=='clang++':
    cflags = '-Weverything -Werror -ferror-limit=1 -Wno-error=padded'
    if int(debug):
        cflags += ' -O0 -g -pg'
    else:
        cflags += ' -O3'
else:
    raise NameError('unsupported compiler')

env = Environment(ENV = os.environ,
                  CXX=compiler,
                  CPPPATH=os.environ['RICH_ROOT']+'/source',
                  LIBPATH=['.',os.environ['HDF5_LIB_PATH']],
                  LIBS=['rich','hdf5','hdf5_cpp'],
                  CXXFLAGS=cflags)
source_file_list = [fname for fname in find_files(os.environ['RICH_ROOT']+'/source','*.cpp')]
lib_file = env.Library('rich',source_file_list)
if None!=source and None!=target:
    env.Program(target,[source])
