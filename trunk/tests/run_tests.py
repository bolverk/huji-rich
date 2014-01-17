#! /usr/bin/python

def gcc_hdf5_include_flags():

    import os
    import sys

    try:
        hdf5_inc_path = os.environ['HDF5_INC_PATH']
        hdf5_inc_path = (hdf5_inc_path+'/').replace('//','/')
    except KeyError:
        print 'You did not specify path to hdf5 include files'
        sys.exit(1)
    return ' -I '+hdf5_inc_path+' '

def gcc_hdf5_library_flags():

    import os
    import sys

    try:
        hdf5_lib_path = os.environ['HDF5_LIB_PATH']
        hdf5_lib_path = (hdf5_lib_path+'/').replace('//','/')
    except KeyError:
        print 'You did not specify path to the hdf5 library'
        sys.exit(1)
    return ' -L '+hdf5_lib_path+'libhdf5_hl_cpp.a '+\
        hdf5_lib_path+'libhdf5_cpp.a '+\
        hdf5_lib_path+'libhdf5_hl.a '+\
        hdf5_lib_path+'libhdf5.a '+\
        '-lpthread -lz -ldl '

def gcc_strict_warning_flags():

    return '-Wfatal-errors '+\
        '-Wmissing-declarations '+\
        '-Wendif-labels '+\
        '-Wabi '+\
        '-Wctor-dtor-privacy '+\
        '-Wstrict-null-sentinel '+\
        '-Winvalid-offsetof '+\
        '-Wsign-promo '+\
        '-Woverloaded-virtual '+\
        '-Wsynth '+\
        '-Wlogical-op '+\
        '-Wno-unused-result '+\
        '-Werror-implicit-function-declaration '+\
        '-Weffc++ '+\
        '-Wno-missing-braces '+\
        '-Wno-missing-field-initializers '+\
        '-Wformat=2 '+\
        '-Wswitch-default '+\
        '-Wswitch-enum '+\
        '-Wpointer-arith '+\
        '-Wshadow '+\
        '-pedantic-errors '+\
        '-Wcast-align '+\
        '-Wcast-qual '+\
        '-Wchar-subscripts '+\
        '-Wcomment '+\
        '-Wconversion '+\
        '-Wdisabled-optimization '+\
        '-Wformat-nonliteral -Wformat-security '+\
        '-Wformat-y2k '+\
        '-Wimport  -Winit-self '+\
        '-Winvalid-pch '+\
        '-Wundef '+\
        '-Wmissing-braces '+\
        '-Wmissing-field-initializers -Wmissing-format-attribute '+\
        '-Wmissing-include-dirs '+\
        '-Wpacked -Wparentheses '+\
        '-Wredundant-decls -Wreturn-type '+\
        '-Wsequence-point -Wsign-compare  -Wstack-protector '+\
        '-Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch  -Wswitch-default '+\
        '-Wswitch-enum -Wtrigraphs  -Wuninitialized '+\
        '-Wunknown-pragmas  -Wunreachable-code -Wunused '+\
        '-Wunused-function  -Wunused-label  -Wunused-parameter '+\
        '-Wunused-value  -Wunused-variable  -Wvariadic-macros '+\
        '-Wvolatile-register-var  -Wwrite-strings '+\
        '-Wall -Wextra -pedantic -Werror'

def make_temp_folder(pref,uname):

    import re

    res = pref+'_'+uname
    res = res.replace('.','')
    res = res.replace('/','_')
    res = res.replace('__','_')
    res = re.sub(r'_$','',res)
    return res

def safe_remove_folder(fpath):
    
    import os
    os.system('touch '+fpath)
    os.system('rm -rf '+fpath)

class Tester:
    """
    A class that manages automatic tests
    """

    def __init__(self, tfname, srcdir):
        """
        Class constructor
        Input:
        tfname - Name of temporary folder
        scrdir - Path to source
        """

        self.__tfname = tfname
        self.__srcdir = srcdir
        return

    def DelTempDirs(self):
        """
        Deletes the temporary folder and all files within
        """

        import os
        os.system('rm -rf '+self.__tfname+'*')
        return

    def RunSingleTest(self, dname,aux_flags=''):
        """
        Runs a single test
        Input:
        dname - Path to directory
        aux_flags - Flags for test script
        """

        import os
        import imp

        temp_folder = make_temp_folder(self.__tfname,dname)
        safe_remove_folder(temp_folder)
        os.mkdir(temp_folder)
        os.system('cp -r ../'+self.__srcdir+' '+temp_folder)
        os.system('cp '+dname+'/* '+temp_folder)
        os.system('cp ./automaker.py '+temp_folder)
        myhome = os.getcwd()
        os.chdir(temp_folder)
        # List of warning flags taken from
        # http://stackoverflow.com/questions/399850/best-compiler-warning-level-for-c-c-compilers
        try:
            hdf5_inc_path = os.environ['HDF5_INC_PATH']
            hdf5_inc_path = (hdf5_inc_path+'/').replace('//','/')
        except:
            print 'You did not specify path to hdf5 include files'
            sys.exit(1)
        try:
            hdf5_lib_path = os.environ['HDF5_LIB_PATH']
            hdf5_lib_path = (hdf5_lib_path+'/').replace('//','/')
        except KeyError:
            print 'You did not specify path to the hdf5 library'
            sys.exit(1)
        automaker = imp.load_source('automaker',os.getcwd()+'/automaker.py')
        intrp = automaker.CppInterpreter(\
            'g++',gcc_hdf5_include_flags()+' '+
            gcc_strict_warning_flags()+' '+
            '-std=c++0x -O3',
            gcc_hdf5_library_flags())
        os.system('g++ -c source/treecode/*.cpp')
        intrp.Interpret('test.cpp')
        mod = imp.load_source('mod',os.getcwd()+'/test.py')
        if aux_flags=='':
            res = mod.main()
        else:
            res = mod.main(aux_flags)
        os.chdir(myhome)
        return res

    def BuildSingleTest(self, dname):
        """
        Runs a single test
        Input:
        dname - Path to directory
        aux_flags - Flags for test script
        """

        import os
        import imp
        import sys

        temp_folder = make_temp_folder(self.__tfname,dname)
        safe_remove_folder(temp_folder)
        os.mkdir(temp_folder)
        os.system('cp -r ../'+self.__srcdir+' '+temp_folder)
        os.system('cp '+dname+'/* '+temp_folder)
        os.system('cp ./automaker.py '+temp_folder)
        myhome = os.getcwd()
        os.chdir(temp_folder)
        # List of warning flags taken from
        # http://stackoverflow.com/questions/399850/best-compiler-warning-level-for-c-c-compilers
        try:
            hdf5_inc_path = os.environ['HDF5_INC_PATH']
            hdf5_inc_path = (hdf5_inc_path+'/').replace('//','/')
        except:
            print 'You did not specify path to hdf5 include files'
            sys.exit(1)
        try:
            hdf5_lib_path = os.environ['HDF5_LIB_PATH']
            hdf5_lib_path = (hdf5_lib_path+'/').replace('//','/')
        except KeyError:
            print 'You did not specify path to the hdf5 library'
            sys.exit(1)
        automaker = imp.load_source('automaker',os.getcwd()+'/automaker.py')
        intrp = automaker.CppInterpreter(\
            'g++',gcc_hdf5_include_flags()+' '+
            gcc_strict_warning_flags()+' '+
            '-std=c++0x -O3',
            gcc_hdf5_library_flags())
        os.system('g++ -c source/treecode/*.cpp')
        intrp.RecursiveCompile('test.cpp')
        intrp.Link('test.cpp')

def ListTestFolders(dname):
    """
    Returns a list of folders with test files
    Input:
    dname - Root folder
    """

    import os
    res  = []
    for roots, dirs, files in os.walk(dname):
        if roots.find('.svn')==-1 and len(roots)>1:
            if os.path.isfile(roots+'/test.py'):
                res.append(roots)
    return res

def main():

    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("test_folder",
                        help='Path to test folder')
    parser.add_argument("-jb","--just_build",
                        action="store_true",
                        help='Only builds the test, but does not run it')
    parser.add_argument("-af","--aux_flags",
                        help='Auxiliary flags for the analysis')
    args = parser.parse_args()

    tst = Tester('temp','source')

    if args.just_build:
        tst.BuildSingleTest(args.test_folder)
    else:
        aux_flags = ''
        if args.aux_flags:
            aux_flags = args.aux_flags
        try:
            tf = tst.RunSingleTest(args.test_folder,
                                   aux_flags)
            if tf:
                print args.test_folder+' ... passed'
            else:
                print args.test_folder+' ... failed'
        except Exception, err:
            print args.test_folder+' ... failed'
            print err

if __name__=='__main__':

    main()
