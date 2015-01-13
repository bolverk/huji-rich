def build_library():

    import os

    os.chdir('..')
    os.system('python ./homebrew_makefile.py')
    os.chdir('convergence')

def temp_folder_name(test_path):

    import re
    
    res = 'temp_'+test_path
    res = res.replace('.','')
    res = res.replace('/','_')
    res = res.replace('__','_')
    res = re.sub(r'_$','',res)
    return res

def new_empty_directory(dirname):

    import os

    os.system('touch '+dirname)
    os.system('rm -rf '+dirname)
    os.system('mkdir '+dirname)

def main():

    import argparse
    import os
    import imp
    homebrew_makefile = imp.load_source(\
        'homebrew_makefile',
        '../homebrew_makefile.py')

    build_library()

    parser = argparse.ArgumentParser()
    parser.add_argument('test_folder',
                        help='path to test folder')
    parser.add_argument('-jb','--just_build',
                        action='store_true',
                        help='Only builds the test, but does not run it')
    args = parser.parse_args()
    
    tfname = temp_folder_name(args.test_folder)
    new_empty_directory(tfname)
    os.system('cp '+
              os.path.join((args.test_folder+'/*').replace('//','/'))+' '+
              tfname)

    # Compile
    os.system('g++ -o '+tfname+'/test.o '+
              '-c -O3 '+
              homebrew_makefile.gcc_strict_warning_flags()+' '+
              (args.test_folder+'/test.cpp').replace('//','/')+' '+
              '-I ..')
    
    # Link
    hdf5_lib_path = os.environ['HDF5_LIB_PATH']
    os.system('g++ -o '+tfname+'/test.exe '+
              tfname+'/test.o '+
              '../library/librich.a '+
              (hdf5_lib_path+'/libhdf5_cpp.a').replace('//','/')+' '+
              (hdf5_lib_path+'/libhdf5_hl.a').replace('//','/')+' '+
              (hdf5_lib_path+'/libhdf5.a').replace('//','/')+' '+
              '-lpthread -lz -ldl ')

    if not args.just_build:
        os.chdir(tfname)
        os.system('./test.exe')
        import imp
        mod = imp.load_source('mod',os.getcwd()+'/test.py')
        try:
            if mod.main():
                print args.test_folder+' ... passed'
            else:
                print args.test_folder+' ... ran and failed'
        except:
            print args.test_folder+' ... failed to run'
        os.chdir('..')        

if __name__ == '__main__':

    main()
