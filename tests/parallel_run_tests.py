def run_single_test(fpath):

    import os
    import run_tests as nrt
    from importlib.machinery import SourceFileLoader
    homebrew_makefile = SourceFileLoader(\
        'homebrew_makefile',
        '../homebrew_makefile.py').load_module()

    tfname = nrt.temp_folder_name(fpath)
    nrt.new_empty_directory(tfname)
    os.system('cp '+
              os.path.join(fpath+'/*').replace('//','/')+' '+
              tfname)

    # Compile
    os.system('g++ -o '+tfname+'/test.o '+
              '-c -O3 '+
              homebrew_makefile.gcc_strict_warning_flags()+' '+
              (fpath+'/test.cpp').replace('//','/')+' '+
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

    # Run
    os.chdir(tfname)
    os.system('./test.exe')
    from importlib.machinery import SourceFileLoader
    mod = SourceFileLoader('mod',os.getcwd()+'/test.py').load_module()
    res = mod.main()
    os.chdir('..')
    return res

def run_single_test_weh(fpath, home_dir):

    import os
    
    try:
        tf = run_single_test(fpath)
        if tf:
            return fpath+' ... passed'
        else:
            return fpath+' ... ran and failed'
    except Exception, err:
        print err
        os.chdir(home_dir)
        return fpath+' ... failed due to error'

def get_all_file_paths():

    import os

    res = []
    for root, dirs, files in os.walk('.'):
        for name in files:
            res.append(os.path.join(root, name))
    return res

def get_source_paths():

    return [fname for fname in get_all_file_paths()
            if 'test.cpp' in fname and not '.svn' in fname]

def get_test_paths():

    return [fpath.replace('/test.cpp','') for fpath in get_source_paths() if not '~' in fpath]

def main():

    import multiprocessing
    import os
    import run_tests as nrt

    nrt.build_library()

    msgs = []

    pool = multiprocessing.Pool()
    [pool.apply_async(run_single_test_weh,
                      (l,os.getcwd()),
                      callback=msgs.append)
     for l in get_test_paths()]
    pool.close()
    pool.join()

    f = open('test_report.txt','w')
    for i in msgs:
        f.write(i+'\n')
    f.close()

if __name__ == '__main__':
    
    main()
