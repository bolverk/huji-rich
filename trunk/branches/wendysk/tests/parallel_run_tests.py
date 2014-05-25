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

def run_single_test(fpath):

    import os
    import imp

    temp_folder = make_temp_folder('temp',
                                   fpath)
    safe_remove_folder(temp_folder)
    os.mkdir(temp_folder)
    os.system('cp -r ../source/ '+temp_folder)
    os.system('cp '+fpath+'/* '+temp_folder)
    os.system('cp automaker.py ./'+temp_folder)
    my_home = os.getcwd()
    os.chdir(temp_folder)
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
    automaker = imp.load_source(
        'automaker',
        os.getcwd()+'/automaker.py')
    intrp = automaker.CppInterpreter(
        'g++','-I '+hdf5_inc_path+' '+
        '-Wfatal-errors '+
        '-Wmissing-declarations '+
        '-Wendif-labels '+
        '-Wabi '+
        '-Wctor-dtor-privacy '+
        '-Wstrict-null-sentinel '+
        '-Winvalid-offsetof '+
        '-Wsign-promo '+
        '-Woverloaded-virtual '+
        '-Wsynth '+
        '-Wlogical-op '+
        '-Wno-unused-result '+
        '-Werror-implicit-function-declaration '+
        '-Weffc++ '+
        '-Wno-missing-braces '+
        '-Wno-missing-field-initializers '+
        '-Wformat=2 '+
        '-Wswitch-default '+
        '-Wswitch-enum '+
        '-Wpointer-arith '+
        '-Wshadow '+
        '-pedantic-errors '+
        '-Wcast-align '+
        '-Wcast-qual '+
        '-Wchar-subscripts '+
        '-Wcomment '+
        '-Wconversion '+
        '-Wdisabled-optimization '+
        '-Wformat-nonliteral -Wformat-security '+
        '-Wformat-y2k '+
        '-Wimport  -Winit-self '+
        '-Winvalid-pch '+
        '-Wundef '+
        '-Wmissing-braces '+
        '-Wmissing-field-initializers -Wmissing-format-attribute '+
        '-Wmissing-include-dirs '+
        '-Wpacked -Wparentheses '+
        '-Wredundant-decls -Wreturn-type '+
        '-Wsequence-point -Wsign-compare  -Wstack-protector '+
        '-Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch  -Wswitch-default '+
        '-Wswitch-enum -Wtrigraphs  -Wuninitialized '+
        '-Wunknown-pragmas  -Wunreachable-code -Wunused '+
        '-Wunused-function  -Wunused-label  -Wunused-parameter '+
        '-Wunused-value  -Wunused-variable  -Wvariadic-macros '+
        '-Wvolatile-register-var  -Wwrite-strings '+
        '-Wall -Wextra -pedantic -Werror -std=c++0x -O3',
        '-L '+hdf5_lib_path+'libhdf5_hl_cpp.a '+
        hdf5_lib_path+'libhdf5_cpp.a '+
        hdf5_lib_path+'libhdf5_hl.a '+
        hdf5_lib_path+'libhdf5.a '+
        '-lpthread -lz -ldl')
    os.system('g++ -c source/treecode/*.cpp')
    intrp.Interpret('test.cpp')

    mod = imp.load_source(
        'mod',
        os.getcwd()+'/test.py')
    res = mod.main()
    os.chdir(my_home)
    return res
        
def run_single_test_weh(fpath,home_dir):

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
        return fpath+' ... failed due to errors'

def main():

    import multiprocessing
    import os

    os.system('rm -rf temp_*')
    ll = ListTestFolders('.')

    pool = multiprocessing.Pool()

    msgs = []
    [pool.apply_async(run_single_test_weh,
                      (l,os.getcwd()),
                      callback=msgs.append)
     for l in ll]
    pool.close()
    pool.join()

    f = open('test_report.txt','w')
    for i in msgs:
        f.write(i+'\n')
    f.close()

    os.system('rm -rf temp_*')

if __name__=='__main__':

    main()
