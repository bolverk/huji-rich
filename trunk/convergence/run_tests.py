#! /usr/bin/python

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

    def RunSingleTest(self, dname):
        """
        Runs a single test
        Input:
        dname - Path to directory
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
        automaker = imp.load_source('automaker',os.getcwd()+'/automaker.py')
        intrp = automaker.CppInterpreter(\
            'g++','-Wfatal-errors -Wall -Wextra -pedantic -Werror -std=c++0x -O3')
        os.system('g++ -c source/treecode/*.cpp')
        intrp.Interpret('test.cpp')
        mod = imp.load_source('mod',os.getcwd()+'/test.py')
        res = mod.main()
        os.chdir(myhome)
        return res

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

import sys
import os

tst = Tester('temp','source')

if len(sys.argv)==1:
    # Run all tests
    tst.DelTempDirs()
    ll = ListTestFolders('.')
    myhome = os.getcwd()
    for i in ll:
        try:
            tf = tst.RunSingleTest(i)
            if tf:
                print i+' ... passed'
            else:
                print i+' ... failed'
        except Exception, err:
            print i+' ... failed'
            print err
        os.chdir(myhome)
        tst.DelTempDirs()
elif len(sys.argv)==2:
    # Run a single test
    ll = sys.argv[1]
    try:
        tf = tst.RunSingleTest(ll)
        if tf:
            print ll+' ... passed'
        else:
            print ll+' ... failed'
    except Exception, err:
        print ll+' ... failed'
        print err
else:
    print 'unknown opion'
    
