#! /usr/bin/python

def IsPureVirtual(fname):
    """
    Returns true if a header file is a pure virtual (no implementation)
    Pure virtual header files are designated as such by the phrase "Pure virtual"
    at the first line
    Input:
    fname - File name
    """

    f = open(fname)
    fl = f.readline()
    f.close()
    return fl.find('Pure virtual')>-1

def RemoveFileFromPath(fpath):
    """
    Removes a file name from the path
    Input:
    fpath - Path to file
    """

    if fpath.find('/')>-1:
        rfpath = fpath[::-1]
        rfpathnf = rfpath[rfpath.find('/'):]
        fpathnf = rfpathnf[::-1]
        res = fpathnf

        return res
    else:
        return ''

def JustFileName(fpath):
    """
    Extracts the file name from the path
    Input:
    fpath - Full path
    """

    if(fpath.find('/')>-1):
        rfpath = fpath[::-1]
        rfname = rfpath[:rfpath.find('/')]
        fname = rfname[::-1]
        return fname
    else:
        return fpath

def JoinPaths(apath, rpath):
    """
    Joins a relative path to an absolute path
    Input:
    apath - Absolute path
    rpath - Relative path
    """

    if len(rpath)==0:
        res = apath
    elif len(apath)==0:
        res = './'+rpath
    elif rpath[0]=='/' or apath[-1]=='/':
        res = apath+rpath
    else:
        res = apath+'/'+rpath
    while(res.find('..')>-1):
        aux1 = res.split('/')
        aux2 = aux1.index('..')
        aux3 = aux1[aux2-1]
        res = res.replace(aux3+'/../','')
    return res

class CppInterpreter:
    """
    Interprets c++ code
    """

    def __init__(self, cname, cflags):
        """
        Class constructor
        Input:
        cname - Compiler name, i.e. string used to invoke compiler
        cflag - Compilation flags
        """
        self.__cname = cname
        self.__cflags = cflags

    def SimpleCompile(self, sfname):
        """
        Creates an object file from a source file
        Input:
        sfname - Name of the source file
        """

        import os
        ef = os.system(self.__cname+' '+self.__cflags+' -c '+sfname)
        if ef!=0:
            raise NameError('Failed to compile')

    def NeedCompile(self,hfname):
        """
        Checkes whether a certain file needs to be compiled
        Returns true if so, and false if the file has already
        been compile (if a newer object file exists)
        """

        if 'treecode' in hfname:
            return False

        res = 1
        if hfname[-4:]!='.hpp':
            print 'NeedCompile only accepts files with the extension .hpp'
            print 'Input file is '+hfname
            raise NameError('NotHeader')
        import os
        if os.path.exists(JustFileName(hfname).replace('.hpp','.o')):
            t_hpp = os.path.getmtime(hfname)
            t_o = os.path.getmtime(JustFileName(hfname).replace('.hpp','.o'))
            t_cpp = os.path.getmtime(hfname.replace('.hpp','.cpp'))
            if (t_o>t_cpp) and (t_o>t_hpp):
                res = 0
        return res                                   

    def RecursiveCompile(self, sfname):
        """
        Checks a source file for dependencies and compiles them recursively
        """

        f = open(sfname)
        for i in f:
            if i.find('#include')>-1 and i.find('"')>i.find('#define'):
                hname = i.split('"')[1]
                fpath = JoinPaths(RemoveFileFromPath(sfname),
                                  RemoveFileFromPath(hname))
                fname = hname.replace('hpp','cpp')
                if self.NeedCompile(fpath+JustFileName(hname)) and \
                        JustFileName(sfname)!=fname:
                    self.RecursiveCompile(fpath+JustFileName(hname))
                    if not IsPureVirtual(fpath+JustFileName(hname)):
                        self.RecursiveCompile(fpath+JustFileName(fname))
        f.close()
        if(sfname[-4:]=='.cpp'):
            self.SimpleCompile(sfname)

    def Link(self, sfname):
        """
        Links object files to an executable
        Input:
        sfname - Name of main file
        To - Do: The function assumes the input file has the extension .cpp
        At the moment I don't think that in this project I would need other
        extensions, so I decided to raise an error in case a file with a
        different extension is entered
        """

        if sfname.find('.cpp')==-1:
            raise NameError('automaker currently supports only .cpp files')

        import os
        ef = os.system(self.__cname+' -o '+sfname.replace('.cpp','.exe')+' *.o')
        if ef!=0:
            raise NameError('Failed to link '+sfname.replace('.cpp','.exe'))

    def Run(self, sfname):
        """
        Runs the executable
        This function assumes the executable was created using the Link function
        Input:
        sfname - Name of the main source file
        """

        if sfname.find('.cpp')==-1:
            raise NameError('automaker currently supports only .cpp files')

        import os
        os.system('./'+sfname.replace('.cpp','.exe'))

    def Interpret(self, sfname):
        """
        Compiles, links and executes
        Input:
        sfname - Name of the main source file
        """

        if sfname.find('.cpp')==-1:
            raise NameError('automaker currently supports only .cpp files')

        self.RecursiveCompile(sfname)
        self.Link(sfname)
        self.Run(sfname)

import sys
if __name__=="__main__":
    fname = sys.argv[1]

    import os
    os.system('rm *.o *.exe')
    myinterp = CppInterpreter('g++','')
    myinterp.Interpret(fname)
