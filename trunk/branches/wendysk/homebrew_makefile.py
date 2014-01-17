def gcc_strict_warning_flags():

    return '-Wfatal-errors '+\
        '-std=c++0x '+\
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
        '-Wunknown-pragmas '+\
        '-Wunused '+\
        '-Wunused-function  -Wunused-label  -Wunused-parameter '+\
        '-Wunused-value  -Wunused-variable  -Wvariadic-macros '+\
        '-Wvolatile-register-var  -Wwrite-strings '+\
        '-Wall -Wextra -pedantic -Werror'

def safe_create(folder_name):

    import os

    if not os.path.isdir(folder_name):
        os.system('mkdir '+folder_name)

def get_all_files(absroot):

    import os

    res = []
    for root, dirs, files in os.walk(absroot):
        for name in files:
            res.append(os.path.join(root,name))
    return res

def get_source_files(absroot):

    return [fname for fname in get_all_files(absroot)
            if fname[-4:]=='.cpp']

def source2lib(fname, lname):

    return lname+'/'+fname.split('/')[-1].replace('.cpp','.o')

def strict_compile(fname, lname):

    import os

    oname = source2lib(fname, lname)
    os.system('g++ -c '+fname+' '+
              '-O3 '+
              gcc_strict_warning_flags()+' '+
              '-o '+oname)

def simple_compile(fname, lname):

    import os

    oname = source2lib(fname, lname)
    os.system('g++ -c '+fname+' '+
              '-O3 '+
              '-o '+oname)

def sensitive_compile(fname, lname):

    if '/treecode/' in fname:
        simple_compile(fname, lname)
    elif 'ConstantPrimitiveEvolution' in fname:
        print 'skipped ConstantPrimitiveEvolution'
    elif 'RT.cpp' in fname:
        simple_compile(fname, lname)
    else:
        strict_compile(fname, lname)

def efficient_compile(fname, lname):

    import os
    import re

    oname = source2lib(fname,lname)
    if not os.path.isfile(oname):
        sensitive_compile(fname,lname)
    else:
        if os.path.getmtime(fname)>os.path.getmtime(oname):
            sensitive_compile(fname, lname)
        elif os.path.getmtime('./homebrew_makefile.py')>os.path.getmtime(oname):
            sensitive_compile(fname, lname)

        rel_dep_list = []
        for line in open(fname).readlines():
            if '#include' in line and not '<' in line:
                rel_dep_list.append(line)
        for dep_file in rel_dep_list:
            rel_path = dep_file.split('"')[-2]
            root_path = fname.replace(fname.split('/')[-1],'')
            header_path = os.path.normpath(os.path.join(root_path, rel_path))
            if os.path.getmtime(header_path)>os.path.getmtime(oname):
                sensitive_compile(fname, lname)
                break

def main():

    import os

    lname = 'library'
    safe_create(lname)

    for fname in get_source_files('source'):
        efficient_compile(fname,lname)

    # Archive
    os.system('ar cr ./'+lname+'/librich.a ./'+lname+'/*.o')

if __name__ == '__main__':

    main()

    
