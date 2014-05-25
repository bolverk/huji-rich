def execute_single(folder_path):

    """
    Executes a single simulation run
    Input:
    folder_path - Path to the run
    """

    import os
    home = os.getcwd()
    os.chdir(folder_path)
    os.system('./test.exe')
    os.system('python ./test.py')
    os.chdir(home)

def execute_all_runs(folder_list, pflag):

    """
    Executes all the simulations
    ptrn - Pattern of folder names
    pflag - Parallel flag (toggles using parallelisation)
    """

    if pflag:

        import multiprocessing
        pool = multiprocessing.Pool()
        pool.map(execute_single,folder_list)
    else:
        [execute_single(folder) for folder in folder_list]

def inalid_syntax_announcement():

    print 'Invalid syntax. Proper syntax:'
    print 'python ./execute_all_runs.py [folder list] [serial/parallel]'
    print 'for example:'
    print 'python ./execute_all_runs.py temp_acoustic_2d_res_* parallel'
    print 'will run all acoustic_2d runs in parallel'

if __name__=='__main__':
    
    import sys

    folder_list = sys.argv[1:-1]
    if sys.argv[-1]=='parallel':
        parallel_flag = True
    elif sys.argv[-1]=='serial':
        parallel_flag = False
    else:
        raise NameError('Unknown parallelisation option '+sys.argv[-1])

    execute_all_runs(folder_list,parallel_flag)
