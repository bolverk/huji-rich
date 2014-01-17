def make_temp_folder(pref,uname):

    import re

    res = pref+'_'+uname
    res = res.replace('.','')
    res = res.replace('/','_')
    res = res.replace('__','_')
    res = re.sub(r'_$','',res)
    return res

def invalid_syntax_announcement():
    print 'Invalid syntax'
    print 'proper syntax:'
    print 'python ./convergence_curve.py [test folder] [low res] [high res] [res step] [output file]'
    print 'exapmple use:'
    print 'python ./convergence_curve.py ./acoustic_2d/base/ 20 101 10 l1_vs_res.txt'

def main():

    import argparse
    import os

    parser = argparse.ArgumentParser(description='Generates convergence curves')
    parser.add_argument('test_folder',
                        help='path to test folder')
    parser.add_argument('lowest_res',
                        type=int,
                        help='lowest bound on res (power of 2)')
    parser.add_argument('highest_res',
                        type=int,
                        help='upper bound on res (power of 2)')
    parser.add_argument('output_file',
                        help='file to write the data to')
    args = parser.parse_args()

    os.system('python ./new_run_tests.py '+args.test_folder)
    exec_folder = make_temp_folder('temp',args.test_folder)

    os.system('python ./prep_dif_res_runs.py '+exec_folder+' '+
              str(args.lowest_res)+' '+str(args.highest_res))

    os.system('python ./execute_all_runs.py '+exec_folder+'_res_* parallel')

    os.system('python ./collect_data.py '+exec_folder+'_res_* '+args.output_file)
    
if __name__=='__main__':
    
    main()
