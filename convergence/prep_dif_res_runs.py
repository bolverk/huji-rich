def prep_dif_res_runs(path_to_exec,res_list,res_file_name='resolution.txt'):
    """
    Prepares a series of the same simulation run with different resolutions
    It assumes the resolution is read off from a text file "resolutions.txt"
    Input:
    path_to_exec - Path to folder with executable
    res_list - List of resolutions to be scanned
    res_file_name - Name of the file were the resolutions is to be written
    """

    if '/' in path_to_exec:
        raise NameError('No slashes expected in path_to_exec')
    
    import os
    for res in res_list:
        folder_name = path_to_exec+'_res_'+str(res)
        os.system('mkdir '+folder_name)
        os.system('cp '+path_to_exec+'/test.exe '+folder_name)
        os.system('cp '+path_to_exec+'/test.py '+folder_name)
        os.system('cp '+path_to_exec+'/*.txt '+folder_name)
        
        f = open(folder_name+'/'+res_file_name,'w')
        f.write(str(res)+'\n')
        f.close()

def check_valid_input(smallest_resolusion, largest_resolution):

    if largest_resolution<=smallest_resolusion:
        raise NameError('smallest resolution is larger than largest resolusion')
    if largest_resolution >= 1e5:
        raise NameError('Resolution is too high')

def main():

    import argparse
    import numpy
    import re

    parser = argparse.ArgumentParser(description='Creates simulation folders of the same calculation with different resolutions')
    parser.add_argument('test_folder',
                        help='path to test folder')
    parser.add_argument('lowest_res',
                        type=int,
                        help='lowest res (power of two)')
    parser.add_argument('highest_res',
                        type=int,
                        help='highest res (power of two)')
    args = parser.parse_args()

    check_valid_input(args.lowest_res, args.highest_res)

    resolution_list = [2**x for x in range(args.lowest_res,
                                           args.highest_res+1)]

    folder_name = re.sub(r'/$','',args.test_folder)
    if '/' in folder_name:
        raise NameError('Expected folder name to be in the same level as the script')

    prep_dif_res_runs(folder_name, resolution_list)

if __name__=='__main__':

    main()
