def main():

    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument('source_path',
                        help='path to file with main function')
    args = parser.parse_args()
    
    # Build
    hdf5_lib_path = os.environ['HDF5_LIB_PATH']
    rich_root = os.environ['RICH_ROOT']
    os.system('g++ -o ./rich.exe '+
              '-O3 '+
              args.source_path+' '+
              '-I '+rich_root+' '+
              (rich_root+'/library/*.o ').replace('//','/')+' '+
              (hdf5_lib_path+'/libhdf5_cpp.a').replace('//','/')+' '+
              (hdf5_lib_path+'/libhdf5_hl.a').replace('//','/')+' '+
              (hdf5_lib_path+'/libhdf5.a').replace('//','/')+' '+
              '-lpthread -lz -ldl ')

if __name__ == '__main__':

    main()
