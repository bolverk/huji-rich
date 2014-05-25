def collect_data(folder_list):

    unsorted_lists = {'resolution':[],
                      'density l1':[],
                      'pressure l1':[],
                      'velocity l1':[]}

    import numpy
    for folder in folder_list:
        res = int(numpy.loadtxt(folder+'/resolution.txt'))
        rawd = numpy.loadtxt(folder+'/gradesheet.txt')
        unsorted_lists['resolution'].append(res)
        unsorted_lists['density l1'].append(rawd[0])
        unsorted_lists['pressure l1'].append(rawd[1])
        unsorted_lists['velocity l1'].append(rawd[2])

    mask = sorted(range(len(unsorted_lists['resolution'])),
                  key=lambda i:unsorted_lists['resolution'][i])
    sorted_lists = {}
    for field in unsorted_lists:
        sorted_lists[field] = [unsorted_lists[field][i] for i in mask]
    
    return sorted_lists

if __name__=='__main__':
    
    import sys

    folder_list = sys.argv[1:-1]
    output_file = sys.argv[-1]
    if not ".txt" in output_file:
        raise NameError("Expected .txt extension in output file. "+
                        "Was an output file provided?")
    temp = collect_data(folder_list)
    
    f = open(output_file,'w')
    for i in range(len(temp['resolution'])):
        for field in ['resolution','density l1','pressure l1','velocity l1']:
            f.write(str(temp[field][i])+' ')
        f.write('\n')
    f.close()
