def get_file_name_from_user():

    from argparse import ArgumentParser

    parser = ArgumentParser(description='replaces old style cast with new style cast')
    parser.add_argument('fname',help='file name')
    return parser.parse_args().fname

def main(fname):

    import re

    raw_data = open(fname).readlines()
    f = open(fname,'w')
    for line in raw_data:
        new_line = re.sub(r'([^a-zA-Z0-9_>])\(([a-zA-Z0-9_]+)\)([a-zA-Z0-9_.]+)\(\);',
                          r'\1static_cast<\2>(\3());',line)
        new_line = re.sub(r'([^a-zA-Z0-9_>])\(([a-zA-Z0-9_]+)\)([a-zA-Z0-9_.]+)\(\),',
                          r'\1static_cast<\2>(\3()),',line)
        new_line = re.sub(r'([^a-zA-Z0-9_>])\(([a-zA-Z0-9_]+)\)([a-zA-Z0-9_.]+)\(\)\)',
                          r'\1static_cast<\2>(\3()))',line)
        new_line = re.sub(r'([^a-zA-Z0-9_>])\(([a-zA-Z0-9_]+)\)([a-zA-Z0-9_.]+);',
                          r'\1static_cast<\2>(\3);',new_line)
        new_line = re.sub(r'([^a-zA-Z0-9_>])\(([a-zA-Z0-9_]+)\)([a-zA-Z0-9_.]+)\)',
                          r'\1static_cast<\2>(\3))',new_line)
        new_line = re.sub(r'([^a-zA-Z0-9_>])\(([a-zA-Z0-9_]+)\)([a-zA-Z0-9_.]+)\]',
                          r'\1static_cast<\2>(\3)]',new_line)
        new_line = re.sub(r'([^a-zA-Z0-9_>])\(([a-zA-Z0-9_]+)\)([a-zA-Z0-9_.]+)-([^>])',
                          r'\1static_cast<\2>(\3)-\4',new_line)
        f.write(new_line)
    f.close()        

if __name__ == '__main__':

    main(get_file_name_from_user())
