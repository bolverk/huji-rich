def main():

    import os
    if os.path.isfile('serial_ignore.txt'):
        return True

    return False

if __name__ == '__main__':

    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
