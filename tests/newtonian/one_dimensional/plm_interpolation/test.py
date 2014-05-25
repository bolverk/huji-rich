#! /usr/bin/python

def main():

    import numpy

    np_list = numpy.loadtxt('np_list.txt')
    cc_pcm = numpy.loadtxt('cc_pcm.txt')
    cc_plm = numpy.loadtxt('cc_plm.txt')

    graphic_flag = False
    if graphic_flag:

        import pylab

        rawd_pcm = numpy.loadtxt('pcm_20_profile.txt')
        rawd_plm = numpy.loadtxt('plm_20_profile.txt')
        pylab.plot(rawd_pcm[:,0],rawd_pcm[:,1],
                   rawd_pcm[:,0],rawd_pcm[:,2],
                   rawd_plm[:,0],rawd_plm[:,2])
        pylab.show()

        pylab.loglog(np_list,cc_pcm,label='PCM')
        pylab.loglog(np_list,cc_plm,label='PLM')
        pylab.show()

    pcm_power = numpy.polyfit(numpy.log(np_list),numpy.log(cc_pcm),1)[0]
    plm_power = numpy.polyfit(numpy.log(np_list),numpy.log(cc_plm),1)[0]

    f = open('gradesheet.txt','w')
    f.write(str(pcm_power)+'\n'+
            str(plm_power)+'\n')
    f.close()

    return abs(pcm_power-(-1))<0.01 and \
        abs(plm_power-(-2))<0.01

if __name__=='__main__':
    import os

    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

