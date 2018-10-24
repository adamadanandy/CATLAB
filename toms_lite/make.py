import subprocess
import sys
import os

gen_name = "toms_lite"

compile_defines_f = {
    'num'             : 0,
    'definevariables' : [],
    'definevalues'    : [],
}

compile_profile_f = {
    'compiler' : 'mpif90',
    'compile_option'   : [],
    'link_option' : [],
    'require' : [],
}

compile_objects_f = [
    'd1mach',
    'i1mach',
    'r1mach',
    'dlamch',
    'slamch',
    'lsame',
    'cacai',
    'cbesj',
    'cbuni',
    'cseri',
    'cunik',
    'gamln',
    'xzsqrt',
    'zbesi',
    'zbknu',
    'zmlri',
    'zuchk',
    'zunk2',
    'cacon',
    'cbesk',
    'cbunk',
    'cshch',
    'cunk1',
    'zacai',
    'zbesj',
    'zmlt',
    'zunhj',
    'zuoik',
    'cairy',
    'cbesy',
    'ckscl',
    'cuchk',
    'cunk2',
    'zacon',
    'zbesk',
    'zbuni',
    'zrati',
    'zuni1',
    'zwrsk',
    'casyi',
    'cbinu',
    'cmlri',
    'cunhj',
    'cuoik',
    'xzabs',
    'zairy',
    'zbesy',
    'zbunk',
    'zs1s2',
    'zuni2',
    'cbesh',
    'cbiry',
    'crati',
    'cuni1',
    'cwrsk',
    'xzexp',
    'zasyi',
    'zbinu',
    'zdiv',
    'zseri',
    'zunik',
    'cbesi',
    'cbknu',
    'cs1s2',
    'cuni2',
    'dgamln',
    'xzlog',
    'zbesh',
    'zbiry',
    'zkscl',
    'zshch',
    'zunk1',
    #'machcon',
    #'zbsubs',
]

generate_lib = {
    'cmd'       : 'ar',
    'opt'       : 'rv',
    'ranlib'    : 'ranlib',
}

def docompile(objects,profile,defines):
    cwdir = os.getcwd()
    for lib in profile['require']:
        libdir = cwdir+'/'+lib
        runcmd = 'python '+libdir+'/make.py'
        print(runcmd)
        subprocess.call(runcmd, cwd = libdir,shell=True)
    cc = profile['compiler']
    opt = ''
    for option in profile['compile_option']:
        opt = opt + ' ' + option
    variables = defines['definevariables']
    defvalues = defines['definevalues']
    dopt = ''
    for itr in range(defines['num']):
        dopt = dopt + ' -D'+variables[itr]+'='+defvalues[itr]
    for obj in objects:
        runcmd = cc+' '+obj+'.f'+' -c'+opt+dopt
        print (runcmd)
        subprocess.call(runcmd,shell=True)

def genlib(libname,objects,genlib,profile):
    libfilename = 'lib'+libname+'.a'
    runcmd = genlib['cmd']+' '+genlib['opt']+' '+libfilename+' '
    for obj in objects:
        runcmd += obj+'.o '
    for lib in profile['require']:
        runcmd += 'lib'+lib+'.a '
    print (runcmd)
    subprocess.call(runcmd,shell=True)
    runcmd = genlib['ranlib']+' '+libfilename
    print (runcmd)
    subprocess.call(runcmd,shell=True)
    runcmd = 'mv '+libfilename+' ../'
    subprocess.call(runcmd,shell=True)

def dolink(output,objects,profile):
    cc = profile['compiler']
    opt = ''
    for option in profile['link_option']:
        opt = opt + ' ' + option
    objs = ''
    for obj in objects:
        objs = objs+' '+obj+'.o'
    runcmd = cc+objs+opt+' -o '+output
    print (runcmd)
    subprocess.call(runcmd,shell=True)

if __name__ == '__main__':
    docompile(compile_objects_f,compile_profile_f,compile_defines_f)
    genlib(gen_name,compile_objects_f,generate_lib,compile_profile_f)
    subprocess.call('rm *.o',shell=True)

