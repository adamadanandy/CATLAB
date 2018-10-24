import subprocess
import sys
import os

gen_name = "cephes_lite"

compile_defines_c = {
    'num'             : 0,
    'definevariables' : [],
    'definevalues'    : [],
}

compile_profile_c = {
    'compiler' : 'mpicc',
    'compile_option'   : ['-I.','-Wall'],
    'link_option' : [],
    'require' : [],
}

compile_objects_c = [
    'airy',
    'chbevl',
    'const',
    'gamma',
    'hyperg',
    'i0',
    'i1',
    'iv',
    'j0',
    'j1',
    'jn',
    'jv',
    'mtherr',
    'polevl',
    'unity',
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
        runcmd = cc+' '+obj+'.c'+' -c'+opt+dopt
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
    docompile(compile_objects_c,compile_profile_c,compile_defines_c)
    genlib(gen_name,compile_objects_c,generate_lib,compile_profile_c)
    subprocess.call('rm *.o',shell=True)

