import subprocess
import sys
import os

gen_name = "catlab"

compile_defines_c = {
    'num'             : 0,
    'definevariables' : [],
    'definevalues'    : [],
}

compile_profile_c = {
    'compiler' : 'mpicc',
    'compile_option'   : ['-I.','-I/home/yyli/include','-Wall'],
    'link_option' : [],
    'require' : ['toms_lite','cephes_lite'],
}

compile_objects_c = [
    'catlab_specfunc',
    'catlab_matrix',
    'catlab_eig',
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
    print (runcmd)
    subprocess.call(runcmd,shell=True)
    runcmd = genlib['ranlib']+' '+libfilename
    print (runcmd)
    subprocess.call(runcmd,shell=True)
    # if there are required libraries, we need to combine them into one lib
    if (profile['require']!=[]):
        for lib in profile['require']:
            thislibruncmd = 'ar x lib'+lib+'.a '
            print(thislibruncmd)
            subprocess.call(thislibruncmd,shell=True)
        runcmd = 'ar x '+libfilename
        print(runcmd)
        subprocess.call(runcmd,shell=True)
        runcmd = 'ar rv '+libfilename+' *.o'
        print(runcmd)
        subprocess.call(runcmd,shell=True)
        runcmd = 'ranlib '+libfilename
        print(runcmd)
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
    subprocess.call('rm *.a',shell=True)

