import os
import subprocess

# c++ gsl include paths. Only add if gsl includes are not located in the search path of gcc
gsl_include_path = None # Set to None if gsl_include_path is already visible to gcc 
gsl_lib_path     = None # Set to None if gsl_lib_path is already visible to gcc
# gsl_include_path = "/Users/ewaldvandyk/bin/gsl/include" 
# gsl_lib_path = "/Users/ewaldvandyk/bin/gsl/lib"

gcc_bin_name = "gcc"


# Path directories
this_path = os.path.realpath(__file__)
this_path = os.path.dirname(this_path)
data_transform_path = os.path.abspath(os.path.join(this_path, "data_transform"))

# Files
gmm_c_file  = os.path.join(data_transform_path, 'gmm.c')
gmm_o_file  = os.path.join(data_transform_path, 'gmm.o')
gmm_so_file = os.path.join(data_transform_path, 'gmm.so')

if gsl_include_path == None:
    compile_args = [gcc_bin_name, "-c", "-fPIC", gmm_c_file, "-o", gmm_o_file]
else:
    compile_args = [gcc_bin_name, "-c", "-fPIC", "-I", gsl_include_path, gmm_c_file, "-o", gmm_o_file]

if gsl_lib_path == None:
    lib_gsl_arg = "-lgsl"
    lib_cblas_arg = "-lgslcblas"
else:
    lib_gsl_arg     = os.path.join(gsl_lib_path, "libgsl.a")
    lib_cblas_arg   = os.path.join(gsl_lib_path, "libgslcblas.a")
link_args = [gcc_bin_name, "-shared", gmm_o_file, lib_gsl_arg, lib_cblas_arg, "-o", gmm_so_file]

proc = subprocess.Popen(compile_args)
proc.wait()
proc = subprocess.Popen(link_args)
proc.wait()


print gsl_include_path


# gcc -c -fPIC c_ctypes_posterior.c
# gcc -I /usr/gsl/lib -shared c_ctypes_posterior.o  -lgsl -lgslcblas -o c_ctypes_posterior.so

# gcc -c -fPIC -I/Users/ewaldvandyk/bin/gsl/include/ gmm.c -o gmm.o
# gcc -shared gmm.o /Users/ewaldvandyk/bin/gsl/lib/libgsl.a /Users/ewaldvandyk/bin/gsl/lib/libgslcblas.a -o gmm.so