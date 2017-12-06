import os

# c++ gsl include paths. Only add if gsl includes are not located in the search path of gcc
gsl_include_path = None 
gsl_lib_path     = None
# gsl_include_path = os.path.join("/", "usr", "local", "include") #Default for OSX installation
# gsl_include_path = os.path.join("/", "usr", "include") #Default for Linux installation
# gsl_lib_path = os.path.join("/", "usr", "local", "lib") #Default for OSX


# Path directories
this_path = os.path.realpath(__file__)
this_path = os.path.dirname(this_path)
data_transform_path = os.path.abspath(os.path.join(this_path, "data_transform"))

# Files
gmm_c_file  = os.path.join(data_transform_path, 'gmm.c')
gmm_o_file  = os.path.join(data_transform_path, 'gmm.o')
gmm_so_file = os.path.join(data_transform_path, 'gmm.so')



print gsl_include_path


# gcc -c -fPIC c_ctypes_posterior.c
# gcc -I /usr/gsl/lib -shared c_ctypes_posterior.o  -lgsl -lgslcblas -o c_ctypes_posterior.so

# gcc -c -fPIC -I/Users/ewaldvandyk/bin/gsl/include/ gmm.c -o gmm.o
# gcc -shared gmm.o /Users/ewaldvandyk/bin/gsl/lib/libgsl.a /Users/ewaldvandyk/bin/gsl/lib/libgslcblas.a -o gmm.so