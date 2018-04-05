#!/usr/bin/env python3

# Search the code for function signatures wrapped in
#
#   void DEVICE_LAUNCHABLE(ca_func(const int* lo, ...));
#
# This should be expanded to two function signatures
#
#   void ca_func(const int* lo, ...);
#
#   __global__ void cuda_ca_func(const int* lo, ...);
#
# The cuda_ version is new
#
# We would then need to write the cuda_ca_func() function

import os
import re
import sys

TEMPLATE = """
__global__ cuda_{}
{{

   int blo[3];
   int bhi[3];
   get_loop_bounds(blo, bhi, lo, hi);
   {};
}}
"""

# for finding a function signature that starts with DEVICE_LAUNCHABLE
sig_re = re.compile("(DEVICE_LAUNCHABLE)(\\()(.*)(\\))(;)", re.IGNORECASE|re.DOTALL)

# for finding just the variable definitions in the function signature (between the ())
decls_re = re.compile("(.*?)(\\()(.*)(\\))", re.IGNORECASE|re.DOTALL)

def doit(headers, cuda_file):

    # keep track of the functions we need to CUDA-ize
    needs_cuda = []

    for hdr in headers:
        # open the header file
        try:
            hin = open(hdr, "r")
        except IOError:
            sys.exit("Cannot open header {}".format(hdr))

        # open the CUDA header for output
        head, tail = os.path.split(hdr)
        ofile = os.path.join(head, "cuda_" + tail)
        try:
            hout = open(ofile, "w")
        except IOError:
            sys.exit("Cannot open output file {}".format(ofile))

        line = hin.readline()
        while line:

            # if the line doesn't have DEVICE_LAUNCHABLE, then output it to the new header
            # otherwise, we need to capture the function signature
            if "DEVICE_LAUNCHABLE" in line:
                launch_sig = "" + line
                sig_end = False
                while not sig_end:
                    line = hin.readline()
                    launch_sig += line
                    if line.strip().endswith(";"):
                        sig_end = True

                # now get just the actual signature
                m = sig_re.search(launch_sig)
                func_sig = m.group(3)

                needs_cuda.append(func_sig)

                # duplicate the signature, once for device and once for global and write these out
                orig_sig = "void {};\n\n".format(func_sig)
                global_sig = "__global__ void cuda_{};\n\n".format(func_sig)

                hout.write(orig_sig)
                hout.write(global_sig)

            else:
                # we just need to output the line unmodified
                hout.write(line)

            line = hin.readline()

        # done with this header
        hin.close()
        hout.close()

    # we are done with all the headers now -- the only thing left is
    # to write out the CUDA version of the file
    with open(cuda_file, "w") as of:
        for sig in needs_cuda:

            # we need to call the cuda_ version of the routine, which
            # means getting rid of the data type definitions and also
            # replacing the _first_ lo and hi with blo and bhi
            dd = decls_re.search(sig)
            vars = []
            for n, v in enumerate(dd.group(3).split(",")):

                # we will assume that our function signatures _always_ include
                # the name of the variable
                _tmp = v.split()
                var = _tmp[-1].replace("*","").replace("&","").strip()

                if n == 0:
                    if var == "lo":
                        var = "blo"
                    else:
                        sys.exit("ERROR: function signatures need to start with lo")

                if n == 1:
                    if var == "hi":
                        var = "bhi"
                    else:
                        sys.exit("ERROR: function signatures need hi as the second argument")

                vars.append(var)

            # reassemble the function sig
            all_vars = ", ".join(vars)
            new_call = "{}({})".format(dd.group(1), all_vars)

            of.write(TEMPLATE.format(sig, new_call))

if __name__ == "__main__":
    HEADERS = ["Castro_F.H"]
    CUDA_FILE = "cuda_interfaces.cpp"
    doit(HEADERS, CUDA_FILE)
