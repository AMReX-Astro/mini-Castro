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
__global__ static void cuda_{}
{{
   int blo[3];
   int bhi[3];
   for (int k = lo[2] + blockIdx.z * blockDim.z + threadIdx.z; k <= hi[2]; k += blockDim.z * gridDim.z) {{
     blo[2] = k;
     bhi[2] = k;
     for (int j = lo[1] + blockIdx.y * blockDim.y + threadIdx.y; j <= hi[1]; j += blockDim.y * gridDim.y) {{
       blo[1] = j;
       bhi[1] = j;
       for (int i = lo[0] + blockIdx.x * blockDim.x + threadIdx.x; i <= hi[0]; i += blockDim.x * gridDim.x) {{
         blo[0] = i;
         bhi[0] = i;
         {};
       }}
     }}
   }}
}}
"""

# for finding a function signature that starts with DEVICE_LAUNCHABLE
sig_re = re.compile("(DEVICE_LAUNCHABLE)(\\()(.*)(\\))(;)", re.IGNORECASE|re.DOTALL)

# for finding just the variable definitions in the function signature (between the ())
decls_re = re.compile("(.*?)(\\()(.*)(\\))", re.IGNORECASE|re.DOTALL)

def doit(headers):

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

        # Now write out the CUDA kernels
        hout.write("\n")
        hout.write("#include <AMReX_BLFort.H>\n")
        hout.write("#include <AMReX_Device.H>\n")
        hout.write("\n")

        hdrmh = hdr.strip(".H")

        # Add an include guard
        hout.write("#ifndef _cuda_" + hdrmh + "_\n")
        hout.write("#define _cuda_" + hdrmh + "_\n\n")

        # Wrap the device declarations in extern "C"
        hout.write("#ifdef AMREX_USE_CUDA\n")
        hout.write("extern \"C\" {\n\n")

        line = hin.readline()
        while line:

            # if the line doesn't have DEVICE_LAUNCHABLE, then skip it.
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

                # First write out the device signature
                device_sig = "__device__ void {};\n\n".format(func_sig)
                hout.write(device_sig)

                # Now write out the global signature. This involves getting
                # rid of the data type definitions and also replacing the
                # lo and hi (which must be in the function definition) with blo and bhi.
                dd = decls_re.search(func_sig)
                vars = []

                has_lo = False
                has_hi = False

                for n, v in enumerate(dd.group(3).split(",")):

                    # we will assume that our function signatures _always_ include
                    # the name of the variable
                    _tmp = v.split()
                    var = _tmp[-1].replace("*","").replace("&","").strip()

                    # Replace AMReX Fortran macros
                    var = var.replace("BL_FORT_FAB_ARG_3D", "BL_FORT_FAB_VAL_3D")

                    if var == "lo":
                        var = "blo"
                        has_lo = True

                    if var == "hi":
                        var = "bhi"
                        has_hi = True

                    vars.append(var)

                if not has_lo or not has_hi:
                    sys.exit("ERROR: function signature must have variables lo and hi defined.")

                # reassemble the function sig
                all_vars = ", ".join(vars)
                new_call = "{}({})".format(dd.group(1), all_vars)

                hout.write(TEMPLATE.format(func_sig, new_call))
                hout.write("\n")

            line = hin.readline()

        # Close out the extern "C" region
        hout.write("\n}\n")
        hout.write("#endif\n")

        # Close out the include guard
        hout.write("\n")
        hout.write("#endif\n")

        hin.close()
        hout.close()

if __name__ == "__main__":
    HEADERS = ["Castro_F.H"]
    doit(HEADERS)
