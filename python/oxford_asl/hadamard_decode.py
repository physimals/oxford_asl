#!/usr/bin/env python
"""
Does naive decoding of Hadamard time-encoded data

Assumes that the Hadamard blocks are in adjacent volumes in the data
and simply decodes each block, outputting the same order
"""
import sys
import argparse

from fsl.data.image import Image

import numpy as np
import scipy.linalg

class ArgumentParser(argparse.ArgumentParser):
    def __init__(self, **kwargs):
        argparse.ArgumentParser.__init__(self, prog="oxford_asl_hadamard_decode", add_help=False, **kwargs)
        self.add_argument("--input", required=True, help="Input file name")
        self.add_argument("--output", required=True, help="Output file name")
        self.add_argument("--hadamard-size", type=int, help="Hadamard matrix size")

def main():
    options = ArgumentParser().parse_args()
    if options.input is None:
        sys.stderr.write("Input file name must be specified")
        sys.exit(1)
    if options.output is None:
        sys.stderr.write("Output file name must be specified")
        sys.exit(1)

    print("Hadamard decoding\n")
    print(" - Input image: %s (assuming Hadamard cycles are in adjacent volumes)" % options.input)
    print(" - Hadamard matrix size: %i" % options.hadamard_size)
    print(" - Output image: %s" % options.output)

    input_img = Image(options.input)
    had_matrix = scipy.linalg.hadamard(options.hadamard_size)
    shape = input_img.shape
    if len(shape) != 4:
        sys.stderr.write("Input image must be 4D")
        sys.exit(1)
    elif shape[3] % options.hadamard_size != 0:
        sys.stderr.write("Input image has %i volumes, inconsistent with Hadamard matrix size of %i" % (shape[3], options.hadamard_size))
        sys.exit(1)

    nvols = shape[3]
    had_size = options.hadamard_size
    input_data = input_img.data
    output_data = np.zeros(list(shape[:3]) + [(had_size-1)*nvols//had_size])
    for vol in range(shape[3]):
        rpt_idx = int(vol / had_size)
        had_idx = vol % had_size
        for sub_bolus in range(had_size-1):
            output_data[..., rpt_idx*(had_size-1)+sub_bolus] += had_matrix[had_idx, sub_bolus+1]*input_data[..., vol]

    output_img = Image(output_data, header=input_img.header)
    output_img.save(options.output)

    print("\nDONE - Output in %s" % options.output)

if __name__ == "__main__":
    main()
