#!/usr/bin/env python3

import sys


def main(input_path, output_path):
    with open(input_path) as input_handle, open(output_path, "w") as output_handle:
        gc_count = 0
        sequence_length = 0
        seen_record = False

        for line in input_handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if seen_record:
                    output_handle.write(f"{gc_count / sequence_length:.3f}\n")
                gc_count = 0
                sequence_length = 0
                seen_record = True
            else:
                gc_count += sum(base in "GgCc" for base in line)
                sequence_length += len(line)

        if seen_record:
            output_handle.write(f"{gc_count / sequence_length:.3f}\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
