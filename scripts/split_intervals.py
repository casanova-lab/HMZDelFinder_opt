#! /usr/bin/env python

import argparse
import math
import sys


def split_interval(interval, window, overlap):
    """Split a large interval into overlapping windows of fixed size.

    >>> split_interval([151, 300], 100, 50)
    [(151, 250), (201, 300)]
    >>> split_interval([101, 201], 100, 50)
    [(101, 200), (102, 201)]
    >>> split_interval([231, 530], 100, 50)
    [(231, 330), (281, 380), (331, 430), (381, 480), (431, 530)]
    >>> split_interval([1, 53], 100, 50)
    [(1, 53)]
    >>> split_interval([1, 200], 100, 0)
    [(1, 100), (101, 200)]
    >>> split_interval([1, 200], 100, 20)
    [(1, 100), (81, 180), (101, 200)]
    """

    length = interval[1] - interval[0] + 1
    n_windows = math.ceil((length - window) / (window - overlap))

    new_intervals = []
    for i in range(0, n_windows):
        new_intervals.append((interval[0] + i * (window - overlap),
                              interval[0] + i * (window - overlap) + window - 1))

    last_interval = (max(interval[1] - window + 1, interval[0]), interval[1])

    if last_interval not in new_intervals:
        new_intervals.append(last_interval)

    return new_intervals


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Split intervals into overlapping windows of fixed size.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', metavar="FILE", type=str, nargs=1,
                        help='Input file (BED format)')
    parser.add_argument('--window', dest='window', default=100, type=int,
                        help='Window size')
    parser.add_argument('--overlap', dest='overlap', default=50, type=int,
                        help='Overlap size')

    args = parser.parse_args()

    if args.window < 1:
        sys.stderr.write("ERROR: Size of sliding windows should be at least 1!")
        sys.exit(1)

    if args.overlap < 0 or args.overlap > args.window:
        sys.stderr.write("ERROR: Size of overlap between sliding windows should be >= 0 and < window_size!")
        sys.exit(1)


    with open(args.input[0]) as fh:
        content = fh.readlines()

    # Split intervals
    for line in content:
        chr, start, stop, *additional_info = line.rstrip().split()
        length = int(stop) - int(start) + 1

        intervals = split_interval((int(start), int(stop)), args.window, args.overlap)
        for interval in intervals:
            sys.stdout.write('\t'.join([chr, str(interval[0]), str(interval[1])] + additional_info) + '\n')
