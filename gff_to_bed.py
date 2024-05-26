import sys
import os

def gff_to_bed():
    requiredVersion = (3, 6)  # Ensure it's set to a Python 3 version now
    check_installation(requiredVersion)

    output_lines = []
    for line in sys.stdin:
        chomped_line = line.rstrip(os.linesep)
        if chomped_line.startswith('##'):
            continue
        else:
            elems = chomped_line.split('\t')
            cols = {
                'seqid': elems[0],
                'source': elems[1],
                'type': elems[2],
                'start': int(elems[3]) - 1,  # BED format is 0-based
                'end': int(elems[4]),
                'score': elems[5],
                'strand': elems[6],
                'phase': elems[7],
                'attributes': elems[8]
            }
            output_lines.append('\t'.join([cols['seqid'],
                                           str(cols['start']),
                                           str(cols['end']),
                                           cols['type'],
                                           cols['score'],
                                           cols['strand']]))
    return output_lines

def check_installation(rv):
    currentVersion = sys.version_info
    if currentVersion < rv:
        sys.stderr.write("[gff2bed.py] - Error: Your Python interpreter must be %d.%d or greater.\n" % (rv[0], rv[1]))
        sys.exit(-1)

if __name__ == '__main__':
    # Adapt to work as a script and as an importable module
    if sys.stdin.isatty():
        print("No data provided on standard input.", file=sys.stderr)
        sys.exit(1)
    else:
        for line in gff_to_bed():
            print(line)