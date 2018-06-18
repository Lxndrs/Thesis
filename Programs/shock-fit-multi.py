"""Program to fit the shock shapes multiple times, with random
missing points so that we can estimate the uncertainties in the
derived radii 

"""
import numpy as np
import os

BASEDIR = 'Multi-Fit'
FRAC_KEEP = 0.6666
DEBUG = True

def find_source_in_line(line):
    return line.split('{')[-1][:-2]


def random_selection_from_regions(regfile):
    """Read a region file and return a list of the lines but with only a
fraction of the points for each source

    """
    with open(regfile) as f:
        reg_lines = f.readlines()

    shock_lines = [line for line in reg_lines if "shock outer" in line]
    other_lines = [line for line in reg_lines if "shock outer" not in line]
    sources = set([find_source_in_line(line) for line in shock_lines])

    source_dict = {source: [] for source in sources}
    for line in shock_lines:
        source_dict[find_source_in_line(line)].append(line) 

    keep_lines = []
    for source, source_lines in source_dict.items():
        npoints = len(source_lines)
        # Number to keep - fraction FRAC_KEEP but at least 4
        nkeep = max(4, int(FRAC_KEEP*npoints))
        if DEBUG:
            print(source, 'contains', npoints, 'points.  Keeping', nkeep, 'of them')
        np.random.shuffle(source_lines)
        sample_lines = source_lines[:nkeep]
        keep_lines.extend(sample_lines)

    return other_lines + keep_lines


if __name__ == '__main__':
    sample_filename = 'LV-positions-2018-will.reg'
    nsamples = 10
    for isample in range(nsamples):
        sample_id = 'samp{:02d}'.format(isample)
        sample_folder = BASEDIR + '/' + sample_id
        new_lines = random_selection_from_regions(sample_filename)
        new_sample_path = sample_folder + '/' + sample_filename.replace('.reg', '-{:s}.reg'.format(sample_id))
        if not os.path.isdir(sample_folder):
            os.makedirs(sample_folder)
        with open(new_sample_path, 'w') as f:
            f.write(''.join(new_lines))
        if DEBUG:
            print('New region written to', new_sample_path)

            

    
