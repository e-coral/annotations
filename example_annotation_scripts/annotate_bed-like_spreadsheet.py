# maintain other excel fields when annotating boundaries in the format:
# chr A, pos A, chr B, pos B
# based on annotate_bed-like_boundaries, annotate_spreadsheet and annotate.annotate_standard_excel_input_file

import os
from pathlib import Path
from src.annotate import annotate
import pandas


# outdir = (Path(__file__).parent.parent / 'example_output_files').resolve()
outdir = Path('/home/emmon/OneDrive/Kate/GI2/results/')
indir = Path('/home/emmon/OneDrive/Kate/GI2/from_kate/')
# indir = (Path(__file__).parent.parent / 'example_input_files').resolve()

def annotate_spreadsheet_by_boundaries():
    """
    For all input spreadsheets, identify the boundaries and add required annotations
    :return: annotated spreadsheet
    """
    # for all input files
    for filename in indir.glob(r'*'):
        outname = f"{filename.stem}_tel-dist"
        infl = os.path.join(indir, filename)
        outfile = os.path.join(outdir, outname)

        with pandas.ExcelWriter(f"{outfile}.xlsx") as writer:
            # read in the input file
            infile = pandas.read_excel(infl, sheet_name=None)
            # for each sheet in the input file (probably only one for each in this case)
            for sname, s in infile.items():
                # if the sheet is not empty
                if not s.empty:
                    # store the original column names to be able to add them back in later
                    orig_columns = [col.strip() for col in s.columns]
                    print(filename)
                    print(orig_columns)

                    # s['end'] = s['posB'] if

                    if s['chrA'] == s['chrB']:
                        print(f"{filename} only has single-crhomosome events")
                    else:
                        print(f"{filename} has some different chrs")
                else:
                    print(f"{filename} {sname} is empty")



if __name__ == '__main__':
    annotate_spreadsheet_by_boundaries()

