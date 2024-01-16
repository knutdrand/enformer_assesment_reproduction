import pandas as pd


def convert_bed_to_csv(bed_filename, csv_filename='tmp.csv'):
    intervals = bnp.read(bed_filename, buffer_type=bnp.io.Bed6Buffer)
    df = intervals.topandas()
    new_df = pd.DataFrame(
        {'ensg': df['name'],
         'chr': [s[3:] for s in df['chromosome']],
         'winS': df['start'],
         'winE': df['stop']})
    new_df.to_csv(csv_filename, index=False, header=False, sep='\t')
