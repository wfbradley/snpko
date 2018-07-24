import numpy as np
import pandas as pd

def syn_df(df, num_rows):

    df_fake = pd.DataFrame(columns=df.columns, index=np.arange(num_rows))
    for field in df_fake.columns:
        fake = np.random.choice(df[field].values, size=num_rows)
        df_fake[field].values[:] = fake

    return(df_fake)


def synthesize(real_data_file, fake_data_file, num_rows):
    '''
    Given a (correctly formatted) input file, synthesize a
    new (CSV) data file by independently scrambling the columns.
    This allows us to produce an anonymized file from real
    data, using the same file format.
    '''

    if real_data_file.endswith('.xlsx'):
        df = pd.read_excel(real_data_file)
    else:
        df = pd.read_csv(real_data_file)

    if num_rows == None:
        num_rows = len(df)
    df_fake = syn_df(df, num_rows)
    if fake_data_file:
        df_fake.to_csv(fake_data_file, index=False)
    else:
        s=df_fake.to_csv(fake_data_file, index=False)
        print (s)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Generate fake data from real data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', type=str, required=True,
                        help='Input file of real data')
    parser.add_argument('--output', type=str,
                        help='Output for synthetic data')
    parser.add_argument('--num_rows', type=int,
                        help='Number of rows; default=same as input')
    args = parser.parse_args()
    synthesize(args.input, args.output, args.num_rows)
