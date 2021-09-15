import pandas as pd
import numpy as np

from .secondlayer import secondlayer

from pathlib import Path

DATA_PATH = str(Path(__file__).absolute().parent) + '/data'

with open(DATA_PATH+'/Database_2nd.txt',mode='r') as handle:
    db2 = pd.read_csv(handle,sep='\t',header=None, names=['0', '1', '2','charge'])
with open(DATA_PATH+'/Database_1st.txt',mode='r') as handle:
    db1 = pd.read_csv(handle,sep='\t',header=None, names=['0', '1','charge'])
with open(DATA_PATH+'/Database_0th.txt',mode='r') as handle:
    db0 = pd.read_csv(handle,sep='\t',header=None, names=['0','charge'])

def mCBAC(f):
    """
    Inputs:
    f = filename of P 1 CIF
       
    Outputs:
    charges_zeroed = numpy array of charges, following neutralization
       
    """
    # --------------------------------------------------
    # Run the python version of secondlayer.f90
    # Return connectivity as a pandas dataframe
    connectivity_df = secondlayer(f)
    # --------------------------------------------------
            
    # Cross-reference the m-CBAC databases
    
    # This section is ripe for a rewrite because I don't like nesting if/thens to this extent
    zero_charges = []
    for index, row in connectivity_df.iterrows():
        match2 = db2.loc[(db2['0'] == row['0']) & (db2['1'] == row['1']) & (db2['2'] == row['2'])]['charge'].values
        if len(match2) > 1:
            raise Exception('ERROR: multiple matches in database (level 2)')
        elif len(match2) == 1:
            #print(row['0'],row['1'],row['2'],match2[0])
            zero_charges.append(match2[0])
        else:
            match1 = db1.loc[(db1['0'] == row['0']) & (db1['1'] == row['1'])]['charge'].values
            if len(match1) > 1:
                raise Exception('ERROR: multiple matches in database (level 1)')
            elif len(match1) == 1:
                #print(row['0'],row['1'],row['2'],match1[0])
                zero_charges.append(match1[0])
            else:
                match0 = db0.loc[(db0['0'] == row['0'])]['charge'].values
                if len(match0) > 1:
                    raise Exception('ERROR: multiple matches in database (level 0)')
                elif len(match0) == 1:
                    #print(row['0'],row['1'],row['2'],match0[0])
                    zero_charges.append(match0[0])
                else:
                    #print(row['0'],row['1'],row['2'],row['0'])
                    zero_charges.append(row['0'])
    
    # Identify sum of charges and number of unknown charges
    Sum_charges = sum([x for x in zero_charges if isinstance(x,float)])
    Num_unknown = sum([1 for x in zero_charges if not isinstance(x,float)])
    print(Sum_charges)
    print(Num_unknown, 'unknown charges')

    if Num_unknown > 0:
        unknown_elements = []
        Unknown_charge = -Sum_charges/Num_unknown
        for i,charge in enumerate(zero_charges):
            if not isinstance(charge,float):
                unknown_elements.append(charge)
                zero_charges[i] = Unknown_charge
        # Add some diagnostics about the missing charges
        with open('unknown_element.list',mode='w') as handle:
            handle.write(f+'\n')
            handle.write(' ')
            handle.write(' '.join(list(set(unknown_elements)))+'\n')
            
    # Neutralize the total charge; weight the shift
    Abs_sum = sum([np.abs(x) for x in zero_charges])
    Real_sum = sum([x for x in zero_charges])
    zero_charges = np.array(zero_charges)
    charges_zeroed = zero_charges - Real_sum*np.abs(zero_charges)/Abs_sum

    return charges_zeroed
