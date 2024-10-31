import pandas as pd
from rdkit.Chem import PandasTools
import rdkit as rdkit
from rdkit import DataStructs
import numpy as np


#generate SDF file
data = pd.read_csv('~/GHDDI_DataScience_lib/HTS_20230317/forGeneralUse/HTS_forGeneralUse_446664.csv', index_col=0) 
PandasTools.AddMoleculeColumnToFrame(data,'Cleaned_SMILES','Molecule',includeFingerprints=True)
PandasTools.WriteSDF(data, '~/HTS_forGeneralUse_446664.sdf', molColName='Molecule', properties=list(data.columns))

# Calculate molecular weight
molwt = [rdkit.Chem.rdMolDescriptors.CalcExactMolWt(rdkit.Chem.MolFromSmiles(data['SMILES'].iloc[i])) for i in range(len(data))]





