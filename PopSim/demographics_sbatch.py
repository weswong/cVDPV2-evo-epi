from models.python.PythonModelManager import PythonModelManager
from simtools.AssetManager.FileList import FileList
from simtools.SetupParser import SetupParser
import argparse

user_model = "Assets\\Demographics_Comps.py"
asset_collection_id = None
asset_files = FileList(root="Assets", recursive=True, max_depth=10)

parser = argparse.ArgumentParser(description = 'Combined Demographics script, used to setup all demography and projection files')
parser.add_argument('-f', default = 1,  action='store_true', 
    help="Simulate demographic transition in Matlab, Bangladesh. Set to 1 for yes, set to 0 for a high fertility simulation")
parser.add_argument('-b', default = 0,  action='store_true', 
    help="Simulate sustained S2, X-reactive immuntiy from a bOPV that extends for the projection simulation")
args = parser.parse_args()
print(args)

configs = []
for n in range(20):
    configs.append({'sim_iter': n+1,
                    'fertility_decline': args.f,
                    'strain_type': 'S2', #strain type for outbreak individuals, only S2 at this point
                    'bOPV': args.b}) #continue bOPV treatment up to 5 years later, boolean, 0 = False DO NOT USE FOR S3
exp_name = "Matlab Demographics:fertility_decline:{f}, outbreak_strain: S2, bOPV_Ximmunity: {b}".format(f = args.f, b = args.b)

if __name__ == "__main__":
    SetupParser.default_block = 'HPC'
    SetupParser.init()

    pmm = PythonModelManager(user_model=user_model, asset_collection_id=asset_collection_id, asset_files=asset_files, exp_name=exp_name,configs=configs)
    pmm.execute(True)
