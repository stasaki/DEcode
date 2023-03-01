# Import any necessary libraries/modules
import sys
import pandas as pd
from scipy import sparse

# Define any necessary functions or classes
def to_sparse(file_name):
    ds = pd.read_csv(file_name,sep="\t")
    ds['feature_vec']=""
    for i in range(ds.shape[0]):
        ds.values[i,5] =sparse.csr_matrix(([float(x) for x in ds.values[i,4].split(",")],
                           ([int(x) for x in ds.values[i,2].split(",")],
                            [int(x) for x in ds.values[i,3].split(",")])), 
                          shape=[int(x) for x in ds.values[i,1].split(",")])
    ds = ds[["Name", "feature_vec"]]
    ds.to_pickle(file_name.replace(".txt",".pkl"))

# Define the main function that will run the script
def main(file_name):
    to_sparse(file_name)

if __name__ == '__main__':
    # Check if the file_name argument was provided
    if len(sys.argv) == 2:
        file_name = sys.argv[1]
        main(file_name)
    else:
        print("Please provide a file name as an argument.")