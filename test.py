import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from scipy.sparse import csc_matrix

# Load RDS file using rpy2
readRDS = robjects.r['readRDS']
sparse_matrix_r = readRDS('saves/peak-calling-output/sparsePeakMatrix.rds')

# Convert R's dgCMatrix (sparse matrix) to Python's CSC (Compressed Sparse Column) format
with localconverter(robjects.default_converter + pandas2ri.converter):
    data = sparse_matrix_r.rx2('x')
    i = sparse_matrix_r.rx2('i')
    p = sparse_matrix_r.rx2('p')
    
    sparse_matrix_python = csc_matrix((data, i, p[:-1]))

print(sparse_matrix_python)
