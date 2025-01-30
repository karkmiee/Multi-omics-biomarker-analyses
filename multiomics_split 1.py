import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.impute import KNNImputer

# Split 5-times into train and test sets
# Save both imputed and non-imputed splits

train_set_list = []
test_set_list = []
train_set_nan_list = []
test_set_nan_list = []


multiomics = pd.read_table("multi_omics_Filtered2.txt")

random_state_list = [1,2,3,4,5]

for i, rs in enumerate(random_state_list):
    test_set, train_set = train_test_split(multiomics,test_size=0.5,stratify=multiomics["Cancer_during_surveillance"], random_state=rs)
    train_nan, test_nan = train_set.copy(), test_set.copy()
    knn_imputer = KNNImputer()
    train_set = knn_imputer.fit_transform(train_set)
    test_set = knn_imputer.transform(test_set)
    train_set = pd.DataFrame(train_set,columns=multiomics.columns)
    test_set = pd.DataFrame(test_set,columns=multiomics.columns)

    train_set_list.append(train_set)
    test_set_list.append(test_set)

    train_set_nan_list.append(train_nan)
    test_set_nan_list.append(test_nan)

    train_set.to_csv(f"train_set_{i}.csv")
    test_set.to_csv(f"test_set_{i}.csv")