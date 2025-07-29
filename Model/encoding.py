# Packages we're using
import os
import numpy as np



def load_dataset(clean, noisy):
    X, Y = [], []

    clean_dir = os.listdir(clean)

    count = 0

    for file_name in  clean_dir:
        clean_file = os.path.join(clean, file_name)
        noisy_file = os.path.join(noisy, file_name)


        clean_f = np.load(clean_file)
        noisy_f = np.load(noisy_file)




        X.append(noisy_f)
        Y.append(clean_f)
        count +=1

        print("File no:", count, "done")

    X = np.array(X, dtype=np.float32)
    Y = np.array(Y, dtype=np.float32)

    return X, Y








