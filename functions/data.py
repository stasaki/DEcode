def load_ml_data(deg_data_file,
                 mRNA_data_loc,
                 mRNA_annotation_data,
                 promoter_data_loc,
                 promoter_annotation_data):
    
    # deg_data_file: a full path to transcriptome data
    # mRNA_data_loc: a full path to a directory where mRNA annotation data is stored
    # mRNA_annotation_data: a file name of mRNA annotation data
    # promoter_data_loc:  a full path to a directory where promoter annotation data is stored
    # promoter_annotation_data: a file name of promoter annotation data
                    
    import pandas as pd
    import numpy as np
    from scipy.sparse import vstack
    
    # Subset DEG data
    import os.path
    root, ext = os.path.splitext(deg_data_file)
    if ext==".txt":
        deg_data = pd.read_csv(deg_data_file,sep="\t")

        if deg_data.columns[0] != 'Name':
            print('The first column must be Name')
            sys.exit(1) 

        med_exp = np.median(deg_data.values[:,1:],axis=1)
        for i in range(deg_data.shape[0]):
            deg_data.iloc[i,1:]  = deg_data.values[i,1:] - med_exp[i]
        deg_data['MedianExp'] = med_exp
    elif ext==".pkl":
        deg_data = pd.read_pickle(deg_data_file)
        if deg_data.columns[-1] != "MedianExp":
            print('The first column must be MedianExp')
            sys.exit(1) 
    else:
        sys.exit(1) 

    # Cocatenate mRNA data
    mRNA_data = pd.read_pickle(mRNA_data_loc+"/"+mRNA_annotation_data[0]+".pkl")
    mRNA_data = pd.merge(deg_data[['Name']], mRNA_data,
        how='inner', on='Name')
    mRNA_feature_name = pd.read_csv(mRNA_data_loc+"/"+mRNA_annotation_data[0]+"_feature_name.txt.gz",
                                      header=None,compression='gzip')
    mRNA_feature_name = mRNA_feature_name.values

    for mRNA_Annot in mRNA_annotation_data[1:]:
        range_data_temp = pd.read_pickle(mRNA_data_loc+"/"+mRNA_Annot+".pkl")
        
        # Filter genes
        range_data_temp = pd.merge(deg_data[['Name']], range_data_temp,
            how='inner', on='Name')
        # Read feature names
        feature_name_temp = pd.read_csv(mRNA_data_loc+"/"+mRNA_Annot+"_feature_name.txt.gz",
                                          header=None,compression='gzip')
        feature_name_temp = feature_name_temp.values

        # Check gene order is identical
        if sum(mRNA_data.Name != range_data_temp.Name) >0:
            raise Exception("gene name does not math!")

        # Remove exon
        indx=feature_name_temp[:,0]!="exon"
        for i in range(mRNA_data.shape[0]):
            mRNA_data.values[i,1]=vstack((mRNA_data.values[i,1],range_data_temp.values[i,1][indx,]))

        # Update feature name
        mRNA_feature_name=np.concatenate((mRNA_feature_name,feature_name_temp[indx]))

    # Convert to dense matrix
    for i in range(mRNA_data.shape[0]):
        mRNA_data.values[i,1]=np.array(mRNA_data.values[i,1].todense())
    
    # Filter features
    has_annot=map(lambda x: np.sum(x,axis=1), mRNA_data.values[:,1])
    has_annot=np.vstack(has_annot)
    indx = np.sum(has_annot>0,axis=0) < 30
    for i in range(mRNA_data.shape[0]):
        mRNA_data.values[i,1]=mRNA_data.values[i,1][~indx,:]
    mRNA_feature_name=mRNA_feature_name[~indx,0]


    # Cocatenate promoter data
    promoter_data = pd.read_pickle(promoter_data_loc+"/"+promoter_annotation_data[0]+".pkl")
    promoter_data = pd.merge(deg_data[['Name']], promoter_data,
        how='inner', on='Name')
    promoter_feature_name = pd.read_csv(promoter_data_loc+"/"+promoter_annotation_data[0]+"_feature_name.txt.gz",
                                          header=None,compression='gzip')
    promoter_feature_name = promoter_feature_name.values

    for mRNA_Annot in promoter_annotation_data[1:]:
        range_data_temp = pd.read_pickle(promoter_data_loc+"/"+mRNA_Annot+".pkl")
        
        # Filter genes
        range_data_temp = pd.merge(deg_data[['Name']], range_data_temp,
            how='inner', on='Name')
        # Read feature names
        feature_name_temp = pd.read_table(promoter_data_loc+"/"+mRNA_Annot+"_feature_name.txt.gz",
                                          header=None,compression='gzip')
        feature_name_temp = feature_name_temp.values

        # Check gene order is identical
        if sum(promoter_data.Name != range_data_temp.Name) >0:
            raise Exception("gene name does not math!")

        # Remove mask
        indx=feature_name_temp[:,0]!="promotor_annot_mask.rds"
        for i in range(promoter_data.shape[0]):
            promoter_data.values[i,1]=np.vstack((promoter_data.values[i,1],range_data_temp.values[i,1][indx,]))

        # update feature name
        promoter_feature_name=np.concatenate((promoter_feature_name,feature_name_temp[indx]))
    
    # Convert to dense matrix
    for i in range(promoter_data.shape[0]):
        promoter_data.values[i,1]=np.array(promoter_data.values[i,1].todense())
    
    # Filter features
    has_annot=map(lambda x: np.sum(x,axis=1), promoter_data.values[:,1])
    has_annot=np.vstack(has_annot)
    indx = np.sum(has_annot>0,axis=0) < 30
    for i in range(promoter_data.shape[0]):
        promoter_data.values[i,1]=promoter_data.values[i,1][~indx,:]
    promoter_feature_name=promoter_feature_name[~indx,0]
    

    # check all matched
    if (sum(deg_data.values[:,0] != mRNA_data.values[:,0]))>0:
        raise Exception("gene name does not math!")
    if (sum(deg_data.values[:,0] != promoter_data.values[:,0]))>0:
        raise Exception("gene name does not math!")

    return deg_data,mRNA_data,promoter_data, mRNA_feature_name, promoter_feature_name


def prep_ml_data_split(deg_data_file,
                       mRNA_data_loc,
                       mRNA_annotation_data,
                       promoter_data_loc,
                       promoter_annotation_data,
                       train_genes,
                       validate_genes,
                       test_genes,
                       outloc,
                       shuffle="None"):
    
    # deg_data_file: a full path to transcriptome data
    # mRNA_data_loc: a full path to a directory where mRNA annotation data is stored
    # mRNA_annotation_data: a file name of mRNA annotation data
    # promoter_data_loc:  a full path to a directory where promoter annotation data is stored
    # promoter_annotation_data: a file name of promoter annotation data
    # train_genes: a full path to a file contating gene ids for training
    # validate_genes: a full path to a file contating gene ids for validation
    # test_genes: a full path to a file contating gene ids for testing
    # outloc: a full path where scaling factors are saved
    # shuffle: which features are shuffled or not
    
    import pandas as pd
    import numpy as np
    
    deg_data, mRNA_data, promoter_data, mRNA_feature_name, promoter_feature_name=load_ml_data(deg_data_file,
                                                                                                         mRNA_data_loc,
                                                                                                         mRNA_annotation_data,
                                                                                                         promoter_data_loc,
                                                                                                         promoter_annotation_data)
    
    
    # Split data into training, validating, and testing subsets
    test = pd.read_csv(test_genes,sep="\t",header=None,compression='gzip').values[:,0]
    validate = pd.read_csv(validate_genes,sep="\t",header=None,compression='gzip').values[:,0]
    train = pd.read_csv(train_genes,sep="\t",header=None,compression='gzip').values[:,0]

    # Split mRNA feature
    X_mRNA_train=mRNA_data.query("Name in @train")
    X_mRNA_validate=mRNA_data.query("Name in @validate")
    X_mRNA_test=mRNA_data.query("Name in @test")

    # Split promoter feature
    X_promoter_train=promoter_data.query("Name in @train")
    X_promoter_validate=promoter_data.query("Name in @validate")
    X_promoter_test=promoter_data.query("Name in @test")


    # Split target data
    Y_train=deg_data.query("Name in @train")
    Y_validate=deg_data.query("Name in @validate")
    Y_test=deg_data.query("Name in @test")

    # Scale fold changes
    std=Y_train.values[:,1:(Y_train.shape[1]-1)].std()
    for i in range(1,Y_train.shape[1]-1):
        Y_train.iloc[:,i]=Y_train.values[:,i]/std
        Y_validate.iloc[:,i]=Y_validate.values[:,i]/std
        Y_test.iloc[:,i]=Y_test.values[:,i]/std
    
    # Store normalization stats 
    DEG_scale_stats = pd.DataFrame({'feature_name' : np.array(Y_train.columns)[1:(Y_train.shape[1]-1)],
                                    'row_indx': range(0,Y_train.shape[1]-2),
                                    'feature_type' : 'deg_stat',
                                    'std': std})
    
    # Scale log2-TPM (assume it is located at the last column)
    std=Y_train.iloc[:,(Y_train.shape[1]-1)].std()
    Y_train.iloc[:,(Y_train.shape[1]-1)]=Y_train.values[:,(Y_train.shape[1]-1)]/std
    Y_validate.iloc[:,(Y_train.shape[1]-1)]=Y_validate.values[:,(Y_train.shape[1]-1)]/std
    Y_test.iloc[:,(Y_train.shape[1]-1)]=Y_test.values[:,(Y_train.shape[1]-1)]/std

    # Store normalization stats 
    DEG_scale_stats=pd.concat((DEG_scale_stats,
                               pd.DataFrame({'feature_name' : np.array(Y_train.columns)[(Y_train.shape[1]-1):],
                                             'row_indx': Y_train.shape[1]-2,
                                             'feature_type' : 'deg_stat',
                                             'std': std})))  

    # Normalize mRNA features
    std = np.hstack(X_mRNA_train.values[:,1]).max(axis=1)

    # Special treatment for STD=0
    std[std==0]=1
    
    for i in range(len(X_mRNA_train.values[:,1])):
        X_mRNA_train.values[i,1]=((X_mRNA_train.values[i,1].transpose())/std).transpose()
    for i in range(len(X_mRNA_validate.values[:,1])):
        X_mRNA_validate.values[i,1]=((X_mRNA_validate.values[i,1].transpose())/std).transpose()
    for i in range(len(X_mRNA_test.values[:,1])):
        X_mRNA_test.values[i,1]=((X_mRNA_test.values[i,1].transpose())/std).transpose()

    # Store normalization stats 
    mRNA_feature_norm_stats = pd.DataFrame({'feature_name' : mRNA_feature_name,
                                            'row_indx': range(len(mRNA_feature_name)),
                                            'feature_type' : 'mRNA_range',
                                            'max':std})

    # Normalizing promoter features
    std = np.hstack(X_promoter_train.values[:,1]).max(axis=1)

    # Special treatment for STD=0
    std[std==0]=1

    # Scaling
    for i in range(len(X_promoter_train.values[:,1])):
        X_promoter_train.values[i,1]=((X_promoter_train.values[i,1].transpose())/std).transpose()
    for i in range(len(X_promoter_validate.values[:,1])):
        X_promoter_validate.values[i,1]=((X_promoter_validate.values[i,1].transpose())/std).transpose()
    for i in range(len(X_promoter_test.values[:,1])):
        X_promoter_test.values[i,1]=((X_promoter_test.values[i,1].transpose())/std).transpose()

    # Store normalization stats 
    promoter_feature_norm_stats = pd.DataFrame({'feature_name' : promoter_feature_name,
                                                'row_indx': range(len(promoter_feature_name)),
                                                'feature_type' : 'promoter_range',
                                                'max':std})

    # conbine all stats 
    feature_norm_stats=pd.concat((DEG_scale_stats,
                                  mRNA_feature_norm_stats,
                                  promoter_feature_norm_stats),sort=False)  

    feature_norm_stats.to_csv(outloc+'feature_norm_stats.txt',sep=",")
    
    # Shuffle features
    if shuffle=="all":
        np.random.seed(1234)
        shuffle_indices = np.random.permutation(np.arange(X_mRNA_train.shape[0]))
        X_mRNA_train.values[:,1] =X_mRNA_train.values[shuffle_indices,1]
        shuffle_indices = np.random.permutation(np.arange(X_mRNA_validate.shape[0]))
        X_mRNA_validate.values[:,1] =X_mRNA_validate.values[shuffle_indices,1]
        shuffle_indices = np.random.permutation(np.arange(X_mRNA_test.shape[0]))
        X_mRNA_test.values[:,1] =X_mRNA_test.values[shuffle_indices,1]

        shuffle_indices = np.random.permutation(np.arange(X_promoter_train.shape[0]))
        X_promoter_train.values[:,1] =X_promoter_train.values[shuffle_indices,1]
        shuffle_indices = np.random.permutation(np.arange(X_promoter_validate.shape[0]))
        X_promoter_validate.values[:,1] =X_promoter_validate.values[shuffle_indices,1]
        shuffle_indices = np.random.permutation(np.arange(X_promoter_test.shape[0]))
        X_promoter_test.values[:,1] =X_promoter_test.values[shuffle_indices,1]
    elif shuffle=="DNA":
        np.random.seed(1234)
        shuffle_indices = np.random.permutation(np.arange(X_promoter_train.shape[0]))
        X_promoter_train.values[:,1] =X_promoter_train.values[shuffle_indices,1]
        shuffle_indices = np.random.permutation(np.arange(X_promoter_validate.shape[0]))
        X_promoter_validate.values[:,1] =X_promoter_validate.values[shuffle_indices,1]
        shuffle_indices = np.random.permutation(np.arange(X_promoter_test.shape[0]))
        X_promoter_test.values[:,1] =X_promoter_test.values[shuffle_indices,1]
    elif shuffle=="RNA":
        np.random.seed(1234)
        shuffle_indices = np.random.permutation(np.arange(X_mRNA_train.shape[0]))
        X_mRNA_train.values[:,1] =X_mRNA_train.values[shuffle_indices,1]
        shuffle_indices = np.random.permutation(np.arange(X_mRNA_validate.shape[0]))
        X_mRNA_validate.values[:,1] =X_mRNA_validate.values[shuffle_indices,1]
        shuffle_indices = np.random.permutation(np.arange(X_mRNA_test.shape[0]))
        X_mRNA_test.values[:,1] =X_mRNA_test.values[shuffle_indices,1]

    return Y_train, Y_validate, Y_test, X_mRNA_train, X_mRNA_validate, X_mRNA_test, X_promoter_train, X_promoter_validate, X_promoter_test

def batch_iter(X_mRNA_train,
               X_promoter_train,
               Y_train,
               batch_size,
               shuffle=True):
    
    import numpy as np
    
    # X_mRNA:_train mRNA features
    # X_promoter_train: promoter features
    # Y_train: target data
    # batch_size: The number of data in each batch
    # shuffle=True: Shuffling the order of data in each epoch
    
    # Number of data in each batch
    data_size = len(Y_train)
    num_batches_per_epoch = int((data_size - 1) / batch_size) + 1
    
    # Obtain size of mRNA feature matrix
    n_feature_mRNA=X_mRNA_train[0].shape[0]
    
    # Obtain size of promoter feature matrix
    n_feature_promoter=X_promoter_train[0].shape[0]
    
    def data_generator():
        while True:
            
            # Shuffle the data at each epoch
            if shuffle:
                shuffle_indices = np.random.permutation(np.arange(data_size))
                shuffled_data_mRNA = X_mRNA_train[shuffle_indices]
                shuffled_data_promoter = X_promoter_train[shuffle_indices]
                shuffled_labels = Y_train[shuffle_indices]
            else:
                shuffled_data_mRNA = X_mRNA_train
                shuffled_data_promoter = X_promoter_train
                shuffled_labels = Y_train
                
            # Generate data for each batch
            for batch_num in range(num_batches_per_epoch):
                start_index = batch_num * batch_size
                end_index = min((batch_num + 1) * batch_size, data_size)
                
                # Prepare mRNA feature matrix
                # Obtain max length of mRNA sequence in this batch
                seq_len_mRNA=map(lambda x: np.shape(x)[1], shuffled_data_mRNA[start_index: end_index])
                max_seq_mRNA=np.max(seq_len_mRNA) 
                
                # Initialize mRNA feature matrix
                X_mRNA=np.zeros((end_index-start_index, max_seq_mRNA, n_feature_mRNA))
                k=0
                for i in range(start_index,end_index):
                    X_mRNA[k,0:seq_len_mRNA[k],:]=shuffled_data_mRNA[i].transpose()
                    k=k+1
                    
                # Prepare promoter feature matrix
                # Obtain max length of promoter sequence in this batch
                seq_len_promoter=map(lambda x: np.shape(x)[1], shuffled_data_promoter[start_index: end_index])
                max_seq_promoter=np.max(seq_len_promoter)
                
                # Initialize promoter feature matrix
                X_promoter=np.zeros((end_index-start_index, max_seq_promoter, n_feature_promoter))
                k=0
                for i in range(start_index,end_index):
                    X_promoter[k,0:seq_len_promoter[k],:]=shuffled_data_promoter[i].transpose()
                    k=k+1                   
                
                # Prepare target data
                Y = shuffled_labels[start_index: end_index]
                
                yield [X_mRNA,X_promoter], Y

    return num_batches_per_epoch, data_generator()

def batch_iter_DeepLIFT(X_mRNA_train,
                        X_promoter_train,
                        Y_train,
                        batch_size, 
                        med_mRNA_len,
                        med_promoter_len,
                        shuffle=True):
    import numpy as np
    
    # X_mRNA_train: mRNA features
    # X_promoter_train: promoter features
    # Y_train: target data
    # batch_size 
    # shuffle=True
    
    # Number of data in each batch
    data_size = len(Y_train)
    num_batches_per_epoch = int((data_size - 1) / batch_size) + 1
    
    # Obtain the size of mRNA feature matrix
    n_feature_mRNA=X_mRNA_train[0].shape[0]
    
    # Obtain the size of promotor feature matrix
    n_feature_promoter=X_promoter_train[0].shape[0]
    
    def data_generator():
        while True:
            # Shuffle the data at each epoch
            if shuffle:
                shuffle_indices = np.random.permutation(np.arange(data_size))
                shuffled_data_mRNA = X_mRNA_train[shuffle_indices]
                shuffled_data_promoter = X_promoter_train[shuffle_indices]
                shuffled_labels = Y_train[shuffle_indices]
            else:
                shuffled_data_mRNA = X_mRNA_train
                shuffled_data_promoter = X_promoter_train
                shuffled_labels = Y_train

            # Generate data for each batch
            for batch_num in range(num_batches_per_epoch):
                start_index = batch_num * batch_size
                end_index = min((batch_num + 1) * batch_size, data_size)
                
                # Prepare mRNA feature matrix
                # Bbtain max length of mRNA sequence in this batch
                seq_len_mRNA=map(lambda x: np.shape(x)[1], shuffled_data_mRNA[start_index: end_index])
                max_seq_mRNA=np.max(seq_len_mRNA) 
                max_seq_mRNA=np.max([max_seq_mRNA,med_mRNA_len])
                
                # Initializing mRNA feature matrix
                X_mRNA=np.zeros((end_index-start_index, max_seq_mRNA, n_feature_mRNA))
                k=0
                for i in range(start_index,end_index):
                    X_mRNA[k,0:seq_len_mRNA[k],:]=shuffled_data_mRNA[i].transpose()
                    k=k+1
                    
                # Prepare promoter feature matrix
                # Obtain max length of promoter sequence in this batch
                seq_len_promoter=map(lambda x: np.shape(x)[1], shuffled_data_promoter[start_index: end_index])
                max_seq_promoter=np.max(seq_len_promoter)
                max_seq_promoter=np.max([max_seq_promoter,med_promoter_len])
                
                # Initializing promoter feature matrix
                X_promoter=np.zeros((end_index-start_index, max_seq_promoter, n_feature_promoter))
                k=0
                for i in range(start_index,end_index):
                    X_promoter[k,0:seq_len_promoter[k],:]=shuffled_data_promoter[i].transpose()
                    k=k+1                   
                
                
                # Prepare target data
                Y = shuffled_labels[start_index: end_index]
                
                yield [X_mRNA,X_promoter], Y

    return num_batches_per_epoch, data_generator()