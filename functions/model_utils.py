def test_prediction(outloc,
                    best_model,
                    X_mRNA_test,
                    X_promoter_test,
                    Y_test):
    
    # outloc: a full path to a directory of the best model
    # best_model: name of the best model
    # X_mRNA_test: mRNA annoation data
    # X_promoter_test: Promoter annoatation data
    # Y_test: Transcriptome data     
    
    import layer_utils
    import metrics
    from keras.models import load_model
    import data
    import numpy as np
    import os 
    
    if not os.path.exists(outloc+best_model+'/test_data/'):
        os.makedirs(outloc+"/"+best_model+'/test_data/')
    
    # Load best model
    model = load_model(outloc+best_model+'_model.h5',
                       custom_objects={'pcor': metrics.pcor,
                                      'GlobalSumPooling1D': layer_utils.GlobalSumPooling1D})

    # Batching testing data
    batch_size=128
    test_steps, test_batches = data.batch_iter(X_mRNA_test.values[:,1],
                                               X_promoter_test.values[:,1],
                                               Y_test.values[:,1:],
                                               batch_size,
                                               shuffle=False)

    # Making prediction
    pred=[]
    actu=[]
    for i in range(test_steps):
        a=test_batches.next()
        b=model.predict(a[0])
        pred.append(b)
        actu.append(np.vstack(a[1]))

    pred=np.vstack(pred)
    actu=np.vstack(actu)
    
    # Save actual and predicted gene expression
    np.savetxt(outloc+best_model+'/test_data/actual.txt',
               actu, delimiter='\t')
    np.savetxt(outloc+best_model+'/test_data/prediction.txt',
               pred, delimiter='\t')
    X_mRNA_test['Name'].to_csv(outloc+best_model+'/test_data/geneid.txt', header=False, index=False, sep='\t')
    
    # gzip text files
    os.system("gzip "+outloc+best_model+'/test_data/actual.txt')
    os.system("gzip "+outloc+best_model+'/test_data/prediction.txt')    
    os.system("gzip "+outloc+best_model+'/test_data/geneid.txt')    
    
def compute_DeepLIFT(outloc,
             best_model,
             X_mRNA_test,
             X_promoter_test,
             Y_test):
    
    # outloc: a full path to a directory of the best model
    # best_model: name of the best model
    # X_mRNA_test: mRNA annoation data
    # X_promoter_test: Promoter annoatation data
    # Y_test: Transcriptome data    
    
    from keras.layers import Dense
    from keras.models import Model
    import layer_utils
    import metrics
    from keras.models import load_model
    import data
    import numpy as np
    import os
    
    # Load model
    model = load_model(outloc+best_model+'_model.h5',
                       custom_objects={'pcor': metrics.pcor,
                                       'GlobalSumPooling1D': layer_utils.GlobalSumPooling1D})

    # Get parameters for the last dens layer
    dens_parameter = model.layers[-1].get_weights()
    
    # Construct a single output model 
    # Adding new layers
    fc = Dense(1,activation='linear')(model.layers[-2].output)
    new_model = Model(inputs=model.input, outputs=fc)
    
    # Paramter for background distribution
    med_mRNA_len=int(np.median(map(lambda x:x.shape[1],X_mRNA_test.values[:,1])))
    med_promoter_len=int(np.median(map(lambda x:x.shape[1],X_promoter_test.values[:,1])))
    gene_names_test = X_mRNA_test.values[:,0]
    
    # DeepLIFT score for each sample
    for out_indx in range(dens_parameter[1].shape[0]):
        print(out_indx)

        # Set dens parameters
        new_model.layers[-1].set_weights([dens_parameter[0][:,out_indx:(out_indx+1)],dens_parameter[1][out_indx:(out_indx+1)]])

        import shap
        method='deepexplainer' 
        if not os.path.exists(outloc+best_model+'/DeepLIFT/'):
            os.makedirs(outloc+"/"+best_model+'/DeepLIFT/')
        
        # Speficy output file names
        outfile_name_at1=outloc+"/"+best_model+'/DeepLIFT/RNA_'+str(out_indx)+'.txt'
        outfile_name_at2=outloc+"/"+best_model+'/DeepLIFT/DNA_'+str(out_indx)+'.txt'

        # Batching testing data
        batch_size=256*4
        test_steps, test_batches = data.batch_iter_DeepLIFT(X_mRNA_test.values[:,1],
                                                            X_promoter_test.values[:,1],
                                                            Y_test.values[:,1:],
                                                            batch_size,
                                                            med_mRNA_len,
                                                            med_promoter_len,
                                                            shuffle=False)

        for i in range(test_steps):
            
            xs_test,ys_test=test_batches.next()

            # Reshape background
            xs_background=[]
            xs_background.append(np.zeros((1,xs_test[0].shape[1],xs_test[0].shape[2])))
            xs_background.append(np.zeros((1,xs_test[1].shape[1],xs_test[1].shape[2])))
            xs_background[0][0,0:med_mRNA_len,0]=1
            xs_background[1][0,0:med_promoter_len,0]=1

            # Compute DeepLIFT scores
            e = shap.DeepExplainer(new_model, xs_background)
            shap_values = e.shap_values(xs_test)
            shap_values=shap_values[0]

            # Export DeepLIFT scores to text
            with open(outfile_name_at1, 'a') as f_handle:
                input_k=0
                for j in range(shap_values[input_k].shape[0]):
                    seq_indx=xs_test[input_k][j,:,0]>0
                    feature_vector=map(str, np.sum(shap_values[input_k][j,seq_indx,:],axis=0))
                    out_txt=str(out_indx)+'\t'+gene_names_test[j+i*batch_size]+'\t'+','.join(feature_vector)+'\n'
                    f_handle.write(out_txt)

            with open(outfile_name_at2, 'a') as f_handle:
                input_k=1
                for j in range(shap_values[input_k].shape[0]):
                    seq_indx=xs_test[input_k][j,:,0]>0
                    feature_vector=map(str, np.sum(shap_values[input_k][j,seq_indx,:],axis=0))
                    out_txt=str(out_indx)+'\t'+gene_names_test[j+i*batch_size]+'\t'+','.join(feature_vector)+'\n'
                    f_handle.write(out_txt)

        # gzip text files
        os.system("gzip "+outfile_name_at1)
        os.system("gzip "+outfile_name_at2)

def coexpression_with_binding_site_removal(genes,
                                           outloc,
                                           best_model,
                                           X_mRNA_test,
                                           X_promoter_test,
                                           Y_test):
    # genes: genes for testing coexpression
    # outloc: a full path to a directory of the best model
    # best_model: name of the best model
    # X_mRNA_test: mRNA annoation data
    # X_promoter_test: Promoter annoatation data
    # Y_test: Transcriptome data     
    
    import layer_utils
    import metrics
    from keras.models import load_model
    import data
    import copy
    from scipy.stats import spearmanr
    import numpy as np
    import os
    
    if not os.path.exists(outloc+best_model+'/binding_site_removal/'):
        os.makedirs(outloc+best_model+'/binding_site_removal/')

    outfile_name=outloc+best_model+'/binding_site_removal/coexpression.txt'
    
    # Load best model
    model = load_model(outloc+best_model+'_model.h5',
                       custom_objects={'pcor': metrics.pcor,
                                      'GlobalSumPooling1D': layer_utils.GlobalSumPooling1D})
    
    # Batch data
    batch_size=128
    test_steps, test_batches = data.batch_iter(X_mRNA_test.query("Name in @genes").values[:,1],
                                               X_promoter_test.query("Name in @genes").values[:,1],
                                               Y_test.query("Name in @genes").values[:,1:],
                                               batch_size,
                                               shuffle=False)
    X=test_batches.next()    
    X_copy = copy.deepcopy(X) 
    
    # None mutated result
    np.random.seed(seed=1234)
    indx0=np.random.uniform(0,1,X[0][0].shape[1]) < -0.5
    indx1=np.random.uniform(0,1,X[0][1].shape[1]) < -0.5
    
    # Predict expression
    Y_predicted=model.predict(X_copy[0])
    
    # Compute correlation
    res=spearmanr(Y_predicted[0,0:-1],Y_predicted[1,0:-1])
    
    #  Save correlation
    with open(outfile_name, 'a') as f_handle:
        out_txt="no_mutation"+'\t'+str(res.correlation)+'\t'+str(res.pvalue)+'\t'+','.join(map(str,np.where(indx0)[0]))+'\t'+','.join(map(str,np.where(indx1)[0]))+'\n'
        f_handle.write(out_txt)
    
    # Simulate random removal of binding sites
    for i in range(10000):
        for tr_id in range(2):
            gene_name = X_mRNA_test.query("Name in @genes").values[tr_id,0]
            X_copy = copy.deepcopy(X) 
            
            indx0=np.random.uniform(0,1,X[0][0].shape[1]) > 0.5
            indx1=np.random.uniform(0,1,X[0][1].shape[1]) > 0.5
            
            # Remove binding sites
            X_copy[0][0][tr_id,indx0,1:]=0
            X_copy[0][1][tr_id,indx1,1:]=0
            
            # Predict expression
            Y_predicted=model.predict(X_copy[0])
            
            # Compute correlation
            res=spearmanr(Y_predicted[0,0:-1],Y_predicted[1,0:-1])
            
            #  Save correlation
            with open(outfile_name, 'a') as f_handle:
                out_txt=gene_name+'\t'+str(res.correlation)+'\t'+str(res.pvalue)+'\t'+','.join(map(str,np.where(indx0)[0]))+'\t'+','.join(map(str,np.where(indx1)[0]))+'\n'
                f_handle.write(out_txt)
    
    os.system("gzip "+outfile_name)

def coexpression_with_KO(genes,
                         outloc,
                         best_model,
                         X_mRNA_test,
                         X_promoter_test,
                         Y_test):
    # genes: genes for testing coexpression
    # outloc: a full path to a directory of the best model
    # best_model: name of the best model
    # X_mRNA_test: mRNA annoation data
    # X_promoter_test: Promoter annoatation data
    # Y_test: Transcriptome data     
    
    import layer_utils
    import metrics
    from keras.models import load_model
    import data
    import copy
    from scipy.stats import spearmanr
    import numpy as np
    import os
    
    if not os.path.exists(outloc+best_model+'/regulator_KO/'):
        os.makedirs(outloc+best_model+'/regulator_KO/')

    outfile_name=outloc+best_model+'/regulator_KO/coexpression.txt'

    # Load best model
    model = load_model(outloc+best_model+'_model.h5',
                       custom_objects={'pcor': metrics.pcor,
                                      'GlobalSumPooling1D': layer_utils.GlobalSumPooling1D})

    # Batch data
    batch_size=128
    test_steps, test_batches = data.batch_iter(X_mRNA_test.query("Name in @genes").values[:,1],
                                               X_promoter_test.query("Name in @genes").values[:,1],
                                               Y_test.query("Name in @genes").values[:,1:],
                                               batch_size,
                                               shuffle=False)
    X=test_batches.next()
    X_copy = copy.deepcopy(X)     
    
    none_zero_indx0=np.where(np.sum(np.sum(X[0][0],axis=1),axis=0)>0)[0][1:]
    none_zero_indx1=np.where(np.sum(np.sum(X[0][1],axis=1),axis=0)>0)[0][1:]
    
    # None mutated result
    np.random.seed(seed=1234)
    indx0=none_zero_indx0[np.random.uniform(0,1,len(none_zero_indx0)) < -0.5]
    indx1=none_zero_indx1[np.random.uniform(0,1,len(none_zero_indx1)) < -0.5]

    # Predict expression
    Y_predicted=model.predict(X_copy[0])
    
    # Compute correlation
    res=spearmanr(Y_predicted[0,0:-1],Y_predicted[1,0:-1])
    
    #  Save correlation
    with open(outfile_name, 'a') as f_handle:
        out_txt=str(res.correlation)+'\t'+str(res.pvalue)+'\t'+','.join(map(str,indx0))+'\t'+','.join(map(str,indx1))+'\n'
        f_handle.write(out_txt)

    # Simulate random removal of binding sites
    for i in range(10000):
        X_copy = copy.deepcopy(X) 

        indx0=none_zero_indx0[np.random.uniform(0,1,len(none_zero_indx0)) > 0.5]
        indx1=none_zero_indx1[np.random.uniform(0,1,len(none_zero_indx1)) > 0.5]

        # Remove binding of regulators
        X_copy[0][0][0,:,indx0]=0
        X_copy[0][1][0,:,indx1]=0

        # Predict expression
        Y_predicted=model.predict(X_copy[0])

        # Compute correlation
        res=spearmanr(Y_predicted[0,0:-1],Y_predicted[1,0:-1])
        
        #  Save correlation
        with open(outfile_name, 'a') as f_handle:
            out_txt=str(res.correlation)+'\t'+str(res.pvalue)+'\t'+','.join(map(str,indx0))+'\t'+','.join(map(str,indx1))+'\n'
            f_handle.write(out_txt)
    
    os.system("gzip "+outfile_name)