def define_network(params):
    import layer_utils
    import metrics
    
    from keras.models import Sequential, Model
    from keras.layers import Dense, Conv1D, Input, BatchNormalization, Activation, concatenate
    from keras.utils import to_categorical
    from keras import optimizers
    
    # Defining mRNA model
    RNA_model = Sequential()
    RNA_model.add(Conv1D(int(params['RNA_n_channel_1st']),
                                int(params['RNA_conv_strides']),
                                input_shape=(None, int(params['n_feature_mRNA'])),
                        use_bias=False))
    if params['ConvRelu']=='Yes':
        RNA_model.add(Activation('relu'))
    for k in range(int(params['RNA_n_ConvLayer'])):
        RNA_model.add(Conv1D(int(params['RNA_n_channel_1st']),
                                    int(params['RNA_conv_strides']),
                             use_bias=False))
        if params['ConvRelu']=='Yes':
            RNA_model.add(Activation('relu'))
    RNA_model.add(layer_utils.GlobalSumPooling1D())
    # Getting a tensor with the output of our RNA model:
    RNA_feature_input = Input(shape=(None, int(params['n_feature_mRNA'])))
    encoded_RNA = RNA_model(RNA_feature_input)
    
    
    # Defining promoter model
    DNA_model = Sequential()
    DNA_model.add(Conv1D(int(params['DNA_n_channel_1st']),
                                int(params['DNA_conv_strides']),
                                input_shape=(None, int(params['n_feature_promoter'])),
                        use_bias=False))
    if params['ConvRelu']=='Yes':
        DNA_model.add(Activation('relu'))
    for k in range(int(params['DNA_n_ConvLayer'])):
        DNA_model.add(Conv1D(int(params['DNA_n_channel_1st']),
                                    int(params['DNA_conv_strides']),
                        use_bias=False))
        if params['ConvRelu']=='Yes':
            DNA_model.add(Activation('relu'))
    DNA_model.add(layer_utils.GlobalSumPooling1D())
    # Getting a tensor with the output of our promoter model:
    DNA_feature_input = Input(shape=(None, int(params['n_feature_promoter'])))
    encoded_DNA = DNA_model(DNA_feature_input)
    
    # Concatenating the mRNA vector and the promoter vector:
    merged = concatenate([encoded_RNA, encoded_DNA])        
    
    # Stacking a deep densely-connected network on top
    for k in range(int(params['Last_fullConLayer'])):
        merged = Dense(int(params['Last_n_channel']))(merged)
        if params['FullRelu']=='Yes':
            merged = Activation('relu')(merged)
        merged = BatchNormalization()(merged)
    
    # Making preditions
    output=Dense(params['n_out'], activation='linear')(merged)
    
    # Putting all inputs and outputs together 
    model = Model(inputs=[RNA_feature_input, DNA_feature_input],
                          outputs=output)        
    
    # Compiling model
    model.compile(loss='mean_squared_error',
                  optimizer=optimizers.Adam(lr=params['lr']),
                  metrics=[metrics.pcor])
    
    return model

