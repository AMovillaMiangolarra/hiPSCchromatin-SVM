import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

# Function to binarise transcriptional output into 0 and 1. 
#Uses the best performing between Kmeans in log space and in linear space (based on inertia changes).
def binarize(df, lines):
    col_id=df.columns[0]

    for line in lines:
        df[f'{line}_binary'] = 0
        
    # Calculate Kmeans binarization in linear and logscale.
    for i in range(df.shape[0]):
        RNA_across_lines = np.array(df[lines].iloc[i], dtype=float)
        kmeans_linear = KMeans(n_clusters=2, init="random").fit(RNA_across_lines.reshape(-1,1))
        kmeans_linear1 = KMeans(n_clusters=1, init="random").fit(RNA_across_lines.reshape(-1,1))
        if sum(RNA_across_lines==0)==0:
            kmeans_log = KMeans(n_clusters=2, init="random").fit(np.log(RNA_across_lines.reshape(-1,1)))
            kmeans_log1 = KMeans(n_clusters=1, init="random").fit(np.log(RNA_across_lines.reshape(-1,1)))
            loginertia=kmeans_log.inertia_/kmeans_log1.inertia_
        else:
            loginertia=1
        #Choose linear or logarithmic. 
        if kmeans_linear.inertia_/kmeans_linear1.inertia_<loginertia:
            binar_gene=kmeans_linear.labels_
            if kmeans_linear.cluster_centers_[0]>kmeans_linear.cluster_centers_[1]:
                binar_gene=(binar_gene-1)*(-1)
        else:
            binar_gene=kmeans_log.labels_
            if kmeans_log.cluster_centers_[0]>kmeans_log.cluster_centers_[1]:
                binar_gene=(binar_gene-1)*(-1)
        # set binary value for each line
        linesp=0
        for line in lines:
            df.loc[i, f'{line}_binary'] = binar_gene[linesp]
            linesp+=1
    return(df)

from sklearn import pipeline
from sklearn import svm
from sklearn import preprocessing

# Function to train an SVM for each gene and process the outputs to find the underlying regulatory mechanisms.  
def SVM_analysis(DEGs, vars, epig_vars):
    #Pipeline for SVM training
    _pipeline = pipeline.Pipeline([
        ('normalize', preprocessing.MinMaxScaler()), 
        ('svc', svm.SVC(kernel='linear', C = 100)) #C=100 to highly penalise errors.
    ])
    
    cols=[['geneid'], epig_vars,['intercept'], [col+'_scale' for col in epig_vars]]
    flat_cols = [x
        for xs in cols
        for x in xs]
    coefficients = pd.DataFrame(columns=flat_cols)

    #Performing SVM training on a per DEG basis
    for gene in DEGs:
        gene_instance = vars[vars['geneid'] == gene]
        result = _pipeline.fit(gene_instance[epig_vars], gene_instance['expression'].values)
        coefficients.loc[len(coefficients)] = [gene] + list(_pipeline.named_steps['svc'].coef_[0]) + [_pipeline.named_steps['svc'].intercept_[0]] + list(_pipeline.named_steps['normalize'].scale_)

    #Post-processing of the SVM coefficients
    norm=np.zeros(len(DEGs))
    for epivar in epig_vars:
        norm+=coefficients[epivar]**2
    coefficients["Normaliz"]=norm
    for index in range(coefficients.shape[0]):
        if coefficients.loc[index,"Normaliz"]>0:
            for epivar in epig_vars:
                coefficients.loc[index,"n"+epivar]=coefficients.loc[index,epivar]/np.sqrt(coefficients.loc[index,'Normaliz'])
        else:
            for epivar in epig_vars:
                coefficients.loc[index,"n"+epivar]==0
                
    for epivar in epig_vars:
        coefficients["n"+epivar+"_scale"]=[-np.log(x) if x<1 else 0 for x in coefficients[epivar+'_scale']]

    for epivar in epig_vars:
        maxnorm=max(abs(coefficients['n'+epivar]*coefficients["n"+epivar+"_scale"]))
        coefficients["n2"+epivar+"_scale"]=coefficients['n'+epivar]*coefficients["n"+epivar+"_scale"]/maxnorm
    
    return(coefficients)

# Function to binarise transcriptional outputs (as above) but for training lines separate from a test line.
# Produces also the binarisation of the test line which can be used as a ground truth of the test.
def binarize_test(df, training_lines, test_line):
    col_id=df.columns[0]

    for line in training_lines:
        df[f'{line}_binary'] = 0
    
    test_line_transcr=[]
    test_line_df=pd.DataFrame(columns=['geneid', 'transcription'])
    # Calculate Kmeans binarization in linear and logscale.
    for i in range(df.shape[0]):
        RNA_across_lines = np.array(df[training_lines].iloc[i], dtype=float)
        kmeans_linear = KMeans(n_clusters=2, init="random").fit(RNA_across_lines.reshape(-1,1))
        kmeans_linear1 = KMeans(n_clusters=1, init="random").fit(RNA_across_lines.reshape(-1,1))
        if sum(RNA_across_lines==0)==0:
            kmeans_log = KMeans(n_clusters=2, init="random").fit(np.log(RNA_across_lines.reshape(-1,1)))
            kmeans_log1 = KMeans(n_clusters=1, init="random").fit(np.log(RNA_across_lines.reshape(-1,1)))
            loginertia=kmeans_log.inertia_/kmeans_log1.inertia_
        else:
            loginertia=1
        #Choose linear or logarithmic. 
        if kmeans_linear.inertia_/kmeans_linear1.inertia_<loginertia:
            binar_gene=kmeans_linear.labels_
            if abs(df[test_line].iloc[i]-kmeans_linear.cluster_centers_[0])>abs(df[test_line].iloc[i]-kmeans_linear.cluster_centers_[1]):
                test_line_transcr.append(1)
            else:
                test_line_transcr.append(0)
            if kmeans_linear.cluster_centers_[0]>kmeans_linear.cluster_centers_[1]:
                binar_gene=(binar_gene-1)*(-1)
                test_line_transcr[i]=-(test_line_transcr[i]-1)
        else:
            binar_gene=kmeans_log.labels_
            if (df[test_line].iloc[i]!=0):
                logtest=np.log(df[test_line].iloc[i])
            else:
                logtest=-10
            if abs(logtest-kmeans_log.cluster_centers_[0])>abs(logtest-kmeans_log.cluster_centers_[1]):
                test_line_transcr.append(1)
            else:
                test_line_transcr.append(0)
            if kmeans_log.cluster_centers_[0]>kmeans_log.cluster_centers_[1]:
                binar_gene=(binar_gene-1)*(-1)
                test_line_transcr[i]=-(test_line_transcr[i]-1)

                
        # set binary value for each line
        linesp=0
        for line in training_lines:
            df.loc[i, f'{line}_binary'] = binar_gene[linesp]
            linesp+=1

        test_line_df.loc[i, 'geneid'] = df['geneid'].iloc[i]
        test_line_df.loc[i, 'transcription'] = test_line_transcr[i]
            
            
    return(df, test_line_df)


# Function to train the SVM with the aim of predicting additional datapoints.
# Provides scores for the prediction via Platt's scaling.
def SVM_prediction(DEGs, vars, genes_test, epig_vars):
    _pipeline = pipeline.Pipeline([
        ('normalize', preprocessing.MinMaxScaler()), 
        ('svc', svm.SVC(kernel='linear', C = 100, probability=True)) #C=100 so that it penalises highly errors
    ])

    probs_prediction=[]
    prediction=[]
    i=0
    for gene in DEGs:
        gene_instance = vars[vars['geneid'] == gene]
        result = _pipeline.fit(gene_instance[epig_vars], gene_instance['expression'].values)
        test_vect=genes_test.loc[genes_test['geneid']==gene,epig_vars]
        probprec =_pipeline.predict_proba(test_vect)
        probs_prediction.append(probprec[0][1])
        prediction.append(_pipeline.predict(test_vect)[0])
    
    return(probs_prediction, prediction)