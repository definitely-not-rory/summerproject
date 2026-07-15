from imports import *

def get_iom_scales(sample,path='/cosma8/data/dp262/dc-dodd1/',save_dir='/cosma/apps/durham/dc-coll7/auriga/iom_scales'):

    os.makedirs(save_dir,exist_ok=True)

    if os.path.exists(f'{path}{sample}.hdf5')!=True:
        print(f'\nNo data for {sample} in {path} located.')
        return
    else:
        print(f'\nImporting data from {path}{sample}.')
        full_path=f'{path}{sample}.hdf5'

    if os.path.exists(f'{save_dir}/{sample}.json')==True:
        valid_input=False 
        valid_inputs=['y','Y','n','N']
        
        while valid_input==False: 
            overwrite=input('\n.json IOM file for %s already exists in %s. Do you want to overwrite it (y/n)?:' % (sample, save_dir)) #Request permission to overwrite pre-existing file.
            if overwrite in valid_inputs: 
                valid_input=True
                if overwrite in valid_inputs[2:]: 
                    print(f'\nIOM scales for {sample} located in {save_dir}, loading from .json file.')

                    with open (f'{save_dir}/{sample}.json', 'r') as f:
                        scales = json.load(f)
                        f.close()

                    return scales
                else:
                    print('\nOverwriting %s.json.' % sample)
    else:
        print(f'\nNo existing IOM scale file for {sample} detected in {save_dir}, generating.')
    
    df = vaex.open(full_path)
    
    iom_axes=['En','Lperp','Lz']

    scales={iom: [np.floor(np.min(df.min(iom))),np.ceil(np.max(df.max(iom)))] for iom in iom_axes}
    
    scales['En'][1]=0
    scales['Lperp'][0]=0
    
    Lz_extent=max(np.abs(scales['Lz']))
    scales['Lz'][0]=-1*Lz_extent
    scales['Lz'][1]=Lz_extent

    with open (f'{save_dir}/{sample}.json', 'w') as f:
        json.dump(scales,f)
        f.close()
    
    with open (f'{save_dir}/{sample}.json', 'r') as f:
        scales = json.load(f)
        f.close()
    
    print(f'\nSaved IOM scales at {save_dir}{sample}.json')

    return scales

def mahalanobis_distance(m1,m2,cov1,cov2):
    a=m2-m1
    b=np.linalg.inv((cov1+cov2))
    return a.dot(b).dot(a.T)

def cluster_distance_matrix(df,halo,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',distance_metric='euclidean'):
    
    if vtoomre==False:
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    features=['scaled_En', 'scaled_Lperp', 'scaled_Lz']

    run_name=f'{halo}_{lsr_def}_{selection_type}'
    save_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/results'

    if os.path.exists(f'{save_dir}/{run_name}_DistanceMatrix_{distance_metric}.npy')==True:
        print(f'Distance matrix for {run_name} located, loading.')
        dist_matrix=np.load(f'{save_dir}/{run_name}_DistanceMatrix_{distance_metric}.npy')
        return dist_matrix
    else:
        print(f'No distance matrix for {run_name} found in {save_dir}, generating.')

    jit_md=jit(mahalanobis_distance,nopython=True,fastmath=True,cache=True)

    unique_labels= np.unique(df.evaluate('label'))
    N_unique=len(unique_labels)

    dist_matrix=np.empty((N_unique*(N_unique-1))//2,dtype=np.double)

    counter=0

    for i in range(0,N_unique-1):
        for j in range(i+1,N_unique):
            cluster1 = df.filter('label==%s'%unique_labels[i]).extract()
            cluster2 = df.filter('label==%s'%unique_labels[j]).extract()

            U = cluster1[features].values
            V = cluster2[features].values

            if distance_metric == 'euclidean':
                min_dist = 1000000
                for u in U:
                    for v in V:
                        dist = np.linalg.norm(u-v)
                        if(dist<min_dist):
                            min_dist = dist
                
            elif distance_metric == 'mahalanobis':
                min_dist = np.sqrt(jit_md(np.mean(U, axis=0), np.mean(V, axis=0),
                                         np.cov(U.T), np.cov(V.T)))

            dist_matrix[counter]=min_dist 
            counter+=1
    
    np.save(f'{save_dir}/{run_name}_DistanceMatrix_{distance_metric}.npy',dist_matrix)

    return dist_matrix
        
def small_num_KS(feh1,feh2,Nstars_lim=20,NMC=100):
    
    main_feh,test_feh=feh1,feh2
    if len(main_feh)<len(test_feh): 
        main_feh, test_feh = feh2,feh1

    KS, pval = stats.ks_2samp(main_feh,test_feh,mode='exact')

    if len(test_feh)<Nstars_lim: 
        orig_pval = pval
    else:
        orig_pval = pval
        comp_percentage = 999.9
        less_percentage = 999.9

    if len(test_feh)<Nstars_lim:            
        NNs,pvals=[],[]
        count=0
        while count<NMC:
            random_sample_feh=np.random.choice(main_feh,len(test_feh))
            KS,pval=stats.ks_2samp(main_feh, random_sample_feh, mode='exact')
            NNs.append(count)
            pvals.append(pval)
            count+=1
        
        less_coords=np.where((np.array(pvals)<orig_pval))
        less_percentage=float(len(np.array(pvals)[less_coords]))/float(len(np.array(pvals)))
        
        same_coords=np.where((np.array(pvals)>0.05))
        same_distr=100.0*float(len(np.array(pvals)[same_coords]))/float(len(np.array(pvals)))        
        
        diff_coords=np.where((np.array(pvals)<0.05))
        diff_distr=100.0*float(len(np.array(pvals)[diff_coords]))/float(len(np.array(pvals)))        
        
        if orig_pval<0.05:
            comp_percentage=diff_distr
        else:
            comp_percentage=same_distr
    
        return orig_pval,comp_percentage,less_percentage
    
    