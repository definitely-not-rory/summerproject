from imports import *

#Function dynamically that generates .json file containing scales for IOM features.
def get_iom_scales(sample,path='/cosma8/data/dp262/dc-dodd1/',save_dir='/cosma/apps/durham/dc-coll7/auriga/iom_scales'):

    os.makedirs(save_dir,exist_ok=True) #Creates directory to save external files to.

    if os.path.exists(f'{path}{sample}.hdf5')!=True: #Verifies that raw particle data exists.
        print(f'\nNo data for {sample} in {path} located.')
        return
    else:
        print(f'\nImporting data from {path}{sample}.')
        full_path=f'{path}{sample}.hdf5' #Initialises quick access variable for loading data.

    if os.path.exists(f'{save_dir}/{sample}.json')==True: #Detects if IOM scaling file has been created before.
        valid_input=False #Initialises valid input variables.
        valid_inputs=['y','Y','n','N']
        
        while valid_input==False: #Continue iterating until valid input is entered.
            overwrite=input('\n.json IOM file for %s already exists in %s. Do you want to overwrite it (y/n)?:' % (sample, save_dir)) #Request permission to overwrite pre-existing file.
            if overwrite in valid_inputs: #Detect if overwrite request input is a valid input.
                valid_input=True #Breaks iteration if input is valid.
                if overwrite in valid_inputs[2:]: #Halts function and loads pre-existing IOM scales if overwrite is denied.
                    print(f'\nIOM scales for {sample} located in {save_dir}, loading from .json file.')

                    with open (f'{save_dir}/{sample}.json', 'r') as f:
                        scales = json.load(f)
                        f.close()

                    return scales
                else:
                    print('\nOverwriting %s.json.' % sample)
    else:
        print(f'\nNo existing IOM scale file for {sample} detected in {save_dir}, generating.')
    
    df = vaex.open(full_path) #Loads sample data to calculate scaling.
    
    iom_axes=['En','Lperp','Lz'] #List of IOM features to calculate scales for.

    scales={iom: [int(np.floor(np.min(df.min(iom)))),int(np.ceil(np.max(df.max(iom))))] for iom in iom_axes} #Generates scale for each IOM as most extreme integer of each feature.
    
    scales['En'][1]=0 #Assigns relevant bounds to 0.
    scales['Lperp'][0]=0
    
    Lz_extent=max(np.abs(scales['Lz'])) #Makes Lz scale symmetric about 0 out to more extreme upper or lower bound.
    scales['Lz'][0]=int(-1*Lz_extent)
    scales['Lz'][1]=int(Lz_extent)

    with open (f'{save_dir}/{sample}.json', 'w') as f: #Saves generated IOM scales externally.
        json.dump(scales,f)
        f.close()
    
    with open (f'{save_dir}/{sample}.json', 'r') as f: #Loads IOM scales from file for exporting.
        scales = json.load(f)
        f.close()
    
    print(f'\nSaved IOM scales at {save_dir}{sample}.json')

    return scales

#Function that calculates the Mahalanobis distance between 2 clusters.
def mahalanobis_distance(m1,m2,cov1,cov2):
    a=m2-m1
    b=np.linalg.inv((cov1+cov2))
    return a.dot(b).dot(a.T)

#Function that generates the IOM distance metric of a dataset for a given distance metric (Euclidean or Mahalanobis).
def cluster_distance_matrix(df,halo,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',distance_metric='euclidean',regen=False):
    
    if vtoomre==False: #Allocates relevant label to particle subselection method ('accreted' -> particles tagged as accreted by AURIGA, 'vtoomre' -> selection over all particles using Toomre velocity threshold).
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    if regen==False:
        regen_select=' not '
    else:
        regen_select=' '

    features=['scaled_En', 'scaled_Lperp', 'scaled_Lz'] #List of features to use when determining distances.

    run_name=f'{halo}_{lsr_def}_{selection_type}' #Initialises quick access strings for clustering run label and file paths.
    save_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/results'

    if os.path.exists(f'{save_dir}/{run_name}_DistanceMatrix_{distance_metric}.npy')==True: #Detects if distance matrix has been generated previously.
        print(f'Distance matrix for {run_name} located. Regeneration{regen_select}selected.')
        if regen==False: #Checks if regeneration of matrix has been requested.
            dist_matrix=np.load(f'{save_dir}/{run_name}_DistanceMatrix_{distance_metric}.npy') #Loads pre-existing distance matrix if it already exists.
            return dist_matrix
    else:
        print(f'No distance matrix for {run_name} found in {save_dir}, generating.')

    jit_md=jit(mahalanobis_distance,nopython=True,fastmath=True,cache=True) #Assigns mahalanobis calculation 'jit' object for faster computation.

    unique_labels= np.unique(df.evaluate('label')) #Retrieves list and number of unique clusters.
    N_unique=len(unique_labels)

    dist_matrix=np.empty((N_unique*(N_unique-1))//2,dtype=np.double) #Initialises empty distance matrix for storage.

    counter=0 #Creates variable for iteration step tracking.

    for i in tqdm(range(0,N_unique-1)):#Iterates across one dimension of distance matrix.
        for j in range(i+1,N_unique): #Iterates over other dimension of distance matrix.
            cluster1 = df.filter('label==%s'%unique_labels[i]).extract() #Selects all stars in clusters associated with current selected distance matrix element. 
            cluster2 = df.filter('label==%s'%unique_labels[j]).extract()

            U = cluster1[features].values #Extracts values of IOM features of each cluster's stars.
            V = cluster2[features].values

            if distance_metric == 'euclidean': #Calculates, if selected, Euclidiean distance between clusters in IOM space and compares to threshold value.
                min_dist = 1000000
                for u in U:
                    for v in V:
                        dist = np.linalg.norm(u-v)
                        if(dist<min_dist):
                            min_dist = dist
                
            elif distance_metric == 'mahalanobis': #Calculates, if selected, Mahalanobis distance between clusters.
                min_dist = np.sqrt(jit_md(np.mean(U, axis=0), np.mean(V, axis=0),
                                         np.cov(U.T), np.cov(V.T)))

            dist_matrix[counter]=min_dist #Assigns determined distance to relevant distance matrix entry.
            counter+=1 #Iterates to next element.
    
    np.save(f'{save_dir}/{run_name}_DistanceMatrix_{distance_metric}.npy',dist_matrix) #Exports generated distance matrix as array stored in .npy file.

    return dist_matrix

#Function that completes KS testing on small populations        
def small_num_KS(feh1,feh2,Nstars_lim=20,NMC=100):
    
    main_feh,test_feh=feh1,feh2 #Loads and orders metallicity data of each cluster by size.
    if len(main_feh)<len(test_feh): 
        main_feh, test_feh = feh2,feh1

    KS, pval = stats.ks_2samp(main_feh,test_feh,mode='exact') #Conducts 2-sample KS test on clusters' data.

    if len(test_feh)<Nstars_lim: #Detects if smaller cluster is under the stellar numerical threshold.
        orig_pval = pval #Assigns probability value only for clusters below threshold.
    else:
        orig_pval = pval #Assigns probability value and 'greater than' state percentage values to clusters above threshold.
        comp_percentage = 999.9
        less_percentage = 999.9

    if len(test_feh)<Nstars_lim: #Re-run KS testing on random samples of data if cluster is below threshold.            
        NNs,pvals=[],[] #Initialises arrays for storing KS test results.
        count=0 #Iteration tracker to ensure correct number of random sample test are conducted.
        while count<NMC: #Iterate until correct number of random sample tests are complete.
            random_sample_feh=np.random.choice(main_feh,len(test_feh)) #Select random sample from larger dataset
            KS,pval=stats.ks_2samp(main_feh, random_sample_feh, mode='exact') #Conduct KS test on random sample.
            NNs.append(count) #Store KS test data in external arrays.
            pvals.append(pval)
            count+=1 #Iterate to next random sample.
        
        less_coords=np.where((np.array(pvals)<orig_pval)) #Calculate percentage of random sample KS tests with probability values lower than that of the test cluster.
        less_percentage=float(len(np.array(pvals)[less_coords]))/float(len(np.array(pvals)))
        
        same_coords=np.where((np.array(pvals)>0.05)) #Calculate percentage of random sample KS tests with probability values higher than probability threshold.
        same_distr=100.0*float(len(np.array(pvals)[same_coords]))/float(len(np.array(pvals)))        
        
        diff_coords=np.where((np.array(pvals)<0.05)) #Calculate percentage of random sample KS tests with probability values lower than probability threshold.
        diff_distr=100.0*float(len(np.array(pvals)[diff_coords]))/float(len(np.array(pvals)))        
        
        if orig_pval<0.05: #Return relevant percentage of random tests dependent on test cluster KS test probability value.
            comp_percentage=diff_distr
        else:
            comp_percentage=same_distr
    
        return orig_pval,comp_percentage,less_percentage
    
    