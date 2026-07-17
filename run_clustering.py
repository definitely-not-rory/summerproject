from imports import *
import aux_functions as calc
import plot_generation as plot

#Function that generates the .yaml file to use Run_Script.py for a given clustering run
def generate_yaml_file(filename, base_data, result_folder, potential, run_name, sample, 
                       path='yaml_files', 
                       data_path='/cosma8/data/dp262/dc-dodd1/',
                       save_path='/cosma8/data/dp262/dc-dodd1/',
                       vtoomre=False,
                       find_iom_scales=True, 
                       skip_steps=[],
                       N_art=100,
                       scales={'En':[-200000,0],'Lperp':[0,5000],'Lz':[-5000,5000]},
                       min_members=10,
                       N_process=4,
                       min_sig=3,
                       dendo_cut=3.2):
    
    if path!='yaml_files': #Locate alternative .yaml save location if specified, if not specified, default to a local 'yaml_files' directory.
        if os.path.isdir(path)!=True: #Verify alternative path exists.
            print('%s does not exist.' % path)
            return
        else:
            print('Alternative path selected: %s' % path)
    else:
        if os.path.isdir(path)!=True: #If default directory selected, detect if it exists.
            print('\nLocal YAML file directory does not exist. Creating \'yaml_files\' directory for storage.')
            os.mkdir(path) #Create 'yaml_files' directory if it does not exist locally.
            print('\nLocal YAML directory \'yaml_files\' created.')
        else:
            print('\nLocal YAML file directory detected.')
        
    current_files=[file for file in os.listdir(path) if file[-5:]=='.yaml'] #Read all existing .yaml files stored in specified directory.

    if filename+'.yaml' in current_files: #Detect if requested .yaml filename is already in use.
        valid_input=False #Yes/No input auxilliary variables.
        valid_inputs=['y','Y','n','N']
        
        while valid_input==False: #Only progress function until Yes/No receives input of 'Y', 'y', 'N' or 'n'.
            overwrite=input('%s.yaml already exists in %s. Do you want to overwrite it (y/n)?:' % (filename, path)) #Request permission to overwrite pre-existing file.
            if overwrite in valid_inputs: #Break loop if a valid input is received.
                valid_input=True
                if overwrite in valid_inputs[2:]: #Check if overwrite permission is 'No', and halt function if it is.
                    return
                else:
                    print('\nOverwriting %s.yaml' % filename)

    full_path='%s/%s.yaml' % (path,filename) #Generate full file path for saving final .yaml file.

    steps=['artificial_vtoomre_sim','linkage','significance','clusters','cluster_fits', 'cluster_dendogram'] #Default lists of steps and parameters used in clustering.
    features=['En','Lperp','Lz']
    
    if find_iom_scales==True:
        scales=calc.get_iom_scales(sample,path=data_path)

    if vtoomre==False: #Check whether or not clustering is to be run with vtoomre cut, and select correct step.
        steps[0]='artificial_no_vtoomre'

    class FlowList(list): #Class and YAML representers to ensure lists and strings appear correctly in final YAML file in the form [a,b,c,d] and "string".
        pass

    def flow_list_representer(dumper, data):
        return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)

    yaml.add_representer(FlowList, flow_list_representer)

    class QuotedString(str):
        pass

    def quoted_str_representer(dumper, data):
        return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='"')

    yaml.add_representer(QuotedString, quoted_str_representer)

    #Default .yaml save data assuming no steps are skipped.
    yaml_data={'run':{
                    'steps':FlowList(steps)
                    },
                
                'data':{
                    'base_data': QuotedString('%s%s.hdf5' % (data_path,base_data)),
                    'result_folder': QuotedString('%s%s/' % (save_path,result_folder)),
                    'pot_name': QuotedString('%s%s.ini' % (data_path,potential)),
                    'run_name': QuotedString('%s' % run_name),
                    'sample': QuotedString('%s%s.hdf5' % (data_path,sample)),
                    
                },
                'art':{'N_art':N_art},
                
                'linkage':{'features':FlowList(features),
                           'scales':{feature: FlowList(scales[feature]) for feature in features},
                           'min_members':min_members,
                           'N_process':N_process
                           },
                
                'cluster':{'min_sig':min_sig},

                'group':{'dendo_cut':dendo_cut}
    }

    step_labels=['art','linkage','sig','clusters','gaussian_fits','labelled_sample'] #User input and yaml_data[data] labels used to select which steps to skip.
    
    step_appends={'art': '+"Art.hdf5"', 
                 'linkage': '+"Linkage.hdf5"', 
                 'sig': '+"Significance.hdf5"', 
                 'clusters': '+"Clusters.hdf5"', 
                 'gaussian_fits': '+"GaussianFits.hdf5"', 
                 'labelled_sample': '+"LabelledSample.hdf5"'} #Amendments made to yaml_data[data] if a step is selected to be skipped.
     
    for to_skip in skip_steps: #Modify yaml_data to skip requested steps.
        step_index=step_labels.index(to_skip) #Locate requested step.
        yaml_data['run']['steps'].pop(step_index) #Remove step from yaml_data[run][steps] and global comparison lists.
        step_labels.pop(step_index)
        yaml_data['data'][to_skip]=step_appends[to_skip] #Amend yaml_data[data] entries to include skipped step's data.
            
    with open (full_path,'w') as yaml_file: #Save data into .yaml file at requested path.
        yaml.dump(yaml_data,yaml_file, default_flow_style=False, sort_keys=False)

    print('\n.yaml file saved for %s at %s.' % (run_name, full_path))

#Function that runs the KapteynClustering algorithm for a specified dataset.
def cluster(halo,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',iom_scaling=True,yaml_params={}): 
    
    lsr_defs=['8kpc','scalelength'] #Labels of permitted definitions of Local Solar Neighbourhood ('8kpc' -> proper radial distance from galactic centre, 'scalelength' -> distance that scales based on properties of individual halo).

    if os.path.isdir(f'{home_dir}{halo}')!=True: #Verifies that data directories for clustering algorithm exist.
        print(f'\nNo directory for {halo} in home directory {home_dir}, creating data structure at {home_dir}{halo}, please populate with correct data files.')
        os.makedirs(f'{home_dir}{halo}',exist_ok=True) #Generates halo directory and relevant subdirectories.
        for lsr in lsr_defs:
            os.makedirs(f'{home_dir}{halo}/{lsr}',exist_ok=True)
            os.makedirs(f'{home_dir}{halo}/{lsr}/accreted',exist_ok=True)
            os.makedirs(f'{home_dir}{halo}/{lsr}/vtoomre',exist_ok=True)
        return
    else:
        print(f'\n{halo} directory located.')
    
    if lsr_def not in lsr_defs:
        print('\nInvalid LSR definition selected, please select \'8kpc\' or \'scalelength\'') #Ensure selected Local Solar Neighbourhood is valid.
        return
    
    if vtoomre==False: #Allocates relevant label to particle subselection method ('accreted' -> particles tagged as accreted by AURIGA, 'vtoomre' -> selection over all particles using Toomre velocity threshold).
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    run_name=f'{halo}_{lsr_def}_{selection_type}' #Initialises quick access strings for clustering run label and file paths.
    run_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/'
    
    req_data=['base_data.hdf5','potential.ini','sample.hdf5'] #Valid suffixes of data files needed to run clustering algorithm.

    check_data_files=[os.path.exists(f'{run_dir}/{run_name}_{data}') for data in req_data] #Checks if each necessary data file for a given run is present.

    missing_data=[data for data in req_data if check_data_files[req_data.index(data)]!=True] #Generates list of any missing data in chosen sample directory.

    if check_data_files!=[True,True,True]: #Checks if all three required files are present.
        print(f'\n{halo} data is incomplete, the following files are missing:\n') #Indicates which, if any, data files are missing.
        print(*[' - '+data for data in missing_data],sep='\n')
        print(f'\n Please move/rename files to haloXX_[LSR]_[SELECTION]_[DATA].hdf5/ini in {home_dir}.')
    
    run_data={file.split('.')[0]: run_name+'_'+file.split('.')[0] for file in req_data} #Dictionary of each type of data file (base_data, sample, potential) and their names.

    for file in ['base_data','sample']: #Verifies all necessary columns are present in each particle data file.
        df=vaex.open(f'{run_dir}{run_data[file]}.hdf5')

        columns=df.get_column_names() #Accesses all current columns in a given data file.
        req_cols=['pos','vel','phi','vT','vR','R'] #List of necessary columns for whom who's presence in each file needs verification.
        missing_cols=[col for col in req_cols if col not in columns] #Generates list of any missing columns in file.

        if 'pos' in missing_cols: #Generates, if necessary, column for Cartesian position as a 3-element array from individual axis component columns.
            df['pos']=np.column_stack((df.evaluate('x'),df.evaluate('y'),df.evaluate('z')))

        if 'vel' in missing_cols: #Generates, if necessary, column for Cartesian velocity as a 3-element array from individual axis component columns.
            df['vel']=np.column_stack((df.evaluate('vx'),df.evaluate('vy'),df.evaluate('vz')))

        if 'phi' in missing_cols: #Generates, if necessary required cylindrical coordinates and velocities from Cartesian coordinates.
            x=df.evaluate('x')
            y=df.evaluate('y')

            vx=df.evaluate('vx')
            vy=df.evaluate('vy')

            df['R']=np.sqrt((x*x)+(y*y))
            df['phi']=np.arctan2(y,x)

            df['vR']=((x*vx)+(y*vy))/df.evaluate('R')
            df['vT']=-((x*vy)-(y*vx))/df.evaluate('R')

        df.export(f'{run_dir}{run_data[file]}_updated.hdf5') #Exports updated dataset (with added missing columns) to a temporary file.
        del df #Releases control of the original .hdf5 file so it can be overwritten.
 
        os.replace(f'{run_dir}{run_data[file]}_updated.hdf5', f'{run_dir}{run_data[file]}.hdf5') #Overwrites original .hdf5 with updated version and removes temporary file.

    generate_yaml_file(run_name,run_data['base_data'],'results',run_data['potential'],run_name,run_data['sample'],data_path=run_dir,save_path=run_dir,vtoomre=vtoomre,find_iom_scales=iom_scaling,**yaml_params) #Generates .yaml file for selected dataset using provided directory paths, whether or not to regenerate the IOM scales, and ncludes any additional requested parameters within the YAML file
 
    kc_main(f'yaml_files/{run_name}.yaml') #Runs clustering algorithm for requested dataset using selected .yaml file.


def chemistry_grouping(halo,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',distance_metric='mahalanobis',plots={},dcut=3.2,p_threshold=0.05):

    lsr_defs=['8kpc','scalelength']

    if lsr_def not in lsr_defs:
        print('\nInvalid LSR definition selected, please select \'8kpc\' or \'scalelength\'')
        return

    if vtoomre==False:
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    if os.path.exists(f'{home_dir}{halo}/{lsr_def}/{selection_type}/results/')!=True:
        print(f'No clustering data detected for {halo} in {home_dir}, generating with default parameters.')

    else:
        print(f'{halo} clustering data located.')

    results_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/results'
    run_name=f'{halo}_{lsr_def}_{selection_type}'

    df=vaex.open(f'{results_dir}/{run_name}_LabelledSample.hdf5')

    with open (f'{home_dir}/iom_scales/{run_name}_sample.json', 'r') as f:
        iom_scales = json.load(f)
        f.close()

    features=[feature for feature in iom_scales]
    scaled_features=[f'scaled_{feature}' for feature in features]

    minmax_values=vaex.from_arrays(**iom_scales)
    
    scaler=vaex.ml.MinMaxScaler(feature_range=[-1,1],features=features,prefix='scaled_')
    scaler.fit(minmax_values)
    df=scaler.transform(df)

    sig_df = df.filter('label!=-1').extract()

    sig_df.export(f'{results_dir}/{run_name}_SignificantSample.hdf5')
    
    unique_labels= np.unique(sig_df.evaluate('label'))
    N_unique=len(unique_labels)

    cmap=plt.get_cmap('gist_ncar',N_unique)

    clusters_cmap, clusters_norm = colors.from_levels_and_colors(unique_labels,[cmap(i) for i in range(cmap.N)],extend='max')

    os.makedirs(f'{results_dir}/plotting',exist_ok=True)

    with open(f'{results_dir}/plotting/clusters_cmap.pkl','wb') as f:
        pickle.dump({'cmap':clusters_cmap,'norm':clusters_norm},f)

    if 'raw' in plots and plots['raw'] == True:
        plot.clusters(halo,lsr_def=lsr_def,vtoomre=vtoomre,home_dir=home_dir,cluster_by='raw')
    
    distance_matrix=calc.cluster_distance_matrix(sig_df,halo,lsr_def=lsr_def,vtoomre=vtoomre,home_dir=home_dir,distance_metric=distance_metric)

    single_linkage=linkage(distance_matrix, 'single')
    np.save(f'{results_dir}/{run_name}_SingleLinkage_{distance_metric}.npy',single_linkage)

    if 'cluster_dendrogram' in plots and plots['cluster_dendrogram']==True:
        plot.cluster_dendrogram(halo,lsr_def=lsr_def,vtoomre=vtoomre,home_dir=home_dir,distance_metric=distance_metric,dcut=dcut)

    if os.path.exists(f'{results_dir}/{run_name}_Leaves.npy') !=True or os.path.exists(f'{results_dir}/plotting/leaf_colours.npy')!=True:
        plot.cluster_dendrogram(halo,lsr_def=lsr_def,vtoomre=vtoomre,home_dir=home_dir,distance_metric=distance_metric,dcut=dcut)

    grouped_clusters=fcluster(single_linkage,dcut,criterion='distance')

    leaves=np.load(f'{results_dir}/{run_name}_Leaves.npy')
    leaf_colours=np.load(f'{results_dir}/plotting/leaf_colours.npy')
    
    cluster_mapping = {i: cluster for i, cluster in enumerate(grouped_clusters)}
    
    groups=np.array([cluster_mapping[leaf] for leaf in leaves])
    
    unique_groups=np.unique(groups)
    group_ids=np.ones(df.count())*-1
    colour_ids = np.empty(df.count(),dtype=object)

    for i in range(len(leaves)):
        leaf_index = np.where(df.evaluate('label')== int(leaves[i])) 
        group_ids[leaf_index] = int(groups[i])
        colour_ids[leaf_index] =leaf_colours[i]
    
    df['groups'] = group_ids
    df['colours'] = colour_ids

    df.export(f'{results_dir}/{run_name}_GroupedSample.hdf5')

    if 'groups' in plots and plots['groups'] == True:
        plot.clusters(halo,lsr_def=lsr_def,vtoomre=vtoomre,home_dir=home_dir,cluster_by='groups')
    
    clusters_per_group=np.array([len(np.unique(df.evaluate('label',selection='groups==%s' % group))) for group in unique_groups])
    
    sorted_indexes=np.flip(np.argsort(clusters_per_group))
    sorted_groups=unique_groups[sorted_indexes]
    sorted_clusters_per_group=clusters_per_group[sorted_indexes]

    N_original_clusters=len(single_linkage)+1
    max_original_clusters=len(single_linkage)

    clusters={i: [i] for i in range(N_original_clusters)}

    KStest_data={}
    test_index=1

    stopped_branches=set()

    cluster_index=N_original_clusters

    for i, (cluster1,cluster2,dist,size) in enumerate(single_linkage):
        
        cluster1,cluster2=int(cluster1),int(cluster2)
        new_index=cluster_index

        if cluster1 in stopped_branches or cluster2 in stopped_branches:
            stopped_branches.add(new_index)
            cluster_index+=1
            continue
        
        labels1=clusters[cluster1]
        labels2=clusters[cluster2]   

        feh1 = df['feh'].values[np.isin(df.evaluate('label'),labels1)]
        feh1 = feh1[~np.isnan(feh1)]
        
        feh2 = df['feh'].values[np.isin(df.evaluate('label'),labels2)]
        feh2 = feh2[~np.isnan(feh2)]

        if (len(feh1)>5)&(len(feh2)>5): 
            if (len(feh1)<20)|(len(feh2)<20): 
                orig_pval, comp_percentage, less_percentage=calc.small_num_KS(feh1,feh2,20,NMC=100)
                if (orig_pval < 0.05) & (comp_percentage>80): 
                    KS =-9999
                    pval = -999
                else: 
                    res= stats.ks_2samp(feh1, feh2, mode='auto')
                    KS, pval,statistic_location, statistic_sign = getattr(res, 'statistic'),getattr(res, 'pvalue'),getattr(res, 'statistic_location'),getattr(res, 'statistic_sign')
            
            else:    
                res= stats.ks_2samp(feh1, feh2, mode='auto')
                KS, pval,statistic_location, statistic_sign = getattr(res, 'statistic'),getattr(res, 'pvalue'),getattr(res, 'statistic_location'),getattr(res, 'statistic_sign')
        else: 
            KS =-9999 
            pval = 9999

        KStest_data[f'test{test_index}']={'cluster1':cluster1,'cluster2':cluster2,'labels1':labels1,'labels2':labels2,'KS':float(KS),'pval':float(pval),'stat_loc':float(statistic_location)}
        test_index+=1

        if KS == -9999:
            if (len(feh1) < 5)|(len(feh1)<20):
                clusters[new_index] = clusters[cluster2]
            elif (len(feh2) < 5)|(len(feh2)<20):
                clusters[new_index] = clusters[cluster1]
    
            cluster_index=cluster_index+1
            continue  

        if pval < p_threshold:
            stopped_branches.add(new_index)
            cluster_index=cluster_index+1
            continue 

        clusters[new_index] = clusters[cluster1] + clusters[cluster2]
        cluster_index=cluster_index+1

    original_clusters=np.arange(0,len(single_linkage)+1,1)
    chemical_groups=original_clusters.copy()

    for cluster in original_clusters:
        for index, step in enumerate(dict(reversed(list(clusters.items())))):
            for f in clusters[step]: 
                if cluster == f: 
                    chemical_groups[cluster] = index
                    break
            if cluster==f:
                break
    
    unique_chemical_groups=np.unique(chemical_groups)
    N_unique_chemgroups=len(unique_chemical_groups)

    label_mapping = {old: new for new, old in enumerate(unique_chemical_groups, start=0)}
    relabelled_chemical_groups = np.array([label_mapping[label] for label in chemical_groups])

    KS_groups=np.zeros(df.count())+-1.0

    labels = df.evaluate('label')

    sorted_unique_labels = np.sort(unique_labels)
    
    for cluster in range(len(sorted_unique_labels)):
        KS_groups[labels==sorted_unique_labels[cluster]]=relabelled_chemical_groups[cluster]

    df['KS_groups']=np.array(KS_groups)
    
    df.export(f'{results_dir}/{run_name}_ChemistryGroups.hdf5')

    KS_df = df.filter('KS_groups!=-1').extract()
    
    unique_chem_labels = np.unique(KS_df.evaluate('KS_groups'))
    N_unique_chem=len(unique_chem_labels)

    chem_cmap=plt.get_cmap('gist_ncar',N_unique_chem)

    chemistry_cmap, chemistry_norm = colors.from_levels_and_colors(unique_chem_labels,[chem_cmap(i) for i in range(chem_cmap.N)],extend='max')

    with open(f'{results_dir}/plotting/chemistry_cmap.pkl','wb') as f:
        pickle.dump({'cmap':chemistry_cmap,'norm':chemistry_norm},f)

    with open (f'{results_dir}/{run_name}_KSTests.json', 'w') as f:
        json.dump(KStest_data,f)
        f.close()
    