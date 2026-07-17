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

#Function that groups individual clusters from clustering algorithm by both Mahalanobis distance, and chemistry via 2 sample KS tests.
def chemistry_grouping(halo,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',distance_metric='mahalanobis',plots={},dcut=3.2,regen_dist_matrix=False,p_threshold=0.05):

    lsr_defs=['8kpc','scalelength'] #Labels of permitted definitions of Local Solar Neighbourhood ('8kpc' -> proper radial distance from galactic centre, 'scalelength' -> distance that scales based on properties of individual halo).

    if lsr_def not in lsr_defs:
        print('\nInvalid LSR definition selected, please select \'8kpc\' or \'scalelength\'') #Ensure selected Local Solar Neighbourhood is valid.
        return

    if vtoomre==False: #Allocates relevant label to particle subselection method ('accreted' -> particles tagged as accreted by AURIGA, 'vtoomre' -> selection over all particles using Toomre velocity threshold).
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    if os.path.exists(f'{home_dir}{halo}/{lsr_def}/{selection_type}/results/')!=True: #Checks if clustering algorithm has been run on selected dataset.
        print(f'No clustering data detected for {halo} in {home_dir}, generating with default parameters.')
        return
    else:
        print(f'{halo} clustering data located.')

    results_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/results' #Initialises quick access strings for clustering run label and file paths.
    run_name=f'{halo}_{lsr_def}_{selection_type}'

    df=vaex.open(f'{results_dir}/{run_name}_LabelledSample.hdf5') #Loads clustered data set from .hdf5 file.

    with open (f'{home_dir}/iom_scales/{run_name}_sample.json', 'r') as f: #Loads generated IOM scales from external .json file.
        iom_scales = json.load(f)
        f.close()

    features=[feature for feature in iom_scales] #Initialises arrays for each IOM feature and their scaled equivalents.
    scaled_features=[f'scaled_{feature}' for feature in features]

    minmax_values=vaex.from_arrays(**iom_scales) #Generates vaex DataFrame of minimum and maximum extents of each IOM feature for scaling.
    
    scaler=vaex.ml.MinMaxScaler(feature_range=[-1,1],features=features,prefix='scaled_') #Generates 'scaler' object to scale IOM features.
    scaler.fit(minmax_values) #Calibrates scaler to fit selected IOM scales.
    df=scaler.transform(df)  #Applies scaling to each IOM feature across whole dataset.

    sig_df = df.filter('label!=-1').extract() #Subselects all stars in significant clusters.

    sig_df.export(f'{results_dir}/{run_name}_SignificantSample.hdf5') #Saves significant subselection as external file.
    
    unique_labels= np.unique(sig_df.evaluate('label')) #Retrieves all unique significant cluster labels.
    N_unique=len(unique_labels) #Total number of significant clusters.

    cmap=plt.get_cmap('gist_ncar',N_unique) #Generates colour map with distinct inidividual colour for each cluster.

    clusters_cmap, clusters_norm = colors.from_levels_and_colors(unique_labels,[cmap(i) for i in range(cmap.N)],extend='max') #Maps and normalises colour map to unique cluster labels.

    os.makedirs(f'{results_dir}/plotting',exist_ok=True) #Creates sub-directory for data files associated with plotting (exported colour maps etc.).

    with open(f'{results_dir}/plotting/clusters_cmap.pkl','wb') as f: #Exports generated unique cluster colour map to external .pkl file.
        pickle.dump({'cmap':clusters_cmap,'norm':clusters_norm},f)

    if 'raw' in plots and plots['raw'] == True: #Plots raw IOM cluster plot if requested.
        plot.clusters(halo,lsr_def=lsr_def,vtoomre=vtoomre,home_dir=home_dir,cluster_by='raw')
    
    distance_matrix=calc.cluster_distance_matrix(sig_df,halo,lsr_def=lsr_def,vtoomre=vtoomre,home_dir=home_dir,distance_metric=distance_metric,regen=regen_dist_matrix) #Generates, if requested, distance matrix using chosen distance metric (Euclidean or Mahalanobis) for dataset.

    single_linkage=linkage(distance_matrix, 'single') #Runs single linkage on resulting clusters from single linkage.
    np.save(f'{results_dir}/{run_name}_SingleLinkage_{distance_metric}.npy',single_linkage) #Saves single linkage output as external .npy file.

    if 'cluster_dendrogram' in plots and plots['cluster_dendrogram']==True: #Plots grouped dendrogram if requested.
        plot.cluster_dendrogram(halo,lsr_def=lsr_def,vtoomre=vtoomre,home_dir=home_dir,distance_metric=distance_metric,dcut=dcut)

    if os.path.exists(f'{results_dir}/{run_name}_Leaves.npy') !=True or os.path.exists(f'{results_dir}/plotting/leaf_colours.npy')!=True: #Generates (via plotting) necessary dendrogram data if not currently present.
        plot.cluster_dendrogram(halo,lsr_def=lsr_def,vtoomre=vtoomre,home_dir=home_dir,distance_metric=distance_metric,dcut=dcut)

    grouped_clusters=fcluster(single_linkage,dcut,criterion='distance') #Separates single-linked clusters into groups using requested cutoff distance threshold.

    leaves=np.load(f'{results_dir}/{run_name}_Leaves.npy') #Loads dendrogram leaf data and colour values.
    leaf_colours=np.load(f'{results_dir}/plotting/leaf_colours.npy')
    
    cluster_mapping = {i: cluster for i, cluster in enumerate(grouped_clusters)} #Assigns each cluster to a dictionary entry that stores it's allocated group.
    
    groups=np.array([cluster_mapping[leaf] for leaf in leaves]) #Assigns each clusters' allocated group to entry in numpy array for all clusters.
    
    unique_groups=np.unique(groups) #Retrieves unique list of groups from dendrogram across clusters.
    group_ids=np.ones(df.count())*-1 #Initialises group labels array such that all stars start as 'ungrouped'.
    colour_ids = np.empty(df.count(),dtype=object) #Initialises equivalent array to store the allocated colour of each star's dendrogram group.
    for i in range(len(leaves)):
        leaf_index = np.where(df.evaluate('label')== int(leaves[i])) #Finds all stars in a given cluster.
        group_ids[leaf_index] = int(groups[i]) #Assigns correct dendrogram group ID to all stars in a given gluster.
        colour_ids[leaf_index] =leaf_colours[i] #Assigns relevant group colour to all stars in a given cluster.
    
    df['groups'] = group_ids #Stores groups and group colours as columns in DataFrame.
    df['colours'] = colour_ids

    df.export(f'{results_dir}/{run_name}_GroupedSample.hdf5') #Exports updated dataset with dendrogram groups to external .hdf5 file.

    if 'groups' in plots and plots['groups'] == True: #Plots dendrogram-grouped IOM plot if requested.
        plot.clusters(halo,lsr_def=lsr_def,vtoomre=vtoomre,home_dir=home_dir,cluster_by='groups')
    
    clusters_per_group=np.array([len(np.unique(df.evaluate('label',selection='groups==%s' % group))) for group in unique_groups]) #Calculates number of distinct clusters in each dendrogram group.
    
    sorted_indexes=np.flip(np.argsort(clusters_per_group)) #Sorts groups and clusters per group by number of clusters per group.
    sorted_groups=unique_groups[sorted_indexes]
    sorted_clusters_per_group=clusters_per_group[sorted_indexes]

    N_original_clusters=len(single_linkage)+1 #Calculates number of original clusters before chemistry KS testing and merging
    max_original_clusters=len(single_linkage) #Stores index of last original cluster.

    clusters={i: [i] for i in range(N_original_clusters)} #Initialises dictionary of unmerged clusters to store which original cluster labels are in a KS-merged cluster.

    KStest_data={} #Initialised dictionary to store necessary KS test data for future exporting.
    test_index=1 #Tracker index to distinguish KS test when storing data.

    stopped_branches=set() #Initialises empty set of branches that are no longer considered in KS-merging.

    cluster_index=N_original_clusters #Assigns index of new merged cluster to next available non-original cluster.

    for i, (cluster1,cluster2,dist,size) in enumerate(single_linkage): #Traverses dendrogram to apply KS tests.
        
        cluster1,cluster2=int(cluster1),int(cluster2) #Loads clusters to test and index of future merged cluster.
        new_index=cluster_index 

        if cluster1 in stopped_branches or cluster2 in stopped_branches: #Detects if merging has been stopped for either of selected clusters, and skips step if so.
            stopped_branches.add(new_index)
            cluster_index+=1
            continue
        
        labels1=clusters[cluster1] #Retrieves labels of all original clusters in selected test clusters.
        labels2=clusters[cluster2]   

        feh1 = df['feh'].values[np.isin(df.evaluate('label'),labels1)] #Retrieves and cleans metallicity data for selected clusters.
        feh1 = feh1[~np.isnan(feh1)]
        
        feh2 = df['feh'].values[np.isin(df.evaluate('label'),labels2)]
        feh2 = feh2[~np.isnan(feh2)]

        if (len(feh1)>5)&(len(feh2)>5):  #Ensures that there are sufficient stars in each cluster to perform KS test.
            if (len(feh1)<20)|(len(feh2)<20):  #Detects if small population KS test is necessary.
                orig_pval, comp_percentage, less_percentage=calc.small_num_KS(feh1,feh2,20,NMC=100) #Completes small populations KS test on data.
                if (orig_pval < 0.05) & (comp_percentage>80): #Checks for failure case in small population KS test.
                    KS =-9999 #Assigns failure state variables.
                    pval = -999
                else:  #Runs full scale KS test on data if small number KS test is successful.
                    res= stats.ks_2samp(feh1, feh2, mode='auto')
                    KS, pval,statistic_location, statistic_sign = getattr(res, 'statistic'),getattr(res, 'pvalue'),getattr(res, 'statistic_location'),getattr(res, 'statistic_sign') #Assigns relevant KS test values to variables.
            
            else: #Runs full scale KS test on dataset if populations are sufficiently large.   
                res= stats.ks_2samp(feh1, feh2, mode='auto')
                KS, pval,statistic_location, statistic_sign = getattr(res, 'statistic'),getattr(res, 'pvalue'),getattr(res, 'statistic_location'),getattr(res, 'statistic_sign') #Assigns relevant KS test values to variables.
        else: #Assign failure state variables if there are an insufficient number of stars to perform any KS testing.
            KS =-9999 
            pval = 9999

        KStest_data[f'test{test_index}']={'cluster1':cluster1,'cluster2':cluster2,'labels1':labels1,'labels2':labels2,'KS':float(KS),'pval':float(pval),'stat_loc':float(statistic_location)} #Stores tested clusters and KS values of current test in dictionary entry.
        test_index+=1 #Iterates to next KS test index.

        if KS == -9999: #Detects KS test failure state
            if (len(feh1) < 5)|(len(feh1)<20): #Detects if failure was due to cluster size of either selected cluster and progresses opposite cluster.
                clusters[new_index] = clusters[cluster2]
            elif (len(feh2) < 5)|(len(feh2)<20):
                clusters[new_index] = clusters[cluster1]
    
            cluster_index=cluster_index+1 #Iterates to next cluster.
            continue  

        if pval < p_threshold: #Stops merging along cluster branch should KS test fail due to low probability and iterates to next cluster.
            stopped_branches.add(new_index)
            cluster_index=cluster_index+1
            continue 

        clusters[new_index] = clusters[cluster1] + clusters[cluster2] #Assigns merged cluster to new index if KS test passes.
        cluster_index=cluster_index+1

    original_clusters=np.arange(0,len(single_linkage)+1,1) #Generates list of original cluster labels.
    chemical_groups=original_clusters.copy() #Initialises chemical group labels based on original cluster labels.

    for cluster in original_clusters: #Iterates over original clusters to assign chemical groups.
        for index, step in enumerate(dict(reversed(list(clusters.items())))): #Iterates over merged clusters to detect final merged clusters.
            for f in clusters[step]: #Searches for original cluster in merged clusters.
                if cluster == f: #Assigns chemical group label to cluster if detected in a given merged cluster.
                    chemical_groups[cluster] = index
                    break
            if cluster==f: #Stops loop if an unmerged cluster is detected.
                break
    
    unique_chemical_groups=np.unique(chemical_groups) #Generates list of unique chemical groups.
    N_unique_chemgroups=len(unique_chemical_groups) #Calculates total number of unique chemical groups.

    label_mapping = {old: new for new, old in enumerate(unique_chemical_groups, start=0)} #Assigns reordered labels for each chemical group.
    relabelled_chemical_groups = np.array([label_mapping[label] for label in chemical_groups]) #Generates array of new reordered chemical group labels for each cluster.

    KS_groups=np.zeros(df.count())+-1.0 #Initialises array of each star's chemical group with stars starting as ungrouped.

    labels = df.evaluate('label') #Loads array of original cluster labels for every star.

    sorted_unique_labels = np.sort(unique_labels) #Generates array of ordered unique cluster labels.
    
    for cluster in range(len(sorted_unique_labels)): #Iterates over every original cluster to assign KS groups.
        KS_groups[labels==sorted_unique_labels[cluster]]=relabelled_chemical_groups[cluster] #Assigns cluster's KS group to all stars in a given original cluster.

    df['KS_groups']=np.array(KS_groups) #Stores KS groups in DataFrame
    
    df.export(f'{results_dir}/{run_name}_ChemistryGroups.hdf5') #Exports chemistry group data to separate external .hdf5 file.

    KS_df = df.filter('KS_groups!=-1').extract() #Generates subsample of all stars assigned to chemical groups.
    
    unique_chem_labels = np.unique(KS_df.evaluate('KS_groups')) #Generates list of unique chemical groups.
    N_unique_chem=len(unique_chem_labels) #Calculates total number of unique chemical groups.

    chem_cmap=plt.get_cmap('gist_ncar',N_unique_chem) #Generates colour map with distinct inidividual colour for each chemical group.

    chemistry_cmap, chemistry_norm = colors.from_levels_and_colors(unique_chem_labels,[chem_cmap(i) for i in range(chem_cmap.N)],extend='max') #Maps and normalises colour map to unique KS group labels.

    with open(f'{results_dir}/plotting/chemistry_cmap.pkl','wb') as f: #Exports generated unique KS group colour map to external .pkl file.
        pickle.dump({'cmap':chemistry_cmap,'norm':chemistry_norm},f)

    with open (f'{results_dir}/{run_name}_KSTests.json', 'w') as f: #Exports all KS test data to external .json file.
        json.dump(KStest_data,f)
        f.close()

    if 'KS_groups' in plots and plots['KS_groups'] == True: #Plots chemically-grouped IOM cluster plot if requested.
        plot.clusters(halo,lsr_def=lsr_def,vtoomre=vtoomre,home_dir=home_dir,cluster_by='chemistry')

    if 'chem_dendrogram' in plots and plots['chem_dendrogram'] == True: #Plots chemically-grouped dendrogram if requested.
        plot.cluster_dendrogram(halo,lsr_def=lsr_def,vtoomre=vtoomre,home_dir=home_dir,distance_metric=distance_metric,dcut=dcut, show_chem=True)