from imports import *
import aux_functions as calc

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
            print('Local YAML file directory does not exist. Creating \'yaml_files\' directory for storage.')
            os.mkdir(path) #Create 'yaml_files' directory if it does not exist locally.
            print('Local YAML directory \'yaml_files\' created.')
        else:
            print('Local YAML file directory detected.')
        
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
        yaml_data['data'][to_skip]=step_appends[to_skip] #Ammend yaml_data[data] entries to include skipped step's data.
            
    with open (full_path,'w') as yaml_file: #Save data into .yaml file at requested path.
        yaml.dump(yaml_data,yaml_file, default_flow_style=False, sort_keys=False)

    print('\n.yaml file saved for %s at %s.' % (run_name, full_path))

def cluster(halo,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',iom_scaling=True,yaml_params={}):
    
    lsr_defs=['8kpc','scalelength']

    if os.path.isdir(f'{home_dir}{halo}')!=True:
        print(f'\nNo directory for {halo} in home directory {home_dir}, creating data structure at {home_dir}{halo}.')
        os.makedirs(f'{home_dir}{halo}',exist_ok=True)
        for lsr in lsr_defs:
            os.makedirs(f'{home_dir}{halo}/{lsr}',exist_ok=True)
            os.makedirs(f'{home_dir}{halo}/{lsr}/accreted',exist_ok=True)
            os.makedirs(f'{home_dir}{halo}/{lsr}/vtoomre',exist_ok=True)
    else:
        print(f'\n{halo} directory located.')
    
    if lsr_def not in lsr_defs:
        print('\nInvalid LSR definition selected, please select \'8kpc\' or \'scalelength\'')
        return
    
    if vtoomre==False:
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    run_name=f'{halo}_{lsr_def}_{selection_type}'
    run_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/'
    
    req_data=['base_data.hdf5','potential.ini','sample.hdf5']

    check_data_files=[os.path.exists(f'{run_dir}/{run_name}_{data}') for data in req_data]

    missing_data=[data for data in req_data if check_data_files[req_data.index(data)]!=True]

    if check_data_files!=[True,True,True]:
        print(f'\n{halo} data is incomplete, the following files are missing:\n')
        print(*[' - '+data for data in missing_data],sep='\n')
        print(f'\n Please move/rename files to haloXX_[LSR]_[SELECTION]_[DATA].hdf5/ini in {home_dir}.')
    
    run_data={file.split('.')[0]: run_name+'_'+file.split('.')[0] for file in req_data}

    generate_yaml_file(run_name,run_data['base_data'],'results',run_data['potential'],run_name,run_data['sample'],data_path=run_dir,save_path=run_dir,vtoomre=vtoomre,find_iom_scales=iom_scaling,**yaml_params)
    
    #kc_main(f'yaml_files/{run_name}.yaml')
    