import os
import yaml

#Function that generates the .yaml file to use Run_Script.py for a given clustering run
def generate_yaml_file(filename, base_data, result_folder, potential, run_name, 
                       path='yaml_files', 
                       vtoomre=True, 
                       skip_steps=[],N_art=100,
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
                    'base_data': QuotedString('/cosma8/data/dp262/dc-dodd1/%s.hdf5' % base_data),
                    'result_folder': QuotedString('/cosma8/data/dp262/dc-dodd1/%s/' % result_folder),
                    'pot_name': QuotedString('/cosma8/data/dp262/dc-dodd1/%s.ini' % potential),
                    'run_name': QuotedString('%s' % run_name),
                    'sample': QuotedString('/cosma8/data/dp262/dc-dodd1/%s.hdf5' % base_data),
                    
                },
                'art':{'N_art':N_art},
                
                'linkage':{'features':FlowList(features),
                           'scales':{feature: FlowList(scales[feature]) for feature in features},
                           'min_members':min_members,
                           'N_processs':N_process
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