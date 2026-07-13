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