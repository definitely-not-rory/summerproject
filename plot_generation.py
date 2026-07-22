from imports import *

def clusters(halo,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',save_dir='figures',cluster_by='raw',flip_Lz=False,mcps='cluster'):
    
    lsr_defs=['8kpc','scalelength']

    if lsr_def not in lsr_defs:
        print('\nInvalid LSR definition selected, please select \'8kpc\' or \'scalelength\'')
        return

    if vtoomre==False:
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    if os.path.exists(f'{home_dir}{halo}/{lsr_def}/{selection_type}/results/')!=True:
        print(f'No clustering data detected for {halo} in {home_dir}, please run clustering.cluster().')
        return
    else:
        print(f'{halo} clustering data located.')

    
    save_path=f'{save_dir}/{halo}/{lsr_def}/{selection_type}/clusters'
    os.makedirs(save_path,exist_ok=True)

    results_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/results'
    run_name=f'{halo}_{lsr_def}_{selection_type}'
    
    x_axes = ['Lz/10e2', 'Lz/10e2', 'Lperp/10e2', 'vx', 'vx', 'vz']
    y_axes = ['En/10e4', 'Lperp/10e2','En/10e4', 'vy', 'vz', 'vy']

    xlabels = ['$L_z$ [$10^3$ kpc km/s]', '$L_z$ [$10^3$ kpc km/s]', '$L_{\perp}$ [$10^3$ kpc km/s]','$v_x$ [km/s]', '$v_x$ [km/s]', '$v_{z}$ [km/s]']
    ylabels = ['$E$ [$10^5$ km$^2$/s$^2$]', '$L_{\perp}$ [$10^3$ kpc km/s]','$E$ [$10^5$ km$^2$/s$^2$]', '$v_{y}$ [km/s]', '$v_z$ [km/s]', '$v_{y}$ [km/s]']

    df=vaex.open(f'{results_dir}/{run_name}_LabelledSample.hdf5')

    if flip_Lz==True:
        df['Lz']=df.evaluate('Lz')*-1
    
    if cluster_by=='raw':
        clusters_df=vaex.open(f'{results_dir}/{run_name}_SignificantSample.hdf5')

        if flip_Lz==True:
            clusters_df['Lz']=clusters_df.evaluate('Lz')*-1

        with open(f'{results_dir}/plotting/clusters_cmap.pkl','rb') as f:
            cmap_data=pickle.load(f)

        clusters_cmap=cmap_data['cmap']
        clusters_norm=cmap_data['norm']

        fig, ax = plt.subplots(2,3,figsize=[15,10])
        plt.tight_layout()

        for i in range(6):
            plt.sca(ax[int(i/3),i%3])
            
            df.scatter(x_axes[i], y_axes[i], s=0.5, c='silver',alpha=0.7, length_limit=6000000)
            clusters_df.scatter(x_axes[i], y_axes[i], s=1, c=clusters_df.evaluate('label'),cmap=clusters_cmap, norm=clusters_norm, alpha=0.6,length_check=False)

            plt.xlabel(xlabels[i])
            plt.ylabel(ylabels[i])

            plt.tight_layout(w_pad=1)
    
    elif cluster_by=='groups':
        grouped_df=vaex.open(f'{results_dir}/{run_name}_GroupedSample.hdf5')

        if flip_Lz==True:
            grouped_df['Lz']=grouped_df.evaluate('Lz')*-1

        groups_only = grouped_df.filter('label>-1').extract()
        
        colour_list=groups_only.evaluate('colours').to_pylist()

        fig, ax = plt.subplots(2,3,figsize=[15,10])
        plt.tight_layout()

        for i in range(6):
            plt.sca(ax[int(i/3),i%3])
            
            df.scatter(x_axes[i], y_axes[i], s=0.5, c='silver',alpha=0.7, length_limit=6000000)
            groups_only.scatter(x_axes[i], y_axes[i], s=1, c=colour_list, alpha=0.6,length_check=False)

            plt.xlabel(xlabels[i])
            plt.ylabel(ylabels[i])

            plt.tight_layout(w_pad=1)
    
    elif cluster_by=='chemistry':
        grouped_df=vaex.open(f'{results_dir}/{run_name}_ChemistryGroups.hdf5')

        if flip_Lz==True:
            grouped_df['Lz']=grouped_df.evaluate('Lz')*-1

        groups_only = grouped_df.filter('KS_groups>-1').extract()

        with open(f'{results_dir}/plotting/chemistry_cmap.pkl','rb') as f:
            cmap_data=pickle.load(f)

        chemistry_cmap=cmap_data['cmap']
        chemistry_norm=cmap_data['norm']

        fig, ax = plt.subplots(2,3,figsize=[15,10])
        plt.tight_layout()

        for i in range(6):
            plt.sca(ax[int(i/3),i%3])
            
            df.scatter(x_axes[i], y_axes[i], s=0.5, c='silver',alpha=0.7, length_limit=6000000)
            groups_only.scatter(x_axes[i], y_axes[i], s=2, c=groups_only.evaluate('KS_groups'),cmap=chemistry_cmap, norm=chemistry_norm, alpha=0.6,length_check=False)

            plt.xlabel(xlabels[i])
            plt.ylabel(ylabels[i])

            plt.tight_layout(w_pad=1)

    elif cluster_by=='progenitor':
        grouped_df=vaex.open(f'{results_dir}/{run_name}_ChemistryGroups.hdf5')

        if flip_Lz==True:
            grouped_df['Lz']=grouped_df.evaluate('Lz')*-1

        groups_only = grouped_df.filter('label>-1').extract()
        
        with open(f'{results_dir}/plotting/prog_cmap.pkl','rb') as f:
            cmap_data=pickle.load(f)

        prog_cmap=cmap_data['cmap']
        prog_norm=cmap_data['norm']

        prog_cmap.set_under('slategrey')

        fig, ax = plt.subplots(2,3,figsize=[15,10])
        plt.tight_layout()

        for i in range(6):
            plt.sca(ax[int(i/3),i%3])
            
            df.scatter(x_axes[i], y_axes[i], s=0.5, c='silver',alpha=0.7, length_limit=6000000)
            groups_only.scatter(x_axes[i], y_axes[i], s=2, c=(groups_only.evaluate(f'{mcps}_mcp')),cmap=prog_cmap, norm=prog_norm,length_check=False,alpha=0.7)

            plt.xlabel(xlabels[i])
            plt.ylabel(ylabels[i])

            plt.tight_layout(w_pad=1)


    plt.savefig(f'{save_path}/{cluster_by}.pdf')
    plt.savefig(f'{save_path}/{cluster_by}.png',dpi=250,bbox_inches='tight')

    

def cluster_dendrogram(halo,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',save_dir='figures',distance_metric='mahalanobis',dcut=3.2,show_chem=False):
    lsr_defs=['8kpc','scalelength']

    if lsr_def not in lsr_defs:
        print('\nInvalid LSR definition selected, please select \'8kpc\' or \'scalelength\'')
        return

    if vtoomre==False:
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    if os.path.exists(f'{home_dir}{halo}/{lsr_def}/{selection_type}/results/')!=True:
        print(f'No clustering data detected for {halo} in {home_dir}, please run clustering.cluster().')
        return
    else:
        print(f'{halo} clustering data located.')

    results_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/results'
    run_name=f'{halo}_{lsr_def}_{selection_type}'

    save_path=f'{save_dir}/{halo}/{lsr_def}/{selection_type}/dendograms'
    os.makedirs(save_path,exist_ok=True)

    sig_df=vaex.open(f'{results_dir}/{run_name}_SignificantSample.hdf5')    
    single_linkage=np.load(f'{results_dir}/{run_name}_SingleLinkage_{distance_metric}.npy')

    if show_chem==True:
        if os.path.exists(f'{home_dir}{halo}/{lsr_def}/{selection_type}/results/{run_name}_ChemistryGroups.hdf5')!=True:
            print(f'No chemistry data detected for {halo} in {home_dir}, please run clustering.chemistry_grouping().')
            return
        else:
            chem_df=vaex.open(f'{results_dir}/{run_name}_ChemistryGroups.hdf5')

            with open(f'{results_dir}/plotting/chemistry_cmap.pkl','rb') as f:
                cmap_data=pickle.load(f)

            chemistry_cmap=cmap_data['cmap']
            chemistry_norm=cmap_data['norm']


    unique_labels= np.unique(sig_df.evaluate('label'))
    N_unique=len(unique_labels)

    def llf(id):
        return str(int(unique_labels[id]))
    
    fig, ax = plt.subplots(figsize=(30,5))
    dendro_result=dendrogram(single_linkage,leaf_label_func=llf,leaf_rotation=90,color_threshold=dcut)

    leaves = dendro_result["leaves"]
    leaf_colours = dendro_result['leaves_color_list']
    
    np.save(f'{results_dir}/{run_name}_Leaves.npy',leaves)
    np.save(f'{results_dir}/plotting/leaf_colours.npy',leaf_colours)

    icoord = np.array(dendro_result["icoord"])
    dcoord = np.array(dendro_result["dcoord"])

    xstick = np.concatenate((icoord[:,:2],icoord[:,2:]))
    ystick = np.concatenate((dcoord[:,:2],dcoord[:,2:][:,::-1]))
    
    base = ystick[:,0] ==0
    base_xstick = xstick[base,0]
    base_ystick = ystick[base,1]
    ix_sort = np.argsort(base_xstick)
    xbase, ybase = base_xstick[ix_sort], base_ystick[ix_sort]

    with open(f'{results_dir}/plotting/clusters_cmap.pkl','rb') as f:
        cmap_data=pickle.load(f)

    clusters_cmap=cmap_data['cmap']
    clusters_norm=cmap_data['norm']

    for (l,x,y) in zip(leaves,xbase,ybase):
        c = clusters_cmap(clusters_norm(l))
        plt.plot([x,x],[0,0.5],c=c,zorder=10, linewidth=3)

    if show_chem==True:
        KSgrouped_df=chem_df.filter('KS_groups!=-1').extract()

        unique_KSgroups=np.unique(KSgrouped_df.evaluate('KS_groups'))

        for group in unique_KSgroups:
            clusters_in_group=np.unique(KSgrouped_df.evaluate('label',selection='KS_groups==%s'%group))
            
            min_border=np.min([xbase[leaves.index(cluster)] for cluster in clusters_in_group])-3
            max_border=np.max([xbase[leaves.index(cluster)] for cluster in clusters_in_group])+3

            colour=chemistry_cmap(chemistry_norm(group))

            ax.axvline(x=min_border,ls='dashed',alpha=0.5,color=colour,zorder=0)
            ax.axvline(x=max_border,ls='dashed',alpha=0.5,color=colour,zorder=0)
            ax.axvspan(min_border,max_border,alpha=0.2,facecolor=colour,zorder=0)
            lo_y_lim,hi_y_lim=ax.get_ylim()
            
            ax.text((min_border+max_border)/2,lo_y_lim+(hi_y_lim-lo_y_lim)*2/3,f'Group {int(group)}',color=colour,ha='center',va='center_baseline',rotation=-90,size=16,zorder=0)

    plt.ylabel(r'$\rm Mahalanobis \; distance$', fontsize=15)
    plt.xlabel('cluster',size=15)

    if dcut==None:
        filename='uncut'
    else:
        ax.axhline(y=dcut, c='grey', lw=1, linestyle='dashed')
        filename=f'cutoff_{dcut}'
    
    if show_chem==True:
        filename+='_chem'

    plt.savefig(f'{save_path}/{filename}.pdf')
    plt.savefig(f'{save_path}/{filename}.png',dpi=250,bbox_inches='tight')

def KS_tests(halo,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',save_dir='figures',selection=None,distance_metric='mahalanobis',p_threshold=0.05,flip_Lz=False):
    lsr_defs=['8kpc','scalelength']

    if lsr_def not in lsr_defs:
        print('\nInvalid LSR definition selected, please select \'8kpc\' or \'scalelength\'')
        return

    if vtoomre==False:
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    if os.path.exists(f'{home_dir}{halo}/{lsr_def}/{selection_type}/results/')!=True:
        print(f'No clustering data detected for {halo} in {home_dir}, please run clustering.cluster().')
        return
    else:
        print(f'{halo} clustering data located.')

    results_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/results'
    run_name=f'{halo}_{lsr_def}_{selection_type}'

    save_path=f'{save_dir}/{halo}/{lsr_def}/{selection_type}/KS_tests'
    os.makedirs(save_path,exist_ok=True)

    df=vaex.open(f'{results_dir}/{run_name}_ChemistryGroups.hdf5')

    if flip_Lz==True:
        df['Lz']=df.evaluate('Lz')*-1

    with open (f'{results_dir}/{run_name}_KSTests.json', 'r') as f:
        KStests_data = json.load(f)
        f.close()

    single_linkage=np.load(f'{results_dir}/{run_name}_SingleLinkage_{distance_metric}.npy')
    max_original_clusters=len(single_linkage)

    with open(f'{results_dir}/plotting/clusters_cmap.pkl','rb') as f:
            cmap_data=pickle.load(f)

    clusters_cmap=cmap_data['cmap']
    clusters_norm=cmap_data['norm']

    to_load={'feh':'feh','Lz':'Lz/10e2','En':'En/10e4','prog':'progenitor_id'}

    if selection!=None:
        if '=' in selection:
            grouping=selection.split('=')[0]
            group_label=int(selection.split('=')[1])
        else:
            grouping=selection
        
        groupings=['group','KS_group','cluster','failed']

        if grouping not in groupings:
            print('\nChosen sub-selection of data is not available, please select from dendrogram: \'group\', KS Test: \'KS_group\', raw clusters: \'cluster\', or failed tests \'failed\'.')
            return

        if grouping=='cluster':
            req_tests=[]
            
            for test in KStests_data:
                test_data=KStests_data[test]
                
                if group_label in test_data['labels1'] or group_label in test_data['labels2']:
                    req_tests.append(test)
        
        elif grouping=='group':
            req_tests=[]

            clusters_in_group=np.unique(df.evaluate('label',selection='groups==%s'%group_label))
            
            for test in KStests_data:
                test_data=KStests_data[test]
                for cluster in clusters_in_group:
                    if cluster in test_data['labels1'] or cluster in test_data['labels2']:
                        if test not in req_tests:
                            req_tests.append(test)
        
        elif grouping=='KS_group':
            req_tests=[]

            clusters_in_group=np.unique(df.evaluate('label',selection='KS_groups==%s'%group_label))
            
            for test in KStests_data:
                test_data=KStests_data[test]
                for cluster in clusters_in_group:
                    if cluster in test_data['labels1'] or cluster in test_data['labels2']:
                        if test not in req_tests:
                            req_tests.append(test)

        elif grouping=='failed':
            req_tests=[]

            for test in KStests_data:
                test_data=KStests_data[test]

                if test_data['pval']<p_threshold:
                    req_tests.append(test)

        for test in req_tests:
            test_data=KStests_data[test]

            cluster1=test_data['cluster1']
            cluster2=test_data['cluster2']

            labels1=test_data['labels1']
            labels2=test_data['labels2']

            KS=test_data['KS']
            pval=test_data['pval']
            stat_loc=test_data['stat_loc']

            cluster1_data={data: df[to_load[data]].values[np.isin(df.evaluate('label'),labels1)] for data in to_load}
            cluster2_data={data: df[to_load[data]].values[np.isin(df.evaluate('label'),labels2)] for data in to_load}

            for data in [cluster1_data,cluster2_data]:
                data['feh']=data['feh'][~np.isnan(data['feh'])]

            mcp1 = stats.mode(np.array(cluster1_data['prog']), keepdims=True).mode[0]
            mcp2 = stats.mode(np.array(cluster2_data['prog']), keepdims=True).mode[0]
            
            fraction1 = np.sum(np.array(cluster1_data['prog']) == mcp1)/len(cluster1_data['prog'])
            fraction2 = np.sum(np.array(cluster2_data['prog']) == mcp2) / len(cluster2_data['prog'])
            
            cluster1_data['colour'] = clusters_cmap(clusters_norm(cluster1)) if cluster1 <= max_original_clusters else 'b'
            cluster2_data['colour'] = clusters_cmap(clusters_norm(cluster2)) if cluster2 <= max_original_clusters else 'r'
            
            fig, axs = plt.subplots(1,3, figsize=[20,5])
            
            plt.sca(axs[0]) 
            
            df.scatter('Lz/10e2','En/10e4',s=0.5, c='silver',alpha=0.7, length_limit=6000000,label='_None')
            plt.scatter(cluster1_data['Lz'],cluster1_data['En'], s=2, color=cluster1_data['colour'], alpha=0.6,label='Cluster %s'%cluster1)
            plt.scatter(cluster2_data['Lz'],cluster2_data['En'], s=2, color=cluster2_data['colour'], alpha=0.6,label='Cluster %s'%cluster2)
            plt.legend()
            
            plt.sca(axs[1]) 
            
            plt.hist(cluster1_data['feh'],bins=np.arange(-3,1,0.1),histtype='step',color=cluster1_data['colour'],density=True,linestyle='-',label='Clusters %s, N = %s'%(labels1,len(cluster1_data['feh'])))
            plt.plot(np.sort(cluster1_data['feh']),np.linspace(0,1.,len(cluster1_data['feh'])),linestyle='-',color=cluster1_data['colour'])
        
            plt.hist(cluster2_data['feh'],bins=np.arange(-3,1,0.1),histtype='step',color=cluster2_data['colour'],density=True,linestyle='--',label='Clusters %s, N = %s'%(labels2,len(cluster2_data['feh'])))
            plt.plot(np.sort(cluster2_data['feh']),np.linspace(0,1.,len(cluster2_data['feh'])),linestyle='--',color=cluster2_data['colour'])
            
            plt.vlines(stat_loc,0,1.2,linestyle=':')
            plt.legend()

            if pval<0.05:
                plt.title('Test %s - KS: %s, pval: %s'%(test[4:],np.round(KS,2),np.round(pval,3)),fontsize=15, color= 'red', fontweight='bold')
            elif KS==-9999:
                plt.title('Test %s - KS: %s, pval: %s'%(test[4:],np.round(KS,2),np.round(pval,3)),fontsize=15, color= 'b', fontweight='bold')
            else:
                plt.title('Test %s - KS: %s, pval: %s'%(test[4:],np.round(KS,2),np.round(pval,3)),fontsize=15, color= 'green', fontweight='bold')
            plt.xlabel('feh')
            plt.xlim(-3,1)
            plt.ylabel('N')

            plt.sca(axs[2])
    
            plt.hist(cluster1_data['prog'],bins=np.arange(-2,20,1),histtype='step',color=cluster1_data['colour'],density=True,linestyle='-',label = f'MCP = {mcp1}, Fraction = {fraction1:.2f}')
            plt.hist(cluster2_data['prog'],bins=np.arange(-2,20,1),histtype='step',color=cluster2_data['colour'],density=True,linestyle='--',label = f'MCP = {mcp2}, Fraction = {fraction2:.2f}')
            plt.legend()

            plt.savefig(f'{save_path}/{test}_clusters_{cluster1}_{cluster2}.pdf')
            plt.savefig(f'{save_path}/{test}_clusters_{cluster1}_{cluster2}.png',dpi=250,bbox_inches='tight')    
    else:
        for test in KStests_data:
            test_data=KStests_data[test]

            cluster1=test_data['cluster1']
            cluster2=test_data['cluster2']

            labels1=test_data['labels1']
            labels2=test_data['labels2']

            KS=test_data['KS']
            pval=test_data['pval']
            stat_loc=test_data['stat_loc']

            cluster1_data={data: df[to_load[data]].values[np.isin(df.evaluate('label'),labels1)] for data in to_load}
            cluster2_data={data: df[to_load[data]].values[np.isin(df.evaluate('label'),labels2)] for data in to_load}

            for data in [cluster1_data,cluster2_data]:
                data['feh']=data['feh'][~np.isnan(data['feh'])]

            mcp1 = stats.mode(np.array(cluster1_data['prog']), keepdims=True).mode[0]
            mcp2 = stats.mode(np.array(cluster2_data['prog']), keepdims=True).mode[0]
            
            fraction1 = np.sum(np.array(cluster1_data['prog']) == mcp1)/len(cluster1_data['prog'])
            fraction2 = np.sum(np.array(cluster2_data['prog']) == mcp2) / len(cluster2_data['prog'])
            
            cluster1_data['colour'] = clusters_cmap(clusters_norm(cluster1)) if cluster1 <= max_original_clusters else 'b'
            cluster2_data['colour'] = clusters_cmap(clusters_norm(cluster2)) if cluster2 <= max_original_clusters else 'r'
            
            fig, axs = plt.subplots(1,3, figsize=[20,5])
            
            plt.sca(axs[0]) 
            
            df.scatter('Lz/10e2','En/10e4',s=0.5, c='silver',alpha=0.7, length_limit=6000000,label='_None')
            plt.scatter(cluster1_data['Lz'],cluster1_data['En'], s=2, color=cluster1_data['colour'], alpha=0.6,label='Cluster %s'%cluster1)
            plt.scatter(cluster2_data['Lz'],cluster2_data['En'], s=2, color=cluster2_data['colour'], alpha=0.6,label='Cluster %s'%cluster2)
            plt.legend()
            
            plt.sca(axs[1]) 
            
            plt.hist(cluster1_data['feh'],bins=np.arange(-3,1,0.1),histtype='step',color=cluster1_data['colour'],density=True,linestyle='-',label='Clusters %s, N = %s'%(labels1,len(cluster1_data['feh'])))
            plt.plot(np.sort(cluster1_data['feh']),np.linspace(0,1.,len(cluster1_data['feh'])),linestyle='-',color=cluster1_data['colour'])
        
            plt.hist(cluster2_data['feh'],bins=np.arange(-3,1,0.1),histtype='step',color=cluster2_data['colour'],density=True,linestyle='--',label='Clusters %s, N = %s'%(labels2,len(cluster2_data['feh'])))
            plt.plot(np.sort(cluster2_data['feh']),np.linspace(0,1.,len(cluster2_data['feh'])),linestyle='--',color=cluster2_data['colour'])
            
            plt.vlines(stat_loc,0,1.2,linestyle=':')
            plt.legend()

            if pval<0.05:
                plt.title('Test %s - KS: %s, pval: %s'%(test,np.round(KS,2),np.round(pval,3)),fontsize=15, color= 'red', fontweight='bold')
            elif KS==-9999:
                plt.title('Test %s - KS: %s, pval: %s'%(test,np.round(KS,2),np.round(pval,3)),fontsize=15, color= 'b', fontweight='bold')
            else:
                plt.title('Test %s - KS: %s, pval: %s'%(test,np.round(KS,2),np.round(pval,3)),fontsize=15, color= 'green', fontweight='bold')
            plt.xlabel('feh')
            plt.xlim(-3,1)
            plt.ylabel('N')

            plt.sca(axs[2])
    
            plt.hist(cluster1_data['prog'],bins=np.arange(-2,20,1),histtype='step',color=cluster1_data['colour'],density=True,linestyle='-',label = f'MCP = {mcp1}, Fraction = {fraction1:.2f}')
            plt.hist(cluster2_data['prog'],bins=np.arange(-2,20,1),histtype='step',color=cluster2_data['colour'],density=True,linestyle='--',label = f'MCP = {mcp2}, Fraction = {fraction2:.2f}')
            plt.legend()

            plt.savefig(f'{save_path}/{test}_clusters_{cluster1}_{cluster2}.pdf')
            plt.savefig(f'{save_path}/{test}_clusters_{cluster1}_{cluster2}.png',dpi=250,bbox_inches='tight')    

def cluster_only(halo,num,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',save_dir='figures',group='cluster',show_clusters=True,show_KS=False,flip_Lz=False):
    lsr_defs=['8kpc','scalelength']

    if lsr_def not in lsr_defs:
        print('\nInvalid LSR definition selected, please select \'8kpc\' or \'scalelength\'')
        return

    if vtoomre==False:
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    if os.path.exists(f'{home_dir}{halo}/{lsr_def}/{selection_type}/results/')!=True:
        print(f'No clustering data detected for {halo} in {home_dir}, please run clustering.cluster().')
        return
    else:
        print(f'{halo} clustering data located.')

    results_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/results'
    run_name=f'{halo}_{lsr_def}_{selection_type}'

    save_path=f'{save_dir}/{halo}/{lsr_def}/{selection_type}/clusters/'
    
    group_types=['cluster','group','KS_group','prog']

    if group not in group_types:
        print('Please select a valid group type:\n - cluster\n - group\n - KS_group\n - prog')
        return
    
    save_path+=group+'s'

    if group!='cluster':
        if show_clusters==True:
            save_path+='/clusters'
        elif group=='progenitor' and show_KS==True:
            save_path+='/KS_groups'
        else:
            save_path+='/groups'
    
    os.makedirs(save_path,exist_ok=True)

    x_axes = ['Lz/10e2', 'Lz/10e2', 'Lperp/10e2', 'vx', 'vx', 'vz']
    y_axes = ['En/10e4', 'Lperp/10e2','En/10e4', 'vy', 'vz', 'vy']

    xlabels = ['$L_z$ [$10^3$ kpc km/s]', '$L_z$ [$10^3$ kpc km/s]', '$L_{\perp}$ [$10^3$ kpc km/s]','$v_x$ [km/s]', '$v_x$ [km/s]', '$v_{z}$ [km/s]']
    ylabels = ['$E$ [$10^5$ km$^2$/s$^2$]', '$L_{\perp}$ [$10^3$ kpc km/s]','$E$ [$10^5$ km$^2$/s$^2$]', '$v_{y}$ [km/s]', '$v_z$ [km/s]', '$v_{y}$ [km/s]']

    df=vaex.open(f'{results_dir}/{run_name}_ChemistryGroups.hdf5')

    if flip_Lz==True:
        df['Lz']=df.evaluate('Lz')*-1

    with open(f'{results_dir}/plotting/clusters_cmap.pkl','rb') as f:
        cmap_data=pickle.load(f)

    clusters_cmap=cmap_data['cmap']
    clusters_norm=cmap_data['norm']

    with open(f'{results_dir}/plotting/chemistry_cmap.pkl','rb') as f:
        cmap_data=pickle.load(f)

    chemistry_cmap=cmap_data['cmap']
    chemistry_norm=cmap_data['norm']

    with open(f'{results_dir}/plotting/prog_cmap.pkl','rb') as f:
        cmap_data=pickle.load(f)

    prog_cmap=cmap_data['cmap']
    prog_norm=cmap_data['norm']


    if group=='cluster':
        unique_clusters=np.unique(df.evaluate('label'))[1:]

        if num not in unique_clusters:
            print('Cluster number not available.')
            return

        clusters_df=df.filter('label==%s'%num).extract()
    
        fig, ax = plt.subplots(2,3,figsize=[15,10])
        plt.tight_layout()

        for i in range(6):
            plt.sca(ax[int(i/3),i%3])
            
            df.scatter(x_axes[i], y_axes[i], s=0.5, c='silver',alpha=0.7, length_limit=6000000)
            clusters_df.scatter(x_axes[i], y_axes[i], s=1, c=clusters_df.evaluate('label'),cmap=clusters_cmap, norm=clusters_norm, alpha=0.4,length_check=False)

            if int(i/3)<1:
                alpha_shape = alphashape.alphashape(list(zip(clusters_df.evaluate(x_axes[i]), clusters_df.evaluate(y_axes[i]))), alpha=0.2)
                ax[int(i/3),i%3].add_patch(plt.Polygon(list(alpha_shape.exterior.coords),fill=False,edgecolor=clusters_cmap(num),linewidth=1.5,ls='dashed'))

            plt.xlabel(xlabels[i])
            plt.ylabel(ylabels[i])

            plt.tight_layout(w_pad=1)
    
    elif group=='group':
        grouped_df=df.filter('groups!=-1').extract()

        unique_groups=np.unique(grouped_df.evaluate('groups'))

        if num not in unique_groups:
            print('Group number not available.')
            return
        
        if show_clusters==False:
            selected_df=df.filter('groups==%s'%num).extract()

            colour_list=selected_df.evaluate('colours').to_pylist()

            fig, ax = plt.subplots(2,3,figsize=[15,10])
            plt.tight_layout()

            for i in range(6):
                plt.sca(ax[int(i/3),i%3])
                
                df.scatter(x_axes[i], y_axes[i], s=0.5, c='silver',alpha=0.7, length_limit=6000000)
                selected_df.scatter(x_axes[i], y_axes[i], s=1,c=colour_list, alpha=0.4,length_check=False)

                if int(i/3)<1:
                    alpha_shape = alphashape.alphashape(list(zip(selected_df.evaluate(x_axes[i]), selected_df.evaluate(y_axes[i]))), alpha=0.2)
                    ax[int(i/3),i%3].add_patch(plt.Polygon(list(alpha_shape.exterior.coords),fill=False,edgecolor=colour_list[0],linewidth=1.5,ls='dashed'))

                plt.xlabel(xlabels[i])
                plt.ylabel(ylabels[i])

                plt.tight_layout(w_pad=1)
        else:
            selected_df=df.filter('groups==%s'%num).extract()

            colour_list=selected_df.evaluate('colours').to_pylist()

            fig, ax = plt.subplots(2,3,figsize=[15,10])
            plt.tight_layout()

            for i in range(6):
                plt.sca(ax[int(i/3),i%3])
                
                df.scatter(x_axes[i], y_axes[i], s=0.5, c='silver',alpha=0.7, length_limit=6000000)
                selected_df.scatter(x_axes[i], y_axes[i], s=1, c=selected_df.evaluate('label'),cmap=clusters_cmap, norm=clusters_norm, alpha=0.4,length_check=False)
                
                if int(i/3)<1:
                    alpha_shape = alphashape.alphashape(list(zip(selected_df.evaluate(x_axes[i]), selected_df.evaluate(y_axes[i]))), alpha=0.2)
                    ax[int(i/3),i%3].add_patch(plt.Polygon(list(alpha_shape.exterior.coords),fill=False,edgecolor=colour_list[0],linewidth=1.5,ls='dashed'))

                plt.xlabel(xlabels[i])
                plt.ylabel(ylabels[i])

                plt.tight_layout(w_pad=1)
    
    elif group=='KS_group':
        KSgrouped_df=df.filter('KS_groups!=-1').extract()

        unique_KSgroups=np.unique(KSgrouped_df.evaluate('KS_groups'))

        if num not in unique_KSgroups:
            print('Group number not available.')
            return
        
        if show_clusters==False:
            selected_df=df.filter('KS_groups==%s'%num).extract()

            fig, ax = plt.subplots(2,3,figsize=[15,10])
            plt.tight_layout()

            for i in range(6):
                plt.sca(ax[int(i/3),i%3])
                
                df.scatter(x_axes[i], y_axes[i], s=0.5, c='silver',alpha=0.7, length_limit=6000000)
                selected_df.scatter(x_axes[i], y_axes[i], s=2, c=selected_df.evaluate('KS_groups'),cmap=chemistry_cmap, norm=chemistry_norm, alpha=0.4,length_check=False)

                if int(i/3)<1:
                    alpha_shape = alphashape.alphashape(list(zip(selected_df.evaluate(x_axes[i]), selected_df.evaluate(y_axes[i]))), alpha=0.2)
                    ax[int(i/3),i%3].add_patch(plt.Polygon(list(alpha_shape.exterior.coords),fill=False,edgecolor=chemistry_cmap(num),linewidth=1.5,ls='dashed'))

                plt.xlabel(xlabels[i])
                plt.ylabel(ylabels[i])

                plt.tight_layout(w_pad=1)
        else:
            selected_df=df.filter('KS_groups==%s'%num).extract()

            fig, ax = plt.subplots(2,3,figsize=[15,10])
            plt.tight_layout()

            for i in range(6):
                plt.sca(ax[int(i/3),i%3])
                
                df.scatter(x_axes[i], y_axes[i], s=0.5, c='silver',alpha=0.7, length_limit=6000000)
                selected_df.scatter(x_axes[i], y_axes[i], s=1, c=selected_df.evaluate('label'),cmap=clusters_cmap, norm=clusters_norm, alpha=0.4,length_check=False)
                
                if int(i/3)<1:
                    alpha_shape = alphashape.alphashape(list(zip(selected_df.evaluate(x_axes[i]), selected_df.evaluate(y_axes[i]))), alpha=0.2)
                    ax[int(i/3),i%3].add_patch(plt.Polygon(list(alpha_shape.exterior.coords),fill=False,edgecolor=chemistry_cmap(num),linewidth=1.5,ls='dashed'))

                plt.xlabel(xlabels[i])
                plt.ylabel(ylabels[i])

                plt.tight_layout(w_pad=1)

    elif group=='prog':
        prog_df=df.filter('progenitor_id!=-1').extract()

        unique_progs=np.unique(prog_df.evaluate('progenitor_id'))


        if num not in unique_progs:
            print('Group number not available.')
            return
        
        if show_clusters==False:
            if show_KS==False:
                selected_df=df.filter('progenitor_id==%s'%num).extract()

                fig, ax = plt.subplots(2,3,figsize=[15,10])
                plt.tight_layout()

                for i in range(6):
                    plt.sca(ax[int(i/3),i%3])
                    
                    df.scatter(x_axes[i], y_axes[i], s=0.5, c='silver',alpha=0.7, length_limit=6000000)
                    selected_df.scatter(x_axes[i], y_axes[i], s=2, c=selected_df.evaluate('progenitor_id'),cmap=prog_cmap, norm=prog_norm, alpha=0.4,length_check=False)

                    if int(i/3)<1:
                        alpha_shape = alphashape.alphashape(list(zip(selected_df.evaluate(x_axes[i]), selected_df.evaluate(y_axes[i]))), alpha=0.2)
                        ax[int(i/3),i%3].add_patch(plt.Polygon(list(alpha_shape.exterior.coords),fill=False,edgecolor=prog_cmap(prog_norm(num)),linewidth=1.5,ls='dashed'))

                    plt.xlabel(xlabels[i])
                    plt.ylabel(ylabels[i])

                    plt.tight_layout(w_pad=1)
            else:
                selected_df=df.filter('KS_mcp==%s'%num).extract()
                prog_df=df.filter('progenitor_id==%s'%num).extract()

                fig, ax = plt.subplots(2,3,figsize=[15,10])
                plt.tight_layout()

                for i in range(6):
                    plt.sca(ax[int(i/3),i%3])
                    
                    df.scatter(x_axes[i], y_axes[i], s=0.5, c='silver',alpha=0.7, length_limit=6000000)
                    selected_df.scatter(x_axes[i], y_axes[i], s=1, c=selected_df.evaluate('KS_groups'),cmap=chemistry_cmap, norm=chemistry_norm, alpha=0.4,length_check=False)
                    
                    if int(i/3)<1:
                        alpha_shape = alphashape.alphashape(list(zip(prog_df.evaluate(x_axes[i]), prog_df.evaluate(y_axes[i]))), alpha=0.2)
                        ax[int(i/3),i%3].add_patch(plt.Polygon(list(alpha_shape.exterior.coords),fill=False,edgecolor=prog_cmap(prog_norm(num)),linewidth=1.5,ls='dashed'))

                    plt.xlabel(xlabels[i])
                    plt.ylabel(ylabels[i])

                    plt.tight_layout(w_pad=1)
        else:
            selected_df=df.filter('cluster_mcp==%s'%num).extract()
            prog_df=df.filter('progenitor_id==%s'%num).extract()

            fig, ax = plt.subplots(2,3,figsize=[15,10])
            plt.tight_layout()

            for i in range(6):
                plt.sca(ax[int(i/3),i%3])
                
                df.scatter(x_axes[i], y_axes[i], s=0.5, c='silver',alpha=0.7, length_limit=6000000)
                selected_df.scatter(x_axes[i], y_axes[i], s=1, c=selected_df.evaluate('label'),cmap=clusters_cmap, norm=clusters_norm, alpha=0.4,length_check=False)
                
                if int(i/3)<1:
                    alpha_shape = alphashape.alphashape(list(zip(prog_df.evaluate(x_axes[i]), prog_df.evaluate(y_axes[i]))), alpha=0.2)
                    ax[int(i/3),i%3].add_patch(plt.Polygon(list(alpha_shape.exterior.coords),fill=False,edgecolor=prog_cmap(num),linewidth=1.5,ls='dashed'))

                plt.xlabel(xlabels[i])
                plt.ylabel(ylabels[i])

                plt.tight_layout(w_pad=1)

    filename=f'{group}_{num}'

    if show_clusters==True:
        filename+='_clusters'
    elif show_KS==True:
        filename+='_KSgroups'

    plt.savefig(f'{save_path}/{filename}.pdf')
    plt.savefig(f'{save_path}/{filename}.png',dpi=250,bbox_inches='tight')

def prog_FeH_hist(halo,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',save_dir='figures',show_only='all',normalised=True,N_cutoff=10):
    lsr_defs=['8kpc','scalelength']

    if lsr_def not in lsr_defs:
        print('\nInvalid LSR definition selected, please select \'8kpc\' or \'scalelength\'')
        return

    if vtoomre==False:
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    if os.path.exists(f'{home_dir}{halo}/{lsr_def}/{selection_type}/results/')!=True:
        print(f'No clustering data detected for {halo} in {home_dir}, please run clustering.cluster().')
        return
    else:
        print(f'{halo} clustering data located.')

    results_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/results'
    run_name=f'{halo}_{lsr_def}_{selection_type}'

    save_path=f'{save_dir}/{halo}/{lsr_def}/{selection_type}/FeH_hists/progenitors'

    os.makedirs(save_path, exist_ok=True)
    os.makedirs(f'{save_path}/individual', exist_ok=True)
    os.makedirs(f'{save_path}/comparison', exist_ok=True)

    df=vaex.open(f'{results_dir}/{run_name}_ChemistryGroups.hdf5')

    with open(f'{results_dir}/plotting/prog_cmap.pkl','rb') as f:
        cmap_data=pickle.load(f)

    prog_cmap=cmap_data['cmap']
    prog_norm=cmap_data['norm']

    if isinstance(show_only, str)!=True:
        if isinstance(show_only, int)==True:
            progenitors=[show_only]
            
            filename=f'/individual/prog_{show_only}'
        else:
            progenitors=show_only
            
            filename='/comparison/prog'
            for prog in show_only:
                filename+=f'_{prog}'

    else:
        all_progenitors,N_progenitor=np.unique(df.evaluate('progenitor_id'),return_counts=True)
        progenitors=all_progenitors[N_progenitor>=N_cutoff]

        filename='/all_progenitors'

    if len(progenitors)>1:
        threshold_progenitor=np.unique([prog_norm(progenitor) for progenitor in progenitors])[-1]
        alphas=np.ones_like(progenitors)*0.4
        for i in range(threshold_progenitor):
            alphas[i]+=0.3
    else:
        alpha=1  

    plt.figure()
    for progenitor in progenitors:
        prog_df=df.filter('progenitor_id==%s'%progenitor)
        plt.hist(prog_df.filter('feh>-3').extract().evaluate('feh'),label=int(progenitor),histtype='step',bins=np.arange(-3,1,0.1),color=prog_cmap(prog_norm(progenitor)),density=normalised,alpha=alphas[int(progenitor-1)])

    plt.xlabel('FeH')
    if normalised==True:
        plt.ylabel('Normalised N')
    else:
        plt.ylabel('N')
        plt.yscale('log')
    
    plt.savefig(f'{save_path}{filename}.pdf')
    plt.savefig(f'{save_path}{filename}.png',dpi=250,bbox_inches='tight')

def group_FeH_hist(halo,num,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',save_dir='figures',group='cluster',normalised=True):
    lsr_defs=['8kpc','scalelength']

    if lsr_def not in lsr_defs:
        print('\nInvalid LSR definition selected, please select \'8kpc\' or \'scalelength\'')
        return

    if vtoomre==False:
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    if os.path.exists(f'{home_dir}{halo}/{lsr_def}/{selection_type}/results/')!=True:
        print(f'No clustering data detected for {halo} in {home_dir}, please run clustering.cluster().')
        return
    else:
        print(f'{halo} clustering data located.')

    results_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/results'
    run_name=f'{halo}_{lsr_def}_{selection_type}'

    save_path=f'{save_dir}/{halo}/{lsr_def}/{selection_type}/FeH_hists/'
    
    group_types=['cluster','group','KS_group']

    if group not in group_types:
        print('Please select a valid group type:\n - cluster\n - group\n - KS_group')
        return
    
    save_path+=group+'s'
    
    os.makedirs(save_path,exist_ok=True)

    df=vaex.open(f'{results_dir}/{run_name}_ChemistryGroups.hdf5')

    with open(f'{results_dir}/plotting/clusters_cmap.pkl','rb') as f:
        cmap_data=pickle.load(f)

    clusters_cmap=cmap_data['cmap']
    clusters_norm=cmap_data['norm']

    group_colours=df.evaluate('colours').to_pylist()

    with open(f'{results_dir}/plotting/chemistry_cmap.pkl','rb') as f:
        cmap_data=pickle.load(f)

    chemistry_cmap=cmap_data['cmap']
    chemistry_norm=cmap_data['norm']

    columns={'cluster':'label','group':'groups','KS_group':'KS_groups'}
    cmaps={'cluster':{'cmap':clusters_cmap,'norm':clusters_norm},'group':group_colours,'KS_group':{'cmap':chemistry_cmap,'norm':chemistry_norm}}

    group_df=df.filter('%s==%s'%(columns[group],num)).extract()

    if isinstance(cmaps[group],dict)==True:
        colour=cmaps[group]['cmap'](cmaps[group]['norm'](num))
    else:
        colour=cmaps[group][num]

    plt.hist(group_df.filter('feh>-3').extract().evaluate('feh'),histtype='step',bins=np.arange(-3,1,0.1),color=colour,density=normalised)

    plt.xlabel('FeH')
    if normalised==True:
        plt.ylabel('Normalised N')
    else:
        plt.ylabel('N')
        plt.yscale('log')

    plt.savefig(f'{save_path}/{group}_{num}.pdf')
    plt.savefig(f'{save_path}/{group}_{num}.png',dpi=250,bbox_inches='tight')
    
def show_progenitors(halo,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',save_dir='figures'):
    lsr_defs=['8kpc','scalelength']

    if lsr_def not in lsr_defs:
        print('\nInvalid LSR definition selected, please select \'8kpc\' or \'scalelength\'')
        return

    if vtoomre==False:
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    if os.path.exists(f'{home_dir}{halo}/{lsr_def}/{selection_type}/results/')!=True:
        print(f'No clustering data detected for {halo} in {home_dir}, please run clustering.cluster().')
        return
    else:
        print(f'{halo} clustering data located.')

    results_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/results'
    run_name=f'{halo}_{lsr_def}_{selection_type}'

    save_path=f'{save_dir}/{halo}/{lsr_def}/{selection_type}/progenitors/'

    os.makedirs(save_path, exist_ok=True)

    df=vaex.open(f'{results_dir}/{run_name}_ChemistryGroups.hdf5')

    with open(f'{results_dir}/plotting/prog_cmap.pkl','rb') as f:
        cmap_data=pickle.load(f)

    prog_cmap=cmap_data['cmap']
    prog_norm=cmap_data['norm']
    
    total_stars=df.count()
    accreted_stars=df.count(selection='progenitor_id!=-1')
    insitu_stars=df.count(selection='progenitor_id==-1')

    progenitors=df.filter('progenitor_id!=-1').extract().evaluate('progenitor_id')

    unique_progenitors=np.sort(np.unique(progenitors))
    progenitor_labels=[str(progenitor) for progenitor in unique_progenitors]

    stars_per_progenitor=[np.count_nonzero(progenitors==progenitor) for progenitor in unique_progenitors]

    sorted_prog_counts=np.array(sorted([[progenitor,int(stars_per_progenitor[progenitor_labels.index(progenitor)])] for progenitor in progenitor_labels],key=lambda x:x[1],reverse=True))

    plt.figure(figsize=(27,9))
    plt.bar('In-Situ',insitu_stars,color='silver')
    for prog in sorted_prog_counts:
        plt.bar(prog[0],int(prog[1]),color=prog_cmap(prog_norm(int(prog[0]))))
    plt.yscale('log')
    plt.xlabel('Progenitor')
    plt.ylabel('Count')
    plt.xticks(rotation=-90)
    plt.gca().axhline(20,c='red',linestyle='dashed')

    plt.savefig(f'{save_path}counts.pdf')
    plt.savefig(f'{save_path}counts.png',dpi=250,bbox_inches='tight')

def dominance_diagram(halo,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',save_dir='figures',group_by='cluster',N_bins=128,insitu_colour='silver',exclude=[],flip_Lz=False):
    lsr_defs=['8kpc','scalelength']

    if lsr_def not in lsr_defs:
        print('\nInvalid LSR definition selected, please select \'8kpc\' or \'scalelength\'')
        return

    if vtoomre==False:
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    if os.path.exists(f'{home_dir}{halo}/{lsr_def}/{selection_type}/results/')!=True:
        print(f'No clustering data detected for {halo} in {home_dir}, please run clustering.cluster().')
        return
    else:
        print(f'{halo} clustering data located.')

    group_types=['cluster','group','KS_group','prog']

    if group_by not in group_types:
        print('Please select a valid group type:\n - cluster\n - group\n - KS_group\n - prog')
        return

    results_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/results'
    run_name=f'{halo}_{lsr_def}_{selection_type}'

    save_path=f'{save_dir}/{halo}/{lsr_def}/{selection_type}/dominance/'

    os.makedirs(save_path, exist_ok=True)

    df=vaex.open(f'{results_dir}/{run_name}_ChemistryGroups.hdf5')

    if flip_Lz==True:
        df['Lz']=df.evaluate('Lz')*-1

    with open(f'{results_dir}/plotting/clusters_cmap.pkl','rb') as f:
        cmap_data=pickle.load(f)

    clusters_cmap=cmap_data['cmap']
    clusters_norm=cmap_data['norm']

    group_colours=df.evaluate('colours').to_pylist()

    with open(f'{results_dir}/plotting/chemistry_cmap.pkl','rb') as f:
        cmap_data=pickle.load(f)

    chemistry_cmap=cmap_data['cmap']
    chemistry_norm=cmap_data['norm']

    with open(f'{results_dir}/plotting/prog_cmap.pkl','rb') as f:
        cmap_data=pickle.load(f)

    prog_cmap=cmap_data['cmap']
    prog_norm=cmap_data['norm']

    with open (f'{home_dir}/iom_scales/{run_name}_sample.json', 'r') as f: #Loads IOM scales from file for exporting.
        scales = json.load(f)
        f.close()
    
    columns={'cluster':'label','group':'groups','KS_group':'KS_groups','prog':'progenitor_id'}
    cmaps={'cluster':{'cmap':clusters_cmap,'norm':clusters_norm},'group':group_colours,'KS_group':{'cmap':chemistry_cmap,'norm':chemistry_norm},'prog':{'cmap':prog_cmap,'norm':prog_norm}}

    x_axes = ['Lz/10e2', 'Lz/10e2', 'Lperp/10e2']
    y_axes = ['En/10e4', 'Lperp/10e2','En/10e4']

    xlabels = ['$L_z$ [$10^3$ kpc km/s]', '$L_z$ [$10^3$ kpc km/s]', '$L_{\perp}$ [$10^3$ kpc km/s]']
    ylabels = ['$E$ [$10^5$ km$^2$/s$^2$]', '$L_{\perp}$ [$10^3$ kpc km/s]','$E$ [$10^5$ km$^2$/s$^2$]']

    xlims=[np.array(scales[iom.split('/')[0]])/(10**(int(iom.split('/')[1].split('e')[1])+1)) for iom in x_axes]
    ylims=[np.array(scales[iom.split('/')[0]])/(10**(int(iom.split('/')[1].split('e')[1])+1)) for iom in y_axes]

    xbins=[np.linspace(xlims[i][0],xlims[i][1],N_bins) for i in range(3)]
    ybins=[np.linspace(ylims[i][0],ylims[i][1],N_bins) for i in range(3)]

    all_groups=np.unique(df.evaluate('%s'%columns[group_by]))
    groups=[group for group in all_groups if group not in exclude]
    
    N_total=np.zeros((3,N_bins-1,N_bins-1))
    
    for i in range(3):
        x_data=df.evaluate(x_axes[i])
        y_data=df.evaluate(y_axes[i])

        N_total[i-1],_,_=np.histogram2d(x_data,y_data,bins=[xbins[i],ybins[i]])

    N_total[(N_total==0.0)]=1.0

    fig, ax = plt.subplots(1,3,figsize=[15,5])
    plt.tight_layout()

    for group in groups:
        if isinstance(cmaps[group_by],dict)!=True:
            if group==-1:
                group_colour=insitu_colour
            else:
                group_colour=cmaps[group_by][int(group)]
        else:
            group_cmap=cmaps[group_by]['cmap']
            group_norm=cmaps[group_by]['norm']

            group_cmap.set_under(insitu_colour)

            group_colour=group_cmap(group_norm(group))
        
        group_df=df.filter('%s==%s'%(columns[group_by],group)).extract()
        
        for i in range(3):
            plt.sca(ax[i])
            x_data=group_df.evaluate(x_axes[i])
            y_data=group_df.evaluate(y_axes[i])

            
            N,_,_=np.histogram2d(x_data,y_data,bins=[xbins[i],ybins[i]])
            X, Y = np.meshgrid(xbins[i]+0.5*(xbins[i][1]-xbins[i][0]), ybins[i]+0.5*(ybins[i][1]-ybins[i][0]))
            
            alpha=N*1.0/N_total[i-1]
            
            N[(N==0)]=np.nan
            N[(N>0)]=1.0
            
            alpha_bias=0.3
            alpha=1.0*alpha.T*(np.log10(N_total[i-1].T)+alpha_bias)/np.max(np.log10(N_total[i-1].T)+alpha_bias)
            
            alpha[(alpha>=0.0)]+=0.15
            alpha[(alpha>=1.0)]=1.0
            
            ax[i].pcolormesh(X,Y,N.T,cmap=colors.ListedColormap(group_colour),alpha=alpha,rasterized=True)

    for i in range(3):
        plt.sca(ax[i])
        plt.xlabel(xlabels[i])
        plt.ylabel(ylabels[i])

    plt.tight_layout(w_pad=1)

    plt.savefig(f'{save_path}{group_by}s.pdf')
    plt.savefig(f'{save_path}{group_by}s.png',dpi=250,bbox_inches='tight')
              
        

