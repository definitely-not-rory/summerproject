from imports import *
import aux_functions as calc
import plot_generation as plot

def progenitor_recovery(halo,lsr_def='8kpc',vtoomre=False,home_dir='/cosma/apps/durham/dc-coll7/auriga/',N_cutoff=20):

    lsr_defs=['8kpc','scalelength'] #Labels of permitted definitions of Local Solar Neighbourhood ('8kpc' -> proper radial distance from galactic centre, 'scalelength' -> distance that scales based on properties of individual halo).

    if lsr_def not in lsr_defs:
        print('\nInvalid LSR definition selected, please select \'8kpc\' or \'scalelength\'') #Ensure selected Local Solar Neighbourhood is valid.
        return

    if vtoomre==False: #Allocates relevant label to particle subselection method ('accreted' -> particles tagged as accreted by AURIGA, 'vtoomre' -> selection over all particles using Toomre velocity threshold).
        selection_type='accreted'
    else:
        selection_type='vtoomre'

    if os.path.exists(f'{home_dir}{halo}/{lsr_def}/{selection_type}/results')!=True: #Checks if clustering algorithm has been run on selected dataset.
        print(f'No clustering data detected for {halo} in {home_dir}, generating with default parameters.')
        return
    else:
        print(f'{halo} clustering data located.')

    results_dir=f'{home_dir}{halo}/{lsr_def}/{selection_type}/results/progenitors' #Initialises quick access strings for clustering run label and file paths.
    run_name=f'{halo}_{lsr_def}_{selection_type}'

    os.makedirs(results_dir,exist_ok=True)

    df=vaex.open(f'{home_dir}{halo}/{lsr_def}/{selection_type}/results/{run_name}_ChemistryGroups.hdf5')

    progenitors=df.filter('progenitor_id!=-1').extract().evaluate('progenitor_id')

    unique_progenitors,N_progenitor=np.unique(progenitors,return_counts=True)
    threshold_progenitors=unique_progenitors[np.where(N_progenitor>=20)]
 
    dominant_predicted_clusters=[]
    dominant_predicted_KS_groups=[]

    cluster_completenesses=[]
    KS_group_completenesses=[]

    for progenitor in threshold_progenitors:
        progenitor_df=df.filter('progenitor_id==%s'%progenitor).extract()
        N_total_prog_stars=progenitor_df.count()

        prog_clusters, prog_stars_per_cluster=np.unique(progenitor_df.filter('label!=-1').extract().evaluate('label'),return_counts=True)
        prog_KS_groups, prog_stars_per_KS_group=np.unique(progenitor_df.filter('progenitor_id==%s'%progenitor).extract().filter('KS_groups!=-1').extract().evaluate('KS_groups'),return_counts=True)

        dominant_predicted_cluster=int(prog_clusters[np.argmax(prog_stars_per_cluster)] if len(prog_clusters)!=0 else -1)
        N_dominant_predicted_cluster=int(prog_stars_per_cluster[np.argmax(prog_stars_per_cluster)] if len(prog_clusters)!=0 else 0)

        dominant_predicted_KS_group=int(prog_KS_groups[np.argmax(prog_stars_per_KS_group)] if len(prog_KS_groups)!=0 else -1)
        N_dominant_predicted_KS_group=int(prog_stars_per_KS_group[np.argmax(prog_stars_per_KS_group)] if len(prog_KS_groups)!=0 else 0)

        dominant_predicted_clusters.append(dominant_predicted_cluster)
        dominant_predicted_KS_groups.append(dominant_predicted_KS_group)

        cluster_completenesses.append(N_dominant_predicted_cluster/N_total_prog_stars)
        KS_group_completenesses.append(N_dominant_predicted_KS_group/N_total_prog_stars)

    clusters=np.unique(df.filter('label!=-1').extract().evaluate('label'))
    cluster_purities=[]

    for cluster in clusters:
        cluster_df=df.filter('label==%s'%cluster).extract()
        N_cluster=cluster_df.count()

        cluster_mcp=np.unique(cluster_df.filter('cluster_mcp!=-1').extract().evaluate('cluster_mcp'))[0]
        N_cluster_mcp=cluster_df.filter('progenitor_id==%s'%cluster_mcp).count()

        cluster_purities.append(N_cluster_mcp/N_cluster)
    
    KS_groups=np.unique(df.filter('KS_groups!=-1').extract().evaluate('KS_groups'))
    KS_group_purities=[]

    for KS_group in KS_groups:
        KS_group_df=df.filter('KS_groups==%s'%KS_group).extract()
        N_KS_group=KS_group_df.count()

        KS_group_mcp=np.unique(KS_group_df.filter('KS_mcp!=-1').extract().evaluate('KS_mcp'))[0]
        N_KS_group_mcp=KS_group_df.filter('progenitor_id==%s'%KS_group_mcp).count()

        KS_group_purities.append(N_KS_group_mcp/N_KS_group)

    cluster_realness=len(np.unique(np.array(cluster_purities)[np.array(cluster_purities)>=(2/3)]))/len(clusters)
    KS_realness=len(np.unique(np.array(KS_group_purities)[np.array(KS_group_purities)>=(2/3)]))/len(KS_groups)

    cluster_recovered_progs=[]
    KS_recovered_progs=[]

    for index,progenitor in enumerate(threshold_progenitors):
        dominant_predicted_cluster=dominant_predicted_clusters[index]
        dominant_predicted_KS_group=dominant_predicted_KS_groups[index]

        cluster_completeness=cluster_completenesses[index]
        KS_completeness=KS_group_completenesses[index]

        cluster_purity=np.array(cluster_purities)[np.where(clusters==dominant_predicted_cluster)][0] if len(np.array(cluster_purities)[np.where(clusters==dominant_predicted_cluster)])!=0 else 0
        KS_purity=np.array(KS_group_purities)[np.where(KS_groups==dominant_predicted_KS_group)][0] if len(np.array(KS_group_purities)[np.where(KS_groups==dominant_predicted_KS_group)])!=0 else 0

        if cluster_purity>=2/3 and cluster_completeness>=0.5:
            cluster_recovered_progs.append(progenitor)
        
        if KS_purity>=2/3 and KS_completeness>=0.5:
            KS_recovered_progs.append(progenitor)

    cluster_recovery_rate=len(cluster_recovered_progs)/len(threshold_progenitors)
    KS_recovery_rate=len(KS_recovered_progs)/len(threshold_progenitors)

    data={'cluster_completeness':cluster_completenesses,'KS_completeness':KS_group_completenesses,'cluster_purity':cluster_purities,'KS_purity':KS_group_purities,'cluster_realness':cluster_realness,'KS_realness':KS_realness,'cluster_recovery':cluster_recovery_rate,'KS_recovery':KS_recovery_rate}
    return data

        
        

  

    


    

        


                

        
        
        

    

    
