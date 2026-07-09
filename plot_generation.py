from imports import *

def show_progenitors(df):
    total_stars=df.count()
    accreted_stars=df.count(selection='prog!=-2')
    insitu_stars=df.count(selection='prog==-2')

    print(f'Total Stars: {total_stars}\nAccreted: {accreted_stars}\nIn-Situ: {insitu_stars}')

    progenitors=df.filter('prog!=-2').extract().evaluate('prog')

    unique_progenitors=np.sort(np.unique(progenitors))
    progenitor_labels=[str(progenitor) for progenitor in unique_progenitors]

    stars_per_progenitor=[np.count_nonzero(progenitors==progenitor) for progenitor in unique_progenitors]

    sorted_prog_counts=np.array(sorted([[progenitor,int(stars_per_progenitor[progenitor_labels.index(progenitor)])] for progenitor in progenitor_labels],key=lambda x:x[1],reverse=True))

    plt.figure(figsize=(18,8))
    plt.bar('In-Situ',insitu_stars,color='teal')
    plt.bar(sorted_prog_counts[:,0],[int(count) for count in sorted_prog_counts[:,1]],color='dodgerblue')
    plt.yscale('log')
    plt.xlabel('Progenitor')
    plt.ylabel('Count')
    plt.xticks(rotation=-45)
    plt.show()

def iom_scatters(df,x,y):
    dimensions={'En':[-150000,-60000],'Lz':[-4000,1000],'Fe_H':[-3,0.5]}

    axes=[x,y]
    cbar_dim=[dimension for dimension in dimensions if dimension not in axes][0]

    df.filter('Fe_H>-3').extract().scatter(x,y,c='silver',s=1)
    df.filter('(Fe_H>-3)&(prog!=-2)').extract().scatter(x,y,c=df.filter('(Fe_H>-3)&(prog!=-2)').extract().evaluate(cbar_dim),cmap='magma',s=1)

    plt.colorbar(label=cbar_dim)
    plt.clim(dimensions[cbar_dim][0],dimensions[cbar_dim][1])
    plt.show()

def stellar_FeH_hists(df, **kwargs):
    if 'normalised' in kwargs:
        normalised=kwargs['normalised']
    else:
        normalised=False
    plt.hist(df.filter('(Fe_H>-3)').extract().evaluate('Fe_H'),label='all stars',histtype='step',bins=np.arange(-3,1,0.2),linestyle='dashed',color='midnightblue',density=normalised)
    plt.hist(df.filter('(Fe_H>-3)&(prog==-2)').extract().evaluate('Fe_H'),label='in-situ',histtype='step',bins=np.arange(-3,1,0.2),color='teal',density=normalised)
    plt.hist(df.filter('(Fe_H>-3)&(prog!=-2)').extract().evaluate('Fe_H'),label='accreted',histtype='step',bins=np.arange(-3,1,0.2),color='dodgerblue',density=normalised)


    plt.xlabel('FeH')
    if normalised==True:
        plt.ylabel('Normalised N')
    else:
        plt.ylabel('N')
        plt.yscale('log')
    plt.legend()
    plt.title('Halo stars')
    plt.show()