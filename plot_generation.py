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
    plt.bar('In-Situ',insitu_stars,color='silver')
    plt.bar(sorted_prog_counts[:,0],[int(count) for count in sorted_prog_counts[:,1]])
    plt.yscale('log')
    plt.xlabel('Progenitor')
    plt.ylabel('Count')
    plt.xticks(rotation=-45)
    plt.show()