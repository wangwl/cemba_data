import seaborn as sns
from .cell import categorical_scatter


def _plot_kmode_overall_single_r(ax, overall_df):
    sns.lineplot(data=overall_df, color='orange',
                 x='k', y='pureness', ax=ax, label='Pureness',
                 markers=True, dashes=True, legend=None)
    sns.scatterplot(data=overall_df, color='orange',
                 x='k', y='pureness', ax=ax, legend=None)
    ax_twin = ax.twinx()
    sns.lineplot(data=overall_df, color='steelblue',
                 x='k', y='completeness', ax=ax_twin, label='Completeness',
                 markers=True, dashes=True, legend=None, )
    sns.scatterplot(data=overall_df, color='steelblue',
                 x='k', y='completeness', ax=ax_twin, legend=None)
    ax.grid()
    ax.set(ylim=(0.98, 1))
    ax_twin.set(ylim=(0.8, 1))
    return ax


def plot_kmode_overall(axes, overall_df):
    for ax, (resolution, sub_df) in zip(axes.flat, overall_df.groupby('resolution')):
        _plot_kmode_overall_single_r(ax, sub_df)
        ax.set_title(f'Resolution = {resolution}')
    return axes


def plot_kmode_stats(axes, coord_data, result_dict, coord_base='umap'):
    cluster_pureness = result_dict['cluster_pureness']
    cluster_completeness = result_dict['cluster_completeness']
    cell_ambiguity = result_dict['cell_ambiguity']
    cluster = result_dict['cluster']

    plot_data = coord_data.copy()
    plot_data['cell_ambiguity'] = cell_ambiguity
    plot_data['cluster'] = cluster
    plot_data['cluster_pureness'] = cluster.map(cluster_pureness)
    plot_data['cluster_completeness'] = cluster.map(cluster_completeness)

    if axes.size != 4:
        raise ValueError('Number of axes is not 4.')
    ax1, ax2, ax3, ax4 = axes.flatten()
    categorical_scatter(data=plot_data, hue='cluster', ax=ax1,
                        coord_base=coord_base, palette='tab20')
    if ax1.is_first_row():
        ax1.set_title('Cluster')
    sns.scatterplot(data=plot_data, x=f'{coord_base}_0', y=f'{coord_base}_1',
                    hue=1 - plot_data['cell_ambiguity'], ax=ax2, hue_norm=(0.95, 1),
                    palette='viridis', s=3, linewidth=0, legend=None)
    if ax2.is_first_row():
        ax2.set_title('1 - Cell Ambiguity')
    sns.scatterplot(data=plot_data, x=f'{coord_base}_0', y=f'{coord_base}_1',
                    hue='cluster_pureness', ax=ax3, hue_norm=(0.7, 1),
                    palette='viridis', s=3, linewidth=0, legend=None)
    if ax3.is_first_row():
        ax3.set_title('Cluster Pureness')
    sns.scatterplot(data=plot_data, x=f'{coord_base}_0', y=f'{coord_base}_1',
                    hue='cluster_completeness', ax=ax4, hue_norm=(0.7, 1),
                    palette='viridis', s=3, linewidth=0, legend=None)
    if ax4.is_first_row():
        ax4.set_title('Cluster Completeness')
    for ax in axes.flat:
        ax.axis('off')
    return axes