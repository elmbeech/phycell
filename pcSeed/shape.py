#########
# title: shape.py
#
# language: python3
# date: 2022-09-18
# license: BSD-3-Clause
# author: Elmar Bucher
#
# description:
#     foo figther the color and the shape.
#########


# libraries
import itertools
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.patches as mpatches
import numpy as np
import os
import pandas as pd
import random
import sys


# functions
def df_label_to_color(df_abc, s_label, s_cmap='viridis', b_shuffle=False):
    '''
    input:
        df_abc: dataframe
        s_label: column name for which a color column should be generated
    output:
        df_abc: updated with color column
    description:
      for the selected label column add a color column to the dataframe
    '''
    ls_label = sorted(df_abc.loc[:,s_label].unique())
    if b_shuffle:
       random.shuffle(ls_label)
    a_color = plt.get_cmap(s_cmap)(np.linspace(0, 1, len(ls_label)))
    d_color = dict(zip(ls_label, a_color))
    df_abc[f'{s_label}_color'] = [colors.to_hex(d_color[s_scene]) for s_scene in df_abc.loc[:,s_label]]

def ax_colorlegend(ax, df_abc, s_label, s_color, r_x_figure2legend_space=0.01, s_fontsize='small'):
    '''
    description:
      add color legend to figure
    '''
    d_color = df_abc.loc[:,[s_label,s_color]].drop_duplicates().set_index(s_label).loc[:,s_color].to_dict()
    lo_patch = []
    for s_label, s_color in sorted(d_color.items()):
        o_patch = mpatches.Patch(color=s_color, label=s_label)
        lo_patch.append(o_patch)
    ax.legend(
        handles = lo_patch, 
        bbox_to_anchor = (1+r_x_figure2legend_space, 0, 0, 0), 
        loc = 'lower left', 
        borderaxespad = 0.00, 
        fontsize = s_fontsize
    )

def z_stack(df_coor, show=True, plot=None, png=None, movie=False, facecolor='white', frame_rate=24):
    """
    plot z stack from shape df
    png None or name
    movie None or name (png)
    """
    # handle input
    i_min_axis = min({df_coor['x'].min() - 1, df_coor['y'].min() - 1})
    i_max_axis = max({df_coor['x'].max() + 1, df_coor['y'].max() + 1})
    li_z = sorted(df_coor['z'].unique(), reverse=True)

    # show and plot mesh or agenst
    if show or not (plot is None):
        # generate plot
        i_z_total = len(li_z)
        fig, ax = plt.subplots(
            nrows=i_z_total, ncols=1, sharex=True, figsize=(4, 4*i_z_total)
        )
        ax = np.ravel(ax)
        for i , i_z in enumerate(li_z):
            df_z = df_coor.loc[df_coor['z']==i_z, :]
            df_z.plot(
                kind = 'scatter',
                x = 'x',
                y = 'y',
                color = df_z.loc[:, 'type_color'],
                xlim = (i_min_axis, i_max_axis),
                ylim = (i_min_axis, i_max_axis),
                grid = True,
                title = f'z_stack: {i_z}',
                ax = ax[i]
            )
        if show:
            # show plot
            fig.show()
        if not (plot is None):
            # handle path and filename
            ls_pathfile = png.split('/')
            s_file = ls_pathfile.pop(-1)
            s_path = '/'.join(ls_pathfile)
            if (len(s_path) > 0):
               s_path = s_path + '/'
            os.makedirs(s_path, exist_ok=True)
            # save plot
            fig.savefig(f'{s_path}{s_file}', facecolor=facecolor)

    # mesh or agnets to png
    if not (png is None):
        # handle path and filename
        ls_pathfile = png.split('/')
        s_shape = ls_pathfile.pop(-1).replace('.png','')
        s_path = '/'.join(ls_pathfile)
        if (len(s_path) > 0):
           s_path = s_path + '/'
        os.makedirs(s_path, exist_ok=True)
        # generate png
        for i , i_z in enumerate(li_z):
            s_title = f'{s_shape}_z_layer'
            fig, ax = plt.subplots(figsize=(4, 4))
            df_z = df_coor.loc[df_coor['z']==i_z, :]
            df_z.plot(
                kind = 'scatter',
                x = 'x',
                y = 'y',
                color = df_z.loc[:, 'type_color'],
                xlim = (i_min_axis, i_max_axis),
                ylim = (i_min_axis, i_max_axis),
                grid = True,
                title = f'{s_title}: {i_z}',
                ax = ax,
            )
            ax_colorlegend(
                ax = ax, 
                df_abc = df_z, 
                s_label = 'type', 
                s_color = 'type_color', 
                r_x_figure2legend_space = -0.99, 
                s_fontsize = 'small'
            )
            # save png
            plt.tight_layout()
            fig.savefig(f'{s_path}{s_title}_{str(i).zfill(6)}.png', facecolor=facecolor)
            plt.close()

    # mesh or agents to movie
    if movie and not (png is None):
        s_pathfile_movie = f'{s_path}{s_title}.mp4'
        os.system(f'ffmpeg -r {frame_rate} -f image2 -i {s_path}{s_title}_%06d.png -vcodec libx264 -pix_fmt yuv420p -strict -2 -tune animation -crf 15 -acodec none {s_pathfile_movie}')


##################
# class definion #
##################
class Shape:
    """
    generic class
    """
    def __init__(self):
        self.density = None
        self.agent = None # x,y,z,type dataframe
        self.coor = None  # one pixel is equal to one um!

    def volume(self):
        """
        output:
           volume: integer

        description:
            function returns volume in pixels.
            a pixel is thought to have 1 micrometer side length.
        """
        return(self.coor.shape[0])

    def union(self, shape):
        """
        # coor and agents
        """
        df_union_coor = pd.merge(
            self.coor,
            shape.coor,
            on=['m','n','p'],
            how='outer',
        )
        df_union_agent = pd.merge(
            self.agent,
            shape.agent,
            on=['m','n','p'],
            how='outer',
        )
        # output
        o_union = Shape()
        o_union.coor = df_union_coor
        o_union.agent = df_union_agent
        return(o_union)

    def intersection(self, shape):
        """
        # coor and agents
        """
        df_intersect = pd.merge(
            self.coor,
            shape.coor,
            on=['m','n','p'],
            how='inner',
        )
        # output
        o_intersection = Shape()
        o_intersection.coor = df_intersect
        return(o_intersection)

    def difference(self, shape):
        """
        # coor and agents
        """
        df_diff = self.coor.copy()
        df_diff['coor'] = df_diff.loc[:,'m'].astype(str) + '_' + df_diff.loc[:,'n'].astype(str) + '_' + df_diff.loc[:,'p'].astype(str)
        es_shape = set(shape.coor.loc[:,'m'].astype(str) + '_' + shape.coor.loc[:,'n'].astype(str) + '_' + shape.coor.loc[:,'p'].astype(str))
        df_diff['intersection'] = df_diff['coor'].isin(es_shape)
        df_diff = df_diff.loc[~ df_diff.loc[:,'intersection'],['m','n','p']]
        # output
        o_diff = Shape()
        o_diff.coor = df_diff
        return(o_diff)

    def symetric_difference(self, shape):
        """
        # coor and agents
`       """
        # left
        df_diff_left = self.difference(shape).coor
        # right
        df_diff_right = shape.coor.copy()
        df_diff_right['coor'] = df_diff_right.loc[:,'m'].astype(str) + '_' + df_diff_right.loc[:,'n'].astype(str) + '_' + df_diff_right.loc[:,'p'].astype(str)
        es_self = set(self.coor.loc[:,'m'].astype(str) + '_' + self.coor.loc[:,'n'].astype(str) + '_' + self.coor.loc[:,'p'].astype(str))
        df_diff_right['intersection'] = df_diff_right['coor'].isin(es_self)
        df_diff_right = df_diff_right.loc[~ df_diff_right.loc[:,'intersection'],['m','n','p']]
        # fusion
        df_diff_sym = pd.merge(
            df_diff_left,
            df_diff_right,
            on=['m','n','p'],
            how='outer'
        )
        # output
        o_diff_sym = Shape()
        o_diff_sym.coor = df_diff_sym
        return(o_diff_sym)

    def set_seeding_density_by_agent_count(self, total):
        """
        input:
            total: integer
                number of agents that about should be seeded.

        output:
            shape.density: float

        description:
            set agent seeding density (agent / volume)
            for that mesh over the number of agents.
        """
        self.density = total / self.volume()

    def set_seeding_density(self, density):
        """
        input:
            density: float
                agent density that should be seeded.

        output:
            shape.density: float

        description:
            set agent seeding density (agens / volume)
            for that mesh explicit.
        """

    def agent_seeding(self, agent_type_fraction):
        """
        input:
            agent_type_fraction: dictionary
                dictionay with agent_type label as key and fraction as value.
                all fractions iagent_type_fractionn the dictionary should sum up to 1.
                e.g. {'a': 1} or {'a': 2/3,'b': 1/3} or {'a': 0.75, 'b': 0.25}

        description:

        """
        # check if agent seeding density is set
        if (self.density is None):
            sys.exit(f'Error @ Shape agent_seeding: no agent seeding density set.\nplese provide density eiter by shape.set_seeding_density or shape.set_seeding_density_by_agent_count function!')

        # check if type fractions sum up to 1
        if (sum(agent_type_fraction.values()) != 1):
            sys.exit(f'Error @ Shape agent_seeding: the agent type fractions do not sum up to 1 {sorted(agent_type_fraction.items())} = {sum(agent_type_fraction.values())}.')

        # save already seeded agents
        df_agent = None
        if not (self.agent is None):
            df_agent = self.agent.copy()
            # check existig types
            es_exit = set(df_agent.loc[:,'type'])
            es_seed = set(agent_type_fraction.keys())
            if (len(es_exit.intesetcion(es_seed)) > 0):
                sys.exit(f'Error @ Shape agent_seeding: provided agent_type set {sorted(es_seed)} overlaps with already seeded agents {sorted(es_exit)}.\nplease adjust the agent_type_fraction dictionary.')

        # random seed agents
        r_total = self.density * self.volume()
        for s_type, r_fract in agent_type_fraction.items():
            # seed agents
            i_seed = int(np.ceil(r_fract * r_total)) # better more then less
            li_index = random.sample(list(self.coor.index), i_seed)
            df_seed = self.coor.loc[li_index,:].copy()
            df_seed.rename({'m':'x', 'n':'y', 'p':'z'}, axis=1, inplace=True)
            df_seed['type'] = s_type
            # jitter position
            df_seed.loc[:,'x'] = df_seed.loc[:,'x'] + np.random.uniform(low=-0.5, high=0.5, size=i_seed)
            df_seed.loc[:,'y'] = df_seed.loc[:,'y'] + np.random.uniform(low=-0.5, high=0.5, size=i_seed)
            df_seed.loc[:,'z'] = df_seed.loc[:,'z'] + np.random.uniform(low=-0.5, high=0.5, size=i_seed)
            # add to results
            if (df_agent is None):
                df_agent = df_seed
            else:
                df_agent = pd.concat([df_agent, df_seed])

        # reset agent_density to None
        self.agent_density = None

        # stor result
        self.agent = df_agent

    def df_mesh(self):
        """
        """
        return(self.coor)

    def df_agent(self):
        """
        """
        return(self.agent)

    def z_stack_mesh(self, show=True, png=None, movie=False, frame_rate=24):
        """
        plot z stack from shape df
        png None or name
        movie None or name (png)
        """
        df_coor = self.df_mesh().rename({'m':'x', 'n':'y', 'p':'z'}, axis=1)
        df_coor['type'] = 'mesh'
        df_label_to_color(df_abc=df_coor, s_label='type', s_cmap='turbo', b_shuffle=False)
        z_stack(
            df_coor = df_coor,
            show = show,
            png = png,
            movie  = movie,
            frame_rate = frame_rate,
        )

    def z_stack_agent(self, show=True, png=None, movie=False, frame_rate=24):
        """
        plot z stack from agent df
        png None or name
        movi None or name (png)
        """
        # plot agents
        df_coor = self.df_agent().copy()
        df_coor.z = df_coor.z.round().astype(int)
        df_label_to_color(df_abc=df_coor, s_label='type', s_cmap='turbo', b_shuffle=False)
        z_stack(
            df_coor = df_coor,
            show = show,
            png = png,
            movie  = movie,
            frame_rate = frame_rate,
        )


class Brick(Shape):
    """
    """
    def __init__(self, x, y, z=1, origin=(0,0,0)):
        """
        """
        # inhert
        super(Brick, self).__init__()

        # x keeping the zero center in mind
        i_x = int(round(x))
        if (i_x % 2 != 0):
            i_x += 1
        i_rx = int(i_x / 2)

        # y keeping the zero center in mind
        i_y = int(round(y))
        if (i_y % 2 != 0):
            i_y += 1
        i_ry = int(i_y / 2)

        # z keeping the zero center in mind
        i_z = int(round(z))
        if (i_z % 2 != 0):
            i_z += 1
        i_rz = int(i_z / 2)

        # handle origin input
        ti_origin = (int(round(origin[0])), int(round(origin[1])), int(round(origin[2])))

        # compute all coordinate combinations
        iti_m = range(ti_origin[0] - i_rx + 1, ti_origin[0] + i_rx)
        iti_n = range(ti_origin[1] - i_ry + 1, ti_origin[1] + i_ry)
        iti_p = range(ti_origin[2] - i_rz + 1, ti_origin[2] + i_rz)
        df_coor = pd.DataFrame(itertools.product(iti_m, iti_n, iti_p), columns=['m','n','p'])

        # store result
        self.coor = df_coor


class Sphere(Shape):
    """
    """
    def __init__(self, d, z=1, origin=(0,0,0)):
        """
        """
        # inhert
        super(Sphere, self).__init__()

        # compute brick coordinates
        if (z is None):
            i_z = d
        else:
            i_z = z
        o_brick = Brick(x=d, y=d, z=i_z, origin=origin)
        df_coor = o_brick.coor

        # comput r keeping the zero center in mind
        i_d = int(round(d))
        if (i_d % 2 != 0):
            i_d += 1
        i_r = int(i_d / 2)

        # handle origin input
        ti_origin = (int(round(origin[0])), int(round(origin[1])), int(round(origin[2])))

        # check with simple pythagoras, if coordinat is inside the sphere
        df_coor['r2'] = i_r**2
        df_coor['m2'] = (df_coor.loc[:,'m'] - ti_origin[0])**2
        df_coor['n2'] = (df_coor.loc[:,'n'] - ti_origin[1])**2
        df_coor['p2'] = (df_coor.loc[:,'p'] - ti_origin[2])**2
        df_coor['m2n2p2'] = df_coor.loc[:,['m2','n2','p2']].sum(axis=1)
        df_coor['inside'] = df_coor.loc[:,'m2n2p2'] <= df_coor.loc[:,'r2']
        df_coor = df_coor.loc[df_coor.inside, ['m','n','p']]

        # store result
        self.coor = df_coor

