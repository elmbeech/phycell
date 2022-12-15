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
from pcSeed import pdplt  # from pwn biotransistor prj
import random
import sys


# functions
def axis(r_min, r_max, r_step, r_origin=0):
    """
    """
    # negative from origin
    er_neg = set()
    if (r_min <= r_origin):
        er_neg = set(np.arange(r_origin, (r_min - r_step), -r_step))
        er_neg = {n for n in er_neg if (n <= r_max) and (n >= r_min)}
    # positive from origin
    er_pos = set()
    if (r_max >= r_origin):
        er_pos = set(np.arange(r_origin, (r_max + r_step), +r_step))
        er_pos = {n for n in er_pos if (n >= r_min) and (n <= r_max)}
    # output
    er_axis = er_neg.union(er_pos)
    return(er_axis)

def z_stack(df_coor, show=False, plot=None, png=None, movie=False, facecolor='white', s=1, figsize=(4, 4), frame_rate=24):
    """
    plot z stack from shape df
    png None or name
    movie None or name (png)
    """
    # handle input
    i_min_axis = min({df_coor['x'].min() - 1, df_coor['y'].min() - 1})
    i_max_axis = max({df_coor['x'].max() + 1, df_coor['y'].max() + 1})
    li_z = sorted(df_coor['z'].unique(), reverse=True)

    # show and plot mesh or agents
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
        if show:
            # show plot
            fig.show()

    # mesh or agents to png
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
            fig, ax = plt.subplots(figsize=figsize)
            df_z = df_coor.loc[df_coor['z']==i_z, :]
            df_z.plot(
                kind = 'scatter',
                x = 'x',
                y = 'y',
                s = s,
                color = df_z.loc[:, 'type_color'],
                xlim = (i_min_axis, i_max_axis),
                ylim = (i_min_axis, i_max_axis),
                grid = True,
                title = f'{s_title}: {i_z}',
                ax = ax,
            )
            pdplt.ax_colorlegend(
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
        self.agent = None # x,y,z,type dataframe
        self.coor = None  # one pixel is equal to one um!
        self.um_p_px = None # 1 is default

    def volume(self):
        """
        output:
           volume: integer in um**2 or um**3

        description:
            function returns volume in pixels.
            a pixel is thought to have 1 micrometer side length.
        """
        #return(self.coor.shape[0])
        i_volumne = self.coor.shape[0]
        i_p = len(set(self.coor.p))
        if (i_p == 1):
            i_volumne = i_volumne * self.um_p_px * self.um_p_px
        else:
            i_volumne = i_volumne * self.um_p_px * self.um_p_px * self.um_p_px
        return(i_volumne)

    def union(self, shape):
        """
        # coor and agents
        """
        if (self.um_p_px != shape.um_p_px):
            sys.exit(f'Error @ Shape.union : the self and shape object are not the same um per pixel scale ({self.um_p_px} != {shape.um_p_px}).\nobjects can not be fused!')
        df_union_coor = pd.merge(
            self.coor,
            shape.coor,
            on=['m','n','p'],
            how='outer',
        )
        df_union_agent = pd.merge(
            self.agent,
            shape.agent,
            on=['x','y','z','type'],
            how='outer',
        )
        df_union_agent = df_union_agent.astype({'x': np.float32, 'y': np.float32, 'z': np.float32, 'type': str})
        # output
        o_union = Shape()
        o_union.coor = df_union_coor
        o_union.agent = df_union_agent
        o_union.um_p_px = self.um_p_px
        return(o_union)

    def intersection(self, shape):
        """
        # coor and agents
        """
        if (self.um_p_px != shape.um_p_px):
            sys.exit(f'Error @ Shape.intersection : the self and shape object are not the same um per pixel scale ({self.um_p_px} != {shape.um_p_px}).\nobjects can not be fused!')
        df_intersect = pd.merge(
            self.coor,
            shape.coor,
            on=['m','n','p'],
            how='inner',
        )
        # output
        o_intersection = Shape()
        o_intersection.coor = df_intersect
        o_intersection.um_p_px = self.um_p_px
        return(o_intersection)

    def difference(self, shape):
        """
        # coor and agents
        """
        if (self.um_p_px != shape.um_p_px):
            sys.exit(f'Error @ Shape.difference : the self and shape object are not the same um per pixel scale ({self.um_p_px} != {shape.um_p_px}).\nobjects can not be fused!')
        df_diff = self.coor.copy()
        df_diff['coor'] = df_diff.loc[:,'m'].astype(str) + '_' + df_diff.loc[:,'n'].astype(str) + '_' + df_diff.loc[:,'p'].astype(str)
        es_shape = set(shape.coor.loc[:,'m'].astype(str) + '_' + shape.coor.loc[:,'n'].astype(str) + '_' + shape.coor.loc[:,'p'].astype(str))
        df_diff['intersection'] = df_diff['coor'].isin(es_shape)
        df_diff = df_diff.loc[~ df_diff.loc[:,'intersection'],['m','n','p']]
        # output
        o_diff = Shape()
        o_diff.coor = df_diff
        o_diff.um_p_px = self.um_p_px
        return(o_diff)

    def symetric_difference(self, shape):
        """
        # coor and agents
        """
        if (self.um_p_px != shape.um_p_px):
            sys.exit(f'Error @ Shape.symetric_difference : the self and shape object are not the same um per pixel scale ({self.um_p_px} != {shape.um_p_px}).\nobjects can not be fused!')
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
        o_diff_sym.um_p_px = self.um_p_px
        return(o_diff_sym)

    def seeding_random(self, agent_type_fraction, agent_count=None, agent_density=None):
        """
        input:
            agent_type_fraction: dictionary
                dictionay with agent_type label as key and fraction as value.
                all fractions iagent_type_fractionn the dictionary should sum up to 1.
                e.g. {'a': 1} or {'a': 2/3,'b': 1/3} or {'a': 0.75, 'b': 0.25}
            agent_density in agent/um**2 or agent/um**3

        description:

        """
        # check if type fractions sum up to 1
        if (sum(agent_type_fraction.values()) != 1):
            sys.exit(f'Error @ Shape seeding_random: the agent type fractions do not sum up to 1 {sorted(agent_type_fraction.items())} = {sum(agent_type_fraction.values())}.')

        # check if agent seeding density is set
        if ((agent_count is None) and (agent_density is None)) or (not (agent_count is None) and not (agent_density is None)):
            sys.exit(f'Error @ Shape seeding_random: either agent_density ({agent_density}) or agent_count ({agent_count}) have to be set.\nNone or both are set!')
        if (agent_density is None):
            agent_density = agent_count / self.volume()

        # save already seeded agents
        df_agent = None
        if not (self.agent is None):
            df_agent = self.agent.copy()
            # check existig types
            es_exist = set(df_agent.loc[:,'type'])
            es_seed = set(agent_type_fraction.keys())
            if (len(es_exist.intesetcion(es_seed)) > 0):
                sys.exit(f'Error @ Shape agent_seeding: provided agent_type set {sorted(es_seed)} overlaps with already seeded agents {sorted(es_exist)}.\nplease adjust the agent_type_fraction dictionary.')

        # random seed agents
        r_total = agent_density * self.volume()
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

        # stor result
        self.agent = df_agent

    def seeding_hexagonal(self, agent_type_fraction, agent_diameter_um, lattice='HPC'):
        """
        agent_diameter_um in um
        lattice='HPC' 'FCC'
        """
        # handle lattice
        if (lattice.upper() == 'HPC'):
            i_layer = 2
        elif (lattice.upper() == 'FCC'):
            i_layer = 3
        else:
            sys.exit(f'Error @ Shape seeding_hexagonal: unknowen hexagonal lattice packing type {lattice}.\nknown are HPC and FCC.')

        # check if type fractions sum up to 1
        if (sum(agent_type_fraction.values()) != 1):
            sys.exit(f'Error @ Shape seeding_hexagonal: the agent type fractions do not sum up to 1 {sorted(agent_type_fraction.items())} = {sum(agent_type_fraction.values())}.')

        # save already seeded agents
        df_agent = None
        if not (self.agent is None):
            df_agent = self.agent.copy()
            # check existig types
            es_exist = set(df_agent.loc[:,'type'])
            es_seed = set(agent_type_fraction.keys())
            if (len(es_exist.intersection(es_seed)) > 0):
                sys.exit(f'Error @ Shape agent_seeding: provided agent_type set {sorted(es_seed)} overlaps with already seeded agents {sorted(es_exist)}.\nplease adjust the agent_type_fraction dictionary.')

        # bue 2022-12-01: code to keep you sane while developing!
        #import itertools
        #import numpy as np
        #import pandas as pd
        #agent_type_fraction = {'abc': 1}
        #agent_diameter_um = 10
        #lattice='HPC'
        #r_d_major = (agent_diameter_um * 2)

        # calculate hexagon radius and such
        r_d_major = (agent_diameter_um * 2) / self.um_p_px  # diameter in pixel
        r_3r_major = r_d_major * 1.5
        r_r_major = r_d_major / 2
        r_r_major_half = r_d_major / 4
        r_r_minor = np.sqrt(r_r_major**2 - r_r_major_half**2)
        r_d_minor = r_r_minor * 2
        r_r_minor_half = r_r_minor / 2
        r_r_minor_third = r_r_minor / 3

        # get min and max
        ai_m = np.array(sorted(set(self.coor.loc[:,'m'])))
        ai_n = np.array(sorted(set(self.coor.loc[:,'n'])))
        ai_p = np.array(sorted(set(self.coor.loc[:,'p'])))
        i_m_min = ai_m.min()
        i_m_max = ai_m.max()
        i_n_min = ai_n.min()
        i_n_max = ai_n.max()
        i_p_min = ai_p.min()
        i_p_max = ai_p.max()

        # bue 2022-12-01: code to keep you sane while developing!
        #i_m_min = -100
        #i_m_max = 100
        #i_n_min= -100
        #i_n_max = 100
        #i_p_min = 0
        #i_p_max = 1

        # layer 1a axis n and m
        lr_axis_l1a_m = sorted(axis(r_min=i_m_min, r_max=i_m_max, r_step=r_d_minor, r_origin=0))
        lr_axis_l1a_n = sorted(axis(r_min=i_n_min, r_max=i_n_max, r_step=r_3r_major, r_origin=0))
        # layer 1b axis n and m
        lr_axis_l1b_m = sorted(axis(r_min=i_m_min, r_max=i_m_max, r_step=r_d_minor, r_origin=0))
        lr_axis_l1b_n = sorted(axis(r_min=i_n_min, r_max=i_n_max, r_step=r_3r_major, r_origin=r_d_major))
        # layer 1c axis n and m
        lr_axis_l1c_m = sorted(axis(r_min=i_m_min, r_max=i_m_max, r_step=r_d_minor, r_origin=r_r_minor))
        lr_axis_l1c_n = sorted(axis(r_min=i_n_min, r_max=i_n_max, r_step=r_3r_major, r_origin=r_r_major_half))
        # layer 1d axis n and m
        lr_axis_l1d_m = sorted(axis(r_min=i_m_min, r_max=i_m_max, r_step=r_d_minor, r_origin=r_r_minor))
        lr_axis_l1d_n = sorted(axis(r_min=i_n_min, r_max=i_n_max, r_step=r_3r_major, r_origin=(r_r_major_half-r_d_major)))
        # layer 1 axis p
        lr_axis_l1_p = sorted(axis(r_min=i_p_min, r_max=i_p_max, r_step=r_d_minor * i_layer, r_origin=0))
        # get layer 1 coordinates
        df_hexcoor_layer1a = pd.DataFrame(itertools.product(lr_axis_l1a_m, lr_axis_l1a_n, lr_axis_l1_p), columns=['x','y','z'])
        df_hexcoor_layer1b = pd.DataFrame(itertools.product(lr_axis_l1b_m, lr_axis_l1b_n, lr_axis_l1_p), columns=['x','y','z'])
        df_hexcoor_layer1c = pd.DataFrame(itertools.product(lr_axis_l1c_m, lr_axis_l1c_n, lr_axis_l1_p), columns=['x','y','z'])
        df_hexcoor_layer1d = pd.DataFrame(itertools.product(lr_axis_l1d_m, lr_axis_l1d_n, lr_axis_l1_p), columns=['x','y','z'])
        df_hexcoor_layer1 = pd.concat([df_hexcoor_layer1a, df_hexcoor_layer1b, df_hexcoor_layer1c, df_hexcoor_layer1d])

        # bue 2022-12-01: code to keep you sane while developing!
        #import matplotlib.pyplot as plt
        #%matplotlib
        #fig, ax = plt.subplots()
        #df_hexcoor_layer1a.plot(kind='scatter', x='x', y='y', grid=True, ax=ax, c='maroon')
        #df_hexcoor_layer1b.plot(kind='scatter', x='x', y='y', grid=True, ax=ax, c='red')
        #df_hexcoor_layer1c.plot(kind='scatter', x='x', y='y', grid=True, ax=ax, c='orange')
        #df_hexcoor_layer1d.plot(kind='scatter', x='x', y='y', grid=True, ax=ax, c='yellow')

        if (ai_p.shape[0] == 1):
            # get layer 2 and 3 coordinates
            df_hexcoor_layer2 = pd.DataFrame()
            df_hexcoor_layer3 = pd.DataFrame()
        else:
            # layer 2a axis n and m
            lr_axis_l2a_m = sorted(axis(r_min=i_m_min, r_max=i_m_max, r_step=r_d_minor, r_origin=(0 + r_r_minor_third)))
            lr_axis_l2a_n = sorted(axis(r_min=i_n_min, r_max=i_n_max, r_step=r_3r_major, r_origin=(0 + r_r_major_half)))
            # layer 2b axis n and m
            lr_axis_l2b_m = sorted(axis(r_min=i_m_min, r_max=i_m_max, r_step=r_d_minor, r_origin=(0 + r_r_minor_third)))
            lr_axis_l2b_n = sorted(axis(r_min=i_n_min, r_max=i_n_max, r_step=r_3r_major, r_origin=(r_d_major + r_r_major_half)))
            # layer 2c axis n and m
            lr_axis_l2c_m = sorted(axis(r_min=i_m_min, r_max=i_m_max, r_step=r_d_minor, r_origin=(r_r_minor + r_r_minor_third)))
            lr_axis_l2c_n = sorted(axis(r_min=i_n_min, r_max=i_n_max, r_step=r_3r_major, r_origin=(r_r_major_half + r_r_major_half)))
            # layer 2d axis n and m
            lr_axis_l2d_m = sorted(axis(r_min=i_m_min, r_max=i_m_max, r_step=r_d_minor, r_origin=(r_r_minor + r_r_minor_third)))
            lr_axis_l2d_n = sorted(axis(r_min=i_n_min, r_max=i_n_max, r_step=r_3r_major, r_origin=(r_r_major_half-r_d_major + r_r_major_half)))
            # layer 2 axis p
            lr_axis_l2_p = sorted(axis(r_min=i_p_min, r_max=i_p_max, r_step=i_layer*r_d_minor, r_origin=r_d_minor))
            # get layer 2 coordinates
            df_hexcoor_layer2a = pd.DataFrame(itertools.product(lr_axis_l2a_m, lr_axis_l2a_n, lr_axis_l2_p), columns=['x','y','z'])
            df_hexcoor_layer2b = pd.DataFrame(itertools.product(lr_axis_l2b_m, lr_axis_l2b_n, lr_axis_l2_p), columns=['x','y','z'])
            df_hexcoor_layer2c = pd.DataFrame(itertools.product(lr_axis_l2c_m, lr_axis_l2c_n, lr_axis_l2_p), columns=['x','y','z'])
            df_hexcoor_layer2d = pd.DataFrame(itertools.product(lr_axis_l2d_m, lr_axis_l2d_n, lr_axis_l2_p), columns=['x','y','z'])
            df_hexcoor_layer2 = pd.concat([df_hexcoor_layer2a, df_hexcoor_layer2b, df_hexcoor_layer2c, df_hexcoor_layer2d])

            if (i_layer < 3):
                # get layer 3 coordinates
                df_hexcoor_layer3 = pd.DataFrame()
            else:
                # layer 3a axis n and m
                lr_axis_l3a_m = sorted(axis(r_min=i_m_min, r_max=i_m_max, r_step=r_d_minor, r_origin=(0 + 2*r_r_minor_third)))
                lr_axis_l3a_n = sorted(axis(r_min=i_n_min, r_max=i_n_max, r_step=r_3r_major, r_origin=(0 + 2*r_r_major_half)))
                # layer 3b axis n and m
                lr_axis_l3b_m = sorted(axis(r_min=i_m_min, r_max=i_m_max, r_step=r_d_minor, r_origin=(0 + 2*r_r_minor_third)))
                lr_axis_l3b_n = sorted(axis(r_min=i_n_min, r_max=i_n_max, r_step=r_3r_major, r_origin=(r_d_major + 2*r_r_major_half)))
                # layer 3c axis n and m
                lr_axis_l3c_m = sorted(axis(r_min=i_m_min, r_max=i_m_max, r_step=r_d_minor, r_origin=(r_r_minor + 2*r_r_minor_third)))
                lr_axis_l3c_n = sorted(axis(r_min=i_n_min, r_max=i_n_max, r_step=r_3r_major, r_origin=(r_r_major_half + 2*r_r_major_half)))
                # layer 3d axis n and m
                lr_axis_l3d_m = sorted(axis(r_min=i_m_min, r_max=i_m_max, r_step=r_d_minor, r_origin=(r_r_minor + 2*r_r_minor_third)))
                lr_axis_l3d_n = sorted(axis(r_min=i_n_min, r_max=i_n_max, r_step=r_3r_major, r_origin=(r_r_major_half-r_d_major + 2*r_r_major_half)))
                # layer 3 axis p
                lr_axis_l3_p = sorted(axis(r_min=i_p_min, r_max=i_p_max, r_step=i_layer*r_d_minor, r_origin=2*r_d_minor))
                # get layer 3 coordinates
                df_hexcoor_layer3a = pd.DataFrame(itertools.product(lr_axis_l3a_m, lr_axis_l3a_n, lr_axis_l3_p), columns=['x','y','z'])
                df_hexcoor_layer3b = pd.DataFrame(itertools.product(lr_axis_l3b_m, lr_axis_l3b_n, lr_axis_l3_p), columns=['x','y','z'])
                df_hexcoor_layer3c = pd.DataFrame(itertools.product(lr_axis_l3c_m, lr_axis_l3c_n, lr_axis_l3_p), columns=['x','y','z'])
                df_hexcoor_layer3d = pd.DataFrame(itertools.product(lr_axis_l3d_m, lr_axis_l3d_n, lr_axis_l3_p), columns=['x','y','z'])
                df_hexcoor_layer3 = pd.concat([df_hexcoor_layer3a, df_hexcoor_layer3b, df_hexcoor_layer3c, df_hexcoor_layer3d])

        # agent seeding coordinates
        df_seed = pd.concat([df_hexcoor_layer1, df_hexcoor_layer2, df_hexcoor_layer3])

        # filter by coordinate magic
        di_x = {}
        for r_x in sorted(set(df_seed.loc[:,'x'])):
            i_m = ai_m[(np.abs(ai_m - r_x)).argmin()]
            di_x.update({r_x: i_m})
        df_seed['m'] = [di_x[r_x] for r_x in df_seed.loc[:,'x']]
        di_y = {}
        for r_y in sorted(set(df_seed.loc[:,'y'])):
            i_n = ai_n[(np.abs(ai_n - r_y)).argmin()]
            di_y.update({r_y: i_n})
        df_seed['n'] = [di_y[r_y] for r_y in df_seed.loc[:,'y']]
        di_z = {}
        for r_z in sorted(set(df_seed.loc[:,'z'])):
            i_p = ai_p[(np.abs(ai_p - r_z)).argmin()]
            di_z.update({r_z: i_p})
        df_seed['p'] = [di_z[r_z] for r_z in df_seed.loc[:,'z']]
        # coor
        df_seed['coor'] = df_seed.loc[:,'m'].astype(str) + '_' + df_seed.loc[:,'n'].astype(str) + '_' + df_seed.loc[:,'p'].astype(str)
        es_coor = set(self.coor.loc[:,'m'].astype(str) + '_' + self.coor.loc[:,'n'].astype(str) + '_' + self.coor.loc[:,'p'].astype(str))
        df_seed = df_seed.loc[df_seed.loc[:,'coor'].isin(es_coor),['x','y','z']]

        # agent celltypeing
        i_total = df_seed.shape[0]
        ls_type = []
        lr_fract = []
        for s_type, r_fract in agent_type_fraction.items():
            ls_type.append(s_type)
            lr_fract.append(r_fract)
        df_seed['type'] = np.random.choice(ls_type, size=i_total, p=lr_fract)

        # add to results
        if (df_agent is None):
            df_agent = df_seed
        else:
            df_agent = pd.concat([df_agent, df_seed])

        # save result
        self.agent = df_agent

    def df_mesh(self):
        """
        """
        return(self.coor)

    def df_agent(self):
        """
        """
        return(self.agent)

    def z_stack_mesh(self, show=False, plot=None, png=None, movie=False, facecolor='white', cmap='inferno', s=1, figsize=(4, 4), frame_rate=24):
        """
        plot z stack from shape df
        png None or name
        movie None or name (png)
        """
        df_coor = self.df_mesh().rename({'m':'x', 'n':'y', 'p':'z'}, axis=1)
        df_coor['type'] = 'mesh'
        pdplt.df_label_to_color(df_abc=df_coor, s_label='type', s_cmap=cmap, b_shuffle=False)
        z_stack(
            df_coor = df_coor,
            show = show,
            plot = plot,
            png = png,
            movie  = movie,
            facecolor = facecolor,
            s = s,
            figsize = figsize,
            frame_rate = frame_rate,
        )

    def z_stack_agent(self, show=False, plot=None, png=None, movie=False, facecolor='white', cmap='inferno', s=1, figsize=(4, 4), frame_rate=24):
        """
        plot z stack from agent df
        png None or name
        movi None or name (png)
        """
        # plot agents
        df_coor = self.df_agent().copy()
        df_coor.z = df_coor.z.round().astype(int)
        if type(cmap) == dict:
            df_coor['type_color'] = [cmap[s_type] for s_type in df_coor.loc[:,'type']]
        else:
            pdplt.df_label_to_color(df_abc=df_coor, s_label='type', s_cmap=cmap, b_shuffle=False)
        z_stack(
            df_coor = df_coor,
            show = show,
            plot = plot,
            png = png,
            movie  = movie,
            facecolor = facecolor,
            s = s,
            figsize = figsize,
            frame_rate = frame_rate,
        )


class Brick(Shape):
    """
    """
    def __init__(self, x, y, z=1, origin=(0,0,0), um_p_px=1):
        """
        """
        # inhert
        super(Brick, self).__init__()

        # handle input
        i_um_p_px = int(round(um_p_px))
        ti_origin = (int(round(origin[0])), int(round(origin[1])), int(round(origin[2])))

        # x keeping the zero in mind
        i_x = int(round(x / i_um_p_px))
        if (i_x % 2 != 0):
            i_x += 1
        i_rx = int((i_x * i_um_p_px) / 2)

        # y keeping the zero in mind
        i_y = int(round(y / i_um_p_px))
        if (i_y % 2 != 0):
            i_y += 1
        i_ry = int((i_y * i_um_p_px) / 2)

        # z keeping the zero in mind
        if (z > i_um_p_px):
            i_z = int(round(z / i_um_p_px))
        else:
            i_z = int(round(z))
        if (i_z % 2 != 0):
            i_z += 1
        i_rz = int((i_z * i_um_p_px) / 2)

        # compute all coordinate combinations
        iti_m = range(ti_origin[0] - i_rx, ti_origin[0] + i_rx + 1, i_um_p_px)
        iti_n = range(ti_origin[1] - i_ry, ti_origin[1] + i_ry + 1, i_um_p_px)
        if (z == 1):
            iti_p = range(1)
        else:
            iti_p = range(ti_origin[2] - i_rz, ti_origin[2] + i_rz + 1, i_um_p_px)
        df_coor = pd.DataFrame(itertools.product(iti_m, iti_n, iti_p), columns=['m','n','p'], dtype=np.int32)

        # store result
        self.um_p_px = i_um_p_px
        self.coor = df_coor


class Sphere(Shape):
    """
    """
    def __init__(self, d, z=1, origin=(0,0,0), um_p_px=1):
        """
        """
        # inhert
        super(Sphere, self).__init__()

        # compute brick coordinates
        if (z is None):
            i_z = d
        else:
            i_z = z
        o_brick = Brick(x=d, y=d, z=i_z, origin=origin, um_p_px=um_p_px)
        df_coor = o_brick.coor

        # comput r keeping the zero in mind
        i_d = int(round(d))
        if (i_d % 2 != 0):
            i_d += 1
        i_r = int(i_d / 2)

        # handle origin input
        ti_origin = (int(round(origin[0])), int(round(origin[1])), int(round(origin[2])))

        # carve with simple pythagoras sphere out of the bick.
        df_coor['r2'] = i_r**2
        df_coor['m2'] = (df_coor.loc[:,'m'] - ti_origin[0])**2
        df_coor['n2'] = (df_coor.loc[:,'n'] - ti_origin[1])**2
        df_coor['p2'] = (df_coor.loc[:,'p'] - ti_origin[2])**2
        df_coor['m2n2p2'] = df_coor.loc[:,['m2','n2','p2']].sum(axis=1)
        df_coor['inside'] = df_coor.loc[:,'m2n2p2'] <= df_coor.loc[:,'r2']
        df_coor = df_coor.loc[df_coor.inside, ['m','n','p']]

        # store result
        self.um_p_px = o_brick.um_p_px
        self.coor = df_coor

