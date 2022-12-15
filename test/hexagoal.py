# const
#cell_radius = 8.412710547954228 # MCF7 core/Physicell_phenotype.cpp
#cell_spacing = 0.95 * 2.0 * cell_radius  # why not 1*2* cell_radius?
#tumor_radius = 250.0

import numpy as np
import pandas as pd
import itertools


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


# set domain 3[mm] * 3[mm]
i_m_min = -100
i_m_max = 100
i_n_min = -100
i_n_max = 100
i_p_min = 0
i_p_max = 1
i_layer = 2


# calculate radius and such
#agent_diameter_um = 2 * 8.412710547954228
agent_diameter_um = 2 * 10
r_d_major = (agent_diameter_um * 2)
r_3r_major = r_d_major * 1.5
r_r_major = r_d_major / 2
r_r_major_half = r_d_major / 4
r_r_minor = np.sqrt(r_r_major**2 - r_r_major_half**2)
r_d_minor = r_r_minor * 2
r_r_minor_half = r_r_minor / 2

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

# plot
%matplotlib
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
df_hexcoor_layer1a.plot(kind='scatter', x='x', y='y', grid=True, ylim=(i_n_min,i_n_max), xlim=(i_m_min,i_m_max), c='maroon', ax=ax)
df_hexcoor_layer1b.plot(kind='scatter', x='x', y='y', grid=True, ylim=(i_n_min,i_n_max), xlim=(i_m_min,i_m_max), c='red', ax=ax)
df_hexcoor_layer1c.plot(kind='scatter', x='x', y='y', grid=True, ylim=(i_n_min,i_n_max), xlim=(i_m_min,i_m_max), c='orange', ax=ax)
df_hexcoor_layer1d.plot(kind='scatter', x='x', y='y', grid=True, ylim=(i_n_min,i_n_max), xlim=(i_m_min,i_m_max), c='yellow', ax=ax)

# agent seeding coordinates
#df_seed = pd.concat([df_hexcoor_layer1a, df_hexcoor_layer1b])

# bue 2022-12-01: code to keep you sane while developing!
#import matplotlib.pyplot as plt
#%matplotlib
#fig, ax = plt.subplots(nrows=1, ncols=3)
#li_ax = [0,0,1,1,2]
#ls_color = ['purple','navy','lime','red','maroon']
#for i, r_z in enumerate(sorted(set(brick.agent.z), reverse=True)):
#    brick.agent.loc[brick.agent.z == r_z,:].plot(kind='scatter', x='x', y='y', grid=True, ylim=(-4,4), xlim=(-4,4), c=ls_color[i], ax=ax[li_ax[i]])


# mesh and hex
%matplotlib
from pcSeed import shape
import matplotlib.pyplot as plt

# 2d
brick = shape.Brick(x=8, y=8, z=1, origin=(0, 0, 0), um_p_px=2.0)
print(brick.coor.shape) # 25
brick.seeding_hexagonal({'a':1/3, 'b':2/3}, agent_diameter_um=2.0)
print(brick.agent.shape) # 50

# point
fig, ax = plt.subplots()
brick.coor.plot(kind='scatter', x='m', y='n', grid=True, ylim=(-4,4), xlim=(-4,4), c='maroon', s=80, ax=ax)
brick.agent.plot(kind='scatter', x='x', y='y', grid=True, ylim=(-4,4), xlim=(-4,4), c='orange', s=40, ax=ax)


# mesh and hex
%matplotlib
from pcSeed import shape
import matplotlib.pyplot as plt

# 3d
brick = shape.Brick(x=8, y=8, z=2, origin=(0, 0, 0), um_p_px=2.0)
print(brick.coor.shape) # 25
#brick.seeding_hexagonal({'a':1/3, 'b':2/3}, agent_diameter_um=2.0, lattice='HPC')
brick.seeding_hexagonal({'a':1/3, 'b':2/3}, agent_diameter_um=2.0, lattice='FCC')
print(brick.agent.shape) # 50
print(sorted(set(brick.coor.p)))

# point
fig, ax = plt.subplots()
brick.coor.plot(kind='scatter', x='m', y='n', grid=True, ylim=(-4,4), xlim=(-4,4), c='maroon', s=80, ax=ax)
brick.agent.loc[brick.agent.z > 0,:].plot(kind='scatter', x='x', y='y', grid=True, ylim=(-4,4), xlim=(-4,4), c='red', s=40, ax=ax)
brick.agent.loc[brick.agent.z == 0,:].plot(kind='scatter', x='x', y='y', grid=True, ylim=(-4,4), xlim=(-4,4), c='orange', s=40, ax=ax)
brick.agent.loc[brick.agent.z < 0,:].plot(kind='scatter', x='x', y='y', grid=True, ylim=(-4,4), xlim=(-4,4), c='yellow', s=40, ax=ax)

