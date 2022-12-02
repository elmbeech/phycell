#####
# title: test_seeding.py
#
# language: python3
# author: bue
# date: 2022-11-30
# license: BSD 3-Clause
#
# description:
#   pytest unit test library for shape.py
#   + https://docs.pytest.org/
#
#   note:
#   assert actual == expected, message
#   == value equality
#   is reference equality
#   pytest.approx for real values
#####

# load library
from pcSeed import shape

# physicell seed function
class TestPcSeedShapeAxis(object):
    ''' test for shape axis '''

    def test_axis_zero_balance(self):
        er_axis = shape.axis(r_min=-2, r_max=+2, r_step=1.1, r_origin=0)
        assert er_axis == {-1.1, 0.0, 1.1}

    def test_axis_zero_negative(self):
        er_axis = shape.axis(r_min=-2, r_max=0, r_step=1, r_origin=0)
        assert er_axis == {-2, -1, 0}

    def test_axis_negative_negative(self):
        er_axis = shape.axis(r_min=-4, r_max=-2, r_step=1, r_origin=0)
        assert er_axis == {-4, -3, -2}

    def test_axis_zero_positive(self):
        er_axis = shape.axis(r_min=0, r_max=+2, r_step=1, r_origin=0)
        assert er_axis == {0, 1, 2}

    def test_axis_positive_positive(self):
        er_axis = shape.axis(r_min=+2, r_max=+4, r_step=1, r_origin=0)
        assert er_axis == {2, 3, 4}

    def test_axis_zero_balance_hugestep(self):
        er_axis = shape.axis(r_min=-2, r_max=+2, r_step=3, r_origin=0)
        assert er_axis == {0}


# physicell seed Shape class
class TestPcSeedShapeShape(object):
    ''' test for shape Shape '''
    skeleton = shape.Shape()

    def test_shape_init(self, mcds=mcds):
        assert (skeleton.agent is None) and \
               (skeleton.coor is None) and \
               (skeleton.um_p_px is None)


# physicell seed Brick class
class TestPcSeedShapeBrick2D(object):
    brick = shape.Brick(x=2, y=2, z=1, origin=(0, 0, 0), um_p_px=1)

    # init
    def test_brick_init(self, brick=brick):
        assert (brick.agent is None) and \
               (brick.um_p_px == 1)
               (brick.coor.shape == (9, 3)) and \
               (all(brick.coor.columns == ['m', 'n', 'p']) and \
               (str(type(brick.coor)) == "<class 'pandas.core.frame.DataFrame'>")

    # seeding
    def test_brick_seeding_hexagonal(self, brick=brick):
        brick.seeding_hexagonal({'a':1/3, 'b':2/3}, agent_diameter_um=1.0, lattice='HPC')
        # bue 2022-12-01: code to keep you sane while developing!
        #brick.agent.plot(kind='scatter', x='x', y='y', grid=True)
        assert (brick.agent.shape == (23, 4)) and \
               (all(brick.agent.columns == ['x', 'y', 'z', 'type']) and \
               (str(type(brick.agent) == "<class 'pandas.core.frame.DataFrame'>") and\
               (set(round(brick.agent.x).astype(int)) == {-1, 0, 1}) and \
               (set(round(brick.agent.y).astype(int)) == {-1, 0, 1}) and \
               (set(round(brick.agent.z).astype(int)) == {0}) and \
               (set(brick.agent.type) == {'a','b'})


class TestPcSeedShapeBrick3D(object):
    brick = shape.Brick(x=8, y=8, z=8, origin=(0, 0, 0), um_p_px=1)

    # init
    def test_brick_init(self):
        assert (brick.agent is None) and \
               (brick.um_p_px == 1)
               (brick.coor.shape == (729, 3)) and \
               (all(brick.coor.columns == ['m', 'n', 'p']) and \
               (str(type(brick.coor)) == "<class 'pandas.core.frame.DataFrame'>")

    # seeding
    def test_brick_seeding_hexagonal_hpc(self, brick=brick):
        brick.seeding_hexagonal({'a':1/3, 'b':2/3}, agent_diameter_um=2.0, lattice='HPC')
        # bue 2022-12-01: code to keep you sane while developing!
        #import matplotlib.pyplot as plt
        #%matplotlib
        #fig, ax = plt.subplots(nrows=1, ncols=3)
        #li_ax = [0,0,1,1,2]
        #ls_color = ['purple','navy','lime','red','maroon']
        #for i, r_z in enumerate(sorted(set(brick.agent.z), reverse=True)):
        #    brick.agent.loc[brick.agent.z == r_z,:].plot(kind='scatter', x='x', y='y', grid=True, ylim=(-4,4), xlim=(-4,4), c=ls_color[i], ax=ax[li_ax[i]])
        assert (brick.agent.shape == (260, 4)) and \
               (all(brick.agent.columns == ['x', 'y', 'z', 'type']) and \
               (str(type(brick.agent) == "<class 'pandas.core.frame.DataFrame'>") and\
               (set(round(brick.agent.x).astype(int)) == {-1, 0, 1}) and \
               (set(round(brick.agent.y).astype(int)) == {-1, 0, 1}) and \
               (set(round(brick.agent.z).astype(int)) == {-1, 0, 1}) and \
               (set(brick.agent.type) == {'a','b'})


    def test_brick_seeding_hexagonal_fcc(self, brick=brick):
        brick.seeding_hexagonal({'a':1/3, 'b':2/3}, agent_diameter_um=2.0, lattice='FCC')
        # bue 2022-12-01: code to keep you sane while developing!
        #import matplotlib.pyplot as plt
        #%matplotlib
        #fig, ax = plt.subplots(nrows=1, ncols=2)
        #li_ax = [0,0,0,1,1]
        #ls_color = ['purple','navy','lime','red','maroon']
        #for i, r_z in enumerate(sorted(set(brick.agent.z), reverse=True)):
        #    brick.agent.loc[brick.agent.z == r_z,:].plot(kind='scatter', x='x', y='y', grid=True, ylim=(-4,4), xlim=(-4,4), c=ls_color[i], ax=ax[li_ax[i]])
        assert (brick.agent.shape == (258, 4)) and \
               (all(brick.agent.columns == ['x', 'y', 'z', 'type']) and \
               (str(type(brick.agent) == "<class 'pandas.core.frame.DataFrame'>") and\
               (set(round(brick.agent.x).astype(int)) == {-1, 0, 1}) and \
               (set(round(brick.agent.y).astype(int)) == {-1, 0, 1}) and \
               (set(round(brick.agent.z).astype(int)) == {-1, 0, 1}) and \
               (set(brick.agent.type) == {'a','b'})


# physicell seed Brick class
#class TestPcSeedShapeShpere(object):
#    sphere = shape.Sphere()

