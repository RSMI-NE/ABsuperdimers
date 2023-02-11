from sqlite3 import enable_shared_cache
from rsmine.coarsegrainer.cg_utils import *
from rsmine.coarsegrainer.analysis_utils import *
import polygons

examples_dir = '.'  # '/home/cluster/mkochj/GitHub/RSMI-NE/examples'
regions_data_dir = examples_dir + '/quasiperiodic_data/regions/'

############# Load the definition of the graph and construct the reference graph

edges = np.loadtxt(os.path.join(
    examples_dir, 'quasiperiodic_data', 'edgedata.dat'))
nodes = np.loadtxt(os.path.join(
    examples_dir, 'quasiperiodic_data', 'nodedata.dat'), usecols=[0])
nodepos = np.loadtxt(os.path.join(
    examples_dir, 'quasiperiodic_data', 'nodedata.dat'), usecols=[1, 2])
order_8 = np.loadtxt(os.path.join(
    examples_dir, 'quasiperiodic_data', 'order_of_8vtx.dat'))

G = construct_reference_graph(edges, nodes)


def rel_angle(vec, factor, vertex_type=None, ref_vec=None):
    """Returns the angle between vec and the reference vector.
    Used for rotating scaled elementary polygons.

    vertex_type is needed to determine the reference vector from the list
    """

    if ref_vec is None and vertex_type is not None:
        ref_vec = factor * nodepos[polygons.super_vertices[vertex_type][1][0]]

    dir = 'cw'
    rot = np.cross(vec, ref_vec)
    vdir = 'cw' if rot > 0 else 'ccw'

    r = (vec[0]*ref_vec[0]+vec[1]*ref_vec[1])/(np.linalg.norm(vec)*np.linalg.norm(ref_vec))

    angle = np.arccos(r)

    if vdir != dir:
        angle = 2*np.pi - angle

    #angle = np.arccos(np.dot(vec, ref_vec) /
    #                  (np.linalg.norm(vec)*np.linalg.norm(ref_vec)))

    if np.isnan(angle):
        return 0
    else:
        return angle


def prepare_data(V_index: int, samples: list = [10705, 12121], examples_dir='.', regions_data_dir=None,):

    if regions_data_dir is None:
        regions_data_dir = examples_dir + '/quasiperiodic_data/regions/'

    #V_order = order_8[V_index][1]
    V_order = order_8[order_8[:, 0] == V_index][0, 1].astype(int)

    # concatenated sample string. An auxiliary variable.
    sample_no = int("".join([str(i) for i in samples]))
    disc = 1

    #We store these in a dictionary for experiment logging purposes only
    EV_params = {
        'env_size': None,  # Not used
        'buffer_size': None,  # Not used
        'block_size': None,  # Not used
        'sample_no': sample_no,
        'sample_seeds': samples,
        'V_index': V_index,
        'V_order': V_order,
    }

    # These are needed by the generator
    data_params = {
        'model': 'dimer_graph',
        'lattice_type': 'networkx',
        'L': 2048,  # unused currently
        # 1 is default. Set to actual value after generating or loading data.
        'N_samples': 1,
        'dimension': 2,
    }

    if V_order != -1:
        mem_edges0_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + 'reg_mem_edgelist_svtx%i_sord%i_mem0.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
    elif V_order == -1:
        mem_edges0_ids = []
        mem_edges1_ids = []
        mem_edges2_ids = []
        mem_edges3_ids = []
        mem_edges4_ids = []

    if V_order == 0:
        mem_edges1_ids = []
        mem_edges2_ids = []
        mem_edges3_ids = []
        mem_edges4_ids = []
    elif V_order == 1:
        mem_edges1_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + 'reg_mem_edgelist_svtx%i_sord%i_mem1.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
        mem_edges2_ids = []
        mem_edges3_ids = []
        mem_edges4_ids = []
    elif V_order == 2:
        mem_edges1_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + 'reg_mem_edgelist_svtx%i_sord%i_mem1.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
        mem_edges2_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + 'reg_mem_edgelist_svtx%i_sord%i_mem2.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
        mem_edges3_ids = []
        mem_edges4_ids = []
    elif V_order > 2:
        mem_edges1_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + 'reg_mem_edgelist_svtx%i_sord%i_mem1.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
        mem_edges2_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + 'reg_mem_edgelist_svtx%i_sord%i_mem2.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
        mem_edges3_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + 'reg_mem_edgelist_svtx%i_sord%i_mem3.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
        mem_edges4_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + 'reg_mem_edgelist_svtx%i_sord%i_mem4.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
        

    if V_order != -1:
        reg_edges0_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + 'reg_edgelist_svtx%i_sord%i_mem0.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
    elif V_order == -1:
        reg_edges0_ids = []
        reg_edges1_ids = []
        reg_edges2_ids = []
        reg_edges3_ids = []
        reg_edges4_ids = []

    if V_order == 0:
        reg_edges1_ids = []
        reg_edges2_ids = []
        reg_edges3_ids = []
        reg_edges4_ids = []
    elif V_order == 1:
        reg_edges1_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + 'reg_edgelist_svtx%i_sord%i_mem1.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
        reg_edges2_ids = []
        reg_edges3_ids = []
        reg_edges4_ids = []
    elif V_order == 2:
        reg_edges1_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + 'reg_edgelist_svtx%i_sord%i_mem1.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
        reg_edges2_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + '/reg_edgelist_svtx%i_sord%i_mem2.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
        reg_edges3_ids = []
        reg_edges4_ids = []
    elif V_order > 2:
        reg_edges1_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + 'reg_edgelist_svtx%i_sord%i_mem1.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
        reg_edges2_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + '/reg_edgelist_svtx%i_sord%i_mem2.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
        reg_edges3_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + '/reg_edgelist_svtx%i_sord%i_mem3.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
        reg_edges4_ids = eids_from_edges(edges, np.loadtxt(
            regions_data_dir + '/reg_edgelist_svtx%i_sord%i_mem4.dat' % (EV_params['V_index'], EV_params['V_order']), dtype=int))
    
    ball4, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=4, cap=8)
    ball5, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=5, cap=8)
    ball6, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=6, cap=8)
    ball8, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=8, cap=11)
    ball10, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=10, cap=15)
    ball12, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=12, cap=15)
    ball16, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=16, cap=19)
    ball32, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=32, cap=35)
    ball40, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=40, cap=50)
    ball48, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=48, cap=50)
    ball54, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=54, cap=60)
    ball64, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=64, cap=70)
    ball68, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=68, cap=70)
    ball72, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=72, cap=84)
    ball84, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=84, cap=90)
    ball24, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=24, cap=26)
    ball20, _ = construct_VE_edgelists(G, EV_params['V_index'], L_B=1, ll=20, cap=26)


    # Construct the dictionary of cases
    #### Be super careful!!! IF the variables  reg_.... is initialized from some previous vertex,
    # there won't be an error, ven of lower order vertex doesn't have this membrane.
    cGV_edges = {}
    cGE_edges = {}

    # Case #1
    auxGV_edges = sorted(set(reg_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(reg_edges1_ids) -
                         set(reg_edges0_ids) - set(mem_edges0_ids))
    cGV_edges[1] = auxGV_edges
    cGE_edges[1] = auxGE_edges

    # Case #2
    auxGV_edges = sorted(set(reg_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(reg_edges1_ids) | set(
        mem_edges1_ids) - set(reg_edges0_ids) - set(mem_edges0_ids))
    cGV_edges[2] = auxGV_edges
    cGE_edges[2] = auxGE_edges

    # Case #3 # here the environment shell is defined by some membrane, so it doesn't work for order-0 vertices
    auxGV_edges = sorted(set(reg_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(reg_edges2_ids) -
                         set(reg_edges0_ids) - set(mem_edges0_ids))
    cGV_edges[3] = auxGV_edges
    cGE_edges[3] = auxGE_edges

    # Case #4
    auxGV_edges = sorted(set(reg_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(reg_edges2_ids) - set(reg_edges1_ids))
    cGV_edges[4] = auxGV_edges
    cGE_edges[4] = auxGE_edges

    # Case #5
    auxGV_edges = sorted(set(reg_edges0_ids) | set(mem_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(reg_edges2_ids) - set(reg_edges1_ids))
    cGV_edges[5] = auxGV_edges
    cGE_edges[5] = auxGE_edges

    # Case #6
    auxGV_edges = sorted(set(reg_edges0_ids) | set(mem_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(reg_edges2_ids) -
                         set(reg_edges1_ids) - set(mem_edges1_ids))
    cGV_edges[6] = auxGV_edges
    cGE_edges[6] = auxGE_edges

    # Case #7
    auxGV_edges = sorted(set(reg_edges0_ids) | set(mem_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball32) - set(reg_edges2_ids))
    cGV_edges[7] = auxGV_edges
    cGE_edges[7] = auxGE_edges

    # Case #8
    auxGV_edges = sorted(set(reg_edges0_ids) | set(mem_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(
        set(ball32) - set(reg_edges2_ids) - set(mem_edges2_ids))
    cGV_edges[8] = auxGV_edges
    cGE_edges[8] = auxGE_edges

    # Case #9 # this case does not include the nearest neighbours in the 1-superlattice
    auxGV_edges = sorted(set(reg_edges1_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(reg_edges2_ids) -
                         set(reg_edges1_ids) - set(mem_edges1_ids))
    cGV_edges[9] = auxGV_edges
    cGE_edges[9] = auxGE_edges

    # Case #10 # this does not have any buffer
    auxGV_edges = sorted(set(reg_edges1_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(reg_edges2_ids) | set(
        mem_edges2_ids) - set(reg_edges1_ids) - set(mem_edges1_ids))
    cGV_edges[10] = auxGV_edges
    cGE_edges[10] = auxGE_edges

    # Case #11
    auxGV_edges = sorted(set(reg_edges1_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(
        set(ball32) - set(reg_edges1_ids) - set(mem_edges1_ids))
    cGV_edges[11] = auxGV_edges
    cGE_edges[11] = auxGE_edges

    # Case #12 # the buffer is terminated exactly in the boundary of the 1-superlattice neighbour blocks
    auxGV_edges = sorted(set(reg_edges1_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(
        set(ball32) - set(reg_edges2_ids) - set(mem_edges2_ids))
    cGV_edges[12] = auxGV_edges
    cGE_edges[12] = auxGE_edges

    # Case #13
    auxGV_edges = sorted(set(reg_edges1_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball32) - set(reg_edges2_ids))
    cGV_edges[13] = auxGV_edges
    cGE_edges[13] = auxGE_edges

    # Case #14
    auxGV_edges = sorted(set(reg_edges1_ids) | set(mem_edges1_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball32) - set(reg_edges2_ids))
    cGV_edges[14] = auxGV_edges
    cGE_edges[14] = auxGE_edges

    # Case #15
    auxGV_edges = sorted(set(reg_edges1_ids) | set(mem_edges1_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(
        set(ball32) - set(reg_edges2_ids) - set(mem_edges2_ids))
    cGV_edges[15] = auxGV_edges
    cGE_edges[15] = auxGE_edges

    # Case #16
    auxGV_edges = sorted(set(reg_edges2_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(
        set(ball32) - set(reg_edges2_ids) - set(mem_edges2_ids))
    cGV_edges[16] = auxGV_edges
    cGE_edges[16] = auxGE_edges

    # Case #17
    auxGV_edges = sorted(set(reg_edges2_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(
        set(ball48) - set(reg_edges2_ids) - set(mem_edges2_ids))
    cGV_edges[17] = auxGV_edges
    cGE_edges[17] = auxGE_edges

    # Case #18 was is ball 20 or ball 24???? # this does not contain the entire nerest-neighbour blocks in the 2-superlattice
    auxGV_edges = sorted(set(reg_edges2_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball48) - set(ball24))
    cGV_edges[18] = auxGV_edges
    cGE_edges[18] = auxGE_edges

    # Case #19
    auxGV_edges = sorted(set(reg_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball16) - set(reg_edges1_ids))
    cGV_edges[19] = auxGV_edges
    cGE_edges[19] = auxGE_edges

    # Case #20
    auxGV_edges = sorted(set(reg_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(
        set(ball16) - set(reg_edges1_ids) - set(mem_edges1_ids))
    cGV_edges[20] = auxGV_edges
    cGE_edges[20] = auxGE_edges

    # Case #21
    auxGV_edges = sorted(set(reg_edges0_ids) | set(mem_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball16) - set(reg_edges1_ids))
    cGV_edges[21] = auxGV_edges
    cGE_edges[21] = auxGE_edges

    # Case #22
    auxGV_edges = sorted(set(reg_edges0_ids) | set(mem_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(
        set(ball16) - set(reg_edges1_ids) - set(mem_edges1_ids))
    cGV_edges[22] = auxGV_edges
    cGE_edges[22] = auxGE_edges

    # Case #23
    auxGV_edges = sorted(set(reg_edges1_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(
        set(ball16) - set(reg_edges1_ids) - set(mem_edges1_ids))
    cGV_edges[23] = auxGV_edges
    cGE_edges[23] = auxGE_edges

    # Case #24
    auxGV_edges = sorted(set(reg_edges1_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball24) - set(ball16))
    cGV_edges[24] = auxGV_edges
    cGE_edges[24] = auxGE_edges

    #Case #25
    auxGV_edges = sorted(set(reg_edges1_ids) | set(mem_edges1_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball24) - set(ball16))
    cGV_edges[25] = auxGV_edges
    cGE_edges[25] = auxGE_edges

    # Case #26
    auxGV_edges = sorted(set(reg_edges1_ids) | set(mem_edges1_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball16) - set(ball12))
    cGV_edges[26] = auxGV_edges
    cGE_edges[26] = auxGE_edges

    # Case #27
    auxGV_edges = sorted(set(reg_edges2_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball32) - set(ball24))
    cGV_edges[27] = auxGV_edges
    cGE_edges[27] = auxGE_edges

    # Case #28
    auxGV_edges = sorted(set(reg_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(
        set(ball8) - set(reg_edges0_ids) - set(mem_edges0_ids))
    cGV_edges[28] = auxGV_edges
    cGE_edges[28] = auxGE_edges

    # Case #29
    auxGV_edges = sorted(set(reg_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(
        set(ball12) - set(reg_edges0_ids) - set(mem_edges0_ids))
    cGV_edges[29] = auxGV_edges
    cGE_edges[29] = auxGE_edges

    #Case #30
    auxGV_edges = sorted(set(reg_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball12) - set(ball5))
    cGV_edges[30] = auxGV_edges
    cGE_edges[30] = auxGE_edges

    # Case #31
    auxGV_edges = sorted(set(reg_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball16) - set(ball5))
    cGV_edges[31] = auxGV_edges
    cGE_edges[31] = auxGE_edges

    # Case #32
    auxGV_edges = sorted(set(reg_edges0_ids) | set(mem_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball16) - set(ball5))
    cGV_edges[32] = auxGV_edges
    cGE_edges[32] = auxGE_edges


    # Case #33
    auxGV_edges = sorted(set(reg_edges2_ids))
    auxGE_edges = sorted(
        set(ball48) - set(reg_edges3_ids) - set(mem_edges3_ids))
    cGV_edges[33] = auxGV_edges
    cGE_edges[33] = auxGE_edges

    # Case #34
    auxGV_edges = sorted(set(reg_edges2_ids))
    auxGE_edges = sorted(
        set(reg_edges4_ids) - set(reg_edges3_ids) - set(mem_edges3_ids))
    cGV_edges[34] = auxGV_edges
    cGE_edges[34] = auxGE_edges


    # Case #35
    auxGV_edges = sorted(set(reg_edges1_ids))
    auxGE_edges = sorted(
        set(reg_edges3_ids) - set(reg_edges2_ids) - set(mem_edges2_ids))
    cGV_edges[35] = auxGV_edges
    cGE_edges[35] = auxGE_edges


    # Case #36
    auxGV_edges = sorted(set(reg_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball12) -
                         set(reg_edges0_ids) - set(mem_edges0_ids))
    cGV_edges[36] = auxGV_edges
    cGE_edges[36] = auxGE_edges


    # Case #37
    # buffer width 2 case for level 0
    # can use for 8-vertices at sigma^1 (doesn't correspond to any level) (indices: 389)
    auxGV_edges = sorted(set(reg_edges0_ids))
    auxGE_edges = sorted(set(ball12) - set(ball4))
    cGV_edges[37] = auxGV_edges
    cGE_edges[37] = auxGE_edges

    # Case #38
    # buffer width 2-3 case for level 1
    # used for 8-vertices at sigma^2=level-0 (indices: 389)
    auxGV_edges = sorted(set(reg_edges1_ids))
    auxGE_edges = sorted(
        set(ball32) - set(ball8))
    cGV_edges[38] = auxGV_edges
    cGE_edges[38] = auxGE_edges

    # Case #39
    # buffer width 3-4 case for level 1
    auxGV_edges = sorted(set(reg_edges1_ids))
    auxGE_edges = sorted(
        set(ball32) - set(ball12))
    cGV_edges[39] = auxGV_edges
    cGE_edges[39] = auxGE_edges

    # Case #40
    # buffer width 3-4 case for level 1
    auxGV_edges = sorted(set(reg_edges1_ids))
    auxGE_edges = sorted(
        set(ball24) - set(ball12))
    cGV_edges[40] = auxGV_edges
    cGE_edges[40] = auxGE_edges

    # Case #41 
    # buffer width 2
    # used for 3-vertices at sigma^2=level-0 (indices: 399, 404)
    auxGV_edges = sorted(set(reg_edges0_ids) | set(mem_edges0_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball12) - set(ball5))
    cGV_edges[41] = auxGV_edges
    cGE_edges[41] = auxGE_edges

    # Case #42
    # buffer width 3-4
    # used for 3-vertices at sigma^3=level-1 (indices: 429, 447)
    auxGV_edges = sorted(set(reg_edges1_ids) | set(mem_edges1_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball32) - set(ball10))
    cGV_edges[42] = auxGV_edges
    cGE_edges[42] = auxGE_edges

    # Case #43
    # buffer width
    # used for 8-vertices at sigma^3=level-1 (indices: 389)
    auxGV_edges = sorted(set(reg_edges2_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball48) - set(ball20))
    cGV_edges[43] = auxGV_edges
    cGE_edges[43] = auxGE_edges

    # Case #44
    auxGV_edges = sorted(set(ball5))
    auxGE_edges = sorted(set(ball32) - set(ball8))
    cGV_edges[44] = auxGV_edges
    cGE_edges[44] = auxGE_edges

    # Case #45
    auxGV_edges = sorted(set(ball8))
    auxGE_edges = sorted(set(ball48) - set(ball16))
    cGV_edges[45] = auxGV_edges
    cGE_edges[45] = auxGE_edges

    # Case #46
    # used for 3-vertices at sigma^1 (doesn't correspond to any level) (indices: 393)
    auxGV_edges, auxGE_edges = construct_VE_edgelists(
        G, EV_params['V_index'], L_B=2, ll=2, cap=16)
    auxGV_edges = list(set(auxGV_edges) - set([807, 877, 879, 46284, 46289, 793]))
    cGV_edges[46] = auxGV_edges
    cGE_edges[46] = auxGE_edges

    # Case #47 
    # use for 8-vertices at sigma^4=level-2 (index:389)
    auxGV_edges = sorted(set(reg_edges3_ids))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    # BUG: environment too large: does not fit in the memory
    auxGE_edges = sorted(set(ball72) - set(ball48))
    cGV_edges[47] = auxGV_edges
    cGE_edges[47] = auxGE_edges

    # Case #48
    # use for 8-vertices at sigma^4=level-2 (index:389) (not perfect tile)
    auxGV_edges = sorted(set(ball16))
    # sorted(set(reg_edges2_ids)-(set(reg_edges0_ids) | set(mem_edges0_ids)))
    auxGE_edges = sorted(set(ball48) - set(ball24))
    cGV_edges[48] = auxGV_edges
    cGE_edges[48] = auxGE_edges


    # Case #49
    # TODO: rewrite this using the `polygons`` module as below

    # use for 8-vertices at sigma^2=level-0 (index:389) (perfect tile)
    auxGV_edges = sorted(set(reg_edges1_ids))
    auxGV_edges = list(set(auxGV_edges) -
                        set([13980, 13982,
                           26928, 26930,
                           39876, 39878,
                           7506, 7508,
                           868, 870,
                           20454, 20456,
                           876, 878,
                           33402, 33404]) -
                        set([13979, 26927, 39875, 7505, 867, 20453, 875, 33401]))
    auxGE_edges = construct_VE_edgelists(G, EV_params['V_index'], L_B=5, ll=4, cap=24)[1]
    cGV_edges[49] = auxGV_edges
    cGE_edges[49] = auxGE_edges


    # Case #50
    # use for 3-vertices at sigma^2=level-0 (index:404, 399, ...) (perfect tile)
    vertex_type = '3fold'
    sigma_scale = 2
    factor = polygons.silver_ratio ** (sigma_scale-1)
    #angle = np.arccos(
    #    np.dot(nodepos[EV_params['V_index']], factor*nodepos[polygons.super_vertices[vertex_type][1][0]]) /
    #    (np.linalg.norm(nodepos[EV_params['V_index']])*np.linalg.norm(factor*nodepos[polygons.super_vertices[vertex_type][1][0]])))
    angle = rel_angle(nodepos[EV_params['V_index']], factor, vertex_type=vertex_type)

    auxGV_edges = polygons.rotate_and_scale_Vgraph(angle, sigma_scale, vertex_type)[0]
    # TODO: can try increasing the buffer size here
    auxGE_edges = construct_VE_edgelists(G, V_index, L_B=4, ll=4, cap=20)[1]
    cGV_edges[50] = auxGV_edges
    cGE_edges[50] = auxGE_edges


    # Case 51
    # Use for 8-vertices at sigma^4=level-2 (index:389) (perfect tile)
    vertex_type = '8fold'
    sigma_scale = 4
    #factor = polygons.silver_ratio ** (sigma_scale-1)

    auxGV_edges = polygons.scale_Vgraph(sigma_scale, vertex_type)[0]
    # TODO: determine the proper buffer size here
    #auxGE_edges = construct_VE_edgelists(G, V_index, L_B=8, ll=16, cap=48)[1]
    #auxGE_edges = sorted(set(ball72) - set(ball54))
    auxGE_edges = sorted(set(ball64) - set(ball40))

    cGV_edges[51] = auxGV_edges
    cGE_edges[51] = auxGE_edges


    # Case 52
    # Use for 3-vertices at sigma^4=level-2 (index:196, 274, ...) (perfect tile)
    vertex_type = '3fold'
    sigma_scale = 4
    factor = polygons.silver_ratio ** (sigma_scale-1)
    angle = rel_angle(nodepos[EV_params['V_index']], factor, vertex_type=vertex_type)

    auxGV_edges = polygons.rotate_and_scale_Vgraph(angle, sigma_scale, vertex_type)[0]
    # TODO: determine the proper buffer size here
    #auxGE_edges = construct_VE_edgelists(G, V_index, L_B=8, ll=16, cap=48)[1]
    auxGE_edges = sorted(set(ball64) - set(ball40))
    cGV_edges[52] = auxGV_edges
    cGE_edges[52] = auxGE_edges


    # Case 53
    # Use for 8-vertices at sigma^3=level-1 (index:389) (perfect tile)
    vertex_type = '8fold'
    sigma_scale = 3
    factor = polygons.silver_ratio ** (sigma_scale-1)

    auxGV_edges = polygons.scale_Vgraph(sigma_scale, vertex_type)[0]
    auxGE_edges = sorted(set(ball48) - set(ball24))
    cGV_edges[53] = auxGV_edges
    cGE_edges[53] = auxGE_edges



    # Case 54
    # Use for 3-vertices at sigma^3=level-1 (index:196, 274, ...) (perfect tile)
    vertex_type = '3fold'
    sigma_scale = 3
    factor = polygons.silver_ratio ** (sigma_scale-1)
    angle = rel_angle(nodepos[EV_params['V_index']],
                    factor, vertex_type=vertex_type)

    auxGV_edges = polygons.rotate_and_scale_Vgraph(
        angle, sigma_scale, vertex_type)[0]
    # TODO: determine the proper buffer size here
    #auxGE_edges = construct_VE_edgelists(G, V_index, L_B=8, ll=16, cap=48)[1]
    auxGE_edges = sorted(set(ball48) - set(ball20))
    cGV_edges[54] = auxGV_edges
    cGE_edges[54] = auxGE_edges


    return EV_params, data_params, cGV_edges, cGE_edges




