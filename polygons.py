from rsmine.coarsegrainer.cg_utils import *
from rsmine.coarsegrainer.analysis_utils import *

from AB_prepare_dataset import *

examples_dir = '../examples' 
regions_data_dir = examples_dir + '/quasiperiodic_data/regions/'

from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import shapely.plotting as splt


elementary_polygon_periphery_coordinates = {}

silver_ratio = 1 + np.sqrt(2)

super_vertices = {'8fold': {1: [389], 2:[389], 3:[389], 4:[389]}, 
                  '3fold': {1: [393], 
                            2: [404, 399, 3742, 7010, 10286, 13554, 16830, 20098], 
                            3: [429, 447], 
                            4: [274, 196]},
                    '4fold': {1: [395]},
                    '5Afold': {1: []},
                    '5Bfold': {1: [553]},
                    '6fold': {1: []},
                    '7fold': {1: []}}


def plot_region(edge_list: list, node_size: float =10, figsize: tuple=(4, 4),):
    weights = np.ones((len(edge_list), 1))
    small_V_visualize(edge_list, weights, edges=edges, nodes=nodes,
                      nodepos=nodepos, node_size=node_size, figsize=figsize)


def sorted_polygon(delV: nx.classes.graph.Graph, edge_type: str, sigma_scale: int = 1, polygon_id: int = 0):

    origin_pos = nodepos[super_vertices[edge_type][sigma_scale][polygon_id]]
    delV_pos = np.array([nodepos[int(node_id)] for node_id in delV.nodes])

    shifted_delV_pos = delV_pos - origin_pos

    sorted_delV_pos = np.array(
        sorted(shifted_delV_pos, key=lambda x: np.arctan2(x[0], x[1]))) + origin_pos

    return sorted_delV_pos


## Determining the elementary polygons

# 8_3-vertex at sigma1 (8_4 at sigma0)

V_id_8sigma1 = super_vertices['8fold'][1][0]

V_edgelist_8sigma1, E_edgelist_8sigma1 = construct_VE_edgelists(
    G, V_id_8sigma1, L_B=2, ll=2, cap=10)

V_8sigma1 = G.edge_subgraph([tuple(edges[id]) for id in V_edgelist_8sigma1])

bulkV_8sigma1 = construct_VE_edgelists(G, V_id_8sigma1, L_B=1, ll=1, cap=8)
delV_edgelist_8sigma1 = list(set(V_edgelist_8sigma1) - set(bulkV_8sigma1[0]))

delV_8sigma1 = V_8sigma1.edge_subgraph(
    [tuple(edges[id]) for id in delV_edgelist_8sigma1])

sorted_delV_pos_8sigma1 = sorted_polygon(delV_8sigma1, '8fold')

elementary_polygon_periphery_coordinates['8fold'] = sorted_delV_pos_8sigma1



# 5B-vertex at sigma1 (7 at sigma0)

V_id_5Bsigma1 = super_vertices['5Bfold'][1][0]

V_edgelist_5Bsigma1 = set(construct_VE_edgelists(G, V_id_5Bsigma1,
                                                 L_B=1, ll=2, cap=20)[0]) - set([46340, 858])
bulkV_5Bsigma1 = set(construct_VE_edgelists(
    G, V_id_5Bsigma1, L_B=1, ll=1, cap=20)[0])

delV_edgelist_5Bsigma1 = list(V_edgelist_5Bsigma1 - bulkV_5Bsigma1)


# Construct the networkx graphs for the bulk and delV regions

V_5Bsigma1 = G.edge_subgraph([tuple(edges[id]) for id in V_edgelist_5Bsigma1])
delV_5Bsigma1 = V_5Bsigma1.edge_subgraph(
    [tuple(edges[id]) for id in delV_edgelist_5Bsigma1])


sorted_delV_pos_5Bsigma1 = sorted_polygon(delV_5Bsigma1, '5Bfold')

elementary_polygon_periphery_coordinates['5Bfold'] = sorted_delV_pos_5Bsigma1



# 4-vertex at sigma1 (6 at sigma0)

V_id_4sigma1 = super_vertices['4fold'][1][0]

V_edgelist_4sigma1 = set(construct_VE_edgelists(G, V_id_4sigma1,
                                                L_B=1, ll=2, cap=20)[0]) - set([865, 873, 817, 829])
bulkV_4sigma1 = set(construct_VE_edgelists(
    G, V_id_4sigma1, L_B=1, ll=1, cap=20)[0])

delV_edgelist_4sigma1 = list(V_edgelist_4sigma1 - bulkV_4sigma1)

# Construct the networkx graphs for the bulk and delV regions

V_4sigma1 = G.edge_subgraph([tuple(edges[id]) for id in V_edgelist_4sigma1])
delV_4sigma1 = V_4sigma1.edge_subgraph(
    [tuple(edges[id]) for id in delV_edgelist_4sigma1])

sorted_delV_pos_4sigma1 = sorted_polygon(delV_4sigma1, '4fold')

elementary_polygon_periphery_coordinates['4fold'] = sorted_delV_pos_4sigma1



# 3-vertex at sigma1 (5A at sigma0)

V_id_3sigma1 = super_vertices['3fold'][1][0]
V_id_3sigma2 = super_vertices['3fold'][2][0]


V_edgelist_3sigma1, _ = construct_VE_edgelists(G, V_id_3sigma1, L_B=2, ll=2, cap=9)

V_edgelist_cleaned_3sigma1 = list(set(V_edgelist_3sigma1) - set([807, 877, 879, 46284, 46289, 793]))

V_3sigma1 = G.edge_subgraph([tuple(edges[id])
                            for id in V_edgelist_cleaned_3sigma1])

bulkV_3sigma1 = construct_VE_edgelists(G, V_id_3sigma1, L_B=1, ll=1, cap=8)
delV_edgelist_3sigma1 = list(
    set(V_edgelist_cleaned_3sigma1) - set(bulkV_3sigma1[0]))


# Construct the networkx graph for delV
delV_3sigma1 = V_3sigma1.edge_subgraph(
    [tuple(edges[id]) for id in delV_edgelist_3sigma1])

sorted_delV_pos_3sigma1 = sorted_polygon(delV_3sigma1, '3fold')

elementary_polygon_periphery_coordinates['3fold'] = sorted_delV_pos_3sigma1




def scale_Vgraph(sigma_scale: int, vertex_type: str, sorted_delV_pos=None, 
                delta: float=1+np.sqrt(2), eps: float=1e-5, central_pos_sigma1=None,
                plot: bool=False, figsize=None):
    """
    Scale the coarse-grained graph block according to the AB inflation rule.
    """

    factor = delta ** (sigma_scale-1) + eps

    #central_pos_sigma1 = nodepos[super_vertices[vertex_type][1][0]]
    if central_pos_sigma1 is None:
        central_pos_sigma1 = nodepos[super_vertices[vertex_type][1][0]]

    central_pos_scaled = delta ** (sigma_scale-1) * central_pos_sigma1
    #central_pos_sigma2 = nodepos[super_vertices[vertex_type][sigma_scale][0]]
    

    if sorted_delV_pos is None:
        sorted_delV_pos = elementary_polygon_periphery_coordinates[vertex_type]

    # scale the graph boundary
    polygon_sigma1 = Polygon(sorted_delV_pos)
    polygon_scaled = Polygon(
        factor * (sorted_delV_pos-central_pos_sigma1) + central_pos_scaled)
    
    # scale the graph
    node_filter = np.array([polygon_scaled.contains(Point(pos)) for pos in nodepos])
    pos_scaled = nodepos[node_filter]
    nodes_scaled = nodes[node_filter]
    pos_dict = dict(zip(nodes_scaled, pos_scaled))

    V_scaled = nx.subgraph(G, nodes_scaled)

    edgelist_scaled = construct_edgelist_from_nodes(G, nodes_scaled)

    if plot:
        if figsize is None:
            figsize = (factor, factor)

        fig, ax = plt.subplots(1, figsize=figsize)
        ax.set_aspect('equal')
        splt.plot_polygon(polygon_sigma1, ax=ax, color='b')
        splt.plot_polygon(polygon_scaled, ax=ax, color='r')
        ax.scatter(pos_scaled[:, 0], pos_scaled[:, 1], s=10, c='r')

    return edgelist_scaled, V_scaled, pos_dict, [polygon_sigma1, polygon_scaled]


def rotate_and_scale_Vgraph(angle:float, sigma_scale: int, vertex_type: str,):

    rot = np.array([[np.cos(angle), -np.sin(angle)],
               [np.sin(angle), np.cos(angle)]])

    sorted_delV_pos = elementary_polygon_periphery_coordinates[vertex_type]

    central_pos_sigma1_rot = np.dot(
        nodepos[super_vertices[vertex_type][1][0]], rot)

    sorted_delV_pos_rot = np.dot(sorted_delV_pos, rot) 

    edgelist_scaled_rot, V_scaled_rot, pos_dict_rot, polygons_rot = scale_Vgraph(sigma_scale, vertex_type, 
                                        sorted_delV_pos=sorted_delV_pos_rot,
                                        central_pos_sigma1=central_pos_sigma1_rot)

    return edgelist_scaled_rot, V_scaled_rot, pos_dict_rot, polygons_rot
    