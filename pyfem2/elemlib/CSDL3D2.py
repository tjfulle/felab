from ._csd_ndlink import ND2NodeLinkElement
class ElasticLink3D2(ND2NodeLinkElement):
    dimensions = 3
    signature = [(1,1,1,0,0,0,0),
                 (1,1,1,0,0,0,0)]  # 2 NODE 3D LINE LINK
