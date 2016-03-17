from ._csd_ndlink import ND2NodeLinkElement
class ElasticLink2D2(ND2NodeLinkElement):
    dimensions = 2
    signature = [(1,1,0,0,0,0,0),  # 2 NODE 2D LINE LINK
                 (1,1,0,0,0,0,0)]
