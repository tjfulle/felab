from ._csd_ndlink import ND2NodeLinkElement
class ElasticLink1D2(ND2NodeLinkElement):
    dimensions = 1
    signature = [(1,0,0,0,0,0,0),  # 2 NODE 1D LINE LINK
                 (1,0,0,0,0,0,0)]
