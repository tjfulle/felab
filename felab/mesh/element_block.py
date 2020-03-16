class element_block:
    def __init__(self, name, id, labels, elefam, elecon):
        self.name = name.upper()
        self.id = id
        self.labels = labels
        self.elefam = elefam
        self.numele = len(labels)
        self.elecon = elecon

    @property
    def eletyp(self):
        return self.elefam

    @eletyp.setter
    def eletyp(self, arg):
        self.elefam = arg
