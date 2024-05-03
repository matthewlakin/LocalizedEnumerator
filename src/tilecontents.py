
class TileContent:

    def __init__(self, list_of_strands):
        # List of strand objects
        self.tilecontents = list_of_strands
        self.strands = list_of_strands
        self.tethered_strands = [strands for strands in list_of_strands]

