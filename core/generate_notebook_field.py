class ItemDictGenerator:
    def __init__(self, item):
        self.item = item

    def generate_dict(self):
        return {self.item: (self.item, None)}


