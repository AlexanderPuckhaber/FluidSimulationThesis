

class MixtureParticle(object):

    def __init__(self, materialList=None):
        if materialList is None:
            self.materialList = []
        else:
            self.materialList = materialList
    
    def addMaterial(self, newMaterial):
        self.materialList.extend([newMaterial])

    def removeMaterial(self, materialString):
        self.materialList.remove(materialString)
    
    def getMaterialList(self):
        return self.materialList

    
    

