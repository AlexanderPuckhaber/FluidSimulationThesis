

class MixtureParticle(object):

    def __init__(self, materialList=None, baseProperties=None):

        if materialList is None:
            self.materials = {}
        else:
            self.materials = materialList

        if baseProperties is None:
            self.baseProperties = {}
        else:
            self.baseProperties = baseProperties
    
    def addMaterial(self, newMaterial=str, materialProperties=None):
        if materialProperties is None:
            self.materials[newMaterial] = {}
        else:
            self.materials[newMaterial] = materialProperties
    
    def addMaterialProperties(self, materialString=str, materialProperties={}):
        self.materials[materialString].update(materialProperties)

    def removeMaterial(self, materialString=str):
        self.materialList.pop(materialString, None)
    
    def getMaterialList(self):
        return self.materials

    def getBaseProperties(self):
        return self.baseProperties

    def getMaterialProperties(self, materialString=str):
        return self.materialList[materialString]

    def generateFullParticleProperties(self):
        self.fullParticleProperties = {}
        # first, add base properties
        self.fullParticleProperties.update(self.baseProperties)

        # next, we need a float for the volume fraction of each material
        # we also need to flatten the dicts containing material properties
        # so 'argon':{rho: .5, m:4} --> {argon_rho: .5, argon_m:}
        for material in self.materials:
            materialName = material
            # insert volume_fraction variable
            # self.fullParticleProperties.update({'frac_' + materialName: {self.materials[material].get('_frac')}})
            tmp = {}
            for materialProperty in self.materials[material]:
                tmp.update({materialName + '_' + materialProperty: self.materials[material].get(materialProperty)})
            self.fullParticleProperties.update(tmp)
        
        return self.fullParticleProperties


        ### TODO: equations to handle the base properties (densities, masses, volumes, etc.) of materials and
        ### multiply them by their volume fractions to get the value for the whole particle
        ### So basically equations of state

        ### also, I need to use/make equations of state for different materials, such as how they respond to heat etc.
        ### and I probably 
    

