import equation

class EquationMixtureParticle(object):
    def __init__(self):
        super(equation, self).__init__()
    
    def get_arrays_used_in_equation(equation):
        """Return two sets, the source and destination arrays used by the equation.
           For the mixture particle, I'm thinking that it should split the mixture particle 
           and make a virtual particle for each phase in the volume fraction. 
           If the equation specified the Oxygen phase, for example, a particle
           with .5 Oxygen and .5 Nitrogen would be split into one 
           Oxygen and another Nitrogen particle (with the correct smaller mass/volume)
           and the Oxygen would be sent to the equation. This method would
           do this for all of the applicable mixture particles...

        """
