from abc import ABC, abstractmethod


class PublicFunction(ABC):

    @abstractmethod
    def energy(self, filename, nstates: int):
        pass
        """

        Args:
            filename: The output file of quantum chemistry software: "gauss.log", "orca.out" ..:
            nstates: int , the current states
        Returns:
            The excited energy list [S0, S1, S2, S3] and  the current total energy(float)
        """

    @abstractmethod
    def grad(self, filename):
        pass
        """

        Args:
            filename: The output file of quantum chemistry software: "gauss.log", "orca.out" ...

        Returns:
           The grad matrix(3* Natom)
        """

    @abstractmethod
    def cico(self, time, filename):
        pass
        """

        Args:
            time: The current dynamics time(default float) or The run time(int)
            filename:  The output file of quantum chemistry software: "gauss.log", "orca.out" ...

        Returns: None
        """

    @abstractmethod
    def check(self, nstates):
        pass
        """

        Args:
            nstates: The current nstates == Inputfile root??

        Returns: None 

        """

    @abstractmethod
    def r_wavefunction(self, filename):
        pass

    @abstractmethod
    def d_wavefunction(self, filename):
        pass

    @abstractmethod
    def renew_calc_states(self, nstates, filename, filename_new=None, **kwargs):
        pass

    # @abstractmethod
    # def keyword(self, filename):
    #     pass